# =============================================================================
# SpaCoDA niche identification -- the Aitchison/PCoA pipeline used in the paper
# =============================================================================
# This file is the niche-finding method of the accompanying manuscript, ported
# into the package so that the package and the paper use one implementation.
#
# It is a faithful port of the study pipeline
#   0.function_auto_select_resoultion_seurat_20260817.r (study code)
# as called by the BRCA and GBM analyses
#   nicheDecon_brca_20260816_01_fixed.r
#   nicheDecon_gbm_20260816_01_fixed.r
# (those two call sites are byte-identical). The numeric behaviour, the
# degenerate-input handling and the returned fields are preserved; only the
# language of messages and the argument defaults differ (see below).
#
# Pipeline, per section:
#   1. Aitchison distance on the cell-type composition of the spots
#      (pseudocount for zeros, never drops spots),
#   2. spatial distance (plain Euclidean by default; Gaussian kernel optional),
#   3. fuse the two distance matrices (default weights 0.7 / 0.3),
#   4. PCoA (cmdscale) on the fused distance -> shared-nearest-neighbour graph,
#   5. scan clustering resolutions and keep the one with maximum modularity,
#   6. niche prototypes = mean composition of each niche.
#
# Defaults follow the PAPER configuration (BRCA/GBM call site), which differs
# from the original helper's own signature in two places:
#   use_rank_norm      paper FALSE      (helper default TRUE)
#   resolution_seq     paper 0.3 .. 3   (helper default 0.3 .. 1.5)
# Both are arguments, so either behaviour is one call away.
# =============================================================================

#' Count cell types that actually occur (column sums > 0)
#' @keywords internal
n_effective_celltypes <- function(prop_mat) {
  sum(colSums(as.matrix(prop_mat)) > 0)
}

#' Aitchison distance between spot compositions
#'
#' Zeros are handled with a multiplicative pseudocount. Spots are never dropped:
#' `zCompositions::cmultRepl` is deliberately not used because it removes rows
#' and columns with more than ~80% zeros, which pure-tumour spots reach in
#' multi-cell-type data, and that silently changes the spot set.
#'
#' @param prop_mat spots x cell types, non-negative proportions (rows need not sum to 1)
#' @param pseudo pseudocount added to every entry before the CLR transform
#' @return n x n distance matrix. If fewer than two cell types carry any signal the
#'   compositions are all identical in the simplex, and an all-zero distance
#'   matrix is returned (a legal input: a homogeneous section).
#' @export
compute_aitchison_distance <- function(prop_mat, pseudo = 1e-6) {
  P <- as.matrix(prop_mat)
  spot_names <- rownames(prop_mat)

  all_zero_cols <- colSums(P) == 0
  if (any(all_zero_cols)) {
    message(sprintf("    dropping %d all-zero cell type(s)", sum(all_zero_cols)))
    P <- P[, !all_zero_cols, drop = FALSE]
  }

  if (ncol(P) < 2) {
    message(sprintf("    degenerate composition: %d effective cell type(s), distances are 0",
                    ncol(P)))
    return(matrix(0, nrow(P), nrow(P), dimnames = list(spot_names, spot_names)))
  }

  P <- P + pseudo
  P <- P / rowSums(P)
  L <- log(P)
  clr <- L - rowMeans(L)

  dist_ait <- as.matrix(dist(clr))
  rownames(dist_ait) <- colnames(dist_ait) <- spot_names
  dist_ait
}

#' Spatial distance with an adaptive-bandwidth Gaussian kernel
#' @keywords internal
compute_spatial_kernel_distance <- function(coords_mat, k_nn = 10) {
  dist_spa <- as.matrix(dist(coords_mat))
  k_use <- min(k_nn, nrow(coords_mat) - 1)
  if (k_use < 1) stop("too few spots to compute a spatial kernel distance")
  nn <- FNN::get.knn(coords_mat, k = k_use)
  sigma <- mean(nn$nn.dist[, k_use])
  if (!is.finite(sigma) || sigma <= 0) sigma <- max(dist_spa) / 10
  1 - exp(-(dist_spa^2) / (2 * sigma^2))
}

#' Plain Euclidean spatial distance
#' @keywords internal
compute_spatial_distance <- function(coords_mat) {
  as.matrix(dist(coords_mat))
}

#' Rank-normalise a distance matrix to (0, 1]
#' @keywords internal
rank_normalize <- function(d) {
  d_vec <- as.vector(d)
  r <- rank(d_vec, ties.method = "average") / length(d_vec)
  matrix(r, nrow = nrow(d), dimnames = dimnames(d))
}

#' Fuse a composition distance with a spatial distance
#'
#' Both distances are put on a comparable scale first: by rank when
#' `use_rank_norm = TRUE`, otherwise by dividing by their maximum. The fused
#' matrix is `expr_weight * composition + spat_weight * space`.
#'
#' @param dist_expr composition distance (see [compute_aitchison_distance()])
#' @param dist_spatial spatial distance
#' @param use_rank_norm rank-normalise both matrices (paper: FALSE)
#' @param gaussian_kernel_distance TRUE if `dist_spatial` is already a 0-1 kernel
#' @param expr_weight,spat_weight fusion weights (paper: 0.7 / 0.3)
#' @export
fuse_distances <- function(dist_expr, dist_spatial,
                           use_rank_norm = FALSE,
                           gaussian_kernel_distance = FALSE,
                           expr_weight = 0.7, spat_weight = 0.3) {
  if (!identical(dim(dist_expr), dim(dist_spatial))) {
    stop(sprintf("distance matrices differ in shape: composition %dx%d vs space %dx%d",
                 nrow(dist_expr), ncol(dist_expr), nrow(dist_spatial), ncol(dist_spatial)))
  }

  expr_range <- max(dist_expr)
  expr_is_flat <- expr_range < 1e-12

  if (isTRUE(use_rank_norm) && !expr_is_flat) {
    expr_norm    <- rank_normalize(dist_expr)
    spatial_norm <- rank_normalize(dist_spatial)
  } else {
    expr_norm <- if (expr_is_flat) dist_expr else dist_expr / expr_range
    if (isTRUE(gaussian_kernel_distance)) {
      spatial_norm <- dist_spatial
    } else {
      sp_max <- max(dist_spatial)
      spatial_norm <- if (sp_max < 1e-12) dist_spatial else dist_spatial / sp_max
    }
  }

  expr_weight * expr_norm + spat_weight * spatial_norm
}

#' Classical PCoA (cmdscale) restricted to the top-k eigenpairs
#'
#' `stats::cmdscale()` runs a full `eigen()` on the N x N double-centred matrix,
#' O(N^3), although only `k` (30) eigenpairs are used: it was ~70% of the niche
#' runtime at a few thousand spots. This computes the same top-k eigenpairs of
#' the same matrix, B = -1/2 J D^2 J, with a Lanczos solver (RSpectra) applied as
#' an operator, so B itself is never formed. The coordinates equal cmdscale's up
#' to the sign of each axis and solver tolerance; distances between spots -- all
#' that the SNN graph uses -- are unchanged. Like cmdscale, axes with
#' non-positive eigenvalues are dropped with a warning.
#' @keywords internal
pcoa_partial <- function(dist_fused, k) {
  n <- nrow(dist_fused)
  x <- dist_fused^2
  storage.mode(x) <- "double"
  op <- function(v, args) {
    v <- v - mean(v)                  # J v
    u <- as.vector(x %*% v)           # D^2 J v
    -0.5 * (u - mean(u))              # -1/2 J D^2 J v
  }
  e <- RSpectra::eigs_sym(op, k = k, which = "LA", n = n,
                          opts = list(tol = 1e-12, maxitr = 10000, retvec = TRUE))
  if (e$nconv < k) stop(sprintf("only %d of %d eigenpairs converged", e$nconv, k))
  ord <- order(e$values, decreasing = TRUE)
  ev <- e$values[ord]
  evec <- e$vectors[, ord, drop = FALSE]
  if (sum(ev > 0) < k) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", sum(ev > 0), k),
            domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
  }
  evec * rep(sqrt(ev), each = n)
}

#' Top-k PCoA coordinates of a distance matrix
#'
#' `method = "full"` is `stats::cmdscale()`, exactly as in the study pipeline.
#' `"partial"` computes only the k needed eigenpairs (see [pcoa_partial()]), and
#' falls back to `"full"` if the solver does not converge. `"auto"` uses
#' `"partial"` from 1,000 spots up when RSpectra is installed.
#' @keywords internal
pcoa_top_k <- function(dist_fused, k, method = c("auto", "partial", "full")) {
  method <- match.arg(method)
  n <- nrow(dist_fused)
  has_rspectra <- requireNamespace("RSpectra", quietly = TRUE)
  if (method == "auto") method <- if (n >= 1000 && has_rspectra) "partial" else "full"
  if (method == "partial") {
    if (!has_rspectra) stop("pcoa_method = \"partial\" needs the RSpectra package")
    pts <- tryCatch(pcoa_partial(dist_fused, k),
                    error = function(e) {
                      warning("partial eigendecomposition failed (", conditionMessage(e),
                              "); falling back to cmdscale()", call. = FALSE)
                      NULL
                    })
    if (!is.null(pts)) return(pts)
  }
  stats::cmdscale(dist_fused, k = k, eig = TRUE)$points
}

#' Build a Seurat object whose reduction is the PCoA of a fused distance
#' @keywords internal
create_seurat_with_fused_distance <- function(dist_fused, prop_mat, coords_mat,
                                              pcoa_dim = 30, pcoa_method = "auto") {
  ## sparse up front: Seurat coerces a dense matrix to dgCMatrix anyway, with a warning
  obj <- Seurat::CreateSeuratObject(counts = as(t(as.matrix(prop_mat)), "CsparseMatrix"),
                                    assay = "Spatial",
                                    min.cells = 0, min.features = 0)
  obj <- Seurat::AddMetaData(obj, metadata = as.data.frame(coords_mat))

  k <- min(pcoa_dim, nrow(dist_fused) - 1)
  if (k < 2) stop("too few spots for PCoA")

  embeddings <- pcoa_top_k(dist_fused, k = k, method = pcoa_method)

  if (is.null(dim(embeddings)) || ncol(embeddings) < 2) {
    stop("PCoA returned fewer than 2 dimensions; the distance matrix is degenerate")
  }
  colnames(embeddings) <- paste0("FUSED_", seq_len(ncol(embeddings)))
  rownames(embeddings) <- rownames(dist_fused)

  obj[["fused"]] <- Seurat::CreateDimReducObject(embeddings = embeddings,
                                                 key = "FUSED_", assay = "Spatial")
  obj
}

#' Modularity of one clustering of a per-section object
#' @keywords internal
compute_modularity <- function(obj, resolution) {
  obj <- Seurat::FindClusters(obj, resolution = resolution, algorithm = 1, verbose = FALSE)
  clusters <- obj$seurat_clusters
  k <- length(unique(clusters))
  n_cells <- ncol(obj)

  if (k <= 1 || k >= n_cells) return(list(modularity = NA_real_, k = k))

  ## Take whichever SNN graph the object actually carries: the assay is "Spatial"
  ## in the study pipeline, but hard-coding that name makes every scan silently
  ## produce NA (igraph errors on a zero-vertex graph) for any other assay.
  gname <- paste0(Seurat::DefaultAssay(obj), "_snn")
  adj <- obj@graphs[[gname]]
  if (is.null(adj)) {
    cand <- grep("_snn$", names(obj@graphs), value = TRUE)
    if (length(cand) == 0) return(list(modularity = NA_real_, k = k))
    adj <- obj@graphs[[cand[1]]]
  }
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
  list(modularity = igraph::modularity(g, as.numeric(clusters), weights = igraph::E(g)$weight),
       k = k)
}

#' Scan clustering resolutions and keep the one with the highest modularity
#'
#' This is the resolution criterion used in the paper (maximum modularity over
#' the scan). It deliberately has **no** null test; [select_resolution()] is the
#' variant that additionally compares against degree-preserving rewired graphs
#' and is used by the CLR path of [integrate_slices()].
#'
#' @param obj a Seurat object carrying an SNN graph (built by
#'   [cluster_single_sample()] or by `FindNeighbors()`)
#' @param resolution_seq candidate resolutions (paper: `seq(0.3, 3, by = 0.1)`)
#' @return list with `best_resolution`, `best_k`, `best_modularity`, `scan_results`
#' @export
select_resolution_by_modularity <- function(obj, resolution_seq = seq(0.3, 3, by = 0.1)) {
  results <- data.frame(resolution = resolution_seq, k = NA_integer_, modularity = NA_real_)

  for (i in seq_along(resolution_seq)) {
    out <- tryCatch(compute_modularity(obj, resolution_seq[i]),
                    error = function(e) list(modularity = NA_real_, k = NA_integer_))
    results$modularity[i] <- out$modularity
    results$k[i]          <- out$k
  }

  valid_results <- results[!is.na(results$modularity), ]
  if (nrow(valid_results) == 0) {
    stop("no resolution produced a valid clustering; widen resolution_seq or check the input")
  }

  best_idx <- which.max(valid_results$modularity)
  list(best_resolution = valid_results$resolution[best_idx],
       best_k          = valid_results$k[best_idx],
       best_modularity = valid_results$modularity[best_idx],
       scan_results    = valid_results)
}

#' A section that is one single niche
#'
#' Used for degenerate input (a homogeneous section, e.g. pure tumour). This is
#' a legal niche state rather than a failure.
#' @keywords internal
make_single_niche_result <- function(prop_mat, coords_mat, reason) {
  n_spot <- nrow(prop_mat)
  spot_clusters <- factor(rep("0", n_spot), levels = "0")
  names(spot_clusters) <- rownames(prop_mat)

  prototypes <- as.data.frame(t(colMeans(as.matrix(prop_mat))))
  rownames(prototypes) <- "0"

  list(spot_clusters = spot_clusters, prototypes = prototypes, seurat_obj = NULL,
       best_resolution = NA_real_, best_k = 1L, best_modularity = NA_real_,
       scan_results = NULL, dist_fused = NULL, n_spots = n_spot,
       homogeneous = TRUE, homogeneous_reason = reason)
}

#' Find the niches of one spatial section
#'
#' The paper's per-section niche clustering: Aitchison distance on the
#' composition, fused with spatial distance, PCoA, shared-nearest-neighbour
#' graph, then the clustering resolution with maximum modularity.
#'
#' A section whose composition carries fewer than `min_effective_ct` cell types,
#' or that has fewer than 10 spots, or whose fused distance is entirely flat, is
#' returned as a **single niche** instead of failing -- a pure-tumour section is
#' a niche state, not an error.
#'
#' @param prop_mat spots x cell types composition (e.g. `t(W_cell)` from
#'   [deconvolve_spots()]); rows must be named
#' @param coords_mat spots x 2 spatial coordinates, same row order as `prop_mat`
#' @param expr_weight,spat_weight weights of the composition and spatial distance
#'   when fused (paper: 0.7 / 0.3)
#' @param k_nn_spatial neighbours used for the adaptive Gaussian bandwidth
#' @param gaussian_kernel_distance use a Gaussian kernel spatial distance instead
#'   of plain Euclidean distance (paper: FALSE)
#' @param use_rank_norm rank-normalise both distances before fusing
#'   (paper: FALSE, the original helper defaulted to TRUE)
#' @param k_nn_snn neighbours for the shared-nearest-neighbour graph (paper: 20)
#' @param resolution_seq resolutions scanned when picking the number of niches
#'   (paper: `seq(0.3, 3, by = 0.1)`)
#' @param pcoa_dim PCoA dimensions requested (paper: 30)
#' @param seed random seed (paper: 42)
#' @param min_effective_ct sections with fewer effective cell types are returned
#'   as a single niche
#' @param verbose print progress
#' @param pcoa_method `"auto"` (default), `"partial"` or `"full"`; see
#'   [pcoa_top_k()]. `"full"` is `cmdscale()` as in the study; `"partial"` gives
#'   the same niches faster. `"auto"` = `"partial"` from 1,000 spots.
#' @param spatial_test run [niche_spatial_test()] on the composition (default
#'   `TRUE`); it runs after the clustering and does not change it
#' @param n_perm permutations for the spatial test
#' @return list with `spot_clusters` (one niche label per spot), `prototypes`
#'   (niches x cell types, the niche signature), `seurat_obj`, `best_resolution`,
#'   `best_k`, `best_modularity`, `scan_results`, `dist_fused`, `n_spots`,
#'   `homogeneous`, `spatial_test`
#' @export
cluster_single_sample <- function(prop_mat, coords_mat,
                                  expr_weight = 0.7, spat_weight = 0.3,
                                  k_nn_spatial = 10,
                                  gaussian_kernel_distance = FALSE,
                                  use_rank_norm = FALSE,
                                  k_nn_snn = 20,
                                  resolution_seq = seq(0.3, 3, by = 0.1),
                                  pcoa_dim = 30,
                                  seed = 42,
                                  min_effective_ct = 2,
                                  verbose = TRUE,
                                  pcoa_method = c("auto", "partial", "full"),
                                  spatial_test = TRUE, n_perm = 999) {
  pcoa_method <- match.arg(pcoa_method)
  set.seed(seed)
  say <- function(...) if (isTRUE(verbose)) cat(...)
  ## T1 runs last, on its own RNG stream, so the clustering is exactly what it
  ## would be without it
  with_test <- function(res) {
    if (!isTRUE(spatial_test)) { res["spatial_test"] <- list(NULL); return(res) }
    st <- niche_spatial_test(prop_mat, coords_mat, n_perm = n_perm, seed = seed)
    res$spatial_test <- st
    if (is.finite(st$p_value)) {
      say(sprintf("  spatial organisation of the composition: ratio = %.2f, p = %.3g\n",
                  st$ratio, st$p_value))
      if (st$p_value > 0.05)
        say("  [!] the composition is not spatially organised (p > 0.05): spatially\n",
            "      contiguous niches here come from the spatial distance term, not from\n",
            "      the composition\n")
    }
    res
  }

  prop_mat   <- as.data.frame(prop_mat)
  coords_mat <- as.matrix(coords_mat)
  stopifnot(all(rownames(prop_mat) == rownames(coords_mat)))
  n_spot_in <- nrow(prop_mat)

  ## ---- 1. degenerate composition: the whole section is one niche -----------
  n_eff_ct <- n_effective_celltypes(prop_mat)
  if (n_eff_ct < min_effective_ct) {
    dominant <- colnames(prop_mat)[which.max(colMeans(as.matrix(prop_mat)))]
    msg <- sprintf("homogeneous composition: only %d effective cell type(s) (dominant = %s)",
                   n_eff_ct, dominant)
    say(sprintf("  [!] %s\n", msg),
        "      -> the section is recorded as a single niche (a legal state, not a failure)\n")
    return(with_test(make_single_niche_result(prop_mat, coords_mat, reason = msg)))
  }

  ## ---- 2. too few spots: also a single niche -------------------------------
  if (n_spot_in < 10) {
    msg <- sprintf("too few spots (%d) for reliable clustering", n_spot_in)
    say(sprintf("  [!] %s -> recorded as a single niche\n", msg))
    return(with_test(make_single_niche_result(prop_mat, coords_mat, reason = msg)))
  }

  ## ---- 3. Aitchison distance ----------------------------------------------
  say("  [1/5] Aitchison distance ...\n")
  dist_ait <- compute_aitchison_distance(prop_mat)
  stopifnot("the Aitchison distance lost spots" = nrow(dist_ait) == n_spot_in)

  ## ---- 4. spatial distance ------------------------------------------------
  if (isTRUE(gaussian_kernel_distance)) {
    say("  [2/5] spatial Gaussian-kernel distance ...\n")
    dist_spatial <- compute_spatial_kernel_distance(coords_mat, k_nn = k_nn_spatial)
  } else {
    say("  [2/5] spatial distance ...\n")
    dist_spatial <- compute_spatial_distance(coords_mat)
  }
  stopifnot("spatial distance has the wrong shape" = nrow(dist_spatial) == n_spot_in)

  ## ---- 5. fuse ------------------------------------------------------------
  say("  [3/5] fusing distance matrices ...\n")
  dist_fused <- fuse_distances(dist_ait, dist_spatial,
                               use_rank_norm = use_rank_norm,
                               gaussian_kernel_distance = gaussian_kernel_distance,
                               expr_weight = expr_weight, spat_weight = spat_weight)

  if (max(dist_fused) < 1e-12) {
    msg <- "fused distance is entirely 0 (no composition and no spatial difference)"
    say(sprintf("  [!] %s -> recorded as a single niche\n", msg))
    return(with_test(make_single_niche_result(prop_mat, coords_mat, reason = msg)))
  }

  ## ---- 6. PCoA + SNN graph ------------------------------------------------
  say("  [4/5] Seurat object + FindNeighbors ...\n")
  obj <- create_seurat_with_fused_distance(dist_fused, prop_mat, coords_mat,
                                           pcoa_dim = pcoa_dim, pcoa_method = pcoa_method)
  stopifnot("the Seurat object lost spots" = ncol(obj) == n_spot_in)

  actual_dims <- ncol(Seurat::Embeddings(obj, "fused"))
  if (actual_dims < pcoa_dim) {
    say(sprintf("    note: PCoA returned %d dims < requested %d, using the actual ones\n",
                actual_dims, pcoa_dim))
  }
  k_use <- min(k_nn_snn, n_spot_in - 1)
  if (k_use < k_nn_snn) {
    say(sprintf("    note: k.param lowered from %d to %d (limited by spot count)\n",
                k_nn_snn, k_use))
  }

  obj <- Seurat::FindNeighbors(obj, reduction = "fused", dims = 1:actual_dims,
                              k.param = k_use, verbose = FALSE)

  ## ---- 7. number of niches = maximum modularity ---------------------------
  say("  [5/5] scanning resolutions, keeping the maximum modularity ...\n")
  optimal <- tryCatch(select_resolution_by_modularity(obj, resolution_seq),
                      error = function(e) NULL)

  if (is.null(optimal)) {
    msg <- "every resolution gives a single cluster (no separable structure)"
    say(sprintf("  [!] %s -> recorded as a single niche\n", msg))
    return(with_test(make_single_niche_result(prop_mat, coords_mat, reason = msg)))
  }

  say(sprintf("  best resolution = %.2f (k = %d, modularity = %.4f)\n",
              optimal$best_resolution, optimal$best_k, optimal$best_modularity))

  obj <- Seurat::FindClusters(obj, resolution = optimal$best_resolution,
                              algorithm = 1, verbose = FALSE)
  spot_clusters <- obj$seurat_clusters
  names(spot_clusters) <- colnames(obj)
  stopifnot("the final clustering lost spots" = length(spot_clusters) == n_spot_in)

  cl_chr <- as.character(spot_clusters)
  prototypes <- stats::aggregate(prop_mat, by = list(cluster = cl_chr), FUN = mean)
  rownames(prototypes) <- as.character(prototypes$cluster)
  prototypes$cluster <- NULL

  with_test(list(spot_clusters = spot_clusters, prototypes = prototypes, seurat_obj = obj,
       best_resolution = optimal$best_resolution, best_k = optimal$best_k,
       best_modularity = optimal$best_modularity, scan_results = optimal$scan_results,
       dist_fused = dist_fused, n_spots = n_spot_in, homogeneous = FALSE))
}
