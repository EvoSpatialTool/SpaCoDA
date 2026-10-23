# SpaCoDA: cell-type deconvolution of spatial spots
# ---------------------------------------------------
#   min_{W>=0} sum_i ||y_i - B w_i||^2_{V_i} + lambda * tr(W L W')
# (lambda = lam, 0 by default: per-spot variance-weighted NNLS)
#
# Plain R for everything except the hot loop (per-spot Gram matrices and
# NNLS solves), which is in src/coda_core.cpp.

# ---------------------------------------------------------------------------
# Signature construction
# ---------------------------------------------------------------------------
build_signature <- function(X, labels, gene_names = NULL, target = 1e4) {
  # X: (n_cells, n_genes) counts
  if (is.null(gene_names)) gene_names <- paste0("g", seq_len(ncol(X)))
  types <- sort(unique(labels))

  tot <- pmax(rowSums(X), 1)
  Xn <- X / tot * target

  B <- matrix(0, ncol(X), length(types))
  cell_sizes <- numeric(length(types))
  for (k in seq_along(types)) {
    m <- labels == types[k]
    B[, k] <- colMeans(Xn[m, , drop = FALSE])
    cell_sizes[k] <- mean(rowSums(X[m, , drop = FALSE]))
  }
  B <- sweep(B, 2, pmax(colSums(B), 1e-12), "/")
  # 0.1.1: label the signature (genes x cell types) so that downstream code and
  # user code can select by name instead of by position.
  dimnames(B) <- list(gene_names, types)

  list(B = B, cell_sizes = cell_sizes, types = types, gene_names = gene_names)
}

select_markers <- function(sig, n_per_type = 150, min_spec = 0.15,
                           min_expr_quantile = 0.5, restrict_to = NULL) {
  B <- sig$B
  spec <- B / pmax(rowSums(B), 1e-12)
  n_genes <- nrow(B); K <- ncol(B)

  eligible <- rep(TRUE, n_genes)
  if (!is.null(restrict_to)) eligible <- sig$gene_names %in% restrict_to

  selected <- integer(0)
  for (k in seq_len(K)) {
    ek <- B[, k]
    pos <- ek[ek > 0]
    thr <- if (length(pos)) quantile(pos, min_expr_quantile) else 0
    ok <- eligible & (spec[, k] >= min_spec) & (ek >= thr)
    if (!any(ok)) ok <- eligible & (ek > 0)
    score <- spec[, k] * log1p(ek * 1e4)
    score[!ok] <- -Inf
    top <- order(score, decreasing = TRUE)[seq_len(min(n_per_type, sum(is.finite(score))))]
    top <- top[is.finite(score[top])]
    selected <- union(selected, top)
  }
  sort(selected)
}

mrna_to_cell <- function(W, cell_sizes) {
  C <- W / pmax(cell_sizes, 1e-12)
  sweep(C, 2, pmax(colSums(C), 1e-12), "/")
}

# ---------------------------------------------------------------------------
# Spatial graph
# ---------------------------------------------------------------------------
# NOTE: the obvious R idiom here is `as.matrix(dist(coords))`, and it is a trap.
# It is O(N^2) in both time and memory: at 3600 spots that is a 104 MB dense
# matrix and 3.3 s -- which was 76% of the total runtime before this fix. At
# 20k spots it would be 3.2 GB and simply fall over. FNN::get.knn is a kd-tree,
# O(N log N), and is the direct equivalent of scipy's cKDTree.
build_spatial_graph <- function(coords, k = 6) {
  n <- nrow(coords)
  kn <- FNN::get.knn(coords, k = k)
  rows <- rep(seq_len(n), each = k)
  cols <- as.vector(t(kn$nn.index))
  vals <- as.vector(t(kn$nn.dist))

  sigma <- median(vals[vals > 0])
  w <- exp(-(vals^2) / (2 * sigma^2))
  A <- sparseMatrix(i = rows, j = cols, x = w, dims = c(n, n))
  A <- pmax(A, t(A))
  diag(A) <- 0
  drop0(A)
}

# ---------------------------------------------------------------------------
# Weights
# ---------------------------------------------------------------------------
weights_from_mu <- function(mu, phi = 0.1, scheme = "nb", max_ratio = 1e4) {
  mu <- pmax(mu, 1e-10)
  V <- switch(scheme,
    none    = matrix(1, nrow(mu), ncol(mu)),
    poisson = 1 / mu,
    nb      = 1 / (mu + phi * mu^2),
    dwls    = 1 / (mu^2),
    stop("unknown scheme")
  )
  if (scheme != "none") {
    vmin <- apply(V, 2, min)
    V <- pmin(V, rep(vmin, each = nrow(V)) * max_ratio)
  }
  sweep(V, 2, pmax(colMeans(V), 1e-12), "/")
}

# ---------------------------------------------------------------------------
# Platform correction
# ---------------------------------------------------------------------------
estimate_platform <- function(Y_mean, B, w_mean, clip = c(0.2, 5), shrink = 0.5) {
  pred <- as.vector(B %*% w_mean)
  Ys <- Y_mean / max(sum(Y_mean), 1e-8)
  Ps <- pred / max(sum(pred), 1e-8)
  ok <- (Ps > 1e-8) & (Ys > 1e-8)
  lg <- rep(0, length(Ys))
  lg[ok] <- log(Ys[ok] + 1e-8) - log(Ps[ok] + 1e-8)
  if (any(ok)) lg[ok] <- lg[ok] - median(lg[ok])
  lg <- lg * (1 - shrink)
  g <- pmin(pmax(exp(lg), clip[1]), clip[2])
  g[!ok] <- 1
  g
}

# ---------------------------------------------------------------------------
# Reference implementation of the sweeps in plain R (engine = "pureR"). Same
# result as gs_sweeps_cpp(), much slower; kept to cross-check the C++ engine.
# ---------------------------------------------------------------------------
gs_sweeps_pureR <- function(Qs_flat, Cs, W, A, deg, lam_eff, n_sweeps, tol,
                            use_spatial) {
  N <- nrow(Cs); K <- ncol(Cs)
  Rs <- vector("list", N)
  for (i in seq_len(N)) {
    Q <- matrix(Qs_flat[i, ], K, K)
    Q <- 0.5 * (Q + t(Q))
    if (use_spatial) diag(Q) <- diag(Q) + lam_eff * deg[i]
    sc <- max(sum(diag(Q)) / K, 1e-12)
    R <- tryCatch(chol(Q), error = function(e) chol(Q + 1e-6 * sc * diag(K)))
    Rs[[i]] <- R
  }
  Ai <- as(A, "TsparseMatrix")
  nb_list <- split(seq_along(Ai@i), Ai@j + 1)
  for (s in seq_len(n_sweeps)) {
    W_old <- W
    for (i in seq_len(N)) {
      cc <- Cs[i, ]
      if (use_spatial) {
        e <- nb_list[[as.character(i)]]
        if (!is.null(e)) {
          m <- W[, Ai@i[e] + 1, drop = FALSE] %*% Ai@x[e]
          cc <- cc + lam_eff * as.vector(m)
        }
      }
      z <- backsolve(Rs[[i]], cc, transpose = TRUE)
      W[, i] <- nnls::nnls(Rs[[i]], z)$x
    }
    if (mean(abs(W - W_old)) < tol) break
  }
  W
}

# ---------------------------------------------------------------------------
# Task A: spot deconvolution
# ---------------------------------------------------------------------------
deconvolve_spots <- function(ST, sig, gene_idx = NULL, coords = NULL,
                             lam = 0, phi = 0.1, weights = "nb",
                             n_outer = 3, n_sweeps = 6, knn = 6,
                             platform_correction = TRUE, tol = 1e-4,
                             engine = c("rcpp", "pureR"), verbose = TRUE,
                             cell_counts = NULL, alpha_scale = 0.1) {
  engine <- match.arg(engine)
  t0 <- Sys.time()

  if (is.null(gene_idx)) {
    gene_idx <- seq_len(nrow(sig$B))
  } else if (is.character(gene_idx)) {
    # gene names are matched against sig$gene_names
    m <- match(gene_idx, sig$gene_names)
    if (anyNA(m)) stop(sum(is.na(m)), " gene(s) in gene_idx are not in sig$gene_names")
    gene_idx <- m
  }


  # SUBTLE: the marker set captures a different fraction f_k of each cell
  # type's mRNA (1.2-1.5x spread). Renormalising B's columns after subsetting
  # rescales column k by 1/f_k, so w_k ends up proportional to n_k*s_k*f_k.
  # The correct divisor for the cell-proportion conversion is s_k * f_k.
  f_capture <- pmax(colSums(sig$B[gene_idx, , drop = FALSE]), 1e-12)
  eff_cell_sizes <- sig$cell_sizes * f_capture

  B <- sig$B[gene_idx, , drop = FALSE]
  B <- sweep(B, 2, pmax(colSums(B), 1e-12), "/")

  Yr <- ST[, gene_idx, drop = FALSE]
  Yr <- Yr / pmax(rowSums(Yr), 1e-12)
  Y <- t(Yr)                      # (G, N)
  G <- nrow(B); K <- ncol(B); N <- ncol(Y)

  use_spatial <- !is.null(coords) && lam > 0
  if (use_spatial) {
    A <- build_spatial_graph(coords, k = knn)
    deg <- rowSums(A)
  } else { A <- NULL; deg <- rep(0, N) }

  # init: plain NNLS
  W <- matrix(0, K, N)
  Q0 <- t(B) %*% B
  C0 <- t(B) %*% Y
  R0 <- chol(Q0 + 1e-8 * sum(diag(Q0)) / K * diag(K))
  for (i in seq_len(N)) {
    z <- backsolve(R0, C0[, i], transpose = TRUE)
    W[, i] <- nnls::nnls(R0, z)$x
  }

  for (outer in seq_len(n_outer)) {
    if (platform_correction && outer > 1) {
      g <- estimate_platform(rowMeans(Y), B, rowMeans(W))
      B <- B * g
      B <- sweep(B, 2, pmax(colSums(B), 1e-12), "/")
    }
    mu <- B %*% W
    V <- weights_from_mu(mu, phi = phi, scheme = weights)

    if (engine == "rcpp") {
      Qs_flat <- all_grams_cpp(B, V)
    } else {
      P <- matrix(0, G, K * K)
      for (k in seq_len(K)) for (l in seq_len(K)) P[, (k-1)*K + l] <- B[, k] * B[, l]
      Qs_flat <- t(V) %*% P
    }
    Cs <- t(V * Y) %*% B

    lam_eff <- 0
    if (use_spatial) {
      tr <- rowSums(Qs_flat[, seq(1, K*K, by = K + 1), drop = FALSE])
      q_scale <- mean(tr) / K
      lam_eff <- lam * q_scale / max(mean(deg), 1e-12)
    }

    if (engine == "rcpp") {
      Asp <- if (use_spatial) as(as(A, "generalMatrix"), "CsparseMatrix")
             else sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(N, N))
      W <- gs_sweeps_cpp(as.matrix(Qs_flat), as.matrix(Cs), W, Asp, deg,
                         lam_eff, n_sweeps, tol, use_spatial)
    } else {
      W <- gs_sweeps_pureR(as.matrix(Qs_flat), as.matrix(Cs), W, A, deg,
                           lam_eff, n_sweeps, tol, use_spatial)
    }
  }

  W <- pmax(W, 0)
  W_mrna <- sweep(W, 2, pmax(colSums(W), 1e-12), "/")
  W_cell <- mrna_to_cell(W_mrna, eff_cell_sizes)

  # 0.1.1: keep the identity of spots and cell types on the returned matrices.
  # The solver works on a bare matrix; without dimnames the caller had no way to
  # align the result to its own metadata, and the cell_counts lookup below
  # silently matched nothing (see NEWS.md).
  spot_ids <- rownames(ST)
  if (is.null(spot_ids)) spot_ids <- colnames(W_mrna)
  dimnames(W_mrna) <- list(sig$types, spot_ids)
  dimnames(W_cell) <- list(sig$types, spot_ids)

  # ---- Optional: cell-count-informed Dirichlet shrinkage ----
  W_cell_refined <- NULL
  if (!is.null(cell_counts)) {
    K <- nrow(W_cell)
    # 0.1.1: accept either a plain vector in spot order or a vector named by
    # spot id, and fail loudly if it lines up with nothing. Before 0.1.1 a
    # mismatch returned the unrefined result with no warning at all.
    if (is.null(names(cell_counts))) {
      if (length(cell_counts) != ncol(W_cell)) {
        stop("cell_counts must be named by spot, or have one value per spot: ",
             length(cell_counts), " values for ", ncol(W_cell), " spots")
      }
      names(cell_counts) <- colnames(W_cell)
    }
    # Empirical Bayes prior: global mean proportion of each cell type
    global_mean <- rowMeans(W_cell, na.rm = TRUE)
    alpha_prior <- global_mean * alpha_scale * K
    # Dirichlet posterior: (N*p + alpha) / (N + sum(alpha))
    # Fully vectorized for speed
    common_spots <- intersect(colnames(W_cell), names(cell_counts))
    if (length(common_spots) == 0) {
      stop("no names(cell_counts) match the spot identifiers in ST; pass ",
           "cell_counts keyed by the same spot ids (e.g. rownames(ST))")
    }
    W_cell_refined <- W_cell
    N_vec <- cell_counts[common_spots]
    W_sub  <- W_cell[, common_spots, drop = FALSE]
    W_cell_refined[, common_spots] <- sweep(
      sweep(W_sub, 2, N_vec, "*") + alpha_prior,
      2, N_vec + sum(alpha_prior), "/")
  }

  list(W_mrna = W_mrna,
         W_cell = if (is.null(W_cell_refined)) W_cell else W_cell_refined,
         W_cell_raw = if (!is.null(W_cell_refined)) W_cell else NULL,
         f_capture = f_capture,
         types = sig$types,
         runtime = as.numeric(difftime(Sys.time(), t0, units = "secs")),
         n_genes_used = G)
}

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
rmse <- function(a, b) sqrt(mean((a - b)^2))
pcc  <- function(a, b) {
  a <- as.vector(a); b <- as.vector(b)
  if (sd(a) < 1e-12 || sd(b) < 1e-12) return(0)
  cor(a, b)
}
jsd <- function(P, Q, eps = 1e-12) {
  P <- P + eps; Q <- Q + eps
  P <- sweep(P, 2, colSums(P), "/"); Q <- sweep(Q, 2, colSums(Q), "/")
  M <- 0.5 * (P + Q)
  kl <- function(a, b) colSums(a * log2(a / b))
  mean(0.5 * kl(P, M) + 0.5 * kl(Q, M))
}

# ===========================================================================
# Spatial smoothing and per-section niches (paper pipeline)
# ===========================================================================

smooth_by_graph <- function(W, A, n_iter = 2, alpha = 0.6) {
  d <- pmax(rowSums(A), 1e-12)
  out <- W
  for (i in seq_len(n_iter)) {
    nb <- as.matrix(out %*% A) / rep(d, each = nrow(out))
    out <- (1 - alpha) * out + alpha * nb
  }
  out
}

#' Find the niches of one spatial section
#'
#' `method = "aitchison"` (default, from SpaCoDA 0.2.0) is the pipeline used in the
#' paper: Aitchison distance on the composition fused with spatial distance,
#' PCoA, shared-nearest-neighbour graph, and the clustering resolution with
#' maximum modularity -- see [cluster_single_sample()], to which the arguments
#' `expr_weight` ... `pcoa_dim` are forwarded. `n_niches`, `knn`, `smooth_iter`
#' and `alpha` belong to the legacy `method = "kmeans"` path and are unused by
#' the default.
#'
#' @param W_cell cell types x spots composition, e.g. `res$W_cell` from
#'   [deconvolve_spots()]; the column names are the spot identifiers
#' @param coords spots x 2 coordinates, in the same order as the columns of
#'   `W_cell` (rownames are matched when present)
#' @param method `"aitchison"` (paper pipeline) or `"kmeans"` (pre-0.2.0
#'   behaviour: spatial smoothing of the proportions followed by k-means)
#' @return list with `labels` (one niche per spot), `profiles` (niches x cell
#'   types, the niche signature), `smoothed` (legacy path only), plus for the
#'   default method `best_resolution`, `scan`, `homogeneous` and `seurat_obj`.
#' @export
define_niches <- function(W_cell, coords,
                          method = c("aitchison", "kmeans"),
                          n_niches = 6, knn = 6, smooth_iter = 2, alpha = 0.6,
                          expr_weight = 0.7, spat_weight = 0.3,
                          k_nn_spatial = 10, gaussian_kernel_distance = FALSE,
                          use_rank_norm = FALSE, k_nn_snn = 20,
                          resolution_seq = seq(0.3, 3, by = 0.1),
                          pcoa_dim = 30, min_effective_ct = 2,
                          seed = 42, verbose = TRUE,
                          pcoa_method = c("auto", "partial", "full"),
                          spatial_test = TRUE, n_perm = 999) {
  method <- match.arg(method)
  pcoa_method <- match.arg(pcoa_method)

  if (method == "aitchison") {
    ## Someone porting a pre-0.2.0 script will pass n_niches/knn/smooth_iter/alpha
    ## and would otherwise silently get a different number of niches.
    if (!missing(n_niches) || !missing(knn) || !missing(smooth_iter) || !missing(alpha)) {
      warning("n_niches / knn / smooth_iter / alpha belong to method = \"kmeans\" and are ",
              "ignored by the default method = \"aitchison\", which chooses the number of ",
              "niches from resolution_seq. Pass method = \"kmeans\" for the pre-0.2.0 behaviour.",
              call. = FALSE)
    }
    prop <- t(as.matrix(W_cell))
    if (is.null(rownames(prop))) rownames(prop) <- paste0("spot", seq_len(nrow(prop)))
    if (is.null(colnames(prop))) colnames(prop) <- paste0("type", seq_len(ncol(prop)))
    co <- as.matrix(coords)
    if (is.null(rownames(co)) || !all(rownames(prop) %in% rownames(co))) {
      if (nrow(co) != nrow(prop))
        stop("coords must have one row per spot (or rownames matching the spots)")
      rownames(co) <- rownames(prop)
    } else {
      co <- co[rownames(prop), , drop = FALSE]
    }
    cs <- cluster_single_sample(prop_mat = prop, coords_mat = co,
                                expr_weight = expr_weight, spat_weight = spat_weight,
                                k_nn_spatial = k_nn_spatial,
                                gaussian_kernel_distance = gaussian_kernel_distance,
                                use_rank_norm = use_rank_norm, k_nn_snn = k_nn_snn,
                                resolution_seq = resolution_seq, pcoa_dim = pcoa_dim,
                                seed = seed, min_effective_ct = min_effective_ct,
                                verbose = verbose, pcoa_method = pcoa_method,
                                spatial_test = spatial_test, n_perm = n_perm)
    profiles <- as.matrix(cs$prototypes)
    return(list(labels = cs$spot_clusters, profiles = profiles, smoothed = NULL,
                best_resolution = cs$best_resolution, best_k = cs$best_k,
                best_modularity = cs$best_modularity, scan = cs$scan_results,
                homogeneous = cs$homogeneous, seurat_obj = cs$seurat_obj,
                spatial_test = cs$spatial_test))
  }

  ## ---- legacy path: spatial smoothing + k-means (pre-0.2.0) ---------------
  if (verbose) message("define_niches(): legacy k-means path (method = 'kmeans')")
  A <- build_spatial_graph(coords, k = knn)
  Wn <- smooth_by_graph(W_cell, A, smooth_iter, alpha)
  Wn <- sweep(Wn, 2, pmax(colSums(Wn), 1e-12), "/")
  set.seed(seed)
  km <- kmeans(t(Wn), centers = n_niches, nstart = 10)
  profiles <- matrix(0, n_niches, nrow(W_cell))
  for (j in seq_len(n_niches)) {
    m <- km$cluster == j
    if (any(m)) profiles[j, ] <- rowMeans(W_cell[, m, drop = FALSE])
  }
  # 0.1.1: carry the cell-type (and niche) identities through, so that
  # profiles / labels can be matched to metadata instead of by position only.
  # 0.2.2: profile rownames are the label values ("1".."J", as in the default
  # path where they are "0".."J-1"), so profiles and labels match by name.
  # They were "niche1".."nicheJ" in 0.1.1 - 0.2.1.
  if (!is.null(rownames(W_cell))) colnames(profiles) <- rownames(W_cell)
  rownames(profiles) <- as.character(seq_len(n_niches))
  labels <- km$cluster
  if (!is.null(colnames(W_cell))) names(labels) <- colnames(W_cell)
  list(labels = labels, profiles = profiles, smoothed = Wn)
}
