# =============================================================================
# Cross-sample integration: super-niches
# =============================================================================
# Two implementations live here.
#
#   method = "prototype"  (default, from SpaCoDA 0.2.0) -- the pipeline used for the
#     BRCA and GBM analyses. Per section, cluster_single_sample() finds the local
#     niches and their prototypes; the prototypes of all sections are then put
#     into one matrix (cell type x prototype), normalised with LogNormalize,
#     scaled, reduced with PCA, batch-corrected with Harmony on the section
#     label, and clustered into super-niches. The clustering resolution is an
#     argument (`resolution`); it was hard-coded to 6 in the study scripts.
#
#   method = "clr"  -- the pre-0.2.0 implementation (kept available): CLR per
#     section, within-section smoothing in CLR space, joint PCA over all spots,
#     Harmony, then a resolution chosen by modularity against a rewired-graph
#     null test (select_resolution()).
#
# Super-niche = a composition state that recurs across sections.
# =============================================================================

integrate_slices_clr <- function(prop_list, coord_list, slice_names,
                             spatial_smooth = 0.5, knn_spatial = 6,
                             n_pcs = 20, do_harmony = TRUE,
                             k_nn_snn = 20,
                             resolution_seq = seq(0.2, 1.5, by = 0.1),
                             seed = 42, verbose = TRUE) {
  set.seed(seed)
  stopifnot(length(prop_list) == length(coord_list),
            length(prop_list) == length(slice_names))
  n_slices <- length(prop_list)

  # ---- 1. Within-slice spatial smoothing, done in CLR space ---------------
  clr_list <- list(); meta_list <- list()
  for (i in seq_len(n_slices)) {
    P <- as.matrix(prop_list[[i]]); co <- as.matrix(coord_list[[i]])
    stopifnot(nrow(P) == nrow(co))
    L <- clr_transform(P)

    if (spatial_smooth > 0) {
      kn <- FNN::get.knn(co, k = knn_spatial)
      # Smoothing in CLR space is meaningful because CLR is Euclidean;
      # averaging raw proportions and renormalizing uses the wrong geometry.
      nb <- matrix(0, nrow(L), ncol(L))
      for (j in seq_len(knn_spatial)) nb <- nb + L[kn$nn.index[, j], , drop = FALSE]
      nb <- nb / knn_spatial
      L <- (1 - spatial_smooth) * L + spatial_smooth * nb
    }
    rn <- paste0(slice_names[i], "|", rownames(P))
    rownames(L) <- rn
    clr_list[[i]] <- L
    meta_list[[i]] <- data.frame(slice = slice_names[i],
                                 x = co[, 1], y = co[, 2],
                                 row.names = rn, stringsAsFactors = FALSE)
  }

  L_all <- do.call(rbind, clr_list)
  meta_all <- do.call(rbind, meta_list)
  if (verbose) message("merged: ", nrow(L_all), " spots x ", ncol(L_all), " cell types")

  # ---- 2. Joint PCA (== joint Aitchison PCoA, O(N*K^2) not O(N^3)) --------
  # This is the core bug fix: all slices must land on the same basis. With
  # independent PCoA the axes rotate/flip per slice (measured axis-2
  # correlation -0.459 and axis-3 only +0.112 when slice abundances differ),
  # so stacking them and feeding Harmony was meaningless.
  n_pcs <- min(n_pcs, ncol(L_all) - 1, nrow(L_all) - 1)
  Lc <- scale(L_all, center = TRUE, scale = FALSE)
  pca <- irlba::prcomp_irlba(Lc, n = n_pcs, center = FALSE, scale. = FALSE)
  emb <- pca$x
  rownames(emb) <- rownames(L_all)
  colnames(emb) <- paste0("PC_", seq_len(ncol(emb)))
  if (verbose) message("joint PCA: ", n_pcs, " PCs, variance explained ",
                       round(100 * sum(pca$sdev^2) / sum(apply(Lc, 2, stats::var)), 1), "%")

  # ---- 3. Seurat object ----------------------------------------------------
  # Note: row names already carry the slice prefix, so merge() must NOT add
  # cell.ids again (would produce slice1_slice1_spot1 and break the embedding
  # alignment). We build a single object directly and avoid merge entirely.
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Package 'Seurat' is required for integrate_slices().")
  obj <- Seurat::CreateSeuratObject(counts = t(exp(L_all)), assay = "Niche",
                                    meta.data = meta_all)
  obj[["pca"]] <- Seurat::CreateDimReducObject(embeddings = emb, key = "PC_",
                                               assay = "Niche")

  red <- "pca"
  if (do_harmony && n_slices > 1) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      warning("Package 'harmony' not installed; skipping batch correction.")
    } else {
      if (verbose) message("Harmony batch correction across slices ...")
      # A per-slice multiplicative scale of the composition scores becomes a
      # slice-specific constant shift in CLR space, which Harmony removes.
      # 0.1.1: project.dim = FALSE. Harmony >= 1.2 defaults it to TRUE, which
      # makes RunHarmony call Seurat::ProjectDim(); ProjectDim requires a
      # scale.data slot, while this assay deliberately holds CLR values, so the
      # call aborted with "data.use %*% cell.embeddings : non-conformable
      # arguments" *after* Harmony had already converged (see NEWS.md).
      nclust_h <- max(2L, min(as.integer(round(ncol(obj) / 30)), 100L))
      obj <- harmony::RunHarmony(obj, group.by.vars = "slice",
                                 reduction.use = "pca",
                                 dims.use = 1:n_pcs,
                                 nclust = nclust_h,
                                 project.dim = FALSE,
                                 verbose = FALSE)
      red <- "harmony"
    }
  }

  # ---- 4. Graph clustering -------------------------------------------------
  obj <- Seurat::FindNeighbors(obj, reduction = red, dims = 1:n_pcs,
                               k.param = k_nn_snn, verbose = FALSE)
  opt <- select_resolution(obj, resolution_seq, graph_name = "Niche_snn",
                           verbose = verbose)
  obj <- Seurat::FindClusters(obj, resolution = opt$best_resolution,
                              verbose = FALSE)

  # ---- 5. Niche composition statistics ------------------------------------
  tab <- table(obj$seurat_clusters, obj$slice)
  tab_prop <- prop.table(tab, margin = 2)
  is_common <- apply(tab_prop, 1, function(x) all(x > 0.05))

  # 0.1.1: also return the null-test statistics. They were only ever printed as
  # messages, so the caller could not report or threshold them.
  list(object = obj, embeddings = emb, reduction = red,
       best_resolution = opt$best_resolution, scan = opt$scan,
       best_modularity = opt$best_modularity,
       null_modularity = opt$null_modularity, z = opt$z,
       niche_summary = tab, niche_summary_prop = tab_prop,
       common_niches = rownames(tab)[is_common],
       clr = L_all)
}

#' Integrate sections and find cross-section super-niches
#'
#' @param prop_list list of composition matrices, one per section, spots x cell
#'   types (e.g. `t(W_cell)` from [deconvolve_spots()]); rows must be named
#' @param coord_list list of spot coordinate matrices, same order as `prop_list`
#' @param slice_names section names
#' @param method `"prototype"` (default) runs the study pipeline: local niches
#'   and their prototypes per section via [cluster_single_sample()], then one
#'   clustering over all prototypes. `"clr"` is the pre-0.2.0 implementation
#'   (CLR over all spots, joint PCA, modularity scan with a rewired-graph null
#'   test) and ignores `resolution`.
#' @param resolution clustering resolution for the super-niches. The study
#'   scripts hard-coded 6; it is an argument here. `NULL` selects the resolution
#'   with maximum modularity from `resolution_seq` instead.
#' @param resolution_seq candidate resolutions when `resolution = NULL`
#' @param pca_dims PCA dimensions for the prototype embedding (study: 15)
#' @param do_harmony batch-correct the section label with Harmony (needs >= 2 sections)
#' @param per_sample_args named list of extra arguments forwarded to
#'   [cluster_single_sample()], e.g. `list(pcoa_dim = 30, k_nn_snn = 20)`
#' @param spatial_smooth,knn_spatial,n_pcs,k_nn_snn,seed,verbose arguments of the
#'   legacy `method = "clr"` path (with the exception of `seed` and `verbose`,
#'   which the prototype path also uses)
#' @return list with
#'   * `object` -- Seurat object of the prototypes; `Idents()`/`seurat_clusters`
#'     are the super-niche of each prototype,
#'   * `super_niche` -- super-niche label per prototype,
#'   * `spot_super_niche` -- super-niche label per spot (named),
#'   * `prototypes` -- prototype x cell-type matrix used as input,
#'   * `per_sample` -- the [cluster_single_sample()] result of each section,
#'   * `best_resolution`, `scan`, `z`, `reduction`,
#'   * `niche_summary`, `niche_summary_prop`, `common_niches` -- super-niche by
#'     section tables (spot level), as consumed by [niche_qc()].
#' @export
integrate_slices <- function(prop_list, coord_list, slice_names,
                             method = c("prototype", "clr"),
                             resolution = 6,
                             resolution_seq = seq(0.3, 3, by = 0.1),
                             pca_dims = 15,
                             do_harmony = TRUE,
                             per_sample_args = list(),
                             spatial_smooth = 0.5, knn_spatial = 6,
                             n_pcs = 20, k_nn_snn = 20,
                             seed = 42, verbose = TRUE) {
  method <- match.arg(method)

  if (method == "clr") {
    return(integrate_slices_clr(prop_list = prop_list, coord_list = coord_list,
                                slice_names = slice_names, spatial_smooth = spatial_smooth,
                                knn_spatial = knn_spatial, n_pcs = n_pcs,
                                do_harmony = do_harmony, k_nn_snn = k_nn_snn,
                                resolution_seq = resolution_seq, seed = seed,
                                verbose = verbose))
  }

  stopifnot(length(prop_list) == length(coord_list),
            length(prop_list) == length(slice_names))
  n_slices <- length(prop_list)
  set.seed(seed)
  say <- function(...) if (isTRUE(verbose)) message(...)

  ## ---- 1. local niches + prototypes, one section at a time ----------------
  per_sample <- vector("list", n_slices)
  for (i in seq_len(n_slices)) {
    say(sprintf("section %d/%d: %s", i, n_slices, slice_names[i]))
    P <- as.matrix(prop_list[[i]])
    args <- c(list(prop_mat = P, coords_mat = as.matrix(coord_list[[i]]),
                   seed = seed, verbose = verbose), per_sample_args)
    per_sample[[i]] <- do.call(cluster_single_sample, args)
  }
  names(per_sample) <- slice_names

  ## ---- 2. stack the prototypes into one matrix ----------------------------
  proto_list <- lapply(seq_len(n_slices), function(i) {
    pt <- as.matrix(per_sample[[i]]$prototypes)
    rownames(pt) <- paste0(slice_names[i], "_", rownames(pt))
    pt
  })
  all_types <- unique(unlist(lapply(proto_list, colnames)))
  proto_list <- lapply(proto_list, function(pt) {
    miss <- setdiff(all_types, colnames(pt))
    if (length(miss)) {
      pt <- cbind(pt, matrix(0, nrow(pt), length(miss),
                             dimnames = list(rownames(pt), miss)))
    }
    pt[, all_types, drop = FALSE]
  })
  prototypes <- do.call(rbind, proto_list)
  sample_of  <- rep(slice_names, vapply(proto_list, nrow, integer(1)))
  say(sprintf("super-niche input: %d prototypes x %d cell types from %d sections",
              nrow(prototypes), ncol(prototypes), n_slices))

  ## ---- 3. cluster the prototypes ------------------------------------------
  obj <- Seurat::CreateSeuratObject(counts = as(t(prototypes), "CsparseMatrix"),
                                    project = "SpaCoDA_superniche")
  obj$sample_id <- sample_of
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                               scale.factor = 1e4, verbose = FALSE)
  Seurat::VariableFeatures(obj) <- rownames(obj)
  obj <- Seurat::ScaleData(obj, features = rownames(obj), verbose = FALSE)

  npcs <- max(2, min(pca_dims, nrow(obj) - 1, ncol(obj) - 1))
  obj <- Seurat::RunPCA(obj, features = rownames(obj), npcs = npcs, verbose = FALSE)

  red <- "pca"
  if (do_harmony && length(unique(sample_of)) > 1) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      warning("Package 'harmony' not installed; skipping batch correction.")
    } else {
      say("Harmony batch correction across sections ...")
      ## Harmony's default nclust = round(N/30) collapses to 0 for a handful of
      ## prototypes; pass its own rule with a floor of 2.
      nclust_h <- max(2L, min(as.integer(round(ncol(obj) / 30)), 100L))
      obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id",
                                 reduction.use = "pca", dims.use = 1:npcs,
                                 nclust = nclust_h,
                                 project.dim = FALSE, verbose = FALSE)
      red <- "harmony"
    }
  }
  ## The study's k.param = 20 assumes hundreds to thousands of prototypes. On a
  ## handful of prototypes, 20 neighbours is nearly the whole cloud, the SNN graph
  ## loses its structure, and every candidate resolution collapses to k = 1 or
  ## k = N -- so no resolution is valid. Cap k at a third of the prototypes, which
  ## is exactly the study value once there are 60 or more of them (BRCA/GBM had
  ## ~2,500 and ~1,800), leaving study-scale runs bit-identical.
  n_proto <- ncol(obj)
  k_proto <- max(2L, min(as.integer(k_nn_snn), as.integer(n_proto %/% 3)))
  if (k_proto < as.integer(k_nn_snn)) {
    say(sprintf("only %d prototypes: SNN k.param capped at %d (study value %d)",
                n_proto, k_proto, as.integer(k_nn_snn)))
  }
  obj <- Seurat::FindNeighbors(obj, reduction = red, dims = 1:npcs,
                               k.param = k_proto, verbose = FALSE)

  if (is.null(resolution)) {
    ## With very few prototypes every candidate resolution degenerates to k = 1
    ## or k = N, so nothing is valid; fall back to resolution 1 with a warning
    ## instead of aborting.
    opt <- tryCatch(select_resolution_by_modularity(obj, resolution_seq = resolution_seq),
                    error = function(e) {
                      warning("no resolution produced a valid clustering (",
                              conditionMessage(e), "); falling back to resolution = 1")
                      NULL
                    })
    if (is.null(opt)) {
      best_res <- 1
      scan     <- data.frame(resolution = 1, k = NA_integer_, modularity = NA_real_)
    } else {
      best_res <- opt$best_resolution
      scan     <- opt$scan_results
      say(sprintf("resolution selected by maximum modularity: %.2f (k = %d)",
                  best_res, opt$best_k))
    }
  } else {
    best_res <- resolution
    scan     <- data.frame(resolution = resolution, k = NA_integer_, modularity = NA_real_)
  }

  ## Seurat's singleton re-assignment can fail when every cluster is a singleton
  ## (only reachable with very few prototypes). The study scripts used the default
  ## group.singletons = TRUE; fall back rather than abort.
  obj <- tryCatch(Seurat::FindClusters(obj, resolution = best_res, verbose = FALSE),
                  error = function(e) {
                    warning("FindClusters with group.singletons = TRUE failed (",
                            conditionMessage(e), "); retrying with group.singletons = FALSE")
                    Seurat::FindClusters(obj, resolution = best_res,
                                         group.singletons = FALSE, verbose = FALSE)
                  })
  super_of_proto <- obj$seurat_clusters
  names(super_of_proto) <- colnames(obj)

  ## A requested resolution can be far too high for a small prototype set: Louvain
  ## then puts every prototype in its own community and Seurat labels them all
  ## "singleton". The study value 6 assumes the study scale (BRCA ~2,500 and GBM
  ## ~1,800 prototypes); say so instead of returning a meaningless partition.
  if (length(unique(super_of_proto)) == 1L && all(as.character(super_of_proto) == "singleton")) {
    warning(sprintf(paste0("resolution = %s leaves all %d prototypes as singletons: the ",
                           "resolution is too high for this number of prototypes. Use ",
                           "resolution = NULL (select by maximum modularity) or a smaller value."),
                    format(best_res), ncol(obj)), call. = FALSE)
  }

  ## Report the achieved k / modularity even when the resolution was supplied, so
  ## that `scan` is informative in both modes.
  if (is.null(resolution)) {
    ## opt is NULL when the scan found no valid resolution and fell back to 1
    best_modularity <- if (is.null(opt)) NA_real_ else opt$best_modularity
  } else {
    best_modularity <- tryCatch({
      gname <- paste0(Seurat::DefaultAssay(obj), "_snn")
      adj <- obj@graphs[[gname]]
      g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)
      igraph::modularity(g, as.numeric(super_of_proto), weights = igraph::E(g)$weight)
    }, error = function(e) NA_real_)
    scan <- data.frame(resolution = best_res, k = length(unique(super_of_proto)),
                       modularity = best_modularity)
  }

  ## ---- 4. carry the super-niche label back to every spot ------------------
  spot_super <- unlist(lapply(seq_len(n_slices), function(i) {
    ## 0.2.2: take the spot names from the factor; as.character() drops them,
    ## so spot_super_niche used to come back unnamed despite the documentation
    sc  <- per_sample[[i]]$spot_clusters
    ids <- paste0(slice_names[i], "_", as.character(sc))
    stats::setNames(as.character(super_of_proto[ids]), names(sc))
  }), use.names = TRUE)

  slice_of_spot <- rep(slice_names,
                       vapply(prop_list, function(x) nrow(as.matrix(x)), integer(1)))
  tab      <- table(spot_super, slice_of_spot)
  tab_prop <- prop.table(tab, margin = 2)
  is_common <- apply(tab_prop, 1, function(x) all(x > 0.05))

  ## per-section T1 results side by side (NULL entries when spatial_test = FALSE)
  st_row <- function(st) data.frame(ratio = st$ratio, z = st$z, p_value = st$p_value)
  spatial_tests <- do.call(rbind, lapply(per_sample, function(r)
    if (is.null(r$spatial_test)) st_row(list(ratio = NA_real_, z = NA_real_, p_value = NA_real_))
    else st_row(r$spatial_test)))
  rownames(spatial_tests) <- slice_names

  list(object = obj, embeddings = Seurat::Embeddings(obj, red), reduction = red,
       spatial_tests = spatial_tests,
       best_resolution = best_res, best_k = length(unique(super_of_proto)),
       best_modularity = best_modularity, null_modularity = NA_real_,
       scan = scan, z = NA_real_,
       super_niche = super_of_proto, spot_super_niche = spot_super,
       prototypes = prototypes, per_sample = per_sample,
       niche_summary = tab, niche_summary_prop = tab_prop,
       common_niches = rownames(tab)[is_common])
}
