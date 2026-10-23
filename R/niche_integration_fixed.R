# =============================================================================
# CLR transform, resolution selection with a null test, and niche QC
# =============================================================================
# Helpers of the legacy CLR integration path (integrate_slices(method = "clr")),
# which was written to avoid the following problems of a naive multi-slice
# PCoA pipeline:
#   1. per-slice PCoA -> slice axes not comparable -> Harmony cannot fix them
#   2. merge(add.cell.ids = ) adds prefixes twice -> embedding row names mismatch
#   3. cross-slice spatial distance set to 0 -> claims "identical location" ->
#      artificially glues slices together
#   4. cmdscale is O(N^3): 96 s at N = 4000, ~43 min at N = 12000
#   5. modularity finds 9-11 "niches" in pure noise (ARI 0.04) without a null test
#
# Aitchison (CLR) distance is invariant to a per-cell-type scaling of all spots
# (max diff 1.5e-05, corr 1.000), so an arbitrary scale of the composition
# scores does not affect the clustering.

# -----------------------------------------------------------------------------
# CLR transform, computed directly (no vegdist/cmdscale detour)
# -----------------------------------------------------------------------------
# Key fact: Aitchison distance == Euclidean distance on CLR-transformed data.
# So there is no reason to build an N x N distance matrix and run cmdscale.
# PCA on the CLR matrix is mathematically equivalent (PCoA of Euclidean == PCA)
# but costs O(N*K) memory instead of O(N^2) and O(N*K^2) time instead of O(N^3).
# This lifts the practical limit from ~6000 spots to hundreds of thousands.
clr_transform <- function(prop_mat, pseudo = NULL) {
  P <- as.matrix(prop_mat)
  if (any(P < 0)) {
    warning("Negative values clipped to 0. If the input comes from z-scored ",
            "values, you are using the wrong function: the Aitchison pipeline ",
            "requires non-negative, unstandardized scores.")
    P[P < 0] <- 0
  }
  # Zero handling: cmultRepl expects COUNTS, not proportions; feeding
  # proportions makes it infer a detection limit from row sums of 1, which is
  # meaningless. We use a simple multiplicative replacement instead.
  if (is.null(pseudo)) {
    pos <- P[P > 0]
    pseudo <- if (length(pos)) min(pos) * 0.65 else 1e-6
  }
  zero_frac <- mean(P == 0)
  if (zero_frac > 0) {
    P[P == 0] <- pseudo
    P <- P / rowSums(P)
  }
  L <- log(P)
  L - rowMeans(L)          # CLR: subtract the per-row log mean
}

# -----------------------------------------------------------------------------
# Resolution selection: modularity + null test
# -----------------------------------------------------------------------------
# When structure exists, max-modularity is near-optimal (with true k = 5 it
# selected k = 5, ARI 0.931 vs a best possible 0.935). When no structure exists
# it still returns 9-11 "niches" (ARI 0.04). Hence the null test: rewire the
# graph keeping the degree sequence, re-cluster, and compare modularity; if the
# real modularity is not clearly higher, the niches are noise. (For a test that
# does not depend on the clustering at all, see niche_spatial_test().)
select_resolution <- function(obj, resolution_seq, graph_name = "Niche_snn",
                              n_null = 10, verbose = TRUE) {
  g <- igraph::graph_from_adjacency_matrix(obj@graphs[[graph_name]],
                                           mode = "undirected", weighted = TRUE)
  scan <- data.frame(resolution = resolution_seq, k = NA_integer_,
                     modularity = NA_real_)
  for (i in seq_along(resolution_seq)) {
    # Run FindClusters only once per resolution. The original code ran it
    # twice (once for modularity, once for k), doubling the cost for nothing.
    o <- Seurat::FindClusters(obj, resolution = resolution_seq[i], verbose = FALSE)
    cl <- o$seurat_clusters
    scan$k[i] <- length(unique(cl))
    if (scan$k[i] > 1 && scan$k[i] < ncol(obj)) {
      scan$modularity[i] <- igraph::modularity(g, as.numeric(cl),
                                               weights = igraph::E(g)$weight)
    }
  }
  valid <- scan[!is.na(scan$modularity), ]
  if (!nrow(valid)) stop("no resolution produced a valid clustering")
  best <- valid[which.max(valid$modularity), ]

  # Null: rewire the graph keeping the degree sequence, measure the baseline
  null_mod <- replicate(n_null, {
    gr <- igraph::rewire(g, igraph::keeping_degseq(niter = igraph::ecount(g) * 5))
    cl <- igraph::cluster_louvain(gr, resolution = best$resolution)
    igraph::modularity(gr, igraph::membership(cl))
  })
  z <- (best$modularity - mean(null_mod)) / max(stats::sd(null_mod), 1e-9)

  if (verbose) {
    message(sprintf("best resolution=%.2f  k=%d  modularity=%.4f", best$resolution,
                    best$k, best$modularity))
    message(sprintf("null modularity=%.4f (sd %.4f)  z=%.1f", mean(null_mod),
                    stats::sd(null_mod), z))
    if (z < 3) message("!! z < 3: modularity is not significantly above random graphs. Your niches are likely noise.")
  }
  list(best_resolution = best$resolution, best_k = best$k,
       best_modularity = best$modularity, null_modularity = mean(null_mod),
       z = z, scan = scan)
}

# -----------------------------------------------------------------------------
# Diagnostics: run this before trusting the result
# -----------------------------------------------------------------------------
niche_qc <- function(res) {
  obj <- res$object
  tab <- res$niche_summary_prop
  cat("=== niche QC ===\n")
  cat("  slices:", ncol(tab), "  niches:", nrow(tab), "\n")
  cat("  shared niches:", length(res$common_niches), "/", nrow(tab), "\n")
  # A niche dominated by a single slice is likely residual batch, not biology
  dom <- apply(tab, 1, max)
  bad <- names(dom)[dom > 0.8]
  if (length(bad)) {
    cat("  !! these niches are dominated (>80%) by one slice, probably uncorrected batch:\n     ",
        paste(bad, collapse = ", "), "\n")
  }
  # Slice mixing: entropy of the slice distribution inside each niche
  ent <- apply(tab, 1, function(x) { p <- x / sum(x); p <- p[p > 0]; -sum(p * log(p)) })
  cat("  slice-mixing entropy (higher = better, max ", round(log(ncol(tab)), 2), "):",
      paste(round(ent, 2), collapse = " "), "\n")
  invisible(list(dominance = dom, entropy = ent))
}
