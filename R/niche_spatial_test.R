# =============================================================================
# Is the composition spatially organised at all?  (permutation test, T1)
# =============================================================================
# The fused composition + space distance of cluster_single_sample() tiles even a
# section with structure-free compositions into spatially contiguous "niches"
# (measured: 7 niches, modularity 0.61, 72% of neighbouring spots in the same
# niche vs 15% by chance). A test on the clustering cannot catch this: shuffling
# compositions keeps the composition clusters, so a modularity statistic under
# that null is *lower* than the observed one exactly when structure is clear.
#
# This test does not look at the clustering. Its statistic is the mean squared
# Aitchison (CLR) distance between each spot and its k nearest spatial
# neighbours -- a multivariate Geary's C. The null shuffles the compositions
# across spots with the coordinates fixed. Small values mean neighbours are more
# alike than random spots, i.e. the composition is spatially organised.
#
# Measured on 800-spot simulations (999 permutations): pure noise ratio 1.003,
# p = 0.68; weak, clear and smooth-gradient structure ratio 0.79 / 0.36 / 0.37,
# all p = 0.001 (the minimum). On 200 independent noise sections the rejection
# rate was 5.5% at 0.05 and 10.0% at 0.10, i.e. calibrated. 20,000 spots with
# 999 permutations take ~6 s.
# A significant result says the composition has spatial structure. It does NOT
# say the structure is made of discrete niches -- a smooth gradient is also
# significant.

#' Permutation test for spatial organisation of the composition
#'
#' @param prop_mat spots x cell types, non-negative proportions
#' @param coords spots x 2 coordinates, same row order as `prop_mat`
#' @param k spatial nearest neighbours per spot
#' @param n_perm number of permutations
#' @param pseudo pseudocount before the CLR transform (as in
#'   [compute_aitchison_distance()])
#' @param seed random seed. The caller's random-number stream is left untouched.
#' @return list with `statistic` (mean squared CLR distance between spatial
#'   neighbours), `null_mean`, `null_sd`, `ratio` (`statistic / null_mean`, < 1
#'   means neighbours are alike), `z`, `p_value` (one-sided, small = organised),
#'   `n_perm`, `k`, `n_spots`, and `note` when the test is not applicable.
#' @export
niche_spatial_test <- function(prop_mat, coords, k = 6, n_perm = 999,
                               pseudo = 1e-6, seed = 42) {
  P  <- as.matrix(prop_mat)
  co <- as.matrix(coords)
  if (nrow(P) != nrow(co)) stop("prop_mat and coords must have one row per spot")
  n <- nrow(P)
  na_result <- function(note)
    list(statistic = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
         ratio = NA_real_, z = NA_real_, p_value = NA_real_, n_perm = 0L,
         k = as.integer(k), n_spots = n, note = note)

  P <- P[, colSums(P) > 0, drop = FALSE]
  if (ncol(P) < 2) return(na_result("fewer than two effective cell types"))
  if (n < k + 2)   return(na_result("too few spots for the neighbour graph"))
  if (n_perm < 1)  stop("n_perm must be >= 1")

  P <- P + pseudo
  P <- P / rowSums(P)
  L <- log(P)
  L <- L - rowMeans(L)

  # sum over directed neighbour pairs of ||L_i - L_j||^2
  #   = k * sum_i ||L_i||^2 + sum_j indeg_j ||L_j||^2 - 2 * sum(L * (A %*% L))
  # evaluated with a sparse A, so no (n*k) x K pair matrix is ever built
  nn <- FNN::get.knn(co, k = k)$nn.index
  A <- Matrix::sparseMatrix(i = rep(seq_len(n), k), j = as.vector(nn), x = 1,
                            dims = c(n, n))
  indeg <- Matrix::colSums(A)
  n_pairs <- n * k
  stat_of <- function(Lp) {
    q <- rowSums(Lp^2)
    (k * sum(q) + sum(indeg * q) - 2 * sum(Lp * as.matrix(A %*% Lp))) / n_pairs
  }
  obs <- stat_of(L)

  # local RNG: leave the caller's stream as it was (downstream Harmony / k-means
  # results must not depend on whether this test ran)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(b) stat_of(L[sample.int(n), , drop = FALSE]),
                 numeric(1))

  nsd <- stats::sd(null)
  list(statistic = obs, null_mean = mean(null), null_sd = nsd,
       ratio = obs / mean(null),
       z = if (is.finite(nsd) && nsd > 0) (obs - mean(null)) / nsd else NA_real_,
       p_value = (1 + sum(null <= obs)) / (n_perm + 1),
       n_perm = as.integer(n_perm), k = as.integer(k), n_spots = n, note = NA_character_)
}
