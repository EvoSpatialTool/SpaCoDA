# 0.2.2 niche additions
#   1. the partial-eigendecomposition PCoA reproduces cmdscale() and the niches
#   2. niche_spatial_test(): calibrated on noise, powered on structure, and it
#      changes neither the niches nor the caller's random-number stream

.sec <- function(n = 400, K = 6, signal = 3, seed = 1) {
  set.seed(seed)
  co <- cbind(x = runif(n, 0, 100), y = runif(n, 0, 100))
  ctr <- cbind(runif(4, 0, 100), runif(4, 0, 100))
  dom <- apply(co, 1, function(p) which.min(colSums((t(ctr) - p)^2)))
  prof <- matrix(rexp(4 * K)^2, 4, K)
  P <- signal * prof[dom, ] + matrix(rexp(n * K), n, K); P <- P / rowSums(P)
  dimnames(P) <- list(paste0("s", seq_len(n)), paste0("t", seq_len(K)))
  rownames(co) <- rownames(P)
  list(P = P, co = co)
}

test_that("partial PCoA gives cmdscale's coordinates up to the sign of each axis", {
  skip_if_not_installed("RSpectra")
  d <- .sec(300)
  D <- fuse_distances(compute_aitchison_distance(d$P), as.matrix(dist(d$co)))
  full <- stats::cmdscale(D, k = 10, eig = TRUE)$points
  part <- SpaCoDA:::pcoa_partial(D, 10)
  expect_equal(dim(part), dim(full))
  sgn <- sign(colSums(full * part))
  expect_equal(unname(sweep(part, 2, sgn, "*")), unname(full), tolerance = 1e-7)
  expect_equal(as.vector(dist(part)), as.vector(dist(full)), tolerance = 1e-8)
})

test_that("pcoa_method does not change the niches", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("RSpectra")
  d <- .sec(600, seed = 2)
  rs <- seq(0.3, 1.5, by = 0.3)
  full <- suppressWarnings(cluster_single_sample(d$P, d$co, resolution_seq = rs,
                                                 pcoa_method = "full", spatial_test = FALSE,
                                                 verbose = FALSE))
  part <- suppressWarnings(cluster_single_sample(d$P, d$co, resolution_seq = rs,
                                                 pcoa_method = "partial", spatial_test = FALSE,
                                                 verbose = FALSE))
  expect_identical(as.character(part$spot_clusters), as.character(full$spot_clusters))
  expect_equal(part$prototypes, full$prototypes)
  expect_equal(part$best_resolution, full$best_resolution)
  expect_equal(part$scan_results, full$scan_results)
})

test_that("niche_spatial_test is calibrated on noise and powered on structure", {
  set.seed(3)
  n <- 400
  co <- cbind(runif(n), runif(n))
  noise <- matrix(rexp(n * 5), n, 5); noise <- noise / rowSums(noise)
  grad  <- cbind(co[, 1], 1 - co[, 1], 0.3, 0.3) + matrix(rexp(n * 4, 8), n, 4)
  grad  <- grad / rowSums(grad)
  tn <- niche_spatial_test(noise, co, n_perm = 199)
  tg <- niche_spatial_test(grad,  co, n_perm = 199)
  expect_gt(tn$p_value, 0.05)
  expect_lt(tg$p_value, 0.01)
  expect_lt(tg$ratio, 1)
  expect_true(abs(tn$ratio - 1) < 0.05)
  ## the statistic is the mean squared CLR distance between kNN pairs
  L <- log(sweep(noise + 1e-6, 1, rowSums(noise + 1e-6), "/"))
  L <- L - rowMeans(L)
  nn <- FNN::get.knn(co, k = 6)$nn.index
  direct <- mean(rowSums((L[rep(seq_len(n), 6), ] - L[as.vector(nn), ])^2))
  expect_equal(tn$statistic, direct, tolerance = 1e-10)
  ## deterministic for a seed, and not applicable to a one-type section
  expect_identical(niche_spatial_test(grad, co, n_perm = 49, seed = 7)$p_value,
                   niche_spatial_test(grad, co, n_perm = 49, seed = 7)$p_value)
  one <- niche_spatial_test(cbind(rep(1, n), 0), co, n_perm = 9)
  expect_true(is.na(one$p_value))
  expect_match(one$note, "effective cell types")
})

test_that("the spatial test leaves the caller's RNG stream untouched", {
  d <- .sec(150)
  set.seed(99); before <- .Random.seed
  niche_spatial_test(d$P, d$co, n_perm = 19)
  expect_identical(.Random.seed, before)
  set.seed(99); a <- runif(3)
  set.seed(99); niche_spatial_test(d$P, d$co, n_perm = 19); b <- runif(3)
  expect_identical(a, b)
})

test_that("running the spatial test does not change the niches", {
  skip_if_not_installed("Seurat")
  d <- .sec(300, seed = 4)
  rs <- seq(0.3, 1.5, by = 0.3)
  on  <- suppressWarnings(cluster_single_sample(d$P, d$co, resolution_seq = rs,
                                                n_perm = 99, verbose = FALSE))
  off <- suppressWarnings(cluster_single_sample(d$P, d$co, resolution_seq = rs,
                                                spatial_test = FALSE, verbose = FALSE))
  expect_identical(as.character(on$spot_clusters), as.character(off$spot_clusters))
  expect_equal(on$prototypes, off$prototypes)
  expect_null(off$spatial_test)
  expect_true(is.numeric(on$spatial_test$p_value))
  ## define_niches passes it through
  dn <- suppressWarnings(define_niches(t(d$P), d$co, resolution_seq = rs, n_perm = 99,
                                       verbose = FALSE))
  expect_equal(dn$spatial_test$p_value, on$spatial_test$p_value)
})

test_that("integrate_slices reports one spatial test per section", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")
  secs <- lapply(1:3, function(i) .sec(200, seed = 10 + i))
  prop <- lapply(secs, `[[`, "P"); coord <- lapply(secs, `[[`, "co")
  sup <- suppressWarnings(integrate_slices(prop, coord, c("a", "b", "c"), resolution = 1,
                                           per_sample_args = list(
                                             resolution_seq = seq(0.3, 1.2, by = 0.3),
                                             n_perm = 49),
                                           verbose = FALSE))
  expect_identical(rownames(sup$spatial_tests), c("a", "b", "c"))
  ## spot_super_niche is named by spot (it came back unnamed before 0.2.2)
  expect_identical(names(sup$spot_super_niche), unlist(lapply(prop, rownames), use.names = FALSE))
  expect_true(all(sup$spatial_tests$p_value < 0.05))
})
