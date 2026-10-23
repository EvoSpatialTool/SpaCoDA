
## Record warnings instead of letting them escape: these toy inputs are small
## enough that Seurat's own SNN builder / GroupSingletons chatter as well, and the
## suite should stay clean while still asserting SpaCoDA's own warnings.
.rec_warn <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}
# Regression tests for the 0.2.0 niche port
#   1. define_niches(method = "aitchison") is the paper pipeline and equals
#      cluster_single_sample() on the same input
#   2. the legacy k-means path still works (method = "kmeans")
#   3. integrate_slices() honours the `resolution` argument (it used to be
#      hard-coded to 6 in the study scripts) and resolution = NULL selects one
#   4. degenerate sections become a single niche instead of an error
#   5. method = "clr" still runs the pre-0.2.0 integration
#
# The full numerical equivalence with the study implementation is established by
# tests/verify_port_vs_original_20261006_01.r on real data (identical labels and
# prototypes, max |difference| = 0).

.make_section <- function(n_spots = 240, n_types = 6, seed = 1, homogeneous = FALSE,
                          grid = 16) {
  set.seed(seed)
  if (homogeneous) {
    P <- matrix(0, n_spots, n_types,
                dimnames = list(paste0("s", seq_len(n_spots)), paste0("t", seq_len(n_types))))
    P[, 1] <- 1
  } else {
    side <- ceiling(sqrt(n_spots))
    gx <- ((seq_len(n_spots) - 1) %% side) + 1
    gy <- ((seq_len(n_spots) - 1) %/% side) + 1
    block <- ((gx > side / 2) + 1) + 2 * (gy > side / 2)          # 4 spatial blocks
    base <- matrix(rexp(n_spots * n_types), n_spots, n_types)
    base[cbind(seq_len(n_spots), ((block - 1) %% n_types) + 1)] <- 8   # block-specific type
    P <- base / rowSums(base)
    rownames(P) <- paste0("s", seq_len(n_spots))
    colnames(P) <- paste0("t", seq_len(n_types))
  }
  side <- ceiling(sqrt(n_spots))
  co <- cbind(x = ((seq_len(n_spots) - 1) %% side) * 100,
              y = ((seq_len(n_spots) - 1) %/% side) * 100)
  rownames(co) <- rownames(P)
  list(P = P, co = co)
}

test_that("cluster_single_sample finds spatially coherent niches", {
  skip_if_not_installed("Seurat")
  sec <- .make_section()
  res <- cluster_single_sample(sec$P, sec$co, resolution_seq = seq(0.3, 1.5, by = 0.2),
                               seed = 1, verbose = FALSE)
  expect_false(res$homogeneous)
  expect_equal(length(res$spot_clusters), nrow(sec$P))
  expect_identical(names(res$spot_clusters), rownames(sec$P))
  expect_equal(nrow(res$prototypes), length(unique(res$spot_clusters)))
  expect_equal(ncol(res$prototypes), ncol(sec$P))
  expect_true(all(c("best_resolution", "best_k", "scan_results", "dist_fused") %in% names(res)))
})

test_that("define_niches defaults to the paper pipeline and matches it exactly", {
  skip_if_not_installed("Seurat")
  sec <- .make_section(seed = 2)
  W <- t(sec$P)                                   # cell types x spots
  a <- define_niches(W, sec$co, resolution_seq = seq(0.3, 1.5, by = 0.2), seed = 1,
                     verbose = FALSE)
  b <- cluster_single_sample(sec$P, sec$co, resolution_seq = seq(0.3, 1.5, by = 0.2),
                             seed = 1, verbose = FALSE)
  expect_identical(as.character(a$labels), as.character(b$spot_clusters))
  expect_equal(unname(as.matrix(a$profiles)), unname(as.matrix(b$prototypes)))
  expect_null(a$smoothed)                         # no spatial smoothing in this path
  expect_equal(ncol(a$profiles), nrow(W))
})

test_that("the legacy k-means path is still available", {
  sec <- .make_section(seed = 3)
  W <- t(sec$P)
  k <- define_niches(W, sec$co, method = "kmeans", n_niches = 3, seed = 0, verbose = FALSE)
  expect_equal(length(unique(k$labels)), 3)
  expect_equal(nrow(k$profiles), 3)
  expect_false(is.null(k$smoothed))
})

test_that("degenerate sections become a single niche, not an error", {
  hom <- .make_section(homogeneous = TRUE)
  r <- cluster_single_sample(hom$P, hom$co, verbose = FALSE)
  expect_true(r$homogeneous)
  expect_equal(length(unique(r$spot_clusters)), 1L)
  expect_equal(nrow(r$prototypes), 1L)

  tiny <- .make_section(n_spots = 6)
  r2 <- cluster_single_sample(tiny$P, tiny$co, verbose = FALSE)
  expect_true(r2$homogeneous)
})

test_that("integrate_slices honours the super-niche resolution argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")
  secs <- lapply(1:4, function(i) .make_section(seed = 10 + i))
  prop <- lapply(secs, `[[`, "P"); coord <- lapply(secs, `[[`, "co")
  names(prop) <- names(coord) <- paste0("sec", 1:4)
  ps <- list(resolution_seq = seq(0.3, 1.5, by = 0.2))

  a <- suppressWarnings(integrate_slices(prop, coord, names(prop), resolution = 1,
                                         per_sample_args = ps, verbose = FALSE))
  b <- suppressWarnings(integrate_slices(prop, coord, names(prop), resolution = 2,
                                         per_sample_args = ps, verbose = FALSE))
  ## the value the study scripts hard-coded: it must now be honoured as an argument
  ## These toy sections are small enough that Seurat's own GroupSingletons
  ## (max(connectivity) with no SNN edge) and the SNN builder warn; the behaviour
  ## under test is SpaCoDA's resolution handling, not Seurat's chatter.
  c6 <- suppressWarnings(integrate_slices(prop, coord, names(prop), resolution = 6,
                                          per_sample_args = ps, verbose = FALSE))
  auto <- suppressWarnings(integrate_slices(prop, coord, names(prop), resolution = NULL,
                                            resolution_seq = seq(0.3, 1.2, by = 0.3),
                                            per_sample_args = ps, verbose = FALSE))

  expect_equal(a$best_resolution, 1)
  expect_equal(b$best_resolution, 2)
  expect_equal(c6$best_resolution, 6)                    # was hard-coded before 0.2.0
  expect_equal(a$scan$resolution, 1)
  expect_equal(c6$scan$resolution, 6)
  ## more resolution can only ever give at least as many super-niches
  expect_gte(length(unique(b$super_niche)), length(unique(a$super_niche)))
  expect_true(is.numeric(auto$best_resolution))          # NULL -> selected by modularity
  expect_gt(length(auto$super_niche), 0)

  ## every spot gets a super-niche and the tables are consistent
  expect_equal(length(a$spot_super_niche), sum(vapply(prop, nrow, integer(1))))
  expect_false(any(is.na(a$spot_super_niche)))
  expect_equal(nrow(a$niche_summary_prop), length(unique(a$super_niche)))
  expect_equal(ncol(a$niche_summary_prop), length(prop))
})

test_that("a resolution too high for the prototype set is reported, not silently returned", {
  skip_if_not_installed("Seurat")
  ## with few prototypes a high resolution puts every prototype in its own
  ## community; Seurat then labels them all "singleton" and modularity is 0
  secs <- lapply(1:2, function(i) .make_section(seed = 30 + i))
  prop <- lapply(secs, `[[`, "P"); coord <- lapply(secs, `[[`, "co")
  out <- .rec_warn(integrate_slices(prop, coord, c("a", "b"), resolution = 20,
                                    per_sample_args = list(
                                      resolution_seq = seq(0.3, 0.6, by = 0.3)),
                                    verbose = FALSE))
  expect_true(any(grepl("singletons: the resolution is too high", out$warnings)))
  expect_true(all(as.character(out$value$super_niche) == "singleton"))
})

test_that("the pre-0.2.0 CLR integration is still reachable", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")
  secs <- lapply(1:3, function(i) .make_section(seed = 20 + i))
  prop <- lapply(secs, `[[`, "P"); coord <- lapply(secs, `[[`, "co")
  r <- suppressWarnings(integrate_slices(prop, coord, paste0("s", 1:3), method = "clr",
                                         n_pcs = 5, resolution_seq = c(0.2, 0.4),
                                         verbose = FALSE))
  expect_true(r$reduction %in% c("pca", "harmony"))
  expect_true(all(c("niche_summary", "common_niches", "clr") %in% names(r)))
  expect_equal(length(r$object$seurat_clusters), sum(vapply(prop, nrow, integer(1))))
})
