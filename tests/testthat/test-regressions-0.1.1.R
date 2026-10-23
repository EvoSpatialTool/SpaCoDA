
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
# Regression tests for the 0.1.1 fixes
#   1. deconvolve_spots() must label its output and must actually apply the
#      cell-count Dirichlet shrinkage
#   2. define_niches() must label profiles / labels
#   3. integrate_slices() must survive harmony >= 1.2 (project.dim = FALSE) and
#      must return the null-test statistics

.make_sim_labelled <- function(G = 120, K = 3, N = 40, n_ref = 300, seed = 7) {
  set.seed(seed)
  B <- matrix(0.001, G, K)
  for (k in seq_len(K)) B[((k - 1) * 20 + 1):(k * 20), k] <- 0.2
  rest <- setdiff(seq_len(G), seq_len(20 * K))
  B[rest, ] <- rexp(length(rest) * K, 1) * 0.002
  B <- sweep(B, 2, colSums(B), "/")
  colnames(B) <- paste0("type", seq_len(K))
  labels <- rep(paste0("type", seq_len(K)), each = n_ref / K)
  sc <- t(vapply(seq_len(n_ref), function(i) B[, labels[i]] * rgamma(G, 20, 20) * 1e4,
                 numeric(G)))
  colnames(sc) <- paste0("g", seq_len(G))
  true_W <- matrix(rexp(N * K, 1), N, K); true_W <- true_W / rowSums(true_W)
  Y <- true_W %*% t(B) * 1e4 + matrix(rnorm(N * G, 0, 0.05), N, G); Y[Y < 0] <- 0
  rownames(Y) <- paste0("spot", seq_len(N))
  colnames(Y) <- paste0("g", seq_len(G))
  list(sc = sc, labels = labels, Y = Y, true_W = true_W,
       coords = cbind(runif(N, 0, 10), runif(N, 0, 10)), K = K, N = N)
}

test_that("deconvolve_spots labels spots and cell types", {
  sim <- .make_sim_labelled()
  sig <- build_signature(sim$sc, sim$labels)
  res <- deconvolve_spots(sim$Y, sig, gene_idx = select_markers(sig, n_per_type = 10),
                          coords = sim$coords, lam = 0, verbose = FALSE)
  expect_identical(rownames(res$W_cell), sig$types)
  expect_identical(colnames(res$W_cell), rownames(sim$Y))
  expect_identical(dimnames(res$W_mrna), dimnames(res$W_cell))
})

test_that("cell_counts shrinkage is applied and is keyed by spot", {
  sim <- .make_sim_labelled()
  sig <- build_signature(sim$sc, sim$labels)
  mk  <- select_markers(sig, n_per_type = 10)
  base <- deconvolve_spots(sim$Y, sig, gene_idx = mk, coords = sim$coords, lam = 0,
                           verbose = FALSE)
  ncell <- setNames(rep(6, sim$N), rownames(sim$Y))
  shrunk <- deconvolve_spots(sim$Y, sig, gene_idx = mk, coords = sim$coords, lam = 0,
                             cell_counts = ncell, alpha_scale = 0.5, verbose = FALSE)
  expect_false(is.null(shrunk$W_cell_raw))                 # refinement happened
  expect_gt(max(abs(shrunk$W_cell - base$W_cell)), 0)      # and changed the answer
  expect_equal(shrunk$W_cell_raw, base$W_cell, tolerance = 1e-12)
  expect_equal(unname(colSums(shrunk$W_cell)), rep(1, sim$N), tolerance = 1e-8)

  # unnamed vector of the right length is accepted (positional)
  pos <- deconvolve_spots(sim$Y, sig, gene_idx = mk, coords = sim$coords, lam = 0,
                          cell_counts = rep(6, sim$N), alpha_scale = 0.5, verbose = FALSE)
  expect_equal(pos$W_cell, shrunk$W_cell, tolerance = 1e-12)

  # wrong length / unmatched names must fail loudly, not silently do nothing
  expect_error(deconvolve_spots(sim$Y, sig, gene_idx = mk, coords = sim$coords, lam = 0,
                                cell_counts = rep(6, sim$N - 1), verbose = FALSE),
               "one value per spot")
  expect_error(deconvolve_spots(sim$Y, sig, gene_idx = mk, coords = sim$coords, lam = 0,
                                cell_counts = setNames(rep(6, sim$N), paste0("zz", 1:sim$N)),
                                verbose = FALSE),
               "no names")
})

test_that("define_niches labels profiles and labels", {
  sim <- .make_sim_labelled()
  sig <- build_signature(sim$sc, sim$labels)
  res <- deconvolve_spots(sim$Y, sig, gene_idx = select_markers(sig, n_per_type = 10),
                          coords = sim$coords, lam = 0, verbose = FALSE)
  ## the default is the paper's Aitchison pipeline; the k-means behaviour belongs
  ## to the legacy path, which is still tested. Since 0.2.2 the profile rownames
  ## are the label values (they were "niche1..J" before).
  nch <- define_niches(res$W_cell, sim$coords, method = "kmeans", n_niches = 2, seed = 1)
  expect_identical(colnames(nch$profiles), rownames(res$W_cell))
  expect_identical(rownames(nch$profiles), as.character(1:2))
  expect_setequal(rownames(nch$profiles), as.character(unique(nch$labels)))
  expect_identical(names(nch$labels), colnames(res$W_cell))
})

test_that("define_niches warns that the legacy arguments are ignored", {
  skip_if_not_installed("Seurat")
  sim <- .make_sim_labelled()
  sig <- build_signature(sim$sc, sim$labels)
  res <- deconvolve_spots(sim$Y, sig, gene_idx = select_markers(sig, n_per_type = 10),
                          coords = sim$coords, lam = 0, verbose = FALSE)
  ## the legacy k-means arguments are accepted but ignored by the 0.2.0 default
  dn <- .rec_warn(define_niches(res$W_cell, sim$coords, n_niches = 3, seed = 1))
  nch <- dn$value
  expect_true(any(grepl("ignored by the default method", dn$warnings)))
  expect_identical(names(nch$labels), colnames(res$W_cell))
})

test_that("integrate_slices works with harmony >= 1.2 and returns the null test", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")
  sim <- .make_sim_labelled(K = 4, N = 160, n_ref = 400)
  sig <- build_signature(sim$sc, sim$labels)
  res <- deconvolve_spots(sim$Y, sig, gene_idx = select_markers(sig, n_per_type = 10),
                          coords = sim$coords, lam = 0, verbose = FALSE)
  W <- t(res$W_cell)
  half <- rep(c("s1", "s2"), each = nrow(W) / 2)
  ## the CLR path is what 0.1.1 shipped; the default since 0.2.0 is the prototype path
  sup <- suppressWarnings(integrate_slices(prop_list = list(W[half == "s1", , drop = FALSE],
                                           W[half == "s2", , drop = FALSE]),
                          coord_list = list(sim$coords[half == "s1", , drop = FALSE],
                                            sim$coords[half == "s2", , drop = FALSE]),
                          slice_names = c("s1", "s2"), method = "clr",
                          spatial_smooth = 0.5, n_pcs = 3, do_harmony = TRUE,
                          resolution_seq = c(0.2, 0.5), verbose = FALSE))
  expect_true(sup$reduction %in% c("harmony", "pca"))
  expect_true(all(c("z", "best_modularity", "null_modularity") %in% names(sup)))
  expect_true(is.finite(sup$z))
  expect_equal(nrow(sup$niche_summary), nrow(sup$niche_summary_prop))
})
