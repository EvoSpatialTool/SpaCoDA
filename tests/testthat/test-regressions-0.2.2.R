# Regression tests for the 0.2.2 fixes: the Seurat wrappers stay sparse, and
# sparse input gives the dense answer

.make_sim_sparse <- function(G = 150, K = 3, N = 50, n_ref = 240, seed = 11) {
  set.seed(seed)
  B <- matrix(rgamma(G * K, 0.3, 1), G, K)
  for (k in seq_len(K)) B[((k - 1) * 20 + 1):(k * 20), k] <- B[((k - 1) * 20 + 1):(k * 20), k] + 5
  B <- sweep(B, 2, colSums(B), "/")
  labels <- rep(paste0("type", seq_len(K)), each = n_ref / K)
  sc <- t(vapply(seq_len(n_ref), function(i)
    as.numeric(rpois(G, B[, match(labels[i], paste0("type", seq_len(K)))] * 800)), numeric(G)))
  colnames(sc) <- paste0("g", seq_len(G)); rownames(sc) <- paste0("cell", seq_len(n_ref))
  W <- matrix(rexp(N * K), N, K); W <- W / rowSums(W)
  Y <- matrix(as.numeric(rpois(N * G, W %*% t(B) * 600)), N, G,
              dimnames = list(paste0("spot", seq_len(N)), colnames(sc)))
  list(sc = sc, labels = labels, Y = Y, coords = cbind(runif(N), runif(N)), K = K)
}

test_that("sparse input gives the dense answer", {
  sim <- .make_sim_sparse()
  sc_sp <- Matrix::Matrix(sim$sc, sparse = TRUE)
  Y_sp  <- Matrix::Matrix(sim$Y, sparse = TRUE)
  sig_d <- build_signature(sim$sc, sim$labels, colnames(sim$sc))
  sig_s <- build_signature(sc_sp, sim$labels, colnames(sim$sc))
  expect_equal(sig_s$B, sig_d$B, tolerance = 1e-12)
  expect_equal(sig_s$cell_sizes, sig_d$cell_sizes, tolerance = 1e-12)
  mk <- select_markers(sig_d, 20)
  rd <- deconvolve_spots(sim$Y, sig_d, gene_idx = mk, verbose = FALSE)
  rs <- deconvolve_spots(Y_sp, sig_s, gene_idx = mk, verbose = FALSE)
  expect_equal(rs$W_cell, rd$W_cell, tolerance = 1e-10)
  co <- sim$coords
  rd2 <- deconvolve_spots(sim$Y, sig_d, gene_idx = mk, coords = co, lam = 0.5, verbose = FALSE)
  rs2 <- deconvolve_spots(Y_sp, sig_s, gene_idx = mk, coords = co, lam = 0.5, verbose = FALSE)
  expect_equal(rs2$W_cell, rd2$W_cell, tolerance = 1e-10)
})

test_that("signature_from_seurat stays sparse and checks the label column", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  sim <- .make_sim_sparse()
  counts <- Matrix::Matrix(t(sim$sc), sparse = TRUE)              # genes x cells
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  obj$celltype <- sim$labels
  sig <- signature_from_seurat(obj, "celltype", min_cells = 1, verbose = FALSE)
  ref <- build_signature(sim$sc, sim$labels, colnames(sim$sc))
  expect_equal(sig$B[rownames(ref$B), ], ref$B, tolerance = 1e-12)
  expect_error(signature_from_seurat(obj, "no_such_column", verbose = FALSE),
               "not found in meta.data")
})
