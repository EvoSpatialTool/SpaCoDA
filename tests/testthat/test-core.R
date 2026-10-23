# 核心功能测试: 签名 / 空间图 / 解卷积精度

.make_sim <- function(G = 200, K = 4, N = 80, n_ref = 400, seed = 42) {
  set.seed(seed)
  B <- matrix(0.001, G, K)
  for (k in 1:K) B[((k-1)*30+1):(k*30), k] <- 0.2
  rest <- setdiff(seq_len(G), seq_len(30*K))
  B[rest, ] <- rexp(length(rest)*K, 1) * 0.002
  B <- sweep(B, 2, colSums(B), "/")
  ref_labels <- rep(1:K, each = n_ref/K)
  sc <- t(vapply(seq_len(n_ref), function(i) B[, ref_labels[i]] * rgamma(G, 20, 20) * 1e4,
                 numeric(G)))
  colnames(sc) <- paste0("g", seq_len(G))
  true_W <- matrix(rexp(N*K, 1), N, K); true_W <- true_W / rowSums(true_W)
  Y <- true_W %*% t(B) * 1e4 + matrix(rnorm(N*G, 0, 0.05), N, G); Y[Y < 0] <- 0
  colnames(Y) <- paste0("g", seq_len(G))
  coords <- cbind(runif(N, 0, 10), runif(N, 0, 10))
  list(sc = sc, ref_labels = ref_labels, Y = Y, true_W = true_W, coords = coords, K = K)
}

test_that("build_signature returns normalized signature", {
  sim <- .make_sim()
  sig <- build_signature(sim$sc, sim$ref_labels)
  expect_equal(ncol(sig$B), sim$K)
  expect_equal(unname(colSums(sig$B)), rep(1, sim$K), tolerance = 1e-8)
  expect_equal(length(sig$cell_sizes), sim$K)
})

test_that("build_spatial_graph is symmetric and sparse", {
  coords <- cbind(runif(50, 0, 10), runif(50, 0, 10))
  A <- build_spatial_graph(coords, k = 6)
  expect_true(inherits(A, "sparseMatrix"))
  expect_equal(sum(abs(A - Matrix::t(A)) > 1e-12), 0)
  expect_equal(sum(Matrix::diag(A) != 0), 0)
})

test_that("deconvolve_spots recovers proportions on simulated data", {
  sim <- .make_sim()
  sig <- build_signature(sim$sc, sim$ref_labels)
  mk  <- select_markers(sig, n_per_type = 25)
  res <- deconvolve_spots(sim$Y, sig, gene_idx = mk, coords = sim$coords, n_sweeps = 10)
  Wt <- t(res$W_cell)
  expect_true(all(Wt >= 0))
  expect_equal(rowSums(Wt), rep(1, nrow(Wt)), tolerance = 1e-6)
  cors <- sapply(seq_len(sim$K), function(k) cor(Wt[, k], sim$true_W[, k]))
  expect_gt(mean(cors), 0.85)
})

test_that("metrics behave", {
  expect_equal(rmse(c(1,2,3), c(1,2,3)), 0)
  expect_gt(jsd(cbind(c(0.5, 0.5)), cbind(c(0.9, 0.1))), 0)
  expect_equal(jsd(cbind(c(0.5, 0.5)), cbind(c(0.5, 0.5))), 0, tolerance = 1e-10)
})
