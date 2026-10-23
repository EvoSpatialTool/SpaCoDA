// SpaCoDA C++ core (RcppArmadillo)
//
// Only the hot loop lives here; everything else is plain R.
//
// What is hot:
//   1. Building the Gram matrix Q_i = B' diag(v_i) B for every spot
//   2. The Gauss-Seidel sweeps: N exact K x K non-negative QPs per sweep
//
//   min_{w>=0} ||y_i - B w_i||^2_{V_i} + lambda * tr(W L W')
// solved by exact block coordinate descent, each block an NNLS-Cholesky solve.
// With lambda = 0 (the default) the blocks are independent and one sweep is an
// exact per-spot weighted NNLS.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ---------------------------------------------------------------------------
// Lawson-Hanson NNLS (active set). Solves min_{x>=0} ||A x - b||^2.
// A is K x K here (tiny), so this is microseconds.
// ---------------------------------------------------------------------------
static vec nnls_lh(const mat& A, const vec& b, int maxit = 100,
                   double tol = 1e-10) {
  const uword n = A.n_cols;
  vec x(n, fill::zeros);
  uvec P(n, fill::zeros);              // passive (free) set indicator
  vec w = A.t() * (b - A * x);

  int iter = 0;
  while (any(P == 0) && any(w.elem(find(P == 0)) > tol) && iter < maxit * (int)n) {
    // pick index with max gradient among the active (zeroed) set
    uvec act = find(P == 0);
    uword j_rel = index_max(w.elem(act));
    uword j = act(j_rel);
    P(j) = 1;

    vec s(n, fill::zeros);
    uvec pas = find(P == 1);
    mat Ap = A.cols(pas);
    vec sp;
    bool ok = solve(sp, Ap, b, solve_opts::no_approx + solve_opts::fast);
    if (!ok) { sp = pinv(Ap) * b; }
    s.elem(pas) = sp;

    int inner = 0;
    while (pas.n_elem > 0 && s.elem(pas).min() <= 0 && inner < maxit) {
      // move as far as we can toward s while staying feasible
      double alpha = datum::inf;
      for (uword t = 0; t < pas.n_elem; ++t) {
        uword idx = pas(t);
        if (s(idx) <= 0) {
          double denom = x(idx) - s(idx);
          if (denom > tol) alpha = std::min(alpha, x(idx) / denom);
        }
      }
      if (!std::isfinite(alpha)) break;
      x = x + alpha * (s - x);
      // drop variables that hit zero
      for (uword t = 0; t < n; ++t) {
        if (P(t) == 1 && std::abs(x(t)) < tol) { P(t) = 0; x(t) = 0.0; }
      }
      pas = find(P == 1);
      s.zeros();
      if (pas.n_elem > 0) {
        mat Ap2 = A.cols(pas);
        vec sp2;
        bool ok2 = solve(sp2, Ap2, b, solve_opts::no_approx + solve_opts::fast);
        if (!ok2) { sp2 = pinv(Ap2) * b; }
        s.elem(pas) = sp2;
      }
      inner++;
    }
    x = s;
    w = A.t() * (b - A * x);
    iter++;
  }
  x.elem(find(x < 0)).zeros();
  return x;
}

// ---------------------------------------------------------------------------
// Exact non-negative QP:  min_{w>=0} w'Qw - 2w'c
// via Cholesky Q = R'R  =>  ||Rw - z||^2 with R'z = c
// ---------------------------------------------------------------------------
static vec nnqp_cpp(mat Q, const vec& c, double ridge = 0.0) {
  const uword K = Q.n_rows;
  Q = 0.5 * (Q + Q.t());
  if (ridge > 0) Q.diag() += ridge;

  double sc = std::max(trace(Q) / (double)K, 1e-12);
  mat R;
  bool ok = false;
  for (int j = 0; j < 8; ++j) {
    mat Qj = Q;
    if (j > 0) Qj.diag() += 1e-8 * std::pow(10.0, (double)j) * sc;
    ok = chol(R, Qj);          // upper: Q = R' R
    if (ok) break;
  }
  if (!ok) {
    mat Qj = Q; Qj.diag() += 1e-3 * sc;
    chol(R, Qj);
  }
  vec z = solve(trimatl(R.t()), c);   // R' z = c
  return nnls_lh(R, z);
}

// ---------------------------------------------------------------------------
// Batched Gram matrices for ALL spots, in ONE BLAS call.
//
//   Q_i[k,l] = sum_g v_ig B_gk B_gl = <v_i, B_:k .* B_:l>
//
// so stacking the K^2 elementwise column products of B into P (G x K^2) turns
// "all Gram matrices" into a single matrix product V' P (one BLAS call).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat all_grams_cpp(const arma::mat& B, const arma::mat& V) {
  const uword G = B.n_rows, K = B.n_cols;
  mat P(G, K * K);
  for (uword k = 0; k < K; ++k)
    for (uword l = 0; l < K; ++l)
      P.col(k * K + l) = B.col(k) % B.col(l);
  return V.t() * P;                    // (N, K^2)
}

// ---------------------------------------------------------------------------
// Gauss-Seidel sweeps with the spatial Laplacian prior.
//
// Local problem at spot i:
//   Q_i = B'V_i B + lam*d_i*I ,  c_i = B'V_i y_i + lam * sum_j A_ij w_j
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat gs_sweeps_cpp(const arma::mat& Qs_flat,   // (N, K^2)
                        const arma::mat& Cs,        // (N, K)
                        arma::mat W,                // (K, N) warm start
                        const arma::sp_mat& A,      // (N, N) adjacency
                        const arma::vec& deg,       // (N)
                        double lam_eff,
                        int n_sweeps,
                        double tol,
                        bool use_spatial) {
  const uword N = Cs.n_rows, K = Cs.n_cols;

  // Pre-factorise: Q_i + lam*d_i*I is constant across sweeps
  std::vector<mat> Rs(N);
  for (uword i = 0; i < N; ++i) {
    mat Q(K, K);
    for (uword k = 0; k < K; ++k)
      for (uword l = 0; l < K; ++l)
        Q(k, l) = Qs_flat(i, k * K + l);
    Q = 0.5 * (Q + Q.t());
    if (use_spatial) Q.diag() += lam_eff * deg(i);

    double sc = std::max(trace(Q) / (double)K, 1e-12);
    mat R; bool ok = false;
    for (int j = 0; j < 8; ++j) {
      mat Qj = Q;
      if (j > 0) Qj.diag() += 1e-8 * std::pow(10.0, (double)j) * sc;
      ok = chol(R, Qj);
      if (ok) break;
    }
    if (!ok) { mat Qj = Q; Qj.diag() += 1e-3 * sc; chol(R, Qj); }
    Rs[i] = R;
  }

  for (int s = 0; s < n_sweeps; ++s) {
    mat W_old = W;
    for (uword i = 0; i < N; ++i) {
      vec c = Cs.row(i).t();
      if (use_spatial) {
        vec m(K, fill::zeros);
        for (sp_mat::const_col_iterator it = A.begin_col(i);
             it != A.end_col(i); ++it) {
          m += (*it) * W.col(it.row());
        }
        c += lam_eff * m;
      }
      vec z = solve(trimatl(Rs[i].t()), c);
      W.col(i) = nnls_lh(Rs[i], z);
    }
    if (mean(mean(abs(W - W_old))) < tol) break;
  }
  return W;
}

// ---------------------------------------------------------------------------
// M independent weighted non-negative ridge regressions (used by the
// count model, deconvolve_spots(model = "nb")).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat wnnls_cpp(const arma::mat& B,   // (G, J)
                    const arma::mat& Y,   // (G, M)
                    const arma::mat& V,   // (G, M)
                    double ridge) {
  const uword J = B.n_cols, M = Y.n_cols;
  mat W(J, M, fill::zeros);
  for (uword m = 0; m < M; ++m) {
    mat Bv = B.each_col() % V.col(m);
    mat Q = B.t() * Bv;
    vec c = Bv.t() * Y.col(m);
    double r = ridge * trace(Q) / (double)J;
    W.col(m) = nnqp_cpp(Q, c, r);
  }
  return W;
}
