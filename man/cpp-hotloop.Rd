\name{cpp-hotloop}
\alias{all_grams_cpp}
\alias{gs_sweeps_cpp}
\title{C++ hot-loop functions (RcppArmadillo)}
\usage{
all_grams_cpp(B, V)
gs_sweeps_cpp(Qs_flat, Cs, W, A, deg, lam_eff, n_sweeps, tol, use_spatial)
}
\arguments{
  \item{B}{signature matrix (genes x cell types).}
  \item{V}{variance-weight matrix matching the expression matrix.}
  \item{Qs_flat}{flattened per-spot Gram matrices (N x K^2).}
  \item{Cs}{per-spot linear terms (N x K).}
  \item{W}{warm-start weight matrix (K x N).}
  \item{A}{spatial adjacency (sparse).}
  \item{deg}{node degrees of the spatial graph.}
  \item{lam_eff}{effective spatial regularization strength.}
  \item{n_sweeps}{number of Gauss-Seidel sweeps.}
  \item{tol}{convergence tolerance.}
  \item{use_spatial}{whether the spatial term is active.}
}
\value{
Matrices used or produced by the solver: the stacked Gram matrices, the
optimized weight matrix.
}
\description{
Direct bindings to the RcppArmadillo hot loop. These are exported for
completeness and testing; most users should call \code{deconvolve_spots} instead,
which prepares the inputs.
}
\keyword{internal}
