\name{niche_spatial_test}
\alias{niche_spatial_test}
\title{Permutation test: is the composition spatially organised?}
\description{
Tests whether neighbouring spots have more similar cell-type compositions than
random pairs of spots. The statistic is the mean squared Aitchison (CLR)
distance between each spot and its \code{k} nearest spatial neighbours (a
multivariate Geary's C); the null distribution shuffles the compositions across
spots while keeping the coordinates fixed.

The test does not use any clustering. This matters because the fused
composition + spatial distance of \code{\link{cluster_single_sample}} also tiles
a section whose compositions carry no spatial structure into contiguous
"niches"; a non-significant result here says that such niches come from the
spatial distance term, not from the composition. A significant result says the
composition is spatially organised, which includes smooth gradients: it does
not by itself show that the structure consists of discrete niches.
}
\usage{
niche_spatial_test(prop_mat, coords, k = 6, n_perm = 999, pseudo = 1e-06,
                   seed = 42)
}
\arguments{
  \item{prop_mat}{spots x cell types, non-negative proportions (e.g. \code{t(W_cell)}).}
  \item{coords}{spots x 2 coordinates, same row order as \code{prop_mat}.}
  \item{k}{spatial nearest neighbours per spot.}
  \item{n_perm}{number of permutations.}
  \item{pseudo}{pseudocount added before the CLR transform, as in
  \code{\link{compute_aitchison_distance}}.}
  \item{seed}{random seed for the permutations. The caller's random-number
  stream is restored afterwards.}
}
\value{
A list with \code{statistic}, \code{null_mean}, \code{null_sd}, \code{ratio}
(\code{statistic / null_mean}; values below 1 mean neighbours are alike),
\code{z}, \code{p_value} (one-sided, \code{(1 + #null <= observed) / (n_perm + 1)}),
\code{n_perm}, \code{k}, \code{n_spots} and \code{note}. When the test does not
apply (fewer than two effective cell types, or fewer than \code{k + 2} spots) the
numbers are \code{NA} and \code{note} says why.
}
\seealso{\code{\link{cluster_single_sample}}, \code{\link{define_niches}}}
\examples{
set.seed(1)
n <- 300
co <- cbind(runif(n), runif(n))
## composition that changes with x: spatially organised
P <- cbind(co[, 1], 1 - co[, 1], 0.2) + matrix(rexp(3 * n, 10), n, 3)
P <- P / rowSums(P)
niche_spatial_test(P, co, n_perm = 99)$p_value
## the same compositions in random places: not organised
niche_spatial_test(P[sample(n), ], co, n_perm = 99)$p_value
}
\keyword{htest}
