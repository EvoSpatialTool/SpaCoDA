\name{fuse_distances}
\alias{fuse_distances}
\title{Fuse a composition distance with a spatial distance}
\description{Both inputs are scaled to a comparable range (by rank, or by their maximum) before the weighted sum.}
\usage{fuse_distances(dist_expr, dist_spatial, use_rank_norm = FALSE,
  gaussian_kernel_distance = FALSE, expr_weight = 0.7, spat_weight = 0.3)}
\arguments{  \item{dist_expr}{composition distance, see \code{\link{compute_aitchison_distance}}}
  \item{dist_spatial}{spatial distance}
  \item{use_rank_norm}{rank-normalise both matrices first (paper: FALSE)}
  \item{gaussian_kernel_distance}{set when \code{dist_spatial} is already a 0-1 kernel}
  \item{expr_weight, spat_weight}{fusion weights (paper: 0.7 / 0.3)}}
\value{The fused distance matrix.}
\seealso{\code{\link{define_niches}}, \code{\link{integrate_slices}}}
\keyword{methods}
