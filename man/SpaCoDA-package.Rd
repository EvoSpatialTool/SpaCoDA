\name{SpaCoDA-package}
\alias{SpaCoDA-package}
\alias{SpaCoDA}
\docType{package}
\title{Fast Spatial Transcriptomics Deconvolution and Niche Discovery}
\description{
Estimates spot-level cell-type proportions from spatial transcriptomics by a
variance-weighted, platform-corrected non-negative least-squares fit against a
single-cell reference (\code{\link{deconvolve_spots}}); an optional
graph-Laplacian smoothness term,
\eqn{\min_{W \ge 0} \sum_i \|y_i - B w_i\|^2_{V_i} + \lambda \mathrm{tr}(W L W^\top)},
is off by default (\eqn{\lambda = 0}). The hot loop is implemented in
RcppArmadillo.

Spatial niches are found per section with \code{\link{define_niches}} /
\code{\link{cluster_single_sample}} (Aitchison distance between spot
compositions fused with the spatial distance, principal coordinates, a
shared-nearest-neighbour graph and the resolution with the maximum modularity)
and integrated across sections into meta-niches with
\code{\link{integrate_slices}}. \code{\link{niche_spatial_test}} tests whether
the composition is spatially organised at all. Seurat wrappers are in
\code{\link{deconvolve_seurat}}.

This version contains the methods used in the accompanying paper.
}
\author{Kai}
\keyword{package}
