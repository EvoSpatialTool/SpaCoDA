\name{select_resolution_by_modularity}
\alias{select_resolution_by_modularity}
\title{Pick a clustering resolution by maximum modularity}
\description{This is the resolution criterion used in the paper. Unlike \code{\link{select_resolution}} it has no null test; that variant compares modularity against degree-preserving rewired graphs and is used by the CLR path of \code{\link{integrate_slices}}.}
\usage{select_resolution_by_modularity(obj, resolution_seq = seq(0.3, 3, by = 0.1))}
\arguments{  \item{obj}{Seurat object with an SNN graph, e.g. from \code{\link{cluster_single_sample}}}
  \item{resolution_seq}{candidate resolutions (paper: 0.3 to 3 by 0.1)}}
\value{A list with \code{best_resolution}, \code{best_k}, \code{best_modularity} and \code{scan_results}.}
\seealso{\code{\link{define_niches}}, \code{\link{integrate_slices}}}
\keyword{methods}
