\name{compute_aitchison_distance}
\alias{compute_aitchison_distance}
\title{Aitchison distance between spot compositions}
\description{Zeros are handled with a multiplicative pseudocount and no spot is ever dropped, unlike \code{zCompositions::cmultRepl}, which removes rows and columns dominated by zeros and would silently change the spot set.}
\usage{compute_aitchison_distance(prop_mat, pseudo = 1e-6)}
\arguments{  \item{prop_mat}{spots x cell types, non-negative proportions}
  \item{pseudo}{pseudocount added before the CLR transform}}
\value{An n x n distance matrix. Fewer than two informative cell types gives an all-zero matrix (a homogeneous section).}
\seealso{\code{\link{define_niches}}, \code{\link{integrate_slices}}}
\keyword{methods}
