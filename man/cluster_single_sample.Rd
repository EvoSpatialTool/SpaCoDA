\name{cluster_single_sample}
\alias{cluster_single_sample}
\title{Find the niches of one spatial section (Aitchison/PCoA pipeline)}
\description{This is the niche-finding method of the accompanying paper: Aitchison distance on the spot composition, fused with the spatial distance, PCoA, a shared-nearest-neighbour graph, and the resolution with maximum modularity. A homogeneous section (fewer than two effective cell types, fewer than ten spots, or a flat fused distance) is returned as a single niche rather than an error.}
\usage{cluster_single_sample(prop_mat, coords_mat, expr_weight = 0.7, spat_weight = 0.3,
  k_nn_spatial = 10, gaussian_kernel_distance = FALSE, use_rank_norm = FALSE,
  k_nn_snn = 20, resolution_seq = seq(0.3, 3, by = 0.1), pcoa_dim = 30,
  seed = 42, min_effective_ct = 2, verbose = TRUE,
  pcoa_method = c("auto", "partial", "full"), spatial_test = TRUE, n_perm = 999)}
\arguments{  \item{prop_mat}{spots x cell types composition; rows must be named}
  \item{coords_mat}{spots x 2 coordinates, same row order as \code{prop_mat}}
  \item{expr_weight, spat_weight}{weights when fusing the composition and spatial distances (paper: 0.7 / 0.3)}
  \item{k_nn_spatial}{neighbours used for the adaptive Gaussian bandwidth}
  \item{gaussian_kernel_distance}{use a Gaussian-kernel spatial distance instead of plain Euclidean (paper: FALSE)}
  \item{use_rank_norm}{rank-normalise both distances before fusing (paper: FALSE)}
  \item{k_nn_snn}{neighbours of the shared-nearest-neighbour graph (paper: 20)}
  \item{resolution_seq}{resolutions scanned when choosing the number of niches (paper: 0.3 to 3 by 0.1)}
  \item{pcoa_dim}{PCoA dimensions requested (paper: 30)}
  \item{seed}{random seed (paper: 42)}
  \item{min_effective_ct}{sections with fewer effective cell types are returned as one niche}
  \item{verbose}{print progress}
  \item{pcoa_method}{how the PCoA is computed. \code{"full"} is \code{cmdscale()}, as in the study. \code{"partial"} computes only the \code{pcoa_dim} eigenpairs that are used (RSpectra Lanczos solver), which gives the same niches much faster on large sections, and falls back to \code{cmdscale()} if the solver does not converge. \code{"auto"} (default) uses \code{"partial"} from 1,000 spots when RSpectra is installed.}
  \item{spatial_test}{run \code{\link{niche_spatial_test}} on the composition (default \code{TRUE}). It runs after the clustering, on its own random-number stream, and does not change the niches.}
  \item{n_perm}{permutations of the spatial test.}}
\value{A list with \code{spot_clusters}, \code{prototypes} (niches x cell types), \code{seurat_obj}, \code{best_resolution}, \code{best_k}, \code{best_modularity}, \code{scan_results}, \code{dist_fused}, \code{n_spots}, \code{homogeneous} and \code{spatial_test} (the result of \code{\link{niche_spatial_test}}, or \code{NULL} when \code{spatial_test = FALSE}).}
\seealso{\code{\link{define_niches}}, \code{\link{integrate_slices}}}
\keyword{methods}
