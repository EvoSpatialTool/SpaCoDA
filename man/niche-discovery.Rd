\name{niche-discovery}
\alias{define_niches}
\alias{clr_transform}
\alias{integrate_slices}
\alias{select_resolution}
\alias{niche_qc}
\title{Spatial niche discovery (paper pipeline)}
\usage{
define_niches(W_cell, coords, method = c("aitchison", "kmeans"),
              n_niches = 6, knn = 6, smooth_iter = 2, alpha = 0.6,
              expr_weight = 0.7, spat_weight = 0.3, k_nn_spatial = 10,
              gaussian_kernel_distance = FALSE, use_rank_norm = FALSE,
              k_nn_snn = 20, resolution_seq = seq(0.3, 3, by = 0.1),
              pcoa_dim = 30, min_effective_ct = 2, seed = 42, verbose = TRUE,
              pcoa_method = c("auto", "partial", "full"),
              spatial_test = TRUE, n_perm = 999)
clr_transform(prop_mat, pseudo = NULL)
integrate_slices(prop_list, coord_list, slice_names,
                 method = c("prototype", "clr"), resolution = 6,
                 resolution_seq = seq(0.3, 3, by = 0.1), pca_dims = 15,
                 do_harmony = TRUE, per_sample_args = list(),
                 spatial_smooth = 0.5, knn_spatial = 6, n_pcs = 20,
                 k_nn_snn = 20, seed = 42, verbose = TRUE)
select_resolution(obj, resolution_seq, graph_name = "Niche_snn",
                  n_null = 10, verbose = TRUE)
niche_qc(res)
}
\arguments{
  \item{W_cell}{spot x cell-type proportion matrix.}
  \item{coords}{spot coordinates.}
  \item{method}{\code{"aitchison"} (default) runs the study pipeline: Aitchison
  distance fused with spatial distance, PCoA, a shared-nearest-neighbour graph and
  the resolution with the maximum modularity. \code{"kmeans"} keeps the pre-0.2.0
  behaviour (spatial smoothing of the proportions followed by k-means) and is what
  uses \code{n_niches}, \code{knn}, \code{smooth_iter} and \code{alpha}.}
  \item{n_niches}{number of niches (\code{method = "kmeans"} only).}
  \item{expr_weight, spat_weight}{weights of the expression and spatial distance in
  the fused distance (study values 0.7 and 0.3).}
  \item{k_nn_spatial}{neighbours used for the spatial distance.}
  \item{gaussian_kernel_distance}{use a Gaussian kernel instead of rank-based
  spatial distances.}
  \item{use_rank_norm}{rank-normalise the expression distance before fusing.}
  \item{pcoa_dim}{PCoA dimensions kept from the fused distance.}
  \item{pcoa_method}{\code{"auto"}, \code{"partial"} or \code{"full"}; see
  \code{\link{cluster_single_sample}}. All three give the same niches;
  \code{"full"} is the study's \code{cmdscale()}.}
  \item{spatial_test, n_perm}{run \code{\link{niche_spatial_test}} with
  \code{n_perm} permutations (default \code{TRUE}, 999); it does not change the
  niches.}
  \item{min_effective_ct}{minimum number of cell types a spot must have to take
  part in the niche clustering without the fallback single-niche result.}
  \item{knn}{neighbourhood size used when defining niches.}
  \item{smooth_iter}{smoothing iterations when defining niches.}
  \item{alpha}{smoothing strength when defining niches.}
  \item{verbose}{print progress messages.}
  \item{prop_mat, prop_list}{proportion matrix, or a list of them per slice.}
  \item{pseudo}{pseudo-count for zeros in the CLR transform (default:
  0.65 x smallest positive entry).}
  \item{coord_list}{list of coordinate matrices per slice.}
  \item{slice_names}{slice names.}
  \item{spatial_smooth}{within-slice smoothing strength (0 = off).}
  \item{knn_spatial}{neighbours for within-slice smoothing.}
  \item{n_pcs}{number of joint PCs.}
  \item{do_harmony}{run Harmony batch correction.}
  \item{k_nn_snn}{k of the shared-nearest-neighbour graph (study value 20). At the
  super-niche stage it is capped at one third of the number of prototypes, which is
  the study value again from 60 prototypes upwards.}
  \item{resolution}{super-niche resolution. The study scripts hard-coded 6, which
  assumes hundreds to thousands of prototypes; \code{NULL} selects it by maximum
  modularity over \code{resolution_seq} instead.}
  \item{pca_dims}{dimensions kept before the super-niche clustering.}
  \item{per_sample_args}{named list of arguments forwarded to the per-section
  clustering, e.g. \code{list(pcoa_dim = 30, k_nn_snn = 20)}.}
  \item{resolution_seq}{clustering resolutions to scan.}
  \item{seed}{random seed.}
  \item{obj}{Seurat object with an SNN graph.}
  \item{graph_name}{name of the graph in the object.}
  \item{n_null}{number of null-graph rewirings for the z-score.}
  \item{res}{result of \code{integrate_slices}.}
}
\description{
Discover spatial niches from cell-type proportions. The default
\code{method = "aitchison"} of \code{define_niches} is the pipeline of the paper:
Aitchison distance between spot compositions fused with spatial distance
(0.7/0.3), PCoA, a shared-nearest-neighbour graph, and the clustering resolution
with the maximum modularity. \code{integrate_slices(method = "prototype")} then
clusters the niche prototypes of all sections to give the cross-sample
super-niches, again at the maximum-modularity resolution. The pre-0.2.0
behaviour remains reachable: \code{method = "kmeans"} for the per-section niches
and \code{method = "clr"} for the CLR / joint-PCA / Harmony integration with a
modularity null test.
}
\value{
\code{define_niches} returns the niche label of every spot, the per-niche
cell-type profiles, the scan over \code{resolution_seq}, the chosen resolution
and, for the default method, \code{spatial_test}. \code{integrate_slices}
(prototype path) adds \code{spatial_tests}, one row per section.
\code{integrate_slices} returns the joint Seurat object, the embeddings, the chosen
resolution (\code{best_resolution}), the number of super-niches (\code{best_k}),
their modularity (\code{best_modularity}), the super-niche label of every spot
(\code{spot_super_niche}), the prototypes and the per-slice composition tables. The
prototype path performs no null test, so \code{z} and \code{null_modularity} are
\code{NA} there; \code{select_resolution} is the variant that adds the null test.
\code{clr_transform} returns the CLR-transformed matrix. \code{niche_qc} returns the dominance and the
slice-mixing entropy of every super-niche.
}
\keyword{cluster}
