\name{seurat-wrappers}
\alias{get_counts}
\alias{signature_from_seurat}
\alias{spatial_from_seurat}
\alias{align_genes}
\alias{deconvolve_seurat}
\title{Seurat wrappers}
\usage{
get_counts(obj, assay = NULL, layer = "counts")
signature_from_seurat(sc_obj, celltype_col, assay = NULL, min_cells = 25,
                      verbose = TRUE)
spatial_from_seurat(st_obj, assay = NULL, image = NULL)
align_genes(sig, target_counts, target_genes, verbose = TRUE)
deconvolve_seurat(sc_obj, st_obj, celltype_col, sc_assay = NULL,
                  st_assay = NULL, n_per_type = 150, lam = NULL,
                  add_to_object = TRUE, verbose = TRUE, assay_name = "SpaCoDA")
}
\arguments{
  \item{obj, sc_obj, st_obj}{Seurat objects (single-cell reference and spatial).}
  \item{assay}{assay name (default: active assay).}
  \item{layer}{data layer (default \code{"counts"}).}
  \item{celltype_col}{metadata column with cell-type labels.}
  \item{image}{spatial image name.}
  \item{sig}{signature list.}
  \item{target_counts, target_genes}{matrix/gene set to align the signature to.}
  \item{verbose}{print progress messages.}
  \item{min_cells}{drop cell types with fewer than this many cells when building the signature.}
  \item{sc_assay, st_assay}{assays of the single-cell and spatial objects (default: active assay).}
  \item{n_per_type}{marker genes per cell type.}
  \item{lam}{spatial smoothing strength. \code{NULL} (default) means 0, i.e. no
  spatial smoothing, the setting used in the paper; the median UMI per spot is
  printed as a diagnostic only.}
  \item{add_to_object}{store the result as a new assay in the spatial object.}
  \item{assay_name}{name of that assay (default \code{"SpaCoDA"}; it was
  \code{"spadecon"} before 0.2.2).}
}
\value{
\code{deconvolve_seurat} returns the \code{deconvolve_spots} result (with
\code{W_cell} as cell types x spots) plus \code{spots}, \code{coords} and, when
\code{add_to_object = TRUE}, \code{object}: the spatial Seurat object with an
added assay (\code{assay_name}) holding the proportions.
Other functions return matrices or signature lists used by the core
pipeline.
}
\description{
Convenience wrappers that run the SpaCoDA pipeline directly on Seurat
objects: extract counts, build the signature from an scRNA-seq object, align
genes, run deconvolution, and store the result back as a new assay.
}
\keyword{manip}
