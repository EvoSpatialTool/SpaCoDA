\name{core-deconvolution}
\alias{build_signature}
\alias{select_markers}
\alias{mrna_to_cell}
\alias{build_spatial_graph}
\alias{weights_from_mu}
\alias{estimate_platform}
\alias{deconvolve_spots}
\alias{smooth_by_graph}
\alias{rmse}
\alias{jsd}
\title{Core deconvolution pipeline}
\usage{
build_signature(X, labels, gene_names = NULL, target = 10000)
select_markers(sig, n_per_type = 150, min_spec = 0.15,
               min_expr_quantile = 0.5, restrict_to = NULL)
mrna_to_cell(W, cell_sizes)
build_spatial_graph(coords, k = 6)
weights_from_mu(mu, phi = 0.1, scheme = "nb", max_ratio = 10000)
estimate_platform(Y_mean, B, w_mean, clip = c(0.2, 5), shrink = 0.5)
deconvolve_spots(ST, sig, gene_idx = NULL, coords = NULL, lam = 0,
                 phi = 0.1, weights = "nb", n_outer = 3, n_sweeps = 6,
                 knn = 6, platform_correction = TRUE, tol = 1e-04,
                 engine = c("rcpp", "pureR"), verbose = TRUE,
                 cell_counts = NULL, alpha_scale = 0.1)
smooth_by_graph(W, A, n_iter = 2, alpha = 0.6)
rmse(a, b)
jsd(P, Q, eps = 1e-12)
}
\arguments{
  \item{X}{cells x genes count matrix for the single-cell reference.}
  \item{labels}{cell-type label per cell.}
  \item{gene_names}{optional gene names.}
  \item{target}{library-size normalization target (default 1e4).}
  \item{sig}{signature list returned by \code{build_signature}.}
  \item{n_per_type}{number of marker genes selected per cell type.}
  \item{min_spec}{minimum within-type specificity of a marker.}
  \item{min_expr_quantile}{expression quantile threshold for markers.}
  \item{restrict_to}{optional gene set to restrict markers to.}
  \item{W}{proportion matrix (cell types x spots or spots x cell types).}
  \item{cell_sizes}{mean library size per cell type, from the signature.}
  \item{coords}{spot coordinates (two columns).}
  \item{k}{number of nearest neighbours for the spatial graph.}
  \item{mu}{mean expression used to build variance weights.}
  \item{phi}{overdispersion parameter of the NB variance scheme.}
  \item{scheme}{variance scheme: \code{"nb"}, \code{"poisson"}, \code{"dwls"} or \code{"none"}.}
  \item{max_ratio}{clipping ratio for extreme weights.}
  \item{Y_mean}{mean spatial expression per gene.}
  \item{B}{signature matrix (genes x cell types).}
  \item{w_mean}{mean estimated proportions.}
  \item{clip}{range of the platform-correction factors.}
  \item{shrink}{shrinkage applied to the log correction factors.}
  \item{ST}{spots x genes spatial count matrix (raw counts), columns in the
  order of \code{sig$gene_names}; dense or sparse.}
  \item{gene_idx}{genes used in the fit: integer indices into
  \code{sig$gene_names} or a character vector of gene names. Default: all genes;
  markers from \code{select_markers} are the usual choice.}
  \item{lam}{spatial smoothing strength (graph-Laplacian lambda). The default 0
  turns the smoothness term off, which is the setting used in the paper; on
  simulations with 4, 6 and 15 cells per spot \code{lam = 0.5} lowered the mean
  Pearson correlation with the truth by 0.13, 0.09 and 0.13. \code{lam > 0} needs
  \code{coords}.}
  \item{weights}{variance scheme passed to \code{weights_from_mu}.}
  \item{n_outer}{re-weighting rounds.}
  \item{n_sweeps}{Gauss-Seidel sweeps per round; only matter when \code{lam > 0}.}
  \item{knn}{neighbours of the spatial graph.}
  \item{platform_correction}{estimate a gene-level platform factor
  (\code{estimate_platform}) between the re-weighting rounds.}
  \item{tol}{convergence tolerance (mean absolute change of the composition).}
  \item{engine}{\code{"rcpp"} (default) or \code{"pureR"} reference implementation.}
  \item{verbose}{print progress messages.}
  \item{cell_counts}{optional number of cells per spot (named by spot, or one
  value per spot in order). When given, the cell proportions are shrunk towards
  the global mean composition, more strongly in spots with few cells; the
  unrefined estimate is returned as \code{W_cell_raw}.}
  \item{alpha_scale}{strength of that shrinkage (prior concentration per cell type).}
  \item{A}{spatial adjacency matrix.}
  \item{n_iter, alpha}{smoothing iterations and strength.}
  \item{a, b}{numeric vectors for \code{rmse}.}
  \item{P, Q}{proportion matrices for Jensen-Shannon divergence.}
  \item{eps}{small constant for numerical stability.}
}
\details{
Each spot's counts on the fitted genes are scaled to proportions and fitted by
non-negative least squares, \eqn{\min_{w \ge 0} \sum_g v_g (y_g - \sum_k B_{gk} w_k)^2},
with weights \eqn{v_g} from the negative-binomial mean-variance relation. The
weights and a gene-level platform factor are re-estimated over \code{n_outer}
rounds. The mRNA shares are converted to cell proportions with the reference cell
sizes and the share of each type's mRNA captured by the fitted genes. This is the
method used in the accompanying paper (with \code{lam = 0}).
}
\value{
\code{build_signature} returns a list with the signature matrix \code{B},
\code{cell_sizes} and \code{types}. \code{deconvolve_spots} returns a list with \code{W_cell} (cell
proportions, cell types x spots), \code{W_mrna} (mRNA shares on the fitted genes),
\code{W_cell_raw} (with \code{cell_counts}), \code{f_capture}, \code{types},
\code{runtime} and \code{n_genes_used}. \code{rmse}/\code{jsd} return non-negative numeric scalars.
}
\description{
Core functions of the SpaCoDA pipeline: build a cell-type signature from a
single-cell reference, select markers, correct platform effects, build the kNN
spatial graph, and estimate spot-level cell-type proportions by variance-weighted
non-negative least squares.
}
\examples{
set.seed(1)
X <- matrix(rpois(200 * 20, 2), 200, 20)
sig <- build_signature(X, rep(1:4, each = 50))
mk <- select_markers(sig, n_per_type = 10)
length(mk) > 0
ST <- matrix(rpois(30 * 20, 3), 30, 20)
res <- deconvolve_spots(ST, sig, gene_idx = mk, verbose = FALSE)
dim(res$W_cell)
}
\keyword{models}
