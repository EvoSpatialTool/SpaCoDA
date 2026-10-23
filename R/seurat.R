# SpaCoDA: Seurat / real-data adapter
# --------------------------------
# Nothing here is exotic: Seurat objects are just count matrices plus metadata,
# and every function in SpaCoDA takes plain matrices. The value of this file is
# that it handles the things that silently corrupt results: gene-axis alignment,
# counts-vs-data slot confusion, Seurat v4/v5 differences, and coordinate
# extraction.

# ---------------------------------------------------------------------------
# Counts extraction that survives Seurat v4 and v5
# ---------------------------------------------------------------------------
get_counts <- function(obj, assay = NULL, layer = "counts") {
  assay <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  # v5 uses layers; v4 uses slots. Try both.
  # Static-call-safe: LayerData only exists in SeuratObject >= 5.
  layer_fn <- if ("LayerData" %in% getNamespaceExports("SeuratObject"))
    getExportedValue("SeuratObject", "LayerData") else NULL
  m <- tryCatch(
    if (!is.null(layer_fn)) layer_fn(obj, assay = assay, layer = layer)
    else stop("LayerData not available"),
    error = function(e) tryCatch(
      Seurat::GetAssayData(obj, assay = assay, slot = layer),
      error = function(e2) stop("Could not extract '", layer, "' from assay '",
                                assay, "': ", conditionMessage(e2))
    )
  )
  if (is.null(m) || nrow(m) == 0) stop("Empty count matrix from assay ", assay)
  m
}

# ---------------------------------------------------------------------------
# Build a SpaCoDA signature straight from a Seurat scRNA object
# ---------------------------------------------------------------------------
#' @param sc_obj      Seurat object, scRNA-seq reference
#' @param celltype_col metadata column holding cell type labels
#' @param assay       assay to use (default: DefaultAssay)
#' @param min_cells   drop cell types with fewer than this many cells. Rare
#'                    types with a handful of cells give noisy signatures and
#'                    are the usual cause of "why is this type everywhere".
signature_from_seurat <- function(sc_obj, celltype_col, assay = NULL,
                                  min_cells = 25, verbose = TRUE) {
  X <- get_counts(sc_obj, assay = assay, layer = "counts")
  # 0.2.2: test the column before as.character(): as.character(NULL) is
  # character(0), not NULL, so a misspelt column used to slip through.
  if (!celltype_col %in% colnames(sc_obj@meta.data))
    stop("Column '", celltype_col, "' not found in meta.data")
  labels <- as.character(sc_obj@meta.data[[celltype_col]])

  keep <- !is.na(labels)
  tab <- table(labels[keep])
  drop <- names(tab)[tab < min_cells]
  if (length(drop) && verbose) {
    message("Dropping ", length(drop), " cell type(s) with < ", min_cells,
            " cells: ", paste(drop, collapse = ", "))
  }
  keep <- keep & !(labels %in% drop)

  X <- X[, keep, drop = FALSE]
  labels <- labels[keep]

  # build_signature() wants cells x genes. 0.2.2: kept sparse -- as.matrix() here made a
  # dense cells x genes copy of the whole reference (~12 GB at 50k x 30k), and
  # build_signature() works on a sparse matrix as it is.
  Xt <- Matrix::t(X)
  if (verbose) message("Signature from ", nrow(Xt), " cells, ", ncol(Xt),
                       " genes, ", length(unique(labels)), " cell types")
  build_signature(Xt, labels, gene_names = rownames(X))
}

# ---------------------------------------------------------------------------
# Extract spatial counts + coordinates from a Seurat spatial object
# ---------------------------------------------------------------------------
spatial_from_seurat <- function(st_obj, assay = NULL, image = NULL) {
  M <- get_counts(st_obj, assay = assay, layer = "counts")
  co <- Seurat::GetTissueCoordinates(st_obj, image = image)
  # Seurat returns different column names by version/platform
  cn <- colnames(co)
  xy <- if (all(c("x", "y") %in% cn)) co[, c("x", "y")]
        else if (all(c("imagerow", "imagecol") %in% cn)) co[, c("imagerow", "imagecol")]
        else co[, 1:2]
  xy <- as.matrix(xy)
  storage.mode(xy) <- "double"

  # align spots between counts and coordinates -- they can differ
  common <- intersect(colnames(M), rownames(xy))
  if (length(common) == 0) {
    if (nrow(xy) == ncol(M)) {
      warning("Spot names do not match between counts and coordinates; ",
              "assuming they are in the same order.")
      common <- colnames(M)
    } else stop("Cannot align spots between counts and coordinates.")
  } else {
    M <- M[, common, drop = FALSE]
    xy <- xy[common, , drop = FALSE]
  }
  # 0.2.2: counts stay sparse (spots x genes); deconvolve_spots() accepts them
  list(counts = Matrix::t(M), coords = xy, spots = common,
       gene_names = rownames(M))
}

# ---------------------------------------------------------------------------
# THE most important function in this file
# ---------------------------------------------------------------------------
#' Align the gene axis between a signature and a target matrix.
#'
#' SpaCoDA indexes genes positionally. If the scRNA reference and the spatial
#' matrix have different gene orders -- and they always do -- then
#' passing them in unaligned silently deconvolves noise and produces
#' confident, wrong answers with no error message. This is by far the most
#' common way to get garbage out of any deconvolution tool.
align_genes <- function(sig, target_counts, target_genes, verbose = TRUE) {
  common <- intersect(sig$gene_names, target_genes)
  if (length(common) < 100) {
    stop("Only ", length(common), " genes in common. Check that both use the ",
         "same identifier type (SYMBOL vs ENSEMBL) and the same species.")
  }
  if (verbose) {
    message("Gene alignment: ", length(common), " common / ",
            length(sig$gene_names), " reference / ", length(target_genes),
            " target")
  }
  si <- match(common, sig$gene_names)
  ti <- match(common, target_genes)
  sig2 <- sig
  sig2$B <- sig$B[si, , drop = FALSE]
  sig2$B <- sweep(sig2$B, 2, pmax(colSums(sig2$B), 1e-12), "/")
  sig2$gene_names <- common
  list(sig = sig2, counts = target_counts[, ti, drop = FALSE],
       gene_names = common)
}

# ---------------------------------------------------------------------------
# Task A, end to end from Seurat objects
# ---------------------------------------------------------------------------
#' @param lam spatial prior strength. Defaults to `0` (no spatial smoothing),
#'   the setting used in the SpaCoDA paper. `lam > 0` enables the graph-Laplacian
#'   smoothness term; we benchmarked it against ground truth on simulated data
#'   spanning 4-15 cells per spot and it reduced accuracy at every setting
#'   tested, so it is off by default. The function still prints your median UMI
#'   per spot as a diagnostic.
#' @param assay_name name of the assay that receives the proportions when
#'   `add_to_object = TRUE` (default `"SpaCoDA"`; `"spadecon"` before 0.2.2)
deconvolve_seurat <- function(sc_obj, st_obj, celltype_col,
                              sc_assay = NULL, st_assay = NULL,
                              n_per_type = 150, lam = NULL,
                              add_to_object = TRUE, verbose = TRUE,
                              assay_name = "SpaCoDA") {
  sig <- signature_from_seurat(sc_obj, celltype_col, assay = sc_assay,
                               verbose = verbose)
  sp <- spatial_from_seurat(st_obj, assay = st_assay)
  al <- align_genes(sig, sp$counts, sp$gene_names, verbose = verbose)

  umi <- median(rowSums(al$counts))
  if (is.null(lam)) {
    # Spatial smoothing is off by default: it is the paper's setting, and the
    # benchmark found lam > 0 harmful across the tested range of cellularity.
    # Set lam explicitly to explore it (see ?deconvolve_spots).
    lam <- 0
    if (verbose) message("median UMI/spot = ", round(umi),
                         " -> using lam = 0 (no spatial smoothing); ",
                         "set lam explicitly to enable it")
  }

  idx <- select_markers(al$sig, n_per_type = n_per_type,
                        restrict_to = al$gene_names)
  if (verbose) message("Using ", length(idx), " marker genes")

  res <- deconvolve_spots(al$counts, al$sig, gene_idx = idx,
                          coords = sp$coords,
                          lam = lam, weights = "nb",
                          platform_correction = TRUE,
                          engine = "rcpp", verbose = FALSE)
  colnames(res$W_cell) <- sp$spots
  rownames(res$W_cell) <- al$sig$types
  res$spots <- sp$spots
  res$coords <- sp$coords

  if (add_to_object) {
    # proportions as an assay -> plottable with SpatialFeaturePlot
    # 0.2.2: the assay is called "SpaCoDA" (was "spadecon"; assay_name = "spadecon"
    # keeps the old name)
    st_obj[[assay_name]] <- SeuratObject::CreateAssayObject(data = res$W_cell)
    res$object <- st_obj
    if (verbose) message("Added '", assay_name, "' assay. Plot with:\n",
                         "  SpatialFeaturePlot(obj, features = '",
                         al$sig$types[1], "')  # after DefaultAssay(obj) <- '", assay_name, "'")
  }
  res
}
