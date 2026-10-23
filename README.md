# SpaCoDA <img src="man/figures/logo.svg" align="right" height="140"/>

Fast compositional deconvolution maps recurrent cellular niches across hundreds of spatial transcriptomics sections

---

## The method in brief

```
scRNA-seq reference ──► build_signature() ──► cell-type signature ──► select_markers()
                                                                          │
spatial counts ──────────────────────────────────────────────────► deconvolve_spots()
                                                                          │  cell-type proportions per spot
                                                                          ▼
coordinates ─────────────────────────────────► define_niches() / cluster_single_sample()   (per section)
                                                                          │  local niches + prototypes
                                                                          ▼
                                                              integrate_slices()           (across sections)
                                                                          │
                                                                          ▼
                                                         meta-niches of every spot
```

**1. Deconvolution.** Each spot is a non-negative mixture of the reference cell-type
expression profiles. The marker-gene counts of each spot are scaled to proportions and
fitted by non-negative least squares, weighted by the negative-binomial mean–variance
relation. The weights and a gene-level platform factor, which absorbs systematic
differences between the reference and the spatial assay, are re-estimated over three
rounds. The fitted mRNA shares are converted to cell proportions using the cell size of
each type. 

**2. Niches within a section.** The compositions of the spots are compared by the
Aitchison distance (log-ratio geometry for proportions). This distance is fused with the
normalised Euclidean distance between spot coordinates (weights 0.7 / 0.3). The fused
distance is embedded by principal coordinates and used to build a shared-nearest-neighbour
graph. Among the candidate resolutions, the Louvain partition with the highest modularity
is kept. Each local niche is summarised by its mean composition (prototype).

**3. Meta-niches across sections.** The prototypes of all sections are normalised,
reduced by PCA, batch-corrected with Harmony and clustered together. The meta-niche label
is then mapped back to every spot.

---

## Installation

```r
install.packages(c("Matrix", "nnls", "FNN", "Rcpp", "RcppArmadillo", "igraph", "irlba"))
install.packages(c("Seurat", "harmony", "RSpectra"))   # niche pipeline

# from the source tarball
devtools::install("EvoSpatialTool/SpaCoDA")
```
---

## Tutorial

### 0. Example data

A single-cell reference with 6 cell types and two spatial sections on a Visium-like grid.
Four niche types occupy regions of the tissue and recur in several places. Section B has no
myeloid niche.

```r
library(SpaCoDA)
set.seed(2026)

K <- 6; G <- 300
types <- c("Tumor", "CAF", "Tcell", "Bcell", "Myeloid", "Endo")

## expression profile of each cell type (40 specific genes per type)
B_true <- matrix(rgamma(G * K, 0.4, 1), G, K)
for (k in 1:K) B_true[(k - 1) * 40 + 1:40, k] <- B_true[(k - 1) * 40 + 1:40, k] + rgamma(40, 3, 0.3)
B_true <- sweep(B_true, 2, colSums(B_true), "/")
genes <- paste0("gene", 1:G)

## single-cell reference: cells x genes counts and cell-type labels
sc_labels <- rep(types, each = 120)
sc_counts <- t(sapply(match(sc_labels, types), function(k) rnbinom(G, mu = B_true[, k] * 2500, size = 5)))
dimnames(sc_counts) <- list(paste0("cell", seq_len(nrow(sc_counts))), genes)

## cell-type composition of four niche types
niche_comp <- rbind(
  tumor_core = c(.75, .10, .03, .02, .07, .03),
  stroma     = c(.10, .55, .05, .02, .13, .15),
  immune_TLS = c(.05, .08, .40, .35, .10, .02),
  myeloid    = c(.25, .10, .05, .02, .50, .08))
colnames(niche_comp) <- types

## one section: spots x genes counts, coordinates and the true niche of each spot
simulate_sample <- function(name, region_types, seed) {
  set.seed(seed)
  g  <- expand.grid(i = 0:31, j = 0:31)
  co <- cbind(x = g$i * 100 + (g$j %% 2) * 50, y = g$j * 86.6)
  n  <- nrow(co)
  centres <- co[sample(n, length(region_types)), ]
  region  <- max.col(-sapply(seq_along(region_types),
                             function(r) colSums((t(co) - centres[r, ])^2)), "first")
  truth   <- rownames(niche_comp)[region_types[region]]
  counts  <- matrix(0, n, G)
  for (i in 1:n) {
    cells <- tabulate(sample(K, rpois(1, 6) + 1, TRUE, niche_comp[truth[i], ]), K)
    counts[i, ] <- rpois(G, B_true %*% cells * 300)
  }
  dimnames(counts) <- list(paste0(name, "_spot", 1:n), genes)
  rownames(co) <- rownames(counts)
  list(counts = counts, coords = co, truth = setNames(truth, rownames(counts)))
}
sA <- simulate_sample("A", c(1, 2, 3, 4, 1, 2, 3, 1), seed = 1)
sB <- simulate_sample("B", c(1, 2, 3, 1, 2, 3, 2, 1), seed = 2)
```

### 1. Cell-type deconvolution

```r
## reference signature (cells x genes in) and marker genes
sig <- build_signature(sc_counts, sc_labels, gene_names = colnames(sc_counts))
mk  <- select_markers(sig, n_per_type = 40)

## deconvolve each section (spots x genes in, genes in the order of the signature)
resA <- deconvolve_spots(sA$counts, sig, gene_idx = mk, verbose = FALSE)
resB <- deconvolve_spots(sB$counts, sig, gene_idx = mk, verbose = FALSE)

round(head(t(resA$W_cell)), 2)   # spots x cell types, rows sum to 1
```

`res$W_cell` holds cell proportions (cell types × spots) and `res$W_mrna` the mRNA shares
before conversion to cells.

> **Real data: align genes first.** `deconvolve_spots()` matches genes by position. Use
> `al <- align_genes(sig, st_counts, colnames(st_counts))` and pass `al$counts` and `al$sig`.

### 2. Niches in one section

```r
nA <- define_niches(resA$W_cell, sA$coords, verbose = FALSE)   # cell types x spots in
table(niche = nA$labels, truth = sA$truth)    # local niches against the simulated truth
round(nA$profiles, 2)                         # prototype (mean composition) of each niche
nA$best_resolution                            # resolution with the highest modularity
nA$spatial_test$p_value                       # is the composition spatially organised?
```

Because the spatial distance is part of the clustering distance, a niche type that occurs
in several places of a section is split into several local niches; they are brought
together again at the meta-niche step.

### 3. Meta-niches across sections

```r
props  <- list(A = t(resA$W_cell), B = t(resB$W_cell))   # spots x cell types per section
coords <- list(A = sA$coords, B = sB$coords)
meta <- integrate_slices(props, coords, names(props), resolution = NULL, verbose = FALSE)

truth <- c(sA$truth, sB$truth)
table(meta = meta$spot_super_niche, truth = truth[names(meta$spot_super_niche)])
round(meta$niche_summary_prop, 3)    # share of each meta-niche per section
meta$spatial_tests                   # spatial test per section
```

`resolution = NULL` picks the meta-niche resolution with the highest modularity. 

### 4. Seurat objects

```r
## sc_obj: Seurat reference with a "celltype" column; st_obj: spatial Seurat object
# res <- deconvolve_seurat(sc_obj, st_obj, celltype_col = "celltype")
# nch <- define_niches(res$W_cell, res$coords)
# res$object$niche <- nch$labels[colnames(res$object)]
# Seurat::SpatialDimPlot(res$object, group.by = "niche")
```

---

## Inputs

| Input | Shape | Notes |
|---|---|---|
| reference counts | cells × genes | raw counts, at least 25 cells per type |
| reference labels | one per cell | |
| spatial counts | spots × genes | raw counts, same genes in the same order as the signature (`align_genes()`) |
| coordinates | spots × 2 | same spot order as the counts |

| Function | Takes | Returns |
|---|---|---|
| `deconvolve_spots` | spots × genes | `W_cell`: cell types × spots |
| `define_niches` | cell types × spots (`res$W_cell`) | `labels`, `profiles`, `spatial_test` |
| `integrate_slices` | list of spots × cell types | `spot_super_niche`, `niche_summary_prop`, `spatial_tests` |

## Main settings

| Setting | Default | Meaning |
|---|---|---|
| `select_markers(n_per_type =)` | 150 | marker genes per cell type |
| `deconvolve_spots(lam =)` | 0 | spatial smoothing; lowered accuracy in every simulation we tried |
| `define_niches(expr_weight =, spat_weight =)` | 0.7, 0.3 | weights of the composition and spatial distance |
| `define_niches(resolution_seq =)` | 0.3–3 | candidate resolutions; the one with maximum modularity is kept |
| `integrate_slices(resolution =)` | numeric | meta-niche resolution; `NULL` selects it by modularity |

