# SpaCoDA 0.3.0 (paper methods)

This version keeps only the methods used in the accompanying paper:

* deconvolution: `deconvolve_spots()` -- variance-weighted NNLS with the
  gene-level platform correction (the former `model = "wls"`, which was the
  default);
* niches: `define_niches()` / `cluster_single_sample()` per section and
  `integrate_slices()` across sections, plus `niche_spatial_test()`.

Removed from 0.2.2: the count model (`deconvolve_spots(model = "nb")`,
`informative_genes()`) and niche method v2 (`find_niches()`, `plot_niches()`).
They are in SpaCoDA 1.0.0. Results of the remaining functions are unchanged
from 0.2.2.

# SpaCoDA 0.2.2

**The package is renamed from `CoDA` to `SpaCoDA`.** `CoDA` differs from CRAN's
`coda` only by case, which CRAN and Bioconductor do not accept. Function names
and arguments are unchanged: replace `library(CoDA)` by `library(SpaCoDA)`. The
`find_niches()` result class is `spacoda_niches` and the default assay of
`deconvolve_seurat()` is `"SpaCoDA"` (both new in this version).

Niche method v2, bug fixes, packaging, a faster niche PCoA and a
spatial-organisation test. **Paper-pipeline niche labels are unchanged** (see below). **Deconvolution numbers are unchanged**: built
against the same BLAS, `deconvolve_spots()` returns bit-identical proportions to
0.2.1 (checked on a 600-spot simulation, `lam = 0` and `lam = 0.5`). Because the
package now links R's own BLAS instead of OpenBLAS (see Packaging), a given
installation may differ from a 0.2.1 build in the last bit (measured max |diff|
3e-16). Sparse input gives the dense answer to the same rounding.

## Removed

* **Niche-level deconvolution of bulk samples is no longer part of the package.**
  `deconvolve_bulk()`, `build_niche_signature()`, `select_niche_markers()`,
  `niche_identifiability_report()`, `niche_to_celltype()` and
  `niche_bulk_from_seurat()` are removed. The C++ solver they shared with the
  count model, `bulk_solve_cpp()`, is renamed `wnnls_cpp()`. Entries below that
  mention these functions describe earlier development versions.

## New

* **`deconvolve_spots(model = "nb")`: a count model, optional.** The default
  `model = "wls"` (the paper's method) is unchanged and returns bit-identical
  results. `model = "nb"` fits the counts as negative binomial with the spot's
  depth as offset and a per-gene dispersion estimated from the residuals; each
  Fisher-scoring step is still one K x K weighted NNLS per spot, and the steps are
  damped by one half (undamped, the iteration cycled with period ~4). The
  platform factor is estimated once from the pseudo-bulk of all expressed genes:
  the bulk composition by least squares in log space, then a per-gene factor
  shrunk towards 1 (empirical Bayes); it is returned as `$gamma`. Spots are
  processed in chunks, so all genes can be used on large data. `lam > 0` is not
  available with `model = "nb"`. `deconvolve_seurat()` passes `model` through.
* **`gene_idx = "informative"` / `informative_genes()`**: every gene expressed in
  the reference and detected in the spatial data; `gene_idx` also accepts gene
  names now.
* Benchmark on simulations (10 cell types incl. two near-identical pairs,
  cell-to-cell heterogeneity, 1 + Poisson(5) cells per spot, gene-level platform
  effect sd 0.7 on the log scale, 600 spots, 2 seeds; RMSE of cell proportions,
  same platform / cross-platform / cross at ~700 UMI per spot / cross at ~15,000):
  default 0.0200 / 0.0259 / 0.0401 / 0.0214; `wls` with all informative genes
  0.0172 / 0.0222 / 0.0354 / 0.0181; `nb` with all informative genes 0.0163 /
  0.0224 / 0.0352 / 0.0181 (about 4 s); RCTD 2.2.1 full mode 0.0164 / 0.0229 /
  0.0356 / 0.0186 (about 400 s on 2 cores). Most of the gain comes from the larger
  gene set.
* **Niche method v2: `find_niches()`, `plot_niches()`.** One call for one or many
  samples. Each spot is described by the CLR composition of its spatial
  neighbourhood at two scales (6 and 18 nearest spots; 20 and 50 for one-hot
  single-cell input); all spots of all samples are clustered together (Leiden on
  a kNN graph), over-split and merged until every pair of niches differs by at
  least `min_diff` (default 0.1) in the proportion of some cell type; clusters that
  are a mixture of two spatially adjacent niches are recorded as their interface
  (`$interface`) instead of a niche. Before that, a per-sample multiplicative
  bias per cell type is estimated (mutual nearest neighbours between samples,
  log ratios where the type is detected in both spots, median polish) and
  removed: K - 1 numbers per sample, which cannot remove a sample-specific niche
  or equalise niche abundances. Returns labels, profiles, per-sample abundance,
  the removed bias, and `niche_spatial_test()` per sample. No Seurat / harmony
  needed. Measured on simulations (Visium-like lattice, Poisson cells per spot,
  boundary mixing, deconvolution noise; v1 = the paper pipeline):
  5 niche types recurring in 12 regions -> v2 5-6 niches, ARI 0.86-0.92, 0.5 s
  (v1: 26-30 niches, ARI 0.15-0.23, 11-18 s); structure-free composition -> v2 1
  niche (v1 19-24); four samples with per-sample cell-type bias (log-scale sd up
  to 0.7) and a niche present in one sample only -> v2 6-7 niches (truth 6), ARI
  0.84-0.89, that niche kept (recall >= 0.90; without the batch step ARI drops to
  0.39-0.55 at sd 0.7; v1 `integrate_slices(resolution = NULL)` ARI 0.14-0.29, the
  sample-specific niche merged with others, precision <= 0.20); one-hot
  single-cell data -> 4-6 niches (truth 4), ARI 0.80-0.88; 53,000 spots in 19 s.
  `min_diff` 0.07-0.15 gives mean ARI 0.87-0.89.
  The paper pipeline (`define_niches()`, `integrate_slices()`) is unchanged.
* `tutorials/00_tutorial_zh.md`: a Chinese tutorial that runs end to end on
  simulated data (deconvolution, v2 niches in one and two samples, v1, bulk).
  `tutorials/01_quickstart.md` and `02_niche_discovery.md` were rewritten: 01
  simulated non-count data, aligned genes without renormalising the signature,
  had a plotting example that failed (`aes(x, y)` on columns `X1`, `X2`) and
  described `deconvolve_seurat()` as returning the Seurat object; 02 still
  described the pre-0.2.0 k-means / CLR methods as the defaults. Both now run as
  written. README, `?deconvolve_spots`, `?deconvolve_seurat` and the package
  help were brought up to date (default `lam = 0`, no shipped files that do not
  exist, v2 niches).
* **Faster niche PCoA, same niches.** `cluster_single_sample()` /
  `define_niches()` gain `pcoa_method = c("auto", "partial", "full")`. `"full"` is
  the study's `cmdscale()`, a full O(N^3) `eigen()` of which only the top
  `pcoa_dim` (30) eigenpairs are used. `"partial"` computes just those eigenpairs
  of the same matrix with RSpectra's Lanczos solver, applied as an operator so the
  double-centred matrix is never formed, and falls back to `cmdscale()` if it
  does not converge. `"auto"` (default) uses `"partial"` from 1,000 spots when
  RSpectra (Suggests; a Seurat dependency) is installed. Checked against 0.2.1's
  `cmdscale()` path on 8 sections of 1,000-4,000 spots (random and Visium-like
  hexagonal layouts, structured, weak and pure-noise compositions) and on a
  three-section `integrate_slices()` run at `resolution = 6` and `NULL`: niche
  labels, prototypes, resolution scans and super-niches are identical. Runtime of
  `cluster_single_sample()`: 1.4x faster at 1,000 spots, 3.3x at 3,000, 4.6x at
  4,000 (130 s -> 28 s). The N x N distance matrices remain, so memory is
  still O(N^2).
* **`niche_spatial_test()`: is the composition spatially organised at all?** A
  permutation test whose statistic is the mean squared CLR (Aitchison) distance
  between each spot and its k nearest spatial neighbours (multivariate Geary's C),
  against compositions shuffled across fixed coordinates. It does not look at the
  clustering, which matters: the fused distance tiles even structure-free
  compositions into contiguous "niches" (7 on random compositions, 72% of
  neighbouring spots in the same niche vs 15% by chance), and a test on the
  clustering's modularity with a shuffle null fails in the other direction
  (clearly structured sections come out non-significant). Measured: pure noise
  p = 0.68; weak, clear and smooth-gradient structure p = 0.001; on 200 noise
  sections the rejection rate is 5.5% at 0.05. ~6 s for 20,000 spots. A
  significant result means spatial organisation, which includes gradients -- not
  by itself that the structure is made of discrete niches.
  `cluster_single_sample()` and `define_niches()` run it by default
  (`spatial_test = TRUE`, `n_perm = 999`) and return it as `spatial_test`;
  `integrate_slices()` collects one row per section in `spatial_tests`. It runs
  after the clustering on its own random-number stream, so niches and super-niches
  are exactly what they are without it.
* README: the "fraction of neighbouring spots in the same niche" diagnostic is
  replaced by `niche_spatial_test()`, since the spatial distance term makes that
  fraction high on noise too.

## Bug fixes

* **`niche_to_celltype()` mixed niches up from 11 niches on.** The prototypes of
  `define_niches()` / `cluster_single_sample()` come out of `aggregate()` in
  character order (`0 1 10 11 2 ...`), while `build_niche_signature()` and
  `deconvolve_bulk()` use the factor order (`0 1 2 ... 10 11`). The product was
  positional, so niche 10 was paired with niche 2's profile, and so on, with no
  error. Niches are now matched by row name; differing niche sets are an error.
* **`niche_bulk_from_seurat(n_niches = )` had no effect.** It passed `n_niches` to
  `define_niches()`, whose default method since 0.2.0 ignores it. The default is
  now `n_niches = NULL` (the paper's method, which picks the number itself); a
  number selects the k-means path with exactly that many niches.
* **`integrate_slices()$spot_super_niche` was unnamed**, although documented as
  named by spot: the names were read after `as.character()` had dropped them, so
  matching super-niches to spots by name returned `NA` throughout. It is now named.
* **`signature_from_seurat()` did not notice a misspelt `celltype_col`**
  (`as.character(NULL)` is `character(0)`, not `NULL`); it now errors.

## Behaviour changes

* `deconvolve_seurat()` stores the proportions in an assay called `"SpaCoDA"`
  (was `"spadecon"`, the package's former working name); the new argument
  `assay_name = "spadecon"` keeps the old name for existing scripts.

* `define_niches(method = "kmeans")`: the `profiles` row names are the label values
  `"1".."J"` (were `"niche1".."nicheJ"`), as in the default path, so that they can
  be matched to `build_niche_signature()` / `deconvolve_bulk()` output by name.
* `signature_from_seurat()` and `spatial_from_seurat()` keep the counts sparse.
  They used to call `as.matrix()`, a dense copy of the whole reference (about
  12 GB for 50k cells x 30k genes). `spatial_from_seurat()$counts` is therefore a
  sparse `Matrix` now.

## Packaging

* `R CMD check` passes: `tests/verify_port_vs_original_20261006_01.r` (which needs
  the study's private data) is excluded from the build and not part of the public repository.
* `src/Makevars` links R's own `$(LAPACK_LIBS) $(BLAS_LIBS)` instead of a
  hard-coded `-lopenblas` that fails where OpenBLAS is not installed; the unused
  OpenMP flags are gone. Compiled objects (`*.o`, `*.so`) are no longer shipped.
* Tests that need Seurat are skipped when it is not installed; `test-core.R` no
  longer uses `%||%` (base R >= 4.4 only, while the package declares R >= 4.1).
* `deconvolve_spots()` no longer triggers Matrix's deprecation warning for
  `as(<ddiMatrix>, "dgCMatrix")`, and the niche functions no longer trigger
  Seurat's "Data is of class matrix. Coercing to dgCMatrix." (they pass the sparse
  matrix Seurat converts to anyway; niche labels and prototypes are identical).

# SpaCoDA 0.2.1 (released as CoDA)

**Default change: `deconvolve_spots(lam = )` is now `0`, not `0.5`. No reported
number changes** — every benchmark and atlas run in the paper passed `lam = 0`
explicitly, so this only aligns the shipped default with what was actually
validated. Previously a user calling the package with defaults got spatial
smoothing that the paper never used.

## Default changes

* **`deconvolve_spots(lam = 0)`.** The graph-Laplacian smoothness term is now
  off by default. Measured on the paper's three simulation difficulty tiers
  (15, 6 and 4 cells per spot), adding `lam = 0.5` cost **0.13, 0.09 and 0.13
  Pearson** respectively — larger than the entire gap between SpaCoDA and RCTD
  (0.019-0.046). The harm was not monotonic in cellularity, so no
  depth- or density-based rule justifies enabling it automatically.
* **`deconvolve_seurat(lam = NULL)` now resolves to `0`.** It previously
  auto-selected from the median UMI per spot (`> 3000` -> 0, `> 800` -> 0.25,
  else 0.5), which enabled spatial smoothing on essentially all real Visium
  data — including the paper's own Xenium-to-Visium datasets (median UMI/spot
  1274, which the old rule mapped to 0.25). The function still prints the
  median UMI as a diagnostic, but no longer acts on it.

## Documentation

* `?deconvolve_spots` and the README's "Choosing `lam`" section now state that
  the Laplacian term is off by default and off in the paper, with the numbers
  behind that choice. `tutorials/01_quickstart.md` no longer describes the
  graph term as active.

## Fixes in the prototype (super-niche) stage

Found by running the ported pipeline on inputs smaller than the study's — the
tutorial's two halves of one Xenium section, 19 prototypes. None of these change
study-scale results: the first three are exact no-ops once there are 60 or more
prototypes, and the fourth only replaces a crash or a silent `NA`.

* **The modularity scan looked for a graph called `Spatial_snn` by name.** For any
  object whose assay is not literally `Spatial`, `obj@graphs[["Spatial_snn"]]` was
  `NULL`, igraph then errored on a zero-vertex graph
  (`Membership vector size differs from number of vertices`), `compute_modularity()`
  swallowed it into `NA`, and **every** candidate resolution looked invalid. The SNN
  graph is now resolved from the object's own assay, with a fallback to whichever
  `*_snn` graph the object carries.
* **`k.param` is capped at a third of the prototypes.** The study value 20 assumes
  hundreds to thousands of prototypes (BRCA ~2,500, GBM ~1,800). On 19 prototypes it
  is nearly the whole point cloud, the SNN graph loses its structure, and the scan
  collapses: measured on 19 points with three planted groups, `k.param = 20` left
  **1** usable resolution in the whole `seq(0.3, 3, by = 0.1)` sweep (best
  modularity 0.0198, k = 2), while the capped `k.param = 6` left **all 28** and
  recovered the planted structure (resolution 0.30, k = 3, modularity 0.6041). The
  cap is `N %/% 3`, so it equals 20 — the study value — from 60 prototypes upwards.
* **`best_modularity` is `NA`, not `NULL`, when the scan falls back.** If no
  resolution was valid the code fell back to `resolution = 1` and then read
  `opt$best_modularity` from the `NULL` scan result, so the returned list carried
  `best_modularity = NULL` and any downstream `round(sup$best_modularity, 3)` died
  with `non-numeric argument to mathematical function`.
* **A resolution too high for the prototype set now warns.** Louvain then puts every
  prototype in its own community and Seurat labels them all `"singleton"`, which
  used to be returned as if it were a partition (the super-niche figure showed a
  single colour, modularity 0). `integrate_slices()` now warns and names the size of
  the prototype set, so the caller can switch to `resolution = NULL`.

# SpaCoDA 0.2.0 (released as CoDA)

**The package's niche finding now is the method of the paper.** Until 0.1.1 the
package shipped a different niche implementation (spatial smoothing of the
proportions followed by k-means, and a CLR/PCA/Harmony integration), while the
BRCA and GBM analyses of the manuscript ran a separate study pipeline. The two
produced different niches from the same composition, so the package could not
reproduce the paper's niche results. 0.2.0 ports the study pipeline into the
package, so one implementation serves both.

Ported from `0.function_auto_select_resoultion_seurat_20260817.r`, as called by
`nicheDecon_brca_20260816_01_fixed.r` and `nicheDecon_gbm_20260816_01_fixed.r`
(those two call sites are byte-identical):

* `cluster_single_sample()` -- per-section niches: Aitchison distance on the
  composition, fused with the spatial distance (weights 0.7/0.3), PCoA, a
  shared-nearest-neighbour graph, and the resolution with maximum modularity.
  Homogeneous or too-small sections are returned as a single niche instead of
  being dropped.
* `compute_aitchison_distance()`, `fuse_distances()`,
  `select_resolution_by_modularity()` -- the building blocks, exported.
* `define_niches(..., method = "aitchison")` is now the default and forwards to
  `cluster_single_sample()`. `method = "kmeans"` keeps the pre-0.2.0 behaviour.
* `integrate_slices(..., method = "prototype")` is now the default: local niches
  and prototypes per section, then one clustering over the stacked prototypes
  (LogNormalize, scale, PCA, Harmony on the section label, SNN). The legacy
  CLR implementation is `method = "clr"`.

**Super-niche resolution is an argument.** The study scripts hard-coded
`FindClusters(resolution = 6)` on the prototypes; `integrate_slices()` now takes
`resolution = 6` as a parameter, and `resolution = NULL` selects it by maximum
modularity from `resolution_seq` instead.

Defaults follow the paper call site, which differs from the original helper's own
signature in two places: `use_rank_norm = FALSE` (helper default `TRUE`) and
`resolution_seq = seq(0.3, 3, by = 0.1)` (helper default `seq(0.3, 1.5, by = 0.1)`).

Verification: `tests/testthat/test-niche-port-0.2.0.R` and
`tests/verify_port_vs_original_20261006_01.r` run the ported functions and the
original study functions on the same input and compare labels and prototypes.
Deconvolution is untouched: `deconvolve_spots()` still returns exactly the numbers
of 0.1.1 (and of 0.1.0 on the default path).

# SpaCoDA 0.1.1 (released as CoDA)

Patch release. **No change to the fitted numbers**: with the default call
(`cell_counts = NULL`, `do_harmony = TRUE`) 0.1.1 returns exactly the same
proportions as 0.1.0 — verified bit-identical on the Xenium breast section used
in `tutorials/`. Everything below fixes either a feature that did not work at
all, or information that was computed and then thrown away.

## Bug fixes

* **`deconvolve_spots(cell_counts = )` did nothing.** The Dirichlet refinement
  matched `cell_counts` against `colnames(W_cell)`, but the solver returned
  `W_cell` without dimnames, so `intersect(colnames(W_cell), names(cell_counts))`
  was always empty and the function returned the *unrefined* composition while
  still reporting `W_cell_raw = NULL` (the only visible hint). The returned
  `W_mrna` / `W_cell` now carry `dimnames = list(cell types, spots)`, taken from
  `rownames(ST)`. `cell_counts` may be a named vector (spot id -> number of cells)
  or a plain vector in spot order; it now raises an informative error instead of
  silently returning the unrefined result when nothing matches.

* **`integrate_slices()` aborted with `data.use %*% cell.embeddings :
  non-conformable arguments` under harmony >= 1.2.** Harmony's `RunHarmony()`
  defaults to `project.dim = TRUE`, which calls `Seurat::ProjectDim()`;
  `ProjectDim` needs a `scale.data` slot, but the assay built here deliberately
  contains CLR values (no scaling step), so the call failed *after* Harmony had
  already converged, making it look like a dimensionality problem. `project.dim`
  is now explicitly `FALSE`.

## Improvements

* `define_niches()` now labels its output: `profiles` gets `colnames` = cell
  types and `rownames` = `niche1..nicheJ`, and `labels` is named by spot. The
  niche -> cell-type conversion (`niche_to_celltype()`) therefore returns
  readable matrices instead of unlabelled ones.
* `build_signature()` returns a labelled `B` (rows = genes, columns = cell types)
  and `build_niche_signature()` one labelled by gene and niche, so `B[, "T cell"]`
  works. This also stops `f_capture` / `eff_cell_sizes` from being anonymous vectors
  inside `deconvolve_spots()`.
* `deconvolve_bulk()` returns a labelled `W` (rows = niches, columns = samples)
  instead of a bare matrix.
* `integrate_slices()` returns the null-test statistics (`z`, `best_modularity`,
  `null_modularity`) that `select_resolution()` computes; previously they were
  only printed as messages.
* Regression tests for both bugs and for the new dimnames: `tests/testthat/test-regressions-0.1.1.R`.
* `tests/testthat/test-core.R` used `t(A)` on a sparse matrix, which only
  dispatches when `Matrix` is attached; it is now `Matrix::t(A)` so the suite
  passes under `R CMD check`.
* `tutorials/02_niche_discovery.md` corrected: `deconvolve_bulk()` takes samples x
  genes (the old snippet passed genes x samples), and it now points at the full worked
  tutorial on real data.
* `README.md` corrected: several tutorials previously showed function signatures
  that do not exist (e.g. `select_resolution(W_cell, coords, n_niches_seq)`,
  `niche_qc(W_cell, coords, labels)`, `integrate_slices(profiles = , sample = )`,
  `select_niche_markers(profiles, sig, n_top)`,
  `deconvolve_seurat(spatial_seurat, sig)`, `niche_to_celltype(profiles)`).

## Compatibility

Nothing was removed. `cell_counts` callers who passed an unnamed vector of the
wrong length used to get a silent no-op and now get an error; passing the right
thing (as documented) is unaffected.
