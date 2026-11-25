# Normalise read depth against a reference

Expected counts per bin and cell, and the log2 ratio of observed to
expected (centred per cell). Every reference of a cell leaves that cell
out.

- `loo`: the pseudobulk of all other reference cells.

- `celltype`: the pseudobulk of the other reference cells of the same
  cell type; cells without a type, and types with fewer than
  `min_cells_type` reference cells, use `loo`. The normal-tissue
  reference: it removes the cell-type-specific chromatin signal in the
  depth.

- `model`: a robust Poisson GLM of the bin counts on GC content,
  mappability proxies and peak coverage, plus a per-cell GC trend
  (`per_cell_gc`, estimated as `gc_method` says). The tumour reference:
  it needs no copy-number-normal cells.

## Usage

``` r
normalise_depth(
  scn,
  method = "loo",
  use = NULL,
  cell_type = NULL,
  ref_cells = NULL,
  min_cells_type = 20,
  per_cell_gc = TRUE,
  gc_method = "lowess",
  ...
)
```

## Arguments

- scn:

  scn object

- method:

  "loo", "celltype" or "model"

- use:

  logical bin filter (default
  [`bin_filter()`](https://ljwharbers.github.io/scphasedcn/reference/bin_filter.md)
  on bins and counts)

- cell_type:

  cell type per cell (default `scn$cells$cell_type`)

- ref_cells:

  optional selector of cells allowed as references (e.g. clean singlets)

- min_cells_type:

  celltype: types with fewer reference cells fall back to the global
  reference

- per_cell_gc:

  model: subtract a per-cell robust GC trend after the shared GLM

- gc_method:

  model with `per_cell_gc`: estimator of the per-cell GC trend.
  `"lowess"` (default): a robust lowess of the bin log2 ratio
  `log2((x + 0.5) / (e + 0.5))` on bin GC. At low depth (most bins with
  0 or 1 fragments) the mean of that log ratio depends nonlinearly on
  `e`, so this fit finds a GC trend that depends on the cell's depth
  even in a cell without one. `"ratio"`: the ratio of summed observed to
  summed expected counts in 20 GC quantile groups of the bins, smoothed
  across the groups (weighted local-linear loess, span 0.75) and
  evaluated at each bin's GC; unbiased at every depth. Bins without GC
  keep the shared model; a cell without fragments gets no trend.
  ([`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  `gc_method`)

- ...:

  unused

## Value

`scn_norm`: list(method, ratio (bins x cells log2, NA on unused bins),
expected (bins x cells counts, NA on unused bins), use, size, params)
