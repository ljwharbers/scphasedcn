# 1 Mb window depth of a precomputed features object

Rebuilds the `win` component of an `scn_features` object saved without
it (e.g. to save disk space): the same cells and bin filter
(`features$use`),
[`normalise_depth()`](https://ljwharbers.github.io/scphasedcn/reference/normalise_depth.md)
with `method` (`model` for tumours, `celltype` with the features' cell
types and the per-cell GC estimator of the features' parameters,
`features$params$gc_method`), then
[`window_matrices()`](https://ljwharbers.github.io/scphasedcn/reference/window_matrices.md).
Needed by the tumour mode (clone pseudobulks).

## Usage

``` r
feature_windows(
  dir,
  features,
  method = "model",
  size = 1e+06,
  min_cells_type = 20
)
```

## Arguments

- dir:

  extraction directory

- features:

  `scn_features` (cells, use)

- method:

  depth normalisation ("model" tumour, "celltype" normal)

- size:

  window size (bp)

- min_cells_type:

  celltype: minimal reference cells of a type

## Value

list `x`, `e`, `windows`
([`window_matrices()`](https://ljwharbers.github.io/scphasedcn/reference/window_matrices.md))
