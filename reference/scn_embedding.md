# Chromatin embedding of the cells of a sample

[`chromatin_embedding()`](https://ljwharbers.github.io/scphasedcn/reference/chromatin_embedding.md)
(high-pass filtered log ratios to the all-cell leave-one-out reference,
top principal components) of the given cells on the used bins of an
extraction directory, with the cell ids as row names, ready for
[`knn_pools()`](https://ljwharbers.github.io/scphasedcn/reference/knn_pools.md).
Built from depth only (no copy-number features), as kNN pools require.

## Usage

``` r
scn_embedding(
  dir,
  cells = NULL,
  use = NULL,
  n_pc = 20,
  hp_window = 21,
  seed = 1
)
```

## Arguments

- dir:

  extraction directory (bins_100kb.h5ad)

- cells:

  cells to embed (default all)

- use:

  logical bin filter over all bins of the sample (default
  [`bin_filter()`](https://ljwharbers.github.io/scphasedcn/reference/bin_filter.md)
  on chr1-22 + chrX)

- n_pc:

  components

- hp_window:

  high-pass window (bins)

- seed:

  RNG seed of the truncated SVD (irlba's random start)

## Value

cells x components matrix
