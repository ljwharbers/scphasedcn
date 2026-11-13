# Subset a features object to a set of cells

Keeps the rows of `long` and `cells` of the given cells and re-indexes
the votes (whose `cell` is the row of `cells`) to the kept rows, so that
[`update_allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/update_allelic_features.md)
and the tumour mode still give each cell its own votes; the 1 Mb window
depth `win` keeps the kept cells' columns. Every other component is kept
as it is.

## Usage

``` r
subset_features(features, cells)
```

## Arguments

- features:

  `scn_features`
  ([`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md))

- cells:

  cells to keep (others are ignored)

## Value

the subset `scn_features` (the object itself when every cell is kept)
