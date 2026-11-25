# Data for repeated clone contrasts

Pre-computes what the clone contrast (`contrast_stat`) needs for many
splits: the per-cell log2 ratio matrix (cells x regions) and the
molecule votes per cell x region x phase block as sparse cells x blocks
matrices.

## Usage

``` r
contrast_data(
  feat_long,
  votes,
  atoms,
  regions,
  cells,
  classes = c("arm", "20Mb"),
  allelic_chroms = NULL
)
```

## Arguments

- feat_long:

  features long table

- votes, atoms, regions:

  votes, atoms and regions tables

- cells:

  cell names in vote index order (rows of the result)

- classes:

  region classes used

- allelic_chroms:

  chromosomes with an allele channel (default all)

## Value

list(L, regions, MA, MB, blk_region, N, cells)
