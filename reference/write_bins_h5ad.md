# Write an scn object as a bins file (`bins_<tag>.h5ad`, anndata layout, gzip only)

Write an scn object as a bins file (`bins_<tag>.h5ad`, anndata layout,
gzip only)

## Usage

``` r
write_bins_h5ad(scn, file, layers = names(scn$layers))
```

## Arguments

- scn:

  scn object

- file:

  output path

- layers:

  which layers to write (default all)

## Value

file path (invisibly)
