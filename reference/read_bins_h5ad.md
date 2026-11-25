# Read a bins file (`bins_<tag>.h5ad`, by default `bins_100kb.h5ad`; [`bins_file()`](https://ljwharbers.github.io/scphasedcn/reference/bins_file.md)) into an scn object (missing layers tolerated)

Read a bins file (`bins_<tag>.h5ad`, by default `bins_100kb.h5ad`;
[`bins_file()`](https://ljwharbers.github.io/scphasedcn/reference/bins_file.md))
into an scn object (missing layers tolerated)

## Usage

``` r
read_bins_h5ad(file, layers = NULL, cells = NULL)
```

## Arguments

- file:

  path

- layers:

  layers to read (default: all present)

- cells:

  optional character vector of cells to keep

## Value

scn object
