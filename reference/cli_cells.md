# Cell subset from a file (one cell per line or a TSV with a `cell` column)

A one-column file has no header unless its first line is `cell` (as
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
reads its cell files); barcodes match with or without a trailing `-1`.

## Usage

``` r
cli_cells(path, all)
```

## Arguments

- path:

  file path or NULL

- all:

  all cell names (returned when `path` is NULL)

## Value

character vector
