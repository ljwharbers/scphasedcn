# Subset an scn object by cells and/or bins

Subset an scn object by cells and/or bins

## Usage

``` r
subset_scn(scn, cells = NULL, bins = NULL)
```

## Arguments

- scn:

  scn object

- cells:

  logical, integer or character (cell names) selector of cells (NULL =
  all). Integer and character selectors give the cells in the requested
  order (negative integers drop cells); a selector that names a cell
  twice, or a cell (name or index) the object does not have, stops

- bins:

  logical/integer selector of bins (NULL = all)

## Value

scn
