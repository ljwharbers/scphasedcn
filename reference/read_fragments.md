# Read fragments (whole file or tabix regions)

Read fragments (whole file or tabix regions)

## Usage

``` r
read_fragments(file, regions = NULL, cells = NULL)
```

## Arguments

- file:

  fragments.tsv.gz

- regions:

  optional data.table (chrom, start, end) from
  [`parse_region()`](https://ljwharbers.github.io/scphasedcn/reference/parse_region.md);
  uses the tabix index

- cells:

  optional character vector of cells to keep

## Value

data.table chrom, start, end, cell, strand, mapq, flags (+ extras),
`hap` (the file's `hap` column as U / A / B / N labels; "U" for every
fragment when the file has none, nothing is decoded from `flags`), `len`
(end - start when absent)
