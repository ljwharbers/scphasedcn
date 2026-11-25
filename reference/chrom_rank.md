# Order chromosome names naturally (chr1..chr22, chrX, chrY, chrM, others)

Other names (e.g. unplaced contigs) come last, in C-locale (byte) order
whatever the session's collation.

## Usage

``` r
chrom_rank(x)
```

## Arguments

- x:

  character vector of chromosome names

## Value

integer ranks usable with [`order()`](https://rdrr.io/r/base/order.html)
