# Event regions of a sample (config `events` map, region strings with 0-based half-open coordinates) as a regions table

Event regions of a sample (config `events` map, region strings with
0-based half-open coordinates) as a regions table

## Usage

``` r
cli_events(s, chrom_sizes = NULL)
```

## Arguments

- s:

  sample entry from
  [`cli_sample()`](https://ljwharbers.github.io/scphasedcn/reference/cli_sample.md)

- chrom_sizes:

  named sizes (closes open ends)

## Value

data.table region, chrom, start, end, event
