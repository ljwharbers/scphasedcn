# The bins file of an extraction directory

The extractor names the bins file after its bin size (`bins_<tag>.h5ad`:
`bins_100kb.h5ad` at the default 100 kb, e.g. `bins_50kb.h5ad` or
`bins_1Mb.h5ad` otherwise, or its `--tag`). Returns `bins_100kb.h5ad`
when the directory has it, otherwise the one `bins_*.h5ad` file in it;
stops when there is none or more than one (naming them).

## Usage

``` r
bins_file(dir)
```

## Arguments

- dir:

  extraction directory

## Value

path of the bins file
