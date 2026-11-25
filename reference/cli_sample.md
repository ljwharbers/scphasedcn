# Sample entry of a sample yaml (with its extraction directory)

The extraction directory is the sample's `extract_dir`, else
`<outputs.root>/<sample>/extract` (where `scphasedcn run` writes it);
NULL when the yaml gives neither.

## Usage

``` r
cli_sample(cfg, sample)
```

## Arguments

- cfg:

  sample yaml read with
  [`read_config()`](https://ljwharbers.github.io/scphasedcn/reference/read_config.md)

- sample:

  sample name

## Value

list with the sample fields plus `name`, `dir` (extraction output dir,
or NULL)
