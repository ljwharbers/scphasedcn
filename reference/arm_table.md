# Chromosome arms (Battenberg coordinates)

p = `[0, cen_left)`, q = `[cen_right, chromosome end)` from the
Battenberg table (`chromosome_coordinates_hg38_chr.txt`: chr, start,
cen.left.base, cen.right.base, end). Battenberg names chrX `chr23` (and
chrY `chr24`); they are renamed. When the table has no chrX row, the
chrX centromere is `params$chrx_centromere`. Acrocentric p arms
(`params$acrocentric_p`) are dropped; chrY arms are kept only for males.

## Usage

``` r
arm_table(centromeres, chrom_sizes, sex = "female", params = caller_defaults())
```

## Arguments

- centromeres:

  path or table of the Battenberg coordinates

- chrom_sizes:

  named chromosome sizes (arms are made for these chromosomes only)

- sex:

  "female" or "male" (case-insensitive, or "F" / "M"; any other value
  stops)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

## Value

data.table `region, chrom, arm, start, end`
