# MEDICC2 input table (MEDICC2 export, step 3)

One row per sample x segment with integer copy numbers capped at `cap`:
`sample_id, chrom, start, end, cn_a, cn_b` (allele-specific) or
`sample_id, chrom, start, end, cn` (total = cn_a + cn_b, for
`--total-copy-numbers --input-allele-columns cn`). Adjacent pieces of a
chromosome whose values are identical in every sample are merged (after
the cap; on cn_a / cn_b, or on cn for the total table), so every sample
has the same segments, as MEDICC2 requires. The normal clone `normal` is
a sample (MEDICC2 `--normal-name N0`); when absent it is added as 1:1 on
diploid and 1:0 on haploid pieces (male chrX), so it is the root of the
tree. Chromosomes are sorted naturally.

## Usage

``` r
medicc2_table(
  x,
  allele_specific = TRUE,
  normal = "N0",
  cap = 8L,
  chroms = paste0("chr", c(1:22, "X"))
)
```

## Arguments

- x:

  [`orient_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/orient_clone_states.md)
  result, or a long table with `clone, chrom, start, end, cn_a, cn_b` (+
  `ploidy`)

- allele_specific:

  TRUE: cn_a / cn_b columns; FALSE: total cn

- normal:

  name of the normal sample

- cap:

  maximal copy number (MEDICC2 caps at 8)

- chroms:

  chromosomes kept

## Value

data.table (MEDICC2 TSV layout); attribute `n_pieces` (pieces before
merging)
