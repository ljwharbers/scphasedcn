# MEDICC2 input table (MEDICC2 export, step 3)

One row per sample x segment with integer copy numbers capped at `cap`:
`sample_id, chrom, start, end, cn_a, cn_b` (allele-specific) or
`sample_id, chrom, start, end, cn` (total = cn_a + cn_b, for
`--total-copy-numbers --input-allele-columns cn`). Adjacent pieces of a
chromosome whose values are identical in every sample are merged (after
the cap; on cn_a / cn_b, or on cn for the total table), so every sample
has the same segments, as MEDICC2 requires; pieces separated by a gap
(e.g. a range dropped by
[`medicc2_pieces()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_pieces.md))
are never merged, so no segment spans a range without states. The normal
clone `normal` is a sample (MEDICC2 `--normal-name N0`); when absent it
is added as 1:1 on diploid and 1:0 on haploid pieces (male chrX), so it
is the root of the tree. Chromosomes are sorted naturally.

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

  maximal copy number (cn_a, cn_b and their total; MEDICC2's default
  maximum is 8). Clone states can exceed it (the state grid of a clone
  at baseline B reaches 2 B + 1 copies, amplicons more): clipped values
  are counted in the attribute `n_capped` and a warning says so; raise
  `cap` (and MEDICC2's maximum) to keep them

- chroms:

  chromosomes kept

## Value

data.table (MEDICC2 TSV layout); attributes `n_pieces` (pieces before
merging) and `n_capped` (sample x piece rows clipped in the table
written: cn_a or cn_b above `cap` (allele-specific), cn_a + cn_b above
it (total))
