# Allelically overdispersed tiles (data-driven allelic mask)

Per tile of class `size`, the method-of-moments beta-binomial
overdispersion of the null rows' allele votes,
`rho = (sum (nA - n/2)^2 - sum n/4) / sum n (n - 1) / 4`, then `z_rho` =
(log rho - median) / MAD of log rho over the tiles with at least
`min_votes` votes per null row. Tiles with `z_rho > k_mad` are masked: a
local stretch whose votes are far more overdispersed than the genome
(e.g. a 5 Mb tile with rho 0.040 vs a genome median of 0.011 over 5 Mb
tiles) and that inflates the chromosome-level allelic LLRs of the cells.
The votes of the atoms in masked tiles
([`atoms_in_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/atoms_in_tiles.md))
are then left out of the allele features
([`update_allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/update_allelic_features.md),
`call_scn(params$mask_allelic = TRUE)`).

## Usage

``` r
allelic_noisy_tiles(
  feat_long,
  null = NULL,
  size = "5Mb",
  k_mad = params$mask_rho_mad,
  min_votes = params$mask_min_votes,
  params = caller_defaults(),
  f0_mad = params$mask_f0_mad,
  f0_min = params$mask_f0_min,
  rules = "rho"
)
```

## Arguments

- feat_long:

  features long table (class, region, chrom, nA, n_hap) with a null
  definition

- null:

  optional logical null flag per row

- size:

  tile class

- k_mad:

  threshold in MADs of log rho (`params$mask_rho_mad`)

- min_votes:

  minimal mean votes per null row of a tile for it to be assessed
  (`params$mask_min_votes`)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- f0_mad:

  bias threshold in MADs of d0 (`params$mask_f0_mad`)

- f0_min:

  minimal deviation d0 of a bias outlier (`params$mask_f0_min`)

- rules:

  the rules that define `masked`: any of `"rho"`, `"f0"`

## Value

data.table tile, chrom, start, end, n_rows, votes, rho, z_rho, assessed,
masked, f0, d0, z_f0, masked_rho, masked_f0 (one row per tile)

## Details

Shared allelic bias (`params$mask_f0`): per tile also the vote-weighted
haplotype fraction of the null rows `f0 = sum nA / sum n_hap`, its
deviation `d0 = |f0 - 0.5|` and `z_f0` = (d0 - median) / MAD of d0 over
the assessed tiles. A tile is a bias outlier (`masked_f0`) when it is
assessed, `z_f0 > f0_mad` and `d0 > f0_min`: a stretch whose allele
fraction is off 1:1 in the null cells for reasons that are not somatic
copy number (germline CNV or segmental duplication, imprinted or
allele-specific accessible regions, mapping bias, false heterozygous
SNPs), against which the 1:1 allele null of every cell is wrong. The
rule sees only the null rows, but it cannot tell such a shared bias from
a real clonal event whose carriers remain among the null cells: an event
frequent enough for that looks the same and can be masked (in a
simulation with a focal event in 35 % of cells, 60 of 96 carriers stayed
in the null set and all four event tiles were masked). The mask is a
diagnostic to be read with the list of masked tiles, not a filter to
apply blindly to samples that may carry clonal events. `|f0 - 0.5|`
covers two orientations of the allele features: under `loo` (block)
orientation the fraction is major-oriented (a biased tile sits above
0.5), under chromosome-anchored orientation it sits at f0 or 1 - f0.
Under `mirrored` orientation each block counts its major allele, so f0
is above 0.5 in every tile by an amount set by the votes per block: the
statistic measures coverage, not shared bias, and
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
does not apply the bias rule (it says so in the run notes and warns).
[`allelic_bias_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_bias_tiles.md)
gives the same table over all cells before calling.

`rules` selects which flags make up `masked` (the tiles whose atoms are
dropped): `"rho"` (the default and the behaviour of earlier versions:
`masked` = `masked_rho`), `"f0"`, or both (union).
