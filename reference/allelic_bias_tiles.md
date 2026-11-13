# Shared allelic bias per tile over all cells (diagnostic)

The per-tile table of
[`allelic_noisy_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_noisy_tiles.md)
computed over every cell (not only null cells) and every eligible row of
a features object, before any caller run: the vote-weighted haplotype
fraction `f0`, its deviation `d0 = |f0 - 0.5|`, the robust z `z_f0` and
the bias flag `masked_f0` (here also `masked`), next to the
overdispersion columns (`rho`, `z_rho`, `masked_rho`). Over all cells a
frequent clonal event also shows up as a bias; the caller's mask
(`params$mask_f0`,
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md))
uses the null cells only, which does not keep an event out of the mask
when its carriers remain among the null cells
([`allelic_noisy_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_noisy_tiles.md)).

## Usage

``` r
allelic_bias_tiles(features, size = "5Mb", params = caller_defaults())
```

## Arguments

- features:

  `scn_features`
  ([`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md))
  or its features long table

- size:

  tile class

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`mask_min_votes`, `mask_rho_mad`, `mask_f0_mad`, `mask_f0_min`)

## Value

data.table as
[`allelic_noisy_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_noisy_tiles.md),
one row per tile, `masked` = `masked_f0`
