# Assign cells to clones by their likelihood under the clone states

Per cell and clone k: the sum over pieces of the depth term
`ll_depth_shift(z, c_k * log2(cn / 2) / pmax(sd_null, sd_floor), tpar)`
and the allele term
`ll_allelic(nA, n, 0.5 + a_k * (f_state - 0.5), rho, two_sided = FALSE)`
on the cell's votes oriented by clone k's pseudobulk (leave-one-out for
k's members). The depth z is taken relative to the reference cells' mean
log2 ratio of each piece (`ref_cells`; removes the shared per-region
bias of the normalisation). With `mixtures`, every pair of clones is
also scored as a mixture (weights `mix_w`; copy numbers mixed; per piece
the votes oriented by the pair's clone that is imbalanced there, as for
the single clones, with the other clone's haplotype agreement taken from
its pooled counts under that orientation). Posterior over singles (prior
`prior`, default uniform) and mixtures (total prior `mix_prior`); the
reported `posterior` is the best clone's posterior among the single
clones.

## Usage

``` r
assign_cells(
  feat_pieces,
  clones,
  segments,
  tpar,
  rho,
  par = NULL,
  prior = NULL,
  mixtures = TRUE,
  min_post = 0.9,
  mix_prior = 0.05,
  ref_cells = NULL,
  mix_w = 0.5,
  mix_margin = 10,
  mix_rel = 0.1,
  params = caller_defaults(),
  mix_pairs = NULL,
  mult_prior = NULL,
  trunks = NULL,
  max_mix_prior = 0.5
)
```

## Arguments

- feat_pieces:

  list(features = per-cell depth features on the pieces
  ([`region_features()`](https://ljwharbers.github.io/scphasedcn/reference/region_features.md)),
  pieces = regions table, votes, atoms, cells = names in vote index
  order)

- clones:

  data.table cell, clone: the clone members used for orientation
  (confident singlets)

- segments:

  segments table

- tpar:

  depth null t parameters
  ([`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md))
  or NULL (Gaussian)

- rho:

  per-cell beta-binomial overdispersion

- par:

  data.table clone, c, a (from
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md));
  default c = a = 1

- prior:

  named clone prior (NULL = uniform; a size-proportional prior lets
  large diffuse clones absorb low-depth cells)

- mixtures:

  score pairwise 50:50 mixtures

- min_post:

  posterior needed for `assigned`

- mix_prior:

  total prior of the mixture hypotheses

- ref_cells:

  reference cells for the depth centring (NULL = none); the clone states
  must be relative to the same reference
  ([`clone_profiles()`](https://ljwharbers.github.io/scphasedcn/reference/clone_profiles.md)
  `ref_cells`)

- mix_w:

  mixture weights of the first clone of a pair (0.5 = 50:50; e.g. c(1/3,
  1/2, 2/3) for unequal nuclei)

- mix_margin:

  a cell is `doublet_like` when the mixtures' posterior is \> 0.5 and
  the best mixture beats the best single clone by more than this
  (log-likelihood)

- mix_rel:

  ... and by more than this fraction of the cell's best single-clone
  log-likelihood ratio (deep cells fit any single clone less perfectly
  in absolute terms)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- mix_pairs:

  optional 2-column matrix of the clone pairs scored as mixtures
  (default all pairs);
  [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
  passes the pairs whose states differ allelically on \>= `mix_min_dist`
  bp
  ([`clone_state_distance()`](https://ljwharbers.github.io/scphasedcn/reference/clone_state_distance.md)):
  a mixture of near-identical (sibling) clones fits a cell between them
  better than either clone, which is not evidence of a doublet

- mult_prior:

  optional named per-cell multiplet prior probability (e.g. a soft
  multiplet score): the cell's total mixture prior is max(`mix_prior`,
  prior), capped at `max_mix_prior`

- trunks:

  optional named clone -\> trunk vector
  ([`clone_trunks()`](https://ljwharbers.github.io/scphasedcn/reference/clone_trunks.md)):
  adds `trunk` and `p_trunk` (posterior summed over the clones of the
  best clone's trunk)

- max_mix_prior:

  cap of the per-cell mixture prior

## Value

data.table cell, clone (the likelihood-best single clone, for every
cell), posterior, ll_margin, doublet_like, pair, ref_pool, p_mix,
ll_mix_gain, ll_best, mix_prior, assigned / confident (posterior \>=
`min_post` and not doublet-like) (plus trunk, p_trunk when `trunks` is
given; attribute `ll`: cells x hypotheses log-likelihood matrix)

## Details

The hypotheses are the clones with members plus the fixed normal clone
`N0` whenever `segments` has it
([`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
`normal_hypothesis`), even when no cell is labelled `N0` yet: its state
is all-neutral, so it needs no orientation, and without it normal cells
would be forced into a tumour clone. Its `c` is the median `c` of `par`
(attenuation 1) when `par` has no `N0` row.
