# Per-cell copy-number states on the MEDICC2 pieces

For each cell x piece, the candidates of
[`cell_piece_candidates()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_candidates.md)
(1:1 plus every oriented `(cn_a, cn_b)` a kept clone has on the piece;
haploid pieces: one-copy states) are scored with the cell's own evidence
([`cell_piece_data()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_data.md));
the cell's clone assignment is not used.

- Depth: `ll_depth_shift(z, mu, tpar)` with
  `mu = c * log2(max(cn, 0.2) / base) / pmax(sd, sd_floor)` (the shift
  model of the clone states,
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  /
  [`assign_cells()`](https://ljwharbers.github.io/scphasedcn/reference/assign_cells.md):
  `base` = the piece ploidy when the depth is centred on the reference
  pool, else 2), relative to the neutral state's shift (the neutral
  state's own depth term is `ll_depth_shift(z, 0) = 0`), Student-t null
  per sd_class x stratum
  ([`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md)
  parameters, fallbacks as in
  [`assign_cells()`](https://ljwharbers.github.io/scphasedcn/reference/assign_cells.md))
  or Gaussian when `tpar` is NULL.

- Alleles: one-sided beta-binomial LLR vs 1:1 of the oriented votes
  `(oA, oN)` at the expected cn_a fraction `0.5 + a * (cn_a / cn - 0.5)`
  (`a` = orientation attenuation of the clone model; 0.5 for cn = 0),
  overdispersion `rho`; haploid pieces have no allele term.

- MAP over the candidates with the prior (`neutral`, `uniform` or
  `clone_freq`); `posterior` is the MAP's posterior and `llr_neutral`
  the log-likelihood (no prior) of the MAP minus that of the neutral
  state.

## Usage

``` r
cell_piece_states(
  features,
  pieces,
  blocks = NULL,
  cells = NULL,
  params = caller_defaults(),
  prior = c("neutral", "uniform", "clone_freq"),
  p_neutral = 0.5,
  clones = NULL,
  c = NULL,
  a = NULL,
  rho = "ref",
  tpar = NULL,
  ref_cells = NULL,
  loo = NULL,
  normal = "N0",
  sd_scale = TRUE,
  ref_center = c("pseudobulk", "mean"),
  center = TRUE,
  s_sys = NULL,
  haploid_sd = sqrt(2),
  data = NULL,
  keep_ll = FALSE,
  rho_grid = c(0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25)
)
```

## Arguments

- features:

  scn_features list (see
  [`cell_piece_data()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_data.md))

- pieces:

  clone x piece table with `cn_a, cn_b` (orient_clone_states()\$pieces
  or the export's pieces.tsv), or a `medicc2_orientation`

- blocks:

  block orientation (orient_clone_states()\$blocks); NULL takes it from
  a `medicc2_orientation`

- cells:

  cells to score (default all)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`sd_floor`, `eps`, `tile_sizes`)

- prior:

  `neutral` (default), `uniform` or `clone_freq` (see
  [`cell_piece_candidates()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_candidates.md)).
  Under the uniform prior a short piece with many candidates and little
  evidence goes to whichever non-neutral state fits marginally better,
  so many rows of normal cells get non-neutral states; `clone_freq`
  gives uninformative cells the majority clone state (tumour states in
  normal cells)

- p_neutral:

  neutral prior mass of `prior = "neutral"`

- clones:

  clones table (`clone, c, a, n_assigned`); used for the default `c` /
  `a` and the `clone_freq` prior

- c, a:

  global depth contrast and allelic attenuation (NULL = cell-weighted
  median of the tumour clones; 1 when `clones` lacks them)

- rho:

  beta-binomial overdispersion of a cell's oriented votes on a piece: a
  number; `"ref"` = fitted
  ([`fit_bb_rho()`](https://ljwharbers.github.io/scphasedcn/reference/fit_bb_rho.md),
  f = 0.5) on the reference cells' oriented votes (diploid pieces; the
  reference pool is balanced everywhere, so it is the null of the allele
  channel, as it is the per-cell null of the caller); `"auto"` = the
  value of `rho_grid` that maximises the marginal likelihood of all
  scored cell x piece rows over their candidates (with the prior); NULL
  = the features' / run's allelic null rho. The run's null rho is fitted
  on votes oriented by the all-cell pseudobulk, which on event-free
  regions is a random orientation; in the true haplotype frame of the
  block orientation a single cell's votes on a piece are far more
  dispersed (often by an order of magnitude).

- tpar:

  depth null t parameters (class, stratum, t_m, t_s, t_df) or NULL
  (Gaussian)

- ref_cells:

  reference pool for the depth centring

- loo:

  NULL or data.table cell, clone (see
  [`cell_piece_data()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_data.md))

- normal:

  name of the normal clone (excluded from the default c / a)

- sd_scale, ref_center, center, haploid_sd:

  see
  [`cell_piece_data()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_data.md)

- s_sys:

  systematic depth SD added to the per-cell SD (NULL = the cell-weighted
  median of the tumour clones' `s_sys`, 0 when `clones` lacks it); see
  [`cell_piece_data()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_data.md)

- data:

  precomputed
  [`cell_piece_data()`](https://ljwharbers.github.io/scphasedcn/reference/cell_piece_data.md)
  output (skips the depth / vote aggregation)

- keep_ll:

  return the per-candidate log-likelihoods

- rho_grid:

  grid of `rho = "auto"`

## Value

list of class `cell_piece_states`: `states` (one row per cell x piece:
cell, piece, chrom, start, end, ploidy, cn_a, cn_b, state (MAP),
posterior, llr_neutral, ll_depth, ll_allele (MAP terms vs neutral),
n_cand, z, log2ratio, offset, sd, n_win, oA, oN (oriented votes), f_obs
(oA / oN), n_blocks), `candidates`, `ll` (cell x piece x candidate, when
`keep_ll`), `params` (c, a, rho, rho_mode, rho_rows (reference rows of
the rho fit), rho_ml (the `auto` grid), s_sys, ref_center, center,
prior, n_ref, loo, sd_scale, depth_null)

## Details

`c`, `a` and `s_sys` are global (one value for every cell): by default
the median over the tumour clones (kept clones of `pieces` except
`normal`), weighted by their assigned cells, of the clones' fitted `c` /
`a` / `s_sys` (`clones` table). A per-cell value would have to come from
the cell's clone, which this function must not use.
