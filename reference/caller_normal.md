# Normal-tissue copy-number caller

1.  null cells: the features' `null_cell` column, or
    [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
    (`null_cells = "auto"` selects when the column is missing, all NA or
    all TRUE), never `multiplets`; 2. Student-t depth null per class x
    stratum and allelic rho on the null arm rows; 3.
    [`score_cells()`](https://ljwharbers.github.io/scphasedcn/reference/score_cells.md)
    with the normal state set; 4. fold-CV budget thresholds
    ([`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md)); 5.
    per-cell calls by
    [`call_policy()`](https://ljwharbers.github.io/scphasedcn/reference/call_policy.md); 6.
    pools within cell type
    ([`make_pools()`](https://ljwharbers.github.io/scphasedcn/reference/make_pools.md))
    and, with an `embedding`, kNN pools
    ([`knn_pools()`](https://ljwharbers.github.io/scphasedcn/reference/knn_pools.md)),
    pooled null thresholds
    ([`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md))
    and pooled calls; 7. events
    ([`consolidate_calls()`](https://ljwharbers.github.io/scphasedcn/reference/consolidate_calls.md)).
    Rows with `eligible = FALSE` are dropped first.

## Usage

``` r
caller_normal(
  feat,
  params = caller_defaults(),
  cell_meta = NULL,
  multiplets = NULL,
  regions = NULL,
  embedding = NULL,
  null_cells = c("auto", "features", "select"),
  depth_model = c("t", "gauss"),
  n_draws = 20000,
  pool_classes = c("chrom", "arm", "50Mb", "20Mb", "10Mb", "5Mb", "event"),
  select_args = list(),
  sex = NULL,
  resolve_gains = TRUE,
  gain_min_dl = 0,
  gain_gate = "dl",
  gain_gate_q = 0.95,
  pool_loo = TRUE,
  pool_standardize = TRUE,
  null_chroms = paste0("chr", 1:22)
)
```

## Arguments

- feat:

  features long table
  ([`vignette("scphasedcn-outputs", package = "scphasedcn")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-outputs.md))

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- cell_meta:

  optional table `cell` (+ `cell_type`) overriding the features' cell
  types

- multiplets:

  cells flagged as multiplets (never null cells, never pooled)

- regions:

  regions table (default
  [`regions_from_features()`](https://ljwharbers.github.io/scphasedcn/reference/regions_from_features.md))

- embedding:

  optional chromatin embedding (cells x components, row names = cells)
  for kNN pools

- null_cells:

  "auto", "features" or "select"

- depth_model:

  "t" or "gauss"

- n_draws:

  null pool-regions per class for the pooled thresholds

- pool_classes:

  region classes that are pooled

- select_args:

  extra arguments for
  [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)

- sex:

  "female", "male" (case-insensitive, or "F" / "M"; any other value
  stops) or NULL (inferred: female when chrX rows carry allele votes);
  sets the chrX rule of
  [`call_policy()`](https://ljwharbers.github.io/scphasedcn/reference/call_policy.md)

- resolve_gains:

  apply
  [`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md)
  to the per-cell calls

- gain_min_dl:

  `min_dl_gain` of
  [`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md)
  (`gain_gate = "dl"`)

- gain_gate:

  depth gate of the per-cell gains, `gate` of
  [`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md):
  `"dl"` (default: gain depth LLR above `gain_min_dl`, i.e. z above the
  midpoint of neutral and the gain shift) or `"z"` (opt-in: z also above
  the `gain_gate_q` quantile of the depth null of the row's class x
  stratum. About `1 - gain_gate_q` of all null-depth rows pass; among
  rows already called as gains with weak allelic support, and on female
  chrX, up to about two to three times that in simulation. Not
  benchmarked on real data as a default)

- gain_gate_q:

  null quantile of the `"z"` gain gate (0.95)

- pool_loo:

  call pools on the leave-one-out robust score `s_rob`
  ([`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md))

- pool_standardize:

  standardise the pooled best LLR per region
  ([`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)
  of the null pools of the same construction) before the per-class
  pooled thresholds

- null_chroms:

  chromosomes whose rows can be null rows (default the autosomes, the
  chromosomes
  [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
  checks); NULL = every `null_region`

## Value

`scphasedcn` (mode normal) with `scores`, `thresholds` (column `level`:
cell, pool_type, pool_knn), `calls` (cell and pool level; `resolved` for
gains changed by
[`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md);
with `gain_gate = "z"` also its `gate_z`, `gate_thr`), `events`,
`pools`,
`null_fit = list(tpar, rho, null_cells, selection, pool_scale, sex, vote_bins, rho_chrom)`
(`vote_bins`: the
[`vote_bins()`](https://ljwharbers.github.io/scphasedcn/reference/vote_bins.md)
table with `params$calib = "votes"`; `rho_chrom`: the
[`rho_by_chrom()`](https://ljwharbers.github.io/scphasedcn/reference/rho_by_chrom.md)
table with `params$calib_rho = "chrom"`)

## Details

Calibration policy (defaults): null cells are selected on autosomal rows
and null rows are autosomal
([`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
`chroms`, `null_chroms`); per-cell chromosome gains need depth support
and lose to a stronger nested loss / CN-LOH
([`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md),
`resolve_gains`; opt-in `gain_gate = "z"` gates on a null quantile of z
instead of the sign of the gain depth LLR); pooled state scores are
standardised per region by the median / MAD of the null pools of the
same construction
([`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md),
`pool_standardize`), and a pooled call must survive dropping its most
influential cell (`pool_loo`: pools are called on `s_rob`, thresholded
on null pools of the same statistic); chrX pools are not reportable in
females
([`call_policy()`](https://ljwharbers.github.io/scphasedcn/reference/call_policy.md)
`sex`).

Opt-in calibrations: `params$calib = "votes"` keys the per-cell
thresholds on class x stratum x vote stratum (`hstratum`, bins of n_hap
of the null rows,
[`vote_bins()`](https://ljwharbers.github.io/scphasedcn/reference/vote_bins.md));
`params$calib_rho = "chrom"` scores the allele channel with a
per-chromosome rho shrunk to the global
([`rho_by_chrom()`](https://ljwharbers.github.io/scphasedcn/reference/rho_by_chrom.md)).
`params$pool_tail_q > 0` (the default) adds a per-region tail factor to
the pooled region scale
([`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)),
cross-fitted with `params$pool_tail_xfit`
([`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)).
Named `event` regions (no null rows) borrow the per-cell thresholds of
their `sd_class`
([`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md));
their pools are called only with `params$pool_events` (experimental, off
by default:
[`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)).
Removed options stop with an error: `calib` "region" / "region_votes" /
"votes_cont", `calib_rho` "cell" / "chrom_cell".
