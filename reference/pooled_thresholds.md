# Null thresholds of pooled scores

Null pools are built with the same construction as the test pools
(`pools_fun`, e.g. random pools within cell type of null cells, or kNN
pools of null cells), scored per pool x region, and the type-7 `q`
quantile is taken per `by` (class). Random constructions are redrawn
(seed + draw) until every class has at least `n_draws` null
pool-regions; a deterministic construction (same pools twice, e.g. kNN)
is used once. With `cv` one threshold per fold from the null cells of
the other folds, plus fold 0 from all null cells; the fold-0 rows and
scale do not depend on `cv`.
[`call_pools()`](https://ljwharbers.github.io/scphasedcn/reference/call_pools.md)
applies fold 0 only (pools have no fold), so
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
computes fold 0 only (`cv = FALSE`). With `loo` the leave-one-out scores
(`s_rob`) are computed and can be thresholded, exactly as
[`call_pools()`](https://ljwharbers.github.io/scphasedcn/reference/call_pools.md)
scores the test pools.

## Usage

``` r
pooled_thresholds(
  d_null,
  pools_fun,
  n_draws = 20000,
  by = "class",
  q = params$null_q,
  cv = TRUE,
  states = caller_states("normal", params),
  scores = c("s_depth", "s_allelic", "s_comb", if (loo) "s_rob"),
  params = caller_defaults(),
  max_draws = 200,
  loo = FALSE,
  standardize = FALSE,
  min_region_pools = 1000
)
```

## Arguments

- d_null:

  null rows scored with `keep_rows = TRUE` (cell, region, class, fold)

- pools_fun:

  function(cells, seed) returning a pools table (pool, cell)

- n_draws:

  minimal null pool-regions per class

- by:

  grouping of the thresholds; with `sexgrp` (autosome / X / Y,
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
  for males with `pool_sex_thr`), a group with fewer than
  `min_region_pools` null pool-regions takes the class threshold
  (`pooled_fallback`), and redraws are counted per class

- q:

  null quantile

- cv:

  per-fold thresholds

- states:

  state set used for scoring

- scores:

  pooled scores to threshold

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- max_draws:

  maximal redraws

- loo:

  leave-one-out scores (adds `s_rob` to the default `scores`)

- standardize:

  standardise the pooled best LLR per region
  ([`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)
  of the fold's null pools) before thresholding; the fold-0 scale is
  returned as attribute `pool_scale` (pass it to
  [`call_pools()`](https://ljwharbers.github.io/scphasedcn/reference/call_pools.md)).
  With `params$pool_tail_q > 0` the scale carries the per-region tail
  factor (`tail_k = params$pool_tail_k`, `n_eff` = the fold's null cells
  / `params$n_pool`). With `params$pool_tail_xfit` (default) the tail is
  cross-fitted: the fold's null cells are split in two halves (seeded),
  pools are built within each half (`n_draws` / 2 and `min_region_pools`
  / 2 each), each half's pools are standardised by the other half's
  scale before the quantile, and the returned scale is the mean of the
  two half scales (geometric mean of the tails)

- min_region_pools:

  with `standardize`, random constructions are redrawn until every
  region also has this many null pools (a stable median / MAD per
  region; fewer pools make the in-sample standardised threshold
  anti-conservative)

## Value

thresholds table (stratum "all"; `n_null` = null pool-regions); with
`standardize`, attribute `pool_scale`
