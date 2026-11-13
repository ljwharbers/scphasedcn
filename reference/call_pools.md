# Pooled calls

Scores the pools
([`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md)
per pool x region), applies the fold-0 pooled thresholds per class and
returns the policy's pool-level rows. With `loo` the pools are scored as
[`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md)
describes; use `score = "s_rob"` (thresholded by
[`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)
with `loo = TRUE`) so a call must survive dropping its most influential
cell.

## Usage

``` r
call_pools(
  d,
  pools,
  thr,
  policy = call_policy("normal"),
  states = caller_states("normal"),
  regions = NULL,
  score = "s_comb",
  loo = score %in% c("s_rob", "s_loo"),
  scale = NULL
)
```

## Arguments

- d:

  rows scored with `keep_rows = TRUE`

- pools:

  pools table (pool, cell)

- thr:

  [`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)
  output

- policy:

  [`call_policy()`](https://ljwharbers.github.io/scphasedcn/reference/call_policy.md)

- states:

  state set used for scoring

- regions:

  optional regions table (coordinates)

- score:

  pooled score

- loo:

  leave-one-out scores (needed for `score = "s_rob"`)

- scale:

  NULL or the
  [`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)
  table the thresholds were computed with (attribute `pool_scale` of
  [`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)
  with `standardize = TRUE`). A region of a class without null pools
  (named `event` regions) is standardised with the class medians of its
  `sd_class` when `thr` has borrowed rows

## Value

calls table (`level = "pool"`, `cell` = pool id) plus `n` (cells in the
pool), `s_comb` and, with `loo`, `s_loo`; with `scale`, `s_comb_llr`
(unstandardised best LLR)
