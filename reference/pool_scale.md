# Per-region scale of pooled null scores

Median and MAD ([`mad()`](https://rdrr.io/r/stats/mad.html),
normal-consistent) of the pooled best LLR (the max over the states'
`cl_<state>` of
[`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md),
i.e. the raw `s_comb`) over the null pools of each region. Pooled scores
standardised with it, `(s - med) / mad`, put every region of a class on
one scale, so a per-class pooled threshold neither over-calls regions
whose null pools are shifted or wide (molecule-poor regions where the
non-carriers of a pool add little negative LLR, shared depth / phasing
artefacts, cell-type specific depth in kNN neighbourhoods) nor ignores
quiet ones. The state stays the likelihood argmax. A region with fewer
than `min_pools` null pools takes the class medians; MADs are floored at
`mad_floor` x the class median MAD.

## Usage

``` r
pool_scale(
  ps,
  states = caller_states("normal"),
  min_pools = 30,
  mad_floor = 0.1,
  tail_q = 0,
  tail_k = 30,
  n_eff = Inf
)
```

## Arguments

- ps:

  pooled null scores
  ([`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md)
  output with `region`, `class`, `cl_<state>`)

- states:

  state set

- min_pools:

  minimal null pools of a region

- mad_floor:

  floor of a region's MAD relative to its class median

- tail_q:

  quantile of the per-region tail factor (0 = no tail factor)

- tail_k:

  shrinkage pseudo-count of the tail factor (effective pools; 30,
  `caller_defaults()$pool_tail_k`)

- n_eff:

  effective number of independent null pools per region (caps `n_pools`
  in the shrinkage weight)

## Value

data.table region, class, n_pools, med, mad (+ `tail`, `u_tail` when
`tail_q > 0`)

## Details

Tail factor (`tail_q > 0`, off by default): the median / MAD standardise
the body of a region's null pools, but a heavy-tailed region keeps
collecting calls at a per-class threshold. Each region's null pools are
standardised by its `med` / `mad`, `u_r` = their type-7 `tail_q`
quantile and `u_c` = the median of `u_r` over the class's regions with
at least `min_pools` pools. With `ne = min(n_pools, n_eff)` (overlapping
kNN pools: one per cell, so the effective number of independent pools is
~ cells / pool size) and `w = ne / (ne + tail_k)`,
`tail = exp(w log(u_r / u_c))`; `tail = 1` when `u_r` or `u_c` is not
positive or the region has fewer than `min_pools` pools.
[`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md)
/ `pool_combine()` then divide by `mad x tail`.
