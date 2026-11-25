# Per-cell calls table

One row per scored cell x region whose class and best state
(`state_comb`) the policy lists at `level`;
`called = score > max(thr, 0)` (the threshold, and a positive
log-likelihood ratio: a negative null threshold never lets a row be
called while neutral is at least as likely). Rows with
`eligible = FALSE` (regions that are not scored) are left out.

## Usage

``` r
call_table(
  d_scored,
  thr_applied = NULL,
  regions = NULL,
  policy = call_policy("normal"),
  score = "s_comb",
  level = "cell",
  keep_cols = c("z", "f", "n_hap")
)
```

## Arguments

- d_scored:

  scored rows
  ([`score_cells()`](https://ljwharbers.github.io/scphasedcn/reference/score_cells.md));
  with `thr_<score>` columns when `thr_applied` is NULL

- thr_applied:

  a thresholds table to apply
  ([`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md)),
  or NULL when `d_scored` already went through
  [`apply_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/apply_thresholds.md)

- regions:

  regions table (coordinates); NULL leaves them NA

- policy:

  [`call_policy()`](https://ljwharbers.github.io/scphasedcn/reference/call_policy.md)

- score:

  score column

- level:

  policy level

- keep_cols:

  further columns of `d_scored` copied to the calls when present

## Value

calls table (check_calls schema) plus `stratum`, `fold` and `keep_cols`
when present
