# Resolve chromosome-level gains against depth and nested calls

A diluted sub-chromosomal imbalance (a 1p loss, a distal CN-LOH) makes
the whole chromosome look like a 2:1 gain: at chromosome scale the
allelic term of the nested event dominates while its depth contrast is
diluted, and the gain state (f 2/3) fits a diluted haplotype fraction
better than 1:0 / 2:0. Per cell, a reportable `gain` call of a `classes`
region therefore needs

1.  depth support, so that the allelic term alone cannot carry a gain
    label:

    - `gate = "dl"` (default): the gain depth LLR of the row itself is
      positive (`dl_gain > min_dl_gain`; for the Student-t null, z above
      the midpoint between neutral and the gain shift, `t_m + mu / 2`
      with `mu = shifts["gain"] / pmax(sd_null, sd_floor)`). For a noisy
      cell (large `sd_null`) `mu` is small and the gate is weak;

    - `gate = "z"` (opt-in): the row's z exceeds
      `max(t_m + mu / 2, t_m + t_s * qt(gate_q, t_df))` with the
      Student-t null of the row's class x stratum (`tpar`,
      [`match_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/match_depth_null.md);
      with `tpar = NULL` the Gaussian null,
      `max(mu / 2, qnorm(gate_q))`): the midpoint stays as a floor and
      the null quantile is added on top. Over all null-depth rows the
      gate passes about `1 - gate_q` of them whatever their `sd_null`;
      this level is unconditional. Among rows already called as gains
      (weak allelic support), and on female chrX (no null rows of its
      own), the pass rate of null-depth rows is higher: up to about two
      to three times `1 - gate_q` in simulation. The gate is opt-in and
      has not been benchmarked on real data as a default. A row with NA
      z has no depth support. `min_dl_gain` is not used. The gate is on
      z, not on `dl_gain`: the Student-t LLR falls back towards 0 for
      very large z
      ([`ll_depth_shift()`](https://ljwharbers.github.io/scphasedcn/reference/ll_depth_shift.md)),
      so a quantile on `dl_gain` would reject the strongest gains; and

2.  (class `chrom`) no nested explanation: no called loss / CN-LOH row
    of the same cell and chromosome in `nested` classes (arm, 50 Mb, 20
    Mb) with a higher score. A gain failing either rule stays in the
    table with `called = TRUE` but `reportable = FALSE` and `resolved`
    naming the reason (`no_depth`, `nested <region> <state>`). With
    `promote`, the best (highest-score) nested loss / CN-LOH call that
    explains a dropped chromosome gain is made reportable
    (`resolved = "replaces <region> gain"`; arm and 50 Mb losses are
    reportable anyway, a 20 Mb one is promoted). Losses are not gated:
    the loss vs CN-LOH label is decided by the scores
    ([`score_cells()`](https://ljwharbers.github.io/scphasedcn/reference/score_cells.md)).

## Usage

``` r
resolve_nested(
  calls,
  scores = NULL,
  min_dl_gain = 0,
  classes = c("chrom", "event"),
  nested = c("arm", "50Mb", "20Mb"),
  promote = TRUE,
  gate = c("dl", "z"),
  gate_q = 0.95,
  tpar = NULL,
  params = caller_defaults()
)
```

## Arguments

- calls:

  cell-level calls table
  ([`call_table()`](https://ljwharbers.github.io/scphasedcn/reference/call_table.md))

- scores:

  scored rows (`cell`, `region`, `dl_gain`; with `gate = "z"` also `z`,
  `sd_null` and, with `tpar`, `class`, `stratum`); NULL when `calls` has
  these columns

- min_dl_gain:

  minimal gain depth LLR of a reportable gain (`gate = "dl"`)

- classes:

  classes whose reportable gains need depth support

- nested:

  classes of the nested calls that can explain a chromosome gain

- promote:

  make the explaining nested call reportable

- gate:

  depth gate of a gain: `"dl"` (gain depth LLR above `min_dl_gain`) or
  `"z"` (z above the midpoint and above the `gate_q` quantile of the
  row's depth null)

- gate_q:

  null quantile of the `"z"` gate (0.95: about 5 % of all null-depth
  rows pass; more among rows already called as gains and on female chrX,
  see Details)

- tpar:

  `"z"` gate: the Student-t fits per class x stratum
  ([`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md),
  completed as in
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md));
  NULL = the Gaussian depth null N(0, 1). A gated row whose class x
  stratum has no fit (a named `event`: no null rows) takes the fit of
  its `sd_class` (a `scores` column) in that stratum; without one it is
  an error
  ([`match_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/match_depth_null.md));
  a row with NA class or stratum has no depth support

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`shifts`, `sd_floor`; the `"z"` gate)

## Value

`calls` with column `resolved` (NA = unchanged) and updated
`reportable`; with `gate = "z"` also `gate_z` (the z the gate used) and
`gate_thr` (its z threshold), set on the gated gain rows and NA
elsewhere
