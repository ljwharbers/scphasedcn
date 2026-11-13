# Null thresholds at the false-call budget

Type-7 quantile `q` of each score on the null rows
([`budget_threshold()`](https://ljwharbers.github.io/scphasedcn/reference/budget_threshold.md)
compatible), per `by` group (class x stratum). With `cv = TRUE` one
threshold per fold k, trained on the null rows of all other folds
(applied to fold k), plus fold 0 trained on all null rows. A group with
fewer than `min_null` null rows uses the group pooled over strata (the
class-wide threshold), flagged `pooled_fallback`.

## Usage

``` r
budget_thresholds(
  d,
  scores = c("s_depth", "s_allelic", "s_comb"),
  by = c("class", "stratum", if ("hstratum" %in% names(d)) "hstratum"),
  q = params$null_q,
  cv = TRUE,
  min_null = params$min_null_rows,
  params = caller_defaults(),
  null = NULL
)
```

## Arguments

- d:

  scored rows with the `by` columns, `fold` (cv) and a null definition

- scores:

  score columns

- by:

  grouping columns (default class x stratum, plus `hstratum` when `d`
  has that column)

- q:

  null quantile

- cv:

  per-fold thresholds

- min_null:

  minimal null rows of a group before falling back

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- null:

  optional logical null flag per row (default `is_null`, else
  `null_cell & null_region & eligible`)

## Value

thresholds table (check_thresholds schema; `stratum` = "all" when not in
`by`; `borrowed_from` when a class borrows its sd_class's thresholds)

## Details

With a vote stratum (`hstratum` column,
[`assign_hstratum()`](https://ljwharbers.github.io/scphasedcn/reference/assign_hstratum.md);
`calib = "votes"`) the groups are class x stratum x hstratum, and a
group short of `min_null` null rows falls back to its class x stratum,
then to its class (the first with `min_null` rows, else the class).

A class without any null row whose rows have an `sd_class` with null
rows (named `event` regions: never null rows) takes the thresholds of
that `sd_class`: its groups are keyed by class x sd_class x the other
`by` columns and use the null rows (and the fallback levels) of the
`sd_class`, so every fold / stratum gets the same threshold as the
`sd_class` rows of that group. Those threshold rows carry
`borrowed_from` (the sd_class; the column is only added when a class
borrows) and
[`apply_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/apply_thresholds.md)
matches them on the rows' `sd_class`.
