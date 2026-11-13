# Apply null thresholds to scored rows

Adds `thr_<score>` and `call_<score>` (`score > thr`, strict). A row
whose `fold` has a threshold in the table gets that fold's
(cross-validated) threshold; other rows (no `fold` column, fold 0, or a
fold the table does not have) get fold 0. A thresholds table whose
`stratum` is all "all" (class-level, e.g. pooled) is matched on class
only. A table with a non-NA `hstratum` column (vote calibration) is also
matched on the rows' `hstratum`. Threshold rows with a `borrowed_from`
sd_class (a class without null rows,
[`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md))
apply to the rows of that class whose `sd_class` is `borrowed_from`.

## Usage

``` r
apply_thresholds(d, thr, scores = unique(thr$score))
```

## Arguments

- d:

  scored rows

- thr:

  [`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md)
  /
  [`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)
  output

- scores:

  score columns to threshold

## Value

copy of `d` with the threshold and call columns
