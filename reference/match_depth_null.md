# Per-row Student-t null parameters

Looks up the
[`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md)
parameters of each row's group (default class x stratum), as the scorer
does; stops when a group has no fit. A completed table whose rows carry
an `sd_class` (a class x stratum completed for several sd_classes, e.g.
named events of different sizes in one stratum) gives a row of `d` with
that `sd_class` the fit of its own sd_class.

## Usage

``` r
match_depth_null(d, tpar, by = c("class", "stratum"))
```

## Arguments

- d:

  rows to score (with the `by` columns; `sd_class` used when `tpar` has
  that column)

- tpar:

  output of
  [`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md)
  (completed for the classes without null rows by the callers)

- by:

  grouping columns

## Value

data.table `t_m`, `t_s`, `t_df` with one row per row of `d` (pass to
[`ll_depth_shift()`](https://ljwharbers.github.io/scphasedcn/reference/ll_depth_shift.md))
