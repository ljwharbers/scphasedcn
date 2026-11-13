# Student-t fit of the null depth z per group

Location-scale Student-t (`MASS::fitdistr(x, "t")`) of the finite null z
per group, on a seeded subsample of at most `n_max` rows per group;
fallback median / mad / df 30 when the fit fails (the scale floored at
the sd and at 1e-3).

## Usage

``` r
fit_depth_null(d, by = c("class", "stratum"), n_max = 20000, seed = 1)
```

## Arguments

- d:

  null rows with `z` and the `by` columns

- by:

  grouping columns (region class x depth stratum)

- n_max:

  maximal rows per group used for the fit

- seed:

  RNG seed (re-set per group)

## Value

data.table with the `by` columns, `t_m`, `t_s`, `t_df`, `n` (finite null
rows)
