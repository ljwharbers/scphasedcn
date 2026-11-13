# Oriented allelic log-likelihood ratio vs balanced

`bb_loglik(nA, n, f, rho, eps) - bb_loglik(nA, n, 0.5, rho, eps)` for
molecule votes `nA` of `n` (one vote per molecule). Two-sided (default)
takes the better of `f` and `1 - f` (the orientation of the counts is
unknown or may be wrong). With several `f` the best LLR over the grid is
returned (the attenuation-aware 1:0 / 2:0 term over `f_att`). Rows with
`n = 0` (or NA) give 0; without any vote `rho` is not used (it may be
NA). The scorer uses llr_20 at `f = 1 - eps`, llr_21 at 2/3, llr_31 at
3/4 and the grid form over `f_att`.

## Usage

``` r
ll_allelic(nA, n, f, rho, eps = 0.02, two_sided = TRUE, per_row = FALSE)
```

## Arguments

- nA:

  votes for haplotype A

- n:

  total votes

- f:

  major-haplotype fraction of the state (scalar or a grid)

- rho:

  beta-binomial overdispersion (null fit,
  [`fit_allelic_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_allelic_null.md))

- eps:

  allele error floor

- two_sided:

  max over f and 1 - f

- per_row:

  `f` is one expected fraction per row (length 1 or `length(n)`), not a
  grid: the LLR of each row at its own `f` (NA results also give 0).
  Used by the tumour mode (`0.5 + a * (f_state - 0.5)` per row).

## Value

numeric vector (length of n)
