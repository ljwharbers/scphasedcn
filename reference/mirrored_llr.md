# Mirrored beta-binomial LLR of one unit

`log(0.5 BB(a; n, f) + 0.5 BB(a; n, 1 - f)) - log BB(a; n, 1/2)`: the
evidence for an imbalance f in either orientation against 1:1. A unit
with n = 1 gives exactly 0.

## Usage

``` r
mirrored_llr(a, n, f, rho, eps = 0.02)
```

## Arguments

- a:

  hap-A votes

- n:

  votes

- f:

  imbalance (major-haplotype fraction)

- rho:

  beta-binomial overdispersion

- eps:

  allele error floor

## Value

numeric vector
