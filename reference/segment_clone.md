# Segment clone profiles (Viterbi)

Viterbi per clone and chromosome over the states of
[`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
(stay probability `1 - p_switch`, switches weighted by the grid prior);
runs shorter than `min_win` windows are absorbed into the adjacent run
whose state explains them better. With the allele part of the emissions
(`ll_a`, from
[`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md))
a run in an imbalanced state whose allelic log-likelihood ratio vs 1:1
is not positive is given the best balanced state instead. Clone profiles
of clones under 50 cells are expected on the 5 Mb tile grid
([`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
does that).

## Usage

``` r
segment_clone(
  ll,
  p_switch = 0.001,
  min_win = 3,
  grid = NULL,
  min_llr = 10,
  min_excess = 0.5,
  bal_frac = 0.75,
  min_bal_len = 2e+07,
  max_contra = 0.1,
  gain_support = 0.3
)
```

## Arguments

- ll:

  a
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  result, or a long table clone, chrom, start, end, state, ll (optional
  window, ploidy; then summaries that need the profile are NA)

- p_switch:

  switch probability between adjacent windows

- min_win:

  minimal windows per segment

- grid:

  state grid (default the fit's grid or
  [`cn_grid()`](https://ljwharbers.github.io/scphasedcn/reference/cn_grid.md))

- min_llr:

  segments with llr above this are `confident` (neutral segments always
  are)

- min_excess:

  with a
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  input: an imbalanced run also needs its pooled split-half oriented
  fraction to exceed 0.5 by at least this share of the expected excess
  `a * (f_state - 0.5)` (0 = off)

- bal_frac, min_bal_len:

  with a
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  input: a balanced non-neutral run (2:2, depth only) needs a mean log2
  of at least `bal_frac` of its expected shift `c * log2(cn / base)` and
  `min_bal_len` bp, else 1:1 (0 = off)

- max_contra:

  with a
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  input: an imbalanced gain / loss whose mean log2 has the opposite sign
  of its expected shift by more than this takes the best supported
  copy-neutral state (CN-LOH or 1:1; Inf = off)

- gain_support:

  with a
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  input: an imbalanced gain (diploid window, copy number above the
  window's base) needs a mean log2 of at least this share of its
  expected shift `c * log2(cn / base)`, measured against the region's
  depth offset (per window, the median log2 of the other clones that are
  1:1 there with= 50 cells). A run that fails takes the best CN-LOH
  state if its pooled oriented fraction supports one (asmin_excess);
  otherwise it keeps its state with label "imbalance" and depth_ok FALSE
  (allelic imbalancewithout the extra copy: a clone mixing LOH and
  heterozygous cells, or a genome whose baseline is not 2). 0 = off(no
  depth_ok column; the behaviour up to 0.1.1). Only applies to a
  fit_clone_states() input.

## Value

segments table
([`check_segments()`](https://ljwharbers.github.io/scphasedcn/reference/check_tables.md))
plus `label`, `ploidy` (and `depth_ok` when `gain_support` \> 0)
