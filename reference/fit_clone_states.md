# Fit clone copy-number states

Per clone, a grid search over the depth contrast `c` (state mean
`c * log2(cn / base)`, `base` = the copy number of the clone's median
depth: the clone baseline B, times the reference ploidy / 2 for
pool-relative profiles; haploid windows (male chrX / chrY, `ploidy` 1)
allow only nB = 0 states and their labels are relative to one copy), the
orientation attenuation `a` (expected major fraction
`0.5 + a * (f_state - 0.5)`) and the window-level overdispersion `rho`.
Emissions: Student-t (`df`) on the window log2 with scale = the window
SE combined with the clone's systematic window noise and inflated for
its autocorrelation (see `clone_noise`), and a two-sided robust
beta-binomial on (nA, n_hap) (a 5 % component with 10x the
overdispersion). The objective is the per-window mixture likelihood over
the states with the grid `prior` weights.

## Usage

``` r
fit_clone_states(
  profile,
  grid = cn_grid(),
  contrast = seq(0.6, 1, 0.05),
  atten = seq(0.6, 1, 0.05),
  rho = c(0.001, 0.003, 0.01, 0.03),
  df = 5,
  eps = 0.02,
  refine = TRUE,
  p_switch = 0.001,
  min_win = 3,
  baseline = 2,
  baselines = c(1.5, 2, 2.5, 3, 3.25, 3.5, 4, 4.5, 5),
  lam = 1,
  min_df = 0.1,
  ploidy_tol = 0.01
)
```

## Arguments

- profile:

  output of
  [`clone_profiles()`](https://ljwharbers.github.io/scphasedcn/reference/clone_profiles.md)

- grid:

  state grid
  ([`cn_grid()`](https://ljwharbers.github.io/scphasedcn/reference/cn_grid.md))

- contrast:

  grid of depth contrasts

- atten:

  grid of orientation attenuations

- rho:

  grid of window-level overdispersions

- df:

  degrees of freedom of the depth t

- eps:

  allele error floor

- refine:

  refine (c, a, rho) on the segment level: segment at the current fit
  (Viterbi, `p_switch`, `min_win`), then maximise the sum over segments
  of the best state's segment likelihood (depth on the
  precision-weighted segment mean, alleles summed over its windows); two
  passes

- p_switch, min_win:

  segmentation parameters of the refinement (as
  [`segment_clone()`](https://ljwharbers.github.io/scphasedcn/reference/segment_clone.md))

- baseline:

  the clone baseline B: one number for every clone (2 = diploid), or
  "fit"

- baselines:

  candidate baselines of `baseline = "fit"`

- lam:

  prior fall-off per copy away from the baseline
  ([`cn_grid()`](https://ljwharbers.github.io/scphasedcn/reference/cn_grid.md))
  of the baseline grids

- min_df:

  imbalanced states of the baseline grids need \|f - 0.5\| \>= this

- ploidy_tol:

  log-likelihood tolerance per window within which half the best
  baseline is preferred

## Value

list(ll = data.table clone, window, chrom, start, end, ploidy, state,
label, ll (emission log-likelihood at the chosen parameters), windows =
the profile with `se_eff` and `base`, par = data.table clone, c, a, rho,
s_sys, loglik, n_win, base (the baseline B), grid (the B = 2 grid),
grids (per clone), baseline_scan (with "fit": clone, B, ll, c, a per
candidate), baseline_sample, baseline_informative (with "fit": TRUE when
at least one clone had allelic information on the scale, i.e. the sample
baseline was fitted, not the default 2), contrast_floor (the lower bound
of `contrast`: a clone whose `c` sits there has a poorly determined
scale, and its states can be off by one copy))

## Details

The baseline B is 2 by default (the diploid baseline of 0.1.x: a
near-triploid clone is then fitted on a compressed scale). With
`baseline = "fit"` every clone is fitted at each of `baselines` on the
grid `baseline_grid(B)` (all states up to 2 B + 1 copies without the
near-balanced imbalanced ones, prior centred on B and normalised, so the
baselines compare on one likelihood scale; a coarse (c, a) grid, no
refinement). A clone whose best attenuation `a` is above the floor of
`atten` has allelic information on the scale and takes its best
baseline, or half of it when that is within `ploidy_tol` x windows (a
whole-genome doubling of the same profile fits as well); the sample
baseline is the cell-weighted median of these clones; the other clones
take the better of 2 and the sample baseline. Each clone is then fitted
at its baseline (B = 2: `grid`, exactly as before).
