# Allele-specific copy-number state grid

nA:nB states relative to a diploid clone baseline: neutral 1:1, loss
1:0, cnloh 2:0, gain 2:1, gain_as 3:0 (allele-specific gain), amp 3:1 /
2:2 / 4:1 / 4:0, del 0:0. `f` is the major-haplotype fraction, `prior` a
log weight used by
[`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
(mixture over states) and by the Viterbi transitions of
[`segment_clone()`](https://ljwharbers.github.io/scphasedcn/reference/segment_clone.md).

## Usage

``` r
cn_grid(max_cn = 5, baseline = NULL, lam = 1)
```

## Arguments

- max_cn:

  maximal total copy number

- baseline:

  NULL (the diploid grid above) or the clone baseline the prior is
  centred on

- lam:

  prior fall-off per copy away from the baseline (with `baseline`)

## Value

data.table state, nA, nB, cn, f, label, prior

## Details

With `baseline` (the copy number of the clone's median depth, e.g. 3 for
a near-triploid clone) the grid holds every state nA \>= nB with nA + nB
\<= `max_cn`, and the prior falls off with the distance of the total
copy number from the baseline: log weight `-lam * |cn - baseline|`,
normalised over the grid (so grids of different baselines compare on one
likelihood scale). Labels stay relative to the diploid normal (see
`state_label()`).
