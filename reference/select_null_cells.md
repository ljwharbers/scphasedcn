# Self-calibrated null cell set (normal tissue)

Starts from all cells of `feat_long` minus `exclude` (e.g.
multiplet-flagged barcodes) and minus the cells of the windowed guard
(see Details). Each iteration fits the depth t null per class x stratum
and the allelic rho on the null rows of the current set, scores every
candidate cell on the `classes` regions (chromosomes and arms) with
`s_comb`, thresholds at the selection quantile `q` and keeps the cells
with at most `max_calls` calls (it stops early, with no null cell, when
no candidate is left: every candidate guarded, or with more calls).
Every iteration re-evaluates all candidates (a cell dropped earlier can
come back), so the set does not ratchet.

## Usage

``` r
select_null_cells(
  feat_long,
  states = caller_states("normal", params),
  params = caller_defaults(),
  n_iter = 2,
  max_calls = 0,
  exclude = NULL,
  classes = c("chrom", "arm"),
  depth_model = "t",
  threshold = c("parametric", "empirical"),
  q = 0.9999,
  n_sim = NULL,
  chroms = paste0("chr", 1:22)
)
```

## Arguments

- feat_long:

  features long table (cell, region, class, stratum, z, sd_null, nA,
  n_hap; fold, null_region, eligible optional; chrom optional, else
  parsed from the region name)

- states:

  state set
  ([`caller_states()`](https://ljwharbers.github.io/scphasedcn/reference/caller_states.md))

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- n_iter:

  iterations

- max_calls:

  cells with more calls than this are not null cells

- exclude:

  cells never in the null set

- classes:

  region classes used for the selection

- depth_model:

  "t" or "gauss"

- threshold:

  "parametric" or "empirical"

- q:

  selection quantile per row

- n_sim:

  simulated rows per class x stratum (parametric; default max(2e4, 20 /
  (1 - q)))

- chroms:

  chromosomes whose rows are scored (default autosomes chr1-chr22); NULL
  = all

## Value

list(null_cells, previous (the set the last iteration started from), log
(per iteration), calls (per candidate cell: n_calls, win_guard, null),
tpar, rho, thresholds)

## Details

Thresholds (`threshold`): "parametric" (default) re-draws the null rows
under the fitted model (t depth, beta-binomial alleles with
[`robust_rho()`](https://ljwharbers.github.io/scphasedcn/reference/robust_rho.md))
and takes the quantile of the simulated scores; it is not inflated by
the events it is meant to find (an empirical 99.9 % quantile of a set in
which more than 0.1 % of the rows carry events sits inside the events).
"empirical" uses fold-cross-validated quantiles of the current null rows
([`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md)),
which is only safe when the starting set is nearly clean. The default
`q` (1e-4 tail per row) keeps the selection from discarding null cells
for their expected 1 / 1,000 false calls (~6 % of cells at 63
chromosome + arm rows); the calling thresholds of
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
are empirical at `null_q` on the selected set.

Only rows on `chroms` (default the autosomes) are scored: chrX in
females carries X-inactivation driven allelic imbalance and male chrX /
chrY are depth-only, so a chrX call says little about whether a cell's
autosomes are a clean null.

Guard (cells never in the set, whatever their calls): with window
columns in `feat_long`
([`window_allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/window_allelic_features.md))
and `params$null_win_guard` (default TRUE), cells with a scored row
(`classes`, i.e. chromosome and arm rows on `chroms`) whose largest
windowed allele LLR over the fraction grid exceeds
`log(1 / (1 - q)) + log(#fractions)` (about 9.2 at the default `q`;
under the null model `P(W > t) <= exp(-t)`, Bonferroni over the
fractions): an allelic imbalance private to the cell, which the oriented
scores can miss. The bound is per row, without a correction over the
about 60 scored rows of a cell: the guard errs towards leaving a cell
out of the null set, which costs little (6 of 4,105 cells of a normal
thyroid sample). It sees only chromosome and arm rows, so a null cell
can still carry a windowed call on a smaller tile. The selection's
scores themselves stay oriented, so that they match what the parametric
thresholds simulate.
