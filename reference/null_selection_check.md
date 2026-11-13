# Selection-free check of the per-cell null

The null cells of normal mode are selected for having no chromosome /
arm call
([`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)),
the thresholds are trained on their rows and the budget is verified on
the same rows by fold cross-validation, so a selection that truncates
the null tail (a too thin parametric null, cells with a noisier than
modelled genome) goes unseen. This check splits the autosomes into two
halves (odd / even chromosome number), selects the null cells on the
rows of one half only, trains the per-cell thresholds as usual on those
rows
([`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
with the selected cells as null cells and the half as `null_chroms`:
depth t null, rho, fold-cross-validated budget thresholds, the run's
`calib` options), and counts the exceedances (`s_comb` above its
threshold) on the null regions of the OTHER half of the same cells: rows
that neither the selection nor the thresholds saw. The two halves are
swapped and the counts combined.

## Usage

``` r
null_selection_check(
  feat,
  params = caller_defaults(),
  multiplets = NULL,
  regions = NULL,
  depth_model = c("t", "gauss"),
  select_args = list(),
  null_chroms = paste0("chr", 1:22),
  false_calls = NULL,
  conf = 0.95
)
```

## Arguments

- feat:

  features long table (as for
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md))

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- multiplets:

  cells never selected as null cells

- regions:

  regions table (default
  [`regions_from_features()`](https://ljwharbers.github.io/scphasedcn/reference/regions_from_features.md))

- depth_model:

  "t" or "gauss"

- select_args:

  extra arguments for
  [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
  (its `chroms` is set per half)

- null_chroms:

  the chromosomes whose rows can be null rows in the run; their numbered
  autosomes are split

- false_calls:

  optional
  [`null_false_calls()`](https://ljwharbers.github.io/scphasedcn/reference/null_false_calls.md)
  `rates` of the run, added as `cv_per_1000` (the run's in-sample
  cross-validated rate on all null rows)

- conf:

  confidence level of the intervals

## Value

data.table class, stratum (`all` = over the strata), n_rows (measured
rows), n_exceed, per_1000 (null and real-event exceedances per 1,000
measured rows: an upper bound of the null rate), lo_per_1000,
hi_per_1000, nominal_per_1000, insample_rows, insample_exceed,
insample_per_1000 (the cross-validated rate on the selection half's null
rows of the same passes), cv_per_1000 (with `false_calls`); attribute
`halves` (half, select_chroms, measure_chroms, n_cells, n_null_cells)

## Details

Cells are selected on one half only, so a cell that carries a real event
on the measured half stays in and the rows of that event count as
exceedances: `per_1000` is the null exceedance rate plus the rate of
real events on the measured half, an upper bound of the null rate. A
rate clearly above the nominal `1 - null_q` (and above the in-sample
cross-validated rate on the selection half, `insample_per_1000`) has two
possible causes: a selected null that is cleaner than the cells it is
applied to (the per-cell false load is then underestimated), or real
events in the cells. Read it next to the share of cells with a called
event: in a sample without real events the excess measures the
selection, in a sample with them it is an upper bound. The halves carry
different chromosomes, so a class whose regions differ in their null
spread (e.g. chromosomes of very different size) also shows up.
Confidence intervals are exact binomial (Clopper-Pearson) over rows;
rows of one cell overlap (a chromosome and its arms), so they are
somewhat too narrow. Costs about one more
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
scoring pass (two half passes without pools). With `calib_rho = "chrom"`
the measured half has no null rows and takes the global rho.
