# Input null of the tumour mode: depth t null and allelic rho

The depth Student-t null
([`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md),
completed by class) and the beta-binomial rho
([`fit_allelic_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_allelic_null.md)
on the arm rows with votes) that
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
uses for discovery, assignment and per-cell confirmation, fitted on the
null rows (`null_region` and `eligible`) of a set of cells:

- `"flagged"`: the cells with `null_cell` TRUE (the `flat` column of the
  cell metadata, or given null cells), whenever there is one and
  `select` is FALSE;

- `"selected"`: otherwise, the cells without a call on the null rows, as
  [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
  finds them (chromosome and arm rows, the tumour state set, `n_iter`
  iterations from all cells minus `exclude` and minus the cells of its
  windowed allele guard, `params$null_win_guard`): the cells kept by
  both of its last two iterations. They are used when there are at least
  `min_cells` of them and the last iteration is stable (Jaccard index of
  the sets it started from and ended with \>= `min_jaccard`). In a
  tumour sample with normal cells that the null rows separate from the
  tumour, these are mostly the normal cells; tumour cells whose events
  all lie outside the null rows or are too weak to call per cell pass as
  well; without such cells the selection usually does not settle;

- `"all"`: otherwise all cells (no stable selection). The fit then
  absorbs the tumour's events on the null rows (inflated rho, heavy
  depth tails), as every tumour-mode run without flagged cells did
  before.

## Usage

``` r
tumour_null_fit(
  long,
  params = caller_defaults(),
  exclude = NULL,
  select = FALSE,
  n_iter = 5,
  min_cells = 5,
  min_jaccard = 0.7,
  max_share = 2/3
)
```

## Arguments

- long:

  features long table
  ([`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
  `long`: cell, region, class, sd_class, chrom, stratum, z, sd_null, nA,
  n_hap, null_region, eligible, null_cell)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- exclude:

  cells never selected (e.g. multiplets)

- select:

  TRUE: select the cells even when some are flagged

- n_iter:

  iterations of the selection (\>= 2)

- min_cells:

  minimal number of selected cells

- min_jaccard:

  minimal Jaccard index of the selection's last iteration

- max_share:

  share of the candidate cells above which a selected set gets a
  `caution` (2 / 3); the fit is not changed

## Value

NULL without null rows (of the flagged cells); else list(tpar, rho,
source = "flagged" / "selected" / "all", null_cells (the cells fitted),
selection (the
[`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
log; NULL for flagged cells), share (selected cells / candidate cells;
NA for flagged cells), caution (NULL, or one line when a selected set
holds more than `max_share` of the candidates), note (one line for the
run notes))

## Details

A fit on all cells is wrong whenever tumour cells carry events on the
null rows: rho and the t tails take up the events, which weakens
assignment and per-cell confirmation. Flagged copy-number-normal cells
(e.g. from a matched assay) remain the better null: a stable selection
does not show that the selected cells are copy-number neutral. In a
tumour in which almost every cell carries the same clonal events on the
null rows, the selection fits its null to those events and then keeps
most of the cells. When it keeps more than `max_share` of the candidate
cells (all cells minus `exclude`, guarded cells included), the result
carries a `caution` and the note says so;
[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
then warns. The note also gives the rho of all cells next to the fitted
one: a selected null whose rho is close to the all-cell rho is a set of
tumour cells that carry the same events.
