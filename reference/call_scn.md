# Single-cell copy-number calls of a sample

Normal tissue (`mode = "normal"`):
[`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
with `celltype` (cell-type matched) depth and the oriented allele
channel -\> optional chromatin embedding
([`scn_embedding()`](https://ljwharbers.github.io/scphasedcn/reference/scn_embedding.md))
for kNN pools -\>
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
(null cells, Student-t / beta-binomial nulls, fold-cross-validated
budget thresholds, per-cell calls, pooled calls, events). Tumour
(`mode = "tumour"`):
[`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
with `model` depth, automatic phasing and the pseudobulk groups
`groups`, event chromosomes excluded from the null -\>
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
(clones, clone segments, cell assignment, per-cell inherited / confirmed
calls).

## Usage

``` r
call_scn(
  dir,
  mode = c("normal", "tumour"),
  cell_meta = NULL,
  cells = NULL,
  sex = NULL,
  events = NULL,
  null_cells = NULL,
  multiplets = NULL,
  control_regions = NULL,
  groups = NULL,
  params = caller_defaults(),
  features = NULL,
  embedding = TRUE,
  verbose = TRUE,
  keep_scores = TRUE,
  caller_args = list(),
  ...
)
```

## Arguments

- dir:

  extraction directory (bins_100kb.h5ad, fragments.tsv.gz + .tbi,
  cells_qc.tsv)

- mode:

  "normal" or "tumour"

- cell_meta:

  NULL, a table or a TSV path with `cell` and optional
  `qc_pass, cell_type, flat, fold, group, stratum, n_frag` (passed to
  [`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md):
  without `cells` only its cells that pass
  [`cell_qc()`](https://ljwharbers.github.io/scphasedcn/reference/cell_qc.md)
  are called, and with a `qc_pass` column only those that also have
  `qc_pass` TRUE; `flat` becomes the features' `null_cell`)

- cells:

  NULL (QC-pass cells) or the cells to call (vector or file; no cell
  QC). The barcodes of `cells`, `cell_meta`, `null_cells`, `multiplets`
  and named `groups` lose a trailing `-1` (as the extractor writes
  them); the run notes give how many of them match the data's cells, and
  a warning says when none does

- sex:

  "female" / "male" (case-insensitive, or "F" / "M"; any other value
  stops; required for real data; NULL uses a `sex` column of
  `cells_qc.tsv` if there is one, as in simulated data; its non-missing
  values must agree)

- events:

  named events: a table `region, chrom, start, end` (0-based half-open),
  [`cli_events()`](https://ljwharbers.github.io/scphasedcn/reference/cli_events.md)
  output or a named vector of region strings
  (`20q_gain = "chr20:28000000-"`; unnamed ones are named after the
  region, see
  [`region_set()`](https://ljwharbers.github.io/scphasedcn/reference/region_set.md));
  scored as `event` regions; in tumour mode their chromosomes are
  excluded from the null rows

- null_cells:

  NULL / "auto" (the features' `null_cell` when informative, else
  selected), "features", "select"
  ([`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)),
  or the null cells themselves (vector or file). Tumour mode: given null
  cells (like `flat` cells of `cell_meta`) are the cells whose null rows
  fit the input depth t null and allelic rho when the features carry no
  `null_fit`; they are not the reference pool (the caller finds and
  cleans the normal clone itself); "select" is ignored with a warning

- multiplets:

  NULL or multiplet flags
  ([`multiplet_table()`](https://ljwharbers.github.io/scphasedcn/reference/multiplet_table.md)
  forms); normal mode: cells with prob \> 0.5 are never null cells and
  never pooled; tumour mode: passed to
  [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md),
  which sets aside the cells above its `multiplet_prob` option (default
  0.9). The `multiplet` column of the cells table flags the cells each
  mode set aside

- control_regions:

  tumour mode: control regions for the null (regions table)

- groups:

  tumour mode: pseudobulk groups of the block orientation (a column of
  `cell_meta`, a vector named by cell, or one label per cell); ignored
  with a warning in normal mode

- params:

  [`caller_params()`](https://ljwharbers.github.io/scphasedcn/reference/caller_params.md)
  /
  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md);
  its optional `normal` / `tumour` blocks
  ([`caller_params()`](https://ljwharbers.github.io/scphasedcn/reference/caller_params.md))
  are the default arguments of
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
  /
  [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
  (`opts`)

- features:

  NULL (computed) or a precomputed `scn_features` (object or .rds path)
  to skip the feature step; the arguments of the feature step (`sex`,
  `events`, `...`; in tumour mode also `cells`, `groups`, `cell_meta`;
  in normal mode `cell_meta` columns other than `cell_type`) are then
  not used, with a warning

- embedding:

  TRUE (compute
  [`scn_embedding()`](https://ljwharbers.github.io/scphasedcn/reference/scn_embedding.md)
  from `dir`; skipped with a warning when that fails), FALSE, or a cells
  x components matrix with cell row names

- verbose:

  log progress (also of
  [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md))

- keep_scores:

  keep the per cell x region scores in the result (normal mode; large)

- caller_args:

  extra arguments for
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
  /
  [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
  (e.g. `n_draws`, `opts`); they take precedence over `params$normal` /
  `params$tumour` (tumour: `opts` merged key by key)

- ...:

  extra arguments for
  [`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
  (e.g. `centromeres`, `classes`, `qc`)

## Value

`scphasedcn` with `features` (the long table), `feature_info` (path,
cells, phasing, notes), `scores` (normal, `keep_scores`), `thresholds`,
`calls`, `events`, `pools`, `null_fit`, `false_calls` (normal:
[`null_false_calls()`](https://ljwharbers.github.io/scphasedcn/reference/null_false_calls.md)),
`embedding`, clone components (tumour), `params` (the caller parameters
plus the arguments used: `normal` = every
[`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
argument of the run, `tumour` = the full tumour options) and `notes`
(run time, host, inputs, git revision)
