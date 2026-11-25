# Changelog

## scphasedcn (development version)

- **Tumour mode fits the clone baseline (ploidy) by default**
  ([`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md)
  `ploidy = "fit"`). Up to now the median depth of every clone was taken
  to be 2 copies, so a near-triploid or whole-genome-doubled clone got
  states on a compressed scale: 2:2 as 1:1, 2:0 as 1:0, 3:0 as 2:0, and
  depth-neutral 2:1 often as 1:1. `fit_clone_states(baseline = "fit")`
  now fits per clone the copy number of its median depth (the baseline:
  1.5 to 5) from depth and alleles together, so states are absolute copy
  numbers.
  - The fit compares baselines on one likelihood scale. Above a baseline
    of 2 the state grid drops the near-balanced states such as 4:3,
    which would only absorb orientation bias.
  - A clone whose allele attenuation is at the grid floor carries no
    allelic information on the scale. It takes 2 or the sample baseline
    (the cell-weighted median of the informative clones).
  - Half the best baseline is preferred within `ploidy_tol`, because a
    whole-genome doubling of the same profile fits equally well.
  - Labels stay relative to the diploid normal: 2:1 is a gain at any
    baseline, 1:1 neutral.
  - `clones` gains `baseline` and `ploidy` (the mean autosomal copy
    number), and the run notes list the baselines.
  - A clone fitted at baseline 2 has exactly the states of 0.1.x.
    `ploidy = 2` restores the old behaviour for every clone; a number
    fixes the baseline.
  - On COLO829 (bulk ploidy 3.1) the clones fit baseline 3. Near-diploid
    tumours stay at 2.
- **Normal mode: windowed, orientation-free allele term**
  ([`window_allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/window_allelic_features.md),
  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  `win_*`). With chromosome-scale (hybrid) phasing the per-cell allele
  term summed oriented votes over the whole region, so an LOH carried by
  one cell only was cancelled at population-phase switches: a cell with
  several whole-chromosome losses and CN-LOH could go uncalled and even
  become a null cell.
  [`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
  now also sums, per cell and region, the mirrored beta-binomial LLR of
  near-complete LOH (f = 0.98) over 1 Mb windows x phase sets, with the
  phase trusted only inside a window. These are the new feature columns
  `wal_f980` and `n_win`. Windows with a shared allelic bias over all
  cells are left out, and only autosomes are scored.
  - [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
    scores the per-cell loss / CN-LOH rows with the better of the
    oriented and the windowed LLR where the windowed LLR is above
    `win_min` (15) (`win_score = "max"`, under the same cross-validated
    thresholds). Gains and pooled scores are unchanged.
  - [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)
    never takes a cell whose windowed LLR on a chromosome or arm row
    exceeds the `exp(-t)` bound at the selection quantile (about 9.2;
    `null_win_guard`).
  - **How the windowed calls are controlled.** The guard bound is below
    `win_min`. The windowed term therefore never raises a chromosome or
    arm row of a selected null cell, and the null false-call rates per
    1,000 are those of the oriented scores: they cannot show false calls
    of the windowed term. Those calls are controlled by `win_min` and
    the model bound `P(W > 15) <= 3e-7` per row. The selection-free null
    check (`call_scn(null_check = TRUE)`) counts them.
  - On a normal thyroid sample (4,105 cells), 2 of 4.0 million
    cell-region rows have W \> 15 (two 5-10 Mb tiles, none on a
    chromosome or arm; 5 and 10 Mb tiles are not called per cell). The
    cell calls number 132 rows, as without the term, and the null rates
    and the null check stay where they were (1.62 / 1.78 per 1,000,
    chromosome / arm rows).
  - On a COLO829BL (lymphoblastoid) sample with hybrid phasing, a cell
    with monosomies 1/3/9/12/21 and CN-LOH 8/18 now gets exactly those
    calls. Of 484 LOH calls, 429 (89 %) lie on ASCAT.sc-HP LOH segments.
    - The null false-call rates stay at about 1 per 1,000.
    - The selection-free null check rises from 2.16 to 2.77 (chromosome
      rows) and from 2.88 to 3.30 (arm rows). This comes from
      window-supported LOH rows of cells selected on the other half of
      the autosomes (85 of 112 on ASCAT.sc-HP LOH segments). Exceedances
      without window support do not rise (230 -\> 221).
    - The guard takes 111 cells out of the null set, which also makes
      the depth null tighter: +39 depth-only rows, mostly male chrX (69
      -\> 105 chrX rows).
  - Under block-only phasing the term is not computed (noted), and
    features computed before this version run as before.
  - New simulator event type `"nucleus:<i>"`, which puts several events
    in one nucleus.
- Tumour mode, `gain_support`
  ([`segment_clone()`](https://ljwharbers.github.io/scphasedcn/reference/segment_clone.md),
  [`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md)):
  an imbalanced gain state needs depth support against the region’s
  offset in the other clones. A run without it becomes CN-LOH when its
  allelic excess supports that, otherwise it keeps its state with
  `label` `imbalance` and a new `depth_ok` column FALSE. MEDICC2 input
  treats such runs as 1:1. It addresses clones that mix LOH and
  heterozygous cells being fitted with false gain states (issue
  [\#4](https://github.com/ljwharbers/scphasedcn/issues/4)).
- **Tumour-mode defaults changed to coarse, reproducible clones.**
  `disc_allele = "fraction"`, `test_allele = "centred"`,
  `drift_test = "modes_event"`, `same_post = "group"` and
  `gain_support = 0.3` are now the defaults of
  [`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md)
  and `caller.yaml`. On the lab’s tumour benchmarks this gives fewer
  clones and doublet clones and much better seed agreement. Allele-only
  subclones smaller than the cell-level haplotype drift are merged. Set
  the five options back to `"z"`, `"z"`, `"none"`, `"clone"` and 0 for
  the behaviour of 0.1.1. Normal mode is unchanged.

## scphasedcn 0.1.1

- Extractor (Python): the `read_id` of a read whose name does not follow
  the SPLONGGET convention (`<CB>_#<uuid>_<sign><piece>of<npieces>`) is
  now a stable 64-bit blake2b digest of the name. It was Python’s
  `hash()`, which is salted per process, so for BAMs from other
  pipelines the `read_id` column of `fragments.tsv.gz` changed from run
  to run. Reads named in the SPLONGGET convention are unaffected, and so
  are depth, alleles and the caller.
- No changes to the R package.

## scphasedcn 0.1.0

First public release of scphasedcn.

- [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md):
  copy-number calls per cell and region from an extraction directory
  (`scphasedcn run`), combining read depth and read-phased alleles.
  - Normal-tissue mode
    ([`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)):
    null cells selected from the data, Student-t / beta-binomial nulls,
    fold-cross-validated thresholds at 1 false call per 1,000 null
    cell-regions, per-cell and pooled calls, events.
  - Tumour mode
    ([`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)):
    clone discovery, allele-specific clone copy-number states,
    likelihood assignment of every cell, per-cell inherited and
    confirmed calls, MEDICC2 input
    ([`medicc2_pieces()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_pieces.md),
    [`medicc2_table()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_table.md)).
- [`write_scphasedcn()`](https://ljwharbers.github.io/scphasedcn/reference/write_scphasedcn.md)
  /
  [`read_scphasedcn()`](https://ljwharbers.github.io/scphasedcn/reference/read_scphasedcn.md):
  flat result tables, the parameters used (`params.yaml`) and run notes.
- Command line:
  `system.file("scripts", "call_scn.R", package = "scphasedcn")`, on
  `--dir` + `--out`, or on a sample yaml
  (`inst/extdata/sample_template.yaml`) with `--config` + `--sample`.
- [`make_cell_meta()`](https://ljwharbers.github.io/scphasedcn/reference/make_cell_meta.md):
  the cell metadata table (cell types, known copy-number-flat cells,
  groups) from the extractor’s `cells_qc.tsv` and your own labels.
- The hg38 centromere table ships with the package
  ([`hg38_centromeres()`](https://ljwharbers.github.io/scphasedcn/reference/hg38_centromeres.md));
  the chrX centromere is masked like the autosomal ones.
- [`simulate_scn()`](https://ljwharbers.github.io/scphasedcn/reference/simulate_scn.md)
  /
  [`write_sim_dataset()`](https://ljwharbers.github.io/scphasedcn/reference/write_sim_dataset.md):
  simulated SPLONGGET-like data in the extractor’s layout.
- Vignettes: getting started, inputs (phasing recipe, extraction, cell
  metadata), outputs.
