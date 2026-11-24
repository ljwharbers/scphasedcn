# scphasedcn (development version)

* **Normal mode: windowed, orientation-free allele term** (`window_allelic_features()`, `caller_defaults()` `win_*`).
  With chromosome-scale (hybrid) phasing the per-cell allele term summed oriented votes over the whole region, so an
  LOH carried by one cell only was cancelled at population-phase switches: a cell with several whole-chromosome
  losses and CN-LOH could go uncalled and even become a null cell. `cell_features()` now also sums, per cell and
  region, the mirrored beta-binomial LLR of near-complete LOH (f = 0.98) over 1 Mb windows x phase sets, with the
  phase trusted only inside a window. These are the new feature columns `wal_f980` and `n_win`. Windows with a shared
  allelic bias over all cells are left out, and only autosomes are scored.
  * `caller_normal()` scores the per-cell loss / CN-LOH rows with the better of the oriented and the windowed LLR
    where the windowed LLR is above `win_min` (15) (`win_score = "max"`, under the same cross-validated thresholds).
    Gains and pooled scores are unchanged.
  * On a COLO829BL (lymphoblastoid) sample with hybrid phasing, a cell with monosomies 1/3/9/12/21 and CN-LOH 8/18
    now gets exactly those calls, and 89 % of the LOH calls lie on ASCAT.sc-HP LOH segments. On a normal thyroid
    sample the calls are practically unchanged (133 vs 132 rows). The null false-call rates are unchanged in both.
  * `event_sd` (experimental, off): rescores already-called cells with a depth scale computed without their called
    chromosomes.
  * `select_null_cells()` never takes a cell whose windowed LLR exceeds the `exp(-t)` bound at the selection quantile
    (`null_win_guard`).
  * Under block-only phasing the term is not computed (noted), and features computed before this version run as
    before.
  * New simulator event type `"nucleus:<i>"`, which puts several events in one nucleus.
* Tumour mode, `gain_support` (`segment_clone()`, `clone_defaults()`): an imbalanced gain state needs depth
  support against the region's offset in the other clones. A run without it becomes CN-LOH when its allelic excess
  supports that, otherwise it keeps its state with `label` `imbalance` and a new `depth_ok` column FALSE. MEDICC2
  input treats such runs as 1:1. It addresses clones that mix LOH and heterozygous cells being fitted with false
  gain states (issue #4).
* **Tumour-mode defaults changed to coarse, reproducible clones.** `disc_allele = "fraction"`, `test_allele = "centred"`,
  `drift_test = "modes_event"`, `same_post = "group"` and `gain_support = 0.3` are now the defaults of
  `clone_defaults()` and `caller.yaml`. On the lab's tumour benchmarks this gives fewer clones and doublet clones and
  much better seed agreement. Allele-only subclones smaller than the cell-level haplotype drift are merged.
  Set the five options back to `"z"`, `"z"`, `"none"`, `"clone"` and 0 for the behaviour of 0.1.1. Normal mode is
  unchanged.

# scphasedcn 0.1.1

* Extractor (Python): the `read_id` of a read whose name does not follow the SPLONGGET convention
  (`<CB>_#<uuid>_<sign><piece>of<npieces>`) is now a stable 64-bit blake2b digest of the name. It was Python's `hash()`,
  which is salted per process, so for BAMs from other pipelines the `read_id` column of `fragments.tsv.gz` changed
  from run to run. Reads named in the SPLONGGET convention are unaffected, and so are depth, alleles and the caller.
* No changes to the R package.

# scphasedcn 0.1.0

First public release of scphasedcn.

* `call_scn()`: copy-number calls per cell and region from an extraction directory (`scphasedcn run`), combining read depth
  and read-phased alleles.
  * Normal-tissue mode (`caller_normal()`): null cells selected from the data, Student-t / beta-binomial nulls,
    fold-cross-validated thresholds at 1 false call per 1,000 null cell-regions, per-cell and pooled calls, events.
  * Tumour mode (`caller_tumour()`): clone discovery, allele-specific clone copy-number states, likelihood
    assignment of every cell, per-cell inherited and confirmed calls, MEDICC2 input (`medicc2_pieces()`,
    `medicc2_table()`).
* `write_scphasedcn()` / `read_scphasedcn()`: flat result tables, the parameters used (`params.yaml`) and run notes.
* Command line: `system.file("scripts", "call_scn.R", package = "scphasedcn")`, on `--dir` + `--out`, or on a sample
  yaml (`inst/extdata/sample_template.yaml`) with `--config` + `--sample`.
* `make_cell_meta()`: the cell metadata table (cell types, known copy-number-flat cells, groups) from the extractor's
  `cells_qc.tsv` and your own labels.
* The hg38 centromere table ships with the package (`hg38_centromeres()`); the chrX centromere is masked like the
  autosomal ones.
* `simulate_scn()` / `write_sim_dataset()`: simulated SPLONGGET-like data in the extractor's layout.
* Vignettes: getting started, inputs (phasing recipe, extraction, cell metadata), outputs.
