# scphasedcn (development version)

* **Tumour mode fits the clone baseline (ploidy) by default** (`clone_defaults()` `ploidy = "fit"`). Up to now the
  median depth of every clone was taken to be 2 copies, so a near-triploid or whole-genome-doubled clone got
  states on a compressed scale: 2:2 as 1:1, 2:0 as 1:0, 3:0 as 2:0, and depth-neutral 2:1 often as 1:1.
  `fit_clone_states(baseline = "fit")` now fits per clone the copy number of its median depth (the baseline:
  1.5 to 5) from depth and alleles together, so states are absolute copy numbers.
  * The fit compares baselines on one likelihood scale. Above a baseline of 2 the state grid drops the near-balanced
    states such as 4:3, which would only absorb orientation bias.
  * A clone whose allele attenuation is at the grid floor carries no allelic information on the scale. It takes 2 or
    the sample baseline (the cell-weighted median of the informative clones).
  * Half the best baseline is preferred within `ploidy_tol`, because a whole-genome doubling of the same profile
    fits equally well.
  * Labels stay relative to the diploid normal: 2:1 is a gain at any baseline, 1:1 neutral.
  * `clones` gains `baseline` and `ploidy` (the mean autosomal copy number), and the run notes list the baselines.
  * A clone fitted at baseline 2 has exactly the states of 0.1.x. `ploidy = 2` restores the old behaviour for every
    clone; a number fixes the baseline.
  * On COLO829 (bulk ploidy 3.1) the clones fit baseline 3, and on CCS15M (bulk 3.46) 3 to 3.5. Near-diploid tumours
    stay at 2.

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
