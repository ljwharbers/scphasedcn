# scphasedcn (development version)

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
