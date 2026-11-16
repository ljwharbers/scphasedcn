# Changelog

## scphasedcn (development version)

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
