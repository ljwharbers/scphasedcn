# scphasedcn (development version)

* New single-cell track plots. `cell_tracks()` builds, for each cell, the depth log2 ratio (1 Mb and aggregated
  windows, on the scale of the clone segments) and the phased B-allele fraction (phase blocks oriented leave-one-out
  by the cell's clone pseudobulk in tumour mode), together with the fitted segments: the clone's segments, inherited
  or confirmed in the cell, or the cell's own events in normal mode. `plot_cell_tracks()` plots one cell and
  `plot_cells_pdf()` writes many cells to a multi-page PDF. ggplot2 and patchwork are suggested packages, and ggrastr
  is used when installed.
* `inst/scripts/plot_cells.R`: command-line front end that plots every cell of a caller run, one PDF per clone (or per
  stratum in normal mode). It recomputes the feature step when no features .rds is given.

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
