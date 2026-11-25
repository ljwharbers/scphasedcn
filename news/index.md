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

- New single-cell track plots.
  [`cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/cell_tracks.md)
  builds, for each cell, the depth log2 ratio (1 Mb and aggregated
  windows, on the scale of the clone segments) and the phased B-allele
  fraction (phase blocks oriented leave-one-out by the cell’s clone
  pseudobulk in tumour mode), together with the fitted segments: the
  clone’s segments, inherited or confirmed in the cell, or the cell’s
  own events in normal mode.
  [`plot_cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/plot_cell_tracks.md)
  plots one cell and
  [`plot_cells_pdf()`](https://ljwharbers.github.io/scphasedcn/reference/plot_cells_pdf.md)
  writes many cells to a multi-page PDF. ggplot2 and patchwork are
  suggested packages, and ggrastr is used when installed.

- `inst/scripts/plot_cells.R`: command-line front end that plots every
  cell of a caller run, one PDF per clone (or per stratum in normal
  mode). It recomputes the feature step when no features .rds is given.

- Tumour mode with a fitted baseline, details and known limits:

  - Per-cell confirmation scores every null row relative to the neutral
    state of its own cell’s clone (its baseline and its depth contrast
    `c`), so a diploid reference cell is not scored as a one-copy gain
    of a near-triploid clone. A piece without a state of a clone is
    taken at that clone’s baseline.
  - Without a reference pool, the neutral of a haploid row (male chrX /
    chrY) is one copy. A haploid gain could not be confirmed by depth
    before. This changes calls on those rows in tumour runs without a
    reference pool, also at the diploid baseline.
  - Baseline candidates are narrowed only after a fit in which at least
    one clone carried allelic information.
  - `ploidy` is “fit” or a number in \[1, 8\], also from a config file.
  - A clone whose depth contrast `c` is at the lower bound of its grid
    (0.6) has a poorly determined scale and can be off by one copy; the
    run notes name such clones. Baselines 0.25 apart are near-ties when
    `c` is free.
  - Labels stay relative to the diploid normal, so per-cell call counts
    rise for clones with a baseline above 2. On a whole-genome-doubled
    tumour the false confirmations on control pieces stayed below 1 per
    1,000 on average, with one small clone at 3 per 1,000.
  - With an odd baseline, a region without allele votes (a female chrX
    is excluded from the allele channel by default) can only take a
    balanced state, one copy off the baseline.
  - In simulations at four seeds the fit put a clone at a quarter-step
    baseline next to the true one (2.25 for 2, 2.75 for 3) in 3 of 12
    tumour runs; in one of them a 2:1 chromosome of a diploid clone came
    out as 3:1. `ploidy = 2` avoids this for tumours known to be
    near-diploid.

- Tumour mode memory and run time: the per-cell confirmation and the
  reference cleaning build their vote tables per clone, the matched null
  is scored one hypothesis at a time, and the baseline fit’s segment
  refinement no longer recomputes the allele terms per depth contrast.
  Results are unchanged. On a whole-genome-doubled tumour of about
  10,000 cells the fitted baseline needs less memory than 0.1.1 did (34
  GiB against 42 GiB).

- Tumour mode without flagged null cells no longer fits its input depth
  t null and allelic rho on all cells. It fits them on the cells without
  a call on the null rows (new
  [`tumour_null_fit()`](https://ljwharbers.github.io/scphasedcn/reference/tumour_null_fit.md))
  when that selection is stable, and falls back to all cells otherwise.
  The run warns, the run notes record which null was used, and the
  warning also says when the selection keeps more than two thirds of the
  cells: a stable selection does not show that the cells are copy-number
  neutral. Results change for such runs; runs with flagged null cells
  give the same calls.

  - `call_scn(mode = "tumour", null_cells = "select")` selects the cells
    for the input null instead of being ignored, and no longer warns
    that “the features carry their own null fit” when it has just fitted
    that null.
  - Without an input rho,
    [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
    and the MEDICC2 helpers use
    [`tumour_null_fit()`](https://ljwharbers.github.io/scphasedcn/reference/tumour_null_fit.md)
    with the caller’s parameters and without multiplets, instead of a
    rho on all null rows of the flagged cells.

- Calls require a positive score, `score > max(thr, 0)`, in every
  comparison of a score with its null threshold: per-cell calls, pooled
  calls, tumour confirmation, null-cell selection, the null check,
  [`event_fdr()`](https://ljwharbers.github.io/scphasedcn/reference/event_fdr.md)
  and the reference-pool cleaning. Results change only where a null
  threshold is negative (vote-rich classes, the opt-in
  `confirm_sided = "one"`). The thresholds tables still report the null
  quantile itself.

- Pooled thresholds are computed for fold 0 only: `thresholds.tsv` no
  longer has pooled fold 1-4 rows, and the pooled step is several times
  faster. Calls are unchanged.

- Normal mode, opt-in attenuation-aware state set for chromosome-scale
  phase sets: `states_normal` is `"normal"` (default, unchanged),
  `"normal_att"` (loss and CN-LOH also scored over `f_att_normal`,
  default 0.70-0.98; `1 - eps` is always included) or `"auto"`
  (`normal_att` when the phasing regime is chromosome). Calls carry
  `f_state`, the best-fitting fraction of the called state. With the
  windowed allele term only the rows at `1 - eps` take the windowed LLR.
  It trades recall of attenuated losses and CN-LOH against more false
  loss / CN-LOH labels and some gain recall; in the author’s benchmarks
  no attenuated grid kept false arm losses within the target, so it is
  not the default.

- `model` depth normalisation, opt-in `gc_method = "ratio"`: the
  per-cell GC trend from summed counts over GC quantile groups instead
  of a lowess of per-bin log ratios, which adds a depth-dependent GC
  bias at low depth. It is applied when features are built, so existing
  feature files must be recomputed. The default `"lowess"` is unchanged.

- Normal mode fixes:

  - A named event in a depth stratum without a depth-null fit of its
    `sd_class` takes that `sd_class`’s fit over strata; before, it got
    no depth null and was never called. Events of different sizes in one
    stratum each use their own `sd_class` fit.
  - [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md):
    a class x stratum whose cells were all dropped is thresholded at
    class level, and a cell with a row that cannot be thresholded is not
    selected; before, such cells could re-enter the null set. It returns
    an empty set instead of an error when every candidate is guarded or
    called.
  - A cross-validation fold whose threshold is NA (empty training set)
    uses the fold-0 threshold.
  - [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md):
    cells without a cell type in `cell_meta` keep the features’ cell
    type.
  - [`make_pools()`](https://ljwharbers.github.io/scphasedcn/reference/make_pools.md)
    and
    [`chrom_rank()`](https://ljwharbers.github.io/scphasedcn/reference/chrom_rank.md)
    order names in C order, so a seed gives the same pools under any
    collation, also for non-ASCII names.
  - [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md),
    [`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md),
    [`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)
    and
    [`null_selection_check()`](https://ljwharbers.github.io/scphasedcn/reference/null_selection_check.md)
    complete a partial `params` list with
    [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md).
  - `win_allelic` values other than TRUE / FALSE / “auto” stop with an
    error;
    [`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)
    default `tail_k` is 30, as
    [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md);
    [`subset_scn()`](https://ljwharbers.github.io/scphasedcn/reference/subset_scn.md)
    keeps the requested order for integer selectors and stops on
    duplicate or unknown cells;
    [`make_cell_meta()`](https://ljwharbers.github.io/scphasedcn/reference/make_cell_meta.md)
    warns when `qc_args` are not applied.
  - Missing categorical values in h5ad `obs` / `var` tables are read as
    `NA` in place, no longer shifting later values.

- Tumour mode fixes:

  - Each EM round removes only its own doublet clones (ids were kept
    across the renumbering of the split pass), so `clones` has unique
    ids and every doublet flag is in `null_fit$clones$doublet_log`.
    Doublet clones flagged in the global round no longer appear as
    `clones` rows.
  - [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
    accepts every multiplet form of
    [`multiplet_table()`](https://ljwharbers.github.io/scphasedcn/reference/multiplet_table.md)
    and stops clearly when `allelic_exclude` leaves no allele votes.
  - `segment_clone(gain_support)`: a gain is above the clone’s integer
    baseline and a run is never replaced by its own state.
    [`split_disagree()`](https://ljwharbers.github.io/scphasedcn/reference/split_disagree.md)
    tests clones on the state grid of their baseline.
  - Documented that the confirmation null refits rho only, never the
    depth t.

- MEDICC2 export:

  - No segment spans a dropped range (`medicc2_pieces(max_gap = 0)` by
    default;
    [`medicc2_table()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_table.md)
    never merges across gaps).
  - Per-cell states keep the chaining flip in the leave-one-out
    orientation: cells of the reference clone were scored on the wrong
    haplotype on flipped pieces.
  - `imbalance` segments enter at the least-imbalanced state of the
    clone’s integer baseline;
    [`medicc2_table()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_table.md)
    warns when it clips at `cap`; `cell_piece_states(baseline = )`
    scores cells of clones with a baseline other than 2; numeric clone
    ids are handled.

- Functions that seed the random-number generator restore the caller’s
  random state (normal mode, tumour mode and the simulator’s genome).
  Results of
  [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
  are unchanged; a script that draws random numbers after calling one of
  these functions directly sees another stream.

- The caller finds the extractor’s bins file under its bin-size name
  (`bins_<tag>.h5ad`) through the new
  [`bins_file()`](https://ljwharbers.github.io/scphasedcn/reference/bins_file.md);
  `bins_100kb.h5ad` is still preferred when present. Bin sizes other
  than 100 kb are readable but not validated.

- Single-cell track plots: the depth baseline follows the features’
  depth reference and the clone’s ploidy baseline, normal-mode events on
  haploid chromosomes get haploid states, a clone’s reportable segments
  that the cell does not confirm are drawn as inherited, and every state
  label has a colour.

- Simulator:

  - [`simulate_scn()`](https://ljwharbers.github.io/scphasedcn/reference/simulate_scn.md):
    `seed` now controls the simulated data when the genome is the
    default or an inline
    [`sim_genome()`](https://ljwharbers.github.io/scphasedcn/reference/sim_genome.md);
    before, every seed gave the same data. Seed 1 with the default
    genome gives the same data as before.
  - Het SNPs are at unique positions; a draw on an occupied position
    moves to the next free one.
  - [`write_bins_h5ad()`](https://ljwharbers.github.io/scphasedcn/reference/write_bins_h5ad.md)
    /
    [`write_alleles_h5ad()`](https://ljwharbers.github.io/scphasedcn/reference/write_alleles_h5ad.md):
    missing values in character columns are written as anndata
    categoricals and read back as `NA`, instead of as the string `"NA"`.

- Tests: the main normal-mode and tumour-mode simulations are also
  checked at three further seeds. Before the seed fix every simulated
  fixture was drawn from one random stream, whatever seed it named.

- Extractor (Python):

  - `extract` / `run` check each chromosome against the BAM header, the
    FASTA and the chromosome sizes before the BAM pass. A mismatch
    (another assembly or naming) stops with a message naming the file,
    instead of binning on the wrong lengths or writing bins without GC.
  - QC-fail reads (flag 0x200) are skipped and counted in a new
    `n_qcfail` column; primary records without SEQ or CIGAR are skipped
    (they crashed the chromosome) and counted in `n_no_seq_cigar`.
  - FRiP (`frip`) counts peak fragments among the clean
    (non-blacklisted) fragments, matching its denominator; values drop
    slightly where blacklisted fragments overlap peaks.
  - Biallelic SNVs require REF and ALT in A/C/G/T everywhere (`extract`,
    `select-hets`, `hybrid-phase`). A SNP inside a reference skip of a
    read is no longer counted as an `other` allele.
  - The bins file is named after the bin size by default
    (`bins_100kb.h5ad` for 100 kb as before); `--tag` still overrides.
  - `merge` without `--chroms` merges the main chromosomes that have a
    finished partial; without an output directory it gives a usage
    error. `index-check` exits with status 1 on a stale or missing
    index.
  - `run` reports a worker killed for memory with the unfinished
    chromosomes and a hint to lower `--jobs`.
  - `select-hets --loci` matches loci contigs that differ by the `chr`
    prefix and exits with an error when no autosomal site is selected.
  - The per-bin layers are built sparse, with identical output and far
    less memory at small bin sizes.
  - Since 0.1.1: `--vote-ps` (one vote per phase set), input metadata in
    each partial with checks at `merge`, `--allow-missing`, an
    interruption-safe `merge`, and a report of barcode-map matches.

- Development versions are marked as such (R 0.1.1.9000, Python
  0.1.2.dev0), and the extractor records its version in its output.
  Python packaging uses an SPDX license expression and needs setuptools
  \>= 77.

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
