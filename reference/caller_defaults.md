# Default parameters of the copy-number caller

One named list shared by all parts of the caller. The calibrated values
(`shifts`, `null_q`, the depth strata) were set on SPLONGGET
normal-tissue and tumour data; see the notes per entry.

## Usage

``` r
caller_defaults()
```

## Value

named list:

- `null_q` threshold quantile of the null score (0.999, type 7: 1 false
  call per 1,000 null cell-regions); calls need score `>` threshold
  (strict, not `>=`);

- `eps` allele error floor of the beta-binomial (0.02);

- `mol_gap` fragments of one cell closer than this (bp) are one molecule
  (2000; the default of
  [`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
  `mol_gap`);

- `mol_split_ps` split molecules by phase set (FALSE): a molecule whose
  fragments lie in different phase sets (across a phase-block junction)
  gives one vote per phase set instead of one vote under the phase set
  of its first fragment
  ([`collapse_molecules()`](https://ljwharbers.github.io/scphasedcn/reference/collapse_molecules.md)
  `split_ps`);

- `votes_clean` votes from clean fragments only (FALSE): fragments with
  the blacklist bit in `flags` (the fragments the extractor leaves out
  of the depth layers) do not vote
  ([`read_votes()`](https://ljwharbers.github.io/scphasedcn/reference/read_votes.md)
  `clean`);

- `sd_floor` floor of the per-cell null SD of the window log2 ratio
  (0.05);

- `strata_breaks` / `strata_labels` depth strata on clean fragments per
  cell (`<20k`, `20-50k`, `>=50k`);

- `shifts` log2 depth shifts of loss (-0.98) and gain (+0.57)
  (calibrated on depth-only spike-ins in normal tissue);

- `k_amp` uncalibrated depth shift factor for cn \>= 4: shift =
  `k_amp * log2(cn / 2)`;

- `f_att` attenuated major-haplotype fractions of a 1:0 / 2:0 state;

- `tile_sizes` named tile sizes (bp);

- `n_folds` cross-validation folds of the null thresholds (4);

- `n_pool` cells per pool (20);

- `chrx_centromere` chrX centromere (bp, hg38) for female chrX arms;

- `acrocentric_p` arms without unique sequence, never regions;

- `min_null_rows` fewer null rows than this in a class x stratum fall
  back to the pooled threshold (5000);

- `tabix` path of the tabix binary ("" = the `tabix` on the PATH; none =
  the R fallback); `seed` RNG seed (1);

- `calib` per-cell calibration of normal mode: `"class"` (one Student-t
  null and one threshold per class x stratum, the default design) or
  `"votes"` (opt-in: thresholds per class x stratum x vote stratum
  `hstratum`, bins of the allele votes n_hap;
  [`vote_bins()`](https://ljwharbers.github.io/scphasedcn/reference/vote_bins.md));

- `vote_bins` maximal number of n_hap bins per class x stratum of the
  vote calibration (4; fewer when the class x stratum has fewer than
  `vote_bins x min_null_rows` null rows with votes);

- `calib_rho` allelic overdispersion of the normal-mode scores:
  `"global"` (one rho on the null arm rows) or `"chrom"` (opt-in: per
  chromosome, shrunk to the global;
  [`rho_by_chrom()`](https://ljwharbers.github.io/scphasedcn/reference/rho_by_chrom.md));

- `rho_chrom_k` shrinkage pseudo-count of the per-chromosome rho (100
  null arm rows; a chromosome with ~4,000 null rows keeps ~98 % of its
  own estimate);

- `mask_allelic` data-driven allelic tile mask in normal mode
  ([`allelic_noisy_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_noisy_tiles.md),
  [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)):
  the votes of 5 Mb tiles whose null overdispersion is an outlier are
  left out of the allele features (FALSE);

- `mask_rho_mad` mask threshold in MADs of log tile rho (2);

- `mask_min_votes` minimal mean votes per null row for a tile to be
  assessed (3);

- `mask_f0` shared-bias allelic tile mask in normal mode
  ([`allelic_noisy_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_noisy_tiles.md),
  [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)):
  the votes of 5 Mb tiles whose vote-weighted null-cell haplotype
  fraction f0 is off 0.5 in an outlying way are left out of the allele
  features (FALSE; independent of `mask_allelic`, both on = the union of
  the two masks;
  [`allelic_bias_tiles()`](https://ljwharbers.github.io/scphasedcn/reference/allelic_bias_tiles.md)
  shows the per-tile fractions before calling); a real clonal event
  whose carriers remain among the null cells can be masked too, so read
  the masked tiles in the run notes; not applied under `mirrored`
  orientation (noted and warned);

- `mask_f0_mad` bias threshold in MADs of `|f0 - 0.5|` over the assessed
  tiles (3);

- `mask_f0_min` minimal `|f0 - 0.5|` of a masked tile (0.03);

- `pool_sex_thr` males: the sex chromosomes get their own pooled
  threshold group (class x autosome / X / Y) with the class threshold as
  fallback (TRUE; no effect in females);

- `pool_events` normal mode: call named `event` regions (no null rows of
  their own) also in pools (FALSE). Off, their pools are scored but get
  no threshold and are never called; their per-cell calls are not
  affected. On (experimental): their pools are standardised with the
  class-median scale of their `sd_class` and thresholded with its pooled
  threshold, which is not calibrated per event (the null pooled score of
  a region depends on its own vote count, which the per-region scale of
  [`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)
  corrects for);

- `pool_tail_q` pooled per-region tail calibration
  ([`pool_scale()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scale.md)):
  quantile of the region-standardised null pools compared with the class
  (0 = off; e.g. 0.99);

- `pool_tail_k` shrinkage pseudo-count of the pooled tail factor (30
  effective pools);

- `pool_tail_xfit` cross-fit the pooled tail
  ([`pooled_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/pooled_thresholds.md)):
  the region scale and the threshold are estimated on pools of two
  disjoint halves of the null cells (TRUE);

- `win_allelic` windowed, orientation-free allele term of
  [`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md)
  ([`window_allelic_features()`](https://ljwharbers.github.io/scphasedcn/reference/window_allelic_features.md)):
  `"auto"` (computed with chromosome-scale, i.e. hybrid, phasing), TRUE
  (always) or FALSE (never);

- `win_size` window size of that term (1e6 bp): the phase is trusted
  only inside one window x phase set;

- `win_f` fractions of the windowed term in addition to `1 - eps` (the
  loss / CN-LOH rows, always included); 0.98 by default (= `1 - eps`);
  gain fractions are not specific at window level
  ([`win_f_grid()`](https://ljwharbers.github.io/scphasedcn/reference/win_f_grid.md));

- `win_bias_max`, `win_bias_min_n` windows whose haplotype fraction
  pooled over all cells is off 0.5 by more than `win_bias_max` (0.1) on
  at least `win_bias_min_n` (50) votes are left out (shared allelic
  bias);

- `win_chrx` also score chrX with the windowed term (FALSE: female
  X-inactivation gives per-cell imbalance);

- `win_score` per-cell allele term of normal mode
  ([`score_cells()`](https://ljwharbers.github.io/scphasedcn/reference/score_cells.md)
  `window`): `"none"` (oriented), `"max"` (the better of the oriented
  and the windowed LLR) or `"only"` (windowed); without window columns
  in the features the oriented term is used ("max");

- `win_min` with `win_score = "max"`, the windowed LLR replaces the
  oriented one only above this value (15; under the null model
  `P(W > 15) <= 3e-7` per row): weak window values of normal cells (W
  5-11 on balanced rows of shallow cells of a normal thyroid sample)
  otherwise displace the per-row null budget, while LOH gives W in the
  tens to hundreds;

- `null_win_guard` null-cell selection
  ([`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)):
  cells with a windowed allele LLR above the `exp(-t)` bound at the
  selection quantile on a chromosome or arm row are never null cells
  (TRUE).

The guard bound (about 9.2) is below `win_min` (15), so the windowed
term never raises a selected null cell's chromosome or arm row: the
per-row null false-call rates are those of the oriented scores, and the
windowed calls are controlled by `win_min` and the model bound
`P(W > t) <= exp(-t)`, not by the empirical null budget. The
selection-free null check
([`null_selection_check()`](https://ljwharbers.github.io/scphasedcn/reference/null_selection_check.md),
`call_scn(null_check = TRUE)`) does count them.
