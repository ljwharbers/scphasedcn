# Outputs: the result object and tables

[`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
returns a `scphasedcn` object; `write_scphasedcn(res, out)` writes it as
flat tables (the command line does both). `read_scphasedcn(out)` reads a
written run back.

## Files

| file | content |
|----|----|
| `cells.tsv` | one row per called cell: cell type, depth stratum, null cell flag, `n_calls`, `n_events`, and the `multiplet` flag (normal mode always; tumour mode when `multiplets` are given: the cells above its `multiplet_prob`); tumour mode adds `clone`, `posterior`, `doublet_like`, `ref_pool` |
| `calls.tsv` | per cell x region calls (normal mode: the called rows; tumour mode: every clone row and every per-cell inherited / confirmed row) |
| `pooled_calls.tsv` | normal mode: called rows of the cell pools (named `event` regions only with `pool_events`, off by default: experimental, not calibrated per event) |
| `events.tsv` | a cell’s overlapping called regions merged into events |
| `thresholds.tsv` | the null thresholds per score, region class and depth stratum (and cross-validation fold); named `event` regions have no null rows of their own and use the thresholds of their `sd_class` (column `borrowed_from`, present when a run has events; the pooled ones only with `pool_events`) |
| `regions.tsv` | the regions scored (arms, chromosomes, 5-50 Mb tiles, events; tumour mode also clone pieces) |
| `null_fit.tsv` | the depth null (Student-t per class x stratum) and the allelic overdispersion `rho` |
| `false_calls.tsv`, `false_events.tsv` | normal mode: the realised false-call rate on the null cells |
| `clones.tsv`, `assignments.tsv`, `clone_segments.tsv`, `clone_profiles.tsv.gz`, `reference_pool.txt` | tumour mode: clones, per-cell assignment, clone copy-number segments, clone window profiles, the normal reference cells |
| `summary.tsv` | headline numbers ([`scphasedcn_summary()`](https://ljwharbers.github.io/scphasedcn/reference/scphasedcn_summary.md)) |
| `params.yaml` | every parameter of the run; reads back with `caller_params("params.yaml")` |
| `NOTES.md` | the run notes (inputs, steps, run time, package version) |
| `result.rds` | the object without the large `scores` / `features` tables |
| `scores.tsv.gz`, `features.tsv.gz` | optional (`write_scores`, `write_features`; large) |

## Calls

| column | meaning |
|----|----|
| `cell` | the cell (normal mode pools: the pool id; tumour mode clone rows: the clone) |
| `region`, `class`, `chrom`, `start`, `end` | the region (0-based, half-open) and its class: `arm`, `chrom`, `5Mb` … `50Mb`, `event` |
| `state` | `loss` (1:0), `gain` (2:1), `cnloh` (2:0, copy-neutral LOH), `gain_as`, `amp`, `del`, `neutral`, or a clone state label |
| `score`, `thr`, `called` | the region score, its null threshold and `score > thr` |
| `level` | `cell`, `pool` or `clone` |
| `source` | `cell` (the cell’s own evidence), `inherited` (from its clone or pool) or `confirmed` (a clone call the cell’s own evidence confirms) |
| `reportable` | the call survives the nesting rules (e.g. a chromosome gain is not reported inside a stronger nested arm loss) |
| `resolved` | normal mode, per cell: why [`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md) changed the row (NA = unchanged): `no_depth` (a chromosome / event gain without depth support), `nested <region> <state>` (a gain explained by a stronger nested loss / CN-LOH), `replaces <region> gain` (the nested call made reportable in its place) |
| `gate_z`, `gate_thr` | normal mode with `gain_gate: "z"`: on the gated chromosome / event gains, the depth z and its threshold (the larger of the midpoint between neutral and the gain shift and the `gain_gate_q` quantile of the depth null); a gain has depth support when `gate_z > gate_thr`. With the default `gain_gate: "dl"` the gate is the sign of the gain depth log-likelihood ratio (z above the midpoint) and these columns are absent |
| `z`, `f`, `n_hap` | normal mode, cell rows: depth z, major-haplotype fraction and molecule votes of the cell in the region (not in tumour mode) |
| `stratum` | normal mode, cell rows: the cell’s depth stratum (clean fragments; `<20k`, `20-50k`, `>=50k`; not in tumour mode) |
| `confident` | tumour mode: the cell’s clone assignment is confident |

## Clone segments (tumour mode)

`clone_segments.tsv`: per clone and segment `chrom`, `start`, `end`,
`n_win` (1 Mb windows), the pooled oriented votes `nA`, `nB`, the total
copy number `cn`, the allele-specific state (`1:1`, `1:0`, `2:1`, `2:0`,
`3:1`, …), `log2` (depth relative to the clone baseline), `f`
(major-haplotype fraction), `llr`, `confident`, `label` (neutral, loss,
cnloh, gain, gain_as, amp, del) and `ploidy` (1 on a male’s chrX /
chrY). The clone `N0` is the fixed normal clone.

`assignments.tsv`: per cell the likelihood-best clone, its `posterior`,
the margin to the second clone, whether the cell looks like a doublet of
two clones (`doublet_like`, `pair`, `p_mix`), and `confident` (posterior
\>= 0.9, not doublet-like).

## The object

``` r

library(scphasedcn)
library(data.table)
sim = simulate_scn(n_cells = 60, sex = "female", cell_types = c(a = 0.5, b = 0.5), events = NULL,
                   multiplet_rate = 0, depth = c(meanlog = log(12000), sdlog = 0.3), phase_block_mean = 1e9,
                   switch_rate = 2e-6, seed = 3)
d = file.path(tempdir(), "out_extract")
invisible(write_sim_dataset(sim, d))
cen = sim$scn$bins[, .(chrom, cstart = 1, len = max(end)), by = chrom][, .(chrom, cstart, cen_left = 0.4 * len,
                                                                             cen_right = 0.4 * len, cend = len)]
res = call_scn(d, "normal", sex = "female", centromeres = cen, null_cells = "select", verbose = FALSE)
names(res)
#>  [1] "mode"         "cells"        "regions"      "features"     "scores"      
#>  [6] "thresholds"   "calls"        "events"       "pools"        "clones"      
#> [11] "assignments"  "segments"     "profiles"     "null_fit"     "params"      
#> [16] "notes"        "false_calls"  "false_events" "embedding"    "feature_info"
#> [21] "run"
res$thresholds[level == "cell" & fold == 0][1:5, .(score, class, stratum, q, thr, n_null)]
#>        score  class stratum     q        thr n_null
#>       <char> <char>  <char> <num>      <num>  <int>
#> 1:   s_depth   10Mb    <20k 0.999  4.0226146    780
#> 2: s_allelic   10Mb    <20k 0.999  3.1379661    780
#> 3:    s_comb   10Mb    <20k 0.999 -0.8174322    780
#> 4:   s_depth   20Mb    <20k 0.999  3.6923150    480
#> 5: s_allelic   20Mb    <20k 0.999  2.8782689    480
out = file.path(tempdir(), "out_caller")
invisible(write_scphasedcn(res, out))
res2 = read_scphasedcn(out)
identical(res2$calls, res$calls)
#> [1] FALSE
```

The components of the object:

| component | normal mode | tumour mode |
|----|----|----|
| `cells` | cell, cell_type, stratum, fold, null_cell, multiplet | the cells with clone, posterior, … |
| `regions` | the scored regions | the scored regions plus clone pieces |
| `scores` | per cell x region scores (with the null flag) | NULL |
| `thresholds` | per-cell (cross-validated) and pool thresholds (`level`) | per depth stratum and clone hypothesis |
| `calls`, `events` | see above | see above |
| `pools` | pool membership (`kind` = `type` or `knn`) | NULL |
| `clones`, `assignments`, `segments`, `profiles` | NULL | see above |
| `null_fit` | depth t parameters, rho, the null cells and how they were selected | depth t parameters, rho, reference pool, clone discovery details |
| `params` | all parameters, including the `normal:` caller arguments | all parameters, including the `tumour:` options |
| `features`, `feature_info`, `embedding`, `notes` | the per-cell region features and the run notes | the same |
