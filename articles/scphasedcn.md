# Getting started with scphasedcn

scphasedcn calls copy-number alterations in single cells from SPLONGGET
data: single-cell long-read DNA sequencing in which every read carries a
10x Multiome cell barcode and is sequenced on Oxford Nanopore. Per cell
and region it combines two channels:

- **depth**: fragments per region against an expected count from a
  reference (cells of the same cell type in normal tissue, a GC /
  mappability model in tumours), scored as a z against a Student-t null;
- **read-phased alleles**: every molecule that covers a phased
  heterozygous SNP votes for one haplotype; allelic imbalance is scored
  with beta-binomial likelihood ratios.

Two modes:

- **normal tissue** (`mode = "normal"`): rare, cell-private events.
  Thresholds come from null cells (cells the caller selects as
  copy-number neutral, or that you give), cross-validated at 1 false
  call per 1,000 null cell-regions; cells are also called in pools
  (random pools within a cell type, chromatin-kNN pools). Named events
  (`events`) are called per cell; their pooled calls are off by default
  (`pool_events`: experimental, not calibrated per event).
- **tumour** (`mode = "tumour"`): clones. The caller discovers clones,
  fits their allele-specific copy-number states, re-assigns every cell
  by likelihood and reports per-cell calls inherited from the clone or
  confirmed by the cell.

## The workflow

1.  **Phase** heterozygous SNPs of the sample (e.g. Clair3 + longphase
    on the pooled BAM) and select them with `scphasedcn select-hets`.
2.  **Extract** with the Python package `scphasedcn` (`scphasedcn run`):
    per-cell 100 kb bin counts, a fragment table with the haplotype of
    every fragment, and per-cell QC.
3.  **Call** with
    [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
    in R, or the command-line script shipped with the package.

[`vignette("scphasedcn-inputs")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-inputs.md)
describes steps 1 and 2, the input files and the cell metadata;
[`vignette("scphasedcn-outputs")`](https://ljwharbers.github.io/scphasedcn/articles/scphasedcn-outputs.md)
the result tables.

## A run on simulated data

The package simulates SPLONGGET-like data and writes it in the layout
the extractor produces. Here: 120 cells of normal tissue (female) in
three cell types, a few cells carrying private events (a chr1p loss, a
chr8 gain), and one phase block per chromosome.

``` r

library(scphasedcn)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
ev = data.table(chrom = c("chr1", "chr8"), start = c(0, 0), end = c(24e6, NA),
                cnA = c(1L, 2L), cnB = c(0L, 1L), type = "private")
sim = simulate_scn(n_cells = 120, sex = "female", cell_types = c(epi = 0.4, fib = 0.35, imm = 0.25),
                   depth = c(meanlog = log(15000), sdlog = 0.3), events = ev, private_frac = 0.08,
                   multiplet_rate = 0, phase_block_mean = 1e9, switch_rate = 2e-6, seed = 7)
dir = file.path(tempdir(), "extract")
invisible(write_sim_dataset(sim, dir))
list.files(dir)
#> [1] "alleles_dna.h5ad"     "bins_100kb.h5ad"      "cells_qc.tsv"        
#> [4] "fragments.tsv.gz"     "fragments.tsv.gz.tbi" "truth_cells.tsv"
```

The simulated genome is small (chr1, chr2, chr8, chr21, chr22, chrX), so
the run gets a centromere table for it; on real hg38 data the default
([`hg38_centromeres()`](https://ljwharbers.github.io/scphasedcn/reference/hg38_centromeres.md))
applies.

``` r

cen = sim$scn$bins[, .(cstart = 1, len = max(end)), by = chrom]
cen = cen[, .(chrom, cstart, cen_left = 0.4 * len, cen_right = 0.4 * len, cend = len)]
```

The cell types are in the simulated `cells_qc.tsv`; on real data they
come from your own annotation through
[`make_cell_meta()`](https://ljwharbers.github.io/scphasedcn/reference/make_cell_meta.md).

``` r

res = call_scn(dir, mode = "normal", sex = "female", centromeres = cen, null_cells = "select", verbose = FALSE)
res
#> <scphasedcn> mode: normal
#>   cells: 120
#>   regions: 97
#>   calls: 13446 rows, 746 called, 653 reportable
#>   clones: 0
#>   notes: 14
```

Calls are per cell and region; `events` merges a cell’s overlapping
calls into events.

``` r

head(res$events[, .(cell, chrom, start, end, state, top_region, score)])
#>         cell  chrom start     end  state top_region     score
#>       <char> <char> <num>   <num> <char>     <char>     <num>
#> 1: cell_0002   chr1     0 2.4e+07   loss      chr1p 144.56874
#> 2: cell_0010   chr1     0 2.4e+07   loss      chr1p 131.50942
#> 3: cell_0018   chr1     0 2.4e+07   loss      chr1p 140.64205
#> 4: cell_0030   chr8     0 3.0e+07   gain       chr8  61.25552
#> 5: cell_0042   chr8     0 3.0e+07   gain       chr8  25.00741
#> 6: cell_0044   chr8     0 3.0e+07   gain       chr8  38.90160
```

Compare with the truth of the simulation (the cells carrying the chr1p
loss):

``` r

carriers = colMeans(sim$truth$cn[sim$scn$bins$chrom == "chr1" & sim$scn$bins$end <= 24e6, ]) < 1.5
called = unique(res$events[chrom == "chr1" & state == "loss", cell])
table(carrier = sim$scn$cells$cell %in% sim$scn$cells$cell[carriers], called = sim$scn$cells$cell %in% called)
#>        called
#> carrier FALSE TRUE
#>   FALSE   104    0
#>   TRUE      0   16
```

Write all tables, the parameters and the run notes:

``` r

out = file.path(tempdir(), "caller")
invisible(write_scphasedcn(res, out))
list.files(out)
#>  [1] "calls.tsv"        "cells.tsv"        "events.tsv"       "false_calls.tsv" 
#>  [5] "false_events.tsv" "NOTES.md"         "null_fit.tsv"     "params.yaml"     
#>  [9] "pooled_calls.tsv" "regions.tsv"      "result.rds"       "summary.tsv"     
#> [13] "thresholds.tsv"
scphasedcn_summary(res)
#>          section                                 item        value     n
#>           <char>                               <char>        <num> <num>
#>  1:        cells                                cells 120.00000000    NA
#>  2:        cells                           null_cells  97.00000000   120
#>  3:        cells                           multiplets   0.00000000   120
#>  4:        cells           cells_with_reportable_call  23.00000000   120
#>  5:        cells                     cells_with_event  23.00000000   120
#>  6:        calls                        cell|arm|gain  23.00000000     0
#>  7:        calls                        cell|arm|loss  16.00000000    16
#>  8:        calls                      cell|chrom|gain  31.00000000    11
#>  9:        calls                       cell|20Mb|gain  13.00000000     0
#> 10:        calls                       cell|20Mb|loss  16.00000000     0
#> 11:        calls                       cell|50Mb|gain  17.00000000     0
#> 12:        calls                        pool|arm|gain 121.00000000   119
#> 13:        calls                      pool|chrom|gain  34.00000000    34
#> 14:        calls                        pool|5Mb|gain  79.00000000    77
#> 15:        calls                       pool|10Mb|gain 194.00000000   194
#> 16:        calls                       pool|20Mb|gain  95.00000000    95
#> 17:        calls                       pool|50Mb|gain 107.00000000   107
#> 18:  false_calls                         arm|per_1000   2.57731959   776
#> 19:  false_calls                       chrom|per_1000   4.12371134   485
#> 20:  false_calls                         5Mb|per_1000   1.06647707  2813
#> 21:  false_calls                        10Mb|per_1000   1.37457045  1455
#> 22:  false_calls                        20Mb|per_1000   2.57731959   776
#> 23:  false_calls                        50Mb|per_1000  10.30927835    97
#> 24:  false_calls       class_x_stratum_range_per_1000   0.00000000    12
#> 25: false_events          arm_calls_per_null_cell|any   0.02061856    97
#> 26: false_events         arm_calls_per_null_cell|loss   0.00000000    97
#> 27: false_events         arm_calls_per_null_cell|gain   0.02061856    97
#> 28: false_events        arm_calls_per_null_cell|cnloh   0.00000000    97
#> 29:       events      reportable_cell_events_per_cell   0.22500000    27
#> 30:       events reportable_cell_events_per_null_cell   0.00000000    97
#> 31:       events                         pooled_calls 630.00000000 10206
#> 32:       events                     pooled_calls|knn 617.00000000  9720
#> 33:       events                    pooled_calls|type  13.00000000   486
#>          section                                 item        value     n
#>           <char>                               <char>        <num> <num>
#>                                    note
#>                                  <char>
#>  1:                                <NA>
#>  2:                                <NA>
#>  3:                                <NA>
#>  4:                                <NA>
#>  5:                                <NA>
#>  6: value = called rows, n = reportable
#>  7: value = called rows, n = reportable
#>  8: value = called rows, n = reportable
#>  9: value = called rows, n = reportable
#> 10: value = called rows, n = reportable
#> 11: value = called rows, n = reportable
#> 12: value = called rows, n = reportable
#> 13: value = called rows, n = reportable
#> 14: value = called rows, n = reportable
#> 15: value = called rows, n = reportable
#> 16: value = called rows, n = reportable
#> 17: value = called rows, n = reportable
#> 18:            null rows, CV thresholds
#> 19:            null rows, CV thresholds
#> 20:            null rows, CV thresholds
#> 21:            null rows, CV thresholds
#> 22:            null rows, CV thresholds
#> 23:            null rows, CV thresholds
#> 24:                          max 10.526
#> 25:                                <NA>
#> 26:                                <NA>
#> 27:                                <NA>
#> 28:                                <NA>
#> 29:                                <NA>
#> 30:                                <NA>
#> 31:        called pool rows / pool rows
#> 32:            kNN pools (one per cell)
#> 33:       random pools within cell type
#>                                    note
#>                                  <char>
```

## From the command line

The package ships a command-line front end; its path is

``` r

system.file("scripts", "call_scn.R", package = "scphasedcn")
#> [1] "/home/runner/work/_temp/Library/scphasedcn/scripts/call_scn.R"
```

``` bash
CALL=$(Rscript -e 'cat(system.file("scripts", "call_scn.R", package = "scphasedcn"))')
Rscript $CALL --dir results/sample/extract --out results/sample/caller --sex female --mode normal \
    --cell-meta cell_meta.tsv --null-cells select
Rscript $CALL --dir results/tumour/extract --out results/tumour/caller --sex male --mode tumour \
    --events "20q_gain=chr20:28000000-;3p_loss=chr3:1-90000000"
Rscript $CALL --help
```

## Parameters

[`caller_params()`](https://ljwharbers.github.io/scphasedcn/reference/caller_params.md)
reads the parameter file shipped with the package
(`inst/config/caller.yaml`, one key per
[`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
entry plus `normal:` / `tumour:` blocks). To change a value, copy the
file, edit it and pass it with `--caller-config` (or
`caller_params(file = ...)` in R). Every run writes the values it used
to `params.yaml`, which reads back with
[`caller_params()`](https://ljwharbers.github.io/scphasedcn/reference/caller_params.md).

``` r

p = caller_params()
p[c("null_q", "eps", "mol_gap", "strata_breaks")]
#> $null_q
#> [1] 0.999
#> 
#> $eps
#> [1] 0.02
#> 
#> $mol_gap
#> [1] 2000
#> 
#> $strata_breaks
#> [1] 20000 50000
```

The defaults were calibrated on SPLONGGET normal-tissue and tumour
samples: the one-copy depth shifts (`shifts`) on depth-only spike-ins
into normal cells, `null_q` at 1 false call per 1,000 null cell-regions.
The depth strata (`strata_breaks`, clean fragments per cell) separate
cells whose power differs; with much deeper or shallower data, move
them.

## Tumour mode

Tumour mode needs read-phased alleles (the allele votes drive the clone
discovery). Phase blocks of longphase are not oriented between blocks;
the caller orients them per clone. Optional inputs: cells known to be
normal (`flat` in the cell metadata: they set the input depth and
allelic null; the reference pool is always the normal clone the caller
finds), groups for the block orientation (`groups`, e.g. spatial
regions), known events (`events`, reported as their own regions) and
multiplet probabilities (`multiplets`).

``` r

res_t = call_scn("results/tumour/extract", mode = "tumour", sex = "male",
                 events = c(`20q_gain` = "chr20:28000000-"))
res_t$clones
head(res_t$segments)
```

MEDICC2 trees of the clones:
[`medicc2_pieces()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_pieces.md)
cuts the clone segments into common pieces,
[`orient_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/orient_clone_states.md)
assigns each clone’s allele-specific copy numbers (`cn_a`, `cn_b`) from
the molecule votes of the features, and
[`medicc2_table()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_table.md)
makes MEDICC2’s input table (write it yourself):

``` r

feat = cell_features("results/tumour/extract", method = "model", sex = "male")
res_t = call_scn("results/tumour/extract", mode = "tumour", features = feat)
pcs = medicc2_pieces(res_t$segments, res_t$clones)
ori = orient_clone_states(pcs, feat, res_t$assignments, res_t$clones)
data.table::fwrite(medicc2_table(ori), "medicc2_input.tsv", sep = "\t")
```
