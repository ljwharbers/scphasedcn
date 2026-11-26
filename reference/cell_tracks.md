# Single-cell depth and B-allele fraction tracks

Builds, for a set of cells, the data behind
[`plot_cell_tracks()`](https://ljwharbers.github.io/scphasedcn/reference/plot_cell_tracks.md):

- `depth`: per cell and window the log2 ratio
  `log2((x + 0.5) / (e + 0.5))` of the features' 1 Mb windows and of
  `agg` windows. As in
  [`clone_profiles()`](https://ljwharbers.github.io/scphasedcn/reference/clone_profiles.md),
  relative to the reference pool when it has at least `min_ref` cells
  and centred on the cell's autosomal median, so the track is on the
  scale of the clone segments. The copy-number baseline (`base`, the
  copy number at log2 0) follows the caller's expected depth. Tumour
  mode: the chromosome's ploidy with a reference pool, else 2 whatever
  the features' depth method (as
  [`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
  scores the cells: without a pool a male chrX sits at log2 -1). Normal
  mode: the chromosome's ploidy when the depth is relative to other
  cells of the sample (a reference pool, or the features' `celltype` /
  `loo` depth), else 2 (`model` depth); without a depth method in the
  features' notes, the ploidy. When the result's `clones` has a
  `baseline` column (the clone's ploidy baseline B, with absolute
  copy-number states), `base` is scaled by B / 2 for the cell's clone
  (NA or no column: 2, i.e. unchanged).

- `baf`: per cell and `agg` window the oriented haplotype-A fraction
  `nA / (nA + nB)` of the cell's molecule votes. Phase blocks are
  oriented with
  [`orient_blocks()`](https://ljwharbers.github.io/scphasedcn/reference/orient_blocks.md)
  (ties dropped) by a leave-one-out pseudobulk: of the cell's clone
  (`orient = "clone"`, the confirmation's orientation) or of all cells
  (`"pseudobulk"`); `"none"` keeps the phased haplotypes as they are
  (meaningful with chromosome-wide phasing). `"auto"`: `"clone"` in
  tumour mode, `"none"` in normal mode with chromosome-wide phasing,
  else `"pseudobulk"`.

- `fit`: the fitted segments of each cell. Tumour mode: the segments of
  the cell's clone, with `source` `confirmed` (a reportable call with
  `source` confirmed: one of the cell's events in the result's
  `events`), `inherited` (a reportable call with `source` inherited: the
  clone's segment carried to the cell without its own confirmation), or
  `clone` (neutral, or a clone segment without a reportable call in the
  cell). Normal mode: the cell's events (nominal states of their labels:
  loss 1:0, cnloh 2:0, gain 2:1, gain_as 3:0, amp 3:1, del 0:0; on a
  haploid chromosome loss 0:0, gain 2:0, amp 3:0) on a neutral baseline;
  events without coordinates take those of their `top_region`, else the
  whole chromosome. Expected log2 `c * log2(max(cn, 0.2) / base)` (`c`
  the clone's depth contrast, 1 in normal mode or when the result has
  none) and expected major fraction `0.5 + a * (f - 0.5)` (`a` the
  clone's orientation attenuation, likewise 1). Haploid chromosomes are
  the chrX / chrY rows with `allelic_ok` FALSE in the result's regions;
  a result without that column is treated as diploid throughout.

- `meta`: one row per cell (clone, posterior, assigned, n_frag, number
  of confirmed / inherited segments).

- `genome`: chromosome offsets for a concatenated genome axis.

## Usage

``` r
cell_tracks(
  features,
  result,
  cells = NULL,
  agg = 5e+06,
  orient = c("auto", "clone", "pseudobulk", "none"),
  min_ref = 30
)
```

## Arguments

- features:

  features object
  ([`cell_features()`](https://ljwharbers.github.io/scphasedcn/reference/cell_features.md))
  with `win`, `votes` and `atoms`

- result:

  caller result
  ([`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
  or
  [`read_scphasedcn()`](https://ljwharbers.github.io/scphasedcn/reference/read_scphasedcn.md))

- cells:

  cells to return (default: the result's cells that are in the features)

- agg:

  size of the aggregated windows (bp; depth points and BAF)

- orient:

  phase-block orientation of the BAF (see Details)

- min_ref:

  minimal reference pool size to reference the depth to the pool

## Value

list of data.tables `depth` (cell, chrom, start, end, size, log2), `baf`
(cell, chrom, start, end, nA, n_hap, baf, se), `fit` (cell, clone,
chrom, start, end, state, cn, label, source, ploidy, base, log2_exp,
f_exp), `meta` and `genome`; attributes `orient` and `mode`
