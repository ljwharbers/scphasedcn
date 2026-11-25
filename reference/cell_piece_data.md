# Per-cell evidence on the MEDICC2 pieces: depth and oriented molecule votes

Depth (the per-cell analogue of
[`clone_profiles()`](https://ljwharbers.github.io/scphasedcn/reference/clone_profiles.md)):
the cell's 1 Mb window counts summed over the piece (windows fully
inside it;
[`region_features()`](https://ljwharbers.github.io/scphasedcn/reference/region_features.md)),
`log2ratio = log2((x + 0.5) / (e + 0.5))` relative to the reference
pool's pseudobulk of the piece (`ref_center = "pseudobulk"`,
`log2((sum x + 0.5) / (sum e + 0.5))` over the reference cells, as the
clone profiles; `"mean"` = the mean of the reference cells' log2 ratios,
as
[`assign_cells()`](https://ljwharbers.github.io/scphasedcn/reference/assign_cells.md),
which a few reference cells with a large negative offset pull down),
then (`center`) minus the cell's own baseline: the window-weighted
median of its log2 ratios over the autosomal pieces (the clone profiles
are centred on the clone's autosomal median in the same way, so states
are relative to the baseline). SD: the cell's null SD of the piece's
size class (`sd_null` of the piece's `sd_class`); with `sd_scale`, a
piece shorter than its tile class (e.g. a 2 Mb piece scored with the 5
Mb SD) gets `sd_null * sqrt(class size / covered window length)` (the
independent-window limit); pieces on haploid chromosomes get
`haploid_sd` x that; `s_sys` (systematic depth bias of a region between
the scored cells and the reference pool, the clone model's `s_sys`) is
added in quadrature: `sd = sqrt(sd_scaled^2 + s_sys^2)`.
`z = log2ratio / pmax(sd, sd_floor)`; pieces without a window have NA
depth (no evidence). Alleles: the cell's molecule votes on the piece
(atom start in `[start, end)`) per phase set, oriented by the block
table: `oA` = votes of the cn_a haplotype (`nA` when `sign > 0`, `nB`
when `sign < 0`), `oN` = oriented votes. With `loo`, the reference
clone's pooled votes are first reduced by the cell's own votes when the
cell was one of the cells that oriented the block (`loo` maps cells to
the clone whose pool they joined; only blocks whose `ref` is that clone
are affected): the orientation is then the one the other cells give, and
the cell's own votes cannot confirm themselves. This is the only use of
a cell's clone membership, and it only removes the cell's own data from
the orientation input.

## Usage

``` r
cell_piece_data(
  features,
  pieces,
  blocks,
  cells = NULL,
  ref_cells = NULL,
  loo = NULL,
  params = caller_defaults(),
  sd_scale = TRUE,
  ref_center = c("pseudobulk", "mean"),
  center = TRUE,
  s_sys = 0,
  haploid_sd = sqrt(2)
)
```

## Arguments

- features:

  scn_features list (`win`, `long`, `regions`, `votes`, `atoms`,
  `cells`; as the tumour mode reads it)

- pieces:

  piece table (`piece, chrom, start, end`, one row per piece or per
  clone x piece)

- blocks:

  block orientation (`piece, ps, sign` and, for `loo`,
  `ref, nA_ref, nB_ref`)

- cells:

  cells to score (default all cells of `features`)

- ref_cells:

  reference cells of the depth centring (\>= 10 used; NULL = no
  centring)

- loo:

  NULL, or data.table `cell, clone`: the cells whose votes formed each
  reference clone's pool

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`sd_floor`)

- sd_scale:

  scale the null SD of pieces shorter than their tile class

- ref_center:

  reference centring: `pseudobulk` or `mean`

- center:

  subtract the cell's autosomal median (window-weighted over pieces)

- s_sys:

  systematic SD (log2) added in quadrature to the per-cell SD

- haploid_sd:

  factor on the null SD of pieces on haploid chromosomes (male chrX /
  chrY): `sd_null` comes from diploid regions, and one copy gives half
  the fragments, so z on haploid pieces is overdispersed

## Value

data.table cell, piece, x, e, n_win (1 Mb windows), log2ratio (centred),
offset (the cell's subtracted baseline), sd_null, sd, z, stratum,
sd_class, n_frag, oA, oN, n_blocks (oriented phase sets with votes)
