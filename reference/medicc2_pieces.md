# Common copy-number pieces of the clones of a tumour run (MEDICC2 export, step 1)

Drops doublet clones, clones with fewer than `min_assigned` assigned
cells and clones without segments; builds the union of the remaining
clones' segment breakpoints
([`segment_pieces()`](https://ljwharbers.github.io/scphasedcn/reference/segment_pieces.md))
on `chroms`, maps every clone's state to every piece (segment containing
the piece midpoint, as `piece_states()`) and keeps the pieces covered by
every kept clone (gaps of any clone, e.g. missing windows or a
clone-specific chromosome start, are dropped). Pieces shorter than
`min_len` are then absorbed into an adjacent piece of the same
chromosome (gap \<= `max_gap`): the neighbour whose states differ in the
fewest clones, the longer one on ties; the neighbour keeps its states
and extends over the short piece. A short piece without an eligible
neighbour is dropped. Adjacent pieces with identical states are not
merged here: merging needs the common haplotype orientation and is done
by
[`medicc2_table()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_table.md).
Segments labelled `imbalance` (allelic imbalance without depth support,
`segment_clone(gain_support = )`) enter as 1:1: MEDICC2 sees no extra
copy where the clone's depth does not show one.

## Usage

``` r
medicc2_pieces(
  segments,
  clones = NULL,
  min_len = 1e+06,
  max_gap = Inf,
  min_assigned = 1,
  chroms = paste0("chr", c(1:22, "X")),
  normal = "N0"
)
```

## Arguments

- segments:

  segments table (`clone, chrom, start, end, state, cn, label, ploidy`;
  `nA, nB, f` optional)

- clones:

  clones table (`clone`, `doublet`, `n_assigned`); NULL keeps every
  clone of `segments`

- min_len:

  pieces shorter than this (bp) are absorbed into a neighbour

- max_gap:

  maximal gap (bp) between a short piece and the neighbour it is
  absorbed into

- min_assigned:

  clones with fewer assigned cells are dropped (the normal clone
  `normal` is always kept)

- chroms:

  chromosomes exported (MEDICC2: chr1-22 and chrX)

- normal:

  name of the normal clone

## Value

long data.table, one row per clone x piece: `piece`
(`seg|all|chrom:start-end`), `chrom`, `start`, `end`, `clone`, `state`,
`major`, `minor`, `cn`, `label`, `ploidy`, `imbalanced` (major != minor
on a diploid piece; these need orientation), `seg_nA`, `seg_nB`, `seg_f`
(the containing segment's pooled votes in the clone's own orientation,
NA when absent). Attributes `clones_dropped` (clone, reason) and
`pieces_dropped` (chrom, start, end, reason: `not_all_clones` /
`short_isolated`) and `n_pieces_raw` (pieces before the coverage
filter).
