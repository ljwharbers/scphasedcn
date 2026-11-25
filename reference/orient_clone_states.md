# Common haplotype orientation of clone states per piece (MEDICC2 export, step 2)

Per piece, a reference clone is chosen among the clones with an
imbalanced state (major != minor on a diploid piece): the clone with the
largest expected orientation signal `a * |major / cn - 0.5| * sqrt(n)`,
where `a` is the clone's allelic attenuation (`clones$a`, 1 when
missing) and `n` its pooled votes on the piece. Each phase block (piece
x phase set) is oriented by the sign of the reference clone's pooled
votes (A - B; ties and blocks without reference votes are left
unoriented and do not count). Every clone's votes on the oriented blocks
are pooled into an oriented fraction `f`: the fraction of the reference
clone's major haplotype. The convention is cn_a = the reference clone's
major haplotype; a clone with an imbalanced state gets
`cn_a = major, cn_b = minor` when `f >= 0.5`, else the swap. Balanced
states (n:n, 0:0) and haploid pieces (male chrX, ploidy 1) pass through
as `cn_a = major, cn_b = minor`.

## Usage

``` r
orient_clone_states(
  pieces,
  features,
  assignments,
  clones = NULL,
  use = c("confident", "assigned"),
  near = 0.05,
  chain = c("parsimony", "phase", "none"),
  min_link = 50
)
```

## Arguments

- pieces:

  [`medicc2_pieces()`](https://ljwharbers.github.io/scphasedcn/reference/medicc2_pieces.md)
  output

- features:

  list with `votes` (cell index, atom, ps, nA, nB), `atoms` and `cells`
  (the cells the vote index refers to: data.table with a `cell` column
  or a character vector), e.g. an `scn_features` object

- assignments:

  assignments table (`cell`, `clone` and the `use` column)

- clones:

  clones table (for the attenuation `a`); NULL = 1

- use:

  cells pooled per clone: `confident` or `assigned`

- near:

  half-width around 0.5 of the low-confidence flag

- chain:

  label chaining between adjacent pieces: `parsimony` (phase-set
  fallback), `phase` (parsimony fallback) or `none` (each piece's cn_a
  is its reference clone's major haplotype)

- min_link:

  minimal summed block weight (votes) of a phase-set link

## Value

list of class `medicc2_orientation`:

- `pieces`: the input rows plus `n_votes` (clone votes on the piece),
  `n_oriented` (on oriented blocks), `f` (oriented fraction of the cn_a
  haplotype), `dev` (\|f - 0.5\|), `is_ref`, `swapped` (cn_a = minor),
  `cn_a`, `cn_b` and `flag` (`ref`, `ok`, `low_conf`, `no_votes`,
  `no_ref`, `balanced`, `haploid`).

- `qc`: one row per piece:
  `piece, chrom, start, end, orient_ref, ref_state, ref_score, flipped`
  (labels flipped by chaining), `link` (`first`, `phase`, `parsimony`,
  `none`), `link_score` (phase: weighted sign agreement between -1 and 1
  before the flip; parsimony: the weighted direction agreement),
  `n_imbalanced, n_blocks` (phase blocks with votes of any kept clone),
  `n_blocks_oriented`, `votes` (all kept clones), `votes_oriented`,
  `ref_votes`, `n_swapped`, `n_low_conf`.

- `blocks`: the reference orientation per block, reusable by a per-cell
  module: one row per piece x phase set with reference votes:
  `piece, chrom, start, end` (the piece; a molecule belongs to it when
  its atom start is in `[start, end)`), `ps` (phase set of the votes
  table), `ref` (reference clone), `nA_ref, nB_ref` and `sign` (+1:
  phase-set haplotype A is cn_a; -1: haplotype B is cn_a; 0: tie,
  unoriented). Blocks absent from the table are unoriented. A cell's
  oriented votes on a piece are `a = sum(sign > 0 ? nA : nB)`,
  `b = sum(sign > 0 ? nB : nA)` over its blocks with `sign != 0`;
  `a / (a + b)` estimates `cn_a / (cn_a + cn_b)` (attenuated).

- `params`: `use`, `near`, `chain`, `min_link`.

## Details

Chaining (`chain`): which haplotype of a piece is called "a" is
arbitrary (only the relation between clones is measured), but MEDICC2
counts an event over adjacent pieces once only when the label is
continuous. Walking along each chromosome, a piece's labels are flipped
as a whole (cn_a \<-\> cn_b in every clone, `f` -\> 1 - `f`, block signs
negated) when, relative to the previous piece, (`parsimony`) more
confidently oriented clones (`ok` / `ref`) change direction (sign of
cn_a - cn_b) than keep it (clones in the same state in both pieces count
1, others 0.5), or (`phase`) the phase sets shared by the two pieces say
that the previous piece's cn_a haplotype is this piece's cn_b (sum over
shared blocks of the product of the two reference signs, weighted by
min(\|nA_ref - nB_ref\|); used when the weights sum to \>= `min_link`
votes). The other rule is the fallback when the first has no evidence.
`parsimony` is the default: a single long phase block straddling a piece
boundary with a phase switch inside it can link the two halves of one
clonal loss with opposite labels. Pieces after a chromosome start or
next to a piece without imbalanced clones are not chained. Flipping
never changes which clones share a haplotype within a piece.

Cells pooled per clone (`use`): `confident` (default) takes the cells
confidently assigned to the clone (posterior at least `min_post`, not
doublet-like), close to the pseudobulk members that produced the clone
states, without the mixed / uncertain cells whose votes would pull `f`
towards 0.5; `assigned` takes the `assigned` column. The reference clone
orients its blocks with its own votes, so its own `f` is biased upwards
(sign reliable, magnitude not; `is_ref` marks these rows). Orientation
errors of the reference on blocks with few votes attenuate the other
clones' `f` towards 0.5 but do not flip it; a clone whose oriented `f`
is within `near` of 0.5 despite an imbalanced state is flagged
`low_conf` (its cn_a / cn_b assignment is a guess).
