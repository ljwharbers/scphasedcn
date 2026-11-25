# Split clones whose cells disagree with the clone state on a segment

Targeted split pass of
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
(`split_pass = "disagree"`). Allele-driven discovery can put cells that
share a haplotype imbalance but differ by one copy in one clone (e.g. a
chr21 loss subclone inside a chr21 CN-LOH clone; a chr21 3:0 / 4:0
LOH-gain subclone inside a 2:1 clone, fitted 3:1). For each clone with
\>= 2 `min_size` members and each of its non-neutral segments (ploidy 2,
\>= `split_min_len`; the copy-neutral 2:0 included), every member is
scored on the discovery halves (`win_disc`, `votes_disc`) under the
clone state and under each alternative of
[`split_alt_states()`](https://ljwharbers.github.io/scphasedcn/reference/split_alt_states.md).
Depth: the cell's segment log2 ratio, centred on the members' median
(the clone state; clone profiles are relative to the clone baseline), in
units of its `sd_null` scaled to the half depth and inflated by the
members' own robust spread, t density at `c * log2(cn / cn_state)`
(clone c from `par`). Alleles: the cell's votes oriented by the clone's
leave-one-out pseudobulk, beta-binomial at `0.5 + a_eff (f - 0.5)`,
a_eff = the attenuation of the members' pooled oriented fraction. The
contradiction of a segment is the fraction of members whose best state
is not the clone state; its noise floor is the same fraction on data
simulated under the clone state at each cell's SD and vote count (5
replicates). The alternative is the state with the largest excess
preference over its noise floor. A segment is a candidate when the
contradiction is \>= 10 % and exceeds its noise floor by \>=
`split_min_excess`, and a two-means split of the per-cell log-likelihood
difference (alternative - clone state) has both parts \>=
max(`min_size`, `small_clone`, 10 %) and a separation \>= `split_sep`.
Candidates are ranked by the excess contradiction (not by separation:
the 8q amplicon's depth spread is unimodal and wide) over all clones; up
to `split_cand` per clone are tested, at most one split per clone and
`split_max` splits in all. A candidate is applied when (1) the
permutation split test (`split_test()`, budget `q`) on the test halves,
restricted to the `cdata` regions the segment overlaps (\>= 80 %, else
\>= 50 % of the shorter of the two; depth + allele), says different, and
(2) the two parts re-fitted on the test halves
([`clone_profiles()`](https://ljwharbers.github.io/scphasedcn/reference/clone_profiles.md),
`fit_states`,
[`segment_clone()`](https://ljwharbers.github.io/scphasedcn/reference/segment_clone.md))
differ on the segment (copy number or major fraction,
[`clone_state_distance()`](https://ljwharbers.github.io/scphasedcn/reference/clone_state_distance.md)
with `min_df = 0.01, min_dcn = 1` as the caller's same-state merge) over
\>= max(`same_dist`, half the segment), with the re-fitted states at \>=
1 copy, both within 2 and one within 1 copies (total) of the clone state
(the clone state of a mixture is an average: e.g. a fitted 4:1 / 3:1
clone that is a mix of 2:1 and 3:0 / 4:0 cells).

## Usage

``` r
split_disagree(
  members,
  segments,
  par,
  inp,
  win_disc,
  votes_disc,
  cdata,
  win_test = NULL,
  votes_test = NULL,
  ref_cells = character(),
  fit_states = NULL,
  opts = clone_defaults(),
  params = caller_defaults(),
  seed = 1,
  frac = 0.5
)
```

## Arguments

- members:

  data.table cell, clone (confident members)

- segments:

  segments table of the members' clones

- par:

  data.table clone, c, a
  ([`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  `par`)

- inp:

  scn_features list or `clone_input` (regions, atoms, cells, features
  for `sd_null`, tpar, rho)

- win_disc, votes_disc:

  discovery halves of the window depth and the votes
  ([`split_win()`](https://ljwharbers.github.io/scphasedcn/reference/split_win.md),
  [`split_votes()`](https://ljwharbers.github.io/scphasedcn/reference/split_votes.md))

- cdata:

  [`contrast_data()`](https://ljwharbers.github.io/scphasedcn/reference/contrast_data.md)
  of the test halves

- win_test, votes_test:

  test halves for the re-fit (default: the full `inp$win` / `inp$votes`)

- ref_cells:

  depth reference cells of the re-fit (used with \>= `min_ref` cells)

- fit_states:

  optional replacement of
  [`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
  (`function(profile, grid)`)

- opts:

  [`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md)
  (`min_size`, `split_sep`, `split_max`, `split_cand`, `split_min_len`,
  `split_min_df`, `split_min_excess`, `q`, `n_perm`, `max_cells`,
  `depth_delta`, `same_dist`, `small_clone`, `size_small`, `min_ref`,
  `p_switch`, `min_win`, `min_llr`)

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`sd_floor`, `eps`)

- seed:

  RNG seed

- frac:

  fraction of the counts in the discovery half

## Value

list(members = relabelled members (split clones get suffixes `a` / `b`,
`b` = the part that prefers the alternative), tests = one row per
candidate: clone, region, chrom, start, end, state, alt, n, contra
(fraction of members whose best state differs from the clone's),
contra_null (its noise floor), excess, contra_alt / contra_alt_null (the
same for the alternative), a_eff, infl (depth SD inflation), sep, n_a,
n_b, d_a, d_b (mean log-likelihood difference per part), outcome, rank,
stat, thr, n_perm, p (permutation p of the test-half split test),
differ, state_a, state_b (re-fitted states at the segment midpoint),
dist_bp, split), cells_b = named list of the moved cells per applied
split (names = the `b` clones), sets = per candidate (rows of `tests`)
the proposed parts list(a, b)
