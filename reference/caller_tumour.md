# Tumour-mode copy-number caller

Full loop:
[`split_votes()`](https://ljwharbers.github.io/scphasedcn/reference/split_votes.md)
(discovery / test halves) -\> clone feature matrix (discovery half,
re-oriented by the all-cell pseudobulk) -\>
[`discover_clones()`](https://ljwharbers.github.io/scphasedcn/reference/discover_clones.md)
-\>
[`merge_clones()`](https://ljwharbers.github.io/scphasedcn/reference/merge_clones.md)
(test half) -\> global EM round -\> within-trunk discovery
([`discover_within()`](https://ljwharbers.github.io/scphasedcn/reference/discover_within.md))
-\>
[`clone_profiles()`](https://ljwharbers.github.io/scphasedcn/reference/clone_profiles.md)
-\>
[`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
-\>
[`segment_clone()`](https://ljwharbers.github.io/scphasedcn/reference/segment_clone.md)
-\> doublet clones -\>
[`assign_cells()`](https://ljwharbers.github.io/scphasedcn/reference/assign_cells.md)
(`em_rounds` hard-EM rounds: the clone states are re-profiled on the
members of the previous round, the singlets with posterior \>=
`member_post`; clones left with fewer than `min_size` members leave the
member set; the assignment's `clone` is every cell's likelihood-best
clone, `confident` = posterior \>= `min_post` and not doublet-like) -\>
[`clean_reference()`](https://ljwharbers.github.io/scphasedcn/reference/clean_reference.md)
-\> optional depth re-pass relative to the clean pool -\> per-cell
calls. Each cell inherits its clone's non-neutral segments
(`source = "inherited"`) and is `confirmed` when its own score under the
clone state (depth + allele LLR vs neutral; allele term two-sided, or
one-sided on the clone-oriented votes with `opts$confirm_sided = "one"`)
is above the threshold of the tumour null: reference-pool cells x all
pieces + all cells x control pieces, scored as the best of the loss /
gain / CN-LOH / amp states (type-7 `null_q` quantile per stratum, strict
`>`; pooled over strata below `min_null_rows`).

## Usage

``` r
caller_tumour(
  feat,
  params = caller_defaults(),
  multiplets = NULL,
  control_regions = NULL,
  opts = clone_defaults(),
  verbose = TRUE
)
```

## Arguments

- feat:

  scn_features list

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)

- multiplets:

  optional data.table cell, prob
  ([`multiplet_table()`](https://ljwharbers.github.io/scphasedcn/reference/multiplet_table.md));
  cells above `opts$multiplet_prob` are excluded from discovery and
  pseudobulks but assigned and flagged

- control_regions:

  optional user control regions (regions table) replacing the automatic
  ones (of the confirmation null and of the reference cleaning,
  [`clean_reference()`](https://ljwharbers.github.io/scphasedcn/reference/clean_reference.md))

- opts:

  tumour-mode options
  ([`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md))

- verbose:

  log progress

## Value

`scphasedcn` (mode tumour) with clones, assignments, segments, profiles,
calls (plus `extra`: intermediate results in `null_fit$clones`)

## Details

Input `feat` (the feature step's `scn_features`, or an equivalent list):
`features` (long table), `regions`, `votes`, `atoms`, `cells`
(data.table with `cell`, in vote index order; `n_frag` / `stratum`),
`win` (1 Mb window depth: list(windows, x, e) with windows x cells
matrices, or a long table cell, chrom, start, end, x, e) and `null_fit`
(list with `tpar` =
[`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md)
table and `rho`).
