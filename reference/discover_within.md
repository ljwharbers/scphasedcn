# Clone discovery within a trunk

Second level of the hierarchical discovery of
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md).
The phase blocks are oriented by the pseudobulk of the trunk's cells
only (leave-one-out, discovery half of the votes): a subclonal event
carried by 5-25 % of a trunk is oriented correctly in most blocks, while
in the all-cell pseudobulk it is outweighed by balanced cells. Features:
per region (classes of `feat_disc`) the signed allelic z and the depth
z, kept only for regions whose allelic z varies across the trunk's cells
beyond the null (variance \> `disp_min`; unit variance expected) or is
imbalanced on average (\|mean\| \> `imb_min`): allelic features select
the regions, depth is used only where the allele channel shows a
clone-level event (a loss vs a gain of the same haplotype differ only in
depth there). With `opts$disc_allele = "fraction"` the allelic features
are the coverage-free oriented fractions of
[`clone_feature_matrix()`](https://ljwharbers.github.io/scphasedcn/reference/clone_feature_matrix.md)
(regions still selected on z); `opts$test_allele` and `opts$drift_test`
set the merge tests
([`merge_clones()`](https://ljwharbers.github.io/scphasedcn/reference/merge_clones.md)).
Discovery
([`discover_clones()`](https://ljwharbers.github.io/scphasedcn/reference/discover_clones.md))
runs on the trunk's deeper cells (`n_frag` \>= the `deep_q` quantile, at
least `4 * min_size` cells), whose per-cell profiles separate subclones;
[`merge_clones()`](https://ljwharbers.github.io/scphasedcn/reference/merge_clones.md)
then tests every pair of leaves (`n_nb = Inf`) on the test halves
(`cdata`). Shallow cells are left to the likelihood assignment.

## Usage

``` r
discover_within(
  inp,
  tcells,
  votes_disc,
  feat_disc,
  cdata,
  n_frag,
  opts = clone_defaults(),
  seed = 1
)
```

## Arguments

- inp:

  a `clone_input`

- tcells:

  trunk cells

- votes_disc:

  discovery half of the votes
  ([`split_votes()`](https://ljwharbers.github.io/scphasedcn/reference/split_votes.md))

- feat_disc:

  per-cell depth features of the discovery half (cell, region, z)

- cdata:

  [`contrast_data()`](https://ljwharbers.github.io/scphasedcn/reference/contrast_data.md)
  of the test halves (all cells)

- n_frag:

  named clean fragments per cell

- opts:

  [`clone_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/clone_defaults.md)

- seed:

  RNG seed

## Value

list(clones = named integer labels of the discovery cells, tests,
regions = selected regions, n_leaves)
