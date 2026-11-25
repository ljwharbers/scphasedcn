# Defaults of the tumour mode

Discovery: `classes`, `depth_delta`, `cap`, `min_n`, `disc_allele`,
`test_allele`, `n_pc`, `k`, `max_leaves`, `min_size`, `q`, `n_perm`,
`drift_test`, `max_cells`, `split_depth` (discovery / test halves of the
window depth as well as of the votes). Merge / split tests
([`merge_clones()`](https://ljwharbers.github.io/scphasedcn/reference/merge_clones.md),
[`discover_within()`](https://ljwharbers.github.io/scphasedcn/reference/discover_within.md),
the split pass; `perm_decision`): `"quantile"` (default) = the observed
contrast above the type-7 `q`-quantile of `n_perm` permutations (with
200 permutations at `q = 0.999` this is "beats about the maximum": ~6
false splits per 1,000 tests of one clone, not 1); `"pvalue"` =
permutation p-value
`(1 + #null >= observed) / (1 + permutations) <= 1 - q`, over
`max(n_perm, ceiling(1 / (1 - q)) - 1)` permutations (999 at
`q = 0.999`) with an exact early stop for "same" (as soon as p can no
longer reach `1 - q`); a "different" verdict runs every permutation. The
run notes give the number of tests and the expected false splits
`n_tests x (1 - q)`.

## Usage

``` r
clone_defaults()
```

## Value

named list (see
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md))

## Details

Allele scale of discovery and tests (coverage-free by default since the
development version): `disc_allele` = the allelic columns of the
discovery feature matrices
([`clone_feature_matrix()`](https://ljwharbers.github.io/scphasedcn/reference/clone_feature_matrix.md),
[`discover_within()`](https://ljwharbers.github.io/scphasedcn/reference/discover_within.md));
`test_allele` = the per-cell allele values of the merge / split
contrast. `"z"` (the default up to 0.1.1) = the standardised count
`(nA - n / 2) / sd(n)`, whose expectation in an imbalanced region grows
like `(f - 0.5) sqrt(n)`: deep and shallow cells of one clone form
separate leaves and test different. `"fraction"` (default of
`disc_allele`) = the cell's oriented fraction `nA / n`, whose
expectation does not depend on the cell's depth (discovery: in units of
the null SD at the region's median vote count, the same scale as `"z"`
for a typical cell; a cell below `min_n` votes gets the region's mean
instead of 0; the regions of
[`discover_within()`](https://ljwharbers.github.io/scphasedcn/reference/discover_within.md)
are still selected on `"z"`); `test_allele = "centred"` (default) = the
standardised count about the compared cells' pooled fraction instead of
1 / 2 (mean 0 at any depth when the two sets share the fraction; weights
the cells like `"z"`). The imbalance term of the contrast (depth
support) always uses `"z"`.

Drift-aware split test (`drift_test`, default `"modes_event"`): the
discovery and test halves of the votes are thinned from the same
molecules, so a cell-level drift of the haplotype share beyond binomial
noise is in both: leaves that discovery cut along the drift of one clone
differ on the test half too, and the permutation null of random splits
does not model that selection. `"none"` (the default up to 0.1.1) = the
permutation contrast decides. `"modes"` = a pair the contrast calls
different stays different when a depth term alone beats the null
threshold (the drift is allelic); an allele-driven difference must in
addition be discrete on at least one of the (up to 3) regions whose
allele term alone beats the threshold: the per-cell oriented fractions
there of the population under discovery (all leaf cells of the merge:
the discovery cells, or a trunk's deep cells; the clone itself in the
split pass; test half, phase blocks oriented leave-one-out by the
population, cells with \>= 10 votes, at most 1,500 cells) are multimodal
by Silverman's critical bandwidth test (smoothed bootstrap, 199 draws,
level 0.01; a continuous drift, also a skewed one, has one mode, a
subclone adds one), and the two clones lie on different sides of the
antimode (their fractions above it differ by \>= 0.5). The population,
not the union of the two clones, is tested because two leaves cut from
one continuum can make a union with two modes. A subclone whose allelic
difference is smaller than the drift has no mode of its own and is
merged. The tests tables gain `differ_contrast` (the contrast's
verdict), `by_depth`, `discrete`, `mode_p`, `mode_sep`. `"modes_event"`
(default) = `"modes"`, and a pair also stays different when the events
of one set are absent from the other: the other set is allelically
balanced (oriented by its own other half, \|t\| \<= 2) on at least 75 %
of the chromosomes, and on at least 2 of them, where the first set is
imbalanced beyond the null threshold (two lineages), or the other set
has no imbalanced region at all and is balanced on one of them (normal
cells against cells with an event; column `by_event`). A cut along a
drift of one clone can bring a set to 1/2 on one arm, not on most of the
other set's events on independent chromosomes, nor make it event-free.
With `"modes"`, a group of normal cells with some tumour cells can merge
with the normal cells (allele-only difference, no mode of its own), and
the merged clone, which carries weak events, is then not replaced by the
normal clone `N0` (no reference pool).

Hierarchy: trunks = average-linkage groups of
[`clone_allelic_distance()`](https://ljwharbers.github.io/scphasedcn/reference/clone_allelic_distance.md)
cut at `trunk_dist` (bp x haplotype fraction; e.g. 100 Mb differing by
0.3); trunk cells = trunk posterior \>= `trunk_post`;
[`discover_within()`](https://ljwharbers.github.io/scphasedcn/reference/discover_within.md)
with `disp_min`, `imb_min`, `deep_q`, `depth_w`, `max_leaves_trunk`,
`n_nb_trunk`; final split pass after the EM rounds (`split_pass`):
`FALSE` (default) or `"disagree"` =
[`split_disagree()`](https://ljwharbers.github.io/scphasedcn/reference/split_disagree.md)
(opt-in: members whose own segment state disagrees with the clone state;
`split_sep`, `split_max` splits per pass, `split_cand` candidates tested
per clone, `split_min_len` minimal segment length, `split_min_df`
major-fraction step of an alternative between amplified states, see
[`split_alt_states()`](https://ljwharbers.github.io/scphasedcn/reference/split_alt_states.md),
`split_min_excess` minimal contradiction above its noise floor). Normal
hypothesis (`normal_hypothesis`): a fixed all-neutral clone `N0`
replaces every fitted clone without a non-neutral segment \>=
`normal_min_event` bp in each round. States: `size_small`,
`small_clone`, `p_switch`, `min_win`, `min_llr`; segment rules of
[`segment_clone()`](https://ljwharbers.github.io/scphasedcn/reference/segment_clone.md)
(all segmentations of the run, including the split pass): `min_excess`
(pooled allelic excess of imbalanced runs), `max_contra` (gains / losses
contradicted by depth become copy-neutral), `bal_frac` and `min_bal_len`
(support of balanced 2:2 runs), `gain_support` (an imbalanced gain needs
this share of its expected depth shift against the region's offset, else
CN-LOH or label "imbalance" with `depth_ok` FALSE; default 0.3, 0 =
off). Clone baseline (`ploidy`): `"fit"` (default) fits per clone the
copy number of its median depth (the baseline B: 3 for a near-triploid
clone) from depth and alleles together
([`fit_clone_states()`](https://ljwharbers.github.io/scphasedcn/reference/fit_clone_states.md)
`baseline = "fit"`), so the states are absolute copy numbers; a number
in \[1, 8\] fixes B for every clone; `2` is the diploid baseline of
0.1.x (a near-triploid or whole-genome-doubled clone is then fitted on a
compressed scale: 2:2 as 1:1, 2:0 as 1:0). `ploidy_tol` (log-likelihood
per window, finite, \>= 0) is the tolerance within which a baseline of
half the best one (a whole-genome doubling of the same profile) is
preferred. Labels stay relative to the diploid normal (2:1 is a gain at
any baseline, 1:1 neutral), so a clone with a baseline above 2 has
non-neutral segments, and per-cell calls, over most of its genome.
`clones` reports `baseline` and `ploidy` (the mean autosomal copy
number); the run notes list the fitted baselines and name every clone
whose depth contrast `c` sits at the lower bound of its grid (its states
can be off by one copy). Assignment: `mix_prior`, `mix_w`, `mix_margin`,
`mix_rel`, `mix_min_dist` (mixture hypotheses only for clone pairs whose
states differ allelically on \>= this many bp), `max_mix_prior`,
`min_post` (confident), `member_post` (pseudobulk members), `em_rounds`,
`same_dist` (clones whose states differ on \< this many bp within a
trunk are merged after each round). `same_post` = the membership
posterior of
[`assign_cells()`](https://ljwharbers.github.io/scphasedcn/reference/assign_cells.md)
that `min_post` and `member_post` cut: `"clone"` (the default up to
0.1.1) = the best clone's posterior; `"group"` (default) = summed over
the clones of the best clone's trunk whose states differ on \<
`same_dist` bp (the clones the same-state merge joins; `assignments`
gains `post_clone`, the single clone's posterior). Two clones with the
same states otherwise share each cell's posterior, few cells reach
`member_post`, both can fall below `min_size` and be dropped before the
merge, and their cells end in another clone (e.g. the normal clone).
Reference: `ref_post`, `min_ref`. Other: `doublet_delta`,
`multiplet_prob` (cells above it are excluded from discovery and
pseudobulks), `control_min_len`, `repass`, `allelic_exclude`.

`checkpoint` (resume cache, "" = off): an .rds path. When the file does
not exist the global discovery (leaves, merged clones) and the first EM
round are saved there; when it exists they are read from it and the run
continues from the within-trunk discovery (useful when rerunning the
later steps on a large sample). The cache is not keyed: reuse it only
with the same features, seed and discovery options.

Confirmation null (per-cell `confirmed` calls, after assignment;
assignment itself always uses the input fits): `ref_null` = "none" (the
input rho / depth t of `feat$null_fit`) or "rho" (opt-in: beta-binomial
rho refitted on the clean reference pool's piece votes as the
confirmation scores them, oriented leave-one-out by a clone pseudobulk;
the `call_scn` recipe on the pool's feature votes is kept as
`null_fit$rho_pool_features`); needs \>= `min_ref` pool cells. The depth
t of the confirmation is always the input one (`null_fit$tpar_confirm` =
`null_fit$tpar`; `depth_t_confirm` rows of `null_fit.tsv`); only rho can
be refitted. `confirm_null` = "generic" (null rows scored as the best of
the tumour state set, one threshold per stratum) or "matched" (each
clone hypothesis clone \| label \| shift \| f scored on the null rows
oriented by that clone's pseudobulk under the same two-state set as the
cells; one threshold per hypothesis x stratum, hypothesis-wide below
`min_null_rows`; held-out fold-CV false-confirmation rates per
hypothesis in `null_fit$confirm_null`; the default). The refitted value
is `null_fit$rho_confirm`. `confirm_sided` = the allele term of the
confirmation scores, cells and null rows alike
([`score_cells()`](https://ljwharbers.github.io/scphasedcn/reference/score_cells.md)
`sided`): "two" (default; a cell whose opposite homologue is imbalanced
also scores high) or "one" (the votes are oriented by the clone's
pseudobulk, so only the clone's major haplotype counts). The fallback
without null rows (type-7 quantile of `cell_state_scores()`) is always
one-sided.

A key that is not one of these defaults stops
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md)
with an error, and so does a value outside its range: `q`, `p_switch`,
`trunk_post`, `min_post`, `ref_post`, `max_mix_prior` and every `mix_w`
strictly between 0 and 1, `mix_prior` in \[0, 1), `multiplet_prob`,
`deep_q` and `member_post` in \[0, 1\], counts (`em_rounds`, `min_size`,
`n_perm`, `min_win`, `min_ref`, ...) whole numbers \>= 1, the caps
`max_cells`, `max_leaves` and `max_leaves_trunk` whole numbers \>= 1 or
`Inf` (no cap), `small_clone` a whole number \>= 0, `trunk_dist` and
`control_min_len` \>= 0, the other sizes / lengths \> 0, the other
tolerances \>= 0.
