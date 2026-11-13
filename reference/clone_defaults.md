# Defaults of the tumour mode

Discovery: `classes`, `depth_delta`, `cap`, `min_n`, `n_pc`, `k`,
`max_leaves`, `min_size`, `q`, `n_perm`, `max_cells`, `split_depth`
(discovery / test halves of the window depth as well as of the votes).
Merge / split tests
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
`n_tests x (1 - q)`. Hierarchy: trunks = average-linkage groups of
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
(support of balanced 2:2 runs). Assignment: `mix_prior`, `mix_w`,
`mix_margin`, `mix_rel`, `mix_min_dist` (mixture hypotheses only for
clone pairs whose states differ allelically on \>= this many bp),
`max_mix_prior`, `min_post` (confident), `member_post` (pseudobulk
members), `em_rounds`, `same_dist` (clones whose states differ on \<
this many bp within a trunk are merged). Reference: `ref_post`,
`min_ref`. Other: `doublet_delta`, `multiplet_prob` (cells above it are
excluded from discovery and pseudobulks), `control_min_len`, `repass`,
`allelic_exclude`.

## Usage

``` r
clone_defaults()
```

## Value

named list (see
[`caller_tumour()`](https://ljwharbers.github.io/scphasedcn/reference/caller_tumour.md))

## Details

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
`null_fit$rho_pool_features`); needs \>= `min_ref` pool cells.
`confirm_null` = "generic" (null rows scored as the best of the tumour
state set, one threshold per stratum) or "matched" (each clone
hypothesis clone \| label \| shift \| f scored on the null rows oriented
by that clone's pseudobulk under the same two-state set as the cells;
one threshold per hypothesis x stratum, hypothesis-wide below
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
