# Per-state log-likelihood ratios and scores per cell x region

Per state (vs neutral): depth term =
[`ll_depth_shift()`](https://ljwharbers.github.io/scphasedcn/reference/ll_depth_shift.md)
with `mu = shift / pmax(ifelse(is.na(sd_null), 1, sd_null), sd_floor)`
under the Student-t null of the row's class x stratum
([`match_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/match_depth_null.md))
or N(0, 1); allelic term = best of the state's rows of
[`ll_allelic()`](https://ljwharbers.github.io/scphasedcn/reference/ll_allelic.md)
(two-sided, or one-sided with `sided = "one"`); combined = sum (missing
channels contribute 0). A state whose rows have different shifts (tumour
`amp`) takes the best row of depth + alleles and reports that row's two
terms.

## Usage

``` r
score_cells(
  d,
  states = caller_states("normal"),
  tpar = NULL,
  rho,
  eps = params$eps,
  depth_model = c("t", "gauss"),
  params = caller_defaults(),
  keep_rows = FALSE,
  sided = c("two", "one")
)
```

## Arguments

- d:

  table with `z`, `sd_null`, `nA`, `n_hap` (or `nB`), and `class`,
  `stratum` for the t model

- states:

  [`caller_states()`](https://ljwharbers.github.io/scphasedcn/reference/caller_states.md)
  table

- tpar:

  [`fit_depth_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_depth_null.md)
  output (t model)

- rho:

  beta-binomial overdispersion of the null allele votes
  ([`fit_allelic_null()`](https://ljwharbers.github.io/scphasedcn/reference/fit_allelic_null.md)):
  one value, or one per row (e.g. per chromosome,
  [`rho_by_chrom()`](https://ljwharbers.github.io/scphasedcn/reference/rho_by_chrom.md));
  rows are then scored in groups of equal rho

- eps:

  allele error floor

- depth_model:

  "t" or "gauss"

- params:

  [`caller_defaults()`](https://ljwharbers.github.io/scphasedcn/reference/caller_defaults.md)
  (`sd_floor`, `eps`)

- keep_rows:

  also return per-row allelic LLRs `ar_<row_id>` (and per-row depth LLRs
  `dr_<row_id>` of multi-shift states); needed by
  [`pool_scores()`](https://ljwharbers.github.io/scphasedcn/reference/pool_scores.md)

- sided:

  allelic term: `"two"` (default) = the better of `f` and `1 - f` per
  allelic row (votes in an arbitrary orientation); `"one"` = `f` on
  haplotype A only (votes already oriented to the state's major
  haplotype, e.g. by a clone pseudobulk in the tumour-mode
  confirmation), so imbalance of the opposite homologue scores negative.
  Applies to every allelic row, including the per-row columns of
  `keep_rows`

## Value

copy of `d` with `dl_<state>`, `al_<state>`, `cl_<state>` per
non-neutral state, `s_comb` (best combined LLR), `state_comb` (argmax,
ties to the first state in table order), `s_depth` (\|z\|),
`state_depth` (sign), `s_allelic` (best allelic LLR), `state_allelic`
(best allelic class, e.g. `loss_or_cnloh`)
