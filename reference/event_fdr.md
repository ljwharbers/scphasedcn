# False-event accounting of a normal-mode result: p- and q-values, expected false calls and events

Output only: no call, threshold, score or existing column changes. The
`q` column is descriptive: it is not a false discovery rate and must not
be used as a filter (see **q-value**).

## Usage

``` r
event_fdr(
  obj,
  tail = c("exponential", "none"),
  tail_n = 50,
  resolve_gains = NULL,
  gain_min_dl = NULL,
  gain_gate = NULL,
  gain_gate_q = NULL,
  alpha = 0.05
)
```

## Arguments

- obj:

  normal-mode `scphasedcn` with `scores`
  ([`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md);
  [`call_scn()`](https://ljwharbers.github.io/scphasedcn/reference/call_scn.md)
  with `keep_scores`)

- tail:

  "exponential" (tail model beyond the top `tail_n` null scores) or
  "none" (the conformal p-value only)

- tail_n:

  null scores above the tail anchor (per calibration group)

- resolve_gains, gain_min_dl:

  the gain gate of the run
  ([`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
  arguments; default the run's recorded `params$normal`, else the
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md)
  defaults)

- gain_gate, gain_gate_q:

  the gain depth gate of the run (`"dl"` or `"z"`,
  [`caller_normal()`](https://ljwharbers.github.io/scphasedcn/reference/caller_normal.md);
  default the run's recorded `params$normal`, else `"z"` when the calls
  carry the gate columns of
  [`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md)
  (`gate_thr`), else `"dl"`; `gain_gate_q` default 0.95). The `"z"` gate
  uses the run's depth null (`null_fit$tpar`)

- alpha:

  q-value level of the descriptive summary count (`n_q_below` in `info`;
  not a false discovery rate)

## Value

list: `calls` (the result's calls with `p`, `q` (descriptive, not a
false discovery rate), `null_borrowed` on the tested per-cell rows, NA
elsewhere), `events` (the result's events with `null_cell`, `p`, `q` (of
the top region) and `null_borrowed`; NA for pool events),
`expected_false` (class, state, chrom_group, null_borrowed, n_rows,
n_tested, n_called, n_events, n_null, n_null_called, n_null_events,
null_per_1000, expected_false_calls (+ `_lo`, `_hi`),
expected_false_events (+ `_lo`, `_hi`), false_frac_calls,
false_frac_events), `groups` (per calibration group: class, stratum,
n_null, p_floor = 1 / (n + 1), borrowed) and `info` (m, tail, tail_n,
alpha, the gain gate)

## Details

**Tested rows.** The per-cell rows of the classes that are reportable
per cell
([`call_policy()`](https://ljwharbers.github.io/scphasedcn/reference/call_policy.md):
chromosomes, arms, 50 Mb tiles, named events) with a threshold, plus the
20 Mb rows that
[`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md)
promoted to reportable. Each row is a hypothesis "this cell-region is
null". Its test statistic is the score of a *reportable* call,
`s' = s_comb` when the row would be reportable if called (its best state
is reportable for its class, and a chromosome / event gain passes the
depth gate of the run's
[`resolve_nested()`](https://ljwharbers.github.io/scphasedcn/reference/resolve_nested.md):
`dl_gain > gain_min_dl` with `gain_gate = "dl"`, z above its
null-quantile threshold with the opt-in `gain_gate = "z"`; a called row
uses its final `reportable` flag, so a gain that lost to a nested loss
or failed the gate counts as not reportable and a promoted row as
reportable), else `-Inf`. A gain the gate removes thus stays in the
tested family with `p = 1`, and the null rows are scored by the same
rule, so the gate is applied to the calls, their p-values and the
expected-false accounting alike.

**p-value.** Against the null rows (`is_null`) of the row's threshold
group (class x depth stratum, x vote stratum with `calib = "votes"`,
with the same pooled fallback as
[`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md))
in the cross-validation folds the row's threshold was trained on (all
folds but its own): `p = (1 + #{null s' >= s'}) / (n + 1)` (a conformal
p-value; `p = 1` for `s' = -Inf`). With `n` null rows p cannot go below
`1 / (n + 1)` (the resolution limit; ~1e-4 to 1e-5 per group on real
data), and reportable null calls are rare, so the empirical p of every
strong call sits at that floor. Under BH with m tests a single event
then cannot reach a q below about `m / (n + 1)`, which is several times
1 (m is larger than n: all cells x all tested rows against the null rows
of one class x stratum), and the q of a strong event depends only on how
many other events are as strong. `tail = "exponential"` (default)
therefore extends the p-value beyond the anchor u, the (`tail_n` + 1)-th
largest null `s_comb` of the group (whatever the state), with an
exponential tail fitted to the top `tail_n` null scores (peaks over
threshold; beta = their mean excess over u, the maximum-likelihood fit)
times the reportable share of that tail:
`p = (a + 1) / (n + 1) exp(-(s' - u) / beta)`, a = reportable rows among
the top `tail_n`; the two pieces meet at u. The exponential tail is
conservative when the null is correctly specified and light-tailed
(near-Gaussian allelic log-likelihood ratios, the bounded Student-t
depth term); it is a model, not a guarantee against artefacts the null
cells do not carry.

**Relative to the selected null.** The p-values, and the tail beyond the
largest null scores, are relative to the *selected* null cells
([`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)).
The selection removes the cells with the heaviest tails (e.g. cells with
a higher allelic overdispersion than modelled), so for cells like those
the p-values are too small, and the tail model extrapolates a truncated
distribution: in such a sample many rows can reach a small q with the
exponential tail and none with `tail = "none"`. Compare the two when the
null is in doubt. Rows on chromosomes without null rows (chrX, whose
rows are never null rows in normal mode) are compared with the autosomal
null of their class x stratum, as their thresholds are
(`null_borrowed = TRUE`).

**Named events.** A class without null rows (named `event` regions)
takes the thresholds of its `sd_class`
([`budget_thresholds()`](https://ljwharbers.github.io/scphasedcn/reference/budget_thresholds.md)),
and its p-values are calibrated on the same rows: the null rows of the
sd_class group (and folds) its threshold was trained on, each scored as
an event row would be (`s_comb` when its best state is reportable for an
event and a gain passes the depth gate, else `-Inf`; the sd_class's own
calling policy differs, e.g. arm gains are not reportable per cell).
`null_borrowed = TRUE`; in `groups` such a group has the sd_class as
`class` and `borrowed = TRUE`. The p-values are as good as the
assumption behind the borrowed thresholds: that the event's scores under
the null look like those of its sd_class (similar size and allele
votes).

**q-value.** Benjamini-Hochberg over all tested rows (m = every tested
row, called or not; rows with `s' = -Inf` have p = 1 and count in m).
Within a class x stratum and fold, q is monotone in the score. `q` is a
descriptive column: it is not a false discovery rate and must not be
used as a filter. It is too small in two independent ways. (1) The unit
is the row, not the event: one real event typically rejects several
nested rows (its chromosome, arms, tiles), which loosens the cut-off for
single-row false events, so the false share among events (an event
carries the q of its top region) with `q` below a cut can be well above
the cut (in a simulation with 300 planted chromosome losses in 3,000
cells, 12 % of the events with q \< 0.05 were false). (2) The p-values
are relative to the selected null cells and the default tail is fitted
to them (see above; on a simulated null in which 10 % of the cells had a
higher allelic overdispersion, 28 events reached q \< 0.05, all false,
and none with `tail = "none"`). No event-level q-value is computed.

**Expected false calls and events** (`expected_false`, one row per class
x state x chromosome group `A` (autosomes) / `X` / `Y`): per depth
stratum, the realised rate of reportable calls of that state among the
null rows of the class (fold-cross-validated thresholds, the calling
policy and the gain gate applied, as for every cell), times the number
of per-cell rows of the class, summed over the strata
(`expected_false_calls`); the same with the null cells' events whose top
region is a null row of that class (`expected_false_events`, so nested
calls count once). Every row of the class counts as a null hypothesis
(no estimate of the null fraction), and the null rate is measured on
cells *selected* for having no chromosome / arm call at a stricter
threshold
([`select_null_cells()`](https://ljwharbers.github.io/scphasedcn/reference/select_null_cells.md)),
whose own rows are cleaner than those of an unselected normal cell: both
estimates are lower bounds of the false load. chrX rows borrow the
autosomal rate (`null_borrowed`), which ignores X-inactivation driven
allelic imbalance. A class without null rows (named events) takes the
rate of its sd_class's null rows called as event rows (per sd_class x
stratum; that call rate also stands for its events).
[`null_selection_check()`](https://ljwharbers.github.io/scphasedcn/reference/null_selection_check.md)
measures how much the selection lowers the rate.
