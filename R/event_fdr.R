# event_fdr.R
# What: the false-event accounting of normal mode, next to null_false_calls() (call_scn.R): per-cell empirical
#       (conformal) p-values of the reportable per-cell rows against the cross-validation null rows their threshold
#       was trained on, row-level Benjamini-Hochberg q-values (descriptive, not a false discovery rate), the expected
#       false calls and events per class x state
#       (event_fdr), and a selection-free check of the per-cell null (null_selection_check).
# Why: normal mode calls at a per cell-region budget (1 / 1,000 null rows); a cell has ~100 reportable rows, so a
#      run's per-cell events carry a false-event load that the counts of calls and events do not show. The null cells
#      are selected for having no calls and the budget is verified on the same cells, so a check on rows that the
#      selection did not see is needed as well. Output only: no call, threshold or score changes.
# Inputs: a normal-mode `scphasedcn` with its scores ([caller_normal()], [call_scn()] before the scores are dropped);
#         a features long table for the selection check.
# Outputs: calls / events with p, q, null_cell, null_borrowed; the expected-false table (`expected_false.tsv`); the
#          selection check table (`null_selection_check.tsv`).

utils::globalVariables(c("p", "rep_", "s_", "key_", "grp_", "borrowed_", "n_rows", "n_tested", "n_events", "n_null_called",
                         "n_null_events", "expected_false_calls", "expected_false_events", "chrom_group", "n_exceed",
                         "insample_rows", "insample_exceed", "null_borrowed", "dl_gain", "top_region", "promoted_", "gkey"))

# ---- calibration groups ----------------------------------------------------------------------------------

# The null rows each scored row's threshold was trained on: the group class x stratum (x hstratum) of the row and
# the null rows [budget_thresholds()] trained that group's threshold on ([threshold_null_rows()]: its own null rows,
# or, with fewer than `min_null`, those of the first coarser level (drop hstratum, then stratum) with `min_null`
# rows, else of the coarsest). A row whose group has no null rows at all (a class without null rows, e.g. named
# events) takes the group of its `sd_class` and is flagged `borrowed`.
# Returns list(key = calibration key per row ("<level>\r<group key at that level>"), borrowed = logical per row,
#              rows = per key the indices (rows of sc) of its null rows)
calibration_keys = function(sc, min_null) {
  by = c("class", "stratum", if ("hstratum" %in% names(sc)) "hstratum")
  isn = sc$is_null %in% TRUE
  keys_for = function(cls) {
    x = copy(sc[, by, with = FALSE])
    x[, class := cls]
    tg = threshold_null_rows(x, by, isn, min_null)
    gk = vapply(seq_len(nrow(tg$groups)), function(g)
      paste0(tg$level[g], "\r", group_key(tg$groups[g], tg$levels[[tg$level[g]]])), character(1))
    gi = match(group_key(x, by), tg$groups$key_)
    rows = tg$use[!duplicated(gk)]
    names(rows) = gk[!duplicated(gk)]
    list(key = gk[gi], rows = rows)
  }
  k1 = keys_for(as.character(sc$class))
  key = k1$key
  rows = k1$rows
  has = lengths(rows)[key] > 0
  borrowed = rep(FALSE, nrow(sc))
  if (any(!has) && "sd_class" %in% names(sc)) {
    k2 = keys_for(as.character(sc$sd_class))
    b = !has & lengths(k2$rows)[k2$key] > 0 & sc$sd_class != sc$class
    key[b] = k2$key[b]
    borrowed[b] = TRUE
    add = setdiff(unique(k2$key[b]), names(rows))
    rows[add] = k2$rows[add]
  }
  list(key = key, borrowed = borrowed, rows = rows)
}

# Empirical (conformal) p-value of reportable scores `s` (-Inf = not reportable) against the calibration rows: `cal`
# their reportable scores (-Inf = not reportable), `plain` their scores whatever the state. Up to the tail anchor
# p = (1 + #{cal >= s}) / (n + 1), n = calibration rows; s = -Inf gives 1. With `tail_n` > 0 and more than `tail_n`
# finite plain scores, a score above the anchor u (the (tail_n + 1)-th largest plain score) takes an exponential tail
# fitted to the top `tail_n` plain scores, times the reportable share of that tail:
# p = (a + 1) / (n + 1) x exp(-(s - u) / beta), beta = their mean excess over u (the maximum-likelihood exponential
# fit), a = how many of them are reportable. Continuous at u (both sides (a + 1) / (n + 1)) and decreasing in s.
# Both pieces are relative to the calibration rows, i.e. to the selected null cells (see event_fdr).
conformal_p = function(s, cal, plain, tail_n = 0) {
  n = length(cal)
  p = rep(1, length(s))
  ok = is.finite(s)
  if (!any(ok)) return(p)
  xs = sort(cal[is.finite(cal)])
  ge = length(xs) - findInterval(s[ok], xs, left.open = TRUE)
  p[ok] = (1 + ge) / (n + 1)
  pl = sort(plain[is.finite(plain)], decreasing = TRUE)
  if (tail_n > 0 && length(pl) > tail_n) {
    u = pl[tail_n + 1L]
    beta = mean(pl[seq_len(tail_n)] - u)
    if (is.finite(beta) && beta > 0) {
      a = sum(is.finite(cal) & cal > u)
      t = which(ok)[s[ok] > u]
      if (length(t)) p[t] = (a + 1) / (n + 1) * exp(-(s[t] - u) / beta)
    }
  }
  p
}

# ---- false-event accounting ----------------------------------------------------------------------------

#' False-event accounting of a normal-mode result: p- and q-values, expected false calls and events
#'
#' Output only: no call, threshold, score or existing column changes. The `q` column is descriptive: it is not a
#' false discovery rate and must not be used as a filter (see **q-value**).
#'
#' **Tested rows.** The per-cell rows of the classes that are reportable per cell ([call_policy()]: chromosomes, arms,
#' 50 Mb tiles, named events) with a threshold, plus the 20 Mb rows that [resolve_nested()] promoted to reportable.
#' Each row is a hypothesis "this cell-region is null". Its test statistic is the score of a *reportable* call,
#' `s' = s_comb` when the row would be reportable if called (its best state is reportable for its class, and a
#' chromosome / event gain passes the depth gate `dl_gain > gain_min_dl` of [resolve_nested()]; a called row uses its
#' final `reportable` flag, so a gain that lost to a nested loss counts as not reportable and a promoted row as
#' reportable), else `-Inf`.
#'
#' **p-value.** Against the null rows (`is_null`) of the row's threshold group (class x depth stratum, x vote stratum
#' with `calib = "votes"`, with the same pooled fallback as [budget_thresholds()]) in the cross-validation folds the
#' row's threshold was trained on (all folds but its own): `p = (1 + #{null s' >= s'}) / (n + 1)` (a conformal
#' p-value; `p = 1` for `s' = -Inf`). With `n` null rows p cannot go below `1 / (n + 1)` (the resolution limit; ~1e-4
#' to 1e-5 per group on real data), and reportable null calls are rare, so the empirical p of every strong call sits
#' at that floor. Under BH with m tests a single event then cannot reach a q below about `m / (n + 1)`, which is
#' several times 1 (m is larger than n: all cells x all tested rows against the null rows of one class x stratum),
#' and the q of a strong event depends only on how many other events are as strong. `tail = "exponential"` (default)
#' therefore extends the p-value beyond the anchor u, the (`tail_n` + 1)-th largest null `s_comb` of the group
#' (whatever the state), with an exponential tail fitted to the top `tail_n` null scores (peaks over threshold; beta
#' = their mean excess over u, the maximum-likelihood fit) times the reportable share of that tail:
#' `p = (a + 1) / (n + 1) exp(-(s' - u) / beta)`, a = reportable rows among the top `tail_n`; the two pieces meet at u.
#' The exponential tail is conservative when the null is correctly specified and light-tailed (near-Gaussian
#' allelic log-likelihood ratios, the bounded Student-t depth term); it is a model, not a guarantee against
#' artefacts the null cells do not carry.
#'
#' **Relative to the selected null.** The p-values, and the tail beyond the largest null scores, are relative to the
#' *selected* null cells ([select_null_cells()]). The selection removes the cells with the heaviest tails (e.g.
#' cells with a higher allelic overdispersion than modelled), so for cells like those the p-values are too small, and
#' the tail model extrapolates a truncated distribution: in such a sample many rows can reach a small q with the
#' exponential tail and none with `tail = "none"`. Compare the two when the null is in doubt.
#' Rows on chromosomes without null rows (chrX, whose rows are never null rows in normal mode) are compared with the
#' autosomal null of their class x stratum, as their thresholds are (`null_borrowed = TRUE`); so are rows of a class
#' without null rows that take the null of their `sd_class`.
#'
#' **q-value.** Benjamini-Hochberg over all tested rows (m = every tested row, called or not; rows with `s' = -Inf`
#' have p = 1 and count in m). Within a class x stratum and fold, q is monotone in the score. `q` is a descriptive
#' column: it is not a false discovery rate and must not be used as a filter. It is too small in two independent
#' ways. (1) The unit is the row, not the event: one real event typically rejects several nested rows (its
#' chromosome, arms, tiles), which loosens the cut-off for single-row false events, so the false share among events
#' (an event carries the q of its top region) with `q` below a cut can be well above the cut (in a simulation with
#' 300 planted chromosome losses in 3,000 cells, 12 % of the events with q < 0.05 were false). (2) The p-values are
#' relative to the selected null cells and the default tail is fitted to them (see above; on a simulated null in
#' which 10 % of the cells had a higher allelic overdispersion, 28 events reached q < 0.05, all false, and none
#' with `tail = "none"`). No event-level q-value is computed.
#'
#' **Expected false calls and events** (`expected_false`, one row per class x state x chromosome group `A`
#' (autosomes) / `X` / `Y`): per depth stratum, the realised rate of reportable calls of that state among the null
#' rows of the class (fold-cross-validated thresholds, the calling policy and the gain gate applied, as for every
#' cell), times the number of per-cell rows of the class, summed over the strata (`expected_false_calls`); the same
#' with the null cells' events whose top region is a null row of that class (`expected_false_events`, so nested
#' calls count once). Every row of the class counts as a null hypothesis (no estimate of the null fraction), and
#' the null rate is measured on cells *selected* for having no chromosome / arm call at a stricter threshold
#' ([select_null_cells()]), whose own rows are cleaner than those of an unselected normal cell: both estimates are
#' lower bounds of the false load. chrX rows borrow the autosomal rate (`null_borrowed`), which ignores
#' X-inactivation driven allelic imbalance. [null_selection_check()] measures how much the selection lowers the rate.
#' @param obj normal-mode `scphasedcn` with `scores` ([caller_normal()]; [call_scn()] with `keep_scores`)
#' @param tail "exponential" (tail model beyond the top `tail_n` null scores) or "none" (the conformal p-value only)
#' @param tail_n null scores above the tail anchor (per calibration group)
#' @param resolve_gains,gain_min_dl the gain gate of the run ([caller_normal()] arguments; default the run's
#'   recorded `params$normal`, else the [caller_normal()] defaults)
#' @param alpha q-value level of the descriptive summary count (`n_q_below` in `info`; not a false discovery rate)
#' @return list: `calls` (the result's calls with `p`, `q` (descriptive, not a false discovery rate), `null_borrowed`
#'   on the tested per-cell rows, NA
#'   elsewhere), `events` (the result's events with `null_cell`, `p`, `q` (of the top region) and `null_borrowed`;
#'   NA for pool events), `expected_false` (class, state, chrom_group, null_borrowed, n_rows, n_tested, n_called,
#'   n_events, n_null, n_null_called, n_null_events, null_per_1000, expected_false_calls (+ `_lo`, `_hi`),
#'   expected_false_events (+ `_lo`, `_hi`), false_frac_calls, false_frac_events), `groups` (per calibration group:
#'   class, stratum, n_null, p_floor = 1 / (n + 1), borrowed) and `info` (m, tail, tail_n, alpha, the gain gate)
#' @export
event_fdr = function(obj, tail = c("exponential", "none"), tail_n = 50, resolve_gains = NULL, gain_min_dl = NULL, alpha = 0.05) {
  tail = match.arg(tail)
  if (!inherits(obj, "scphasedcn") || !identical(obj$mode, "normal")) stop("event_fdr: needs a normal-mode scphasedcn result")
  if (is.null(obj$scores)) stop("event_fdr: the result has no scores (call_scn(keep_scores = TRUE) or caller_normal())")
  if (is.null(obj$calls)) stop("event_fdr: the result has no calls")
  na = obj$params$normal %||% list()
  resolve_gains = resolve_gains %||% na$resolve_gains %||% TRUE
  gain_min_dl = gain_min_dl %||% na$gain_min_dl %||% 0
  params = obj$params %||% caller_defaults()
  min_null = params$min_null_rows %||% caller_defaults()$min_null_rows
  sc = as.data.table(obj$scores)
  for (cn in c("cell", "region", "class", "stratum", "s_comb", "is_null"))
    if (!cn %in% names(sc)) stop("event_fdr: scores lack ", cn)
  if (!"chrom" %in% names(sc)) sc[, chrom := region_chrom(region, class)]
  if (!"fold" %in% names(sc)) sc[, fold := 0L]
  calls = as.data.table(obj$calls)
  ci = which(calls$level == "cell")
  x = calls[ci, .(cell, region, class, state, score, thr, called, reportable,
                  resolved = if ("resolved" %in% names(calls)) resolved else NA_character_)]
  k = match(paste(x$cell, x$region, sep = "\r"), paste(sc$cell, sc$region, sep = "\r"))
  if (anyNA(k)) stop("event_fdr: ", sum(is.na(k)), " cell-level calls have no scored row")
  x[, `:=`(is_null = sc$is_null[k] %in% TRUE, fold = as.integer(sc$fold[k]), stratum = as.character(sc$stratum[k]),
           chrom = sc$chrom[k], dl_gain = if ("dl_gain" %in% names(sc)) sc$dl_gain[k] else NA_real_,
           sd_class = if ("sd_class" %in% names(sc)) as.character(sc$sd_class[k]) else class)]
  # reportable if called: the final flag of called rows; policy flag + gain depth gate of uncalled rows
  x[, rep_ := reportable %in% TRUE]
  if (isTRUE(resolve_gains))
    x[!(called %in% TRUE) & class %in% c("chrom", "event") & state == "gain" & !is.na(dl_gain) & dl_gain <= gain_min_dl, rep_ := FALSE]
  x[, s_ := fifelse(rep_ & !is.na(score), score, -Inf)]
  x[, promoted_ := !is.na(resolved) & startsWith(resolved, "replaces")]
  sex = obj$null_fit$sex
  pol = call_policy("normal", sex)
  rep_classes = unique(pol[level == "cell" & reportable == TRUE, class])
  tested = (x$class %in% rep_classes | x$promoted_) & !is.na(x$thr)
  # calibration groups on the scored rows of the calls (the thresholds' groups, from all null rows of the run)
  cols = intersect(c("class", "sd_class", "stratum", "hstratum"), names(sc))
  ck = calibration_keys(sc[, c(cols, "is_null"), with = FALSE], min_null)
  x[, `:=`(key_ = ck$key[k], borrowed_ = ck$borrowed[k])]
  null_chroms = unique(sc$chrom[sc$is_null %in% TRUE])
  x[, null_borrowed := borrowed_ | !chrom %in% null_chroms]
  # p-values per calibration group x fold: the group's null rows (those its threshold was trained on) of the other
  # folds (all folds for rows without one)
  folds = sort(unique(x$fold[!is.na(x$fold) & x$fold != 0L]))
  x[, p := NA_real_]
  tn = if (tail == "exponential") tail_n else 0
  xi = rep(NA_integer_, nrow(sc)); xi[k] = seq_along(k)      # the call row of each scored row
  grp = list()
  for (g in unique(x$key_[tested])) {
    nul = xi[ck$rows[[g]]]
    nul = nul[!is.na(nul)]
    tg = which(x$key_ == g & tested)
    grp[[length(grp) + 1]] = data.table(gkey = g, n_null = length(nul))
    for (fo in unique(x$fold[tg])) {
      t_ = tg[x$fold[tg] %in% fo]
      cal = if (!is.na(fo) && fo %in% folds) nul[!is.na(x$fold[nul]) & x$fold[nul] != fo] else nul
      set(x, i = t_, j = "p", value = conformal_p(x$s_[t_], x$s_[cal], x$score[cal], tn))
    }
  }
  m = sum(tested & !is.na(x$p))
  x[, q := NA_real_]
  it = which(tested & !is.na(x$p))
  if (length(it)) set(x, i = it, j = "q", value = stats::p.adjust(x$p[it], "BH", n = m))
  x[!tested, null_borrowed := NA]
  # calls: new columns on the cell rows
  out_calls = copy(calls)
  for (cn in c("p", "q", "null_borrowed")) {
    v = if (cn == "null_borrowed") rep(NA, nrow(calls)) else rep(NA_real_, nrow(calls))
    v[ci] = x[[cn]]
    set(out_calls, j = cn, value = v)
  }
  # events: null_cell, p / q of the top region, null_borrowed
  ev = if (is.null(obj$events)) data.table() else copy(as.data.table(obj$events))
  nullc = as.data.table(obj$cells)[null_cell %in% TRUE, as.character(cell)]
  if (nrow(ev)) {
    isc = ev$level == "cell"
    ke = match(paste(ev$cell, ev$top_region, sep = "\r"), paste(x$cell, x$region, sep = "\r"))
    ke[!isc] = NA_integer_
    ev[, null_cell := fifelse(isc, cell %in% nullc, NA)]
    ev[, `:=`(p = x$p[ke], q = x$q[ke], null_borrowed = x$null_borrowed[ke])]
  } else ev = data.table(cell = character(), chrom = character(), start = numeric(), end = numeric(), state = character(),
                         top_region = character(), score = numeric(), level = character(), top_class = character(),
                         n_regions = integer(), null_cell = logical(), p = numeric(), q = numeric(), null_borrowed = logical())
  # expected false calls / events per class x state x chromosome group
  ef = expected_false_table(x, tested, ev, pol)
  # calibration groups: class / stratum of each key, n_null, resolution
  groups = rbindlist(grp)
  if (nrow(groups)) {
    kk = match(groups$gkey, x$key_)
    lev = as.integer(sub("\r.*$", "", groups$gkey))
    hs = "hstratum" %in% names(sc)
    bw = x$borrowed_[kk]
    cls = x$class[kk]
    if (any(bw)) cls[bw] = as.character(sc$sd_class[k[kk[bw]]])
    groups[, `:=`(class = cls, stratum = fifelse(lev == 1L | (lev == 2L & hs), x$stratum[kk], "all"),
                  borrowed = bw, p_floor = 1 / (n_null + 1))]
    if (hs) groups[, hstratum := fifelse(lev == 1L, as.character(sc$hstratum[k[kk]]), NA_character_)]
    groups[, gkey := NULL]
    setcolorder(groups, intersect(c("class", "stratum", "hstratum", "n_null", "p_floor", "borrowed"), names(groups)))
    groups = groups[order(match(class, region_classes()), stratum)]
  }
  ce = ev[level == "cell"]
  info = list(m = m, tail = tail, tail_n = tail_n, alpha = alpha, resolve_gains = resolve_gains, gain_min_dl = gain_min_dl,
              n_q_below = sum(ce$q < alpha, na.rm = TRUE), n_cell_events = nrow(ce), n_cell_events_null_cells = sum(ce$null_cell %in% TRUE),
              n_null_cells = length(nullc), p_floor_max = if (nrow(groups)) max(groups$p_floor) else NA_real_)
  list(calls = out_calls[], events = ev[], expected_false = ef, groups = groups, info = info)
}

# expected false calls / events per class x state x chromosome group (see event_fdr): x = the cell-level calls with
# is_null, stratum, chrom, rep_, key_; ev = the events with null_cell
expected_false_table = function(x, tested, ev, pol) {
  x = copy(x)
  x[, grp_ := sex_group(chrom)]
  x[, fin_ := called %in% TRUE & reportable %in% TRUE]
  combos = unique(rbind(pol[level == "cell" & reportable == TRUE & class %in% unique(x$class[tested]), .(class, state)],
                        x[fin_ == TRUE, .(class, state)]))
  ce = ev[level == "cell"]
  ke = match(paste(ce$cell, ce$top_region, sep = "\r"), paste(x$cell, x$region, sep = "\r"))
  ce = ce[, .(class = x$class[ke], state, stratum = x$stratum[ke], grp_ = x$grp_[ke], top_null = x$is_null[ke] %in% TRUE)]
  # rate class: the class itself, or its sd_class when the class has no null rows (borrowed)
  ncls = unique(x$class[x$is_null])
  rows = x[!is.na(thr), .(n_rows = .N), by = .(class, stratum, grp_)]
  nul = x[is_null == TRUE, .(n_null = .N), by = .(class, stratum)]
  out = list()
  for (i in seq_len(nrow(combos))) {
    cl = combos$class[i]; st = combos$state[i]
    rc = cl
    if (!cl %in% ncls) {
      sdc = x$sd_class[match(cl, x$class)]
      rc = if (!is.na(sdc) && sdc %in% ncls) sdc else NA_character_
    }
    nr = rows[class == cl]
    kn = x[is_null == TRUE & class %in% rc & fin_ & state == st, .N, by = stratum]
    ke_ = ce[top_null & class %in% rc & state == st, .N, by = stratum]
    nn = nul[class %in% rc]
    for (g in sort(unique(nr$grp_))) {
      r = nr[grp_ == g]
      tot_null = sum(nn$n_null); tot_k = sum(kn$N); tot_e = sum(ke_$N)
      rate = function(cnt, s) {
        n0 = nn$n_null[match(s, nn$stratum)]
        c0 = cnt$N[match(s, cnt$stratum)]; c0[is.na(c0)] = 0
        pooled = if (tot_null > 0) sum(cnt$N) / tot_null else NA_real_
        fifelse(!is.na(n0) & n0 > 0, c0 / n0, pooled)
      }
      e_calls = sum(r$n_rows * rate(kn, r$stratum))
      e_events = sum(r$n_rows * rate(ke_, r$stratum))
      scale0 = if (tot_null > 0) sum(r$n_rows) / tot_null else NA_real_
      ci = function(e, cnt) {
        sc_ = if (cnt > 0) e / cnt else scale0
        c(if (cnt > 0) sc_ * stats::qgamma(0.025, cnt) else 0, sc_ * stats::qgamma(0.975, cnt + 1))
      }
      cc = ci(e_calls, tot_k); cv = ci(e_events, tot_e)
      xg = x[class == cl & grp_ == g]
      n_called = sum(xg$fin_ & xg$state == st)
      n_ev = nrow(ce[class == cl & state == st & grp_ == g])
      out[[length(out) + 1]] = data.table(
        class = cl, state = st, chrom_group = g, null_borrowed = g != "A" | !identical(rc, cl),
        n_rows = sum(r$n_rows), n_tested = sum(xg$rep_ & xg$state == st & !is.na(xg$thr)), n_called = n_called, n_events = n_ev,
        n_null = as.integer(tot_null), n_null_called = as.integer(tot_k), n_null_events = as.integer(tot_e),
        null_per_1000 = if (tot_null > 0) 1000 * tot_k / tot_null else NA_real_,
        expected_false_calls = e_calls, expected_false_calls_lo = cc[1], expected_false_calls_hi = cc[2],
        expected_false_events = e_events, expected_false_events_lo = cv[1], expected_false_events_hi = cv[2],
        false_frac_calls = if (n_called > 0) min(1, e_calls / n_called) else NA_real_,
        false_frac_events = if (n_ev > 0) min(1, e_events / n_ev) else NA_real_)
    }
  }
  res = rbindlist(out)
  if (!nrow(res)) return(res)
  res[order(match(class, region_classes()), match(state, c("loss", "gain", "cnloh")), chrom_group)]
}
utils::globalVariables(c("fin_", "top_null", "N", "score", "called", "thr"))

# ---- selection-free null check ---------------------------------------------------------------------------

#' Selection-free check of the per-cell null
#'
#' The null cells of normal mode are selected for having no chromosome / arm call ([select_null_cells()]), the
#' thresholds are trained on their rows and the budget is verified on the same rows by fold cross-validation, so a
#' selection that truncates the null tail (a too thin parametric null, cells with a noisier than modelled genome)
#' goes unseen. This check splits the autosomes into two halves (odd / even chromosome number), selects the null
#' cells on the rows of one half only, trains the per-cell thresholds as usual on those rows ([caller_normal()] with
#' the selected cells as null cells and the half as `null_chroms`: depth t null, rho, fold-cross-validated budget
#' thresholds, the run's `calib` options), and counts the exceedances (`s_comb` above its threshold) on the null
#' regions of the OTHER half of the same cells: rows that neither the selection nor the thresholds saw. The two halves
#' are swapped and the counts combined.
#'
#' Cells are selected on one half only, so a cell that carries a real event on the measured half stays in and the
#' rows of that event count as exceedances: `per_1000` is the null exceedance rate plus the rate of real events on
#' the measured half, an upper bound of the null rate. A rate clearly above the nominal `1 - null_q` (and above the
#' in-sample cross-validated rate on the selection half, `insample_per_1000`) has two possible causes: a selected
#' null that is cleaner than the cells it is applied to (the per-cell false load is then underestimated), or real
#' events in the cells. Read it next to the share of cells with a called event: in a sample without real events
#' the excess measures the selection, in a sample with them it is an upper bound. The halves carry different
#' chromosomes, so a class whose regions differ in their null spread (e.g. chromosomes of very different size) also
#' shows up.
#' Confidence intervals are exact binomial (Clopper-Pearson) over rows; rows of one cell overlap (a chromosome and its
#' arms), so they are somewhat too narrow. Costs about one more [caller_normal()] scoring pass (two half passes
#' without pools). With `calib_rho = "chrom"` the measured half has no null rows and takes the global rho.
#' @param feat features long table (as for [caller_normal()])
#' @param params [caller_defaults()]
#' @param multiplets cells never selected as null cells
#' @param regions regions table (default [regions_from_features()])
#' @param depth_model "t" or "gauss"
#' @param select_args extra arguments for [select_null_cells()] (its `chroms` is set per half)
#' @param null_chroms the chromosomes whose rows can be null rows in the run; their numbered autosomes are split
#' @param false_calls optional [null_false_calls()] `rates` of the run, added as `cv_per_1000` (the run's in-sample
#'   cross-validated rate on all null rows)
#' @param conf confidence level of the intervals
#' @return data.table class, stratum (`all` = over the strata), n_rows (measured rows), n_exceed, per_1000 (null and
#'   real-event exceedances per 1,000 measured rows: an upper bound of the null rate),
#'   lo_per_1000, hi_per_1000, nominal_per_1000, insample_rows, insample_exceed, insample_per_1000 (the
#'   cross-validated rate on the selection half's null rows of the same passes), cv_per_1000 (with `false_calls`);
#'   attribute `halves` (half, select_chroms, measure_chroms, n_cells, n_null_cells)
#' @export
null_selection_check = function(feat, params = caller_defaults(), multiplets = NULL, regions = NULL, depth_model = c("t", "gauss"),
                                select_args = list(), null_chroms = paste0("chr", 1:22), false_calls = NULL, conf = 0.95) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  depth_model = match.arg(depth_model)
  d = copy(as.data.table(feat))
  d[, cell := as.character(cell)]
  if (!"chrom" %in% names(d)) d[, chrom := region_chrom(region, class)]
  if (!"null_region" %in% names(d)) d[, null_region := TRUE]
  if (!"eligible" %in% names(d)) d[, eligible := TRUE]
  if (!"fold" %in% names(d) || anyNA(d$fold)) d[, fold := assign_folds(cell, params$n_folds, params$seed, if ("origin" %in% names(d)) origin else NULL)]
  d = d[eligible %in% TRUE]
  multiplets = intersect(as.character(multiplets), d$cell)
  auto = paste0("chr", 1:22)
  auto = auto[auto %in% (null_chroms %||% auto) & auto %in% d$chrom]
  num = as.integer(sub("^chr", "", auto))
  halves = list(odd = auto[num %% 2L == 1L], even = auto[num %% 2L == 0L])
  if (min(lengths(halves)) == 0) stop("null_selection_check: needs autosomes with odd and even numbers")
  sa0 = select_args[setdiff(names(select_args), "chroms")]
  meas = list(); ins = list(); hv = list()
  for (h in names(halves)) {
    sel_ch = halves[[h]]; oth = setdiff(auto, sel_ch)
    sel = do.call(select_null_cells, c(list(feat_long = d, params = params, exclude = multiplets, depth_model = depth_model,
                                            chroms = sel_ch), sa0))
    d[, null_cell := cell %in% sel$null_cells]
    r = caller_normal(d, params = params, multiplets = multiplets, regions = regions, null_cells = "features",
                      depth_model = depth_model, pool_classes = character(), resolve_gains = FALSE, null_chroms = sel_ch)
    sc = r$scores
    k = !is.na(sc$call_s_comb)
    mrow = k & sc$null_cell %in% TRUE & sc$null_region %in% TRUE & sc$chrom %in% oth
    irow = k & sc$is_null %in% TRUE
    meas[[h]] = sc[mrow, .(n_rows = .N, n_exceed = sum(call_s_comb)), by = .(class, stratum = as.character(stratum))]
    ins[[h]] = sc[irow, .(insample_rows = .N, insample_exceed = sum(call_s_comb)), by = .(class, stratum = as.character(stratum))]
    hv[[h]] = data.table(half = h, select_chroms = paste(sel_ch, collapse = ","), measure_chroms = paste(oth, collapse = ","),
                         n_cells = uniqueN(d$cell), n_null_cells = length(sel$null_cells))
  }
  agg = function(l, cols) {
    a = rbindlist(l)[, lapply(.SD, sum), by = .(class, stratum), .SDcols = cols]
    rbind(a, a[, c(list(stratum = "all"), lapply(.SD, sum)), by = class, .SDcols = cols])
  }
  m = agg(meas, c("n_rows", "n_exceed"))
  i = agg(ins, c("insample_rows", "insample_exceed"))
  out = merge(m, i, by = c("class", "stratum"), all = TRUE)
  for (cn in c("n_rows", "n_exceed", "insample_rows", "insample_exceed")) set(out, i = which(is.na(out[[cn]])), j = cn, value = 0L)
  a = (1 - conf) / 2
  out[, `:=`(per_1000 = fifelse(n_rows > 0, 1000 * n_exceed / n_rows, NA_real_),
             lo_per_1000 = fifelse(n_exceed > 0, 1000 * stats::qbeta(a, n_exceed, n_rows - n_exceed + 1), 0),
             hi_per_1000 = fifelse(n_rows > n_exceed, 1000 * stats::qbeta(1 - a, n_exceed + 1, n_rows - n_exceed), 1000),
             nominal_per_1000 = 1000 * (1 - params$null_q),
             insample_per_1000 = fifelse(insample_rows > 0, 1000 * insample_exceed / insample_rows, NA_real_))]
  out[n_rows == 0, `:=`(lo_per_1000 = NA_real_, hi_per_1000 = NA_real_)]
  if (!is.null(false_calls) && nrow(false_calls)) {
    fc = as.data.table(false_calls)
    out[, cv_per_1000 := fc$per_1000[match(paste(class, stratum), paste(fc$class, fc$stratum))]]
  } else out[, cv_per_1000 := NA_real_]
  setcolorder(out, c("class", "stratum", "n_rows", "n_exceed", "per_1000", "lo_per_1000", "hi_per_1000", "nominal_per_1000",
                     "insample_rows", "insample_exceed", "insample_per_1000", "cv_per_1000"))
  out = out[order(match(class, region_classes()), stratum != "all", stratum)]
  setattr(out, "halves", rbindlist(hv))
  out[]
}
utils::globalVariables(c("lo_per_1000", "hi_per_1000", "nominal_per_1000", "insample_per_1000", "cv_per_1000", "per_1000"))
