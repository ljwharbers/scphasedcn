# caller.R
# What: the scorer of the copy-number caller. Per cell x region state
#       log-likelihood ratios (depth location shift under a Student-t / Gaussian null + oriented
#       beta-binomial allele term, max over each state's allelic rows), null thresholds at the
#       false-call budget per class x depth stratum (fold-cross-validated, type-7 quantile, strict
#       `>`), per-cell calls under a calling policy, consolidation of nested calls into events,
#       pooled scores and pooled thresholds (random pools within a cell type, or chromatin kNN
#       pools), self-calibration of the null cell set in normal tissue, and the normal-tissue
#       wrapper caller_normal().
# Why: one scorer for normal tissue and tumours. With the `normal` state set its scores, cross-validated
#      thresholds, calls and pooled scores are tested against an independent reference implementation
#      (tests/testthat/helper-refscorer.R); the tumour state set adds attenuated 1:0 / 2:0, allele-specific gains
#      and amplifications.
# Inputs: features long table (vignette("scphasedcn-outputs", package = "scphasedcn")) or equivalent cell x region
#         tables.
# Outputs: scored tables, thresholds / calls tables (check_thresholds, check_calls), events,
#          pools, `scphasedcn` objects (mode normal).

utils::globalVariables(c("row_id", "shift", "fold", "is_null", "null_cell", "null_region", "sd_class",
                         "stratum", "n_hap", "z", "sd_null", "origin", "call_s_comb", "n_calls", "ev",
                         "top_region", "reportable", "thr", "q", "n_null", "pooled_fallback", "center",
                         "level", "source", "s_comb", "state_comb", "n_called", "per_1000", "rate", "key_",
                         "fold_", "lab_", "n_regions", "chrom_", "start_", "end_", "cls_", "state_",
                         "region_", "score_", "thr_", "called_", "rep_", "cell_type"))

# ---- state sets ----------------------------------------------------------------------------------

#' Copy-number state sets of the scorer
#'
#' One row per (state, allelic configuration). A state may have several rows
#' (e.g. gain = 2:1 at f 2/3 and 3:1 at f 3/4); its allelic LLR is the best
#' row (two-sided), its depth LLR the row's log2 shift.
#'   * `normal`: neutral 1:1; loss 1:0 at f = 1 - eps (shift loss);
#'     gain 2:1 (2/3) and 3:1 (3/4), both at the gain shift; cnloh 2:0 at
#'     1 - eps (shift 0). The normal-tissue state set.
#'   * `normal_att`: as `normal` with the loss and cnloh rows over `f_att` (an
#'     attenuation-aware sensitivity variant).
#'   * `tumour` (per-cell detection and assignment): loss 1:0 at f in
#'     c(0.80, f_att); cnloh 2:0 over f_att; gain 2:1 (2/3); gain_as 3:0 over
#'     f_att (gain shift); amp 3:1 (0.75), 2:2 (0.5), 4:1 (0.80) with shift
#'     `k_amp * log2(cn / 2)`.
#' @param set "normal", "normal_att" or "tumour"
#' @param params [caller_defaults()]
#' @return data.table state, row_id, nA, nB, cn, shift (log2), f (major-haplotype fraction);
#'   the first row is `neutral`
#' @export
caller_states = function(set = c("normal", "normal_att", "tumour"), params = caller_defaults()) {
  if (length(set) == 1 && set %in% c("e4", "e4_att")) set = sub("e4", "normal", set, fixed = TRUE)   # older names
  set = match.arg(set)
  eps = params$eps; sh = params$shifts; fa = params$f_att
  mk = function(state, a, b, shift, f) data.table(state = state, nA = as.integer(a), nB = as.integer(b),
                                                  cn = as.integer(a + b), shift = shift, f = f)
  amp = function(a, b, f) mk("amp", a, b, params$k_amp * log2((a + b) / 2), f)
  st = if (set %in% c("normal", "normal_att")) {
    fl = if (set == "normal_att") fa else 1 - eps
    rbind(mk("neutral", 1, 1, 0, 0.5), mk("loss", 1, 0, sh[["loss"]], fl), mk("gain", 2, 1, sh[["gain"]], 2 / 3),
          mk("gain", 3, 1, sh[["gain"]], 3 / 4), mk("cnloh", 2, 0, 0, fl))
  } else {
    rbind(mk("neutral", 1, 1, 0, 0.5), mk("loss", 1, 0, sh[["loss"]], c(0.80, fa)), mk("cnloh", 2, 0, 0, fa),
          mk("gain", 2, 1, sh[["gain"]], 2 / 3), mk("gain_as", 3, 0, sh[["gain"]], fa),
          amp(3, 1, 0.75), amp(2, 2, 0.5), amp(4, 1, 0.80))
  }
  st[, row_id := paste0(state, "_", seq_len(.N)), by = state]
  setcolorder(st, c("state", "row_id", "nA", "nB", "cn", "shift", "f"))
  st[]
}

# states without any depth shift (CN-LOH): purely allelic, undefined on depth-only (allelic_ok FALSE) rows
pure_allelic_states = function(ns) {
  sn = unique(ns$state)
  sn[vapply(sn, function(s) all(ns$shift[ns$state == s] == 0), logical(1))]
}

# states that share an identical allelic row set are indistinguishable by alleles: one label
# ("loss_or_cnloh" in the normal set), in table order
allelic_classes = function(ns) {
  sn = unique(ns$state)
  fkey = vapply(sn, function(s) paste(sprintf("%.15g", sort(unique(ns$f[ns$state == s]))), collapse = ","), character(1))
  u = unique(fkey)
  list(label = vapply(u, function(k) paste(sn[fkey == k], collapse = "_or_"), character(1)),
       states = lapply(u, function(k) sn[fkey == k]))
}

# ---- scoring ---------------------------------------------------------------------------------------

#' Per-state log-likelihood ratios and scores per cell x region
#'
#' Per state (vs neutral): depth term = [ll_depth_shift()] with
#' `mu = shift / pmax(ifelse(is.na(sd_null), 1, sd_null), sd_floor)` under the
#' Student-t null of the row's class x stratum ([match_depth_null()]) or
#' N(0, 1); allelic term = best of the state's rows of [ll_allelic()]
#' (two-sided, or one-sided with `sided = "one"`); combined = sum (missing channels contribute 0). A state whose
#' rows have different shifts (tumour `amp`) takes the best row of depth +
#' alleles and reports that row's two terms.
#' @param d table with `z`, `sd_null`, `nA`, `n_hap` (or `nB`), and `class`,
#'   `stratum` for the t model
#' @param states [caller_states()] table
#' @param tpar [fit_depth_null()] output (t model)
#' @param rho beta-binomial overdispersion of the null allele votes ([fit_allelic_null()]): one value, or one
#'   per row (e.g. per chromosome, [rho_by_chrom()]); rows are then scored in groups of equal rho
#' @param eps allele error floor
#' @param depth_model "t" or "gauss"
#' @param params [caller_defaults()] (`sd_floor`, `eps`)
#' @param keep_rows also return per-row allelic LLRs `ar_<row_id>` (and
#'   per-row depth LLRs `dr_<row_id>` of multi-shift states); needed by [pool_scores()]
#' @param sided allelic term: `"two"` (default) = the better of `f` and `1 - f` per allelic row (votes in an arbitrary
#'   orientation); `"one"` = `f` on haplotype A only (votes already oriented to the state's major haplotype, e.g. by a
#'   clone pseudobulk in the tumour-mode confirmation), so imbalance of the opposite homologue scores negative. Applies to
#'   every allelic row, including the per-row columns of `keep_rows`
#' @return copy of `d` with `dl_<state>`, `al_<state>`, `cl_<state>` per
#'   non-neutral state, `s_comb` (best combined LLR), `state_comb` (argmax,
#'   ties to the first state in table order), `s_depth` (|z|), `state_depth`
#'   (sign), `s_allelic` (best allelic LLR), `state_allelic` (best allelic
#'   class, e.g. `loss_or_cnloh`)
#' @export
score_cells = function(d, states = caller_states("normal"), tpar = NULL, rho, eps = params$eps,
                       depth_model = c("t", "gauss"), params = caller_defaults(), keep_rows = FALSE, sided = c("two", "one")) {
  depth_model = match.arg(depth_model)
  two = match.arg(sided) == "two"
  d = copy(as.data.table(d))
  n = nrow(d)
  if (!"z" %in% names(d)) d[, z := NA_real_]
  if (!"sd_null" %in% names(d)) d[, sd_null := NA_real_]
  nA = if ("nA" %in% names(d)) as.numeric(d$nA) else rep(0, n)
  nh = if ("n_hap" %in% names(d)) as.numeric(d$n_hap) else if ("nB" %in% names(d)) nA + d$nB else rep(0, n)
  ns = as.data.table(states)[state != "neutral"]
  if (nrow(ns) == 0) stop("score_cells: no non-neutral state")
  sn = unique(ns$state)
  sdp = pmax(fifelse(is.na(d$sd_null), 1, d$sd_null), params$sd_floor)
  tp = NULL
  if (depth_model == "t") {
    if (is.null(tpar)) stop("score_cells: depth_model 't' needs `tpar` (fit_depth_null)")
    tp = match_depth_null(d, tpar)
  }
  fk = function(x) sprintf("%.17g", x)
  dcache = list()
  for (s in unique(ns$shift)) dcache[[fk(s)]] = if (s == 0) rep(0, n) else ll_depth_shift(d$z, s / sdp, tp)
  acache = list()
  if (length(rho) != 1 && length(rho) != n) stop("score_cells: `rho` must be one value or one per row")
  rgrp = if (length(rho) > 1 && length(unique(rho)) > 1) split(seq_len(n), rho) else NULL   # per-row rho: score by group
  if (length(rho) > 1 && is.null(rgrp)) rho = rho[1]
  for (f in unique(ns$f)) {
    if (is.null(rgrp)) a = ll_allelic(nA, nh, f, rho, eps, two_sided = two) else {
      a = numeric(n)
      for (g in rgrp) a[g] = ll_allelic(nA[g], nh[g], f, rho[g[1]], eps, two_sided = two)
    }
    acache[[fk(f)]] = fifelse(is.na(a), 0, a)
  }
  ar = lapply(seq_len(nrow(ns)), function(i) acache[[fk(ns$f[i])]])
  names(ar) = ns$row_id
  cl = list()
  for (s in sn) {
    r = which(ns$state == s)
    if (length(unique(ns$shift[r])) == 1) {
      dl = dcache[[fk(ns$shift[r[1]])]]
      al = if (length(r) == 1) ar[[r]] else do.call(pmax, unname(ar[r]))
      cs = dl + al
    } else {
      comb = lapply(r, function(i) dcache[[fk(ns$shift[i])]] + ar[[i]])
      best = max.col(do.call(cbind, comb), ties.method = "first")
      D = do.call(cbind, lapply(r, function(i) dcache[[fk(ns$shift[i])]]))
      A = do.call(cbind, ar[r])
      idx = cbind(seq_len(n), best)
      dl = D[idx]; al = A[idx]; cs = do.call(pmax, comb)
    }
    set(d, j = paste0("dl_", s), value = dl)
    set(d, j = paste0("al_", s), value = al)
    set(d, j = paste0("cl_", s), value = cs)
    cl[[s]] = cs
  }
  # depth-only rows (allelic_ok FALSE: male chrX / chrY): a purely allelic state (no depth shift, e.g. CN-LOH) has
  # no evidence there; its LLR would be 0 and could win the argmax over negative depth states, so it is excluded
  noal = if ("allelic_ok" %in% names(d)) which(d$allelic_ok %in% FALSE) else integer()
  if (length(noal)) for (s in pure_allelic_states(ns)) {
    cl[[s]][noal] = -Inf
    set(d, i = noal, j = paste0("cl_", s), value = -Inf)
  }
  M = do.call(cbind, cl)
  set(d, j = "s_comb", value = if (length(cl) == 1) cl[[1]] else do.call(pmax, unname(cl)))
  set(d, j = "state_comb", value = sn[max.col(M, ties.method = "first")])
  zz = fifelse(is.na(d$z), 0, d$z)
  set(d, j = "s_depth", value = abs(zz))
  set(d, j = "state_depth", value = fifelse(zz < 0, "loss", "gain"))
  set(d, j = "s_allelic", value = if (length(ar) == 1) ar[[1]] else do.call(pmax, unname(ar)))
  ac = allelic_classes(ns)
  A = do.call(cbind, lapply(ac$states, function(ss) {
    r = which(ns$state %in% ss)
    if (length(r) == 1) ar[[r]] else do.call(pmax, unname(ar[r]))
  }))
  set(d, j = "state_allelic", value = unname(ac$label)[max.col(A, ties.method = "first")])
  if (keep_rows) {
    for (i in seq_len(nrow(ns))) set(d, j = paste0("ar_", ns$row_id[i]), value = ar[[i]])
    for (s in sn) {
      r = which(ns$state == s)
      if (length(unique(ns$shift[r])) > 1) for (i in r) set(d, j = paste0("dr_", ns$row_id[i]), value = dcache[[fk(ns$shift[i])]])
    }
  }
  d[]
}

# ---- folds and thresholds -------------------------------------------------------------------------

#' Cross-validation folds of cells
#'
#' Random, balanced folds over cells (or over their origins: spiked copies
#' inherit the fold of the cell they were made from, so a copy and its
#' original are never in different folds).
#' @param cells cell ids (may repeat, e.g. one per row)
#' @param n_folds number of folds
#' @param seed RNG seed
#' @param origin optional origin id per entry of `cells`
#' @return integer fold (1..n_folds) per entry of `cells`
#' @export
assign_folds = function(cells, n_folds = 4, seed = 1, origin = NULL) {
  base = if (is.null(origin)) as.character(cells) else fifelse(is.na(origin), as.character(cells), as.character(origin))
  u = unique(base)
  set.seed(seed)
  fu = sample(rep_len(seq_len(n_folds), length(u)))
  as.integer(fu[match(base, u)])
}

# null rows: `null` if given, else `is_null`, else null_cell & null_region & eligible
null_flag = function(d, null = NULL) {
  if (!is.null(null)) {
    if (length(null) != nrow(d)) stop("null flag must have one value per row")
    return(!is.na(null) & as.logical(null))
  }
  if ("is_null" %in% names(d)) return(!is.na(d$is_null) & d$is_null)
  if (!"null_cell" %in% names(d)) stop("no null definition: need `is_null` or `null_cell` (+ null_region, eligible)")
  v = d$null_cell
  if ("null_region" %in% names(d)) v = v & d$null_region
  if ("eligible" %in% names(d)) v = v & d$eligible
  !is.na(v) & v
}

q7 = function(x, q) {
  x = x[is.finite(x)]
  if (length(x) == 0) NA_real_ else unname(quantile(x, q, type = 7))
}

group_key = function(d, by) {
  if (length(by) == 0) return(rep("all", nrow(d)))
  do.call(paste, c(lapply(by, function(b) as.character(d[[b]])), sep = "\r"))
}

# The null rows a group's threshold is trained on ([budget_thresholds()]; also the calibration rows of [event_fdr()]).
# Levels: 1 = the `by` group, then `by` without the vote stratum (`hstratum`), then without the depth stratum
# (`stratum`). A group with fewer than `min_null` null rows at level 1 takes the first coarser level with `min_null`
# rows, else the coarsest (by = class x stratum: the class). `isn` = null flag per row of `d`. A class without null
# rows whose rows have an `sd_class` with null rows (named events; when `d` has `sd_class`) borrows the null rows
# (and the fallback levels) of that sd_class: its groups carry `borrowed_from`.
# Returns list(levels = the `by` columns per level, groups = unique `by` rows of d (+ `borrowed_from` when a class
#              borrows) with key_ (their level-1 key; that of the sd_class for borrowed groups), level = the level
#              used per group, use = per group the indices (rows of d) of its null rows, borrowed = a class borrows)
threshold_null_rows = function(d, by, isn, min_null) {
  levs = list(by); cur = by
  for (drop in c("hstratum", "stratum")) if (drop %in% cur) { cur = setdiff(cur, drop); levs[[length(levs) + 1]] = cur }
  # classes without null rows (named events) borrow the null rows of their sd_class
  src = NULL
  if ("class" %in% by && "sd_class" %in% names(d)) {
    cls = as.character(d$class); sdc = as.character(d$sd_class)
    has = unique(cls[isn])
    b = !cls %in% has & !is.na(sdc) & sdc != cls & sdc %in% has
    if (any(b)) src = fifelse(b, sdc, NA_character_)
  }
  groups = if (is.null(src)) unique(d[, by, with = FALSE]) else unique(cbind(d[, by, with = FALSE], borrowed_from = src))
  gsrc = groups
  if (!is.null(src)) gsrc = copy(groups)[!is.na(borrowed_from), class := borrowed_from]
  groups[, key_ := group_key(gsrc, by)]
  ni = which(isn)
  idx_l = lapply(levs, function(l) split(ni, group_key(d, l)[ni]))
  key_l = lapply(levs, function(l) group_key(gsrc, l))
  level = integer(nrow(groups)); use = vector("list", nrow(groups))
  for (g in seq_len(nrow(groups))) {
    l = 1L
    u = idx_l[[1]][[key_l[[1]][g]]] %||% integer()
    while (length(u) < min_null && l < length(levs)) {
      l = l + 1L
      u = idx_l[[l]][[key_l[[l]][g]]] %||% integer()
    }
    level[g] = l; use[[g]] = u
  }
  list(levels = levs, groups = groups, level = level, use = use, borrowed = !is.null(src))
}

#' Null thresholds at the false-call budget
#'
#' Type-7 quantile `q` of each score on the null rows ([budget_threshold()]
#' compatible), per `by` group (class x stratum). With `cv = TRUE` one
#' threshold per fold k, trained on the null rows of all other folds (applied
#' to fold k), plus fold 0 trained on all null rows. A group with fewer than
#' `min_null` null rows uses the group pooled over strata (the class-wide
#' threshold), flagged `pooled_fallback`.
#'
#' With a vote stratum (`hstratum` column, [assign_hstratum()]; `calib = "votes"`) the groups are class x
#' stratum x hstratum, and a group short of `min_null` null rows falls back to its class x stratum, then to
#' its class (the first with `min_null` rows, else the class).
#'
#' A class without any null row whose rows have an `sd_class` with null rows (named `event` regions: never null
#' rows) takes the thresholds of that `sd_class`: its groups are keyed by class x sd_class x the other `by`
#' columns and use the null rows (and the fallback levels) of the `sd_class`, so every fold / stratum gets the
#' same threshold as the `sd_class` rows of that group. Those threshold rows carry `borrowed_from` (the
#' sd_class; the column is only added when a class borrows) and [apply_thresholds()] matches them on the rows'
#' `sd_class`.
#' @param d scored rows with the `by` columns, `fold` (cv) and a null definition
#' @param scores score columns
#' @param by grouping columns (default class x stratum, plus `hstratum` when `d` has that column)
#' @param q null quantile
#' @param cv per-fold thresholds
#' @param min_null minimal null rows of a group before falling back
#' @param params [caller_defaults()]
#' @param null optional logical null flag per row (default `is_null`, else
#'   `null_cell & null_region & eligible`)
#' @return thresholds table (check_thresholds schema; `stratum` = "all" when not in `by`; `borrowed_from` when a
#'   class borrows its sd_class's thresholds)
#' @export
budget_thresholds = function(d, scores = c("s_depth", "s_allelic", "s_comb"), by = c("class", "stratum", if ("hstratum" %in% names(d)) "hstratum"),
                             q = params$null_q, cv = TRUE, min_null = params$min_null_rows,
                             params = caller_defaults(), null = NULL) {
  d = as.data.table(d)
  isn = null_flag(d, null)
  fold = if (cv && "fold" %in% names(d)) as.integer(d$fold) else rep(0L, nrow(d))
  folds = setdiff(sort(unique(fold[!is.na(fold)])), 0L)
  tg = threshold_null_rows(d, by, isn, min_null)
  groups = tg$groups
  out = list()
  for (g in seq_len(nrow(groups))) {
    use = tg$use[[g]]
    fb = tg$level[g] > 1L
    fu = fold[use]
    for (fo in c(0L, folds)) {
      tr = if (fo == 0L) rep(TRUE, length(use)) else !is.na(fu) & fu != fo
      for (s in scores) {
        out[[length(out) + 1]] = list(g = g, score = s, fold = fo, thr = q7(d[[s]][use[tr]], q), n_null = sum(tr), pooled_fallback = fb)
      }
    }
  }
  res = rbindlist(out)
  res = cbind(groups[res$g, c(by, if (tg$borrowed) "borrowed_from"), with = FALSE], res[, !"g"])
  if (!"stratum" %in% by) res[, stratum := "all"]
  if (!"class" %in% by) res[, class := "all"]
  res[, `:=`(q = q, n_null = as.integer(n_null), fold = as.integer(fold))]
  for (b in intersect(c("class", "stratum", "hstratum"), names(res))) set(res, j = b, value = as.character(res[[b]]))
  setcolorder(res, c("score", "class", "stratum", "fold", "q", "thr", "n_null", "pooled_fallback"))
  res[]
}

#' Apply null thresholds to scored rows
#'
#' Adds `thr_<score>` and `call_<score>` (`score > thr`, strict). A row
#' whose `fold` has a threshold in the table gets that fold's (cross-validated)
#' threshold; other rows (no `fold` column, fold 0, or a fold the table does
#' not have) get fold 0. A thresholds table whose `stratum` is all "all"
#' (class-level, e.g. pooled) is matched on class only. A table with a non-NA `hstratum` column (vote
#' calibration) is also matched on the rows' `hstratum`. Threshold rows with a `borrowed_from` sd_class (a class
#' without null rows, [budget_thresholds()]) apply to the rows of that class whose `sd_class` is `borrowed_from`.
#' @param d scored rows
#' @param thr [budget_thresholds()] / [pooled_thresholds()] output
#' @param scores score columns to threshold
#' @return copy of `d` with the threshold and call columns
#' @export
apply_thresholds = function(d, thr, scores = unique(thr$score)) {
  d = copy(as.data.table(d))
  thr = as.data.table(thr)
  fold = if ("fold" %in% names(d)) as.integer(d$fold) else rep(0L, nrow(d))
  fold[is.na(fold)] = 0L
  for (s in scores) {
    tt = thr[score == s]
    if (nrow(tt) == 0) stop("apply_thresholds: no thresholds for score ", s)
    by = c(if (!all(tt$class == "all")) "class", if (!all(tt$stratum == "all")) "stratum",
           if ("hstratum" %in% names(tt) && !all(is.na(tt$hstratum))) "hstratum",
           if ("sexgrp" %in% names(tt) && !all(is.na(tt$sexgrp))) "sexgrp")
    if ("sexgrp" %in% by && !"sexgrp" %in% names(d)) stop("apply_thresholds: thresholds are keyed on sexgrp but the rows have no sexgrp column")
    if ("hstratum" %in% by && !"hstratum" %in% names(d)) stop("apply_thresholds: thresholds are keyed on hstratum but the rows have no hstratum column")
    kd = group_key(d, by); kt = group_key(tt, by)
    if ("class" %in% by && "borrowed_from" %in% names(tt) && any(!is.na(tt$borrowed_from))) {
      # borrowed thresholds: keyed on class x sd_class (rows of other classes keep an empty prefix)
      bp = unique(paste(tt$class, tt$borrowed_from, sep = "\r")[!is.na(tt$borrowed_from)])
      sdc = if ("sd_class" %in% names(d)) as.character(d$sd_class) else rep(NA_character_, nrow(d))
      kd = paste(fifelse(paste(d$class, sdc, sep = "\r") %in% bp, sdc, ""), kd, sep = "\n")
      kt = paste(fifelse(is.na(tt$borrowed_from), "", as.character(tt$borrowed_from)), kt, sep = "\n")
    }
    k1 = match(paste(kd, fold, sep = "\t"), paste(kt, tt$fold, sep = "\t"))
    k0 = match(paste(kd, 0L, sep = "\t"), paste(kt, tt$fold, sep = "\t"))
    k = fifelse(fold != 0L & !is.na(k1), k1, k0)
    th = tt$thr[k]
    set(d, j = paste0("thr_", s), value = th)
    set(d, j = paste0("call_", s), value = d[[s]] > th)
  }
  d[]
}

#' Realised false-call rates on the null rows
#'
#' @param d rows with `call_<score>` columns ([apply_thresholds()]) and a null definition
#' @param scores scores (without the `call_` prefix)
#' @param by grouping of the per-1,000 rates
#' @param null optional null flag per row
#' @return list: `rates` (by, score, n_null, n_called, per_1000) and
#'   `per_cell` (score, n_cells, events_per_cell = called autosomal arm null
#'   rows per null cell averaged over null cells (the D1' analogue),
#'   frac_cells = fraction of null cells with at least one)
#' @export
realised_false_calls = function(d, scores = "s_comb", by = c("class", "stratum"), null = NULL) {
  d = as.data.table(d)
  isn = null_flag(d, null)
  x = d[isn]
  chrom = if ("chrom" %in% names(x)) x$chrom else sub("[pq]$", "", sub(":.*$", "", sub("^[^|]*\\|", "", x$region)))
  arm_rows = x$class == "arm" & chrom %in% paste0("chr", 1:22)
  rates = list(); pc = list()
  for (s in scores) {
    cc = x[[paste0("call_", s)]]
    r = cbind(x[, by, with = FALSE], cc_ = cc)[, .(n_null = sum(!is.na(cc_)), n_called = sum(cc_, na.rm = TRUE)), by = by]
    r[, `:=`(score = s, per_1000 = 1000 * n_called / n_null)]
    rates[[s]] = r
    a = data.table(cell = x$cell[arm_rows], k = cc[arm_rows])[, .(k = sum(k, na.rm = TRUE)), by = cell]
    pc[[s]] = data.table(score = s, n_cells = nrow(a), events_per_cell = if (nrow(a)) mean(a$k) else NA_real_,
                         frac_cells = if (nrow(a)) mean(a$k > 0) else NA_real_)
  }
  rates = rbindlist(rates)
  setcolorder(rates, c(by, "score"))
  list(rates = rates[], per_cell = rbindlist(pc))
}

# ---- vote strata and per-chromosome allelic rho ---------------------------------------------------

#' Vote strata: bins of the allele votes per class x depth stratum
#'
#' Under the null the allelic LLR of a fixed gain fraction f is linear in nA: its mean is about
#' `n_hap x 0.5 log(4 f (1 - f))` and its spread grows with sqrt(n_hap), so one null quantile per class is set
#' by the vote-rich regions (large chromosomes) and is too loose for vote-poor ones (small chromosomes such as
#' chr21 carry several times fewer votes than chr1). Thresholds keyed on class x stratum x `hstratum` (bins of n_hap) pool regions with similar
#' vote counts instead of estimating one threshold per region.
#'
#' Per class x stratum, the null rows with votes (`n_hap > 0`) are cut at their type-7 quantiles into
#' `k = max(1, min(n_bins, floor(n / min_null)))` bins of (nearly) equal size, so a bin keeps `min_null`
#' null rows. Bin i covers `lo < n_hap <= hi`. [assign_hstratum()] labels rows `h1`..`hk`, `none` (no votes)
#' or `all` (a class x stratum without bins).
#' @param d rows with `class`, `stratum`, `n_hap` and a null definition ([budget_thresholds()])
#' @param n_bins maximal bins per class x stratum (`params$vote_bins`)
#' @param min_null minimal null rows per bin (`params$min_null_rows`)
#' @param null optional logical null flag per row
#' @param params [caller_defaults()]
#' @return data.table class, stratum, hstratum, lo, hi, n_null (null rows with votes in the class x stratum)
#' @export
vote_bins = function(d, n_bins = params$vote_bins, min_null = params$min_null_rows, null = NULL, params = caller_defaults()) {
  d = as.data.table(d)
  miss = setdiff(c("class", "stratum", "n_hap"), names(d))
  if (length(miss)) stop("vote_bins: d lacks ", paste(miss, collapse = ", "))
  isn = null_flag(d, null)
  x = data.table(class = as.character(d$class[isn]), stratum = as.character(d$stratum[isn]), n = as.numeric(d$n_hap[isn]))[is.finite(n) & n > 0]
  if (!nrow(x)) return(data.table(class = character(), stratum = character(), hstratum = character(), lo = numeric(), hi = numeric(), n_null = integer()))
  out = x[, {
    k = max(1L, min(as.integer(n_bins), as.integer(floor(.N / max(min_null, 1)))))
    br = if (k > 1) unique(stats::quantile(n, seq_len(k - 1) / k, type = 7, names = FALSE)) else numeric()
    .(hstratum = paste0("h", seq_len(length(br) + 1L)), lo = c(0, br), hi = c(br, Inf), n_null = .N)
  }, by = .(class, stratum)]
  out[, n_null := as.integer(n_null)]
  out[]
}

#' Vote stratum of each row
#'
#' @param d rows with `class`, `stratum`, `n_hap`
#' @param bins [vote_bins()] output
#' @return character vector: `h<i>` (bin of the row's class x stratum), `none` (n_hap 0 or NA) or `all`
#'   (class x stratum without bins). Rows of a class without bins (named `event` regions) use the bins of their
#'   `sd_class` when `d` has that column (they borrow its thresholds, [budget_thresholds()]).
#' @export
assign_hstratum = function(d, bins) {
  bins = as.data.table(bins)
  n = as.numeric(d$n_hap)
  none = is.na(n) | n <= 0
  h = rep("all", length(n))
  cls = as.character(d$class)
  if ("sd_class" %in% names(d)) {
    b = !cls %in% bins$class & as.character(d$sd_class) %in% bins$class
    cls[b] = as.character(d$sd_class)[b]
  }
  key = paste(cls, d$stratum, sep = "\r")
  bk = paste(bins$class, bins$stratum, sep = "\r")
  for (k in unique(bk)) {
    i = which(key == k & !none)
    if (!length(i)) next
    hi = bins$hi[bk == k]
    br = hi[is.finite(hi)]
    h[i] = paste0("h", findInterval(n[i], br, left.open = TRUE) + 1L)
  }
  h[none] = "none"
  h
}
utils::globalVariables(c("hstratum", "lo", "hi"))

#' Per-chromosome allelic overdispersion, shrunk to the global
#'
#' The beta-binomial rho of the null allele votes can differ between chromosomes (small chromosomes up to
#' about twice the genome-wide value), so one global rho under-penalises the allelic tails of the more overdispersed ones.
#' Per chromosome, rho is fitted ([fit_allelic_null()]) on the null rows of `classes` (default the arms, as the
#' global rho) and shrunk to the global: `rho_c = (n x rho_raw + k x rho_global) / (n + k)`, n = null rows with
#' votes. Chromosomes without null rows (chrX in normal mode) keep the global rho ([rho_rows()]).
#' @param d rows with `nA`, `n_hap`, `class`, `chrom` (or a parsable `region`) and a null definition
#' @param rho global rho (NULL: fitted on the same rows)
#' @param k shrinkage pseudo-count (null rows; `params$rho_chrom_k`)
#' @param null optional logical null flag per row
#' @param classes region classes of the fit
#' @param params [caller_defaults()]
#' @return data.table chrom, n_rows, rho_raw, rho_global, w (= n / (n + k)), rho
#' @export
rho_by_chrom = function(d, rho = NULL, k = params$rho_chrom_k, null = NULL, classes = "arm", params = caller_defaults()) {
  d = as.data.table(d)
  isn = null_flag(d, null)
  ch = if ("chrom" %in% names(d)) d$chrom else region_chrom(d$region, d$class)
  ok = isn & d$class %in% classes & !is.na(d$n_hap) & d$n_hap > 0
  if (!any(ok)) ok = isn & !is.na(d$n_hap) & d$n_hap > 0
  if (!any(ok)) stop("rho_by_chrom: no null rows with votes")
  if (is.null(rho)) rho = fit_allelic_null(d$nA[ok], d$n_hap[ok])
  f = fit_allelic_null(d$nA[ok], d$n_hap[ok], by = ch[ok])
  setnames(f, c("group", "rho"), c("chrom", "rho_raw"))
  f[, `:=`(rho_global = rho, w = n_rows / (n_rows + k), rho = (n_rows * rho_raw + k * rho) / (n_rows + k))]
  setcolorder(f, c("chrom", "n_rows", "rho_raw", "rho_global", "w", "rho"))
  f[order(match(chrom, paste0("chr", c(1:22, "X", "Y"))))]
}
utils::globalVariables(c("rho_raw", "rho_global"))

#' Per-row rho from a [rho_by_chrom()] table (the global rho where the chromosome has none)
#' @param d rows with `chrom` (or a parsable `region`)
#' @param tab [rho_by_chrom()] output
#' @param rho global rho
#' @return numeric vector, one per row
#' @export
rho_rows = function(d, tab, rho) {
  ch = if ("chrom" %in% names(d)) d$chrom else region_chrom(d$region, d$class)
  r = tab$rho[match(ch, tab$chrom)]
  r[is.na(r)] = rho
  r
}

# pooled threshold group of a chromosome: X, Y or A (autosome)
sex_group = function(chrom) fifelse(chrom %in% "chrX", "X", fifelse(chrom %in% "chrY", "Y", "A"))
utils::globalVariables("sexgrp")

#' Allelically overdispersed tiles (data-driven allelic mask)
#'
#' Per tile of class `size`, the method-of-moments beta-binomial overdispersion of the null rows' allele votes,
#' `rho = (sum (nA - n/2)^2 - sum n/4) / sum n (n - 1) / 4`, then `z_rho` = (log rho - median) / MAD of log rho
#' over the tiles with at least `min_votes` votes per null row. Tiles with `z_rho > k_mad` are masked: a local
#' stretch whose votes are far more overdispersed than the genome (e.g. a 5 Mb tile with rho 0.040 vs a
#' genome median of 0.011 over 5 Mb tiles) and that inflates the chromosome-level allelic LLRs of the cells. The votes of the
#' atoms in masked tiles ([atoms_in_tiles()]) are then left out of the allele features
#' ([update_allelic_features()], `call_scn(params$mask_allelic = TRUE)`).
#'
#' Shared allelic bias (`params$mask_f0`): per tile also the vote-weighted haplotype fraction of the null rows
#' `f0 = sum nA / sum n_hap`, its deviation `d0 = |f0 - 0.5|` and `z_f0` = (d0 - median) / MAD of d0 over the
#' assessed tiles. A tile is a bias outlier (`masked_f0`) when it is assessed, `z_f0 > f0_mad` and
#' `d0 > f0_min`: a stretch whose allele fraction is off 1:1 in the null cells for reasons that are not somatic
#' copy number (germline CNV or segmental duplication, imprinted or allele-specific accessible regions, mapping
#' bias, false heterozygous SNPs), against which the 1:1 allele null of every cell is wrong. The rule sees only
#' the null rows, but it cannot tell such a shared bias from a real clonal event whose carriers remain among the
#' null cells: an event frequent enough for that looks the same and can be masked (in a simulation with a focal
#' event in 35 % of cells, 60 of 96 carriers stayed in the null set and all four event tiles were masked). The
#' mask is a diagnostic to be read with the list of masked tiles, not a filter to apply blindly to samples that
#' may carry clonal events.
#' `|f0 - 0.5|` covers two orientations of the allele features: under `loo` (block) orientation the fraction is
#' major-oriented (a biased tile sits above 0.5), under chromosome-anchored orientation it sits at f0 or 1 - f0.
#' Under `mirrored` orientation each block counts its major allele, so f0 is above 0.5 in every tile by an amount
#' set by the votes per block: the statistic measures coverage, not shared bias, and `call_scn()` does not apply
#' the bias rule (it says so in the run notes and warns).
#' [allelic_bias_tiles()] gives the same table over all cells before calling.
#'
#' `rules` selects which flags make up `masked` (the tiles whose atoms are dropped): `"rho"` (the default and
#' the behaviour of earlier versions: `masked` = `masked_rho`), `"f0"`, or both (union).
#' @param feat_long features long table (class, region, chrom, nA, n_hap) with a null definition
#' @param null optional logical null flag per row
#' @param size tile class
#' @param k_mad threshold in MADs of log rho (`params$mask_rho_mad`)
#' @param min_votes minimal mean votes per null row of a tile for it to be assessed (`params$mask_min_votes`)
#' @param params [caller_defaults()]
#' @param f0_mad bias threshold in MADs of d0 (`params$mask_f0_mad`)
#' @param f0_min minimal deviation d0 of a bias outlier (`params$mask_f0_min`)
#' @param rules the rules that define `masked`: any of `"rho"`, `"f0"`
#' @return data.table tile, chrom, start, end, n_rows, votes, rho, z_rho, assessed, masked, f0, d0, z_f0,
#'   masked_rho, masked_f0 (one row per tile)
#' @export
allelic_noisy_tiles = function(feat_long, null = NULL, size = "5Mb", k_mad = params$mask_rho_mad, min_votes = params$mask_min_votes,
                               params = caller_defaults(), f0_mad = params$mask_f0_mad, f0_min = params$mask_f0_min,
                               rules = "rho") {
  rules = unique(as.character(rules))
  if (!length(rules) || !all(rules %in% c("rho", "f0"))) stop("allelic_noisy_tiles: `rules` must be \"rho\", \"f0\" or both")
  if ("f0" %in% rules && (is.null(f0_mad) || is.null(f0_min))) stop("allelic_noisy_tiles: the f0 rule needs f0_mad and f0_min (params$mask_f0_mad, mask_f0_min)")
  if (is.null(f0_mad)) f0_mad = 3
  if (is.null(f0_min)) f0_min = 0.03
  d = as.data.table(feat_long)
  isn = null_flag(d, null)
  ch = if ("chrom" %in% names(d)) d$chrom else region_chrom(d$region, d$class)
  i = which(isn & d$class == size)
  x = data.table(tile = d$region[i], chrom = ch[i], a = as.numeric(d$nA[i]), n = as.numeric(d$n_hap[i]))
  x[is.na(n), `:=`(n = 0, a = 0)]
  t = x[, .(n_rows = .N, votes = sum(n), S = sum((a - n / 2)^2), N4 = sum(n / 4), P = sum(n * (n - 1) / 4), A = sum(a)), by = .(tile, chrom)]
  t[, rho := fifelse(P > 0, pmax((S - N4) / P, 1e-6), NA_real_)]
  t[, assessed := !is.na(rho) & votes / n_rows >= min_votes]
  lr = log(t$rho[t$assessed])
  m = if (length(lr)) stats::median(lr) else NA_real_
  s = if (length(lr) > 1) stats::mad(lr) else NA_real_
  t[, z_rho := (log(rho) - m) / s]
  t[, masked_rho := assessed & is.finite(z_rho) & z_rho > k_mad]
  # shared allelic bias: vote-weighted null fraction, robust z of its deviation from 1:1 over the assessed tiles
  t[, f0 := fifelse(votes > 0, A / votes, NA_real_)]
  t[, d0 := abs(f0 - 0.5)]
  dv = t$d0[t$assessed]
  m0 = if (length(dv)) stats::median(dv) else NA_real_
  s0 = if (length(dv) > 1) stats::mad(dv) else NA_real_
  t[, z_f0 := (d0 - m0) / s0]
  t[, masked_f0 := assessed & !is.na(z_f0) & z_f0 > f0_mad & d0 > f0_min]
  t[, masked := ("rho" %in% rules & masked_rho) | ("f0" %in% rules & masked_f0)]
  t[, start := suppressWarnings(as.numeric(sub("^.*:", "", tile)))]
  t[, end := start + unname(params$tile_sizes[size])]
  t = t[, .(tile, chrom, start, end, n_rows, votes, rho, z_rho, assessed, masked, f0, d0, z_f0, masked_rho, masked_f0)]
  t[order(match(chrom, paste0("chr", c(1:22, "X", "Y"))), start)]
}
utils::globalVariables(c("assessed", "masked", "z_rho", "N4", "P", "tile", "votes", "A", "f0", "d0", "z_f0", "masked_rho", "masked_f0"))

#' Shared allelic bias per tile over all cells (diagnostic)
#'
#' The per-tile table of [allelic_noisy_tiles()] computed over every cell (not only null cells) and every
#' eligible row of a features object, before any caller run: the vote-weighted haplotype fraction `f0`, its
#' deviation `d0 = |f0 - 0.5|`, the robust z `z_f0` and the bias flag `masked_f0` (here also `masked`), next to
#' the overdispersion columns (`rho`, `z_rho`, `masked_rho`). Over all cells a frequent clonal event also shows up
#' as a bias; the caller's mask (`params$mask_f0`, `call_scn()`) uses the null cells only, which does not keep
#' an event out of the mask when its carriers remain among the null cells ([allelic_noisy_tiles()]).
#' @param features `scn_features` ([cell_features()]) or its features long table
#' @param size tile class
#' @param params [caller_defaults()] (`mask_min_votes`, `mask_rho_mad`, `mask_f0_mad`, `mask_f0_min`)
#' @return data.table as [allelic_noisy_tiles()], one row per tile, `masked` = `masked_f0`
#' @export
allelic_bias_tiles = function(features, size = "5Mb", params = caller_defaults()) {
  d = if (is.data.frame(features)) as.data.table(features) else as.data.table(features$long)
  if (is.null(d) || !nrow(d)) stop("allelic_bias_tiles: `features` must be an scn_features object or a features long table")
  use = if ("eligible" %in% names(d)) d$eligible %in% TRUE else rep(TRUE, nrow(d))
  allelic_noisy_tiles(d, null = use, size = size, params = params, rules = "f0")
}

#' Atoms whose start lies in the masked tiles
#' @param atoms atoms table (atom, chrom, start, end)
#' @param tiles [allelic_noisy_tiles()] output (rows with `masked`), or any table with chrom, start, end
#' @return integer atom ids
#' @export
atoms_in_tiles = function(atoms, tiles) {
  atoms = as.data.table(atoms); tiles = as.data.table(tiles)
  if ("masked" %in% names(tiles)) tiles = tiles[masked %in% TRUE]
  if (!nrow(tiles) || !nrow(atoms)) return(integer())
  hit = rep(FALSE, nrow(atoms))
  for (i in seq_len(nrow(tiles))) hit = hit | (atoms$chrom == tiles$chrom[i] & atoms$start >= tiles$start[i] & atoms$start < tiles$end[i])
  as.integer(atoms$atom[hit])
}

# t fits for every class x stratum of d: missing groups take the fit of the row's sd_class, then the
# class pooled over strata (fit on the pooled null rows is not available here, so the median of the
# class's strata), else stop
complete_tpar = function(d, tpar) {
  tpar = as.data.table(tpar)
  need = unique(d[, intersect(c("class", "sd_class", "stratum"), names(d)), with = FALSE])
  if (!"sd_class" %in% names(need)) need[, sd_class := class]
  miss = need[!paste(class, stratum) %in% paste(tpar$class, tpar$stratum)]
  if (nrow(miss) == 0) return(tpar)
  add = list()
  for (i in seq_len(nrow(miss))) {
    m = miss[i]
    src = tpar[class == m$sd_class & stratum == m$stratum]
    if (nrow(src) == 0) src = tpar[class == m$class][, .(t_m = median(t_m), t_s = median(t_s), t_df = median(t_df), n = sum(n))]
    if (nrow(src) == 0) src = tpar[class == m$sd_class][, .(t_m = median(t_m), t_s = median(t_s), t_df = median(t_df), n = sum(n))]
    if (nrow(src) == 0) stop(sprintf("complete_tpar: no t fit for class %s (sd_class %s)", m$class, m$sd_class))
    add[[i]] = data.table(class = m$class, stratum = m$stratum, t_m = src$t_m[1], t_s = src$t_s[1], t_df = src$t_df[1], n = 0L)
  }
  rbind(tpar, rbindlist(add), fill = TRUE)
}

#' Robust overdispersion of null allele votes
#'
#' Median-based estimate of the beta-binomial rho that tolerates a minority of
#' imbalanced rows (events in a not yet cleaned null set): with
#' `zj = (nA + U - 1/2 - n/2) / sqrt(n/4)` (U uniform jitter, seeded) the
#' variance is `1 + (n - 1) rho + 1/(3n)`, and rho solves
#' `median(|zj| / sqrt(1 + (n - 1) rho + 1/(3n))) = qnorm(0.75)`. The MLE
#' ([fit_allelic_null()]) is used for the final null; this one only seeds the
#' null-cell selection.
#' @param nA votes for haplotype A
#' @param n total votes
#' @param min_n rows with fewer votes are ignored (normal approximation)
#' @param seed jitter seed
#' @return rho in `[1e-4, 0.5]`; the MLE when fewer than 50 rows have `n >= min_n`
#' @export
robust_rho = function(nA, n, min_n = 30, seed = 1) {
  ok = !is.na(n) & !is.na(nA) & n >= min_n
  if (!any(!is.na(n) & n > 0)) return(1e-4)
  if (sum(ok) < 50) return(fit_allelic_null(nA, n))
  set.seed(seed)
  a = nA[ok]; m = n[ok]
  zj = abs(a + runif(length(a)) - 0.5 - m / 2) / sqrt(m / 4)
  f = function(r) median(zj / sqrt(1 + (m - 1) * r + 1 / (3 * m))) - qnorm(0.75)
  if (f(1e-4) <= 0) return(1e-4)
  if (f(0.5) >= 0) return(0.5)
  stats::uniroot(f, c(1e-4, 0.5), tol = 1e-9)$root
}

# Parametric null thresholds: every template row (the current null rows: class, stratum, sd_null, n_hap) is
# re-drawn under the fitted null (z ~ location-scale t of its class x stratum, or N(0, 1); nA ~ BB(n_hap, 1/2,
# rho)), scored with the same states, and the type-7 q quantile of s_comb is taken per class x stratum. Rows are
# replicated until each group has at least n_sim draws.
param_thresholds = function(tmpl, states, tpar, rho, q, params, depth_model = "t", n_sim = NULL, seed = 1,
                            by = c("class", "stratum")) {
  n_sim = n_sim %||% max(2e4, ceiling(20 / (1 - q)))
  tm = as.data.table(tmpl)[, c(by, "sd_null", "n_hap"), with = FALSE]
  tm[, n_hap := fifelse(is.na(n_hap), 0, n_hap)]
  gk = group_key(tm, by)
  reps = ceiling(n_sim / table(gk)[gk])
  sim = tm[rep(seq_len(nrow(tm)), as.integer(reps))]
  set.seed(seed)
  n = nrow(sim)
  if (depth_model == "t") {
    tp = match_depth_null(sim, tpar, by)
    sim[, z := tp$t_m + tp$t_s * stats::rt(n, tp$t_df)]
  } else sim[, z := rnorm(n)]
  p = if (rho > 0) stats::rbeta(n, 0.5 * (1 - rho) / rho, 0.5 * (1 - rho) / rho) else rep(0.5, n)
  sim[, nA := as.numeric(rbinom(n, n_hap, p))]
  sc = score_cells(sim, states, tpar, rho, depth_model = depth_model, params = params)
  res = sc[, .(thr = q7(s_comb, q), n_null = .N), by = by]
  res[, `:=`(score = "s_comb", fold = 0L, q = q, pooled_fallback = FALSE)]
  if (!"stratum" %in% by) res[, stratum := "all"]
  res[, n_null := as.integer(n_null)]
  setcolorder(res, c("score", "class", "stratum", "fold", "q", "thr", "n_null", "pooled_fallback"))
  res[]
}

#' Self-calibrated null cell set (normal tissue)
#'
#' Starts from all cells of `feat_long` minus `exclude` (e.g. multiplet-flagged
#' barcodes). Each iteration fits the depth t null per class x stratum and the
#' allelic rho on the null rows of the current set, scores every candidate
#' cell on the `classes` regions (chromosomes and arms) with `s_comb`,
#' thresholds at the selection quantile `q` and keeps the cells with at most
#' `max_calls` calls. Every iteration re-evaluates all candidates (a cell
#' dropped earlier can come back), so the set does not ratchet.
#'
#' Thresholds (`threshold`): "parametric" (default) re-draws the null rows
#' under the fitted model (t depth, beta-binomial alleles with [robust_rho()])
#' and takes the quantile of the simulated scores; it is not inflated by the
#' events it is meant to find (an empirical 99.9 % quantile of a set in which
#' more than 0.1 % of the rows carry events sits inside the events).
#' "empirical" uses fold-cross-validated quantiles of the current null rows
#' ([budget_thresholds()]), which is only safe when the starting set is nearly
#' clean. The default `q` (1e-4 tail per row) keeps the selection from
#' discarding null cells for their expected 1 / 1,000 false calls (~6 % of
#' cells at 63 chromosome + arm rows); the calling thresholds of
#' [caller_normal()] are empirical at `null_q` on the selected set.
#'
#' Only rows on `chroms` (default the autosomes) are scored: chrX in females
#' carries X-inactivation driven allelic imbalance and male chrX / chrY are
#' depth-only, so a chrX call says little about whether a cell's autosomes are
#' a clean null.
#' @param feat_long features long table (cell, region, class, stratum, z, sd_null, nA, n_hap; fold,
#'   null_region, eligible optional; chrom optional, else parsed from the region name)
#' @param states state set ([caller_states()])
#' @param params [caller_defaults()]
#' @param n_iter iterations
#' @param max_calls cells with more calls than this are not null cells
#' @param exclude cells never in the null set
#' @param classes region classes used for the selection
#' @param depth_model "t" or "gauss"
#' @param threshold "parametric" or "empirical"
#' @param q selection quantile per row
#' @param n_sim simulated rows per class x stratum (parametric; default max(2e4, 20 / (1 - q)))
#' @param chroms chromosomes whose rows are scored (default autosomes chr1-chr22); NULL = all
#' @return list(null_cells, log (per iteration), calls (per candidate cell: n_calls, null), tpar, rho, thresholds)
#' @export
select_null_cells = function(feat_long, states = caller_states("normal", params), params = caller_defaults(), n_iter = 2,
                             max_calls = 0, exclude = NULL, classes = c("chrom", "arm"), depth_model = "t",
                             threshold = c("parametric", "empirical"), q = 0.9999, n_sim = NULL, chroms = paste0("chr", 1:22)) {
  threshold = match.arg(threshold)
  d = as.data.table(feat_long)[class %in% classes & !cell %in% exclude]
  if (!is.null(chroms)) {
    ch = if ("chrom" %in% names(d)) d$chrom else region_chrom(d$region, d$class)
    d = d[ch %in% chroms]
  }
  if (nrow(d) == 0) stop("select_null_cells: no rows of classes ", paste(classes, collapse = ", "), if (!is.null(chroms)) " on the given chromosomes")
  if (!"null_region" %in% names(d)) d[, null_region := TRUE]
  if (!"eligible" %in% names(d)) d[, eligible := TRUE]
  if (!"fold" %in% names(d) || anyNA(d$fold)) d[, fold := assign_folds(cell, params$n_folds, params$seed, if ("origin" %in% names(d)) origin else NULL)]
  cand = unique(d$cell)
  cur = cand
  ok_row = !is.na(d$null_region) & d$null_region & !is.na(d$eligible) & d$eligible
  log = list(); nc = NULL; tpar = NULL; rho = NA_real_; thr = NULL
  for (it in seq_len(n_iter)) {
    isn = ok_row & d$cell %in% cur
    tpar = if (depth_model == "t") complete_tpar(d, fit_depth_null(d[isn])) else NULL
    arm = isn & d$class == "arm"
    if (!any(arm & d$n_hap > 0, na.rm = TRUE)) arm = isn
    rho = if (threshold == "parametric") robust_rho(d$nA[arm], d$n_hap[arm], seed = params$seed) else fit_allelic_null(d$nA[arm], d$n_hap[arm])
    sc = score_cells(d, states, tpar, rho, depth_model = depth_model, params = params)
    thr = if (threshold == "parametric") {
      param_thresholds(d[isn], states, tpar, rho, q, params, depth_model, n_sim, seed = params$seed + it)
    } else budget_thresholds(sc, "s_comb", q = q, cv = TRUE, min_null = params$min_null_rows, params = params, null = isn)
    sc = apply_thresholds(sc, thr, "s_comb")
    nc = sc[, .(n_calls = sum(call_s_comb & eligible, na.rm = TRUE)), by = cell]
    new = nc[n_calls <= max_calls, cell]
    fr = 1000 * mean(sc$call_s_comb[isn], na.rm = TRUE)
    log[[it]] = data.table(iter = it, n_in = length(cur), n_out = length(new), n_dropped = length(setdiff(cur, new)),
                           n_added = length(setdiff(new, cur)), rho = rho, call_rate_per_1000 = fr)
    cur = new
  }
  nc[, null := cell %in% cur]
  list(null_cells = cur, log = rbindlist(log), calls = nc[], tpar = tpar, rho = rho, thresholds = thr)
}

# ---- calling policy, calls, events -------------------------------------------------------------------

#' Calling policy: which class x state is called at which level
#'
#' Normal tissue: per cell, chromosomes in every state and
#' arm / 50 Mb losses and CN-LOH are reportable; 20 Mb losses / CN-LOH and
#' sub-chromosome gains are computed per cell but not reportable
#' (`reportable = FALSE`, recall 0.24-0.74 per cell); 5 / 10 Mb tiles and
#' sub-chromosome gains are pooled (`level = "pool"`, every class x state
#' reportable). Named `event` regions are reportable per cell. Tumour:
#' per-cell rows (detection / clone assignment) are not reportable; the state
#' is reported at clone level.
#'
#' Column `chrx` says whether rows on chrX are reportable at that level. In
#' normal tissue of a female (or unknown sex) chrX pools are not reportable:
#' X-inactivation makes neighbouring cells share the active X, so allelic
#' imbalance of chrX at neighbourhood (pool) scale is expected biology, not a
#' copy-number event; per-cell chrX calls stay reportable. Male chrX is
#' depth-only and its pools are reportable. Pass a policy with `chrx = TRUE`
#' to override.
#' @param mode "normal" or "tumour"
#' @param sex "female", "male" or NULL (unknown, treated as female)
#' @return data.table class, state, level (cell | pool | clone), reportable, chrx
#' @export
call_policy = function(mode = c("normal", "tumour"), sex = NULL) {
  mode = match.arg(mode)
  if (mode == "normal") {
    st = c("loss", "gain", "cnloh")
    cellp = rbind(
      CJ(class = "chrom", state = st, reportable = TRUE),
      CJ(class = c("arm", "50Mb"), state = c("loss", "cnloh"), reportable = TRUE),
      CJ(class = c("arm", "50Mb"), state = "gain", reportable = FALSE),
      CJ(class = "20Mb", state = st, reportable = FALSE),
      CJ(class = "event", state = st, reportable = TRUE))
    cellp[, level := "cell"]
    poolp = CJ(class = c("chrom", "arm", "50Mb", "20Mb", "10Mb", "5Mb", "event"), state = st, reportable = TRUE)
    poolp[, level := "pool"]
    out = rbind(cellp, poolp)
    out[, chrx := level == "cell" | identical(sex, "male")]
  } else {
    st = c("loss", "cnloh", "gain", "gain_as", "amp")
    cl = c("chrom", "arm", "50Mb", "20Mb", "10Mb", "5Mb", "event", "segment")
    out = rbind(CJ(class = cl, state = st, reportable = FALSE)[, level := "cell"],
                CJ(class = cl, state = st, reportable = TRUE)[, level := "clone"])
    out[, chrx := TRUE]
  }
  setcolorder(out, c("class", "state", "level", "reportable", "chrx"))
  out[]
}

# reportable flag of call rows under a policy: the class x state entry, and FALSE on chrX when the policy's `chrx`
# column (if any) says chrX is not reportable at that level
policy_reportable = function(x, pol) {
  k = match(paste(x$class, x$state), paste(pol$class, pol$state))
  r = pol$reportable[k]
  if ("chrx" %in% names(pol)) r = r & (pol$chrx[k] %in% TRUE | is.na(x$chrom) | x$chrom != "chrX")
  r
}

# chromosome of a region name (tiles "5Mb|chr1:0", arms "chr1p", chromosomes "chr1"); NA otherwise
region_chrom = function(region, class) {
  ch = fifelse(grepl("\\|", region), sub(":.*$", "", sub("^[^|]*\\|", "", region)), NA_character_)
  ch = fifelse(class == "arm", sub("[pq]$", "", region), ch)
  fifelse(class == "chrom", region, ch)
}

#' Minimal regions table from a features table
#'
#' For callers given only features: `region`, `class`, `sd_class`, `chrom`
#' and tile coordinates parsed from tile names (`5Mb|chr1:0` -> 0, 5e6);
#' arms / chromosomes / events get NA coordinates.
#' @param feat features long table
#' @param params [caller_defaults()] (tile sizes)
#' @return regions table (check_regions schema)
#' @export
regions_from_features = function(feat, params = caller_defaults()) {
  f = as.data.table(feat)
  if (!"sd_class" %in% names(f)) f = copy(f)[, sd_class := class]
  if (!"chrom" %in% names(f)) f = copy(f)[, chrom := region_chrom(region, class)]
  el = if ("eligible" %in% names(f)) f$eligible else rep(TRUE, nrow(f))
  nr = if ("null_region" %in% names(f)) f$null_region else rep(TRUE, nrow(f))
  r = data.table(region = f$region, class = f$class, sd_class = f$sd_class, chrom = f$chrom, el = el, nr = nr)[
    , .(class = class[1], sd_class = sd_class[1], chrom = chrom[1], eligible = any(el, na.rm = TRUE),
        null_region = all(nr, na.rm = TRUE)), by = region]
  tile = grepl("\\|", r$region) & r$class %in% names(params$tile_sizes)
  r[, start := fifelse(tile, suppressWarnings(as.numeric(sub("^.*:", "", region))), NA_real_)]
  r[, end := fifelse(tile, start + unname(params$tile_sizes[class]), NA_real_)]
  r[, `:=`(n_bins = NA_integer_, n_use = NA_integer_, allelic_ok = TRUE)]
  setcolorder(r, c("region", "class", "sd_class", "chrom", "start", "end", "n_bins", "n_use", "eligible", "null_region", "allelic_ok"))
  r[]
}

add_coords = function(x, regions) {
  if (is.null(regions)) {
    if (!"chrom" %in% names(x)) x[, chrom := region_chrom(region, class)]
    if (!"start" %in% names(x)) x[, start := NA_real_]
    if (!"end" %in% names(x)) x[, end := NA_real_]
    return(x)
  }
  rg = as.data.table(regions)
  k = match(x$region, rg$region)
  x[, `:=`(chrom = fifelse(is.na(k), if ("chrom" %in% names(x)) chrom else NA_character_, rg$chrom[k]),
           start = as.numeric(rg$start[k]), end = as.numeric(rg$end[k]))]
  x
}

#' Per-cell calls table
#'
#' One row per scored cell x region whose class and best state (`state_comb`)
#' the policy lists at `level`; `called = score > thr`. Rows with
#' `eligible = FALSE` (regions that are not scored) are left out.
#' @param d_scored scored rows ([score_cells()]); with `thr_<score>` columns
#'   when `thr_applied` is NULL
#' @param thr_applied a thresholds table to apply ([budget_thresholds()]), or
#'   NULL when `d_scored` already went through [apply_thresholds()]
#' @param regions regions table (coordinates); NULL leaves them NA
#' @param policy [call_policy()]
#' @param score score column
#' @param level policy level
#' @param keep_cols further columns of `d_scored` copied to the calls when present
#' @return calls table (check_calls schema) plus `stratum`, `fold` and `keep_cols` when present
#' @export
call_table = function(d_scored, thr_applied = NULL, regions = NULL, policy = call_policy("normal"), score = "s_comb", level = "cell",
                      keep_cols = c("z", "f", "n_hap")) {
  d = as.data.table(d_scored)
  if (!is.null(thr_applied)) d = apply_thresholds(d, thr_applied, score)
  tc = paste0("thr_", score); sc = sub("^s_", "state_", score)
  if (!tc %in% names(d)) stop("call_table: no ", tc, " column; pass a thresholds table")
  lv = level
  pol = as.data.table(policy)
  pol = pol[pol$level == lv]
  keep = paste(d$class, d[[sc]]) %in% paste(pol$class, pol$state)
  if ("eligible" %in% names(d)) keep = keep & d$eligible %in% TRUE
  extra = intersect(c("chrom", "stratum", "fold", keep_cols), names(d))
  x = data.table(cell = as.character(d$cell[keep]), region = d$region[keep], class = d$class[keep],
                 state = d[[sc]][keep], score = d[[score]][keep], thr = d[[tc]][keep])
  for (e in extra) set(x, j = e, value = d[[e]][keep])
  x = add_coords(x, regions)
  x[, called := !is.na(score) & !is.na(thr) & score > thr]
  x[, `:=`(level = level, source = "cell")]
  x[, reportable := policy_reportable(x, pol)]
  setcolorder(x, c("cell", "region", "class", "chrom", "start", "end", "state", "score", "thr", "called", "level", "source", "reportable"))
  check_calls(x)
  x[]
}

#' Consolidate nested calls into events
#'
#' Per cell (or pool / clone) x chromosome x state, called regions nested in a
#' larger called region of the same state are merged into it (the largest
#' called region, ties to the higher score); overlapping but not nested calls
#' stay separate events.
#' @param calls calls table ([call_table()], [call_pools()])
#' @param regions optional regions table to fill missing coordinates
#' @param reportable_only only reportable calls
#' @return data.table cell, chrom, start, end, state, top_region, score, level, top_class, n_regions
#' @export
consolidate_calls = function(calls, regions = NULL, reportable_only = TRUE) {
  x = as.data.table(calls)[called == TRUE]
  if (reportable_only) x = x[reportable == TRUE]
  empty = data.table(cell = character(), chrom = character(), start = numeric(), end = numeric(), state = character(),
                     top_region = character(), score = numeric(), level = character(), top_class = character(), n_regions = integer())
  if (nrow(x) == 0) return(empty)
  if (!is.null(regions)) {
    rg = as.data.table(regions); k = match(x$region, rg$region)
    fs = which(is.na(x$start) & !is.na(k)); fe = which(is.na(x$end) & !is.na(k))
    set(x, i = fs, j = "start", value = as.numeric(rg$start[k[fs]]))
    set(x, i = fe, j = "end", value = as.numeric(rg$end[k[fe]]))
  }
  # regions without coordinates (arms / chromosomes of a features-only run) nest by name: chrom > arm > tiles
  x[, `:=`(start_ = fifelse(is.na(start), 0, start), end_ = fifelse(is.na(end), fifelse(class == "chrom", Inf, 1e12), end))]
  x[, size_ := end_ - start_]
  setorder(x, level, cell, chrom, state, -size_, -score)
  nest = function(s, e) {
    id = integer(length(s)); es = numeric(); ee = numeric()
    for (i in seq_along(s)) {
      h = which(s[i] >= es & e[i] <= ee)
      if (length(h)) id[i] = h[1] else { es = c(es, s[i]); ee = c(ee, e[i]); id[i] = length(es) }
    }
    id
  }
  x[, ev := nest(start_, end_), by = .(level, cell, chrom, state)]
  out = x[, .(start = start[1], end = end[1], top_region = region[1], score = score[1], top_class = class[1], n_regions = .N),
          by = .(level, cell, chrom, state, ev)]
  out[, ev := NULL]
  setcolorder(out, c("cell", "chrom", "start", "end", "state", "top_region", "score", "level", "top_class", "n_regions"))
  out[]
}
utils::globalVariables("size_")

#' Resolve chromosome-level gains against depth and nested calls
#'
#' A diluted sub-chromosomal imbalance (a 1p loss, a distal CN-LOH) makes the
#' whole chromosome look like a 2:1 gain: at chromosome scale the allelic term
#' of the nested event dominates while its depth contrast is diluted, and the
#' gain state (f 2/3) fits a diluted haplotype fraction better than 1:0 / 2:0.
#' Per cell, a reportable `gain` call of a `classes` region therefore needs
#'   1. depth support, so that the allelic term alone cannot carry a gain label:
#'      * `gate = "dl"` (default): the gain depth LLR of the row itself is
#'        positive (`dl_gain > min_dl_gain`; for the Student-t null, z above
#'        the midpoint between neutral and the gain shift, `t_m + mu / 2` with
#'        `mu = shifts["gain"] / pmax(sd_null, sd_floor)`). For a noisy cell
#'        (large `sd_null`) `mu` is small and the gate is weak;
#'      * `gate = "z"` (opt-in): the row's z exceeds
#'        `max(t_m + mu / 2, t_m + t_s * qt(gate_q, t_df))` with the Student-t
#'        null of the row's class x stratum (`tpar`, [match_depth_null()]; with
#'        `tpar = NULL` the Gaussian null, `max(mu / 2, qnorm(gate_q))`): the
#'        midpoint stays as a floor and the null quantile is added on top.
#'        Over all null-depth rows the gate passes about `1 - gate_q` of them
#'        whatever their `sd_null`; this level is unconditional. Among rows
#'        already called as gains (weak allelic support), and on female chrX
#'        (no null rows of its own), the pass rate of null-depth rows is
#'        higher: up to about two to three times `1 - gate_q` in simulation.
#'        The gate is opt-in and has not been benchmarked on real data as a
#'        default. A row with NA z has no depth support.
#'        `min_dl_gain` is not used. The gate is on z, not on `dl_gain`: the
#'        Student-t LLR falls back towards 0 for very large z
#'        ([ll_depth_shift()]), so a quantile on `dl_gain` would reject the
#'        strongest gains; and
#'   2. (class `chrom`) no nested explanation: no called loss / CN-LOH row of
#'      the same cell and chromosome in `nested` classes (arm, 50 Mb, 20 Mb)
#'      with a higher score.
#' A gain failing either rule stays in the table with `called = TRUE` but
#' `reportable = FALSE` and `resolved` naming the reason (`no_depth`,
#' `nested <region> <state>`). With `promote`, the best (highest-score) nested
#' loss / CN-LOH call that explains a dropped chromosome gain is made
#' reportable (`resolved = "replaces <region> gain"`; arm and 50 Mb losses are
#' reportable anyway, a 20 Mb one is promoted). Losses are not gated: the
#' loss vs CN-LOH label is decided by the scores ([score_cells()]).
#' @param calls cell-level calls table ([call_table()])
#' @param scores scored rows (`cell`, `region`, `dl_gain`; with `gate = "z"` also `z`, `sd_null` and, with
#'   `tpar`, `class`, `stratum`); NULL when `calls` has these columns
#' @param min_dl_gain minimal gain depth LLR of a reportable gain (`gate = "dl"`)
#' @param classes classes whose reportable gains need depth support
#' @param nested classes of the nested calls that can explain a chromosome gain
#' @param promote make the explaining nested call reportable
#' @param gate depth gate of a gain: `"dl"` (gain depth LLR above `min_dl_gain`) or `"z"` (z above the
#'   midpoint and above the `gate_q` quantile of the row's depth null)
#' @param gate_q null quantile of the `"z"` gate (0.95: about 5 % of all null-depth rows pass; more among rows
#'   already called as gains and on female chrX, see Details)
#' @param tpar `"z"` gate: the Student-t fits per class x stratum ([fit_depth_null()], completed as in
#'   [caller_normal()]); NULL = the Gaussian depth null N(0, 1). A gated row whose class x stratum has no fit (a named
#'   `event`: no null rows) takes the fit of its `sd_class` (a `scores` column) in that stratum; without one it is an
#'   error ([match_depth_null()]); a row with NA class or stratum has no depth support
#' @param params [caller_defaults()] (`shifts`, `sd_floor`; the `"z"` gate)
#' @return `calls` with column `resolved` (NA = unchanged) and updated `reportable`; with `gate = "z"` also
#'   `gate_z` (the z the gate used) and `gate_thr` (its z threshold), set on the gated gain rows and NA elsewhere
#' @export
resolve_nested = function(calls, scores = NULL, min_dl_gain = 0, classes = c("chrom", "event"), nested = c("arm", "50Mb", "20Mb"),
                          promote = TRUE, gate = c("dl", "z"), gate_q = 0.95, tpar = NULL, params = caller_defaults()) {
  gate = match.arg(gate)
  if (gate == "z" && (!is.numeric(gate_q) || length(gate_q) != 1 || is.na(gate_q) || gate_q <= 0 || gate_q >= 1))
    stop("resolve_nested: gate_q must be one number in (0, 1)")
  x = copy(as.data.table(calls))
  if (!"resolved" %in% names(x)) x[, resolved := NA_character_]
  if (gate == "z") x[, `:=`(gate_z = NA_real_, gate_thr = NA_real_)]
  lv = if ("level" %in% names(x)) x$level else rep("cell", nrow(x))
  g = which(lv == "cell" & x$class %in% classes & x$state == "gain" & x$called %in% TRUE & x$reportable %in% TRUE)
  if (!length(g)) return(x[])
  key2 = function(a, b) paste(a, b, sep = "|")
  if (gate == "dl") {
    dl = if (!is.null(scores)) {
      sc = as.data.table(scores)
      if (!"dl_gain" %in% names(sc)) stop("resolve_nested: `scores` has no dl_gain column (score_cells)")
      sc$dl_gain[match(key2(x$cell[g], x$region[g]), key2(sc$cell, sc$region))]
    } else if ("dl_gain" %in% names(x)) x$dl_gain[g] else stop("resolve_nested: need `scores` or a dl_gain column")
    no_depth = !is.na(dl) & dl <= min_dl_gain
  } else {
    zg = gain_gate_z(x[g], scores, gate_q, tpar, params)
    no_depth = !(zg$z > zg$thr) %in% TRUE                       # NA z (or no null fit): no depth support
    set(x, i = g, j = "gate_z", value = zg$z)
    set(x, i = g, j = "gate_thr", value = zg$thr)
  }
  # best nested loss / CN-LOH call of the same cell x chromosome
  ni = which(lv == "cell" & x$class %in% nested & x$state %in% c("loss", "cnloh") & x$called %in% TRUE)
  nst = x[ni, .(cell, chrom, region, state, score)][, row_ := ni]
  setorder(nst, cell, chrom, -score)
  best = nst[!duplicated(nst[, .(cell, chrom)])]
  kb = match(key2(x$cell[g], x$chrom[g]), key2(best$cell, best$chrom))
  is_chr = x$class[g] == "chrom"
  nested_win = is_chr & !is.na(kb) & best$score[kb] > x$score[g]
  nested_win[is.na(nested_win)] = FALSE
  drop = no_depth | nested_win
  why = ifelse(nested_win, paste("nested", best$region[kb], best$state[kb]), "no_depth")
  why = ifelse(no_depth & nested_win, paste("no_depth;", why), why)
  if (any(drop)) {
    set(x, i = g[drop], j = "reportable", value = FALSE)
    set(x, i = g[drop], j = "resolved", value = why[drop])
    pr = which(drop & is_chr & !is.na(kb))
    if (promote && length(pr)) {
      rows = best$row_[kb[pr]]
      set(x, i = rows, j = "resolved", value = paste("replaces", x$region[g[pr]], "gain"))
      set(x, i = rows, j = "reportable", value = TRUE)
    }
  }
  x[]
}
utils::globalVariables(c("resolved", "row_", "gate_z", "gate_thr"))

# the "z" gain gate of resolve_nested(): per gain row its z and the z threshold max(t_m + mu / 2, t_m + t_s * q) with
# mu = shift_gain / pmax(sd_null, sd_floor) (sd_null NA -> 1, as in score_cells) and q the gate_q quantile of the row's
# null (Student-t of its class x stratum, or N(0, 1) without tpar). NA z, or NA class / stratum with tpar, gives NA (no
# depth support); NA sd_null counts as 1. A class x stratum that tpar has no fit for (a named event: no null rows) takes
# the fit of the row's sd_class in that stratum, as complete_tpar() does for score_cells; without one it stops
# (match_depth_null)
gain_gate_z = function(rows, scores, gate_q, tpar, params) {
  rows = as.data.table(rows)
  need = c("z", "sd_null", if (!is.null(tpar)) c("class", "stratum"))
  src = if (!is.null(scores)) {
    sc = as.data.table(scores)
    miss = setdiff(need, names(sc))
    if (length(miss)) stop("resolve_nested: gate 'z' needs `scores` column(s) ", paste(miss, collapse = ", "))
    k = match(paste(rows$cell, rows$region, sep = "|"), paste(sc$cell, sc$region, sep = "|"))
    sc[k, intersect(c(need, "sd_class"), names(sc)), with = FALSE]
  } else {
    miss = setdiff(need, names(rows))
    if (length(miss)) stop("resolve_nested: gate 'z' needs `scores` or calls column(s) ", paste(miss, collapse = ", "))
    rows[, intersect(c(need, "sd_class"), names(rows)), with = FALSE]
  }
  z = as.numeric(src$z)
  mu = params$shifts[["gain"]] / pmax(fifelse(is.na(src$sd_null), 1, as.numeric(src$sd_null)), params$sd_floor)
  m = rep(0, nrow(rows)); s = rep(1, nrow(rows)); q = rep(stats::qnorm(gate_q), nrow(rows))
  if (!is.null(tpar)) {
    m[] = NA_real_; s[] = NA_real_; q[] = NA_real_
    ok = which(!is.na(src$class) & !is.na(src$stratum))
    if (length(ok)) {
      key = src[ok, .(class = as.character(class), stratum)]
      if ("sd_class" %in% names(src)) {                          # no fit of the class: the fit of its sd_class
        tk = paste(tpar$class, tpar$stratum, sep = "\r")
        b = !paste(key$class, key$stratum, sep = "\r") %in% tk & !is.na(src$sd_class[ok]) &
          paste(src$sd_class[ok], key$stratum, sep = "\r") %in% tk
        set(key, i = which(b), j = "class", value = as.character(src$sd_class[ok][b]))
      }
      tp = match_depth_null(key, tpar)
      m[ok] = tp$t_m; s[ok] = tp$t_s; q[ok] = stats::qt(gate_q, tp$t_df)
    }
  }
  list(z = z, thr = pmax(m + mu / 2, m + s * q))
}

# ---- pools -----------------------------------------------------------------------------------------

#' Pooled state scores
#'
#' Per pool (x `by` group, e.g. region): the per-state depth LLRs and the
#' per-row allelic LLRs are summed over the pool's rows, then each state takes
#' its best row (depth + alleles); scores as in [score_cells()] (`s_depth` =
#' best pooled depth LLR over the states with a depth shift).
#'
#' Leave-one-out (`loo = TRUE`): for every state row (depth + one allelic
#' row) the largest single-cell contribution is removed from the pooled sum;
#' `s_loo` is the best such row of the called state (`state_comb`) and
#' `s_rob = pmin(s_comb, s_loo)`: a pool passes on `s_rob` only if it still
#' passes after dropping its most influential cell, so one strong carrier
#' cannot make ~20 overlapping kNN pools call. `state_rob` = `state_comb`;
#' `lo_<state>` holds every state's leave-one-out score.
#'
#' Region scale (`scale`, a [pool_scale()] table of null pools of the same
#' construction): the best state is chosen on the LLRs, then `s_comb` and
#' `s_loo` are standardised per region, `(s - med) / mad` of the region's null
#' pools; `s_comb_llr` keeps the raw best LLR.
#' @param d rows scored with `keep_rows = TRUE` ([score_cells()]); with a
#'   `pool` column when `pools` is NULL
#' @param pools NULL or a table `pool`, `cell` (a cell may be in several pools)
#' @param states the state set used for scoring
#' @param by grouping within a pool (default region and class when present)
#' @param loo add the leave-one-out scores `s_loo`, `s_rob`, `state_rob`, `lo_<state>`
#' @param scale NULL or a [pool_scale()] table (per-region null median / MAD of the pooled best LLR)
#' @return data.table pool, by, n, cl_<state>, s_depth, state_depth, s_allelic, state_allelic, s_comb, state_comb
#'   (+ s_loo, s_rob, state_rob, lo_<state>; + s_comb_llr with `scale`)
#' @export
pool_scores = function(d, pools = NULL, states = caller_states("normal"), by = intersect(c("region", "class"), names(d)),
                       loo = FALSE, scale = NULL) {
  if (!is.null(scale) && !all(c("region", "class") %in% by)) stop("pool_scores: `scale` needs `by` to contain region and class")
  d = as.data.table(d)
  ns = as.data.table(states)[state != "neutral"]
  sn = unique(ns$state)
  rt = pool_row_terms(ns)
  multi = rt$multi
  dcols = rt$dcols; acols = rt$acols
  cols = c(dcols, acols)
  miss = setdiff(cols, names(d))
  if (length(miss)) stop("pool_scores: missing ", paste(miss, collapse = ", "), " (score_cells(keep_rows = TRUE) with the same states)")
  rcols = character()
  if (loo) {
    d = d[, unique(c(intersect(c("cell", "pool", "allelic_ok"), names(d)), by, cols)), with = FALSE]
    rcols = paste0("rc_", rt$rows$row_id)
    for (i in seq_len(nrow(rt$rows))) {
      dc = rt$rows$dcol[i]
      set(d, j = rcols[i], value = (if (is.na(dc)) 0 else d[[dc]]) + d[[rt$rows$acol[i]]])
    }
  }
  if (!is.null(pools)) {
    p = unique(as.data.table(pools)[, .(pool, cell)])
    d = d[, c("cell", intersect("allelic_ok", names(d)), by, cols, rcols), with = FALSE][p, on = "cell", allow.cartesian = TRUE, nomatch = NULL]
  } else if (!"pool" %in% names(d)) stop("pool_scores: need a `pool` column or `pools`")
  P = d[, c(lapply(.SD, sum), list(n = .N)), by = c("pool", by), .SDcols = cols]
  # pools of depth-only rows (allelic_ok FALSE everywhere, male chrX / chrY): no purely allelic state
  noal = if ("allelic_ok" %in% names(d) && !"allelic_ok" %in% by) which(!d[, .(a = any(allelic_ok %in% TRUE)), by = c("pool", by)]$a) else integer()
  if (loo) {
    M = d[, lapply(.SD, max), by = c("pool", by), .SDcols = rcols]
    for (r in rcols) set(P, j = r, value = M[[r]])
    rm(M)
  }
  cl = list(); dep = list(); lo = list()
  for (s in sn) {
    r = which(ns$state == s)
    A = lapply(r, function(i) P[[paste0("ar_", ns$row_id[i])]])
    if (!multi[[s]]) {
      sh = ns$shift[r[1]]
      al = if (length(A) == 1) A[[1]] else do.call(pmax, A)
      if (sh == 0) cl[[s]] = al else {
        D = P[[paste0("dl_", s)]]
        cl[[s]] = D + al
        dep[[s]] = D
      }
      if (loo) {
        D0 = if (sh == 0) 0 else P[[paste0("dl_", s)]]
        L = lapply(seq_along(r), function(j) D0 + A[[j]] - P[[paste0("rc_", ns$row_id[r[j]])]])
        lo[[s]] = if (length(L) == 1) L[[1]] else do.call(pmax, L)
      }
    } else {
      Dr = lapply(r, function(i) P[[paste0("dr_", ns$row_id[i])]])
      cl[[s]] = do.call(pmax, Map(`+`, Dr, A))
      dep[[s]] = do.call(pmax, Dr)
      if (loo) lo[[s]] = do.call(pmax, lapply(seq_along(r), function(j) Dr[[j]] + A[[j]] - P[[paste0("rc_", ns$row_id[r[j]])]]))
    }
    if (length(noal) && s %in% pure_allelic_states(ns)) {
      cl[[s]][noal] = -Inf
      if (loo) lo[[s]][noal] = -Inf
    }
    set(P, j = paste0("cl_", s), value = cl[[s]])
  }
  if (length(dep)) {
    set(P, j = "s_depth", value = if (length(dep) == 1) dep[[1]] else do.call(pmax, unname(dep)))
    set(P, j = "state_depth", value = names(dep)[max.col(do.call(cbind, dep), ties.method = "first")])
  } else P[, `:=`(s_depth = 0, state_depth = NA_character_)]
  allA = lapply(ns$row_id, function(r) P[[paste0("ar_", r)]])
  set(P, j = "s_allelic", value = if (length(allA) == 1) allA[[1]] else do.call(pmax, allA))
  ac = allelic_classes(ns)
  AC = do.call(cbind, lapply(ac$states, function(ss) {
    r = which(ns$state %in% ss)
    if (length(r) == 1) allA[[r]] else do.call(pmax, allA[r])
  }))
  set(P, j = "state_allelic", value = unname(ac$label)[max.col(AC, ties.method = "first")])
  if (loo) {
    for (s in sn) set(P, j = paste0("lo_", s), value = lo[[s]])
    P[, (rcols) := NULL]
  }
  pool_combine(P, sn, loo, scale)
}

# best state of pooled state scores (cl_<state>, and lo_<state> with loo): s_comb / state_comb, s_loo (leave-one-out
# score of the called state), s_rob = pmin(s_comb, s_loo), state_rob. The state is always the likelihood argmax. With
# `scale` ([pool_scale()]) s_comb and s_loo are then standardised per region, (score - med) / mad of the region's null
# pools' best LLR (mad x tail when the scale has a `tail` column), and s_comb_llr keeps the unstandardised best LLR.
# A region without null pools takes the class medians of med / mad (tail 1), those of its `sd_class` (when `P` has
# that column) if its class has no null pools at all.
pool_combine = function(P, sn, loo, scale = NULL) {
  cl = lapply(sn, function(s) P[[paste0("cl_", s)]]); names(cl) = sn
  set(P, j = "s_comb", value = if (length(cl) == 1) cl[[1]] else do.call(pmax, unname(cl)))
  set(P, j = "state_comb", value = sn[max.col(do.call(cbind, cl), ties.method = "first")])
  if (loo) {
    LM = do.call(cbind, lapply(sn, function(s) P[[paste0("lo_", s)]]))
    set(P, j = "s_loo", value = LM[cbind(seq_len(nrow(P)), match(P$state_comb, sn))])
  }
  if (!is.null(scale)) {
    sca = as.data.table(scale)
    set(P, j = "s_comb_llr", value = P$s_comb)
    k = match(P$region, sca$region)
    med = sca$med[k]; mad_ = sca$mad[k]
    if ("tail" %in% names(sca)) {                               # per-region tail factor (pool_scale tail_q > 0)
      tl = sca$tail[k]
      tl[is.na(tl) | !is.finite(tl) | tl <= 0] = 1
      mad_ = mad_ * tl
    }
    fb = which(is.na(k))
    if (length(fb)) {                                           # regions without null pools: class medians (tail 1)
      cm = sca[, .(m = median(med), d = median(mad)), by = class]
      kc = match(P$class[fb], cm$class)
      if ("sd_class" %in% names(P)) {                           # a class without null pools (events): its sd_class
        ks = match(P$sd_class[fb], cm$class)
        kc = fifelse(is.na(kc), ks, kc)
      }
      med[fb] = cm$m[kc]; mad_[fb] = cm$d[kc]
      med[is.na(med)] = 0; mad_[is.na(mad_)] = 1
    }
    set(P, j = "s_comb", value = (P$s_comb - med) / mad_)
    if (loo) set(P, j = "s_loo", value = (P$s_loo - med) / mad_)
  }
  if (loo) {
    set(P, j = "s_rob", value = pmin(P$s_comb, P$s_loo))
    set(P, j = "state_rob", value = P$state_comb)
  }
  P[]
}

#' Per-region scale of pooled null scores
#'
#' Median and MAD (`mad()`, normal-consistent) of the pooled best LLR (the max
#' over the states' `cl_<state>` of [pool_scores()], i.e. the raw `s_comb`)
#' over the null pools of each region. Pooled scores standardised with it,
#' `(s - med) / mad`, put every region of a class on one scale, so a per-class
#' pooled threshold neither over-calls regions whose null pools are shifted or
#' wide (molecule-poor regions where the non-carriers of a pool add little
#' negative LLR, shared depth / phasing artefacts, cell-type specific depth in
#' kNN neighbourhoods) nor ignores quiet ones. The state stays the likelihood
#' argmax. A region with fewer than `min_pools` null pools takes the class
#' medians; MADs are floored at `mad_floor` x the class median MAD.
#'
#' Tail factor (`tail_q > 0`, off by default): the median / MAD standardise the body of a region's null pools,
#' but a heavy-tailed region keeps collecting calls at a per-class threshold. Each region's null pools are
#' standardised by its `med` / `mad`, `u_r` = their type-7 `tail_q` quantile and `u_c` = the median of `u_r`
#' over the class's regions with at least `min_pools` pools. With `ne = min(n_pools, n_eff)` (overlapping kNN
#' pools: one per cell, so the effective number of independent pools is ~ cells / pool size) and
#' `w = ne / (ne + tail_k)`, `tail = exp(w log(u_r / u_c))`; `tail = 1` when `u_r` or `u_c` is not positive
#' or the region has fewer than `min_pools` pools. [pool_scores()] / `pool_combine()` then divide by
#' `mad x tail`.
#' @param ps pooled null scores ([pool_scores()] output with `region`, `class`, `cl_<state>`)
#' @param states state set
#' @param min_pools minimal null pools of a region
#' @param mad_floor floor of a region's MAD relative to its class median
#' @param tail_q quantile of the per-region tail factor (0 = no tail factor)
#' @param tail_k shrinkage pseudo-count of the tail factor (effective pools)
#' @param n_eff effective number of independent null pools per region (caps `n_pools` in the shrinkage weight)
#' @return data.table region, class, n_pools, med, mad (+ `tail`, `u_tail` when `tail_q > 0`)
#' @export
pool_scale = function(ps, states = caller_states("normal"), min_pools = 30, mad_floor = 0.1, tail_q = 0, tail_k = 100, n_eff = Inf) {
  ps = as.data.table(ps)
  sn = unique(as.data.table(states)[state != "neutral"]$state)
  s = do.call(pmax, unname(lapply(sn, function(x) ps[[paste0("cl_", x)]])))
  out = data.table(region = ps$region, class = ps$class, s_ = s)[, .(n_pools = .N, med = median(s_), mad = stats::mad(s_)), by = .(region, class)]
  cm = out[n_pools >= min_pools, .(m = median(med), d = median(mad)), by = class]
  k = match(out$class, cm$class)
  few = which(out$n_pools < min_pools & !is.na(k))
  set(out, i = few, j = "med", value = cm$m[k[few]])
  set(out, i = few, j = "mad", value = cm$d[k[few]])
  fl = mad_floor * cm$d[k]
  v = pmax(out$mad, fifelse(is.na(fl), 0, fl))
  v[!is.finite(v) | v <= 0] = 1
  set(out, j = "mad", value = v)
  if (!is.null(tail_q) && is.finite(tail_q) && tail_q > 0) {
    if (tail_q >= 1) stop("pool_scale: tail_q must be in [0, 1)")
    kr = match(ps$region, out$region)
    u = data.table(region = ps$region, u_ = (s - out$med[kr]) / out$mad[kr])[, .(u_tail = q7(u_, tail_q)), by = region]
    set(out, j = "u_tail", value = u$u_tail[match(out$region, u$region)])
    uc = out[n_pools >= min_pools & is.finite(u_tail), .(u_c = median(u_tail)), by = class]
    u_c = uc$u_c[match(out$class, uc$class)]
    ne = pmin(out$n_pools, n_eff)
    w = ne / (ne + tail_k)
    ok = out$n_pools >= min_pools & is.finite(out$u_tail) & is.finite(u_c) & out$u_tail > 0 & u_c > 0
    tl = rep(1, nrow(out))
    tl[ok] = exp(w[ok] * log(out$u_tail[ok] / u_c[ok]))
    set(out, j = "tail", value = tl)
  }
  out[]
}
utils::globalVariables(c("n_pools", "s_comb_llr", "s_", "med", "u_", "u_tail", "u_c"))

# summed columns of the pooled scores and the depth / allelic column of every state row: dl_<state> for a
# single-shift state with a depth shift (NA for shift 0), dr_<row> for multi-shift states, ar_<row>
pool_row_terms = function(ns) {
  sn = unique(ns$state)
  multi = vapply(sn, function(s) length(unique(ns$shift[ns$state == s])) > 1, logical(1))
  pre = function(p, x) if (length(x)) paste0(p, x) else character()
  dcols = c(pre("dl_", sn[!multi & vapply(sn, function(s) ns$shift[ns$state == s][1] != 0, logical(1))]),
            pre("dr_", ns$row_id[ns$state %in% sn[multi]]))
  rows = data.table(state = ns$state, row_id = ns$row_id, acol = paste0("ar_", ns$row_id))
  rows[, dcol := fifelse(multi[state], paste0("dr_", row_id), fifelse(ns$shift == 0, NA_character_, paste0("dl_", state)))]
  list(multi = multi, dcols = dcols, acols = rows$acol, rows = rows[])
}
utils::globalVariables(c("dcol", "acol"))

#' Random pools of cells within groups
#'
#' Shuffles the cells of each group (cell type) and splits them into
#' `round(n / n_pool)` pools of (nearly) equal size; groups with fewer than
#' `ceiling(n_pool / 2)` cells give no pool.
#' @param cells cell ids
#' @param group group per cell (NULL = one group)
#' @param n_pool target pool size
#' @param seed RNG seed
#' @return data.table pool, cell, group
#' @export
make_pools = function(cells, group = NULL, n_pool = 20, seed = 1) {
  set.seed(seed)
  g = if (is.null(group)) rep("all", length(cells)) else as.character(group)
  g[is.na(g)] = "NA"
  out = list()
  for (gg in sort(unique(g))) {
    cc = cells[g == gg]
    n = length(cc)
    if (n < ceiling(n_pool / 2)) next
    k = max(1L, round(n / n_pool))
    cc = cc[sample.int(n)]
    out[[gg]] = data.table(pool = paste0(gg, "|pool", rep_len(seq_len(k), n)), cell = cc, group = gg)
  }
  if (!length(out)) return(data.table(pool = character(), cell = character(), group = character()))
  rbindlist(out)
}

#' Chromatin-neighbourhood pools
#'
#' One pool per cell: the cell and its `k - 1` nearest neighbours (Euclidean,
#' [FNN::get.knn()]) in a chromatin embedding ([chromatin_embedding()]). Never
#' build the embedding from copy-number features: pools of cells that share
#' an event would make the pooled null anti-conservative.
#' @param emb cells x components matrix with cell ids as row names
#' @param k pool size (including the centre cell)
#' @return data.table pool, cell, group ("knn"), center
#' @export
knn_pools = function(emb, k = 20) {
  emb = as.matrix(emb)
  if (is.null(rownames(emb))) stop("knn_pools: `emb` needs cell ids as row names")
  n = nrow(emb)
  kk = min(k - 1L, n - 1L)
  if (kk < 1) stop("knn_pools: need at least 2 cells")
  nn = FNN::get.knn(emb, k = kk)$nn.index
  cells = rownames(emb)
  data.table(pool = paste0("knn|", rep(cells, kk + 1L)), cell = c(cells, cells[as.vector(nn)]),
             group = "knn", center = rep(cells, kk + 1L))[order(pool)]
}

# Thresholds of classes without null rows: a class of `d` that has no threshold row with null rows (n_null 0
# everywhere, or no rows: named `event` regions are never null rows) and whose rows' sd_class has thresholds gets a copy
# of every threshold row of that sd_class (folds, strata, sex groups), class set to its own and `borrowed_from` = the
# sd_class; apply_thresholds() matches them on the rows' sd_class. Unchanged (same object) when nothing is borrowed.
borrow_thresholds = function(thr, d) {
  if (!is.data.table(thr)) thr = as.data.table(thr)
  if (is.null(d) || !nrow(thr) || !all(c("class", "sd_class") %in% names(d)) || all(thr$class == "all")) return(thr)
  pairs = unique(data.table(class = as.character(d$class), sd_class = as.character(d$sd_class)))
  has = unique(thr$class[thr$n_null > 0])
  pairs = pairs[!class %in% has & !is.na(sd_class) & sd_class != class & sd_class %in% has]
  if (!nrow(pairs)) return(thr)
  thr = copy(thr)
  if (!"borrowed_from" %in% names(thr)) thr[, borrowed_from := NA_character_]
  add = rbindlist(lapply(seq_len(nrow(pairs)), function(i)
    copy(thr[class == pairs$sd_class[i] & is.na(borrowed_from)])[, `:=`(class = pairs$class[i], borrowed_from = pairs$sd_class[i])]))
  rbind(thr[!class %in% pairs$class], add)
}
utils::globalVariables("borrowed_from")

#' Null thresholds of pooled scores
#'
#' Null pools are built with the same construction as the test pools
#' (`pools_fun`, e.g. random pools within cell type of null cells, or kNN
#' pools of null cells), scored per pool x region, and the type-7 `q`
#' quantile is taken per `by` (class). Random constructions are redrawn
#' (seed + draw) until every class has at least `n_draws` null pool-regions;
#' a deterministic construction (same pools twice, e.g. kNN) is used once.
#' With `cv` one threshold per fold from the null cells of the other folds,
#' plus fold 0 from all null cells. With `loo` the leave-one-out scores (`s_rob`) are computed and can be
#' thresholded, exactly as [call_pools()] scores the test pools.
#' @param d_null null rows scored with `keep_rows = TRUE` (cell, region, class, fold)
#' @param pools_fun function(cells, seed) returning a pools table (pool, cell)
#' @param n_draws minimal null pool-regions per class
#' @param by grouping of the thresholds; with `sexgrp` (autosome / X / Y, [caller_normal()] for males with
#'   `pool_sex_thr`), a group with fewer than `min_region_pools` null pool-regions takes the class threshold
#'   (`pooled_fallback`), and redraws are counted per class
#' @param q null quantile
#' @param cv per-fold thresholds
#' @param states state set used for scoring
#' @param scores pooled scores to threshold
#' @param params [caller_defaults()]
#' @param max_draws maximal redraws
#' @param loo leave-one-out scores (adds `s_rob` to the default `scores`)
#' @param standardize standardise the pooled best LLR per region ([pool_scale()] of the fold's null pools)
#'   before thresholding; the fold-0 scale is returned as attribute `pool_scale` (pass it to [call_pools()]).
#'   With `params$pool_tail_q > 0` the scale carries the per-region tail factor (`tail_k = params$pool_tail_k`,
#'   `n_eff` = the fold's null cells / `params$n_pool`). With
#'   `params$pool_tail_xfit` (default) the tail is cross-fitted: the fold's null cells are split in two halves
#'   (seeded), pools are built within each half (`n_draws` / 2 and `min_region_pools` / 2 each), each half's
#'   pools are standardised by the other half's scale before the quantile, and the returned scale is the mean of
#'   the two half scales (geometric mean of the tails)
#' @return thresholds table (stratum "all"; `n_null` = null pool-regions); with `standardize`, attribute
#'   `pool_scale`
#' @param min_region_pools with `standardize`, random constructions are redrawn until every region also has
#'   this many null pools (a stable median / MAD per region; fewer pools make the in-sample standardised
#'   threshold anti-conservative)
#' @export
pooled_thresholds = function(d_null, pools_fun, n_draws = 20000, by = "class", q = params$null_q, cv = TRUE,
                             states = caller_states("normal", params), scores = c("s_depth", "s_allelic", "s_comb", if (loo) "s_rob"),
                             params = caller_defaults(), max_draws = 200, loo = FALSE, standardize = FALSE, min_region_pools = 1000) {
  d_null = as.data.table(d_null)
  sn = unique(as.data.table(states)[state != "neutral"]$state)
  keep = unique(c(by, scores, if (standardize) c("region", "class", paste0("cl_", sn), if (loo) paste0("lo_", sn))))
  scale0 = NULL
  folds = if (cv && "fold" %in% names(d_null)) setdiff(sort(unique(d_null$fold[!is.na(d_null$fold)])), 0L) else integer()
  tail_q = params$pool_tail_q %||% 0
  xfit = standardize && tail_q > 0 && isTRUE(params$pool_tail_xfit %||% TRUE)
  scale_fun = function(A, nc) pool_scale(A, states, tail_q = tail_q, tail_k = params$pool_tail_k %||% 30, n_eff = nc / params$n_pool)
  draw = function(tr, cells, seed0, n_draws, min_region_pools) {
    acc = list(); prev = NULL
    for (it in seq_len(max_draws)) {
      pl = as.data.table(pools_fun(cells, seed0 + it))
      if (nrow(pl) == 0) break
      sig = paste(pl$pool, pl$cell)
      if (!is.null(prev) && identical(sig, prev)) break
      prev = sig
      ps = pool_scores(tr, pl, states, by = unique(c("region", "class", by)), loo = loo)
      acc[[it]] = ps[, keep, with = FALSE]
      cnt = rbindlist(acc)[, .N, by = setdiff(by, "sexgrp")]   # draws are counted per class (sex groups fall back)
      enough = min(cnt$N) >= n_draws
      if (standardize) enough = enough && min(rbindlist(lapply(acc, function(a) a[, .(region)]))[, .N, by = region]$N) >= min_region_pools
      if (enough) break
    }
    rbindlist(acc)
  }
  out = list()
  for (fo in c(0L, folds)) {
    tr = if (fo == 0L) d_null else d_null[!is.na(fold) & fold != fo]
    cells = unique(tr$cell)
    if (!xfit) {
      A = draw(tr, cells, params$seed + 1000L * fo, n_draws, min_region_pools)
      if (nrow(A) == 0) next
      if (standardize) {
        sca = scale_fun(A, length(cells))
        A = pool_combine(A, sn, loo, sca)
        if (fo == 0L) scale0 = sca
      }
    } else {
      # cross-fitted tail: two disjoint halves of the training null cells, pools built within each half; each
      # half's pools are standardised by the other half's scale (median / MAD / tail), so the tail factor and
      # the threshold quantile never see the same pools. The fold's scale (fold 0: returned for calling) is the
      # average of the two half scales (med, mad: mean; tail: geometric mean).
      set.seed(params$seed + 7919L + fo)
      sh = sample(cells)
      halves = split(sh, rep(1:2, length.out = length(sh)))
      Ah = lapply(1:2, function(h) draw(tr[cell %in% halves[[h]]], halves[[h]], params$seed + 1000L * fo + 500L * h,
                                        ceiling(n_draws / 2), ceiling(min_region_pools / 2)))
      if (min(vapply(Ah, nrow, 1L)) == 0) next
      Sh = lapply(1:2, function(h) scale_fun(Ah[[h]], length(halves[[h]])))
      A = rbind(pool_combine(Ah[[1]], sn, loo, Sh[[2]]), pool_combine(Ah[[2]], sn, loo, Sh[[1]]), fill = TRUE)
      sca = rbind(Sh[[1]], Sh[[2]], fill = TRUE)[, .(n_pools = sum(n_pools), med = mean(med), mad = mean(mad),
                                                     u_tail = mean(u_tail), tail = exp(mean(log(tail)))), by = .(region, class)]
      if (fo == 0L) scale0 = sca
    }
    for (s in scores) {
      r = A[, .(thr = q7(get(s), q), n_null = .N), by = by]
      r[, pooled_fallback := FALSE]
      if ("sexgrp" %in% by) {                                   # sex-chromosome groups: class threshold below min_region_pools
        coarse = setdiff(by, "sexgrp")
        rc = A[, .(thr_c = q7(get(s), q)), by = coarse]
        k = match(group_key(r, coarse), group_key(rc, coarse))
        few = r$n_null < min_region_pools
        r[few, `:=`(thr = rc$thr_c[k[few]], pooled_fallback = TRUE)]
      }
      r[, `:=`(score = s, fold = fo)]
      out[[length(out) + 1]] = r
    }
  }
  res = rbindlist(out)
  if (!"class" %in% by) res[, class := "all"]
  if (!"stratum" %in% by) res[, stratum := "all"]
  res[, `:=`(q = q, n_null = as.integer(n_null), fold = as.integer(fold))]
  setcolorder(res, c("score", "class", "stratum", "fold", "q", "thr", "n_null", "pooled_fallback"))
  if (standardize) setattr(res, "pool_scale", scale0)
  res[]
}

#' Pooled calls
#'
#' Scores the pools ([pool_scores()] per pool x region), applies the fold-0
#' pooled thresholds per class and returns the policy's pool-level rows. With
#' `loo` the pools are scored as [pool_scores()] describes; use
#' `score = "s_rob"` (thresholded by [pooled_thresholds()] with `loo = TRUE`) so
#' a call must survive dropping its most influential cell.
#' @param d rows scored with `keep_rows = TRUE`
#' @param pools pools table (pool, cell)
#' @param thr [pooled_thresholds()] output
#' @param policy [call_policy()]
#' @param states state set used for scoring
#' @param regions optional regions table (coordinates)
#' @param score pooled score
#' @param loo leave-one-out scores (needed for `score = "s_rob"`)
#' @param scale NULL or the [pool_scale()] table the thresholds were computed with (attribute `pool_scale` of
#'   [pooled_thresholds()] with `standardize = TRUE`). A region of a class without null pools (named `event`
#'   regions) is standardised with the class medians of its `sd_class` when `thr` has borrowed rows
#' @return calls table (`level = "pool"`, `cell` = pool id) plus `n` (cells in the pool), `s_comb` and,
#'   with `loo`, `s_loo`; with `scale`, `s_comb_llr` (unstandardised best LLR)
#' @export
call_pools = function(d, pools, thr, policy = call_policy("normal"), states = caller_states("normal"), regions = NULL, score = "s_comb",
                      loo = score %in% c("s_rob", "s_loo"), scale = NULL) {
  d = as.data.table(d)
  if ("eligible" %in% names(d)) d = d[eligible %in% TRUE]
  # borrowed thresholds (a class without null pools, borrow_thresholds()) are matched on the rows' sd_class
  brw = "borrowed_from" %in% names(thr) && any(!is.na(thr$borrowed_from)) && "sd_class" %in% names(d)
  ps = pool_scores(d, pools, states, by = c("region", "class", if (brw) "sd_class"), loo = loo, scale = scale)
  if ("sexgrp" %in% names(thr) && !all(is.na(thr$sexgrp))) {
    # the rows' chromosome (a named event has none in its name)
    ch = if ("chrom" %in% names(d)) d$chrom[match(ps$region, d$region)] else region_chrom(ps$region, ps$class)
    ps[, sexgrp := sex_group(ch)]
  }
  ps = apply_thresholds(ps, thr, score)
  pol = as.data.table(policy)[level == "pool"]
  sc = sub("^s_", "state_", score)
  keep = paste(ps$class, ps[[sc]]) %in% paste(pol$class, pol$state)
  x = data.table(cell = as.character(ps$pool[keep]), region = ps$region[keep], class = ps$class[keep], state = ps[[sc]][keep],
                 score = ps[[score]][keep], thr = ps[[paste0("thr_", score)]][keep], n = ps$n[keep], s_comb = ps$s_comb[keep])
  if (loo) x[, s_loo := ps$s_loo[keep]]
  if (!is.null(scale)) x[, s_comb_llr := ps$s_comb_llr[keep]]
  if ("chrom" %in% names(d)) x[, chrom := d$chrom[match(region, d$region)]]
  x = add_coords(x, regions)
  x[, called := !is.na(score) & !is.na(thr) & score > thr]
  x[, `:=`(level = "pool", source = "cell")]
  x[, reportable := policy_reportable(x, pol)]
  setcolorder(x, c("cell", "region", "class", "chrom", "start", "end", "state", "score", "thr", "called", "level", "source", "reportable"))
  check_calls(x)
  x[]
}

# ---- normal-tissue caller ------------------------------------------------------------------------------

#' Normal-tissue copy-number caller
#'
#' 1. null cells: the features' `null_cell` column, or [select_null_cells()]
#'    (`null_cells = "auto"` selects when the column is missing, all NA or all
#'    TRUE), never `multiplets`; 2. Student-t depth null per class x stratum
#'    and allelic rho on the null arm rows; 3. [score_cells()] with the normal
#'    state set; 4. fold-CV budget thresholds ([budget_thresholds()]); 5.
#'    per-cell calls by [call_policy()]; 6. pools within cell type
#'    ([make_pools()]) and, with an `embedding`, kNN pools ([knn_pools()]),
#'    pooled null thresholds ([pooled_thresholds()]) and pooled calls; 7.
#'    events ([consolidate_calls()]). Rows with `eligible = FALSE` are dropped first.
#'
#' Calibration policy (defaults): null cells are selected on autosomal rows
#' and null rows are autosomal ([select_null_cells()] `chroms`, `null_chroms`);
#' per-cell chromosome gains need depth support and lose to a stronger nested
#' loss / CN-LOH ([resolve_nested()], `resolve_gains`; opt-in `gain_gate = "z"` gates on a null quantile of z
#' instead of the sign of the gain depth LLR); pooled state scores are
#' standardised per region by the median / MAD of the null pools of the same
#' construction ([pool_scale()], `pool_standardize`), and a pooled call must
#' survive dropping its most influential cell (`pool_loo`: pools are called on
#' `s_rob`, thresholded on null pools of the same statistic); chrX pools are
#' not reportable in females ([call_policy()] `sex`).
#'
#' Opt-in calibrations: `params$calib = "votes"` keys the per-cell thresholds on class x stratum x vote stratum
#' (`hstratum`, bins of n_hap of the null rows, [vote_bins()]); `params$calib_rho = "chrom"` scores the allele
#' channel with a per-chromosome rho shrunk to the global ([rho_by_chrom()]). `params$pool_tail_q > 0` (the
#' default) adds a per-region tail factor to the pooled region scale ([pool_scale()]), cross-fitted with
#' `params$pool_tail_xfit` ([pooled_thresholds()]). Named `event` regions (no null rows) borrow the per-cell
#' thresholds of their `sd_class` ([budget_thresholds()]); their pools are called only with `params$pool_events`
#' (experimental, off by default: [caller_defaults()]). Removed options stop with
#' an error: `calib` "region" / "region_votes" / "votes_cont", `calib_rho` "cell" / "chrom_cell".
#' @param feat features long table (`vignette("scphasedcn-outputs", package = "scphasedcn")`)
#' @param params [caller_defaults()]
#' @param cell_meta optional table `cell` (+ `cell_type`) overriding the features' cell types
#' @param multiplets cells flagged as multiplets (never null cells, never pooled)
#' @param regions regions table (default [regions_from_features()])
#' @param embedding optional chromatin embedding (cells x components, row names = cells) for kNN pools
#' @param null_cells "auto", "features" or "select"
#' @param depth_model "t" or "gauss"
#' @param n_draws null pool-regions per class for the pooled thresholds
#' @param pool_classes region classes that are pooled
#' @param select_args extra arguments for [select_null_cells()]
#' @param sex "female", "male" (case-insensitive, or "F" / "M"; any other value stops) or NULL (inferred: female when
#'   chrX rows carry allele votes); sets the chrX rule of [call_policy()]
#' @param resolve_gains apply [resolve_nested()] to the per-cell calls
#' @param gain_min_dl `min_dl_gain` of [resolve_nested()] (`gain_gate = "dl"`)
#' @param gain_gate depth gate of the per-cell gains, `gate` of [resolve_nested()]: `"dl"` (default: gain depth
#'   LLR above `gain_min_dl`, i.e. z above the midpoint of neutral and the gain shift) or `"z"` (opt-in: z also
#'   above the `gain_gate_q` quantile of the depth null of the row's class x stratum. About `1 - gain_gate_q` of
#'   all null-depth rows pass; among rows already called as gains with weak allelic support, and on female chrX,
#'   up to about two to three times that in simulation. Not benchmarked on real data as a default)
#' @param gain_gate_q null quantile of the `"z"` gain gate (0.95)
#' @param pool_loo call pools on the leave-one-out robust score `s_rob` ([pool_scores()])
#' @param pool_standardize standardise the pooled best LLR per region ([pool_scale()] of the null pools of the
#'   same construction) before the per-class pooled thresholds
#' @param null_chroms chromosomes whose rows can be null rows (default the autosomes, the chromosomes
#'   [select_null_cells()] checks); NULL = every `null_region`
#' @return `scphasedcn` (mode normal) with `scores`, `thresholds` (column `level`:
#'   cell, pool_type, pool_knn), `calls` (cell and pool level; `resolved` for gains changed by
#'   [resolve_nested()]; with `gain_gate = "z"` also its `gate_z`, `gate_thr`), `events`, `pools`, `null_fit = list(tpar, rho, null_cells, selection, pool_scale, sex,
#'   vote_bins, rho_chrom)` (`vote_bins`: the [vote_bins()] table with `params$calib = "votes"`; `rho_chrom`: the
#'   [rho_by_chrom()] table with `params$calib_rho = "chrom"`)
#' @export
caller_normal = function(feat, params = caller_defaults(), cell_meta = NULL, multiplets = NULL, regions = NULL, embedding = NULL,
                         null_cells = c("auto", "features", "select"), depth_model = c("t", "gauss"), n_draws = 20000,
                         pool_classes = c("chrom", "arm", "50Mb", "20Mb", "10Mb", "5Mb", "event"), select_args = list(),
                         sex = NULL, resolve_gains = TRUE, gain_min_dl = 0, gain_gate = "dl", gain_gate_q = 0.95,
                         pool_loo = TRUE, pool_standardize = TRUE, null_chroms = paste0("chr", 1:22)) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  null_cells = match.arg(null_cells); depth_model = match.arg(depth_model)
  sex = check_sex(sex, "caller_normal")
  if (!identical(length(gain_gate), 1L) || !gain_gate %in% c("dl", "z"))
    stop("caller_normal: gain_gate must be 'dl' or 'z', not '", paste(gain_gate, collapse = ", "), "'", call. = FALSE)
  if (gain_gate == "z" && (!is.numeric(gain_gate_q) || length(gain_gate_q) != 1 || is.na(gain_gate_q) || gain_gate_q <= 0 || gain_gate_q >= 1))
    stop("caller_normal: gain_gate_q must be one number in (0, 1)", call. = FALSE)
  d = copy(as.data.table(feat))
  need = c("cell", "region", "class", "stratum", "z", "sd_null", "nA", "n_hap")
  miss = setdiff(need, names(d))
  if (length(miss)) stop("caller_normal: features lack ", paste(miss, collapse = ", "))
  d[, cell := as.character(cell)]
  if (!"sd_class" %in% names(d)) d[, sd_class := class]
  if (!"null_region" %in% names(d)) d[, null_region := TRUE]
  if (!"eligible" %in% names(d)) d[, eligible := TRUE]
  if (!"chrom" %in% names(d)) d[, chrom := region_chrom(region, class)]
  if (!is.null(cell_meta) && "cell_type" %in% names(cell_meta)) d[, cell_type := as.character(cell_meta$cell_type)[match(cell, cell_meta$cell)]]
  if (!"cell_type" %in% names(d)) d[, cell_type := "all"]
  if (!"fold" %in% names(d) || anyNA(d$fold)) d[, fold := assign_folds(cell, params$n_folds, params$seed, if ("origin" %in% names(d)) origin else NULL)]
  d = d[eligible %in% TRUE]                              # ineligible regions are not scored (see the regions table)
  multiplets = intersect(as.character(multiplets), d$cell)
  if (is.null(regions)) regions = regions_from_features(d, params)
  if (!"allelic_ok" %in% names(d) && "allelic_ok" %in% names(regions)) {   # depth-only regions (male chrX / chrY)
    ao = as.logical(regions$allelic_ok)[match(d$region, regions$region)]
    if (any(ao %in% FALSE)) d[, allelic_ok := !(ao %in% FALSE)]
  }
  # 1. null cells
  nc_col = if ("null_cell" %in% names(d)) d$null_cell else NULL
  use_sel = null_cells == "select" || (null_cells == "auto" && (is.null(nc_col) || all(is.na(nc_col)) || all(nc_col %in% TRUE)))
  sel = NULL
  if (use_sel) {
    sa = c(list(feat_long = d, params = params, exclude = multiplets, depth_model = depth_model), select_args)
    if (!"chroms" %in% names(sa) && !is.null(null_chroms)) sa$chroms = null_chroms
    sel = do.call(select_null_cells, sa)
    d[, null_cell := cell %in% sel$null_cells]
  } else {
    if (is.null(nc_col)) stop("caller_normal: null_cells = 'features' but no null_cell column")
    d[, null_cell := !is.na(null_cell) & null_cell & !cell %in% multiplets]
  }
  d[, is_null := null_cell & !is.na(null_region) & null_region & !is.na(eligible) & eligible]
  if (!is.null(null_chroms)) d[!chrom %in% null_chroms, is_null := FALSE]      # null rows only where the null cells were checked
  # 2. null fits
  calib = params$calib %||% "class"
  if (!identical(length(calib), 1L) || !calib %in% c("class", "votes"))
    stop("caller_normal: params$calib must be 'class' or 'votes', not '", paste(calib, collapse = ", "), "'",
         call. = FALSE)
  calib_rho = params$calib_rho %||% "global"
  if (!identical(length(calib_rho), 1L) || !calib_rho %in% c("global", "chrom"))
    stop("caller_normal: params$calib_rho must be 'global' or 'chrom', not '", paste(calib_rho, collapse = ", "), "'",
         call. = FALSE)
  pool_events = params$pool_events %||% FALSE
  if (!is.logical(pool_events) || length(pool_events) != 1L || is.na(pool_events))
    stop("caller_normal: params$pool_events must be TRUE or FALSE, not '", paste(format(pool_events), collapse = ", "), "'",
         call. = FALSE)
  vbins = NULL
  if (calib == "votes") {                                      # thresholds per class x stratum x vote stratum
    vbins = vote_bins(d, params$vote_bins %||% 4, params$min_null_rows, null = d$is_null, params = params)
    d[, hstratum := assign_hstratum(d, vbins)]
  }
  tpar = if (depth_model == "t") complete_tpar(d, fit_depth_null(d[is_null == TRUE])) else NULL
  arm = d$is_null & d$class == "arm"
  rho = if (any(arm & d$n_hap > 0, na.rm = TRUE)) fit_allelic_null(d$nA[arm], d$n_hap[arm]) else fit_allelic_null(d$nA[d$is_null], d$n_hap[d$is_null])
  rho_tab = NULL; rho_sc = rho
  novotes_note = NULL
  if (is.na(rho)) {
    # no null row has an allele vote: without any vote at all (an unphased extraction, a VCF without phase sets, no
    # allele layer) every allelic term is 0 and the calls are depth-only; votes elsewhere but not on the null rows
    # cannot be scored without a rho
    if (any(d$n_hap > 0, na.rm = TRUE))
      stop("caller_normal: the null rows have no allele votes, so the allelic overdispersion (rho) cannot be fitted",
           call. = FALSE)
    novotes_note = "no allele votes (unphased extraction or no phase sets?): the calls are depth-only"
    warning("caller_normal: ", novotes_note, call. = FALSE)
  }
  if (calib_rho == "chrom" && is.null(novotes_note)) {                                  # per-chromosome rho, shrunk to the global
    rho_tab = rho_by_chrom(d, rho, params$rho_chrom_k %||% 100, null = d$is_null, params = params)
    rho_sc = rho_rows(d, rho_tab, rho)
  }
  # 3-5. scores, thresholds, per-cell calls
  states = caller_states("normal", params)
  scores = c("s_depth", "s_allelic", "s_comb")
  sc = score_cells(d, states, tpar, rho_sc, depth_model = depth_model, params = params, keep_rows = TRUE)
  thr = budget_thresholds(sc, scores, q = params$null_q, cv = TRUE, min_null = params$min_null_rows, params = params)
  sc = apply_thresholds(sc, thr, scores)
  if (is.null(sex)) sex = infer_sex(d)
  policy = call_policy("normal", sex)
  calls = call_table(sc, NULL, regions, policy)
  if (resolve_gains) calls = resolve_nested(calls, sc, min_dl_gain = gain_min_dl, gate = gain_gate, gate_q = gain_gate_q,
                                            tpar = tpar, params = params)
  thr[, level := "cell"]
  # 6. pools
  pc = setdiff(unique(d$cell), multiplets)
  ctype = d[, .(cell_type = cell_type[1]), by = cell]
  ct = setNames(ctype$cell_type, ctype$cell)
  sp = sc[class %in% pool_classes]
  # pooled null: the null rows, plus for a male the chrX / chrY rows of the null cells (one copy is their null state),
  # so that the sex-chromosome pools get their own region scale instead of the autosomal class fallback. Only classes
  # with null rows take part: a class without any (named events) borrows the thresholds of its sd_class in every sex
  # group (borrow_thresholds()) instead of holding sex-chromosome null pools that would stop it from borrowing
  sexnull = identical(sex, "male") & sp$chrom %in% c("chrX", "chrY") & sp$null_cell %in% TRUE & sp$eligible %in% TRUE &
    sp$class %in% sp$class[sp$is_null %in% TRUE]
  pnull = sp[is_null == TRUE | sexnull]
  # male: the sex chromosomes get their own pooled threshold group (class x {autosome, X, Y}); their depth-only,
  # neighbourhood-sensitive null is heavier-tailed than the autosomal one that sets the class threshold
  pby = "class"
  if (isTRUE(params$pool_sex_thr %||% TRUE) && any(sexnull)) {
    pnull[, sexgrp := sex_group(chrom)]
    pby = c("class", "sexgrp")
  }
  pools = list(); pthr = list(); pcalls = list(); pscale = list()
  # classes without null pools (named events) take the pooled thresholds of their sd_class only with pool_events
  # (experimental: their pools are standardised with the sd_class's median scale, not their own); otherwise they get
  # no pooled threshold and their pools are scored but never called
  pborrow = function(th) if (pool_events) borrow_thresholds(th, sp) else th
  preg = unique(sp[, .(region, class)])
  pscore = if (pool_loo) "s_rob" else "s_comb"
  pools$type = make_pools(pc, ct[pc], params$n_pool, params$seed)
  if (nrow(pools$type) && nrow(pnull)) {
    fun = function(cells, seed) make_pools(cells, ct[cells], params$n_pool, seed)
    th_ = pooled_thresholds(pnull, fun, n_draws = n_draws, by = pby, q = params$null_q, states = states, params = params,
                                  loo = pool_loo, standardize = pool_standardize)
    pscale$type = attr(th_, "pool_scale")
    pthr$type = pborrow(th_)[, level := "pool_type"]
    pcalls$type = call_pools(sp, pools$type, pthr$type, policy, states, regions, score = pscore, loo = pool_loo,
                             scale = pscale$type)
  }
  if (!is.null(embedding)) {
    emb = as.matrix(embedding)
    pk = intersect(pc, rownames(emb))
    pools$knn = knn_pools(emb[pk, , drop = FALSE], params$n_pool)
    nk = intersect(unique(pnull$cell), rownames(emb))
    if (length(nk) > params$n_pool) {
      fun = function(cells, seed) knn_pools(emb[intersect(cells, rownames(emb)), , drop = FALSE], params$n_pool)
      th_ = pooled_thresholds(pnull[cell %in% nk], fun, n_draws = n_draws, by = pby, q = params$null_q, states = states, params = params,
                                   loo = pool_loo, standardize = pool_standardize)
      pscale$knn = attr(th_, "pool_scale")
      pthr$knn = pborrow(th_)[, level := "pool_knn"]
      pcalls$knn = call_pools(sp, pools$knn, pthr$knn, policy, states, regions, score = pscore, loo = pool_loo,
                              scale = pscale$knn)
    }
  }
  all_calls = rbindlist(c(list(calls), pcalls), fill = TRUE)
  events = consolidate_calls(all_calls, regions)
  # regions left without any threshold after borrowing (a class without null rows whose sd_class has no thresholds
  # either, e.g. an event whose sd_class is not among the region classes): scored but never called. Without
  # pool_events the pools of classes without null pools are left uncalled by design and not listed
  na_thr = function(x, col) if (is.null(x) || !nrow(x)) character() else x[, .(nothr_ = all(is.na(get(col)))), by = region][nothr_ == TRUE, region]
  nt_cell = na_thr(sc, "thr_s_comb")
  pc_ = rbindlist(pcalls, fill = TRUE)
  if (!pool_events && nrow(pc_)) pc_ = pc_[class %in% pnull$class]
  nt_pool = na_thr(pc_, "thr")
  nothr = unique(c(nt_cell, nt_pool))
  nothr_note = NULL
  if (length(nothr)) {
    lab = vapply(nothr, function(r) paste0(r, " (", paste(c(if (r %in% nt_cell) "per cell", if (r %in% nt_pool) "pooled"),
                                                             collapse = ", "), ")"), character(1))
    nothr_note = sprintf("no threshold for %d region(s), never called: %s (no null rows in their class or sd_class)",
                         length(nothr), paste(lab, collapse = ", "))
    warning("caller_normal: ", nothr_note, call. = FALSE)
  }
  cells = d[, .(cell_type = cell_type[1], stratum = stratum[1], fold = fold[1], null_cell = null_cell[1]), by = cell]
  cells[, multiplet := cell %in% multiplets]
  new_scphasedcn("normal", cells = cells, regions = regions, scores = sc,
              thresholds = rbindlist(c(list(thr), pthr), fill = TRUE), calls = all_calls, events = events,
              pools = rbindlist(lapply(names(pools), function(k) pools[[k]][, kind := k]), fill = TRUE),
              null_fit = list(tpar = tpar, rho = rho, null_cells = cells[null_cell == TRUE, cell], selection = sel,
                              pool_scale = if (length(pscale)) rbindlist(lapply(names(pscale), function(k) copy(pscale[[k]])[, kind := k])) else NULL, sex = sex,
                              vote_bins = vbins, rho_chrom = rho_tab),
              params = params,
              notes = c(sprintf("null cells: %s (%d)", if (use_sel) "select_null_cells" else "features", sum(cells$null_cell)),
                        sprintf("rho %.5g", rho), novotes_note,
                        if (!is.null(rho_tab)) sprintf("allelic rho per chromosome (rho_chrom_k %g): %.3g-%.3g (global %.3g)", params$rho_chrom_k %||% 100,
                                                       min(rho_tab$rho), max(rho_tab$rho), rho),
                        if (!is.null(vbins)) sprintf("cell thresholds per class x stratum x vote stratum (n_hap bins, up to %d per class x stratum; %d bins)",
                                                     as.integer(params$vote_bins %||% 4), nrow(vbins)),
                        if ((params$pool_tail_q %||% 0) > 0) sprintf("pooled tail calibration q %g (pool_tail_k %g, cross-fitted %s)", params$pool_tail_q,
                                                                params$pool_tail_k %||% 30, isTRUE(params$pool_tail_xfit %||% TRUE)),
                        sprintf("sex %s; gains resolved %s (%s); pools: region-standardised %s, called on %s",
                                sex %||% "unknown", resolve_gains,
                                if (gain_gate == "z") sprintf("depth gate z: above the midpoint and the %g null quantile", gain_gate_q)
                                else sprintf("min dl_gain %g", gain_min_dl), pool_standardize, pscore),
                        if (pool_events && any(vapply(pthr, function(x) any(!is.na(x$borrowed_from %||% NA)), logical(1))))
                          paste("pooled calls of named events (pool_events, experimental): standardised with the class-median",
                                "scale of their sd_class and thresholded with its pooled threshold, not calibrated per event"),
                        nothr_note))
}

# sex for the chrX policy: female when chrX rows have allele votes (males are hemizygous: depth-only chrX), NULL without chrX rows
infer_sex = function(d) {
  x = d[chrom == "chrX"]
  if (!nrow(x)) return(NULL)
  if ("allelic_ok" %in% names(x) && !all(is.na(x$allelic_ok))) return(if (any(x$allelic_ok %in% TRUE)) "female" else "male")
  if (mean(x$n_hap > 0, na.rm = TRUE) > 0.5) "female" else "male"
}
utils::globalVariables(c("kind", "multiplet", "cc_", "s_loo", "nothr_"))
