# clones_tumour.R
# What: tumour mode of the copy-number caller, part 4 of 4: the tumour-mode caller
#       [caller_tumour()] and its options (clone_defaults, check_clone_opts), per-cell state scores, the
#       adapters to the scorer for scoring / confirmation (b_*), and per-cell calls inherited from / confirmed
#       against the clone states. The other parts: clones_discovery.R (clone discovery on signed allelic +
#       depth features, kNN-smoothed PCA, Ward tree over-cut, sibling merges at a permutation-calibrated
#       1 / 1,000 split test), clones_states.R (per-clone profiles, nA:nB states, Viterbi segments) and
#       clones_assign.R (likelihood re-assignment with 50:50 mixtures as doublet hypotheses, doublet clones,
#       allelic cleaning of the normal reference pool, automatic control regions).
# Why: in tumours per-cell detection is reliable but per-cell states are not (attenuated haplotype fractions,
#      allele-specific amplicons, a heavy-tailed per-cell depth term); clone-level depth and orientation settle
#      the state, and cells are then scored under the clone states. In a tumour sample the pool of
#      copy-number-flat cells can carry tumour haplotypes, so the normal reference is cleaned with the allele
#      channel.
# Inputs: an `scn_features`-like list (see [caller_tumour()] for the components used) or the
#         individual tables (features, regions, votes, atoms; vignette("scphasedcn-outputs", package = "scphasedcn")).
# Outputs: data.tables (clone feature matrix, clones, profiles, state log-likelihoods, segments,
#          assignments, calls) and a `scphasedcn` object in mode "tumour".
# The globalVariables() block below covers the column names used across all four clones_*.R files.

utils::globalVariables(c(
  "region", "class", "eligible", "cell", "atom", "ps", "nA", "nB", "n_hap", "z", "log2ratio", "sd_null", "chrom",
  "start", "end", "pos", "grp", "gA", "gB", "aA", "aB", "rA", "rB", "s", "oA", "oB", "n_blocks", "clone", "window",
  "state", "ll", "label", "cn", "ploidy", "w", "allelic_ok", "orient", "sd_class", "stratum", "n_frag", "x", "e",
  "posterior", "doublet_like", "assigned", "multiplet", "ref_pool", "i.nA", "i.nB", "blk", "ci", "piece", "len",
  "n_win", "llr", "confident", "f", "log2", "log2_se", "f_se", "n_cells", "pair", "ll_margin", "score", "thr",
  "called", "source", "level", "reportable", "null_cell", "null_region", "w5", "half", "key", "f_state", "mu",
  "shift", "fexp", "n", "seg_start", "seg_end", "is_ref", "flag", "ps_key", "kind", "prior", "nA_s", "nB_s"))

# ---- the tumour-mode caller ----------------------------------------------------------------------------

#' Defaults of the tumour mode
#'
#' Discovery: `classes`, `depth_delta`, `cap`, `min_n`, `n_pc`, `k`, `max_leaves`, `min_size`, `q`, `n_perm`, `max_cells`,
#' `split_depth` (discovery / test halves of the window depth as well as of the votes). Hierarchy: trunks =
#' average-linkage groups of [clone_allelic_distance()] cut at `trunk_dist` (bp x haplotype fraction; e.g. 100 Mb
#' differing by 0.3); trunk cells = trunk posterior >= `trunk_post`; [discover_within()] with `disp_min`, `imb_min`,
#' `deep_q`, `depth_w`, `max_leaves_trunk`, `n_nb_trunk`; final split pass after the EM rounds (`split_pass`): `FALSE`
#' (default) or `"disagree"` = [split_disagree()] (opt-in: members whose own segment state disagrees with the clone
#' state; `split_sep`, `split_max` splits per pass, `split_cand` candidates tested per clone, `split_min_len` minimal
#' segment length, `split_min_df` major-fraction step of an alternative between amplified states, see
#' [split_alt_states()], `split_min_excess` minimal contradiction above its noise floor).
#' Normal hypothesis (`normal_hypothesis`): a fixed all-neutral clone `N0` replaces every fitted clone without a
#' non-neutral segment >= `normal_min_event` bp in each round. States: `size_small`, `small_clone`, `p_switch`, `min_win`, `min_llr`;
#' segment rules of [segment_clone()] (all segmentations of the run, including the split pass): `min_excess` (pooled
#' allelic excess of imbalanced runs), `max_contra` (gains / losses contradicted by depth become copy-neutral), `bal_frac`
#' and `min_bal_len` (support of balanced 2:2 runs).
#' Assignment: `mix_prior`, `mix_w`, `mix_margin`, `mix_rel`,
#' `mix_min_dist` (mixture hypotheses only for clone pairs whose states differ allelically on >= this many bp),
#' `max_mix_prior`, `min_post` (confident), `member_post` (pseudobulk members), `em_rounds`, `same_dist` (clones whose
#' states differ on < this many bp within a trunk are merged). Reference: `ref_post`, `min_ref`. Other:
#' `doublet_delta`, `multiplet_prob` (cells above it are excluded from discovery and pseudobulks), `control_min_len`,
#' `repass`, `allelic_exclude`.
#'
#' `checkpoint` (resume cache, "" = off): an .rds path. When the file does not exist the global discovery
#' (leaves, merged clones) and the first EM round are saved there; when it exists they are read from it and the run
#' continues from the within-trunk discovery (useful when rerunning the later steps on a large sample). The cache is
#' not keyed: reuse it only with the same features, seed and discovery options.
#'
#' Confirmation null (per-cell `confirmed` calls, after assignment; assignment itself always uses the input fits):
#' `ref_null` = "none" (the input rho / depth t of `feat$null_fit`) or "rho" (opt-in: beta-binomial rho refitted on the
#' clean reference pool's piece votes as the confirmation scores them, oriented leave-one-out by a clone pseudobulk;
#' the `call_scn` recipe on the pool's feature votes is kept as `null_fit$rho_pool_features`); needs >= `min_ref`
#' pool cells. `confirm_null` = "generic" (null rows scored as the best of the
#' tumour state set, one threshold per stratum) or "matched" (each clone hypothesis clone | label | shift | f scored on
#' the null rows oriented by that clone's pseudobulk under the same two-state set as the cells; one threshold per
#' hypothesis x stratum, hypothesis-wide below `min_null_rows`; held-out fold-CV false-confirmation rates per
#' hypothesis in `null_fit$confirm_null`; the default). The refitted value is `null_fit$rho_confirm`.
#'
#' A key that is not one of these defaults stops [caller_tumour()] with an error, and so does a value outside its range:
#' `q`, `p_switch`, `trunk_post`, `min_post`, `ref_post`, `max_mix_prior` and every `mix_w` strictly
#' between 0 and 1, `mix_prior` in \[0, 1), `multiplet_prob`, `deep_q` and `member_post` in \[0, 1\], counts (`em_rounds`,
#' `min_size`, `n_perm`, `min_win`, `min_ref`, ...) whole numbers >= 1, the caps `max_cells`, `max_leaves` and
#' `max_leaves_trunk` whole numbers >= 1 or `Inf` (no cap), `small_clone` a whole number >= 0, `trunk_dist` and
#' `control_min_len` >= 0, the other sizes / lengths > 0, the other tolerances >= 0.
#' @return named list (see [caller_tumour()])
#' @export
clone_defaults = function() {
  list(classes = c("arm", "20Mb"), depth_delta = 0.2, cap = 6, min_n = 5, n_pc = 20, k = 20, max_leaves = 40, min_size = 20,
       q = 0.999, n_perm = 200, max_cells = 250, size_small = 5e6, small_clone = 50, split_depth = TRUE,
       trunk_dist = 3e7, trunk_post = 0.9, disp_min = 1.5, imb_min = 1, deep_q = 0.5, depth_w = 1,
       max_leaves_trunk = 20, n_nb_trunk = Inf, split_pass = FALSE, split_sep = 2, split_max = 4, split_cand = 3,
       split_min_len = 1e7, split_min_df = 0.15, split_min_excess = 0.05,
       p_switch = 1e-3, min_win = 3, min_llr = 10, min_excess = 0.5, max_contra = 0.1, bal_frac = 0.75, min_bal_len = 2e7,
       mix_prior = 0.05, mix_w = c(1 / 3, 1 / 2, 2 / 3), mix_margin = 10,
       mix_rel = 0.1, mix_min_dist = 3e7, max_mix_prior = 0.5, min_post = 0.9, member_post = 0.9, em_rounds = 2, same_dist = 5e6,
       ref_post = 0.9, min_ref = 30, doublet_delta = 10, multiplet_prob = 0.9, control_min_len = 2e7,
       repass = TRUE, allelic_exclude = c("chrX", "chrY"), checkpoint = "", normal_hypothesis = TRUE,
       normal_min_event = 1e7, ref_null = "none", confirm_null = "matched")
}

# tumour options: clone_defaults() updated from `opts`; unknown keys stop with an error, option values are checked
# match.arg style, numeric options against their ranges (probabilities strictly inside (0, 1), counts whole numbers
# >= 1, sizes > 0, tolerances >= 0)
check_clone_opts = function(opts) {
  opts = as.list(opts)
  unknown = setdiff(names(opts), names(clone_defaults()))
  if (length(unknown)) stop("caller_tumour: unknown option(s) ", paste(unknown, collapse = ", "), " (see clone_defaults())", call. = FALSE)
  opts = utils::modifyList(clone_defaults(), opts)
  sp = opts$split_pass
  if (!(isFALSE(sp) || identical(sp, "disagree")))
    stop("caller_tumour: opts$split_pass must be FALSE or \"disagree\", not ", paste(format(sp), collapse = ", "),
         call. = FALSE)
  rn = as.character(opts$ref_null)
  if (length(rn) != 1 || !rn %in% c("none", "rho"))
    stop("caller_tumour: opts$ref_null must be \"none\" or \"rho\", not ", paste(rn, collapse = ", "), call. = FALSE)
  opts$confirm_null = match.arg(as.character(opts$confirm_null)[1], c("generic", "matched"))
  for (k in c("min_excess", "max_contra", "bal_frac", "min_bal_len")) {
    v = opts[[k]]
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0)
      stop("caller_tumour: opts$", k, " (segment_clone() rule) must be one number >= 0, not ", paste(format(v), collapse = ", "), call. = FALSE)
  }
  # ranges of the other numeric options: `what` names the allowed range in the error ("one <noun> <what>")
  chk = function(keys, ok, what, n = 1, noun = "number") for (k in keys) {
    v = opts[[k]]
    good = is.numeric(v) && (is.null(n) || length(v) == n) && length(v) >= 1 && !anyNA(v) && all(ok(v))
    if (!good) stop("caller_tumour: opts$", k, " must be ", if (is.null(n)) paste0(noun, "s ") else paste0("one ", noun, " "), what,
                    ", not ", if (length(v)) paste(format(v), collapse = ", ") else "empty", call. = FALSE)
  }
  whole = function(v) is.finite(v) & v == round(v)
  chk(c("q", "p_switch", "trunk_post", "min_post", "ref_post", "max_mix_prior"), function(v) v > 0 & v < 1,
      "strictly between 0 and 1")
  chk("mix_w", function(v) v > 0 & v < 1, "strictly between 0 and 1", n = NULL)
  chk("mix_prior", function(v) v >= 0 & v < 1, "in [0, 1)")
  chk(c("multiplet_prob", "deep_q", "member_post"), function(v) v >= 0 & v <= 1, "in [0, 1]")   # member_post: posterior >= it
  chk(c("em_rounds", "min_size", "min_n", "n_pc", "k", "n_perm", "split_max", "split_cand", "min_win", "min_ref"),
      function(v) whole(v) & v >= 1, ">= 1", noun = "whole number")
  # caps on a count: Inf = no cap (no subsampling / no limit on the leaves of the over-cut)
  chk(c("max_leaves", "max_leaves_trunk", "max_cells"), function(v) (whole(v) & v >= 1) | v == Inf, ">= 1 or Inf",
      noun = "whole number")
  chk("small_clone", function(v) whole(v) & v >= 0, ">= 0", noun = "whole number")   # 0: no clone is small
  chk(c("cap", "size_small", "split_min_len", "n_nb_trunk"), function(v) v > 0, "> 0")
  chk(c("trunk_dist", "control_min_len"), function(v) v >= 0, ">= 0")   # 0: every clone its own trunk / every piece a control
  chk(c("depth_delta", "disp_min", "imb_min", "depth_w", "split_sep", "split_min_df", "split_min_excess", "min_llr", "mix_margin",
        "mix_rel", "mix_min_dist", "same_dist", "doublet_delta", "normal_min_event"), function(v) v >= 0, ">= 0")
  for (k in c("split_depth", "repass", "normal_hypothesis")) {
    v = opts[[k]]
    if (!is.logical(v) || length(v) != 1 || is.na(v))
      stop("caller_tumour: opts$", k, " must be TRUE or FALSE, not ", paste(format(v), collapse = ", "), call. = FALSE)
  }
  if (!is.character(opts$classes) || !length(opts$classes) || !all(opts$classes %in% region_classes()))
    stop("caller_tumour: opts$classes must be region classes (", paste(region_classes(), collapse = ", "), "), not ",
         paste(format(opts$classes), collapse = ", "), call. = FALSE)
  if (!is.null(opts$allelic_exclude) && !is.character(opts$allelic_exclude))
    stop("caller_tumour: opts$allelic_exclude must be chromosome names", call. = FALSE)
  opts
}

# run notes on degenerate nulls of the tumour mode (the calls are not changed): `ref_thr` = the reference-pool
# cleaning threshold of clean_reference() (Inf: no control rows, nothing removed); `n_null` = null rows per
# confirmation hypothesis (named; all 0: no null rows at all, every threshold Inf); a hypothesis with fewer than
# 10 / (1 - q) null rows has a poorly determined type-7 `q` quantile (below 1 / (1 - q) rows: about the largest null
# score, which falls short of the true quantile)
confirm_null_notes = function(ref_thr = NULL, n_null = NULL, q = 0.999) {
  out = character()
  if (length(ref_thr) == 1 && !is.na(ref_thr) && is.infinite(ref_thr))
    out = c(out, paste("reference pool cleaning: no control rows to calibrate the allelic threshold (Inf), so no cell was",
                       "removed from the reference pool"))
  if (length(n_null)) {
    if (all(n_null == 0)) {
      out = c(out, paste("per-cell confirmation: no null rows (no reference-pool cell and no control piece): every",
                         "threshold is Inf and no cell call is confirmed"))
    } else {
      n_min = ceiling(10 / (1 - q))
      few = n_null[n_null < n_min]
      if (length(few))
        out = c(out, sprintf(paste("per-cell confirmation: %d of %d hypotheses have fewer than %d null rows (10 / (1 - q), q %g): %s;",
                                   "their thresholds rest on few null rows (below 1 / (1 - q) rows the quantile is about the",
                                   "largest null score and falls short of the true one, so false confirmations there can exceed",
                                   "1 - q; Inf without null rows)"), length(few), length(n_null), n_min, q,
                             paste(sprintf("%s %d", names(few), as.integer(few)), collapse = ", ")))
    }
  }
  out
}

#' Input null of the tumour mode: depth t null and allelic rho
#'
#' The depth Student-t null ([fit_depth_null()], completed by class) and the beta-binomial rho
#' ([fit_allelic_null()] on the arm rows with votes) that [caller_tumour()] uses for discovery, assignment and
#' per-cell confirmation, fitted on the null rows (`null_region` and `eligible`) of a set of cells:
#'   * `"flagged"`: the cells with `null_cell` TRUE (the `flat` column of the cell metadata, or given null cells),
#'     whenever there is one and `select` is FALSE;
#'   * `"selected"`: otherwise, the copy-number-neutral cells that [select_null_cells()] finds on the null rows
#'     (chromosome and arm rows, the tumour state set, `n_iter` iterations from all cells minus `exclude`): the cells
#'     kept by both of its last two iterations. They are used when there are at least `min_cells` of them and the last
#'     iteration is stable (Jaccard index of the sets it started from and ended with >= `min_jaccard`). In a tumour
#'     sample the cells without a call on the null rows are the normal cells (and tumour cells whose events all lie
#'     outside the null rows); without such cells the selection does not settle;
#'   * `"all"`: otherwise all cells (no neutral cells found). The fit then absorbs the tumour's events on the null
#'     rows (inflated rho, heavy depth tails), as every tumour-mode run without flagged cells did before.
#'
#' A fit on all cells is wrong whenever tumour cells carry events on the null rows: rho and the t tails take up the
#' events, which weakens assignment and per-cell confirmation. Flagged copy-number-normal cells (e.g. from a matched
#' assay) remain the better null: the selection can only find cells that look neutral at the data's resolution. The
#' note gives the rho of all cells next to the fitted one: a selected null whose rho is close to the all-cell rho is
#' a set of tumour cells that carry the same events (in a sample without normal cells).
#' @param long features long table ([cell_features()] `long`: cell, region, class, sd_class, chrom, stratum, z,
#'   sd_null, nA, n_hap, null_region, eligible, null_cell)
#' @param params [caller_defaults()]
#' @param exclude cells never selected (e.g. multiplets)
#' @param select TRUE: select null cells even when some are flagged
#' @param n_iter iterations of the selection (>= 2)
#' @param min_cells minimal number of selected cells
#' @param min_jaccard minimal Jaccard index of the selection's last iteration
#' @return NULL without null rows (of the flagged cells); else list(tpar, rho, source = "flagged" / "selected" /
#'   "all", null_cells (the cells fitted), selection (the [select_null_cells()] log; NULL for flagged cells), note
#'   (one line for the run notes))
#' @export
tumour_null_fit = function(long, params = caller_defaults(), exclude = NULL, select = FALSE, n_iter = 5, min_cells = 5,
                           min_jaccard = 0.7) {
  long = as.data.table(long)
  nr = !is.na(long$null_region) & long$null_region & !is.na(long$eligible) & long$eligible
  flagged = "null_cell" %in% names(long) && any(long$null_cell %in% TRUE)
  sel = NULL
  if (flagged && !select) {
    src = "flagged"
    nr = nr & long$null_cell %in% TRUE
    if (!any(nr)) return(NULL)
    note = sprintf("tumour input null (depth t, allelic rho): fitted on the %d flagged null cells", uniqueN(long$cell[nr]))
  } else {
    if (!any(nr)) return(NULL)
    what = if (flagged) "null-cell selection asked for" else "no flagged null cells"
    sl = select_null_cells(long[nr], states = caller_states("tumour", params), params = params, n_iter = n_iter,
                           exclude = exclude)
    sel = sl$log
    keep = intersect(sl$null_cells, sl$previous)
    jac = length(keep) / length(union(sl$null_cells, sl$previous))
    why = sprintf("%d cells without a call on the null rows in the last two iterations, Jaccard %.2f", length(keep), jac)
    if (length(keep) >= min_cells && is.finite(jac) && jac >= min_jaccard) {
      src = "selected"
      nr = nr & long$cell %in% keep
      note = sprintf("tumour input null (depth t, allelic rho): %s; fitted on %d copy-number-neutral cells selected on the null rows (%s)",
                     what, length(keep), why)
    } else {
      src = "all"
      note = sprintf(paste("tumour input null (depth t, allelic rho): %s; no stable set of at least %d copy-number-neutral cells",
                           "(%s), so fitted on all cells: rho and the depth tails absorb the tumour's events on the null rows"),
                     what, min_cells, why)
    }
  }
  tpar = complete_tpar(long, fit_depth_null(long[nr]))
  arm = nr & long$class == "arm" & long$n_hap > 0
  rho = if (any(arm, na.rm = TRUE)) fit_allelic_null(long$nA[arm %in% TRUE], long$n_hap[arm %in% TRUE]) else NA_real_
  note = paste0(note, sprintf("; rho %s", signif(rho, 3)))
  if (src == "selected") {
    a0 = !is.na(long$null_region) & long$null_region & !is.na(long$eligible) & long$eligible & long$class == "arm" & long$n_hap > 0
    note = paste0(note, sprintf(" (all cells: %s)", signif(fit_allelic_null(long$nA[a0 %in% TRUE], long$n_hap[a0 %in% TRUE]), 3)))
  }
  list(tpar = tpar, rho = rho, source = src, null_cells = unique(long$cell[nr]), selection = sel, note = note)
}

# per-cell call scores on a set of regions for cells of clone k (depth + allele under clone k's state)
cell_state_scores = function(inp, regs, rows, k_members, tpar, rho, params, ref_cells = NULL) {
  # rows: data.table cell, region, cn, ploidy, f_state, c, a, oclone (orientation clone)
  d = region_features(inp, regs, cells = unique(rows$cell))
  if (length(ref_cells) >= 10) {
    rf = region_features(inp, regs, cells = ref_cells)
    rm_ = rf[, .(rmean = mean(log2ratio, na.rm = TRUE)), by = region]
    d[rm_, log2ratio := log2ratio - i.rmean, on = "region"]
    d[, z := log2ratio / pmax(sd_null, params$sd_floor)]
  }
  d = merge(rows, d, by = c("cell", "region"))
  tp = cl_match_tpar(d, tpar)
  d[, mu := c * log2(pmax(cn, 0.2) / if (length(ref_cells) >= 10) ploidy else 2) / pmax(sd_null, params$sd_floor)]
  d[, s_depth := ll_depth_shift(z, mu, tp)]
  cells = inp$cells$cell
  d[, s_allele := 0]
  rr = regs[region %in% d$region]
  tv = trip_votes(inp$votes, inp$atoms, rr, length(cells), inp$allelic_chroms)
  if (nrow(tv$trip)) for (k in unique(d$oclone)) {
    oc = oriented_counts(tv, match(k_members[clone == k, cell], cells))
    oc[, `:=`(cell = cells[ci], region = rr$region[ri])]
    d[oc, `:=`(oA = i.oA, oN = i.oN), on = .(cell, region)]
    sel = d$oclone == k & !is.na(d$oN)
    if (any(sel)) d[sel, `:=`(s_allele = fifelse(abs(f_state - 0.5) > 1e-9,
                                                 ll_allelic(oA, oN, 0.5 + a * (f_state - 0.5), rho, params$eps, two_sided = FALSE, per_row = TRUE), 0),
                              nA_o = oA, n_o = oN)]
    d[, c("oA", "oN") := NULL]
  }
  if (!"nA_o" %in% names(d)) d[, `:=`(nA_o = 0, n_o = 0)]
  d[is.na(nA_o), `:=`(nA_o = 0, n_o = 0)]
  d[, base_ := if (length(ref_cells) >= 10) ploidy else 2L]
  d[, score := s_depth + s_allele]
  d[]
}
utils::globalVariables(c("s_depth", "s_allele", "i.oA", "i.oN", "oclone"))

b_fun = function(f) get(f, envir = asNamespace("scphasedcn"))

# score rows (cell_state_scores output) with the scorer: one non-neutral state per row (the clone state, shift
# c * log2(cn / base), major fraction 0.5 + a (f_state - 0.5)) or, for null rows, a given state set (B's tumour set)
b_score = function(rows, tpar, rho, params, states = NULL) {
  r = copy(rows)
  r[, `:=`(class = sd_class, nA = nA_o, n_hap = n_o)]
  tp = if (is.null(tpar)) NULL else b_fun("complete_tpar")(r, tpar)
  dm = if (is.null(tp)) "gauss" else "t"
  if (!is.null(states)) return(b_fun("score_cells")(r, states, tpar = tp, rho = rho, depth_model = dm, params = params))
  r[, `:=`(shift_ = c * log2(pmax(cn, 0.2) / base_), f_ = 0.5 + a * (f_state - 0.5))]
  r[, grp_ := .GRP, by = .(label, shift_, f_)]
  out = rbindlist(lapply(split(r, by = "grp_"), function(x) {
    b_score_hyp(x, x$label[1], x$shift_[1], x$f_[1], tp, rho, params)[, .(cell, region, s_comb)]
  }))
  out
}
utils::globalVariables(c("nA_o", "n_o", "shift_", "f_", "grp_", "base_", "s_comb", "thr_s_comb"))

# one clone hypothesis scored with the scorer's score_cells under the two-state set {neutral, label at depth shift
# `shift` and major fraction `f`}; x carries class (= sd_class), stratum, z, sd_null, nA, n_hap (oriented votes);
# tpar = completed t fits (NULL: Gaussian depth). Returns score_cells' table (s_comb = the hypothesis LLR).
b_score_hyp = function(x, label, shift, f, tpar, rho, params, cn = NULL) {
  cn = as.integer(cn %||% if ("cn" %in% names(x)) x$cn[1] else 3L)
  st = data.table(state = c("neutral", label), row_id = c("neutral_1", paste0(label, "_1")),
                  nA = c(1L, cn), nB = c(1L, 0L), cn = c(2L, cn), shift = c(0, shift), f = c(0.5, f))
  b_fun("score_cells")(x, st, tpar = tpar, rho = rho, depth_model = if (is.null(tpar)) "gauss" else "t", params = params)
}

# key of a clone hypothesis (matched confirmation null): clone | label | shift | f
hyp_key = function(clone, label, shift, f) paste(clone, label, sprintf("%.6g", shift), sprintf("%.6g", f), sep = "|")

# matched confirmation null: null rows (cell_state_scores output, one copy per event clone with oclone = that clone,
# `kind` pool / control) scored per hypothesis of the cell rows `sc` (clone, label, cn, c, a, f_state, base_) under
# the same two-state set; returns list(sc = sc with hyp, nb = null rows cell, region, stratum, kind, hyp, s_comb)
b_matched_null = function(sc, nrows, tpar, rho, params) {
  sc[, hyp := hyp_key(clone, label, c * log2(pmax(cn, 0.2) / base_), 0.5 + a * (f_state - 0.5))]
  hy = unique(sc[, .(hyp, clone, label, cn, shift = c * log2(pmax(cn, 0.2) / base_), f = 0.5 + a * (f_state - 0.5))], by = "hyp")
  r = copy(nrows)
  r[, `:=`(class = sd_class, nA = nA_o, n_hap = n_o)]
  tp = if (is.null(tpar)) NULL else b_fun("complete_tpar")(r, tpar)
  by_k = split(r, by = "oclone")
  nb = rbindlist(lapply(seq_len(nrow(hy)), function(i) {
    x = by_k[[hy$clone[i]]]
    if (is.null(x) || !nrow(x)) return(NULL)
    b_score_hyp(x, hy$label[i], hy$shift[i], hy$f[i], tp, rho, params, cn = hy$cn[i])[, .(cell, region, stratum, kind, hyp = hy$hyp[i], s_comb)]
  }))
  list(sc = sc, nb = nb, hyps = hy)
}
utils::globalVariables(c("hyp", "kind", "call_s_comb", "n_called"))

#' Tumour-mode copy-number caller
#'
#' Full loop: [split_votes()] (discovery / test halves) -> clone feature matrix (discovery half, re-oriented by the
#' all-cell pseudobulk) -> [discover_clones()] -> [merge_clones()] (test half) -> global EM round -> within-trunk
#' discovery ([discover_within()]) -> [clone_profiles()] -> [fit_clone_states()] -> [segment_clone()] -> doublet clones ->
#' [assign_cells()] (`em_rounds` hard-EM rounds: the clone states are re-profiled on the members of the previous round,
#' the singlets with posterior >= `member_post`; clones left with fewer than `min_size` members leave the member set; the
#' assignment's `clone` is every cell's likelihood-best clone, `confident` = posterior >= `min_post` and not
#' doublet-like) ->
#' [clean_reference()] -> optional depth re-pass relative to the clean pool -> per-cell calls. Each cell
#' inherits its clone's non-neutral segments (`source = "inherited"`) and is `confirmed` when its own
#' score under the clone state (depth + allele LLR vs neutral) is above the threshold of the tumour null:
#' reference-pool cells x all pieces + all cells x control pieces, scored as the best of the loss / gain /
#' CN-LOH / amp states (type-7 `null_q` quantile per stratum, strict `>`; pooled over strata below
#' `min_null_rows`).
#'
#' Input `feat` (the feature step's `scn_features`, or an equivalent list): `features` (long table),
#' `regions`, `votes`, `atoms`, `cells` (data.table with `cell`, in vote index order; `n_frag` / `stratum`),
#' `win` (1 Mb window depth: list(windows, x, e) with windows x cells matrices, or a long table cell, chrom,
#' start, end, x, e) and `null_fit` (list with `tpar` = [fit_depth_null()] table and `rho`; without it, or without its
#' `rho`, the input null is [tumour_null_fit()] of the features, as in [call_scn()]).
#' @param feat scn_features list
#' @param params [caller_defaults()]
#' @param multiplets optional data.table cell, prob ([multiplet_table()]); cells above `opts$multiplet_prob` are excluded from
#'   discovery and pseudobulks but assigned and flagged
#' @param control_regions optional user control regions (regions table) replacing the automatic ones (of the
#'   confirmation null and of the reference cleaning, [clean_reference()])
#' @param opts tumour-mode options ([clone_defaults()])
#' @param verbose log progress
#' @return `scphasedcn` (mode tumour) with clones, assignments, segments, profiles, calls (plus `extra`:
#'   intermediate results in `null_fit$clones`)
#' @export
caller_tumour = function(feat, params = caller_defaults(), multiplets = NULL, control_regions = NULL, opts = clone_defaults(),
                         verbose = TRUE) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  say = function(...) if (verbose) message(format(Sys.time(), "%H:%M:%S"), " [tumour] ", ...)
  opts = check_clone_opts(opts)
  inp = clone_input(feat)
  if (!nrow(inp$votes) || length(inp$rho) != 1 || !is.finite(inp$rho))
    stop("caller_tumour: no allele votes (no molecule votes, or none on the null rows to fit rho); the tumour mode needs ",
         "read-phased allele votes", call. = FALSE)
  cells = inp$cells$cell
  notes = character()
  if (length(opts$allelic_exclude)) {
    # female chrX: Tn5 fragments come mostly from the active X, so X-inactivation patches look like opposite
    # haplotype imbalance; chrX / chrY are depth-only in the clone steps
    ex = inp$atoms[chrom %in% opts$allelic_exclude, atom]
    inp$votes = inp$votes[!atom %in% ex]
    inp$allelic_chroms = setdiff(inp$allelic_chroms, opts$allelic_exclude)
  }
  mult_prob = NULL; mult = character()
  if (!is.null(multiplets)) {
    md = as.data.table(multiplets)
    mult_prob = stats::setNames(pmin(pmax(as.numeric(md$prob), 0), 1), as.character(md$cell))
    mult_prob[is.na(mult_prob)] = 0
    mult = names(mult_prob)[mult_prob > opts$multiplet_prob]
  }
  disc_cells = setdiff(intersect(cells, unique(inp$features$cell)), mult)
  say(length(cells), " cells, ", length(disc_cells), " used for discovery")
  # all-cell oriented votes per arm (trunk distance, clone_allelic_distance); computed once
  arms = inp$regions[class == "arm" & eligible == TRUE & allelic_ok == TRUE & chrom %in% inp$allelic_chroms]
  arm_oc = NULL
  if (nrow(arms) && nrow(inp$votes)) {
    arm_oc = cl_orient(inp$votes, arms, inp$atoms, NULL, loo = TRUE)
    arm_oc[, cell := cells[cell]]
  }
  nfr = if ("n_frag" %in% names(inp$cells)) stats::setNames(as.numeric(inp$cells$n_frag), cells) else {
    u = unique(inp$features[, .(cell, n_frag)], by = "cell"); stats::setNames(as.numeric(u$n_frag), u$cell)
  }

  # 1. discovery and merging
  # discovery on one half of the molecule votes (and of the window depth), the merge / split tests on the other
  # (split_votes, split_win): leaves formed on a difference are never tested on the same molecules / fragments
  vs = split_votes(inp$votes, 0.5, params$seed)
  inp_disc = inp; inp_disc$votes = vs$disc
  ws = if (isTRUE(opts$split_depth) && !is.null(inp$win)) split_win(inp$win, 0.5, params$seed) else NULL
  feat_disc = if (!is.null(ws)) half_depth_features(inp, ws$disc, opts$classes) else inp$features
  feat_test = if (!is.null(ws)) half_depth_features(inp, ws$test, opts$classes) else inp$features
  # split pass mode: FALSE = none, "disagree" = split_disagree() (which needs the window and vote halves themselves:
  # proposals on the discovery halves, test and re-fit on the test halves)
  split_mode = if (identical(opts$split_pass, "disagree")) "disagree" else "none"
  halves = if (split_mode == "disagree") list(win_disc = ws$disc, win_test = ws$test) else NULL
  if (split_mode == "disagree" && is.null(ws)) notes = c(notes, "split_pass disagree without split_depth window halves: proposals use the full window depth")
  rm(ws); invisible(gc(verbose = FALSE))
  cd = contrast_data(feat_test, vs$test, inp$atoms, inp$regions, cells, opts$classes, inp$allelic_chroms)
  rm(feat_test); invisible(gc(verbose = FALSE))
  # development checkpoint (opts$checkpoint = an .rds path): the global discovery and the first round are saved
  # there and re-used by later runs with the same input and discovery options
  ck = if (length(opts$checkpoint) && nzchar(opts$checkpoint) && file.exists(opts$checkpoint)) readRDS(opts$checkpoint) else NULL
  if (!is.null(ck)) {
    disc = ck$disc; mc = ck$mc; X = NULL
    say("checkpoint ", opts$checkpoint, ": ", length(unique(disc$leaves)), " leaves, ", length(unique(mc$clones)), " clones after merging")
  } else {
    X = clone_feature_matrix(feat_disc[cell %in% disc_cells], opts$classes, inp$rho, opts$cap, opts$min_n,
                             vs$disc, inp$atoms, inp$regions, cells,
                             keep_orient = character())
    disc = discover_clones(X, opts$n_pc, opts$k, opts$max_leaves, opts$min_size, params$seed)
    disc$pc = NULL
    say(length(unique(disc$leaves)), " leaves")
    mc = merge_clones(NULL, NULL, inp$atoms, inp$regions, disc$leaves, disc, q = opts$q, n_perm = opts$n_perm,
                      seed = params$seed, rho = inp$rho, classes = opts$classes, cells = cells, max_cells = opts$max_cells, cdata = cd,
                      depth_delta = opts$depth_delta)
    say(length(unique(mc$clones)), " clones after merging (", nrow(mc$tests), " split tests)")
  }
  members = data.table(cell = names(mc$clones), clone = paste0("G", mc$clones))
  rm(X); invisible(gc(verbose = FALSE))

  grid = cn_grid()
  fit_states = function(profile, grid) fit_clone_states(profile, grid, p_switch = opts$p_switch, min_win = opts$min_win)
  ref_pool = character()
  run_states = function(members, ref_cells) {
    sz = members[, .N, by = clone]
    size = stats::setNames(ifelse(sz$N < opts$small_clone, opts$size_small, 1e6), sz$clone)
    # a single clone is never its own depth reference (its profile would be flat by construction): median-centred
    if (nrow(sz) < 2) ref_cells = character()
    else if (length(ref_cells) < opts$min_ref) {
      # provisional reference before the clean pool exists: an allelically event-free clone of the median-centred
      # fit (no imbalanced non-neutral segment >= 10 Mb; the allele channel is unaffected by the normalisation's
      # shared window bias, which the reference removes from the depth); none -> median-centred profiles
      ref_cells = character()
      p0 = clone_profiles(inp$win, inp$votes, inp$atoms, members, size = size, seed = params$seed, cells = cells,
                          rho = inp$rho, allelic_chroms = inp$allelic_chroms, haploid = inp$haploid)
      f0 = fit_states(p0, grid)
      s0 = segment_clone(f0, opts$p_switch, opts$min_win, grid, opts$min_llr, opts$min_excess, opts$bal_frac, opts$min_bal_len,
                        opts$max_contra)
      dirty = s0[label != "neutral" & nA_of(state) != nB_of(state) & (end - start) >= 1e7, unique(clone)]
      bal = p0[n_hap > 0 & ploidy == 2, .(imb = mean((nA - n_hap / 2)^2 / (n_hap / 4))), by = clone]
      bal = bal[clone %in% sz[N >= opts$min_ref, clone] & !clone %in% dirty]
      if (nrow(bal)) {
        k0 = bal[which.min(imb), clone]
        ref_cells = members[clone == k0, cell]
        notes <<- c(notes, sprintf("provisional depth reference: clone %s (allelically event-free; mean window chi2 %.2f)", k0, min(bal$imb)))
      } else {
        notes <<- c(notes, "no allelically event-free clone: median-centred clone profiles")
        return(list(profiles = p0, fit = f0, segments = s0, ref = character()))
      }
    }
    prof = clone_profiles(inp$win, inp$votes, inp$atoms, members, ref_cells = ref_cells, size = size, seed = params$seed,
                          cells = cells, rho = inp$rho, allelic_chroms = inp$allelic_chroms, haploid = inp$haploid, min_ref = opts$min_ref)
    fit = fit_states(prof, grid)
    seg = segment_clone(fit, opts$p_switch, opts$min_win, grid, opts$min_llr, opts$min_excess, opts$bal_frac, opts$min_bal_len,
                       opts$max_contra)
    list(profiles = prof, fit = fit, segments = seg, ref = if (length(ref_cells) >= opts$min_ref) ref_cells else character())
  }
  doublet_clones = character()
  # one hard-EM round: clone states on the members -> doublet clones -> trunks / mixture pairs -> likelihood assignment
  # of every cell (clone = the likelihood-best clone; `confident` = posterior >= min_post and not doublet-like)
  trunks_of = function(members, K_) {
    ad = if (!is.null(arm_oc) && length(K_) >= 2) clone_allelic_distance(arm_oc, members[clone %in% K_], arms) else NULL
    tr = if (is.null(ad)) stats::setNames(rep("T1", length(K_)), K_) else
      stats::setNames(paste0("T", stats::cutree(stats::hclust(stats::as.dist(ad), method = "average"), h = opts$trunk_dist)), rownames(ad))
    list(tr = tr, ad = ad)
  }
  em_round = function(members, ref_pool, label) {
    st = run_states(members, ref_pool)
    tk = trunks_of(members, sort(unique(st$segments$clone)))
    # doublet clones: mixtures of two clones of different trunks only
    xp = outer(tk$tr, tk$tr, "!=") & upper.tri(diag(length(tk$tr)))
    ij = which(xp, arr.ind = TRUE)
    dc = flag_doublet_clones(members[, .(n_cells = .N), by = clone], st$segments, st$profiles, st$fit, opts$doublet_delta,
                             pairs = cbind(names(tk$tr)[ij[, 1]], names(tk$tr)[ij[, 2]]))
    if (any(dc$doublet)) {
      doublet_clones <<- union(doublet_clones, dc[doublet == TRUE, clone])
      say("doublet clone(s): ", paste(dc[doublet == TRUE, paste0(clone, " (", best_pair, ")")], collapse = ", "))
      members = members[!clone %in% doublet_clones]
      st$segments = st$segments[!clone %in% doublet_clones]
      st$profiles = st$profiles[!clone %in% doublet_clones]
      st$fit$par = st$fit$par[!clone %in% doublet_clones]
      tk = trunks_of(members, sort(unique(st$segments$clone)))
    }
    if (isTRUE(opts$normal_hypothesis) && length(unique(st$segments$clone)) >= 1) {
      # the diploid normal state is always a hypothesis in a tumour (normal admixture): a fixed all-neutral clone N0
      # replaces every fitted clone without a non-neutral segment >= normal_min_event (a small fitted normal clone drifts:
      # low-depth tumour cells join it, it acquires events and the reference pool is lost)
      nn = add_normal_clone(st, opts$normal_min_event)
      if (length(nn$replaced)) members = copy(members)[clone %in% nn$replaced, clone := "N0"]
      st = nn$st
      tk = trunks_of(members[clone != "N0"], setdiff(sort(unique(st$segments$clone)), "N0"))
      tk$tr = c(tk$tr, N0 = "TN")
    }
    tr = tk$tr; ad = tk$ad
    dm = clone_state_distance(st$segments)
    # mixture hypotheses: clones of different trunks whose states differ allelically on >= mix_min_dist
    ij = which(dm >= opts$mix_min_dist & upper.tri(dm) & outer(tr[rownames(dm)], tr[colnames(dm)], "!="), arr.ind = TRUE)
    mp = cbind(rownames(dm)[ij[, 1]], colnames(dm)[ij[, 2]])
    pieces = segment_pieces(st$segments, inp$regions)
    fp = list(features = region_features(inp, pieces), pieces = pieces, votes = inp$votes, atoms = inp$atoms, cells = cells)
    asg = assign_cells(fp, members, st$segments, inp$tpar, inp$rho, par = st$fit$par, mixtures = TRUE, min_post = opts$min_post,
                       mix_prior = opts$mix_prior, ref_cells = st$ref, mix_w = opts$mix_w, mix_margin = opts$mix_margin, mix_rel = opts$mix_rel,
                       params = params, mix_pairs = mp, mult_prior = mult_prob, trunks = tr, max_mix_prior = opts$max_mix_prior)
    asg[, multiplet := cell %in% mult]
    say(label, ": ", uniqueN(st$segments$clone), " clones in ", uniqueN(tr), " trunk(s), ", nrow(mp), " mixture pairs; confident ",
        sum(asg$confident), " / ", nrow(asg), "; trunk posterior >= ", opts$trunk_post, " ", sum(asg$p_trunk >= opts$trunk_post),
        "; doublet-like ", sum(asg$doublet_like))
    list(st = st, asg = asg, trunks = tr, dist = dm, trunk_dist = ad)
  }
  next_members = function(asg) {
    m = asg[doublet_like == FALSE & multiplet == FALSE & posterior >= opts$member_post, .(cell, clone)]
    small = m[, .N, by = clone][N < opts$min_size, clone]
    if (length(small)) { m = m[!clone %in% small]; notes <<- c(notes, paste("dropped small clone(s):", paste(small, collapse = ", "))) }
    m
  }
  # clones of one trunk whose segment states differ on < same_dist bp (any state difference) are merged
  merge_same = function(r, members) {
    dm = clone_state_distance(r$st$segments, min_df = 0.01, min_dcn = 1)
    K = rownames(dm)
    if (length(K) < 2) return(list(members = members, n = 0L))
    g = stats::cutree(stats::hclust(stats::as.dist(dm), method = "single"), h = max(opts$same_dist - 1, 0))
    key = paste(r$trunks[K], g)
    rep_ = K[match(key, key)]
    n = sum(rep_ != K)
    if (n) {
      say("merging clones with the same states: ", paste(sprintf("%s -> %s", K[rep_ != K], rep_[rep_ != K]), collapse = ", "))
      members = copy(members)[, clone := rep_[match(clone, K)]]
    }
    list(members = members, n = n)
  }
  trunk_info = NULL
  if (!is.null(ck$r)) {
    r = ck$r; doublet_clones = ck$doublet_clones
    say("checkpoint: global round re-used (", uniqueN(r$st$segments$clone), " clones)")
  } else r = em_round(members, ref_pool, "global round")
  if (length(opts$checkpoint) && nzchar(opts$checkpoint) && is.null(ck))
    saveRDS(list(disc = disc, mc = mc, r = r, doublet_clones = doublet_clones), opts$checkpoint)
  # 2. within-trunk discovery: blocks oriented by each trunk's pseudobulk, allelically selected regions
  tc = r$asg[doublet_like == FALSE & multiplet == FALSE & p_trunk >= opts$trunk_post, .(cell, trunk)]
  newm = list(); trunk_info = list()
  for (t_ in sort(unique(tc$trunk))) {
    tcells = tc[trunk == t_, cell]
    if (length(tcells) < 4 * opts$min_size) {
      newm[[t_]] = data.table(cell = tcells, clone = t_)
      trunk_info[[t_]] = data.table(trunk = t_, n_cells = length(tcells), n_regions = NA_integer_, n_leaves = 1L, n_clones = 1L)
      next
    }
    dw = discover_within(inp, tcells, vs$disc, feat_disc, cd, nfr, opts, params$seed)
    newm[[t_]] = data.table(cell = names(dw$clones), clone = paste0(t_, ".", dw$clones))
    trunk_info[[t_]] = data.table(trunk = t_, n_cells = length(tcells), n_regions = length(dw$regions), n_leaves = dw$n_leaves,
                                  n_clones = uniqueN(dw$clones), regions = paste(dw$regions, collapse = " "))
    say("trunk ", t_, " (", paste(names(r$trunks)[r$trunks == t_], collapse = ", "), "): ", length(tcells), " cells, ",
        length(dw$regions), " regions selected, ", dw$n_leaves, " leaves -> ", uniqueN(dw$clones), " clone(s)")
  }
  trunk_info = rbindlist(trunk_info, fill = TRUE)
  members = rbindlist(newm)
  lev = unique(members$clone)
  trunk_info[, clones := vapply(trunk, function(t_) paste(paste0("C", which(sub("\\..*", "", lev) == t_)), collapse = ","), "")]
  members[, clone := paste0("C", match(clone, lev))]
  if (split_mode == "disagree") { halves$votes_disc = vs$disc; halves$votes_test = vs$test }
  rm(vs, inp_disc); invisible(gc(verbose = FALSE))
  ms = list(n = 0L)
  for (round in seq_len(opts$em_rounds)) {
    r = em_round(members, ref_pool, paste("round", round))
    members = next_members(r$asg)
    ms = merge_same(r, members)
    members = ms$members
  }
  if (ms$n > 0) { r = em_round(members, ref_pool, "after merging"); members = next_members(r$asg) }
  split_tests = NULL; asg_pre_split = r$asg[, .(cell, clone, posterior, doublet_like, confident)]
  if (split_mode == "disagree") {
    # clones whose members disagree with the clone state on a segment (one copy apart on the same haplotype, or a
    # different major fraction): one pass, at most split_max splits (split_disagree)
    sp = split_disagree(members, r$st$segments, r$st$fit$par, inp, halves$win_disc, halves$votes_disc, cd,
                        win_test = halves$win_test, votes_test = halves$votes_test, ref_cells = r$st$ref, fit_states = fit_states,
                        opts = opts, params = params, seed = params$seed)
    split_tests = sp$tests
    n_cand = if (is.null(split_tests)) 0L else nrow(split_tests)
    say("disagreement split pass: ", n_cand, " candidate(s), ", if (n_cand) sum(!is.na(split_tests$differ)) else 0L, " tested, ",
        if (n_cand) sum(split_tests$differ %in% TRUE) else 0L, " different on the test halves, ",
        if (n_cand) sum(split_tests$split) else 0L, " split")
    if (n_cand && any(split_tests$split)) {
      say("disagreement split(s): ", paste(split_tests[split == TRUE, sprintf("%s %s %s -> %s (contradiction %.2f, %d / %d cells; re-fit %s / %s)",
                                                                             clone, region, state, alt, contra, n_a, n_b, state_a, state_b)], collapse = ", "))
      members = copy(sp$members)
      moved = sp$cells_b
      lev = unique(members$clone); members[, clone := paste0("C", match(clone, lev))]
      r = em_round(members, ref_pool, "after splitting"); members = next_members(r$asg)
      ms = merge_same(r, members); members = ms$members
      if (ms$n > 0) { r = em_round(members, ref_pool, "after merging"); members = next_members(r$asg) }
      # does each split survive the EM rounds and the same-state merge? (majority clone of the two parts afterwards)
      for (kb in names(moved)) {
        k0 = sub("b$", "", kb)
        ca_ = sp$members[clone == paste0(k0, "a"), cell]; cb_ = moved[[kb]]
        maj = function(cc) { x = r$asg[cell %in% cc & doublet_like == FALSE, .N, by = clone][order(-N)]
                             if (nrow(x)) sprintf("%s %.2f", x$clone[1], x$N[1] / sum(x$N)) else NA_character_ }
        ma = maj(ca_); mb = maj(cb_)
        split_tests[split == TRUE & clone == k0, `:=`(after_a = ma, after_b = mb,
                                                     separate = !is.na(ma) && !is.na(mb) && sub(" .*", "", ma) != sub(" .*", "", mb))]
      }
    }
  }
  rm(cd, feat_disc, halves); invisible(gc(verbose = FALSE))
  st = r$st; asg = r$asg
  # 3. clean reference and optional depth re-pass
  cr = clean_reference(inp, asg, st$segments, par = st$fit$par, members = members, min_post = opts$ref_post,
                       min_event = opts$normal_min_event,
                       control = if (is.null(control_regions)) NULL else as.data.table(control_regions)[, .(region, chrom, start, end)])
  ref_pool = cr$pool
  notes = c(notes, confirm_null_notes(ref_thr = cr$thr))
  say("reference pool: ", length(ref_pool), " cells (normal clone(s) ", paste(cr$normal, collapse = ", "), "; flagged ",
      sum(cr$table$flag %||% FALSE), ", thr ", signif(cr$thr, 3), ")")
  if (opts$repass && length(ref_pool) >= opts$min_ref && uniqueN(members$clone) >= 2) {
    r = em_round(members, ref_pool, "depth re-pass relative to the clean pool")
    st = r$st; asg = r$asg
    members = next_members(asg)
  }
  pool_cells = ref_pool
  asg[, ref_pool := cell %in% pool_cells]

  # confirmation null (opts$ref_null): the allelic rho refitted on the clean reference pool; assignment above keeps
  # the input fits (options validated by check_clone_opts)
  ref_null = opts$ref_null
  confirm_null = opts$confirm_null
  rho_c = inp$rho; tpar_c = inp$tpar; null_confirm = NULL
  refit = ref_null != "none" && length(ref_pool) >= opts$min_ref
  if (ref_null != "none" && !refit)
    notes = c(notes, sprintf("ref_null %s: reference pool %d < min_ref %d cells, input null kept", ref_null, length(ref_pool), opts$min_ref))
  rho_feat = NA_real_
  if (refit && all(c("null_region", "eligible", "class", "n_hap", "nA", "chrom") %in% names(inp$features))) {
    # diagnostic only: the call_scn recipe of the input rho (null, eligible, autosomal arm rows with votes) on the pool
    # cells. Its votes carry the feature step's orientation (e.g. by the cell's own group pseudobulk: self-orientation inflates
    # rho), not the leave-one-out clone orientation the confirmation scores, so the refit below uses the latter.
    fr = inp$features[cell %in% ref_pool & null_region %in% TRUE & eligible %in% TRUE & class == "arm" & n_hap > 0 &
                      !chrom %in% c("chrX", "chrY")]
    if (length(inp$allelic_chroms)) fr = fr[chrom %in% inp$allelic_chroms]
    rho_feat = if (nrow(fr)) fit_allelic_null(fr$nA, fr$n_hap) else NA_real_
  }
  matched = confirm_null == "matched"

  # 3. per-cell calls
  sg = st$segments
  ctl = if (!is.null(control_regions)) as.data.table(control_regions) else control_regions(sg, opts$control_min_len, inp$regions)
  if (!nrow(ctl)) notes = c(notes, "no automatic control regions (no piece neutral in every clone)")
  ev = sg[label != "neutral"]
  calls = NULL
  thr_tab = NULL
  if (nrow(ev)) {
    ev_reg = seg_regions(ev[, .(clone, chrom, start, end)], inp$regions)
    ev[, region := ev_reg$region]
    par = st$fit$par
    rows = asg[, .(cell, clone)][ev[, .(clone, region, cn, ploidy, f_state = fifelse(cn > 0, nA_of(state) / cn, 0.5), label)],
                                 on = "clone", allow.cartesian = TRUE, nomatch = NULL]
    rows[par, `:=`(c = i.c, a = i.a), on = "clone"]
    rows[, oclone := clone]
    # null: reference-pool cells x all pieces + all cells x control pieces, best of the generic states (generic) or one
    # copy per event clone, oriented by that clone's pseudobulk (matched)
    pieces = segment_pieces(sg, inp$regions)
    gen = data.table(cn = c(1L, 3L, 2L, 4L), f_state = c(1, 2 / 3, 1, 0.75))
    cmed = stats::median(par$c); amed = stats::median(par$a)
    ev_clones = intersect(unique(ev$clone), unique(members$clone))
    mk_null = function(cc, regs) {
      if (!length(cc) || !nrow(regs)) return(NULL)
      regs = regs[!chrom %in% inp$haploid]
      if (!nrow(regs)) return(NULL)
      base = CJ(cell = cc, region = regs$region)
      base[asg, oclone := i.clone, on = "cell"]
      base = base[!is.na(oclone) & oclone %in% members$clone]
      if (matched) {
        if (!nrow(base)) return(NULL)
        base[, oclone := NULL]
        return(rbindlist(lapply(ev_clones, function(k) {
          pk = par[clone == k]
          copy(base)[, `:=`(oclone = k, cn = 3L, ploidy = 2L, f_state = 2 / 3, c = if (nrow(pk)) pk$c[1] else cmed,
                            a = if (nrow(pk)) pk$a[1] else amed)]
        })))
      }
      rbindlist(lapply(seq_len(nrow(gen)), function(j) {
        r = copy(base)[, `:=`(cn = gen$cn[j], ploidy = 2L, f_state = gen$f_state[j], c = cmed, a = amed)]
        r[, ploidy := fifelse(regs$chrom[match(region, regs$region)] %in% inp$haploid, 1L, 2L)]
        r
      }))
    }
    nr = rbindlist(list(mk_null(ref_pool, pieces), mk_null(asg[multiplet == FALSE, cell], ctl)), fill = TRUE)
    nx = NULL; nsc_rows = data.table()
    if (!is.null(nr) && nrow(nr)) {
      nr[, kind := fifelse(region %in% ctl$region, "control", "pool")]
      regs_n = rbindlist(list(pieces, ctl), fill = TRUE)
      nx = cell_state_scores(inp, unique(regs_n, by = "region"), nr, members, inp$tpar, rho_c, params, st$ref)
      nsc_rows = unique(nx, by = c("cell", "region", if (matched) "oclone"))
    }
    if (refit && nrow(nsc_rows)) {
      # rho refitted on the votes the confirmation scores: the pool cells' piece votes, oriented leave-one-out by a clone
      # pseudobulk (one orientation per cell x piece: the cell's clone N0 with the generic null, the first event clone with
      # the matched null). In tumour samples this rho can depend on the orientation clone and be inflated on event
      # pieces (ambient tumour fragments in pool cells); a per-clone fit on control pieces is a possible refinement.
      pa = unique(nsc_rows[cell %in% ref_pool & n_o > 0], by = c("cell", "region"))
      rr = if (nrow(pa)) fit_allelic_null(pa$nA_o, pa$n_o) else NA_real_
      if (length(rr) == 1 && is.finite(rr)) {
        rho_c = rr
        nx[, s_allele := fifelse(abs(f_state - 0.5) > 1e-9, ll_allelic(nA_o, n_o, 0.5 + a * (f_state - 0.5), rho_c, params$eps, two_sided = FALSE, per_row = TRUE), 0)]
        nx[, score := s_depth + s_allele]
      }
      msg = sprintf(paste("confirmation rho refitted on the reference pool (%d cells, %d piece rows, clone-oriented votes): %s",
                          "(input %s; the feature-step recipe on the pool %s)%s"), length(ref_pool), nrow(pa), signif(rr, 3),
                    signif(inp$rho, 3), signif(rho_feat, 3), if (identical(rho_c, rr)) "" else "; not finite, input kept")
      say(msg); notes = c(notes, msg)
    }
    sc = cell_state_scores(inp, ev_reg, rows, members, inp$tpar, rho_c, params, st$ref)
    nsc = if (!is.null(nx)) nx[, .(score = max(score)), by = .(cell, region, stratum)] else
      data.table(cell = character(), region = character(), stratum = character(), score = numeric())
    rm(nx)
    use_b = nrow(nsc_rows) > 0             # the scorer; without null rows the type-7 quantile fallback below
    if (use_b && matched) {
      # the scorer, matched null: each clone hypothesis (clone | label | shift | f) against the null rows oriented by that
      # clone and scored under the same two-state set; one threshold per hypothesis x stratum (hypothesis-wide below
      # min_null_rows and for a stratum without null rows)
      bs = b_score(sc, tpar_c, rho_c, params)
      sc[bs, score := i.s_comb, on = .(cell, region)]
      mn = b_matched_null(sc, nsc_rows, tpar_c, rho_c, params)
      sc = mn$sc; nb = mn$nb
      nb[, class := hyp]
      nb_class = nb$class
      bt = b_fun("budget_thresholds")(nb, scores = "s_comb", by = c("class", "stratum"), cv = FALSE, params = params,
                                      null = rep(TRUE, nrow(nb)))
      sc[, class := hyp]
      sc = b_fun("apply_thresholds")(sc, bt, "s_comb")
      miss = is.na(sc$thr_s_comb)
      if (any(miss)) {
        bw = b_fun("budget_thresholds")(nb, scores = "s_comb", by = "class", cv = FALSE, params = params, null = rep(TRUE, nrow(nb)))
        sc[miss, thr_s_comb := bw$thr[match(class, bw$class)]]
      }
      sc[, thr := thr_s_comb]
      thr_tab = bt[, .(hyp = class, stratum, thr, n_null, pooled_fallback)]
      # held-out false-confirmation rate per hypothesis (fold CV over the null cells), control pieces vs pool rows
      nb[, fold := b_fun("assign_folds")(cell, seed = params$seed)]
      btc = b_fun("budget_thresholds")(nb, scores = "s_comb", by = c("class", "stratum"), cv = TRUE, params = params,
                                       null = rep(TRUE, nrow(nb)))
      nbc = b_fun("apply_thresholds")(nb, btc, "s_comb")
      rt = function(x) x[, .(n_null = sum(!is.na(call_s_comb)), n_called = sum(call_s_comb, na.rm = TRUE))]
      null_confirm = rbind(nbc[, rt(.SD), by = .(hyp, kind)], nbc[, rt(.SD), by = kind][, hyp := "all"],
                           rt(nbc)[, `:=`(hyp = "all", kind = "all")], use.names = TRUE)
      null_confirm[, per_1000 := 1000 * n_called / n_null]
      rm(nb, nbc, btc, mn)
      notes = c(notes, sprintf(paste("per-cell confirmation: the scorer (score_cells / budget_thresholds / apply_thresholds), matched null",
                                     "(%d clone hypotheses; held-out false confirmations %.2f / 1,000 on control pieces, %.2f on pool rows)"),
                               uniqueN(thr_tab$hyp), null_confirm[hyp == "all" & kind == "control", per_1000][1],
                               null_confirm[hyp == "all" & kind == "pool", per_1000][1]))
    } else if (use_b) {
      # the scorer: the cell's s_comb under its clone's state vs the tumour state set on the null rows
      bs = b_score(sc, tpar_c, rho_c, params)
      sc[bs, score := i.s_comb, on = .(cell, region)]
      nb = b_score(nsc_rows, tpar_c, rho_c, params, states = b_fun("caller_states")("tumour", params))
      nb[, class := "segment"]
      nb_n = nrow(nb)
      bt = b_fun("budget_thresholds")(nb, scores = "s_comb", by = c("class", "stratum"), cv = FALSE, params = params,
                                      null = rep(TRUE, nrow(nb)))
      sc[, class := "segment"]
      sc = b_fun("apply_thresholds")(sc, bt, "s_comb")
      sc[, thr := thr_s_comb]
      thr_tab = bt[, .(stratum, thr, n_null)]
      notes = c(notes, "per-cell confirmation: the scorer (score_cells / budget_thresholds / apply_thresholds), caller_states('tumour') null")
    }
    if (!use_b) thr_tab = {
      pooled = if (nrow(nsc)) unname(stats::quantile(nsc$score, params$null_q, type = 7)) else Inf
      tt = nsc[, .(thr = unname(stats::quantile(score, params$null_q, type = 7)), n_null = .N), by = stratum]
      tt[n_null < params$min_null_rows, thr := pooled]
      rbind(tt, data.table(stratum = "pooled", thr = pooled, n_null = nrow(nsc)))
    }
    if (!use_b) {
      sc[thr_tab, thr := i.thr, on = "stratum"]
      sc[is.na(thr), thr := thr_tab[stratum == "pooled", thr][1]]
    }
    # degenerate confirmation nulls (notes only): no null rows at all, or hypotheses with too few null rows
    n_null_h = if (use_b && matched) {
      hc = table(factor(nb_class, levels = unique(c(nb_class, sc$hyp))))
      stats::setNames(as.integer(hc), names(hc))
    } else if (use_b) c(segment = nb_n) else c(segment = nrow(nsc))
    notes = c(notes, confirm_null_notes(n_null = n_null_h, q = params$null_q))
    sc[is.na(thr), thr := Inf]
    calls = sc[, .(cell, region, class = "segment", chrom, start = ev$start[match(region, ev$region)],
                   end = ev$end[match(region, ev$region)], state = label, score, thr, called = score > thr, level = "cell",
                   source = fifelse(score > thr, "confirmed", "inherited"),
                   reportable = ev$confident[match(region, ev$region)])]
    # every cell has calls from its likelihood-best clone; `confident` = the cell's clone membership is confident
    calls[, confident := asg$confident[match(cell, asg$cell)] %in% TRUE]
    clone_rows = ev[, .(cell = clone, region, class = "segment", chrom, start, end, state = label, score = llr, thr = opts$min_llr,
                        called = llr > opts$min_llr, level = "clone", source = "cell", reportable = confident, confident = TRUE)]
    calls = rbind(clone_rows, calls)
    all_regions = rbind(inp$regions, pieces, ctl, ev_reg, fill = TRUE)
  } else all_regions = rbind(inp$regions, segment_pieces(sg, inp$regions), fill = TRUE)
  all_regions = unique(all_regions, by = "region")
  all_regions[is.na(sd_class), sd_class := class]
  clone_tab = members[, .(n_members = .N), by = clone]
  clone_tab = merge(asg[, .(n_cells = .N, n_assigned = sum(assigned)), by = clone], clone_tab, by = "clone", all = TRUE)
  clone_tab = merge(clone_tab, st$fit$par[, .(clone, c, a, rho, s_sys)], by = "clone", all.x = TRUE)
  clone_tab[, normal := clone %in% cr$normal]
  clone_tab[, trunk := r$trunks[clone]]
  clone_tab = rbind(clone_tab, data.table(clone = doublet_clones, doublet = TRUE), fill = TRUE)
  clone_tab[is.na(doublet), doublet := FALSE]
  thresholds = if (!is.null(thr_tab)) thr_tab[, .(score = "segment_state", class = "segment", stratum, fold = 0L,
                                                   q = params$null_q, thr, n_null = as.integer(n_null),
                                                   pooled_fallback = stratum == "pooled" | n_null < params$min_null_rows)] else NULL
  if (!is.null(thr_tab) && "hyp" %in% names(thr_tab)) thresholds[, `:=`(pooled_fallback = thr_tab$pooled_fallback, hyp = thr_tab$hyp)]
  # events: the called (confirmed) reportable calls consolidated per cell / clone, an empty table otherwise
  events = if (exists("consolidate_calls", envir = asNamespace("scphasedcn"), inherits = FALSE)) {
    b_fun("consolidate_calls")(if (is.null(calls)) data.table(cell = character(), region = character(), class = character(),
                                 chrom = character(), start = numeric(), end = numeric(), state = character(), score = numeric(),
                                 called = logical(), level = character(), reportable = logical()) else calls, all_regions)
  } else NULL
  new_scphasedcn("tumour", cells = inp$cells, regions = all_regions, features = NULL, thresholds = thresholds, calls = calls, events = events,
              clones = clone_tab, assignments = asg, segments = sg, profiles = st$profiles,
              null_fit = list(tpar = inp$tpar, rho = inp$rho, ref_thr = cr$thr, rho_confirm = rho_c, tpar_confirm = tpar_c,
                              rho_pool_features = rho_feat, ref_null = ref_null, confirm_null = null_confirm,
                              clones = list(leaves = disc$leaves, merged = mc$clones, tests = mc$tests, par = st$fit$par, reference = cr$table, doublet_clones = doublet_clones,
                                            trunks = r$trunks, trunk_discovery = trunk_info, state_distance = r$dist,
                                            trunk_distance = r$trunk_dist,
                                            split_tests = split_tests, assignments_pre_split = asg_pre_split)),
              params = c(params, list(tumour = opts)), notes = notes)
}
utils::globalVariables(c("pool_cells", "imb", "prob", "N", "i.c", "i.a", "i.clone", "i.thr", "n_null", "n_members", "normal",
                         "trunk", "p_trunk", "clones", "mix_prior"))
