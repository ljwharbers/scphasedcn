# clones_states.R
# What: tumour mode of the copy-number caller, part 2 of 4: clone pseudobulk profiles (1 Mb
#       depth log2 with SE, split-half oriented haplotype fraction) and clone copy-number states on an nA:nB
#       grid with a fitted depth contrast and orientation attenuation (allele-specific gains included),
#       depth / allele emissions, clone noise and Viterbi segments.
# Why: in tumours per-cell states are not
#      reliable (attenuated haplotype fractions, allele-specific amplicons, a heavy-tailed per-cell depth
#      term); clone-level depth and orientation settle the state.
# Inputs: 1 Mb window matrices, votes, atoms and a clone assignment (see [clone_profiles()]).
# Outputs: data.tables (profiles, state log-likelihoods, fit parameters, segments).
# See also: clones_discovery.R, clones_assign.R, clones_tumour.R.

# ---- clone profiles ----------------------------------------------------------------------------------

# aggregate the 1 Mb window matrices to `size` windows (window start = floor(start / size) * size)
win_aggregate = function(win, size) {
  wt = copy(win$windows)
  if (size <= 1e6 + 1) return(win)
  wt[, w5 := floor(start / size) * size]
  key_ = paste0(wt$chrom, ":", format(wt$w5, scientific = FALSE, trim = TRUE))
  ind = group_indicator(key_)
  nw = wt[, .(chrom = chrom[1], start = min(w5), end = max(end)), by = .(window = key_)]
  nw = nw[match(rownames(ind), window)]
  nw = nw[order(chrom_rank(chrom), start)]
  ind = ind[nw$window, , drop = FALSE]
  x = as.matrix(ind %*% win$x); e = as.matrix(ind %*% win$e)
  colnames(x) = colnames(win$x); colnames(e) = colnames(win$e)
  list(windows = nw, x = x, e = e)
}

row_sd = function(M) {
  n = rowSums(!is.na(M)); m = rowSums(M, na.rm = TRUE) / pmax(n, 1)
  sqrt(pmax((rowSums(M^2, na.rm = TRUE) - n * m^2) / pmax(n - 1, 1), 0))
}

#' Clone pseudobulk profiles
#'
#' Per clone x window (1 Mb, or `size`): depth log2 of the pooled counts (`log2((sum x + 0.5) / (sum e +
#' 0.5))`), relative to the clean reference pool's pseudobulk when `ref_cells` has >= `min_ref` cells, then
#' centred on the clone's autosomal median (states are relative to the clone baseline); SE = SD of the
#' per-cell window log2 ratios / sqrt(cells) (combined with the reference SE). Haplotype fraction: the
#' clone's cells are split into two random halves; each phase block is oriented on one half and counted on
#' the other, summed per window: `nA`, `nB`, `n_hap`, `f`, `f_se` (binomial SE x
#' the beta-binomial inflation `1 + (nbar - 1) * rho`, nbar = molecule-weighted molecules per cell), and
#' `n_blocks`.
#' @param win 1 Mb window depth: list(windows (window, chrom, start, end), x, e) with x / e windows x cells
#' @param votes,atoms votes and atoms tables
#' @param assign data.table cell, clone (cells used for each clone's pseudobulk)
#' @param ref_cells reference pool cells (NULL = none)
#' @param size window size (bp; a multiple of 1 Mb), one value or one per clone (named)
#' @param seed RNG seed of the split halves
#' @param cells cell names in vote index order (default the columns of `win$x`)
#' @param rho beta-binomial overdispersion
#' @param allelic_chroms chromosomes with an allele channel (NULL = all)
#' @param haploid chromosomes with one copy in a normal cell (male chrX / chrY): `ploidy` 1
#' @param min_ref minimal reference pool size
#' @return data.table clone, window, chrom, start, end, ploidy, base (copy number of log2 = 0: the ploidy when
#'   relative to the pool, else 2), n_cells, log2, log2_se, nA, nB, n_hap, f, f_se, n_blocks, ref ("pool" or
#'   "median"), size
#' @export
clone_profiles = function(win, votes, atoms, assign, ref_cells = NULL, size = 1e6, seed = 1, cells = colnames(win$x),
                          rho = 0.01, allelic_chroms = NULL, haploid = character(), min_ref = 30) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  assign = as.data.table(assign)[!is.na(clone), .(cell, clone = as.character(clone))]
  cl = sort(unique(assign$clone))
  if (length(size) == 1) size = stats::setNames(rep(size, length(cl)), cl)
  at = as.data.table(atoms)[, .(atom, chrom, astart = start)]
  v_all = as.data.table(votes)[, .(cell, atom, ps, nA, nB)]
  v_all = at[v_all, on = "atom", nomatch = NULL]
  if (!is.null(allelic_chroms)) v_all = v_all[chrom %in% allelic_chroms]
  v_all[, cname := cells[cell]]
  out = list()
  for (sz in unique(size)) {
    w = win_aggregate(win, sz)
    lr = log2((w$x + 0.5) / (w$e + 0.5))
    auto = !w$windows$chrom %in% c(haploid, "chrX", "chrY")
    ref = intersect(ref_cells, colnames(w$x))
    use_ref = length(ref) >= min_ref
    if (use_ref) {
      rl = log2((rowSums(w$x[, ref, drop = FALSE]) + 0.5) / (rowSums(w$e[, ref, drop = FALSE]) + 0.5))
      rse = row_sd(lr[, ref, drop = FALSE]) / sqrt(length(ref))
    }
    for (k in names(size)[size == sz]) {
      cc = intersect(assign[clone == k, cell], colnames(w$x))
      if (!length(cc)) next
      l2 = log2((rowSums(w$x[, cc, drop = FALSE]) + 0.5) / (rowSums(w$e[, cc, drop = FALSE]) + 0.5))
      se = row_sd(lr[, cc, drop = FALSE]) / sqrt(length(cc))
      if (use_ref) { l2 = l2 - rl; se = sqrt(se^2 + rse^2) }
      l2 = l2 - stats::median(l2[auto], na.rm = TRUE)
      pr = data.table(clone = k, window = w$windows$window, chrom = w$windows$chrom, start = w$windows$start,
                      end = w$windows$end, n_cells = length(cc), log2 = l2, log2_se = pmax(se, 1e-3),
                      ref = if (use_ref) "pool" else "median", size = sz)
      # split-half oriented votes
      vk = v_all[cname %in% cc]
      if (nrow(vk)) {
        set.seed(seed)
        hh = stats::setNames(sample(rep(1:2, length.out = length(cc))), cc)
        vk[, half := hh[cname]]
        bt = vk[, .(gA = sum(nA), gB = sum(nB)), by = .(chrom, ps, half)]
        o = dcast(bt, chrom + ps ~ half, value.var = c("gA", "gB"), fill = 0)
        for (cn_ in c("gA_1", "gA_2", "gB_1", "gB_2")) if (!cn_ %in% names(o)) o[, (cn_) := 0]
        o[, `:=`(s1 = sign(gA_1 - gB_1), s2 = sign(gA_2 - gB_2))]
        vk[o, `:=`(s = fifelse(half == 1L, i.s2, i.s1)), on = .(chrom, ps)]
        vk = vk[s != 0]
        vk[, `:=`(oA = fifelse(s > 0, nA, nB), oB = fifelse(s > 0, nB, nA))]
        vk[, wstart := floor(astart / sz) * sz]
        pc = vk[, .(oA = sum(oA), oB = sum(oB)), by = .(cname, chrom, wstart)]
        pw = pc[, .(nA = sum(oA), nB = sum(oB), nbar = sum((oA + oB)^2) / max(sum(oA + oB), 1)), by = .(chrom, wstart)]
        nbk = vk[, .(n_blocks = uniqueN(ps)), by = .(chrom, wstart)]
        pw = merge(pw, nbk, by = c("chrom", "wstart"))
        pr[, wstart := floor(start / sz) * sz]
        pr = merge(pr, pw, by = c("chrom", "wstart"), all.x = TRUE, sort = FALSE)
        pr[, wstart := NULL]
      } else pr[, `:=`(nA = NA_real_, nB = NA_real_, nbar = NA_real_, n_blocks = NA_real_)]
      out[[k]] = pr
    }
  }
  p = rbindlist(out, fill = TRUE)
  p[is.na(nA), `:=`(nA = 0, nB = 0, nbar = 1, n_blocks = 0)]
  p[, n_hap := nA + nB]
  p[, f := fifelse(n_hap > 0, nA / n_hap, NA_real_)]
  p[, f_se := fifelse(n_hap > 0, sqrt(pmax(f * (1 - f), 0.01) / n_hap * (1 + (nbar - 1) * rho)), NA_real_)]
  p[, ploidy := fifelse(chrom %in% haploid, 1L, 2L)]
  # copy-number baseline of the depth: the reference's ploidy when relative to a pool, else the diploid median
  p[, base := fifelse(ref == "pool", ploidy, 2L)]
  p[, nbar := NULL]
  p = p[order(clone, chrom_rank(chrom), start)]
  setcolorder(p, c("clone", "window", "chrom", "start", "end", "ploidy", "base", "n_cells", "log2", "log2_se", "nA", "nB",
                   "n_hap", "f", "f_se", "n_blocks", "ref", "size"))
  p[]
}
utils::globalVariables(c("cname", "wstart", "nbar", "i.s1", "i.s2", "gA_1", "gA_2", "gB_1", "gB_2", "s1", "s2", "astart"))

# ---- states, fit and segmentation -----------------------------------------------------------------

#' Allele-specific copy-number state grid
#'
#' nA:nB states relative to a diploid clone baseline: neutral 1:1, loss 1:0, cnloh 2:0, gain 2:1,
#' gain_as 3:0 (allele-specific gain), amp 3:1 / 2:2 / 4:1 / 4:0, del 0:0. `f` is the major-haplotype
#' fraction, `prior` a log weight used by [fit_clone_states()] (mixture over states) and by the Viterbi
#' transitions of [segment_clone()].
#'
#' With `baseline` (the copy number of the clone's median depth, e.g. 3 for a near-triploid clone) the grid holds
#' every state nA >= nB with nA + nB <= `max_cn`, and the prior falls off with the distance of the total copy number
#' from the baseline: log weight `-lam * |cn - baseline|`, normalised over the grid (so grids of different baselines
#' compare on one likelihood scale). Labels stay relative to the diploid normal (see `state_label()`).
#' @param max_cn maximal total copy number
#' @param baseline NULL (the diploid grid above) or the clone baseline the prior is centred on
#' @param lam prior fall-off per copy away from the baseline (with `baseline`)
#' @return data.table state, nA, nB, cn, f, label, prior
#' @export
cn_grid = function(max_cn = 5, baseline = NULL, lam = 1) {
  if (is.null(baseline)) {
    g = data.table(state = c("1:1", "1:0", "2:0", "2:1", "3:0", "3:1", "2:2", "4:1", "4:0", "0:0"),
                   nA = c(1L, 1L, 2L, 2L, 3L, 3L, 2L, 4L, 4L, 0L), nB = c(1L, 0L, 0L, 1L, 0L, 1L, 2L, 1L, 0L, 0L),
                   label = c("neutral", "loss", "cnloh", "gain", "gain_as", "amp", "amp", "amp", "amp", "del"),
                   prior = log(c(0.5, 0.08, 0.08, 0.08, 0.04, 0.04, 0.04, 0.03, 0.03, 0.02)))
    g[, cn := nA + nB]
    g = g[cn <= max_cn]
  } else {
    g = CJ(nA = 0:max_cn, nB = 0:max_cn)[nA >= nB & nA + nB <= max_cn]
    g[, `:=`(cn = nA + nB, state = paste0(nA, ":", nB))]
    g[, label := state_label(nA, nB, 2L)]
    g[, prior := -lam * abs(cn - baseline)]
    g[, prior := prior - (max(prior) + log(sum(exp(prior - max(prior)))))]
    # the 1:1 state first, as in the diploid grid (ties of the per-window maxima resolve to it)
    g = g[order(state != "1:1", cn, -nA)]
  }
  g[, f := fifelse(cn > 0, nA / cn, 0.5)]
  setcolorder(g, c("state", "nA", "nB", "cn", "f", "label", "prior"))
  g[]
}

# label of a state relative to the diploid normal (haploid chromosomes, ploidy 1: relative to one copy, 1:0 is
# neutral): 0 copies del, 1 loss, 2 neutral (1:1) or cnloh (2:0), 3 gain (2:1 ...) or gain_as (3:0), >= 4 amp. The
# states of the diploid cn_grid() keep their grid labels. `grid` is unused (kept for callers that pass it).
state_label = function(nA, nB, ploidy, grid = NULL) {
  cn = nA + nB
  lab = ifelse(cn == 0, "del", ifelse(cn == 1, "loss", ifelse(cn == 2, ifelse(nB == 0, "cnloh", "neutral"),
               ifelse(cn == 3, ifelse(nB == 0, "gain_as", "gain"), "amp"))))
  hap = ploidy == 1
  lab[hap] = ifelse(cn[hap] == 1, "neutral", ifelse(cn[hap] == 0, "del", ifelse(cn[hap] == 2, "gain", "amp")))
  lab
}

# per-window x state emission matrices for one clone profile: depth part (depends on c) and allele part (a, rho)
depth_emission = function(p, grid, c_, df = 5) {
  W = nrow(p)
  mean_ = c_ * log2(pmax(outer(rep(1, W), grid$cn), 0.2) / p$base)
  D = stats::dt((p$log2 - mean_) / p$se_eff, df, log = TRUE) - log(p$se_eff)
  D[!is.finite(p$log2), ] = 0
  hw = which(p$ploidy == 1)
  if (length(hw)) D[hw, !(grid$nB == 0 & grid$cn <= 3)] = -Inf
  D
}
# robust window likelihood: (1 - out_w) BB(rho) + out_w BB(rho_out) (orientation damage makes a few windows far
# more dispersed than the rest; the heavy component keeps them from overruling the depth)
bb_robust = function(x, n, f, rho_, eps, out_w = 0.05, rho_out = NULL) {
  rho_out = rho_out %||% max(0.05, 10 * rho_)
  l1 = bb_loglik(x, n, f, rho_, eps) + log(1 - out_w); l2 = bb_loglik(x, n, f, rho_out, eps) + log(out_w)
  m = pmax(l1, l2)
  m + log(exp(l1 - m) + exp(l2 - m))
}
allele_emission = function(p, grid, a_, rho_, eps = 0.02, out_w = 0.05) {
  W = nrow(p); K = nrow(grid)
  A = matrix(0, W, K)
  has = which(p$n_hap > 0 & p$ploidy == 2)
  if (length(has)) {
    x = p$nA[has]; n = p$n_hap[has]
    l0 = bb_robust(x, n, 0.5, rho_, eps, out_w)
    for (j in seq_len(K)) {
      fj = 0.5 + a_ * (grid$f[j] - 0.5)
      A[has, j] = if (abs(fj - 0.5) < 1e-9) l0 else pmax(bb_robust(x, n, fj, rho_, eps, out_w), bb_robust(x, n, 1 - fj, rho_, eps, out_w))
    }
  }
  A
}


# systematic window noise of a clone profile: residuals from a running median (21 windows) per chromosome;
# s_sys^2 = MAD(residual)^2 - median(SE^2); lag-1 autocorrelation of the residuals -> evidence inflation
# sqrt((1 + r1) / (1 - r1)) of the window SE (an AR(1) run of n windows carries ~n (1 - r1) / (1 + r1) windows of
# independent evidence; chromatin waves are autocorrelated over several Mb)
clone_noise = function(p) {
  r = p[is.finite(log2), {
    k = min(21L, .N - (.N + 1L) %% 2L)
    rr = if (.N >= 5) log2 - stats::runmed(log2, max(3L, k), endrule = "median") else rep(NA_real_, .N)
    list(r = rr, r_prev = c(NA_real_, rr[-.N]))
  }, by = chrom]
  ok = is.finite(r$r)
  if (sum(ok) < 10) return(list(s_sys = 0, acf1 = 0, infl = 1))
  s_tot = stats::mad(r$r[ok])
  s_sys = sqrt(max(0, s_tot^2 - stats::median(p$log2_se^2, na.rm = TRUE)))
  pp = is.finite(r$r) & is.finite(r$r_prev)
  a1 = if (sum(pp) > 10) suppressWarnings(stats::cor(r$r[pp], r$r_prev[pp], method = "spearman")) else 0
  a1 = min(max(ifelse(is.finite(a1), a1, 0), 0), 0.9)
  list(s_sys = s_sys, acf1 = a1, infl = sqrt((1 + a1) / (1 - a1)))
}

#' Fit clone copy-number states
#'
#' Per clone, a grid search over the depth contrast `c` (state mean `c * log2(cn / base)`, `base` = the
#' copy number of the clone's median depth: the clone baseline B, times the reference ploidy / 2 for pool-relative
#' profiles; haploid windows (male chrX / chrY, `ploidy` 1) allow only nB = 0 states and their labels are relative
#' to one copy),
#' the orientation attenuation `a` (expected major fraction `0.5 + a * (f_state - 0.5)`) and the
#' window-level overdispersion `rho`. Emissions: Student-t (`df`) on the window log2 with scale = the
#' window SE combined with the clone's systematic window noise and inflated for its autocorrelation
#' (see `clone_noise`), and a two-sided robust beta-binomial on (nA, n_hap) (a 5 % component with 10x the
#' overdispersion). The objective is the per-window mixture likelihood over the states with the grid
#' `prior` weights.
#'
#' The baseline B is 2 by default (the diploid baseline of 0.1.x: a near-triploid clone is then fitted on a
#' compressed scale). With `baseline = "fit"` every clone is fitted at each of `baselines` on the grid
#' `baseline_grid(B)` (all states up to 2 B + 1 copies without the near-balanced imbalanced ones, prior centred on
#' B and normalised, so the baselines compare on one likelihood scale; a coarse (c, a) grid, no refinement). A clone
#' whose best attenuation `a` is above the floor of `atten` has allelic information on the scale and takes its best
#' baseline, or half of it when that is within `ploidy_tol` x windows (a whole-genome doubling of the same profile
#' fits as well); the sample baseline is the cell-weighted median of these clones; the other clones take the better
#' of 2 and the sample baseline. Each clone is then fitted at its baseline (B = 2: `grid`, exactly as before).
#' @param profile output of [clone_profiles()]
#' @param grid state grid ([cn_grid()])
#' @param contrast grid of depth contrasts
#' @param atten grid of orientation attenuations
#' @param rho grid of window-level overdispersions
#' @param df degrees of freedom of the depth t
#' @param eps allele error floor
#' @param refine refine (c, a, rho) on the segment level: segment at the current fit (Viterbi, `p_switch`,
#'   `min_win`), then maximise the sum over segments of the best state's segment likelihood (depth on the
#'   precision-weighted segment mean, alleles summed over its windows); two passes
#' @param p_switch,min_win segmentation parameters of the refinement (as [segment_clone()])
#' @param baseline the clone baseline B: one number for every clone (2 = diploid), or "fit"
#' @param baselines candidate baselines of `baseline = "fit"`
#' @param lam prior fall-off per copy away from the baseline ([cn_grid()]) of the baseline grids
#' @param min_df imbalanced states of the baseline grids need |f - 0.5| >= this
#' @param ploidy_tol log-likelihood tolerance per window within which half the best baseline is preferred
#' @return list(ll = data.table clone, window, chrom, start, end, ploidy, state, label, ll (emission
#'   log-likelihood at the chosen parameters), windows = the profile with `se_eff` and `base`, par = data.table clone,
#'   c, a, rho, s_sys, loglik, n_win, base (the baseline B), grid (the B = 2 grid), grids (per clone), baseline_scan
#'   (with "fit": clone, B, ll, c, a per candidate), baseline_sample, baseline_informative (with "fit": TRUE when at
#'   least one clone had allelic information on the scale, i.e. the sample baseline was fitted, not the default 2),
#'   contrast_floor (the lower bound of `contrast`: a
#'   clone whose `c` sits there has a poorly determined scale, and its states can be off by one copy))
#' @export
fit_clone_states = function(profile, grid = cn_grid(), contrast = seq(0.6, 1, 0.05), atten = seq(0.6, 1, 0.05),
                            rho = c(0.001, 0.003, 0.01, 0.03), df = 5, eps = 0.02, refine = TRUE, p_switch = 1e-3, min_win = 3,
                            baseline = 2, baselines = c(1.5, 2, 2.5, 3, 3.25, 3.5, 4, 4.5, 5), lam = 1, min_df = 0.1,
                            ploidy_tol = 0.01) {
  prof = as.data.table(profile)
  if (!"ploidy" %in% names(prof)) prof[, ploidy := 2L]
  if (!"base" %in% names(prof)) prof[, base := 2L]
  # the profile's base is the copy number of log2 = 0 under a diploid clone baseline (the pool ploidy, or 2 when
  # median-centred); a clone baseline B scales it by B / 2
  prof[, base0_ := as.numeric(base)]
  fit_b = !is.numeric(baseline)
  if (fit_b && !identical(baseline, "fit")) stop("fit_clone_states: baseline must be a number or \"fit\"", call. = FALSE)
  if (!fit_b && (length(baseline) != 1 || !is.finite(baseline) || baseline <= 0))
    stop("fit_clone_states: baseline must be one positive number or \"fit\"", call. = FALSE)
  lls = list(); pars = list(); wins = list(); grids = list(); scans = list()
  K = unique(prof$clone)
  P = list()
  for (k in K) {
    p = prof[clone == k]
    ns = clone_noise(p)
    p[, se_eff := sqrt(log2_se^2 + ns$s_sys^2) * ns$infl]
    P[[k]] = list(p = p, ns = ns)
  }
  B_of = stats::setNames(rep(if (fit_b) 2 else baseline, length(K)), K)
  if (fit_b) {
    sc = baseline_scan(P, baselines, contrast, atten, rho, df, eps, lam, min_df)
    B_of = choose_baselines(sc, P, min(atten), ploidy_tol)
    scans = sc
  }
  for (k in K) {
    p = copy(P[[k]]$p); ns = P[[k]]$ns
    B = B_of[[k]]
    g = if (abs(B - 2) < 1e-9) grid else baseline_grid(B, lam, min_df)
    p[, base := base0_ * B / 2]
    best = fit_clone_core(p, g, contrast, atten, rho, df, eps, refine, p_switch, min_win)
    E = best$E
    lls[[k]] = data.table(clone = k, window = rep(p$window, nrow(g)), chrom = rep(p$chrom, nrow(g)),
                          start = rep(p$start, nrow(g)), end = rep(p$end, nrow(g)), ploidy = rep(p$ploidy, nrow(g)),
                          state = rep(g$state, each = nrow(p)), ll = as.vector(E),
                          ll_a = as.vector(allele_emission(p, g, best$a, best$rho, eps)))
    pars[[k]] = data.table(clone = k, c = best$c, a = best$a, rho = best$rho, s_sys = ns$s_sys, acf1 = ns$acf1, infl = ns$infl,
                           loglik = best$ll, n_win = nrow(p), base = B)
    p[, base0_ := NULL]
    wins[[k]] = p
    grids[[k]] = g
  }
  ll = rbindlist(lls)
  ll[, label := state_label(as.integer(sub(":.*", "", state)), as.integer(sub(".*:", "", state)), ploidy)]
  structure(list(ll = ll, windows = rbindlist(wins), par = rbindlist(pars), grid = grid, grids = grids,
                 baseline_scan = if (length(scans)) scans else NULL,
                 baseline_sample = if (fit_b) attr(B_of, "sample") else NULL,
                 baseline_informative = if (fit_b) attr(B_of, "informative") else NULL, contrast_floor = min(contrast)),
            class = "clone_states")
}
utils::globalVariables(c("base0_", "B"))

# state grid of a clone baseline B != 2: every nA >= nB up to max(5, 2 B + 1) copies without the near-balanced
# imbalanced states (0 < |f - 0.5| < min_df: 3:2 stays, 4:3 / 5:4 go; above a baseline of ~3 they only absorb the
# orientation bias of the split-half fractions), prior centred on B ([cn_grid()] `baseline`)
baseline_grid = function(B, lam = 1, min_df = 0.1) {
  g = cn_grid(max(5L, as.integer(ceiling(2 * B)) + 1L), baseline = B, lam = lam)
  g = g[nA == nB | abs(f - 0.5) >= min_df - 1e-9]
  g[, prior := prior - (max(prior) + log(sum(exp(prior - max(prior)))))]
  g[]
}

# per clone x candidate baseline: the window mixture log-likelihood (normalised prior, so baselines compare) of a
# coarse (c, a, rho) grid without refinement, and the attenuation at the best point
baseline_scan = function(P, baselines, contrast, atten, rho, df, eps, lam, min_df) {
  cs = contrast[unique(c(seq(1, length(contrast), by = 2), length(contrast)))]
  as_ = atten[unique(c(seq(1, length(atten), by = 2), length(atten)))]
  out = list()
  for (B in baselines) {
    g = baseline_grid(B, lam, min_df)
    for (k in names(P)) {
      p = copy(P[[k]]$p)[, base := base0_ * B / 2]
      b = fit_clone_core(p, g, cs, as_, rho, df, eps, refine = FALSE)
      out[[length(out) + 1]] = data.table(clone = k, B = B, ll = b$ll + nrow(p) * max(g$prior), c = b$c, a = b$a,
                                          n_cells = suppressWarnings(max(p$n_cells, na.rm = TRUE)), n_win = nrow(p))
    }
  }
  rbindlist(out)
}

# baselines from a baseline_scan(): a clone whose allele attenuation at its best baseline is above the grid floor
# carries allelic information on the scale and takes the smallest baseline within ploidy_tol x windows of its best
# (WGD-type ties go to the smaller genome); the sample baseline is the cell-weighted median of those clones (2 when
# there are none); every other clone takes the better of 2 and the sample baseline
choose_baselines = function(sc, P, a_floor, ploidy_tol) {
  K = names(P)
  # the best baseline, or half of it (a whole-genome doubling of the same profile) when that is within the tolerance
  pick = function(s) {
    b = s[which.max(ll)]
    h = s[abs(B - b$B / 2) <= 0.25 + 1e-9 & ll >= b$ll - ploidy_tol * b$n_win]
    if (nrow(h)) h[which.max(ll)] else b
  }
  bestk = rbindlist(lapply(K, function(k) sc[clone == k][which.max(ll)]))
  pk = rbindlist(lapply(K, function(k) pick(sc[clone == k])))
  inf = bestk$a > a_floor + 1e-9
  Bs = if (any(inf)) {
    w = pk$n_cells[inf]; w[!is.finite(w)] = 1
    b = pk$B[inf]; o = order(b); cw = cumsum(w[o]) / sum(w)
    b[o][which(cw >= 0.5)[1]]
  } else 2
  out = stats::setNames(numeric(length(K)), K)
  for (i in seq_along(K)) {
    s = sc[clone == K[i]]
    out[i] = if (inf[i]) pk$B[i] else { s2 = s[B %in% c(2, Bs)]; if (nrow(s2)) pick(s2)$B else 2 }
  }
  attr(out, "sample") = Bs
  attr(out, "informative") = any(inf)
  out
}

# one clone: grid search over (c, a, rho) of the window mixture likelihood, then (refine) the segment-level refinement
fit_clone_core = function(p, grid, contrast, atten, rho, df, eps, refine = TRUE, p_switch = 1e-3, min_win = 3) {
    lp = grid$prior - max(grid$prior)
    best = list(ll = -Inf)
    Ds = lapply(contrast, function(c_) depth_emission(p, grid, c_, df))
    for (r_ in rho) for (a_ in atten) {
      A = allele_emission(p, grid, a_, r_, eps)
      for (ic in seq_along(contrast)) {
        M = Ds[[ic]] + A + matrix(lp, nrow(p), length(lp), byrow = TRUE)
        mx = do.call(pmax, as.data.frame(M))
        tot = sum(mx + log(rowSums(exp(M - mx))))
        if (tot > best$ll) best = list(ll = tot, c = contrast[ic], a = a_, rho = r_, E = Ds[[ic]] + A)
      }
    }
    if (refine) for (it in seq_len(2)) {
      # segment-level refinement: segment the windows at the current fit, then choose (c, a, rho) by the
      # segment likelihood (depth: the segment's precision-weighted mean log2; alleles: summed window terms),
      # which a few noisy windows cannot move the way they move the per-window mixture
      sid = integer(nrow(p)); off = 0L
      for (ch in unique(p$chrom)) {
        ii = which(p$chrom == ch)
        pth = seg_path(best$E[ii, , drop = FALSE], p_switch, lp, min_win)
        rr = rle(pth); sid[ii] = off + rep(seq_along(rr$lengths), rr$lengths); off = off + length(rr$lengths)
      }
      fin = is.finite(p$log2)
      wv = 1 / p$se_eff^2
      sm = rowsum(cbind(ifelse(fin, p$log2 * wv, 0), ifelse(fin, wv, 0)), sid)
      m_seg = sm[, 1] / pmax(sm[, 2], 1e-12); se_seg = 1 / sqrt(pmax(sm[, 2], 1e-12))
      base_seg = as.vector(rowsum(p$base, sid) / tabulate(sid))
      hap_seg = as.vector(rowsum(as.numeric(p$ploidy == 1), sid)) > 0
      mask = matrix(0, length(m_seg), nrow(grid)); mask[hap_seg, !(grid$nB == 0 & grid$cn <= 3)] = -Inf
      has_d = sm[, 2] > 0
      best2 = list(ll = -Inf)
      for (ic in seq_along(contrast)) {
        mu = contrast[ic] * log2(pmax(outer(rep(1, length(m_seg)), grid$cn), 0.2) / base_seg)
        Dg = stats::dt((m_seg - mu) / se_seg, df, log = TRUE) - log(se_seg)
        Dg[!has_d, ] = 0
        for (r_ in rho) for (a_ in atten) {
          Ag = rowsum(allele_emission(p, grid, a_, r_, eps), sid)
          tot = sum(apply(Dg + Ag + mask, 1, max))
          if (tot > best2$ll) best2 = list(ll = tot, c = contrast[ic], a = a_, rho = r_)
        }
      }
      changed = !(best2$c == best$c && best2$a == best$a && best2$rho == best$rho)
      best$c = best2$c; best$a = best2$a; best$rho = best2$rho
      best$E = depth_emission(p, grid, best$c, df) + allele_emission(p, grid, best$a, best$rho, eps)
      if (!changed) break
    }
    best
}
utils::globalVariables(c("se_eff", "d", "base", "r_prev"))

viterbi_path = function(E, p_switch, lw) {
  W = nrow(E); K = ncol(E)
  if (W == 0) return(integer())
  w = exp(lw); w = w / sum(w)
  # transition from i to j != i: p_switch * w_j / (1 - w_i)
  lt = matrix(0, K, K)
  for (i in seq_len(K)) { lt[i, ] = log(p_switch * w / (1 - w[i])); lt[i, i] = log(1 - p_switch) }
  dp = matrix(-Inf, W, K); bp = matrix(0L, W, K)
  dp[1, ] = lw + E[1, ]
  for (t in seq_len(W)[-1]) {
    cand = dp[t - 1, ] + lt                # K x K: from i (rows) to j (cols)
    bi = max.col(t(cand), ties.method = "first")
    dp[t, ] = cand[cbind(bi, seq_len(K))] + E[t, ]
    bp[t, ] = bi
  }
  path = integer(W); path[W] = which.max(dp[W, ])
  for (t in rev(seq_len(W - 1))) path[t] = bp[t + 1, path[t + 1]]
  path
}

# Viterbi path of one chromosome, runs shorter than min_win absorbed into the better adjacent run
seg_path = function(E, p_switch, lw, min_win) {
  path = viterbi_path(E, p_switch, lw)
  if (nrow(E) >= min_win) repeat {
    r = rle(path)
    if (length(r$lengths) <= 1 || all(r$lengths >= min_win)) break
    j = which(r$lengths < min_win)[which.min(r$lengths[r$lengths < min_win])]
    ends = cumsum(r$lengths); sts = ends - r$lengths + 1
    rows = sts[j]:ends[j]
    cand = c(if (j > 1) r$values[j - 1], if (j < length(r$values)) r$values[j + 1])
    sc = vapply(cand, function(s_) sum(E[rows, s_]), numeric(1))
    path[rows] = cand[which.max(sc)]
  }
  path
}

#' Segment clone profiles (Viterbi)
#'
#' Viterbi per clone and chromosome over the states of [fit_clone_states()] (stay probability
#' `1 - p_switch`, switches weighted by the grid prior); runs shorter than `min_win` windows are absorbed
#' into the adjacent run whose state explains them better. With the allele part of the emissions (`ll_a`, from
#' [fit_clone_states()]) a run in an imbalanced state whose allelic log-likelihood ratio vs 1:1 is not positive is
#' given the best balanced state instead. Clone profiles of clones under 50 cells are
#' expected on the 5 Mb tile grid ([caller_tumour()] does that).
#'
#' The rules below are anchored on the clone's integer baseline nb (its fitted baseline `base` rounded; 2 without a
#' fit or at the diploid baseline): copy-neutral states have nb copies (CN-LOH 2:0 and 1:1 when nb = 2, e.g. 2:1 and
#' 3:0 when nb = 3), and the balanced reference is nb/2:nb/2 (nb even) or the better supported of the balanced states
#' next to nb (nb odd, e.g. 1:1 or 2:2 when nb = 3). At nb = 2 they are the rules up to 0.1.x.
#' @param ll a [fit_clone_states()] result, or a long table clone, chrom, start, end, state, ll (optional
#'   window, ploidy; then summaries that need the profile are NA)
#' @param p_switch switch probability between adjacent windows
#' @param min_win minimal windows per segment
#' @param grid state grid (default the fit's grid or [cn_grid()])
#' @param min_llr segments with llr above this are `confident` (neutral segments always are)
#' @param min_excess with a [fit_clone_states()] input: an imbalanced run also needs its pooled split-half oriented
#'   fraction to exceed 0.5 by at least this share of the expected excess `a * (f_state - 0.5)` (0 = off)
#' @param bal_frac,min_bal_len with a [fit_clone_states()] input: a balanced run other than the balanced reference
#'   (e.g. 2:2 at baseline 2; depth only) needs a mean log2 of at least `bal_frac` of its expected shift
#'   `c * log2(cn / base)` and `min_bal_len` bp, else it takes the balanced reference (1:1 at baseline 2) (0 = off)
#' @param max_contra with a [fit_clone_states()] input: an imbalanced run whose copy number differs from nb and whose
#'   mean log2 has the opposite sign of its expected shift by more than this takes the best supported copy-neutral
#'   imbalanced state (nb copies: CN-LOH at baseline 2), or the balanced reference when none is supported (Inf = off)
#' @param gain_support with a [fit_clone_states()] input: an imbalanced gain (diploid window, copy number above the
#'   clone's integer baseline nb) needs a mean log2 of at least this share of its expected shift `c * log2(cn / base)`, measured
#'   against the region's depth offset (per window, the median over the other clones that are 1:1 there, with at
#'   least 50 cells, of their log2 minus the 1:1 expectation). A run that fails takes the best copy-neutral imbalanced
#'   state (nb copies: CN-LOH at baseline 2) if its pooled oriented fraction supports one (as `min_excess`); otherwise
#'   it keeps its state with `label` "imbalance" and `depth_ok` FALSE (allelic imbalance without the extra copy,
#'   typically a clone mixing LOH and heterozygous cells). 0 = off (no `depth_ok` column; the behaviour up to 0.1.1).
#'   Only applies to a [fit_clone_states()] input.
#' @return segments table ([check_segments()]) plus `label`, `ploidy` (and `depth_ok` when `gain_support` > 0)
#' @export
segment_clone = function(ll, p_switch = 1e-3, min_win = 3, grid = NULL, min_llr = 10, min_excess = 0.5, bal_frac = 0.75,
                         min_bal_len = 2e7, max_contra = 0.1, gain_support = 0.3) {
  win = NULL; par = NULL; grids = list()
  if (inherits(ll, "clone_states")) {
    win = ll$windows; par = ll$par; grids = ll$grids %||% list(); grid = grid %||% ll$grid; ll = ll$ll
  }
  grid = grid %||% cn_grid()
  ll = as.data.table(ll)
  if (!"window" %in% names(ll)) ll[, window := paste0(chrom, ":", format(start, scientific = FALSE, trim = TRUE))]
  if (!"ploidy" %in% names(ll)) ll[, ploidy := 2L]
  out = list()
  keep = list()   # per clone x chromosome: windows, emissions, final path (the gain-support pass needs every clone's path)
  for (k in unique(ll$clone)) for (ch in unique(ll[clone == k, chrom])) {
    # the clone's state grid (a fitted baseline B != 2 has its own grid) and its integer baseline nb: the copy-neutral
    # states are those with nb copies, the balanced reference of the depth-only rules is nb/2:nb/2 (nb even) or the
    # balanced states next to it (nb odd); with B = 2 these are 2 copies and 1:1, the rules up to 0.1.x
    gk = grids[[k]] %||% grid
    st = intersect(gk$state, unique(ll[clone == k]$state))
    g = gk[match(st, state)]
    lw = g$prior - max(g$prior)
    Bk = if (!is.null(par) && "base" %in% names(par) && nrow(par[clone == k]) == 1) par[clone == k, base] else 2
    nb = as.integer(floor(Bk + 0.5))
    bal = which(g$nA == g$nB)
    n0b = if (nb %% 2 == 0) which(g$nA == g$nB & g$cn == nb) else which(g$nA == g$nB & abs(g$cn - nb) == 1)
    wk = NULL; pk = NULL
    d = ll[clone == k & chrom == ch]
    wn = unique(d[, .(window, start, end, ploidy)])[order(start)]
    E = matrix(-Inf, nrow(wn), length(st))
    E[cbind(match(d$window, wn$window), match(d$state, st))] = d$ll
    E[is.na(E)] = -Inf
    path = seg_path(E, p_switch, lw, min_win)
    if ("ll_a" %in% names(d)) {
      # an imbalanced state needs the allele channel's support over its run (allelic LLR vs 1:1 > 0); otherwise the
      # run takes the best balanced state: depth offsets over long stretches must not create 2:1 calls at f ~ 0.5
      A = matrix(0, nrow(wn), length(st))
      A[cbind(match(d$window, wn$window), match(d$state, st))] = d$ll_a
      n11 = which(st == "1:1")
      # best balanced reference state of a run (by its emissions)
      ref_bal = function(rows) n0b[which.max(colSums(E[rows, n0b, drop = FALSE]))]
      r = rle(path); ends = cumsum(r$lengths); sts = ends - r$lengths + 1
      wk = if (!is.null(win)) win[clone == k & chrom == ch][match(wn$window, window)] else NULL
      pk = if (!is.null(par)) par[clone == k] else NULL
      for (j in seq_along(r$values)) {
        rows = sts[j]:ends[j]
        s_ = r$values[j]
        if (wn$ploidy[rows[1]] != 2 || !length(n11) || !length(n0b)) next
        if (g$nA[s_] != g$nB[s_]) {
          ok_ = sum(A[rows, s_]) - sum(A[rows, n11]) > 0
          # ... and, pooled over the run, the split-half oriented fraction must show at least `min_excess` of the
          # expected excess over 0.5 (the per-window two-sided allelic term is biased towards imbalance where the
          # oriented fraction is ~0.5 with noise)
          fo = 0.5
          if (!is.null(wk) && nrow(pk) == 1) {
            nh = sum(wk$n_hap[rows], na.rm = TRUE)
            fo = if (nh > 0) sum(wk$nA[rows], na.rm = TRUE) / nh else 0.5
          }
          if (ok_ && !is.null(wk) && nrow(pk) == 1 && is.finite(min_excess) && min_excess > 0)
            ok_ = (fo - 0.5) >= min_excess * pk$a * (g$f[s_] - 0.5)
          if (ok_ && !is.null(wk) && nrow(pk) == 1 && g$cn[s_] != nb && is.finite(max_contra)) {
            # an imbalanced gain / loss whose run depth clearly contradicts it (opposite sign beyond max_contra log2)
            # takes the best copy-neutral state instead (CN-LOH if its allelic excess is supported, else 1:1): a clone
            # with a minority of imbalanced cells must not be given a gain at a depth below its baseline
            exp_ = pk$c * log2(max(g$cn[s_], 0.2) / (wk$base[rows[1]] %||% 2))
            obs_ = mean(wk$log2[rows], na.rm = TRUE)
            if (is.finite(obs_) && sign(obs_) == -sign(exp_) && abs(obs_) > max_contra) {
              cn2 = which(g$cn == nb & g$nA != g$nB)
              cn2 = cn2[(fo - 0.5) >= min_excess * pk$a * (g$f[cn2] - 0.5)]
              path[rows] = if (length(cn2)) cn2[which.max(colSums(E[rows, cn2, drop = FALSE]))] else ref_bal(rows)
              next
            }
          }
          if (!ok_) path[rows] = bal[which.max(colSums(E[rows, bal, drop = FALSE]))]
        }
      }
      if (!is.null(wk) && !is.null(pk) && nrow(pk) == 1 && is.finite(bal_frac) && bal_frac > 0) {
        # a balanced non-neutral run (depth only: 2:2, 0:0) needs >= min_bal_len and a mean log2 of >= bal_frac of its
        # expected shift; otherwise 1:1 (clone-level depth has chromatin / GC waves that the allele channel cannot check)
        r = rle(path); ends = cumsum(r$lengths); sts = ends - r$lengths + 1
        for (j in seq_along(r$values)) {
          rows = sts[j]:ends[j]
          s_ = r$values[j]
          # (nb even: the nb/2:nb/2 reference is exempt; nb odd: both balanced neighbours are tested and a run that
          # fails takes the better supported of them)
          if (wn$ploidy[rows[1]] != 2 || g$nA[s_] != g$nB[s_] || !length(n0b) || (nb %% 2 == 0 && s_ %in% n0b)) next
          exp_ = pk$c * log2(max(g$cn[s_], 0.2) / (wk$base[rows[1]] %||% 2))
          obs_ = mean(wk$log2[rows], na.rm = TRUE)
          len_ = wn$end[rows[length(rows)]] - wn$start[rows[1]]
          if (len_ < min_bal_len || !is.finite(obs_) || sign(obs_) != sign(exp_) || abs(obs_) < bal_frac * abs(exp_))
            path[rows] = ref_bal(rows)
        }
      }
    }
    keep[[length(keep) + 1]] = list(k = k, ch = ch, wn = wn, E = E, path = path, wk = wk, pk = pk,
                                    nodepth = logical(nrow(wn)), st = st, g = g, nb = nb)
  }
  gs_on = is.finite(gain_support) && gain_support > 0 && !is.null(win) && !is.null(par)
  if (gs_on) keep = gain_support_pass(keep, gain_support, min_excess)
  for (it in keep) {
    k = it$k; ch = it$ch; wn = it$wn; E = it$E; path = it$path; st = it$st
    # a flagged run is its own segment even when a neighbour has the same state with depth support
    key = paste(path, it$nodepth)
    r = rle(key); ends = cumsum(r$lengths); sts = ends - r$lengths + 1
    for (j in seq_along(r$values)) {
      rows = sts[j]:ends[j]
      s_ = path[rows[1]]
      pl = wn$ploidy[rows[1]]
      n0 = if (pl == 1) which(st == "1:0") else which(st == "1:1")
      llr = sum(E[rows, s_]) - sum(E[rows, n0])
      out[[length(out) + 1]] = data.table(clone = k, chrom = ch, start = wn$start[rows[1]], end = wn$end[rows[length(rows)]],
                                          n_win = length(rows), state = st[s_], ploidy = pl, llr = llr,
                                          depth_ok = !it$nodepth[rows[1]],
                                          w_from = wn$window[rows[1]], w_ids = list(wn$window[rows]))
    }
  }
  seg = rbindlist(out)
  seg[, `:=`(nA_s = nA_of(state), nB_s = nB_of(state))]
  seg[, cn := as.integer(nA_s + nB_s)]
  seg[, label := state_label(nA_s, nB_s, ploidy)]
  if (gs_on) seg[depth_ok == FALSE, label := "imbalance"] else seg[, depth_ok := NULL]
  if (!is.null(win)) {
    sm = rbindlist(lapply(seq_len(nrow(seg)), function(i) {
      wv = win[clone == seg$clone[i] & window %in% seg$w_ids[[i]]]
      data.table(nA = as.integer(sum(wv$nA)), nB = as.integer(sum(wv$nB)), log2 = mean(wv$log2, na.rm = TRUE))
    }))
    seg[, c("nA", "nB", "log2") := sm]
  } else seg[, `:=`(nA = NA_integer_, nB = NA_integer_, log2 = NA_real_)]
  seg[, f := fifelse(nA + nB > 0, nA / (nA + nB), NA_real_)]
  seg[, confident := label == "neutral" | llr > min_llr]
  seg[, c("w_from", "w_ids", "nA_s", "nB_s") := NULL]
  seg[, n_win := as.integer(n_win)]
  setcolorder(seg, intersect(c("clone", "chrom", "start", "end", "n_win", "nA", "nB", "cn", "state", "log2", "f", "llr", "confident",
                               "label", "ploidy", "depth_ok"), names(seg)))
  seg = seg[order(clone, chrom_rank(chrom), start)]
  check_segments(seg)
  seg[]
}
utils::globalVariables(c("w_ids", "w_from", "depth_ok"))

# gain-support pass of segment_clone(): `keep` holds per clone x chromosome the windows (wn), emissions (E), final path,
# the clone's profile windows (wk), parameters (pk), states (st, g) and integer baseline (nb). The regional depth offset
# of a window is the median, over the clones that are 1:1 there (diploid window, >= 50 cells), of their log2 minus the
# 1:1 expectation c log2(2 / base) (0 for a diploid baseline); a clone's own window never enters its offset because
# only non-1:1 runs are tested. An imbalanced gain (more copies than the clone's integer baseline nb; at a half-integer
# baseline such as 2.5, nb = 3 and 2:1 is copy-neutral, not a gain) whose run log2 minus the run's mean
# offset is below gain_support x its expected shift takes the best copy-neutral imbalanced state (nb copies; CN-LOH for
# a diploid baseline) when its pooled oriented fraction supports one (as min_excess), else it is flagged
# (nodepth = TRUE) and keeps its state.
gain_support_pass = function(keep, gain_support, min_excess) {
  ref = rbindlist(lapply(keep, function(it) {
    if (is.null(it$wk) || !all(c("log2", "n_cells") %in% names(it$wk))) return(NULL)
    n11 = which(it$st == "1:1")
    if (!length(n11)) return(NULL)
    ok = it$path == n11 & it$wn$ploidy == 2 & it$wk$n_cells >= 50 & is.finite(it$wk$log2)
    if (!any(ok)) return(NULL)
    e11 = if (!is.null(it$pk) && nrow(it$pk) == 1 && "base" %in% names(it$wk)) it$pk$c * log2(2 / it$wk$base[ok]) else 0
    data.table(window = it$wn$window[ok], log2 = it$wk$log2[ok] - e11)
  }))
  off = if (nrow(ref)) ref[, .(off = stats::median(log2)), by = window] else data.table(window = character(), off = numeric())
  for (i in seq_along(keep)) {
    it = keep[[i]]
    wk = it$wk; pk = it$pk; g = it$g
    if (is.null(wk) || is.null(pk) || nrow(pk) != 1 || !any(it$st == "1:1")) next
    path = it$path
    r = rle(path); ends = cumsum(r$lengths); sts = ends - r$lengths + 1
    for (j in seq_along(r$values)) {
      rows = sts[j]:ends[j]
      s_ = r$values[j]
      base_ = wk$base[rows[1]] %||% 2
      # a gain = more copies than the clone's integer baseline nb (the copy-neutral states have nb copies, as in
      # segment_clone()); at an integer baseline nb is the diploid window's base, at B = 2.5 a 2:1 run is copy-neutral
      if (it$wn$ploidy[rows[1]] != 2 || g$nA[s_] == g$nB[s_] || g$cn[s_] <= it$nb) next
      exp_ = pk$c * log2(g$cn[s_] / base_)
      obs_ = mean(wk$log2[rows], na.rm = TRUE)
      o_ = off$off[match(it$wn$window[rows], off$window)]
      o_ = if (any(is.finite(o_))) mean(o_, na.rm = TRUE) else 0
      if (!is.finite(obs_) || !is.finite(exp_) || exp_ <= 0 || obs_ - o_ >= gain_support * exp_) next
      nh = sum(wk$n_hap[rows], na.rm = TRUE)
      fo = if (nh > 0) sum(wk$nA[rows], na.rm = TRUE) / nh else 0.5
      cn2 = setdiff(which(g$cn == it$nb & g$nA != g$nB), s_)     # never the run's own state
      if (is.finite(min_excess) && min_excess > 0) cn2 = cn2[(fo - 0.5) >= min_excess * pk$a * (g$f[cn2] - 0.5)]
      if (length(cn2)) path[rows] = cn2[which.max(colSums(it$E[rows, cn2, drop = FALSE]))]
      else it$nodepth[rows] = TRUE
    }
    it$path = path
    keep[[i]] = it
  }
  keep
}
