# window_allelic.R
# What: the windowed, orientation-free allele term of normal mode. Per cell x unit (one subwindow of at most
#       `win_size` bp x one phase set) the mirrored beta-binomial LLR of an imbalance f vs 1:1, summed per region;
#       the window columns of the long features table (`wal_f<fff>`, `n_win`) and their note.
# Why: on chromosome-scale (hybrid) phasing the oriented per-cell allele sums cancel an imbalance private to one cell
#      at population-phase switches, so a cell with whole-chromosome LOH can go uncalled and even become a null cell
#      (seen in a COLO829BL lymphoblastoid sample). Inside a window the phase is reliable; across windows it is not
#      needed.
# Inputs: votes / atoms / regions of cell_features(). Outputs: per cell x region window LLRs (cell_features() writes
#         them into `long`; score_cells(window = ) and select_null_cells() use them).

utils::globalVariables(c("sw", "masked", "i.nA", "i.n", "n_win"))

#' Haplotype fractions of the windowed allele term
#'
#' The fractions at which [window_allelic_features()] evaluates its windowed LLR: `params$win_f` (default 0.98, the
#' loss / CN-LOH rows of the normal state set, `1 - eps`). Only near-complete LOH is specific at window level: cells
#' of normal tissue show a mild chromosome-scale allelic skew (best window fraction 0.67-0.80 on many chromosomes in
#' the deepest cells of a normal thyroid sample), which a window term at gain fractions would call; at 0.98 no cell of
#' that sample has a chromosome with W > 10, while 99 % of the COLO829BL rows with W > 20 lie on ASCAT.sc-HP LOH
#' segments.
#' @param params [caller_defaults()] (`win_f`, `eps`)
#' @return sorted numeric vector
#' @export
win_f_grid = function(params = caller_defaults()) sort(unique(round(params$win_f %||% (1 - params$eps), 6)))

#' Column name of the windowed allele LLR at fraction f
#' @param f haplotype fraction(s)
#' @return `wal_f<1000 f>` (e.g. `wal_f980` for 0.98, `wal_f667` for 2/3)
#' @export
win_col = function(f) sprintf("wal_f%03d", as.integer(round(1000 * f)))

#' Mirrored beta-binomial LLR of one unit
#'
#' `log(0.5 BB(a; n, f) + 0.5 BB(a; n, 1 - f)) - log BB(a; n, 1/2)`: the evidence for an imbalance f in either
#' orientation against 1:1. A unit with n = 1 gives exactly 0.
#' @param a hap-A votes
#' @param n votes
#' @param f imbalance (major-haplotype fraction)
#' @param rho beta-binomial overdispersion
#' @param eps allele error floor
#' @return numeric vector
#' @export
mirrored_llr = function(a, n, f, rho, eps = 0.02) {
  l1 = bb_loglik(a, n, f, rho, eps); l2 = bb_loglik(a, n, 1 - f, rho, eps)
  log(0.5) + pmax(l1, l2) + log1p(exp(-abs(l1 - l2))) - bb_loglik(a, n, 0.5, rho, eps)
}

#' Window units of the votes
#'
#' Groups the atoms into subwindows (the `size` grid of each chromosome cut at every region start / end, so that each
#' region is a union of subwindows) and sums the votes per cell x subwindow x phase set: one unit is a stretch of at
#' most `size` bp in one phase set, inside which the phase is trusted.
#' @param votes votes table ([read_votes()])
#' @param atoms atoms table
#' @param regions regions table (the subwindows are cut at their boundaries; only their chromosomes are used)
#' @param size window size (bp)
#' @return list `units` (cell, sw, nA, n; units with n >= 1) and `windows` (sw, chrom, start, end)
#' @export
window_units = function(votes, atoms, regions, size = 1e6) {
  at = as.data.table(atoms); regions = as.data.table(regions)
  at = at[chrom %in% unique(regions$chrom)]
  sw_of = rep(NA_integer_, max(c(as.data.table(atoms)$atom, 0L)))
  wins = list(); off = 0L
  for (ch in unique(at$chrom)) {
    a = at[chrom == ch]
    hi = max(a$end)
    rr = regions[chrom == ch]
    cut = sort(unique(c(seq(0, hi, by = size), rr$start, rr$end[!is.na(rr$end)])))
    cut = cut[cut < hi]
    k = findInterval(a$start, cut)
    used = sort(unique(k))
    sw_of[a$atom] = off + match(k, used)
    wins[[ch]] = data.table(sw = off + seq_along(used), chrom = ch, start = cut[used], end = c(cut, hi)[used + 1])
    off = off + length(used)
  }
  w = rbindlist(wins)
  if (!nrow(w)) w = data.table(sw = integer(), chrom = character(), start = numeric(), end = numeric())
  v = as.data.table(votes)
  v = v[atom <= length(sw_of)]
  sw = sw_of[v$atom]
  ok = !is.na(sw)
  u = data.table(cell = v$cell[ok], sw = sw[ok], ps = v$ps[ok], nA = v$nA[ok], n = v$nA[ok] + v$nB[ok])
  u = u[, .(nA = sum(nA), n = sum(n)), by = .(cell, sw, ps)][n > 0]
  u[, ps := NULL]
  list(units = u[], windows = w[])
}

#' Windowed, orientation-free allele features per cell and region
#'
#' For each unit u of [window_units()] (one cell x subwindow of at most `params$win_size` bp x phase set) and each
#' fraction f of [win_f_grid()] (default 0.98: near-complete LOH), the mirrored beta-binomial LLR [mirrored_llr()]
#' `L_u(f) = log(0.5 BB(a_u; n_u, f) + 0.5 BB(a_u; n_u, 1 - f)) - log BB(a_u; n_u, 1/2)` (overdispersion `rho`,
#' allele error floor `eps`), summed over the units of each region: `W(f) = sum_u L_u(f)`. The phase is trusted only
#' inside a unit, so neither the orientation of the cell's haplotypes relative to other cells nor phase switches
#' between windows matter: an allelic imbalance private to one cell is visible on chromosome-scale (hybrid) phasing,
#' where the oriented chromosome sums of [allelic_features()] cancel it at population-phase switches. A unit with one
#' vote contributes exactly 0, and under the null model `E[exp W(f)] = 1`, so `P(W(f) > t) <= exp(-t)` (the null-cell
#' guard of [select_null_cells()]; the bound assumes independent votes inside a window, which holds well at f near 1
#' but not at gain fractions, see [win_f_grid()]).
#'
#' Needs phase sets that are reliable over a window: on block-only phasing a window holds many short unoriented
#' blocks and the units are blocks of one or two molecules (the per-block mirrored statistic, which has no power).
#' Subwindows whose votes over all cells have a shared allelic bias (`|f - 0.5| > params$win_bias_max` on at least
#' `params$win_bias_min_n` votes; germline CNV, mapping bias, imprinting) are left out. Only autosomes are scored
#' unless `params$win_chrx` (female chrX carries X-inactivation driven imbalance).
#' @param votes votes table
#' @param atoms atoms table
#' @param regions regions table (`allelic_ok` regions are used)
#' @param params [caller_defaults()] (`eps`, `win_f`, `win_size`, `win_bias_max`, `win_bias_min_n`, `win_chrx`, `seed`)
#' @param rho NULL (fitted: [fit_bb_rho()] on the units with n >= 2 of the null regions, at most 2e6 seeded units) or
#'   a value
#' @param exclude_atoms NULL or atom ids whose votes are dropped (allelic mask)
#' @return data.table `cell` (index), `region`, `n_win` (units with n >= 2), one `wal_f<fff>` column per fraction of
#'   [win_f_grid()] ([win_col()]); regions without a unit of the cell are absent. Attributes: `rho`, `windows` (the
#'   subwindows with their pooled `nA`, `n`, `f` and `masked`), `f_grid`, `regions` (the regions in scope)
#' @export
window_allelic_features = function(votes, atoms, regions, params = caller_defaults(), rho = NULL, exclude_atoms = NULL) {
  regions = as.data.table(regions)
  if ("allelic_ok" %in% names(regions)) regions = regions[allelic_ok %in% TRUE]
  chroms = c(AUTOSOMES_CALLER, if (isTRUE(params$win_chrx)) "chrX")
  regions = regions[chrom %in% chroms]
  fg = win_f_grid(params)
  cols = win_col(fg)
  finish = function(o, rho, w = NULL) {
    setattr(o, "rho", rho); setattr(o, "windows", w); setattr(o, "f_grid", fg); setattr(o, "regions", regions$region)
    o[]
  }
  empty = function() {
    o = data.table(cell = integer(), region = character(), n_win = numeric())
    for (k in cols) set(o, j = k, value = numeric())
    o
  }
  votes = as.data.table(votes)
  if (length(exclude_atoms)) votes = votes[!atom %in% exclude_atoms]
  if (!nrow(regions) || !nrow(votes)) return(finish(empty(), NA_real_))
  wu = window_units(votes, atoms, regions, size = params$win_size %||% 1e6)
  u = wu$units; w = wu$windows
  pw = u[, .(nA = sum(nA), n = sum(n)), by = sw]
  w[pw, `:=`(nA = i.nA, n = i.n), on = "sw"]
  w[is.na(n), `:=`(nA = 0, n = 0)]
  w[, f := fifelse(n > 0, nA / n, NA_real_)]
  w[, masked := !is.na(f) & n >= (params$win_bias_min_n %||% 50) & abs(f - 0.5) > (params$win_bias_max %||% 0.1)]
  u = u[n >= 2 & !sw %in% w$sw[w$masked]]
  # subwindow -> region (start-inside rule) per layer of non-overlapping regions
  map = atom_region_map(w[, .(atom = sw, chrom, start, end)], regions)
  if (is.null(rho)) {
    nullsw = if ("null_region" %in% names(regions)) unique(map$atom[regions$null_region[map$ri] %in% TRUE]) else w$sw
    s = u[sw %in% nullsw]
    if (!nrow(s)) s = u
    set.seed(params$seed %||% 1)
    s = s[sample.int(.N, min(.N, 2e6))]
    rho = if (nrow(s)) fit_bb_rho(s$nA, s$n) else NA_real_
    rm(s)
  }
  if (!nrow(u) || is.na(rho)) return(finish(empty(), rho, w))
  eps = params$eps
  L = vapply(fg, function(f) mirrored_llr(u$nA, u$n, f, rho, eps), numeric(nrow(u)))
  if (!is.matrix(L)) L = matrix(L, nrow = nrow(u))
  colnames(L) = cols
  out = lapply(sort(unique(map$layer)), function(l) {
    ml = map[layer == l]
    ri = ml$ri[match(u$sw, ml$atom)]
    sel = which(!is.na(ri))
    if (!length(sel)) return(NULL)
    agg = data.table(cell = u$cell[sel], ri = ri[sel], L[sel, , drop = FALSE])
    agg[, c(list(n_win = as.numeric(.N)), lapply(.SD, sum)), by = .(cell, ri), .SDcols = cols]
  })
  o = rbindlist(out)
  if (!nrow(o)) return(finish(empty(), rho, w))
  o[, region := regions$region[ri]][, ri := NULL]
  o[, cell := as.integer(cell)]
  setcolorder(o, c("cell", "region", "n_win", cols))
  finish(o, rho, w)
}

# whether the windowed allele term is computed for a features run: `params$win_allelic` "auto" (= chromosome-scale
# orientation), TRUE (always) or FALSE (never)
win_allelic_on = function(params, orientation) {
  wa = params$win_allelic %||% "auto"
  if (isFALSE(wa) || identical(tolower(as.character(wa)), "false")) return(FALSE)
  if (isTRUE(wa) || identical(tolower(as.character(wa)), "true")) return(TRUE)
  identical(unname(orientation), "chromosome")
}

# writes the window columns into the long features table (by reference): regions in the window scope get 0 where the
# cell has no unit (no evidence), regions out of scope NA (the scorer falls back to the oriented term there)
add_window_columns = function(long, wa, cells) {
  cols = c("n_win", win_col(attr(wa, "f_grid")))
  scope = long$region %in% attr(wa, "regions")
  k = match(paste(long$cell, long$region), paste(cells[wa$cell], wa$region))
  for (cn in cols) {
    v = as.numeric(wa[[cn]])[k]
    set(long, j = cn, value = fifelse(scope, fifelse(is.na(v), 0, v), NA_real_))
  }
  invisible(long)
}

# drops the window columns of a long features table (by reference)
drop_window_columns = function(long) {
  cols = intersect(c("n_win", grep("^wal_f[0-9]+$", names(long), value = TRUE)), names(long))
  if (length(cols)) long[, (cols) := NULL]
  invisible(long)
}

# the note of the windowed allele term of a features run
win_allelic_note = function(wa, params) {
  w = attr(wa, "windows")
  sprintf("windowed allele term: %s bp units, rho %.4g, f %s; %d of %d subwindows masked for shared bias (|f - 0.5| > %s on >= %s votes); %s",
          format(params$win_size %||% 1e6, scientific = FALSE), attr(wa, "rho"), paste(round(attr(wa, "f_grid"), 3), collapse = "/"),
          if (is.null(w)) 0L else sum(w$masked), if (is.null(w)) 0L else nrow(w), format(params$win_bias_max %||% 0.1),
          format(params$win_bias_min_n %||% 50), if (isTRUE(params$win_chrx)) "autosomes and chrX" else "autosomes only")
}

# the guard statistic per row: the largest window LLR over the fraction grid (-Inf where no window column / NA)
win_guard_stat = function(d, params = caller_defaults()) {
  wc = intersect(win_col(win_f_grid(params)), names(d))
  if (!length(wc)) return(NULL)
  W = lapply(wc, function(k) fifelse(is.na(d[[k]]), -Inf, as.numeric(d[[k]])))
  list(stat = if (length(W) == 1) W[[1]] else do.call(pmax, W), n_f = length(wc))
}
