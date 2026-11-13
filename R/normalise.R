# normalise.R
# What: one interface normalise_depth(scn, method = "loo" | "celltype" | "model") returning
#       an `scn_norm` (log2 ratio bins x cells, expected counts, size factors).
#       loo: global pseudobulk leave-one-out (LOO); celltype: cell-type-matched LOO with
#       minority fallback (the normal-tissue reference); model: robust Poisson GLM
#       on GC / mappability proxies / peak coverage with a per-cell GC trend (lowess, or ratio of sums
#       over GC groups) (RIDDLER-style; the tumour reference). Also the high-pass-filtered
#       chromatin embedding used for the caller's kNN neighbourhoods.
# Why: the depth channel's chromatin wave is cell-type specific and the noise
#      is cell-private; every reference for cell i excludes cell i (no leak).
# Inputs: scn (counts, bins, cells with optional cell_type column).
# Outputs: scn_norm list(method, ratio, expected, use, size, params).

#' Normalise read depth against a reference
#'
#' Expected counts per bin and cell, and the log2 ratio of observed to expected (centred per cell). Every
#' reference of a cell leaves that cell out.
#'   * `loo`: the pseudobulk of all other reference cells.
#'   * `celltype`: the pseudobulk of the other reference cells of the same cell type; cells without a type, and
#'     types with fewer than `min_cells_type` reference cells, use `loo`. The normal-tissue reference: it removes
#'     the cell-type-specific chromatin signal in the depth.
#'   * `model`: a robust Poisson GLM of the bin counts on GC content, mappability proxies and peak coverage,
#'     plus a per-cell GC trend (`per_cell_gc`, estimated as `gc_method` says). The tumour reference: it needs no
#'     copy-number-normal cells.
#' @param scn scn object
#' @param method "loo", "celltype" or "model"
#' @param use logical bin filter (default [bin_filter()] on bins and counts)
#' @param cell_type cell type per cell (default `scn$cells$cell_type`)
#' @param ref_cells optional selector of cells allowed as references (e.g. clean singlets)
#' @param min_cells_type celltype: types with fewer reference cells fall back to the global reference
#' @param per_cell_gc model: subtract a per-cell robust GC trend after the shared GLM
#' @param gc_method model with `per_cell_gc`: estimator of the per-cell GC trend. `"lowess"` (default): a robust
#'   lowess of the bin log2 ratio `log2((x + 0.5) / (e + 0.5))` on bin GC. At low depth (most bins with 0 or 1
#'   fragments) the mean of that log ratio depends nonlinearly on `e`, so this fit finds a GC trend that depends on
#'   the cell's depth even in a cell without one. `"ratio"`: the ratio of summed observed to summed expected counts
#'   in 20 GC quantile groups of the bins, smoothed across the groups (weighted local-linear loess, span 0.75) and
#'   evaluated at each bin's GC; unbiased at every depth. Bins without GC keep the shared model; a cell
#'   without fragments gets no trend. ([caller_defaults()] `gc_method`)
#' @param ... unused
#' @return `scn_norm`: list(method, ratio (bins x cells log2, NA on unused bins),
#'   expected (bins x cells counts, NA on unused bins), use, size, params)
#' @export
normalise_depth = function(scn, method = "loo", use = NULL, cell_type = NULL, ref_cells = NULL,
                           min_cells_type = 20, per_cell_gc = TRUE, gc_method = "lowess", ...) {
  method = norm_method(method)
  gc_method = gc_trend_method(gc_method)
  use = use %||% bin_filter(scn$bins, scn$counts)
  types = cell_type %||% scn$cells$cell_type
  refs = if (is.null(ref_cells)) rep(TRUE, ncol(scn$counts)) else seq_len(ncol(scn$counts)) %in% cell_selector(scn, ref_cells)
  X = dense(scn$counts[use, , drop = FALSE])
  s = colSums(X)
  bins_use = scn$bins[use]
  E = switch(method,
    loo = norm_loo(X, s, refs),
    celltype = norm_celltype(X, s, types, refs, min_cells_type),
    model = norm_model(X, s, bins_use, per_cell_gc, gc_method)
  )
  finish_norm(scn, X, E, use, s, method, list(min_cells_type = min_cells_type, per_cell_gc = per_cell_gc,
                                              gc_method = gc_method))
}

# depth normalisation method: "loo", "celltype" or "model" (aliases N1, N2, N4)
norm_method = function(method) {
  alias = c(N1 = "loo", N2 = "celltype", N4 = "model")
  method = method[1]
  if (method %in% names(alias)) method = alias[[method]]
  match.arg(method, c("loo", "celltype", "model"))
}

# estimator of the per-cell GC trend of `model`: "lowess" or "ratio" (NULL = "lowess", e.g. parameters of older runs)
gc_trend_method = function(gc_method) {
  if (is.null(gc_method) || !length(gc_method)) return("lowess")
  match.arg(as.character(gc_method[1]), c("lowess", "ratio"))
}

finish_norm = function(scn, X, E, use, s, method, params) {
  # rescale expected to the cell's total over used bins, then centre the log2 ratio per cell
  E = sweep(E, 2, colSums(E) / s, "/")
  E[!is.finite(E)] = NA
  ratio = log2((X + 0.5) / (E + 0.5))
  ratio = sweep(ratio, 2, apply(ratio, 2, median, na.rm = TRUE), "-")
  full_r = matrix(NA_real_, nrow(scn$counts), ncol(scn$counts)); full_r[use, ] = ratio
  full_e = matrix(NA_real_, nrow(scn$counts), ncol(scn$counts)); full_e[use, ] = E
  colnames(full_r) = colnames(full_e) = scn$cells$cell
  structure(list(method = method, ratio = full_r, expected = full_e, use = use, size = s, params = params,
                 cells = scn$cells$cell), class = "scn_norm")
}

#' @export
print.scn_norm = function(x, ...) {
  cat(sprintf("scn_norm %s: %d used bins x %d cells\n", x$method, sum(x$use), length(x$cells)))
  invisible(x)
}

lowess_fit = function(x, y, f = 0.3) {
  ok = is.finite(x) & is.finite(y)
  out = rep(NA_real_, length(x))
  if (sum(ok) < 20) { out[ok] = median(y[ok]); return(out) }
  lw = lowess(x[ok], y[ok], f = f, iter = 3)
  # lowess returns sorted unique-ish x; interpolate back
  out[ok] = approx(lw$x, lw$y, xout = x[ok], rule = 2, ties = mean)$y
  out
}

# ---- loo: global pseudobulk LOO --------------------------------------------------------------
norm_loo = function(X, s, refs) {
  Tb = rowSums(X[, refs, drop = FALSE]); S = sum(s[refs])
  # for a reference cell: (T - x_i) / (S - s_i); for a non-reference cell: T / S
  own = X * rep(refs, each = nrow(X))
  E = (Tb - own) / rep(S - s * refs, each = nrow(X))
  E = E * rep(s, each = nrow(X))
  E
}

# ---- celltype: cell-type matched LOO with minority fallback -------------------------------------
norm_celltype = function(X, s, types, refs, min_cells_type = 20) {
  if (is.null(types)) return(norm_loo(X, s, refs))
  E = norm_loo(X, s, refs)
  for (t in setdiff(unique(types), NA)) {          # NA type = no matched reference (global fallback)
    inT = !is.na(types) & types == t
    rt = inT & refs
    if (sum(rt) < min_cells_type) next
    Tb = rowSums(X[, rt, drop = FALSE]); S = sum(s[rt])
    cols = which(inT)
    own = X[, cols, drop = FALSE] * rep(refs[cols], each = nrow(X))
    Et = (Tb - own) / rep(S - s[cols] * refs[cols], each = nrow(X))
    E[, cols] = Et * rep(s[cols], each = nrow(X))
  }
  E
}

# ---- high-pass-filtered chromatin embedding ------------------------------------------------
high_pass = function(R, bins, window = 21) {
  if (window %% 2 == 0) window = window + 1
  out = R
  for (ch in unique(bins$chrom)) {
    rows = which(bins$chrom == ch)
    if (length(rows) < window) { out[rows, ] = R[rows, , drop = FALSE] - rep(colMeans(R[rows, , drop = FALSE]), each = length(rows)); next }
    for (i in seq_len(ncol(R))) {
      v = R[rows, i]
      out[rows, i] = v - stats::runmed(v, window, endrule = "median")
    }
  }
  out
}

#' High-pass-filtered chromatin embedding (kNN neighbourhoods of the caller)
#' @param X used bins x cells counts
#' @param s size factors
#' @param bins used bins table
#' @param hp_window running-median window (bins)
#' @param n_pc number of PCs
#' @return cells x n_pc score matrix
#' @export
chromatin_embedding = function(X, s, bins, hp_window = 21, n_pc = 20) {
  E1 = norm_loo(X, s, rep(TRUE, ncol(X)))
  R = log2((X + 0.5) / (E1 + 0.5))
  R = high_pass(R, bins, hp_window)
  R = R - rowMeans(R)
  n_pc = min(n_pc, ncol(X) - 1, nrow(X) - 1)
  sv = irlba::irlba(t(R), nv = n_pc)
  sv$u %*% diag(sv$d, n_pc, n_pc)
}

# ---- model: robust Poisson GLM on bin covariates (RIDDLER-style) ------------------------------
norm_model = function(X, s, bins, per_cell_gc = TRUE, gc_method = "lowess", gc_groups = 20, gc_span = 0.75) {
  Tb = rowSums(X)
  cov = data.frame(y = Tb)
  form = "y ~ ns(gc, 4)"
  cov$gc = bins$gc
  informative = function(v) is.numeric(v) && sum(is.finite(v)) > 10 && sd(v, na.rm = TRUE) > 1e-6
  if ("peak_bp" %in% names(bins) && informative(bins$peak_bp)) { cov$peak = log1p(bins$peak_bp); form = paste(form, "+ peak") }
  if ("frac_trf" %in% names(bins) && informative(bins$frac_trf)) { cov$trf = bins$frac_trf; form = paste(form, "+ trf") }
  if ("frac_blacklist" %in% names(bins) && informative(bins$frac_blacklist) && mean(bins$frac_blacklist > 0) > 0.05) { cov$bl = bins$frac_blacklist; form = paste(form, "+ bl") }
  # glmrob's Mqle often stops at its iteration limit on ~30k bins with a "did not converge" warning; the last iterate
  # is used either way (as before). That one warning is muffled and reported once as a message; other warnings pass.
  nonconv = FALSE
  fit = tryCatch(withCallingHandlers(
                   robustbase::glmrob(stats::as.formula(form), family = stats::poisson(), data = cov, method = "Mqle"),
                   warning = function(w) if (grepl("did not converge", conditionMessage(w), fixed = TRUE)) {
                     nonconv <<- TRUE
                     invokeRestart("muffleWarning")
                   }),
                 error = function(e) NULL)
  if (nonconv && !is.null(fit)) message("norm_model: robustbase::glmrob did not converge (iteration limit); using its last iterate")
  mu = if (is.null(fit)) {
    fit2 = stats::glm(stats::as.formula(form), family = stats::poisson(), data = cov)
    predict(fit2, type = "response")
  } else predict(fit, type = "response")
  mu = pmax(mu, 1e-9)
  E = outer(mu / sum(mu), s)
  if (per_cell_gc && gc_trend_method(gc_method) == "ratio") {
    # per-cell GC trend from aggregated counts (gc_group_trend()), applied column by column so that no second
    # bins x cells matrix is allocated
    tr = gc_group_trend(X, E, bins$gc, mu, n_groups = gc_groups, span = gc_span)
    if (!is.null(tr)) {
      for (i in seq_len(ncol(X))) {
        f = tr$F[tr$lo, i] * (1 - tr$w) + tr$F[tr$hi, i] * tr$w
        E[tr$rows, i] = E[tr$rows, i] * 2^f
      }
    }
  } else if (per_cell_gc) {
    # per-cell robust GC trend of the residual: lowess (C, ~20 ms per cell) in the parent process; no forks, so the
    # dense matrices are never duplicated (the rlm + mclapply version was OOM-killed at 120 GB on 10k cells)
    for (i in seq_len(ncol(X))) {
      r = log2((X[, i] + 0.5) / (E[, i] + 0.5))
      fit = lowess_fit(bins$gc, r, f = 0.3)
      ok = is.finite(fit)
      E[ok, i] = E[ok, i] * 2^fit[ok]
    }
  }
  E
}

# Per-cell GC trend by ratio of sums (norm_model(gc_method = "ratio")). The bins with GC are split into `n_groups` GC
# quantile groups, shared by all cells. Per group and cell: summed counts over summed expected counts (0.5 added to
# both sums, so a cell without fragments gets 1). Sums are unbiased at any depth, unlike the mean of per-bin log
# ratios of sparse counts; the ratios are smoothed before the log is taken, because the log of a single group's ratio
# is itself biased at low depth by an amount that depends on the group's expected count, i.e. on GC. Each group sits
# at its `mu`-weighted mean GC. The group ratios are smoothed by a
# local-linear loess (`span`, surface "direct") weighted by the group's share of the shared model `mu`, evaluated on
# a grid of `n_grid` GC values over the range of the bins (so the outer bins get the local-linear extrapolation, not
# the outer group's value), and interpolated linearly in GC from the grid to the bins. The weights are the same for
# every cell, so the smoother is one n_grid x G linear operator applied to all cells at once (loess is linear in the
# response: column j is the fit to the j-th unit vector). With fewer than 8 groups the group values are only
# interpolated (constant beyond the outer groups). The smoothed ratio is floored at 1/16 before the log (an
# extrapolated ratio of a nearly empty cell can reach 0). Returns NULL when there are fewer than 2 groups (constant
# GC), else list(rows = the bins with GC, lo / hi / w = grid points and weight per such bin, F = grid x cells log2
# factors).
gc_group_trend = function(X, E, gc, mu, n_groups = 20, span = 0.75, n_grid = 101) {
  rows = which(is.finite(gc))
  if (length(rows) < 2) return(NULL)
  br = unique(stats::quantile(gc[rows], seq(0, 1, length.out = n_groups + 1), names = FALSE))
  if (length(br) < 2) return(NULL)
  g = as.integer(cut(gc[rows], br, include.lowest = TRUE))
  g = match(g, sort(unique(g)))                                   # drop empty groups (ties)
  G = max(g)
  if (G < 2) return(NULL)
  grp = rep(G + 1L, nrow(X)); grp[rows] = g                       # G + 1: bins without GC, dropped below
  sx = rowsum(X, grp, reorder = TRUE)[seq_len(G), , drop = FALSE] # groups x cells
  se = rowsum(E, grp, reorder = TRUE)[seq_len(G), , drop = FALSE]
  Q = (sx + 0.5) / (se + 0.5)
  m = pmax(mu[rows], 0)
  wg = as.vector(tapply(m, g, sum))
  mg = as.vector(tapply(m * gc[rows], g, sum)) / wg
  mg[!is.finite(mg)] = as.vector(tapply(gc[rows], g, mean))[!is.finite(mg)]
  wg = wg / sum(wg)
  grid = seq(min(gc[rows]), max(gc[rows]), length.out = n_grid)
  S = vapply(seq_len(G), function(j) {
    y = as.numeric(seq_len(G) == j)
    if (G >= 8) {
      fit = stats::loess(y ~ mg, weights = wg, span = span, degree = 1, control = stats::loess.control(surface = "direct"))
      as.numeric(stats::predict(fit, newdata = data.frame(mg = grid)))
    } else stats::approx(mg, y, xout = grid, rule = 2)$y
  }, numeric(n_grid))
  F = log2(pmax(S %*% Q, 1 / 16))
  lo =findInterval(gc[rows], grid, all.inside = TRUE)
  w = (gc[rows] - grid[lo]) / (grid[lo + 1L] - grid[lo])
  list(rows = rows, lo = lo, hi = lo + 1L, w = pmin(pmax(w, 0), 1), F = F)
}
