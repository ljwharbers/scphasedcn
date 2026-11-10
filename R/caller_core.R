# caller_core.R
# What: the shared foundation of the copy-number caller: default parameters, the depth stratum, the two channel
#       log-likelihood primitives (depth location shift under a Student-t or
#       Gaussian null; oriented beta-binomial allele LLR vs 1:1), the null fits
#       (location-scale t per class x stratum; beta-binomial rho), schema
#       validators for the tables the caller modules exchange, the `scphasedcn`
#       result object and the tile naming rule.
# Why: the feature step (features.R), the scorer (caller.R), the tumour mode (clones_*.R) and the
#      pipeline (call_scn.R) share these definitions so that every part scores, thresholds and names
#      regions the same way. The primitives are tested for exact parity with an independent reference
#      implementation of the scorer in the tests (tests/testthat/helper-refscorer.R).
# Inputs: vectors / data.tables (schemas in vignette("scphasedcn-outputs", package = "scphasedcn")).
# Outputs: vectors, data.tables, `scphasedcn` objects.

utils::globalVariables(c("z", "t_m", "t_s", "t_df", "snp1", "rho", "n_rows", "nA_", "n_", "ok_", "weight", "stat", "stat_min"))

# ---- defaults --------------------------------------------------------------------------------------

#' Default parameters of the copy-number caller
#'
#' One named list shared by all parts of the caller. The calibrated values (`shifts`, `null_q`, the depth
#' strata) were set on SPLONGGET normal-tissue and tumour data; see the notes per entry.
#' @return named list:
#'   * `null_q` threshold quantile of the null score (0.999, type 7: 1 false call per 1,000 null cell-regions);
#'     calls need score `>` threshold (strict, not `>=`);
#'   * `eps` allele error floor of the beta-binomial (0.02);
#'   * `mol_gap` fragments of one cell closer than this (bp) are one molecule (2000; the default of
#'     [cell_features()] `mol_gap`);
#'   * `sd_floor` floor of the per-cell null SD of the window log2 ratio (0.05);
#'   * `strata_breaks` / `strata_labels` depth strata on clean fragments per cell
#'     (`<20k`, `20-50k`, `>=50k`);
#'   * `shifts` log2 depth shifts of loss (-0.98) and gain (+0.57) (calibrated on depth-only spike-ins in normal tissue);
#'   * `k_amp` uncalibrated depth shift factor for cn >= 4: shift = `k_amp * log2(cn / 2)`;
#'   * `f_att` attenuated major-haplotype fractions of a 1:0 / 2:0 state;
#'   * `tile_sizes` named tile sizes (bp);
#'   * `n_folds` cross-validation folds of the null thresholds (4);
#'   * `n_pool` cells per pool (20);
#'   * `chrx_centromere` chrX centromere (bp, hg38) for female chrX arms;
#'   * `acrocentric_p` arms without unique sequence, never regions;
#'   * `min_null_rows` fewer null rows than this in a class x stratum fall back to the pooled threshold (5000);
#'   * `tabix` path of the tabix binary ("" = the `tabix` on the PATH; none = the R fallback); `seed` RNG seed (1);
#'   * `calib` per-cell calibration of normal mode: `"class"` (one Student-t null and one threshold per
#'     class x stratum, the default design) or `"votes"` (opt-in: thresholds per class x stratum x vote stratum
#'     `hstratum`, bins of the allele votes n_hap; [vote_bins()]);
#'   * `vote_bins` maximal number of n_hap bins per class x stratum of the vote calibration (4; fewer when
#'     the class x stratum has fewer than `vote_bins x min_null_rows` null rows with votes);
#'   * `calib_rho` allelic overdispersion of the normal-mode scores: `"global"` (one rho on the null arm
#'     rows) or `"chrom"` (opt-in: per chromosome, shrunk to the global; [rho_by_chrom()]);
#'   * `rho_chrom_k` shrinkage pseudo-count of the per-chromosome rho (100 null arm rows; a chromosome with ~4,000
#'     null rows keeps ~98 % of its own estimate);
#'   * `mask_allelic` data-driven allelic tile mask in normal mode ([allelic_noisy_tiles()], `call_scn()`): the
#'     votes of 5 Mb tiles whose null overdispersion is an outlier are left out of the allele features (FALSE);
#'   * `mask_rho_mad` mask threshold in MADs of log tile rho (2);
#'   * `mask_min_votes` minimal mean votes per null row for a tile to be assessed (3);
#'   * `pool_sex_thr` males: the sex chromosomes get their own pooled threshold group (class x autosome / X / Y)
#'     with the class threshold as fallback (TRUE; no effect in females);
#'   * `chrx_allelic` normal mode, females: score chrX rows with the allele channel (TRUE, the default) or on
#'     depth only (FALSE: [caller_normal()] drops the chrX allele votes before scoring and treats the rows like
#'     male chrX, `allelic_ok = FALSE`, so purely allelic states such as CN-LOH are not scored there, and
#'     thresholds them with depth-only thresholds of the autosomal null rows, [depth_only_thresholds()]).
#'     X-inactivation makes every female cell allelically imbalanced on chrX, so with TRUE chrX collects
#'     allele-driven "gains" / CN-LOH without depth support. The null cells, the null fits, the default
#'     thresholds and every autosomal score are the same either way; chrX pools of a female stay non-reportable
#'     ([call_policy()]). No effect in males, on the features, or in tumour mode (which drops the chrX votes
#'     itself, [clone_defaults()] `allelic_exclude`);
#'   * `pool_tail_q` pooled per-region tail calibration ([pool_scale()]): quantile of the region-standardised
#'     null pools compared with the class (0 = off; e.g. 0.99);
#'   * `pool_tail_k` shrinkage pseudo-count of the pooled tail factor (30 effective pools);
#'   * `pool_tail_xfit` cross-fit the pooled tail ([pooled_thresholds()]): the region scale and the
#'     threshold are estimated on pools of two disjoint halves of the null cells (TRUE).
#' @export
caller_defaults = function() {
  list(
    null_q = 0.999,
    eps = 0.02,
    mol_gap = 2000,
    sd_floor = 0.05,
    strata_breaks = c(2e4, 5e4),
    strata_labels = c("<20k", "20-50k", ">=50k"),
    shifts = c(loss = -0.98, gain = 0.57),
    k_amp = 0.975,
    f_att = c(0.85, 0.90, 0.95, 0.98),
    tile_sizes = c("5Mb" = 5e6, "10Mb" = 1e7, "20Mb" = 2e7, "50Mb" = 5e7),
    n_folds = 4,
    n_pool = 20,
    chrx_centromere = c(58605580, 62412542),
    acrocentric_p = paste0("chr", c(13, 14, 15, 21, 22), "p"),
    min_null_rows = 5000,
    tabix = "",
    seed = 1,
    calib = "class",
    vote_bins = 4,
    calib_rho = "global",
    rho_chrom_k = 100,
    mask_allelic = FALSE,
    mask_rho_mad = 2,
    mask_min_votes = 3,
    pool_sex_thr = TRUE,
    chrx_allelic = TRUE,
    pool_tail_q = 0.99,
    pool_tail_k = 30,
    pool_tail_xfit = TRUE
  )
}

#' Depth stratum of a cell
#'
#' Strata on clean fragments per cell:
#' `cut(n, c(-Inf, breaks, Inf), right = FALSE)`, so a cell with exactly
#' 20,000 fragments is in `20-50k`.
#' @param n clean fragments per cell
#' @param breaks interior breaks (default `caller_defaults()$strata_breaks`)
#' @param labels stratum labels (length(breaks) + 1); generated as `<a`,
#'   `a-b`, `>=b` (in thousands) when NULL
#' @return factor
#' @export
depth_stratum = function(n, breaks = caller_defaults()$strata_breaks, labels = NULL) {
  if (is.null(labels)) {
    k = as.character(breaks / 1e3)
    labels = c(paste0("<", k[1], "k"),
               if (length(k) > 1) paste0(k[-length(k)], "-", k[-1], "k"),
               paste0(">=", k[length(k)], "k"))
  }
  if (length(labels) != length(breaks) + 1) stop("depth_stratum: need length(breaks) + 1 labels")
  cut(n, c(-Inf, breaks, Inf), labels = labels, right = FALSE)
}

# ---- channel log-likelihoods -------------------------------------------------------------------------

#' Depth log-likelihood ratio of a location shift
#'
#' `log f0(z - mu) - log f0(z)`: the evidence that a cell-region's window z
#' (`log2ratio / pmax(sd_null, sd_floor)`) is shifted by `mu` (a state's log2
#' shift divided by the same floored `sd_null`). `f0` is the location-scale
#' Student-t fitted to the null z of the row's class x stratum when `tpar` is
#' given (`dt((x - t_m) / t_s, t_df)`; [fit_depth_null()]), else N(0, 1)
#' (`mu * z - mu^2 / 2`). Missing z contributes 0 (no evidence). This is
#' the depth part (`dl_loss`, `dl_gain`) of the region score of [score_cells()].
#' @param z window z per row (NA allowed)
#' @param mu shift in z units (scalar or one per row)
#' @param tpar NULL (Gaussian) or the t parameters: a list / 1-row data.table
#'   with `t_m`, `t_s`, `t_df`, or the same columns with one value per row
#' @return numeric vector (length of z)
#' @export
ll_depth_shift = function(z, mu, tpar = NULL) {
  zz = fifelse(is.na(z), 0, z)
  if (is.null(tpar)) {
    out = mu * zz - mu^2 / 2
  } else {
    if (!all(c("t_m", "t_s", "t_df") %in% names(tpar))) stop("ll_depth_shift: tpar needs t_m, t_s, t_df")
    m = tpar[["t_m"]]; s = tpar[["t_s"]]; df = tpar[["t_df"]]
    for (v in list(m, s, df)) if (!length(v) %in% c(1L, length(z))) stop("ll_depth_shift: tpar must have 1 row or one row per z")
    lt = function(x) dt((x - m) / s, df, log = TRUE)
    out = lt(zz - mu) - lt(zz)
  }
  fifelse(is.na(z), 0, out)
}

#' Oriented allelic log-likelihood ratio vs balanced
#'
#' `bb_loglik(nA, n, f, rho, eps) - bb_loglik(nA, n, 0.5, rho, eps)` for
#' molecule votes `nA` of `n` (one vote per molecule). Two-sided (default)
#' takes the better of `f` and `1 - f` (the orientation of the counts is
#' unknown or may be wrong). With several `f` the best LLR over the grid is
#' returned (the attenuation-aware 1:0 / 2:0 term over `f_att`). Rows with
#' `n = 0` (or NA) give 0. The scorer uses llr_20 at
#' `f = 1 - eps`, llr_21 at 2/3, llr_31 at 3/4 and the grid form over `f_att`.
#' @param nA votes for haplotype A
#' @param n total votes
#' @param f major-haplotype fraction of the state (scalar or a grid)
#' @param rho beta-binomial overdispersion (null fit, [fit_allelic_null()])
#' @param eps allele error floor
#' @param two_sided max over f and 1 - f
#' @param per_row `f` is one expected fraction per row (length 1 or `length(n)`), not a grid: the LLR of
#'   each row at its own `f` (NA results also give 0). Used by the tumour mode (`0.5 + a * (f_state - 0.5)` per row).
#' @return numeric vector (length of n)
#' @export
ll_allelic = function(nA, n, f, rho, eps = 0.02, two_sided = TRUE, per_row = FALSE) {
  ll0 = bb_loglik(nA, n, 0.5, rho, eps)
  if (per_row) {
    if (!length(f) %in% c(1L, length(n))) stop("ll_allelic: with per_row = TRUE, f needs length 1 or length(n)")
    l = if (two_sided) pmax(bb_loglik(nA, n, f, rho, eps), bb_loglik(nA, n, 1 - f, rho, eps)) else bb_loglik(nA, n, f, rho, eps)
    out = l - ll0
    out[is.na(n) | n == 0 | is.na(out)] = 0
    return(out)
  }
  best = rep(-Inf, length(n))
  for (fk in f) {
    l = if (two_sided) pmax(bb_loglik(nA, n, fk, rho, eps), bb_loglik(nA, n, 1 - fk, rho, eps)) else bb_loglik(nA, n, fk, rho, eps)
    best = pmax(best, l - ll0)
  }
  best[is.na(n) | n == 0] = 0
  best
}

# ---- null fits -------------------------------------------------------------------------------------------

#' Student-t fit of the null depth z per group
#'
#' Location-scale Student-t (`MASS::fitdistr(x, "t")`) of the finite null z
#' per group, on a seeded subsample of at most `n_max` rows per group;
#' fallback median / mad / df 30 when the fit fails.
#' @param d null rows with `z` and the `by` columns
#' @param by grouping columns (region class x depth stratum)
#' @param n_max maximal rows per group used for the fit
#' @param seed RNG seed (re-set per group)
#' @return data.table with the `by` columns, `t_m`, `t_s`, `t_df`, `n` (finite null rows)
#' @export
fit_depth_null = function(d, by = c("class", "stratum"), n_max = 20000, seed = 1) {
  d = as.data.table(d)
  d[is.finite(z), {
    set.seed(seed)
    x = if (.N > n_max) z[sample.int(.N, n_max)] else z
    f = tryCatch(suppressWarnings(MASS::fitdistr(x, "t")$estimate), error = function(e) c(m = median(x), s = mad(x), df = 30))
    list(t_m = f[["m"]], t_s = f[["s"]], t_df = f[["df"]], n = .N)
  }, by = by]
}

#' Per-row Student-t null parameters
#'
#' Looks up the [fit_depth_null()] parameters of each row's group (default
#' class x stratum), as the scorer does; stops when a group has no fit.
#' @param d rows to score (with the `by` columns)
#' @param tpar output of [fit_depth_null()]
#' @param by grouping columns
#' @return data.table `t_m`, `t_s`, `t_df` with one row per row of `d` (pass to [ll_depth_shift()])
#' @export
match_depth_null = function(d, tpar, by = c("class", "stratum")) {
  key = function(x) do.call(paste, c(lapply(by, function(b) as.character(x[[b]])), sep = "\r"))
  k = match(key(d), key(tpar))
  if (anyNA(k)) stop(sprintf("match_depth_null: no t fit for %s", paste(head(unique(key(d)[is.na(k)]), 5), collapse = ", ")))
  data.table(t_m = tpar$t_m[k], t_s = tpar$t_s[k], t_df = tpar$t_df[k])
}

#' Beta-binomial overdispersion of the null allele votes
#'
#' [fit_bb_rho()] (balanced, f = 0.5) on the rows with `n > 0`.
#' @param nA votes for haplotype A (oriented)
#' @param n total votes
#' @param by NULL (one rho) or grouping: a vector with one label per row, or
#'   a list / data.table of such vectors
#' @return rho (numeric) when `by` is NULL; else data.table with the group
#'   columns (`group` for an unnamed vector), `rho` and `n_rows` (rows with n > 0;
#'   rho NA for a group without any)
#' @export
fit_allelic_null = function(nA, n, by = NULL) {
  ok = !is.na(n) & n > 0
  if (is.null(by)) {
    if (!any(ok)) return(NA_real_)
    return(fit_bb_rho(nA[ok], n[ok]))
  }
  g = if (is.atomic(by)) data.table(group = by) else as.data.table(by)
  if (nrow(g) != length(n)) stop("fit_allelic_null: `by` must have one value per row")
  gcols = copy(names(g))
  g[, `:=`(nA_ = nA, n_ = n, ok_ = ok)]
  out = g[, list(rho = if (any(ok_)) fit_bb_rho(nA_[ok_], n_[ok_]) else NA_real_, n_rows = sum(ok_)), by = gcols]
  out[]
}

# ---- schemas and validators --------------------------------------------------------------------------

#' Region classes of the caller
#' @return character vector: arm, chrom, 5Mb, 10Mb, 20Mb, 50Mb, event, segment
#' @export
region_classes = function() c("arm", "chrom", "5Mb", "10Mb", "20Mb", "50Mb", "event", "segment")

caller_schemas = function() {
  list(
    regions = c(region = "character", class = "character", sd_class = "character", chrom = "character",
                start = "numeric", end = "numeric", n_bins = "integer", n_use = "integer",
                eligible = "logical", null_region = "logical", allelic_ok = "logical"),
    votes = c(cell = "integer", atom = "integer", ps = "numeric", nA = "integer", nB = "integer"),
    atoms = c(atom = "integer", chrom = "character", start = "numeric", end = "numeric", bin = "integer"),
    features = c(cell = "character", region = "character", class = "character", sd_class = "character",
                 chrom = "character", x = "numeric", e = "numeric", log2ratio = "numeric", sd_null = "numeric",
                 z = "numeric", nA = "numeric", nB = "numeric", n_hap = "numeric", n_blocks = "numeric",
                 f = "numeric", orient = "character", stratum = "character", n_frag = "numeric",
                 cell_type = "character", fold = "integer", null_cell = "logical", null_region = "logical",
                 eligible = "logical"),
    thresholds = c(score = "character", class = "character", stratum = "character", fold = "integer",
                   q = "numeric", thr = "numeric", n_null = "integer", pooled_fallback = "logical"),
    calls = c(cell = "character", region = "character", class = "character", chrom = "character",
              start = "numeric", end = "numeric", state = "character", score = "numeric", thr = "numeric",
              called = "logical", level = "character", source = "character", reportable = "logical"),
    segments = c(clone = "character", chrom = "character", start = "numeric", end = "numeric",
                 n_win = "integer", nA = "integer", nB = "integer", cn = "integer", state = "character",
                 log2 = "numeric", f = "numeric", llr = "numeric", confident = "logical")
  )
}

# does a column satisfy a schema type? An all-NA column passes any type; "integer" also accepts whole doubles.
col_has_type = function(v, type) {
  if (all(is.na(v)) && (is.logical(v) || is.numeric(v) || is.character(v))) return(TRUE)
  switch(type,
         character = is.character(v),
         numeric = is.numeric(v),
         integer = is.integer(v) || (is.double(v) && all(is.na(v) | v == round(v))),
         logical = is.logical(v),
         stop("unknown schema type ", type))
}

check_schema = function(x, what, allowed = list()) {
  if (!is.data.frame(x)) stop(sprintf("check_%s: expected a data.frame / data.table, got %s", what, class(x)[1]), call. = FALSE)
  sch = caller_schemas()[[what]]
  miss = setdiff(names(sch), names(x))
  if (length(miss)) stop(sprintf("check_%s: missing column(s): %s", what, paste(miss, collapse = ", ")), call. = FALSE)
  bad = names(sch)[!vapply(names(sch), function(cn) col_has_type(x[[cn]], sch[[cn]]), logical(1))]
  if (length(bad)) {
    stop(sprintf("check_%s: wrong type for %s", what,
                 paste(sprintf("`%s` (expected %s, got %s)", bad, sch[bad],
                               vapply(bad, function(cn) class(x[[cn]])[1], character(1))), collapse = ", ")), call. = FALSE)
  }
  for (cn in names(allowed)) {
    v = x[[cn]]
    v = unique(v[!is.na(v)])
    out = setdiff(v, allowed[[cn]])
    if (length(out)) stop(sprintf("check_%s: `%s` has value(s) %s; allowed: %s", what, cn,
                                  paste(head(out, 5), collapse = ", "), paste(allowed[[cn]], collapse = "|")), call. = FALSE)
  }
  invisible(x)
}

#' Validate the caller's tables against their schemas
#'
#' Each validator stops with a message naming the missing column(s), the
#' column(s) of the wrong type, or values outside the allowed set, and
#' returns its input invisibly. Extra columns are allowed. Types: character
#' (factors are rejected), numeric (integer or double), integer (integer or
#' whole-number double), logical; an all-NA column passes any type. Schemas
#' are listed in `vignette("scphasedcn-outputs", package = "scphasedcn")`.
#' @param x table to check
#' @return `x`, invisibly
#' @name check_tables
NULL

#' @describeIn check_tables regions: `region` unique, `class` in [region_classes()]
#' @export
check_regions = function(x) {
  check_schema(x, "regions", list(class = region_classes()))
  if (anyDuplicated(x$region)) stop(sprintf("check_regions: duplicated region(s): %s",
                                            paste(head(unique(x$region[duplicated(x$region)]), 5), collapse = ", ")), call. = FALSE)
  if (any(x$end <= x$start, na.rm = TRUE)) stop("check_regions: end <= start (coordinates are 0-based half-open)", call. = FALSE)
  invisible(x)
}

#' @describeIn check_tables molecule votes: `cell` (index into the cells), `atom`, `ps`, `nA`, `nB`
#' @export
check_votes = function(x) {
  check_schema(x, "votes")
  if (any(x$nA < 0 | x$nB < 0, na.rm = TRUE)) stop("check_votes: negative vote counts", call. = FALSE)
  invisible(x)
}

#' @describeIn check_tables allele atoms: `atom`, `chrom`, `start`, `end`, `bin`
#' @export
check_atoms = function(x) {
  check_schema(x, "atoms")
  if (anyDuplicated(x$atom)) stop("check_atoms: duplicated atom ids", call. = FALSE)
  invisible(x)
}

#' @describeIn check_tables features long table (one row per cell x region)
#' @export
check_features = function(x) {
  check_schema(x, "features", list(class = region_classes()))
  invisible(x)
}

#' @describeIn check_tables null thresholds per score x class x stratum x fold
#' @export
check_thresholds = function(x) {
  check_schema(x, "thresholds")
  invisible(x)
}

#' @describeIn check_tables calls: `level` cell|pool|clone, `source` cell|inherited|confirmed
#' @export
check_calls = function(x) {
  check_schema(x, "calls", list(level = c("cell", "pool", "clone"), source = c("cell", "inherited", "confirmed")))
  invisible(x)
}

#' @describeIn check_tables clone segments
#' @export
check_segments = function(x) {
  check_schema(x, "segments")
  invisible(x)
}

# ---- result object ------------------------------------------------------------------------------------

#' The caller's result object
#'
#' An S3 list of class `scphasedcn` holding everything one run of the caller
#' produces. Components that are not NULL are validated with the matching
#' `check_*` function.
#' @param mode "normal" (normal tissue: per-cell calls, pooled small events)
#'   or "tumour" (clone discovery, clone states, cell re-assignment)
#' @param cells character vector of cell barcodes or a data.table with a
#'   `cell` column (row order defines the integer `cell` index of the votes)
#' @param regions regions table ([check_regions()])
#' @param features features long table ([check_features()])
#' @param scores per cell x region x state scores
#' @param thresholds thresholds table ([check_thresholds()])
#' @param calls calls table ([check_calls()])
#' @param events called events (the scorer / pipeline)
#' @param pools pool definitions
#' @param clones clone table
#' @param assignments cell-to-clone assignments
#' @param segments clone segments ([check_segments()])
#' @param profiles per-clone or per-cell profiles (the tumour mode / pipeline)
#' @param null_fit null fits (depth t parameters, allelic rho)
#' @param params parameter list ([caller_defaults()])
#' @param notes character vector of run notes
#' @return object of class `scphasedcn`
#' @export
new_scphasedcn = function(mode, cells, regions, features = NULL, scores = NULL, thresholds = NULL, calls = NULL,
                       events = NULL, pools = NULL, clones = NULL, assignments = NULL, segments = NULL,
                       profiles = NULL, null_fit = NULL, params = caller_defaults(), notes = character()) {
  mode = match.arg(mode, c("normal", "tumour"))
  if (is.character(cells)) cells = data.table(cell = cells)
  if (!is.data.frame(cells) || !"cell" %in% names(cells)) stop("new_scphasedcn: `cells` must be a character vector or a table with a `cell` column")
  cells = as.data.table(cells)
  if (anyDuplicated(cells$cell)) stop("new_scphasedcn: duplicated cells")
  check_regions(regions)
  if (!is.null(features)) check_features(features)
  if (!is.null(thresholds)) check_thresholds(thresholds)
  if (!is.null(calls)) check_calls(calls)
  if (!is.null(segments)) check_segments(segments)
  structure(list(mode = mode, cells = cells, regions = as.data.table(regions), features = features, scores = scores,
                 thresholds = thresholds, calls = calls, events = events, pools = pools, clones = clones,
                 assignments = assignments, segments = segments, profiles = profiles, null_fit = null_fit,
                 params = params, notes = as.character(notes)),
            class = "scphasedcn")
}

#' @export
print.scphasedcn = function(x, ...) {
  n_clones = if (is.null(x$clones)) 0L else if (is.data.frame(x$clones)) {
    if ("clone" %in% names(x$clones)) uniqueN(x$clones$clone) else nrow(x$clones)
  } else length(x$clones)
  cat(sprintf("<scphasedcn> mode: %s\n", x$mode))
  cat(sprintf("  cells: %d\n", nrow(x$cells)))
  cat(sprintf("  regions: %d\n", nrow(x$regions)))
  if (is.null(x$calls)) cat("  calls: none\n") else
    cat(sprintf("  calls: %d rows, %d called, %d reportable\n", nrow(x$calls),
                sum(x$calls$called, na.rm = TRUE), sum(x$calls$called & x$calls$reportable, na.rm = TRUE)))
  cat(sprintf("  clones: %d\n", n_clones))
  if (length(x$notes)) cat(sprintf("  notes: %d\n", length(x$notes)))
  invisible(x)
}

# objects of the working-name class "scncall" (0.1.0-pre builds, e.g. a result.rds loaded with readRDS()) print the same
#' @export
print.scncall = function(x, ...) {
  print.scphasedcn(x, ...)
  invisible(x)
}

# ---- naming --------------------------------------------------------------------------------------------

#' Tile region name
#'
#' `<label>|<window id>` with the window id of [window_ids()]
#' (`chrom:start`, start as a plain integer), e.g. `"5Mb|chr1:0"`.
#' @param label tile size label ("5Mb", "10Mb", "20Mb", "50Mb")
#' @param chrom chromosome
#' @param start tile start (bp, 0-based)
#' @return character vector
#' @export
tile_name = function(label, chrom, start) {
  paste0(label, "|", chrom, ":", formatC(start, format = "d", big.mark = ""))
}
