# call_scn.R
# What: the pipeline of the copy-number caller: the top-level call_scn() (features ->
#       optional chromatin embedding -> normal-tissue caller, or tumour-mode clone caller), the caller
#       parameters from inst/config/caller.yaml (caller_params), the output writer / reader
#       (write_scphasedcn, read_scphasedcn), a headline summary (scphasedcn_summary) and the null false-call
#       accounting on the null cells (null_false_calls).
# Why: the feature step (features.R), the scorer (caller.R) and the tumour mode (clones_*.R) are separate building
#      blocks; a user runs one
#      function (or inst/scripts/call_scn.R) on an extraction directory and gets one result object and one
#      output directory with flat tables, a parameter file and run notes.
# Inputs: an extraction directory (read_sample_depth() layout) or a precomputed `scn_features` (cell_features()),
#         optional cell metadata, null cells, multiplets, events.
# Outputs: `scphasedcn` objects; <out>/{cells,calls,events,pooled_calls,thresholds,null_fit,false_calls,...}.tsv,
#          result.rds, params.yaml, NOTES.md.

utils::globalVariables(c("called", "reportable", "level", "cell", "n_calls", "n_events", "prob", "qc_pass",
                         "state_comb", "call_s_comb", "is_null", "class", "stratum", "chrom", "n_null", "n_called",
                         "per_1000", "kind", "value", "ref_pool", "null_cell", "section", "item", "n"))

# ---- parameters -------------------------------------------------------------------------------------

#' Caller parameters from the caller config
#'
#' [caller_defaults()] updated from a YAML file (default `inst/config/caller.yaml` of the package: one key per
#' default, named vectors as maps) and then from `overrides`. Values are coerced to the type of the default
#' (numeric, logical, character; names kept for named vectors). Empty strings and nulls keep the default
#' (e.g. `tabix: ""` = the `tabix` on the PATH). Unknown keys in the file are ignored with a
#' warning; unknown keys in `overrides` are added as-is.
#'
#' Two optional blocks hold the arguments of the two callers: `tumour:` (tumour-mode options, any
#' [clone_defaults()] key) and `normal:` (arguments of [caller_normal()] other than the inputs [call_scn()] sets,
#' e.g. `n_draws`, `pool_loo`). When a block is present the result carries it as `$tumour` / `$normal`, filled
#' up to the full defaults ([clone_defaults()], the defaults of [caller_normal()]); [call_scn()] passes it to the
#' caller of its mode (its `caller_args` take precedence) and records the values used in the result's `params`,
#' so the `params.yaml` of [write_scphasedcn()] reads back with this function.
#' @param file YAML file (NULL = the package's `caller.yaml`; a missing file gives the defaults)
#' @param overrides named list of values that replace the file's / default values
#' @return named list (the keys of [caller_defaults()], plus `tumour` / `normal` when the file has those blocks)
#' @export
caller_params = function(file = NULL, overrides = list()) {
  p = caller_defaults()
  if (is.null(file)) file = system.file("config", "caller.yaml", package = "scphasedcn")
  if (length(file) == 1 && nzchar(file) && file.exists(file)) {
    y = read_config(file)
    if (!is.null(y$caller) && is.list(y$caller)) y = y$caller
    blocks = c("tumour", "normal")
    unknown = setdiff(names(y), c(names(p), blocks))
    if (length(unknown)) warning("caller_params: unknown key(s) in ", file, ": ", paste(unknown, collapse = ", "), call. = FALSE)
    for (k in intersect(names(y), names(p))) {
      v = coerce_param(y[[k]], p[[k]])
      if (!is.null(v)) p[[k]] = v
    }
    for (b in intersect(blocks, names(y)))
      p[[b]] = block_params(y[[b]], if (b == "tumour") clone_defaults() else normal_arg_defaults(), b, file)
  }
  for (k in names(overrides)) p[[k]] = if (k %in% names(caller_defaults())) coerce_param(overrides[[k]], caller_defaults()[[k]]) %||% p[[k]] else overrides[[k]]
  p
}

# defaults of the caller_normal() arguments that call_scn() does not set itself (the `normal:` block)
normal_arg_defaults = function() {
  f = formals(caller_normal)
  keep = setdiff(names(f), c("feat", "params", "cell_meta", "multiplets", "regions", "embedding", "null_cells"))
  out = lapply(keep, function(k) eval(f[[k]]))
  names(out) = keep
  out$depth_model = out$depth_model[1]
  out
}

# a `tumour:` / `normal:` block: the defaults updated from the yaml values (coerced to the defaults' types; a logical
# default also takes a string, e.g. split_pass "disagree"); unknown keys are dropped with a warning
block_params = function(y, def, block, file = "") {
  if (is.null(y) || !is.list(y)) return(def)
  unknown = setdiff(names(y), names(def))
  if (length(unknown)) warning("caller_params: unknown key(s) in the ", block, " block of ", file, ": ", paste(unknown, collapse = ", "), call. = FALSE)
  for (k in intersect(names(y), names(def))) {
    d = def[[k]]
    v = if (is.null(d) || is.list(d)) y[[k]] else coerce_param(y[[k]], d)
    if (is.logical(d) && length(v) && all(is.na(v))) {
      s = unlist(y[[k]], use.names = FALSE)
      if (is.character(s) && !all(s %in% c("", "null", "NULL", "~"))) v = s
    }
    if (!is.null(v)) def[k] = list(v)
  }
  def
}

# split call_scn's params into the caller parameters and the argument blocks of the two callers
split_params = function(params) {
  list(params = params[setdiff(names(params), c("tumour", "normal"))], tumour = params$tumour, normal = params$normal)
}

# coerce a config value to the type (and names) of the default; NULL / "" / "null" / "~" keep the default
coerce_param = function(v, def) {
  if (is.null(v) || length(v) == 0) return(NULL)
  nm = names(v)
  if (is.list(v)) { nm = names(v); v = unlist(v, use.names = FALSE) }
  if (is.character(v) && length(v) == 1 && v %in% c("", "null", "NULL", "~")) return(NULL)
  out = if (is.logical(def)) as.logical(v) else if (is.numeric(def)) as.numeric(v) else as.character(v)
  if (!is.null(nm)) names(out) = gsub("^[\"']|[\"']$", "", nm) else if (!is.null(names(def)) && length(out) == length(def)) names(out) = names(def)
  out
}

# minimal YAML writer for a parameter list (scalars, vectors as inline lists, named vectors as block maps)
yaml_lines = function(x, indent = 0) {
  pad = strrep(" ", indent)
  q = function(v) {
    if (is.logical(v)) return(ifelse(is.na(v), "null", ifelse(v, "true", "false")))
    if (is.numeric(v)) return(ifelse(is.na(v), "null", format(v, scientific = FALSE, trim = TRUE, digits = 15)))
    sprintf("\"%s\"", gsub("\"", "\\\\\"", as.character(v)))
  }
  out = character()
  for (k in names(x)) {
    v = x[[k]]
    if (is.null(v)) { out = c(out, sprintf("%s%s: null", pad, k)); next }
    if (is.list(v)) { out = c(out, sprintf("%s%s:", pad, k), yaml_lines(v, indent + 2)); next }
    if (!is.null(names(v)) && length(v) > 0) {
      out = c(out, sprintf("%s%s:", pad, k), sprintf("%s  %s: %s", pad, names(v), q(v)))
    } else if (length(v) == 1) {
      out = c(out, sprintf("%s%s: %s", pad, k, q(v)))
    } else out = c(out, sprintf("%s%s: [%s]", pad, k, paste(q(v), collapse = ", ")))
  }
  out
}

# ---- inputs -------------------------------------------------------------------------------------------

# a table from a path or a data.frame (NULL stays NULL)
as_input_table = function(x, what) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1 && file.exists(x)) x = fread(x)
  if (!is.data.frame(x)) stop(sprintf("call_scn: `%s` must be a table or a path to one", what))
  x = as.data.table(x)
  if (!"cell" %in% names(x)) setnames(x, 1, "cell")
  x[, cell := as.character(cell)]
  x[]
}

# a cell vector from a vector, a file (one cell per line or a `cell` column) or a table
as_cell_vector = function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(as.character(if ("cell" %in% names(x)) x$cell else x[[1]]))
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    d = read_cell_file(x)
    return(as.character(if ("cell" %in% names(d)) d$cell else d[[1]]))
  }
  as.character(x)
}

# a one-column file of cells has no header unless its first line is "cell"
read_cell_file = function(path) {
  d = fread(path, header = "auto")
  if (ncol(d) == 1 && names(d)[1] != "cell") d = fread(path, header = FALSE)
  d
}

#' Multiplet table of the caller
#'
#' Normalises the forms a multiplet flag comes in to a table `cell, prob`: a character vector (or a file with
#' one cell per line) of flagged cells (prob 1); a table with `prob` (kept); with a logical `multiplet` column
#' (prob 1 / 0); with a logical `clean_singlet` column (e.g. a multiplet scorer's table: prob 0 for clean
#' singlets, 1 otherwise); any other table (prob 1 for every listed cell).
#' @param x NULL, a vector, a path or a table
#' @return data.table cell, prob (NULL for NULL input)
#' @export
multiplet_table = function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && !(length(x) == 1 && file.exists(x))) return(data.table(cell = x, prob = 1))
  if (is.character(x)) {
    d = read_cell_file(x)
    if (ncol(d) == 1) return(data.table(cell = as.character(d[[1]]), prob = 1))
    x = d
  }
  d = as.data.table(x)
  if (!"cell" %in% names(d)) setnames(d, 1, "cell")
  p = if ("prob" %in% names(d)) as.numeric(d$prob) else if ("multiplet" %in% names(d)) as.numeric(as.logical(d$multiplet)) else
    if ("clean_singlet" %in% names(d)) 1 - as.numeric(as.logical(d$clean_singlet)) else rep(1, nrow(d))
  p[is.na(p)] = 0
  data.table(cell = as.character(d$cell), prob = p)
}

#' Chromatin embedding of the cells of a sample
#'
#' [chromatin_embedding()] (high-pass filtered log ratios to the all-cell leave-one-out reference, top principal
#' components) of the given cells on the used bins of an extraction directory, with the cell ids as row
#' names, ready for [knn_pools()]. Built from depth only (no copy-number features), as kNN pools require.
#' @param dir extraction directory (bins_100kb.h5ad)
#' @param cells cells to embed (default all)
#' @param use logical bin filter over all bins of the sample (default [bin_filter()] on chr1-22 + chrX)
#' @param n_pc components
#' @param hp_window high-pass window (bins)
#' @param seed RNG seed of the truncated SVD (irlba's random start)
#' @return cells x components matrix
#' @export
scn_embedding = function(dir, cells = NULL, use = NULL, n_pc = 20, hp_window = 21, seed = 1) {
  scn = read_sample_depth(dir)
  if (!is.null(cells)) {
    miss = setdiff(cells, scn$cells$cell)
    if (length(miss)) stop(sprintf("scn_embedding: %d cells not in %s", length(miss), dir))
    scn = subset_scn(scn, cells = cells)
  }
  if (is.null(use)) use = bin_filter(scn$bins, scn$counts) & scn$bins$chrom %in% c(AUTOSOMES_CALLER, "chrX")
  if (length(use) != nrow(scn$bins)) stop("scn_embedding: `use` must have one value per bin")
  X = dense(scn$counts[use, , drop = FALSE])
  s = colSums(X)
  set.seed(seed)                                   # irlba starts from a random vector: seed it for reproducible kNN pools
  emb = chromatin_embedding(X, s, scn$bins[use], hp_window = hp_window, n_pc = n_pc)
  rownames(emb) = scn$cells$cell
  emb
}

#' 1 Mb window depth of a precomputed features object
#'
#' Rebuilds the `win` component of an `scn_features` object saved without it (e.g. to save disk space):
#' the same cells and bin filter (`features$use`), [normalise_depth()] with `method` (`model` for tumours, `celltype` with the
#' features' cell types), then [window_matrices()]. Needed by the tumour mode (clone pseudobulks).
#' @param dir extraction directory
#' @param features `scn_features` (cells, use)
#' @param method depth normalisation ("model" tumour, "celltype" normal)
#' @param size window size (bp)
#' @param min_cells_type celltype: minimal reference cells of a type
#' @return list `x`, `e`, `windows` ([window_matrices()])
#' @export
feature_windows = function(dir, features, method = "model", size = 1e6, min_cells_type = 20) {
  scn = read_sample_depth(dir)
  cells = features$cells$cell
  miss = setdiff(cells, scn$cells$cell)
  if (length(miss)) stop(sprintf("feature_windows: %d cells not in %s", length(miss), dir))
  scn = subset_scn(scn, cells = cells)
  use = features$use %||% (bin_filter(scn$bins, scn$counts) & scn$bins$chrom %in% c(AUTOSOMES_CALLER, "chrX"))
  norm = normalise_depth(scn, method = method, use = use, cell_type = features$cells$cell_type, min_cells_type = min_cells_type)
  window_matrices(scn, norm, size)
}

# code revision of a run: SCPHASEDCN_REVISION (env), the commit an install_github / pak install records (RemoteSha), or
# `git rev-parse` only when the working directory is a scphasedcn source checkout (a run from the repository
# root), so a run never records the revision of an unrelated repository
git_rev = function() {
  env = Sys.getenv("SCPHASEDCN_REVISION")
  if (nzchar(env)) return(env)
  sha = tryCatch(utils::packageDescription("scphasedcn")$RemoteSha, error = function(e) NULL)
  if (is.character(sha) && length(sha) == 1 && nzchar(sha)) return(substr(sha, 1, 7))
  src = file.exists("DESCRIPTION") &&
    isTRUE(tryCatch(read.dcf("DESCRIPTION", fields = "Package")[1, 1] == "scphasedcn", error = function(e) FALSE))
  if (src) {
    r = tryCatch(suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
                 error = function(e) character())
    if (length(r) == 1 && nzchar(r) && is.null(attr(r, "status"))) return(r)
  }
  if (file.exists("REVISION")) return(paste0(readLines("REVISION", n = 1), " (snapshot)"))
  NA_character_
}

# ---- null false calls -----------------------------------------------------------------------------

#' Realised false calls on the null rows of a scored table
#'
#' Rates per class x depth stratum (and per class over all strata) of `call_s_comb` on the null rows
#' (`is_null`, else `null_cell & null_region & eligible`), and the per-null-cell false arm events (called
#' autosomal arm null rows per null cell), overall and per called state.
#' @param scores scored table with `call_s_comb` ([caller_normal()]'s `scores`)
#' @return list(rates = data.table class, stratum, n_null, n_called, per_1000; events = data.table state,
#'   n_cells, arm_calls_per_cell, frac_cells)
#' @export
null_false_calls = function(scores) {
  d = as.data.table(scores)
  if (!"call_s_comb" %in% names(d)) stop("null_false_calls: no call_s_comb column (apply_thresholds first)")
  isn = null_flag(d)
  x = d[isn, .(cell, class, stratum, chrom = if ("chrom" %in% names(d)) chrom else region_chrom(region, class),
               st = state_comb, k = !is.na(call_s_comb) & call_s_comb)]
  r1 = x[, .(n_null = .N, n_called = sum(k)), by = .(class, stratum)]
  r2 = x[, .(stratum = "all", n_null = .N, n_called = sum(k)), by = class]
  rates = rbind(r1, r2)
  rates[, per_1000 := 1000 * n_called / n_null]
  rates = rates[order(match(class, region_classes()), stratum != "all", stratum)]
  arm = x[class == "arm" & chrom %in% AUTOSOMES_CALLER]
  ev = rbindlist(lapply(c("any", "loss", "gain", "cnloh"), function(s) {
    a = arm[, .(k = sum(k & (s == "any" | st == s))), by = cell]
    data.table(state = s, n_cells = nrow(a), arm_calls_per_cell = if (nrow(a)) mean(a$k) else NA_real_,
               frac_cells = if (nrow(a)) mean(a$k > 0) else NA_real_)
  }))
  list(rates = rates[], events = ev)
}
utils::globalVariables(c("st", "k", "region"))

# ---- top level ----------------------------------------------------------------------------------------

#' Single-cell copy-number calls of a sample
#'
#' Normal tissue (`mode = "normal"`): [cell_features()] with `celltype` (cell-type matched) depth and the oriented allele
#' channel -> optional chromatin embedding ([scn_embedding()]) for kNN pools -> [caller_normal()] (null cells,
#' Student-t / beta-binomial nulls, fold-cross-validated budget thresholds, per-cell calls, pooled calls,
#' events). Tumour (`mode = "tumour"`): [cell_features()] with `model` depth, automatic phasing and the pseudobulk
#' groups `groups`, event chromosomes excluded from the null -> [caller_tumour()] (clones, clone
#' segments, cell assignment, per-cell inherited / confirmed calls).
#' @param dir extraction directory (bins_100kb.h5ad, fragments.tsv.gz + .tbi, cells_qc.tsv)
#' @param mode "normal" or "tumour"
#' @param cell_meta NULL, a table or a TSV path with `cell` and optional `cell_type, flat, fold, group,
#'   stratum, n_frag` (passed to [cell_features()]; `flat` becomes the features' `null_cell`)
#' @param cells NULL (QC-pass cells) or the cells to call (vector or file)
#' @param sex "female" / "male" (required for real data; NULL uses a `sex` column of `cells_qc.tsv` if there is one, as in
#'   simulated data)
#' @param events named events: a table `region, chrom, start, end`, [cli_events()] output or a named vector of
#'   region strings (`20q_gain = "chr20:28000000-"`); scored as `event` regions; in tumour mode their chromosomes
#'   are excluded from the null rows
#' @param null_cells NULL / "auto" (the features' `null_cell` when informative, else selected), "features",
#'   "select" ([select_null_cells()]), or the null cells themselves (vector or file)
#' @param multiplets NULL or multiplet flags ([multiplet_table()] forms); normal mode: cells with prob > 0.5
#'   are never null cells and never pooled; tumour mode: passed to `caller_tumour()`
#' @param control_regions tumour mode: control regions for the null (regions table)
#' @param groups tumour mode: pseudobulk groups of the block orientation (a column of `cell_meta`, a vector named
#'   by cell, or one label per cell); ignored with a warning in normal mode
#' @param params [caller_params()] / [caller_defaults()]; its optional `normal` / `tumour` blocks
#'   ([caller_params()]) are the default arguments of [caller_normal()] / `caller_tumour()` (`opts`)
#' @param features NULL (computed) or a precomputed `scn_features` (object or .rds path) to skip the feature step
#' @param embedding TRUE (compute [scn_embedding()] from `dir`; skipped with a warning when that fails), FALSE,
#'   or a cells x components matrix with cell row names
#' @param verbose log progress
#' @param keep_scores keep the per cell x region scores in the result (normal mode; large)
#' @param caller_args extra arguments for [caller_normal()] / `caller_tumour()` (e.g. `n_draws`, `opts`); they take
#'   precedence over `params$normal` / `params$tumour` (tumour: `opts` merged key by key)
#' @param phasing phasing / block orientation of the allele features ([cell_features()]: "auto", "chromosome",
#'   "block", "loo" or "mirrored"; both modes; ignored with precomputed `features`). Normal mode with `loo` on
#'   block phasing (read-backed blocks without a population phase) gives a warning and a run note: per-cell allelic
#'   calls of events carried by few cells need chromosome-scale phasing (`vignette("scphasedcn-inputs")`,
#'   `scphasedcn hybrid-phase`)
#' @param ... extra arguments for [cell_features()] (e.g. `centromeres`, `classes`, `qc`)
#' @return `scphasedcn` with `features` (the long table), `feature_info` (path, cells, phasing, phasing_chrom (the
#'   per-chromosome phasing regimes, [phasing_by_chrom()]), notes), `scores`
#'   (normal, `keep_scores`), `thresholds`, `calls`, `events`, `pools`, `null_fit`, `false_calls` (normal:
#'   [null_false_calls()]), `embedding`, clone components (tumour), `params` (the caller parameters plus the
#'   arguments used: `normal` = every [caller_normal()] argument of the run, `tumour` = the full tumour options) and
#'   `notes` (run time, host, inputs, git revision)
#' @export
call_scn = function(dir, mode = c("normal", "tumour"), cell_meta = NULL, cells = NULL, sex = NULL, events = NULL,
                    null_cells = NULL, multiplets = NULL, control_regions = NULL, groups = NULL,
                    params = caller_defaults(), features = NULL, embedding = TRUE, verbose = TRUE,
                    keep_scores = TRUE, caller_args = list(), phasing = "auto", ...) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  mode = match.arg(mode)
  phasing = match.arg(phasing, c("auto", "chromosome", "block", "loo", "mirrored"))
  # caller argument blocks of the parameters (caller_params `normal:` / `tumour:`); caller_args take precedence
  sp = split_params(params)
  params = sp$params
  if (mode == "normal") {
    caller_args = utils::modifyList(sp$normal %||% list(), as.list(caller_args))
  } else {
    topts = utils::modifyList(as.list(sp$tumour %||% list()), as.list(caller_args$opts %||% list()))
    if (length(topts)) caller_args$opts = topts
  }
  if (mode == "normal" && !is.null(groups)) {
    warning("call_scn: `groups` (pseudobulk orientation groups) is used in tumour mode only; ignored", call. = FALSE)
    groups = NULL
  }
  t0 = Sys.time()
  mins = function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say = function(...) if (verbose) cli_log("[call_scn] ", ...)
  cm = as_input_table(cell_meta, "cell_meta")
  cells = as_cell_vector(cells)
  mult = multiplet_table(multiplets)
  nc_mode = "auto"; nc_cells = NULL
  if (!is.null(null_cells)) {
    if (is.character(null_cells) && length(null_cells) == 1 && null_cells %in% c("auto", "features", "select")) nc_mode = null_cells
    else { nc_cells = as_cell_vector(null_cells); nc_mode = "features" }
  }
  if (mode == "tumour" && !is.null(nc_cells) && !is.null(cm)) cm[, flat := cell %in% nc_cells]
  # pseudobulk groups as a vector named by cell (cell_features only knows its own meta columns)
  grp = groups
  if (!is.null(grp) && is.character(grp) && length(grp) == 1 && is.null(names(grp))) {
    if (is.null(cm) || !grp %in% names(cm)) stop("call_scn: `groups` names a column that `cell_meta` does not have: ", grp)
    grp = setNames(as.character(cm[[grp]]), cm$cell)
  }
  inputs = c(sprintf("dir: %s", dir), sprintf("mode: %s", mode),
             sprintf("cell_meta: %s", if (is.character(cell_meta)) cell_meta else if (is.null(cm)) "none" else sprintf("table (%d cells)", nrow(cm))),
             sprintf("cells: %s", if (is.null(cells)) "QC pass" else length(cells)),
             sprintf("sex: %s", sex %||% "from cells"), sprintf("null_cells: %s", if (is.null(nc_cells)) nc_mode else sprintf("given (%d)", length(nc_cells))),
             sprintf("multiplets: %s", if (is.null(mult)) "none" else sprintf("%d flagged of %d", sum(mult$prob > 0.5), nrow(mult))))

  # 1. features
  feat_path = NULL
  if (is.character(features)) { feat_path = features; say("reading features ", features); features = readRDS(features) }
  if (is.null(features)) {
    say("cell_features (", if (mode == "normal") "celltype" else "model", " depth)")
    if (mode == "normal") {
      features = cell_features(dir, cells = cells, cell_meta = cm, method = "celltype", sex = sex, events = events, params = params,
                               phasing = phasing, verbose = verbose, ...)
    } else {
      # chromosomes of the events (null_exclude): from the region strings, so open ends ("chr20:28000000-") and bare
      # chromosomes need no chromosome sizes here (cell_features builds the event table with them)
      nex = if (is.null(events)) character() else
        if (is.data.frame(events)) unique(as.character(events$chrom)) else unique(sub(":.*$", "", unlist(events)))
      features = cell_features(dir, cells = cells, cell_meta = cm, method = "model", sex = sex, events = events, null_exclude = nex,
                               phasing = phasing, groups = grp, params = params, verbose = verbose, ...)
    }
    t_feat = mins()
  } else {
    if (!is.list(features) || is.null(features$long)) stop("call_scn: `features` must be an scn_features object (cell_features()) or its .rds")
    if (phasing != "auto") warning("call_scn: `phasing` is ignored with precomputed `features`", call. = FALSE)
    t_feat = 0
    if (!is.null(cells) && mode == "normal") {
      features$long = features$long[cell %in% cells]
      features$cells = features$cells[cell %in% cells]
    }
  }
  long = features$long
  if (!is.null(nc_cells)) { long = copy(long); long[, null_cell := cell %in% nc_cells] }
  say(sprintf("features: %d cells x %d regions (%.1f min)", nrow(features$cells), nrow(features$regions), mins()))
  ftab = features$cells
  phase_note = phasing_run_note(features$phasing, mode)
  if (!is.null(phase_note)) warning("call_scn: ", phase_note, call. = FALSE)
  # chromosomes whose phasing regime differs from the run's (a note of cell_features() that the run notes repeat)
  mixed_note = if (phasing_ok(features$phasing))
    phasing_mixed_note(features$phasing_chrom, features$phasing[["detected"]], features$phasing[["orientation"]])

  if (mode == "normal") {
    if (!is.null(control_regions)) warning("call_scn: `control_regions` is used in tumour mode only", call. = FALSE)
    # 2. chromatin embedding for kNN pools
    emb = NULL
    if (isTRUE(embedding)) {
      t_e = Sys.time()
      emb = tryCatch(scn_embedding(dir, ftab$cell, features$use, seed = params$seed), error = function(e) {
        warning("call_scn: no chromatin embedding (", conditionMessage(e), "); kNN pools skipped", call. = FALSE)
        NULL
      })
      if (!is.null(emb)) say(sprintf("embedding: %d cells x %d components (%.1f min)", nrow(emb), ncol(emb), as.numeric(difftime(Sys.time(), t_e, units = "mins"))))
    } else if (!isFALSE(embedding) && !is.null(embedding)) {
      emb = as.matrix(embedding)
      if (is.null(rownames(emb))) stop("call_scn: `embedding` needs cell ids as row names")
    }
    # 3. normal-tissue caller
    mcells = if (is.null(mult)) NULL else mult[prob > 0.5, cell]
    ctm = if (!is.null(cm) && "cell_type" %in% names(cm)) cm[, .(cell, cell_type = as.character(cell_type))] else NULL
    # 2b. optional allelic tile mask: null cells of a first selection pass -> overdispersed 5 Mb tiles -> allele
    #     features re-aggregated without their votes (depth untouched); caller_normal then runs as usual
    amask = NULL; mask_note = NULL
    if (isTRUE(params$mask_allelic)) {
      t_m = Sys.time()
      if (is.null(features$votes) || is.null(features$atoms)) {
        warning("call_scn: mask_allelic needs features$votes and features$atoms; no mask applied", call. = FALSE)
      } else {
        sel0 = select_null_cells(long[eligible %in% TRUE], params = params, exclude = mcells)
        isn0 = long$cell %in% sel0$null_cells & long$null_region %in% TRUE & long$eligible %in% TRUE & long$chrom %in% AUTOSOMES_CALLER
        amask = allelic_noisy_tiles(long, null = isn0, params = params)
        nm = sum(amask$masked)
        if (nm) {
          ex = atoms_in_tiles(features$atoms, amask)
          features = update_allelic_features(features, ex, params = params)
          long = features$long
          if (!is.null(nc_cells)) { long = copy(long); long[, null_cell := cell %in% nc_cells] }
        }
        mask_note = sprintf("allelic mask: %d of %d assessed 5 Mb tiles masked (%s MADs of log rho; %s) from %d first-pass null cells; %.1f min",
                            nm, sum(amask$assessed), format(params$mask_rho_mad), if (nm) paste(amask[masked == TRUE, tile], collapse = ", ") else "none",
                            length(sel0$null_cells), as.numeric(difftime(Sys.time(), t_m, units = "mins")))
        say(mask_note)
      }
    }
    say("caller_normal (null cells: ", nc_mode, ")")
    res = do.call(caller_normal, c(list(feat = long, params = params, cell_meta = ctm, multiplets = mcells, regions = features$regions,
                                        embedding = emb, null_cells = nc_mode), caller_args))
    if (!is.null(amask)) {
      res$null_fit$allelic_mask = amask
      res$notes = c(res$notes, mask_note)
    }
    fc = null_false_calls(res$scores)
    res$false_calls = fc$rates
    res$false_events = fc$events
    if (!keep_scores) res["scores"] = list(NULL)
    res$embedding = emb
  } else {
    if (is.null(features$votes)) stop("call_scn: tumour mode needs the molecule votes (features$votes); compute the features without `features` or save them with their votes")
    if (is.null(features$win)) {
      say("rebuilding 1 Mb window depth (model) for the clone pseudobulks")
      features$win = feature_windows(dir, features, "model")
    }
    if (is.null(features$null_fit)) {
      # depth t null and allelic rho on the null rows (null cells if any, else all cells) for the tumour mode
      nr = !is.na(long$null_region) & long$null_region & !is.na(long$eligible) & long$eligible
      if (any(long$null_cell %in% TRUE)) nr = nr & long$null_cell %in% TRUE
      if (any(nr)) {
        tpar = complete_tpar(long, fit_depth_null(long[nr]))
        arm = nr & long$class == "arm" & long$n_hap > 0
        rho = if (any(arm, na.rm = TRUE)) fit_allelic_null(long$nA[arm %in% TRUE], long$n_hap[arm %in% TRUE]) else NA_real_
        features$null_fit = list(tpar = tpar, rho = rho)
      }
    }
    args = list(feat = features, params = params, multiplets = mult, control_regions = control_regions)
    say("caller_tumour")
    res = do.call(caller_tumour, c(args, caller_args))
    res$embedding = NULL
  }
  res$notes = c(res$notes, mixed_note, phase_note)

  # cells table: QC, stratum, fold, null cell, number of reportable cell calls, clone columns
  ct = as.data.table(res$cells)
  ct[, cell := as.character(cell)]
  if (!is.null(ftab)) {
    for (cn in c("qc_pass", "n_frag", "stratum", "fold", "cell_type")) if (cn %in% names(ftab) && !cn %in% names(ct))
      set(ct, j = cn, value = ftab[[cn]][match(ct$cell, ftab$cell)])
    if (!"null_cell" %in% names(ct)) ct[, null_cell := cell %in% unique(long[null_cell %in% TRUE, cell])]
  }
  setnames(ct, "qc_pass", "qc", skip_absent = TRUE)
  if (!is.null(res$calls)) {
    cc = as.data.table(res$calls)[level == "cell" & called %in% TRUE & reportable %in% TRUE, .(n_calls = .N), by = cell]
    ct[, n_calls := cc$n_calls[match(cell, cc$cell)]]
  } else ct[, n_calls := NA_integer_]
  ct[is.na(n_calls), n_calls := 0L]
  if (!is.null(res$events) && nrow(res$events)) {
    ee = as.data.table(res$events)[level == "cell", .(n_events = .N), by = cell]
    ct[, n_events := ee$n_events[match(cell, ee$cell)]]
    ct[is.na(n_events), n_events := 0L]
  } else ct[, n_events := 0L]
  if (!is.null(res$assignments)) {
    asg = as.data.table(res$assignments)
    for (cn in intersect(c("clone", "posterior", "doublet_like", "ref_pool", "assigned"), names(asg))) if (!cn %in% names(ct))
      set(ct, j = cn, value = asg[[cn]][match(ct$cell, as.character(asg$cell))])
  }
  if (!is.null(mult) && !"multiplet" %in% names(ct)) ct[, multiplet := cell %in% mult[prob > 0.5, cell]]
  lead = intersect(c("cell", "qc", "stratum", "fold", "null_cell", "n_calls", "n_events", "clone", "posterior", "doublet_like", "ref_pool"), names(ct))
  setcolorder(ct, lead)
  res$cells = ct

  res$features = long
  res$feature_info = list(path = feat_path, cells = ftab, phasing = features$phasing, phasing_chrom = features$phasing_chrom,
                          notes = features$notes, use = features$use)
  # parameters of the run: the caller's own record (tumour: + the full `tumour` options) over the input parameters,
  # plus every caller_normal() argument (defaults filled in) in normal mode
  rp = res$params %||% list()
  pr = params
  for (k in names(rp)) pr[k] = list(rp[[k]])
  if (mode == "normal") {
    na = normal_arg_defaults()
    for (k in names(caller_args)) na[k] = list(caller_args[[k]])
    pr$normal = na
  }
  res$params = pr
  host = Sys.info()[["nodename"]]
  res$notes = c(res$notes,
                sprintf("call_scn %s: %.1f min (features %.1f min) on %s, %s", mode, mins(), t_feat, host, format(Sys.time(), "%F %T")),
                sprintf("scphasedcn %s from %s; git %s", as.character(utils::packageVersion("scphasedcn")),
                        tryCatch(find.package("scphasedcn"), error = function(e) "?"), git_rev()),
                if (!is.null(feat_path)) sprintf("features: %s", feat_path) else "features: computed",
                inputs)
  res$run = list(mode = mode, minutes = mins(), host = host, time = Sys.time(), git = git_rev(), inputs = inputs)
  say(sprintf("done (%.1f min)", mins()))
  res
}
utils::globalVariables(c("flat", "cell_type", "multiplet"))

# the run note (and warning text) of normal mode on block phasing: allele features oriented by the other cells'
# pseudobulk (`loo`) on unoriented read-backed blocks; NULL otherwise (tumour mode, chromosome-wide phasing, other
# orientations, or features without a phasing record)
phasing_run_note = function(phasing, mode) {
  if (mode != "normal" || !phasing_ok(phasing)) return(NULL)
  if (!(phasing[["detected"]] == "block" && phasing[["orientation"]] %in% c("loo", "block"))) return(NULL)
  paste("normal mode on block phasing: the allele features were oriented by the pseudobulk of the other cells (loo),",
        "which is blind to events carried by few cells and aligns allelic bias shared by all cells. Per-cell allelic",
        "calls for events carried by few cells are not supported with block-only phasing; use chromosome-scale",
        "(hybrid) phasing: see vignette(\"scphasedcn-inputs\") and `scphasedcn hybrid-phase`.")
}

# a phasing record with `detected` and `orientation`
phasing_ok = function(phasing) !is.null(phasing) && all(c("detected", "orientation") %in% names(phasing))

# ---- output ------------------------------------------------------------------------------------------------

# null fit as one table: t rows (class, stratum, t_m, t_s, t_df, n) + a rho row
null_fit_table = function(nf) {
  if (is.null(nf)) return(data.table(kind = character()))
  t = if (!is.null(nf$tpar)) copy(as.data.table(nf$tpar))[, kind := "depth_t"] else NULL
  r = if (!is.null(nf$rho) && length(nf$rho) == 1) data.table(kind = "allelic_rho", rho = as.numeric(nf$rho)) else
    if (is.data.frame(nf$rho)) copy(as.data.table(nf$rho))[, kind := "allelic_rho"] else NULL
  extra = if (!is.null(nf$ref_thr) && length(nf$ref_thr) == 1 && is.numeric(nf$ref_thr)) data.table(kind = "reference_thr", thr = nf$ref_thr) else NULL
  # tumour mode: the confirmation null (rho / depth t refitted on the reference pool, or the input fits)
  tc = if (!is.null(nf$tpar_confirm)) copy(as.data.table(nf$tpar_confirm))[, kind := "depth_t_confirm"] else NULL
  rc = if (!is.null(nf$rho_confirm) && length(nf$rho_confirm) == 1) data.table(kind = "allelic_rho_confirm", rho = as.numeric(nf$rho_confirm)) else NULL
  out = rbindlist(list(t, r, extra, tc, rc), fill = TRUE)
  if (nrow(out)) setcolorder(out, "kind")
  out[]
}

flatten_profiles = function(p) {
  if (is.null(p)) return(NULL)
  if (is.data.frame(p)) return(as.data.table(p))
  if (is.list(p)) {
    parts = Filter(is.data.frame, p)
    if (length(parts)) return(rbindlist(lapply(names(parts), function(k) copy(as.data.table(parts[[k]]))[, part := k]), fill = TRUE))
  }
  NULL
}
utils::globalVariables("part")

#' Write a caller result to a directory
#'
#' Flat tables: `cells.tsv` (cell, qc, stratum, fold, null_cell, n_calls, n_events, clone / posterior /
#' doublet_like / ref_pool when present), `calls.tsv` (cell- and clone-level calls: called rows in normal mode,
#' all inherited / confirmed rows in tumour mode), `pooled_calls.tsv` (called pool rows), `events.tsv`,
#' `thresholds.tsv`, `null_fit.tsv` (depth t rows + rho), `false_calls.tsv` + `false_events.tsv` (normal mode),
#' tumour mode: `clones.tsv`, `assignments.tsv`, `clone_segments.tsv`, `clone_profiles.tsv.gz`,
#' `reference_pool.txt`; `regions.tsv`; `result.rds` (the object without scores and features), `params.yaml`,
#' `summary.tsv` ([scphasedcn_summary()]) and `NOTES.md`. Optional: `scores.tsv.gz`, `features.tsv.gz`.
#' @param obj `scphasedcn` from [call_scn()]
#' @param out output directory (created)
#' @param write_scores also write the per cell x region scores
#' @param write_features also write the features long table
#' @return named character vector of the written paths, invisibly
#' @export
write_scphasedcn = function(obj, out, write_scores = FALSE, write_features = FALSE) {
  obj = relabel_legacy(obj)
  if (!inherits(obj, "scphasedcn")) stop("write_scphasedcn: not a scphasedcn object")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  files = character()
  w = function(d, name) {
    if (is.null(d)) d = data.table()
    p = file.path(out, name)
    fwrite(as.data.table(d), p, sep = "\t")
    files[[name]] <<- p
  }
  calls = if (is.null(obj$calls)) NULL else as.data.table(obj$calls)
  w(obj$cells, "cells.tsv")
  if (!is.null(calls)) {
    cl = calls[level != "pool"]
    if (obj$mode == "normal") cl = cl[called %in% TRUE]
    w(cl, "calls.tsv")
    w(calls[level == "pool" & called %in% TRUE], "pooled_calls.tsv")
  } else { w(NULL, "calls.tsv"); w(NULL, "pooled_calls.tsv") }
  w(obj$events, "events.tsv")
  w(obj$thresholds, "thresholds.tsv")
  w(null_fit_table(obj$null_fit), "null_fit.tsv")
  w(obj$regions, "regions.tsv")
  if (obj$mode == "normal") {
    w(obj$false_calls, "false_calls.tsv")
    w(obj$false_events, "false_events.tsv")
  }
  if (obj$mode == "tumour" || !is.null(obj$clones)) {
    w(obj$clones, "clones.tsv")
    w(obj$assignments, "assignments.tsv")
    w(obj$segments, "clone_segments.tsv")
    pr = flatten_profiles(obj$profiles)
    p = file.path(out, "clone_profiles.tsv.gz")
    fwrite(if (is.null(pr)) data.table() else pr, p, sep = "\t", compress = "gzip")
    files[["clone_profiles.tsv.gz"]] = p
    asg = as.data.table(obj$assignments %||% data.table())
    ref = if ("ref_pool" %in% names(asg)) asg[ref_pool %in% TRUE, as.character(cell)] else character()
    p = file.path(out, "reference_pool.txt")
    writeLines(ref, p)
    files[["reference_pool.txt"]] = p
  }
  if (write_scores && !is.null(obj$scores)) {
    p = file.path(out, "scores.tsv.gz"); fwrite(obj$scores, p, sep = "\t", compress = "gzip"); files[["scores.tsv.gz"]] = p
  }
  if (write_features && !is.null(obj$features)) {
    p = file.path(out, "features.tsv.gz"); fwrite(obj$features, p, sep = "\t", compress = "gzip"); files[["features.tsv.gz"]] = p
  }
  sm = scphasedcn_summary(obj)
  w(sm, "summary.tsv")
  writeLines(c("# scphasedcn caller parameters", yaml_lines(obj$params)), file.path(out, "params.yaml"))
  files[["params.yaml"]] = file.path(out, "params.yaml")
  slim = obj
  slim["scores"] = list(NULL)
  slim["features"] = list(NULL)
  slim$files = c(files, result.rds = file.path(out, "result.rds"), NOTES.md = file.path(out, "NOTES.md"))
  saveRDS(slim, file.path(out, "result.rds"))
  files[["result.rds"]] = file.path(out, "result.rds")
  md_tab = function(d) {
    d = as.data.table(d)
    fmt = function(v) if (is.numeric(v)) ifelse(is.na(v), "", formatC(v, digits = 4, format = "g")) else ifelse(is.na(v), "", as.character(v))
    c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", strrep("---|", ncol(d))),
      if (nrow(d)) paste0("| ", do.call(paste, c(lapply(d, fmt), sep = " | ")), " |"))
  }
  notes = c(sprintf("# scphasedcn result (%s mode)", obj$mode), "",
            sprintf("Written %s on %s.", format(Sys.time(), "%F %T"), Sys.info()[["nodename"]]), "",
            "## Run notes", "", paste0("- ", obj$notes), "",
            "## Summary", "", md_tab(sm), "",
            "## Files", "", paste0("- `", names(files), "`"), "- `NOTES.md` (this file)", "",
            "`result.rds` holds the full object without `scores` and `features` (read it with `read_scphasedcn()`).",
            "`calls.tsv` holds the cell- and clone-level calls (normal mode: called rows only; the rows that were",
            "scored but not called are in `result.rds`); `pooled_calls.tsv` the called pool rows.")
  writeLines(notes, file.path(out, "NOTES.md"))
  files[["NOTES.md"]] = file.path(out, "NOTES.md")
  invisible(files)
}

# results written before the rename (working name scncall, 0.1.0-pre builds) carry the class "scncall":
# relabel it so that they read, write and print as scphasedcn objects
relabel_legacy = function(x) {
  if (!inherits(x, "scphasedcn") && inherits(x, "scncall")) class(x)[class(x) == "scncall"] = "scphasedcn"
  x
}

#' Read a caller result written by [write_scphasedcn()]
#' @param out output directory (or the path of its `result.rds`)
#' @return `scphasedcn`
#' @export
read_scphasedcn = function(out) {
  p = if (dir.exists(out)) file.path(out, "result.rds") else out
  if (!file.exists(p)) stop("read_scphasedcn: no result.rds in ", out)
  x = relabel_legacy(readRDS(p))
  if (!inherits(x, "scphasedcn")) stop("read_scphasedcn: ", p, " is not a scphasedcn object")
  x
}

#' Headline numbers of a caller result
#'
#' @param obj `scphasedcn`
#' @return data.table section, item, value, n (denominator or count), note. Sections: `cells` (cells, null
#'   cells, multiplets, cells with a reportable call / event), `calls` (called cell / pool / clone rows per
#'   level, class and state; `n` = reportable among them), `false_calls` (null false calls per 1,000 per class
#'   over all strata; `n` = null rows), `false_events` (false autosomal arm calls per null cell, any state and
#'   per state), `events` (reportable events per cell and per null cell, pooled calls), `clones` (clones,
#'   assigned cells, reference pool).
#' @export
scphasedcn_summary = function(obj) {
  out = list()
  add = function(section, item, value, n = NA_real_, note = NA_character_)
    out[[length(out) + 1]] <<- data.table(section = section, item = item, value = as.numeric(value), n = as.numeric(n), note = as.character(note))
  ct = as.data.table(obj$cells)
  nc = nrow(ct)
  add("cells", "cells", nc)
  if ("null_cell" %in% names(ct)) add("cells", "null_cells", sum(ct$null_cell %in% TRUE), nc)
  if ("multiplet" %in% names(ct)) add("cells", "multiplets", sum(ct$multiplet %in% TRUE), nc)
  if ("n_calls" %in% names(ct)) add("cells", "cells_with_reportable_call", sum(ct$n_calls > 0, na.rm = TRUE), nc)
  if ("n_events" %in% names(ct)) add("cells", "cells_with_event", sum(ct$n_events > 0, na.rm = TRUE), nc)
  calls = if (is.null(obj$calls)) data.table() else as.data.table(obj$calls)
  if (nrow(calls)) {
    s = calls[called %in% TRUE, .(called = .N, reportable = sum(reportable %in% TRUE)), by = .(level, class, state)]
    s = s[order(level, match(class, region_classes()), state)]
    for (i in seq_len(nrow(s))) add("calls", paste(s$level[i], s$class[i], s$state[i], sep = "|"), s$called[i], s$reportable[i], "value = called rows, n = reportable")
  }
  if (!is.null(obj$false_calls) && nrow(obj$false_calls)) {
    f = as.data.table(obj$false_calls)[stratum == "all"]
    for (i in seq_len(nrow(f))) add("false_calls", paste0(f$class[i], "|per_1000"), f$per_1000[i], f$n_null[i], "null rows, CV thresholds")
    fs = as.data.table(obj$false_calls)[stratum != "all"]
    if (nrow(fs)) add("false_calls", "class_x_stratum_range_per_1000", min(fs$per_1000), nrow(fs), sprintf("max %.3f", max(fs$per_1000)))
  }
  if (!is.null(obj$false_events) && nrow(obj$false_events)) {
    fe = as.data.table(obj$false_events)
    for (i in seq_len(nrow(fe))) add("false_events", paste0("arm_calls_per_null_cell|", fe$state[i]), fe$arm_calls_per_cell[i], fe$n_cells[i])
  }
  ev = if (is.null(obj$events)) data.table() else as.data.table(obj$events)
  if (nc > 0) {
    ce = if (nrow(ev)) ev[level == "cell"] else data.table(cell = character())
    add("events", "reportable_cell_events_per_cell", nrow(ce) / nc, nrow(ce))
    if ("null_cell" %in% names(ct) && any(ct$null_cell %in% TRUE)) {
      nn = ct[null_cell %in% TRUE, cell]
      add("events", "reportable_cell_events_per_null_cell", sum(ce$cell %in% nn) / length(nn), length(nn))
    }
  }
  if (nrow(calls)) {
    pc = calls[level == "pool"]
    if (nrow(pc)) {
      kind_of = ifelse(startsWith(pc$cell, "knn|"), "knn", "type")
      add("events", "pooled_calls", sum(pc$called %in% TRUE), nrow(pc), "called pool rows / pool rows")
      for (k in sort(unique(kind_of))) add("events", paste0("pooled_calls|", k), sum(pc$called[kind_of == k] %in% TRUE),
                                         sum(kind_of == k), if (k == "knn") "kNN pools (one per cell)" else "random pools within cell type")
    }
  }
  if (!is.null(obj$clones)) {
    cl = as.data.table(obj$clones)
    add("clones", "clones", if ("clone" %in% names(cl)) uniqueN(cl$clone) else nrow(cl))
    if ("assigned" %in% names(ct)) add("clones", "assigned_cells", sum(ct$assigned %in% TRUE), nc)
    if ("ref_pool" %in% names(ct)) add("clones", "reference_pool", sum(ct$ref_pool %in% TRUE), nc)
  }
  rbindlist(out)
}
utils::globalVariables(c("assigned", "state"))
