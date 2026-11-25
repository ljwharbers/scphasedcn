#!/usr/bin/env Rscript
# plot_cells.R
# What: command-line front end of the single-cell track plots: scphasedcn::cell_tracks() on a caller run and its
#       features, then plot_cells_pdf() for every cell (one PDF per clone in tumour mode, per stratum in normal mode,
#       in chunks of --chunk cells).
# Why: one command per caller run, to look at the data behind each cell's calls.
# Usage: Rscript plot_cells.R --run <caller output dir> [--features <rds>] [--out <dir>] [options]
#        Without --features the feature step is recomputed (cell_features()) on the run's extraction directory
#        (--dir, default the `dir:` line of the run's NOTES.md), cells, sex (NOTES.md) and params.yaml.
#        RLIB (env) puts a package library first on the library path.
#        The script's path: system.file("scripts", "plot_cells.R", package = "scphasedcn").
# Inputs: the caller run (result.rds, NOTES.md, params.yaml); a features .rds or the extraction directory.
# Outputs: <out>/<group>[_<k>].pdf, <out>/cells.tsv (cell -> file, page), optionally <out>/tracks.rds.

defaults = list(run = "", tracks = "", features = "", dir = "", sex = "", out = "", cells = "", per_page = 4,
                chunk = 500, agg = 5e6, orient = "auto", cores = 8, save_tracks = FALSE, min_hap = 3)
usage = paste(
  "Single-cell depth / BAF track plots (scphasedcn::cell_tracks, plot_cells_pdf) of one caller run.",
  "Needs the R packages ggplot2 and patchwork (ggrastr, when installed, rasterises the 1 Mb points).",
  "",
  "  --run <dir>          caller output directory (result.rds)",
  "  --tracks <rds>       tracks.rds of an earlier run of this script (--save-tracks): re-plot without features",
  "  --features <rds>     the run's scn_features .rds (default: recompute the feature step)",
  "  --dir <path>         extraction directory for the recompute (default: the `dir:` line of <run>/NOTES.md)",
  "  --sex female|male    sample sex for the recompute (default: the `sex:` line of <run>/NOTES.md)",
  "  --out <dir>          output directory (default <run>/cell_plots)",
  "  --cells <file|list>  cells to plot, a file (one per line) or a ';' list (default all cells of the run)",
  "  --per-page <n>       cells per page",
  "  --chunk <n>          cells per PDF",
  "  --agg <bp>           aggregated window size (depth points, BAF)",
  "  --orient <x>         auto | clone | pseudobulk | none (BAF phase-block orientation)",
  "  --min-hap <n>        minimal molecule votes of a BAF point",
  "  --cores <n>          data.table threads",
  "  --save-tracks        also write tracks.rds",
  sep = "\n")
if (any(commandArgs(trailingOnly = TRUE) %in% c("-h", "--help"))) {
  cat(usage, "\n\noptions (defaults):\n", sep = "")
  for (k in names(defaults)) cat(sprintf("  --%s (default: %s)\n", gsub("_", "-", k), paste(format(defaults[[k]]), collapse = ",")))
  quit(save = "no", status = 0)
}

local({ l = Sys.getenv("RLIB"); if (nzchar(l)) .libPaths(c(l, .libPaths())) })
suppressPackageStartupMessages({ library(scphasedcn); library(data.table) })
fail = function(...) { message("plot_cells.R: ", ..., "\nUsage: plot_cells.R --run <caller output dir> [options]; see --help"); quit(save = "no", status = 1) }
a = tryCatch(cli_args(defaults = defaults, usage = usage), error = function(e) fail(conditionMessage(e)))
unknown = setdiff(names(a), names(defaults))
if (length(unknown)) fail("unknown option(s) ", paste0("--", gsub("_", "-", unknown), collapse = ", "))
if (!is.character(a$run) || !nzchar(a$run)) fail("--run is required")
if (!a$orient %in% c("auto", "clone", "pseudobulk", "none")) fail("--orient must be auto, clone, pseudobulk or none")
for (k in c("per_page", "chunk", "agg", "cores", "min_hap")) {
  v = a[[k]]
  ok = is.numeric(v) && length(v) == 1 && is.finite(v) && (if (k == "min_hap") v >= 0 else v > 0)
  if (!ok) fail(sprintf("--%s must be a %s number", gsub("_", "-", k), if (k == "min_hap") "non-negative" else "positive"))
}
if (!is.logical(a$save_tracks) || is.na(a$save_tracks)) fail("--save-tracks takes no value (or TRUE / FALSE)")
for (k in c("tracks", "features", "dir", "sex", "out", "cells"))
  if (!is.character(a[[k]])) fail(sprintf("--%s needs a value", k))
if (nzchar(a$tracks) && !file.exists(a$tracks)) fail("--tracks file not found: ", a$tracks)
if (nzchar(a$features) && !file.exists(a$features)) fail("--features file not found: ", a$features)
setDTthreads(as.integer(a$cores))
res = tryCatch(read_scphasedcn(a$run), error = function(e) fail(conditionMessage(e)))
out = if (nzchar(a$out)) a$out else file.path(a$run, "cell_plots")
notes = if (file.exists(file.path(a$run, "NOTES.md"))) readLines(file.path(a$run, "NOTES.md")) else character()
note_val = function(key) { l = grep(sprintf("^- %s: ", key), notes, value = TRUE); if (length(l)) sub(sprintf("^- %s: ", key), "", l[1]) else "" }
for (p in c("ggplot2", "patchwork")) if (!requireNamespace(p, quietly = TRUE)) fail("needs the R package ", p)

if (nzchar(a$tracks)) {
  cli_log("reading tracks ", a$tracks)
  feat = NULL
} else if (nzchar(a$features)) {
  cli_log("reading features ", a$features)
  feat = readRDS(a$features)
} else {
  dir = if (nzchar(a$dir)) a$dir else note_val("dir")
  sex = if (nzchar(a$sex)) a$sex else note_val("sex")
  if (identical(sex, "from cells")) sex = ""   # the run took the sex from the cells table: cell_features() does too
  if (!nzchar(dir) || !dir.exists(dir)) fail("no extraction directory: give --dir (NOTES.md has none, or it is gone: ", dir, ")")
  pf = file.path(a$run, "params.yaml")
  params = scphasedcn:::split_params(caller_params(if (file.exists(pf)) pf else NULL))$params
  ct = as.data.table(res$cells)
  cm = ct[, intersect(c("cell", "cell_type", "group", "flat", "fold", "stratum", "n_frag"), names(ct)), with = FALSE]
  method = if (res$mode == "normal") "celltype" else "model"
  if (method == "celltype" && (!"cell_type" %in% names(cm) || all(is.na(cm$cell_type)))) cm[, cell_type := "all"]
  cli_log(sprintf("cell_features (%s depth) on %s: %d cells, sex %s", method, dir, nrow(ct), sex))
  feat = cell_features(dir, cells = ct$cell, cell_meta = cm, method = method, sex = if (nzchar(sex)) sex else NULL,
                       params = params, verbose = TRUE)
}

cells = if (nzchar(a$cells)) {
  if (file.exists(a$cells)) readLines(a$cells) else strsplit(a$cells, ";", fixed = TRUE)[[1]]
} else NULL
if (is.null(feat)) {
  tr = readRDS(a$tracks)
  if (!inherits(tr, "cell_tracks")) fail("--tracks is not a cell_tracks .rds: ", a$tracks)
  fcells = tr$meta$cell
} else fcells = colnames(feat$win$x)
if (is.null(cells)) cells = intersect(as.character(res$cells$cell), fcells)
drop = setdiff(cells, fcells)
if (length(drop)) { cli_log(sprintf("%d cells are not in the features / tracks; skipped", length(drop))); cells = setdiff(cells, drop) }
if (!length(cells)) fail("no cells to plot")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (!is.null(feat)) {
  cli_log(sprintf("cell_tracks: %d cells", length(cells)))
  tr = cell_tracks(feat, res, cells = cells, agg = a$agg, orient = a$orient)
  rm(feat); invisible(gc())
}
tr$meta = tr$meta[cell %in% cells]
cli_log(sprintf("tracks: %d depth, %d BAF, %d segment rows; orientation %s", nrow(tr$depth), nrow(tr$baf), nrow(tr$fit),
                attr(tr, "orient")))
if (isTRUE(a$save_tracks) && !nzchar(a$tracks)) saveRDS(tr, file.path(out, "tracks.rds"))

# groups: clone (tumour) or stratum (normal); cells by decreasing fragment count
m = copy(tr$meta)
m[, grp := if (res$mode == "tumour") fifelse(is.na(clone), "noclone", as.character(clone)) else
  fifelse(is.na(stratum), "all", as.character(stratum))]
m[, grp := gsub("[^A-Za-z0-9_.-]", "_", sub("^>=", "ge", sub("^<", "lt", grp)))]   # strata "<20k" -> lt20k, ">=50k" -> ge50k
m = m[order(grp, -n_frag)]
m[, k := (seq_len(.N) - 1L) %/% a$chunk + 1L, by = grp]
m[, n_k := max(k), by = grp]
m[, file := fifelse(n_k > 1, sprintf("%s_%02d.pdf", grp, k), sprintf("%s.pdf", grp))]
m[, page := (seq_len(.N) - 1L) %/% a$per_page + 1L, by = file]
for (f in unique(m$file)) {
  t0 = Sys.time()
  cc = m[file == f, cell]
  plot_cells_pdf(tr, file.path(out, f), cells = cc, per_page = a$per_page, min_hap = a$min_hap)
  cli_log(sprintf("%s: %d cells, %.1f MB, %.1f min", f, length(cc), file.size(file.path(out, f)) / 2^20,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
fwrite(m[, .(cell, clone, stratum, n_frag, n_confirmed, n_inherited, file, page)], file.path(out, "cells.tsv"), sep = "\t")
cli_log("done: ", out)
