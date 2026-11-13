#!/usr/bin/env Rscript
# call_scn.R
# What: command-line front end of the copy-number caller: scphasedcn::call_scn() on one extraction directory (written by
#       `scphasedcn run`), then write_scphasedcn() into the output directory.
# Why: one command per sample; the R API stays the single place where the caller is defined.
# Usage: Rscript call_scn.R --dir <extract dir> --out <dir> --sex female|male [--mode normal|tumour] [options]
#        or, with a sample yaml (see system.file("extdata", "sample_template.yaml", package = "scphasedcn")):
#        Rscript call_scn.R --config <sample yaml> --sample <key> [--tag <tag>] [options]
#        Caller parameters: --caller-config (default the package's caller.yaml, incl. its normal: / tumour: blocks);
#        --seed, --n-draws, --confirm-null and --split-pass override single values of it.
#        Lists (--cells, --null-cells, --multiplets given inline, --events) use ';' as separator (so that values
#        survive tools that split at commas, e.g. sbatch --export). RLIB (env) puts a package library first on the library path.
#        The script's path: system.file("scripts", "call_scn.R", package = "scphasedcn").
# Inputs: the extraction directory (bins_100kb.h5ad, fragments.tsv.gz + .tbi, cells_qc.tsv), optional cell metadata,
#         an optional sample yaml and caller yaml.
# Outputs: <out>/ cells.tsv, calls.tsv, pooled_calls.tsv, events.tsv, thresholds.tsv, null_fit.tsv, false_calls.tsv,
#          (tumour: clones.tsv, assignments.tsv, clone_segments.tsv, clone_profiles.tsv.gz, reference_pool.txt),
#          result.rds, params.yaml, summary.tsv, NOTES.md.

defaults = list(sample = "", mode = "", dir = "", sex = "", events = "", config = "", caller_config = "", cell_meta = "",
                cells = "", null_cells = "auto", multiplets = "", groups = "", control_regions = "", features = "",
                centromeres = "", out = "", tag = "default", cores = 8, write_scores = FALSE, write_features = FALSE,
                no_embedding = FALSE, seed = "", n_draws = "", confirm_null = "", split_pass = "", phasing = "auto")
usage = paste(
  "Single-cell copy-number caller (scphasedcn::call_scn) on one sample.",
  "",
  "Inputs: either --dir + --out (+ --sex), or --config + --sample.",
  "  --dir <path>            extraction directory written by `scphasedcn run` (bins_100kb.h5ad, fragments.tsv.gz + .tbi,",
  "                          cells_qc.tsv)",
  "  --out <dir>             output directory (with --config / --sample: default <outputs.root>/<sample>/caller/<tag>)",
  "  --sex female|male       sample sex (required unless the sample yaml gives it; also F / M, any case)",
  "  --config <yaml>         sample yaml (layout: system.file('extdata', 'sample_template.yaml', package = 'scphasedcn');",
  "                          default: the SCPHASEDCN_SAMPLE_CONFIG environment variable)",
  "  --sample <key>          sample of the --config yaml (extraction dir, sex, events, type)",
  "  --mode normal|tumour    default: the sample's `type` (tumour / tumor; normal when missing)",
  "  --events <list>         known events, name=region;name=region (e.g. chr5_gain=chr5;20q_gain=chr20:28000000-;",
  "                          0-based half-open, open end = chromosome end); a bare region is named after itself",
  "  --centromeres <tsv>     centromere table (default: the sample yaml's reference.centromeres, else the hg38 table",
  "                          shipped with the package)",
  "  --caller-config <yaml>  caller parameters (default: the package's caller.yaml; its normal: / tumour: blocks are",
  "                          the caller_normal() arguments / tumour-mode options clone_defaults())",
  "  --cell-meta <tsv>       cell + optional qc_pass/cell_type/flat/fold/group/stratum/n_frag (see make_cell_meta();",
  "                          only its cells that pass the cell QC are called, and qc_pass TRUE when it has that column;",
  "                          cell_type 'unknown' = no type)",
  "  --cells <file|list>     cells to call, without cell QC (default: QC pass, and in cell-meta when given)",
  "  --null-cells <x>        auto | features | select | a file or ';' list of null cells",
  "  --multiplets <x>        file (one cell per line, or a table with prob / multiplet / clean_singlet) or ';' list",
  "  --groups <column>       tumour mode: cell-meta column of the pseudobulk orientation groups (e.g. a spatial region)",
  "  --control-regions <tsv> tumour mode: control regions table",
  "  --phasing <x>           auto | chromosome | block | loo | mirrored: block orientation of the allele features",
  "                          (default auto: chromosome-wide phase sets as they are, else loo; normal mode:",
  "                          chromosome-scale phasing recommended, see `scphasedcn hybrid-phase`)",
  "  --features <rds>        precomputed scn_features (skips the feature step)",
  "  --tag <tag>             output tag (with --config / --sample)",
  "  --cores <n>             data.table threads",
  "  --write-scores          also write scores.tsv.gz (large)",
  "  --write-features        also write features.tsv.gz (large)",
  "  --no-embedding          no chromatin embedding (no kNN pools)",
  "  --seed <n>              RNG seed (overrides the caller config's seed)",
  "  --n-draws <n>           normal mode: null pool-regions per class for the pooled thresholds (overrides normal: n_draws)",
  "  --confirm-null <x>      tumour mode: matched | generic per-cell confirmation null (overrides tumour: confirm_null)",
  "  --split-pass <x>        tumour mode: FALSE | disagree final clone split pass (overrides tumour: split_pass)",
  sep = "\n")
if (any(commandArgs(trailingOnly = TRUE) %in% c("-h", "--help"))) {
  cat(usage, "\n\noptions (defaults):\n", sep = "")
  for (k in names(defaults)) cat(sprintf("  --%s (default: %s)\n", gsub("_", "-", k), paste(format(defaults[[k]]), collapse = ",")))
  quit(save = "no", status = 0)
}

RLIB = Sys.getenv("RLIB")
if (nzchar(RLIB)) .libPaths(c(RLIB, .libPaths()))
suppressPackageStartupMessages({ library(scphasedcn); library(data.table) })
fail = function(...) { message("call_scn.R: ", ...); quit(save = "no", status = 1) }
a = cli_args(defaults = defaults, usage = usage)
unknown = setdiff(names(a), names(defaults))
if (length(unknown)) fail("unknown option(s) ", paste0("--", gsub("_", "-", unknown), collapse = ", "), "; see --help")
t0 = Sys.time()
setDTthreads(as.integer(a$cores))
split_list = function(x) if (!nzchar(x)) NULL else if (file.exists(x)) x else trimws(strsplit(x, ";", fixed = TRUE)[[1]])

# sample yaml: --config, else the SCPHASEDCN_SAMPLE_CONFIG environment variable (a default yaml for batch wrappers)
cfg = NULL
cfg_file = if (nzchar(a$config)) a$config else Sys.getenv("SCPHASEDCN_SAMPLE_CONFIG")
if (nzchar(cfg_file)) {
  if (!file.exists(cfg_file)) fail("sample yaml not found: ", cfg_file)
  cfg = read_config(cfg_file)
} else if (nzchar(a$sample)) fail("--sample needs --config <sample yaml>")
s = if (nzchar(a$sample)) cli_sample(cfg, a$sample) else list(name = "sample")
dir = if (nzchar(a$dir)) a$dir else s$dir
if (is.null(dir) || !nzchar(dir)) fail("give --dir <extraction dir> (or --config + --sample)")
if (!nzchar(a$features) && !file.exists(file.path(dir, "bins_100kb.h5ad"))) fail("no bins_100kb.h5ad in ", dir)
mode = if (nzchar(a$mode)) a$mode else if (isTRUE(tolower(s$type) %in% c("tumour", "tumor"))) "tumour" else "normal"
if (identical(tolower(mode), "tumor")) mode = "tumour"
if (!mode %in% c("normal", "tumour")) fail("--mode must be normal or tumour")
if (!a$phasing %in% c("auto", "chromosome", "block", "loo", "mirrored")) fail("--phasing must be auto, chromosome, block, loo or mirrored")
sex = if (nzchar(a$sex)) a$sex else s$sex
if (!is.null(sex)) sex = tryCatch(scphasedcn:::check_sex(sex, "--sex"), error = function(e) fail(conditionMessage(e)))
# events: name=region or a bare region (named after the region, e.g. chr7 -> event|chr7); 0-based half-open coordinates
events = if (nzchar(a$events)) scphasedcn:::parse_event_arg(a$events) else if (!is.null(s$events)) unlist(s$events) else NULL
cen = if (nzchar(a$centromeres)) a$centromeres else cfg$reference$centromeres
if (!is.null(cen) && !file.exists(cen)) fail("centromere table not found: ", cen)
params = caller_params(if (nzchar(a$caller_config)) a$caller_config else NULL)
if (nzchar(a$seed)) params$seed = as.numeric(a$seed)
# caller arguments from the flags (over the caller config's normal: / tumour: blocks, merged in call_scn)
caller_args = if (mode == "normal") {
  if (nzchar(a$n_draws)) list(n_draws = as.numeric(a$n_draws)) else list()
} else {
  o = list()
  if (nzchar(a$confirm_null)) o$confirm_null = a$confirm_null
  if (nzchar(a$split_pass)) o$split_pass = if (toupper(a$split_pass) %in% c("TRUE", "FALSE")) as.logical(a$split_pass) else a$split_pass
  if (length(o)) list(opts = o) else list()
}
out = if (nzchar(a$out)) a$out else if (nzchar(a$sample) && length(cfg$outputs$root)) {
  file.path(cfg$outputs$root, s$name, "caller", a$tag)
} else fail("give --out <dir>")

cm = NULL
if (nzchar(a$cell_meta)) {
  cm = fread(a$cell_meta)
  if (!"cell" %in% names(cm)) fail("--cell-meta needs a `cell` column")
  if ("cell_type" %in% names(cm)) cm[cell_type %in% c("unknown", ""), cell_type := NA]
  if ("stratum" %in% names(cm)) cm[, stratum := as.character(stratum)]
  keep = intersect(c("cell", "qc_pass", "cell_type", "flat", "fold", "group", "stratum", "n_frag", a$groups), names(cm))
  cm = cm[, keep, with = FALSE]
}
# --cells bypasses the cell QC; otherwise call_scn() selects the cells of --cell-meta that pass cell_qc() on
# cells_qc.tsv and, when the table has a qc_pass column, have qc_pass TRUE, exactly as the R API does with the same table
cells = split_list(a$cells)
nc = a$null_cells
null_cells = if (nc %in% c("auto", "features", "select")) nc else split_list(nc)
mult = split_list(a$multiplets)
ctrl = if (nzchar(a$control_regions)) fread(a$control_regions) else NULL
cli_log("call_scn: sample ", s$name, ", mode ", mode, ", dir ", dir, ", out ", out)

# with --features the feature step is skipped: its inputs (sex, centromeres) are not passed
pre = nzchar(a$features)
res = call_scn(dir, mode = mode, cell_meta = cm, cells = cells, sex = if (pre) NULL else sex, events = events, null_cells = null_cells,
               multiplets = mult, control_regions = ctrl, groups = if (nzchar(a$groups)) a$groups else NULL, params = params,
               features = if (pre) a$features else NULL, embedding = !isTRUE(a$no_embedding),
               keep_scores = isTRUE(a$write_scores), caller_args = caller_args, phasing = a$phasing,
               centromeres = if (pre) NULL else cen)
print(res)
files = write_scphasedcn(res, out, write_scores = isTRUE(a$write_scores), write_features = isTRUE(a$write_features))
print(scphasedcn_summary(res))
cli_log("wrote ", length(files), " files to ", out, " (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
