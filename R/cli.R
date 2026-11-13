# cli.R
# What: helpers for the thin command-line wrappers (inst/scripts/call_scn.R and user scripts): argument parsing
#       without argparse, sample lookup in a sample yaml, logging, writing tables, loading cell subsets.
# Why: every wrapper needs the same ten lines; keeping them here keeps the scripts to the analysis logic.
# Inputs: commandArgs(), a sample yaml (inst/extdata/sample_template.yaml layout). Outputs: lists / files.

#' Parse `--key value` and `--flag` arguments
#' @param args character vector (default `commandArgs(trailingOnly = TRUE)`)
#' @param defaults named list of defaults; values are coerced to the default's type
#' @param usage text printed with `--help`
#' @return named list. With `--help` the usage is printed and a script (`Rscript`) quits; an interactive session gets
#'   the defaults back with attribute `help = TRUE`
#' @examples
#' cli_args(c("--dir", "extract", "--cores", "4", "--no-embedding"),
#'          defaults = list(dir = "", cores = 8))
#' @export
cli_args = function(args = commandArgs(trailingOnly = TRUE), defaults = list(), usage = NULL) {
  out = defaults
  if (any(args %in% c("-h", "--help"))) {
    if (!is.null(usage)) cat(usage, "\n")
    cat("options:\n"); for (k in names(defaults)) cat(sprintf("  --%s (default: %s)\n", k, paste(format(defaults[[k]]), collapse = ",")))
    if (interactive()) return(invisible(structure(defaults, help = TRUE)))
    quit(save = "no", status = 0)
  }
  i = 1
  while (i <= length(args)) {
    a = args[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a)
    key = sub("^--", "", a)
    key = gsub("-", "_", key)
    if (i == length(args) || startsWith(args[i + 1], "--")) { out[[key]] = TRUE; i = i + 1; next }
    val = args[i + 1]
    def = defaults[[key]]
    out[[key]] = if (is.numeric(def)) as.numeric(val) else if (is.logical(def)) as.logical(val) else val
    i = i + 2
  }
  out
}

#' Sample entry of a sample yaml (with its extraction directory)
#'
#' The extraction directory is the sample's `extract_dir`, else `<outputs.root>/<sample>/extract` (where
#' `scphasedcn run` writes it).
#' @param cfg sample yaml read with [read_config()]
#' @param sample sample name
#' @return list with the sample fields plus `name`, `dir` (extraction output dir)
#' @export
cli_sample = function(cfg, sample) {
  s = cfg$samples[[sample]]
  if (is.null(s)) stop("sample not in config: ", sample)
  s$name = sample
  s$dir = s$extract_dir %||% file.path(cfg$outputs$root, sample, "extract")
  s
}

#' Timestamped message
#' @param ... pasted
#' @export
cli_log = function(...) message(format(Sys.time(), "%H:%M:%S"), " ", ...)

#' Write a data.table as TSV into the output directory
#' @param d data.table
#' @param out output directory (created)
#' @param name file name
#' @return path
#' @export
cli_write = function(d, out, name) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  p = file.path(out, name)
  fwrite(d, p, sep = "\t")
  cli_log("wrote ", p, " (", nrow(d), " rows)")
  invisible(p)
}

#' Cell subset from a file (one cell per line or a TSV with a `cell` column)
#'
#' A one-column file has no header unless its first line is `cell` (as [call_scn()] reads its cell files); barcodes
#' match with or without a trailing `-1`.
#' @param path file path or NULL
#' @param all all cell names (returned when `path` is NULL)
#' @return character vector
#' @export
cli_cells = function(path, all) {
  if (is.null(path) || isFALSE(path) || !nzchar(path)) return(all)
  d = read_cell_file(path)
  v = if ("cell" %in% names(d)) d$cell else d[[1]]
  all[strip_barcode(all) %in% strip_barcode(v)]                    # a trailing -1 on either side does not matter
}

#' Event regions of a sample (config `events` map, region strings with 0-based half-open coordinates) as a regions table
#' @param s sample entry from [cli_sample()]
#' @param chrom_sizes named sizes (closes open ends)
#' @return data.table region, chrom, start, end, event
#' @export
cli_events = function(s, chrom_sizes = NULL) {
  if (is.null(s$events) || length(s$events) == 0) return(NULL)
  r = parse_region(unlist(s$events), chrom_sizes)
  r[, event := names(s$events)]
  r[]
}

# --events of call_scn.R: "name=region;region;..." -> region strings named by event ("" for a bare region, which
# event_table() names after the region)
parse_event_arg = function(x) {
  ev = trimws(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE)[[1]])
  ev = ev[nzchar(ev)]
  named = grepl("=", ev, fixed = TRUE)
  stats::setNames(trimws(ifelse(named, sub("^[^=]*=", "", ev), ev)), trimws(ifelse(named, sub("=.*$", "", ev), "")))
}

#' Arm regions from an annotated bins table
#' @param bins bins with chrom and arm
#' @return regions table (region = chrom + arm)
#' @export
arm_regions = function(bins) {
  b = bins[!is.na(arm) & arm %in% c("p", "q")]
  r = b[, .(start = min(start), end = max(end)), by = .(chrom, arm)]
  r[, region := paste0(chrom, ":", start, "-", end)]
  r[, name := paste0(chrom, arm)]
  r[, .(region, chrom, start, end, name)]
}
