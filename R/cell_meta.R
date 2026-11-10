# cell_meta.R
# What: make_cell_meta(), the cell metadata table of call_scn() / `call_scn.R --cell-meta` from the extractor's
#       cells_qc.tsv plus optional user tables (cell types, known copy-number-flat cells, pseudobulk groups).
# Why: the caller needs only a few per-cell labels; this builds them from what a user has (e.g. RNA cell-type
#      annotations of the same barcodes) without any sample-specific script.
# Inputs: cells_qc.tsv (or an extraction directory), optional tables / vectors keyed by cell barcode.
# Outputs: data.table cell, qc_pass, cell_type, flat, group (written with data.table::fwrite for the CLI).

#' Cell metadata for the caller
#'
#' Builds the `cell_meta` table of [call_scn()] (the `--cell-meta` TSV of `call_scn.R`). Every column but `cell` is
#' optional for the caller:
#' * `qc_pass`: the cells to call. [call_scn()] (and `call_scn.R --cell-meta`) without `cells` calls the cells of the
#'   table with `qc_pass` TRUE; without the column, those of the table that pass [cell_qc()].
#' * `cell_type`: the depth reference of normal-tissue mode is matched by cell type (types with fewer than 20 cells,
#'   and cells with `NA`, use the reference of all cells). Without cell types the reference is the all-cell one.
#' * `flat`: cells known to be copy-number neutral (e.g. normal cells from an orthogonal assay). Normal mode uses them
#'   as the null cells with `null_cells = "features"`; without `flat`, use `null_cells = "select"` (`--null-cells
#'   select`): the caller selects its own null cells. Tumour mode fits its input depth / allelic null on them (when the
#'   features carry no null fit); its reference pool is not taken from them but is the normal clone the caller finds
#'   and cleans itself.
#' * `group`: tumour mode only, pseudobulk groups of the phase-block orientation (`groups = "group"`, `--groups
#'   group`), e.g. spatial regions or RNA clusters; NA = all-cell pseudobulk.
#'
#' Barcodes are matched as given (a trailing `-1` is stripped from every input, as the extractor does).
#' @param cells_qc `cells_qc.tsv` written by the extractor (`scphasedcn run`; path), its table, or the extraction directory; NULL = the
#'   cells of the other inputs
#' @param cell_types NULL, a path or table with columns `cell` and `cell_type` (else the first two columns), or a
#'   character vector named by cell
#' @param flat NULL, cell barcodes of known copy-number-flat cells (vector, or a file with one cell per line, without a
#'   header unless the first line is `cell`, or with a `cell` column)
#' @param groups NULL, a path or table with `cell` and `group` (else the first two columns), or a vector named by cell
#' @param qc_args named list of [cell_qc()] arguments for `qc_pass` (e.g. `list(min_frag = 10000)`)
#' @return data.table `cell`, `cell_type`, `flat`, `group`, plus `qc_pass` ([cell_qc()] on `cells_qc`) when
#'   `cells_qc` is given
#' @export
#' @examples
#' qc = data.table::data.table(cell = c("AAAC", "AAAG", "AAAT"), n_clean = c(40000, 25000, 900))
#' make_cell_meta(qc, cell_types = c(AAAC = "epithelial", AAAG = "immune"), flat = "AAAG")
make_cell_meta = function(cells_qc = NULL, cell_types = NULL, flat = NULL, groups = NULL, qc_args = list()) {
  strip = strip_barcode
  as_map = function(x, col) {
    if (is.null(x)) return(NULL)
    if (is.character(x) && length(x) == 1 && is.null(names(x)) && file.exists(x)) x = fread(x)
    if (is.data.frame(x)) {
      x = as.data.table(x)
      cc = if (all(c("cell", col) %in% names(x))) c("cell", col) else names(x)[1:2]
      return(setNames(as.character(x[[cc[2]]]), strip(x[[cc[1]]])))
    }
    if (is.null(names(x))) stop("make_cell_meta: `", col, "` needs cell names")
    setNames(as.character(x), strip(names(x)))
  }
  qc = NULL
  if (!is.null(cells_qc)) {
    if (is.character(cells_qc) && dir.exists(cells_qc)) cells_qc = file.path(cells_qc, "cells_qc.tsv")
    qc = if (is.character(cells_qc)) read_cells_qc(cells_qc) else as.data.table(cells_qc)
    if (!"cell" %in% names(qc)) setnames(qc, 1, "cell")
    qc[, cell := strip(cell)]
  }
  ct = as_map(cell_types, "cell_type")
  gr = as_map(groups, "group")
  fl = NULL
  if (!is.null(flat)) {
    if (is.character(flat) && length(flat) == 1 && file.exists(flat)) {
      d = read_cell_file(flat)                    # a one-column file has no header unless its first line is "cell"
      flat = if ("cell" %in% names(d)) d$cell else d[[1]]
    }
    fl = strip(flat)
  }
  cells = if (!is.null(qc)) qc$cell else unique(c(names(ct), names(gr), fl))
  if (!length(cells)) stop("make_cell_meta: no cells (give cells_qc or one of the tables)")
  out = data.table(cell = cells)
  if (!is.null(qc)) out[, qc_pass := if ("qc_pass" %in% names(qc)) as.logical(qc$qc_pass) else do.call(cell_qc, c(list(qc), qc_args))$qc_pass]
  out[, cell_type := if (is.null(ct)) NA_character_ else unname(ct[cell])]
  out[, flat := if (is.null(fl)) NA else cell %in% fl]
  out[, group := if (is.null(gr)) NA_character_ else unname(gr[cell])]
  miss = function(m) if (is.null(m)) 0L else sum(!names(m) %in% cells)
  n_miss = miss(ct) + miss(gr) + (if (is.null(fl)) 0L else sum(!fl %in% cells))
  if (n_miss > 0) message("make_cell_meta: ", n_miss, " labelled barcodes are not in cells_qc (dropped)")
  out[]
}
