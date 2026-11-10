# utils.R
# What: small shared helpers (null-coalescing, chromosome ordering, region
#       parsing, bin lookup, sparse/dense coercions, group indicators, cell
#       selectors, fragment haplotype codes, budget thresholds, the sample sex check).
# Why: every module needs the same handful of primitives; keeping them in one
#      file avoids divergent re-implementations.
# Inputs: plain vectors / data.tables / Matrix objects. Outputs: same.

`%||%` = function(a, b) if (is.null(a)) b else a

# sample sex of the caller's entry points: "female" / "male", case-insensitive, or "F" / "M"; NULL stays NULL; any
# other value stops (a "Male" read as female would make chrX allelic, drop chrY and assume ploidy 2)
check_sex = function(sex, where = "scphasedcn") {
  if (is.null(sex)) return(NULL)
  s = tolower(trimws(as.character(unlist(sex))))
  out = if (length(s) == 1 && !is.na(s)) c(female = "female", f = "female", male = "male", m = "male")[s] else NA
  if (is.na(out)) stop(sprintf("%s: `sex` must be \"female\" or \"male\" (or \"F\" / \"M\"), not %s", where,
                               if (length(sex) == 1) sprintf("\"%s\"", as.character(unlist(sex))) else sprintf("a vector of length %d", length(sex))),
                       call. = FALSE)
  unname(out)
}

#' Order chromosome names naturally (chr1..chr22, chrX, chrY, chrM, others)
#' @param x character vector of chromosome names
#' @return integer ranks usable with `order()`
#' @export
chrom_rank = function(x) {
  base = sub("^chr", "", x)
  num = suppressWarnings(as.integer(base))
  r = ifelse(!is.na(num), num,
             ifelse(base == "X", 100L, ifelse(base == "Y", 101L, ifelse(base == "M" | base == "MT", 102L, 200L))))
  r + match(x, sort(unique(x))) / 1e6
}

#' Parse a region string "chr20:28000000-" or "chr5" into a data.table
#' @param x character vector of regions (`chrom`, `chrom:start-end`, `chrom:start-`)
#' @param chrom_sizes optional named integer vector to fill open ends
#' @return data.table with columns region, chrom, start, end (0-based half-open; NA end = chromosome end)
#' @export
parse_region = function(x, chrom_sizes = NULL) {
  out = lapply(x, function(r) {
    parts = strsplit(r, ":", fixed = TRUE)[[1]]
    chrom = parts[1]
    s = 0; e = NA_real_
    if (length(parts) > 1) {
      se = strsplit(parts[2], "-", fixed = TRUE)[[1]]
      s = as.numeric(gsub("[^0-9]", "", se[1]))
      if (length(se) > 1 && nzchar(se[2])) e = as.numeric(gsub("[^0-9]", "", se[2]))
    }
    if (is.na(e) && !is.null(chrom_sizes) && chrom %in% names(chrom_sizes)) e = as.numeric(chrom_sizes[[chrom]])
    data.table(region = r, chrom = chrom, start = s, end = e)
  })
  rbindlist(out)
}

#' Index of the bin containing each position
#'
#' Bins must tile each chromosome in order (as produced by [make_bins()] or
#' the Python extraction layer). Positions outside any bin give NA.
#' @param bins data.table with chrom, start, end (sorted by chrom then start)
#' @param chrom character vector
#' @param pos numeric vector (0-based)
#' @return integer vector of row indices into `bins`
#' @export
locate_bins = function(bins, chrom, pos) {
  idx = rep(NA_integer_, length(pos))
  starts_by = split(seq_len(nrow(bins)), bins$chrom)
  for (ch in unique(chrom)) {
    rows = starts_by[[ch]]
    if (is.null(rows)) next
    sel = which(chrom == ch)
    k = findInterval(pos[sel], bins$start[rows])
    ok = k >= 1 & pos[sel] < bins$end[rows][pmax(k, 1)]
    idx[sel[ok]] = rows[k[ok]]
  }
  idx
}

#' Bins overlapping a set of regions
#' @param bins data.table with chrom, start, end
#' @param regions data.table from [parse_region()] (NA end = to chromosome end)
#' @return list (one per region) of integer bin indices
#' @export
region_bins = function(bins, regions) {
  lapply(seq_len(nrow(regions)), function(r) {
    e = regions$end[r]; if (is.na(e)) e = Inf
    which(bins$chrom == regions$chrom[r] & bins$end > regions$start[r] & bins$start < e)
  })
}

#' Coerce to a bins x cells dgCMatrix
#' @param m matrix-like
#' @return dgCMatrix
#' @export
as_dgc = function(m) {
  if (is(m, "dgCMatrix")) return(m)
  as(as(as(Matrix::Matrix(m, sparse = TRUE), "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

#' Dense double matrix from any matrix-like
#' @param m matrix-like
#' @return base matrix
#' @export
dense = function(m) {
  if (is.matrix(m)) return(m)
  as.matrix(m)
}

#' Aggregation matrix mapping rows (bins) to groups
#' @param groups integer/character group id per row (NA = dropped)
#' @return sparse groups x rows indicator matrix (dgCMatrix) with group labels as rownames
#' @export
group_indicator = function(groups) {
  ok = !is.na(groups)
  g = factor(groups[ok], levels = unique(groups[ok]))
  m = sparseMatrix(i = as.integer(g), j = which(ok), x = 1, dims = c(nlevels(g), length(groups)))
  rownames(m) = levels(g)
  m
}

chrom_sizes_from_bins = function(bins) {
  s = bins[, .(len = max(end)), by = chrom]
  setNames(s$len, s$chrom)
}

# column indices of an scn object for a cell selector (NULL = all, logical, character names or integer indices)
cell_selector = function(scn, cells) {
  n = ncol(scn$counts)
  if (is.null(cells)) return(seq_len(n))
  if (is.logical(cells)) return(which(cells))
  if (is.character(cells)) return(match(cells, scn$cells$cell))
  as.integer(cells)
}

#' Haplotype codes used in the fragments `hap` column
#'
#' The Python extraction writes `hap` as 0 = no SNP (unphased), 1 = hap A,
#' 2 = hap B, 3 = ambiguous. [hap_code()] / [hap_label()] convert between the
#' integer codes and the labels "U" / "A" / "B" / "N" used in R.
#' @param hap character labels
#' @return integer codes
#' @export
hap_code = function(hap) {
  code = c(U = 0L, A = 1L, B = 2L, N = 3L)[as.character(hap)]
  code[is.na(code)] = 0L
  as.integer(code)
}

#' @rdname hap_code
#' @param code integer codes (or already labels, returned unchanged)
#' @export
hap_label = function(code) {
  if (is.character(code)) return(code)
  c("U", "A", "B", "N")[pmin(pmax(as.integer(code), 0L), 3L) + 1L]
}

#' Threshold at a false-call budget from null scores
#' @param null_scores scores of null (event-free) cell-regions
#' @param budget_per_1000 allowed false calls per 1,000 cell-regions
#' @return threshold (type-7 quantile at 1 - budget / 1000)
#' @export
budget_threshold = function(null_scores, budget_per_1000 = 1) {
  x = null_scores[is.finite(null_scores)]
  unname(quantile(x, 1 - budget_per_1000 / 1000))
}

# ---- RNG state -------------------------------------------------------------------------------------------
# The entry points (call_scn, cell_features, caller_normal, caller_tumour, simulate_scn) seed the RNG internally
# for reproducibility. rng_save() / rng_restore() put the caller's global RNG state back when the outermost entry
# point returns; nested entry points leave the stream alone, so a run draws exactly the numbers it always drew.
.scphasedcn_rng = new.env(parent = emptyenv())

rng_save = function() {
  if (isTRUE(.scphasedcn_rng$active)) return(NULL)
  .scphasedcn_rng$active = TRUE
  had = exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  list(had = had, seed = if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL)
}

rng_restore = function(state) {
  if (is.null(state)) return(invisible())
  .scphasedcn_rng$active = FALSE
  if (state$had) assign(".Random.seed", state$seed, envir = globalenv())
  else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  invisible()
}
