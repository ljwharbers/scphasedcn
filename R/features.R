# features.R
# What: the feature step of the copy-number caller: per-cell per-region depth and
#       allele features from an extraction directory (bins_100kb.h5ad + fragments.tsv.gz). Cell QC,
#       Battenberg arms, the regions table (arms, chromosomes, 5-50 Mb
#       tiles, named events; containment rule for depth, start-inside rule for alleles, the sd_class
#       rule), depth x / e / log2ratio / per-cell null SD / z by sparse region x bin products, molecule
#       votes (fragments of one cell within `mol_gap` collapsed, majority haplotype) counted per atom
#       and phase set, phase-block orientation (chromosome-wide anchored phase, leave-one-out
#       pseudobulk / clone orientation, or mirrored), and the long features table.
# Why: one function, cell_features(), builds every feature table that the scorer and the
#      tumour mode consume (schemas in vignette("scphasedcn-outputs", package = "scphasedcn")).
# Inputs: an extraction directory (bins_100kb.h5ad, fragments.tsv.gz, cells_qc.tsv; read_sample_depth()), a centromere table,
#         cell metadata.
# Outputs: an `scn_features` list: long (features), regions, cells, votes, atoms, win, params, notes.

utils::globalVariables(c("qc_pass", "stratum", "cen_left", "cen_right", "cstart", "cend", "csz", "sd_class",
                         "null_region", "allelic_ok", "n_use", "atom", "ri", "layer", "isA", "isB", "mol", "new",
                         "gA", "gB", "rA", "rB", "grp", "n_blocks", "n_hap", "f_mirror", "fold", "null_cell",
                         "nA_maj", "region_name", "ci", "hA", "hB", "i.gA", "i.gB", "i.hA", "i.hB",
                         "i.f_mirror", "i.n_blocks", "flags"))

AUTOSOMES_CALLER = paste0("chr", 1:22)

# ---- cells ---------------------------------------------------------------------------------------

#' Cell QC of the caller
#'
#' A cell passes when it has at least `min_frag` clean fragments, a fraction of MAPQ-60 fragments of
#' at least `min_mapq60` and at most `max_chrM` chrM fragments.
#' Clean fragments are `n_clean`, falling back to `n_pass` and then `n_frag`. A QC column that is
#' missing from the cell table is not applied; an NA value fails the cell (`qc_pass == TRUE` is required),
#' except NA `frac_chrM`, which passes.
#' @param scn scn object (uses `scn$cells`) or a cell table
#' @param min_frag minimal clean fragments
#' @param min_mapq60 minimal fraction of MAPQ-60 fragments
#' @param max_chrM maximal fraction of chrM fragments
#' @param params [caller_defaults()] (depth strata)
#' @return data.table `cell, qc_pass, stratum, n_frag`
#' @export
cell_qc = function(scn, min_frag = 5000, min_mapq60 = 0.5, max_chrM = 0.2, params = caller_defaults()) {
  cells = if (inherits(scn, "scn")) scn$cells else as.data.table(scn)
  nm = names(cells)
  ncol_ = intersect(c("n_clean", "n_pass", "n_frag"), nm)
  if (!length(ncol_)) stop("cell_qc: cells need n_clean, n_pass or n_frag")
  n = as.numeric(cells[[ncol_[1]]])
  pass = n >= min_frag
  if ("frac_mapq60" %in% nm) pass = pass & cells$frac_mapq60 >= min_mapq60
  if ("frac_chrM" %in% nm) pass = pass & (is.na(cells$frac_chrM) | cells$frac_chrM <= max_chrM)
  pass[is.na(pass)] = FALSE
  data.table(cell = as.character(cells$cell), qc_pass = pass,
             stratum = as.character(depth_stratum(n, params$strata_breaks, params$strata_labels)), n_frag = n)
}

# ---- arms and regions --------------------------------------------------------------------------------

#' Chromosome arms (Battenberg coordinates)
#'
#' p = `[0, cen_left)`, q = `[cen_right, chromosome end)` from the Battenberg table
#' (`chromosome_coordinates_hg38_chr.txt`: chr, start, cen.left.base, cen.right.base, end). Battenberg names chrX `chr23` (and chrY `chr24`); they are renamed. When the
#' table has no chrX row, the chrX centromere is `params$chrx_centromere`. Acrocentric p arms
#' (`params$acrocentric_p`) are dropped; chrY arms are kept only for males.
#' @param centromeres path or table of the Battenberg coordinates
#' @param chrom_sizes named chromosome sizes (arms are made for these chromosomes only)
#' @param sex "female" or "male" (case-insensitive, or "F" / "M"; any other value stops)
#' @param params [caller_defaults()]
#' @return data.table `region, chrom, arm, start, end`
#' @export
arm_table = function(centromeres, chrom_sizes, sex = "female", params = caller_defaults()) {
  sex = check_sex(sex, "arm_table")
  ce = if (is.character(centromeres)) fread(centromeres) else copy(as.data.table(centromeres))
  setnames(ce, 1:5, c("chrom", "cstart", "cen_left", "cen_right", "cend"))
  ce[, chrom := as.character(chrom)]
  ce[chrom %in% c("chr23", "23"), chrom := "chrX"]
  ce[chrom %in% c("chr24", "24"), chrom := "chrY"]
  ce = ce[chrom %in% names(chrom_sizes)]
  if (!"chrX" %in% ce$chrom && "chrX" %in% names(chrom_sizes)) {
    ce = rbind(ce[, .(chrom, cstart, cen_left, cen_right, cend)],
               data.table(chrom = "chrX", cstart = 1, cen_left = params$chrx_centromere[1],
                          cen_right = params$chrx_centromere[2], cend = chrom_sizes[["chrX"]]), fill = TRUE)
  }
  if (!identical(sex, "male")) ce = ce[chrom != "chrY"]
  ce[, csz := as.numeric(chrom_sizes[chrom])]
  arms = rbind(ce[, .(region = paste0(chrom, "p"), chrom, arm = "p", start = 0, end = as.numeric(cen_left))],
               ce[, .(region = paste0(chrom, "q"), chrom, arm = "q", start = as.numeric(cen_right), end = csz)])
  arms = arms[!region %in% params$acrocentric_p & end > start]
  arms[order(chrom_rank(chrom), arm)]
}

# events: data.table (region, chrom, start, end), cli_events() output (event = name), or a named
# character vector of region strings (the sample config `events` map; 0-based half-open, an open end = chromosome
# end). An event without a name is named after its region string ("chr20:28000000-"), a bare chromosome "chr7"
# "event|chr7" (the chromosome region is already called "chr7").
event_table = function(events, chrom_sizes) {
  if (is.null(events)) return(NULL)
  if (is.character(events) || is.list(events) && !is.data.frame(events)) {
    ev = parse_region(unlist(events), chrom_sizes)
    nm = names(events)
    if (!is.null(nm)) {
      miss = is.na(nm) | !nzchar(trimws(nm))
      ev[!miss, region := nm[!miss]]
    }
    if (is.null(nm) || any(miss)) {
      un = if (is.null(nm)) rep(TRUE, nrow(ev)) else miss
      ev[un & !grepl(":", region, fixed = TRUE), region := paste0("event|", region)]
    }
  } else {
    ev = copy(as.data.table(events))
    if ("event" %in% names(ev)) ev[, region := event]
    if (!all(c("region", "chrom", "start", "end") %in% names(ev))) stop("events need region, chrom, start, end")
  }
  ev[, start := as.numeric(start)][, end := as.numeric(end)]
  ev[is.na(end), end := as.numeric(chrom_sizes[chrom])]
  ev[, .(region = as.character(region), chrom = as.character(chrom), start, end)]
}

#' sd_class of a region (which null regions give its per-cell null SD)
#'
#' Arm, chrom and tile regions use their own class. Any other region (event, segment) uses, in order:
#' `chrom` if it contains >= 90 % of its chromosome's bins; `arm` if it contains >= 60 % of the bins of
#' one arm; else the tile size nearest to its length on a log scale. Bins are counted under the depth
#' (containment) rule over all bins (see [region_set()]).
#' @param regions table with `chrom, start, end` and optionally `class`
#' @param bins bins table
#' @param arms [arm_table()] output
#' @param params [caller_defaults()] (`tile_sizes`)
#' @return character vector
#' @export
sd_class_of = function(regions, bins, arms, params = caller_defaults()) {
  r = as.data.table(regions)
  ts = params$tile_sizes
  out = if ("class" %in% names(r)) as.character(r$class) else rep("event", nrow(r))
  todo = which(!out %in% c("arm", "chrom", names(ts)))
  if (!length(todo)) return(out)
  bins = as.data.table(bins); arms = as.data.table(arms)
  for (i in todo) {
    ch = r$chrom[i]; e = r$end[i]; if (is.na(e)) e = Inf
    bc = bins$chrom == ch
    inr = bc & bins$start >= r$start[i] & bins$end <= e
    if (any(bc) && sum(inr) >= 0.9 * sum(bc)) { out[i] = "chrom"; next }
    a = arms[chrom == ch]
    fa = vapply(seq_len(nrow(a)), function(k) {
      ina = bc & bins$start >= a$start[k] & bins$end <= a$end[k]
      if (!any(ina)) 0 else sum(inr & ina) / sum(ina)
    }, numeric(1))
    if (length(fa) && max(fa) >= 0.6) { out[i] = "arm"; next }
    len = (if (any(bc)) min(e, max(bins$end[bc])) else e) - r$start[i]
    out[i] = names(ts)[which.min(abs(log(len) - log(ts)))]
  }
  out
}

#' Region x bin indicator under the depth (containment) rule
#'
#' Row r marks the bins of region r's chromosome with `bin.start >= start & bin.end <= end`
#' (NA end = chromosome end), restricted to `use`.
#' @param bins bins table
#' @param regions regions table (`region, chrom, start, end`)
#' @param use logical bin filter (NULL = all bins)
#' @return sparse regions x bins dgCMatrix (row names = region names)
#' @export
region_bin_matrix = function(bins, regions, use = NULL) {
  ok = if (is.null(use)) rep(TRUE, nrow(bins)) else as.logical(use)
  by_chr = split(seq_len(nrow(bins)), bins$chrom)
  ii = vector("list", nrow(regions)); jj = vector("list", nrow(regions))
  bs = bins$start; be = bins$end
  for (r in seq_len(nrow(regions))) {
    rows = by_chr[[regions$chrom[r]]]
    if (is.null(rows)) next
    e = regions$end[r]; if (is.na(e)) e = Inf
    b = rows[bs[rows] >= regions$start[r] & be[rows] <= e & ok[rows]]
    ii[[r]] = rep.int(r, length(b)); jj[[r]] = b
  }
  i = as.integer(unlist(ii)); j = as.integer(unlist(jj))
  sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = c(nrow(regions), nrow(bins)),
               dimnames = list(regions$region, NULL))
}

#' The regions table of the caller
#'
#' Arms ([arm_table()]), whole chromosomes (first to last bin), tiles of `params$tile_sizes`
#' (one per [window_ids()] window, named by [tile_name()], `end = min(start + size, chromosome end)`)
#' and named events, on the chromosomes that have usable bins (chrY only for males). Regions without a
#' usable bin are dropped. Columns follow `vignette("scphasedcn-outputs", package = "scphasedcn")`: `n_bins` / `n_use` count all / usable
#' bins under the containment rule; tiles are `eligible` when `n_bins >= 0.9 * size / bin width` and
#' `n_use >= 0.7 * n_bins` (the eligibility rule of the spike-in calibration), other regions when they
#' have a usable bin; `null_region` = arm, chrom and tile regions on autosomes not in `null_exclude`
#' (chromosome or region names); `allelic_ok` is FALSE on male chrX / chrY; `sd_class` by [sd_class_of()].
#' @param bins bins table (chrom, start, end)
#' @param centromeres Battenberg centromere table or path ([arm_table()])
#' @param classes region classes to build (`arm`, `chrom`, tile labels); `event` is added when `events` is given
#' @param events NULL, a data.table `region, chrom, start, end` (0-based half-open; NA end = chromosome end),
#'   [cli_events()] output or a named vector of region strings (`"chr20:28000000-"`, `"chr7"`; [parse_region()]). An
#'   event without a name is named after its region string (a bare chromosome `"chr7"`: `"event|chr7"`). An event
#'   that yields no region (a chromosome without usable bins, or no whole usable bin inside) is dropped with a
#'   warning.
#' @param sex "female" or "male" (case-insensitive, or "F" / "M"; any other value stops)
#' @param null_exclude chromosomes (or region names) whose regions never contribute null rows
#' @param use logical bin filter (NULL = all bins)
#' @param params [caller_defaults()]
#' @return regions data.table (passes [check_regions()])
#' @export
region_set = function(bins, centromeres, classes = c("arm", "chrom", "5Mb", "10Mb", "20Mb", "50Mb"), events = NULL,
                      sex = "female", null_exclude = character(), use = NULL, params = caller_defaults()) {
  sex = check_sex(sex, "region_set")
  bins = as.data.table(bins)
  use = if (is.null(use)) rep(TRUE, nrow(bins)) else as.logical(use)
  bad = setdiff(classes, c("arm", "chrom", names(params$tile_sizes), "event"))
  if (length(bad)) stop("region_set: unknown class(es) ", paste(bad, collapse = ", "))
  cs = chrom_sizes_from_bins(bins)
  chroms = unique(bins$chrom[use])
  if (!identical(sex, "male")) chroms = setdiff(chroms, "chrY")
  arms = arm_table(centromeres, cs, sex, params)
  out = list()
  if ("arm" %in% classes) out$arm = arms[chrom %in% chroms, .(region, class = "arm", chrom, start, end)]
  if ("chrom" %in% classes)
    out$chrom = bins[chrom %in% chroms, .(start = min(start), end = max(end)), by = chrom][, .(region = chrom, class = "chrom", chrom, start, end)]
  b = bins[chrom %in% chroms]
  for (lab in intersect(names(params$tile_sizes), classes)) {
    size = params$tile_sizes[[lab]]
    ids = unique(window_ids(b, size))
    t0 = data.table(chrom = sub(":[^:]*$", "", ids), start = as.numeric(sub("^.*:", "", ids)))
    t0[, end := pmin(start + size, as.numeric(cs[chrom]))]
    out[[lab]] = t0[, .(region = tile_name(lab, chrom, start), class = lab, chrom, start, end)]
  }
  ev = event_table(events, cs)
  if (!is.null(ev) && nrow(ev)) out$event = ev[chrom %in% chroms, .(region, class = "event", chrom, start, end)]
  r = rbindlist(out)
  M = region_bin_matrix(bins, r)
  r[, n_bins := as.integer(round(Matrix::rowSums(M)))]
  r[, n_use := as.integer(round(as.vector(M %*% as.numeric(use))))]
  if (!is.null(ev) && nrow(ev)) {
    # events that yield no region: a chromosome without usable bins (unknown, chrY of a female, not in `use`), or no
    # whole usable bin inside the event (containment rule)
    no_chr = ev[!chrom %in% chroms, region]
    no_bin = r[class == "event" & n_use == 0, region]
    if (length(no_chr) || length(no_bin))
      warning("region_set: event(s) dropped without a region: ",
              paste(c(if (length(no_chr)) sprintf("%s (chromosome %s not in the usable bins)", no_chr, ev[!chrom %in% chroms, chrom]),
                      if (length(no_bin)) sprintf("%s (no whole usable bin inside)", no_bin)), collapse = ", "), call. = FALSE)
  }
  r = r[n_use > 0]
  bw = stats::median(bins$end - bins$start)
  tsz = params$tile_sizes[r$class]
  r[, eligible := ifelse(class %in% names(params$tile_sizes), n_bins >= 0.9 * tsz / bw & n_use >= 0.7 * n_bins, n_use > 0)]
  r[, null_region := class %in% c("arm", "chrom", names(params$tile_sizes)) & chrom %in% AUTOSOMES_CALLER &
        !chrom %in% null_exclude & !region %in% null_exclude]
  r[, allelic_ok := !(identical(sex, "male") & chrom %in% c("chrX", "chrY"))]
  r[, sd_class := sd_class_of(r, bins, arms, params)]
  ord = order(match(r$class, region_classes()), chrom_rank(r$chrom), r$start, r$end)
  r = r[ord, .(region, class, sd_class, chrom, start, end, n_bins, n_use, eligible, null_region, allelic_ok)]
  if (anyDuplicated(r$region)) stop("region_set: duplicated region names: ", paste(head(r$region[duplicated(r$region)]), collapse = ", "))
  check_regions(r)
  r[]
}

# ---- depth ------------------------------------------------------------------------------------------------

mad_sd = function(m) apply(m, 2, function(v) 1.4826 * stats::mad(v, center = stats::median(v, na.rm = TRUE), constant = 1, na.rm = TRUE))

# null reference regions of an sd_class that the regions table does not contain (e.g. an event with
# sd_class 10Mb when no 10Mb tiles were requested): built from the bins on the null chromosomes
sd_reference_regions = function(bins, use, cls, null_chroms, params) {
  b = bins[use & bins$chrom %in% null_chroms]
  if (cls == "chrom") return(b[, .(region = chrom, chrom = chrom[1], start = min(start), end = max(end)), by = chrom][, .(region, chrom, start, end)])
  if (cls == "arm") {
    if (!"arm" %in% names(b)) return(NULL)
    return(b[arm %in% c("p", "q"), .(start = min(start), end = max(end)), by = .(chrom, arm)][, .(region = paste0(chrom, arm), chrom, start, end)])
  }
  if (cls %in% names(params$tile_sizes)) {
    size = params$tile_sizes[[cls]]
    ids = unique(window_ids(b, size))
    t0 = data.table(region = ids, chrom = sub(":[^:]*$", "", ids), start = as.numeric(sub("^.*:", "", ids)))
    t0[, end := start + size]
    return(t0)
  }
  NULL
}

# regions x cells matrices x, e, log2ratio, sd_null, z
depth_mats = function(scn, norm, regions, params = caller_defaults()) {
  regions = as.data.table(regions)
  use = norm$use
  M = region_bin_matrix(scn$bins, regions, use)
  E = norm$expected; E[is.na(E)] = 0
  x = as.matrix(M %*% scn$counts); e = as.matrix(M %*% E)
  dimnames(x) = dimnames(e) = NULL
  lr = log2((x + 0.5) / (e + 0.5))
  sdm = matrix(NA_real_, nrow(regions), ncol(x))
  null_chroms = unique(regions$chrom[regions$null_region])
  if (!length(null_chroms)) null_chroms = AUTOSOMES_CALLER
  for (k in unique(regions$sd_class)) {
    ref = which(regions$class == k & regions$null_region & regions$n_use > 0)
    s = if (length(ref)) mad_sd(lr[ref, , drop = FALSE]) else {
      rr = sd_reference_regions(as.data.table(scn$bins), use, k, null_chroms, params)
      if (is.null(rr) || !nrow(rr)) rep(NA_real_, ncol(x)) else {
        Mr = region_bin_matrix(scn$bins, rr, use)
        Mr = Mr[Matrix::rowSums(Mr) > 0, , drop = FALSE]
        mad_sd(log2((as.matrix(Mr %*% scn$counts) + 0.5) / (as.matrix(Mr %*% E) + 0.5)))
      }
    }
    rows = which(regions$sd_class == k)
    sdm[rows, ] = rep(s, each = length(rows))
  }
  z = lr / pmax(sdm, params$sd_floor)
  list(x = x, e = e, log2ratio = lr, sd_null = sdm, z = z)
}

#' Depth features per cell and region
#'
#' Observed (`x`) and expected (`e`) fragments of each region's usable bins (containment rule) by sparse
#' region x bin products, `log2ratio = log2((x + 0.5) / (e + 0.5))`,
#' `sd_null` = 1.4826 x MAD (`mad(constant = 1)`) of the cell's `log2ratio` over the null regions
#' (`null_region`, with a usable bin) whose class is the row's `sd_class` (tiles: all null tiles of that
#' size), and `z = log2ratio / pmax(sd_null, sd_floor)`. When the regions table
#' has no null region of a needed sd_class, the reference regions are built from the bins of the null
#' chromosomes (arms from `bins$arm`, chromosomes, or tiles).
#' @param scn scn object (counts, bins, cells)
#' @param norm `scn_norm` from [normalise_depth()] (expected counts, `use`)
#' @param regions regions table ([region_set()])
#' @param params [caller_defaults()] (`sd_floor`)
#' @return data.table keyed by `cell, region`: cell, region, class, sd_class, chrom, x, e, log2ratio, sd_null, z
#' @export
depth_features = function(scn, norm, regions, params = caller_defaults()) {
  regions = as.data.table(regions)
  dm = depth_mats(scn, norm, regions, params)
  nr = nrow(regions); cells = scn$cells$cell; nc = length(cells)
  d = data.table(cell = rep(cells, each = nr), region = rep(regions$region, nc), class = rep(regions$class, nc),
                 sd_class = rep(regions$sd_class, nc), chrom = rep(regions$chrom, nc),
                 x = as.vector(dm$x), e = as.vector(dm$e), log2ratio = as.vector(dm$log2ratio),
                 sd_null = as.vector(dm$sd_null), z = as.vector(dm$z))
  setkeyv(d, c("cell", "region"))
  d[]
}

#' Window count matrices for clone pseudobulks
#'
#' Observed and expected fragments summed over the usable bins of [window_ids()] windows (default 1 Mb).
#' @param scn scn object
#' @param norm `scn_norm` ([normalise_depth()])
#' @param size window size (bp)
#' @return list `x`, `e` (windows x cells matrices, row names = window ids, column names = cells) and
#'   `windows` (window, chrom, start, end, n_use)
#' @export
window_matrices = function(scn, norm, size = 1e6) {
  ids = window_ids(scn$bins, size)
  ids[!norm$use] = NA
  ind = group_indicator(ids)
  E = norm$expected; E[is.na(E)] = 0
  x = as.matrix(ind %*% scn$counts); e = as.matrix(ind %*% E)
  dimnames(x) = dimnames(e) = list(rownames(ind), scn$cells$cell)
  wt = window_table(scn$bins[!is.na(ids)], ids[!is.na(ids)])[match(rownames(ind), window)]
  setnames(wt, "n_bins", "n_use")
  list(x = x, e = e, windows = wt[])
}

#' Per-cell depth features on arbitrary regions
#'
#' Sums the 1 Mb window counts (`win`) of the windows inside each region (window start >= region start
#' and window end <= region end) and returns the feature columns of the long features table: x, e, log2ratio, sd_null (the
#' cell's `sd_null` of the region's `sd_class` from the run's features; the nearest available class when
#' missing), z, stratum, n_frag. For regions that are not in the feature table (clone pieces and segments,
#' MEDICC2 pieces). No allelic columns: the callers score those with a clone-specific orientation.
#' @param feat scn_features list (or `clone_input`) with the 1 Mb window depth `win`
#' @param regs regions table (region, chrom, start, end, class, sd_class)
#' @param cells cells to compute (default all)
#' @param sd_floor floor of sd_null
#' @return features long table: cell, region, x, e, class, sd_class, chrom, log2ratio, sd_null, z, n_frag, stratum
#' @export
region_features = function(feat, regs, cells = NULL, sd_floor = caller_defaults()$sd_floor) {
  inp = clone_input(feat)
  win = inp$win
  if (is.null(win)) stop("region_features: input lacks the 1 Mb window depth `win`")
  rg = as.data.table(regs)
  cc = cells %||% colnames(win$x)
  cc = intersect(cc, colnames(win$x))
  wt = win$windows
  hit = wt[, .(wi = .I, chrom, start, end)][rg[, .(ri = .I, chrom, rs = start, re = end)],
                                             on = .(chrom, start >= rs, end <= re), .(wi = x.wi, ri = i.ri), nomatch = NULL, allow.cartesian = TRUE]
  W = Matrix::sparseMatrix(i = hit$ri, j = hit$wi, x = 1, dims = c(nrow(rg), nrow(wt)))
  xr = as.matrix(W %*% win$x[, cc, drop = FALSE]); er = as.matrix(W %*% win$e[, cc, drop = FALSE])
  d = data.table(cell = rep(cc, each = nrow(rg)), region = rep(rg$region, length(cc)), x = as.vector(xr), e = as.vector(er))
  d[rg, `:=`(class = i.class, sd_class = i.sd_class, chrom = i.chrom), on = "region"]
  d[, log2ratio := log2((x + 0.5) / (e + 0.5))]
  sn = unique(inp$features[, .(cell, sd_class, sd_null)], by = c("cell", "sd_class"))
  avail = unique(sn$sd_class)
  order_cls = c("chrom", "arm", "50Mb", "20Mb", "10Mb", "5Mb")
  d[, sdc := sd_class]
  miss = setdiff(unique(d$sdc), avail)
  for (m in miss) {
    pos = match(m, order_cls)
    alt = avail[order(abs(match(avail, order_cls) - pos))][1]
    d[sdc == m, sdc := alt]
  }
  d[sn, sd_null := i.sd_null, on = .(cell, sdc = sd_class)]
  d[, z := log2ratio / pmax(sd_null, sd_floor)]
  cm = inp$cells
  if (!"n_frag" %in% names(cm)) cm = unique(inp$features[, .(cell, n_frag, stratum)], by = "cell")
  if (!"stratum" %in% names(cm)) cm[, stratum := as.character(depth_stratum(n_frag))]
  d[cm, `:=`(n_frag = i.n_frag, stratum = i.stratum), on = "cell"]
  d[, sdc := NULL]
  d[]
}
utils::globalVariables(c("wi", "ri", "rs", "re", "x.wi", "i.ri", "i.class", "i.sd_class", "i.chrom", "sdc", "i.sd_null",
                         "i.n_frag", "i.stratum"))

# ---- molecules and votes -------------------------------------------------------------------------------

#' Collapse fragments into molecules
#'
#' Fragments with haplotype A or B, sorted by cell,
#' chromosome and start; a new molecule starts when `start > cummax(end) + gap` of the previous
#' fragments of the same cell and chromosome; one vote per molecule for the majority haplotype
#' (ties dropped); `ps` is the first fragment's phase set. The running maximum is computed per group
#' with an integer offset (exact for coordinates below 2^53 / number of groups).
#'
#' Phase sets are not oriented relative to each other, so a molecule whose fragments lie in different phase sets
#' (a molecule across a block junction) mixes unrelated haplotype labels. With `split_ps = TRUE` the molecules are
#' formed as above and then split by phase set: the fragments of each phase set of a molecule give their own vote
#' (majority haplotype of those fragments, ties dropped) in that phase set, with `start` / `end` / `n_frag` of
#' those fragments. Molecules within one phase set are the same under both settings.
#' @param fr fragments with `cell, chrom, start, end, hap, ps`; `hap` "A"/"B" (other labels dropped) or
#'   integer codes 1 / 2 (= A / B)
#' @param gap maximal gap (bp) between fragments of one molecule
#' @param split_ps one vote per molecule x phase set instead of one per molecule (default FALSE; the caller
#'   parameter `mol_split_ps`)
#' @return data.table `cell, chrom, start, end, ps, n_frag, hap` (hap "A"/"B"), one row per molecule (per molecule x
#'   phase set with `split_ps`)
#' @export
collapse_molecules = function(fr, gap = 2000, split_ps = FALSE) {
  h = fr$hap
  isa = if (is.character(h)) h == "A" else h == 1L
  isb = if (is.character(h)) h == "B" else h == 2L
  keep = which(isa | isb)
  d = data.table(cell = fr$cell[keep], chrom = fr$chrom[keep], start = fr$start[keep], end = fr$end[keep],
                 isA = isa[keep], ps = fr$ps[keep])
  if (!nrow(d)) return(data.table(cell = d$cell, chrom = character(), start = numeric(), end = numeric(), ps = d$ps,
                                  n_frag = integer(), hap = character()))
  setorderv(d, c("cell", "chrom", "start"))
  n = nrow(d)
  first = c(TRUE, d$cell[-1] != d$cell[-n] | d$chrom[-1] != d$chrom[-n])
  g = cumsum(first)
  big = max(as.numeric(d$end)) + 1
  off = (g - 1) * big
  cm = cummax(as.numeric(d$end) + off) - off
  d[, mol := cumsum(first | c(TRUE, d$start[-1] > cm[-n] + gap))]
  d[, isB := !isA]
  if (isTRUE(split_ps)) {
    m = d[, .(cell = cell[1], chrom = chrom[1], start = min(start), end = max(end),
              nA = sum(isA), nB = sum(isB), n_frag = .N), by = .(mol, ps)]
    setcolorder(m, c("mol", "cell", "chrom", "start", "end", "ps"))
  } else {
    m = d[, .(cell = cell[1], chrom = chrom[1], start = min(start), end = max(end), ps = ps[1],
              nA = sum(isA), nB = sum(isB), n_frag = .N), by = mol]
  }
  m = m[nA != nB]
  m[, hap := ifelse(nA > nB, "A", "B")]
  m[, c("nA", "nB", "mol") := NULL]
  m[]
}

#' Allele atoms
#'
#' The smallest genomic units votes are counted in: the bins of the regions' chromosomes, split at
#' every region start / end that is not a bin boundary, so that "molecule start inside a region" is
#' "atom inside the region" for every region.
#' @param bins bins table
#' @param regions regions table (`chrom, start, end`)
#' @return data.table `atom, chrom, start, end, bin` (bin = row of `bins`; passes [check_atoms()])
#' @export
atoms_from_regions = function(bins, regions) {
  bins = as.data.table(bins); regions = as.data.table(regions)
  chroms = intersect(unique(bins$chrom), unique(regions$chrom))
  out = lapply(chroms, function(ch) {
    rows = which(bins$chrom == ch)
    bs = bins$start[rows]; lo = min(bs); hi = max(bins$end[rows])
    rr = regions[chrom == ch]
    cut = c(bs, rr$start, rr$end)
    cut = sort(unique(cut[!is.na(cut) & cut >= lo & cut < hi]))
    data.table(chrom = ch, start = cut, end = c(cut[-1], hi), bin = rows[findInterval(cut, bs)])
  })
  a = rbindlist(out)
  a[, atom := seq_len(.N)]
  setcolorder(a, c("atom", "chrom", "start", "end", "bin"))
  a[, bin := as.integer(bin)]
  a[]
}

# blacklist bit of the fragments' `flags` column (F_BLACKLIST in the extractor's bam.py). The extractor's depth
# layers count the "clean" fragments, `(flags & F_BLACKLIST) == 0` (cli.py, cmd_extract), i.e. those whose start is
# outside the blacklist; `votes_clean` applies the same definition to the votes.
FRAG_FLAG_BLACKLIST = 8L

# hap-assigned fragments with a phase set of one chromosome: start, end, cell, hap (1/2), ps; `clean = TRUE` keeps only
# clean fragments (blacklist bit unset). The tabix pipeline runs under `set -o pipefail` with its output in a temporary
# file: no record on the contig (exit 0, empty output) gives no rows, a failing command (missing / corrupt file or
# index, tabix error) stops with the command's error
read_hap_fragments = function(frag_file, chrom, tabix = NULL, clean = FALSE) {
  empty = data.table(start = numeric(), end = numeric(), cell = character(), hap = integer(), ps = numeric())
  if (!is.null(tabix) && nzchar(tabix) && file.exists(tabix) && file.exists(paste0(frag_file, ".tbi"))) {
    flt = if (isTRUE(clean)) sprintf(" && int($7/%d)%%2==0", FRAG_FLAG_BLACKLIST) else ""
    pipe = sprintf("set -o pipefail; %s %s %s | awk -F'\\t' '($11==1||$11==2) && $13>0%s' | cut -f2,3,4,7,11,13",
                   shQuote(tabix), shQuote(frag_file), shQuote(chrom), flt)
    out = tempfile(fileext = ".tsv"); err = tempfile(fileext = ".err")
    on.exit(unlink(c(out, err)), add = TRUE)
    status = suppressWarnings(system2("bash", c("-c", shQuote(pipe)), stdout = out, stderr = err))
    if (!identical(as.integer(status), 0L)) {
      msg = if (file.exists(err)) paste(utils::head(readLines(err, warn = FALSE), 5), collapse = "; ") else ""
      stop(sprintf("read_hap_fragments: reading %s (%s) with tabix failed (exit status %s)%s", frag_file, chrom, status,
                   if (nzchar(msg)) paste0(": ", msg) else ""), call. = FALSE)
    }
    if (!file.exists(out) || file.size(out) == 0) return(empty)
    fr = fread(out, header = FALSE, sep = "\t", col.names = c("start", "end", "cell", "flags", "hap", "ps"),
               colClasses = list(character = 3))
    if (ncol(fr) != 6) stop(sprintf("read_hap_fragments: %s (%s) does not have the fragment table's columns", frag_file, chrom),
                            call. = FALSE)
    fr[, flags := NULL]
    return(fr)
  }
  fr = read_fragments(frag_file, regions = data.table(chrom = chrom, start = 0, end = NA_real_))
  if (!nrow(fr) || !all(c("hap", "ps") %in% names(fr))) return(empty)
  chr0 = chrom
  fr = fr[fr$chrom == chr0]
  fr[, hap := hap_code(hap)]
  if (isTRUE(clean)) fr = fr[bitwAnd(as.integer(flags), FRAG_FLAG_BLACKLIST) == 0L]
  fr[hap %in% 1:2 & ps > 0, .(start, end, cell = as.character(cell), hap, ps)]
}

#' Molecule votes per cell, atom and phase set
#'
#' Per chromosome: `tabix <file> <chrom> | cut -f2,3,4,7,11,13` read with `fread(cmd = )` (an awk filter keeps hap 1 / 2 with `ps > 0` before R sees the rows), or
#' [read_fragments()] when the tabix binary or index is missing; fragments of `cells` with haplotype
#' A / B and a phase set are collapsed into molecules ([collapse_molecules()]), each molecule is assigned
#' to the atom containing its start and the votes are aggregated per cell x atom x phase set. Memory
#' stays per chromosome (integer cell indices, `gc()` between chromosomes). A failing tabix command stops with its
#' error; a chromosome with atoms but no haplotype-assigned fragment of `cells` gives a warning.
#'
#' `clean = TRUE` keeps only the fragments the extractor counts in the depth layers (no blacklist bit in `flags`,
#' the extractor's `clean = (flags & F_BLACKLIST) == 0`), in the tabix path (awk) and in the fallback;
#' `split_ps = TRUE` gives one vote per molecule x phase set ([collapse_molecules()]).
#' @param frag_file fragments.tsv.gz (+ .tbi)
#' @param cells cell barcodes; the votes' `cell` is the index into this vector
#' @param atoms [atoms_from_regions()] output
#' @param chroms chromosomes to read (default: the atoms' chromosomes)
#' @param gap molecule gap (bp)
#' @param tabix path of the tabix binary ("" = the `tabix` on the PATH; NULL / not found = the R fallback through
#'   Rsamtools)
#' @param split_ps split molecules by phase set (caller parameter `mol_split_ps`; default FALSE)
#' @param clean votes from clean (not blacklisted) fragments only (caller parameter `votes_clean`; default FALSE)
#' @return votes data.table `cell, atom, ps, nA, nB` (passes [check_votes()])
#' @export
read_votes = function(frag_file, cells, atoms, chroms = NULL, gap = 2000, tabix = caller_defaults()$tabix,
                      split_ps = FALSE, clean = FALSE) {
  atoms = as.data.table(atoms)
  chroms = chroms %||% unique(atoms$chrom)
  if (identical(tabix, "")) tabix = unname(Sys.which("tabix"))
  out = vector("list", length(chroms))
  none = character()                                 # chromosomes with atoms but no haplotype fragment of the cells
  for (k in seq_along(chroms)) {
    ch = chroms[k]
    a = atoms[chrom == ch][order(start)]
    if (!nrow(a)) next
    fr = read_hap_fragments(frag_file, ch, tabix, clean = clean)
    if (!nrow(fr)) { none = c(none, ch); next }
    ci = match(fr$cell, cells)
    fr = data.table(cell = ci, chrom = ch, start = as.numeric(fr$start), end = as.numeric(fr$end), hap = as.integer(fr$hap), ps = fr$ps)[!is.na(cell)]
    if (!nrow(fr)) { none = c(none, ch); next }
    big = nrow(fr) > 5e6
    m = collapse_molecules(fr, gap = gap, split_ps = split_ps)
    rm(fr)
    k_at = findInterval(m$start, a$start)
    ok = k_at > 0 & m$start < a$end[pmax(k_at, 1)]
    m = m[ok][, atom := a$atom[k_at[ok]]]
    m[, isA := hap == "A"][, isB := !isA]
    out[[k]] = m[, .(nA = as.integer(sum(isA)), nB = as.integer(sum(isB))), by = .(cell, atom, ps)]
    rm(m)
    if (big) gc(verbose = FALSE)
  }
  if (length(none))
    warning(sprintf("read_votes: no haplotype-assigned fragment with a phase set of the given cells on %d chromosome(s) (%s): their allele features are empty (unphased extraction, or barcodes that do not match?)",
                    length(none), paste(utils::head(none, 10), collapse = ", ")), call. = FALSE)
  v = rbindlist(out)
  if (!nrow(v)) v = data.table(cell = integer(), atom = integer(), ps = numeric(), nA = integer(), nB = integer())
  v[, cell := as.integer(cell)][, atom := as.integer(atom)][, ps := as.numeric(ps)]
  setorderv(v, c("cell", "atom", "ps"))
  v[]
}

#' Phasing regime of a votes table
#'
#' `"chromosome"` (one phase set per chromosome, e.g. a hybrid longphase + population phase) when the
#' median over chromosomes of the number of phase sets needed to hold `frac` of the chromosome's votes is
#' at most 2, else `"block"` (unoriented blocks). Vote-weighted so that the residual small phase sets of a
#' hybrid phase (many small phase sets next to one that holds nearly all votes) do not count; without residual sets this is the "median number of distinct phase sets <= 2" rule.
#' @param votes votes table
#' @param atoms atoms table
#' @param frac vote fraction the counted phase sets must hold
#' @return "chromosome" or "block"
#' @export
detect_phasing = function(votes, atoms, frac = 0.9) {
  if (!nrow(votes)) return("chromosome")
  ch = atoms$chrom[match(votes$atom, atoms$atom)]
  t = data.table(chrom = ch, ps = votes$ps, n = votes$nA + votes$nB)[, .(n = sum(n)), by = .(chrom, ps)]
  setorderv(t, c("chrom", "n"), c(1L, -1L))
  k = t[, .(k = which(cumsum(n) >= frac * sum(n))[1]), by = chrom]$k
  if (stats::median(k) <= 2) "chromosome" else "block"
}

# ---- region votes ----------------------------------------------------------------------------------------

# greedy interval colouring: regions in one layer never overlap (each atom maps to <= 1 region per layer)
region_layers = function(regions) {
  o = order(regions$chrom, regions$start, regions$end)
  layer = integer(nrow(regions))
  last = list()                                    # per layer: named vector chrom -> last end
  for (r in o) {
    ch = regions$chrom[r]; s = regions$start[r]; e = regions$end[r]; if (is.na(e)) e = Inf
    l = 0L
    for (k in seq_along(last)) {
      le = last[[k]][ch]
      if (is.na(le) || le <= s) { l = k; break }
    }
    if (l == 0L) { last[[length(last) + 1]] = numeric(); l = length(last) }
    last[[l]][ch] = e
    layer[r] = l
  }
  layer
}

# atom -> region index (start-inside rule) with the region's layer
atom_region_map = function(atoms, regions) {
  at = as.data.table(atoms)
  by_chr = split(seq_len(nrow(at)), at$chrom)
  lay = region_layers(regions)
  m = rbindlist(lapply(seq_len(nrow(regions)), function(r) {
    idx = by_chr[[regions$chrom[r]]]
    if (is.null(idx)) return(NULL)
    e = regions$end[r]; if (is.na(e)) e = Inf
    k = idx[at$start[idx] >= regions$start[r] & at$start[idx] < e]
    if (!length(k)) return(NULL)
    data.table(atom = at$atom[k], ri = r)
  }))
  if (!nrow(m)) m = data.table(atom = integer(), ri = integer())
  m[, layer := lay[ri]]
  m[]
}

# per layer: votes summed per cell x region x phase set
layer_blocks = function(votes, map, l) {
  ml = map[layer == l]
  ri = ml$ri[match(votes$atom, ml$atom)]
  sel = which(!is.na(ri))
  if (!length(sel)) return(data.table(cell = integer(), ri = integer(), ps = numeric(), nA = integer(), nB = integer()))
  b = data.table(cell = votes$cell[sel], ri = ri[sel], ps = votes$ps[sel], nA = votes$nA[sel], nB = votes$nB[sel])
  b[, .(nA = sum(nA), nB = sum(nB)), by = .(cell, ri, ps)]
}

# leave-one-out orientation of one layer's blocks (see orient_blocks()); tied blocks are dropped
orient_layer = function(b, groups, loo) {
  if (!nrow(b)) return(data.table(cell = integer(), ri = integer(), nA = numeric(), nB = numeric(), n_blocks = integer()))
  b[, grp := if (is.null(groups)) NA_character_ else as.character(groups[cell])]
  all_tot = b[, .(gA = sum(nA), gB = sum(nB)), by = .(ri, ps)]
  b[all_tot, `:=`(gA = i.gA, gB = i.gB), on = .(ri, ps)]
  if (any(!is.na(b$grp))) {
    grp_tot = b[!is.na(grp), .(hA = sum(nA), hB = sum(nB)), by = .(grp, ri, ps)]
    b[grp_tot, `:=`(gA = i.hA, gB = i.hB), on = .(grp, ri, ps)]
  }
  if (loo) b[, `:=`(rA = gA - nA, rB = gB - nB)] else b[, `:=`(rA = gA, rB = gB)]
  b = b[rA != rB]
  b[, flip := rA < rB]
  b[, .(nA = sum(ifelse(flip, nB, nA)), nB = sum(ifelse(flip, nA, nB)), n_blocks = .N), by = .(cell, ri)]
}

check_groups = function(groups, votes) {
  if (is.null(groups)) return(NULL)
  if (length(groups) < max(votes$cell, 0L)) stop("groups: need one value per cell (index of the votes' `cell`)")
  g = as.character(groups)
  g[!is.na(g) & !nzchar(g)] = NA
  g
}

#' Orient phase blocks by a leave-one-out pseudobulk
#'
#' Per region and phase set, the pseudobulk votes of the
#' cell's group (all cells for a cell with NA group or when `groups` is NULL) minus the cell's own votes
#' (`loo = TRUE`) decide the block's orientation: tied blocks are dropped, blocks with `rA < rB` are
#' flipped, and the oriented votes are summed per cell and region. Region membership: molecule (atom)
#' start inside `[start, end)`.
#' @param votes votes table ([read_votes()])
#' @param regions regions table
#' @param atoms atoms table
#' @param groups NULL or one group label per cell (indexed by the votes' `cell`; NA = all-cell pseudobulk)
#' @param loo leave the cell's own votes out of the pseudobulk
#' @return data.table `cell, region, nA, nB, n_blocks` (rows with at least one non-tied block)
#' @export
orient_blocks = function(votes, regions, atoms, groups = NULL, loo = TRUE) {
  regions = as.data.table(regions)
  g = check_groups(groups, votes)
  map = atom_region_map(atoms, regions)
  out = lapply(sort(unique(map$layer)), function(l) orient_layer(layer_blocks(votes, map, l), g, loo))
  o = rbindlist(out)
  if (!nrow(o)) return(data.table(cell = integer(), region = character(), nA = numeric(), nB = numeric(), n_blocks = integer()))
  o[, .(cell, region = regions$region[ri], nA, nB, n_blocks)]
}

# mirrored allele statistics (per-block orientation marginalised)
mirrored_stats = function(hc, rho, eps = 0.02, states = imbalance_states(), grid = seq(0.5, 1 - eps, by = 0.01)) {
  mir = function(a, n, f) {
    l1 = bb_loglik(a, n, f, rho, eps); l2 = bb_loglik(a, n, 1 - f, rho, eps)
    log(0.5) + pmax(l1, l2) + log1p(exp(-abs(l1 - l2)))
  }
  tab = hc[n > 0, .(cell, region, nA, n)]
  for (k in seq_len(nrow(states))) tab[, (paste0("ll_", states$state[k])) := mir(nA, n, states$f[k])]
  for (g in seq_along(grid)) tab[, (paste0("g", g)) := mir(nA, n, grid[g])]
  llcols = paste0("ll_", states$state); gcols = paste0("g", seq_along(grid))
  out = tab[, c(list(n = sum(n), n_blocks = .N), lapply(.SD, sum)), by = .(cell, region), .SDcols = c(llcols, gcols)]
  G = as.matrix(out[, gcols, with = FALSE])
  out[, f_mirror := grid[max.col(G, ties.method = "first")]]
  out[, (gcols) := NULL]
  out[]
}

#' Allele features per cell and region
#'
#' Votes summed per cell and region (molecule start inside the region) after orienting the phase
#' blocks: `chromosome` = plain sums (chromosome-wide anchored phase; `orient = "anchored"`); `loo` =
#' [orient_blocks()] (`orient = "pseudobulk"` without groups, `"clone"` with groups); `mirrored` = block
#' orientation marginalised ( `f` = the grid MLE `f_mirror`, `nA` / `nB` = the
#' block-wise majority / minority votes; `orient = "none"`); `auto` = `chromosome` when
#' [detect_phasing()] finds chromosome-wide phase sets, else `loo`. Regions with `allelic_ok = FALSE` are
#' skipped.
#' @param votes votes table
#' @param atoms atoms table
#' @param regions regions table
#' @param orientation "auto", "chromosome", "loo" or "mirrored"
#' @param groups NULL or one group label per cell (see [orient_blocks()])
#' @param rho beta-binomial overdispersion for `mirrored` (NULL: fitted on the blocks of null regions with
#'   at least 2 votes, at most 2e6 seeded rows)
#' @param params [caller_defaults()] (`eps`, `seed`)
#' @param exclude_atoms NULL or atom ids whose votes are dropped before the region sums (an allelic mask of
#'   overdispersed tiles, [allelic_noisy_tiles()]); orientation and block logic are otherwise unchanged
#' @return data.table `cell` (index), `region, nA, nB, n_hap, n_blocks, f, orient`; attribute `rho` for mirrored
#' @export
allelic_features = function(votes, atoms, regions, orientation = c("auto", "chromosome", "loo", "mirrored"), groups = NULL,
                            rho = NULL, params = caller_defaults(), exclude_atoms = NULL) {
  orientation = match.arg(orientation)
  if (length(exclude_atoms)) votes = votes[!atom %in% exclude_atoms]
  regions = as.data.table(regions)
  if ("allelic_ok" %in% names(regions)) regions = regions[allelic_ok == TRUE]
  if (orientation == "auto") orientation = if (detect_phasing(votes, atoms) == "chromosome") "chromosome" else "loo"
  g = check_groups(groups, votes)
  label = switch(orientation, chromosome = "anchored", loo = if (is.null(g) || all(is.na(g))) "pseudobulk" else "clone", mirrored = "none")
  map = atom_region_map(atoms, regions)
  layers = sort(unique(map$layer))
  rho_used = NULL
  if (orientation == "mirrored" && is.null(rho)) {
    nb = rbindlist(lapply(layers, function(l) layer_blocks(votes, map, l)))
    nb = nb[(nA + nB) >= 2]
    if ("null_region" %in% names(regions) && any(regions$null_region[nb$ri])) nb = nb[regions$null_region[ri]]
    set.seed(params$seed); nb = nb[sample.int(.N, min(.N, 2e6))]
    rho = fit_bb_rho(nb$nA, nb$nA + nb$nB)
    rm(nb)
  }
  out = lapply(layers, function(l) {
    b = layer_blocks(votes, map, l)
    if (!nrow(b)) return(NULL)
    switch(orientation,
      chromosome = b[, .(nA = sum(nA), nB = sum(nB), n_blocks = .N), by = .(cell, ri)],
      loo = orient_layer(b, g, TRUE),
      mirrored = {
        b[, nA_maj := pmax(nA, nB)]
        mx = b[, .(nA = sum(nA_maj), nB = sum(nA + nB - nA_maj)), by = .(cell, ri)]
        ms = mirrored_stats(b[, .(cell, region = ri, nA, n = nA + nB)], rho = rho, eps = params$eps)
        mx[ms, `:=`(f = i.f_mirror, n_blocks = i.n_blocks), on = .(cell, ri = region)]
        mx
      })
  })
  o = rbindlist(out, fill = TRUE)
  if (!nrow(o)) o = data.table(cell = integer(), ri = integer(), nA = numeric(), nB = numeric(), n_blocks = integer(), f = numeric())
  o[, n_hap := as.numeric(nA + nB)]
  if (orientation != "mirrored") o[, f := ifelse(n_hap > 0, nA / n_hap, NA_real_)]
  res = o[, .(cell = as.integer(cell), region = regions$region[ri], nA = as.numeric(nA), nB = as.numeric(nB), n_hap,
              n_blocks = as.numeric(n_blocks), f, orient = label)]
  if (orientation == "mirrored") setattr(res, "rho", rho)
  res[]
}

# ---- top level ---------------------------------------------------------------------------------------------

# the sample reader of the caller: bins_100kb.h5ad (requested layers only, no obsm), the fragments path and the
# cells_qc.tsv columns merged into `cells`; no alleles file (the caller reads its allele votes from the fragment table)
read_sample_depth = function(dir, layers = character()) {
  scn = read_bins_h5ad(file.path(dir, "bins_100kb.h5ad"), layers = layers)
  scn$obsm = list()
  ff = file.path(dir, "fragments.tsv.gz")
  if (file.exists(ff)) scn$fragments = ff
  qf = file.path(dir, "cells_qc.tsv")
  if (file.exists(qf)) {
    qc = read_cells_qc(qf)
    newcols = setdiff(names(qc), names(scn$cells))
    if (length(newcols)) scn$cells = merge(scn$cells, qc[, c("cell", newcols), with = FALSE], by = "cell", all.x = TRUE, sort = FALSE)
    scn$cells = scn$cells[match(colnames(scn$counts), cell)]
  }
  scn
}

# cells of a feature run, in the order of `all` (the sample's cells): `cells` as given (no QC); else the cells that pass
# cell_qc() (`qct`, with the run's `qc` arguments), restricted to the cells of `cell_meta` and, when it has a `qc_pass`
# column, to those with `qc_pass` TRUE (both must pass: the table's QC never overrides the run's `qc`).
# The R API and the CLI (which passes --cell-meta as `cell_meta`) select the same cells from the same table.
feature_cells = function(all, qct, cells = NULL, cell_meta = NULL) {
  if (!is.null(cells)) return(list(keep = intersect(all, cells), how = "given"))
  meta_qc = !is.null(cell_meta) && "qc_pass" %in% names(cell_meta)
  k = qct[qc_pass == TRUE]$cell
  if (meta_qc) k = intersect(k, as.character(cell_meta$cell)[as.logical(cell_meta$qc_pass) %in% TRUE])
  if (!is.null(cell_meta)) k = intersect(k, cell_meta$cell)
  list(keep = intersect(all, k), how = if (meta_qc) "QC pass and cell_meta qc_pass" else "QC pass")
}

# sex of a sample from its cells' `sex` column: the one value of the non-missing entries ([check_sex()] spellings);
# NA when there is none (missing or empty everywhere); stops when the cells disagree
cells_sex = function(x, where = "cell_features") {
  v = as.character(x)
  v = unique(v[!is.na(v) & nzchar(trimws(v))])
  if (!length(v)) return(NA_character_)
  s = unique(vapply(v, check_sex, character(1), where = where, USE.NAMES = FALSE))
  if (length(s) > 1)
    stop(sprintf("%s: the cells' `sex` column mixes both sexes (%s); pass `sex` for the sample", where,
                 paste(sprintf("\"%s\"", v), collapse = ", ")), call. = FALSE)
  s
}

# cell metadata of a feature run: `meta` (one row per kept cell: the QC n_frag / stratum and default labels) with the
# columns of `cell_meta` used as-is, except that a cell without a value for the QC-derived stratum / n_frag (not in
# cell_meta, e.g. with `cells` given, or NA there) keeps its QC value
apply_cell_meta = function(meta, cell_meta) {
  if (!is.null(cell_meta)) {
    cm = as.data.table(cell_meta)[match(meta$cell, cell)]
    for (cn in intersect(c("cell_type", "group", "flat", "fold", "stratum", "n_frag"), names(cm))) {
      v = cm[[cn]]
      if (cn == "stratum") v = fifelse(is.na(v), as.character(meta$stratum), as.character(v))
      if (cn == "n_frag") v = fifelse(is.na(v), as.numeric(meta$n_frag), as.numeric(v))
      meta[, (cn) := v]
    }
  }
  meta[, `:=`(cell_type = as.character(cell_type), group = as.character(group), flat = as.logical(flat),
              fold = as.integer(fold), stratum = as.character(stratum), n_frag = as.numeric(n_frag))]
  meta[]
}

#' Per-cell per-region features of a sample
#'
#' Reads an extraction directory, applies the cell QC ([cell_qc()]; with a `qc_pass` column of `cell_meta` a cell must
#' also have `qc_pass` TRUE; `cells` restricts the run to a set of cells instead, without QC), annotates arms, filters bins (default: [bin_filter()] on the kept cells, chr1-22
#' and chrX), normalises depth ([normalise_depth()]: `celltype` for normal tissue, `model` for tumours), builds the
#' regions ([region_set()]), the depth features ([depth_features()]), the molecule votes
#' ([read_votes()]) and the oriented allele features ([allelic_features()]), and joins them into the
#' long features table (`vignette("scphasedcn-outputs", package = "scphasedcn")`).
#' @param dir extraction directory (bins_100kb.h5ad, fragments.tsv.gz + .tbi, cells_qc.tsv)
#' @param cells NULL (QC-pass cells) or the cells to use (no QC)
#' @param cell_meta NULL or a data.table with `cell` and optional `qc_pass, cell_type, flat, fold, group, stratum,
#'   n_frag` (used as-is; NA cell_type = no type-matched reference; NA or missing `stratum` / `n_frag` keep the QC
#'   values). With `cells` NULL only cells of `cell_meta` are
#'   used, and with a `qc_pass` column (e.g. [make_cell_meta()]) only those with `qc_pass` TRUE that also pass
#'   [cell_qc()] with the `qc` arguments
#' @param method depth normalisation, "celltype" or "model" (any [normalise_depth()] method is accepted)
#' @param ref_cells reference cells of the normalisation (selector; NULL = all)
#' @param min_cells_type celltype: minimal reference cells of a type
#' @param use NULL or a logical bin filter
#' @param sex "female" / "male" (case-insensitive, or "F" / "M"; any other value stops; NULL: the cells' `sex` column, whose
#'   non-missing values must all give the same sex)
#' @param centromeres centromere table or path (Battenberg layout, see [arm_table()]; NULL: the hg38 table shipped with
#'   the package, [hg38_centromeres()])
#' @param classes region classes
#' @param events named events (see [region_set()])
#' @param null_exclude chromosomes / regions without null rows (known events)
#' @param phasing "auto", "chromosome", "block" (= "loo"), "loo" or "mirrored"
#' @param groups NULL, the name of a `cell_meta` column, a vector named by cell, or one label per kept
#'   cell: the pseudobulk groups of the block orientation (NA = all-cell pseudobulk)
#' @param mol_gap molecule gap (bp; default `params$mol_gap`, 2000); `params$mol_split_ps` and
#'   `params$votes_clean` set the molecule split by phase set and the blacklist filter of the votes ([read_votes()])
#' @param keep_votes keep the votes table in the result
#' @param layers bin layers to read (none are needed by the depth normalisation)
#' @param qc named list of [cell_qc()] arguments
#' @param params [caller_defaults()]
#' @param verbose log each step (the run notes) as it happens
#' @param exclude_atoms NULL or atom ids whose votes are left out of the allele features (see [allelic_features()])
#' @return `scn_features` list: `long` (features, keyed by cell and region), `regions`, `cells` (row order =
#'   the votes' cell index), `votes`, `atoms`, `win` ([window_matrices()]), `params`, `notes`, `use`, `phasing`
#' @export
cell_features = function(dir, cells = NULL, cell_meta = NULL, method = c("celltype", "model"), ref_cells = NULL, min_cells_type = 20,
                         use = NULL, sex = NULL, centromeres = NULL, classes = c("arm", "chrom", "5Mb", "10Mb", "20Mb", "50Mb"),
                         events = NULL, null_exclude = character(), phasing = "auto", groups = NULL, mol_gap = params$mol_gap,
                         keep_votes = TRUE, layers = character(), qc = list(), params = caller_defaults(),
                         verbose = FALSE, exclude_atoms = NULL) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  method = norm_method(method)
  phasing = match.arg(phasing, c("auto", "chromosome", "block", "loo", "mirrored"))
  sex = check_sex(sex, "cell_features")
  t0 = Sys.time()
  notes = character()
  note = function(...) {
    notes <<- c(notes, paste0(...))
    if (verbose) cli_log(sprintf("[%.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...)
    invisible()
  }
  centromeres = centromeres %||% hg38_centromeres()

  scn = read_sample_depth(dir, layers = layers)
  qct = do.call(cell_qc, c(list(scn), qc, list(params = params)))
  sel = feature_cells(scn$cells$cell, qct, cells, cell_meta)
  keep = sel$keep
  if (!length(keep)) stop("cell_features: no cells")
  note(sprintf("cells: %d of %d (%s)", length(keep), nrow(scn$cells), sel$how))
  sex = sex %||% (if ("sex" %in% names(scn$cells)) cells_sex(scn$cells$sex[match(keep, scn$cells$cell)]) else NULL)
  if (is.null(sex) || is.na(sex)) stop("cell_features: `sex` is needed")
  sex = check_sex(sex, "cell_features")
  scn = subset_scn(scn, cells = keep)

  # cell metadata: QC defaults, cell_meta used as-is
  meta = qct[match(keep, cell), .(cell, n_frag, stratum)]
  meta[, `:=`(cell_type = NA_character_, group = NA_character_, flat = NA, fold = NA_integer_)]
  if ("cell_type" %in% names(scn$cells)) meta[, cell_type := as.character(scn$cells$cell_type)]
  meta = apply_cell_meta(meta, cell_meta)

  # bins, usable bins, normalisation
  scn$bins = annotate_arms(scn$bins, centromeres)
  if (is.null(use)) use = bin_filter(scn$bins, scn$counts) & scn$bins$chrom %in% c(AUTOSOMES_CALLER, "chrX")
  note(sprintf("usable bins: %d of %d", sum(use), length(use)))
  norm = normalise_depth(scn, method = method, use = use, cell_type = meta$cell_type, ref_cells = ref_cells,
                         min_cells_type = min_cells_type)
  norm$ratio = NULL
  note(sprintf("depth: %s (ref cells: %s, min_cells_type %d)", method,
               if (is.null(ref_cells)) "all" else length(cell_selector(scn, ref_cells)), as.integer(min_cells_type)))

  regions = region_set(scn$bins, centromeres, classes = classes, events = events, sex = sex, null_exclude = null_exclude,
                       use = use, params = params)
  note(sprintf("regions: %s", paste(regions[, .N, by = class][, paste0(class, " ", N)], collapse = ", ")))
  dm = depth_mats(scn, norm, regions, params)
  win = window_matrices(scn, norm, 1e6)
  rm(norm); gc(verbose = FALSE)

  # alleles
  ra = regions[allelic_ok == TRUE]
  atoms = atoms_from_regions(scn$bins, ra)
  votes = read_votes(scn$fragments, keep, atoms, chroms = unique(ra$chrom), gap = mol_gap, tabix = params$tabix,
                     split_ps = isTRUE(params$mol_split_ps), clean = isTRUE(params$votes_clean))
  detected = detect_phasing(votes, atoms)
  orientation = switch(phasing, auto = if (detected == "chromosome") "chromosome" else "loo", block = "loo", phasing)
  note(sprintf("votes: %d rows (%d molecules); phasing detected %s, orientation %s", nrow(votes), sum(votes$nA) + sum(votes$nB),
               detected, orientation))
  g = groups
  if (!is.null(g)) {
    if (is.character(g) && length(g) == 1 && g %in% names(meta)) g = meta[[g]]
    else if (!is.null(names(g))) g = unname(g[keep])
    if (length(g) != length(keep)) stop("cell_features: `groups` must give one label per cell")
  }
  al = allelic_features(votes, atoms, ra, orientation = orientation, groups = g, params = params, exclude_atoms = exclude_atoms)
  if (length(exclude_atoms)) note(sprintf("allelic mask: votes of %d atoms left out", length(exclude_atoms)))

  # long table
  nr = nrow(regions); nc = length(keep)
  ri = match(al$region, regions$region)
  A = matrix(0, nr, nc); B = matrix(0, nr, nc); NB = matrix(0, nr, nc); Fm = matrix(NA_real_, nr, nc)
  idx = cbind(ri, al$cell)
  A[idx] = al$nA; B[idx] = al$nB; NB[idx] = al$n_blocks; Fm[idx] = al$f
  lab = if (nrow(al)) al$orient[1] else switch(orientation, chromosome = "anchored", loo = if (is.null(g)) "pseudobulk" else "clone", mirrored = "none")
  rep_r = function(v) rep(v, nc)
  rep_c = function(v) rep(v, each = nr)
  long = data.table(cell = rep_c(keep), region = rep_r(regions$region), class = rep_r(regions$class),
                    sd_class = rep_r(regions$sd_class), chrom = rep_r(regions$chrom),
                    x = as.vector(dm$x), e = as.vector(dm$e), log2ratio = as.vector(dm$log2ratio),
                    sd_null = as.vector(dm$sd_null), z = as.vector(dm$z),
                    nA = as.vector(A), nB = as.vector(B))
  rm(dm, A, B); gc(verbose = FALSE)
  long[, n_hap := nA + nB]
  long[, n_blocks := as.vector(NB)]
  long[, f := ifelse(n_hap > 0, as.vector(Fm), NA_real_)]
  rm(NB, Fm)
  long[, orient := ifelse(rep_r(regions$allelic_ok), lab, "none")]
  long[, `:=`(stratum = rep_c(meta$stratum), n_frag = rep_c(meta$n_frag), cell_type = rep_c(meta$cell_type),
              group = rep_c(meta$group), fold = rep_c(meta$fold), null_cell = rep_c(meta$flat),
              null_region = rep_r(regions$null_region), eligible = rep_r(regions$eligible),
              allelic_ok = rep_r(regions$allelic_ok))]
  setkeyv(long, c("cell", "region"))
  note(sprintf("long: %d rows (%d cells x %d regions); %.1f min", nrow(long), nc, nr,
               as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cells_out = merge(meta, qct[, .(cell, qc_pass)], by = "cell", sort = FALSE)[match(keep, cell)]
  check_regions(regions); check_atoms(atoms); check_votes(votes); check_features(long)
  structure(list(long = long, regions = regions, cells = cells_out, votes = if (keep_votes) votes else NULL, atoms = atoms,
                 win = win, params = params, notes = notes, use = use, phasing = c(detected = detected, orientation = orientation)),
            class = "scn_features")
}

#' Re-aggregate the allele features of a features object without some atoms
#'
#' Recomputes `nA`, `nB`, `n_hap`, `n_blocks` and `f` of the `allelic_ok` regions of `features$long` from
#' `features$votes` after dropping the votes of `exclude_atoms` (an allelic mask, [allelic_noisy_tiles()]), with
#' the features' own orientation (`features$phasing["orientation"]`). The depth columns are untouched.
#' @param features `scn_features` ([cell_features()]) with `votes` and `atoms`
#' @param exclude_atoms atom ids to leave out
#' @param groups pseudobulk groups of the block orientation (as in [cell_features()]; NULL in normal mode)
#' @param params [caller_defaults()]
#' @return the features object with updated allele columns and a note
#' @export
update_allelic_features = function(features, exclude_atoms, groups = NULL, params = caller_defaults()) {
  if (is.null(features$votes) || is.null(features$atoms)) stop("update_allelic_features: the features object needs votes and atoms")
  ra = as.data.table(features$regions)
  if ("allelic_ok" %in% names(ra)) ra = ra[allelic_ok %in% TRUE]
  orient = if (!is.null(features$phasing) && "orientation" %in% names(features$phasing)) unname(features$phasing[["orientation"]]) else "auto"
  if (orient == "block") orient = "loo"
  al = allelic_features(features$votes, features$atoms, ra, orientation = orient, groups = groups, params = params,
                        exclude_atoms = exclude_atoms)
  long = copy(as.data.table(features$long))
  cells = as.character(features$cells$cell)
  k = match(paste(long$cell, long$region), paste(cells[al$cell], al$region))
  aok = long$region %in% ra$region
  pick = function(v, old) fifelse(aok, fifelse(is.na(k), 0, as.numeric(v[k])), as.numeric(old))
  long[, `:=`(nA = pick(al$nA, nA), nB = pick(al$nB, nB))]
  long[, n_hap := nA + nB]
  long[, n_blocks := pick(al$n_blocks, n_blocks)]
  long[, f := fifelse(aok, fifelse(is.na(k) | n_hap <= 0, NA_real_, as.numeric(al$f[k])), f)]
  setkeyv(long, c("cell", "region"))
  features$long = long
  features$notes = c(features$notes, sprintf("allele features re-aggregated without %d masked atoms", length(exclude_atoms)))
  features
}

#' Subset a features object to a set of cells
#'
#' Keeps the rows of `long` and `cells` of the given cells and re-indexes the votes (whose `cell` is the row of
#' `cells`) to the kept rows, so that [update_allelic_features()] and the tumour mode still give each cell its own
#' votes; the 1 Mb window depth `win` keeps the kept cells' columns. Every other component is kept as it is.
#' @param features `scn_features` ([cell_features()])
#' @param cells cells to keep (others are ignored)
#' @return the subset `scn_features` (the object itself when every cell is kept)
#' @export
subset_features = function(features, cells) {
  ids = as.character(features$cells$cell)
  keep = ids %in% as.character(cells)
  if (all(keep)) return(features)
  idx = which(keep)
  features$cells = features$cells[keep]
  features$long = features$long[cell %in% ids[idx]]
  if (!is.null(features$votes)) {
    v = features$votes[cell %in% idx]
    v[, cell := as.integer(match(cell, idx))]
    features$votes = v
  }
  w = features$win
  if (is.list(w) && is.matrix(w$x) && !is.null(colnames(w$x))) {
    k = colnames(w$x) %in% ids[idx]
    w$x = w$x[, k, drop = FALSE]; w$e = w$e[, k, drop = FALSE]
    features$win = w
  }
  features
}

#' @export
print.scn_features = function(x, ...) {
  cat(sprintf("<scn_features> %d cells x %d regions (%d rows)\n", nrow(x$cells), nrow(x$regions), nrow(x$long)))
  cat(sprintf("  phasing: %s\n", paste(names(x$phasing), x$phasing, sep = " = ", collapse = ", ")))
  if (!is.null(x$votes)) cat(sprintf("  votes: %d rows, atoms: %d\n", nrow(x$votes), nrow(x$atoms)))
  for (n in x$notes) cat("  ", n, "\n", sep = "")
  invisible(x)
}
