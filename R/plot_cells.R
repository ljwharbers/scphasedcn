# plot_cells.R
# What: single-cell genome tracks: per-cell depth log2 ratio and phased B-allele fraction with the fitted copy number
#       (clone segments inherited or confirmed in the cell, or the cell's own events in normal mode), and their plots.
# Why: the caller's tables say which states a cell carries; the tracks show the data behind them, one cell at a time.
# Inputs: a features object ([cell_features()]: `win`, `votes`, `atoms`, `cells`) and a caller result
#         ([call_scn()] / [read_scphasedcn()]).
# Outputs: [cell_tracks()] long tables; [plot_cell_tracks()] a patchwork plot; [plot_cells_pdf()] a multi-page PDF.
# See also: clones_states.R (clone_profiles(), the same depth scale for clone pseudobulks).

# nominal allele-specific state of a call label (normal mode events; relative to a diploid baseline)
label_states = function() {
  data.table(label = c("neutral", "loss", "cnloh", "gain", "gain_as", "amp", "del"),
             state = c("1:1", "1:0", "2:0", "2:1", "3:0", "3:1", "0:0"))
}

# colours of the state labels (ggsci npg)
state_colours = function() {
  c(neutral = "grey30", loss = "#4DBBD5", del = "#3C5488", cnloh = "#00A087", gain = "#F39B7F",
    gain_as = "#7E6148", amp = "#DC0000")
}

# chromosome offsets on a concatenated genome axis
genome_layout = function(windows) {
  g = as.data.table(windows)[, .(len = max(end)), by = chrom]
  g = g[order(chrom_rank(chrom))]
  g[, offset := c(0, cumsum(as.numeric(len))[-.N])]
  g[, `:=`(mid = offset + len / 2, shade = seq_len(.N) %% 2L == 0L)]
  g[]
}

#' Single-cell depth and B-allele fraction tracks
#'
#' Builds, for a set of cells, the data behind [plot_cell_tracks()]:
#' - `depth`: per cell and window the log2 ratio `log2((x + 0.5) / (e + 0.5))` of the features' 1 Mb windows and of
#'   `agg` windows. As in [clone_profiles()], relative to the reference pool when it has at least `min_ref` cells
#'   (then the copy-number baseline is the chromosome's ploidy) and centred on the cell's autosomal median, so the
#'   track is on the scale of the clone segments.
#' - `baf`: per cell and `agg` window the oriented haplotype-A fraction `nA / (nA + nB)` of the cell's molecule votes.
#'   Phase blocks are oriented with [orient_blocks()] (ties dropped) by a leave-one-out pseudobulk: of the cell's
#'   clone (`orient = "clone"`, the confirmation's orientation) or of all cells (`"pseudobulk"`); `"none"` keeps
#'   the phased haplotypes as they are (meaningful with chromosome-wide phasing). `"auto"`: `"clone"` in tumour mode,
#'   `"none"` in normal mode with chromosome-wide phasing, else `"pseudobulk"`.
#' - `fit`: the fitted segments of each cell. Tumour mode: the segments of the cell's clone, with `source`
#'   `confirmed` (called in the cell), `inherited` (the clone's call carried to the cell), or `clone` (neutral, or a
#'   clone call not reported in the cell). Normal mode: the cell's events (nominal states of their labels: loss 1:0,
#'   cnloh 2:0, gain 2:1, gain_as 3:0, amp 3:1, del 0:0) on a neutral baseline. Expected log2 `c * log2(cn / base)`
#'   (`c` the clone's depth contrast, 1 in normal mode) and expected major fraction `0.5 + a * (f - 0.5)` (`a` the
#'   clone's orientation attenuation, 1 in normal mode).
#' - `meta`: one row per cell (clone, posterior, assigned, n_frag, number of confirmed / inherited segments).
#' - `genome`: chromosome offsets for a concatenated genome axis.
#' @param features features object ([cell_features()]) with `win`, `votes` and `atoms`
#' @param result caller result ([call_scn()] or [read_scphasedcn()])
#' @param cells cells to return (default: the result's cells that are in the features)
#' @param agg size of the aggregated windows (bp; depth points and BAF)
#' @param orient phase-block orientation of the BAF (see Details)
#' @param min_ref minimal reference pool size to reference the depth to the pool
#' @return list of data.tables `depth` (cell, chrom, start, end, size, log2), `baf` (cell, chrom, start, end, nA,
#'   n_hap, baf, se), `fit` (cell, clone, chrom, start, end, state, cn, label, source, ploidy, base, log2_exp,
#'   f_exp), `meta` and `genome`; attributes `orient` and `mode`
#' @export
cell_tracks = function(features, result, cells = NULL, agg = 5e6, orient = c("auto", "clone", "pseudobulk", "none"),
                       min_ref = 30) {
  orient = match.arg(orient)
  result = relabel_legacy(result)
  win = features$win
  if (is.null(win) || is.null(features$votes) || is.null(features$atoms))
    stop("cell_tracks: `features` needs `win`, `votes` and `atoms` (cell_features(keep_votes = TRUE))")
  win = as_clone_win(win, colnames(win$x))
  fcells = colnames(win$x)
  if (is.null(cells)) cells = intersect(as.character(result$cells$cell), fcells)
  miss = setdiff(cells, fcells)
  if (length(miss)) stop(sprintf("cell_tracks: %d cells are not in the features (e.g. %s)", length(miss), miss[1]))
  mode = result$mode
  regs = as.data.table(result$regions)
  haploid = character()
  if ("allelic_ok" %in% names(regs)) {
    ao = regs[class %in% c("chrom", "arm") & allelic_ok == FALSE, unique(chrom)]
    haploid = intersect(ao, c("chrX", "chrY"))
  }
  # cell -> clone, reference pool, clone parameters
  asg = if (!is.null(result$assignments)) as.data.table(result$assignments) else data.table(cell = character(), clone = character())
  ref = if ("ref_pool" %in% names(asg)) intersect(asg[ref_pool %in% TRUE, cell], fcells) else character()
  use_ref = length(ref) >= min_ref
  cl = if (!is.null(result$clones)) as.data.table(result$clones) else data.table(clone = character(), c = numeric(), a = numeric())
  if (orient == "auto") {
    chromwide = !is.null(features$phasing) && identical(unname(features$phasing[["detected"]]), "chromosome")
    orient = if (mode == "tumour" && nrow(asg)) "clone" else if (chromwide) "none" else "pseudobulk"
  }

  # ---- depth ----
  depth_size = function(w, size) {
    lr = log2((w$x[, cells, drop = FALSE] + 0.5) / (w$e[, cells, drop = FALSE] + 0.5))
    if (use_ref) {
      rl = log2((rowSums(w$x[, ref, drop = FALSE]) + 0.5) / (rowSums(w$e[, ref, drop = FALSE]) + 0.5))
      lr = lr - rl
    }
    auto = !w$windows$chrom %in% c(haploid, "chrX", "chrY")
    lr = sweep(lr, 2, apply(lr[auto, , drop = FALSE], 2, stats::median, na.rm = TRUE))
    ok = rowSums(w$e[, cells, drop = FALSE]) > 0
    wt = as.data.table(w$windows)[ok]
    data.table(cell = rep(cells, each = nrow(wt)), chrom = rep(wt$chrom, length(cells)), start = rep(wt$start, length(cells)),
               end = rep(wt$end, length(cells)), size = size, log2 = as.vector(lr[ok, , drop = FALSE]))
  }
  wa = win_aggregate(win, agg)
  depth = rbind(depth_size(win, 1e6), depth_size(wa, agg))

  # ---- BAF ----
  aw = as.data.table(wa$windows)[, .(region = window, chrom, start, end = start + agg)]
  votes = as.data.table(features$votes)
  groups = NULL
  if (orient == "clone") groups = asg$clone[match(fcells, asg$cell)]
  if (orient == "none") {
    at = as.data.table(features$atoms)[, .(atom, chrom, astart = start)]
    v = at[votes[cell %in% match(cells, fcells), .(cell, atom, nA, nB)], on = "atom", nomatch = NULL]
    v[, start := floor(astart / agg) * agg]
    ob = v[, .(nA = sum(nA), nB = sum(nB)), by = .(cell, chrom, start)]
  } else {
    ob = as.data.table(orient_blocks(votes, aw, features$atoms, groups = groups, loo = TRUE))
    ob = ob[cell %in% match(cells, fcells)]
    ob = aw[ob, on = "region"][, .(cell, chrom, start, nA, nB)]
  }
  ob[, cell := fcells[cell]]
  ob[, end := start + agg]
  ob[, n_hap := nA + nB]
  ob = ob[n_hap > 0]
  ob[, baf := nA / n_hap]
  ob[, se := sqrt(pmax(baf * (1 - baf), 0.01) / n_hap)]
  baf = ob[, .(cell, chrom, start, end, nA, n_hap, baf, se)]

  # ---- fitted segments ----
  cl_of = asg$clone[match(cells, asg$cell)]
  base_of = function(ploidy) if (use_ref) ploidy else rep(2L, length(ploidy))
  if (mode == "tumour" && !is.null(result$segments)) {
    seg = as.data.table(result$segments)[, .(clone, chrom, start, end, state, label, ploidy)]
    fit = data.table(cell = cells, clone = cl_of)[!is.na(clone)][seg, on = "clone", allow.cartesian = TRUE, nomatch = NULL]
    calls = as.data.table(result$calls)
    calls = calls[level == "cell" & called %in% TRUE & cell %in% cells]
    if (nrow(calls) && "source" %in% names(calls)) {
      calls = calls[, .(src = fifelse(any(source == "confirmed"), "confirmed",
                                      fifelse(any(reportable %in% TRUE), "inherited", "clone"))),
                    by = .(cell, chrom, start, end)]
      fit[calls, src := i.src, on = .(cell, chrom, start, end)]
    } else fit[, src := NA_character_]
    fit[, source := fifelse(is.na(src), "clone", src)][, src := NULL]
    # cells without a clone: the neutral baseline
    no_cl = cells[is.na(cl_of)]
    if (length(no_cl)) fit = rbind(fit, neutral_baseline(no_cl, win$windows, haploid), fill = TRUE)
  } else {
    fit = neutral_baseline(cells, win$windows, haploid)
    ev = if (!is.null(result$events)) as.data.table(result$events)[level == "cell" & cell %in% cells] else data.table()
    if (nrow(ev)) {
      ls_ = label_states()
      ev = ev[, .(cell, clone = NA_character_, chrom, start, end, label = state,
                  state = ls_$state[match(state, ls_$label)], source = "confirmed")]
      ev[, ploidy := fifelse(chrom %in% haploid, 1L, 2L)]
      fit = rbind(fit, ev, fill = TRUE)
    }
  }
  fit[, cn := nA_of(state) + nB_of(state)]
  fit[, f := fifelse(cn > 0, nA_of(state) / cn, 0.5)]
  fit[, base := base_of(ploidy)]
  fit[, c_ := cl$c[match(clone, cl$clone)]][is.na(c_), c_ := 1]
  fit[, a_ := cl$a[match(clone, cl$clone)]][is.na(a_), a_ := 1]
  fit[, log2_exp := c_ * log2(pmax(cn, 0.25) / base)]
  fit[, f_exp := fifelse(ploidy == 1L, NA_real_, 0.5 + a_ * (f - 0.5))]   # haploid: no heterozygous SNPs
  fit[, c("c_", "a_", "f") := NULL]
  fit = fit[order(match(cell, cells), chrom_rank(chrom), start)]
  setcolorder(fit, c("cell", "clone", "chrom", "start", "end", "state", "cn", "label", "source", "ploidy", "base",
                     "log2_exp", "f_exp"))

  # ---- meta ----
  ct = as.data.table(result$cells)[match(cells, cell)]
  meta = data.table(cell = cells, clone = cl_of,
                    posterior = if ("posterior" %in% names(asg)) asg$posterior[match(cells, asg$cell)] else NA_real_,
                    assigned = if ("assigned" %in% names(asg)) asg$assigned[match(cells, asg$cell)] else NA,
                    n_frag = if ("n_frag" %in% names(ct)) ct$n_frag else NA_real_,
                    stratum = if ("stratum" %in% names(ct)) ct$stratum else NA_character_)
  nsrc = fit[, .(n_confirmed = sum(source == "confirmed"), n_inherited = sum(source == "inherited")), by = cell]
  meta = nsrc[meta, on = "cell"]
  setcolorder(meta, c("cell", "clone", "posterior", "assigned", "n_frag", "stratum"))
  meta[is.na(n_confirmed), n_confirmed := 0L][is.na(n_inherited), n_inherited := 0L]

  structure(list(depth = depth[], baf = baf[], fit = fit[], meta = meta[], genome = genome_layout(win$windows)),
            orient = orient, mode = mode, class = c("cell_tracks", "list"))
}
utils::globalVariables(c("src", "i.src", "c_", "a_", "log2_exp", "f_exp", "n_confirmed", "n_inherited", "baf", "se",
                         "shade", "offset", "len", "mid", "gx", "gx0", "gx1", "lt", "size"))

# one neutral segment per chromosome and cell
neutral_baseline = function(cells, windows, haploid) {
  g = as.data.table(windows)[, .(start = min(start), end = max(end)), by = chrom]
  out = data.table(cell = rep(cells, each = nrow(g)), clone = NA_character_, chrom = rep(g$chrom, length(cells)),
                   start = rep(g$start, length(cells)), end = rep(g$end, length(cells)))
  out[, ploidy := fifelse(chrom %in% haploid, 1L, 2L)]
  out[, `:=`(state = fifelse(ploidy == 1L, "1:0", "1:1"), label = "neutral", source = "clone")]
  out[]
}

#' @export
print.cell_tracks = function(x, ...) {
  cat(sprintf("<cell_tracks> %s mode, %d cells, BAF orientation: %s\n", attr(x, "mode"), nrow(x$meta), attr(x, "orient")))
  invisible(x)
}

# ---- plots -------------------------------------------------------------------------------------------

need_plot_pkgs = function() {
  for (p in c("ggplot2", "patchwork"))
    if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("plot_cell_tracks: needs the `%s` package", p), call. = FALSE)
}

theme_tracks = function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(), legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold", size = base_size),
                   plot.margin = ggplot2::margin(2, 5, 2, 5))
}

# add genome-axis coordinates
to_genome = function(d, g, cols = c("start", "end")) {
  d = copy(d)
  off = g$offset[match(d$chrom, g$chrom)]
  d[, gx0 := get(cols[1]) + off][, gx1 := get(cols[2]) + off][, gx := (gx0 + gx1) / 2]
  d[!is.na(gx)]
}

chrom_shading = function(g) {
  sh = g[shade == TRUE]
  ggplot2::geom_rect(data = sh, ggplot2::aes(xmin = offset, xmax = offset + len, ymin = -Inf, ymax = Inf),
                     inherit.aes = FALSE, fill = "grey93", colour = NA)
}

chrom_axis = function(g, labels = TRUE) {
  lab = sub("^chr", "", g$chrom)
  ggplot2::scale_x_continuous(breaks = g$mid, labels = if (labels) lab else NULL, expand = c(0, 0),
                              limits = c(0, max(g$offset + g$len)))
}

points_layer = function(..., raster = TRUE) {
  if (raster && requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(ggplot2::geom_point(...), dpi = 200)
  else ggplot2::geom_point(...)
}

#' Plot the tracks of one cell
#'
#' Three stacked tracks on a concatenated genome axis:
#' - depth log2 ratio (grey: 1 Mb windows; black: `agg` windows) with the expected log2 of the fitted segments;
#' - oriented B-allele fraction (`agg` windows, point size by the number of molecule votes) with the expected major /
#'   minor fractions (none on haploid chromosomes);
#' - a ribbon of the fitted nA:nB states (labelled when longer than 20 Mb), with segments confirmed in the cell
#'   outlined.
#'
#' Segments are coloured by state label; inherited segments are dashed, confirmed and clone (neutral) segments solid.
#' Needs ggplot2 and patchwork (ggrastr, when installed, rasterises the 1 Mb points).
#' @param tracks [cell_tracks()] output
#' @param cell cell barcode
#' @param ylim_log2 limits of the log2 track (points outside are squished onto the border)
#' @param min_hap minimal number of molecule votes of a BAF point
#' @param raster rasterise the 1 Mb points (ggrastr)
#' @param base_size base font size
#' @return a patchwork object
#' @export
plot_cell_tracks = function(tracks, cell, ylim_log2 = c(-2, 2), min_hap = 3, raster = TRUE, base_size = 11) {
  need_plot_pkgs()
  cc = cell
  g = tracks$genome
  m = tracks$meta[cell == cc]
  if (!nrow(m)) stop("plot_cell_tracks: cell not in the tracks: ", cc)
  d = to_genome(tracks$depth[cell == cc], g)
  b = to_genome(tracks$baf[cell == cc & n_hap >= min_hap], g)
  f = to_genome(tracks$fit[cell == cc], g)
  f[, lt := fifelse(source == "inherited", "inherited", fifelse(source == "confirmed", "confirmed", "clone"))]
  cols = state_colours()
  aes = ggplot2::aes
  # identical scales (limits, key widths) in every panel, so that patchwork collects one legend
  lt_lim = intersect(c("clone", "confirmed", "inherited"), unique(f$lt))
  lab_lim = intersect(names(cols), unique(f$label))
  lt_scale = ggplot2::scale_linetype_manual(values = c(clone = "solid", confirmed = "solid", inherited = "22"),
                                            limits = lt_lim, name = "segment")
  col_scale = ggplot2::scale_colour_manual(values = cols, limits = lab_lim, name = "state")
  title = sprintf("%s  |  clone %s%s  |  %s fragments  |  %d confirmed, %d inherited segments", cc,
                  ifelse(is.na(m$clone), "-", m$clone),
                  ifelse(is.na(m$posterior), "", sprintf(" (p %.2f%s)", m$posterior, ifelse(m$assigned %in% FALSE, ", unassigned", ""))),
                  ifelse(is.na(m$n_frag), "?", format(round(m$n_frag), big.mark = ",")), m$n_confirmed, m$n_inherited)

  # depth
  d1 = d[size == 1e6]; d5 = d[size != 1e6]
  sq = function(v) pmin(pmax(v, ylim_log2[1]), ylim_log2[2])
  pd = ggplot2::ggplot() + chrom_shading(g) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    points_layer(data = d1, aes(gx, sq(log2)), colour = "grey60", size = 0.15, alpha = 0.5, raster = raster) +
    ggplot2::geom_point(data = d5, aes(gx, sq(log2)), colour = "black", size = 0.5) +
    ggplot2::geom_segment(data = f, aes(x = gx0, xend = gx1, y = sq(log2_exp), yend = sq(log2_exp), colour = label,
                                        linetype = lt), linewidth = 1) +
    ggplot2::scale_y_continuous(limits = ylim_log2, name = "log2 ratio") +
    col_scale + lt_scale + chrom_axis(g, FALSE) + ggplot2::labs(x = NULL, title = title) + theme_tracks(base_size)

  # BAF
  fb = rbind(f[!is.na(f_exp), .(gx0, gx1, y = f_exp, label, lt)], f[!is.na(f_exp) & f_exp != 0.5, .(gx0, gx1, y = 1 - f_exp, label, lt)])
  pb = ggplot2::ggplot() + chrom_shading(g) +
    ggplot2::geom_hline(yintercept = 0.5, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(data = b, aes(gx, baf, size = n_hap), colour = "black", alpha = 0.6, stroke = 0) +
    ggplot2::geom_segment(data = fb, aes(x = gx0, xend = gx1, y = y, yend = y, colour = label, linetype = lt),
                          linewidth = 1) +
    ggplot2::scale_size_area(max_size = 1.6, name = "votes", guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), name = "BAF") +
    col_scale + lt_scale + chrom_axis(g, FALSE) + ggplot2::labs(x = NULL) + theme_tracks(base_size)

  # nA:nB ribbon
  pr =ggplot2::ggplot(f) +
    ggplot2::geom_rect(aes(xmin = gx0, xmax = gx1, ymin = 0, ymax = 1, fill = label),
                       colour = ifelse(f$lt == "confirmed", "black", NA), linewidth = 0.4) +
    ggplot2::geom_text(data = f[label != "neutral" & (gx1 - gx0) > 2e7], aes(x = gx, y = 0.5, label = state),
                       size = base_size / 4, colour = "white") +
    ggplot2::scale_fill_manual(values = cols, limits = lab_lim, name = "state", guide = "none") +
    ggplot2::scale_y_continuous(breaks = NULL, name = "nA:nB", expand = c(0, 0)) +
    chrom_axis(g, TRUE) + ggplot2::labs(x = NULL) + theme_tracks(base_size) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank())
  patchwork::wrap_plots(pd, pb, pr, ncol = 1, heights = c(1.3, 1, 0.25), guides = "collect")
}
utils::globalVariables("y")

#' Write the tracks of many cells to a multi-page PDF
#'
#' `per_page` cells per page ([plot_cell_tracks()]), with `grDevices::cairo_pdf(onefile = TRUE)`.
#' @param tracks [cell_tracks()] output
#' @param file PDF path
#' @param cells cells (default: all cells of the tracks)
#' @param per_page cells per page
#' @param width,height page size (inches)
#' @param ... passed to [plot_cell_tracks()]
#' @return `file`, invisibly
#' @export
plot_cells_pdf = function(tracks, file, cells = tracks$meta$cell, per_page = 4, width = 14, height = 4 * per_page, ...) {
  need_plot_pkgs()
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  pages = split(cells, ceiling(seq_along(cells) / per_page))
  for (pg in pages) {
    ps = lapply(pg, function(c_) plot_cell_tracks(tracks, c_, ...))
    if (length(ps) < per_page) ps = c(ps, rep(list(patchwork::plot_spacer()), per_page - length(ps)))
    print(patchwork::wrap_plots(ps, ncol = 1))
  }
  invisible(file)
}
