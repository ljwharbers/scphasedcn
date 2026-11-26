# plot_cells.R
# What: single-cell genome tracks: per-cell depth log2 ratio and phased B-allele fraction with the fitted copy number
#       (clone segments inherited or confirmed in the cell, or the cell's own events in normal mode), and their plots.
# Why: the caller's tables say which states a cell carries; the tracks show the data behind them, one cell at a time.
# Inputs: a features object ([cell_features()]: `win`, `votes`, `atoms`, `cells`) and a caller result
#         ([call_scn()] / [read_scphasedcn()]).
# Outputs: [cell_tracks()] long tables; [plot_cell_tracks()] a patchwork plot; [plot_cells_pdf()] a multi-page PDF.
# See also: clones_states.R (clone_profiles(), the same depth scale for clone pseudobulks).

# nominal allele-specific state of a call label (normal mode events) relative to the chromosome's ploidy: diploid
# (loss 1:0, cnloh 2:0, gain 2:1, ...) or haploid (male chrX / chrY: loss / del 0:0, gain 2:0, amp 3:0, as
# state_label() names haploid states); NA for a label without a nominal state
label_states = function() {
  data.table(label = c("neutral", "loss", "cnloh", "gain", "gain_as", "amp", "del"),
             state = c("1:1", "1:0", "2:0", "2:1", "3:0", "3:1", "0:0"),
             state_hap = c("1:0", "0:0", NA, "2:0", "2:0", "3:0", "0:0"))
}
event_state = function(label, ploidy) {
  ls_ = label_states()
  k = match(label, ls_$label)
  fifelse(ploidy == 1L, ls_$state_hap[k], ls_$state[k])
}

# colours of the state labels (ggsci npg); `labels`: also a colour for every other label present (e.g. the tumour
# label "imbalance", or labels of a newer caller), so that no segment falls on the scale's NA colour
state_colours = function(labels = NULL) {
  cols = c(neutral = "grey30", loss = "#4DBBD5", del = "#3C5488", cnloh = "#00A087", gain = "#F39B7F",
           gain_as = "#7E6148", amp = "#DC0000", imbalance = "#E64B35")
  extra = setdiff(unique(stats::na.omit(as.character(labels))), names(cols))
  if (length(extra)) {
    pal = c("#8491B4", "#91D1C2", "#B09C85", "#F0B400", "#6F2DA8")
    cols = c(cols, stats::setNames(rep_len(pal, length(extra)), extra))
  }
  cols
}

# depth normalisation of a features object, from its notes ("depth: <method> ..."); NA when not noted
features_depth_method = function(features) {
  l = grep("^depth: ", features$notes %||% character(), value = TRUE)
  if (length(l)) sub("^depth: ([A-Za-z0-9_]+).*$", "\\1", l[1]) else NA_character_
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
#'   and centred on the cell's autosomal median, so the track is on the scale of the clone segments. The
#'   copy-number baseline (`base`, the copy number at log2 0) follows the caller's expected depth. Tumour mode: the
#'   chromosome's ploidy with a reference pool, else 2 whatever the features' depth method (as [caller_tumour()]
#'   scores the cells: without a pool a male chrX sits at log2 -1). Normal mode: the chromosome's ploidy when the
#'   depth is relative to other cells of the sample (a reference pool, or the features' `celltype` / `loo` depth),
#'   else 2 (`model` depth); without a depth method in the features' notes, the ploidy.
#'   When the result's `clones` has a `baseline` column (the clone's ploidy baseline B, with absolute copy-number
#'   states), `base` is scaled by B / 2 for the cell's clone (NA or no column: 2, i.e. unchanged).
#' - `baf`: per cell and `agg` window the oriented haplotype-A fraction `nA / (nA + nB)` of the cell's molecule votes.
#'   Phase blocks are oriented with [orient_blocks()] (ties dropped) by a leave-one-out pseudobulk: of the cell's
#'   clone (`orient = "clone"`, the confirmation's orientation) or of all cells (`"pseudobulk"`); `"none"` keeps
#'   the phased haplotypes as they are (meaningful with chromosome-wide phasing). `"auto"`: `"clone"` in tumour mode,
#'   `"none"` in normal mode with chromosome-wide phasing, else `"pseudobulk"`.
#' - `fit`: the fitted segments of each cell. Tumour mode: the segments of the cell's clone, with `source`
#'   `confirmed` (a reportable call with `source` confirmed: one of the cell's events in the result's `events`),
#'   `inherited` (a reportable call with `source` inherited: the clone's segment carried to the cell without its own
#'   confirmation), or `clone` (neutral, or a clone segment without a reportable call in the cell). Normal mode: the cell's events (nominal states of their labels: loss 1:0,
#'   cnloh 2:0, gain 2:1, gain_as 3:0, amp 3:1, del 0:0; on a haploid chromosome loss 0:0, gain 2:0, amp 3:0) on a
#'   neutral baseline; events without coordinates take those of their `top_region`, else the whole chromosome.
#'   Expected log2 `c * log2(max(cn, 0.2) / base)` (`c` the clone's depth contrast, 1 in normal mode or when the
#'   result has none) and expected major fraction `0.5 + a * (f - 0.5)` (`a` the clone's orientation attenuation,
#'   likewise 1). Haploid chromosomes are the chrX / chrY rows with `allelic_ok` FALSE in the result's regions; a
#'   result without that column is treated as diploid throughout.
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
  cells = unique(as.character(cells))
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
  cl = if (!is.null(result$clones)) as.data.table(result$clones) else data.table(clone = character())
  for (k in c("c", "a", "baseline")) if (!k %in% names(cl)) cl[, (k) := rep(NA_real_, .N)]
  # copy-number baseline of the depth: the ploidy when the depth is relative to cells of the sample. Tumour mode
  # follows caller_tumour(), whose expected depth is pool-relative only with a pool, whatever the depth method
  dm = features_depth_method(features)
  rel_depth = use_ref || (mode != "tumour" && (if (is.na(dm)) mode == "normal" else dm %in% c("celltype", "loo")))
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
  base_of = function(ploidy) if (rel_depth) as.integer(ploidy) else rep(2L, length(ploidy))
  if (mode == "tumour" && !is.null(result$segments)) {
    seg = as.data.table(result$segments)
    if (!"label" %in% names(seg)) seg[, label := state_label(nA_of(state), nB_of(state), fifelse(chrom %in% haploid, 1L, 2L))]
    if (!"ploidy" %in% names(seg)) seg[, ploidy := fifelse(chrom %in% haploid, 1L, 2L)]
    seg = seg[, .(clone, chrom, start, end, state, label, ploidy)]
    fit = data.table(cell = cells, clone = cl_of)[!is.na(clone)][seg, on = "clone", allow.cartesian = TRUE, nomatch = NULL]
    calls = if (is.null(result$calls)) data.table() else as.data.table(result$calls)
    # caller_tumour() gives every cell its clone's segments: `source` "confirmed" (and `called`) when the cell's own
    # score passes, else "inherited" with `called` FALSE; a row counts when it is reportable (the result's events are
    # the reportable confirmed calls)
    need = c("cell", "chrom", "start", "end", "level", "source")
    calls = if (nrow(calls) && all(need %in% names(calls))) calls[level == "cell" & cell %in% cells] else data.table()
    if (nrow(calls)) {
      if (!"reportable" %in% names(calls)) calls[, reportable := TRUE]
      calls = calls[, .(src = fifelse(any(source == "confirmed" & reportable %in% TRUE), "confirmed",
                                      fifelse(any(source == "inherited" & reportable %in% TRUE), "inherited", NA_character_))),
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
      # events without coordinates (regions of a features-only run): those of the top region, else the chromosome
      ev[, `:=`(start = as.numeric(start), end = as.numeric(end))]
      if ("top_region" %in% names(ev) && all(c("region", "start", "end") %in% names(regs))) {
        k = match(ev$top_region, regs$region)
        ev[is.na(start), start := as.numeric(regs$start[k[is.na(start)]])]
        ev[is.na(end), end := as.numeric(regs$end[k[is.na(end)]])]
      }
      ce = as.data.table(win$windows)[, .(cs = min(start), cend = max(end)), by = chrom]
      ev[ce, `:=`(start = fifelse(is.na(start), as.numeric(i.cs), start), end = fifelse(is.na(end), as.numeric(i.cend), end)),
         on = "chrom"]
      ev = ev[, .(cell, clone = NA_character_, chrom, start, end, label = state, ploidy = fifelse(chrom %in% haploid, 1L, 2L),
                  source = "confirmed")]
      ev[, state := event_state(label, ploidy)]
      fit = rbind(fit, ev, fill = TRUE)
    }
  }
  fit[, cn := nA_of(state) + nB_of(state)]
  fit[, f := fifelse(cn > 0, nA_of(state) / cn, 0.5)]
  # clone ploidy baseline B (clones$baseline; NA / absent: 2): the cell's median depth is B copies, so every baseline
  # scales by B / 2 (relative depth: ploidy x B / 2, e.g. a male chrX against a same-sex reference; model depth: B)
  fit[, b_ := as.numeric(cl$baseline[match(clone, cl$clone)])][is.na(b_), b_ := 2]
  fit[, base := base_of(ploidy) * b_ / 2][, b_ := NULL]
  fit[, c_ := as.numeric(cl$c[match(clone, cl$clone)])][is.na(c_), c_ := 1]
  fit[, a_ := as.numeric(cl$a[match(clone, cl$clone)])][is.na(a_), a_ := 1]
  fit[, log2_exp := c_ * log2(pmax(cn, 0.2) / base)]   # 0.2: the floor of the clone-state depth emission
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
                         "shade", "offset", "len", "mid", "gx", "gx0", "gx1", "lt", "size", "cs", "cend", "i.cs", "i.cend", "b_"))

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

need_plot_pkgs = function(where = "plot_cell_tracks") {
  for (p in c("ggplot2", "patchwork"))
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("%s: needs the `%s` package (install.packages(\"%s\"))", where, p, p), call. = FALSE)
  if (utils::packageVersion("ggplot2") < "3.4.0")
    stop(sprintf("%s: needs ggplot2 >= 3.4.0 (`linewidth`), found %s", where, utils::packageVersion("ggplot2")), call. = FALSE)
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
  if (raster) ggrastr::rasterise(ggplot2::geom_point(...), dpi = 200) else ggplot2::geom_point(...)
}

# ---- graphics devices ----
# thin wrappers of the base devices (mocked in the tests)
has_cairo = function() isTRUE(capabilities("cairo"))
pdf_cairo = function(...) grDevices::cairo_pdf(...)
pdf_plain = function(...) grDevices::pdf(...)

# call `open` (a device-opening function) and return the number of the device it really opened (NULL: none). A
# warning or an error counts as a failure: a device it left open is closed again
open_device_checked = function(open) {
  before = grDevices::dev.list()
  ok = tryCatch({ open(); TRUE }, warning = function(w) FALSE, error = function(e) FALSE)
  new = setdiff(grDevices::dev.list(), before)
  if (!ok) {
    for (d in new) close_device(d)
    return(NULL)
  }
  if (length(new)) new[length(new)] else NULL
}

# close device `dev` if it is open (never the null device) and make `prev` current again if it is still open
close_device = function(dev, prev = NULL) {
  if (length(dev) == 1 && !is.na(dev) && dev > 1 && dev %in% grDevices::dev.list()) grDevices::dev.off(dev)
  if (length(prev) == 1 && !is.na(prev) && prev > 1 && prev %in% grDevices::dev.list()) grDevices::dev.set(prev)
  invisible()
}

# open a PDF device on `file`: cairo_pdf() when R reports cairo and it really opens a device (on some systems it
# reports the capability but cannot load cairo and opens nothing), else pdf(); the number of the opened device, now
# current. Stops when neither opens a device
open_pdf_device = function(file, width, height) {
  dev = if (has_cairo()) open_device_checked(function() pdf_cairo(file, width = width, height = height, onefile = TRUE))
  if (is.null(dev)) dev = open_device_checked(function() pdf_plain(file, width = width, height = height, onefile = TRUE))
  if (is.null(dev)) stop(sprintf("plot_cells_pdf: could not open a PDF device on %s (neither grDevices::cairo_pdf() nor grDevices::pdf())",
                                 file), call. = FALSE)
  grDevices::dev.set(dev)
  dev
}

# can ggrastr rasterise here? A rasterised layer is drawn once on a throw-away pdf(NULL) device (ggrastr renders with
# a Cairo-based device, which can be missing or fail to load; that shows only when the plot is drawn). Checked once
# per session; FALSE without ggrastr
.scphasedcn_raster = new.env(parent = emptyenv())
raster_probe = function() {
  p = ggplot2::ggplot(data.frame(x = 1:2, y = 1:2), ggplot2::aes(x, y)) + ggrastr::rasterise(ggplot2::geom_point(), dpi = 20)
  prev = grDevices::dev.cur()
  dev = open_device_checked(function() pdf_plain(NULL))
  if (is.null(dev)) stop("no device to draw on")
  on.exit(close_device(dev, prev), add = TRUE)
  grDevices::dev.set(dev)
  print(p)
  TRUE
}
raster_ok = function() {
  if (!requireNamespace("ggrastr", quietly = TRUE)) return(FALSE)
  if (is.null(.scphasedcn_raster$ok)) {
    why = tryCatch({ raster_probe(); NULL }, warning = function(w) conditionMessage(w), error = function(e) conditionMessage(e))
    .scphasedcn_raster$ok = is.null(why)
    if (!is.null(why)) message("raster = TRUE: ggrastr cannot rasterise here (", why, "); the 1 Mb points are drawn unrasterised")
  }
  .scphasedcn_raster$ok
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
#' Needs ggplot2 and patchwork (ggrastr, when installed, rasterises the 1 Mb points; when it cannot rasterise here,
#' e.g. without a working Cairo, the points are drawn unrasterised, with one message per session).
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
  raster = isTRUE(raster) && raster_ok()
  cc = cell
  g = tracks$genome
  m = tracks$meta[cell == cc]
  if (!nrow(m)) stop("plot_cell_tracks: cell not in the tracks: ", cc)
  d = to_genome(tracks$depth[cell == cc], g)
  b = to_genome(tracks$baf[cell == cc & n_hap >= min_hap], g)
  f = to_genome(tracks$fit[cell == cc], g)
  f[, lt := fifelse(source == "inherited", "inherited", fifelse(source == "confirmed", "confirmed", "clone"))]
  if (length(cc) != 1) stop("plot_cell_tracks: `cell` must be one cell barcode")
  cols = state_colours(f$label)
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
#' `per_page` cells per page ([plot_cell_tracks()]), with `grDevices::cairo_pdf(onefile = TRUE)`, or
#' `grDevices::pdf()` when R has no cairo or `cairo_pdf()` opens no device (a warning or an error counts as that). The
#' cells are checked before the file is opened; only `file` (and its directory) is written. Only the device opened
#' here is closed, also when plotting fails, and the device that was current before is made current again.
#' @param tracks [cell_tracks()] output
#' @param file PDF path
#' @param cells cells (default: all cells of the tracks)
#' @param per_page cells per page
#' @param width,height page size (inches)
#' @param ... passed to [plot_cell_tracks()]
#' @return `file`, invisibly
#' @export
plot_cells_pdf = function(tracks, file, cells = tracks$meta$cell, per_page = 4, width = 14, height = 4 * per_page, ...) {
  need_plot_pkgs("plot_cells_pdf")
  cells = unique(as.character(cells))
  if (!length(cells)) stop("plot_cells_pdf: no cells to plot")
  miss = setdiff(cells, tracks$meta$cell)
  if (length(miss)) stop(sprintf("plot_cells_pdf: %d cells are not in the tracks (e.g. %s)", length(miss), miss[1]))
  if (length(per_page) != 1 || !is.finite(per_page) || per_page < 1) stop("plot_cells_pdf: `per_page` must be a number >= 1")
  per_page = as.integer(per_page)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  prev = grDevices::dev.cur()
  dev = open_pdf_device(file, width, height)
  on.exit(close_device(dev, prev), add = TRUE)
  pages = split(cells, ceiling(seq_along(cells) / per_page))
  for (pg in pages) {
    ps = lapply(pg, function(c_) plot_cell_tracks(tracks, c_, ...))
    if (length(ps) < per_page) ps = c(ps, rep(list(patchwork::plot_spacer()), per_page - length(ps)))
    grDevices::dev.set(dev)            # the page goes to this file, whatever plotting did to the current device
    print(patchwork::wrap_plots(ps, ncol = 1))
  }
  invisible(file)
}
