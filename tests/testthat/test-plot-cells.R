# test-plot-cells.R: single-cell tracks (R/plot_cells.R) on the sim_tum2 tumour simulation, with a result built from
# the simulation's truth (clones = cell types; cloneA chr21 2:1, cloneB chr1 2:0 CN-LOH).

pc_cache = new.env()
pc_features = function() {
  if (is.null(pc_cache$feat)) {
    s = sim_tum2
    d = file.path(tempdir(), "plotcells_sim_tum2")
    dir.create(d, showWarnings = FALSE)
    sc = s$scn; sc$layers = list(); sc$obsm = list()
    invisible(capture.output(suppressMessages(write_bins_h5ad(sc, file.path(d, "bins_100kb.h5ad")))))
    fr = copy(s$fragments)[, cell := s$scn$cells$cell[cell]]
    write_fragments(fr, file.path(d, "fragments.tsv.gz"))
    fwrite(s$scn$cells, file.path(d, "cells_qc.tsv"), sep = "\t")
    cs = s$scn$bins[, .(len = max(end)), by = chrom]
    ce = cs[, .(chrom, cstart = 1, cen_left = 0.4 * len, cen_right = 0.4 * len, cend = len)]
    pc_cache$feat = suppressMessages(cell_features(d, cells = s$scn$cells$cell, centromeres = ce, sex = "male",
                                                   method = "model"))
  }
  pc_cache$feat
}

pc_result = function(feat) {
  ct = sim_tum2$scn$cells[, .(cell, clone = c(normal = "N0", cloneA = "CA", cloneB = "CB")[cell_type])]
  chroms = feat$win$windows[, .(start = min(start), end = max(end)), by = chrom]
  seg = rbindlist(lapply(c("N0", "CA", "CB"), function(k) chroms[, .(clone = k, chrom, start, end)]))
  seg[, state := fifelse(chrom == "chrX", "1:0", "1:1")]
  seg[clone == "CA" & chrom == "chr21", state := "2:1"]
  seg[clone == "CB" & chrom == "chr1", state := "2:0"]
  seg[, `:=`(n_win = 1L, nA = 0L, nB = 0L, cn = nA_of(state) + nB_of(state), log2 = 0, f = 0.5, llr = 0,
             confident = TRUE, label = fifelse(state == "2:1", "gain", fifelse(state == "2:0", "cnloh", "neutral")),
             ploidy = fifelse(chrom == "chrX", 1L, 2L))]
  cb = ct[clone == "CB", cell]
  c1 = chroms[chrom == "chr1"]
  # as caller_tumour() writes them: every cell carries its clone's segments; `called` (score above the confirmation
  # threshold) exactly when `source` is "confirmed", the others are "inherited" (reportable when the clone is confident)
  sc = c(10, 0.5)[seq_along(cb) %% 2 + 1]
  calls = data.table(cell = cb, region = "seg|CB|chr1", class = "segment", chrom = "chr1", start = c1$start,
                     end = c1$end, state = "cnloh", score = sc, thr = 1, called = sc > 1, level = "cell",
                     source = fifelse(sc > 1, "confirmed", "inherited"), reportable = TRUE)
  asg = ct[, .(cell, clone, posterior = 1, assigned = TRUE, ref_pool = clone == "N0")]
  clones = data.table(clone = c("N0", "CA", "CB"), c = c(1, 0.8, 0.8), a = c(1, 0.9, 0.9))
  regs = feat$regions
  new_scphasedcn("tumour", cells = sim_tum2$scn$cells[, .(cell, n_frag)], regions = regs, calls = calls,
                 assignments = asg, segments = seg, clones = clones)
}

test_that("cell_tracks returns depth, oriented BAF and fitted segments that follow the clone states", {
  feat = pc_features()
  res = pc_result(feat)
  tr = cell_tracks(feat, res)
  expect_s3_class(tr, "cell_tracks")
  expect_identical(attr(tr, "orient"), "clone")
  expect_setequal(names(tr), c("depth", "baf", "fit", "meta", "genome"))
  expect_setequal(unique(tr$depth$size), c(1e6, 5e6))
  asg = res$assignments
  ct = function(k) asg[clone == k, cell]
  # cloneB chr1 CN-LOH: clone-oriented BAF near 1 (major haplotype); normal cells near 0.5
  b = tr$baf[chrom == "chr1", .(f = sum(nA) / sum(n_hap)), by = .(cb = cell %in% ct("CB"))]
  expect_gt(b[cb == TRUE, f], 0.8)
  expect_lt(abs(b[cb == FALSE, f] - 0.5), 0.1)
  # cloneA chr21 gain: higher depth than normal cells
  d = tr$depth[size == 5e6 & chrom == "chr21", .(m = stats::median(log2)), by = .(ca = cell %in% ct("CA"))]
  expect_gt(d[ca == TRUE, m] - d[ca == FALSE, m], 0.3)
  # fitted segments: sources from the calls, expected values from the clone parameters
  f = tr$fit[chrom == "chr1" & cell %in% ct("CB")]
  expect_setequal(unique(f$source), c("confirmed", "inherited"))
  expect_equal(tr$meta[cell %in% ct("CB"), sum(n_inherited)], sum(!res$calls$called))
  # an inherited call that is not reportable (clone not confident) stays a plain clone segment
  r_nr = res; r_nr$calls = copy(res$calls)[called == FALSE, reportable := FALSE]
  expect_false("inherited" %in% cell_tracks(feat, r_nr, cells = ct("CB"))$fit$source)
  expect_equal(unique(f$f_exp), 0.5 + 0.9 * 0.5)
  expect_equal(unique(f$log2_exp), 0)   # 2:0 keeps the copy number
  g = tr$fit[chrom == "chr21" & cell %in% ct("CA")]
  expect_equal(unique(g$log2_exp), 0.8 * log2(3 / 2))
  expect_true(all(tr$fit[cell %in% ct("N0"), source] == "clone"))
  expect_equal(nrow(tr$meta), nrow(res$cells))
  # subset of cells
  tr2 = cell_tracks(feat, res, cells = ct("CB")[1:3])
  expect_setequal(unique(tr2$depth$cell), ct("CB")[1:3])
  expect_error(cell_tracks(feat, res, cells = "nope"), "not in the features")
})

test_that("plot_cell_tracks builds a plot and plot_cells_pdf writes a multi-page PDF", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  feat = pc_features()
  res = pc_result(feat)
  cb = res$assignments[clone == "CB", cell][1:3]
  tr = cell_tracks(feat, res, cells = cb)
  p = plot_cell_tracks(tr, cb[1], raster = FALSE)
  expect_s3_class(p, "patchwork")
  f = tempfile(fileext = ".png")
  grDevices::png(f, width = 1400, height = 600)
  print(p)
  grDevices::dev.off()
  expect_gt(file.size(f), 0)
  expect_error(plot_cell_tracks(tr, "nope"), "not in the tracks")
  pdf_f = tempfile(fileext = ".pdf")
  plot_cells_pdf(tr, pdf_f, per_page = 2, raster = FALSE)
  expect_gt(file.size(pdf_f), 0)
})

test_that("normal-mode tracks use the cell's events on a neutral baseline", {
  feat = pc_features()
  res = pc_result(feat)
  cells = res$cells$cell[1:4]
  ev = data.table(cell = cells[1], chrom = "chr21", start = 0, end = 2e7, state = "gain", top_region = "chr21",
                  score = 5, level = "cell", top_class = "chrom", n_regions = 1L)
  nres = new_scphasedcn("normal", cells = res$cells, regions = res$regions, events = ev)
  tr = cell_tracks(feat, nres, cells = cells, orient = "pseudobulk")
  expect_identical(attr(tr, "orient"), "pseudobulk")
  f = tr$fit[cell == cells[1] & source == "confirmed"]
  expect_identical(f$state, "2:1")
  expect_equal(f$log2_exp, log2(3 / 2))
  expect_true(all(tr$fit[cell != cells[1], label] == "neutral"))
  expect_identical(tr$fit[cell == cells[2] & chrom == "chrX", state], "1:0")
})

# ---- edge inputs and results of the current caller --------------------------------------------------

# the sim_tum2 extraction directory of pc_features(), and features with cell-type matched depth (normal mode's)
pc_dir = function() { pc_features(); file.path(tempdir(), "plotcells_sim_tum2") }
pc_features_celltype = function() {
  if (is.null(pc_cache$feat_ct)) {
    s = sim_tum2
    cs = s$scn$bins[, .(len = max(end)), by = chrom]
    ce = cs[, .(chrom, cstart = 1, cen_left = 0.4 * len, cen_right = 0.4 * len, cend = len)]
    pc_cache$feat_ct = suppressMessages(cell_features(pc_dir(), cells = s$scn$cells$cell, centromeres = ce, sex = "male",
                                                      method = "celltype"))
  }
  pc_cache$feat_ct
}
# the colours of the segments of a plot's first (depth) panel
pc_segment_colours = function(p) {
  ld = ggplot2::layer_data(p[[1]], 5)   # shading, hline, 1 Mb points, agg points, segments
  unique(ld$colour)
}

test_that("male chrX: the neutral segment follows the depth reference (cell-type matched depth sits at 0)", {
  s = sim_tum2
  n0 = s$scn$cells[cell_type == "normal", cell][1:20]
  for (m in c("celltype", "model")) {
    feat = if (m == "celltype") pc_features_celltype() else pc_features()
    nres = new_scphasedcn("normal", cells = s$scn$cells[, .(cell, n_frag)], regions = feat$regions)
    tr = cell_tracks(feat, nres, cells = n0, orient = "pseudobulk")
    obs = tr$depth[size == 5e6 & chrom == "chrX", stats::median(log2)]
    ex = unique(tr$fit[chrom == "chrX", log2_exp])
    expect_length(ex, 1L)
    expect_lt(abs(obs - ex), 0.25, label = sprintf("%s depth: chrX data %.2f vs fitted %.2f", m, obs, ex))
  }
})

test_that("normal-mode events on a haploid chromosome get haploid states", {
  feat = pc_features()
  res = pc_result(feat)
  cells = res$cells$cell[1:3]
  cx = feat$win$windows[chrom == "chrX", .(start = min(start), end = max(end))]
  ev = data.table(cell = cells[1:2], chrom = "chrX", start = cx$start, end = cx$end, state = c("loss", "gain"),
                  top_region = "chrX", score = 5, level = "cell", top_class = "chrom", n_regions = 1L)
  nres = new_scphasedcn("normal", cells = res$cells, regions = res$regions, events = ev)
  tr = cell_tracks(feat, nres, cells = cells, orient = "pseudobulk")
  f = tr$fit[chrom == "chrX" & source == "confirmed"]
  expect_identical(f[cell == cells[1], state], "0:0")
  expect_identical(f[cell == cells[2], state], "2:0")
  expect_true(all(is.na(f$f_exp)))
  expect_true(all(f$ploidy == 1L))
  # an event without coordinates takes those of its top region
  r21 = res$regions[region == "chr21"]
  ev2 = data.table(cell = cells[3], chrom = "chr21", start = NA_real_, end = NA_real_, state = "loss",
                   top_region = "chr21", score = 5, level = "cell", top_class = "chrom", n_regions = 1L)
  tr2 = cell_tracks(feat, new_scphasedcn("normal", cells = res$cells, regions = res$regions, events = ev2), cells = cells[3])
  g = tr2$fit[source == "confirmed"]
  expect_equal(c(g$start, g$end), c(r21$start, r21$end))
})

test_that("cell_tracks takes partial tumour results (no calls, clones without c / a, no allelic_ok, no ploidy)", {
  feat = pc_features()
  res = pc_result(feat)
  cb = res$assignments[clone == "CB", cell][1:2]
  r1 = res; r1$calls = NULL
  tr = cell_tracks(feat, r1, cells = cb)
  expect_true(all(tr$fit$source == "clone"))
  r2 = res; r2$clones = r2$clones[, .(clone)]
  r2$segments = copy(r2$segments)[clone == "CB" & chrom == "chr2", state := "0:0"]
  tr2 = cell_tracks(feat, r2, cells = cb)
  expect_equal(unique(tr2$fit[chrom == "chr1", f_exp]), 1)          # a = 1: 2:0 at its full major fraction
  expect_equal(unique(tr2$fit[chrom == "chr2", log2_exp]), log2(0.2 / 2))   # the clone-state emission's floor
  r3 = res; r3$regions = r3$regions[, !"allelic_ok"]; r3$segments = r3$segments[, !"ploidy"]
  tr3 = cell_tracks(feat, r3, cells = cb)
  expect_true(all(tr3$fit[chrom == "chrX", ploidy] == 2L))          # sex unknown without allelic_ok: diploid
  expect_identical(cell_tracks(feat, res, cells = c(cb, cb))$meta$cell, cb)
})

test_that("expected log2 lines are relative to the clone's ploidy baseline (clones$baseline)", {
  feat = pc_features()
  res = pc_result(feat)
  ca = res$assignments[clone == "CA", cell][1:2]
  # CA: baseline 3 copies; its chr21 2:1 then sits at the cell's median, a 2:2 on chr22 above it
  res$clones = copy(res$clones)[, baseline := c(2, 3, NA)]
  res$segments = copy(res$segments)[clone == "CA" & chrom == "chr22", state := "2:2"]
  for (pool in c(TRUE, FALSE)) {
    r = res
    if (!pool) r$assignments = copy(r$assignments)[, ref_pool := FALSE]   # model depth, no reference pool
    tr = cell_tracks(feat, r, cells = c(ca, res$assignments[clone == "CB", cell][1]))
    f = tr$fit[cell %in% ca]
    expect_equal(unique(f[chrom == "chr21", log2_exp]), 0, info = paste("pool", pool))
    expect_equal(unique(f[chrom == "chr22", log2_exp]), 0.8 * log2(4 / 3), info = paste("pool", pool))
    expect_equal(unique(f[chrom == "chr1", log2_exp]), 0.8 * log2(2 / 3), info = paste("pool", pool))
    # male chrX (one copy): against a pool of the same sex the baseline is 1 x 3 / 2 copies, on model depth 3
    expect_equal(unique(f[chrom == "chrX", log2_exp]), 0.8 * log2(1 / if (pool) 1.5 else 3), info = paste("pool", pool))
    # a clone with baseline NA keeps the diploid baseline
    expect_equal(unique(tr$fit[!cell %in% ca & chrom == "chr2", log2_exp]), 0)
  }
  # without the column: as before
  r0 = pc_result(feat)
  expect_equal(unique(cell_tracks(feat, r0, cells = ca)$fit[chrom == "chr21", log2_exp]), 0.8 * log2(3 / 2))
})

test_that("tumour tracks without a reference pool follow the caller's expected depth, whatever the depth method", {
  # a male tumour result without a pool on cell-type matched features: caller_tumour() scores the cells without a pool
  # against B copies (a male chrX 1:0 at c log2(1 / B)), so the tracks do too; with a pool the ploidy x B / 2
  feat = pc_features_celltype()
  expect_identical(scphasedcn:::features_depth_method(feat), "celltype")
  res = pc_result(feat)
  res$clones = copy(res$clones)[, baseline := c(2, 3, NA)]
  ca = res$assignments[clone == "CA", cell][1:2]
  n0 = res$assignments[clone == "N0", cell][1]
  r0 = res; r0$assignments = copy(r0$assignments)[, ref_pool := FALSE]
  f = cell_tracks(feat, r0, cells = c(ca, n0))$fit
  expect_equal(unique(f[cell %in% ca & chrom == "chrX", log2_exp]), 0.8 * log2(1 / 3))
  expect_equal(unique(f[cell %in% ca & chrom == "chrX", base]), 3)
  expect_equal(unique(f[cell %in% ca & chrom == "chr1", log2_exp]), 0.8 * log2(2 / 3))
  expect_equal(unique(f[cell == n0 & chrom == "chrX", log2_exp]), log2(1 / 2))
  expect_equal(unique(f[cell == n0 & chrom == "chr1", log2_exp]), 0)
  # with the pool: unchanged (ploidy x B / 2)
  fp = cell_tracks(feat, res, cells = c(ca, n0))$fit
  expect_equal(unique(fp[cell %in% ca & chrom == "chrX", log2_exp]), 0.8 * log2(1 / 1.5))
  expect_equal(unique(fp[cell == n0 & chrom == "chrX", log2_exp]), 0)
  # normal mode keeps the features' depth method: cell-type matched depth puts a male chrX at the ploidy
  nres = new_scphasedcn("normal", cells = res$cells, regions = res$regions)
  expect_equal(unique(cell_tracks(feat, nres, cells = n0, orient = "pseudobulk")$fit[chrom == "chrX", log2_exp]), 0)
})

test_that("plots colour every state label, including imbalance and labels the palette does not know", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  feat = pc_features()
  res = pc_result(feat)
  res$segments[clone == "CB" & chrom == "chr1", label := "imbalance"]
  res$segments[clone == "CB" & chrom == "chr2", `:=`(state = "2:1", label = "some_new_label")]
  cb = res$assignments[clone == "CB", cell][1]
  tr = cell_tracks(feat, res, cells = cb)
  p = plot_cell_tracks(tr, cb, raster = FALSE)
  cols = pc_segment_colours(p)
  expect_false(anyNA(cols))
  expect_false("grey50" %in% cols)                                   # the scale's NA colour
  expect_length(cols, 3L)                                            # neutral, imbalance, the unknown label
  # no known label at all
  tr$fit[, label := "unknown_state"]
  p2 = plot_cell_tracks(tr, cb, raster = FALSE)
  expect_false("grey50" %in% pc_segment_colours(p2))
})

test_that("plots build for a cell without BAF points and on one chromosome", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  feat = pc_features()
  res = pc_result(feat)
  cb = res$assignments[clone == "CB", cell][1]
  tr = cell_tracks(feat, res, cells = cb)
  p = plot_cell_tracks(tr, cb, raster = FALSE, min_hap = 1e9)
  expect_s3_class(ggplot2::ggplot_build(p[[2]]), "ggplot_built")
  t1 = tr
  for (k in c("depth", "baf", "fit")) t1[[k]] = t1[[k]][chrom == "chr21"]
  t1$genome = genome_layout(feat$win$windows[chrom == "chr21"])
  p1 = plot_cell_tracks(t1, cb, raster = FALSE)
  for (i in 1:3) expect_s3_class(ggplot2::ggplot_build(p1[[i]]), "ggplot_built")
})

test_that("plot_cells_pdf checks the cells before opening the file and closes only its own device", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  feat = pc_features()
  res = pc_result(feat)
  cb = res$assignments[clone == "CB", cell][1:2]
  tr = cell_tracks(feat, res, cells = cb)
  f = tempfile(fileext = ".pdf")
  expect_error(plot_cells_pdf(tr, f, cells = c(cb[1], "nope")), "not in the tracks")
  expect_false(file.exists(f))
  expect_error(plot_cells_pdf(tr, f, cells = character()), "no cells")
  expect_false(file.exists(f))
  png_f = tempfile(fileext = ".png")
  grDevices::png(png_f)
  dev0 = grDevices::dev.cur()
  plot_cells_pdf(tr, f, per_page = 1, raster = FALSE)
  expect_identical(grDevices::dev.cur(), dev0)
  grDevices::dev.off(dev0)
  expect_gt(file.size(f), 0)
})

test_that("cell_tracks and the plots run on call_scn() results, the same after write / read", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  feat = pc_features()
  res = suppressWarnings(call_scn(pc_dir(), "tumour", features = feat, cells = sim_tum2$scn$cells$cell, verbose = FALSE))
  tr = cell_tracks(feat, res)
  expect_identical(attr(tr, "orient"), "clone")
  expect_true(all(tr$fit$source %in% c("confirmed", "inherited", "clone")))
  cc = as.data.table(res$calls)[level == "cell" & called %in% TRUE]
  if (nrow(cc) && "source" %in% names(cc)) {
    # every called cell-level row of an assigned cell is found among the cell's fitted segments
    cc = cc[cell %in% tr$meta[!is.na(clone), cell]]
    hit = tr$fit[source != "clone"][cc, on = .(cell, chrom, start, end), nomatch = NULL]
    expect_equal(nrow(unique(hit[, .(cell, chrom, start, end)])), nrow(unique(cc[, .(cell, chrom, start, end)])))
  }
  # reportable inherited calls (called FALSE by construction) are drawn as inherited segments
  inh = as.data.table(res$calls)[level == "cell" & source == "inherited" & reportable %in% TRUE &
                                   cell %in% tr$meta[!is.na(clone), cell]]
  inh = unique(inh[, .(cell, chrom, start, end)])
  inh = inh[!tr$fit[source == "confirmed"], on = .(cell, chrom, start, end)]
  expect_gt(nrow(inh), 0)
  expect_equal(nrow(tr$fit[source == "inherited"]), nrow(inh))
  expect_equal(sum(tr$meta$n_inherited), nrow(inh))
  out = file.path(tempdir(), "plotcells_tum_out")
  write_scphasedcn(res, out)
  tr2 = cell_tracks(feat, read_scphasedcn(out))
  expect_equal(tr2, tr)
  cl = tr$meta[!is.na(clone)][order(-n_confirmed)]$cell[1]
  expect_s3_class(ggplot2::ggplot_build(plot_cell_tracks(tr, cl, raster = FALSE)[[1]]), "ggplot_built")
  # normal mode on cell-type matched depth
  fct = pc_features_celltype()
  nres = suppressWarnings(call_scn(pc_dir(), "normal", features = fct, cells = sim_tum2$scn$cells$cell, verbose = FALSE,
                                   embedding = FALSE))
  ntr = cell_tracks(fct, nres)
  expect_identical(attr(ntr, "mode"), "normal")
  ev = if (is.null(nres$events)) data.table(level = character()) else as.data.table(nres$events)
  if (nrow(ev[level == "cell"])) {
    expect_gt(nrow(ntr$fit[source == "confirmed"]), 0)
    expect_false(anyNA(ntr$fit$state))
  }
  c1 = ntr$meta[order(-n_confirmed)]$cell[1]
  expect_s3_class(ggplot2::ggplot_build(plot_cell_tracks(ntr, c1, raster = FALSE)[[3]]), "ggplot_built")
})

test_that("inst/scripts/plot_cells.R --help exits 0 without data", {
  script = system.file("scripts", "plot_cells.R", package = "scphasedcn")
  expect_true(file.exists(script))
  rs = file.path(R.home("bin"), "Rscript")
  out = suppressWarnings(system2(rs, c(script, "--help"), stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"))
  expect_true(any(grepl("--run", out)))
  expect_true(any(grepl("--orient", out)))
})

test_that("plot_cells.R fails with usage on bad input and plots a written run", {
  # the script loads the installed package: run it under R CMD check (or SCPHASEDCN_TEST_CLI=true after an install)
  skip_if_not(testthat::is_checking() || isTRUE(as.logical(Sys.getenv("SCPHASEDCN_TEST_CLI", "FALSE"))))
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  script = system.file("scripts", "plot_cells.R", package = "scphasedcn")
  rs = file.path(R.home("bin"), "Rscript")
  run = function(...) suppressWarnings(system2(rs, c(script, ...), stdout = TRUE, stderr = TRUE))
  bad = list(list(args = character(), msg = "--run is required"),
             list(args = c("--run", file.path(tempdir(), "no_such_run")), msg = "no result.rds"),
             list(args = c("--run", tempdir(), "--bogus", "1"), msg = "unknown option"),
             list(args = c("--run", tempdir(), "--orient", "sideways"), msg = "--orient must be"),
             list(args = c("--run", tempdir(), "--per-page", "0"), msg = "--per-page must be"))
  for (b in bad) {
    log = run(b$args)
    expect_false(is.null(attr(log, "status")), info = b$msg)
    expect_true(any(grepl(b$msg, log, fixed = TRUE)), info = paste(b$msg, paste(log, collapse = "\n")))
    expect_true(any(grepl("--help", log, fixed = TRUE)), info = b$msg)
  }
  feat = pc_features()
  res = pc_result(feat)
  rd = file.path(tempdir(), "plotcells_cli_run")
  unlink(rd, recursive = TRUE)
  write_scphasedcn(res, rd)
  ff = file.path(tempdir(), "plotcells_cli_features.rds")
  saveRDS(feat, ff)
  od = file.path(tempdir(), "plotcells_cli_out")
  cb = shQuote(paste(res$assignments[clone == "CB", cell][1:2], collapse = ";"))   # system2 runs a shell: quote the ";"
  log = run("--run", rd, "--features", ff, "--out", od, "--cells", cb, "--cores", "1")
  expect_null(attr(log, "status"))
  expect_true(file.exists(file.path(od, "CB.pdf")))
  expect_identical(nrow(fread(file.path(od, "cells.tsv"))), 2L)
})
