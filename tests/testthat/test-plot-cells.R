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
  calls = data.table(cell = cb, region = "seg|CB|chr1", class = "segment", chrom = "chr1", start = c1$start,
                     end = c1$end, state = "cnloh", score = 10, thr = 1, called = TRUE, level = "cell",
                     source = c("confirmed", "inherited")[seq_along(cb) %% 2 + 1], reportable = TRUE)
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
