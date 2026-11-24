# test-window-allelic.R: the windowed, orientation-free allele term (R/window_allelic.R): the unit LLR, the window
# units, the per-region sums with the shared-bias mask, its use in score_cells() and in the null-cell guard of
# select_null_cells(), and end to end on sim_norm_aneu (three cells with private whole-chromosome LOH) and its
# block-phased twin.

wa_centromeres = function(bins) {
  cs = bins[, .(len = max(end)), by = chrom]
  cs[, .(chrom, cstart = 1, cen_left = 0.4 * len, cen_right = 0.4 * len, cend = len)]
}
wa_cache = new.env()
wa_dir = function(name) {
  key = paste0(name, "_dir")
  if (is.null(wa_cache[[key]])) {
    d = file.path(tempdir(), paste0("winall_", name))
    unlink(d, recursive = TRUE)
    invisible(capture.output(suppressMessages(write_sim_dataset(get(name), d))))
    wa_cache[[key]] = d
  }
  wa_cache[[key]]
}
# small chr1 of 3 Mb: 100 kb atoms, arms p = [0, 1.5 Mb), q = [1.5, 3 Mb) and the chromosome
wa_toy = function() {
  atoms = data.table(atom = 1:30, chrom = "chr1", start = (0:29) * 1e5, end = (1:30) * 1e5, bin = 1:30)
  regions = data.table(region = c("chr1p", "chr1q", "chr1"), class = c("arm", "arm", "chrom"), sd_class = c("arm", "arm", "chrom"),
                       chrom = "chr1", start = c(0, 1.5e6, 0), end = c(1.5e6, 3e6, 3e6), n_bins = c(15, 15, 30),
                       n_use = c(15, 15, 30), eligible = TRUE, null_region = TRUE, allelic_ok = TRUE)
  list(atoms = atoms, regions = regions)
}

test_that("mirrored_llr equals the direct mixture formula, is 0 at n = 1 and has expectation 1 (as a ratio) under 1:1", {
  rho = 0.03; eps = 0.02
  bb = function(a, n, f) bb_loglik(a, n, f, rho, eps)
  a = c(0, 3, 7, 10); n = c(10, 10, 10, 10)
  ref = log(0.5 * exp(bb(a, n, 0.9)) + 0.5 * exp(bb(a, n, 0.1))) - bb(a, n, 0.5)
  expect_equal(mirrored_llr(a, n, 0.9, rho, eps), ref, tolerance = 1e-10)
  expect_equal(mirrored_llr(c(0, 1), c(1, 1), 0.98, rho, eps), c(0, 0), tolerance = 1e-12)
  for (k in c(2, 6, 15)) {
    x = 0:k
    expect_equal(sum(exp(bb(x, k, 0.5)) * exp(mirrored_llr(x, rep(k, k + 1), 0.98, rho, eps))), 1, tolerance = 1e-10)
  }
  expect_equal(win_col(c(0.98, 2 / 3)), c("wal_f980", "wal_f667"))
  expect_equal(win_f_grid(), 0.98)
  expect_equal(win_f_grid(modifyList(caller_defaults(), list(win_f = c(0.9, 0.98)))), c(0.9, 0.98))
  expect_equal(win_f_grid(modifyList(caller_defaults(), list(eps = 0.01))), c(0.98, 0.99))   # 1 - eps always in
})

test_that("window_units cuts the window grid at region boundaries and sums votes per cell x subwindow x phase set", {
  tt = wa_toy()
  votes = data.table(cell = c(1L, 1L, 1L, 1L, 2L), atom = c(1L, 10L, 12L, 16L, 25L), ps = c(5, 5, 5, 5, 5),
                     nA = c(1L, 2L, 1L, 3L, 0L), nB = c(0L, 1L, 0L, 0L, 2L))
  wu = window_units(votes, tt$atoms, tt$regions, size = 1e6)
  expect_equal(wu$windows$start, c(0, 1e6, 1.5e6, 2e6))
  expect_equal(wu$windows$end, c(1e6, 1.5e6, 2e6, 3e6))
  u = wu$units[order(cell, sw)]
  expect_equal(u$sw, c(1L, 2L, 3L, 4L))                       # atom 10 (0.9 Mb) in window 1; 12 (1.1 Mb) in 2; 16 (1.5 Mb) in 3
  expect_equal(u$n, c(4, 1, 3, 2))
  expect_equal(u$nA, c(3, 1, 3, 0))
  # two phase sets in one window are two units
  v2 = rbind(votes, data.table(cell = 1L, atom = 2L, ps = 9, nA = 2L, nB = 0L))
  expect_equal(nrow(window_units(v2, tt$atoms, tt$regions)$units[cell == 1 & sw == 1]), 2L)
})

test_that("window_allelic_features: private LOH scores high, balanced cells negative, shared-bias windows are masked", {
  tt = wa_toy()
  set.seed(3)
  mk = function(cell, f, atoms = 1:30, n = 4) {
    a = rbinom(length(atoms), n, f)
    data.table(cell = as.integer(cell), atom = atoms, ps = 1, nA = a, nB = n - a)
  }
  # cell 1: LOH on chr1q (all votes on one haplotype, either one per window: the orientation is irrelevant)
  votes = rbind(mk(1, 0.5, 1:15), mk(1, 1, 16:20), mk(1, 0, 21:30), rbindlist(lapply(2:20, function(i) mk(i, 0.5))))
  p = caller_defaults()
  wa = window_allelic_features(votes, tt$atoms, tt$regions, params = p, rho = 0.01)
  expect_true(all(c("cell", "region", "n_win", "wal_f980") %in% names(wa)))
  expect_equal(attr(wa, "rho"), 0.01)
  w1 = wa[cell == 1 & region == "chr1q", wal_f980]
  expect_gt(w1, 20)
  expect_lt(wa[cell == 1 & region == "chr1p", wal_f980], 0)
  expect_true(all(wa[cell > 1 & region == "chr1q", wal_f980] < 0))
  expect_equal(wa[cell == 1 & region == "chr1", wal_f980], wa[cell == 1 & region == "chr1p", wal_f980] + w1, tolerance = 1e-8)
  # a window where every cell leans to haplotype A (shared bias) is masked and contributes nothing
  vb = rbind(votes[!(atom %in% 21:30)], rbindlist(lapply(1:20, function(i) mk(i, 0.9, 21:30))))
  wb = window_allelic_features(vb, tt$atoms, tt$regions, params = p, rho = 0.01)
  expect_true(any(attr(wb, "windows")$masked))
  expect_equal(attr(wb, "windows")[masked == TRUE, start], 2e6)
  expect_true(all(wb[region == "chr1q", n_win] <= 1))          # only the 1.5-2 Mb subwindow is left on chr1q
  # rho fitted when not given; excluded atoms leave the sums
  wf = window_allelic_features(votes, tt$atoms, tt$regions, params = p)
  expect_true(attr(wf, "rho") > 0 && attr(wf, "rho") < 0.2)
  we = window_allelic_features(votes, tt$atoms, tt$regions, params = p, rho = 0.01, exclude_atoms = 16:30)
  expect_false(any(we$region == "chr1q"))
  # chrX is out of scope unless win_chrx
  ax = copy(tt$atoms)[, chrom := "chrX"]; rx = copy(tt$regions)[, `:=`(chrom = "chrX", region = sub("chr1", "chrX", region))]
  expect_equal(nrow(window_allelic_features(votes, ax, rx, params = p, rho = 0.01)), 0L)
  expect_gt(nrow(window_allelic_features(votes, ax, rx, params = modifyList(p, list(win_chrx = TRUE)), rho = 0.01)), 0L)
})

test_that("score_cells(window =) takes the better LLR on rows with a window column; pooled per-row columns stay oriented", {
  d = data.table(cell = c("a", "b", "c"), region = "r", class = "chrom", stratum = "all", z = c(-3, 0, 0.3), sd_null = 0.3,
                 nA = c(50, 50, 52), n_hap = c(100, 100, 100), wal_f980 = c(120, -40, NA))
  p = caller_defaults()
  s0 = score_cells(d, rho = 0.01, depth_model = "gauss", keep_rows = TRUE)
  s1 = score_cells(d, rho = 0.01, depth_model = "gauss", keep_rows = TRUE, window = "max")
  s2 = score_cells(d, rho = 0.01, depth_model = "gauss", keep_rows = TRUE, window = "only")
  expect_equal(s1$al_loss, c(120, s0$al_loss[2], s0$al_loss[3]))     # c: NA window -> oriented
  expect_equal(s1$al_cnloh, s1$al_loss)
  expect_equal(s1$al_gain, s0$al_gain)                                # no wal_f667 column: gains untouched
  expect_equal(s2$al_loss, c(120, -40, s0$al_loss[3]))
  expect_equal(s1$ar_loss_1, s0$ar_loss_1)
  expect_equal(s1$state_comb[1], "loss")
  expect_gt(s1$s_comb[1], s0$s_comb[1])
  expect_equal(s1$s_allelic[1], 120)
  expect_equal(score_cells(d, rho = 0.01, depth_model = "gauss", window = "max", sided = "one")$al_loss,
               score_cells(d, rho = 0.01, depth_model = "gauss", sided = "one")$al_loss)
  expect_equal(score_cells(d[, -"wal_f980"], rho = 0.01, depth_model = "gauss", window = "max")$s_comb, s0$s_comb)
  # gate: a windowed LLR at or below win_min (15) leaves the oriented term; win_min 0 lets it in
  dg = copy(d)[, wal_f980 := c(10, 10, 10)]
  expect_equal(score_cells(dg, rho = 0.01, depth_model = "gauss", window = "max")$al_loss, s0$al_loss)
  expect_equal(score_cells(dg, rho = 0.01, depth_model = "gauss", window = "max",
                           params = modifyList(p, list(win_min = 0)))$al_loss, pmax(s0$al_loss, 10))
  expect_equal(score_cells(dg, rho = 0.01, depth_model = "gauss", window = "only")$al_loss, c(10, 10, 10))
})

test_that("select_null_cells: the windowed guard keeps cells with a private LOH out of the null set", {
  set.seed(5)
  cells = sprintf("c%03d", 1:120)
  regs = c(paste0("chr", 1:6, "p"), paste0("chr", 1:6, "q"), paste0("chr", 1:6))
  d = CJ(cell = cells, region = regs)
  d[, class := ifelse(grepl("[pq]$", region), "arm", "chrom")]
  d[, chrom := sub("[pq]$", "", region)]
  d[, `:=`(stratum = "all", sd_null = 0.2, z = rnorm(.N), n_hap = 200)]
  d[, nA := as.numeric(rbinom(.N, n_hap, 0.5))]
  d[, wal_f980 := -abs(rnorm(.N, 60, 10))]
  d[cell == "c007" & chrom == "chr3", wal_f980 := 80]                  # private LOH of one cell, oriented f near 0.5
  p = caller_defaults()
  sel = select_null_cells(d, params = p, n_sim = 2e4)
  expect_false("c007" %in% sel$null_cells)
  expect_true(sel$calls[cell == "c007", win_guard])
  expect_equal(sum(sel$calls$win_guard), 1L)
  expect_true(all(sel$log$n_win_guard == 1L))
  expect_equal(sel$log$n_in[1], length(cells) - 1L)                    # out from the first fit on
  s1 = select_null_cells(d, params = p, n_sim = 2e4, n_iter = 1)
  expect_false("c007" %in% s1$null_cells)
  expect_equal(s1$log$n_in, length(cells) - 1L)
  sel0 = select_null_cells(d, params = modifyList(p, list(null_win_guard = FALSE)), n_sim = 2e4)
  expect_true("c007" %in% sel0$null_cells)
  expect_false(any(sel0$calls$win_guard))
  # columns of another fraction grid: warned, no guard
  dx = copy(d); setnames(dx, "wal_f980", "wal_f970")
  expect_warning(sx <- select_null_cells(dx, params = p, n_sim = 2e4), "no windowed null-cell guard")
  expect_true("c007" %in% sx$null_cells)
  # below the bound (log 1e4 = 9.2): not guarded
  d[cell == "c007" & chrom == "chr3", wal_f980 := 8]
  expect_true("c007" %in% select_null_cells(d, params = p, n_sim = 2e4)$null_cells)
})

test_that("call_scn on sim_norm_aneu: cells with private whole-chromosome LOH are called and are not null cells", {
  skip_on_cran()
  cen = wa_centromeres(sim_norm_aneu$scn$bins)
  d = wa_dir("sim_norm_aneu")
  cells = sim_norm_aneu$scn$cells$cell
  aneu = unique(cells[sim_norm_aneu$truth$nuclei[nucleus %in% sim_norm_aneu_nuclei, cell]])
  expect_length(aneu, 3L)
  fe = suppressWarnings(cell_features(d, cells = cells, centromeres = cen, sex = "female", method = "celltype"))
  expect_true(all(c("n_win", "wal_f980") %in% names(fe$long)))
  expect_true(is.list(fe$win_allelic) && fe$win_allelic$rho > 0)
  expect_true(any(grepl("^windowed allele term: ", fe$notes)))
  expect_true(all(is.na(fe$long[chrom == "chrX", wal_f980])))
  expect_true(all(!is.na(fe$long[chrom != "chrX" & allelic_ok == TRUE, wal_f980])))
  # the window term separates the aneuploid cells' LOH chromosomes from everything else
  w = fe$long[class == "chrom" & chrom %in% c("chr2", "chr21", "chr22")]
  expect_true(all(w[cell %in% aneu, wal_f980] > 20))
  expect_lt(max(fe$long[class == "chrom" & !cell %in% aneu & chrom == "chr22", wal_f980]), 9)   # no other event on chr22
  run = function(p) suppressWarnings(call_scn(d, "normal", features = fe, params = p, embedding = FALSE, verbose = FALSE))
  r1 = run(caller_defaults())
  r0 = run(modifyList(caller_defaults(), list(win_score = "none", null_win_guard = FALSE)))
  expect_false(any(r1$cells[cell %in% aneu, null_cell]))
  cl = r1$calls[level == "cell" & called & reportable & class == "chrom" & state %in% c("loss", "cnloh")]
  for (a in aneu) expect_gte(length(intersect(cl[cell == a, region], c("chr2", "chr21", "chr22"))), 2L)
  expect_true(any(grepl("^windowed allele term: ", r1$notes)))
  expect_equal(r1$null_fit$win_allelic$rho, fe$win_allelic$rho)
  expect_equal(null_fit_table(r1$null_fit)[kind == "allelic_rho_window", rho], fe$win_allelic$rho)
  expect_true(any(grepl("per-cell allele term: the better", r1$notes)))
  # no extra calls on the other cells' chromosomes without an event: the clean chromosome rows of non-carriers
  tr = sim_norm_aneu$truth
  carr = unique(cells[tr$nuclei$cell[unique(unlist(tr$carriers))]])
  clean = setdiff(cells, c(carr, aneu))
  nc = function(r) r$calls[level == "cell" & called & reportable & cell %in% clean & chrom != "chrX", .N]
  expect_lte(nc(r1), nc(r0) + 1L)
  # window columns that match no allelic row fraction (features of another eps): warned, noted, oriented scores
  fe2 = copy(fe$long); setnames(fe2, "wal_f980", "wal_f990")
  ng = modifyList(caller_defaults(), list(null_win_guard = FALSE))
  expect_warning(r3 <- caller_normal(fe2, params = ng, regions = fe$regions, pool_classes = character()),
                 "match no allelic row fraction")
  expect_true(any(grepl("oriented only \\(window columns wal_f990", r3$notes)))
  r0n = caller_normal(fe$long, params = modifyList(ng, list(win_score = "none")), regions = fe$regions,
                      pool_classes = character())
  expect_equal(r3$scores$s_comb, r0n$scores$s_comb)
})

test_that("update_allelic_features recomputes the window columns without the masked atoms", {
  skip_on_cran()
  cen = wa_centromeres(sim_norm_aneu$scn$bins)
  fe = suppressWarnings(cell_features(wa_dir("sim_norm_aneu"), cells = sim_norm_aneu$scn$cells$cell, centromeres = cen,
                                      sex = "female", method = "celltype"))
  ex = fe$atoms[chrom == "chr22" & start < 10e6, atom]                  # mask the first half of chr22
  fm = update_allelic_features(fe, ex)
  k = c("cell", "region")
  a = fe$long[class == "chrom" & chrom == "chr22", c(k, "n_win", "wal_f980"), with = FALSE]
  b = fm$long[class == "chrom" & chrom == "chr22", c(k, "n_win", "wal_f980"), with = FALSE]
  expect_true(all(b$n_win <= a$n_win) && any(b$n_win < a$n_win))
  aneu = unique(sim_norm_aneu$scn$cells$cell[sim_norm_aneu$truth$nuclei[nucleus %in% sim_norm_aneu_nuclei, cell]])
  expect_true(all(b[cell %in% aneu, wal_f980] < a[cell %in% aneu, wal_f980]))   # less LOH evidence left
  # other chromosomes: same units, only rho is refitted
  expect_gt(stats::cor(fm$long[class == "chrom" & chrom == "chr21", wal_f980], fe$long[class == "chrom" & chrom == "chr21", wal_f980]), 0.999)
  expect_equal(fm$long[class == "chrom" & chrom == "chr21", n_win], fe$long[class == "chrom" & chrom == "chr21", n_win])
  expect_true(is.list(fm$win_allelic) && fm$win_allelic$rho > 0)
})

test_that("block-only phasing: no window columns, a note, and the oriented scores", {
  skip_on_cran()
  cen = wa_centromeres(sim_norm_aneu_blk$scn$bins)
  d = wa_dir("sim_norm_aneu_blk")
  fe = suppressWarnings(cell_features(d, cells = sim_norm_aneu_blk$scn$cells$cell, centromeres = cen, sex = "female",
                                      method = "celltype"))
  expect_false(any(grepl("^wal_f", names(fe$long))))
  expect_null(fe$win_allelic)
  expect_true(any(grepl("windowed allele term: not computed", fe$notes)))
  r = suppressWarnings(call_scn(d, "normal", features = fe, embedding = FALSE, verbose = FALSE))
  expect_true(any(grepl("per-cell allele term: oriented only", r$notes)))
  # forced on (win_allelic TRUE) the columns are computed even on blocks
  fe2 = suppressWarnings(cell_features(d, cells = sim_norm_aneu_blk$scn$cells$cell, centromeres = cen, sex = "female",
                                       method = "celltype", params = modifyList(caller_defaults(), list(win_allelic = TRUE))))
  expect_true("wal_f980" %in% names(fe2$long))
})
