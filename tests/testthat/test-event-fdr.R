# test-event-fdr.R: the false-event accounting of normal mode (R/event_fdr.R): event_fdr() p- / q-values and the
# expected false calls / events on fake features with and without planted events, and null_selection_check() on a pure
# null and on a null that is too thin for a subset of cells.

fdr_cache = new.env()
fdr_run = function(planted) {
  key = if (planted) "planted" else "null"
  if (is.null(fdr_cache[[key]])) {
    rg = fake_regions(tiles = c("20Mb" = 2e7))
    ev = data.table(cell = 1:120, region = rep(c("chr8", "chr1p"), each = 60), state = rep(c("gain", "loss"), each = 60))
    ff = fake_features(1500, rg, events = if (planted) ev else NULL, seed = 61, hap_frac = 0.1)
    p = modifyList(caller_defaults(), list(min_null_rows = 1000))
    r = caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character(), params = p)
    fdr_cache[[key]] = list(r = r, f = event_fdr(r))
  }
  fdr_cache[[key]]
}

test_that("event_fdr adds columns only; its null groups are the thresholds' groups", {
  x = fdr_run(FALSE)
  r = x$r; f = x$f
  # every pre-existing column and row unchanged
  expect_equal(f$calls[, names(r$calls), with = FALSE], r$calls)
  expect_equal(f$events[, names(r$events), with = FALSE], r$events)
  expect_true(all(c("p", "q", "null_borrowed") %in% names(f$calls)))
  expect_true(all(c("null_cell", "p", "q", "null_borrowed") %in% names(f$events)))
  # the calibration groups hold exactly the null rows the s_comb thresholds were trained on (fold 0)
  th = r$thresholds[level == "cell" & score == "s_comb" & fold == 0]
  g = merge(f$groups, th[, .(class, stratum, n_thr = n_null)], by = c("class", "stratum"))
  expect_equal(nrow(g), nrow(f$groups))
  expect_equal(g$n_null, g$n_thr)
  expect_equal(f$groups$p_floor, 1 / (f$groups$n_null + 1))
  # tested rows: chromosomes and arms (the classes reportable per cell), not the 20 Mb tiles; chrX borrows the null
  cc = f$calls[level == "cell" & !is.na(p)]
  expect_setequal(unique(cc$class), c("chrom", "arm"))
  expect_equal(f$info$m, nrow(cc))
  expect_true(all(cc[chrom == "chrX", null_borrowed]) && !any(cc[chrom != "chrX", null_borrowed]))
  # events: null_cell from the cells table, q of the top region
  e = f$events[level == "cell"]
  expect_equal(e$null_cell, e$cell %in% r$cells[null_cell == TRUE, cell])
  expect_equal(e$q, cc$q[match(paste(e$cell, e$top_region), paste(cc$cell, cc$region))])
  r0 = r; r0["scores"] = list(NULL)
  expect_error(event_fdr(r0), "no scores")
})

test_that("event_fdr on a pure null: expected false ~ realised calls and events, p super-uniform, q not small", {
  f = fdr_run(FALSE)$f
  ef = f$expected_false
  expect_true(all(c("class", "state", "chrom_group", "n_rows", "n_tested", "n_called", "n_events", "null_per_1000",
                    "expected_false_calls", "expected_false_events", "false_frac_events") %in% names(ef)))
  expect_true(all(ef[chrom_group == "X", null_borrowed]))
  for (k in list(c("n_called", "expected_false_calls"), c("n_events", "expected_false_events"))) {
    n = sum(ef[[k[1]]]); E = sum(ef[[k[2]]])
    expect_gt(E, 3)
    expect_gte(n, stats::qpois(0.001, E)); expect_lte(n, stats::qpois(0.999, E))
  }
  cc = f$calls[level == "cell" & !is.na(p)]
  for (t in c(0.001, 0.01, 0.05, 0.2)) expect_lte(mean(cc$p <= t), t + 0.005)
  expect_lte(sum(cc$q < 0.05), 1)
  expect_gt(stats::median(cc$q), 0.5)
  expect_lte(f$info$n_q_below, 1)
})

test_that("event_fdr: q monotone in the score per class x stratum x fold; strong planted events get q < 0.05", {
  x = fdr_run(TRUE)
  r = x$r; f = x$f
  cc = f$calls[level == "cell" & !is.na(p) & p < 1]
  mono = cc[order(-score), .(p_ok = !is.unsorted(p), q_ok = !is.unsorted(q)), by = .(class, stratum, fold)]
  expect_true(all(mono$p_ok) && all(mono$q_ok))
  car = function(i) sprintf("cell_%04d", i)
  hit = f$calls[level == "cell" & called & reportable & ((region == "chr8" & state == "gain" & cell %in% car(1:60)) |
                                                           (region == "chr1p" & state == "loss" & cell %in% car(61:120)))]
  expect_gte(nrow(hit), 110)
  expect_true(all(hit$q < 0.05))
  e = f$events[level == "cell" & ((chrom == "chr8" & cell %in% car(1:60)) | (chrom == "chr1" & cell %in% car(61:120)))]
  expect_true(all(e$q < 0.05))
  expect_false(any(e$null_cell))
  # the expected false chromosome gains match the gains called outside the planted carriers
  ef = f$expected_false[class == "chrom" & state == "gain" & chrom_group == "A"]
  fp = nrow(f$events[level == "cell" & top_class == "chrom" & state == "gain" & chrom != "chrX" & !(chrom == "chr8" & cell %in% car(1:60))])
  expect_gte(fp, stats::qpois(0.001, ef$expected_false_events)); expect_lte(fp, stats::qpois(0.999, ef$expected_false_events))
  expect_lt(ef$false_frac_events, 0.5)
  # without the tail model p never goes below the resolution floor 1 / (n + 1) and the planted q-values are coarser
  f0 = event_fdr(r, tail = "none")
  c0 = f0$calls[level == "cell" & !is.na(p)]
  expect_gte(min(c0$p), min(f0$groups$p_floor))
  h0 = c0[paste(cell, region) %in% paste(hit$cell, hit$region)]
  expect_gt(stats::median(h0$q), stats::median(hit$q))
})

# a run in which thresholds fall back: min_null_rows 3000, the >=50k cells are not null cells (that stratum has no null
# rows), chrom <20k has too few; 30 >=50k cells carry a chr8 gain or a chr1p loss
fb_run = function() {
  if (is.null(fdr_cache$fallback)) {
    rg = fake_regions(tiles = c("20Mb" = 2e7))
    f0 = fake_features(1500, rg, seed = 61, hap_frac = 0.1)
    hi = head(unique(f0[stratum == ">=50k", cell]), 30)
    ev = data.table(cell = hi, region = rep(c("chr8", "chr1p"), each = 15), state = rep(c("gain", "loss"), each = 15))
    ff = fake_features(1500, rg, events = ev, seed = 61, hap_frac = 0.1)
    ff[, null_cell := null_cell & stratum != ">=50k"]
    p = modifyList(caller_defaults(), list(min_null_rows = 3000))
    r = caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character(), null_cells = "features", params = p)
    fdr_cache$fallback = list(r = r, f = event_fdr(r), ev = ev)
  }
  fdr_cache$fallback
}

# the calibration rows of every s_comb threshold group and fold: the null rows of the group's calibration key in the
# other folds (all of them for fold 0)
calib_counts = function(r) {
  sc = as.data.table(r$scores)
  cols = intersect(c("class", "sd_class", "stratum", "hstratum"), names(sc))
  ck = calibration_keys(sc[, c(cols, "is_null"), with = FALSE], r$params$min_null_rows)
  th = r$thresholds[level == "cell" & score == "s_comb" & class %in% c("chrom", "arm")]
  i = match(paste(th$class, th$stratum), paste(sc$class, sc$stratum))
  th[, n_cal := mapply(function(g, fo) { cr = ck$rows[[g]]; if (fo == 0L) length(cr) else sum(sc$fold[cr] != fo) },
                       ck$key[i], fold)]
  th[]
}

test_that("event_fdr: in fallback groups the p-values use the null rows the thresholds were trained on", {
  x = fb_run()
  r = x$r; f = x$f
  th = r$thresholds[level == "cell" & score == "s_comb" & class %in% c("chrom", "arm")]
  expect_true(th[stratum == ">=50k" & fold == 0, all(pooled_fallback)])
  expect_true(th[class == "chrom" & stratum == "<20k" & fold == 0, pooled_fallback])
  expect_false(th[class == "arm" & stratum == "<20k" & fold == 0, pooled_fallback])
  # groups: a fallback threshold group is calibrated on its class-wide null rows, as many as its threshold
  t0 = th[fold == 0, .(class, stratum = fifelse(pooled_fallback, "all", stratum), n_thr = n_null)]
  g = merge(f$groups, unique(t0), by = c("class", "stratum"))
  expect_equal(nrow(g), nrow(f$groups))
  expect_equal(g$n_null, g$n_thr)
  # per group and fold: as many calibration rows as the threshold of that group and fold was trained on
  cc = calib_counts(r)
  expect_equal(cc$n_cal, cc$n_null)
})

test_that("event_fdr: a stratum without null rows gets p-values from the class-wide null, never p = 1 when called", {
  x = fb_run()
  f = x$f
  cc = f$calls[level == "cell" & !is.na(p)]
  car = cc[paste(cell, region) %in% paste(x$ev$cell, x$ev$region) & called & reportable]
  expect_gte(nrow(car), 25)
  expect_true(all(car$stratum == ">=50k"))
  expect_true(all(is.finite(car$p) & car$p < 1e-3))
  expect_false(any(cc[called & reportable, p] == 1))
  expect_lt(f$info$p_floor_max, 1e-3)
})

test_that("event_fdr without fallback: each group's calibration rows are its class x stratum null rows; p is conformal", {
  r = fdr_run(FALSE)$r
  expect_false(any(r$thresholds[level == "cell", pooled_fallback]))
  cc = calib_counts(r)
  expect_equal(cc$n_cal, cc$n_null)
  # p (tail none) = (1 + #{null s' >= s'}) / (n + 1) over the null rows of the same class x stratum in the other folds
  f = event_fdr(r, tail = "none")
  sc = as.data.table(r$scores)
  x = f$calls[level == "cell"]
  k = match(paste(x$cell, x$region), paste(sc$cell, sc$region))
  rep = x$reportable %in% TRUE
  rep[!(x$called %in% TRUE) & x$class %in% c("chrom", "event") & x$state == "gain" & (sc$dl_gain[k] <= 0) %in% TRUE] = FALSE
  s = fifelse(rep & !is.na(x$score), x$score, -Inf)
  isn = sc$is_null[k] %in% TRUE
  fo = sc$fold[k]
  st = as.character(sc$stratum[k])
  tst = !is.na(x$p)
  expect_gt(sum(tst & is.finite(s) & x$p < 1), 100)
  pe = rep(NA_real_, nrow(x))
  for (g in unique(paste(x$class, st, fo)[tst])) {
    i = which(tst & paste(x$class, st, fo) == g)
    cal = s[isn & x$class == x$class[i[1]] & st == st[i[1]] & fo != fo[i[1]]]
    xs = sort(cal[is.finite(cal)])
    pe[i] = fifelse(is.finite(s[i]), (1 + length(xs) - findInterval(s[i], xs, left.open = TRUE)) / (length(cal) + 1), 1)
  }
  expect_equal(x$p[tst], pe[tst])
  expect_equal(x$q[tst], stats::p.adjust(pe[tst], "BH"))
})

# null_selection_check on a small genome (8 autosomes + chrX), 1,200 cells: a pure null (cached)
nsc_genome = function() fake_regions(sim_genome(setNames(rep(3e7, 9), c(paste0("chr", 1:8), "chrX"))), tiles = c("10Mb" = 1e7))
nsc_null = function() {
  if (is.null(fdr_cache$nsc0)) {
    rg = nsc_genome()
    f0 = fake_features(1200, rg, seed = 71)
    fdr_cache$nsc0 = list(rg = rg, f0 = f0, n0 = null_selection_check(f0, regions = rg, select_args = list(n_sim = 2e4)))
  }
  fdr_cache$nsc0
}

test_that("null_selection_check: ~1 / 1,000 on a pure null, clearly higher when the null is too thin for a subset of cells", {
  x0 = nsc_null()
  rg = x0$rg; f0 = x0$f0; n0 = x0$n0
  expect_true(all(c("class", "stratum", "n_rows", "n_exceed", "per_1000", "lo_per_1000", "hi_per_1000", "nominal_per_1000",
                    "insample_per_1000", "cv_per_1000") %in% names(n0)))
  hv = attr(n0, "halves")
  expect_equal(hv$select_chroms, c("chr1,chr3,chr5,chr7", "chr2,chr4,chr6,chr8"))
  a0 = n0[stratum == "all"]
  expect_setequal(a0$class, c("chrom", "arm", "10Mb"))
  expect_true(all(a0$per_1000 > 0.3 & a0$per_1000 < 2.5))
  expect_true(all(a0$lo_per_1000 <= a0$insample_per_1000))       # no class above its in-sample (CV) rate
  expect_equal(sum(n0[stratum != "all", n_rows]), sum(a0$n_rows))
  # 150 of 1,200 cells with an allelic overdispersion (rho 0.05) the robust parametric selection does not model: the
  # selection removes the worst of them, the selected half looks clean (in-sample ~1 / 1,000) and the other half does not
  f1 = copy(f0)
  bad = unique(f1$cell)[1:150]
  set.seed(5)
  f1[cell %in% bad, nA := as.numeric(rbinom(.N, n_hap, rbeta(.N, 0.5 * 0.95 / 0.05, 0.5 * 0.95 / 0.05)))]
  f1[, nB := n_hap - nA]
  n1 = null_selection_check(f1, regions = rg, select_args = list(n_sim = 2e4))
  ch0 = n0[class == "chrom" & stratum == "all"]; ch1 = n1[class == "chrom" & stratum == "all"]
  expect_gt(ch1$lo_per_1000, 2)
  expect_gt(ch1$per_1000, 2.5 * ch0$per_1000)
  expect_lt(ch1$insample_per_1000, 2)
  expect_error(null_selection_check(f0[chrom %in% c("chr1", "chr3")], regions = rg), "odd and even")
})

test_that("null_selection_check counts real events on the measured half: per_1000 bounds the null rate from above", {
  # 24 of 1,200 cells carry a whole-chromosome gain of chr2 (even half): selected as null cells on the odd half, their
  # chr2 rows are measured and exceed the thresholds, although the null itself is correct
  x0 = nsc_null()
  car = sprintf("cell_%04d", 1:24)
  ev = CJ(cell = car, region = x0$rg[chrom == "chr2", region])[, state := "gain"]
  f3 = fake_features(1200, x0$rg, events = ev, seed = 71)
  n3 = null_selection_check(f3, regions = x0$rg, select_args = list(n_sim = 2e4))
  ch0 = x0$n0[class == "chrom" & stratum == "all"]; ch3 = n3[class == "chrom" & stratum == "all"]
  expect_gt(ch3$per_1000, 2 * ch0$per_1000)
  expect_gt(ch3$lo_per_1000, ch0$per_1000)
  # the excess is about the carriers' chr2 rows over the measured rows, not a worse null on the selection half
  expect_lt(ch3$insample_per_1000, 2)
  expect_gt(ch3$n_exceed - ch0$n_exceed, 0.5 * length(car))
})

test_that("event_fdr: named events (thresholds borrowed from sd_class) get p-values calibrated like their thresholds", {
  # two named events on a pure null: one borrows from arm (whose per-cell gains are not reportable, an event gain is),
  # one from 10Mb (a class without per-cell call rows). Their p-values come from the sd_class's null rows scored as
  # event rows: about uniform, never all 1 and not far too small
  bins = sim_genome()
  rg = fake_regions(bins, tiles = c("10Mb" = 1e7, "20Mb" = 2e7))
  ev = data.table(region = c("ev_1q", "ev_2_10M"), class = "event", sd_class = c("arm", "10Mb"), chrom = c("chr1", "chr2"),
                  start = c(24e6, 1e7), end = c(6e7, 2e7))
  ev[, n_bins := vapply(seq_len(.N), function(i) sum(bins$chrom == chrom[i] & bins$start >= start[i] & bins$end <= end[i]), integer(1))]
  ev[, `:=`(n_use = n_bins, eligible = TRUE, null_region = FALSE, allelic_ok = TRUE)]
  rg = rbind(rg, ev, use.names = TRUE)
  ff = fake_features(2000, rg, seed = 21, hap_frac = 0.1)
  r = caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character())
  expect_true(all(r$thresholds[level == "cell" & class == "event", borrowed_from] %in% c("arm", "10Mb")))
  f = event_fdr(r)
  ce = f$calls[level == "cell" & class == "event"]
  expect_false(anyNA(ce$p))
  expect_true(all(ce$null_borrowed))
  s = ce[, .(n = .N, le_1e3 = mean(p <= 1e-3), le_0.01 = mean(p <= 0.01), le_0.1 = mean(p <= 0.1), all_1 = all(p == 1)), by = region]
  expect_setequal(s$region, ev$region)
  expect_false(any(s$all_1))
  expect_true(all(s$le_0.1 > 0.02 & s$le_0.1 < 0.15))
  expect_true(all(s$le_0.01 < 0.02))
  expect_true(all(s$le_1e3 < 0.004))
  # the borrowed calibration groups: null rows of the sd_class
  g = f$groups[borrowed == TRUE]
  expect_setequal(unique(g$class), c("arm", "10Mb"))
  expect_true(all(g$n_null > 0))
  # expected false calls of the event class: its sd_classes' null rows called as event rows
  ef = f$expected_false[class == "event"]
  expect_true(all(ef$null_borrowed))
  expect_true(all(ef$n_null == r$scores[is_null == TRUE & class %in% c("arm", "10Mb"), .N]))
})

test_that("event_fdr, default gain gate (dl): an uncalled row is reportable unless it is a chrom / event gain with dl_gain <= 0", {
  x = fdr_run(TRUE)
  r = x$r; f = x$f
  expect_identical(f$info$gain_gate, "dl")
  sc = r$scores
  ce = f$calls[level == "cell" & !is.na(p)]
  ce[, dl := sc$dl_gain[match(paste(cell, region), paste(sc$cell, sc$region))]]
  ce[, rep_old := fifelse(called %in% TRUE, reportable %in% TRUE,
                          reportable %in% TRUE & !(class %in% c("chrom", "event") & state == "gain" & !is.na(dl) & dl <= 0))]
  expect_true(all(ce[rep_old == FALSE, p] == 1))
  expect_identical(event_fdr(r, gain_gate = "dl"), f)
  # n_tested of the expected-false table: the reportable rows by that rule
  nt = ce[rep_old == TRUE & class == "chrom" & state == "gain" & chrom != "chrX", .N]
  expect_equal(f$expected_false[class == "chrom" & state == "gain" & chrom_group == "A", n_tested], nt)
})

test_that("event_fdr with gain_gate = 'z': a gain the z gate removes is not reportable for its p-value or the null accounting", {
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ev = data.table(cell = 1:60, region = "chr8", state = "gain")
  ff = fake_features(1500, rg, events = ev, seed = 61, hap_frac = 0.1)
  # noisy cells (4 x sd_null): the gain shift in z is small, so the null quantile of the z gate is above the dl
  # midpoint and the two gates differ
  ff[, sd_null := 4 * sd_null]
  p = modifyList(caller_defaults(), list(min_null_rows = 1000))
  r = caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character(), params = p, gain_gate = "z")
  f = event_fdr(r)                                            # gate inferred from the calls' gate columns
  expect_identical(f$info$gain_gate, "z")
  expect_identical(event_fdr(r, gain_gate = "z", gain_gate_q = 0.95), f)
  fd = event_fdr(r, gain_gate = "dl")                         # the old behaviour: the dl rule on uncalled rows
  sc = r$scores
  ce = f$calls[level == "cell" & class == "chrom" & state == "gain" & !is.na(p)]
  ks = match(paste(ce$cell, ce$region), paste(sc$cell, sc$region))
  zg = gain_gate_z(sc[ks, .(cell, region)], sc, 0.95, r$null_fit$tpar, p)
  ce[, `:=`(z_pass = (zg$z > zg$thr) %in% TRUE, dl = sc$dl_gain[ks])]
  # rows that pass the dl gate but fail the z gate: counted before, not now
  gated = ce[!(called %in% TRUE) & !z_pass & !(dl <= 0)]
  expect_gt(nrow(gated), 20)
  expect_true(all(gated$p == 1))
  pd = fd$calls[level == "cell" & class == "chrom" & state == "gain" & !is.na(p)]
  expect_true(all(pd$p[match(paste(gated$cell, gated$region), paste(pd$cell, pd$region))] < 1))
  # a called gain the gate dropped is not reportable and has p = 1
  expect_true(all(ce[called %in% TRUE & !(reportable %in% TRUE), p] == 1))
  # the null accounting and n_tested follow the z rule
  efz = f$expected_false[class == "chrom" & state == "gain" & chrom_group == "A"]
  efd = fd$expected_false[class == "chrom" & state == "gain" & chrom_group == "A"]
  expect_equal(efz$n_tested, ce[chrom != "chrX" & ((called %in% TRUE & reportable %in% TRUE) | (!(called %in% TRUE) & z_pass)), .N])
  expect_lt(efz$n_tested, efd$n_tested)
  expect_equal(efz$n_null_called, sum(r$calls[level == "cell" & class == "chrom" & state == "gain" & chrom != "chrX" &
                                               paste(cell, region) %in% paste(sc$cell, sc$region)[sc$is_null %in% TRUE],
                                             called %in% TRUE & reportable %in% TRUE]))
  # the strongest null calibration statistic of a gain uses the z rule too: no gated null gain is reportable
  nz = ce[cell %in% r$cells[null_cell == TRUE, cell] & !(called %in% TRUE) & !z_pass]
  expect_true(all(nz$p == 1))
})
