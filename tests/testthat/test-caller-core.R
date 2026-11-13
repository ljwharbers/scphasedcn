# test-caller-core.R: parity of the caller primitives with the frozen reference scorer (helper-refscorer.R), the
# depth strata, schema validators, the scphasedcn object, tile names and the simulated fragment phase sets.

core_regions = fake_regions()
core_feat = fake_features(240, core_regions,
                          events = data.table(cell = c(1:20, 21:40, 41:60), region = rep(c("chr1p", "chr8", "chr2q"), each = 20),
                                              state = rep(c("loss", "gain", "cnloh"), each = 20)))

test_that("fake tables pass the validators", {
  expect_silent(check_regions(core_regions))
  expect_silent(check_features(core_feat))
  atoms = fake_atoms()
  expect_silent(check_atoms(atoms))
  v = fake_votes(30, atoms, clones = rep(c("A", "B"), 15),
                 events = data.table(clone = "B", chrom = "chr22", start = 0, end = NA, cnA = 1L, cnB = 0L))
  expect_silent(check_votes(v))
  # clone B lost one haplotype of chr22: after undoing the block flips, f ~ 1 in clone B and ~ 0.5 in clone A
  b = attr(v, "blocks")
  vv = merge(v, fake_atoms()[, .(atom, chrom)], by = "atom")
  vv = merge(vv, b[, .(chrom, ps, flip)], by = c("chrom", "ps"))
  expect_equal(nrow(vv), nrow(v))
  expect_gt(mean(b$flip), 0.3); expect_lt(mean(b$flip), 0.7)
  vv[flip == TRUE, c("nA", "nB") := .(nB, nA)]
  vv = vv[chrom == "chr22"]
  vv[, clone := attr(v, "clone")[cell]]
  fr = vv[, .(f = sum(nA) / sum(nA + nB)), by = clone]
  expect_gt(fr[clone == "B", f], 0.95)
  expect_lt(abs(fr[clone == "A", f] - 0.5), 0.05)
  # events shift z and f in the right direction
  m = core_feat[, .(z = median(z), f = median(f, na.rm = TRUE)), by = truth]
  expect_lt(m[truth == "loss", z], -3)
  expect_gt(m[truth == "gain", z], 2)
  expect_gt(m[truth == "cnloh", f], 0.9)
  expect_equal(core_feat$z, core_feat$log2ratio / pmax(core_feat$sd_null, 0.05))
  expect_equal(core_feat$log2ratio, log2((core_feat$x + 0.5) / (core_feat$e + 0.5)))
})

test_that("ll_depth_shift equals the depth part of the reference scorer (t and Gaussian)", {
  d = copy(core_feat)
  d[sample(.N, 50), z := NA]
  d[sample(.N, 20), sd_null := NA]
  tpar = ref_fit_t(d[null_cell == TRUE], n_max = 2000)
  shift = caller_defaults()$shifts
  sdp = pmax(fifelse(is.na(d$sd_null), 1, d$sd_null), 0.05)
  for (model in c("t", "gauss")) {
    ref = ref_score(d, shift, model = model, tpar = tpar)
    tp = if (model == "t") match_depth_null(d, tpar) else NULL
    l_loss = ll_depth_shift(d$z, shift[["loss"]] / sdp, tp)
    l_gain = ll_depth_shift(d$z, shift[["gain"]] / sdp, tp)
    expect_lte(max(abs(l_loss - ref$dl_loss)), 1e-12)
    expect_lte(max(abs(l_gain - ref$dl_gain)), 1e-12)
    expect_true(all(l_loss[is.na(d$z)] == 0))
  }
  # scalar tpar (list) and Gaussian formula
  expect_equal(ll_depth_shift(c(1, NA), 2), c(2 * 1 - 2, 0))
  tp1 = list(t_m = 0.1, t_s = 1.2, t_df = 6)
  expect_equal(ll_depth_shift(1.5, -1, tp1), dt((2.5 - 0.1) / 1.2, 6, log = TRUE) - dt((1.5 - 0.1) / 1.2, 6, log = TRUE))
  expect_error(match_depth_null(data.table(class = "nope", stratum = "<20k"), tpar), "no t fit")
})

test_that("ll_allelic equals the reference oriented LLRs and the attenuation grid", {
  set.seed(3)
  n = c(0, rpois(500, 30), 1, 2)
  nA = rbinom(length(n), n, sample(c(0.5, 0.7, 0.95, 0.05), length(n), replace = TRUE))
  rho = 0.004
  ref = ref_allelic_llr(nA, n, rho)
  eps = caller_defaults()$eps
  expect_identical(ll_allelic(nA, n, 1 - eps, rho), ref$llr_20)
  expect_identical(ll_allelic(nA, n, 2 / 3, rho), ref$llr_21)
  expect_identical(ll_allelic(nA, n, 3 / 4, rho), ref$llr_31)
  fa = caller_defaults()$f_att
  expect_identical(ll_allelic(nA, n, fa, rho), ref_llr_grid(nA, n, rho, fa))
  expect_equal(ll_allelic(0, 0, 0.9, rho), 0)
  # one-sided: no credit for the mirrored fraction
  one = ll_allelic(c(0, 20), c(20, 20), 0.9, rho, two_sided = FALSE)
  expect_lt(one[1], 0); expect_gt(one[2], 0)
  expect_equal(ll_allelic(c(0, 20), c(20, 20), 0.9, rho), rep(one[2], 2))
})

test_that("ll_allelic(per_row = TRUE) scores each row at its own f (the former clones.R ll_allelic_rows)", {
  set.seed(4)
  n = c(0, NA, rpois(300, 25), 3)
  nA = rbinom(length(n), ifelse(is.na(n), 0, n), 0.6)
  f = runif(length(n), 0.3, 0.99)
  rho = 0.01; eps = 0.02
  ref = bb_loglik(nA, n, f, rho, eps) - bb_loglik(nA, n, 0.5, rho, eps)    # ll_allelic_rows' formula
  ref[is.na(n) | n == 0 | is.na(ref)] = 0
  pr = ll_allelic(nA, n, f, rho, eps, two_sided = FALSE, per_row = TRUE)
  expect_identical(pr, ref)
  expect_identical(pr[1:2], c(0, 0))
  # each row equals the grid mode at that row's f alone
  k = 3:12
  expect_identical(pr[k], vapply(k, function(i) ll_allelic(nA[i], n[i], f[i], rho, eps, two_sided = FALSE), numeric(1)))
  # a scalar f is recycled; two-sided per row takes the better of f and 1 - f
  expect_identical(ll_allelic(nA, n, 0.7, rho, eps, two_sided = FALSE, per_row = TRUE),
                   ll_allelic(nA, n, rep(0.7, length(n)), rho, eps, two_sided = FALSE, per_row = TRUE))
  ts = ll_allelic(nA, n, f, rho, eps, per_row = TRUE)
  expect_true(all(ts >= pr))
  expect_identical(ts[k], vapply(k, function(i) ll_allelic(nA[i], n[i], f[i], rho, eps), numeric(1)))
  expect_error(ll_allelic(nA, n, f[1:2], rho, per_row = TRUE), "length 1 or length")
})

test_that("fit_depth_null reproduces the reference t fit exactly", {
  d = core_feat[null_cell == TRUE]
  expect_identical(fit_depth_null(d, n_max = 1500), ref_fit_t(d, n_max = 1500))
  expect_identical(fit_depth_null(d, by = "class", n_max = 3000, seed = 7), ref_fit_t(d, by = "class", n_max = 3000, seed = 7))
  tp = fit_depth_null(d, n_max = 1500)
  expect_true(all(tp$t_df > 2))
  expect_gt(median(tp$t_df), 4); expect_lt(median(tp$t_df), 15)   # fake null z ~ t(7)
})

test_that("fit_allelic_null fits one rho or one per group", {
  set.seed(5)
  n = rpois(4000, 40)
  p = rbeta(4000, 0.5 * (1 - 0.01) / 0.01, 0.5 * (1 - 0.01) / 0.01)
  nA = rbinom(4000, n, p)
  r = fit_allelic_null(c(nA, 0), c(n, 0))
  expect_equal(r, fit_bb_rho(nA, n))
  expect_gt(r, 0.003); expect_lt(r, 0.03)
  g = fit_allelic_null(nA, n, by = rep(c("a", "b"), 2000))
  expect_equal(names(g), c("group", "rho", "n_rows"))
  expect_equal(g[group == "a", rho], fit_bb_rho(nA[c(TRUE, FALSE)], n[c(TRUE, FALSE)]))
  g2 = fit_allelic_null(nA, n, by = data.table(stratum = rep(c("x", "y"), each = 2000)))
  expect_equal(g2$stratum, c("x", "y"))
})

test_that("depth_stratum equals the reference definition", {
  e3_depth_stratum = function(n) cut(n, c(-Inf, 2e4, 5e4, Inf), labels = c("<20k", "20-50k", ">=50k"), right = FALSE)
  n = c(0, 19999, 2e4, 20001, 49999, 5e4, 1e6, NA, runif(100, 0, 1e5))
  expect_identical(depth_stratum(n), e3_depth_stratum(n))
  expect_identical(levels(depth_stratum(1, breaks = c(1e4, 2.5e4, 5e4))), c("<10k", "10-25k", "25-50k", ">=50k"))
})

test_that("validators reject missing columns, wrong types and unknown values", {
  expect_error(check_regions(core_regions[, !"sd_class"]), "missing column\\(s\\): sd_class")
  expect_error(check_features(core_feat[, !c("z", "fold")]), "missing column\\(s\\): z, fold")
  bad = copy(core_regions)[, n_bins := as.character(n_bins)]
  expect_error(check_regions(bad), "wrong type for `n_bins`")
  bad = copy(core_regions)[1, class := "7Mb"]
  expect_error(check_regions(bad), "`class` has value\\(s\\) 7Mb")
  expect_error(check_regions(rbind(core_regions, core_regions[1])), "duplicated region")
  bad = copy(core_feat)[, stratum := factor(stratum)]
  expect_error(check_features(bad), "wrong type for `stratum`")
  expect_error(check_votes(data.table(cell = 1L, atom = 1L, ps = 10, nA = 1.5, nB = 0L)), "wrong type for `nA`")
  thr = data.table(score = "s_comb", class = "arm", stratum = "<20k", fold = 0L, q = 0.999, thr = 3.2,
                   n_null = 12000L, pooled_fallback = FALSE)
  expect_silent(check_thresholds(thr))
  calls = data.table(cell = "cell_0001", region = "chr1p", class = "arm", chrom = "chr1", start = 0, end = 24e6,
                     state = "loss", score = 12, thr = 3.2, called = TRUE, level = "cell", source = "cell", reportable = TRUE)
  expect_silent(check_calls(calls))
  expect_error(check_calls(copy(calls)[, level := "bogus"]), "`level`")
  seg = data.table(clone = "c1", chrom = "chr1", start = 0, end = 24e6, n_win = 24L, nA = 300L, nB = 5L, cn = 1L,
                   state = "1:0", log2 = -0.9, f = 0.98, llr = 150, confident = TRUE)
  expect_silent(check_segments(seg))
  expect_error(check_segments(seg[, !"confident"]), "confident")
  expect_error(check_votes(list(cell = 1L)), "expected a data.frame")
  # all-NA columns pass any type
  expect_silent(check_segments(copy(seg)[, f := NA]))
})

test_that("new_scphasedcn builds and prints a result object", {
  calls = data.table(cell = c("cell_0001", "cell_0002"), region = "chr1p", class = "arm", chrom = "chr1", start = 0,
                     end = 24e6, state = "loss", score = c(12, 1), thr = 3.2, called = c(TRUE, FALSE), level = "cell",
                     source = "cell", reportable = TRUE)
  x = new_scphasedcn("normal", cells = unique(core_feat$cell), regions = core_regions, features = core_feat, calls = calls,
                  clones = data.table(clone = c("c1", "c2")))
  expect_s3_class(x, "scphasedcn")
  expect_equal(x$params, caller_defaults())
  out = capture.output(print(x))
  expect_match(out[1], "mode: normal")
  expect_true(any(grepl("cells: 240", out)))
  expect_true(any(grepl(sprintf("regions: %d", nrow(core_regions)), out)))
  expect_true(any(grepl("calls: 2 rows, 1 called", out)))
  expect_true(any(grepl("clones: 2", out)))
  expect_error(new_scphasedcn("bogus", cells = "a", regions = core_regions))
  expect_error(new_scphasedcn("tumour", cells = "a", regions = core_regions[, !"class"]), "class")
})

test_that("tile names follow the window_ids format", {
  expect_identical(tile_name("5Mb", "chr1", 0), "5Mb|chr1:0")
  expect_identical(tile_name("50Mb", "chr10", 1e8), "50Mb|chr10:100000000")
  b = sim_genome()
  expect_identical(tile_name("10Mb", b$chrom, floor(b$start / 1e7) * 1e7), paste0("10Mb|", window_ids(b, 1e7)))
  d = caller_defaults()
  expect_equal(d$null_q, 0.999); expect_equal(unname(d$shifts), c(-0.98, 0.57))
  expect_identical(levels(depth_stratum(1)), d$strata_labels)
})

test_that("simulated fragments carry the phase set of their first SNP", {
  s = sim_small
  fr = s$fragments
  expect_true(all(fr$ps[fr$nsnp > 0] > 0))
  expect_true(all(fr$ps[fr$nsnp == 0] == -1L))
  # ps of the first covered SNP
  set.seed(9)
  k = sample(which(fr$nsnp > 0), 300)
  first_ps = vapply(k, function(i) {
    sp = s$snps[chrom == fr$chrom[i] & pos >= fr$start[i] & pos < fr$end[i]]
    sp$ps[which.min(sp$pos)]
  }, numeric(1))
  expect_equal(as.numeric(fr$ps[k]), first_ps)
  # hybrid-like phasing: one block per chromosome, anchored
  fn = sim_norm_f$fragments
  main = fn[ps > 0, .N, by = .(chrom, ps)][, .(frac = max(N) / sum(N)), by = chrom]
  expect_gt(min(main$frac), 0.5)
  expect_gt(sum(main$frac == 1), 3)
  expect_gt(fn[hap %in% 1:2, mean(hap_label(hap) == hap_true)], 0.9)
  # unoriented blocks: hap agrees with the truth within a block, up to one flip per block
  ft = sim_tum2$fragments
  agree = ft[hap %in% 1:2, .(a = mean(hap_label(hap) == hap_true), n = .N), by = .(chrom, ps)][n >= 20]
  expect_gt(nrow(agree), 100)
  expect_gt(mean(agree$a > 0.85 | agree$a < 0.15), 0.95)
  expect_gt(mean(agree$a < 0.15), 0.25)
  expect_gt(mean(agree$a > 0.85), 0.25)
  # the clone truth is carried by the event table
  expect_setequal(unique(sim_tum2$truth$events$type), c("cloneA", "cloneB"))
  expect_equal(sim_norm_f$params$sex, "female")
})

test_that("fit_depth_null: the fallback scale is never 0, and ll_depth_shift refuses a zero scale", {
  # mostly tied z: the t fit fails and the MAD is 0; the fallback scale is floored so no depth LLR is NaN
  d = data.table(class = "arm", stratum = "s", z = c(rep(0, 50), 1, -1))
  tp = fit_depth_null(d)
  expect_gt(tp$t_s, 0)
  expect_equal(tp$t_df, 30)
  expect_true(all(is.finite(ll_depth_shift(c(0, 1, -1, 3), -1, tp))))
  sc = score_cells(d[, `:=`(sd_null = 1, nA = 0, n_hap = 0)], tpar = tp, rho = 0.01)
  expect_false(anyNA(sc$s_comb))
  expect_error(ll_depth_shift(0, 1, list(t_m = 0, t_s = 0, t_df = 5)), "t_s must be > 0")
  expect_error(ll_depth_shift(c(0, 1), 1, list(t_m = 0, t_s = c(1, -1), t_df = 5)), "t_s must be > 0")
})
