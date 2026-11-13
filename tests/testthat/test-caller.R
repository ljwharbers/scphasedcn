# test-caller.R: the scorer (R/caller.R). Exact parity with the frozen reference scorer (helper-refscorer.R: scores,
# allelic LLRs, attenuation grid, cross-validated and all-null calls, pooled scores), budget calibration and recall on fake
# features, null-cell self-calibration on the sim_norm_f simulation, policy / consolidation / pools, and the
# normal-tissue wrapper end to end.

shifts = caller_defaults()$shifts

# 50k rows: arm / chrom / 20 Mb tiles x 3 strata, events on three regions, some missing z / sd_null
par_regions = fake_regions(tiles = c("20Mb" = 2e7))
par_feat = fake_features(1725, par_regions,
                         events = data.table(cell = c(1:30, 31:60, 61:90), region = rep(c("chr1p", "chr8", "chr2q"), each = 30),
                                             state = rep(c("loss", "gain", "cnloh"), each = 30)), seed = 11)
set.seed(12)
par_feat[sample(.N, 80), z := NA]
par_feat[sample(.N, 40), sd_null := NA]
par_tpar = fit_depth_null(par_feat[null_cell == TRUE], n_max = 4000)
par_rho = 0.002
par_ref = function(model = "t", att = FALSE) {
  r = cbind(par_feat, ref_allelic_llr(par_feat$nA, par_feat$n_hap, par_rho))
  if (att) r[, llr_20 := ref_llr_grid(nA, n_hap, par_rho, REF_F_ATT)]
  ref_score(r, shifts, model = model, tpar = if (model == "t") par_tpar else NULL)
}

test_that("caller_states defines the base, attenuation and tumour state sets", {
  st = caller_states("normal")
  expect_equal(st$state, c("neutral", "loss", "gain", "gain", "cnloh"))
  expect_equal(st$f, c(0.5, 0.98, 2 / 3, 3 / 4, 0.98))
  expect_equal(st$shift, c(0, -0.98, 0.57, 0.57, 0))
  expect_equal(st$row_id, c("neutral_1", "loss_1", "gain_1", "gain_2", "cnloh_1"))
  expect_identical(caller_states("e4"), st)                     # older name, accepted
  att = caller_states("normal_att")
  expect_equal(att[state == "loss", f], caller_defaults()$f_att_normal)
  expect_equal(att[state == "cnloh", f], caller_defaults()$f_att_normal)
  tu = caller_states("tumour")
  expect_setequal(unique(tu$state), c("neutral", "loss", "cnloh", "gain", "gain_as", "amp"))
  expect_equal(tu[state == "loss", f], c(0.80, caller_defaults()$f_att))
  expect_equal(tu[state == "amp", paste0(nA, ":", nB)], c("3:1", "2:2", "4:1"))
  expect_equal(tu[state == "amp", shift], 0.975 * log2(c(4, 4, 5) / 2))
  expect_equal(tu[state == "gain_as", unique(shift)], 0.57)
})

test_that("score_cells equals the reference scorer column by column (t, Gaussian, attenuation)", {
  expect_gte(nrow(par_feat), 50000)
  expect_equal(uniqueN(par_feat$stratum), 3)
  cols = c("dl_loss", "dl_gain", "al_loss", "al_gain", "al_cnloh", "cl_loss", "cl_gain", "cl_cnloh", "s_comb", "s_depth", "s_allelic")
  for (v in list(list("t", FALSE, "normal"), list("gauss", FALSE, "normal"), list("t", TRUE, "normal_att"))) {
    ref = par_ref(v[[1]], v[[2]])
    # the attenuation grid of the reference: f_att_normal = REF_F_ATT
    me = score_cells(par_feat, caller_states(v[[3]], modifyList(caller_defaults(), list(f_att_normal = REF_F_ATT))),
                     if (v[[1]] == "t") par_tpar else NULL, par_rho, depth_model = v[[1]])
    for (cn in cols) expect_lte(max(abs(me[[cn]] - ref[[cn]])), 1e-12, label = paste(v[[1]], v[[3]], cn))
    expect_identical(me$state_comb, ref$state_comb)
    expect_identical(me$state_depth, ref$state_depth)
    expect_identical(me$state_allelic, ref$state_allelic)
  }
  # per-row allelic LLRs (keep_rows) are the reference llr_20 / 21 / 31
  mk = score_cells(par_feat, caller_states("normal"), par_tpar, par_rho, keep_rows = TRUE)
  ra = ref_allelic_llr(par_feat$nA, par_feat$n_hap, par_rho)
  expect_identical(mk$ar_loss_1, ra$llr_20); expect_identical(mk$ar_gain_1, ra$llr_21); expect_identical(mk$ar_gain_2, ra$llr_31)
  expect_identical(mk$ar_cnloh_1, ra$llr_20)
  # the scorer finds the planted states
  m = mk[truth != "neutral", mean(state_comb == truth), by = truth]
  expect_true(all(m$V1 > 0.8))
})

test_that("budget_thresholds + apply_thresholds reproduce the reference CV and full calls exactly", {
  ref = par_ref("t")[, is_null := null_cell]
  me = score_cells(par_feat, caller_states("normal"), par_tpar, par_rho)[, is_null := null_cell]
  sc = c("s_depth", "s_allelic", "s_comb")
  rc = ref_cv_calls(ref, scores = sc, conc = FALSE)
  thr = budget_thresholds(me, sc, min_null = 0)
  expect_silent(check_thresholds(thr))
  mc = apply_thresholds(me, thr, sc)
  setkey(rc, cell, region); setkey(mc, cell, region)
  expect_equal(nrow(rc), nrow(mc))
  for (s in sc) {
    expect_identical(mc[[paste0("thr_", s)]], rc[[paste0("thr_", s)]])
    expect_identical(mc[[paste0("call_", s)]], rc[[paste0("call_", s)]])
  }
  # thresholds table: 1 + 4 folds per class x stratum x score; fold k trained without fold k
  expect_setequal(unique(thr$fold), 0:4)
  expect_equal(nrow(thr), uniqueN(me[, .(class, stratum)]) * 5 * length(sc))
  x = me[class == "arm" & stratum == "20-50k" & is_null == TRUE]
  expect_identical(thr[score == "s_comb" & class == "arm" & stratum == "20-50k" & fold == 2, thr], unname(quantile(x[fold != 2, s_comb], 0.999, type = 7)))
  expect_equal(thr[score == "s_comb" & class == "arm" & stratum == "20-50k" & fold == 0, thr], budget_threshold(x$s_comb, 1))
  # no CV = ref_calls (all null rows, fold 0)
  rn = ref_calls(ref, scores = sc, conc = FALSE); setkey(rn, cell, region)
  mn = apply_thresholds(me, budget_thresholds(me, sc, cv = FALSE, min_null = 0), sc); setkey(mn, cell, region)
  for (s in sc) {
    expect_identical(mn[[paste0("thr_", s)]], rn[[paste0("thr_", s)]])
    expect_identical(mn[[paste0("call_", s)]], rn[[paste0("call_", s)]])
  }
  # strict ">" on a discrete atom: the threshold sits on the atom and the atom is not called
  at = data.table(cell = paste0("c", 1:2000), class = "arm", stratum = "<20k", fold = rep(1:4, 500), is_null = TRUE,
                  s_allelic = c(rep(0, 1990), rep(5, 10)))
  ta = budget_thresholds(at, "s_allelic", min_null = 0)
  expect_true(all(ta$thr == 5))
  ca = apply_thresholds(rbind(at, data.table(cell = "x", class = "arm", stratum = "<20k", fold = 1L, is_null = FALSE, s_allelic = 5.01)), ta, "s_allelic")
  expect_false(any(ca[s_allelic == 5, call_s_allelic]))
  expect_true(ca[cell == "x", call_s_allelic])
  # pooled fallback: a sparse stratum takes the class-wide threshold
  fb = budget_thresholds(me, "s_comb", min_null = 8000, cv = FALSE)
  small = me[is_null == TRUE, .N, by = .(class, stratum)][N < 8000]
  expect_gt(nrow(small), 0)
  k = fb[class == small$class[1] & stratum == small$stratum[1]]
  expect_true(k$pooled_fallback)
  expect_equal(k$thr, unname(quantile(me[is_null == TRUE & class == small$class[1], s_comb], 0.999, type = 7)))
  # rows without a fold, or with a fold the table lacks, get fold 0
  me2 = copy(me)[1:10, fold := 9L]
  a2 = apply_thresholds(me2, thr, "s_comb")
  k0 = thr[score == "s_comb" & fold == 0][match(paste(me2$class[1:10], me2$stratum[1:10]), paste(class, stratum))]
  expect_identical(a2$thr_s_comb[1:10], k0$thr)
})

test_that("pool_scores equals the reference pool scores", {
  mk = score_cells(par_feat, caller_states("normal"), par_tpar, par_rho, keep_rows = TRUE)
  set.seed(4)
  idx = sample(nrow(mk), 6000)
  pool = rep(1:300, each = 20)
  y = mk[idx, .(dl_loss, dl_gain, al_loss, llr_21 = ar_gain_1, llr_31 = ar_gain_2)][, pool := pool]
  rp = ref_pool_scores(y)
  mp = pool_scores(mk[idx][, pool := pool], states = caller_states("normal"), by = character())
  for (cn in c("s_depth", "s_allelic", "s_comb", "state_depth", "state_allelic", "state_comb", "n")) expect_identical(mp[[cn]], rp[[cn]])
  # with a pools table (cell -> pool), per region
  pl = data.table(pool = rep(c("p1", "p2"), each = 20), cell = sprintf("cell_%04d", c(1:20, 101:120)))
  ps = pool_scores(mk, pl, caller_states("normal"))
  expect_equal(ps[pool == "p1" & region == "chr1p", state_comb], "loss")
  expect_equal(ps[pool == "p1" & region == "chr1p", n], 20L)
  expect_equal(ps[pool == "p1" & region == "chr1p", cl_loss], mk[cell %in% pl[pool == "p1", cell] & region == "chr1p", sum(dl_loss) + sum(ar_loss_1)])
  expect_error(pool_scores(score_cells(par_feat[1:10], caller_states("normal"), par_tpar, par_rho), pl), "keep_rows")
})

test_that("budget calibration: realised CV rate 0.7-1.4 per 1,000 per class; chromosome loss / CN-LOH recall > 0.8", {
  rg = fake_regions(tiles = c("5Mb" = 5e6, "10Mb" = 1e7))
  ev = data.table(cell = 1:160, region = rep(c("chr2", "chr8", "chr22", "chr1"), each = 40), state = rep(c("loss", "cnloh"), each = 80))
  ff = fake_features(2700, rg, events = ev, seed = 21)
  sc = score_cells(ff, caller_states("normal"), fit_depth_null(ff[null_cell == TRUE]), fit_allelic_null(ff[null_cell == TRUE]$nA, ff[null_cell == TRUE]$n_hap))
  thr = budget_thresholds(sc, "s_comb")
  sc = apply_thresholds(sc, thr, "s_comb")
  expect_gte(sc[null_cell == TRUE, .N], 2e5)
  rf = realised_false_calls(sc, "s_comb", by = "class")
  expect_true(all(rf$rates$per_1000 >= 0.7 & rf$rates$per_1000 <= 1.4), label = paste(round(rf$rates$per_1000, 2), collapse = " "))
  expect_lt(rf$per_cell$events_per_cell, 0.05)
  ok = sc[truth != "neutral", .(recall = mean(call_s_comb & state_comb == truth)), by = truth]
  expect_true(all(ok$recall > 0.8), label = paste(ok$truth, round(ok$recall, 2), collapse = " "))
})

# features of a simulation (stand-in for cell_features()): celltype expected counts, per-cell window z with sd_null = 1.4826 MAD over
# the cell's regions of the same class, and anchored molecule votes (hybrid-like phasing: one block per chromosome)
sim_features = function(sim, tiles = c("20Mb" = 2e7)) {
  scn = sim$scn; bins = scn$bins
  rg = fake_regions(bins, tiles = tiles)
  norm = normalise_depth(scn, "celltype")
  X = dense(scn$counts); E = norm$expected; use = norm$use
  cells = scn$cells$cell
  fr = sim$fragments[hap %in% 1:2]
  fr[, ci := if (is.numeric(cell)) as.integer(cell) else match(cell, cells)]
  out = list()
  for (i in seq_len(nrow(rg))) {
    r = rg[i]
    b = which(bins$chrom == r$chrom & bins$start >= r$start & bins$end <= r$end & use)
    h = fr[chrom == r$chrom & start >= r$start & start < r$end, .(nA = sum(hap == 1L), nB = sum(hap == 2L)), by = ci]
    nA = numeric(length(cells)); nB = numeric(length(cells)); nA[h$ci] = h$nA; nB[h$ci] = h$nB
    out[[i]] = data.table(cell = cells, region = r$region, class = r$class, sd_class = r$sd_class, chrom = r$chrom,
                          x = colSums(X[b, , drop = FALSE]), e = colSums(E[b, , drop = FALSE]), nA = nA, nB = nB)
  }
  d = rbindlist(out)
  d[, log2ratio := log2((x + 0.5) / (e + 0.5))]
  d[, sd_null := 1.4826 * mad(log2ratio, constant = 1), by = .(cell, sd_class)]
  d[, z := log2ratio / pmax(sd_null, 0.05)]
  d[, `:=`(n_hap = nA + nB, n_blocks = 1, f = fifelse(nA + nB > 0, nA / (nA + nB), NA_real_), orient = "anchored")]
  d = merge(d, scn$cells[, .(cell, n_frag = as.numeric(n_frag), cell_type)], by = "cell")
  d[, `:=`(stratum = as.character(depth_stratum(n_frag)), fold = assign_folds(cell), null_cell = NA, null_region = TRUE, eligible = TRUE)]
  d[]
}

test_that("select_null_cells drops >= 90 % of sim_norm_f's private-event carriers and keeps the rest", {
  sf = sim_features(sim_norm_f)
  tr = sim_norm_f$truth
  carriers = sim_norm_f$scn$cells$cell[unique(tr$nuclei$cell[unique(unlist(tr$carriers))])]
  mult = tr$cells[multiplet == TRUE, cell]
  sel = select_null_cells(sf, exclude = mult, chroms = NULL)          # all chromosomes: the chrX-loss carriers count too
  expect_equal(nrow(sel$log), 2)
  sc = setdiff(carriers, mult)
  expect_gt(length(sc), 20)
  expect_gte(mean(!sc %in% sel$null_cells), 0.9)
  expect_gte(mean(setdiff(unique(sf$cell), c(carriers, mult)) %in% sel$null_cells), 0.9)
  expect_false(any(mult %in% sel$null_cells))
  # the empirical variant cannot clean a set this contaminated (> 0.1 % event rows): documented limitation
  emp = select_null_cells(sf, exclude = mult, threshold = "empirical", n_iter = 1, chroms = NULL)
  expect_lt(mean(!sc %in% emp$null_cells), 0.5)
})

test_that("assign_folds is balanced and copies follow their origin", {
  f = assign_folds(paste0("c", 1:100), 4, seed = 3)
  expect_equal(as.vector(table(f)), rep(25L, 4))
  cells = c(paste0("c", 1:100), paste0("c", 1:10, "__spk"))
  g = assign_folds(cells, 4, seed = 3, origin = c(paste0("c", 1:100), paste0("c", 1:10)))
  expect_identical(g[1:100], f)
  expect_identical(g[101:110], f[1:10])
})

test_that("call_policy, call_table and consolidate_calls follow the normal-tissue rules", {
  pol = call_policy("normal")
  expect_true(pol[level == "cell" & class == "chrom" & state == "gain", reportable])
  expect_true(pol[level == "cell" & class == "arm" & state == "cnloh", reportable])
  expect_false(pol[level == "cell" & class == "arm" & state == "gain", reportable])
  expect_false(pol[level == "cell" & class == "20Mb" & state == "loss", reportable])
  expect_equal(nrow(pol[level == "cell" & class %in% c("5Mb", "10Mb")]), 0)
  expect_true(all(pol[level == "pool", reportable]))
  expect_true(all(!call_policy("tumour")[level == "cell", reportable]))
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  d = data.table(cell = "c1", region = c("chr1", "chr1p", "20Mb|chr1:0", "chr2q", "chr2", "20Mb|chr8:0"),
                 class = c("chrom", "arm", "20Mb", "arm", "chrom", "20Mb"), stratum = "<20k",
                 s_comb = c(30, 20, 10, 15, 1, 12), state_comb = c("loss", "loss", "loss", "cnloh", "cnloh", "gain"), thr_s_comb = 5)
  ct = call_table(d, NULL, rg, pol)
  expect_silent(check_calls(ct))
  expect_equal(ct[region == "chr2", called], FALSE)
  expect_equal(ct[region == "20Mb|chr8:0", reportable], FALSE)
  expect_equal(ct[region == "chr1p", c(start, end)], c(0, 24e6))
  ev = consolidate_calls(ct, rg)
  expect_equal(nrow(ev), 2)
  expect_equal(ev[chrom == "chr1", top_region], "chr1")
  expect_equal(ev[chrom == "chr1", n_regions], 2L)                  # chr1p nested in chr1; 20 Mb tile not reportable
  expect_equal(ev[chrom == "chr2", top_region], "chr2q")
  expect_equal(nrow(consolidate_calls(ct, rg, reportable_only = FALSE)), 3)
  ce = call_table(copy(d)[, eligible := region != "chr2q"], NULL, rg, pol)       # ineligible regions are not called
  expect_false("chr2q" %in% ce$region)
  thr = data.table(score = "s_comb", class = c("chrom", "arm", "20Mb"), stratum = "<20k", fold = 0L, q = 0.999, thr = 25, n_null = 1L, pooled_fallback = FALSE)
  ct2 = call_table(d[, !"thr_s_comb"], thr, rg, pol)
  expect_equal(ct2[called == TRUE, region], "chr1")
})

test_that("pools: random partitions within type, kNN neighbourhoods, pooled thresholds and calls", {
  cells = sprintf("c%03d", 1:130); type = rep(c("a", "b", "c"), c(60, 61, 9))
  p = make_pools(cells, type, n_pool = 20, seed = 2)
  expect_setequal(unique(p$group), c("a", "b"))                   # 9 cells < n_pool / 2: no pool
  expect_equal(p[, .N, by = pool][, sort(N)], c(20L, 20L, 20L, 20L, 20L, 21L))
  expect_true(all(p[, uniqueN(group), by = pool]$V1 == 1))
  expect_identical(make_pools(cells, type, 20, seed = 2), p)
  set.seed(1)
  emb = matrix(rnorm(130 * 4), 130, 4, dimnames = list(cells, NULL))
  emb[1:30, 1] = emb[1:30, 1] + 20
  kp = knn_pools(emb, k = 10)
  expect_equal(kp[, .N, by = pool][, unique(N)], 10L)
  expect_true(all(kp[center %in% cells[1:30], cell] %in% cells[1:30]))
  # pooled null thresholds on null cells with the same construction, then pooled calls find a shared loss
  rg = fake_regions(tiles = c("10Mb" = 1e7))
  ff = fake_features(400, rg, events = data.table(cell = 1:20, region = "10Mb|chr2:10000000", state = "loss"), seed = 5)
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  sc = score_cells(ff, caller_states("normal"), fit_depth_null(ff[null_cell == TRUE]), 0.002, keep_rows = TRUE)
  ct = setNames(ff[, cell_type[1], by = cell]$V1, ff[, cell_type[1], by = cell]$cell)
  fun = function(cc, seed) make_pools(cc, ct[cc], 20, seed)
  pt = pooled_thresholds(sc[null_cell == TRUE], fun, n_draws = 3000)
  expect_silent(check_thresholds(pt))
  expect_setequal(unique(pt$fold), 0:4)
  expect_true(all(pt[score == "s_comb" & fold == 0, n_null] >= 3000))
  test_pools = rbind(data.table(pool = "carriers", cell = sprintf("cell_%04d", 1:20)), make_pools(sprintf("cell_%04d", 21:400), NULL, 20, 9), fill = TRUE)
  cp = call_pools(sc, test_pools, pt, call_policy("normal"), caller_states("normal"), rg)
  expect_silent(check_calls(cp))
  expect_true(cp[cell == "carriers" & region == "10Mb|chr2:10000000", called & state == "loss"])
  expect_lt(cp[cell != "carriers", mean(called)], 0.01)
})

test_that("tumour states score attenuated losses, allele-specific gains and amplifications", {
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ff = fake_features(300, rg, events = data.table(cell = 1:20, region = "chr8", state = "loss"), seed = 7)
  st = caller_states("tumour")
  sc = score_cells(ff, st, fit_depth_null(ff[null_cell == TRUE]), 0.002, keep_rows = TRUE)
  expect_true(all(c("cl_amp", "dl_amp", "al_amp", "dr_amp_1", "dr_amp_3", "ar_gain_as_4", "s_comb") %in% names(sc)))
  expect_gt(sc[truth == "loss", mean(state_comb == "loss")], 0.8)
  # multi-shift state: the reported depth and allele terms are the best row's and add up
  expect_equal(sc$cl_amp, sc$dl_amp + sc$al_amp)
  expect_equal(sc$cl_amp, pmax(sc$dr_amp_1 + sc$ar_amp_1, sc$dr_amp_2 + sc$ar_amp_2, sc$dr_amp_3 + sc$ar_amp_3))
  ps = pool_scores(sc, data.table(pool = "p", cell = sprintf("cell_%04d", 1:20)), st)
  expect_equal(ps[region == "chr8", state_comb], "loss")
})

test_that("caller_normal runs end to end and returns a valid scphasedcn", {
  rg = fake_regions(tiles = c("5Mb" = 5e6, "10Mb" = 1e7, "20Mb" = 2e7))
  ev = data.table(cell = c(1:10, 11:20), region = rep(c("chr1p", "chr8"), each = 10), state = rep(c("loss", "gain"), each = 10))
  ff = fake_features(300, rg, events = ev, seed = 3)
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  ff[, null_cell := NA]
  set.seed(2)
  emb = matrix(rnorm(300 * 5), 300, 5, dimnames = list(unique(ff$cell), NULL))
  res = caller_normal(ff, regions = rg, embedding = emb, multiplets = "cell_0300", n_draws = 1500)
  expect_s3_class(res, "scphasedcn")
  expect_equal(res$mode, "normal")
  expect_silent(check_calls(res$calls))
  expect_silent(check_thresholds(res$thresholds))
  expect_setequal(unique(res$calls$level), c("cell", "pool"))
  expect_setequal(unique(res$thresholds$level), c("cell", "pool_type", "pool_knn"))
  expect_false("cell_0300" %in% res$null_fit$null_cells)
  expect_false(any(sprintf("cell_%04d", 1:10) %in% res$null_fit$null_cells))
  # the planted chr1p losses are reportable per-cell events
  e = res$events[level == "cell" & chrom == "chr1" & state == "loss"]
  expect_gte(sum(sprintf("cell_%04d", 1:10) %in% e$cell), 8)
  expect_true(all(c("tpar", "rho", "null_cells", "selection") %in% names(res$null_fit)))
  out = capture.output(print(res))
  expect_match(out[1], "mode: normal")
  # a null_cell column is used as given ("features"); "auto" selects when it is all TRUE (no information)
  ff2 = fake_features(200, rg[class %in% c("arm", "chrom")], seed = 4)
  r2 = caller_normal(ff2, regions = rg[class %in% c("arm", "chrom")], n_draws = 500, null_cells = "features")
  expect_equal(length(r2$null_fit$null_cells), 200)
  expect_null(r2$null_fit$selection)
  r3 = caller_normal(ff2, regions = rg[class %in% c("arm", "chrom")], n_draws = 500)
  expect_false(is.null(r3$null_fit$selection))
  expect_gte(length(r3$null_fit$null_cells), 190)
})

# ---- calibration fixes (nested gains, pooled centring / leave-one-out, chrX policy) ------------------------------

test_that("resolve_nested: sim_norm_f carriers of the 1p loss / chr2 CN-LOH get no chromosome gain, true events stay", {
  sf = sim_features(sim_norm_f)
  rg = fake_regions(sim_norm_f$scn$bins, tiles = c("20Mb" = 2e7))
  tr = sim_norm_f$truth
  cells = sim_norm_f$scn$cells$cell
  car = lapply(tr$carriers, function(n) unique(cells[tr$nuclei$cell[n]]))     # 1p loss, chr8 gain, chr2 CN-LOH, chr21 gain, chrX loss
  mult = tr$cells[multiplet == TRUE, cell]
  r0 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, resolve_gains = FALSE, pool_classes = character())
  r1 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, pool_classes = character())
  rep_cells = function(r, reg, st) unique(r$calls[level == "cell" & called & reportable & region == reg & state == st, cell])
  # without resolution the diluted nested events are reported as chromosome gains
  expect_gte(length(intersect(rep_cells(r0, "chr1", "gain"), car[[1]])), 5)
  expect_gte(length(intersect(rep_cells(r0, "chr2", "gain"), car[[3]])), 4)
  # with resolution: no chromosome gain on carriers of the nested events; the nested events and the true gains are called
  expect_length(intersect(rep_cells(r1, "chr1", "gain"), car[[1]]), 0)
  expect_length(intersect(rep_cells(r1, "chr2", "gain"), car[[3]]), 0)
  expect_gte(mean(car[[1]] %in% rep_cells(r1, "chr1p", "loss")), 0.8)
  expect_gte(mean(car[[3]] %in% rep_cells(r1, "chr2q", "cnloh")), 0.8)
  expect_gte(mean(car[[2]] %in% rep_cells(r1, "chr8", "gain")), 0.7)
  expect_gte(mean(car[[4]] %in% rep_cells(r1, "chr21", "gain")), 0.7)
  ev = r1$events[level == "cell"]
  expect_equal(nrow(ev[cell %in% c(car[[1]], car[[3]]) & chrom %in% c("chr1", "chr2") & state == "gain"]), 0)
  x = r1$calls[!is.na(resolved)]
  expect_true(all(x[state == "gain", called & !reportable]))
  expect_true(all(grepl("^no_depth|^nested ", x[state == "gain", resolved])))
  # every dropped gain lacks depth support or lost to a stronger nested loss / CN-LOH
  g = merge(x[state == "gain", .(cell, region, resolved)], r1$scores[, .(cell, region, dl_gain)], by = c("cell", "region"))
  expect_true(all(g$dl_gain <= 0 | grepl("nested", g$resolved)))
})

test_that("resolve_nested rules on a hand-made calls table", {
  ct = data.table(cell = "c1", region = c("chr1", "chr1p", "chr2", "chr2q", "20Mb|chr2:20000000", "chr3"), class = c("chrom", "arm", "chrom", "arm", "20Mb", "chrom"),
                  chrom = c("chr1", "chr1", "chr2", "chr2", "chr2", "chr3"), start = NA_real_, end = NA_real_,
                  state = c("gain", "loss", "gain", "gain", "cnloh", "gain"), score = c(20, 30, 25, 9, 12, 15), thr = 5, called = TRUE,
                  level = "cell", source = "cell", reportable = c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE), dl_gain = c(2, -1, 3, 1, 0, -0.5))
  r = resolve_nested(ct)
  expect_false(r[region == "chr1", reportable])                   # nested 1p loss scores higher
  expect_match(r[region == "chr1", resolved], "nested chr1p loss")
  expect_true(r[region == "chr2", reportable])                    # depth support, nested CN-LOH weaker
  expect_false(r[region == "chr3", reportable])                   # no depth support
  expect_equal(r[region == "chr3", resolved], "no_depth")
  expect_equal(r[region == "chr1p", resolved], "replaces chr1 gain")
  expect_true(is.na(r[region == "20Mb|chr2:20000000", resolved]))
  r2 = resolve_nested(ct, min_dl_gain = 3.5)                      # chr2 then lacks depth support: the 20 Mb CN-LOH is promoted
  expect_false(r2[region == "chr2", reportable])
  expect_true(r2[region == "20Mb|chr2:20000000", reportable])
  expect_equal(nrow(resolve_nested(ct[0])), 0)
  expect_error(resolve_nested(ct[, !"dl_gain"]), "dl_gain")
})

# fake features: a region-biased null (every cell's z at one 10 Mb tile shifted down), a molecule-poor tile with one strong
# single carrier, and optionally a pool-wide chr8 gain in `gain_cells`
pool_fake = function(n = 400, seed = 31, gain_cells = integer()) {
  rg = fake_regions(tiles = c("10Mb" = 1e7))
  ev = rbind(data.table(cell = 1L, region = "10Mb|chr2:0", state = "loss"),
             if (length(gain_cells)) data.table(cell = gain_cells, region = "chr8", state = "gain"))
  ff = fake_features(n, rg, events = ev, seed = seed)
  ff[region == "10Mb|chr8:10000000", z := z - 2]                  # region bias shared by all cells (null and test)
  set.seed(seed + 1)
  ff[region == "10Mb|chr2:0" & cell != "cell_0001", n_hap := as.numeric(rpois(.N, 1) + 1)]
  ff[region == "10Mb|chr2:0" & cell != "cell_0001", nA := as.numeric(rbinom(.N, n_hap, 0.5))]
  ff[cell == "cell_0001" & region == "10Mb|chr2:0", `:=`(n_hap = 400, nA = 392)]
  ff[, nB := n_hap - nA]
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  ff[]
}
null_pool_set = function(cells, ct, seeds) rbindlist(lapply(seeds, function(k) make_pools(cells, ct[cells], 20, k)[, pool := paste0(k, "|", pool)]))

test_that("pool_scale / pool_scores: region standardisation, LOO = dropping the top cell", {
  ff = pool_fake()
  st = caller_states("normal")
  sc = score_cells(ff, st, fit_depth_null(ff[null_cell == TRUE]), 0.002, keep_rows = TRUE)
  nul = sc[null_cell == TRUE]
  # region standardisation: the molecule-poor tile (its non-carriers add little negative LLR) stops collecting the class's null calls
  ct = setNames(ff[, cell_type[1], by = cell]$V1, ff[, cell_type[1], by = cell]$cell)
  fun = function(cc, seed) make_pools(cc, ct[cc], 20, seed)
  p0 = pooled_thresholds(nul, fun, n_draws = 4000, cv = FALSE)
  p1 = pooled_thresholds(nul, fun, n_draws = 4000, cv = FALSE, loo = TRUE, standardize = TRUE, min_region_pools = 300,
                         params = modifyList(caller_defaults(), list(pool_tail_q = 0)))   # plain med / mad (no tail)
  sca = attr(p1, "pool_scale")
  expect_true(all(c("region", "class", "n_pools", "med", "mad") %in% names(sca)))
  expect_gte(min(sca$n_pools), 300)
  expect_gt(sca[region == "10Mb|chr2:0", med], sca[region == "10Mb|chr8:0", med] + 20)
  expect_true("s_rob" %in% p1$score)
  test = null_pool_set(unique(nul$cell), ct, 101:220)
  c0 = call_pools(sc, test, p0, call_policy("normal"), st)
  c1 = call_pools(sc, test, p1, call_policy("normal"), st, score = "s_rob", scale = sca)
  expect_true(all(c("s_comb", "s_loo", "s_comb_llr") %in% names(c1)))
  share = function(x) x[class == "10Mb" & called == TRUE, max(table(region)) / .N]
  expect_gt(c0[class == "10Mb" & called == TRUE, mean(region == "10Mb|chr2:0")], 0.5)   # unstandardised: one tile collects the calls
  expect_lt(share(c1), 0.35)                                      # standardised: spread over the 22 tiles
  expect_lt(c1[, mean(called)], 0.004)
  # leave-one-out score of one pool = best row of the called state after dropping that row's top cell
  pl = data.table(pool = "p", cell = sprintf("cell_%04d", 1:20))
  ps = pool_scores(sc, pl, st, loo = TRUE)
  x = sc[cell %in% pl$cell & region == "10Mb|chr2:0"]
  rc = x$dl_loss + x$ar_loss_1
  y = ps[region == "10Mb|chr2:0"]
  expect_equal(y$cl_loss, sum(rc), tolerance = 1e-9)
  expect_equal(y$lo_loss, sum(rc) - max(rc), tolerance = 1e-9)
  expect_equal(y$s_loo, y[[paste0("lo_", y$state_comb)]])
  expect_equal(y$s_rob, min(y$s_comb, y$s_loo))
  # standardised scores: (cl - med) / mad of the region
  pz = pool_scores(sc, pl, st, loo = TRUE, scale = sca)
  praw = pool_scores(sc, pl, st, loo = TRUE)
  k = match(pz$region, sca$region)
  expect_equal(pz$s_comb, (praw$s_comb - sca$med[k]) / sca$mad[k])
  expect_equal(pz$s_loo, (praw$s_loo - sca$med[k]) / sca$mad[k])
  expect_identical(pz$state_comb, praw$state_comb)                # the state stays the likelihood argmax
  expect_equal(pz$s_comb_llr, praw$s_comb)
  # defaults unchanged: no LOO / standardised columns
  expect_false(any(c("s_loo", "s_rob", "s_comb_llr") %in% names(pool_scores(sc, pl, st))))
})

test_that("leave-one-out: one strong carrier does not make its ~20 overlapping kNN pools call; a pool-wide gain does", {
  cells = sprintf("cell_%04d", 1:400)
  set.seed(5)
  emb = matrix(rnorm(400 * 4), 400, 4, dimnames = list(cells, NULL))
  kp = knn_pools(emb, 20)
  pw = kp[pool == "knn|cell_0200", cell]
  ff = pool_fake(gain_cells = match(pw, cells))
  st = caller_states("normal")
  sc = score_cells(ff, st, fit_depth_null(ff[null_cell == TRUE]), 0.002, keep_rows = TRUE)
  fun = function(cc, seed) knn_pools(emb[cc, , drop = FALSE], 20)
  nul = sc[null_cell == TRUE]
  t0 = pooled_thresholds(nul, fun, cv = FALSE, standardize = TRUE)
  t1 = pooled_thresholds(nul, fun, cv = FALSE, loo = TRUE, standardize = TRUE)
  carrier_pools = setdiff(kp[cell == "cell_0001", unique(pool)], kp[cell %in% pw, unique(pool)])
  expect_gte(length(carrier_pools), 10)
  c0 = call_pools(sc, kp, t0, call_policy("normal"), st, scale = attr(t0, "pool_scale"))
  c1 = call_pools(sc, kp, t1, call_policy("normal"), st, score = "s_rob", scale = attr(t1, "pool_scale"))
  r = "10Mb|chr2:0"
  expect_gte(c0[region == r & cell %in% carrier_pools, mean(called)], 0.8)          # (CN-LOH: the pool's depth does not support a loss)
  expect_lte(c1[region == r & cell %in% carrier_pools, mean(called)], 0.1)
  expect_true(c1[cell == "knn|cell_0200" & region == "chr8", called & state == "gain"])
})

test_that("chrX: pools are not reportable in females, cell calls are; select_null_cells scores autosomes by default", {
  pol = call_policy("normal")
  expect_true(all(c("class", "state", "level", "reportable", "chrx") %in% names(pol)))
  expect_true(all(pol[level == "cell", chrx]))
  expect_false(any(pol[level == "pool", chrx]))
  expect_true(all(call_policy("normal", "male")$chrx))
  expect_true(all(call_policy("tumour")$chrx))
  rg = fake_regions(tiles = c("10Mb" = 1e7))
  d = data.table(cell = "c", region = c("chrX", "chr8"), class = "chrom", stratum = "<20k", s_comb = 10, state_comb = "gain", thr_s_comb = 1)
  expect_true(all(call_table(d, NULL, rg, pol)$reportable))
  # pooled rows on chrX: called but not reportable (female / unknown), reportable under a male policy
  ff = fake_features(60, rg[class %in% c("chrom", "arm")], seed = 8)
  ff[region %in% c("chrX", "chr8"), z := z + 8]
  st = caller_states("normal")
  sc = score_cells(ff, st, NULL, 0.002, depth_model = "gauss", keep_rows = TRUE)
  thr = data.table(score = "s_comb", class = c("chrom", "arm"), stratum = "all", fold = 0L, q = 0.999, thr = 20, n_null = 1L, pooled_fallback = FALSE)
  pl = data.table(pool = "p", cell = unique(ff$cell)[1:20])
  cf = call_pools(sc, pl, thr, call_policy("normal", "female"), st, rg)
  expect_true(all(cf[region %in% c("chrX", "chr8"), called]))
  expect_false(cf[region == "chrX", reportable])
  expect_true(cf[region == "chr8", reportable])
  expect_true(call_pools(sc, pl, thr, call_policy("normal", "male"), st, rg)[region == "chrX", reportable])
  # null-cell selection: a cell whose only event is on chrX stays a null cell by default
  ev = data.table(cell = c(1:15, 16:30), region = rep(c("chrX", "chr2"), each = 15), state = "loss")
  fs = fake_features(600, rg[class %in% c("chrom", "arm")], events = ev, seed = 9)
  fs[, null_cell := NA]
  s_auto = select_null_cells(fs)
  s_all = select_null_cells(fs, chroms = NULL)
  xc = sprintf("cell_%04d", 1:15); ac = sprintf("cell_%04d", 16:30)
  expect_gte(mean(xc %in% s_auto$null_cells), 0.8)
  expect_lte(mean(xc %in% s_all$null_cells), 0.2)
  expect_lte(mean(ac %in% s_auto$null_cells), 0.2)
  expect_error(select_null_cells(fs, chroms = "chrZ"), "no rows")
  # caller_normal infers the sex from chrX allele votes
  expect_equal(infer_sex(data.table(chrom = "chrX", n_hap = c(10, 20), allelic_ok = TRUE)), "female")
  expect_equal(infer_sex(data.table(chrom = "chrX", n_hap = 0, allelic_ok = FALSE)), "male")
  expect_null(infer_sex(data.table(chrom = "chr1", n_hap = 1)))
})

# ---- calibration options and the pooled tail -------------------------------------------------------------------

test_that("caller_normal checks its calibration options", {
  sf = sim_features(sim_norm_f)
  rg = fake_regions(sim_norm_f$scn$bins, tiles = c("20Mb" = 2e7))
  p = caller_defaults()
  expect_equal(p$calib, "class")
  expect_error(caller_normal(sf, regions = rg, n_draws = 500, pool_classes = character(), params = modifyList(p, list(calib = "x"))), "calib")
  expect_error(caller_normal(sf, regions = rg, n_draws = 500, pool_classes = character(), params = modifyList(p, list(calib_rho = "x"))), "calib_rho")
  expect_false(p$pool_events)
  for (v in list("yes", NA, c(TRUE, FALSE), 1))
    expect_error(caller_normal(sf, regions = rg, n_draws = 500, pool_classes = character(), params = modifyList(p, list(pool_events = v))),
                 "pool_events must be TRUE or FALSE")
})

test_that("pool_scale tail_q: a heavier-tailed region gets tail > 1, homogeneous regions ~1; pool_combine divides by it", {
  set.seed(21)
  st = caller_states("normal")
  n = 3000
  ps = rbindlist(lapply(1:20, function(i) {
    s = if (i == 1) stats::rt(n, 3) else rnorm(n)
    data.table(pool = seq_len(n), region = sprintf("10Mb|chr%d:0", i), class = "10Mb", cl_loss = s, cl_gain = -Inf, cl_cnloh = -Inf)
  }))
  s0 = pool_scale(ps, st)
  expect_identical(names(s0), c("region", "class", "n_pools", "med", "mad"))   # default columns unchanged
  expect_identical(pool_scale(ps, st, tail_q = 0), s0)
  s1 = pool_scale(ps, st, tail_q = 0.99)
  expect_true(all(c("region", "class", "n_pools", "med", "mad", "tail") %in% names(s1)))
  expect_equal(s1[, .(region, class, n_pools, med, mad)], s0)
  heavy = s1[region == "10Mb|chr1:0", tail]
  expect_gt(heavy, 1.3)
  expect_true(all(abs(s1[region != "10Mb|chr1:0", tail] - 1) < 0.1))
  # shrinkage: n_eff caps the weight, tail_k = Inf switches it off
  s2 = pool_scale(ps, st, tail_q = 0.99, n_eff = 50)
  expect_lt(s2[region == "10Mb|chr1:0", tail], heavy)
  expect_gt(s2[region == "10Mb|chr1:0", tail], 1)
  w = 50 / 150
  u_c = median(s2[n_pools >= 30, u_tail])
  expect_equal(s2[region == "10Mb|chr1:0", tail], exp(w * log(s2[region == "10Mb|chr1:0", u_tail] / u_c)))
  expect_true(all(pool_scale(ps, st, tail_q = 0.99, tail_k = Inf)$tail == 1))
  # pool_combine: standardised by mad x tail; a region without null pools takes the class medians and tail 1
  P = data.table(pool = 1:3, region = c("10Mb|chr1:0", "10Mb|chr2:0", "10Mb|chr99:0"), class = "10Mb", cl_loss = 10, cl_gain = 0, cl_cnloh = 0)
  c0 = scphasedcn:::pool_combine(copy(P), c("loss", "gain", "cnloh"), FALSE, s0)
  c1 = scphasedcn:::pool_combine(copy(P), c("loss", "gain", "cnloh"), FALSE, s1)
  k = match(P$region[1:2], s1$region)
  expect_equal(c1$s_comb[1:2], (10 - s1$med[k]) / (s1$mad[k] * s1$tail[k]))
  expect_equal(c1$s_comb[3], (10 - median(s1$med)) / median(s1$mad))
  expect_equal(c1$s_comb[3], c0$s_comb[3])
  expect_lt(c1$s_comb[1], c0$s_comb[1])
  # pooled_thresholds passes params$pool_tail_q: the scale attribute carries the tail column
  ff = pool_fake(n = 200)
  sc = score_cells(ff, st, fit_depth_null(ff[null_cell == TRUE]), 0.002, keep_rows = TRUE)
  ct = setNames(ff[, cell_type[1], by = cell]$V1, ff[, cell_type[1], by = cell]$cell)
  fun = function(cc, seed) make_pools(cc, ct[cc], 20, seed)
  pt = pooled_thresholds(sc[null_cell == TRUE], fun, n_draws = 500, cv = FALSE, loo = TRUE, standardize = TRUE, min_region_pools = 100,
                         params = modifyList(caller_defaults(), list(pool_tail_q = 0.99)))
  expect_true("tail" %in% names(attr(pt, "pool_scale")))
  expect_true(all(is.finite(attr(pt, "pool_scale")$tail)))
})

# ---- vote strata (calib = "votes") and the cross-fitted pooled tail ------------------------------------------------

# one class / stratum: 5 vote-poor regions (n_hap ~100) and 15 vote-rich ones (~1,200), null allele votes and z
vote_fake = function(n_cells = 3000, seed = 3, rho = 0.002) {
  set.seed(seed)
  d = CJ(ci = seq_len(n_cells), ri = 1:20)
  d[, `:=`(cell = sprintf("c%05d", ci), region = sprintf("r%02d", ri), class = "chrom", stratum = "20-50k",
           fold = as.integer((ci - 1L) %% 4L + 1L), sd_null = 0.2, is_null = TRUE)]
  d[, n_hap := as.numeric(rpois(.N, fifelse(ri <= 5, 100, 1200)))]
  p = rbeta(nrow(d), 0.5 * (1 - rho) / rho, 0.5 * (1 - rho) / rho)
  d[, nA := as.numeric(rbinom(.N, n_hap, p))]
  d[, z := rnorm(.N)]
  d[, c("ci", "ri") := NULL]
  d[]
}

test_that("vote_bins / assign_hstratum: quantile bins per class x stratum, adaptive count, none, reproducible", {
  d = vote_fake(n_cells = 1000)
  d[1:7, n_hap := 0]
  d[8, n_hap := NA]
  expect_equal(caller_defaults()$vote_bins, 4)
  vb = vote_bins(d, 4, 2000)
  expect_true(all(c("class", "stratum", "hstratum", "lo", "hi", "n_null") %in% names(vb)))
  expect_equal(vb$hstratum, paste0("h", 1:4))
  expect_equal(vb$lo[-1], vb$hi[-4])
  expect_equal(vb$hi[1:3], unname(quantile(d[n_hap > 0, n_hap], 1:3 / 4, type = 7)))
  expect_equal(vb$n_null[1], d[n_hap > 0, .N])
  h = assign_hstratum(d, vb)
  expect_true(all(h[1:8] == "none"))
  tab = table(h[-(1:8)])
  expect_true(all(abs(tab / sum(tab) - 0.25) < 0.02))
  i = which(h == "h2")
  expect_true(all(d$n_hap[i] > vb$lo[2] & d$n_hap[i] <= vb$hi[2]))
  expect_identical(assign_hstratum(copy(d)[, null_cell := FALSE], vb), h)        # depends on n_hap only
  # fewer bins when a bin would have < min_null null rows; one bin = "h1" for everything with votes
  expect_equal(nrow(vote_bins(d, 4, 6000)), 3L)
  expect_equal(nrow(vote_bins(d, 4, 1e6)), 1L)
  expect_true(all(assign_hstratum(d, vote_bins(d, 4, 1e6))[-(1:8)] == "h1"))
  # a class x stratum without bins: "all"
  d2 = copy(d)[9:13][, stratum := ">=50k"]
  expect_true(all(assign_hstratum(d2, vb) == "all"))
  expect_equal(nrow(vote_bins(d[0])), 0L)
})

test_that("budget_thresholds: hstratum keys the thresholds, falls back to class x stratum then class; default unchanged", {
  d = vote_fake()
  st = caller_states("normal")
  sc = score_cells(d, st, fit_depth_null(d), 0.002)
  t0 = budget_thresholds(sc, "s_comb", min_null = 5000)
  expect_false("hstratum" %in% names(t0))
  # a constant vote stratum gives the class x stratum thresholds
  tc = budget_thresholds(copy(sc)[, hstratum := "h1"], "s_comb", min_null = 5000)
  expect_true("hstratum" %in% names(tc))
  expect_identical(tc$thr, t0$thr)
  expect_identical(apply_thresholds(copy(sc)[, hstratum := "h1"], tc, "s_comb")$call_s_comb, apply_thresholds(sc, t0, "s_comb")$call_s_comb)
  # vote strata: per-region false calls back near the budget (the class threshold is set by the vote-rich regions)
  vb = vote_bins(sc, 4, 5000)
  sh = copy(sc)[, hstratum := assign_hstratum(sc, vb)]
  t1 = budget_thresholds(sh, "s_comb", min_null = 5000)
  expect_setequal(unique(t1$hstratum), paste0("h", 1:4))
  expect_false(any(t1$pooled_fallback))
  r0 = realised_false_calls(apply_thresholds(sc, t0, "s_comb"), "s_comb", by = "region")$rates
  r1 = realised_false_calls(apply_thresholds(sh, t1, "s_comb"), "s_comb", by = "region")$rates
  expect_gt(r0[region <= "r05", sum(n_called) / sum(n_null)] * 1000, 3)
  expect_lt(r0[region > "r05", sum(n_called) / sum(n_null)] * 1000, 0.3)
  expect_lt(abs(r1[region <= "r05", sum(n_called) / sum(n_null)] * 1000 - 1), 0.6)
  expect_lt(abs(r1[region > "r05", sum(n_called) / sum(n_null)] * 1000 - 1), 0.6)
  expect_error(apply_thresholds(sc, t1, "s_comb"), "hstratum")
  # fallback: a small vote stratum takes its class x stratum; if that is small too, the class
  fb = copy(sh)[region == "r01" & fold == 1, hstratum := "hx"]
  tf = budget_thresholds(fb, "s_comb", min_null = 1000, cv = FALSE)
  k = tf[hstratum == "hx"]
  expect_true(k$pooled_fallback)
  expect_equal(k$thr, unname(quantile(fb[stratum == "20-50k", s_comb], 0.999, type = 7)))
  fb2 = rbind(fb, copy(fb[1:600])[, `:=`(stratum = "<20k", hstratum = "h1")])
  tf2 = budget_thresholds(fb2, "s_comb", min_null = 1000, cv = FALSE)
  expect_true(tf2[stratum == "<20k", pooled_fallback])
  expect_equal(tf2[stratum == "<20k", thr], unname(quantile(fb2$s_comb, 0.999, type = 7)))
  expect_false(tf2[stratum == "20-50k" & hstratum == "h2", pooled_fallback])
})

test_that("caller_normal: calib = 'votes' on sim_norm_f; the default path has no vote strata", {
  sf = sim_features(sim_norm_f)
  rg = fake_regions(sim_norm_f$scn$bins, tiles = c("20Mb" = 2e7))
  mult = sim_norm_f$truth$cells[multiplet == TRUE, cell]
  p = modifyList(caller_defaults(), list(min_null_rows = 200))
  r0 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, pool_classes = character(), params = p)
  expect_false("hstratum" %in% names(r0$scores))
  expect_false("hstratum" %in% names(r0$thresholds))
  expect_null(r0$null_fit$vote_bins)
  r1 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, pool_classes = character(), params = modifyList(p, list(calib = "votes")))
  vb = r1$null_fit$vote_bins
  expect_false(is.null(vb))
  expect_true(all(c("hstratum") %in% names(r1$scores)))
  expect_true(all(r1$scores$hstratum %in% c(unique(vb$hstratum), "none", "all")))
  expect_gt(uniqueN(r1$thresholds[level == "cell", hstratum]), 1)
  expect_silent(check_thresholds(r1$thresholds))
  expect_true(any(grepl("vote stratum", r1$notes)))
  expect_identical(r1$scores$z, r0$scores$z)                                    # votes leave the depth z alone
  # the planted chromosome gains are still found
  car = lapply(sim_norm_f$truth$carriers, function(n) unique(sim_norm_f$scn$cells$cell[sim_norm_f$truth$nuclei$cell[n]]))
  rep_cells = function(r, reg, st) unique(r$calls[level == "cell" & called & reportable & region == reg & state == st, cell])
  expect_gte(mean(car[[2]] %in% rep_cells(r1, "chr8", "gain")), 0.6)
})

test_that("cross-fitted pooled tail: the tail and the threshold come from disjoint pools; null stays in range", {
  ff = pool_fake(n = 400)
  st = caller_states("normal")
  sc = score_cells(ff, st, fit_depth_null(ff[null_cell == TRUE]), 0.002, keep_rows = TRUE)
  nul = sc[null_cell == TRUE]
  ct = setNames(ff[, cell_type[1], by = cell]$V1, ff[, cell_type[1], by = cell]$cell)
  fun = function(cc, seed) make_pools(cc, ct[cc], 20, seed)
  test = null_pool_set(unique(nul$cell), ct, 1001:1200)
  rate = function(p) {
    pt = pooled_thresholds(nul, fun, n_draws = 4000, cv = FALSE, loo = TRUE, standardize = TRUE, min_region_pools = 300, params = p)
    cp = call_pools(sc, test, pt, call_policy("normal"), st, score = "s_rob", scale = attr(pt, "pool_scale"))
    list(pt = pt, r = 1000 * mean(cp$called))
  }
  p = caller_defaults()
  expect_true(p$pool_tail_xfit)
  a0 = rate(modifyList(p, list(pool_tail_q = 0)))                                # no tail: xfit is not used
  ax = rate(modifyList(p, list(pool_tail_q = 0.99)))
  an = rate(modifyList(p, list(pool_tail_q = 0.99, pool_tail_xfit = FALSE)))
  expect_false("tail" %in% names(attr(a0$pt, "pool_scale")))
  sx = attr(ax$pt, "pool_scale")
  expect_true(all(c("region", "class", "n_pools", "med", "mad", "u_tail", "tail") %in% names(sx)))
  expect_true(all(is.finite(sx$tail) & sx$tail > 0))
  expect_gte(ax$r, 0.1)
  expect_lte(ax$r, 2)
  expect_gte(ax$pt[score == "s_rob", min(thr)], an$pt[score == "s_rob", min(thr)])   # held-out scale: not optimistic
  # caller_normal on sim_norm_f with the cross-fitted tail: pooled null of null-cell pools on the autosomes at or below
  # the budget (chrX is not a null region: female X-inactivation imbalance, not reportable in pools). With 200 cells
  # (~10 pools per draw, halves of ~90 cells) this is a smoke test, not a calibration measurement.
  sf = sim_features(sim_norm_f)
  rg = fake_regions(sim_norm_f$scn$bins, tiles = c("20Mb" = 2e7))
  mult = sim_norm_f$truth$cells[multiplet == TRUE, cell]
  r = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 2000, params = modifyList(p, list(pool_tail_q = 0.99)))
  ps = r$null_fit$pool_scale
  expect_true("tail" %in% names(ps))
  expect_true(any(grepl("cross-fitted TRUE", r$notes)))
  nc = r$null_fit$null_cells
  ctn = setNames(r$cells$cell_type, r$cells$cell)
  tp = null_pool_set(nc, ctn, 2001:2100)
  cn = call_pools(r$scores[cell %in% nc], tp, r$thresholds[level == "pool_type"], call_policy("normal"), caller_states("normal"),
                  score = "s_rob", scale = ps[kind == "type"])
  expect_gt(nrow(cn[chrom != "chrX"]), 10000)
  expect_lte(1000 * cn[chrom != "chrX", mean(called)], 3)
  expect_false(any(cn[chrom == "chrX", reportable]))
})

# ---- per-chromosome rho (calib_rho = "chrom") -----------------------------------------------------------------------

# one class / stratum, 20 regions whose mean n_hap runs log-evenly from 20 to 2,000; region r01..r20 on chromosomes
# chr1..chr20; chr7 with a larger allelic overdispersion (rho7)
cont_fake = function(n_cells = 3000, seed = 5, rho = 0.002, rho7 = rho, class = "chrom") {
  set.seed(seed)
  d = CJ(ci = seq_len(n_cells), ri = 1:20)
  d[, `:=`(cell = sprintf("c%05d", ci), region = sprintf("r%02d", ri), chrom = paste0("chr", ri), class = class, stratum = "20-50k",
           fold = as.integer((ci - 1L) %% 4L + 1L), sd_null = 0.2, is_null = TRUE)]
  d[, n_hap := as.numeric(rpois(.N, exp(log(20) + (ri - 1) / 19 * log(100))))]
  r = fifelse(d$ri == 7L, rho7, rho)
  p = rbeta(nrow(d), 0.5 * (1 - r) / r, 0.5 * (1 - r) / r)
  d[, nA := as.numeric(rbinom(.N, n_hap, p))]
  d[, z := rnorm(.N)]
  d[, c("ci", "ri") := NULL]
  d[]
}

test_that("score_cells with a per-row rho scores each rho group as with a scalar rho; equal values = scalar", {
  d = cont_fake(n_cells = 300)
  st = caller_states("normal")
  tp = fit_depth_null(d)
  s0 = score_cells(d, st, tp, 0.002, keep_rows = TRUE)
  expect_identical(score_cells(d, st, tp, rep(0.002, nrow(d)), keep_rows = TRUE), s0)
  rr = fifelse(d$chrom == "chr7", 0.006, 0.002)
  s1 = score_cells(d, st, tp, rr, keep_rows = TRUE)
  i = which(d$chrom == "chr7")
  s7 = score_cells(d[i], st, tp, 0.006, keep_rows = TRUE)
  expect_identical(s1[i], s7)
  expect_identical(s1[-i], s0[-i])
  expect_error(score_cells(d, st, tp, c(0.002, 0.003)), "one value or one per row")
})

test_that("rho_by_chrom recovers an injected per-chromosome rho and shrinks it to the global", {
  d = cont_fake(n_cells = 3000, rho = 0.002, rho7 = 0.008, class = "arm")
  expect_equal(caller_defaults()$calib_rho, "global")
  expect_equal(caller_defaults()$rho_chrom_k, 100)
  r0 = rho_by_chrom(d, k = 0)
  expect_true(all(c("chrom", "n_rows", "rho_raw", "rho_global", "w", "rho") %in% names(r0)))
  expect_equal(nrow(r0), 20L)
  expect_equal(r0$rho, r0$rho_raw)
  expect_equal(r0[chrom == "chr7", rho], 0.008, tolerance = 0.25)
  expect_true(all(abs(r0[chrom %in% paste0("chr", 10:20), rho] / 0.002 - 1) < 0.5))   # vote-rich chromosomes (n_hap > 130)
  r1 = rho_by_chrom(d, k = 1000)
  x = r1[chrom == "chr7"]
  expect_equal(x$rho, (x$n_rows * x$rho_raw + 1000 * x$rho_global) / (x$n_rows + 1000))
  expect_equal(x$w, x$n_rows / (x$n_rows + 1000))
  expect_lt(x$rho, x$rho_raw)
  expect_gt(x$rho, 2 * x$rho_global)
  # per-row rho: the global where a chromosome has no null rows
  dd = rbind(d[1:3], data.table(chrom = "chrX"), fill = TRUE)
  expect_equal(rho_rows(dd, r1, 0.0015), c(r1$rho[match(d$chrom[1:3], r1$chrom)], 0.0015))
})

test_that("caller_normal: calib = 'votes' with calib_rho = 'chrom' on sim_norm_f; the default records neither", {
  sf = sim_features(sim_norm_f)
  rg = fake_regions(sim_norm_f$scn$bins, tiles = c("20Mb" = 2e7))
  mult = sim_norm_f$truth$cells[multiplet == TRUE, cell]
  p = modifyList(caller_defaults(), list(min_null_rows = 100))
  r0 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, pool_classes = character(), params = p)
  expect_null(r0$null_fit$rho_chrom)
  expect_null(r0$null_fit$vote_bins)
  r1 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, pool_classes = character(),
                     params = modifyList(p, list(calib = "votes", calib_rho = "chrom")))
  rc = r1$null_fit$rho_chrom
  expect_false(is.null(rc))
  expect_true(all(rc$chrom %in% paste0("chr", 1:22)))
  expect_false(is.null(r1$null_fit$vote_bins))
  expect_silent(check_thresholds(r1$thresholds))
  expect_false(anyNA(r1$scores$thr_s_comb))
  expect_true(any(grepl("per chromosome", r1$notes)) && any(grepl("vote stratum", r1$notes)))
  expect_error(caller_normal(sf, regions = rg, n_draws = 500, pool_classes = character(), params = modifyList(p, list(calib_rho = "x"))), "calib_rho")
  car = lapply(sim_norm_f$truth$carriers, function(n) unique(sim_norm_f$scn$cells$cell[sim_norm_f$truth$nuclei$cell[n]]))
  rep_cells = function(r, reg, st) unique(r$calls[level == "cell" & called & reportable & region == reg & state == st, cell])
  expect_gte(mean(car[[2]] %in% rep_cells(r1, "chr8", "gain")), 0.6)
})

test_that("male chrX is depth-only in cell and pooled scores: no CN-LOH state, pooled null at the null level", {
  # male sample: chrX regions allelic_ok = FALSE and no allele votes; a depth bias shared by all cells on chrX (as a
  # neighbourhood-specific X accessibility), which only the chrX null pools can absorb
  rg = fake_regions(tiles = c("10Mb" = 1e7), sex = "male")
  expect_true(all(!rg[chrom == "chrX", allelic_ok]))
  ff = fake_features(400, rg, seed = 41)
  ff[chrom == "chrX", `:=`(nA = 0, nB = 0, n_hap = 0, f = NA_real_, z = z - 0.6)]
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  set.seed(42)
  emb = matrix(rnorm(400 * 5), 400, 5, dimnames = list(unique(ff$cell), NULL))
  res = caller_normal(ff, regions = rg, embedding = emb, sex = "male", null_cells = "features", n_draws = 1500)
  pc = res$calls[level == "pool" & chrom == "chrX"]
  expect_gt(nrow(pc), 100)
  expect_false(any(pc$state == "cnloh"))
  expect_lte(mean(pc$called), 0.01)
  expect_true(all(pc[called == TRUE, reportable]))                  # male chrX pools are reportable
  expect_false(any(res$calls[level == "cell" & chrom == "chrX", state] == "cnloh"))
  expect_true(all(is.infinite(res$scores[chrom == "chrX", cl_cnloh])))
  expect_true(all(is.finite(res$scores[chrom != "chrX", cl_cnloh])))
  # the chrX regions have their own pooled scale (from the null cells' chrX pools), not the class fallback
  ps = res$null_fit$pool_scale
  expect_true(all(rg[chrom == "chrX" & class %in% c("chrom", "arm", "10Mb"), region] %in% ps[kind == "type", region]))
  # pool_scores on a depth-only pool: no purely allelic state even when every other state is negative
  st = caller_states("normal")
  sc = res$scores[chrom == "chrX" & region == "chrX"]
  pl = data.table(pool = "p", cell = sc$cell[1:20])
  pz = pool_scores(sc, pl, st)
  expect_true(is.infinite(pz$cl_cnloh))
  expect_true(pz$state_comb %in% c("loss", "gain"))
  # pool_sex_thr (default): chrX has its own pooled threshold group; off: one threshold per class
  expect_true(caller_defaults()$pool_sex_thr)
  pt = res$thresholds[level %in% c("pool_type", "pool_knn")]
  expect_true("sexgrp" %in% names(pt))
  expect_setequal(unique(pt$sexgrp), c("A", "X"))
  expect_true(all(pt[sexgrp == "X", n_null] > 0))
  expect_true(all(is.na(res$thresholds[level == "cell", sexgrp])))
  r_off = caller_normal(ff, regions = rg, sex = "male", null_cells = "features", n_draws = 1500, params = modifyList(caller_defaults(), list(pool_sex_thr = FALSE)))
  expect_false("sexgrp" %in% names(r_off$thresholds))
  # female: no sex chromosomes in the pooled null, so no sex groups and the old path
  rf = fake_regions(tiles = c("10Mb" = 1e7))
  f2 = fake_features(300, rf, seed = 43)
  f2[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  r_f = caller_normal(f2, regions = rf, sex = "female", null_cells = "features", n_draws = 1500)
  expect_false("sexgrp" %in% names(r_f$thresholds))
  r_f0 = caller_normal(f2, regions = rf, sex = "female", null_cells = "features", n_draws = 1500, params = modifyList(caller_defaults(), list(pool_sex_thr = FALSE)))
  expect_identical(r_f$calls, r_f0$calls)
  expect_identical(r_f$thresholds, r_f0$thresholds)
})

# ---- allelic tile mask (mask_allelic) -------------------------------------------------------------------------------

# null cells x 5 Mb tiles of chr1-chr6 (10 tiles each), ~20 votes per row, beta-binomial at rho 0.01; one tile
# (5Mb|chr3:25000000) at rho 0.3 with ~60 votes per row; chromosome rows = the sums of the cell's tiles
mask_fake = function(n_cells = 600, seed = 13, rho = 0.01, noisy = "5Mb|chr3:25000000", noisy_rho = 0.3, noisy_n = 60) {
  set.seed(seed)
  tl = CJ(chrom = paste0("chr", 1:6), start = seq(0, 45e6, by = 5e6))
  tl[, region := tile_name("5Mb", chrom, start)]
  d = CJ(cell = sprintf("c%04d", seq_len(n_cells)), region = tl$region)
  d = merge(d, tl, by = "region")
  d[, `:=`(class = "5Mb", stratum = "20-50k", is_null = TRUE, fold = (match(cell, unique(cell)) - 1L) %% 4L + 1L)]
  d[, n_hap := as.numeric(rpois(.N, fifelse(region == noisy, noisy_n, 20)))]
  r = fifelse(d$region == noisy, noisy_rho, rho)
  p = rbeta(nrow(d), 0.5 * (1 - r) / r, 0.5 * (1 - r) / r)
  d[, nA := as.numeric(rbinom(.N, n_hap, p))]
  d[]
}

test_that("allelic_noisy_tiles masks exactly an injected overdispersed tile; its chromosome's false allelic calls drop", {
  d = mask_fake()
  expect_false(caller_defaults()$mask_allelic)
  expect_equal(caller_defaults()$mask_rho_mad, 2)
  m = allelic_noisy_tiles(d, k_mad = 4, min_votes = 3)
  expect_true(all(c("tile", "chrom", "start", "end", "n_rows", "votes", "rho", "z_rho", "assessed", "masked") %in% names(m)))
  expect_equal(nrow(m), 60L)
  expect_equal(m[masked == TRUE, tile], "5Mb|chr3:25000000")
  expect_equal(m[masked == TRUE, c(start, end)], c(25e6, 30e6))
  expect_equal(m[tile == "5Mb|chr1:0", rho], 0.01, tolerance = 0.3)
  expect_gt(m[masked == TRUE, rho], 0.08)
  expect_false(any(allelic_noisy_tiles(d, k_mad = 4, min_votes = 1000)$masked))     # too few votes: not assessed
  # chromosome-level allele-only calls with and without the masked tile's votes
  chrom_rows = function(x) x[, .(nA = sum(nA), n_hap = sum(n_hap)), by = .(cell, chrom, fold)][
    , `:=`(region = chrom, class = "chrom", stratum = "20-50k", is_null = TRUE, z = 0, sd_null = 0.2)][]
  st = caller_states("normal")
  rate = function(x) {
    sc = score_cells(x, st, NULL, 0.01, depth_model = "gauss")
    sc = apply_thresholds(sc, budget_thresholds(sc, "s_allelic", min_null = 100), "s_allelic")
    sc[, .(per_1000 = 1000 * mean(call_s_allelic)), by = chrom]
  }
  r0 = rate(chrom_rows(d))
  r1 = rate(chrom_rows(d[!region %in% m[masked == TRUE, tile]]))
  expect_gt(r0[chrom == "chr3", per_1000], 3)
  expect_lt(r1[chrom == "chr3", per_1000], r0[chrom == "chr3", per_1000] / 2)
})

test_that("masked atoms are absent from the allele region sums (allelic_features exclude_atoms, update_allelic_features)", {
  rg = fake_regions(tiles = c("5Mb" = 5e6))
  at = fake_atoms()
  v = fake_votes(40, at, rep("a", 40), NULL, flip_prob = 0, seed = 3, mol_per_atom = 2)
  tiles = data.table(tile = "5Mb|chr2:10000000", chrom = "chr2", start = 10e6, end = 15e6, masked = TRUE)
  ex = atoms_in_tiles(at, tiles)
  expect_equal(sort(ex), at[chrom == "chr2" & start >= 10e6 & start < 15e6, atom])
  expect_length(atoms_in_tiles(at, tiles[0]), 0)
  r = rg[class %in% c("chrom", "5Mb")]
  a0 = allelic_features(v, at, r, orientation = "chromosome")
  a1 = allelic_features(v, at, r, orientation = "chromosome", exclude_atoms = ex)
  expect_identical(a1, allelic_features(v[!atom %in% ex], at, r, orientation = "chromosome"))
  expect_false("5Mb|chr2:10000000" %in% a1$region)
  k = function(a, reg) a[region == reg][order(cell), nA + nB]
  vx = v[atom %in% ex, .(n = sum(nA + nB)), by = cell][order(cell)]
  expect_equal(k(a0, "chr2")[vx$cell] - k(a1, "chr2")[vx$cell], vx$n)
  expect_identical(a1[region == "chr1"], a0[region == "chr1"])
  # update_allelic_features on a features object: only the allele columns change, as allelic_features with the mask
  cells = sprintf("cell_%04d", 1:40)
  long = CJ(cell = cells, region = r$region)
  long = merge(long, r[, .(region, class, chrom, null_region, eligible, allelic_ok)], by = "region")
  long[, `:=`(z = 0.1, sd_null = 0.2)]
  long[a0[, .(cell = cells[cell], region, nA, nB, n_blocks, f)], `:=`(nA = i.nA, nB = i.nB, n_blocks = i.n_blocks, f = i.f), on = .(cell, region)]
  long[is.na(nA), `:=`(nA = 0, nB = 0, n_blocks = 0)]
  long[, n_hap := nA + nB]
  fo = list(long = long, regions = r, cells = data.table(cell = cells), votes = v, atoms = at, phasing = c(detected = "chromosome", orientation = "chromosome"))
  fu = update_allelic_features(fo, ex)
  u = fu$long
  chk = merge(u[, .(cell, region, nA, n_hap, z)], a1[, .(cell = cells[cell], region, nA_ref = nA, n_ref = n_hap)], by = c("cell", "region"), all.x = TRUE)
  chk[is.na(n_ref), `:=`(nA_ref = 0, n_ref = 0)]
  expect_equal(chk$nA, chk$nA_ref)
  expect_equal(chk$n_hap, chk$n_ref)
  expect_true(all(u$z == 0.1))
  expect_true(any(grepl("masked atoms", fu$notes)))
  expect_error(update_allelic_features(fo[setdiff(names(fo), "votes")], ex), "votes")
})

# ---- named events (class `event`: no null rows) ------------------------------------------------------------------

test_that("named event regions borrow the thresholds of their sd_class (per cell; pooled with pool_events); other classes unchanged", {
  rg0 = fake_regions(tiles = c("20Mb" = 2e7))
  ev_row = copy(rg0[region == "chr8"])[, `:=`(region = "chr8_event", class = "event", sd_class = "chrom", null_region = FALSE)]
  rg = rbind(rg0, ev_row)
  ff = fake_features(300, rg, events = data.table(cell = 1:12, region = "chr8_event", state = "loss"), seed = 21)
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  set.seed(5)
  emb = matrix(rnorm(300 * 5), 300, 5, dimnames = list(unique(ff$cell), NULL)); emb[1:12, 1] = emb[1:12, 1] + 30
  run = function(ff, rg, pe = NULL) caller_normal(ff, regions = rg, embedding = emb, null_cells = "features", n_draws = 1500,
                                                  params = if (is.null(pe)) caller_defaults() else modifyList(caller_defaults(), list(pool_events = pe)))
  res = run(ff, rg, pe = TRUE)
  th = res$thresholds
  expect_true("borrowed_from" %in% names(th))
  te = th[level == "cell" & class == "event"]
  expect_gt(nrow(te), 0)
  expect_false(anyNA(te$thr))
  expect_true(all(te$borrowed_from == "chrom"))
  tc = th[level == "cell" & class == "chrom"]
  k = match(paste(te$score, te$stratum, te$fold), paste(tc$score, tc$stratum, tc$fold))
  expect_false(anyNA(k))
  expect_identical(te$thr, tc$thr[k])
  expect_true(all(is.na(th[class != "event", borrowed_from])))
  # per-cell rows carry the borrowed threshold; the spiked carriers are called (reportable events)
  ev_sc = res$scores[region == "chr8_event"]
  expect_false(anyNA(ev_sc$thr_s_comb))
  ce = res$calls[level == "cell" & region == "chr8_event" & called & state == "loss"]
  expect_gte(sum(sprintf("cell_%04d", 1:12) %in% ce$cell), 10)
  expect_true(all(ce$reportable))
  expect_lt(sum(!ce$cell %in% sprintf("cell_%04d", 1:12)), 3)
  # pool_events: a non-NA borrowed pooled threshold, the event's pools are thresholded, and the run notes it
  for (lv in c("pool_type", "pool_knn")) {
    tp = th[level == lv & class == "event"]
    expect_gt(nrow(tp), 0)
    expect_false(anyNA(tp$thr))
    expect_true(all(tp$borrowed_from == "chrom"))
  }
  pe = res$calls[level == "pool" & region == "chr8_event"]
  expect_gt(nrow(pe), 0)
  expect_false(anyNA(pe$thr))
  expect_gt(sum(pe$called), 0)
  expect_length(grep("pooled calls of named events (pool_events, experimental)", res$notes, fixed = TRUE), 1)
  # default (pool_events FALSE): the event's pools are scored but get no borrowed threshold and are never called
  rd = run(ff, rg)
  expect_false(any(rd$thresholds[level != "cell", class] == "event"))
  pd = rd$calls[level == "pool" & region == "chr8_event"]
  expect_identical(nrow(pd), nrow(pe))
  expect_true(all(is.na(pd$thr)))
  expect_false(any(pd$called))
  expect_false(any(grepl("pool_events", rd$notes)))
  expect_false(any(grepl("no threshold", rd$notes)))
  # per-cell event calls, scores and thresholds do not depend on pool_events
  expect_identical(rd$scores, res$scores)
  expect_identical(rd$thresholds[level == "cell"], res$thresholds[level == "cell"])
  expect_identical(rd$calls[level == "cell"], res$calls[level == "cell"])
  # the other classes: identical pooled thresholds and calls with and without pool_events, and without the event row
  r0 = run(ff[class != "event"], rg0)
  expect_false("borrowed_from" %in% names(r0$thresholds))
  th1 = th[class != "event"][, borrowed_from := NULL]
  expect_identical(th1, r0$thresholds)
  expect_identical(rd$thresholds[class != "event"][, borrowed_from := NULL], r0$thresholds)
  expect_identical(res$calls[region != "chr8_event"], r0$calls)
  expect_identical(rd$calls[region != "chr8_event"], r0$calls)
  expect_identical(rd$null_fit$pool_scale, res$null_fit$pool_scale)
  # vote calibration: the event rows take the vote strata of their sd_class
  p = modifyList(caller_defaults(), list(calib = "votes", min_null_rows = 200))
  rv = caller_normal(ff, regions = rg, null_cells = "features", n_draws = 500, pool_classes = character(), params = p)
  sv = rv$scores[region == "chr8_event"]
  expect_true(all(sv$hstratum %in% c(rv$null_fit$vote_bins[class == "chrom", hstratum], "none")))
  expect_false(anyNA(sv$thr_s_comb))
  tv = rv$thresholds[level == "cell" & class == "chrom" & score == "s_comb"]
  expect_identical(sv$thr_s_comb, tv$thr[match(paste(sv$stratum, sv$hstratum, sv$fold), paste(tv$stratum, tv$hstratum, tv$fold))])
})

test_that("male sample: an autosomal and a chrX named event are both pooled-thresholded with pool_events; other classes unchanged", {
  # male: chrX is depth-only and its rows of the null cells join the pooled null in their own sex group. A named
  # event never has null rows of its own, so with pool_events an autosomal event borrows its sd_class threshold of
  # the autosomal group and a chrX event the one of the X group, whether or not the other event is present.
  rg0 = fake_regions(tiles = c("20Mb" = 2e7), sex = "male")
  mk = function(src, name) copy(rg0[region == src])[, `:=`(region = name, class = "event", sd_class = "chrom", null_region = FALSE)]
  rgA = rbind(rg0, mk("chr8", "ev8"))
  rgB = rbind(rgA, mk("chrX", "evX"))
  ffB = fake_features(300, rgB, events = data.table(cell = 1:12, region = "ev8", state = "loss"), seed = 21)
  ffB[chrom == "chrX", `:=`(nA = 0, nB = 0, n_hap = 0, f = NA_real_)]
  ffB[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  run = function(ff, rg, pe = TRUE) caller_normal(ff, regions = rg, null_cells = "features", n_draws = 1500, sex = "male",
                                                 params = modifyList(caller_defaults(), list(pool_events = pe)))
  rB = run(ffB, rgB)
  rA = run(ffB[region != "evX"], rgA)
  r0 = run(ffB[class != "event"], rg0)
  pt = rB$thresholds[level == "pool_type"]
  pe = rB$calls[level == "pool" & class == "event"]
  expect_setequal(unique(pe$region), c("ev8", "evX"))
  expect_false(anyNA(pe$thr))
  # each event is thresholded in the sex group of its chromosome, with its sd_class threshold
  te = pt[class == "event" & fold == 0L & score == "s_rob"]
  tc = pt[class == "chrom" & fold == 0L & score == "s_rob"]
  expect_true(all(te$borrowed_from == "chrom"))
  expect_true(all(pe[region == "ev8", thr] == tc[sexgrp == "A", thr]))
  expect_true(all(pe[region == "evX", thr] == tc[sexgrp == "X", thr]))
  # the autosomal event's pools are called as without the chrX event
  expect_identical(rB$calls[level == "pool" & region == "ev8"], rA$calls[level == "pool" & region == "ev8"])
  expect_gt(sum(rB$calls[level == "pool" & region == "ev8", called]), 0)
  # other classes: identical thresholds and calls without any event
  th1 = rB$thresholds[class != "event"][, borrowed_from := NULL]
  expect_identical(th1, r0$thresholds)
  expect_identical(rB$calls[class != "event"], r0$calls)
  # default (pool_events FALSE): no pooled event threshold in any sex group, no pooled event call; the rest unchanged
  rD = run(ffB, rgB, pe = FALSE)
  expect_false(any(rD$thresholds[level != "cell", class] == "event"))
  pd = rD$calls[level == "pool" & class == "event"]
  expect_setequal(unique(pd$region), c("ev8", "evX"))
  expect_true(all(is.na(pd$thr)))
  expect_false(any(pd$called))
  expect_identical(rD$thresholds[class != "event"][, borrowed_from := NULL], r0$thresholds)
  expect_identical(rD$calls[class != "event"], r0$calls)
  expect_identical(rD$calls[level == "cell"], rB$calls[level == "cell"])
  expect_identical(rD$scores, rB$scores)
})

test_that("a region left without thresholds after borrowing warns once and is noted; results unchanged", {
  # an event whose sd_class (10 Mb tiles) is not among the region classes: nothing to borrow, never called
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ev = copy(rg[region == "chr8"])[, `:=`(region = "ev8_10Mb", class = "event", sd_class = "10Mb", null_region = FALSE)]
  rg2 = rbind(rg, ev)
  ff = fake_features(200, rg2, events = data.table(cell = 1:10, region = "ev8_10Mb", state = "loss"), seed = 3)
  ws = character()
  r = withCallingHandlers(caller_normal(ff, regions = rg2, null_cells = "features", n_draws = 500),
                          warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_true(all(is.na(r$scores[region == "ev8_10Mb", thr_s_comb])))
  expect_length(grep("ev8_10Mb", ws), 1)
  expect_match(ws[grep("ev8_10Mb", ws)], "ev8_10Mb (per cell", fixed = TRUE)
  expect_true(any(grepl("no threshold for 1 region(s), never called: ev8_10Mb (per cell)", r$notes, fixed = TRUE)))
  # pool_events borrows nothing either (no pooled threshold of the sd_class): same calls and notes
  rp = suppressWarnings(caller_normal(ff, regions = rg2, null_cells = "features", n_draws = 500,
                                      params = modifyList(caller_defaults(), list(pool_events = TRUE))))
  expect_identical(rp$notes, r$notes)
  expect_identical(rp$calls, r$calls)
  r0 = caller_normal(ff[region != "ev8_10Mb"], regions = rg, null_cells = "features", n_draws = 500)
  expect_identical(r$calls[region != "ev8_10Mb"], r0$calls)
  expect_false(any(grepl("no threshold", r0$notes)))
})

test_that("without any allele votes caller_normal runs depth-only with one warning, in every null-cell path", {
  # an unphased extraction (or a VCF without phase sets): nA = nB = n_hap = 0 everywhere, rho cannot be fitted
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ff = fake_features(300, rg, seed = 3)
  ff[, `:=`(nA = 0, nB = 0, n_hap = 0, f = NA_real_)]
  run = function(...) {
    ws = character()
    r = withCallingHandlers(caller_normal(ff, regions = rg, n_draws = 500, ...),
                            warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(r = r, ws = ws)
  }
  for (a in list(list(null_cells = "features"), list(null_cells = "select"),
                 list(null_cells = "select", select_args = list(threshold = "empirical")),
                 list(null_cells = "features", params = modifyList(caller_defaults(), list(calib_rho = "chrom"))))) {
    x = do.call(run, a)
    expect_length(grep("no allele votes", x$ws), 1)
    expect_length(x$ws, 1)
    expect_true(any(grepl("no allele votes", x$r$notes)))
    al = grep("^(al|ar)_", names(x$r$scores), value = TRUE)
    expect_gt(length(al), 0)
    expect_true(all(as.matrix(x$r$scores[, al, with = FALSE]) == 0))     # every allelic term is 0: depth-only
    expect_gt(nrow(x$r$calls), 0)
  }
  # ll_allelic without votes does not need rho
  expect_identical(ll_allelic(c(0, 0), c(0, NA), 2 / 3, NA_real_), c(0, 0))
  # votes, but none on the null rows: a clear error
  f2 = copy(ff)[, null_cell := cell != "cell_0001"]
  f2[cell == "cell_0001", `:=`(nA = 3, nB = 1, n_hap = 4)]
  expect_error(caller_normal(f2, regions = rg, null_cells = "features", n_draws = 500), "null rows have no allele votes")
})

# ---- attenuation-aware state set of normal mode (params$normal_states, f_att_normal) ---------------------------------

# fake features with weak per-cell depth (sd_null x depth_k: a one-copy loss shifts z by about -0.7, as in sparse
# single-cell data, so the allele term decides the state) and a switch error: a share `w` of the votes of the event
# regions is flipped in every cell (A <-> B), so the carriers' oriented fraction falls from 0.98 to 0.98 - w x 0.96
att_features = function(n_cells, rg, ev, w, depth_k = 8, seed = 1) {
  ff = fake_features(n_cells, rg, events = ev, seed = seed)
  sh = caller_defaults()$shifts
  set.seed(seed + 100)
  ff[, sd_null := sd_null * depth_k]
  delta = fifelse(ff$truth == "loss", sh[["loss"]], fifelse(ff$truth == "gain", sh[["gain"]], 0))
  ff[, z := rt(.N, 7) + delta / pmax(sd_null, 0.05)]
  ff[, log2ratio := z * pmax(sd_null, 0.05)]
  ff[, x := (e + 0.5) * 2^log2ratio - 0.5]
  i = which(ff$region %in% unique(ev$region))
  xa = rbinom(length(i), ff$nA[i], w); xb = rbinom(length(i), ff$nB[i], w)
  ff[i, nA := nA - xa + xb]
  ff[i, `:=`(nB = n_hap - nA, f = fifelse(n_hap > 0, nA / n_hap, NA_real_))]
  ff[]
}

test_that("normal_att takes its grid from f_att_normal (tumour keeps f_att); normal_state_set resolves 'auto' from the regime", {
  p = modifyList(caller_defaults(), list(f_att_normal = c(0.98, 0.6, 0.75)))
  att = caller_states("normal_att", p)
  expect_equal(att[state == "loss", f], c(0.6, 0.75, 0.98))
  expect_equal(att[state == "cnloh", f], c(0.6, 0.75, 0.98))
  expect_equal(att[state == "loss", row_id], paste0("loss_", 1:3))
  expect_equal(att[state == "gain", f], c(2 / 3, 3 / 4))
  expect_identical(caller_states("tumour", p), caller_states("tumour"))
  expect_identical(caller_states("normal", p), caller_states("normal"))
  expect_equal(caller_states("normal_att")[state == "loss", f], caller_defaults()$f_att_normal)
  expect_lt(min(caller_defaults()$f_att_normal), 0.85)
  p0 = caller_defaults(); p0$f_att_normal = NULL                  # parameter lists without the key: f_att, as before
  expect_equal(caller_states("normal_att", p0)[state == "loss", f], caller_defaults()$f_att)
  # loss and CN-LOH share one allelic row set: one allelic class
  expect_equal(unname(scphasedcn:::allelic_classes(att[state != "neutral"])$label), c("loss_or_cnloh", "gain"))
  for (bad in list(0.4, c(0.7, 1), "0.7", NA_real_, numeric()))
    expect_error(caller_states("normal_att", modifyList(caller_defaults(), list(f_att_normal = bad))), "f_att_normal")
  # normal_state_set
  chr = c(detected = "chromosome", orientation = "chromosome")
  blk = c(detected = "block", orientation = "loo")
  st = function(v, ph = NULL, votes = TRUE) normal_state_set(modifyList(caller_defaults(), list(normal_states = v)), ph, votes)$set
  expect_equal(normal_state_set()$set, "normal")
  expect_equal(st("normal", chr), "normal")
  expect_equal(st("normal_att", blk), "normal_att")
  expect_equal(st("auto", chr), "normal_att")
  expect_equal(st("auto", c(detected = "chromosome", orientation = "loo")), "normal_att")   # the merged phase sets decide
  expect_equal(st("auto", blk), "normal")
  expect_equal(st("auto", NULL), "normal")
  expect_equal(st("auto", chr, votes = FALSE), "normal")
  expect_match(normal_state_set(modifyList(caller_defaults(), list(normal_states = "auto")), chr)$reason, "phasing regime chromosome")
  for (bad in list("x", NA_character_, c("normal", "auto"), 1, character()))
    expect_error(normal_state_set(modifyList(caller_defaults(), list(normal_states = bad))), "normal_states")
  p0 = caller_defaults(); p0$normal_states = NULL
  expect_equal(normal_state_set(p0)$set, "normal")
})

test_that("caller_normal: normal_states unset and 'normal' give identical results; invalid values stop", {
  rg = fake_regions(tiles = c("5Mb" = 5e6, "20Mb" = 2e7))
  ev = data.table(cell = c(1:10, 11:20), region = rep(c("chr1p", "chr8"), each = 10), state = rep(c("loss", "gain"), each = 10))
  ff = fake_features(240, rg, events = ev, seed = 3)
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  ff[, null_cell := NA]                                                  # null cells selected
  set.seed(2)
  emb = matrix(rnorm(240 * 5), 240, 5, dimnames = list(unique(ff$cell), NULL))
  p0 = caller_defaults(); p0$normal_states = NULL; p0$f_att_normal = NULL
  run = function(p, ...) caller_normal(ff, regions = rg, embedding = emb, n_draws = 1000, params = p, ...)
  r0 = run(p0)
  r1 = run(modifyList(caller_defaults(), list(normal_states = "normal")))
  for (k in setdiff(names(r0), "params")) expect_identical(r1[[k]], r0[[k]], label = k)
  expect_equal(r1$null_fit$state_set, "normal")
  expect_false("f_state" %in% names(r1$scores))
  expect_false("f_state" %in% names(r1$calls))
  expect_length(grep("^state set normal \\(normal_states = \"normal\"\\)$", r1$notes), 1)
  # "auto" on block phasing is the normal set: same results
  r2 = run(modifyList(caller_defaults(), list(normal_states = "auto")), phasing = c(detected = "block", orientation = "loo"))
  for (k in setdiff(names(r0), c("params", "notes"))) expect_identical(r2[[k]], r0[[k]], label = k)
  expect_true(any(grepl("normal_states = \"auto\": phasing regime block", r2$notes, fixed = TRUE)))
  for (bad in list(list(normal_states = "att"), list(normal_states = NA_character_), list(f_att_normal = c(0.3, 0.9))))
    expect_error(caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character(), params = modifyList(caller_defaults(), bad)),
                 names(bad))
})

test_that("normal_att calls attenuated losses as losses where the normal set labels them gains; 'auto' follows the regime", {
  rg = fake_regions(tiles = c("20Mb" = 2e7))
  ev = data.table(cell = c(1:20, 21:40), region = rep(c("chr1p", "chr2q"), each = 20), state = rep(c("loss", "cnloh"), each = 20))
  ff = att_features(300, rg, ev, w = 0.3)                                # carriers' fraction about 0.69
  car = ev[, .(cell = sprintf("cell_%04d", cell), region, state)]
  run = function(set, ...) caller_normal(ff, regions = rg, null_cells = "features", n_draws = 500, pool_classes = character(),
                                         params = modifyList(caller_defaults(), list(normal_states = set)), ...)
  ok = function(r, reg) {
    x = r$calls[level == "cell" & region == reg & cell %in% car[region == reg, cell]]
    c(state = sum(x$called & x$reportable & x$state == car[region == reg, state][1]), gain = sum(x$called & x$state == "gain"))
  }
  rn = run("normal")
  ra = run("normal_att")
  # the normal set: the attenuated 1p loss is a (non-reportable arm) gain, never a loss
  expect_lte(ok(rn, "chr1p")[["state"]], 1)
  expect_gte(ok(rn, "chr1p")[["gain"]], 8)
  # normal_att: most carriers are reportable losses, few labelled gains
  expect_gte(ok(ra, "chr1p")[["state"]], 7)
  expect_lt(ok(ra, "chr1p")[["gain"]], ok(rn, "chr1p")[["gain"]] - 5)
  expect_gte(ok(ra, "chr2q")[["state"]], ok(rn, "chr2q")[["state"]])
  # the scores carry one allelic row per grid point and the best-fitting fraction of the called state
  expect_true(all(paste0("ar_loss_", 1:4) %in% names(ra$scores)))
  fl = ra$calls[level == "cell" & called & state == "loss" & region == "chr1p" & cell %in% car$cell]
  expect_true(all(fl$f_state %in% c(0.70, 0.80)))
  sc = ra$scores
  i = which(sc$state_comb == "loss" & sc$n_hap > 0)
  best = caller_defaults()$f_att_normal[max.col(as.matrix(sc[i, paste0("ar_loss_", 1:4), with = FALSE]), ties.method = "first")]
  expect_equal(sc$f_state[i], best)
  expect_true(all(is.na(sc$f_state[sc$n_hap %in% 0])))
  expect_equal(ra$null_fit$state_set, "normal_att")
  expect_true(any(grepl("^state set normal_att \\(normal_states = \"normal_att\"\\); loss and CN-LOH rows over f 0.7, 0.8, 0.9, 0.98$", ra$notes)))
  # auto: chromosome-scale phase sets -> normal_att (identical results), block phasing -> normal
  rc = run("auto", phasing = c(detected = "chromosome", orientation = "chromosome"))
  expect_equal(rc$null_fit$state_set, "normal_att")
  expect_identical(rc$calls, ra$calls)
  expect_identical(rc$thresholds, ra$thresholds)
  expect_true(any(grepl("state set normal_att (normal_states = \"auto\": phasing regime chromosome)", rc$notes, fixed = TRUE)))
  rb = run("auto", phasing = c(detected = "block", orientation = "loo"))
  expect_equal(rb$null_fit$state_set, "normal")
  expect_identical(rb$calls, rn$calls)
})

test_that("normal_att is used by the null-cell selection, the pools and the named events", {
  rg0 = fake_regions(tiles = c("20Mb" = 2e7))
  ev_row = copy(rg0[region == "chr8"])[, `:=`(region = "chr8_event", class = "event", sd_class = "chrom", null_region = FALSE)]
  rg = rbind(rg0, ev_row)
  ev = data.table(cell = c(1:15, 16:30), region = rep(c("chr1p", "chr8_event"), each = 15), state = "loss")
  ff = att_features(300, rg, ev, w = 0.3, seed = 4)
  ff[, cell_type := rep(c("epi", "fib"), length.out = .N)[match(cell, unique(cell))]]
  ff[, null_cell := NA]
  set.seed(5)
  emb = matrix(rnorm(300 * 5), 300, 5, dimnames = list(unique(ff$cell), NULL))
  pa = modifyList(caller_defaults(), list(normal_states = "normal_att"))
  r = caller_normal(ff, regions = rg, embedding = emb, n_draws = 1000, params = pa)
  # null-cell selection under the normal_att set
  sa = select_null_cells(ff, caller_states("normal_att", pa), pa)
  sn = select_null_cells(ff, caller_states("normal", pa), pa)
  expect_identical(r$null_fit$selection$thresholds, sa$thresholds)
  expect_identical(r$null_fit$selection$null_cells, sa$null_cells)
  expect_false(identical(sa$thresholds$thr, sn$thresholds$thr))
  # select_args can still pass their own states
  r_s = caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character(), params = pa,
                      select_args = list(states = caller_states("normal", pa)))
  expect_identical(r_s$null_fit$selection$thresholds, sn$thresholds)
  # pooled scores of the cell-type pools: the normal_att states with the run's region scale
  sca = r$null_fit$pool_scale[kind == "type"][, kind := NULL]
  sp = r$scores[class %in% c("chrom", "arm", "50Mb", "20Mb", "10Mb", "5Mb", "event") & eligible %in% TRUE]
  ps = pool_scores(sp, r$pools[kind == "type", .(pool, cell)], caller_states("normal_att", pa), by = c("region", "class"),
                   loo = TRUE, scale = sca)
  pc = r$calls[level == "pool" & grepl("^(epi|fib)\\|pool", cell)]
  k = match(paste(pc$cell, pc$region), paste(ps$pool, ps$region))
  expect_false(anyNA(k))
  expect_equal(pc$s_comb, ps$s_comb[k])
  expect_equal(pc$state, ps$state_rob[k])
  expect_true(all(c("pool_type", "pool_knn") %in% r$thresholds$level))
  # every scoring step of the run receives the normal_att state table (the pooled scores of null pools hardly differ
  # between the sets, so the states passed are recorded instead)
  seen = list()
  rec = function(nm, f) function(...) {
    a = as.list(match.call(f, sys.call()))[-1]
    seen[[nm]] <<- c(seen[[nm]], list(eval(a$states, parent.frame())))
    f(...)
  }
  ns = asNamespace("scphasedcn")
  fns = c("select_null_cells", "score_cells", "pooled_thresholds", "call_pools")
  orig = mget(fns, envir = ns)
  do.call(local_mocked_bindings, c(lapply(setNames(fns, fns), function(nm) rec(nm, orig[[nm]])), list(.package = "scphasedcn", .env = environment())))
  invisible(caller_normal(ff, regions = rg, embedding = emb, n_draws = 500, params = pa))
  want = caller_states("normal_att", pa)
  expect_setequal(names(seen), fns)
  for (nm in fns) for (x in seen[[nm]]) expect_identical(x, want, label = nm)
  expect_gte(length(seen$pooled_thresholds), 2)                         # cell-type and kNN pools
  expect_gte(length(seen$call_pools), 2)
  # the named event borrows the per-cell thresholds of its sd_class (computed under normal_att); its attenuated
  # carriers are reportable losses, more of them than under the normal set
  th = r$thresholds[level == "cell"]
  te = th[class == "event"]; tc = th[class == "chrom"]
  expect_gt(nrow(te), 0)
  expect_identical(te$thr, tc$thr[match(paste(te$score, te$stratum, te$fold), paste(tc$score, tc$stratum, tc$fold))])
  ce = r$calls[level == "cell" & region == "chr8_event" & called & reportable & state == "loss"]
  expect_gte(sum(sprintf("cell_%04d", 16:30) %in% ce$cell), 5)
  rn = caller_normal(ff, regions = rg, n_draws = 500, pool_classes = character())
  cn = rn$calls[level == "cell" & region == "chr8_event" & called & reportable & state == "loss"]
  expect_gt(sum(sprintf("cell_%04d", 16:30) %in% ce$cell), sum(sprintf("cell_%04d", 16:30) %in% cn$cell))
})
