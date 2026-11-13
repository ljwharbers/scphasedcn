# test-clones.R: tumour mode (R/clones_*.R) on hand-made tables, fake_votes and the cached
# simulations sim_tum2 (two clones + normal, unoriented short phase blocks) and sim_norm_f (normal tissue).
# The caller input is built here from a simulation (model depth, molecule votes, all-cell leave-one-out
# orientation) until the feature step's cell_features() is available; results are cached in caller_test_cache.

clone_cached = function(name, build) {
  if (is.null(caller_test_cache[[name]])) caller_test_cache[[name]] = build()
  caller_test_cache[[name]]
}

# scn_features-like list from a simulation: features (arms, chromosomes, 5 / 10 / 20 Mb tiles), regions, votes,
# atoms (one per bin), cells, 1 Mb window depth (win) and null fits. Depth: model expected counts, each cell centred
# on the median log2 ratio of its autosomal 5 Mb windows. Null cells = non-multiplet cells of `null_types`.
clone_sim_input = function(sim, null_types = "normal", tiles = c("5Mb" = 5e6, "10Mb" = 1e7, "20Mb" = 2e7), mol_gap = 2000) {
  scn = sim$scn
  nm = suppressWarnings(normalise_depth(scn, "model"))
  bins = scn$bins
  sex = scn$cells$sex[1]
  regions = fake_regions(bins, tiles = tiles, sex = sex)
  if (sex == "male") regions[chrom %in% c("chrX", "chrY"), null_region := FALSE]
  cells = scn$cells$cell
  use = nm$use
  X = as.matrix(scn$counts); X[!use, ] = 0
  E = nm$expected; E[!use, ] = 0
  w5 = window_ids(bins, 5e6); w5[!use | bins$chrom %in% c("chrX", "chrY")] = NA
  i5 = group_indicator(w5)
  med = apply(log2(as.matrix(i5 %*% X) / as.matrix(i5 %*% E)), 2, median)
  E = sweep(E, 2, 2^med, "*")
  wid = window_ids(bins, 1e6)
  ind = group_indicator(wid)
  wt = window_table(bins, wid)[match(rownames(ind), window)]
  win = list(windows = wt[, .(window, chrom, start, end)], x = as.matrix(ind %*% X), e = as.matrix(ind %*% E))
  colnames(win$x) = colnames(win$e) = cells
  rb = lapply(seq_len(nrow(regions)), function(r) which(bins$chrom == regions$chrom[r] & bins$start >= regions$start[r] & bins$end <= regions$end[r]))
  RI = Matrix::sparseMatrix(i = rep(seq_along(rb), lengths(rb)), j = unlist(rb), x = 1, dims = c(nrow(regions), nrow(bins)))
  xr = as.matrix(RI %*% X); er = as.matrix(RI %*% E)
  d = data.table(cell = rep(cells, each = nrow(regions)), ci = rep(seq_along(cells), each = nrow(regions)),
                 region = rep(regions$region, length(cells)), x = as.vector(xr), e = as.vector(er))
  d = merge(d, regions[, .(region, class, sd_class, chrom, null_region, eligible, allelic_ok)], by = "region", sort = FALSE)
  d[, log2ratio := log2((x + 0.5) / (e + 0.5))]
  d[, null_cell := (scn$cells$cell_type %in% null_types & !sim$truth$cells$multiplet)[ci]]
  d[, sd_null := mad(log2ratio[null_region]), by = .(cell, sd_class)]
  d[, z := log2ratio / pmax(sd_null, 0.05)]
  fr = sim$fragments[hap %in% c(1L, 2L) & ps > 0, .(cell, chrom, start, end, hap = ifelse(hap == 1L, "A", "B"), ps)]
  setorder(fr, cell, chrom, start)
  fr[, new := c(TRUE, start[-1] > cummax(end)[-.N] + mol_gap), by = .(cell, chrom)]
  fr[, mol := cumsum(new)]
  m = fr[, .(start = min(start), ps = ps[1], nA = sum(hap == "A"), nB = sum(hap == "B")), by = .(cell, chrom, mol)][nA != nB]
  if (sex == "male") m = m[!chrom %in% c("chrX", "chrY")]
  atoms = fake_atoms(bins)
  m[, atom := locate_bins(bins, chrom, start)]
  votes = m[, .(nA = sum(nA > nB), nB = sum(nB > nA)), by = .(cell, atom, ps)]
  votes[, `:=`(cell = as.integer(cell), atom = as.integer(atom), ps = as.numeric(ps), nA = as.integer(nA), nB = as.integer(nB))]
  setorder(votes, cell, atom, ps)
  ob = cl_orient(votes, regions, atoms, groups = NULL, loo = TRUE)
  ob[, cell := cells[cell]]
  d = merge(d, ob, by = c("cell", "region"), all.x = TRUE, sort = FALSE)
  d[is.na(nA) | allelic_ok == FALSE, `:=`(nA = 0, nB = 0, n_blocks = 0)]
  d[, n_hap := nA + nB][, f := fifelse(n_hap > 0, nA / n_hap, NA_real_)]
  nf = scn$cells$n_frag
  d[, `:=`(orient = "pseudobulk", n_frag = as.numeric(nf[ci]), stratum = as.character(depth_stratum(nf[ci])),
           cell_type = scn$cells$cell_type[ci], fold = as.integer((ci - 1L) %% 4L + 1L))]
  d[, c("ci", "allelic_ok") := NULL]
  nr = d[null_cell & null_region & eligible]
  list(features = d[], regions = regions, votes = votes, atoms = atoms,
       cells = data.table(cell = cells, n_frag = as.numeric(nf), stratum = as.character(depth_stratum(nf)), cell_type = scn$cells$cell_type),
       win = win, null_fit = list(tpar = fit_depth_null(nr, by = c("class", "stratum")), rho = fit_allelic_null(nr$nA, nr$n_hap)))
}

sim_truth = function(sim) data.table(cell = sim$scn$cells$cell, type = sim$scn$cells$cell_type,
                                     mult = sim$truth$cells$multiplet, types = sim$truth$cells$types)

# fake clone profile: windows of 1 Mb on chr1-chr3 (40 windows each) with given states; depth c * log2(cn / 2) +
# N(0, sd); votes Binomial(n_hap, 0.5 + a (f - 0.5)) on the major haplotype (clone_profiles' split-half orientation is
# one-sided: a block's major haplotype in one half is its major haplotype in the other; `flip = TRUE` flips windows at
# random). `mix`: a list of two state vectors averaged 50:50 (a doublet profile).
fake_clone_profile = function(states, clone = "K1", c_ = 0.8, a_ = 0.85, sd = 0.05, n_hap = 400, seed = 1, mix = NULL, flip = FALSE) {
  set.seed(seed)
  w = data.table(chrom = rep(c("chr1", "chr2", "chr3"), each = 40), start = rep((0:39) * 1e6, 3))
  w[, end := start + 1e6]
  w[, window := paste0(chrom, ":", format(start, scientific = FALSE, trim = TRUE))]
  if (is.null(mix)) {
    nA_ = as.integer(sub(":.*", "", states)); nB_ = as.integer(sub(".*:", "", states))
    cn = nA_ + nB_; fs = ifelse(cn > 0, nA_ / cn, 0.5)
  } else {
    a1 = as.integer(sub(":.*", "", mix[[1]])); b1 = as.integer(sub(".*:", "", mix[[1]]))
    a2 = as.integer(sub(":.*", "", mix[[2]])); b2 = as.integer(sub(".*:", "", mix[[2]]))
    cn = (a1 + b1 + a2 + b2) / 2; fs = (a1 + a2) / (a1 + b1 + a2 + b2)
  }
  w[, `:=`(clone = clone, ploidy = 2L, base = 2L, n_cells = 100L, log2 = c_ * log2(pmax(cn, 0.2) / 2) + rnorm(.N, 0, sd), log2_se = sd)]
  fe = 0.5 + a_ * (fs - 0.5)
  x = rbinom(nrow(w), n_hap, fe)
  fl = runif(nrow(w)) < 0.5 & flip
  w[, `:=`(nA = ifelse(fl, n_hap - x, x), n_hap = n_hap)]
  w[, `:=`(nB = n_hap - nA, f = nA / n_hap, f_se = sqrt(0.25 / n_hap), n_blocks = 20, ref = "pool", size = 1e6)]
  w[]
}
fake_states = function() {
  s = rep("1:1", 120)
  s[11:25] = "1:0"          # chr1 10-25 Mb loss
  s[40 + 21:40] = "2:1"     # chr2 20-40 Mb gain
  s[80 + 1:15] = "3:0"      # chr3 0-15 Mb allele-specific gain
  s[80 + 31:40] = "2:0"     # chr3 30-40 Mb CN-LOH
  s
}

# ---- unit tests on hand-made inputs ---------------------------------------------------------------

test_that("cn_grid has the states and labels of the tumour model", {
  g = cn_grid()
  expect_setequal(g$state, c("1:1", "1:0", "2:0", "2:1", "3:0", "3:1", "2:2", "4:1", "4:0", "0:0"))
  expect_equal(g[state == "3:0", label], "gain_as")
  expect_equal(g[state == "2:0", label], "cnloh")
  expect_equal(g[state == "2:1", f], 2 / 3)
  expect_true(all(g[state %in% c("3:1", "2:2", "4:1", "4:0"), label] == "amp"))
  expect_false("4:1" %in% cn_grid(4)$state)
})

test_that("clone_ari is 1 for identical partitions and ~0 for unrelated ones", {
  x = rep(1:3, each = 50)
  expect_equal(clone_ari(x, paste0("k", x)), 1)
  set.seed(1)
  expect_lt(abs(clone_ari(x, sample(x))), 0.05)
})

test_that("fit_clone_states recovers the depth contrast and attenuation of a fake profile", {
  p = fake_clone_profile(fake_states(), c_ = 0.8, a_ = 0.85)
  fit = fit_clone_states(p)
  expect_lte(abs(fit$par$c - 0.8), 0.1)
  expect_lte(abs(fit$par$a - 0.85), 0.1)
  seg = segment_clone(fit)
  check_segments(seg)
  get = function(ch, pos) seg[chrom == ch & start <= pos & end > pos, state]
  expect_equal(get("chr1", 17e6), "1:0")
  expect_equal(get("chr2", 30e6), "2:1")
  expect_equal(get("chr3", 7e6), "3:0")
  expect_equal(get("chr3", 35e6), "2:0")
  expect_equal(get("chr1", 35e6), "1:1")
})

test_that("clone_ami is 1 for identical partitions and ~0 for unrelated ones", {
  x = rep(1:4, each = 40)
  expect_equal(clone_ami(x, paste0("k", x)), 1)
  set.seed(2)
  expect_lt(abs(clone_ami(x, sample(x))), 0.05)
  expect_true(is.na(clone_ami(1, 2)))
  y = x; y[1:20] = 2                      # half of cluster 1 merged into cluster 2
  expect_gt(clone_ami(x, y), 0.6); expect_lt(clone_ami(x, y), 1)
})

test_that("split_win thins the window counts into two halves", {
  win = list(windows = data.table(window = c("chr1:0", "chr1:1000000"), chrom = "chr1", start = c(0, 1e6), end = c(1e6, 2e6)),
             x = matrix(c(10, 20, 30, 40), 2, dimnames = list(NULL, c("a", "b"))), e = matrix(c(12, 18, 33, 37), 2, dimnames = list(NULL, c("a", "b"))))
  s = split_win(win, 0.5, seed = 3)
  expect_equal(s$disc$x + s$test$x, win$x)
  expect_equal(s$disc$e + s$test$e, win$e)
  expect_true(all(s$disc$x >= 0 & s$test$x >= 0))
})

test_that("two_means separates a bimodal vector and finds little in a unimodal one", {
  set.seed(4)
  x = c(rnorm(100, -3), rnorm(60, 0))
  tm = two_means(x)
  expect_gt(tm$sep, 2.5)
  expect_equal(sum(tm$g == tm$g[1]), 100, tolerance = 5)
  expect_lt(two_means(rnorm(200))$sep, 3)
})

test_that("clone_state_distance and clone_trunks group clones that share their allelic events", {
  seg = function(k, ch, s, e, st) data.table(clone = k, chrom = ch, start = s, end = e, state = st)
  base = function(k) rbind(seg(k, "chr1", 0, 100e6, "1:1"), seg(k, "chr2", 0, 100e6, "1:1"), seg(k, "chr3", 0, 100e6, "1:1"))
  sg = rbind(base("N"),
             base("R1")[chrom == "chr1", state := "1:0"], base("R2")[chrom == "chr1", state := "1:0"][chrom == "chr3", state := "4:0"],
             base("R3")[chrom == "chr1", state := "1:0"][chrom == "chr3", state := "3:0"],
             base("L1")[chrom == "chr2", state := "2:1"])
  sg[, cn := as.integer(nA_of(state) + nB_of(state))]
  sg[, `:=`(ploidy = 2L, label = state_label(nA_of(state), nB_of(state), 2L))]
  d = clone_state_distance(sg)
  expect_equal(d["R1", "N"], 100e6)
  expect_equal(d["R2", "R3"], 0)                 # 4:0 vs 3:0: depth only, both f = 1
  expect_equal(d["R1", "R2"], 100e6)             # 1:1 vs 4:0
  expect_equal(d["R1", "L1"], 200e6)
  expect_equal(clone_state_distance(sg, min_df = 0.01, min_dcn = 1)["R2", "R3"], 100e6)
  tr = clone_trunks(sg, max_dist = 150e6)
  expect_equal(tr[["R1"]], tr[["R2"]]); expect_equal(tr[["R2"]], tr[["R3"]])
  expect_false(tr[["L1"]] == tr[["R1"]])
})

test_that("segment_clone reproduces a hand-made 3-segment chromosome", {
  g = cn_grid()
  truth = c(rep("1:1", 10), rep("1:0", 10), rep("2:1", 10))
  ll = CJ(w = 1:30, state = g$state)
  ll[, `:=`(clone = "K", chrom = "chr5", start = (w - 1) * 1e6, end = w * 1e6)]
  ll[, ll := ifelse(state == truth[w], 0, -5)]
  ll[w == 15 & state == "2:1", ll := 1]       # one noisy window does not open a segment (min_win 3)
  seg = segment_clone(ll)
  expect_equal(nrow(seg), 3)
  expect_equal(seg$start, c(0, 10e6, 20e6))
  expect_equal(seg$end, c(10e6, 20e6, 30e6))
  expect_equal(seg$state, c("1:1", "1:0", "2:1"))
  expect_equal(seg$label, c("neutral", "loss", "gain"))
  expect_equal(seg$cn, c(2L, 1L, 3L))
})

test_that("cl_orient (orient_blocks) orients blocks by the group pseudobulk with the cell left out", {
  atoms = fake_atoms()
  cl = rep(c("A", "B"), each = 30)
  ev = data.table(clone = "A", chrom = "chr22", start = 0, end = NA, cnA = 1L, cnB = 0L)
  v = fake_votes(60, atoms, cl, ev, mol_per_atom = 1, seed = 3)
  rg = fake_regions()[region %in% c("chr22", "chr21")]
  ob = cl_orient(v, rg, atoms, groups = cl)
  ob[, `:=`(clone = cl[cell], f = nA / (nA + nB))]
  expect_gt(ob[clone == "A" & region == "chr22", median(f)], 0.9)
  expect_lt(ob[clone == "B" & region == "chr22", abs(mean(f) - 0.5)], 0.05)
  expect_lt(ob[region == "chr21", abs(mean(f) - 0.5)], 0.05)       # leave-one-out: no self-orientation bias
  # without groups the all-cell pseudobulk still orients clone A's loss (it is the only imbalance)
  o2 = cl_orient(v, rg, atoms, groups = NULL)
  o2[, `:=`(clone = cl[cell], f = nA / (nA + nB))]
  expect_gt(o2[clone == "A" & region == "chr22", median(f)], 0.8)
  expect_named(ob[, 1:5], c("cell", "region", "nA", "nB", "n_blocks"))
})

test_that("flag_doublet_clones dissolves a clone explained by a 50:50 mixture of two others", {
  s = fake_states()
  pn = fake_clone_profile(rep("1:1", 120), clone = "N", seed = 2)
  pt = fake_clone_profile(s, clone = "T", seed = 3)
  pm = fake_clone_profile(NULL, clone = "M", seed = 4, mix = list(rep("1:1", 120), s))
  fit = fit_clone_states(rbind(pn, pt, pm))
  seg = segment_clone(fit)
  dc = flag_doublet_clones(data.table(clone = c("N", "T", "M"), n_cells = c(300L, 200L, 25L)), seg, fit$windows, fit)
  expect_true(dc[clone == "M", doublet])
  expect_false(any(dc[clone != "M", doublet]))
  expect_equal(dc[clone == "M", best_pair], "N+T")
})

# ---- sim_tum2: the full tumour mode ----------------------------------------------------------------

tum2_input = function() clone_cached("clones_tum2_input", function() clone_sim_input(sim_tum2))
tum2_result = function() clone_cached("clones_tum2_result", function() caller_tumour(tum2_input(), verbose = FALSE))

test_that("clone_feature_matrix orients unoriented input by the all-cell pseudobulk", {
  inp = tum2_input()
  X1 = clone_feature_matrix(inp$features, rho = inp$null_fit$rho)
  f0 = copy(inp$features)[, `:=`(orient = "none", nA = nB, nB = nA)]    # scrambled, unoriented
  X0 = clone_feature_matrix(f0, rho = inp$null_fit$rho, votes = inp$votes, atoms = inp$atoms, regions = inp$regions,
                            cells = inp$cells$cell)
  expect_equal(dim(X0), dim(X1))
  expect_equal(X0[, colnames(X1)], X1)
  expect_true(all(abs(X1) <= 6))
})

test_that("caller_tumour finds the clones of sim_tum2 and their states", {
  res = tum2_result()
  tr = sim_truth(sim_tum2)
  expect_s3_class(res, "scphasedcn")
  expect_equal(res$mode, "tumour")
  check_segments(res$segments)
  check_calls(res$calls)
  a = merge(res$assignments, tr, by = "cell")
  sing = a[mult == FALSE]
  expect_gte(clone_ari(sing$clone, sing$type), 0.9)
  map = sing[, .N, by = .(clone, type)][order(-N)][!duplicated(clone)]
  cl = function(type_) map[type == type_, clone]
  expect_length(unique(map$type), 3)
  sg = res$segments
  at = function(k, ch, pos) sg[clone == k & chrom == ch & start <= pos & end > pos]
  expect_equal(at(cl("cloneB"), "chr22", 10e6)$state, "1:0")
  expect_equal(at(cl("cloneA"), "chr21", 10e6)$state, "2:1")
  expect_equal(at(cl("cloneB"), "chr1", 30e6)$state, "2:0")
  for (k in c(cl("cloneA"), cl("cloneB"))) {
    s8 = at(k, "chr8", 22e6)
    expect_equal(s8$label, "amp")
    expect_lte(abs(s8$start - 15e6), 2e6)
  }
  s2 = at(cl("cloneB"), "chr2", 17e6)
  expect_equal(s2$label, "loss")
  expect_lte(abs(s2$start - 10e6), 2e6)
  expect_lte(abs(s2$end - 25e6), 2e6)
  expect_true(all(sg[clone == cl("normal"), label] == "neutral"))
  expect_true(res$clones[clone == cl("normal"), normal])
  # assignment, doublets, reference pool
  expect_gte(sing[, mean(map$type[match(clone, map$clone)] == type)], 0.95)
  ab = a[grepl("cloneA", types) & grepl("cloneB", types)]
  expect_gt(nrow(ab), 0)
  expect_gte(mean(ab$doublet_like), 0.5)
  pool = a[ref_pool == TRUE]
  expect_gt(nrow(pool), 10)
  expect_gte(mean(pool$type == "normal" & !pool$mult), 0.98)
  # per-cell calls: inherited from the clone segments, mostly confirmed on the cells' own evidence
  cc = res$calls[level == "cell"]
  expect_true(all(cc$source %in% c("inherited", "confirmed")))
  expect_equal(cc$called, cc$score > cc$thr)
  expect_gt(cc[grepl("chr22", region), mean(source == "confirmed")], 0.8)
  expect_true(all(res$calls[level == "clone", called]))
  expect_true(any(grepl("the scorer", res$notes)))
  check_thresholds(res$thresholds)
})

test_that("caller_tumour: the defaults of confirm_sided, perm_decision and profile_orient are today's code paths", {
  # all-default results equal the run with the three options spelled out at their pre-option behaviour (the run before
  # these options existed is bit-identical: assignments, segments, calls; checked against the previous release)
  d = clone_defaults()
  expect_identical(d[c("confirm_sided", "perm_decision", "profile_orient")], list(confirm_sided = "two", perm_decision = "quantile", profile_orient = "ps"))
  r0 = tum2_result()
  r1 = caller_tumour(tum2_input(), opts = list(confirm_sided = "two", perm_decision = "quantile", profile_orient = "ps"), verbose = FALSE)
  for (nm in c("assignments", "segments", "calls", "clones", "thresholds", "profiles", "events", "null_fit", "notes")) expect_identical(r1[[nm]], r0[[nm]])
  expect_false("p" %in% names(r0$null_fit$clones$tests))
})

test_that("caller_tumour stops on unknown or invalid tumour options before any work", {
  inp = tum2_input()
  expect_error(caller_tumour(inp, opts = list(split_pass = TRUE), verbose = FALSE), "split_pass")
  expect_error(caller_tumour(inp, opts = list(ref_null = "x"), verbose = FALSE), "ref_null")
  expect_error(caller_tumour(inp, opts = list(confirm_null = "x"), verbose = FALSE), "should be one of")
  expect_error(caller_tumour(inp, opts = list(confirm_sided = "x"), verbose = FALSE), "should be one of")
  expect_error(caller_tumour(inp, opts = list(perm_decision = "x"), verbose = FALSE), "should be one of")
  expect_error(caller_tumour(inp, opts = list(profile_orient = "x"), verbose = FALSE), "should be one of")
  expect_error(caller_tumour(inp, opts = list(bogus = 1), verbose = FALSE), "unknown option")
  expect_identical(scphasedcn:::check_clone_opts(list()), clone_defaults())
  expect_identical(scphasedcn:::check_clone_opts(list(split_pass = "disagree"))$split_pass, "disagree")
  # segment_clone() rules: clone_defaults() keys with segment_clone()'s defaults, checked as numbers >= 0
  rules = c("min_excess", "max_contra", "bal_frac", "min_bal_len")
  expect_equal(clone_defaults()[rules], lapply(formals(segment_clone)[rules], eval))
  expect_identical(scphasedcn:::check_clone_opts(list(max_contra = Inf, min_excess = 0))[rules],
                   list(min_excess = 0, max_contra = Inf, bal_frac = 0.75, min_bal_len = 2e7))
  expect_error(caller_tumour(inp, opts = list(bal_frac = -1), verbose = FALSE), "bal_frac")
  expect_error(caller_tumour(inp, opts = list(min_bal_len = "x"), verbose = FALSE), "min_bal_len")
})

test_that("control_regions and segment_pieces partition the clone breakpoints", {
  res = tum2_result()
  pcs = segment_pieces(res$segments, tum2_input()$regions)
  check_regions(pcs)
  expect_true(all(grepl("^seg\\|all\\|", pcs$region)))
  expect_equal(pcs[chrom == "chr2", sort(unique(c(start, end)))], c(0, 10e6, 25e6, 40e6))
  ctl = control_regions(res$segments, min_len = 5e6)
  expect_true(all(ctl$end - ctl$start >= 5e6))
  expect_true(all(c("seg|control|chr2:0-10000000", "seg|control|chr8:0-15000000") %in% ctl$region))
  expect_false(any(ctl$chrom %in% c("chr1", "chr21", "chr22")))
})

test_that("region_features sums the 1 Mb windows inside each region and scales by the cell's class sd_null", {
  inp = clone_input(tum2_input())
  rg = inp$regions[class %in% c("5Mb", "20Mb", "arm")]
  cc = inp$cells$cell[c(1, 50, 120)]
  d = region_features(inp, rg, cells = cc, sd_floor = 0.1)
  expect_equal(nrow(d), nrow(rg) * length(cc))
  expect_setequal(names(d), c("cell", "region", "x", "e", "class", "sd_class", "chrom", "log2ratio", "sd_null", "z", "n_frag", "stratum"))
  # windows fully inside the region, by hand
  wt = inp$win$windows
  r = rg[class == "20Mb"][2]
  w = wt[chrom == r$chrom & start >= r$start & end <= r$end, window]
  expect_equal(d[cell == cc[2] & region == r$region, x], sum(inp$win$x[w, cc[2]]))
  expect_equal(d[cell == cc[2] & region == r$region, e], sum(inp$win$e[w, cc[2]]))
  # Mb-aligned tiles: the same counts as the bin-level feature table, the cell's sd_null of the class
  f = inp$features[cell %in% cc & class == "5Mb", .(cell, region, x0 = x, e0 = e, sd0 = sd_null)]
  m = merge(d, f, by = c("cell", "region"))
  expect_equal(nrow(m), nrow(f))
  expect_equal(m$x, m$x0)
  expect_equal(m$e, m$e0)
  expect_equal(m$sd_null, m$sd0)
  expect_equal(d$z, d$log2ratio / pmax(d$sd_null, 0.1))
  expect_equal(d$log2ratio, log2((d$x + 0.5) / (d$e + 0.5)))
  # a class missing from the features takes a nearest available class's sd_null (50Mb: arm and 20Mb are one step away)
  r2 = copy(r)[, `:=`(region = "x|1", class = "50Mb", sd_class = "50Mb")]
  d2 = region_features(inp, r2, cells = cc[1])
  expect_true(d2$sd_null %in% inp$features[cell == cc[1] & sd_class %in% c("arm", "20Mb"), unique(sd_null)])
  expect_error(region_features(within(unclass(inp), win <- NULL), rg), "lacks the 1 Mb window depth")
})

test_that("region_block_votes sums the votes of overlapping regions as the feature step's per-layer layer_blocks", {
  inp = clone_input(tum2_input())
  rg = inp$regions[class %in% c("arm", "20Mb")]
  rb = scphasedcn:::region_block_votes(inp$votes, inp$atoms, rg)
  m = scphasedcn:::atom_region_map(inp$atoms, rg)
  expect_gt(max(m$layer), 1L)                     # arms and 20 Mb tiles overlap
  lb = rbindlist(lapply(sort(unique(m$layer)), function(l) scphasedcn:::layer_blocks(inp$votes, m, l)))
  lb[, region := rg$region[ri]]
  x = merge(rb, lb, by = c("cell", "region", "ps"))
  expect_equal(nrow(x), nrow(rb))
  expect_equal(nrow(x), nrow(lb))
  expect_equal(x$nA.x, x$nA.y)
  expect_equal(x$nB.x, x$nB.y)
  # one block id per region x phase set, numbered in order of first appearance
  expect_equal(uniqueN(rb[, .(region, ps)]), max(rb$blk))
  expect_equal(unique(rb$blk), seq_len(max(rb$blk)))
})

# ---- sim_norm_f: no over-splitting -------------------------------------------------------------------

norm_input = function() clone_cached("clones_norm_input", function() clone_sim_input(sim_norm_f, null_types = c("epi", "fib", "imm")))

test_that("normal tissue (sim_norm_f) is one clone", {
  # as caller_tumour: discovery on one half of the molecule votes, the merge test on the other half
  inp = norm_input()
  cells = inp$cells$cell
  vs = split_votes(inp$votes, 0.5, seed = 1)
  expect_equal(sum(vs$disc$nA) + sum(vs$test$nA), sum(inp$votes$nA))
  X = clone_feature_matrix(inp$features, rho = inp$null_fit$rho, votes = vs$disc, atoms = inp$atoms, regions = inp$regions,
                           cells = cells, keep_orient = character())
  disc = discover_clones(X)
  expect_gt(length(unique(disc$leaves)), 1)
  expect_true(all(table(disc$leaves) >= 20))
  cd = contrast_data(inp$features, vs$test, inp$atoms, inp$regions, cells)
  mc = merge_clones(inp$features, NULL, inp$atoms, inp$regions, disc$leaves, disc, rho = inp$null_fit$rho, cells = cells, cdata = cd)
  expect_length(unique(mc$clones), 1)
})

test_that("the split test keeps its budget on random splits of normal cells (5-run proxy)", {
  # Proxy for "< 1 false split in 20 runs": five random 60 / 60 splits of sim_norm_f singlets, each tested
  # with the full 200-permutation null (no early stop); at q = 0.999 a false split is expected 0.005 times.
  inp = norm_input()
  cd = contrast_data(inp$features, inp$votes, inp$atoms, inp$regions, inp$cells$cell)
  sing = which(!sim_norm_f$truth$cells$multiplet)
  res = rbindlist(lapply(1:5, function(r) {
    set.seed(100 + r)
    a = sample(sing, 60); b = sample(setdiff(sing, a), 60)
    split_test(cd, a, b, inp$null_fit$rho, n_perm = 200, seed = r, early = FALSE)
  }))
  expect_equal(sum(res$differ), 0)
  expect_true(all(res$n_perm == 200))
})

# ---- perm_decision: the merge / split test's decision rule ---------------------------------------------------------

test_that("perm_decide: false-different rate on exchangeable groups ~6 / 1,000 (quantile of 200) vs <= 1 / 1,000 (p-value)", {
  # exchangeable groups: |difference of means| of two random halves of 40 iid values, the null by permuting the labels
  n_rep = 20000
  run = function(dec, seed) {
    set.seed(seed)
    rbindlist(lapply(seq_len(n_rep), function(r) {
      x = stats::rnorm(40)
      st = function(g) abs(sum(x[g]) - sum(x[-g]))
      d = scphasedcn:::perm_decide(st(1:20), function() st(sample.int(40, 20)), q = 0.999, n_perm = 200, perm_decision = dec)
      list(differ = d$differ, n_perm = d$n_perm, p = d$p)
    }))
  }
  rq = run("quantile", 21); rp = run("pvalue", 22)
  ci = function(k) signif(1000 * stats::binom.test(k, n_rep)$conf.int, 3)
  message(sprintf("exchangeable groups, %d tests: quantile (n_perm 200, q 0.999) %d false differ = %.2f / 1,000 (95%% CI %s); ",
                  n_rep, sum(rq$differ), 1000 * mean(rq$differ), paste(ci(sum(rq$differ)), collapse = "-")),
          sprintf("pvalue %d = %.2f / 1,000 (95%% CI %s); mean permutations %.1f / %.1f",
                  sum(rp$differ), 1000 * mean(rp$differ), paste(ci(sum(rp$differ)), collapse = "-"), mean(rq$n_perm), mean(rp$n_perm)))
  # the quantile of 200 permutations: (n + 1 - h) / (n + 1) = 1.199 / 201 ~ 6 / 1,000, outside a binomial CI around 1 / 1,000
  expect_gt(mean(rq$differ), 0.0035); expect_lt(mean(rq$differ), 0.0085)
  expect_lt(stats::binom.test(sum(rq$differ), n_rep, 0.001, alternative = "greater")$p.value, 1e-6)
  # the p-value decision keeps 1 / 1,000 (exact under exchangeability)
  expect_gt(stats::binom.test(sum(rp$differ), n_rep, 0.001, alternative = "greater")$p.value, 0.001)
  expect_lte(mean(rp$differ), 0.002)
  # every "different" ran all 999 permutations at p = 1 / 1,000; "same" stops early (mean ~ log(999) permutations)
  expect_true(all(rp[differ == TRUE, n_perm] == 999))
  expect_true(all(abs(rp[differ == TRUE, p] - 0.001) < 1e-12))
  expect_true(all(rp[differ == FALSE, p] > 0.001))
  expect_lt(mean(rp$n_perm), 30)
  expect_identical(scphasedcn:::perm_min(0.999), 999L)
  expect_identical(scphasedcn:::perm_min(0.99), 99L)
  expect_identical(scphasedcn:::perm_min(0.95), 19L)
})

test_that("perm_decide: the p-value early stop is exact (same verdict as the full run) and the quantile path is unchanged", {
  for (r in 1:200) {
    set.seed(r)
    obs = stats::rnorm(1, 2.5)
    nul = stats::rnorm(1200)
    i = 0; draw = function() { i <<- i + 1; nul[i] }
    e = scphasedcn:::perm_decide(obs, draw, 0.999, 200, TRUE, "pvalue")
    i = 0
    f = scphasedcn:::perm_decide(obs, draw, 0.999, 200, FALSE, "pvalue")
    expect_identical(e$differ, f$differ)
    expect_equal(f$n_perm, 999)
    expect_equal(f$p, (1 + sum(nul[1:999] >= obs)) / 1000)
    expect_identical(f$differ, f$p <= 0.001)
    # n_perm above the minimum: up to max(n_perm, 999) permutations, p <= 1 - q
    i = 0
    g = scphasedcn:::perm_decide(obs, draw, 0.999, 1200, FALSE, "pvalue")
    expect_equal(g$n_perm, 1200)
    expect_identical(g$differ, (1 + sum(nul >= obs)) / 1201 <= 0.001)
  }
  # quantile: the loop of split_test before perm_decision (type-7 quantile of the draws, early stops)
  set.seed(5); nul = stats::rexp(200); obs = 4
  i = 0; draw = function() { i <<- i + 1; nul[i] }
  d = scphasedcn:::perm_decide(obs, draw, 0.999, 200, FALSE, "quantile")
  expect_equal(d$thr, unname(stats::quantile(nul, 0.999, type = 7)))
  expect_identical(d$differ, obs > d$thr)
})

test_that("split_test(perm_decision = \"pvalue\") separates sim_tum2's clones and keeps one clone together", {
  inp = tum2_input()
  cells = inp$cells$cell
  tr = sim_truth(sim_tum2)[mult == FALSE]
  cd = contrast_data(inp$features, inp$votes, inp$atoms, inp$regions, cells)
  ia = match(tr[type == "cloneA", cell], cells); ib = match(tr[type == "cloneB", cell], cells)
  t0 = Sys.time()
  dd = split_test(cd, ia, ib, inp$null_fit$rho, perm_decision = "pvalue")
  tm = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  expect_true(dd$differ)
  expect_equal(dd$n_perm, 999)
  expect_equal(dd$p, 1 / 1000)
  set.seed(3); sa = sample(ia, 30); sb = sample(setdiff(ia, sa), 30)
  ss = split_test(cd, sa, sb, inp$null_fit$rho, perm_decision = "pvalue")
  expect_false(ss$differ)
  expect_lt(ss$n_perm, 999)
  expect_gt(ss$p, 0.001)
  message(sprintf("split_test pvalue: cloneA vs cloneB %d permutations in %.1f s; a random split of cloneA stops after %d", dd$n_perm, tm, ss$n_perm))
  # the default is the quantile decision without a p column
  dq = split_test(cd, ia, ib, inp$null_fit$rho)
  expect_false("p" %in% names(dq))
  expect_identical(dq, split_test(cd, ia, ib, inp$null_fit$rho, perm_decision = "quantile"))
})

test_that("caller_tumour(perm_decision = \"pvalue\") on sim_tum2: same clones, p in the tests, the run notes count the tests", {
  r0 = tum2_result()
  t0 = Sys.time()
  r1 = clone_cached("clones_tum2_pvalue", function() caller_tumour(tum2_input(), opts = list(perm_decision = "pvalue"), verbose = FALSE))
  tm = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr = sim_truth(sim_tum2)
  a = merge(r1$assignments, tr, by = "cell")[mult == FALSE]
  expect_gte(clone_ari(a$clone, a$type), 0.9)
  tt = r1$null_fit$clones$tests
  expect_true("p" %in% names(tt))
  expect_true(all(tt[differ == TRUE, n_perm] >= 999))
  expect_true(all(tt$differ == (tt$p <= 0.001)))
  expect_false("p" %in% names(r0$null_fit$clones$tests))
  nt = grep("^split tests", r1$notes, value = TRUE)
  expect_length(nt, 1)
  expect_match(nt, "perm_decision pvalue")
  expect_match(grep("^split tests", r0$notes, value = TRUE), "perm_decision quantile.*6.0 / 1,000")
  message("sim_tum2 perm_decision pvalue: ", round(tm), " s; ", nt, "\n  default: ", grep("^split tests", r0$notes, value = TRUE),
          "\n  global tests (pvalue): ", nrow(tt), ", permutations ", sum(tt$n_perm), " (quantile: ", nrow(r0$null_fit$clones$tests), ", ",
          sum(r0$null_fit$clones$tests$n_perm), ")")
})

# ---- sim_tum3: the confirmation null refitted on the reference pool, matched null --------------------------------
# sim_tum2's clones at half the depth with an attenuated allele signal (clone a ~ 0.7). The input rho is inflated 10x,
# as in tumour samples, where call_scn fits it on copy-number-flat cells that partly carry tumour haplotypes.
# 5 Mb control pieces: sim_tum3 has two (chr2 0-10 Mb, chr8 0-15 Mb), none at the default 20 Mb.

tum3_input = function() clone_cached("clones_tum3_input", function() clone_sim_input(sim_tum3))
tum3_inflated = function() clone_cached("clones_tum3_inflated", function() {
  x = tum3_input()
  x$null_fit$rho = 10 * x$null_fit$rho
  x
})
tum3_run = function(name, opts = list()) clone_cached(paste0("clones_tum3_", name), function()
  caller_tumour(tum3_inflated(), opts = c(list(control_min_len = 5e6), opts), verbose = FALSE))
tum3_map = function(res) {
  a = merge(res$assignments, sim_truth(sim_tum3), by = "cell")[mult == FALSE]
  a[, .N, by = .(clone, type)][order(-N)][!duplicated(clone)]
}
tum3_conf = function(res, type_, ch) {
  x = merge(res$calls[level == "cell" & chrom == ch], sim_truth(sim_tum3), by = "cell")[type == type_ & mult == FALSE]
  mean(x$source == "confirmed")
}

test_that("sim_tum3: attenuated clone allele signal; ref_null = \"none\" keeps the input fits; the generic null has no hyp rows", {
  inp = tum3_inflated()
  res = tum3_run("generic", list(confirm_null = "generic"))                        # default is now matched
  map = tum3_map(res)
  kA = map[type == "cloneA", clone]
  expect_length(kA, 1)
  expect_gte(res$clones[clone == kA, a], 0.55)
  expect_lte(res$clones[clone == kA, a], 0.85)
  expect_equal(res$segments[clone == kA & chrom == "chr21" & start <= 10e6 & end > 10e6, state], "2:1")
  expect_equal(res$null_fit$rho_confirm, inp$null_fit$rho)
  expect_identical(res$null_fit$tpar_confirm, inp$null_fit$tpar)
  expect_null(res$null_fit$confirm_null)
  expect_false("hyp" %in% names(res$thresholds))
})

test_that("ref_null refits rho on the reference pool; the matched null confirms more clone gains within budget", {
  inp = tum3_input()
  res0 = tum3_run("generic", list(confirm_null = "generic"))
  res = tum3_run("rho_matched", list(ref_null = "rho", confirm_null = "matched"))
  # the refit is after assignment: clones and assignments are those of the default run
  expect_equal(res$assignments[, .(cell, clone, posterior)], res0$assignments[, .(cell, clone, posterior)])
  expect_equal(res$segments, res0$segments)
  # rho refitted on the pool's clone-oriented piece votes: the 10x inflation is gone (the simulator has no allelic
  # overdispersion beyond the binomial, so the clean fit, ~5e-4, and the refit are both ~0 and their ratio is noise)
  rc = res$null_fit$rho_confirm
  expect_true(is.finite(rc))
  expect_lt(rc, 2 * inp$null_fit$rho)
  expect_true(is.finite(res$null_fit$rho_pool_features))
  expect_lt(res$null_fit$rho_pool_features, 2 * inp$null_fit$rho)
  expect_true(any(grepl("confirmation rho refitted", res$notes)))
  # matched thresholds: one per clone hypothesis (x stratum), every cell call has a finite threshold
  check_thresholds(res$thresholds)
  expect_true("hyp" %in% names(res$thresholds))
  cc = res$calls[level == "cell"]
  expect_true(all(is.finite(cc$thr)))
  expect_equal(cc$called, cc$score > cc$thr)
  hy = res$thresholds$hyp
  expect_setequal(sub("\\|.*", "", hy), unique(res$segments[label != "neutral", clone]))
  # chr21 2:1 (f ~ 0.61) in cloneA: more confirmed than with the generic null on the inflated rho; other events kept
  expect_gte(tum3_conf(res, "cloneA", "chr21") - tum3_conf(res0, "cloneA", "chr21"), 0.05)
  for (ch in c("chr22", "chr2", "chr1")) expect_gte(tum3_conf(res, "cloneB", ch), 0.95)
  expect_gte(tum3_conf(res, "cloneA", "chr8"), 0.95)
  # held-out (fold-CV) false confirmations of the matched null
  cn = res$null_fit$confirm_null
  expect_true(all(c("hyp", "kind", "n_null", "n_called", "per_1000") %in% names(cn)))
  expect_setequal(cn[hyp != "all", unique(hyp)], hy)
  expect_gt(cn[hyp == "all" & kind == "control", n_null], 1000)
  expect_lte(cn[hyp == "all" & kind == "all", per_1000], 3)
  expect_lte(cn[hyp == "all" & kind == "control", per_1000], 5)
  # null_fit.tsv carries the confirmation fits
  nf = null_fit_table(res$null_fit)
  expect_equal(nf[kind == "allelic_rho_confirm", rho], rc)
  expect_gt(nrow(nf[kind == "depth_t_confirm"]), 0)
})

# ---- profile_orient: clone profile orientation per window x phase set ------------------------------------------------

truth_profiles = function(sim, inp, orient) {
  tr = sim_truth(sim)[mult == FALSE]
  members = tr[, .(cell, clone = type)]
  clone_profiles(inp$win, inp$votes, inp$atoms, members, ref_cells = members[clone == "normal", cell], cells = inp$cells$cell,
                 rho = inp$null_fit$rho, allelic_chroms = setdiff(unique(inp$regions$chrom), c("chrX", "chrY")),
                 haploid = c("chrX", "chrY"), orient = orient)
}

test_that("clone_profiles(orient = \"window\") keeps two opposite-haplotype events inside one phase set; \"ps\" loses them", {
  inp = clone_cached("clones_cw_input", function() clone_sim_input(sim_tum_cw))
  expect_equal(sim_tum_cw$snps[chrom == "chr2", uniqueN(ps)], 1)        # one phase set over both chr2 events
  p0 = truth_profiles(sim_tum_cw, inp, "ps"); p1 = truth_profiles(sim_tum_cw, inp, "window")
  expect_identical(truth_profiles(sim_tum_cw, inp, c("ps", "window")), p0)          # "ps" is the default
  fa = function(p, lo, hi) p[clone == "cloneA" & chrom == "chr2" & start >= lo & start < hi, sum(nA) / sum(n_hap)]
  # one orientation for the chromosome: the arms point opposite ways (one of them below 0.5)
  expect_lt(min(fa(p0, 0, 16e6), fa(p0, 16e6, 40e6)), 0.45)
  # per window: both arms oriented to their own major haplotype (2:1 at f ~ 2/3)
  expect_gt(fa(p1, 0, 16e6), 0.6); expect_gt(fa(p1, 16e6, 40e6), 0.6)
  expect_equal(p1[, .(clone, window, log2)], p0[, .(clone, window, log2)])         # depth untouched
  s0 = segment_clone(fit_clone_states(p0)); s1 = segment_clone(fit_clone_states(p1))
  at = function(sg, ch, pos) sg[clone == "cloneA" & chrom == ch & start <= pos & end > pos]
  message("sim_tum_cw cloneA chr2 (p 2:1, q 1:2): ps f p / q arm ", round(fa(p0, 0, 16e6), 3), " / ", round(fa(p0, 16e6, 40e6), 3),
          ", state p / q ", at(s0, "chr2", 8e6)$state, " / ", at(s0, "chr2", 30e6)$state, "; window f ", round(fa(p1, 0, 16e6), 3),
          " / ", round(fa(p1, 16e6, 40e6), 3), ", state ", at(s1, "chr2", 8e6)$state, " / ", at(s1, "chr2", 30e6)$state)
  # "ps": the imbalanced run fails min_excess (pooled f ~ 0.5) and is reset to neutral although the depth is +0.58 log2
  expect_true("neutral" %in% c(at(s0, "chr2", 8e6)$label, at(s0, "chr2", 30e6)$label))
  expect_equal(at(s1, "chr2", 8e6)$state, "2:1"); expect_equal(at(s1, "chr2", 30e6)$state, "2:1")
  expect_equal(at(s0, "chr22", 10e6)$state, "1:0"); expect_equal(at(s1, "chr22", 10e6)$state, "1:0")
})

test_that("caller_tumour(profile_orient = \"window\") calls both chr2 gains of sim_tum_cw", {
  inp = clone_cached("clones_cw_input", function() clone_sim_input(sim_tum_cw))
  r1 = caller_tumour(inp, opts = list(profile_orient = "window"), verbose = FALSE)
  r0 = caller_tumour(inp, verbose = FALSE)
  tr = sim_truth(sim_tum_cw)
  a1 = merge(r1$assignments, tr, by = "cell")[mult == FALSE]
  expect_gte(clone_ari(a1$clone, a1$type), 0.9)
  kA = a1[type == "cloneA", .N, by = clone][order(-N), clone][1]
  sg = r1$segments[clone == kA]
  expect_true(all(sg[chrom == "chr2", label] == "gain"))
  expect_equal(sg[chrom == "chr22" & start <= 10e6 & end > 10e6, state], "1:0")
  # per-cell confirmation is not fixed by this option: it orients the votes per piece x phase set (trip_votes /
  # oriented_counts), and the chr2 piece (one 2:1 segment over both arms) holds both opposite-haplotype gains in one
  # phase set, so few cloneA cells confirm chr2 (reported, not asserted)
  cc = merge(r1$calls[level == "cell" & chrom == "chr2"], tr, by = "cell")[mult == FALSE & type == "cloneA"]
  expect_gt(nrow(cc), 0)
  # the default ("ps") calls no chr2 gain in any clone
  expect_false(any(r0$segments[chrom == "chr2", label] != "neutral"))
  a0 = merge(r0$assignments, tr, by = "cell")[mult == FALSE]
  message("sim_tum_cw caller_tumour: ps ARI ", round(clone_ari(a0$clone, a0$type), 2), ", non-neutral segments ",
          paste(r0$segments[label != "neutral", paste(clone, chrom, state)], collapse = ", "), "; window ARI ",
          round(clone_ari(a1$clone, a1$type), 2), ", ", paste(r1$segments[label != "neutral", paste(clone, chrom, start, end, state)], collapse = ", "),
          "; cloneA cells with chr2 confirmed ", round(mean(cc$source == "confirmed"), 3))
})

test_that("profile_orient \"window\" on sim_tum2 (short phase blocks): near-identical profiles and the same clone states", {
  inp = tum2_input()
  p0 = truth_profiles(sim_tum2, inp, "ps"); p1 = truth_profiles(sim_tum2, inp, "window")
  m = merge(p0[, .(clone, window, f0 = f, n0 = n_hap)], p1[, .(clone, window, f1 = f, n1 = n_hap)], by = c("clone", "window"))
  # cost of the smaller orientation units: blocks that cross a window boundary are oriented per window part
  ev = m[clone != "normal" & abs(f0 - 0.5) > 0.1]
  message(sprintf("sim_tum2 profiles ps vs window: mean |f diff| %.4f (all windows), imbalanced windows mean f %.4f vs %.4f (oriented excess ratio %.3f), votes %.0f vs %.0f",
                  mean(abs(m$f0 - m$f1), na.rm = TRUE), mean(ev$f0), mean(ev$f1), mean(ev$f1 - 0.5) / mean(ev$f0 - 0.5), sum(m$n0), sum(m$n1)))
  expect_lt(mean(abs(m$f0 - m$f1), na.rm = TRUE), 0.02)
  expect_gt(mean(ev$f1 - 0.5) / mean(ev$f0 - 0.5), 0.9)
  r1 = clone_cached("clones_tum2_window", function() caller_tumour(inp, opts = list(profile_orient = "window"), verbose = FALSE))
  r0 = tum2_result()
  expect_equal(r1$assignments$clone, r0$assignments$clone)
  expect_equal(r1$segments[, .(clone, chrom, start, end, state)], r0$segments[, .(clone, chrom, start, end, state)])
  expect_identical(r1$params$tumour$profile_orient, "window")
})

# ---- confirm_sided: one-sided per-cell confirmation (clone-oriented votes) ---------------------------------------------

test_that("b_score_hyp: two-sided confirms a cell whose opposite homologue is imbalanced; confirm_sided \"one\" does not", {
  # votes oriented by the clone (haplotype A = the clone's major haplotype): 40 / 40 vs 0 / 40, depth at the loss shift
  x = data.table(cell = c("same", "opposite"), region = "r", class = "arm", stratum = "x", z = -5, sd_null = 0.2,
                 nA = c(40, 0), n_hap = 40)
  p = caller_defaults()
  s2 = scphasedcn:::b_score_hyp(x, "loss", -1, 0.95, NULL, 0.01, p, cn = 1L)
  s1 = scphasedcn:::b_score_hyp(x, "loss", -1, 0.95, NULL, 0.01, p, cn = 1L, sided = "one")
  expect_identical(scphasedcn:::b_score_hyp(x, "loss", -1, 0.95, NULL, 0.01, p, cn = 1L, sided = "two"), s2)
  expect_equal(s2$s_comb[1], s2$s_comb[2])
  expect_gt(s2$s_comb[2], 15)
  expect_equal(s1$s_comb[1], s2$s_comb[1])
  expect_lt(s1$s_comb[2], -20)
  # b_score (cell_state_scores rows) passes it through
  rows = data.table(cell = c("same", "opposite"), region = "r", sd_class = "arm", stratum = "x", z = -5, sd_null = 0.2,
                    nA_o = c(40, 0), n_o = 40, label = "loss", cn = 1L, c = 1, a = 1, f_state = 0.95, base_ = 2L)
  b1 = scphasedcn:::b_score(rows, NULL, 0.01, p, sided = "one")
  expect_equal(b1$s_comb, s1$s_comb)
  expect_equal(scphasedcn:::b_score(rows, NULL, 0.01, p)$s_comb, s2$s_comb)
  expect_identical(scphasedcn:::check_clone_opts(list())$confirm_sided, "two")
  expect_error(scphasedcn:::check_clone_opts(list(confirm_sided = "both")), "should be one of")
})

# per-event confirmed rate of the tumour cells (singlets of the truth type) and false confirmations: confirmed calls on
# singlets whose truth type does not carry the event (normal cells, the other clone)
tum_conf_table = function(res, sim) {
  tr = sim_truth(sim)
  x = merge(res$calls[level == "cell"], tr, by = "cell")[mult == FALSE]
  ev = data.table(chrom = c("chr8", "chr21", "chr22", "chr2", "chr1"), carrier = c("cloneA|cloneB", "cloneA", "cloneB", "cloneB", "cloneB"))
  x = merge(x, ev, by = "chrom")
  x[, carries := mapply(grepl, carrier, type)]
  x[, .(n_true = sum(carries), conf_true = mean(source[carries] == "confirmed"), n_false = sum(!carries),
        conf_false = sum(source[!carries] == "confirmed")), by = chrom][order(chrom)]
}

test_that("confirm_sided \"one\": assignments and segments unchanged, events still confirmed, control false confirmations at the null rate", {
  for (nm in c("tum2", "tum3")) {
    r0 = if (nm == "tum2") tum2_result() else tum3_run("matched")
    r1 = if (nm == "tum2") clone_cached("clones_tum2_one", function() caller_tumour(tum2_input(), opts = list(confirm_sided = "one"), verbose = FALSE)) else
      tum3_run("one", list(confirm_sided = "one"))
    sim = if (nm == "tum2") sim_tum2 else sim_tum3
    # confirmation is after assignment: clones, assignments and segments are those of the default run
    expect_identical(r1$assignments, r0$assignments)
    expect_identical(r1$segments, r0$segments)
    expect_identical(r1$clones, r0$clones)
    expect_identical(r1$params$tumour$confirm_sided, "one")
    expect_true(any(grepl("one-sided", r1$notes)))
    c0 = tum_conf_table(r0, sim); c1 = tum_conf_table(r1, sim)
    message(nm, " confirmed rate per event (two-sided / one-sided):\n",
            paste(capture.output(print(merge(c0, c1, by = c("chrom", "n_true", "n_false"), suffixes = c(".two", ".one")))), collapse = "\n"))
    # true events stay confirmed (within 10 points of the two-sided rate) ...
    m = merge(c0, c1, by = "chrom", suffixes = c(".two", ".one"))
    expect_true(all(m$conf_true.one >= m$conf_true.two - 0.1))
    # ... and the held-out false confirmations of the matched null stay at the budget
    cn = r1$null_fit$confirm_null
    if (!is.null(cn) && nrow(cn[hyp == "all" & kind == "control"])) {
      message(nm, " held-out false confirmations per 1,000 (two / one): control ",
              signif(r0$null_fit$confirm_null[hyp == "all" & kind == "control", per_1000], 3), " / ",
              signif(cn[hyp == "all" & kind == "control", per_1000], 3), ", pool ",
              signif(r0$null_fit$confirm_null[hyp == "all" & kind == "pool", per_1000], 3), " / ",
              signif(cn[hyp == "all" & kind == "pool", per_1000], 3))
      expect_lte(cn[hyp == "all" & kind == "control", per_1000], 5)
      expect_lte(cn[hyp == "all" & kind == "all", per_1000], 3)
    }
    cc = r1$calls[level == "cell"]
    expect_equal(cc$called, cc$score > cc$thr)
    check_calls(r1$calls)
  }
})

# ---- sim_tum4: the disagreement split (split_pass = "disagree") ------------------------------------------------------
# sim_tum2's clones plus cloneC (cloneB with chr21 1:0 instead of 2:0, same haplotype: a depth-only difference: a loss
# subclone inside a CN-LOH clone) and cloneD (cloneA with chr21 3:0 instead of 2:1: an LOH gain). At the defaults the
# allele-driven discovery separates cloneD from cloneA but keeps cloneB and cloneC together.

tum4_input = function() clone_cached("clones_tum4_input", function() clone_sim_input(sim_tum4))
tum4_run = function(name, opts = list()) clone_cached(paste0("clones_tum4_", name), function()
  caller_tumour(tum4_input(), opts = opts, verbose = FALSE))
tum4_sing = function(res) merge(res$assignments, sim_truth(sim_tum4), by = "cell")[mult == FALSE]
tum4_chr21 = function(res) {
  a = tum4_sing(res)
  map = a[, .N, by = .(type, clone)][order(-N)][!duplicated(type)]
  map[, state := vapply(clone, function(k) res$segments[clone == k & chrom == "chr21" & start <= 15e6 & end > 15e6, state][1], "")]
  stats::setNames(map$state, map$type)
}

test_that("split_alt_states: one copy apart; depth-only alternatives only next to <= 2 copies", {
  expect_setequal(split_alt_states("2:0")$state, c("1:0", "3:0", "1:1", "2:1"))
  expect_true(all(c("3:0", "4:0") %in% split_alt_states("3:1")$state))
  expect_false("2:1" %in% split_alt_states("3:1")$state)             # same depth step between amplified states, f 0.67 vs 0.75
  expect_false("3:0" %in% split_alt_states("4:0")$state)             # amplicon depth: 4:0 vs 3:0 is not a subclone
  expect_true("2:0" %in% split_alt_states("1:0")$state)
  expect_equal(nrow(split_alt_states("9:9")), 0)
  expect_identical(clone_defaults()$split_pass, FALSE)
})

test_that("split_disagree splits merged clones on the chr21 state they disagree on, not on the shared amplicon", {
  inp = tum4_input()
  cells = inp$cells$cell
  tr = sim_truth(sim_tum4)[mult == FALSE]
  # truth-merged clones: K1 = cloneB + cloneC (chr21 2:0 / 1:0), K2 = cloneA + cloneD (chr21 2:1 / 3:0), N = normal
  grp = c(cloneB = "K1", cloneC = "K1", cloneA = "K2", cloneD = "K2", normal = "N")
  members = tr[, .(cell, clone = grp[type])]
  ref = members[clone == "N", cell]
  vs = split_votes(inp$votes, 0.5, seed = 1); ws = split_win(inp$win, 0.5, seed = 1)
  cd = contrast_data(half_depth_features(clone_input(inp), ws$test, c("arm", "20Mb")), vs$test, inp$atoms, inp$regions, cells)
  prof = clone_profiles(inp$win, inp$votes, inp$atoms, members, ref_cells = ref, cells = cells, rho = inp$null_fit$rho,
                        allelic_chroms = setdiff(unique(inp$regions$chrom), c("chrX", "chrY")), haploid = c("chrX", "chrY"))
  fit = fit_clone_states(prof)
  sg = segment_clone(fit)
  sp = split_disagree(members[clone != "N"], sg[clone != "N"], fit$par, inp, ws$disc, vs$disc, cd, win_test = ws$test,
                      votes_test = vs$test, ref_cells = ref, seed = 1)
  tt = sp$tests
  expect_true(all(c("clone", "region", "state", "alt", "contra", "sep", "p", "differ", "state_a", "state_b", "split", "outcome") %in% names(tt)))
  expect_equal(sort(tt[split == TRUE, clone]), c("K1", "K2"))
  expect_true(all(tt[split == TRUE, chrom] == "chr21"))
  expect_false(any(tt[chrom == "chr8", split]))                     # the shared 3:1 amplicon is never split
  # each split separates the two truth clones (the proposal is made on the discovery halves: part purity >= 0.8, the
  # caller's EM round after the pass reassigns on all molecules) and the re-fitted parts differ on chr21
  m = merge(sp$members, tr[, .(cell, type)], by = "cell")
  pur = m[clone != "N", .(purity = max(table(type)) / .N, type = names(which.max(table(type)))), by = clone]
  expect_equal(nrow(pur), 4)
  expect_true(all(pur$purity >= 0.8))
  expect_setequal(pur$type, c("cloneA", "cloneB", "cloneC", "cloneD"))
  expect_true(all(tt[split == TRUE, state_a != state_b]))
  expect_setequal(c(tt[clone == "K1" & split == TRUE, c(state_a, state_b)]), c("2:0", "1:0"))
  # ranked by contradiction: the applied candidates are the clones' most contradicted segments
  expect_true(all(tt[split == TRUE, rank] <= 2))
})

test_that("caller_tumour(split_pass = \"disagree\") separates the chr21 subclones of sim_tum4; the default is unchanged", {
  r0 = tum4_run("default")
  r1 = tum4_run("disagree", list(split_pass = "disagree"))
  # the default: cloneB and cloneC (depth-only chr21 difference) share a clone, no split pass ran
  s0 = tum4_sing(r0)
  expect_lt(clone_ari(s0$clone, s0$type), 0.9)
  expect_null(r0$null_fit$clones$split_tests)
  r0f = tum4_run("false", list(split_pass = FALSE))
  for (nm in c("assignments", "segments", "calls", "clones", "thresholds", "null_fit")) expect_identical(r0f[[nm]], r0[[nm]])
  # disagreement split: both pairs separate, with their chr21 states
  s1 = tum4_sing(r1)
  expect_gte(clone_ari(s1$clone, s1$type), 0.9)
  st = tum4_chr21(r1)
  expect_equal(unname(st[c("cloneA", "cloneB", "cloneC", "cloneD")]), c("2:1", "2:0", "1:0", "3:0"))
  expect_equal(length(unique(s1[, .N, by = .(type, clone)][order(-N)][!duplicated(type), clone])), 5)
  tt = r1$null_fit$clones$split_tests
  expect_true(any(tt$split))
  expect_true(all(tt[split == TRUE, chrom] == "chr21"))
  expect_true(all(tt[split == TRUE, separate] %in% TRUE))
  # the doublet-like rate does not rise
  expect_lte(mean(r1$assignments$doublet_like), mean(r0$assignments$doublet_like) + 0.02)
  check_segments(r1$segments)
  check_calls(r1$calls)
  check_thresholds(r1$thresholds)
})

# ---- coverage-free allele scale, same-state posteriors, drift-aware split test -------------------------------------

test_that("assign_cells(same_groups): two clones with the same states no longer split each cell's posterior", {
  inp0 = tum2_input(); inp = clone_input(inp0)
  r0 = tum2_result()
  sg = copy(r0$segments); cl_ = copy(r0$clones); par0 = copy(r0$null_fit$clones$par)
  k = cl_[normal == FALSE & doublet == FALSE][order(-n_members), clone][1]
  sg2 = rbind(sg, copy(sg[clone == k])[, clone := "Kdup"])
  m = copy(r0$assignments)[confident == TRUE, .(cell, clone)]
  set.seed(1)
  half = m[clone == k, sample(cell, floor(.N / 2))]
  m[cell %in% half, clone := "Kdup"]
  par = rbind(par0, copy(par0[clone == k])[, clone := "Kdup"])
  tr = stats::setNames(rep("T1", length(unique(sg2$clone))), unique(sg2$clone))
  grp = stats::setNames(names(tr), names(tr)); grp["Kdup"] = k      # the twins form one group
  pieces = segment_pieces(sg2, inp$regions)
  fp = list(features = region_features(inp, pieces), pieces = pieces, votes = inp$votes, atoms = inp$atoms, cells = inp$cells$cell)
  a0 = assign_cells(fp, m, sg2, inp$tpar, inp$rho, par = par, params = caller_defaults(), trunks = tr)
  a1 = assign_cells(fp, m, sg2, inp$tpar, inp$rho, par = par, params = caller_defaults(), trunks = tr, same_groups = grp)
  kc = a0[clone %in% c(k, "Kdup"), cell]
  expect_lt(mean(a0[cell %in% kc, posterior] >= 0.9), 0.2)        # the posterior is split between the twins
  expect_gt(mean(a1[cell %in% kc, posterior] >= 0.9), 0.8)        # summed over the same-state group
  expect_equal(a1$post_clone, a0$posterior)
  expect_identical(a1$clone, a0$clone)
  expect_false("post_clone" %in% names(a0))
})
