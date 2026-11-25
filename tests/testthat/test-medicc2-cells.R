# test-medicc2-cells.R: per-cell states on the MEDICC2 pieces (R/medicc2_cells.R) on a hand-made male sample: chr1
# (60 Mb), chr2 (40 Mb), chrX (40 Mb). Pieces: chr1 0-30 Mb, chr1 30-60 Mb (P), chr2, chrX. Clones on P: K10 1:0, K20 2:0,
# K21 2:1, K01 0:1 (the opposite haplotype lost), N0 1:1; K20 also gains chrX (2:0, haploid). Cells are simulated with
# known states (depth: Poisson window counts; votes: fake_votes with random block flips; the block orientation table
# holds the true flips).

mc_toy = function(n_each = 6L, depth_mb = 300, mol_per_atom = 0.5, seed = 7) {
  set.seed(seed)
  sizes = c(chr1 = 60e6, chr2 = 40e6, chrX = 40e6)
  bins = make_bins(sizes, 1e5)
  regions = fake_regions(bins, tiles = c("5Mb" = 5e6, "10Mb" = 1e7), sex = "male")
  atoms = data.table(atom = seq_len(nrow(bins)), chrom = bins$chrom, start = bins$start, end = bins$end, bin = seq_len(nrow(bins)))
  types = c("ref", "bal", "s10", "s20", "s21", "s01")
  lab = c(rep("ref", 20L), rep(types[-1], each = n_each))
  cells = sprintf("c%03d_%s", seq_along(lab), lab)
  ev = data.table(clone = c("s10", "s20", "s21", "s01"), chrom = "chr1", start = 30e6, end = 60e6,
                  cnA = c(1L, 2L, 2L, 0L), cnB = c(0L, 0L, 1L, 1L))
  votes = fake_votes(length(lab), atoms[chrom != "chrX"], lab, ev, flip_prob = 0.5, seed = seed, block_atoms = 10L,
                     mol_per_atom = mol_per_atom)
  # depth: 1 Mb windows, x ~ Poisson(e * cn / normal cn); s20 also carries 2 copies of chrX
  wt = rbindlist(lapply(names(sizes), function(ch) data.table(chrom = ch, start = seq(0, sizes[[ch]] - 1e6, by = 1e6))))
  wt[, `:=`(end = start + 1e6, window = paste0(chrom, ":", format(start, scientific = FALSE, trim = TRUE)))]
  cn_true = function(l) {
    cn = fifelse(wt$chrom == "chrX", 1, 2)
    P = wt$chrom == "chr1" & wt$start >= 30e6
    cn[P] = switch(l, s10 = 1, s20 = 2, s21 = 3, s01 = 1, 2)
    if (l == "s20") cn[wt$chrom == "chrX"] = 2
    cn
  }
  e = matrix(depth_mb, nrow(wt), length(lab), dimnames = list(wt$window, cells))
  x = e
  for (j in seq_along(lab)) x[, j] = stats::rpois(nrow(wt), e[, j] * cn_true(lab[j]) / fifelse(wt$chrom == "chrX", 1, 2))
  # null SD per cell x sd_class: Poisson SD of the class size, x 1.5
  sdc = regions[, .(w = stats::median(end - start)), by = sd_class]
  long = CJ(cell = cells, region = regions$region)
  long[regions, sd_class := i.sd_class, on = "region"]
  long[, sd_null := 1.5 / (log(2) * sqrt(depth_mb * sdc$w[match(sd_class, sdc$sd_class)] / 1e6))]
  cm = data.table(cell = cells, n_frag = depth_mb * 140, stratum = "20-50k")
  features = list(long = long, regions = regions, cells = cm, votes = votes, atoms = atoms, win = list(windows = wt, x = x, e = e),
                  rho = 0.005)
  # pieces: clone x piece oriented states
  pc = data.table(piece = c("p1a", "p1b", "p2", "pX"), chrom = c("chr1", "chr1", "chr2", "chrX"), start = c(0, 30e6, 0, 0),
                  end = c(30e6, 60e6, 40e6, 40e6), ploidy = c(2L, 2L, 2L, 1L))
  st = data.table(clone = c("N0", "K10", "K20", "K21", "K01"), a1 = c(1L, 1L, 2L, 2L, 0L), b1 = c(1L, 0L, 0L, 1L, 1L),
                  xa = c(1L, 1L, 2L, 1L, 1L), n_assigned = c(50, 10, 10, 10, 30), c = c(0.7, 1, 1, 1, 1), a = c(1, 1, 1, 1, 0.8))
  pieces = rbindlist(lapply(seq_len(nrow(pc)), function(i) {
    p = pc[i]
    cbind(p, st[, .(clone, cn_a = if (p$piece == "p1b") a1 else if (p$piece == "pX") xa else 1L,
                    cn_b = if (p$piece == "p1b") b1 else if (p$piece == "pX") 0L else 1L)])
  }))
  # block orientation: the true flips (sign +1 = phase-set haplotype A is cn_a = the simulated A)
  tb = copy(attr(votes, "blocks"))
  tb[, bstart := ps - 1]
  bl = pc[ploidy == 2L][tb, on = .(chrom, start <= bstart, end > bstart), nomatch = NULL,
                        .(piece, chrom, start = x.start, end = x.end, ps, ref = "K10", sign = fifelse(flip, -1L, 1L))]
  list(features = features, pieces = pieces, blocks = bl, cells = cells, lab = stats::setNames(lab, cells),
       clones = st[, .(clone, n_assigned, c, a)], ref = cells[lab == "ref"])
}

test_that("cell_piece_candidates: neutral + clone states, haploid pieces, priors", {
  x = mc_toy(n_each = 1L)
  cd = cell_piece_candidates(x$pieces)
  expect_setequal(cd[piece == "p1b", state], c("1:1", "1:0", "2:0", "2:1", "0:1"))
  expect_equal(cd[piece == "p1a", state], "1:1")
  expect_setequal(cd[piece == "pX", state], c("1:0", "2:0"))
  expect_equal(cd[piece == "pX" & neutral == TRUE, state], "1:0")
  expect_true(all(abs(cd[, sum(prior), by = piece]$V1 - 1) < 1e-12))
  expect_equal(cd[piece == "p1b", unique(prior)], 0.2)
  cf = cell_piece_candidates(x$pieces, x$clones, prior = "clone_freq")
  # 0:1 (K01, 30 cells) is more probable a priori than 1:0 (K10, 10 cells); 1:1 has N0 (50)
  expect_gt(cf[piece == "p1b" & state == "0:1", prior], cf[piece == "p1b" & state == "1:0", prior])
  expect_equal(cf[piece == "p1b"][which.max(prior), state], "1:1")
  expect_true(all(abs(cf[, sum(prior), by = piece]$V1 - 1) < 1e-12))
  # haploid pieces drop states with a second haplotype; a piece without 1:1 clones still gets the neutral candidate
  p2 = copy(x$pieces)[piece == "pX" & clone == "K21", cn_b := 1L][piece == "p2", `:=`(cn_a = 2L, cn_b = 1L)]
  c2 = cell_piece_candidates(p2, x$clones, prior = "clone_freq")
  expect_false("1:1" %in% c2[piece == "pX", state])
  expect_setequal(c2[piece == "p2", state], c("1:1", "2:1"))
  expect_gt(c2[piece == "p2" & state == "1:1", prior], 0)
  cn = cell_piece_candidates(x$pieces, prior = "neutral", p_neutral = 0.6)
  expect_equal(cn[piece == "p1b" & neutral == TRUE, prior], 0.6)
  expect_equal(cn[piece == "p1b" & neutral == FALSE, unique(prior)], 0.1)
  expect_equal(cn[piece == "p1a", prior], 1)
})

test_that("cell_piece_states recovers simulated 1:0 / 2:0 / 2:1 / 0:1 cells and keeps balanced cells 1:1", {
  x = mc_toy()
  r = cell_piece_states(x$features, x$pieces, x$blocks, params = caller_defaults(), clones = x$clones, ref_cells = x$ref)
  expect_s3_class(r, "cell_piece_states")
  s = r$states
  expect_equal(nrow(s), length(x$cells) * 4L)
  # global c / a: assigned-cell weighted median of the tumour clones (N0 excluded)
  expect_equal(r$params$c, 1)
  expect_equal(r$params$a, 0.8)
  s[, lab := x$lab[cell]]
  p = s[piece == "p1b"]
  truth = c(ref = "1:1", bal = "1:1", s10 = "1:0", s20 = "2:0", s21 = "2:1", s01 = "0:1")
  expect_true(all(p$state == truth[p$lab]))
  expect_true(all(p[lab %in% c("s10", "s20", "s21", "s01"), posterior] > 0.99))
  expect_true(all(p[lab %in% c("s10", "s20", "s21", "s01"), llr_neutral] > 10))
  expect_true(all(p[lab %in% c("ref", "bal"), llr_neutral] == 0))
  # the allele term distinguishes 1:0 from 0:1 (same depth), the depth term 1:0 from 2:0 (same fraction)
  expect_true(all(p[lab == "s10", f_obs] > 0.8) && all(p[lab == "s01", f_obs] < 0.2))
  expect_true(all(abs(p[lab == "s20", z]) < 2) && all(p[lab == "s10", z] < -2) && all(p[lab == "s21", z] > 2))
  # balanced pieces stay neutral in every cell; chrX: 2:0 only in s20 (haploid, depth only)
  expect_true(all(s[piece %in% c("p1a", "p2"), state] == "1:1"))
  expect_true(all(s[piece == "pX", state] == fifelse(s[piece == "pX", lab] == "s20", "2:0", "1:0")))
  expect_true(all(s[piece == "pX", ll_allele] == 0))
  expect_true(all(s[piece == "p1b" & lab != "ref", oN] > 50))
})

test_that("cell_piece_states: a swapped block orientation swaps cn_a and cn_b", {
  x = mc_toy()
  b2 = copy(x$blocks)[, sign := -sign]
  r1 = cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones, ref_cells = x$ref)$states[piece == "p1b"]
  r2 = cell_piece_states(x$features, x$pieces, b2, clones = x$clones, ref_cells = x$ref)$states[piece == "p1b"]
  m = merge(r1[, .(cell, a1 = cn_a, b1 = cn_b, f1 = f_obs)], r2[, .(cell, a2 = cn_a, b2 = cn_b, f2 = f_obs)], by = "cell")
  m[, lab := x$lab[cell]]
  # 1:0 <-> 0:1 are both candidates; 2:0 / 2:1 have no swapped candidate and fall back to the best remaining state
  expect_true(all(m[lab %in% c("s10", "s01"), a1 == b2 & b1 == a2]))
  expect_equal(m[lab == "s10", unique(paste(a2, b2))], "0 1")
  expect_true(all(abs(m$f1 + m$f2 - 1) < 1e-9 | is.na(m$f1)))
  expect_true(all(m[lab %in% c("ref", "bal"), a2 == 1 & b2 == 1]))
})

test_that("cell_piece_states: few votes and shallow depth fall back to the prior; clone_freq prior and t null", {
  x = mc_toy(depth_mb = 3, mol_per_atom = 0.005, seed = 11)
  r = cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones, ref_cells = x$ref, keep_ll = TRUE)
  p = r$states[piece == "p1b"]
  p[, lab := x$lab[cell]]
  # 1:0 and 0:1 have the same depth; with ~1 vote per cell the haplotype stays uncertain
  expect_lt(median(p[lab %in% c("s10", "s01"), posterior]), 0.95)
  expect_true(all(abs(r$ll[, sum(post), by = .(cell, piece)]$V1 - 1) < 1e-9))
  tp = data.table(class = c("5Mb", "10Mb", "arm", "chrom"), stratum = "20-50k", t_m = 0, t_s = 1, t_df = 4)
  r2 = cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones, ref_cells = x$ref, prior = "clone_freq", tpar = tp)
  expect_equal(r2$params$depth_null, "t")
  expect_equal(r2$params$prior, "clone_freq")
  # with no evidence the clone_freq MAP is the most frequent state (1:1 on p1b)
  expect_gt(mean(r2$states[piece == "p1b", state] == "1:1"), 0.5)
  # c / a overrides
  r3 = cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones, c = 0.6, a = 0.7, rho = NULL)
  expect_equal(c(r3$params$c, r3$params$a), c(0.6, 0.7))
  expect_equal(r3$params$n_ref, 0L)
  expect_equal(r3$params$rho_mode, "null_fit")
  expect_equal(r3$params$rho, 0.005)
  expect_error(cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones), "ref_cells")
})

test_that("cell_piece_states: rho from the reference cells (default) and by marginal likelihood", {
  x = mc_toy(n_each = 4L)
  # overdispersed votes: per cell x phase block, all molecules of a block from one haplotype with probability 0.5
  v = copy(x$features$votes)
  set.seed(5)
  v[, clump := runif(1) < 0.5, by = .(cell, ps)]
  v[clump == TRUE, c("nA", "nB") := if (runif(1) < 0.5) .(nA + nB, 0L) else .(0L, nA + nB), by = .(cell, ps)]
  f2 = x$features; f2$votes = v[, .(cell, atom, ps, nA = as.integer(nA), nB = as.integer(nB))]
  r0 = cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones, ref_cells = x$ref)
  r1 = cell_piece_states(f2, x$pieces, x$blocks, clones = x$clones, ref_cells = x$ref)
  expect_equal(r1$params$rho_mode, "ref")
  expect_gt(r1$params$rho_rows, 10)
  expect_gt(r1$params$rho, 5 * r0$params$rho)
  rc = x$ref
  # the simulated cells are still recovered with the larger rho (depth + strongly imbalanced votes)
  s1 = r1$states[piece == "p1b"][, lab := x$lab[cell]]
  expect_gt(mean(s1$state == c(ref = "1:1", bal = "1:1", s10 = "1:0", s20 = "2:0", s21 = "2:1", s01 = "0:1")[s1$lab]), 0.9)
  ra = cell_piece_states(f2, x$pieces, x$blocks, clones = x$clones, ref_cells = rc, rho = "auto", rho_grid = c(0.001, 0.01, 0.05, 0.1, 0.2))
  expect_equal(ra$params$rho_mode, "auto")
  expect_equal(nrow(ra$params$rho_ml), 5L)
  expect_gte(ra$params$rho, 0.01)
})

test_that("cell_piece_data: leave-one-out removes the cell's own votes from the block orientation", {
  x = mc_toy(n_each = 2L)
  d0 = cell_piece_data(x$features, x$pieces, x$blocks, ref_cells = x$ref)
  cc = x$cells[x$lab == "s10"][1]
  ci = match(cc, x$cells)
  # a reference orientation made of this cell's votes only: with loo every block becomes a tie and is dropped
  own = x$features$votes[cell == ci][x$features$atoms[, .(atom, chrom, astart = start)], on = "atom", nomatch = NULL]
  own = own[chrom == "chr1" & astart >= 30e6, .(nA_ref = sum(nA), nB_ref = sum(nB)), by = ps]
  b = copy(x$blocks)[piece == "p1b"][own, on = "ps", nomatch = NULL]
  b[, sign := as.integer(sign(nA_ref - nB_ref))]
  d1 = cell_piece_data(x$features, x$pieces, b, cells = cc)
  expect_gt(d1[piece == "p1b", oN], 0)
  expect_gt(d1[piece == "p1b", oA / oN], 0.95)
  d2 = cell_piece_data(x$features, x$pieces, b, cells = cc, loo = data.table(cell = cc, clone = "K10"))
  expect_equal(d2[piece == "p1b", oN], 0)
  # membership of another clone does not change anything
  d3 = cell_piece_data(x$features, x$pieces, b, cells = cc, loo = data.table(cell = cc, clone = "K20"))
  expect_equal(d3[piece == "p1b", oN], d1[piece == "p1b", oN])
  # depth: centred on the reference cells, n_win windows per piece, sd scaled for pieces shorter than their class
  expect_equal(d0[piece == "p1b", unique(n_win)], 30L)
  rc = x$ref
  expect_lt(abs(d0[cell %in% rc & piece == "p1a", median(log2ratio)]), 0.05)
  expect_true(all(d0$sd >= d0$sd_null))
  # s_sys is added in quadrature
  d4 = cell_piece_data(x$features, x$pieces, x$blocks, ref_cells = x$ref, s_sys = 0.2)
  expect_equal(d4$sd, sqrt(d0$sd^2 + 0.04))
})

test_that("cell_piece_states: a cell's global depth offset is removed by its own autosomal median", {
  x = mc_toy(n_each = 3L)
  f2 = x$features
  f2$win$x = f2$win$x * 1.4
  # every scored cell 40 % deeper than its expectation on every window, the reference pool unchanged
  rc = x$ref
  f2$win$x[, rc] = x$features$win$x[, rc]
  r1 = cell_piece_states(x$features, x$pieces, x$blocks, clones = x$clones, ref_cells = rc)$states
  r2 = cell_piece_states(f2, x$pieces, x$blocks, clones = x$clones, ref_cells = rc)$states
  m = merge(r1[, .(cell, piece, s1 = state, o1 = offset)], r2[, .(cell, piece, s2 = state, o2 = offset)], by = c("cell", "piece"))
  expect_true(all(m$s1 == m$s2))
  expect_true(all(abs(m[!cell %in% rc, o2 - o1] - log2(1.4)) < 0.02))
  r3 = cell_piece_states(f2, x$pieces, x$blocks, clones = x$clones, ref_cells = rc, center = FALSE)$states
  r3[, lab := x$lab[cell]]
  expect_true(any(r3[lab == "bal", state] != "1:1"))
})

test_that("cell_piece_data / cell_piece_states with loo on a piece flipped by chaining: the reference clone's cells keep its convention", {
  x = mc_toy(n_each = 6L)
  rc = x$cells[x$lab == "s20"]
  ci = match(rc, x$cells)
  # K20's pooled votes orient the blocks of p1b (unflipped: cn_a = K20's major haplotype)
  own = x$features$votes[cell %in% ci][x$features$atoms[, .(atom, chrom, astart = start)], on = "atom", nomatch = NULL]
  own = own[chrom == "chr1" & astart >= 30e6, .(nA_ref = sum(nA), nB_ref = sum(nB)), by = ps]
  b0 = copy(x$blocks)[piece == "p1b"][own, on = "ps", nomatch = NULL]
  b0[, `:=`(ref = "K20", sign = as.integer(sign(nA_ref - nB_ref)))]
  # the same table after a chaining flip of p1b: signs negated, reference counts unchanged, cn_a / cn_b swapped
  b1 = copy(b0)[, sign := -sign]
  b1_in = copy(b1)
  p1 = copy(x$pieces)[piece == "p1b", `:=`(cn_a = cn_b, cn_b = cn_a)]
  loo = data.table(cell = rc, clone = "K20")
  d0 = cell_piece_data(x$features, x$pieces, b0, cells = rc, loo = loo)[piece == "p1b"]
  d1 = cell_piece_data(x$features, p1, b1, cells = rc, loo = loo)[piece == "p1b"]
  expect_identical(b1, b1_in)                                         # the caller's block table is not modified
  expect_gt(min(d0$oN), 0)
  expect_equal(d1$oN, d0$oN)
  expect_equal(d1$oA, d0$oN - d0$oA)
  # the flip flag of orient_clone_states() gives the same; a table without the flip is unchanged
  d2 = cell_piece_data(x$features, p1, copy(b1)[, flipped := TRUE], cells = rc, loo = loo)[piece == "p1b"]
  expect_equal(d2$oA, d1$oA)
  d3 = cell_piece_data(x$features, x$pieces, copy(b0)[, flipped := FALSE], cells = rc, loo = loo)[piece == "p1b"]
  expect_equal(d3$oA, d0$oA)
  # states: the truly 2:0 cells of the reference clone are called as their clone, 2:0 unflipped and 0:2 after the flip
  st = function(p, b) cell_piece_states(x$features, p, b, cells = rc, clones = x$clones, ref_cells = x$ref, loo = loo,
                                         rho = 0.01)$states[piece == "p1b"]
  expect_true(all(st(x$pieces, b0)$state == "2:0"))
  expect_true(all(st(p1, b1)$state == "0:2"))
})
