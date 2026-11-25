# test-medicc2.R: MEDICC2 export (R/medicc2.R) on a hand-made tumour run: clones A, B, C (+ N0) on chr1 (60 Mb),
# chr2 (40 Mb) and a male chrX (40 Mb), a doublet clone D and an empty clone E. Truth: A and C lose haplotype B of
# chr1 30-60 Mb, B loses haplotype A there; B gains haplotype A on chr2 20-40 Mb; A is 2:2 on chr2; B gains chrX.
# C's chr1 breakpoint is 0.5 Mb off (a short piece to absorb), C's chr2 starts at 1 Mb (a gap) and C has an extra
# breakpoint at chr2 10 Mb (identical pieces to merge).

m2_toy = function(c_cnloh = FALSE) {
  seg = function(clone, chrom, start, end, state, ploidy = 2L)
    data.table(clone = clone, chrom = chrom, start = start, end = end, state = state, ploidy = ploidy)
  sg = rbind(
    seg("A", "chr1", 0, 30e6, "1:1"), seg("A", "chr1", 30e6, 60e6, "1:0"), seg("A", "chr2", 0, 40e6, "2:2"),
    seg("A", "chrX", 0, 40e6, "1:0", 1L),
    seg("B", "chr1", 0, 30e6, "1:1"), seg("B", "chr1", 30e6, 60e6, "1:0"), seg("B", "chr2", 0, 20e6, "1:1"),
    seg("B", "chr2", 20e6, 40e6, "2:1"), seg("B", "chrX", 0, 40e6, "2:0", 1L),
    seg("C", "chr1", 0, 30.5e6, "1:1"),
    if (c_cnloh) rbind(seg("C", "chr1", 30.5e6, 45.5e6, "1:0"), seg("C", "chr1", 45.5e6, 60e6, "2:0")) else seg("C", "chr1", 30.5e6, 60e6, "1:0"), seg("C", "chr2", 1e6, 10e6, "1:1"),
    seg("C", "chr2", 10e6, 40e6, "1:1"), seg("C", "chrX", 0, 40e6, "1:0", 1L),
    seg("E", "chr1", 0, 60e6, "1:0"), seg("E", "chr2", 0, 40e6, "1:1"), seg("E", "chrX", 0, 40e6, "1:0", 1L),
    seg("N0", "chr1", 0, 60e6, "1:1"), seg("N0", "chr2", 0, 40e6, "1:1"), seg("N0", "chrX", 0, 40e6, "1:0", 1L))
  sg[, `:=`(cn = nA_of(state) + nB_of(state), label = "x", nA = NA_integer_, nB = NA_integer_, f = NA_real_)]
  clones = data.table(clone = c("A", "B", "C", "N0", "D", "E"), n_assigned = c(40, 40, 40, 40, NA, 0),
                      a = c(1, 1, 1, 1, NA, 1), doublet = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  bins = make_bins(c(chr1 = 60e6, chr2 = 40e6, chrX = 40e6), 1e5)
  atoms = data.table(atom = seq_len(nrow(bins)), chrom = bins$chrom, start = bins$start, end = bins$end, bin = seq_len(nrow(bins)))
  lab = rep(c("A", "B", "C", "N0"), each = 40)
  ev = data.table(clone = c("A", "C", "B", "B", "A"), chrom = c("chr1", "chr1", "chr1", "chr2", "chr2"),
                  start = c(30e6, 30.5e6, 30e6, 20e6, 0), end = c(60e6, 60e6, 60e6, 40e6, 40e6),
                  cnA = c(1L, 1L, 0L, 2L, 2L), cnB = c(0L, 0L, 1L, 1L, 2L))
  if (c_cnloh) ev = rbind(ev, data.table(clone = "C", chrom = "chr1", start = 45.5e6, end = 60e6, cnA = 2L, cnB = 0L))
  votes = fake_votes(length(lab), atoms, lab, ev, flip_prob = 0.5, seed = 3, block_atoms = 10L, mol_per_atom = 0.5)
  cells = sprintf("cell%03d", seq_along(lab))
  asg = data.table(cell = cells, clone = lab, confident = TRUE, assigned = TRUE)
  list(segments = sg, clones = clones, features = list(votes = votes, atoms = atoms, cells = data.table(cell = cells)),
       assignments = asg)
}

test_that("medicc2_pieces: clone filter, coverage, short-piece absorption", {
  x = m2_toy()
  pc = medicc2_pieces(x$segments, x$clones, min_len = 1e6)
  expect_setequal(unique(pc$clone), c("A", "B", "C", "N0"))
  dr = attr(pc, "clones_dropped")
  expect_equal(dr[order(clone), reason], c("doublet", "no_assigned_cells"))
  pd = attr(pc, "pieces_dropped")
  expect_true(pd[chrom == "chr2" & start == 0 & end == 1e6, reason] == "not_all_clones")
  u = unique(pc[, .(chrom, start, end)])
  expect_true(all(u$end - u$start >= 1e6))
  # the 30-30.5 Mb piece (A, B 1:0, C 1:1) joins the 30.5-60 Mb piece (all 1:0)
  expect_equal(u[chrom == "chr1", paste(start, end)], c("0 3e+07", "3e+07 6e+07"))
  expect_equal(pc[chrom == "chr1" & start == 30e6, sort(state)], c("1:0", "1:0", "1:0", "1:1"))
  # every clone on every piece; chr2 keeps the 10 Mb breakpoint until the table merges it
  expect_true(all(pc[, .N, by = piece]$N == 4))
  expect_equal(u[chrom == "chr2", start], c(1e6, 10e6, 20e6))
  expect_equal(pc[chrom == "chr2" & start == 20e6 & clone == "B", imbalanced], TRUE)
  expect_equal(pc[chrom == "chr2" & clone == "A", unique(imbalanced)], FALSE)
  expect_equal(pc[chrom == "chrX", unique(imbalanced)], FALSE)
})

test_that("medicc2_pieces: an imbalance segment without depth support enters as 1:1", {
  x = m2_toy()
  x$segments[clone == "B" & chrom == "chr2" & start == 20e6, label := "imbalance"]
  pc = medicc2_pieces(x$segments, x$clones, min_len = 1e6)
  b = pc[chrom == "chr2" & start == 20e6 & clone == "B"]
  expect_equal(b$state, "1:1"); expect_equal(b$cn, 2); expect_false(b$imbalanced)
  # every other row as without the label
  p0 = medicc2_pieces(m2_toy()$segments, x$clones, min_len = 1e6)
  expect_equal(pc[!(chrom == "chr2" & clone == "B"), .(piece, clone, state)], p0[!(chrom == "chr2" & clone == "B"), .(piece, clone, state)])
})

test_that("orient_clone_states: opposite haplotype losses are swapped, balanced and haploid pass through", {
  x = m2_toy()
  pc = medicc2_pieces(x$segments, x$clones)
  o = orient_clone_states(pc, x$features, x$assignments, x$clones)
  expect_s3_class(o, "medicc2_orientation")
  p = o$pieces[chrom == "chr1" & start == 30e6]
  ab = function(k) p[clone == k, c(cn_a, cn_b)]
  expect_equal(ab("A"), ab("C"))
  expect_equal(ab("B"), rev(ab("A")))
  expect_equal(sort(ab("A")), c(0L, 1L))
  expect_equal(sum(p$is_ref), 1L)
  expect_true(all(p[clone != "N0", dev] > 0.3))
  expect_true(all(p[clone %in% c("A", "B", "C"), flag] %in% c("ref", "ok")))
  # N0 on an oriented piece: balanced, f near 0.5
  expect_equal(p[clone == "N0", c(cn_a, cn_b, flag == "balanced")], c(1L, 1L, 1L))
  expect_lt(p[clone == "N0", dev], 0.05)
  # single imbalanced clone: it is its own reference; balanced 2:2 passes
  q = o$pieces[chrom == "chr2" & start == 20e6]
  expect_equal(q[clone == "B", c(cn_a, cn_b)], c(2L, 1L))
  expect_true(q[clone == "B", is_ref])
  expect_equal(q[clone == "A", c(cn_a, cn_b)], c(2L, 2L))
  expect_equal(q[clone == "A", flag], "balanced")
  # haploid chrX
  r = o$pieces[chrom == "chrX"]
  expect_equal(r[clone == "B", c(cn_a, cn_b)], c(2L, 0L))
  expect_true(all(r$flag == "haploid"))
  expect_true(all(r$cn_b == 0L))
  # QC and block orientation
  expect_equal(o$qc[chrom == "chr1" & start == 30e6, n_imbalanced], 3L)
  expect_true(is.na(o$qc[chrom == "chr1" & start == 0, orient_ref]))
  b = o$blocks[chrom == "chr1" & start == 30e6]
  expect_true(nrow(b) > 0 && all(b$sign %in% c(-1L, 0L, 1L)))
  expect_equal(uniqueN(b$ref), 1L)
  # the reference orientation recovers the true block flips (up to one global sign)
  tb = attr(x$features$votes, "blocks")[chrom == "chr1"]
  m = b[tb, on = "ps", nomatch = NULL][sign != 0]
  agree = mean((m$sign > 0) == !m$flip)
  expect_true(agree > 0.95 || agree < 0.05)
})

test_that("orient_clone_states chains the haplotype label across adjacent pieces", {
  x = m2_toy(c_cnloh = TRUE)
  pc = medicc2_pieces(x$segments, x$clones)
  for (ch in c("phase", "parsimony")) {
    o = orient_clone_states(pc, x$features, x$assignments, x$clones, chain = ch, min_link = 10)
    p = o$pieces[chrom == "chr1" & start >= 30e6]
    expect_equal(uniqueN(p$start), 2L)
    # A, B keep their labels over the C breakpoint; C keeps the haplotype of A (2:0 after 1:0)
    for (k in c("A", "B")) expect_equal(uniqueN(p[clone == k, paste(cn_a, cn_b)]), 1L)
    expect_equal(p[clone == "C", sign(cn_a - cn_b)], p[clone == "A", sign(cn_a - cn_b)])
    expect_equal(o$qc[chrom == "chr1" & start == 45.5e6, link], ch)
    # f is the fraction of the cn_a haplotype after any flip
    q = p[flag %in% c("ok", "ref")]
    expect_true(all((q$f >= 0.5) == (q$cn_a > q$cn_b)))
  }
  # block signs follow the flips: the signs of a phase set spanning the breakpoint agree
  b = o$blocks[chrom == "chr1" & start >= 30e6 & sign != 0]
  sp = b[, .(n = uniqueN(piece), s = uniqueN(sign)), by = ps][n == 2]
  expect_true(nrow(sp) >= 1 && all(sp$s == 1))
  o0 = orient_clone_states(pc, x$features, x$assignments, x$clones, chain = "none")
  expect_true(all(o0$qc$link == "first") && !any(o0$qc$flipped))
  # without chaining the two pieces have different references (B, C) whose majors are opposite haplotypes
  expect_equal(o0$qc[chrom == "chr1" & start >= 30e6, orient_ref], c("B", "C"))
  expect_equal(o0$pieces[chrom == "chr1" & start >= 30e6 & clone == "A", uniqueN(paste(cn_a, cn_b))], 2L)
})

test_that("orient_clone_states flags clones with an imbalanced state but no oriented signal", {
  x = m2_toy()
  # declare C 2:1 on chr2 20-40 Mb although its votes are balanced
  sg = copy(x$segments)
  sg = rbind(sg[!(clone == "C" & chrom == "chr2" & start == 10e6)],
             data.table(clone = "C", chrom = "chr2", start = c(10e6, 20e6), end = c(20e6, 40e6), state = c("1:1", "2:1"),
                        ploidy = 2L, cn = c(2L, 3L), label = "x", nA = NA_integer_, nB = NA_integer_, f = NA_real_))
  o = orient_clone_states(medicc2_pieces(sg, x$clones), x$features, x$assignments, x$clones)
  q = o$pieces[chrom == "chr2" & start == 20e6]
  expect_equal(q[clone == "C", flag], "low_conf")
  expect_equal(o$qc[chrom == "chr2" & start == 20e6, n_low_conf], 1L)
})

test_that("medicc2_table: N0 sample, merging, total = a + b, cap and order", {
  x = m2_toy()
  o = orient_clone_states(medicc2_pieces(x$segments, x$clones), x$features, x$assignments, x$clones)
  as_ = medicc2_table(o)
  expect_equal(names(as_), c("sample_id", "chrom", "start", "end", "cn_a", "cn_b"))
  expect_equal(as_$sample_id[1], "N0")
  expect_setequal(unique(as_$sample_id), c("N0", "A", "B", "C"))
  expect_equal(as_[sample_id == "N0" & chrom == "chrX", c(cn_a, cn_b)], c(1L, 0L))
  # chr2 1-10 and 10-20 Mb are identical in every sample: one segment 1-20 Mb
  expect_equal(as_[sample_id == "A" & chrom == "chr2", paste(start, end)], c("1000000 20000000", "20000000 40000000"))
  # same segments in every sample
  seg_key = as_[, paste(chrom, start, end)]
  expect_true(all(table(seg_key) == 4))
  expect_equal(unique(as_$chrom), c("chr1", "chr2", "chrX"))
  expect_true(is.integer(as_$start) && is.integer(as_$cn_a))
  tot = medicc2_table(o, allele_specific = FALSE)
  expect_equal(names(tot), c("sample_id", "chrom", "start", "end", "cn"))
  # total merges more: chr1 30-60 Mb is 1 copy in A, B and C; chr2 is one piece in A (4) ... still 2 pieces (B 2 vs 3)
  m = merge(as_, tot, by = c("sample_id", "chrom", "start", "end"))
  expect_equal(m$cn, m$cn_a + m$cn_b)
  expect_lte(nrow(tot), nrow(as_))
  # without N0 in the input, the table adds a diploid / haploid-X normal
  d = o$pieces[clone != "N0"]
  t2 = medicc2_table(d)
  expect_equal(t2[sample_id == "N0" & chrom == "chr1", unique(paste(cn_a, cn_b))], "1 1")
  expect_equal(t2[sample_id == "N0" & chrom == "chrX", unique(paste(cn_a, cn_b))], "1 0")
  # cap
  d2 = copy(d)[clone == "A" & chrom == "chr2", cn_a := 12L]
  expect_warning(ta <- medicc2_table(d2), "above cap = 8 and were clipped")
  expect_equal(max(ta$cn_a), 8L)
  expect_warning(tt <- medicc2_table(d2, allele_specific = FALSE), "above cap = 8 and were clipped")
  expect_equal(max(tt$cn), 8L)
})

test_that("medicc2_pieces / orient_clone_states take numeric clone ids and do not modify the caller's tables", {
  x = m2_toy()
  ids = c(A = 3, B = 5, C = 7, N0 = 9, D = 11, E = 13)
  sg = copy(x$segments)[, clone := ids[clone]]
  cl = copy(x$clones)[, clone := ids[clone]]
  cl_in = copy(cl)
  pn = medicc2_pieces(sg, cl, min_len = 1e6, normal = 9)
  pc = medicc2_pieces(x$segments, x$clones, min_len = 1e6)
  expect_identical(cl, cl_in)                                        # no `reason` column added to the caller's table
  expect_type(pn$clone, "character")
  key = function(p) p[, .(piece, clone, state, cn, imbalanced)][order(piece, clone)]
  pc_n = copy(pc)[, clone := as.character(ids[clone])]
  expect_identical(key(pn), key(pc_n))
  expect_identical(attr(pn, "clones_dropped")[order(clone), reason], attr(pc, "clones_dropped")[order(clone), reason])
  asg = copy(x$assignments)[, clone := ids[clone]]
  on = orient_clone_states(pn, x$features, asg, cl)
  oc = orient_clone_states(pc, x$features, x$assignments, x$clones)
  expect_identical(on$pieces[order(piece, clone), .(cn_a, cn_b, flag)], oc$pieces[, clone := as.character(ids[clone])][order(piece, clone), .(cn_a, cn_b, flag)])
})

test_that("MEDICC2 export: no segment and no absorbed piece spans a dropped range", {
  # medicc2_table: chr1 [0, 10 Mb) and [30, 50 Mb) with equal states, [10, 30 Mb) dropped -> two segments
  d = data.table(clone = rep(c("K1", "K2"), each = 2), chrom = "chr1", start = c(0, 30e6), end = c(10e6, 50e6),
                 cn_a = rep(c(1L, 2L), each = 2), cn_b = 0L, ploidy = 2L)
  tb = medicc2_table(d)
  expect_equal(tb[sample_id == "K1", .(start, end)], data.table(start = c(0L, 30000000L), end = c(10000000L, 50000000L)))
  expect_true(all(tb[, .N, by = sample_id]$N == 2L))
  # contiguous pieces with equal states are still merged
  d2 = copy(d)[start == 30e6, start := 10e6]
  expect_equal(nrow(medicc2_table(d2)[sample_id == "K1"]), 1L)
  # medicc2_pieces: a short piece next to a dropped range is absorbed into its contiguous neighbour only
  seg = function(clone, start, end, state) data.table(clone = clone, chrom = "chr1", start = start, end = end, state = state)
  sg = rbind(seg("A", c(0, 30.5e6), c(30.5e6, 50e6), c("1:0", "1:1")),
             seg("B", c(0, 30e6), c(10e6, 50e6), c("1:0", "1:0")))
  pc = medicc2_pieces(sg, min_len = 1e6)
  u = unique(pc[, .(start, end)])[order(start)]
  expect_equal(u, data.table(start = c(0, 30e6), end = c(10e6, 50e6)))
  pd = attr(pc, "pieces_dropped")
  expect_equal(pd[reason == "not_all_clones", .(start, end)], data.table(start = 10e6, end = 30e6))
  # the absorbed short piece takes the contiguous neighbour's states
  expect_equal(pc[start == 30e6][order(clone), state], c("1:1", "1:0"))
  # max_gap opts in to absorption across the gap (the old default): the left piece then claims [10, 30 Mb)
  pg = medicc2_pieces(sg, min_len = 1e6, max_gap = Inf)
  expect_equal(unique(pg[, .(start, end)])[order(start)], data.table(start = c(0, 30.5e6), end = c(30.5e6, 50e6)))
})

test_that("orient_clone_states: the block table carries the piece's chaining flip", {
  x = m2_toy(c_cnloh = TRUE)
  pc = medicc2_pieces(x$segments, x$clones)
  for (ch in c("parsimony", "none")) {
    o = orient_clone_states(pc, x$features, x$assignments, x$clones, chain = ch, min_link = 10)
    b = o$blocks
    expect_true(is.logical(b$flipped))
    expect_identical(b$flipped, o$qc$flipped[match(b$piece, o$qc$piece)])
    expect_true(all(b$sign == fifelse(b$flipped, -1L, 1L) * sign(b$nA_ref - b$nB_ref)))
  }
})

test_that("medicc2_pieces: an imbalance segment enters at the balanced state of its clone's integer baseline", {
  x = m2_toy()
  x$segments[clone == "B" & chrom == "chr2" & start == 20e6, label := "imbalance"]
  at = function(B) {
    cl = copy(x$clones)[, baseline := 2][clone == "B", baseline := B]
    medicc2_pieces(x$segments, cl, min_len = 1e6)[chrom == "chr2" & start == 20e6 & clone == "B"]
  }
  b2 = at(2)
  expect_equal(b2$state, "1:1"); expect_equal(b2$cn, 2); expect_equal(b2$label, "neutral")
  expect_identical(b2, medicc2_pieces(x$segments, x$clones, min_len = 1e6)[chrom == "chr2" & start == 20e6 & clone == "B"])
  b4 = at(4)
  expect_equal(b4$state, "2:2"); expect_equal(b4$cn, 4); expect_false(b4$imbalanced)
  expect_equal(at(3.75)$state, "2:2")                                  # nb = 4
  b3 = at(3)                                                           # odd nb: no balanced state at the baseline, kept
  expect_equal(b3$state, "2:1"); expect_equal(b3$cn, 3); expect_true(b3$imbalanced)
})

test_that("medicc2_table counts and warns about copy numbers above the cap", {
  d = data.table(clone = c("K1", "K2"), chrom = "chr1", start = 0, end = 10e6, cn_a = c(9L, 2L), cn_b = c(1L, 1L), ploidy = 2L)
  expect_warning(tb <- medicc2_table(d), "1 sample x piece rows have a copy number above cap = 8")
  expect_equal(attr(tb, "n_capped"), 1L)
  expect_equal(tb[sample_id == "K1", c(cn_a, cn_b)], c(8L, 1L))
  expect_no_warning(t2 <- medicc2_table(d, cap = 12L))
  expect_equal(attr(t2, "n_capped"), 0L)
})
