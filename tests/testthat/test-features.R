# test-features.R: the feature step (R/features.R) on hand-made cases and the two cached simulations of
# helper-caller.R: molecule collapse parity with a frozen reference, region naming and the depth / allele membership
# rules, sd_class, atoms, depth features vs a manual loop, oriented allele fractions of event carriers,
# the leave-one-out leak test, male chrX and clone-grouped orientation.

# frozen reference implementation of the molecule collapse (the parity oracle of collapse_molecules())
ref_e2_molecules = function(fr, gap = 2000) {
  d = fr[hap %in% c("A", "B"), .(cell, chrom, start, end, hap, ps)]
  setorder(d, cell, chrom, start)
  d[, new := c(TRUE, start[-1] > cummax(end)[-.N] + gap), by = .(cell, chrom)]
  d[, mol := cumsum(new)]
  m = d[, .(start = min(start), end = max(end), ps = ps[1], nA = sum(hap == "A"), nB = sum(hap == "B"), n_frag = .N), by = .(cell, chrom, mol)]
  m = m[nA != nB]
  m[, hap := ifelse(nA > nB, "A", "B")]
  m[, c("nA", "nB", "mol") := NULL]
  m[]
}

# centromere table of sim_genome(): p = [0, 0.4 x length)
sim_centromeres = function(bins, chr1_left = NULL, chr1_right = NULL) {
  cs = bins[, .(len = max(end)), by = chrom]
  ce = cs[, .(chrom, cstart = 1, cen_left = 0.4 * len, cen_right = 0.4 * len, cend = len)]
  if (!is.null(chr1_left)) ce[chrom == "chr1", `:=`(cen_left = chr1_left, cen_right = chr1_right)]
  ce[]
}

feat_cache = new.env()
sim_features = function(name, ...) {
  if (is.null(feat_cache[[name]])) {
    s = get(name)
    d = file.path(tempdir(), paste0("feat_", name))
    # bins (counts only), fragments and cells_qc: all cell_features() reads
    dir.create(d, showWarnings = FALSE)
    sc = s$scn; sc$layers = list(); sc$obsm = list()
    invisible(capture.output(suppressMessages(write_bins_h5ad(sc, file.path(d, "bins_100kb.h5ad")))))
    fr = copy(s$fragments)[, cell := s$scn$cells$cell[cell]]
    write_fragments(fr, file.path(d, "fragments.tsv.gz"))
    fwrite(s$scn$cells, file.path(d, "cells_qc.tsv"), sep = "\t")
    feat_cache[[name]] = cell_features(d, cells = s$scn$cells$cell, centromeres = sim_centromeres(s$scn$bins), ...)
    feat_cache[[paste0(name, "_dir")]] = d
  }
  feat_cache[[name]]
}

test_that("collapse_molecules equals the reference molecule collapse", {
  fr = data.table(
    cell = c("c1", "c1", "c1", "c1", "c1", "c2", "c1", "c2", "c2", "c2", "c1", "c1"),
    chrom = "chr1",
    start = c(100, 600, 2600, 5000, 7001, 150, 9502, 400, 2401, 9000, 20000, 20100),
    end = c(500, 600, 3000, 5000, 7500, 400, 9600, 400, 3000, 9100, 20050, 20200),
    hap = c("A", "A", "B", "A", "B", "B", "A", "A", "B", "B", "A", "B"),
    ps = c(11, 12, 13, 14, 15, 21, 16, 22, 23, 24, 17, 18))
  # c1: 600 -> 2600 and 3000 -> 5000 are exactly gap (one molecule, A A B A -> A, ps 11); 5000 -> 7001 is
  # gap + 1 (new: B, ps 15); 7500 -> 9502 is gap + 2 (new: A, ps 16); 20000 / 20100: tie (dropped).
  # c2 interleaved with c1; 400 -> 2401 is gap + 1 (new molecule)
  fr = rbind(fr, data.table(cell = "c1", chrom = "chr2", start = 50, end = 60, hap = "N", ps = 1))
  ref = ref_e2_molecules(fr)
  got = collapse_molecules(fr)
  expect_equal(got, ref)
  expect_equal(got[cell == "c1" & chrom == "chr1", hap], c("A", "B", "A"))
  expect_equal(got[cell == "c1" & chrom == "chr1", ps], c(11, 15, 16))
  expect_equal(got[cell == "c1" & chrom == "chr1", n_frag], c(4L, 1L, 1L))
  expect_equal(nrow(got[cell == "c1" & start == 20000]), 0L)
  # integer haplotype codes and a shuffled input give the same molecules
  fr2 = copy(fr)[sample(.N)][, hap := hap_code(hap)]
  expect_equal(collapse_molecules(fr2)[order(cell, chrom, start)], ref[order(cell, chrom, start)])
  # random fragments: exact parity
  set.seed(3)
  n = 4000
  rf = data.table(cell = sample(sprintf("b%02d", 1:15), n, TRUE), chrom = sample(c("chr1", "chr2"), n, TRUE),
                  start = round(runif(n, 0, 2e6)))
  rf[, end := start + round(rexp(n, 1 / 800))][, hap := sample(c("A", "B", "U", "N"), n, TRUE, c(0.4, 0.4, 0.15, 0.05))]
  rf[, ps := sample(1:5, n, TRUE)]
  expect_equal(collapse_molecules(rf, gap = 1500), ref_e2_molecules(rf, gap = 1500))
})

test_that("collapse_molecules split_ps: one vote per molecule x phase set", {
  # two fragments 500 bp apart in different phase sets: one molecule under the first PS today, two votes split
  fr = data.table(cell = "c1", chrom = "chr1", start = c(1000, 2000), end = c(1500, 2600), hap = "A", ps = c(50, 4000))
  m0 = collapse_molecules(fr)
  expect_equal(m0[, .(start, end, ps, n_frag, hap)], data.table(start = 1000, end = 2600, ps = 50, n_frag = 2L, hap = "A"))
  m1 = collapse_molecules(fr, split_ps = TRUE)
  expect_equal(names(m1), names(m0))
  expect_equal(m1[, .(start, end, ps, n_frag, hap)],
               data.table(start = c(1000, 2000), end = c(1500, 2600), ps = c(50, 4000), n_frag = 1L, hap = "A"))
  # PS 50 votes A, B, A and PS 4000 votes B: a 2:2 tie (dropped) today, A in PS 50 and B in PS 4000 split
  fr2 = data.table(cell = "c1", chrom = "chr1", start = c(1000, 1600, 2300, 2900), end = c(1400, 2000, 2700, 3200),
                   hap = c("A", "B", "B", "A"), ps = c(50, 50, 4000, 50))
  expect_equal(nrow(collapse_molecules(fr2)), 0L)
  m2 = collapse_molecules(fr2, split_ps = TRUE)
  expect_equal(m2[, .(start, end, ps, n_frag, hap)],
               data.table(start = c(1000, 2300), end = c(3200, 2700), ps = c(50, 4000), n_frag = c(3L, 1L), hap = c("A", "B")))
  # one phase set per cell and chromosome: identical under both settings
  set.seed(4)
  n = 3000
  rf = data.table(cell = sample(sprintf("b%02d", 1:10), n, TRUE), chrom = sample(c("chr1", "chr2"), n, TRUE),
                  start = round(runif(n, 0, 2e6)))
  rf[, end := start + round(rexp(n, 1 / 800))][, hap := sample(c("A", "B", "U"), n, TRUE, c(0.45, 0.45, 0.1))]
  rf[, ps := fifelse(chrom == "chr1", 11, 22)]
  expect_equal(collapse_molecules(rf, split_ps = TRUE), collapse_molecules(rf))
  # mixed phase sets: the split never merges molecules, it only splits them by phase set
  rf[, ps := sample(1:3, n, TRUE)]
  ms = collapse_molecules(rf, split_ps = TRUE)
  expect_gt(nrow(ms), nrow(collapse_molecules(rf)))
  expect_equal(anyDuplicated(ms[, .(cell, chrom, start, ps)]), 0L)
  expect_equal(collapse_molecules(rf[0], split_ps = TRUE), collapse_molecules(rf[0]))
})

test_that("read_votes clean: blacklisted fragments do not vote (tabix and fallback)", {
  # flags: 2 = MAPQ 60, 8 = blacklist start, 16 = TRF start (only the blacklist bit makes a fragment not clean)
  fr = data.table(chrom = "chr1", start = c(1000, 100000, 200000, 300000, 1500000, 1000),
                  end = c(1500, 100500, 200500, 300500, 1500500, 1500),
                  cell = c("c1", "c1", "c1", "c1", "c1", "c2"), strand = "+", mapq = 60L,
                  flags = c(2L, 8L, 10L, 2L, 18L, 8L))
  fr[, `:=`(len = end - start, piece = 1L, npieces = 1L, hap = c(1L, 2L, 1L, 1L, 2L, 2L), nsnp = 1L, ps = 50,
            read_id = seq_len(.N))]
  f = file.path(tempfile("votes_clean"), "fragments.tsv.gz")
  dir.create(dirname(f))
  suppressWarnings(write_fragments(fr, f))
  atoms = data.table(atom = 1:2, chrom = "chr1", start = c(0, 1e6), end = c(1e6, 2e6), bin = 1:2)
  v_all = data.table(cell = c(1L, 1L, 2L), atom = c(1L, 2L, 1L), ps = 50, nA = c(3L, 0L, 0L), nB = c(1L, 1L, 1L))
  v_clean = data.table(cell = c(1L, 1L), atom = 1:2, ps = 50, nA = c(2L, 0L), nB = c(0L, 1L))
  paths = list(fallback = NULL)
  if (nzchar(Sys.which("tabix")) && file.exists(paste0(f, ".tbi"))) paths$tabix = ""
  for (p in names(paths)) {
    tb = paths[[p]]
    v0 = read_votes(f, c("c1", "c2"), atoms, tabix = tb)
    expect_equal(v0, v_all, ignore_attr = TRUE, info = p)
    expect_equal(read_votes(f, c("c1", "c2"), atoms, tabix = tb, clean = FALSE, split_ps = FALSE), v0, info = p)
    v1 = read_votes(f, c("c1", "c2"), atoms, tabix = tb, clean = TRUE)
    expect_equal(v1, v_clean, ignore_attr = TRUE, info = p)
    expect_silent(check_votes(v1))
  }
})

test_that("cell_features passes mol_split_ps and votes_clean to read_votes", {
  ev = data.table(region = "loh2", chrom = "chr2", start = 20e6, end = 40e6)
  sim_features("sim_norm_f", method = "celltype", sex = "female", events = ev, classes = c("arm", "chrom", "10Mb"))
  got = list()
  local_mocked_bindings(read_votes = function(...) { got <<- list(...); stop("read_votes reached") })
  d = feat_cache$sim_norm_f_dir
  cf = function(params) cell_features(d, cells = sim_norm_f$scn$cells$cell, centromeres = sim_centromeres(sim_norm_f$scn$bins),
                                      sex = "female", classes = "arm", params = params)
  expect_error(cf(caller_defaults()), "read_votes reached")
  expect_false(got$split_ps); expect_false(got$clean)
  p = modifyList(caller_defaults(), list(mol_split_ps = TRUE, votes_clean = TRUE))
  expect_error(cf(p), "read_votes reached")
  expect_true(got$split_ps); expect_true(got$clean)
})

test_that("regions: names, arms, containment and start-inside rules, sd_class", {
  bins = sim_genome()
  ce = sim_centromeres(bins, chr1_left = 24.05e6, chr1_right = 24.05e6)
  r = region_set(bins, ce, classes = c("arm", "chrom", "5Mb", "50Mb"), sex = "female", null_exclude = "chr8")
  expect_silent(check_regions(r))
  expect_true(all(c("chr1p", "chr1q", "chrXp", "chrXq", "chr21q", "chr22q", "chr1", "5Mb|chr1:0", "5Mb|chr1:55000000",
                    "50Mb|chr1:50000000") %in% r$region))
  expect_false(any(c("chr21p", "chr22p") %in% r$region))                  # acrocentric p arms
  expect_equal(r[region == "chr1p", c(start, end)], c(0, 24.05e6))
  expect_equal(r[region == "50Mb|chr1:50000000", end], 60e6)
  expect_equal(r[region == "chr1p", n_bins], 240L)                          # 24.0-24.1 Mb straddles the centromere
  expect_equal(r[region == "chr1q", n_bins], 359L)
  expect_false(any(r[chrom == "chr8", null_region]))
  expect_true(all(r[chrom == "chrX", !null_region & allelic_ok]))
  expect_false(r[region == "50Mb|chr1:50000000", eligible])                # 10 Mb of a 50 Mb tile
  expect_true(r[region == "5Mb|chr1:0", eligible])
  rm_ = region_set(bins, ce, classes = c("arm", "chrom"), sex = "male")
  expect_true(all(!rm_[chrom == "chrX", allelic_ok]))
  M = region_bin_matrix(bins, r[region %in% c("chr1p", "chr1q")])
  expect_equal(as.vector(Matrix::rowSums(M)), c(240, 359))
  expect_equal(sum(M[, bins$chrom == "chr1" & bins$start == 24e6]), 0)
  # atoms split at the non-aligned centromere; a molecule starting at 24.05 Mb belongs to q, one before to p
  at = atoms_from_regions(bins, r)
  expect_silent(check_atoms(at))
  expect_true(all(c(24e6, 24.05e6) %in% at[chrom == "chr1", start]))
  a1 = at[chrom == "chr1" & start == 24e6, atom]; a2 = at[chrom == "chr1" & start == 24.05e6, atom]
  v = data.table(cell = 1L, atom = c(a1, a2), ps = 1, nA = c(3L, 5L), nB = c(0L, 1L))
  af = allelic_features(v, at, r[class == "arm"], orientation = "chromosome")
  expect_equal(af[region == "chr1p", c(nA, nB)], c(3, 0))
  expect_equal(af[region == "chr1q", c(nA, nB)], c(5, 1))
  expect_equal(af$orient[1], "anchored")
  # sd_class of events
  ev = data.table(region = c("e_chrom", "e_arm", "e_tile"), class = "event", chrom = c("chr2", "chr1", "chr1"),
                  start = c(0, 30e6, 30e6), end = c(38e6, 55.2e6, 42e6))
  arms = arm_table(ce, chrom_sizes_from_bins(bins), "female")
  expect_equal(sd_class_of(ev, bins, arms), c("chrom", "arm", "10Mb"))
  re = region_set(bins, ce, classes = "arm", events = ev[, .(region, chrom, start, end)])
  expect_equal(re[class == "event", sd_class], c("arm", "10Mb", "chrom")[match(re[class == "event", region], c("e_arm", "e_tile", "e_chrom"))])
  expect_false(any(re[class == "event", null_region]))
})

test_that("depth_features equals a manual computation and z is null-centred with MAD 1", {
  s = sim_norm_f
  scn = s$scn
  use = bin_filter(scn$bins, scn$counts)
  ce = sim_centromeres(scn$bins)
  r = region_set(scn$bins, ce, classes = c("arm", "chrom", "10Mb"), events = data.table(region = "ev", chrom = "chr2", start = 20e6, end = 40e6),
                 null_exclude = "chr21", use = use)
  norm = normalise_depth(scn, "celltype", use = use, cell_type = scn$cells$cell_type)
  d = depth_features(scn, norm, r)
  expect_equal(nrow(d), nrow(r) * ncol(scn$counts))
  E = norm$expected; E[is.na(E)] = 0
  X = as.matrix(scn$counts)
  B = lapply(seq_len(nrow(r)), function(k) which(scn$bins$chrom == r$chrom[k] & scn$bins$start >= r$start[k] & scn$bins$end <= r$end[k] & use))
  xm = t(vapply(B, function(b) colSums(X[b, , drop = FALSE]), numeric(ncol(X))))
  em = t(vapply(B, function(b) colSums(E[b, , drop = FALSE]), numeric(ncol(X))))
  lr = log2((xm + 0.5) / (em + 0.5))
  sdm = matrix(NA_real_, nrow(r), ncol(X))
  for (j in seq_len(ncol(X))) for (k in seq_len(nrow(r))) {
    ref = r$class == r$sd_class[k] & r$null_region
    sdm[k, j] = 1.4826 * mad(lr[ref, j], constant = 1)
  }
  man = data.table(cell = rep(scn$cells$cell, each = nrow(r)), region = rep(r$region, ncol(X)), x = as.vector(xm), e = as.vector(em),
                   log2ratio = as.vector(lr), sd_null = as.vector(sdm), z = as.vector(lr / pmax(sdm, 0.05)))
  setkeyv(man, c("cell", "region"))
  expect_equal(d$x, man$x)
  expect_equal(d$e, man$e, tolerance = 1e-10)
  expect_equal(d$sd_null, man$sd_null, tolerance = 1e-10)
  expect_equal(d$z, man$z, tolerance = 1e-10)
  expect_equal(d[region == "ev", unique(sd_class)], "arm")
  # null-centred, MAD 1 per cell on the null arms (sd above the floor)
  nul = d[region %in% r[null_region == TRUE & class == "arm", region]]
  st = nul[, .(med = median(z), mad1 = 1.4826 * mad(z, constant = 1), sd = sd_null[1]), by = cell]
  expect_lt(abs(median(st$med)), 0.25)
  expect_equal(st[sd > 0.05, mad1], rep(1, sum(st$sd > 0.05)), tolerance = 1e-8)
  w = window_matrices(scn, norm)
  expect_equal(colSums(w$x), colSums(X[use, ]), ignore_attr = TRUE)
  expect_equal(nrow(w$windows), nrow(w$x))
})

test_that("cell_features on sim_norm_f: anchored orientation, carriers vs non-carriers, validators", {
  s = sim_norm_f
  ev = data.table(region = "loh2", chrom = "chr2", start = 20e6, end = 40e6)
  fe = sim_features("sim_norm_f", method = "celltype", sex = "female", events = ev, classes = c("arm", "chrom", "10Mb"))
  expect_s3_class(fe, "scn_features")
  expect_silent(check_features(fe$long)); expect_silent(check_votes(fe$votes)); expect_silent(check_atoms(fe$atoms))
  expect_silent(check_regions(fe$regions))
  expect_equal(nrow(fe$long), nrow(fe$cells) * nrow(fe$regions))
  expect_equal(unname(fe$phasing["detected"]), "chromosome")
  expect_true(all(fe$long$orient == "anchored"))
  expect_true(all(is.na(fe$long[n_hap == 0, f])))
  # atoms reproduce per-bin molecule counts
  fr = read_fragments(file.path(feat_cache$sim_norm_f_dir, "fragments.tsv.gz"))[ps > 0]
  m = collapse_molecules(fr[, .(cell, chrom, start, end, hap, ps)])
  m[, bin := locate_bins(s$scn$bins, chrom, start)]
  per_bin = m[, .N, by = bin][order(bin)]
  vb = fe$votes[, .(n = sum(nA + nB)), by = .(bin = fe$atoms$bin[match(atom, fe$atoms$atom)])][order(bin)]
  expect_equal(vb$n, per_bin$N); expect_equal(vb$bin, per_bin$bin)
  # tabix and the R fallback give the same votes
  v2 = read_votes(file.path(feat_cache$sim_norm_f_dir, "fragments.tsv.gz"), fe$cells$cell, fe$atoms, chroms = c("chr21", "chrX"), tabix = NULL)
  expect_equal(v2, fe$votes[atom %in% fe$atoms[chrom %in% c("chr21", "chrX"), atom]])
  # carriers of the private 1:0 / 2:0 events (truth copy numbers of the barcode; multiplets excluded)
  b = s$scn$bins; cn = s$truth$cn; cnA = s$truth$cnA; single = !s$truth$cells$multiplet
  cells = s$scn$cells$cell
  car = list(chr1p = colMeans(cn[b$chrom == "chr1" & b$end <= 24e6, ]) < 1.5,
             loh2 = colMeans((cn - cnA)[b$chrom == "chr2" & b$start >= 20e6, ]) < 0.5,
             chrX = colMeans(cn[b$chrom == "chrX", ]) < 1.5)
  for (rg in names(car)) {
    x = fe$long[region == rg][match(cells, cell)]
    expect_gt(sum(car[[rg]] & single), 3)
    expect_true(all(x$f[car[[rg]] & single] > 0.85), info = rg)
  }
  # non-carriers at chromosome scale: 0.5 +/- 0.05
  evchrom = list(chr1 = car$chr1p, chr2 = car$loh2, chrX = car$chrX,
                 chr8 = colMeans(cn[b$chrom == "chr8", ]) > 2.5, chr21 = colMeans(cn[b$chrom == "chr21", ]) > 2.5)
  for (ch in c("chr1", "chr2", "chr8", "chr21", "chr22", "chrX")) {
    non = if (is.null(evchrom[[ch]])) rep(TRUE, length(cells)) else !evchrom[[ch]]
    x = fe$long[region == ch][match(cells, cell)]
    expect_lt(abs(median(x$f[non], na.rm = TRUE) - 0.5), 0.05)
  }
})

test_that("orient_blocks: a cell's own votes never orient its blocks (leak test)", {
  atoms = fake_atoms()
  regions = fake_regions(tiles = c("10Mb" = 1e7))
  cl = rep(c("A", "B"), each = 20)
  v = fake_votes(40, atoms, clones = cl, events = data.table(clone = "B", chrom = "chr22", start = 0, end = NA, cnA = 1L, cnB = 0L),
                 mol_per_atom = 1)
  for (g in list(NULL, cl)) {
    o1 = orient_blocks(v, regions, atoms, groups = g)
    for (c0 in c(1L, 25L)) {
      v2 = copy(v)[cell == c0, `:=`(nA = 3L * nA, nB = 3L * nB)]
      o2 = orient_blocks(v2, regions, atoms, groups = g)
      a = o1[cell == c0][order(region)]; b2 = o2[cell == c0][order(region)]
      expect_equal(b2$region, a$region)
      expect_equal(b2$nA, 3 * a$nA); expect_equal(b2$nB, 3 * a$nB); expect_equal(b2$n_blocks, a$n_blocks)
    }
  }
  # clone-grouped LOO orientation recovers the clone B loss
  o = orient_blocks(v, regions, atoms, groups = cl)[region == "chr22"]
  o[, f := nA / (nA + nB)]
  expect_gt(min(o[cell > 20, f]), 0.9)
  expect_lt(abs(median(o[cell <= 20, f]) - 0.5), 0.1)
  expect_equal(detect_phasing(v, atoms), "block")
  # chromosome-wide phase sets plus a few residual small sets (as a real hybrid phase) stay "chromosome"
  vc = copy(v)[, ps := as.numeric(match(atoms$chrom[match(atom, atoms$atom)], unique(atoms$chrom)))]
  vc[sample(.N, 50), ps := 1000 + seq_len(50)]
  expect_equal(detect_phasing(vc, atoms), "chromosome")
  # mirrored: n_hap = all votes of the region, f in [0.5, 1)
  mi = allelic_features(v, atoms, regions[class == "chrom"], orientation = "mirrored", rho = 0.01)
  tot = merge(v, atoms[, .(atom, chrom)], by = "atom")[, .(n = sum(nA + nB)), by = .(cell, region = chrom)]
  mm = merge(mi, tot, by = c("cell", "region"))
  expect_equal(nrow(mm), nrow(tot)); expect_equal(mm$n_hap, as.numeric(mm$n))
  expect_true(all(mi$f >= 0.5 & mi$f < 1)); expect_equal(mi$orient[1], "none")
})

test_that("cell_features on sim_tum2 (male): chrX depth-only, clone-grouped orientation", {
  s = sim_tum2
  fe = sim_features("sim_tum2", method = "celltype", sex = "male", classes = c("arm", "chrom"), groups = "cell_type",
                    null_exclude = c("chr1", "chr2", "chr8", "chr21", "chr22"))
  expect_silent(check_features(fe$long))
  expect_equal(unname(fe$phasing), c("block", "loo"))
  x = fe$long[chrom == "chrX"]
  expect_gt(nrow(x), 0)
  expect_true(all(!x$allelic_ok) && all(x$n_hap == 0) && all(x$orient == "none"))
  expect_true(all(is.finite(x$z)))
  expect_false(any(fe$atoms$chrom == "chrX"))
  expect_equal(unique(fe$long[chrom != "chrX", orient]), "clone")
  single = !s$truth$cells$multiplet
  cells = s$scn$cells$cell; type = s$scn$cells$cell_type
  y = fe$long[region == "chr22"][match(cells, cell)]
  expect_gt(median(y$f[type == "cloneB" & single]), 0.8)
  expect_gt(mean(y$f[type == "cloneB" & single] > 0.7), 0.9)
  expect_lt(abs(median(y$f[type == "normal" & single]) - 0.5), 0.05)
})

# ---- population-phase orientation prior ------------------------------------------------------------------

# frozen copy of orient_layer() (the block orientation must stay bit-identical)
ref_orient_layer = function(b, groups, loo) {
  if (!nrow(b)) return(data.table(cell = integer(), ri = integer(), nA = numeric(), nB = numeric(), n_blocks = integer()))
  b[, grp := if (is.null(groups)) NA_character_ else as.character(groups[cell])]
  all_tot = b[, .(gA = sum(nA), gB = sum(nB)), by = .(ri, ps)]
  b[all_tot, `:=`(gA = i.gA, gB = i.gB), on = .(ri, ps)]
  if (any(!is.na(b$grp))) {
    grp_tot = b[!is.na(grp), .(hA = sum(nA), hB = sum(nB)), by = .(grp, ri, ps)]
    b[grp_tot, `:=`(gA = i.hA, gB = i.hB), on = .(grp, ri, ps)]
  }
  if (loo) b[, `:=`(rA = gA - nA, rB = gB - nB)] else b[, `:=`(rA = gA, rB = gB)]
  b = b[rA != rB]
  b[, flip := rA < rB]
  b[, .(nA = sum(ifelse(flip, nB, nA)), nB = sum(ifelse(flip, nA, nB)), n_blocks = .N), by = .(cell, ri)]
}

test_that("orient_blocks reproduces the reference leave-one-out orientation exactly", {
  atoms = fake_atoms(); regions = fake_regions(tiles = c("10Mb" = 1e7))
  cl = rep(c("A", "B"), each = 20)
  v = fake_votes(40, atoms, clones = cl, events = data.table(clone = "B", chrom = "chr22", start = 0, end = NA, cnA = 1L, cnB = 0L),
                 mol_per_atom = 0.3)
  map = scphasedcn:::atom_region_map(atoms, regions)
  for (g in list(NULL, cl)) {
    gg = scphasedcn:::check_groups(g, v)
    ref = rbindlist(lapply(sort(unique(map$layer)), function(l) ref_orient_layer(scphasedcn:::layer_blocks(v, map, l), gg, TRUE)))
    ref = ref[, .(cell, region = regions$region[ri], nA, nB, n_blocks)]
    o = orient_blocks(v, regions, atoms, groups = g)
    expect_identical(o, ref)
  }
  expect_false(exists("read_phase_prior", envir = asNamespace("scphasedcn"), inherits = FALSE))
})
