# test-normalise.R: the depth normalisations (loo, celltype, model) run; leave-one-out references do not leak the cell's own counts.
test_that("the depth normalisations return a full ratio matrix; the older method names are aliases", {
  scn = sim_small$scn
  clean = clean_cells(sim_small)
  for (m in c("loo", "celltype", "model")) {
    nm = suppressMessages(normalise_depth(scn, m, ref_cells = clean))
    expect_s3_class(nm, "scn_norm")
    expect_equal(dim(nm$ratio), dim(scn$counts))
    expect_true(all(is.finite(nm$ratio[nm$use, ])), info = m)
    expect_true(all(is.na(nm$ratio[!nm$use, ])), info = m)
  }
  expect_identical(normalise_depth(scn, "N2")$ratio, normalise_depth(scn, "celltype")$ratio)
  expect_error(normalise_depth(scn, "N3"))
})

test_that("LOO references of loo / celltype do not depend on the cell's own counts", {
  scn = sim_flat$scn
  scn$cells[, cell_type := rep(c("a", "b"), length.out = .N)]
  i = 5
  b = which(scn$bins$chrom == "chr2")
  scn2 = scn
  scn2$counts[b, i] = scn2$counts[b, i] * 2
  # chromosome-level log2 ratio: summed counts over summed expected on the used bins
  chrom_ratio = function(s, norm) {
    ids = s$bins$chrom; ids[!norm$use] = NA
    ind = group_indicator(ids)
    X = as.matrix(ind %*% s$counts)
    E = as.matrix(ind %*% replace(norm$expected, is.na(norm$expected), 0))
    log2((X + 0.5) / (E + 0.5))
  }
  for (m in c("loo", "celltype")) {
    n1 = normalise_depth(scn, m); n2 = normalise_depth(scn2, m)
    # the reference profile of cell i (expected / its total) must not move when its own counts change
    p1 = n1$expected[, i] / sum(n1$expected[, i], na.rm = TRUE); p2 = n2$expected[, i] / sum(n2$expected[, i], na.rm = TRUE)
    u = n1$use
    expect_lt(sum(abs(p2[u] - p1[u])) / sum(p1[u]), 1e-10, label = m)
    # the doubled chromosome shows as a +1 log2 contrast against the rest of the genome,
    # relative to the cell's own (noise) offset before doubling
    R = chrom_ratio(scn2, n2); R1 = chrom_ratio(scn, n1)
    contrast = function(M) M["chr2", i] - median(M[rownames(M) != "chr2", i])
    expect_lt(abs((contrast(R) - contrast(R1)) - 1), 0.1, label = m)
    expect_lt(abs(contrast(R1)), 0.3, label = m)
  }
})

test_that("normalise_extra_bins adds bins with the reference of the existing normalisation, other bins untouched", {
  s = simulate_scn(n_cells = 40, genome = sim_genome(c(chr1 = 20e6, chrX = 10e6, chrY = 10e6)), sex = "male", events = NULL, seed = 4)
  scn = s$scn
  y = scn$bins$chrom == "chrY"
  use0 = bin_filter(scn$bins, scn$counts) & !y
  for (m in c("loo", "celltype")) {
    n0 = normalise_depth(scn, m, use = use0)
    n1 = normalise_extra_bins(scn, n0, y)
    expect_identical(n1$expected[!y, ], n0$expected[!y, ])
    expect_identical(n1$use, use0 | y)
    # chrY: the same leave-one-out reference with each cell's size factor of the other bins
    X = as.matrix(scn$counts[y, ]); sz = n0$size
    if (m == "loo") {
      ref = (rowSums(X) - X) / rep(sum(sz) - sz, each = nrow(X)) * rep(sz, each = nrow(X))
      expect_equal(unname(n1$expected[y, ]), unname(ref))
    }
    expect_equal(stats::median(colSums(X) / colSums(n1$expected[y, ])), 1, tolerance = 0.1)
  }
  expect_error(normalise_extra_bins(scn, normalise_depth(scn, "model", use = use0), y), "loo / celltype")
})
