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

# ---- per-cell GC trend of `model`: lowess (default) and ratio of sums (gc_method = "ratio") ----------------------
# one-chromosome layout with smooth GC, a shared GC effect on accessibility exp(4 (gc - 0.42)), Poisson cells at a
# given mean count per bin; `eff` an optional per-cell log2 GC effect (+eff in odd cells, -eff in even cells: an
# effect shared by all cells is absorbed by the shared GLM)
gc_cells = function(lambda, n_cells = 20, n_bins = 28000, eff = NULL, seed = 1) {
  set.seed(seed)
  bins = make_bins(c(chrA = n_bins * 1e5), 1e5)
  bins[, gc := pmin(pmax(0.42 + as.numeric(stats::filter(rnorm(.N, 0, 0.09), rep(1 / 5, 5), circular = TRUE)), 0.31), 0.64)]
  a = if (is.null(eff)) rep(0, n_cells) else rep(c(1, -1), length.out = n_cells)
  P = exp(4 * (bins$gc - 0.42)) * 2^outer(if (is.null(eff)) rep(0, n_bins) else eff(bins$gc), a)
  lam = sweep(P, 2, colSums(P), "/") * lambda * n_bins
  X = matrix(as.numeric(rpois(length(lam), lam)), n_bins, n_cells)
  list(bins = bins, X = X, s = colSums(X), lam = lam)
}
# log2(expected / truth) per GC decile, both rescaled to the cell's total (as finish_norm does), mean over `cols`
gc_decile_bias = function(E, sm, cols = seq_len(ncol(E))) {
  E = sweep(E, 2, colSums(E) / sm$s, "/")[, cols, drop = FALSE]
  L = sweep(sm$lam, 2, colSums(sm$lam) / sm$s, "/")[, cols, drop = FALSE]
  q = cut(sm$bins$gc, stats::quantile(sm$bins$gc, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  rowMeans(rowsum(log2(E / L), q) / as.vector(table(q)))
}

test_that("the default per-cell GC step (lowess) is unchanged; gc_method unset = \"lowess\"", {
  sm = gc_cells(0.5, n_cells = 12, n_bins = 3000)
  E1 = suppressMessages(norm_model(sm$X, sm$s, sm$bins))
  E2 = suppressMessages(norm_model(sm$X, sm$s, sm$bins, gc_method = "lowess"))
  expect_identical(E1, E2)
  # the per-cell step as released in 0.1.1: a lowess of the bin log2 ratio on GC, applied to the shared model
  E0 = suppressMessages(norm_model(sm$X, sm$s, sm$bins, per_cell_gc = FALSE))
  for (i in seq_len(ncol(sm$X))) {
    fit = lowess_fit(sm$bins$gc, log2((sm$X[, i] + 0.5) / (E0[, i] + 0.5)), f = 0.3)
    ok = is.finite(fit)
    E0[ok, i] = E0[ok, i] * 2^fit[ok]
  }
  expect_identical(E1, E0)
  # through normalise_depth, and parameters of older runs (no gc_method) mean lowess
  scn = sim_small$scn
  n1 = suppressMessages(normalise_depth(scn, "model"))
  n2 = suppressMessages(normalise_depth(scn, "model", gc_method = "lowess"))
  n3 = suppressMessages(normalise_depth(scn, "model", gc_method = NULL))
  expect_identical(n1$expected, n2$expected); expect_identical(n1$expected, n3$expected)
  expect_identical(n1$params$gc_method, "lowess")
  n4 = suppressMessages(normalise_depth(scn, "model", gc_method = "ratio"))
  expect_identical(n4$params$gc_method, "ratio")
  expect_false(identical(n1$expected, n4$expected))
  expect_equal(colSums(n4$expected, na.rm = TRUE), n4$size, ignore_attr = TRUE)   # cell totals preserved
  expect_error(normalise_depth(scn, "model", gc_method = "spline"))
  expect_identical(caller_defaults()$gc_method, "lowess")
})

test_that("ratio-of-sums GC trend has no depth-dependent GC bias where the lowess has one", {
  # no per-cell GC effect: the expected counts should follow the truth at every depth (mean counts per bin 0.17,
  # 0.27, 1.7 = about 5k, 8k, 50k fragments over 30k bins)
  spread = sapply(c(0.17, 0.27, 1.7), function(l) {
    sm = gc_cells(l, seed = 3)
    b_low = gc_decile_bias(suppressMessages(norm_model(sm$X, sm$s, sm$bins, gc_method = "lowess")), sm)
    b_rat = gc_decile_bias(suppressMessages(norm_model(sm$X, sm$s, sm$bins, gc_method = "ratio")), sm)
    c(lowess = unname(b_low[10] - b_low[1]), ratio = unname(b_rat[10] - b_rat[1]), ratio_max = max(abs(b_rat)))
  })
  expect_true(all(abs(spread["ratio", ]) < 0.03), info = paste(round(spread["ratio", ], 3), collapse = " "))
  expect_true(all(spread["ratio_max", ] < 0.03))
  # the lowess: a GC trend that changes with depth
  expect_gt(max(spread["lowess", ]) - min(spread["lowess", ]), 0.1)
})

test_that("ratio-of-sums GC trend recovers a known per-cell GC effect", {
  eff = function(gc) 3 * (gc - 0.42) + 10 * (gc - 0.42)^2
  for (l in c(0.17, 1.7)) {
    sm = gc_cells(l, eff = eff, seed = 4)
    plus = seq(1, ncol(sm$X), 2)                          # the +eff cells
    b0 = gc_decile_bias(suppressMessages(norm_model(sm$X, sm$s, sm$bins, per_cell_gc = FALSE)), sm, plus)
    b = gc_decile_bias(suppressMessages(norm_model(sm$X, sm$s, sm$bins, gc_method = "ratio")), sm, plus)
    expect_gt(max(abs(b0)), 0.15)                         # the shared model alone misses the effect
    expect_lt(max(abs(b)), 0.05, label = sprintf("lambda %.2f: %s", l, paste(round(b, 3), collapse = " ")))
  }
})

test_that("ratio-of-sums GC trend handles empty cells, bins without GC and constant GC", {
  sm = gc_cells(0.3, n_cells = 5, n_bins = 2000)
  sm$X[, 2] = 0; sm$s = colSums(sm$X)
  E0 = suppressMessages(norm_model(sm$X, sm$s, sm$bins, per_cell_gc = FALSE))
  E = suppressMessages(norm_model(sm$X, sm$s, sm$bins, gc_method = "ratio"))
  expect_true(all(is.finite(E)))
  expect_equal(E[, 2], E0[, 2])                           # no fragments: no trend
  # bins without GC are left out of the groups (they keep the shared model)
  gc = sm$bins$gc; gc[c(5, 500)] = NA
  tr = gc_group_trend(sm$X, E0, gc, rep(1, nrow(sm$X)))
  expect_false(any(c(5, 500) %in% tr$rows))
  expect_length(tr$rows, nrow(sm$X) - 2)
  expect_true(all(is.finite(tr$F)))
  expect_equal(unname(tr$F[, 2]), rep(0, nrow(tr$F)))
  # constant GC: no groups, no trend
  expect_null(gc_group_trend(sm$X, E0, rep(0.4, nrow(sm$X)), rep(1, nrow(sm$X))))
  # few groups (heavy ties): unsmoothed group ratios, finite
  tr = gc_group_trend(sm$X, E0, rep(c(0.4, 0.5, 0.6), length.out = nrow(sm$X)), rep(1, nrow(sm$X)))
  expect_true(all(is.finite(tr$F)))
  expect_true(all(tr$w >= 0 & tr$w <= 1))
})
