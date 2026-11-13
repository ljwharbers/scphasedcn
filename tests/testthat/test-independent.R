# test-independent.R: the scoring primitives against references written out here, independently of the package code
# (the parity oracle of helper-refscorer.R calls the package's own bb_loglik() and copies its depth null fit, so it
# cannot catch an error in either).

# beta-binomial log-pmf with mean f and overdispersion rho: alpha = f (1 - rho) / rho, beta = (1 - f) (1 - rho) / rho,
# P(a | n) = choose(n, a) B(a + alpha, n - a + beta) / B(alpha, beta), written with rising factorials
# (B(x + k, y + m) / B(x, y) = prod_{i<k} (x + i) prod_{i<m} (y + i) / prod_{i<k+m} (x + y + i)); f is first clamped to
# [eps, 1 - eps] (the allele error floor of bb_loglik()); rho = 0 is the binomial
ref_bb = function(a, n, f, rho, eps = 0.02) {
  f = min(max(f, eps), 1 - eps)
  lc = lgamma(n + 1) - lgamma(a + 1) - lgamma(n - a + 1)
  if (rho == 0) return(lc + a * log(f) + (n - a) * log(1 - f))
  al = f * (1 - rho) / rho
  be = (1 - f) * (1 - rho) / rho
  rf = function(x, k) if (k == 0) 0 else sum(log(x + 0:(k - 1)))
  lc + rf(al, a) + rf(be, n - a) - rf(al + be, n)
}

test_that("bb_loglik equals the beta-binomial log-pmf (mean f, overdispersion rho, f clamped to [eps, 1 - eps])", {
  g = expand.grid(n = c(0, 1, 5, 23, 80), f = c(0.5, 2 / 3, 0.75, 0.98, 1, 0.01), rho = c(1e-4, 0.002, 0.02, 0.3, 0.5, 0))
  for (i in seq_len(nrow(g))) {
    n = g$n[i]
    for (a in unique(c(0, floor(n / 3), n))) {
      expect_equal(bb_loglik(a, n, g$f[i], g$rho[i]), ref_bb(a, n, g$f[i], g$rho[i]), tolerance = 1e-10,
                   info = paste(a, n, g$f[i], g$rho[i]))
      expect_equal(bb_loglik(a, n, g$f[i], g$rho[i], eps = 0.001), ref_bb(a, n, g$f[i], g$rho[i], eps = 0.001),
                   tolerance = 1e-10, info = paste("eps 0.001", a, n, g$f[i], g$rho[i]))
    }
  }
  # n = 0: the only outcome has probability 1
  expect_equal(bb_loglik(0, 0, 0.8, 0.01), 0)
  # vectorised over rows
  a = c(0, 3, 7); n = c(4, 9, 7)
  expect_equal(bb_loglik(a, n, 0.75, 0.02), mapply(ref_bb, a, n, 0.75, 0.02), tolerance = 1e-10)
})

test_that("ll_allelic: two-sided and one-sided LLR against 1:1 from the reference log-pmf", {
  rows = data.table(nA = c(0, 3, 12, 2, 40, 5), n = c(0, 4, 15, 9, 50, 5))
  for (rho in c(1e-4, 0.01, 0.2)) for (f in c(2 / 3, 0.75, 0.98)) {
    two = mapply(function(a, n) if (n == 0) 0 else max(ref_bb(a, n, f, rho), ref_bb(a, n, 1 - f, rho)) - ref_bb(a, n, 0.5, rho),
                 rows$nA, rows$n)
    one = mapply(function(a, n) if (n == 0) 0 else ref_bb(a, n, f, rho) - ref_bb(a, n, 0.5, rho), rows$nA, rows$n)
    expect_equal(ll_allelic(rows$nA, rows$n, f, rho), two, tolerance = 1e-10, info = paste(rho, f))
    expect_equal(ll_allelic(rows$nA, rows$n, f, rho, two_sided = FALSE), one, tolerance = 1e-10, info = paste(rho, f))
    expect_equal(ll_allelic(rows$nA, rows$n, f, rho, two_sided = FALSE, per_row = TRUE), one, tolerance = 1e-10)
  }
  # a grid of f: the best LLR over the grid
  grid = c(0.6, 0.8, 0.95)
  best = mapply(function(a, n) if (n == 0) 0 else max(sapply(grid, function(f) max(ref_bb(a, n, f, 0.01), ref_bb(a, n, 1 - f, 0.01)))) -
                  ref_bb(a, n, 0.5, 0.01), rows$nA, rows$n)
  expect_equal(ll_allelic(rows$nA, rows$n, grid, 0.01), best, tolerance = 1e-10)
})

test_that("ll_depth_shift: log t density ratio of a location shift, written out", {
  # location-scale Student-t log density
  ldt = function(x, m, s, df) lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * log(df * pi) - log(s) -
    (df + 1) / 2 * log(1 + ((x - m) / s)^2 / df)
  cases = data.table(z = c(0, -1.3, 2.4, 5, -7), mu = c(-1, 0.6, -2.2, 3, 1.5), m = c(0, 0.1, -0.2, 0.05, 0),
                     s = c(1, 1.3, 0.8, 2, 1.1), df = c(5, 30, 3, 100, 7))
  for (i in seq_len(nrow(cases))) {
    x = cases[i]
    expect_equal(ll_depth_shift(x$z, x$mu, list(t_m = x$m, t_s = x$s, t_df = x$df)),
                 ldt(x$z - x$mu, x$m, x$s, x$df) - ldt(x$z, x$m, x$s, x$df), tolerance = 1e-10, info = i)
  }
  # one row of parameters per z
  expect_equal(ll_depth_shift(cases$z, cases$mu, cases[, .(t_m = m, t_s = s, t_df = df)]),
               ldt(cases$z - cases$mu, cases$m, cases$s, cases$df) - ldt(cases$z, cases$m, cases$s, cases$df), tolerance = 1e-10)
  # Gaussian: log N(z - mu) - log N(z) = mu z - mu^2 / 2; NA z gives 0
  expect_equal(ll_depth_shift(c(1.7, NA), -0.9), c(dnorm(1.7 + 0.9, log = TRUE) - dnorm(1.7, log = TRUE), 0), tolerance = 1e-12)
})

test_that("orient_layer: a hand-computed leave-one-out orientation (3 cells, 2 blocks)", {
  # block 1 (ps 1): cell votes A / B = 5/1, 4/0, 0/3 (totals 9 / 4); block 2 (ps 2): 0/2, 1/3, 2/2 (totals 3 / 7).
  # Leave-one-out pseudobulk of block 1: cell 1 4/3, cell 2 5/4, cell 3 9/1 -> A ahead for all, kept as is.
  # Block 2: cell 1 3/5, cell 2 2/4, cell 3 1/5 -> B ahead for all, flipped (A and B swapped).
  # Oriented sums: cell 1 (5 + 2) / (1 + 0) = 7 / 1, cell 2 (4 + 3) / (0 + 1) = 7 / 1, cell 3 (0 + 2) / (3 + 2) = 2 / 5.
  b = data.table(cell = c(1L, 2L, 3L, 1L, 2L, 3L), ri = 1L, ps = c(1, 1, 1, 2, 2, 2),
                 nA = c(5, 4, 0, 0, 1, 2), nB = c(1, 0, 3, 2, 3, 2))
  o = scphasedcn:::orient_layer(copy(b), NULL, TRUE)[order(cell)]
  expect_equal(o$cell, 1:3)
  expect_equal(o$nA, c(7, 7, 2))
  expect_equal(o$nB, c(1, 1, 5))
  expect_equal(o$n_blocks, c(2L, 2L, 2L))
  # a block whose leave-one-out pseudobulk is tied is dropped: cell 1 alone on block 3 (pseudobulk 0 / 0)
  o3 = scphasedcn:::orient_layer(rbind(copy(b), data.table(cell = 1L, ri = 1L, ps = 3, nA = 4, nB = 0)), NULL, TRUE)[order(cell)]
  expect_equal(o3$nA, c(7, 7, 2))
  expect_equal(o3$n_blocks, c(2L, 2L, 2L))
})

test_that("score_cells on depth-only rows: no purely allelic state, s_comb is the best depth term", {
  # Gaussian depth model: dl = mu z - mu^2 / 2 with mu = shift / max(sd_null, sd_floor); no votes, so every allelic
  # term is 0 and, with allelic_ok FALSE, CN-LOH (no depth shift) is excluded even when both depth terms are negative
  p = caller_defaults()
  d = data.table(cell = "c", region = c("r1", "r2", "r3"), class = "chrom", stratum = "s", z = c(-0.2, -3, 1.8),
                 sd_null = c(0.3, 0.02, 0.4), nA = 0, n_hap = 0, allelic_ok = FALSE)
  sc = score_cells(d, caller_states("normal"), rho = 0.01, depth_model = "gauss")
  sdp = pmax(d$sd_null, p$sd_floor)
  dl = function(shift) { mu = shift / sdp; mu * d$z - mu^2 / 2 }
  dl_loss = dl(p$shifts[["loss"]]); dl_gain = dl(p$shifts[["gain"]])
  expect_true(all(is.infinite(sc$cl_cnloh) & sc$cl_cnloh < 0))
  expect_false(any(sc$state_comb == "cnloh"))
  expect_equal(sc$s_comb, pmax(dl_loss, dl_gain), tolerance = 1e-12)
  expect_equal(sc$state_comb, ifelse(dl_loss >= dl_gain, "loss", "gain"))
  expect_lt(sc$s_comb[1], 0)                                     # both depth terms negative: CN-LOH would have won at 0
})
