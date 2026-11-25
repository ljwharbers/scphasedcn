# test-allelic.R: beta-binomial primitives and switch-error attenuation.
test_that("bb_loglik reduces to the binomial at rho = 0 and fit_bb_rho recovers overdispersion", {
  expect_equal(bb_loglik(3, 10, 0.4, rho = 0), dbinom(3, 10, 0.4, log = TRUE))
  set.seed(5)
  n = rep(60, 4000)
  p = rbeta(length(n), 0.5 * (1 - 0.05) / 0.05, 0.5 * (1 - 0.05) / 0.05)
  a = rbinom(length(n), n, p)
  expect_lt(abs(fit_bb_rho(a, n) - 0.05), 0.01)
  expect_lt(fit_bb_rho(rbinom(length(n), n, 0.5), n), 0.005)
  expect_setequal(imbalance_states()$state, c("1:1", "1:0", "2:1", "3:1", "4:1", "2:0"))
})

test_that("switch errors attenuate the anchored imbalance; theory and simulation agree in direction", {
  expect_lt(switch_attenuation(0.02, 750), switch_attenuation(1e-4, 750))
  expect_gt(switch_attenuation(1e-4, 750), 0.85)
  s_bad = simulate_scn(n_cells = 80, seed = 1, multiplet_rate = 0, switch_rate = 0.02)
  s_good = simulate_scn(n_cells = 80, seed = 1, multiplet_rate = 0, switch_rate = 1e-4)
  # hap-A fraction of the phased fragments on chr22 (monosomy 1:0 in the tumour cells)
  f_hap = function(s) {
    b = s$scn$bins$chrom == "chr22"
    nA = colSums(s$scn$layers$hapA[b, , drop = FALSE]); nB = colSums(s$scn$layers$hapB[b, , drop = FALSE])
    ifelse(nA + nB > 0, nA / (nA + nB), NA_real_)
  }
  tum = s_bad$scn$cells$cell_type == "tumour"
  dev_bad = median(abs(f_hap(s_bad)[tum] - 0.5), na.rm = TRUE); dev_good = median(abs(f_hap(s_good)[tum] - 0.5), na.rm = TRUE)
  expect_lt(dev_bad, dev_good - 0.1)
  expect_gt(dev_good, 0.35)
})
