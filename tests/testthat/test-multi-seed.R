# test-multi-seed.R: the main normal-mode and tumour-mode fixtures (sim_norm_f, sim_tum2) re-simulated at other seeds
# (2, 3, 4; same parameters, only the seed differs), with the key assertions of their seed-1 tests that held at every
# seed when this file was written. The other assertions of those tests were measured at seeds 1-4 too; the ones that
# fail at some seed are not repeated here (multiplet carriers diluting a private event, a sub-arm CN-LOH with neutral
# depth, near-zero rho ratios, a quarter-integer baseline): they are tracked separately. ~50 s per seed; skipped on CRAN
# like the other slow end-to-end tests.

ms_seeds = c(2L, 3L, 4L)

# the fixtures' parameters (helper-caller.R), seed free
ms_sim_norm_f = function(seed) simulate_scn(
  n_cells = 200, sex = "female", cell_types = c(epi = 0.4, fib = 0.35, imm = 0.25),
  depth = c(meanlog = log(15000), sdlog = 0.35), events = sim_norm_f_events(), private_frac = 0.05,
  multiplet_rate = 0.05, phase_block_mean = 1e9, switch_rate = 2e-6, block_flip_prob = 0, seed = seed)
ms_sim_tum2 = function(seed) simulate_scn(
  n_cells = 200, sex = "male", cell_types = c(normal = 0.3, cloneA = 0.35, cloneB = 0.35),
  depth = c(meanlog = log(8000), sdlog = 0.35), events = sim_tum2_events(),
  multiplet_rate = 0.05, phase_block_mean = 2e5, block_flip_prob = 0.5, seed = seed)

test_that("multi-seed: the re-simulation functions use the fixtures' parameters", {
  skip_on_cran()
  a = ms_sim_norm_f(1L)
  expect_identical(a$scn$counts, sim_norm_f$scn$counts)
  expect_identical(a$fragments, sim_norm_f$fragments)
  b = ms_sim_tum2(1L)
  expect_identical(b$scn$counts, sim_tum2$scn$counts)
  expect_identical(b$fragments, sim_tum2$fragments)
})

test_that("multi-seed sim_norm_f: null-cell selection drops the private-event carriers; nested events resolved", {
  skip_on_cran()
  for (s in ms_seeds) {
    sim = ms_sim_norm_f(s)
    sf = sim_features(sim)
    tr = sim$truth
    cells = sim$scn$cells$cell
    car = lapply(tr$carriers, function(n) unique(cells[tr$nuclei$cell[n]]))     # 1p loss, chr8 gain, chr2 CN-LOH, chr21 gain, chrX loss
    carriers = unique(unlist(car))
    mult = tr$cells[multiplet == TRUE, cell]
    # select_null_cells (all chromosomes): >= 90 % of the singlet carriers out, >= 90 % of the other cells in
    sel = select_null_cells(sf, exclude = mult, chroms = NULL)
    sc = setdiff(carriers, mult)
    expect_gt(length(sc), 20, label = paste("seed", s, "singlet carriers"))
    expect_gte(mean(!sc %in% sel$null_cells), 0.9, label = paste("seed", s, "carriers dropped"))
    expect_gte(mean(setdiff(unique(sf$cell), c(carriers, mult)) %in% sel$null_cells), 0.9, label = paste("seed", s, "non-carriers kept"))
    expect_false(any(mult %in% sel$null_cells), label = paste("seed", s, "multiplets in the null"))
    expect_equal(nrow(sel$log), 2, label = paste("seed", s, "selection iterations"))
    # resolve_nested: no chromosome gain on carriers of the nested 1p loss / chr2 CN-LOH; the chr8 gain is found
    rg = fake_regions(sim$scn$bins, tiles = c("20Mb" = 2e7))
    r1 = caller_normal(sf, regions = rg, multiplets = mult, n_draws = 500, pool_classes = character())
    rep_cells = function(r, reg, st) unique(r$calls[level == "cell" & called & reportable & region == reg & state == st, cell])
    expect_length(intersect(rep_cells(r1, "chr1", "gain"), car[[1]]), 0)
    expect_length(intersect(rep_cells(r1, "chr2", "gain"), car[[3]]), 0)
    expect_equal(nrow(r1$events[level == "cell" & cell %in% c(car[[1]], car[[3]]) & chrom %in% c("chr1", "chr2") & state == "gain"]), 0,
                 label = paste("seed", s, "gain events on nested carriers"))
    expect_gte(mean(car[[2]] %in% rep_cells(r1, "chr8", "gain")), 0.7, label = paste("seed", s, "chr8 gain recall"))
  }
})

test_that("multi-seed sim_tum2: caller_tumour finds the clones, their states and the normal clone", {
  skip_on_cran()
  for (s in ms_seeds) {
    sim = ms_sim_tum2(s)
    res = caller_tumour(clone_sim_input(sim), verbose = FALSE)
    lab = function(x) paste("seed", s, x)
    a = merge(res$assignments, sim_truth(sim), by = "cell")
    sing = a[mult == FALSE]
    expect_gte(clone_ari(sing$clone, sing$type), 0.9, label = lab("ARI"))
    map = sing[, .N, by = .(clone, type)][order(-N)][!duplicated(clone)]
    expect_length(unique(map$type), 3)
    cl = function(type_) map[type == type_, clone]
    sg = res$segments
    at = function(k, ch, pos) sg[clone == k & chrom == ch & start <= pos & end > pos]
    expect_equal(at(cl("cloneB"), "chr22", 10e6)$state, "1:0", label = lab("cloneB chr22"))
    expect_equal(at(cl("cloneA"), "chr21", 10e6)$state, "2:1", label = lab("cloneA chr21"))
    expect_equal(at(cl("cloneB"), "chr1", 30e6)$state, "2:0", label = lab("cloneB chr1"))
    for (k in c(cl("cloneA"), cl("cloneB"))) {
      s8 = at(k, "chr8", 22e6)
      expect_equal(s8$label, "amp", label = lab("chr8 amp"))
      expect_lte(abs(s8$start - 15e6), 2e6, label = lab("chr8 amp start"))
    }
    s2 = at(cl("cloneB"), "chr2", 17e6)
    expect_equal(s2$label, "loss", label = lab("cloneB chr2 loss"))
    expect_lte(abs(s2$start - 10e6), 2e6, label = lab("chr2 loss start"))
    expect_lte(abs(s2$end - 25e6), 2e6, label = lab("chr2 loss end"))
    expect_true(all(sg[clone == cl("normal"), label] == "neutral"), label = lab("normal clone neutral"))
    expect_true(res$clones[clone == cl("normal"), normal], label = lab("normal clone flagged"))
    expect_true(all(res$clones$baseline == 2), label = lab("diploid baselines"))
    expect_gte(sing[, mean(map$type[match(clone, map$clone)] == type)], 0.95, label = lab("assignment accuracy"))
    pool = a[ref_pool == TRUE]
    expect_gt(nrow(pool), 10, label = lab("reference pool size"))
    expect_gte(mean(pool$type == "normal" & !pool$mult), 0.98, label = lab("reference pool purity"))
    cc = res$calls[level == "cell"]
    expect_gt(cc[grepl("chr22", region), mean(source == "confirmed")], 0.8, label = lab("chr22 confirmed"))
  }
})
