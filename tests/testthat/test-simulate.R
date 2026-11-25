# test-simulate.R: the simulator's seeding (seed controls the data, the genome its own seed, the caller's stream kept)
#                  and its het-SNP positions (one SNP per position).

# integer-valued fingerprint of a simulation: fragments, depth counts and truth (drawn before the SNPs), then the SNPs
# and the allele / haplotype counts
sim_fingerprint = function(s) {
  fr = s$fragments
  cs = Matrix::summary(as(s$scn$counts, "CsparseMatrix"))
  sp = s$snps
  c(n_frag = nrow(fr), frag_start = sum(as.numeric(fr$start)), frag_end = sum(as.numeric(fr$end)),
    frag_mix = sum((fr$start %% 1009) * (seq_len(nrow(fr)) %% 1013)),
    counts = sum(cs$x), counts_ij = sum(as.numeric(cs$i) * cs$j * cs$x), cn = sum(s$truth$cn_nucleus),
    n_snp = nrow(sp), snp_pos = sum(as.numeric(sp$pos)), snp_hap = sum((sp$hap_of_alt == "A") * (seq_len(nrow(sp)) %% 1013)),
    dp = sum(s$scn$alleles$DP), ad = sum(s$scn$alleles$AD), hapA = sum(s$scn$layers$hapA), hapB = sum(s$scn$layers$hapB))
}
pre_snp = c("n_frag", "frag_start", "frag_end", "frag_mix", "counts", "counts_ij", "cn")

test_that("simulate_scn: seed controls the data with the default genome; the genome does not depend on seed", {
  a = simulate_scn(n_cells = 20, seed = 1)
  b = simulate_scn(n_cells = 20, seed = 2)
  expect_false(identical(a$scn$counts, b$scn$counts))
  expect_false(identical(a$fragments$start, b$fragments$start))
  # the genome is sim_genome() (its own seed 1) whatever the simulation seed
  expect_identical(a$scn$bins, b$scn$bins)
  expect_equal(a$scn$bins, sim_genome(), ignore_attr = TRUE)
  # the same genome value gives the same data, passed precomputed or evaluated in the call
  g = sim_genome()
  expect_identical(simulate_scn(n_cells = 20, genome = g, seed = 2)$scn$counts, b$scn$counts)
  # seed 1 with the default genome: the data it gave before the fix (the part drawn before the SNPs; the SNP part
  # differs only at former duplicate positions, see below)
  expect_equal(sim_fingerprint(a)[pre_snp],
               c(n_frag = 73358, frag_start = 1550714655678, frag_end = 1550776049591, frag_mix = 18653385374,
                 counts = 73358, counts_ij = 725615605, cn = 106800))
})

test_that("sim_genome / simulate_scn leave the caller's random stream alone; an inline genome uses its own seed", {
  set.seed(5); x = runif(3)
  set.seed(5); g2 = sim_genome(seed = 2); expect_identical(runif(3), x)
  set.seed(5); invisible(simulate_scn(n_cells = 5, seed = 3)); expect_identical(runif(3), x)
  s = simulate_scn(n_cells = 5, genome = sim_genome(seed = 2), seed = 3)
  expect_equal(s$scn$bins, g2, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(g2$gc, sim_genome()$gc)))
})

test_that("unique_positions: draws without collisions unchanged; collisions moved to the next free position", {
  up = scphasedcn:::unique_positions
  expect_identical(up(c(0, 5, 9), 10), c(0, 5, 9))
  expect_identical(up(c(3, 3, 3, 4, 8), 100), c(3, 4, 5, 6, 8))
  expect_identical(up(c(0, 7, 9, 9), 10), c(0, 7, 8, 9))      # past the end: moved down from size - 1
  expect_identical(up(c(9, 9, 9), 3), c(0, 1, 2))
  expect_identical(up(numeric(0), 10), numeric(0))
  expect_identical(up(4, 10), 4)
  expect_error(up(c(1, 1, 1), 2), "do not fit")
  set.seed(1)
  p = sort(floor(runif(3000) * 60000))
  expect_true(anyDuplicated(p) > 0)
  q = up(p, 60000)
  expect_false(anyDuplicated(q) > 0)
  expect_true(all(diff(q) > 0) && min(q) >= 0 && max(q) <= 59999)
  expect_true(all(q >= p))                                    # moved up only (no collision near the end here)
  expect_lte(sum(q != p), 2 * sum(duplicated(p)))             # only the collisions (and rarely a neighbour) move
})

test_that("simulate_scn: one SNP per position (unique chrom:pos ids); a draw without duplicates is unchanged", {
  a = simulate_scn(n_cells = 20, seed = 1)
  expect_false(anyDuplicated(a$snps$id) > 0)
  expect_false(anyDuplicated(a$scn$alleles$snps$snp) > 0)
  expect_true(a$snps[, all(diff(pos) > 0), by = chrom][, all(V1)])
  # at a low SNP density the seed-1 draw has no duplicate: the whole simulation is the one from before the change
  b = simulate_scn(n_cells = 20, seed = 1, snp_density = 1 / 2e5)
  expect_equal(sim_fingerprint(b),
               c(n_frag = 73358, frag_start = 1550714655678, frag_end = 1550776049591, frag_mix = 18653385374,
                 counts = 73358, counts_ij = 725615605, cn = 106800, n_snp = 995, snp_pos = 20166708141,
                 snp_hap = 258894, dp = 322, ad = 172, hapA = 186, hapB = 136))
})
