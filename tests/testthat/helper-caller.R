# helper-caller.R
# What: shared test data for the caller modules (caller_core, features, caller, clones, call_scn):
#       two cached simulations and generators of fake caller tables.
#       - sim_norm_f: normal tissue, female, 3 cell types, hybrid-like phasing (one block per chromosome),
#         cell-private events only.
#       - sim_tum2: male tumour, two clones with a shared trunk, unoriented short phase blocks.
#       - sim_tum3: sim_tum2 at half the depth with in-block switch errors (clone allele attenuation a ~ 0.7).
#       - fake_regions(), fake_features(), fake_atoms(), fake_votes(): tables in the package's table
#         schemas (check_tables) with known truth, for unit tests that do not need the full simulator.
# Why: the tests of the feature step, scorer, tumour mode and pipeline share the same data. The simulations are
#      built lazily (active bindings, first use only, then cached), so test files that do not use them
#      pay nothing; both together take ~20 s.
# Inputs: none. Outputs: objects in the test environment.

caller_test_cache = new.env()
lazy_cached = function(name, build, env = parent.frame()) {
  makeActiveBinding(name, function() {
    if (is.null(caller_test_cache[[name]])) caller_test_cache[[name]] = build()
    caller_test_cache[[name]]
  }, env)
}

# ---- simulations ---------------------------------------------------------------------------------

# Cell-private events of the normal-tissue simulation. The simulator draws carriers of a "private"
# event independently per event: each nucleus carries it with probability `private_frac` (0.05),
# so ~10 of 200 barcodes carry each event (a few carry two).
sim_norm_f_events = function() {
  data.table(chrom = c("chr1", "chr8", "chr2", "chr21", "chrX"),
             start = c(0, 0, 20e6, 0, 0), end = c(24e6, NA, NA, NA, NA),
             cnA = c(1L, 2L, 2L, 2L, 1L), cnB = c(0L, 1L, 0L, 1L, 0L), type = "private")
}
# chr1p loss (1:0), chr8 trisomy (2:1), chr2 distal CN-LOH (2:0, 20-40 Mb), chr21 trisomy (2:1),
# chrX loss (1:0; female). Default sim_genome(): chr1 60 Mb (p arm 0-24 Mb), chr2 40, chr8 30,
# chr21 20, chr22 20, chrX 40 Mb; 100 kb bins.
# Hybrid (longphase + population) phasing: phase_block_mean = 1e9 gives one anchored block per chromosome
# (seed 101: chrX has one extra break). The default switch_rate (1e-4 per SNP) would put ~5 switches on a
# 30,000-SNP chromosome and scramble chromosome-wide orientation; 2e-6 (~1 switch per 500,000 SNPs) keeps the
# anchored haplotype fraction intact (hap agrees with hap_true for > 90 % of phased fragments).
# multiplet_rate 0.05: a realistic residual rate after the soft multiplet filter.
lazy_cached("sim_norm_f", function() simulate_scn(
  n_cells = 200, sex = "female", cell_types = c(epi = 0.4, fib = 0.35, imm = 0.25),
  depth = c(meanlog = log(15000), sdlog = 0.35), events = sim_norm_f_events(), private_frac = 0.05,
  multiplet_rate = 0.05, phase_block_mean = 1e9, switch_rate = 2e-6, block_flip_prob = 0, seed = 101))

# Clone truth of the tumour simulation: `type` = the cell type (clone) that carries the event. The
# simulator applies one event row to one cell type, so the trunk event (chr8 distal 3:1 in both clones)
# is written as two identical rows, one per clone. chr2 10-25 Mb is interstitial; chr1 2:0 is
# whole-chromosome CN-LOH.
sim_tum2_events = function() {
  data.table(chrom = c("chr8", "chr8", "chr21", "chr22", "chr2", "chr1"),
             start = c(15e6, 15e6, 0, 0, 10e6, 0), end = c(NA, NA, NA, NA, 25e6, NA),
             cnA = c(3L, 3L, 2L, 1L, 1L, 2L), cnB = c(1L, 1L, 1L, 0L, 0L, 0L),
             type = c("cloneA", "cloneB", "cloneA", "cloneB", "cloneB", "cloneB"))
}
# Unoriented longphase-like blocks: mean block 200 kb, each block flipped with probability 0.5; default
# switch_rate (1e-4 per SNP, ~100 SNPs per block). Depth: meanlog log(8000).
lazy_cached("sim_tum2", function() simulate_scn(
  n_cells = 200, sex = "male", cell_types = c(normal = 0.3, cloneA = 0.35, cloneB = 0.35),
  depth = c(meanlog = log(8000), sdlog = 0.35), events = sim_tum2_events(),
  multiplet_rate = 0.05, phase_block_mean = 2e5, block_flip_prob = 0.5, seed = 102))
# sim_tum3: sim_tum2's clones and events, attenuated allele signal and half the depth (a hard confirmation regime):
# switch_rate 3.5e-3 per SNP in ~100-SNP blocks gives switch_attenuation(3.5e-3, 100) ~ 0.72 (with the 4 % allele
# error the clone attenuation a ~ 0.66, chr21 2:1 at f ~ 0.61); depth meanlog log(4000).
lazy_cached("sim_tum3", function() simulate_scn(
  n_cells = 200, sex = "male", cell_types = c(normal = 0.3, cloneA = 0.35, cloneB = 0.35),
  depth = c(meanlog = log(4000), sdlog = 0.35), events = sim_tum2_events(),
  multiplet_rate = 0.05, phase_block_mean = 2e5, block_flip_prob = 0.5, switch_rate = 3.5e-3, seed = 102))
# sim_tum4: sim_tum2's clones plus two subclones that differ from their parent clone only on chr21 (the loss-inside-CN-LOH
# and LOH-gain cases of the disagreement split): cloneB now carries chr21 CN-LOH (2:0) and cloneC is cloneB
# with a chr21 loss of the same haplotype (1:0: one copy apart, same allele configuration); cloneD is cloneA with chr21
# 3:0 instead of 2:1 (an LOH gain: same total copy number, different major fraction). ~100 cells per tumour clone.
# Three event-free chromosomes (chr3-chr5) are added to sim_genome(): with only chr1 / chr2 / chr8 / chr21 / chr22, every
# autosome carries an event in some clone and a tumour cell's chromosome-class sd_null (MAD over its null chromosomes)
# is several times its depth noise, so a one-copy depth difference on chr21 carries no per-cell evidence.
sim_tum4_events = function() {
  ev = sim_tum2_events()
  b = ev[type == "cloneB"]
  rbind(ev, data.table(chrom = "chr21", start = 0, end = NA, cnA = 2L, cnB = 0L, type = "cloneB"),
        copy(b)[, type := "cloneC"], data.table(chrom = "chr21", start = 0, end = NA, cnA = 1L, cnB = 0L, type = "cloneC"),
        copy(ev[type == "cloneA" & chrom != "chr21"])[, type := "cloneD"],
        data.table(chrom = "chr21", start = 0, end = NA, cnA = 3L, cnB = 0L, type = "cloneD"))
}
lazy_cached("sim_tum4", function() simulate_scn(
  n_cells = 460, sex = "male",
  genome = sim_genome(c(chr1 = 60e6, chr2 = 40e6, chr3 = 40e6, chr4 = 40e6, chr5 = 40e6, chr8 = 30e6, chr21 = 20e6, chr22 = 20e6, chrX = 40e6)), cell_types = c(normal = 0.13, cloneA = 0.2175, cloneB = 0.2175, cloneC = 0.2175, cloneD = 0.2175),
  depth = c(meanlog = log(8000), sdlog = 0.35), events = sim_tum4_events(),
  multiplet_rate = 0.05, phase_block_mean = 2e5, block_flip_prob = 0.5, seed = 104))

# ---- fake tables ---------------------------------------------------------------------------------

#' Regions table (check_regions schema) on a simulated genome: arms (p = [0, cen), q = [cen, end) with
#' cen = 0.4 x length as in sim_genome), whole chromosomes and tiles of the given sizes.
fake_regions = function(bins = sim_genome(), tiles = c("5Mb" = 5e6, "10Mb" = 1e7), sex = "female") {
  cs = bins[, .(len = max(end)), by = chrom]
  arms = rbind(cs[, .(region = paste0(chrom, "p"), class = "arm", chrom, start = 0, end = 0.4 * len)],
               cs[, .(region = paste0(chrom, "q"), class = "arm", chrom, start = 0.4 * len, end = len)])
  chr = cs[, .(region = chrom, class = "chrom", chrom, start = 0, end = len)]
  tl = rbindlist(lapply(names(tiles), function(lab) cs[, {
    s = seq(0, len - 1, by = tiles[[lab]])
    .(region = tile_name(lab, chrom, s), class = lab, start = s, end = pmin(s + tiles[[lab]], len))
  }, by = chrom]))
  r = rbind(arms, chr, tl, use.names = TRUE)
  r[, n_bins := vapply(seq_len(.N), function(i) sum(bins$chrom == chrom[i] & bins$start >= start[i] & bins$end <= end[i]), integer(1))]
  r[, `:=`(sd_class = class, n_use = n_bins, eligible = n_bins > 0, null_region = TRUE,
           allelic_ok = !(sex == "male" & chrom %in% c("chrX", "chrY")))]
  setcolorder(r, c("region", "class", "sd_class", "chrom", "start", "end", "n_bins", "n_use", "eligible", "null_region", "allelic_ok"))
  r[]
}

#' Features long table (check_features schema) with known truth.
#' Null rows: z ~ t(tdf), nA ~ beta-binomial(n_hap, 0.5, rho). Event rows (`events`: cell (name or index),
#' region, state in loss | gain | cnloh): z shifted by the state's log2 shift / max(sd_null, 0.05)
#' (caller_defaults()$shifts; cnloh 0) and nA at f = 0.98 (loss, cnloh) or 2/3 (gain), haplotype A major.
#' Events change only the listed rows (overlapping regions are not propagated). x, e, log2ratio, z and
#' sd_null satisfy the definitions of the features table (vignette("scphasedcn-outputs")) row by row; sd_null is drawn per cell x sd_class (not re-estimated).
#' Extra column `truth` (neutral | loss | gain | cnloh). `hap_frac`: fraction of a cell's fragments that give a
#' molecule vote (n_hap ~ Poisson(hap_frac x n_frag x n_bins / genome bins)). Strata from n_frag ~ lognormal(log 3e4, 0.6).
fake_features = function(n_cells, regions, events = NULL, tdf = 7, rho = 0.002, seed = 1,
                         hap_frac = 0.03, cell_type = "normal", n_folds = 4L) {
  set.seed(seed)
  cells = data.table(cell = sprintf("cell_%04d", seq_len(n_cells)), ci = seq_len(n_cells),
                     n_frag = round(rlnorm(n_cells, log(3e4), 0.6)), s_cell = runif(n_cells, 0.10, 0.20))
  cells[, `:=`(stratum = as.character(depth_stratum(n_frag)), fold = as.integer((ci - 1L) %% n_folds + 1L))]
  rg = as.data.table(regions)
  tot_bins = rg[class == "chrom", sum(n_bins)]
  if (!is.finite(tot_bins) || tot_bins == 0) tot_bins = sum(rg$n_bins)
  sdc = rg[, .(w = median(end - start)), by = sd_class]
  d = CJ(ci = seq_len(n_cells), ri = seq_len(nrow(rg)))
  d = cbind(cells[d$ci], rg[d$ri, .(region, class, sd_class, chrom, start, end, n_bins, eligible, null_region)])
  d[, sd_null := s_cell * sqrt(2e7 / sdc$w[match(sd_class, sdc$sd_class)])]
  d[, truth := "neutral"]
  if (!is.null(events) && nrow(events) > 0) {
    ev = as.data.table(events)
    if (is.numeric(ev$cell)) ev[, cell := sprintf("cell_%04d", as.integer(cell))]
    d[ev, truth := i.state, on = .(cell, region)]
  }
  sh = caller_defaults()$shifts
  delta = fifelse(d$truth == "loss", sh[["loss"]], fifelse(d$truth == "gain", sh[["gain"]], 0))
  d[, z := rt(.N, tdf) + delta / pmax(sd_null, 0.05)]
  d[, log2ratio := z * pmax(sd_null, 0.05)]
  d[, e := n_frag * n_bins / tot_bins]
  d[, x := (e + 0.5) * 2^log2ratio - 0.5]
  d[, n_hap := as.numeric(rpois(.N, hap_frac * n_frag * n_bins / tot_bins))]
  ftrue = fifelse(d$truth %in% c("loss", "cnloh"), 0.98, fifelse(d$truth == "gain", 2 / 3, 0.5))
  p = rbeta(nrow(d), ftrue * (1 - rho) / rho, (1 - ftrue) * (1 - rho) / rho)
  d[, nA := as.numeric(rbinom(.N, n_hap, p))]
  d[, nB := n_hap - nA]
  d[, f := fifelse(n_hap > 0, nA / n_hap, NA_real_)]
  d[, n_blocks := fifelse(n_hap > 0, 1 + as.numeric(rpois(.N, 1)), 0)]
  d[, `:=`(orient = "anchored", cell_type = cell_type, null_cell = !cell %in% d[truth != "neutral", unique(cell)])]
  d[, c("ci", "s_cell", "start", "end", "n_bins") := NULL]
  setcolorder(d, c("cell", "region", "class", "sd_class", "chrom", "x", "e", "log2ratio", "sd_null", "z", "nA", "nB",
                   "n_hap", "n_blocks", "f", "orient", "stratum", "n_frag", "cell_type", "fold", "null_cell",
                   "null_region", "eligible", "truth"))
  d[]
}

#' Allele atoms (check_atoms schema): one atom per bin of a simulated genome.
fake_atoms = function(bins = sim_genome()) {
  data.table(atom = seq_len(nrow(bins)), chrom = bins$chrom, start = bins$start, end = bins$end, bin = seq_len(nrow(bins)))
}

#' Molecule votes (check_votes schema) with clone structure and random per-block orientation.
#' @param clones clone label per cell (length n_cells), or named proportions to sample from
#' @param events data.table clone, chrom, start, end (NA = chromosome end), cnA, cnB; baseline 1:1
#' @param flip_prob probability that a phase block is reported flipped (votes nA / nB swapped)
#' @param block_atoms consecutive atoms per phase block (ps = first atom start + 1)
#' @param mol_per_atom mean molecules per cell and atom at two copies
#' @return votes data.table (cell, atom, ps, nA, nB; rows with nA + nB > 0) with attributes `clone`
#'   (clone per cell), `blocks` (chrom, block, ps, flip) and `truth` (per cell x atom: cnA, cnB)
fake_votes = function(n_cells, atoms, clones, events, flip_prob = 0.5, seed = 1, block_atoms = 10L, mol_per_atom = 0.5) {
  set.seed(seed)
  if (is.numeric(clones) && !is.null(names(clones))) clones = sample(names(clones), n_cells, replace = TRUE, prob = clones)
  if (length(clones) != n_cells) stop("fake_votes: `clones` needs one label per cell")
  at = copy(as.data.table(atoms))
  setorder(at, chrom, start)
  at[, block := (seq_len(.N) - 1L) %/% block_atoms + 1L, by = chrom]
  at[, ps := start[1] + 1, by = .(chrom, block)]
  blocks = unique(at[, .(chrom, block, ps)])
  blocks[, flip := runif(.N) < flip_prob]
  na = nrow(at)
  cnA = matrix(1L, na, n_cells); cnB = matrix(1L, na, n_cells)
  if (!is.null(events) && nrow(events) > 0) for (k in seq_len(nrow(events))) {
    ev = events[k]
    e = if (is.na(ev$end)) Inf else ev$end
    rows = which(at$chrom == ev$chrom & at$start >= ev$start & at$end <= e)
    who = which(clones == ev$clone)
    cnA[rows, who] = as.integer(ev$cnA); cnB[rows, who] = as.integer(ev$cnB)
  }
  cn = cnA + cnB
  nmol = matrix(rpois(length(cn), mol_per_atom * cn / 2), na, n_cells)
  pA = pmin(pmax(cnA / pmax(cn, 1), 0.02), 0.98)
  nA = matrix(rbinom(length(nmol), nmol, pA), na, n_cells)
  nB = nmol - nA
  fl = blocks$flip[match(paste(at$chrom, at$block), paste(blocks$chrom, blocks$block))]
  v = data.table(cell = rep(seq_len(n_cells), each = na), atom = rep(at$atom, n_cells), ps = rep(at$ps, n_cells),
                 nA = as.vector(nA), nB = as.vector(nB), flip = rep(fl, n_cells))
  v[flip == TRUE, c("nA", "nB") := .(nB, nA)]
  v = v[nA + nB > 0, .(cell, atom, ps, nA = as.integer(nA), nB = as.integer(nB))]
  setorder(v, cell, atom)
  setattr(v, "clone", clones)
  setattr(v, "blocks", blocks)
  setattr(v, "truth", list(cnA = cnA, cnB = cnB, atom = at$atom))
  v[]
}
