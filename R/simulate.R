# simulate.R
# What: generative simulator of SPLONGGET single-nucleus DNA data in the same
#       structures the Python extraction layer produces (bins x cells count
#       layers, haplotype layers, per-SNP AD/DP, fragment table).
# Why: the caller tests need a ground truth: the simulator provides cells with
#      depth, per-bin accessibility (with within-bin peak clustering),
#      cell-type chromatin waves, cell-private autocorrelated noise,
#      allele-specific CN events, multiplets as unions of nuclei, phased SNPs
#      with switch errors and ONT allele errors, and fragments with lengths.
# Inputs: parameters (see simulate_scn()). Outputs: list(scn, truth, ...).

#' Mini genome for simulations
#'
#' @param chroms named numeric vector of chromosome lengths (bp)
#' @param bin_size bin width in bp
#' @param centromere_frac position of the centromere as a fraction of length
#' @param seed RNG seed for the bin covariates
#' @return data.table of bins: bin_id, chrom, start, end, gc, frac_blacklist,
#'   frac_trf, is_centromere, peak_bp, arm, n_peaks
#' @export
sim_genome = function(chroms = c(chr1 = 60e6, chr2 = 40e6, chr8 = 30e6, chr21 = 20e6, chr22 = 20e6, chrX = 40e6),
                      bin_size = 1e5, centromere_frac = 0.4, seed = 1) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  set.seed(seed)
  bins = make_bins(chroms, bin_size)
  nb = nrow(bins)
  # GC: smooth random walk around 0.41, clipped to a plausible range
  gc = 0.41 + smooth_noise(nb, sd = 0.05, window = 15)
  bins[, gc := pmin(pmax(gc, 0.3), 0.6)]
  bins[, frac_blacklist := 0]
  bl = sample(nb, max(1, round(0.01 * nb)))
  bins[bl, frac_blacklist := runif(length(bl), 0.3, 1)]
  bins[, frac_trf := pmin(abs(rnorm(nb, 0.01, 0.02)), 0.3)]
  bins[, is_centromere := FALSE]
  for (ch in names(chroms)) {
    cen = centromere_frac * chroms[[ch]]
    bins[chrom == ch & start < cen + 1.5e6 & end > cen - 1.5e6, is_centromere := TRUE]
  }
  bins[, arm := ifelse(start < centromere_frac * chroms[chrom], "p", "q")]
  bins[, n_peaks := rpois(nb, 3)]
  bins[, peak_bp := n_peaks * 500L]
  bins[is_centromere == TRUE, `:=`(n_peaks = 0L, peak_bp = 0L)]
  bins[]
}

smooth_noise = function(n, sd = 1, window = 10) {
  x = rnorm(n + window, 0, sd)
  if (window <= 1) return(x[seq_len(n)])
  k = rep(1 / window, window)
  y = stats::filter(x, k, sides = 2)
  y = y[seq_len(n)]
  y[is.na(y)] = 0
  as.numeric(y) * sqrt(window)
}

#' Default clonal tumour events for the simulator
#'
#' 8q-like amplicon (cnA 3 + cnB 1 = 4 copies on the distal half of chr8),
#' trisomy 21 (2 + 1), monosomy 22 (1 + 0). Carried by cell type "tumour".
#' @return data.table chrom, start, end, cnA, cnB, type
#' @export
sim_events_tumour = function() {
  data.table(chrom = c("chr8", "chr21", "chr22"),
             start = c(15e6, 0, 0), end = c(NA, NA, NA),
             cnA = c(3L, 2L, 1L), cnB = c(1L, 1L, 0L), type = "tumour")
}

#' Simulate a SPLONGGET-like single-nucleus DNA data set
#'
#' @param n_cells number of barcodes (cells) in the output
#' @param genome bins table from [sim_genome()]
#' @param sex "male" (chrX and chrY one copy, hap A) or "female" (no chrY: a chrY in the genome gets 0 copies)
#' @param cell_types named numeric vector of type proportions
#' @param depth lognormal parameters (meanlog, sdlog) of fragments per nucleus
#' @param events data.table (chrom, start, end, cnA, cnB, type); `type` may be a
#'   cell type name or "private" (then `private_frac` of the nuclei carry it)
#' @param private_frac fraction of nuclei carrying a "private" event
#' @param multiplet_rate fraction of barcodes that are multi-nucleus pools
#' @param multiplet_sizes sizes (number of nuclei) and their probabilities
#' @param peak_frac fraction of fragments placed at accessibility peaks (0 =
#'   uniform placement within a bin)
#' @param peak_sd spread (bp) of fragment starts around a peak
#' @param access_sd log-sd of the shared per-bin accessibility
#' @param wave_sd log-sd of the cell-type chromatin wave (2-5 Mb scale)
#' @param noise_phi cell-private multiplicative noise (log-sd), autocorrelated over ~5 bins
#' @param frag_len lognormal parameters of fragment length (bp)
#' @param snp_density het SNPs per bp
#' @param phase_block_mean mean phase-block length (bp)
#' @param switch_rate per-SNP switch-error probability inside a block (1e-4
#'   keeps a 750-SNP block mostly in one parity, attenuation ~0.93 by
#'   [switch_attenuation()]; 0.02 scrambles it)
#' @param block_flip_prob probability that a block's orientation is flipped
#'   relative to the true haplotypes (0 = anchored blocks)
#' @param allele_error per-observation allele error rate (ONT)
#' @param seed RNG seed
#' @return list with `scn` (an `scn` object, see [new_scn()]), `truth` (cells,
#'   nuclei, cn/cnA per nucleus, events), `fragments` (data.table in the Python
#'   `fragments.tsv.gz` layout, `ps` = phase set of the first covered SNP or -1, plus nucleus,
#'   copy, hap_true, hap, n_snps), `snps` (with hap_of_alt_true), `params`
#' @export
simulate_scn = function(n_cells = 200, genome = sim_genome(), sex = c("male", "female"),
                        cell_types = c(normal = 0.6, tumour = 0.4),
                        depth = c(meanlog = log(3000), sdlog = 0.35),
                        events = sim_events_tumour(), private_frac = 0.05,
                        multiplet_rate = 0.15, multiplet_sizes = c(`2` = 0.8, `3` = 0.2),
                        peak_frac = 0.5, peak_sd = 200, access_sd = 0.5, wave_sd = 0.15,
                        noise_phi = 0.25, frag_len = c(meanlog = log(700), sdlog = 0.6),
                        snp_density = 1 / 2000, phase_block_mean = 1.5e6, switch_rate = 1e-4,
                        block_flip_prob = 0, allele_error = 0.04, seed = 1) {
  rng = rng_save(); on.exit(rng_restore(rng), add = TRUE)
  sex = match.arg(sex)
  set.seed(seed)
  bins = copy(genome)
  nb = nrow(bins)
  chrom_sizes = chrom_sizes_from_bins(bins)
  bins[, bin_id := seq_len(nb)]

  # ---- barcodes and nuclei ---------------------------------------------------
  is_mult = runif(n_cells) < multiplet_rate
  k = ifelse(is_mult, as.integer(names(multiplet_sizes))[
    sample.int(length(multiplet_sizes), n_cells, replace = TRUE, prob = multiplet_sizes)], 1L)
  n_nuc = sum(k)
  nuc_cell = rep(seq_len(n_cells), k)
  types = names(cell_types)
  nuc_type = sample(types, n_nuc, replace = TRUE, prob = cell_types)
  nuc_depth = rlnorm(n_nuc, depth[["meanlog"]], depth[["sdlog"]])
  nuclei = data.table(nucleus = seq_len(n_nuc), cell = nuc_cell, type = nuc_type, depth = nuc_depth)

  # ---- truth copy numbers per nucleus ----------------------------------------
  cnA = matrix(1L, nb, n_nuc); cnB = matrix(1L, nb, n_nuc)
  if (sex == "male") cnB[bins$chrom %in% c("chrX", "chrY"), ] = 0L       # one chrX and (when the genome has it) one chrY
  if (sex == "female") { cnA[bins$chrom == "chrY", ] = 0L; cnB[bins$chrom == "chrY", ] = 0L }
  carriers = vector("list", 0)
  if (!is.null(events) && nrow(events) > 0) {
    for (e in seq_len(nrow(events))) {
      ev = events[e]
      rows = region_bins(bins, data.table(chrom = ev$chrom, start = ev$start, end = ev$end))[[1]]
      if (ev$type == "private") {
        who = which(runif(n_nuc) < private_frac)
      } else who = which(nuc_type == ev$type)
      cnA[rows, who] = as.integer(ev$cnA); cnB[rows, who] = as.integer(ev$cnB)
      carriers[[e]] = who
    }
  }
  cn = cnA + cnB

  # ---- accessibility, waves, noise -------------------------------------------
  access = exp(smooth_noise(nb, sd = access_sd, window = 3) - 8 * (bins$gc - 0.41)^2 * 10 +
                 0.2 * log1p(bins$peak_bp / 500))
  access[bins$is_centromere] = access[bins$is_centromere] * 0.05
  access[bins$frac_blacklist > 0.5] = access[bins$frac_blacklist > 0.5] * 5
  wave = sapply(types, function(t) exp(smooth_noise(nb, sd = wave_sd, window = 30)))
  colnames(wave) = types
  eps = matrix(0, nb, n_nuc)
  for (i in seq_len(n_nuc)) eps[, i] = exp(smooth_noise(nb, sd = noise_phi, window = 5))
  w = access * wave[, nuc_type, drop = FALSE] * eps * cn / 2
  w[cn == 0] = 0
  size = nuc_depth / colSums(w)
  mu = sweep(w, 2, size, "*")
  nfr = matrix(rpois(length(mu), mu), nb, n_nuc)

  # ---- fragments ----------------------------------------------------------------
  tot = sum(nfr)
  nz = which(nfr > 0, arr.ind = TRUE)
  fr = data.table(bin = rep(nz[, 1], nfr[nz]), nucleus = rep(nz[, 2], nfr[nz]))
  fr[, copy := 1L + floor(runif(.N) * cn[cbind(bin, nucleus)])]
  fr[, hap_true := ifelse(copy <= cnA[cbind(bin, nucleus)], "A", "B")]
  bsz = bins$end - bins$start
  fr[, chrom := bins$chrom[bin]]
  at_peak = runif(tot) < peak_frac & bins$n_peaks[fr$bin] > 0
  peak_id = floor(runif(tot) * pmax(bins$n_peaks[fr$bin], 1))
  # deterministic peak centres per bin: hash of (bin, peak_id)
  centre = bins$start[fr$bin] + ((peak_id * 7919 + fr$bin * 104729) %% bsz[fr$bin])
  start = ifelse(at_peak, centre + rnorm(tot, 0, peak_sd), bins$start[fr$bin] + runif(tot) * bsz[fr$bin])
  start = pmin(pmax(floor(start), bins$start[fr$bin]), bins$end[fr$bin] - 1)
  len = pmin(pmax(round(rlnorm(tot, frag_len[["meanlog"]], frag_len[["sdlog"]])), 100), 50000)
  fr[, start := as.numeric(start)]
  fr[, end := pmin(start + len, chrom_sizes[chrom])]
  fr[, strand := sample(c("+", "-"), .N, replace = TRUE)]
  fr[, mapq := 60L]
  # one molecule cannot yield overlapping fragments: drop overlaps within (nucleus, chrom, copy)
  setorder(fr, nucleus, chrom, copy, start)
  fr[, tmp := .GRP, by = .(nucleus, chrom, copy)]
  keep = nonoverlap_keep(fr$tmp, fr$start, fr$end)
  fr = fr[keep]
  fr[, tmp := NULL]
  fr[, cell := nuclei$cell[nucleus]]

  # ---- SNPs, phase blocks, switch errors ----------------------------------------
  snps = rbindlist(lapply(names(chrom_sizes), function(ch) {
    n = rpois(1, chrom_sizes[[ch]] * snp_density)
    data.table(chrom = ch, pos = sort(floor(runif(n) * chrom_sizes[[ch]])))
  }))
  snps[, hap_of_alt_true := sample(c("A", "B"), .N, replace = TRUE)]
  snps[, block := cumsum(c(TRUE, runif(.N - 1) < 1 / (phase_block_mean * snp_density))), by = chrom]
  snps[, ps := pos[1] + 1, by = .(chrom, block)]
  snps[, switch := cumsum(runif(.N) < switch_rate) %% 2 == 1, by = .(chrom, block)]
  flipped_blocks = snps[, .(flip = runif(1) < block_flip_prob), by = .(chrom, block)]
  snps = merge(snps, flipped_blocks, by = c("chrom", "block"), sort = FALSE)
  snps[, orient := xor(switch, flip)]
  snps[, hap_of_alt := ifelse(orient, ifelse(hap_of_alt_true == "A", "B", "A"), hap_of_alt_true)]
  setorder(snps, chrom, pos)
  snps[, snp := seq_len(.N)]
  snps[, id := paste0(chrom, ":", pos)]
  # chrX SNPs are kept in a male: hemizygous, every read carries hap A (a real het list would drop them)

  # ---- fragment x SNP observations ---------------------------------------------
  obs = fragment_snp_obs(fr, snps, allele_error)
  fr[, n_snps := 0L]
  fr[, hap := "U"]
  # phase set of the fragment: `ps` of the first SNP it covers, -1 when it covers none (the Python
  # extraction's sentinel, alleles.py `vote()`); deterministic, adds no RNG draws
  fr[, ps := -1L]
  if (nrow(obs) > 0) {
    votes = obs[, .(n_snps = .N,
                    vA = sum(alt_obs == (hap_of_alt == "A")),
                    vB = sum(alt_obs == (hap_of_alt == "B")),
                    snp1 = min(snp)), by = idx]
    fr[votes$idx, n_snps := votes$n_snps]
    fr[votes$idx, hap := ifelse(votes$vA > votes$vB, "A", ifelse(votes$vB > votes$vA, "B", "N"))]
    fr[votes$idx, ps := as.integer(snps$ps[votes$snp1])]
  }
  fr[, flags := 2L]                        # F_MAPQ60 (all simulated fragments are MAPQ 60)

  # ---- barcode-level matrices --------------------------------------------------
  cells = data.table(cell = sprintf("cell_%04d", seq_len(n_cells)))
  built = build_layers_from_fragments(fr, bins, n_cells = n_cells, cell_index = fr$cell)
  AD_DP = alleles_from_obs(obs, fr, snps, n_cells)

  # obs columns and truth
  dom = nuclei[, .(cell_type = type[1], n_nuclei = .N, members = paste(nucleus, collapse = ";"),
                   types = paste(type, collapse = ";"), depth_true = sum(depth)), by = cell]
  setorder(dom, cell)
  obs_dt = data.table(cell = cells$cell,
                      n_frag = as.integer(colSums(built$counts)),
                      span_mean = colSums(built$layers$span) / pmax(colSums(built$counts), 1),
                      frac_chrX = colSums(built$counts[bins$chrom == "chrX", , drop = FALSE]) / pmax(colSums(built$counts), 1),
                      cell_type = dom$cell_type)
  truth_cells = data.table(cell = cells$cell, dom[, .(n_nuclei, multiplet = n_nuclei > 1, members, types, depth_true)])

  # true CN per barcode = member nuclei summed (what the barcode actually contains)
  nuc_by_cell = split(nuclei$nucleus, nuclei$cell)
  cn_cell = sapply(seq_len(n_cells), function(i) rowSums(cn[, nuc_by_cell[[i]], drop = FALSE]))
  cnA_cell = sapply(seq_len(n_cells), function(i) rowSums(cnA[, nuc_by_cell[[i]], drop = FALSE]))
  colnames(cn_cell) = colnames(cnA_cell) = cells$cell

  scn = new_scn(counts = built$counts, bins = bins, cells = obs_dt, layers = built$layers,
                obsm = built$obsm,
                uns = list(params = list(simulated = TRUE, sex = sex, bin_size = bins$end[1] - bins$start[1])),
                alleles = list(AD = AD_DP$AD, DP = AD_DP$DP,
                               snps = snps[, .(snp = id, chrom, pos, ref = "A", alt = "G", hap_of_alt, ps)]))
  scn$cells[, sex := sex]
  colnames(scn$counts) = cells$cell
  for (nm in names(scn$layers)) colnames(scn$layers[[nm]]) = cells$cell
  colnames(scn$alleles$AD) = colnames(scn$alleles$DP) = cells$cell
  fr_out = fr[, .(chrom, start, end, cell, strand, mapq, flags, len = end - start, piece = 1L, npieces = 1L,
                  hap = hap_code(hap), nsnp = n_snps, ps, read_id = seq_len(.N), nucleus, copy, hap_true, bin)]
  setorder(fr_out, chrom, start, cell)
  list(scn = scn,
       truth = list(cells = truth_cells, nuclei = nuclei, cn = cn_cell, cnA = cnA_cell,
                    cn_nucleus = cn, cnA_nucleus = cnA, events = events, carriers = carriers,
                    access = access, wave = wave, size = size),
       fragments = fr_out, snps = snps,
       params = list(n_cells = n_cells, sex = sex, cell_types = cell_types, depth = depth,
                     multiplet_rate = multiplet_rate, peak_frac = peak_frac, frag_len = frag_len,
                     snp_density = snp_density, switch_rate = switch_rate, allele_error = allele_error,
                     block_flip_prob = block_flip_prob, seed = seed))
}

# Drop fragments that overlap an earlier-accepted fragment of the same molecule (same group, sorted by
# group then start): one Tn5-cut molecule yields non-overlapping fragments. Returns a logical keep vector.
# Exact semantics of a sequential scan (keep k when start[k] >= end of the last kept fragment of its group;
# kept fragments are disjoint and sorted, so that end is the running maximum over the kept ones). Fragments
# that start at or after the end of every earlier fragment of the group are always kept and reset the scan;
# only the runs between them are scanned in R.
nonoverlap_keep = function(group, start, end) {
  n = length(start)
  keep = rep(TRUE, n)
  if (n < 2) return(keep)
  first = c(TRUE, group[-1] != group[-n])
  prev_max = c(-Inf, stats::ave(end, cumsum(first), FUN = cummax)[-n])   # max end over earlier fragments
  prev_max[first] = -Inf
  free = start >= prev_max
  run = cumsum(free)
  for (idx in split(which(run %in% run[!free]), run[run %in% run[!free]])) {
    cur = end[idx[1]]
    for (k in idx[-1]) if (start[k] < cur) keep[k] = FALSE else cur = end[k]
  }
  keep
}

# Fragment x SNP observations: one row per (fragment, covered SNP) with the
# observed allele (alt_obs) after ONT error; idx = row of the fragment table.
fragment_snp_obs = function(fr, snps, allele_error) {
  out = vector("list", 0)
  for (ch in unique(fr$chrom)) {
    sp = snps[chrom == ch]
    if (nrow(sp) == 0) next
    sel = which(fr$chrom == ch)
    lo = findInterval(fr$start[sel] - 1, sp$pos)      # SNPs with pos >= start
    hi = findInterval(fr$end[sel] - 1, sp$pos)        # SNPs with pos < end
    nn = hi - lo
    has = which(nn > 0)
    if (length(has) == 0) next
    idx = rep(sel[has], nn[has])
    snp_i = unlist(lapply(has, function(h) seq.int(lo[h] + 1, hi[h])))
    o = data.table(idx = idx, snp = sp$snp[snp_i], hap_of_alt_true = sp$hap_of_alt_true[snp_i],
                   hap_of_alt = sp$hap_of_alt[snp_i])
    out[[ch]] = o
  }
  if (length(out) == 0) return(data.table(idx = integer(), snp = integer(), hap_of_alt_true = character(),
                                          hap_of_alt = character(), alt_obs = logical()))
  o = rbindlist(out)
  frag_hap = fr$hap_true[o$idx]
  truth_alt = frag_hap == o$hap_of_alt_true
  err = runif(nrow(o)) < allele_error
  o[, alt_obs := xor(truth_alt, err)]
  o
}

alleles_from_obs = function(obs, fr, snps, n_cells) {
  ns = nrow(snps)
  if (nrow(obs) == 0) {
    z = sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(ns, n_cells))
    return(list(AD = z, DP = z))
  }
  cell_i = as.integer(fr$cell[obs$idx])
  agg = data.table(snp = obs$snp, cell = cell_i, alt = as.integer(obs$alt_obs))[, .(DP = .N, AD = sum(alt)), by = .(snp, cell)]
  DP = sparseMatrix(i = agg$snp, j = agg$cell, x = as.numeric(agg$DP), dims = c(ns, n_cells))
  AD = sparseMatrix(i = agg$snp, j = agg$cell, x = as.numeric(agg$AD), dims = c(ns, n_cells))
  list(AD = AD, DP = DP)
}

#' Build the bin-level layers from a fragment table
#'
#' Counts fragments by the bin of their start, spans, fragment-length classes,
#' haplotype layers from the fragment `hap` column, and the obsm histograms
#' (fragment length, MAPQ) in the Python extraction layout. Used by the
#' simulator.
#' @param fr data.table with chrom, start, end, cell (integer index or factor),
#'   hap (optional; "U"/"A"/"B"/"N" or the integer codes), mapq (optional),
#'   bin (optional; else located from start)
#' @param bins bins table
#' @param n_cells number of cells (columns)
#' @param cell_index integer cell index per fragment (1..n_cells)
#' @param fraglen_breaks fragment-length histogram breaks (Python layout)
#' @return list(counts, layers, obsm)
#' @export
build_layers_from_fragments = function(fr, bins, n_cells, cell_index = fr$cell,
                                       fraglen_breaks = c(0, 200, 300, 500, 750, 1000, 2000, 5000, 10000, Inf)) {
  nb = nrow(bins)
  ci = as.integer(cell_index)
  b = if ("bin" %in% names(fr)) fr$bin else locate_bins(bins, fr$chrom, fr$start)
  ok = !is.na(b)
  len = fr$end - fr$start
  sm = function(i, j, x) {
    sel = ok & !is.na(i)
    sparseMatrix(i = i[sel], j = j[sel], x = as.numeric(x[sel]), dims = c(nb, n_cells))
  }
  counts = sm(b, ci, rep(1, nrow(fr)))
  layers = list(span = sm(b, ci, len))
  hap = if ("hap" %in% names(fr)) hap_label(fr$hap) else rep("U", nrow(fr))
  layers$frag_short = sm(b, ci, as.numeric(len < 300))
  layers$frag_long = sm(b, ci, as.numeric(len >= 1000))
  layers$hapA = sm(b, ci, as.numeric(hap == "A"))
  layers$hapB = sm(b, ci, as.numeric(hap == "B"))
  layers$hap_ambig = sm(b, ci, as.numeric(hap == "N"))

  # obsm histograms (cells x k), same layouts as the Python writer
  hist_by_cell = function(code, cells, k, labels) {
    m = matrix(0, n_cells, k)
    if (length(code) > 0) {
      ok = !is.na(code) & code >= 1 & code <= k
      t = table(factor(cells[ok], levels = seq_len(n_cells)), factor(code[ok], levels = seq_len(k)))
      m[] = as.numeric(t)
    }
    colnames(m) = labels
    m
  }
  mapq = if ("mapq" %in% names(fr)) fr$mapq else rep(60L, nrow(fr))
  mapq_code = ifelse(mapq == 0, 1L, ifelse(mapq < 30, 2L, ifelse(mapq < 60, 3L, 4L)))
  obsm = list(
    fraglen_hist = hist_by_cell(as.integer(cut(len[ok], fraglen_breaks, right = FALSE)), ci[ok], length(fraglen_breaks) - 1,
                                c("lt200", "200_300", "300_500", "500_750", "750_1000", "1k_2k", "2k_5k", "5k_10k", "ge10k")),
    mapq_hist = hist_by_cell(mapq_code[ok], ci[ok], 4, c("mapq0", "mapq1_29", "mapq30_59", "mapq60"))
  )
  list(counts = counts, layers = layers, obsm = obsm)
}
