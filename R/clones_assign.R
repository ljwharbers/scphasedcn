# clones_assign.R
# What: tumour mode of the copy-number caller, part 3 of 4: pieces (the union of all clone
#       segment breakpoints) and per-cell features on arbitrary regions, oriented counts, clone state /
#       allelic distances and trunks, likelihood re-assignment of cells (with 50:50 mixtures as doublet
#       hypotheses), doublet clones, allelic cleaning of the normal reference pool, automatic control regions
#       and the all-neutral normal clone.
# Why: cells are scored under the clone states (per-cell states are unreliable in tumours). In a tumour
#      sample the pool of copy-number-flat cells can carry tumour haplotypes, so the normal reference is
#      cleaned with the allele channel.
# Inputs: segments and profiles (clones_states.R), per-cell features on the pieces, votes and atoms.
# Outputs: data.tables (pieces, assignments, doublet flags, reference pool, control regions).
# See also: clones_discovery.R, clones_states.R, clones_tumour.R.

# ---- pieces and per-cell features on arbitrary regions ---------------------------------------------

#' Pieces: the union of all clone segment breakpoints
#'
#' Per chromosome, every segment start / end of every clone is a breakpoint; the pieces between adjacent
#' breakpoints are regions of class `segment` named `seg|<tag>|chrom:start-end`.
#' @param segments segments table
#' @param regions regions table of the run (for `sd_class`, n_bins, allelic_ok)
#' @param tag clone field of the region name ("all" for the shared pieces)
#' @return regions table ([check_regions()]) of the pieces
#' @export
segment_pieces = function(segments, regions = NULL, tag = "all") {
  sg = as.data.table(segments)
  pc = sg[, {
    b = sort(unique(c(start, end)))
    list(start = b[-length(b)], end = b[-1])
  }, by = chrom]
  # keep pieces covered by at least one segment
  cov = sg[pc, on = .(chrom, start <= start, end >= end), .N, by = .EACHI]
  pc = pc[cov$N > 0]
  seg_regions(pc, regions, tag)
}

# regions table for segments / pieces (data.table chrom, start, end [, clone])
seg_regions = function(pc, regions = NULL, tag = "all") {
  pc = as.data.table(pc)
  tg = if ("clone" %in% names(pc)) pc$clone else rep(tag, nrow(pc))
  out = data.table(region = sprintf("seg|%s|%s:%s-%s", tg, pc$chrom, format(pc$start, scientific = FALSE, trim = TRUE),
                                    format(pc$end, scientific = FALSE, trim = TRUE)),
                   class = "segment", sd_class = NA_character_, chrom = pc$chrom, start = as.numeric(pc$start),
                   end = as.numeric(pc$end), n_bins = NA_integer_, n_use = NA_integer_, eligible = TRUE,
                   null_region = FALSE, allelic_ok = TRUE)
  if (!is.null(regions)) {
    rg = as.data.table(regions)
    out[, sd_class := seg_sd_class(out, rg)]
    ao = rg[class %in% c("chrom", "arm"), .(ok = all(allelic_ok)), by = chrom]
    out[ao, allelic_ok := i.ok, on = "chrom"]
    bl = rg[class == "chrom", .(chrom, len = end - start, nb = n_bins, nu = n_use)]
    out[bl, `:=`(n_bins = as.integer(round(i.nb * (end - start) / i.len)), n_use = as.integer(round(i.nu * (end - start) / i.len))), on = "chrom"]
  }
  out = unique(out, by = "region")
  out[]
}
utils::globalVariables(c("ok", "i.ok", "i.nb", "i.nu", "i.len", "nb", "nu"))

# sd_class of segment regions (the regions-table rule): chrom if >= 90 % of the chromosome, arm if >= 60 % of an
# arm, else the nearest tile size on a log scale (coverage measured in bp here)
seg_sd_class = function(sg, regions) {
  chr = regions[class == "chrom", .(chrom, cs = start, ce = end)]
  arms = regions[class == "arm", .(chrom, as_ = start, ae = end)]
  tiles = intersect(c("5Mb", "10Mb", "20Mb", "50Mb"), unique(regions$class))
  tsz = caller_defaults()$tile_sizes[tiles]
  vapply(seq_len(nrow(sg)), function(i) {
    ch = chr[chrom == sg$chrom[i]]
    len = sg$end[i] - sg$start[i]
    if (nrow(ch) && len >= 0.9 * (ch$ce[1] - ch$cs[1])) return("chrom")
    am = arms[chrom == sg$chrom[i]]
    if (nrow(am)) {
      ov = pmax(0, pmin(am$ae, sg$end[i]) - pmax(am$as_, sg$start[i])) / (am$ae - am$as_)
      if (any(ov >= 0.6)) return("arm")
    }
    if (!length(tsz)) return("arm")
    names(tsz)[which.min(abs(log(len) - log(tsz)))]
  }, character(1))
}
utils::globalVariables(c("cs", "ce", "as_", "ae"))
utils::globalVariables("i.region")

# per-row t parameters with fallback (sd_class x stratum -> class arm x stratum -> first row)
cl_match_tpar = function(d, tpar) {
  if (is.null(tpar)) return(NULL)
  tp = as.data.table(tpar)
  k1 = match(paste(d$sd_class, d$stratum), paste(tp$class, tp$stratum))
  k2 = match(paste("arm", d$stratum), paste(tp$class, tp$stratum))
  k3 = match(d$stratum, tp$stratum)
  k = ifelse(!is.na(k1), k1, ifelse(!is.na(k2), k2, ifelse(!is.na(k3), k3, 1L)))
  data.table(t_m = tp$t_m[k], t_s = tp$t_s[k], t_df = tp$t_df[k])
}

# oriented counts per cell x region under the pseudobulk of a cell set (leave-one-out for its members)
# tv: list(trip = data.table ci, region index, blk, nA, nB; MA, MB sparse cells x blocks; nr)
oriented_counts = function(tv, members) {
  ind = numeric(tv$n_cells); ind[members] = 1
  GA = as.vector(Matrix::crossprod(tv$MA, ind)); GB = as.vector(Matrix::crossprod(tv$MB, ind))
  nk = tv$n_cells * tv$nr
  if (is.null(tv$XA)) tv = trip_index(tv)
  # every row with the block's pseudobulk orientation (sparse products over blocks) ...
  pos = as.numeric(GA > GB); neg = as.numeric(GA < GB)
  sA = as.vector(tv$XA %*% pos + tv$XB %*% neg); sN = as.vector(tv$XN %*% (pos + neg))
  # ... then the members' own rows re-oriented without their own votes (leave-one-out), where that changes the sign
  tr = tv$trip
  m = which(ind[tr$ci] > 0)
  if (length(m)) {
    b = tr$blk[m]; xa = tr$nA[m]; xb = tr$nB[m]
    d0 = GA[b] - GB[b]; d1 = d0 - xa + xb
    ch = which(sign(d0) != sign(d1))
    if (length(ch)) {
      d0 = d0[ch]; d1 = d1[ch]; xa = xa[ch]; xb = xb[ch]
      dA = ((d1 > 0) * xa + (d1 < 0) * xb) - ((d0 > 0) * xa + (d0 < 0) * xb)
      dN = ((d1 != 0) - (d0 != 0)) * (xa + xb)
      key = tv$key[m[ch]]
      sA = sA + bin_sum(key, dA, nk); sN = sN + bin_sum(key, dN, nk)
    }
  }
  k = which(sN > 0 | sA > 0)
  # keys are column-major: cell + (region - 1) * n_cells
  data.table(ci = as.integer((k - 1) %% tv$n_cells + 1), ri = as.integer((k - 1) %/% tv$n_cells + 1), oA = sA[k], oN = sN[k])
}

# sparse (cell x region) x block vote matrices of a trip_votes() table (keys cell + (region - 1) * n_cells)
trip_index = function(tv) {
  tr = tv$trip
  nk = tv$n_cells * tv$nr; nb = ncol(tv$MA)
  tv$key = as.integer(tr$ci + (tr$ri - 1L) * tv$n_cells)
  tv$XA = Matrix::sparseMatrix(i = tv$key, j = tr$blk, x = tr$nA, dims = c(nk, nb))
  tv$XB = Matrix::sparseMatrix(i = tv$key, j = tr$blk, x = tr$nB, dims = c(nk, nb))
  tv$XN = tv$XA + tv$XB
  tv
}
utils::globalVariables(c("oN"))

trip_votes = function(votes, atoms, regs, n_cells, allelic_chroms = NULL) {
  rg = as.data.table(regs)
  if (!is.null(allelic_chroms)) rg = rg[chrom %in% allelic_chroms]
  if ("allelic_ok" %in% names(rg)) rg = rg[allelic_ok == TRUE]
  v = region_block_votes(votes, atoms, rg)
  v = v[cell <= n_cells]
  nb = max(c(1L, v$blk))
  list(trip = v[, .(ci = cell, ri = match(region, regs$region), blk, nA = as.numeric(nA), nB = as.numeric(nB))],
       MA = Matrix::sparseMatrix(i = v$cell, j = v$blk, x = as.numeric(v$nA), dims = c(n_cells, nb)),
       MB = Matrix::sparseMatrix(i = v$cell, j = v$blk, x = as.numeric(v$nB), dims = c(n_cells, nb)),
       nr = nrow(regs), n_cells = n_cells) |> trip_index()
}

# per clone x piece state (from the clone's segment containing the piece midpoint)
piece_states = function(pieces, segments, grid = cn_grid(10)) {
  sg = as.data.table(segments)
  pc = as.data.table(pieces)[, .(region, chrom, mid = (start + end) / 2)]
  out = sg[pc, on = .(chrom, start <= mid, end > mid), .(clone = x.clone, region = i.region, state = x.state, cn = x.cn,
                                                        ploidy = x.ploidy, label = x.label), nomatch = NULL, allow.cartesian = TRUE]
  out[, `:=`(nA_s = as.integer(sub(":.*", "", state)), nB_s = as.integer(sub(".*:", "", state)))]
  out[, f_state := fifelse(cn > 0, nA_s / cn, 0.5)]
  out[]
}
utils::globalVariables(c("mid", "x.clone", "x.state", "x.cn", "x.ploidy", "x.label"))

#' Allelic state distance between clones
#'
#' Base pairs of the shared pieces ([segment_pieces()]) on which two clones' states differ in their major-haplotype
#' fraction by at least `min_df` (e.g. 1:1 vs 2:1, 1:0 vs 2:1) or in total copy number by at least `min_dcn`. Depth-only
#' differences of one copy between amplified states (4:0 vs 3:0) do not count: a cell between two such clones is
#' fitted better by their mixture than by either clone, which is not evidence of a doublet.
#' @param segments segments table
#' @param min_df minimal difference of the state major fractions
#' @param min_dcn minimal total copy-number difference that counts on its own
#' @param min_len pieces shorter than this are ignored
#' @return symmetric numeric matrix clones x clones (bp)
#' @export
clone_state_distance = function(segments, min_df = 0.15, min_dcn = 2, min_len = 5e6) {
  sg = as.data.table(segments)
  K = sort(unique(sg$clone))
  out = matrix(0, length(K), length(K), dimnames = list(K, K))
  if (length(K) < 2) return(out)
  pcs = segment_pieces(sg)
  pcs = pcs[(end - start) >= min_len]
  if (!nrow(pcs)) return(out)
  st = piece_states(pcs, sg)
  st[ploidy == 1, f_state := fifelse(cn > 0, 1, 0.5)]
  st[, len := pcs$end[match(region, pcs$region)] - pcs$start[match(region, pcs$region)]]
  Fm = dcast(st, region ~ clone, value.var = "f_state")
  Cm = dcast(st, region ~ clone, value.var = "cn")
  ln = st[match(Fm$region, region), len]
  Fm = as.matrix(Fm[, K, with = FALSE]); Cm = as.matrix(Cm[, K, with = FALSE])
  for (i in seq_along(K)) for (j in seq_along(K)) if (i < j) {
    d = (abs(Fm[, i] - Fm[, j]) >= min_df - 1e-9) | (abs(Cm[, i] - Cm[, j]) >= min_dcn)
    out[i, j] = out[j, i] = sum(ln[d %in% TRUE])
  }
  out
}

#' Trunks of a clone set
#'
#' Groups clones that share their clonal events: complete-linkage clustering of [clone_state_distance()] cut at
#' `max_dist` bp (two clones of one trunk differ allelically on less than `max_dist` of the genome). Trunks are the
#' units of block orientation for within-trunk clone discovery (subclonal events of 5-25 % of a trunk are oriented
#' far better by the trunk's pseudobulk than by all cells) and of the doublet hypotheses (mixtures of clones within a
#' trunk are not scored).
#' @param segments segments table
#' @param max_dist cut height (bp)
#' @param ... passed to [clone_state_distance()]
#' @return named character vector clone -> trunk id (`T1`, `T2`, ...)
#' @export
clone_trunks = function(segments, max_dist = 5e7, ...) {
  Dm = clone_state_distance(segments, ...)
  K = rownames(Dm)
  if (length(K) < 2) return(stats::setNames(rep("T1", length(K)), K))
  tr = stats::cutree(stats::hclust(stats::as.dist(Dm), method = "complete"), h = max_dist)
  stats::setNames(paste0("T", tr), K)
}

#' Allelic distance between clones under the all-cell orientation
#'
#' Clone pseudobulk haplotype fractions per arm with every phase block oriented by the all-cell pseudobulk (`oc`: per
#' cell x arm counts from `orient_blocks` without groups, computed once): `D[k, l] = sum over arms of len x max(|f_k -
#' f_l| - floor, 0)` (bp x fraction). Events carried by a large share of the cells (trunk events, e.g. an arm loss in
#' one tumour region, a chromosome lost on opposite haplotypes in two lineages, a gain in all tumour cells) are oriented correctly in nearly every block
#' and give large differences; a subclone's private event is outweighed in the all-cell pseudobulk and contributes
#' little. This is what makes the distance a trunk distance, unlike the clone state distance, where sibling clones
#' with one misfitted segment are as far apart as different lineages.
#' @param oc data.table cell (name), region (arm), nA, nB: all-cell oriented votes
#' @param members data.table cell, clone
#' @param arms regions table of the arms (region, start, end)
#' @param floor noise floor of |f_k - f_l|
#' @param min_mol arms with fewer pooled molecules in a clone are ignored for that clone's pairs
#' @return symmetric numeric matrix clones x clones
#' @export
clone_allelic_distance = function(oc, members, arms, floor = 0.03, min_mol = 200) {
  m = as.data.table(members)[, .(cell, clone = as.character(clone))]
  x = merge(as.data.table(oc)[, .(cell, region, nA, nB)], m, by = "cell")
  p = x[, .(nA = sum(nA), n = sum(nA + nB)), by = .(clone, region)]
  p[, f := fifelse(n >= min_mol, nA / n, NA_real_)]
  K = sort(unique(m$clone))
  ar = as.data.table(arms)[region %in% p$region]
  p = p[region %in% ar$region]
  Fm = matrix(NA_real_, length(K), nrow(ar), dimnames = list(K, ar$region))
  Fm[cbind(match(p$clone, K), match(p$region, ar$region))] = p$f
  ln = ar$end - ar$start
  out = matrix(0, length(K), length(K), dimnames = list(K, K))
  for (i in seq_along(K)) for (j in seq_along(K)) if (i < j) {
    d = pmax(abs(Fm[i, ] - Fm[j, ]) - floor, 0)
    out[i, j] = out[j, i] = sum(ln * d, na.rm = TRUE)
  }
  out
}

#' Assign cells to clones by their likelihood under the clone states
#'
#' Per cell and clone k: the sum over pieces of the depth term `ll_depth_shift(z, c_k * log2(cn / 2) /
#' pmax(sd_null, sd_floor), tpar)` and the allele term `ll_allelic(nA, n, 0.5 + a_k * (f_state - 0.5), rho,
#' two_sided = FALSE)` on the cell's votes oriented by clone k's pseudobulk (leave-one-out for k's members).
#' The depth z is taken relative to the reference cells' mean log2 ratio of each piece (`ref_cells`;
#' removes the shared per-region bias of the normalisation). With `mixtures`, every pair of clones is also
#' scored as a mixture (weights `mix_w`; copy numbers mixed; per piece the votes oriented by the pair's clone
#' that is imbalanced there, as for the single clones, with the other clone's haplotype agreement taken from its
#' pooled counts under that orientation).
#' Posterior over singles (prior `prior`, default uniform) and mixtures (total prior `mix_prior`); the
#' reported `posterior` is the best clone's posterior among the single clones.
#' @param feat_pieces list(features = per-cell depth features on the pieces ([region_features()]),
#'   pieces = regions table, votes, atoms, cells = names in vote index order)
#' @param clones data.table cell, clone: the clone members used for orientation (confident singlets)
#' @param segments segments table
#' @param tpar depth null t parameters ([fit_depth_null()]) or NULL (Gaussian)
#' @param rho per-cell beta-binomial overdispersion
#' @param par data.table clone, c, a (from [fit_clone_states()]); default c = a = 1
#' @param prior named clone prior (NULL = uniform; a size-proportional prior lets large diffuse clones absorb low-depth cells)
#' @param mixtures score pairwise 50:50 mixtures
#' @param min_post posterior needed for `assigned`
#' @param mix_prior total prior of the mixture hypotheses
#' @param ref_cells reference cells for the depth centring (NULL = none); the clone states must be relative to the
#'   same reference ([clone_profiles()] `ref_cells`)
#' @param mix_w mixture weights of the first clone of a pair (0.5 = 50:50; e.g. c(1/3, 1/2, 2/3) for unequal nuclei)
#' @param mix_margin a cell is `doublet_like` when the mixtures' posterior is > 0.5 and the best mixture beats the
#'   best single clone by more than this (log-likelihood)
#' @param mix_rel ... and by more than this fraction of the cell's best single-clone log-likelihood ratio (deep
#'   cells fit any single clone less perfectly in absolute terms)
#' @param params [caller_defaults()]
#' @param mix_pairs optional 2-column matrix of the clone pairs scored as mixtures (default all pairs);
#'   [caller_tumour()] passes the pairs whose states differ allelically on >= `mix_min_dist` bp
#'   ([clone_state_distance()]): a mixture of near-identical (sibling) clones fits a cell between them better than
#'   either clone, which is not evidence of a doublet
#' @param mult_prior optional named per-cell multiplet prior probability (e.g. a soft multiplet score): the cell's total
#'   mixture prior is max(`mix_prior`, prior), capped at `max_mix_prior`
#' @param trunks optional named clone -> trunk vector ([clone_trunks()]): adds `trunk` and `p_trunk` (posterior summed
#'   over the clones of the best clone's trunk)
#' @param max_mix_prior cap of the per-cell mixture prior
#' @param same_groups optional named clone -> group vector (e.g. the caller's same-state groups): `posterior` (and so
#'   `assigned` / `confident`) is then summed over the clones of the best clone's group, and the single best clone's
#'   posterior is kept as `post_clone`. Two clones with the same states otherwise share each cell's posterior and
#'   neither reaches `min_post`
#' @return data.table cell, clone (the likelihood-best single clone, for every cell), posterior, ll_margin,
#'   doublet_like, pair, ref_pool, p_mix, ll_mix_gain, ll_best, mix_prior, assigned / confident (posterior >= `min_post`
#'   and not doublet-like) (plus trunk, p_trunk when `trunks` is given; attribute `ll`: cells x hypotheses log-likelihood matrix)
#' @export
assign_cells = function(feat_pieces, clones, segments, tpar, rho, par = NULL, prior = NULL, mixtures = TRUE, min_post = 0.9,
                        mix_prior = 0.05, ref_cells = NULL, mix_w = 0.5, mix_margin = 10, mix_rel = 0.1, params = caller_defaults(),
                        mix_pairs = NULL, mult_prior = NULL, trunks = NULL, max_mix_prior = 0.5, same_groups = NULL) {
  fp = feat_pieces
  pieces = as.data.table(fp$pieces)
  cells = fp$cells
  d = as.data.table(fp$features)[region %in% pieces$region]
  clones = as.data.table(clones)[!is.na(clone), .(cell, clone = as.character(clone))]
  K = sort(unique(clones$clone))
  par = if (is.null(par)) data.table(clone = K, c = 1, a = 1) else as.data.table(par)[, .(clone = as.character(clone), c, a)]
  ps_ = piece_states(pieces, segments)
  ps_ = ps_[clone %in% K]
  # depth
  if (length(ref_cells) >= 10) {
    rm_ = d[cell %in% ref_cells, .(rmean = mean(log2ratio, na.rm = TRUE)), by = region]
    d[rm_, log2ratio := log2ratio - i.rmean, on = "region"]
  }
  d[, z := log2ratio / pmax(sd_null, params$sd_floor)]
  tp = cl_match_tpar(d, tpar)
  sc = unique(d$cell)
  ci_d = match(d$cell, sc)
  depth_ll = function(cnv, plv, cc_) {
    base = if (length(ref_cells) >= 10) plv else 2
    mu = cc_ * log2(pmax(cnv, 0.2) / base) / pmax(d$sd_null, params$sd_floor)
    l = ll_depth_shift(d$z, mu, tp)
    l[!is.finite(l)] = 0
    rowsum(l, ci_d)[as.character(seq_along(sc)), 1]
  }
  st_mat = function(k, col) { x = ps_[clone == k]; x[[col]][match(d$region, x$region)] }
  LL = matrix(0, length(sc), length(K), dimnames = list(sc, K))
  for (k in K) {
    cnv = st_mat(k, "cn"); plv = st_mat(k, "ploidy")
    cnv[is.na(cnv)] = 2; plv[is.na(plv)] = 2
    LL[, k] = depth_ll(cnv, plv, par[clone == k, c])
  }
  # alleles
  tv = trip_votes(fp$votes, fp$atoms, pieces, length(cells))
  vcell = match(sc, cells)
  member = lapply(K, function(k) match(clones[clone == k, cell], cells)); names(member) = K
  al_ll = function(oc, fexp_by_piece) {
    fe = fexp_by_piece[oc$ri]
    l = ll_allelic(oc$oA, oc$oN, fe, rho, params$eps, two_sided = FALSE, per_row = TRUE)
    l[is.na(fe) | abs(fe - 0.5) < 1e-9] = 0
    out = numeric(length(sc))
    s_ = rowsum(l, oc$ci)
    k_ = match(as.integer(rownames(s_)), vcell)
    out[k_[!is.na(k_)]] = s_[!is.na(k_), 1]
    out
  }
  fexp = function(k, orient_sign = NULL) {
    x = ps_[clone == k]
    a_ = par[clone == k, a]
    f = x$f_state[match(pieces$region, x$region)]
    f[is.na(f)] = 0.5
    if (!is.null(orient_sign)) f = ifelse(orient_sign < 0, 1 - f, f)
    0.5 + a_ * (f - 0.5)
  }
  occ = list()
  if (nrow(tv$trip)) for (k in K) {
    occ[[k]] = oriented_counts(tv, member[[k]])
    LL[, k] = LL[, k] + al_ll(occ[[k]], fexp(k))
  }
  hyp = K
  LM = NULL
  pr = if (length(K) >= 2) utils::combn(K, 2) else matrix(character(), 2, 0)
  if (!is.null(mix_pairs) && ncol(pr)) {
    mp = as.matrix(mix_pairs)
    ok_ = paste(pr[1, ], pr[2, ]) %in% c(paste(mp[, 1], mp[, 2]), paste(mp[, 2], mp[, 1]))
    pr = pr[, ok_, drop = FALSE]
  }
  if (mixtures && ncol(pr) >= 1) {
    wn = function(k, l, w) if (abs(w - 0.5) < 1e-9) paste(k, l, sep = "+") else sprintf("%s+%s@%.2f", k, l, w)
    LM = matrix(0, length(sc), ncol(pr) * length(mix_w),
                dimnames = list(sc, unlist(lapply(seq_len(ncol(pr)), function(j) vapply(mix_w, function(w) wn(pr[1, j], pr[2, j], w), "")))))
    col = 0L
    for (j in seq_len(ncol(pr))) {
      k = pr[1, j]; l = pr[2, j]
      ck = par[clone == k, c]; cl_ = par[clone == l, c]
      cnk = st_mat(k, "cn"); cnl = st_mat(l, "cn"); cnk[is.na(cnk)] = 2; cnl[is.na(cnl)] = 2
      plv = st_mat(k, "ploidy"); plv[is.na(plv)] = 2
      if (nrow(tv$trip)) {
        # each piece is scored under the orientation of the pair's clone that is imbalanced there (the larger one
        # if both are), i.e. the same orientations the single-clone hypotheses use: a pooled pair orientation
        # would be better determined than a single clone's and favour the mixtures systematically
        xk = ps_[clone == k]; xl = ps_[clone == l]
        ck_ = xk$cn[match(pieces$region, xk$region)]; cl2 = xl$cn[match(pieces$region, xl$region)]
        ck_[is.na(ck_)] = 2; cl2[is.na(cl2)] = 2
        fk0 = fexp(k); fl0 = fexp(l)
        imb_k = abs(fk0 - 0.5) > 1e-9; imb_l = abs(fl0 - 0.5) > 1e-9
        ori_k = imb_k & (!imb_l | length(member[[k]]) >= length(member[[l]]))
        ori_l = imb_l & !ori_k
        # haplotype agreement of the other clone under the orienting clone's orientation (pooled over its members)
        agree = function(oc, m) {
          x = oc[ci %in% m, .(oA = sum(oA), oN = sum(oN)), by = ri]
          s_ = rep(1, nrow(pieces)); s_[x$ri] = ifelse(x$oA / pmax(x$oN, 1) < 0.5, -1, 1)
          s_
        }
        sl_k = agree(occ[[k]], member[[l]]); sk_l = agree(occ[[l]], member[[k]])
        fl_k = 0.5 + (fl0 - 0.5) * sl_k        # l's major fraction under k's orientation
        fk_l = 0.5 + (fk0 - 0.5) * sk_l        # k's major fraction under l's orientation
        ock = occ[[k]][ri %in% which(ori_k)]; ocl = occ[[l]][ri %in% which(ori_l)]
      }
      for (w in mix_w) {
        col = col + 1L
        LM[, col] = depth_ll(w * cnk + (1 - w) * cnl, plv, w * ck + (1 - w) * cl_)
        if (nrow(tv$trip)) {
          den = w * ck_ + (1 - w) * cl2
          fm_k = ifelse(den > 0, (w * ck_ * fk0 + (1 - w) * cl2 * fl_k) / pmax(den, 1e-9), 0.5)
          fm_l = ifelse(den > 0, (w * ck_ * fk_l + (1 - w) * cl2 * fl0) / pmax(den, 1e-9), 0.5)
          LM[, col] = LM[, col] + al_ll(ock, fm_k) + al_ll(ocl, fm_l)
        }
      }
    }
    hyp = c(K, colnames(LM))
  }
  sizes = as.vector(table(clones$clone)[K])
  pk = if (is.null(prior)) rep(1 / length(K), length(K)) else as.vector(prior[K] / sum(prior[K]))
  # per-cell total mixture prior: mix_prior, raised to the cell's multiplet prior probability (e.g. a soft multiplet score)
  pm = rep(if (is.null(LM)) 0 else mix_prior, length(sc))
  if (!is.null(LM) && !is.null(mult_prior)) {
    mpv = as.numeric(mult_prior[sc]); mpv[is.na(mpv)] = 0
    pm = pmin(pmax(pm, mpv), max_mix_prior)
  }
  ALL = cbind(LL, LM)
  LP = matrix(log(pk), length(sc), length(K), byrow = TRUE) + log1p(-pm)
  if (!is.null(LM)) LP = cbind(LP, matrix(log(pm / ncol(LM)), length(sc), ncol(LM)))
  M = ALL + LP
  mx = apply(M, 1, max)
  P = exp(M - mx); P = P / rowSums(P)
  nk = length(K)
  best = max.col(LL + matrix(log(pk), nrow(LL), nk, byrow = TRUE), ties.method = "first")
  srt = t(apply(LL, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
  if (nk == 1) srt = cbind(LL[, 1], NA_real_)
  p_mix = if (is.null(LM)) rep(0, nrow(P)) else rowSums(P[, -seq_len(nk), drop = FALSE])
  bp = if (is.null(LM)) rep(NA_character_, nrow(P)) else colnames(LM)[max.col(LM, ties.method = "first")]
  # posterior over the single clones (membership), and doublet-like only when a mixture is both more probable
  # and better by `mix_margin` nats than the best single clone (many near-identical mixtures of sibling clones
  # otherwise soak up probability from cells that fit their clone slightly imperfectly)
  MS = M[, seq_len(nk), drop = FALSE]
  PS = exp(MS - apply(MS, 1, max)); PS = PS / rowSums(PS)
  dll = if (is.null(LM)) rep(-Inf, nrow(P)) else apply(LM, 1, max) - apply(LL, 1, max)
  # the misfit of the single-clone model grows with a cell's depth (e.g. deep cells with ~3x the median fragments), so the
  # mixture gain must also be a fraction `mix_rel` of the cell's best single-clone log-likelihood ratio
  llb = apply(LL, 1, max)
  out = data.table(cell = sc, clone = K[best], posterior = PS[cbind(seq_len(nrow(PS)), best)],
                   ll_margin = srt[, 1] - srt[, 2], doublet_like = p_mix > 0.5 & dll > pmax(mix_margin, mix_rel * pmax(llb, 0)),
                   pair = bp, ref_pool = FALSE,
                   p_mix = p_mix, ll_mix_gain = dll, ll_best = LL[cbind(seq_len(nrow(LL)), best)], mix_prior = pm)
  if (!is.null(same_groups)) {
    # posterior summed over the clones of the best clone's group (clones without a group are their own group)
    gk = as.character(same_groups[K]); gk[is.na(gk)] = K[is.na(gk)]
    GP = PS %*% (outer(gk, unique(gk), "==") * 1)
    out[, post_clone := posterior]
    out[, posterior := GP[cbind(seq_len(nrow(GP)), match(gk[best], unique(gk)))]]
  }
  out[, assigned := posterior >= min_post & !doublet_like]
  out[, confident := assigned]
  if (!is.null(trunks)) {
    # trunk-level membership: posterior summed over the clones of the best clone's trunk
    tk = as.character(trunks[K]); tk[is.na(tk)] = K[is.na(tk)]
    TP = PS %*% (outer(tk, unique(tk), "==") * 1)
    colnames(TP) = unique(tk)
    out[, trunk := tk[best]]
    out[, p_trunk := TP[cbind(seq_len(nrow(TP)), match(tk[best], colnames(TP)))]]
  }
  out[doublet_like == FALSE, pair := NA_character_]
  setattr(out, "ll", ALL)
  out[]
}
utils::globalVariables(c("rmean", "i.rmean", "a", "c", "ll_mix_gain", "post_clone"))

#' Flag clones that are 50:50 mixtures of two other clones
#'
#' For each clone m and each pair (k, l) of larger clones, the profile of m is scored under the mixture
#' of k's and l's segment states (copy number averaged; major fraction mixed with the better of same /
#' opposite haplotype per window, attenuated by m's `a`) with m's depth contrast, and under m's own
#' segments. A clone whose mixture log-likelihood is within `delta` of its own is a doublet clone.
#' @param clones data.table clone, n_cells (or cell, clone)
#' @param segments segments table
#' @param profiles clone profiles ([clone_profiles()])
#' @param fit optional [fit_clone_states()] result (for c, a, rho and the effective SE)
#' @param delta log-likelihood tolerance
#' @param df depth t degrees of freedom
#' @param pairs optional 2-column matrix of the clone pairs that may form a doublet clone (default all pairs);
#'   [caller_tumour()] passes the pairs of different trunks: a clone between two sibling clones of one trunk is a
#'   subclone with intermediate states, not a doublet population
#' @return data.table clone, best_pair, ll_self, ll_mix, delta_ll, doublet
#' @export
flag_doublet_clones = function(clones, segments, profiles, fit = NULL, delta = 10, df = 5, pairs = NULL) {
  cl = as.data.table(clones)
  if (!"n_cells" %in% names(cl)) cl = cl[, .(n_cells = .N), by = clone]
  cl[, clone := as.character(clone)]
  sg = as.data.table(segments)
  pr = as.data.table(if (!is.null(fit)) fit$windows else profiles)
  if (!"se_eff" %in% names(pr)) pr[, se_eff := log2_se]
  if (!"base" %in% names(pr)) pr[, base := 2L]
  par = if (!is.null(fit)) fit$par else data.table(clone = unique(pr$clone), c = 1, a = 1, rho = 0.01)
  if (nrow(cl) < 3) return(data.table(clone = cl$clone, best_pair = NA_character_, ll_self = NA_real_, ll_mix = NA_real_,
                                      delta_ll = NA_real_, doublet = FALSE))
  win_state = function(k, p) {
    s = sg[clone == k]
    j = s[p[, .(chrom, mid = (start + end) / 2)], on = .(chrom, start <= mid, end > mid), which = TRUE, mult = "first"]
    list(nA = as.integer(sub(":.*", "", s$state[j])), nB = as.integer(sub(".*:", "", s$state[j])))
  }
  emis = function(p, cn, fmaj, c_, a_, rho_, pl) {
    m = c_ * log2(pmax(cn, 0.2) / p$base)
    dl = stats::dt((p$log2 - m) / p$se_eff, df, log = TRUE) - log(p$se_eff)
    dl[!is.finite(dl)] = 0
    fj = 0.5 + a_ * (fmaj - 0.5)
    al = ifelse(p$n_hap > 0 & pl == 2, pmax(bb_robust(p$nA, p$n_hap, fj, rho_, 0.02), bb_robust(p$nA, p$n_hap, 1 - fj, rho_, 0.02)), 0)
    dl + al
  }
  res = rbindlist(lapply(cl$clone, function(m) {
    p = pr[clone == m]
    pm = par[clone == as.character(m)]
    sm = win_state(m, p)
    pl = p$ploidy %||% rep(2L, nrow(p))
    cn_m = sm$nA + sm$nB
    self = sum(emis(p, cn_m, ifelse(cn_m > 0, sm$nA / pmax(cn_m, 1), 0.5), pm$c, pm$a, pm$rho, pl), na.rm = TRUE)
    bigger = cl[clone != m & n_cells > cl[clone == m, n_cells], clone]
    if (length(bigger) < 2) return(data.table(clone = m, best_pair = NA_character_, ll_self = self, ll_mix = NA_real_))
    pp = utils::combn(bigger, 2)
    if (!is.null(pairs)) {
      ok_ = paste(pp[1, ], pp[2, ]) %in% c(paste(pairs[, 1], pairs[, 2]), paste(pairs[, 2], pairs[, 1]))
      pp = pp[, ok_, drop = FALSE]
      if (!ncol(pp)) return(data.table(clone = m, best_pair = NA_character_, ll_self = self, ll_mix = NA_real_))
    }
    mix = vapply(seq_len(ncol(pp)), function(j) {
      sk = win_state(pp[1, j], p); sl = win_state(pp[2, j], p)
      cnk = sk$nA + sk$nB; cnl = sl$nA + sl$nB
      tot = cnk + cnl
      same = ifelse(tot > 0, (sk$nA + sl$nA) / pmax(tot, 1), 0.5)
      opp = ifelse(tot > 0, (sk$nA + sl$nB) / pmax(tot, 1), 0.5)
      e1 = emis(p, tot / 2, same, pm$c, pm$a, pm$rho, pl); e2 = emis(p, tot / 2, opp, pm$c, pm$a, pm$rho, pl)
      sum(pmax(e1, e2), na.rm = TRUE)
    }, numeric(1))
    data.table(clone = m, best_pair = paste(pp[1, which.max(mix)], pp[2, which.max(mix)], sep = "+"), ll_self = self, ll_mix = max(mix))
  }), fill = TRUE)
  res[, delta_ll := ll_self - ll_mix]
  res[, doublet := !is.na(delta_ll) & delta_ll < delta]
  # a doublet clone must carry non-neutral segments (a normal clone is never dissolved)
  nn = sg[label != "neutral" & (end - start) >= 1e7, unique(clone)]
  res[!clone %in% nn, doublet := FALSE]
  res[]
}
utils::globalVariables(c("best_pair", "ll_self", "ll_mix", "delta_ll", "doublet"))

#' Clean the normal reference pool with the allele channel
#'
#' Reference pool = cells of the normal clone(s) (no non-neutral segment >= `min_event` bp) with posterior at
#' least `min_post`, not doublet-like with a tumour clone in the mixture (a mixture of two normal clones is not a reason
#' to exclude a cell), and without a per-cell allelic call on the clone events: for every tumour clone k and each of
#' its imbalanced non-neutral segments of >= `min_event` bp, the cell's votes are oriented by k's pseudobulk (so the
#' tumour's major haplotype is "A") and scored with the one-sided `ll_allelic(nA, n, f_grid, rho, two_sided = FALSE)`
#' (best excess of the tumour haplotype over the grid); the cell's event score is the maximum over these clone x
#' event regions, and a call is event score > `thr`. Without `thr` (`null = "max"`), it is the type-7 `q` quantile over
#' the candidate cells of the same maximum over the control regions (neutral in every clone; orientation by the same
#' clones), so the multiplicity of the maximum is in the null as well; `null = "region"` is the quantile of the
#' single-region control scores (the first build's rule, which flags too many cells when there are many events).
#' @param feat scn_features list (or `clone_input`)
#' @param assign assignments ([assign_cells()])
#' @param segments segments table
#' @param thr allelic call threshold (NULL = from the control-region null)
#' @param par data.table clone, c, a (attenuation; default a = 1)
#' @param members data.table cell, clone used for orientation (default the assigned cells)
#' @param min_post minimal posterior
#' @param min_event minimal length of a non-neutral segment that makes a clone non-normal
#' @param q null quantile
#' @param control optional control regions (default [control_regions()] with min_len 5 Mb)
#' @param f_grid major-haplotype fractions of the one-sided allelic score
#' @param null "max" or "region" (see above)
#' @return list(pool = cell names, table = per-cell event / control max score and flag, thr, normal = normal clones)
#' @export
clean_reference = function(feat, assign, segments, thr = NULL, par = NULL, members = NULL, min_post = 0.9,
                           min_event = 1e7, q = 0.999, control = NULL, f_grid = c(0.6, 2 / 3, 0.75, 0.85, 0.95),
                           null = c("max", "region")) {
  null = match.arg(null)
  inp = clone_input(feat)
  sg = as.data.table(segments)
  asg = as.data.table(assign)
  ab = sg[label != "neutral" & (end - start) >= min_event, unique(clone)]
  normal = setdiff(unique(sg$clone), ab)
  # doublet-like cells are excluded only when the best mixture involves a tumour clone
  pair_tumour = vapply(strsplit(sub("@.*", "", ifelse(is.na(asg$pair), "", asg$pair)), "+", fixed = TRUE),
                       function(p) any(p %in% ab), logical(1))
  cand = asg[clone %in% normal & posterior >= min_post & (doublet_like == FALSE | !pair_tumour), cell]
  if (!length(cand)) return(list(pool = character(), table = data.table(cell = character(), score = numeric(), flag = logical()),
                                 thr = thr %||% NA_real_, normal = normal))
  if (is.null(members)) members = asg[assigned == TRUE, .(cell, clone)]
  members = as.data.table(members)
  par = if (is.null(par)) data.table(clone = unique(sg$clone), a = 1) else as.data.table(par)
  ev = sg[!clone %in% normal & label != "neutral" & nA_of(state) != nB_of(state) & cn > 0 & (end - start) >= min_event]
  if (!nrow(ev)) return(list(pool = cand, table = data.table(cell = cand, score = 0, flag = FALSE), thr = thr %||% NA_real_, normal = normal))
  cells = inp$cells$cell
  ev_reg = seg_regions(ev[, .(clone, chrom, start, end)], inp$regions)
  ev_reg[, `:=`(eclone = ev$clone, fs = nA_of(ev$state) / ev$cn)]
  if (is.null(control)) control = control_regions(sg, min_len = 5e6, regions = inp$regions)
  ctl = as.data.table(control)
  ctl = ctl[!chrom %in% inp$haploid]
  all_r = unique(rbindlist(list(ev_reg[, names(ctl), with = FALSE], ctl), fill = TRUE), by = "region")
  tv = trip_votes(inp$votes, inp$atoms, all_r, length(cells), inp$allelic_chroms)
  ci_c = match(cand, cells)
  scores = rbindlist(lapply(unique(ev_reg$eclone), function(k) {
    if (!nrow(tv$trip)) return(data.table())
    oc = oriented_counts(tv, match(members[clone == k, cell], cells))
    oc = oc[ci %in% ci_c]
    oc[, region := all_r$region[ri]]
    oc = oc[region %in% c(ev_reg[eclone == k, region], ctl$region)]
    oc[, `:=`(cell = cells[ci], score = ll_allelic(oA, oN, f_grid, inp$rho, two_sided = FALSE), is_ctl = region %in% ctl$region)]
    oc
  }))
  obs = scores[is_ctl == FALSE]
  nul = scores[is_ctl == TRUE]
  tab = data.table(cell = cand)
  if (nrow(obs)) tab = merge(tab, obs[, .(score = max(score), n_ev = .N), by = cell], by = "cell", all.x = TRUE) else tab[, `:=`(score = 0, n_ev = 0L)]
  if (nrow(nul)) tab = merge(tab, nul[, .(score_ctl = max(score), n_ctl = .N), by = cell], by = "cell", all.x = TRUE) else tab[, `:=`(score_ctl = 0, n_ctl = 0L)]
  tab[is.na(score), `:=`(score = 0, n_ev = 0L)]
  tab[is.na(score_ctl), `:=`(score_ctl = 0, n_ctl = 0L)]
  if (is.null(thr)) {
    thr = if (!nrow(nul)) Inf else if (null == "max") unname(stats::quantile(tab$score_ctl, q, type = 7)) else
      unname(stats::quantile(nul$score, q, type = 7))
  }
  tab[, flag := score > thr]
  list(pool = tab[flag == FALSE, cell], table = tab, thr = thr, normal = normal)
}
utils::globalVariables(c("eclone", "fs", "is_ctl", "score_ctl", "n_ev", "n_ctl"))

nA_of = function(state) as.integer(sub(":.*", "", state))
nB_of = function(state) as.integer(sub(".*:", "", state))

#' Control regions of a tumour: pieces neutral in every clone
#'
#' The automatic tumour null: the pieces ([segment_pieces()]) that are neutral in every clone and at least
#' `min_len` long. User-supplied control regions override these in [caller_tumour()].
#' @param segments segments table
#' @param min_len minimal piece length (bp)
#' @param regions optional regions table (for sd_class / allelic_ok)
#' @return regions table (class `segment`, names `seg|control|...`)
#' @export
control_regions = function(segments, min_len = 2e7, regions = NULL) {
  sg = as.data.table(segments)
  pcs = segment_pieces(sg, regions)
  st = piece_states(pcs, sg)
  bad = st[label != "neutral", unique(region)]
  out = pcs[!region %in% bad & (end - start) >= min_len]
  out[, region := sub("^seg\\|all\\|", "seg|control|", region)]
  out[, null_region := TRUE]
  out[]
}

# fixed normal clone N0: all-neutral segments over the chromosomes of the fit (1:0 on haploid chromosomes), replacing
# every fitted clone without a non-neutral segment >= min_event (returned in `replaced`); depth contrast = the median
# of the fitted clones, attenuation 1. st: run_states() result (profiles, fit, segments, ref)
add_normal_clone = function(st, min_event = 1e7) {
  sg = as.data.table(st$segments)
  ev = sg[label != "neutral" & (end - start) >= min_event, unique(clone)]
  repl = setdiff(unique(sg$clone), c(ev, "N0"))
  chr = sg[, .(start = min(start), end = max(end), ploidy = as.integer(ploidy[1])), by = chrom]
  n0 = chr[, .(clone = "N0", chrom, start, end, n_win = 1L, nA = NA_integer_, nB = NA_integer_,
               cn = as.integer(fifelse(ploidy == 1L, 1L, 2L)), state = fifelse(ploidy == 1L, "1:0", "1:1"), log2 = 0, f = 0.5,
               llr = 0, confident = TRUE, label = "neutral", ploidy)]
  out = rbind(sg[!clone %in% c(repl, "N0")], n0, fill = TRUE)
  st$segments = out[order(clone, chrom_rank(chrom), start)]
  par = as.data.table(st$fit$par)
  p0 = data.table(clone = "N0", c = if (nrow(par)) stats::median(par$c) else 1, a = 1,
                  rho = if (nrow(par)) stats::median(par$rho) else 0.01, s_sys = 0)
  st$fit$par = rbind(par[!clone %in% c(repl, "N0")], p0, fill = TRUE)
  # the profile of the replaced clones (pooled) is kept under N0 for reporting
  if (!is.null(st$profiles) && length(repl)) {
    pr = as.data.table(st$profiles)
    pn = pr[clone %in% repl][order(-n_cells)][!duplicated(window)]
    if (nrow(pn)) pn[, clone := "N0"]
    st$profiles = rbind(pr[!clone %in% c(repl, "N0")], pn, fill = TRUE)
  }
  list(st = st, replaced = repl)
}
