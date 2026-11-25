# medicc2_cells.R
# What: per-cell allele-specific copy-number states on the common MEDICC2 pieces of a tumour run (per-cell MEDICC2
#       input). Every cell x piece is scored with the cell's own evidence only (its window depth over the piece and
#       its molecule votes oriented by the export's block orientation) against the oriented (cn_a, cn_b) states the
#       kept clones have on that piece, plus 1:1; the MAP state, its posterior and the LLR vs 1:1 are returned.
# Why: a per-cell MEDICC2 tree needs one allele-specific state per cell and piece in the same haplotype convention as
#      the clone export (cn_a = the same haplotype in every sample). Inheriting the clone state would make the tree
#      a copy of the clone assignment; scoring each cell independently lets the tree test it.
# Inputs: scn_features (win, long, regions, votes, atoms, cells), the export's pieces (orient_clone_states()$pieces /
#         pieces.tsv) and block orientation (orient_clone_states()$blocks / block_orientation.rds), the run's clones
#         table (c, a, n_assigned) and null fit (tpar, rho), the reference pool. Outputs: data.tables.

# weighted median (weights >= 0; NA values dropped)
wmedian = function(x, w = rep(1, length(x))) {
  ok = is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x = x[ok]; w = w[ok]
  o = order(x); x = x[o]; w = cumsum(w[o]) / sum(w)
  x[which(w >= 0.5)[1]]
}

#' Candidate per-cell states of the MEDICC2 pieces
#'
#' Per piece, the distinct oriented `(cn_a, cn_b)` of the clones plus the neutral state (1:1 on diploid pieces, 1:0
#' on haploid pieces such as male chrX; haploid pieces keep only `cn_b = 0` states). The prior is uniform over a
#' piece's candidates (`uniform`), `p_neutral` on the neutral state and the rest split evenly (`neutral`), or
#' (`clone_freq`) proportional to the summed size (`n_assigned`, else 1 per clone) of the clones carrying the state plus
#' a pseudo-count of `pseudo` x the total size (so that a state no clone has, e.g. the neutral state of an all-gained
#' piece, keeps a small prior).
#' @param pieces long table clone x piece with `piece, chrom, start, end, clone, cn_a, cn_b` (+ `ploidy`, default 2)
#' @param clones optional clones table (`clone`, `n_assigned`) for the `clone_freq` prior
#' @param prior `uniform`, `clone_freq` or `neutral` (`p_neutral` on the neutral state, the rest split evenly over
#'   the other candidates: the clone model's state prior, [cn_grid()], also gives the neutral state 0.5)
#' @param p_neutral neutral prior mass of `prior = "neutral"`
#' @param pseudo pseudo-count (fraction of the total clone size) of the `clone_freq` prior
#' @return data.table piece, chrom, start, end, ploidy, cn_a, cn_b, cn, state (`cn_a:cn_b`), neutral, n_clones
#'   (clones with the state), w_clones (their summed size), prior (sums to 1 per piece)
#' @export
cell_piece_candidates = function(pieces, clones = NULL, prior = c("uniform", "clone_freq", "neutral"), pseudo = 0.01, p_neutral = 0.5) {
  prior = match.arg(prior)
  pl = as.data.table(pieces)
  if (!all(c("piece", "chrom", "start", "end", "clone", "cn_a", "cn_b") %in% names(pl)))
    stop("cell_piece_candidates: pieces need piece, chrom, start, end, clone, cn_a, cn_b")
  pl = pl[, .(piece, chrom, start = as.numeric(start), end = as.numeric(end), clone = as.character(clone),
              cn_a = as.integer(cn_a), cn_b = as.integer(cn_b), ploidy = if ("ploidy" %in% names(pl)) as.integer(ploidy) else 2L)]
  w = rep(1, nrow(pl))
  if (!is.null(clones) && "n_assigned" %in% names(clones)) {
    cs = as.data.table(clones)
    w = as.numeric(cs$n_assigned[match(pl$clone, as.character(cs$clone))])
    w[!is.finite(w)] = 0
  }
  pl[, w := w]
  pp = pl[, .(chrom = chrom[1], start = start[1], end = end[1], ploidy = min(ploidy)), by = piece]
  hap = pp[ploidy == 1L, piece]
  # haploid pieces: one-copy states only (cn_b = 0)
  pl = pl[!(piece %in% hap & cn_b != 0L)]
  st = pl[, .(n_clones = .N, w_clones = sum(w)), by = .(piece, cn_a, cn_b)]
  neu = pp[, .(piece, cn_a = 1L, cn_b = fifelse(ploidy == 1L, 0L, 1L), n_clones = 0L, w_clones = 0)]
  st = rbind(st, neu[!st, on = .(piece, cn_a, cn_b)])
  st = pp[st, on = "piece"]
  st[, `:=`(cn = cn_a + cn_b, state = paste0(cn_a, ":", cn_b), neutral = cn_a == 1L & cn_b == fifelse(ploidy == 1L, 0L, 1L))]
  if (prior == "uniform") st[, prior := 1 / .N, by = piece] else if (prior == "neutral") {
    st[, prior := if (.N == 1L) 1 else fifelse(neutral, p_neutral, (1 - p_neutral) / (.N - 1L)), by = piece]
  } else {
    tot = sum(unique(pl[, .(clone, w)])$w)
    st[, prior := (w_clones + pseudo * max(tot, 1)) / sum(w_clones + pseudo * max(tot, 1)), by = piece]
  }
  st = st[order(chrom_rank(chrom), start, -neutral, cn, cn_a)]
  setcolorder(st, c("piece", "chrom", "start", "end", "ploidy", "cn_a", "cn_b", "cn", "state", "neutral", "n_clones", "w_clones", "prior"))
  st[]
}
utils::globalVariables(c("w", "w_clones", "n_clones", "neutral"))

#' Per-cell evidence on the MEDICC2 pieces: depth and oriented molecule votes
#'
#' Depth (the per-cell analogue of [clone_profiles()]): the cell's 1 Mb window counts summed over the piece (windows
#' fully inside it; [region_features()]), `log2ratio = log2((x + 0.5) / (e + 0.5))` relative to the reference
#' pool's pseudobulk of the piece (`ref_center = "pseudobulk"`, `log2((sum x + 0.5) / (sum e + 0.5))` over the
#' reference cells, as the clone profiles; `"mean"` = the mean of the reference cells' log2 ratios, as
#' [assign_cells()], which a few reference cells with a large negative offset pull down), then (`center`) minus the
#' cell's own baseline: the window-weighted median of its log2 ratios over the autosomal pieces (the clone profiles are
#' centred on the clone's autosomal median in the same way, so states are relative to the baseline). SD: the cell's
#' null SD of the piece's size class (`sd_null` of the piece's `sd_class`); with `sd_scale`, a piece shorter than its
#' tile class (e.g. a 2 Mb piece scored with the 5 Mb SD) gets `sd_null * sqrt(class size / covered window length)`
#' (the independent-window limit); pieces on haploid chromosomes get `haploid_sd` x that; `s_sys` (systematic depth
#' bias of a region between the scored cells and the reference pool, the clone model's `s_sys`) is added in
#' quadrature: `sd = sqrt(sd_scaled^2 + s_sys^2)`.
#' `z = log2ratio / pmax(sd, sd_floor)`; pieces without a window have NA depth (no evidence).
#' Alleles: the cell's molecule votes on the piece (atom start in `[start, end)`) per phase set, oriented by the
#' block table: `oA` = votes of the cn_a haplotype (`nA` when `sign > 0`, `nB` when `sign < 0`), `oN` = oriented votes.
#' With `loo`, the reference clone's pooled votes are first reduced by the cell's own votes when the cell was one of
#' the cells that oriented the block (`loo` maps cells to the clone whose pool they joined; only blocks whose `ref`
#' is that clone are affected): the orientation is then the one the other cells give, in the block table's current
#' label convention (on a piece whose labels the chaining of [orient_clone_states()] flipped, `flipped`, the sign of
#' the remaining votes is negated as `sign` is; a table without `flipped` takes the flip from the piece's non-tied
#' blocks), and the cell's own votes cannot confirm themselves. This is the only use of a cell's clone membership, and it only removes the cell's own data
#' from the orientation input.
#' @param features scn_features list (`win`, `long`, `regions`, `votes`, `atoms`, `cells`; as the tumour mode reads it)
#' @param pieces piece table (`piece, chrom, start, end`, one row per piece or per clone x piece)
#' @param blocks block orientation (`piece, ps, sign` and, for `loo`, `ref, nA_ref, nB_ref`; `flipped` optional)
#' @param cells cells to score (default all cells of `features`)
#' @param ref_cells reference cells of the depth centring (>= 10 used; NULL = no centring)
#' @param loo NULL, or data.table `cell, clone`: the cells whose votes formed each reference clone's pool
#' @param params [caller_defaults()] (`sd_floor`)
#' @param sd_scale scale the null SD of pieces shorter than their tile class
#' @param ref_center reference centring: `pseudobulk` or `mean`
#' @param center subtract the cell's autosomal median (window-weighted over pieces)
#' @param s_sys systematic SD (log2) added in quadrature to the per-cell SD
#' @param haploid_sd factor on the null SD of pieces on haploid chromosomes (male chrX / chrY): `sd_null` comes from
#'   diploid regions, and one copy gives half the fragments, so z on haploid pieces is overdispersed
#' @return data.table cell, piece, x, e, n_win (1 Mb windows), log2ratio (centred), offset (the cell's subtracted
#'   baseline), sd_null, sd, z, stratum, sd_class, n_frag, oA, oN, n_blocks (oriented phase sets with votes)
#' @export
cell_piece_data = function(features, pieces, blocks, cells = NULL, ref_cells = NULL, loo = NULL, params = caller_defaults(),
                           sd_scale = TRUE, ref_center = c("pseudobulk", "mean"), center = TRUE, s_sys = 0, haploid_sd = sqrt(2)) {
  ref_center = match.arg(ref_center)
  inp = clone_input(features)
  pc = unique(as.data.table(pieces)[, .(piece, chrom, start = as.numeric(start), end = as.numeric(end))], by = "piece")
  pc = pc[order(chrom_rank(chrom), start)]
  all_cells = inp$cells$cell
  cells = cells %||% all_cells
  miss = setdiff(cells, all_cells)
  if (length(miss)) stop(sprintf("cell_piece_data: %d cells are not in the features (e.g. %s)", length(miss), miss[1]))
  # ---- depth ----
  regs = seg_regions(pc[, .(chrom, start, end)], inp$regions)
  regs[, region := pc$piece]
  ref = if (length(ref_cells) >= 10) intersect(ref_cells, colnames(inp$win$x)) else character()
  d = region_features(inp, regs, cells = union(cells, ref), sd_floor = params$sd_floor)
  wt = inp$win$windows
  hit = wt[, .(chrom, ws = as.numeric(start), we = as.numeric(end))][regs[, .(region, chrom, rs = start, re = end)],
                                                                     on = .(chrom, ws >= rs, we <= re), .(region = i.region),
                                                                     nomatch = NULL, allow.cartesian = TRUE]
  nw = hit[, .(n_win = .N), by = region]
  regs[, n_win := 0L]
  regs[nw, n_win := i.n_win, on = "region"]
  d[regs, `:=`(n_win = i.n_win, len = i.end - i.start, sd_class = i.sd_class), on = "region"]
  d[n_win == 0L | e <= 0, log2ratio := NA_real_]
  if (length(ref)) {
    rm_ = if (ref_center == "mean") d[cell %in% ref, .(rmean = mean(log2ratio, na.rm = TRUE)), by = region] else
      d[cell %in% ref & !is.na(log2ratio), .(rmean = log2((sum(x) + 0.5) / (sum(e) + 0.5))), by = region]
    d[rm_, log2ratio := log2ratio - i.rmean, on = "region"]
  }
  d[, offset := 0]
  if (center) {
    auto = !regs$chrom %in% c(inp$haploid, "chrX", "chrY")
    ar = regs$region[auto]
    if (length(ar)) {
      cm_ = d[region %in% ar & is.finite(log2ratio) & n_win > 0, .(off = wmedian(log2ratio, n_win)), by = cell]
      d[cm_, offset := i.off, on = "cell"]
      d[, log2ratio := log2ratio - offset]
    }
  }
  ts = params$tile_sizes %||% caller_defaults()$tile_sizes
  d[, sd := sd_null]
  if (sd_scale) {
    csz = ts[d$sd_class]
    k = !is.na(csz) & d$n_win > 0
    d[k, sd := sd_null * sqrt(pmax(1, csz[k] / pmax(n_win * 1e6, 1)))]
  }
  if (length(inp$haploid) && haploid_sd != 1) d[chrom %in% inp$haploid, sd := sd * haploid_sd]
  if (is.finite(s_sys) && s_sys > 0) d[, sd := sqrt(sd^2 + s_sys^2)]
  d[, z := log2ratio / pmax(sd, params$sd_floor)]
  d = d[cell %in% cells, .(cell, piece = region, x, e, n_win, log2ratio, offset, sd_null, sd, z, stratum = as.character(stratum),
                           sd_class, n_frag)]
  # ---- alleles ----
  out_al = data.table(cell = character(), piece = character(), oA = numeric(), oN = numeric(), n_blocks = integer())
  votes = inp$votes
  bl = copy(as.data.table(blocks))                    # columns are added below: never modify the caller's table
  if (!is.null(votes) && nrow(votes) && nrow(bl)) {
    ci = match(cells, all_cells)
    atoms = inp$atoms
    pc[, pix := seq_len(.N)]
    am = atoms[, .(atom, chrom, pos = start)][pc, on = .(chrom, pos >= start, pos < end), .(atom = x.atom, pix = i.pix), nomatch = NULL]
    at_pi = rep(NA_integer_, max(atoms$atom))
    at_pi[am$atom] = am$pix
    sel = votes$cell %in% ci
    v = data.table(cell = votes$cell[sel], pix = at_pi[votes$atom[sel]], ps = as.numeric(votes$ps[sel]),
                   nA = as.numeric(votes$nA[sel]), nB = as.numeric(votes$nB[sel]))
    v = v[!is.na(pix), .(nA = sum(nA), nB = sum(nB)), by = .(cell, pix, ps)]
    bl[, `:=`(pix = match(piece, pc$piece), ps = as.numeric(ps))]
    bl = bl[!is.na(pix)]
    if (!"sign" %in% names(bl)) stop("cell_piece_data: blocks need `sign`")
    if (!is.null(loo) && all(c("nA_ref", "nB_ref") %in% names(bl))) {
      # the piece's label flip (orient_clone_states() chaining): sign = phi_ * sign(nA_ref - nB_ref), phi_ = -1 on flipped
      # pieces; from `flipped`, else from the piece's non-tied oriented blocks (1 when it has none)
      if ("flipped" %in% names(bl)) bl[, phi_ := fifelse(flipped %in% TRUE, -1L, 1L)] else
        bl[, phi_ := { k = which(sign != 0 & nA_ref != nB_ref); if (length(k)) as.integer(sign[k[1]] * sign(nA_ref[k[1]] - nB_ref[k[1]])) else 1L },
           by = pix]
    }
    cols = intersect(c("ref", "nA_ref", "nB_ref", "phi_"), names(bl))
    v = bl[, c("pix", "ps", "sign", cols), with = FALSE][v, on = .(pix, ps), nomatch = NULL]
    if (!is.null(loo) && nrow(v)) {
      if (!all(c("ref", "nA_ref", "nB_ref") %in% names(v))) stop("cell_piece_data: `loo` needs blocks with ref, nA_ref, nB_ref")
      lc = as.data.table(loo)
      own = as.character(lc$clone)[match(all_cells[v$cell], as.character(lc$cell))]
      m = !is.na(own) & own == v$ref
      v[, sign := as.integer(fifelse(m, phi_ * sign((nA_ref - nB_ref) - (nA - nB)), sign))]
    }
    v = v[!is.na(sign) & sign != 0L]
    if (nrow(v)) {
      out_al = v[, .(oA = sum(fifelse(sign > 0, nA, nB)), oN = sum(nA + nB), n_blocks = .N), by = .(cell, pix)]
      out_al[, `:=`(cell = all_cells[cell], piece = pc$piece[pix])]
      out_al = out_al[, .(cell, piece, oA, oN, n_blocks)]
    }
  }
  d[out_al, `:=`(oA = i.oA, oN = i.oN, n_blocks = i.n_blocks), on = .(cell, piece)]
  if (!"oA" %in% names(d)) d[, `:=`(oA = NA_real_, oN = NA_real_, n_blocks = NA_integer_)]
  d[is.na(oN), `:=`(oA = 0, oN = 0, n_blocks = 0L)]
  d[]
}
utils::globalVariables(c("n_win", "i.n_win", "len", "i.end", "i.start", "sd", "pix", "piece", "oA", "i.oA", "i.oN", "i.n_blocks",
                         "rs", "re", "nA_ref", "nB_ref", "ws", "we", "offset", "off", "i.off", "phi_", "flipped"))

#' Per-cell copy-number states on the MEDICC2 pieces
#'
#' For each cell x piece, the candidates of [cell_piece_candidates()] (1:1 plus every oriented `(cn_a, cn_b)` a kept
#' clone has on the piece; haploid pieces: one-copy states) are scored with the cell's own evidence
#' ([cell_piece_data()]); the cell's clone assignment is not used.
#' - Depth: `ll_depth_shift(z, mu, tpar)` with `mu = c * log2(max(cn, 0.2) / base) / pmax(sd, sd_floor)` (the shift
#'   model of the clone states, [fit_clone_states()] / [assign_cells()]: `base` = the piece ploidy when the depth is
#'   centred on the reference pool, else 2), relative to the neutral state's shift (the neutral state's own depth
#'   term is `ll_depth_shift(z, 0) = 0`), Student-t null per sd_class x stratum ([fit_depth_null()] parameters,
#'   fallbacks as in [assign_cells()]) or Gaussian when `tpar` is NULL.
#' - Alleles: one-sided beta-binomial LLR vs 1:1 of the oriented votes `(oA, oN)` at the expected cn_a fraction
#'   `0.5 + a * (cn_a / cn - 0.5)` (`a` = orientation attenuation of the clone model; 0.5 for cn = 0), overdispersion
#'   `rho`; haploid pieces have no allele term.
#' - MAP over the candidates with the prior (`neutral`, `uniform` or `clone_freq`); `posterior` is the MAP's posterior and
#'   `llr_neutral` the log-likelihood (no prior) of the MAP minus that of the neutral state.
#'
#' `c`, `a` and `s_sys` are global (one value for every cell): by default the median over the tumour clones (kept clones
#' of `pieces` except `normal`), weighted by their assigned cells, of the clones' fitted `c` / `a` / `s_sys` (`clones`
#' table). A per-cell value would have to come from the cell's clone, which this function must not use.
#' @param features scn_features list (see [cell_piece_data()])
#' @param pieces clone x piece table with `cn_a, cn_b` (orient_clone_states()$pieces or the export's pieces.tsv), or a
#'   `medicc2_orientation`
#' @param blocks block orientation (orient_clone_states()$blocks); NULL takes it from a `medicc2_orientation`
#' @param cells cells to score (default all)
#' @param params [caller_defaults()] (`sd_floor`, `eps`, `tile_sizes`)
#' @param prior `neutral` (default), `uniform` or `clone_freq` (see [cell_piece_candidates()]). Under the uniform prior
#'   a short piece with many candidates and little evidence goes to whichever non-neutral state fits marginally better,
#'   so many rows of normal cells get non-neutral states; `clone_freq` gives uninformative cells the majority clone
#'   state (tumour states in normal cells)
#' @param p_neutral neutral prior mass of `prior = "neutral"`
#' @param clones clones table (`clone, c, a, n_assigned`); used for the default `c` / `a` and the `clone_freq` prior
#' @param c,a global depth contrast and allelic attenuation (NULL = cell-weighted median of the tumour clones; 1 when
#'   `clones` lacks them)
#' @param rho beta-binomial overdispersion of a cell's oriented votes on a piece: a number; `"ref"` = fitted
#'   ([fit_bb_rho()], f = 0.5) on the reference cells' oriented votes (diploid pieces; the reference pool is balanced
#'   everywhere, so it is the null of the allele channel, as it is the per-cell null of the caller); `"auto"` = the value
#'   of `rho_grid` that maximises the marginal likelihood of all scored cell x piece rows over their candidates (with the
#'   prior); NULL = the features' / run's allelic null rho. The run's null rho is fitted on votes oriented by the
#'   all-cell pseudobulk, which on event-free regions is a random orientation; in the true haplotype frame of the block
#'   orientation a single cell's votes on a piece are far more dispersed (often by an order of magnitude).
#' @param rho_grid grid of `rho = "auto"`
#' @param tpar depth null t parameters (class, stratum, t_m, t_s, t_df) or NULL (Gaussian)
#' @param ref_cells reference pool for the depth centring
#' @param loo NULL or data.table cell, clone (see [cell_piece_data()])
#' @param normal name of the normal clone (excluded from the default c / a)
#' @param sd_scale,ref_center,center,haploid_sd see [cell_piece_data()]
#' @param s_sys systematic depth SD added to the per-cell SD (NULL = the cell-weighted median of the tumour clones'
#'   `s_sys`, 0 when `clones` lacks it); see [cell_piece_data()]
#' @param data precomputed [cell_piece_data()] output (skips the depth / vote aggregation)
#' @param keep_ll return the per-candidate log-likelihoods
#' @return list of class `cell_piece_states`: `states` (one row per cell x piece: cell, piece, chrom, start, end,
#'   ploidy, cn_a, cn_b, state (MAP), posterior, llr_neutral, ll_depth, ll_allele (MAP terms vs neutral), n_cand, z,
#'   log2ratio, offset, sd, n_win, oA, oN (oriented votes), f_obs (oA / oN), n_blocks), `candidates`, `ll` (cell x piece x
#'   candidate, when `keep_ll`), `params` (c, a, rho, rho_mode, rho_rows (reference rows of the rho fit), rho_ml (the
#'   `auto` grid), s_sys, ref_center, center, prior, n_ref, loo, sd_scale, depth_null)
#' @export
cell_piece_states = function(features, pieces, blocks = NULL, cells = NULL, params = caller_defaults(),
                             prior = c("neutral", "uniform", "clone_freq"), p_neutral = 0.5, clones = NULL, c = NULL, a = NULL, rho = "ref", tpar = NULL,
                             ref_cells = NULL, loo = NULL, normal = "N0", sd_scale = TRUE, ref_center = c("pseudobulk", "mean"),
                             center = TRUE, s_sys = NULL, haploid_sd = sqrt(2), data = NULL, keep_ll = FALSE,
                             rho_grid = c(0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25)) {
  ref_center = match.arg(ref_center)
  prior = match.arg(prior)
  if (inherits(pieces, "medicc2_orientation")) {
    blocks = blocks %||% pieces$blocks
    pieces = pieces$pieces
  }
  pl = as.data.table(pieces)
  if (is.null(blocks) && is.null(data)) stop("cell_piece_states: need `blocks` (block orientation)")
  K = setdiff(unique(as.character(pl$clone)), normal)
  if (!is.null(clones)) {
    cs = as.data.table(clones)[as.character(clone) %in% K]
    wt_ = if ("n_assigned" %in% names(cs)) as.numeric(cs$n_assigned) else rep(1, nrow(cs))
    wt_[!is.finite(wt_)] = 0
    if (is.null(c) && "c" %in% names(cs)) c = wmedian(as.numeric(cs$c), wt_)
    if (is.null(a) && "a" %in% names(cs)) a = wmedian(as.numeric(cs$a), wt_)
    if (is.null(s_sys) && "s_sys" %in% names(cs)) s_sys = wmedian(as.numeric(cs$s_sys), wt_)
  }
  if (is.null(c) || !is.finite(c)) c = 1
  if (is.null(a) || !is.finite(a)) a = 1
  if (is.null(s_sys) || !is.finite(s_sys)) s_sys = 0
  c_ = c; a_ = a
  use_ref = length(ref_cells) >= 10
  rho_mode = if (is.character(rho)) match.arg(rho, c("ref", "auto")) else "given"
  if (rho_mode == "ref" && !use_ref) stop("cell_piece_states: rho = \"ref\" needs >= 10 `ref_cells`")
  if (is.null(rho)) {
    rho = if (!is.null(features)) clone_input(features)$rho else NA_real_
    rho_mode = "null_fit"
  }
  if (rho_mode %in% c("given", "null_fit") && !is.finite(rho)) stop("cell_piece_states: need `rho`")
  cand = cell_piece_candidates(pl, clones, prior, p_neutral = p_neutral)
  dcells = if (!is.null(cells) && rho_mode == "ref") union(cells, ref_cells) else cells
  d = as.data.table(data %||% cell_piece_data(features, pl, blocks, dcells, ref_cells, loo, params, sd_scale, ref_center, center, s_sys,
                                              haploid_sd))
  if (!"offset" %in% names(d)) d[, offset := 0]
  pp = unique(cand[, .(piece, chrom, start, end, ploidy)])
  d = d[piece %in% pp$piece]
  d[pp, `:=`(ploidy = i.ploidy, chrom = i.chrom, start = i.start, end = i.end), on = "piece"]
  rho_rows = NA_integer_
  if (rho_mode == "ref") {
    rr = d[cell %in% ref_cells & ploidy == 2L & oN > 0]
    if (nrow(rr) < 50) stop("cell_piece_states: fewer than 50 reference cell x piece rows with oriented votes for rho = \"ref\"")
    rho = fit_bb_rho(rr$oA, rr$oN)
    rho_rows = nrow(rr)
  }
  if (!is.null(cells)) d = d[cell %in% cells]
  # long table cell x piece x candidate
  d[, row := .I]
  L = cand[, .(piece, cn_a, cn_b, cn, state, neutral, prior)][d[, .(row, cell, piece, ploidy, z, sd, sd_class, stratum, oA, oN)],
                                                                  on = "piece", allow.cartesian = TRUE]
  tp = cl_match_tpar(L, tpar)
  base = if (use_ref) L$ploidy else 2
  L[, mu := c_ * log2(pmax(cn, 0.2) / base) / pmax(sd, params$sd_floor)]
  mu0 = c_ * log2(fifelse(L$ploidy == 1L, 1, 2) / base) / pmax(L$sd, params$sd_floor)
  L[, ll_depth := ll_depth_shift(z, mu, tp) - ll_depth_shift(z, mu0, tp)]
  L[!is.finite(ll_depth), ll_depth := 0]
  L[, fexp := fifelse(cn > 0, 0.5 + a_ * (cn_a / pmax(cn, 1) - 0.5), 0.5)]
  al = function(r_) fifelse(L$ploidy == 1L | abs(L$fexp - 0.5) < 1e-9, 0, ll_allelic(L$oA, L$oN, L$fexp, r_, params$eps, two_sided = FALSE, per_row = TRUE))
  rho_ml = NULL
  if (rho_mode == "auto") {
    # marginal log-likelihood of every row over its candidates (allele part absolute: + the neutral's bb_loglik)
    has = d$ploidy == 2L & d$oN > 0
    rho_ml = data.table(rho = rho_grid, loglik = vapply(rho_grid, function(r_) {
      L[, lp_ := ll_depth + al(r_) + log(prior)]
      lse = L[, { m = max(lp_); m + log(sum(exp(lp_ - m))) }, by = row]$V1
      sum(lse) + sum(bb_loglik(d$oA[has], d$oN[has], 0.5, r_, params$eps))
    }, numeric(1)))
    rho = rho_ml$rho[which.max(rho_ml$loglik)]
    L[, lp_ := NULL]
  }
  L[, ll_allele := al(rho)]
  L[, ll := ll_depth + ll_allele]
  L[, lp := ll + log(prior)]
  L[, post := { m = max(lp); p = exp(lp - m); p / sum(p) }, by = row]
  best = L[order(row, -lp, -neutral)][!duplicated(row)]
  neu = L[neutral == TRUE, .(row, ll_neu = ll)]
  best[neu, ll_neu := i.ll_neu, on = "row"]
  nc = L[, .(n_cand = .N), by = row]
  st = d[best[, .(row, cn_a, cn_b, state, posterior = post, llr_neutral = ll - ll_neu, ll_depth, ll_allele)], on = "row"]
  st[nc, n_cand := i.n_cand, on = "row"]
  st[, f_obs := fifelse(oN > 0, oA / oN, NA_real_)]
  st = st[order(cell, chrom_rank(chrom), start), .(cell, piece, chrom, start, end, ploidy, cn_a, cn_b, state, posterior, llr_neutral,
                                                    ll_depth, ll_allele, n_cand, z, log2ratio, offset, sd, n_win, oA, oN, f_obs, n_blocks)]
  ll_out = if (keep_ll) L[d[, .(row, cell)], on = "row"][, .(cell = i.cell, piece, state, cn_a, cn_b, prior, ll_depth, ll_allele, ll, post)] else NULL
  structure(list(states = st[], candidates = cand, ll = ll_out,
                 params = list(c = c_, a = a_, p_neutral = p_neutral, haploid_sd = haploid_sd, rho = rho, rho_mode = rho_mode, rho_rows = rho_rows, rho_ml = rho_ml, s_sys = s_sys, ref_center = ref_center, center = center, prior = prior, n_ref = if (use_ref) length(ref_cells) else 0L,
                               loo = !is.null(loo), sd_scale = sd_scale, depth_null = if (is.null(tpar)) "gauss" else "t")),
            class = "cell_piece_states")
}
utils::globalVariables(c("lp_", "V1", "row", "mu", "ll_depth", "fexp", "ll_allele", "lp", "post", "ll_neu", "i.ll_neu", "n_cand", "i.n_cand",
                         "f_obs", "llr_neutral", "i.ploidy", "i.cell", "posterior"))

#' @export
print.cell_piece_states = function(x, ...) {
  s = x$states
  cat(sprintf("<cell_piece_states> %d cells x %d pieces\n", uniqueN(s$cell), uniqueN(s$piece)))
  cat(sprintf("  c %.3g, a %.3g, rho %.3g (%s), s_sys %.3g, prior %s, depth null %s, reference cells %d\n", x$params$c, x$params$a,
              x$params$rho, x$params$rho_mode, x$params$s_sys, x$params$prior, x$params$depth_null, as.integer(x$params$n_ref)))
  neu = s[, mean(cn_a == 1L & cn_b == fifelse(ploidy == 1L, 0L, 1L))]
  cat(sprintf("  MAP neutral: %.1f %% of cell x piece rows; median posterior %.3f\n", 100 * neu, stats::median(s$posterior)))
  invisible(x)
}
