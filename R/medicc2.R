# medicc2.R
# What: export of the tumour-mode clone copy-number states (clone_segments / clones / assignments of a caller run)
#       as MEDICC2 input: common pieces over all clones, one haplotype convention per piece across clones
#       (block orientation by a reference clone's pooled votes), and the MEDICC2 TSV (allele-specific or total).
# Why: every clone's `state` (major:minor) is oriented by its own pseudobulk, so major / minor do not name the same
#      haplotype across clones. MEDICC2's allele-specific mode treats cn_a / cn_b as fixed haplotypes; two clones that
#      lost opposite haplotypes must get swapped cn_a / cn_b, and two that lost the same one identical values.
# Inputs: segments, clones and assignments tables of a tumour run (vignette("scphasedcn-outputs", package = "scphasedcn")), the run's molecule votes
#         and atoms (scn_features). Outputs: data.tables (pieces, orientation QC, block orientation, MEDICC2 table).

#' Common copy-number pieces of the clones of a tumour run (MEDICC2 export, step 1)
#'
#' Drops doublet clones, clones with fewer than `min_assigned` assigned cells and clones without segments; builds the
#' union of the remaining clones' segment breakpoints ([segment_pieces()]) on `chroms`, maps every clone's state to
#' every piece (segment containing the piece midpoint, as `piece_states()`) and keeps the pieces covered by every kept
#' clone (gaps of any clone, e.g. missing windows or a clone-specific chromosome start, are dropped). Pieces shorter
#' than `min_len` are then absorbed into an adjacent piece of the same chromosome (gap <= `max_gap`, by default 0:
#' only a neighbour that touches it): the neighbour whose states differ in the fewest clones, the longer one on ties;
#' the neighbour keeps its states and extends over the short piece. A short piece without an eligible neighbour is
#' dropped. A neighbour across a gap (a dropped range, see `pieces_dropped`, or a range without any segment) is not
#' eligible by default: the absorbing piece would extend over the gap and claim states where not every clone has one.
#' Adjacent pieces with identical states are not merged here: merging needs the common haplotype orientation and is
#' done by [medicc2_table()].
#' Segments labelled `imbalance` (allelic imbalance without depth support, `segment_clone(gain_support = )`) enter as
#' 1:1: MEDICC2 sees no extra copy where the clone's depth does not show one.
#' @param segments segments table (`clone, chrom, start, end, state, cn, label, ploidy`; `nA, nB, f` optional)
#' @param clones clones table (`clone`, `doublet`, `n_assigned`); NULL keeps every clone of `segments`
#' @param min_len pieces shorter than this (bp) are absorbed into a neighbour
#' @param max_gap maximal gap (bp) between a short piece and the neighbour it is absorbed into (default 0: contiguous
#'   neighbours only; a larger value lets the absorbing piece extend over the gap, i.e. claim states there)
#' @param min_assigned clones with fewer assigned cells are dropped (the normal clone `normal` is always kept)
#' @param chroms chromosomes exported (MEDICC2: chr1-22 and chrX)
#' @param normal name of the normal clone
#' @return long data.table, one row per clone x piece: `piece` (`seg|all|chrom:start-end`), `chrom`, `start`, `end`,
#'   `clone`, `state`, `major`, `minor`, `cn`, `label`, `ploidy`, `imbalanced` (major != minor on a diploid piece; these
#'   need orientation), `seg_nA`, `seg_nB`, `seg_f` (the containing segment's pooled votes in the clone's own
#'   orientation, NA when absent). Attributes `clones_dropped` (clone, reason) and `pieces_dropped` (chrom, start, end,
#'   reason: `not_all_clones` / `short_isolated`) and `n_pieces_raw` (pieces before the coverage filter).
#' @export
medicc2_pieces = function(segments, clones = NULL, min_len = 1e6, max_gap = 0, min_assigned = 1,
                          chroms = paste0("chr", c(1:22, "X")), normal = "N0") {
  sg = as.data.table(segments)[chrom %in% chroms]
  if (!nrow(sg)) stop("medicc2_pieces: no segments on ", paste(chroms, collapse = ", "))
  # clone ids as character (numeric ids would select the state matrix columns by position); the clones table is
  # copied, so the caller's table is not modified
  sg[, clone := as.character(clone)]
  dropped = data.table(clone = character(), reason = character())
  if (!is.null(clones)) {
    cl = copy(as.data.table(clones))[, clone := as.character(clone)]
    if (!"doublet" %in% names(cl)) cl[, doublet := FALSE]
    if (!"n_assigned" %in% names(cl)) cl[, n_assigned := NA_real_]
    cl[, reason := fifelse(doublet %in% TRUE, "doublet",
                           fifelse(clone != normal & (is.na(n_assigned) | n_assigned < min_assigned), "no_assigned_cells", NA_character_))]
    cl[is.na(reason) & !clone %in% sg$clone, reason := "no_segments"]
    dropped = cl[!is.na(reason), .(clone, reason)]
    sg = sg[clone %in% cl[is.na(reason), clone]]
  }
  if (!nrow(sg)) stop("medicc2_pieces: no clone left")
  K = sort(unique(sg$clone))
  sg[, `:=`(major = nA_of(state), minor = nB_of(state))]
  if (!"ploidy" %in% names(sg)) sg[, ploidy := 2L]
  for (v in c("nA", "nB", "f")) if (!v %in% names(sg)) sg[, (v) := NA_real_]
  if (!"label" %in% names(sg)) sg[, label := state_label(major, minor, ploidy)]
  if (!"cn" %in% names(sg)) sg[, cn := major + minor]
  # imbalanced runs without depth support (segment_clone(gain_support = )): copy-neutral for the tree (no extra copy)
  if (any(sg$label == "imbalance"))
    sg[label == "imbalance", `:=`(state = "1:1", major = nA_of("1:1"), minor = nB_of("1:1"), cn = 2 + 0 * cn, label = "neutral")]
  pcs = segment_pieces(sg)[, .(region, chrom, start, end)]
  n_raw = nrow(pcs)
  st = sg[pcs[, .(region, chrom, mid = (start + end) / 2)], on = .(chrom, start <= mid, end > mid),
          .(region = i.region, clone = x.clone, state = x.state, major = x.major, minor = x.minor, cn = x.cn, label = x.label,
            ploidy = as.integer(x.ploidy), seg_nA = as.numeric(x.nA), seg_nB = as.numeric(x.nB), seg_f = as.numeric(x.f)),
          nomatch = NULL, allow.cartesian = TRUE]
  cov = st[, .(n = uniqueN(clone)), by = region]
  full = cov[n == length(K), region]
  pdrop = pcs[!region %in% full, .(chrom, start, end, reason = "not_all_clones")]
  pcs = pcs[region %in% full]
  pcs = pcs[order(chrom_rank(chrom), start)]
  # state matrix pieces x clones for the absorption of short pieces
  S = as.matrix(dcast(st[region %in% pcs$region], region ~ clone, value.var = "state"), rownames = "region")[pcs$region, K, drop = FALSE]
  P = pcs[, .(orig = region, chrom, start, end)]
  repeat {
    len = P$end - P$start
    short = which(len < min_len)
    if (!length(short)) break
    ii = short[which.min(len[short])]
    nb = c(ii - 1L, ii + 1L)
    nb = nb[nb >= 1 & nb <= nrow(P)]
    nb = nb[P$chrom[nb] == P$chrom[ii]]
    gap = vapply(nb, function(j) if (j < ii) P$start[ii] - P$end[j] else P$start[j] - P$end[ii], numeric(1))
    nb = nb[gap <= max_gap]
    if (!length(nb)) {
      pdrop = rbind(pdrop, P[ii, .(chrom, start, end, reason = "short_isolated")])
      P = P[-ii]; S = S[-ii, , drop = FALSE]
      next
    }
    d = vapply(nb, function(j) sum(S[j, ] != S[ii, ]), numeric(1))
    jj = nb[order(d, -(P$end[nb] - P$start[nb]))[1]]
    ns = min(P$start[jj], P$start[ii]); ne = max(P$end[jj], P$end[ii])
    set(P, jj, c("start", "end"), list(ns, ne))
    P = P[-ii]; S = S[-ii, , drop = FALSE]
  }
  out = st[P, on = .(region = orig), nomatch = NULL]
  out[, piece := sprintf("seg|all|%s:%s-%s", chrom, format(start, scientific = FALSE, trim = TRUE), format(end, scientific = FALSE, trim = TRUE))]
  out[, imbalanced := major != minor & ploidy != 1L]
  out = out[, .(piece, chrom, start, end, clone, state, major, minor, cn = as.integer(cn), label, ploidy, imbalanced,
                seg_nA, seg_nB, seg_f)]
  out = out[order(chrom_rank(chrom), start, clone)]
  setattr(out, "clones_dropped", dropped)
  setattr(out, "pieces_dropped", pdrop[order(chrom_rank(chrom), start)])
  setattr(out, "n_pieces_raw", n_raw)
  out[]
}
utils::globalVariables(c("doublet", "n_assigned", "reason", "x.major", "x.minor", "x.nA", "x.nB", "x.f", "i.region",
                         "seg_nA", "seg_nB", "seg_f", "major", "minor", "ploidy", "imbalanced", "orig", "i.chrom", "i.start",
                         "i.end", "x.start", "x.end"))

#' Common haplotype orientation of clone states per piece (MEDICC2 export, step 2)
#'
#' Per piece, a reference clone is chosen among the clones with an imbalanced state (major != minor on a diploid
#' piece): the clone with the largest expected orientation signal `a * |major / cn - 0.5| * sqrt(n)`, where `a` is
#' the clone's allelic attenuation (`clones$a`, 1 when missing) and `n` its pooled votes on the piece. Each phase block
#' (piece x phase set) is oriented by the sign of the reference clone's pooled votes (A - B; ties and blocks without
#' reference votes are left unoriented and do not count). Every clone's votes on the oriented blocks are pooled into
#' an oriented fraction `f`: the fraction of the reference clone's major haplotype. The convention is cn_a = the
#' reference clone's major haplotype; a clone with an imbalanced state gets `cn_a = major, cn_b = minor` when
#' `f >= 0.5`, else the swap. Balanced states (n:n, 0:0) and haploid pieces (male chrX, ploidy 1) pass through as
#' `cn_a = major, cn_b = minor`.
#'
#' Chaining (`chain`): which haplotype of a piece is called "a" is arbitrary (only the relation between clones is
#' measured), but MEDICC2 counts an event over adjacent pieces once only when the label is continuous. Walking along
#' each chromosome, a piece's labels are flipped as a whole (cn_a <-> cn_b in every clone, `f` -> 1 - `f`, block signs
#' negated) when, relative to the previous piece, (`parsimony`) more confidently oriented clones (`ok` / `ref`)
#' change direction (sign of cn_a - cn_b) than keep it (clones in the same state in both pieces count 1, others 0.5),
#' or (`phase`) the phase sets shared by the two pieces say that the previous piece's cn_a haplotype is this piece's
#' cn_b (sum over shared blocks of the product of the two reference signs, weighted by min(|nA_ref - nB_ref|); used
#' when the weights sum to >= `min_link` votes). The other rule is the fallback when the first has no evidence.
#' `parsimony` is the default: a single long phase block straddling a piece boundary with a phase switch inside it
#' can link the two halves of one clonal loss with opposite labels. Pieces after a chromosome
#' start or next to a piece without imbalanced clones are not chained. Flipping never changes which clones share a
#' haplotype within a piece.
#'
#' Cells pooled per clone (`use`): `confident` (default) takes the cells confidently assigned to the clone (posterior at least
#' `min_post`, not doublet-like), close to the pseudobulk members that produced the clone states, without the
#' mixed / uncertain cells whose votes would pull `f` towards 0.5; `assigned` takes the `assigned` column. The
#' reference clone orients its blocks with its own votes, so its own `f` is biased upwards (sign reliable, magnitude
#' not; `is_ref` marks these rows). Orientation errors of the reference on blocks with few votes attenuate the other
#' clones' `f` towards 0.5 but do not flip it; a clone whose oriented `f` is within `near` of 0.5 despite an imbalanced
#' state is flagged `low_conf` (its cn_a / cn_b assignment is a guess).
#' @param pieces [medicc2_pieces()] output
#' @param features list with `votes` (cell index, atom, ps, nA, nB), `atoms` and `cells` (the cells the vote index
#'   refers to: data.table with a `cell` column or a character vector), e.g. an `scn_features` object
#' @param assignments assignments table (`cell`, `clone` and the `use` column)
#' @param clones clones table (for the attenuation `a`); NULL = 1
#' @param use cells pooled per clone: `confident` or `assigned`
#' @param near half-width around 0.5 of the low-confidence flag
#' @param chain label chaining between adjacent pieces: `parsimony` (phase-set fallback), `phase` (parsimony
#'   fallback) or `none` (each piece's cn_a is its reference clone's major haplotype)
#' @param min_link minimal summed block weight (votes) of a phase-set link
#' @return list of class `medicc2_orientation`:
#'   - `pieces`: the input rows plus `n_votes` (clone votes on the piece), `n_oriented` (on oriented blocks), `f`
#'     (oriented fraction of the cn_a haplotype), `dev` (|f - 0.5|), `is_ref`, `swapped` (cn_a = minor), `cn_a`, `cn_b`
#'     and `flag` (`ref`, `ok`, `low_conf`, `no_votes`, `no_ref`, `balanced`, `haploid`).
#'   - `qc`: one row per piece: `piece, chrom, start, end, orient_ref, ref_state, ref_score, flipped` (labels flipped
#'     by chaining), `link` (`first`, `phase`, `parsimony`, `none`), `link_score` (phase: weighted sign agreement between
#'     -1 and 1 before the flip; parsimony: the weighted direction agreement), `n_imbalanced, n_blocks`
#'     (phase blocks with votes of any kept clone), `n_blocks_oriented`, `votes` (all kept clones), `votes_oriented`,
#'     `ref_votes`, `n_swapped`, `n_low_conf`.
#'   - `blocks`: the reference orientation per block, reusable by a per-cell module: one row per piece x phase set
#'     with reference votes: `piece, chrom, start, end` (the piece; a molecule belongs to it when its atom start is in
#'     `[start, end)`), `ps` (phase set of the votes table), `ref` (reference clone), `nA_ref, nB_ref` and `sign`
#'     (+1: phase-set haplotype A is cn_a; -1: haplotype B is cn_a; 0: tie, unoriented). Blocks absent from the table
#'     are unoriented. A cell's oriented votes on a piece are `a = sum(sign > 0 ? nA : nB)`, `b = sum(sign > 0 ? nB :
#'     nA)` over its blocks with `sign != 0`; `a / (a + b)` estimates `cn_a / (cn_a + cn_b)` (attenuated).
#'   - `params`: `use`, `near`, `chain`, `min_link`.
#' @export
orient_clone_states = function(pieces, features, assignments, clones = NULL, use = c("confident", "assigned"), near = 0.05,
                               chain = c("parsimony", "phase", "none"), min_link = 50) {
  use = match.arg(use); chain = match.arg(chain)
  pl = copy(as.data.table(pieces))[, clone := as.character(clone)]
  votes = features$votes; atoms = as.data.table(features$atoms)
  if (is.null(votes) || is.null(atoms)) stop("orient_clone_states: `features` needs votes and atoms")
  cells = features$cells
  cells = if (is.data.frame(cells)) cells$cell else as.character(cells)
  if (!length(cells)) stop("orient_clone_states: `features` needs the cells the votes refer to")
  K = sort(unique(pl$clone))
  pc = unique(pl[, .(piece, chrom, start, end)])[order(chrom_rank(chrom), start)]
  pc[, pix := seq_len(.N)]
  pl[pc, pix := i.pix, on = "piece"]
  pl[, cl := match(clone, K)]
  # cell index -> clone index (cells pooled per clone)
  asg = copy(as.data.table(assignments))
  if (!use %in% names(asg)) stop("orient_clone_states: assignments lack `", use, "`")
  asg[, clone := as.character(clone)]
  asg = asg[get(use) %in% TRUE & clone %in% K]
  ci_cl = rep(NA_integer_, length(cells))
  k = match(asg$cell, cells)
  ci_cl[k[!is.na(k)]] = match(asg$clone[!is.na(k)], K)
  # atom -> piece index (molecule rule: atom start inside [start, end))
  am = atoms[, .(atom, chrom, pos = start)][pc, on = .(chrom, pos >= start, pos < end), .(atom = x.atom, pix = i.pix), nomatch = NULL]
  at_pi = rep(NA_integer_, max(atoms$atom))
  at_pi[am$atom] = am$pix
  cv = ci_cl[votes$cell]; pv = at_pi[votes$atom]
  ok = which(!is.na(cv) & !is.na(pv))
  agg = data.table(pix = pv[ok], ps = votes$ps[ok], cl = cv[ok], nA = as.numeric(votes$nA[ok]), nB = as.numeric(votes$nB[ok]))
  rm(cv, pv, ok)
  agg = agg[, .(nA = sum(nA), nB = sum(nB)), by = .(pix, ps, cl)]
  tot = agg[, .(n_votes = sum(nA + nB)), by = .(pix, cl)]
  pl[tot, n_votes := i.n_votes, on = .(pix, cl)]
  pl[is.na(n_votes), n_votes := 0]
  # reference clone per piece
  att = rep(1, length(K))
  if (!is.null(clones) && "a" %in% names(clones)) {
    ca = as.data.table(clones)[, .(clone = as.character(clone), a = as.numeric(a))]
    att = ca$a[match(K, ca$clone)]
    att[is.na(att)] = 1
  }
  pl[, ref_score := fifelse(imbalanced & n_votes > 0, att[cl] * abs(major / pmax(cn, 1) - 0.5) * sqrt(n_votes), NA_real_)]
  ref = pl[!is.na(ref_score)][order(pix, -ref_score, clone)][!duplicated(pix), .(pix, ref_cl = cl, ref_score, ref_state = state)]
  # block orientation by the reference clone's pooled votes
  blk = agg[ref, on = .(pix, cl = ref_cl), nomatch = NULL][, .(pix, ps, ref = K[cl], nA_ref = nA, nB_ref = nB, sign = as.integer(sign(nA - nB)))]
  ov = blk[sign != 0][agg, on = .(pix, ps), nomatch = NULL]
  ov = ov[, .(n_oriented = sum(nA + nB), oA = sum(fifelse(sign > 0, nA, nB))), by = .(pix, cl)]
  pl[ov, `:=`(n_oriented = i.n_oriented, oA = i.oA), on = .(pix, cl)]
  pl[is.na(n_oriented), `:=`(n_oriented = 0, oA = 0)]
  pl[, f := fifelse(n_oriented > 0, oA / n_oriented, NA_real_)]
  pl[, dev := abs(f - 0.5)]
  pl[ref, `:=`(ref_cl = i.ref_cl), on = "pix"]
  pl[, is_ref := !is.na(ref_cl) & cl == ref_cl]
  pl[, swapped := imbalanced & !is.na(ref_cl) & !is.na(f) & f < 0.5]
  pl[, `:=`(cn_a = fifelse(swapped, minor, major), cn_b = fifelse(swapped, major, minor))]
  pl[, flag := fifelse(ploidy == 1L, "haploid",
               fifelse(!imbalanced, "balanced",
               fifelse(is.na(ref_cl), "no_ref",
               fifelse(is_ref, "ref",
               fifelse(is.na(f), "no_votes",
               fifelse(dev < near, "low_conf", "ok"))))))]
  # chaining: the per-piece convention (cn_a = the reference's major) is arbitrary; flip whole pieces so that cn_a is
  # the same haplotype as in the previous piece of the chromosome (phase-set link, else fewest direction changes)
  n_p = nrow(pc)
  flipped = logical(n_p); link = rep("first", n_p); link_score = rep(NA_real_, n_p)
  if (chain != "none" && n_p > 1) {
    pl = pl[order(pix, cl)]
    D = matrix(pl[, sign(cn_a - cn_b)], n_p, length(K), byrow = TRUE)
    W = matrix(pl[, as.numeric(flag %in% c("ok", "ref"))], n_p, length(K), byrow = TRUE)
    ST = matrix(pl$state, n_p, length(K), byrow = TRUE)
    if (nrow(pl) != n_p * length(K) || any(pl$cl != rep(seq_along(K), n_p))) stop("orient_clone_states: every clone needs a row on every piece")
    has_ref = seq_len(n_p) %in% ref$pix
    bl = lapply(seq_len(n_p), function(k) blk[pix == k & sign != 0, .(ps, s = sign, w = abs(nA_ref - nB_ref))])
    for (k in 2:n_p) {
      if (pc$chrom[k] != pc$chrom[k - 1]) next
      link[k] = "none"
      if (!has_ref[k] || !has_ref[k - 1]) next
      by_phase = function() {
        m = merge(bl[[k - 1]], bl[[k]], by = "ps")
        w = pmin(m$w.x, m$w.y)
        if (!nrow(m) || sum(w) < min_link) return(NULL)
        ph = sum(m$s.x * m$s.y * w)
        if (ph == 0) return(NULL)
        list(flip = ph < 0, link = "phase", score = ph / sum(w))
      }
      by_parsimony = function() {
        S = sum(W[k - 1, ] * W[k, ] * ifelse(ST[k - 1, ] == ST[k, ], 1, 0.5) * D[k - 1, ] * D[k, ])
        if (S == 0) return(NULL)
        list(flip = S < 0, link = "parsimony", score = S)
      }
      dec = if (chain == "phase") by_phase() %||% by_parsimony() else by_parsimony() %||% by_phase()
      flip = FALSE
      if (!is.null(dec)) { flip = dec$flip; link[k] = dec$link; link_score[k] = dec$score }
      if (isTRUE(flip)) {
        flipped[k] = TRUE; D[k, ] = -D[k, ]
        bl[[k]][, s := -s]
      }
    }
    fp = which(flipped)
    if (length(fp)) {
      pl[pix %in% fp, `:=`(cn_a = cn_b, cn_b = cn_a, f = 1 - f)]
      blk[pix %in% fp, sign := -sign]
    }
    pl[, swapped := imbalanced & cn_a < cn_b]
  }
  qc_b = agg[, .(n_blocks = uniqueN(ps), votes = sum(nA + nB)), by = pix]
  qc_o = blk[sign != 0, .(n_blocks_oriented = .N, ref_votes = sum(nA_ref + nB_ref)), by = pix]
  qc_v = pl[, .(n_imbalanced = sum(imbalanced), votes_oriented = sum(n_oriented), n_swapped = sum(swapped),
                n_low_conf = sum(flag == "low_conf")), by = pix]
  qc = pc[, .(pix, piece, chrom, start, end, flipped = flipped, link = link, link_score = link_score)]
  qc[ref, `:=`(orient_ref = K[i.ref_cl], ref_state = i.ref_state, ref_score = i.ref_score), on = "pix"]
  qc[qc_v, `:=`(n_imbalanced = i.n_imbalanced, votes_oriented = i.votes_oriented, n_swapped = i.n_swapped,
                n_low_conf = i.n_low_conf), on = "pix"]
  qc[qc_b, `:=`(n_blocks = i.n_blocks, votes = i.votes), on = "pix"]
  qc[qc_o, `:=`(n_blocks_oriented = i.n_blocks_oriented, ref_votes = i.ref_votes), on = "pix"]
  for (v in c("n_blocks", "votes", "n_blocks_oriented", "ref_votes")) qc[is.na(get(v)), (v) := 0]
  setcolorder(qc, c("pix", "piece", "chrom", "start", "end", "orient_ref", "ref_state", "ref_score", "flipped", "link", "link_score",
                    "n_imbalanced", "n_blocks", "n_blocks_oriented", "votes", "votes_oriented", "ref_votes", "n_swapped", "n_low_conf"))
  blocks = blk[pc, on = "pix", nomatch = NULL][, .(piece, chrom, start, end, ps, ref, nA_ref, nB_ref, sign)]
  blocks = blocks[order(chrom_rank(chrom), start, ps)]
  pl = pl[order(pix, clone)]
  pl[, c("pix", "cl", "ref_cl", "oA") := NULL]
  qc[, pix := NULL]
  structure(list(pieces = pl[], qc = qc[], blocks = blocks[], params = list(use = use, near = near, chain = chain, min_link = min_link)), class = "medicc2_orientation")
}
utils::globalVariables(c("pix", "cl", "i.pix", "x.atom", "n_votes", "i.n_votes", "ref_score", "ref_cl", "ref_state", "nA_ref",
                         "nB_ref", "n_oriented", "oA", "i.n_oriented", "i.oA", "dev", "i.ref_cl", "is_ref", "swapped", "cn_a",
                         "cn_b", "i.ref_state", "i.ref_score", "i.n_imbalanced", "i.votes_oriented", "i.n_swapped",
                         "i.n_low_conf", "i.n_blocks", "i.votes", "i.n_blocks_oriented", "i.ref_votes", "n_blocks_oriented",
                         "ref_votes", "votes_oriented", "n_swapped", "n_low_conf", "n_imbalanced", "votes", "orient_ref", "n_blocks",
                         "s", "flipped", "link", "link_score"))

#' MEDICC2 input table (MEDICC2 export, step 3)
#'
#' One row per sample x segment with integer copy numbers capped at `cap`: `sample_id, chrom, start, end, cn_a, cn_b`
#' (allele-specific) or `sample_id, chrom, start, end, cn` (total = cn_a + cn_b, for `--total-copy-numbers
#' --input-allele-columns cn`). Adjacent pieces of a chromosome whose values are identical in every sample are merged
#' (after the cap; on cn_a / cn_b, or on cn for the total table), so every sample has the same segments, as MEDICC2
#' requires; pieces separated by a gap (e.g. a range dropped by [medicc2_pieces()]) are never merged, so no segment
#' spans a range without states. The normal clone `normal` is a sample (MEDICC2 `--normal-name N0`); when absent it is added as 1:1 on
#' diploid and 1:0 on haploid pieces (male chrX), so it is the root of the tree. Chromosomes are sorted naturally.
#' @param x [orient_clone_states()] result, or a long table with `clone, chrom, start, end, cn_a, cn_b` (+ `ploidy`)
#' @param allele_specific TRUE: cn_a / cn_b columns; FALSE: total cn
#' @param normal name of the normal sample
#' @param cap maximal copy number (MEDICC2 caps at 8)
#' @param chroms chromosomes kept
#' @return data.table (MEDICC2 TSV layout); attribute `n_pieces` (pieces before merging)
#' @export
medicc2_table = function(x, allele_specific = TRUE, normal = "N0", cap = 8L, chroms = paste0("chr", c(1:22, "X"))) {
  d = as.data.table(if (inherits(x, "medicc2_orientation")) x$pieces else x)
  if (!all(c("clone", "chrom", "start", "end", "cn_a", "cn_b") %in% names(d))) stop("medicc2_table: need clone, chrom, start, end, cn_a, cn_b")
  d = d[chrom %in% chroms, .(clone = as.character(clone), chrom, start, end, cn_a, cn_b,
                              ploidy = if ("ploidy" %in% names(d)) as.integer(ploidy) else 2L)]
  if (!normal %in% d$clone) {
    n0 = d[, .(ploidy = min(ploidy)), by = .(chrom, start, end)]
    n0 = n0[, .(clone = normal, chrom, start, end, cn_a = 1L, cn_b = fifelse(ploidy == 1L, 0L, 1L), ploidy)]
    d = rbind(d, n0)
  }
  d[, `:=`(cn_a = pmin(as.integer(round(cn_a)), as.integer(cap)), cn_b = pmin(as.integer(round(cn_b)), as.integer(cap)))]
  d[, cn := pmin(cn_a + cn_b, as.integer(cap))]
  pcs = unique(d[, .(chrom, start, end)])[order(chrom_rank(chrom), start)]
  n_pc = d[, .N, by = .(chrom, start, end)]
  if (any(n_pc$N != uniqueN(d$clone))) stop("medicc2_table: every sample needs a value on every piece")
  val = if (allele_specific) d[, paste(cn_a, cn_b)] else as.character(d$cn)
  d[, v := val]
  kk = d[order(clone), .(kk = paste(v, collapse = "|")), by = .(chrom, start, end)]
  pcs[kk, kk := i.kk, on = .(chrom, start, end)]
  # a run never spans a gap between pieces (a range without states, e.g. dropped by medicc2_pieces())
  pcs[, gap_ := cumsum(c(TRUE, start[-1] != end[-.N])), by = chrom]
  pcs[, run := rleid(chrom, kk, gap_)]
  runs = pcs[, .(start = min(start), end = max(end), chrom = chrom[1]), by = run]
  d[pcs, run := i.run, on = .(chrom, start, end)]
  out = d[, .(cn_a = cn_a[1], cn_b = cn_b[1], cn = cn[1]), by = .(sample_id = clone, run)]
  out[runs, `:=`(chrom = i.chrom, start = i.start, end = i.end), on = "run"]
  smp = c(normal, setdiff(sort(unique(out$sample_id)), normal))
  out = out[order(match(sample_id, smp), chrom_rank(chrom), start)]
  out[, `:=`(start = as.integer(start), end = as.integer(end))]
  out = if (allele_specific) out[, .(sample_id, chrom, start, end, cn_a, cn_b)] else out[, .(sample_id, chrom, start, end, cn)]
  setattr(out, "n_pieces", nrow(pcs))
  out[]
}
utils::globalVariables(c("v", "run", "kk", "i.kk", "i.run", "sample_id", "gap_"))
