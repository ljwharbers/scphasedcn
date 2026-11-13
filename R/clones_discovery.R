# clones_discovery.R
# What: tumour mode of the copy-number caller, part 1 of 4: clone discovery. Small shared
#       utilities (clone agreement indices, input accessors, phase-block orientation, discovery / test
#       halves), the signed allelic + depth clone feature matrix (kNN-smoothed PCA, Ward tree over-cut), the
#       clone contrast and its permutation-calibrated 1 / 1,000 split test, sibling merges, within-trunk
#       discovery and the opt-in disagreement split pass.
# Why: in tumours per-cell detection is reliable but per-cell states are not; clones are found first, and
#      their states settle the calls.
# Inputs: an `scn_features`-like list (see [caller_tumour()]) or the individual tables
#         (features, regions, votes, atoms; vignette("scphasedcn-outputs", package = "scphasedcn")).
# Outputs: clone feature matrices, leaves / merged clones, split-test statistics.
# See also: clones_states.R (profiles, states, segments), clones_assign.R (pieces, assignment, doublets,
#           reference cleaning), clones_tumour.R (the tumour-mode caller and its options).

# ---- small utilities -----------------------------------------------------------------------------

#' Adjusted Rand index of two labelings
#'
#' Hubert-Arabie ARI on the cells labelled in both vectors (NA dropped).
#' @param x,y label vectors of equal length
#' @return numeric ARI (1 = identical partitions)
#' @export
clone_ari = function(x, y) {
  ok = !is.na(x) & !is.na(y)
  x = as.character(x[ok]); y = as.character(y[ok])
  if (length(x) < 2) return(NA_real_)
  tab = table(x, y)
  n = sum(tab)
  sc = sum(choose(tab, 2)); a = sum(choose(rowSums(tab), 2)); b = sum(choose(colSums(tab), 2))
  ex = a * b / choose(n, 2); mx = (a + b) / 2
  if (isTRUE(all.equal(mx, ex))) return(1)
  (sc - ex) / (mx - ex)
}

#' Adjusted mutual information of two labelings
#'
#' AMI (Vinh, Epps & Bailey 2010; arithmetic-mean normalisation, scikit-learn's default) on the cells labelled
#' in both vectors (NA and "" dropped): `(MI - E(MI)) / (mean(H(x), H(y)) - E(MI))`, E(MI) under the
#' hypergeometric (permutation) model.
#' @param x,y label vectors of equal length
#' @return numeric AMI (1 = identical partitions, ~0 = independent)
#' @export
clone_ami = function(x, y) {
  ok = !is.na(x) & !is.na(y) & x != "" & y != ""
  x = as.character(x[ok]); y = as.character(y[ok])
  n = length(x)
  if (n < 2) return(NA_real_)
  tab = table(x, y)
  a = rowSums(tab); b = colSums(tab)
  if (length(a) == 1 && length(b) == 1) return(1)
  h = function(p) { p = p[p > 0] / sum(p); -sum(p * log(p)) }
  ii = which(tab > 0, arr.ind = TRUE)
  nz = tab[ii]
  mi = sum(nz / n * log(n * nz / (a[ii[, 1]] * b[ii[, 2]])))
  emi = 0
  for (i in seq_along(a)) for (j in seq_along(b)) {
    lo = max(1, a[i] + b[j] - n); hi = min(a[i], b[j])
    if (lo > hi) next
    nij = lo:hi
    lp = lfactorial(a[i]) + lfactorial(b[j]) + lfactorial(n - a[i]) + lfactorial(n - b[j]) - lfactorial(n) -
      lfactorial(nij) - lfactorial(a[i] - nij) - lfactorial(b[j] - nij) - lfactorial(n - a[i] - b[j] + nij)
    emi = emi + sum(nij / n * log(n * nij / (a[i] * b[j])) * exp(lp))
  }
  den = (h(a) + h(b)) / 2 - emi
  if (abs(den) < 1e-12) return(if (abs(mi - emi) < 1e-12) 1 else 0)
  unname((mi - emi) / den)
}

# fetch a component of the scn_features list under one of several names
cl_part = function(feat, names, required = TRUE) {
  for (nm in names) if (!is.null(feat[[nm]])) return(feat[[nm]])
  if (required) stop(sprintf("caller_tumour: input lacks `%s`", names[1]), call. = FALSE)
  NULL
}

# normalised view of the scn_features input used by the clone functions
clone_input = function(feat) {
  if (inherits(feat, "clone_input")) return(feat)
  features = as.data.table(cl_part(feat, c("features", "long")))
  regions = as.data.table(cl_part(feat, "regions"))
  cells = cl_part(feat, "cells", required = FALSE)
  if (is.null(cells)) cells = data.table(cell = unique(features$cell))
  if (is.character(cells)) cells = data.table(cell = cells)
  cells = as.data.table(cells)
  nf = cl_part(feat, c("null_fit"), required = FALSE)
  tpar = feat$tpar %||% nf$tpar %||% nf$depth
  rho = feat$rho %||% nf$rho %||% nf$allelic
  if (is.data.frame(rho)) rho = stats::median(rho$rho, na.rm = TRUE)
  if (is.null(rho) || !is.finite(rho)) {
    nr = features[null_cell & null_region & eligible]
    rho = if (nrow(nr)) fit_allelic_null(nr$nA, nr$n_hap) else 0.01
  }
  win = cl_part(feat, c("win", "windows_1mb", "win1mb"), required = FALSE)
  if (!is.null(win)) win = as_clone_win(win, cells$cell)
  haploid = character()
  ao = regions[class %in% c("chrom", "arm") & allelic_ok == FALSE, unique(chrom)]
  if (length(ao)) haploid = intersect(ao, c("chrX", "chrY"))
  structure(list(features = features, regions = regions, votes = as.data.table(cl_part(feat, "votes")),
                 atoms = as.data.table(cl_part(feat, "atoms")), cells = cells, win = win, tpar = tpar, rho = rho,
                 haploid = haploid, allelic_chroms = regions[allelic_ok == TRUE, unique(chrom)]),
            class = "clone_input")
}

# 1 Mb window depth as list(windows, x, e): x / e windows x cells (colnames = cells). Accepts that list or a
# long table cell, chrom, start, end, x, e.
as_clone_win = function(win, cells) {
  if (is.list(win) && !is.data.frame(win) && all(c("windows", "x", "e") %in% names(win))) {
    win$windows = as.data.table(win$windows)
    if (!"window" %in% names(win$windows)) win$windows[, window := paste0(chrom, ":", format(start, scientific = FALSE, trim = TRUE))]
    return(win)
  }
  d = as.data.table(win)
  if (!all(c("cell", "chrom", "start", "end", "x", "e") %in% names(d))) stop("win: need list(windows, x, e) or a long table cell, chrom, start, end, x, e")
  d[, window := paste0(chrom, ":", format(start, scientific = FALSE, trim = TRUE))]
  wt = unique(d[, .(window, chrom, start, end)])
  setorder(wt, chrom, start)
  wt = wt[order(chrom_rank(chrom), start)]
  cc = intersect(cells, unique(d$cell))
  x = matrix(0, nrow(wt), length(cc), dimnames = list(wt$window, cc)); e = x
  ii = cbind(match(d$window, wt$window), match(d$cell, cc))
  ok = !is.na(ii[, 2])
  x[ii[ok, ]] = d$x[ok]; e[ii[ok, ]] = d$e[ok]
  list(windows = wt, x = x, e = e)
}

# votes summed per cell x region x phase block with an integer block id (region x ps): one pass over all regions
# (overlapping classes at once) on the the feature step atom -> region map (atom_region_map(), start-inside rule); the
# block ids follow the votes order, which fixes the column order of the contrast / trip vote matrices
region_block_votes = function(votes, atoms, regions) {
  rg = as.data.table(regions)
  m = atom_region_map(atoms, rg)[, .(atom, region = rg$region[ri])]
  v = as.data.table(votes)[, .(cell, atom, ps, nA, nB)]
  v = m[v, on = "atom", nomatch = NULL, allow.cartesian = TRUE]
  v = v[, .(nA = sum(nA), nB = sum(nB)), by = .(cell, region, ps)]
  v[, blk := .GRP, by = .(region, ps)]
  v[]
}

# group label per cell index (length n_cells; NA = no group)
group_vector = function(groups, n_cells, cells = NULL) {
  if (is.null(groups)) return(rep("all", n_cells))
  if (is.data.frame(groups)) {
    g = rep(NA_character_, n_cells)
    idx = if (is.numeric(groups$cell)) as.integer(groups$cell) else match(groups$cell, cells)
    g[idx[!is.na(idx)]] = as.character(groups[[setdiff(names(groups), "cell")[1]]])[!is.na(idx)]
    return(g)
  }
  if (!is.null(names(groups)) && !is.null(cells)) {
    g = rep(NA_character_, n_cells)
    k = match(names(groups), cells)
    g[k[!is.na(k)]] = as.character(groups)[!is.na(k)]
    return(g)
  }
  if (length(groups) != n_cells) stop("groups: need one label per cell index")
  as.character(groups)
}

# the feature step's orient_blocks() as integer cell index / numeric counts (the clone steps' vote orientation)
cl_orient = function(votes, regions, atoms, groups = NULL, loo = TRUE) {
  out = as.data.table(orient_blocks(votes, regions, atoms, groups, loo = loo))
  out[, .(cell = as.integer(cell), region, nA = as.numeric(nA), nB = as.numeric(nB), n_blocks = as.numeric(n_blocks))]
}

#' Split molecule votes into two independent halves
#'
#' Binomial thinning of each cell x atom x phase-set count (`nA` and `nB` separately, probability `frac`):
#' molecule votes are independent, so the halves are independent samples of the same cells. Clone discovery
#' clusters cells on one half and [merge_clones()] tests the leaves on the other; testing leaves on the data
#' that formed them inflates every contrast (leaves are selected to differ), which a permutation null of random
#' splits does not reproduce.
#' @param votes votes table
#' @param frac fraction of the votes in the first half
#' @param seed RNG seed
#' @return list(disc, test): votes tables (rows with no votes dropped)
#' @export
split_votes = function(votes, frac = 0.5, seed = 1) {
  set.seed(seed)
  v = as.data.table(votes)
  a1 = stats::rbinom(nrow(v), v$nA, frac); b1 = stats::rbinom(nrow(v), v$nB, frac)
  d = v[, .(cell, atom, ps, nA = as.integer(a1), nB = as.integer(b1))]
  t = v[, .(cell, atom, ps, nA = as.integer(nA - a1), nB = as.integer(nB - b1))]
  list(disc = d[nA + nB > 0], test = t[nA + nB > 0])
}

#' Split window depth into two independent halves
#'
#' Binomial thinning of the 1 Mb window counts (`x`; probability `frac`) with the expected counts scaled by `frac`
#' and `1 - frac`: the depth analogue of [split_votes()]. Clone discovery uses one half and the merge / split tests
#' the other, so a depth difference that formed two leaves is not tested on the same counts.
#' @param win 1 Mb window depth list(windows, x, e)
#' @param frac fraction of the counts in the first half
#' @param seed RNG seed
#' @return list(disc, test): window depth lists
#' @export
split_win = function(win, frac = 0.5, seed = 1) {
  set.seed(seed)
  x = round(win$x)
  x1 = matrix(stats::rbinom(length(x), as.vector(x), frac), nrow(x), dimnames = dimnames(x))
  list(disc = list(windows = win$windows, x = x1, e = win$e * frac),
       test = list(windows = win$windows, x = x - x1, e = win$e * (1 - frac)))
}

# per-cell depth features (the long-table columns used by the clone steps) of the eligible regions of `classes` from a
# window depth list (e.g. one half of split_win()); allelic columns 0 (oriented separately)
half_depth_features = function(inp, win, classes) {
  rg = inp$regions[class %in% classes & eligible == TRUE]
  inp2 = inp; inp2$win = win
  d = region_features(inp2, rg)
  d[, `:=`(eligible = TRUE, nA = 0, nB = 0, n_hap = 0, orient = "none")]
  d[]
}

# ---- clone feature matrix ------------------------------------------------------------------------

#' Cells x features matrix for clone discovery
#'
#' One depth column (`d|<region>`: the window z, capped at +/- `cap`, NA -> 0) and one signed allelic
#' column (`a|<region>`) per eligible region of the given classes. Allelic z =
#' `(nA - n / 2) / sqrt(n / 4 * (1 + (n - 1) * rho))` on the oriented votes, capped at +/- `cap`, 0 below
#' `min_n` molecules; the sign keeps which haplotype is in excess, so subclones that are imbalanced on
#' opposite haplotypes separate (e.g. two subclones with a loss of different chr21 homologues). Input whose `orient` is not in `keep_orient`
#' (unoriented block-phased `none`, or `clone` = oriented by known groups such as spatial regions, which would leak the
#' grouping into discovery) is first re-oriented by the all-cell leave-one-out pseudobulk (`orient_blocks`,
#' NULL groups) when `votes`, `atoms`, `regions` and `cells` are given.
#' @param feat_long features long table ([cell_features()] `long`)
#' @param classes region classes used
#' @param rho beta-binomial overdispersion of the molecule votes
#' @param cap cap of |z|
#' @param min_n minimal molecules for an allelic value
#' @param votes,atoms,regions,cells optional: orient unoriented input (cells = names in vote index order)
#' @param keep_orient `orient` labels used as given (chromosome-wide anchored phase, all-cell pseudobulk)
#' @return numeric matrix cells x features (rownames = cells)
#' @export
clone_feature_matrix = function(feat_long, classes = c("arm", "20Mb"), rho, cap = 6, min_n = 5,
                                votes = NULL, atoms = NULL, regions = NULL, cells = NULL, keep_orient = c("anchored", "pseudobulk")) {
  f = as.data.table(feat_long)[class %in% classes & eligible == TRUE,
                               .(cell, region, z, nA, nB, orient = if ("orient" %in% names(feat_long)) orient else "pseudobulk")]
  if (!nrow(f)) stop("clone_feature_matrix: no eligible rows of classes ", paste(classes, collapse = ", "))
  if (any(!f$orient %in% keep_orient) && !is.null(votes)) {
    if (is.null(cells)) stop("clone_feature_matrix: `cells` (names in vote index order) needed to orient")
    rg = as.data.table(regions)[region %in% f$region]
    ob = cl_orient(votes, rg, atoms, NULL, loo = TRUE)
    ob[, cell := cells[cell]]
    f[, c("nA", "nB") := NULL]
    f = merge(f, ob[, .(cell, region, nA, nB)], by = c("cell", "region"), all.x = TRUE)
    f[is.na(nA), `:=`(nA = 0, nB = 0)]
  }
  f[, n := nA + nB]
  f[, az := fifelse(n >= min_n, (nA - n / 2) / sqrt(n / 4 * (1 + (n - 1) * rho)), 0)]
  f[, az := pmin(pmax(az, -cap), cap)]
  f[, dz := pmin(pmax(fifelse(is.finite(z), z, 0), -cap), cap)]
  cc = unique(f$cell); rr = sort(unique(f$region))
  ii = cbind(match(f$cell, cc), match(f$region, rr))
  D = matrix(0, length(cc), length(rr)); A = D
  D[ii] = f$dz; A[ii] = f$az
  colnames(D) = paste0("d|", rr); colnames(A) = paste0("a|", rr)
  X = cbind(D, A)
  rownames(X) = cc
  X[, colSums(X != 0) > 0, drop = FALSE]
}
utils::globalVariables(c("az", "dz", "den", "n_a", "n_b"))

# PCA (irlba) of the centred matrix + kNN smoothing (mean over each cell and its k neighbours)
smooth_embedding = function(X, n_pc = 20, k = 20, seed = 1) {
  set.seed(seed)
  Xc = sweep(X, 2, colMeans(X), "-")
  keep = apply(Xc, 2, stats::var) > 0
  Xc = Xc[, keep, drop = FALSE]
  n_pc = max(1L, min(n_pc, ncol(Xc) - 1L, nrow(Xc) - 2L))
  full_svd = function() { s0 = svd(Xc, nu = n_pc, nv = 0); list(u = s0$u, d = s0$d[seq_len(n_pc)]) }
  # irlba only for a small share of the singular values of a large matrix (it warns and may not converge otherwise;
  # then the exact SVD is used)
  sv = if (n_pc >= 0.2 * min(dim(Xc)) || min(dim(Xc)) < 100) full_svd() else
    tryCatch(irlba::irlba(Xc, nv = n_pc), warning = function(w) full_svd())
  pc = sv$u %*% diag(sv$d, n_pc, n_pc)
  rownames(pc) = rownames(X)
  k = min(k, nrow(pc) - 1L)
  if (k < 1) return(list(pc = pc, smooth = pc))
  nn = FNN::get.knn(pc, k = k)$nn.index
  S = pc
  for (j in seq_len(k)) S = S + pc[nn[, j], , drop = FALSE]
  S = S / (k + 1)
  rownames(S) = rownames(X)
  list(pc = pc, smooth = S)
}

# reassign the cells of clusters smaller than min_size to the nearest remaining centroid
absorb_small = function(lab, S, min_size) {
  repeat {
    sz = table(lab)
    if (length(sz) <= 1 || min(sz) >= min_size) break
    small = names(sz)[which.min(sz)]
    cent = rowsum(S, lab) / as.vector(table(lab)[rownames(rowsum(S, lab))])
    others = setdiff(rownames(cent), small)
    idx = which(lab == small)
    d = as.matrix(FNN::get.knnx(cent[others, , drop = FALSE], S[idx, , drop = FALSE], k = 1)$nn.index)
    lab[idx] = others[d[, 1]]
  }
  lab
}

#' Discover clone leaves: smoothed PCA + Ward tree over-cut
#'
#' `irlba` PCA of the centred feature matrix, kNN smoothing of the PCs (mean over each cell and its `k`
#' neighbours, FNN), Ward.D2 `hclust` on the smoothed PCs, cut at `max_leaves`; leaves smaller than
#' `min_size` are merged into the nearest leaf (centroid distance). The leaves over-split on purpose;
#' [merge_clones()] merges them back with a calibrated test.
#' @param X cells x features matrix ([clone_feature_matrix()])
#' @param n_pc principal components
#' @param k neighbours for the smoothing
#' @param max_leaves leaves of the over-cut
#' @param min_size minimal leaf size (cells)
#' @param seed RNG seed
#' @return list(tree = hclust, leaves = named integer labels (names = cells), pc, smooth)
#' @export
discover_clones = function(X, n_pc = 20, k = 20, max_leaves = 40, min_size = 20, seed = 1) {
  em = smooth_embedding(X, n_pc, k, seed)
  tree = stats::hclust(stats::dist(em$smooth), method = "ward.D2")
  lab = as.character(stats::cutree(tree, k = min(max_leaves, nrow(X))))
  lab = absorb_small(lab, em$smooth, min_size)
  lab = as.integer(factor(lab, levels = unique(lab)))
  names(lab) = rownames(X)
  list(tree = tree, leaves = lab, pc = em$pc, smooth = em$smooth)
}

# ---- clone contrast and the permutation split test ---------------------------------------------------

#' Data for repeated clone contrasts
#'
#' Pre-computes what the clone contrast (`contrast_stat`) needs for many splits: the per-cell log2 ratio matrix (cells x
#' regions) and the molecule votes per cell x region x phase block as sparse cells x blocks matrices.
#' @param feat_long features long table
#' @param votes,atoms,regions votes, atoms and regions tables
#' @param cells cell names in vote index order (rows of the result)
#' @param classes region classes used
#' @param allelic_chroms chromosomes with an allele channel (default all)
#' @return list(L, regions, MA, MB, blk_region, N, cells)
#' @export
contrast_data = function(feat_long, votes, atoms, regions, cells, classes = c("arm", "20Mb"), allelic_chroms = NULL) {
  f = as.data.table(feat_long)[class %in% classes & eligible == TRUE & cell %in% cells, .(cell, region, log2ratio)]
  rr = sort(unique(f$region))
  L = matrix(NA_real_, length(cells), length(rr), dimnames = list(cells, rr))
  L[cbind(match(f$cell, cells), match(f$region, rr))] = f$log2ratio
  rg = as.data.table(regions)[region %in% rr]
  if ("allelic_ok" %in% names(rg)) rg = rg[allelic_ok == TRUE]
  if (!is.null(allelic_chroms)) rg = rg[chrom %in% allelic_chroms]
  out = list(L = L, regions = rr, cells = cells)
  if (!is.null(votes) && nrow(rg)) {
    v = region_block_votes(votes, atoms, rg)
    v = v[cell <= length(cells)]
    nb = max(v$blk)
    out$MA = Matrix::sparseMatrix(i = v$cell, j = v$blk, x = as.numeric(v$nA), dims = c(length(cells), nb))
    out$MB = Matrix::sparseMatrix(i = v$cell, j = v$blk, x = as.numeric(v$nB), dims = c(length(cells), nb))
    bm = unique(v[, .(blk, region)])
    out$blk_region = match(bm$region[order(bm$blk)], rr)
    nn = v[, .(n = sum(nA + nB)), by = .(cell, region)]
    N = matrix(0, length(cells), length(rr)); N[cbind(nn$cell, match(nn$region, rr))] = nn$n
    out$N = N
    out$trip = v[, .(ci = cell, blk, nA = as.numeric(nA), nB = as.numeric(nB))]
  }
  class(out) = "contrast_data"
  out
}

# weighted bin counts: sum of w per integer key in 1..n (duplicates summed in C by the sparse constructor)
bin_sum = function(key, w, n) {
  if (!length(key)) return(numeric(n))
  as.vector(Matrix::sparseMatrix(i = as.integer(key), j = rep.int(1L, length(key)), x = as.numeric(w), dims = c(n, 1L)))
}

# Yuen's trimmed-mean test per column of two matrices (rows = cells, NA = missing): returns t^2 / 2 of the
# trimmed mean difference beyond `delta`, with winsorised variances; 0 when a side has < 3 values
yuen_cols = function(A, B, trim = 0.25, delta = 0) {
  a = trim_moments(A, trim); b = trim_moments(B, trim)
  out = pmax(abs(a$m - b$m) - delta, 0)^2 / (a$d + b$d) / 2
  out[a$n < 3 | b$n < 3 | !is.finite(out)] = 0
  out
}

# one-sample version: t^2 / 2 of the trimmed mean of each column vs 0
yuen_one = function(A, trim = 0.25) {
  a = trim_moments(A, trim)
  out = a$m^2 / a$d / 2
  out[a$n < 3 | !is.finite(out)] = 0
  out
}

# trimmed mean (m), its squared standard error from the winsorised variance (d) and n per column (NA = missing)
trim_moments = function(M, trim = 0.25) {
    n = colSums(!is.na(M))
    S = apply(M, 2, sort, na.last = TRUE)
    if (!is.matrix(S)) S = matrix(S, nrow = nrow(M))
    S[is.na(S)] = 0
    C1 = rbind(0, apply(S, 2, cumsum)); C2 = rbind(0, apply(S^2, 2, cumsum))
    if (!is.matrix(C1)) { C1 = matrix(C1, ncol = ncol(M)); C2 = matrix(C2, ncol = ncol(M)) }
    g = floor(trim * n); h = n - 2 * g
    j = seq_len(ncol(M))
    lo = cbind(pmax(g, 0) + 1, j); hi = cbind(pmax(n - g, 0) + 1, j)
    ts = C1[hi] - C1[lo]; ts2 = C2[hi] - C2[lo]
    m = ts / pmax(h, 1)
    ylo = S[cbind(pmin(pmax(g + 1, 1), nrow(S)), j)]; yhi = S[cbind(pmin(pmax(n - g, 1), nrow(S)), j)]
    sw = ts + g * ylo + g * yhi; sw2 = ts2 + g * ylo^2 + g * yhi^2
    vw = pmax((sw2 - sw^2 / pmax(n, 1)) / pmax(n - 1, 1), 1e-8)
    d = (n - 1) * vw / pmax(h * (h - 1), 1)
    list(n = n, m = m, d = d)
}

# sparse (cell x region) x block vote matrices of contrast data (rows: cell + (region - 1) * n_cells): the
# oriented counts are XA %*% [block oriented A] + XB %*% [block flipped], so a permutation costs a few sparse
# products over blocks instead of passes over all cell x region x block rows
contrast_agg = function(cd) {
  nc = length(cd$cells); nr = length(cd$regions); nb = ncol(cd$MA)
  key = as.integer(cd$trip$ci + (cd$blk_region[cd$trip$blk] - 1L) * nc)
  XA = Matrix::sparseMatrix(i = key, j = cd$trip$blk, x = cd$trip$nA, dims = c(nc * nr, nb))
  XB = Matrix::sparseMatrix(i = key, j = cd$trip$blk, x = cd$trip$nB, dims = c(nc * nr, nb))
  list(XA = XA, XB = XB, XN = XA + XB)
}

# Contrast between two cell sets
#
# Max over regions of max(depth term, allele term), both Yuen trimmed-mean tests (`trim` = 25 % per side)
# between the per-cell values of the two sets, reported as t^2 / 2. Depth: the per-cell log2 ratio, only
# the part of the trimmed-mean difference beyond `depth_delta` (log2) counts (cell types differ in
# chromatin depth waves by a few tenths of a log2 unit at arm / 20 Mb scale; a clonal one-copy change at
# contrast c >= 0.6 by >= 0.35). Alleles: per-cell oriented z `(nA - n / 2) / sqrt(n / 4 * (1 + (n - 1) *
# rho))` (NA below 5 molecules), with every other cell of one set orienting the phase blocks (its
# pseudobulk) and the remaining cells of that set compared with the other set under that orientation (no
# compared cell contributes to the orientation); the larger of the two directions is kept, so imbalance on
# opposite haplotypes and imbalanced vs balanced are both seen. Trimming makes the test about the bulk of each set: a minority of cells with
# private events (normal tissue) does not make two leaves different. A region's depth term counts only when its
# allele term exceeds `allele_support` (default 8, t ~ 4; regions without an allele channel, chrX / chrY, never separate clones): leaves of normal tissue that share
# a chromatin / replication state differ in depth by tenths of a log2 unit at arm scale, while
# every clonal copy-number change except a balanced one (2:2) also shifts the haplotype fraction.
contrast_stat = function(cd, a, b, rho, depth_delta = 0.2, trim = 0.25, min_n = 5, allele_support = 8,
                         depth_support = c("imbalance", "difference")) {
  depth_support = match.arg(depth_support)
  nr = length(cd$regions)
  dep = yuen_cols(cd$L[a, , drop = FALSE], cd$L[b, , drop = FALSE], trim, depth_delta)
  al = rep(0, nr); imb = rep(0, nr)
  if (!is.null(cd$trip)) {
    nc = length(cd$cells)
    # fixed aggregation (triplet -> cell x region, column-major); a permutation only changes the orientation
    agg = cd$agg %||% contrast_agg(cd)
    na_ = length(a); nb_ = length(b)
    for (h in 1:2) {
      # orient by every other cell of one set; compare the rest of that set with the other set, so no compared
      # cell contributes to the orientation (leave-one-out would make the orienting set's own cells
      # systematically different from the other set's)
      own = if (h == 1) a else b
      oth = if (h == 1) b else a
      ori = own[seq(1, length(own), by = 2)]; own_t = setdiff(own, ori)
      if (length(own_t) < 3) next
      ind = numeric(nc); ind[ori] = 1
      GA = as.vector(Matrix::crossprod(cd$MA, ind)); GB = as.vector(Matrix::crossprod(cd$MB, ind))
      pos_ = as.numeric(GA > GB); neg_ = as.numeric(GA < GB)
      OA = matrix(as.vector(agg$XA %*% pos_ + agg$XB %*% neg_), nc, nr); ON = matrix(as.vector(agg$XN %*% (pos_ + neg_)), nc, nr)
      rows = c(own_t, oth)
      OA = OA[rows, , drop = FALSE]; ON = ON[rows, , drop = FALSE]
      Z = (OA - ON / 2) / sqrt(pmax(ON, 1) / 4 * (1 + (ON - 1) * rho))
      Z[ON < min_n] = NA
      al = pmax(al, yuen_cols(Z[seq_along(own_t), , drop = FALSE], Z[length(own_t) + seq_along(oth), , drop = FALSE], trim, 0))
      # imbalance of the orienting set's own compared cells (one-sample, oriented by the other half of the set)
      imb = pmax(imb, yuen_one(Z[seq_along(own_t), , drop = FALSE], trim))
    }
  }
  dep[!is.finite(dep)] = 0; al[!is.finite(al)] = 0; imb[!is.finite(imb)] = 0
  # a depth difference counts only where the allele channel supports a copy-number change: the two sets differ
  # allelically there or (depth_support = "imbalance") either set is allelically imbalanced there (every clonal
  # copy-number change except a balanced one shifts the haplotype fraction; a loss and a gain / CN-LOH of the same
  # haplotype differ only in depth, e.g. a chr21 1:0 vs 2:0 subclone pair). Cells sharing a chromatin / replication state cluster
  # together and differ in depth by tenths of a log2 unit at balanced regions; regions without an allele channel
  # (chrX / chrY) never count
  has_al = if (is.null(cd$blk_region)) rep(FALSE, nr) else seq_len(nr) %in% cd$blk_region
  sup = if (depth_support == "imbalance") pmax(al, imb) else al
  dep_used = if (is.finite(allele_support)) ifelse(has_al & sup > allele_support, dep, 0) else dep
  list(stat = max(c(dep_used, al)), by_region = data.table(region = cd$regions, depth = dep, allelic = al, imbalance = imb,
                                                           depth_used = dep_used))
}

# contrast data restricted to a set of cells (rows renumbered 1..length(idx) in the order of idx)
subset_contrast_data = function(cd, idx) {
  out = cd
  out$cells = cd$cells[idx]
  out$L = cd$L[idx, , drop = FALSE]
  if (!is.null(cd$trip)) {
    pos = integer(length(cd$cells)); pos[idx] = seq_along(idx)
    tr = cd$trip[pos[ci] > 0]
    tr[, ci := pos[ci]]
    out$trip = tr
    out$MA = cd$MA[idx, , drop = FALSE]; out$MB = cd$MB[idx, , drop = FALSE]
    out$N = NULL
    out$agg = contrast_agg(out)
  }
  out
}

# permutation split test: is the contrast of (a, b) above the q-quantile of random splits of a + b at the
# same sizes? Sequential early stop: "same" as soon as enough null values reach the observed one that the
# q-quantile of all n_perm values cannot fall below it (exact), "different" after 20 / 50 / 100 permutations
# when the observed exceeds 3 x the null maximum (heuristic). With an early stop the reported thr is the
# quantile of the permutations run.
split_test = function(cd, a, b, rho, q = 0.999, n_perm = 200, max_cells = 250, seed = 1, early = TRUE, depth_delta = 0.2,
                      with_p = FALSE) {
  set.seed(seed)
  if (length(a) > max_cells) a = sample(a, max_cells)
  if (length(b) > max_cells) b = sample(b, max_cells)
  # restrict the contrast data to the tested cells once (permutations then touch only their votes)
  cd = subset_contrast_data(cd, c(a, b))
  a = seq_along(a); b = length(a) + seq_along(b)
  obs = contrast_stat(cd, a, b, rho, depth_delta)
  pool = c(a, b); na = length(a)
  null = numeric(0)
  # exact stop for "same": the type-7 q-quantile of n_perm values is at least their order statistic
  # floor(q (n_perm - 1)) + 1, so once n_perm - floor(q (n_perm - 1)) null values reach the observed one the
  # pair cannot differ (n_perm 200, q 0.999: two values)
  k_same = n_perm - floor(q * (n_perm - 1))
  for (p in seq_len(n_perm)) {
    pp = sample(pool)
    null[p] = contrast_stat(cd, pp[seq_len(na)], pp[-seq_len(na)], rho, depth_delta)$stat
    if (early && sum(null >= obs$stat) >= k_same) break
    if (early && p %in% c(20, 50, 100) && obs$stat > 3 * max(null)) break
  }
  thr = unname(stats::quantile(null, q, type = 7))
  top = obs$by_region[order(-pmax(depth_used, allelic))][1]
  out = data.table(stat = obs$stat, thr = thr, n_perm = length(null), differ = obs$stat > thr,
                   n_a = length(a), n_b = length(b), top_region = top$region, top_depth = top$depth, top_allelic = top$allelic)
  # permutation p (with an early stop: over the permutations run); only on request, so the merge tests' tables keep their columns
  if (with_p) out[, p := (1 + sum(null >= obs$stat)) / (1 + length(null))]
  out
}

#' Merge over-split clone leaves with a permutation-calibrated split test
#'
#' Bottom-up: each round, every clone's `n_nb` nearest clones (Ward distance between the clones' centroids in
#' the smoothed PC space) are candidate pairs; the pairs are tested in order of distance with the clone contrast (`split_test`)
#' against the `q`-quantile of the same statistic on random splits of the pooled pair at matching sizes
#' (`n_perm` permutations, sequential early stop), and non-different pairs are merged (each clone at most once
#' per round). Stops when every candidate pair of the remaining clones differs at the budget (`q = 0.999`, aimed at 1
#' false split per 1,000 tests; with `n_perm = 200` the type-7 0.999 quantile of the permutations sits near their
#' maximum, so the quantile rule alone splits a null pair at about 1 / 201 per test, and the early "different" stop
#' (observed statistic above 3 x the null maximum after 20 / 50 / 100 permutations) adds about 0.0006-0.001 per test
#' at the minimum leaf size, 20 against 20 cells; there is no correction for the number of tests). Each half is
#' subsampled to `max_cells` cells (observed and null alike).
#' Nearest neighbours rather than Ward siblings only: once two sibling clones differ, a sibling-only rule never
#' compares either of them with any other clone.
#' @param feat_long features long table
#' @param votes,atoms,regions votes, atoms, regions
#' @param leaves named leaf labels (names = cells), from [discover_clones()]
#' @param tree the [discover_clones()] result (its `smooth` coordinates are used) or a cells x dims matrix
#' @param q null quantile (budget)
#' @param n_perm permutations per test
#' @param seed RNG seed
#' @param rho beta-binomial overdispersion
#' @param classes region classes of the contrast
#' @param cells cell names in vote index order
#' @param max_cells subsample per half
#' @param cdata optional pre-computed [contrast_data()]
#' @param depth_delta minimal depth effect (log2) of the contrast
#' @param n_nb nearest clones tested per clone and round
#' @return list(clones = named integer labels, tests = data.table of all tests)
#' @export
merge_clones = function(feat_long, votes, atoms, regions, leaves, tree, q = 0.999, n_perm = 200, seed = 1,
                        rho = 0.01, classes = c("arm", "20Mb"), cells = NULL, max_cells = 250, cdata = NULL, depth_delta = 0.2,
                        n_nb = 3) {
  S = if (is.matrix(tree)) tree else tree$smooth
  S = S[names(leaves), , drop = FALSE]
  if (is.null(cells)) cells = names(leaves)
  if (is.null(cdata)) cdata = contrast_data(feat_long, votes, atoms, regions, cells, classes)
  idx = match(names(leaves), cdata$cells)
  cur = as.character(leaves)
  comp = split(as.character(leaves), cur); comp = lapply(comp, unique)
  cache = list(); tests = list()
  ckey = function(l) paste(sort(as.integer(comp[[l]])), collapse = ",")
  run_test = function(l1, l2) {
    k1 = ckey(l1); k2 = ckey(l2)
    key = if (k1 < k2) paste(k1, k2, sep = "|") else paste(k2, k1, sep = "|")
    if (!is.null(cache[[key]])) return(cache[[key]])
    tt = split_test(cdata, idx[cur == l1], idx[cur == l2], rho, q, n_perm, max_cells, seed = seed + length(cache), depth_delta = depth_delta)
    tt[, `:=`(key = key, l1 = l1, l2 = l2)]
    cache[[key]] <<- tt
    tests[[length(tests) + 1]] <<- tt
    tt
  }
  repeat {
    labs = sort(unique(cur))
    if (length(labs) < 2) break
    sz = as.vector(table(cur)[labs])
    cent = rowsum(S, cur)[labs, , drop = FALSE] / sz
    # Ward distance between clone centroids; candidates = each clone's n_nb nearest clones, closest first
    D = as.matrix(stats::dist(cent))^2 * outer(sz, sz, function(x, y) x * y / (x + y))
    diag(D) = Inf
    cand = unique(rbindlist(lapply(seq_along(labs), function(i) {
      j = order(D[i, ])[seq_len(min(n_nb, length(labs) - 1))]
      data.table(i = pmin(i, j), j = pmax(i, j))
    })))
    cand[, d := D[cbind(i, j)]]
    setorder(cand, d)
    merged = FALSE
    used = character()
    for (r in seq_len(nrow(cand))) {
      l1 = labs[cand$i[r]]; l2 = labs[cand$j[r]]
      if (l1 %in% used || l2 %in% used) next
      tt = run_test(l1, l2)
      if (!tt$differ) {
        cur[cur == l2] = l1
        comp[[l1]] = c(comp[[l1]], comp[[l2]]); comp[[l2]] = NULL
        used = c(used, l1, l2)
        merged = TRUE
      }
    }
    if (!merged) break
  }
  lab = as.integer(factor(cur, levels = unique(cur)))
  names(lab) = names(leaves)
  list(clones = lab, tests = rbindlist(tests, fill = TRUE))
}
utils::globalVariables(c("differ", "l1", "l2", "depth_used"))

#' Clone discovery within a trunk
#'
#' Second level of the hierarchical discovery of [caller_tumour()]. The phase blocks are oriented by the pseudobulk
#' of the trunk's cells only (leave-one-out, discovery half of the votes): a subclonal event carried by 5-25 % of a
#' trunk is oriented correctly in most blocks, while in the all-cell pseudobulk it is outweighed by balanced cells.
#' Features: per region (classes of `feat_disc`) the signed allelic z and the depth z, kept only for regions whose
#' allelic z varies across the trunk's cells beyond the null (variance > `disp_min`; unit variance expected) or is
#' imbalanced on average (|mean| > `imb_min`): allelic features select the regions, depth is used only where the
#' allele channel shows a clone-level event (a loss vs a gain of the same haplotype differ only in depth there).
#' Discovery ([discover_clones()]) runs on the trunk's deeper cells (`n_frag` >= the `deep_q` quantile, at least
#' `4 * min_size` cells), whose per-cell profiles separate subclones; [merge_clones()] then tests every pair of leaves
#' (`n_nb = Inf`) on the test halves (`cdata`). Shallow cells are left to the likelihood assignment.
#' @param inp a `clone_input`
#' @param tcells trunk cells
#' @param votes_disc discovery half of the votes ([split_votes()])
#' @param feat_disc per-cell depth features of the discovery half (cell, region, z)
#' @param cdata [contrast_data()] of the test halves (all cells)
#' @param n_frag named clean fragments per cell
#' @param opts [clone_defaults()]
#' @param seed RNG seed
#' @return list(clones = named integer labels of the discovery cells, tests, regions = selected regions, n_leaves)
#' @export
discover_within = function(inp, tcells, votes_disc, feat_disc, cdata, n_frag, opts = clone_defaults(), seed = 1) {
  cells = inp$cells$cell
  ti = match(tcells, cells)
  rg = inp$regions[class %in% opts$classes & eligible == TRUE & allelic_ok == TRUE & chrom %in% inp$allelic_chroms]
  ob = cl_orient(votes_disc[cell %in% ti], rg, inp$atoms, NULL, loo = TRUE)
  ob[, cell := cells[cell]]
  f = as.data.table(feat_disc)[cell %in% tcells & class %in% opts$classes, .(cell, region, z)]
  f = merge(f, ob[, .(cell, region, nA, nB)], by = c("cell", "region"), all.x = TRUE)
  f[is.na(nA), `:=`(nA = 0, nB = 0)]
  f[, n := nA + nB]
  f[, az := fifelse(n >= opts$min_n, (nA - n / 2) / sqrt(n / 4 * (1 + (n - 1) * inp$rho)), NA_real_)]
  f[, az := pmin(pmax(az, -opts$cap), opts$cap)]
  f[, dz := pmin(pmax(fifelse(is.finite(z), z, 0), -opts$cap), opts$cap)]
  sel = f[, .(v = stats::var(az, na.rm = TRUE), m = mean(az, na.rm = TRUE), n_ok = sum(!is.na(az))), by = region]
  keep = sel[n_ok >= 10 & (v > opts$disp_min | abs(m) > opts$imb_min), region]
  deep = tcells[n_frag[tcells] >= stats::quantile(n_frag[tcells], opts$deep_q, na.rm = TRUE)]
  if (length(deep) < 4 * opts$min_size) deep = tcells[order(-n_frag[tcells])][seq_len(min(length(tcells), 4 * opts$min_size))]
  one = function(why) list(clones = stats::setNames(rep(1L, length(deep)), deep), tests = data.table(), regions = keep,
                           n_leaves = 1L, note = why)
  if (!length(keep) || length(deep) < 2 * opts$min_size) return(one("no selected region or too few cells"))
  f = f[cell %in% deep & region %in% keep]
  f[is.na(az), az := 0]
  rr = sort(keep); cc = deep
  ii = cbind(match(f$cell, cc), match(f$region, rr))
  A = matrix(0, length(cc), length(rr)); Dm = A
  A[ii] = f$az; Dm[ii] = f$dz
  X = cbind(A, opts$depth_w * Dm)
  colnames(X) = c(paste0("a|", rr), paste0("d|", rr)); rownames(X) = cc
  X = X[, colSums(X != 0) > 0, drop = FALSE]
  if (ncol(X) < 2) return(one("fewer than 2 informative features"))
  disc = discover_clones(X, opts$n_pc, opts$k, opts$max_leaves_trunk, opts$min_size, seed)
  if (length(unique(disc$leaves)) < 2) return(one("one leaf"))
  mc = merge_clones(NULL, NULL, inp$atoms, inp$regions, disc$leaves, disc, q = opts$q, n_perm = opts$n_perm, seed = seed,
                    rho = inp$rho, classes = opts$classes, cells = cells, max_cells = opts$max_cells, cdata = cdata,
                    depth_delta = opts$depth_delta, n_nb = opts$n_nb_trunk)
  list(clones = mc$clones, tests = mc$tests, regions = keep, n_leaves = length(unique(disc$leaves)), note = "")
}
utils::globalVariables(c("n_ok", "v", "m"))

# two-means on a numeric vector (1-D Lloyd from the quartiles; split_disagree): list(g = 1 / 2 per value, m = means,
# sep = mean difference over the pooled within-group SD)
two_means = function(x, n_iter = 25) {
  ctr = stats::quantile(x, c(0.25, 0.75), names = FALSE)
  g = rep(1L, length(x))
  for (i in seq_len(n_iter)) {
    g_new = ifelse(abs(x - ctr[1]) <= abs(x - ctr[2]), 1L, 2L)
    if (length(unique(g_new)) < 2) break
    ctr = c(mean(x[g_new == 1]), mean(x[g_new == 2]))
    if (all(g_new == g)) break
    g = g_new
  }
  if (length(unique(g)) < 2) return(list(g = g, m = c(mean(x), mean(x)), sep = 0))
  s = sqrt((stats::var(x[g == 1]) + stats::var(x[g == 2])) / 2)
  list(g = g, m = ctr, sep = abs(diff(ctr)) / max(s, 1e-9))
}

#' Alternative states of a clone segment for the disagreement split
#'
#' The states a sub-population of a clone could carry on a segment in clone state `state`: states of the grid (major
#' haplotype first, at least one copy) within one copy of its total copy number. An alternative whose major fraction
#' differs by less than `min_df` (a depth-only difference) is kept only when one of the two states has <= 2 copies (2:0 vs
#' 1:0 vs 3:0): between amplified states one copy is a small log2 step, and amplicon copy number varies continuously
#' between cells (e.g. an 8q amplicon at cn 6-16), as [clone_state_distance()] also ignores it.
#' @param state clone state (`"nA:nB"`)
#' @param grid state grid ([cn_grid()])
#' @param min_df minimal major-fraction difference of an alternative between amplified states
#' @return grid rows of the alternatives (state, nA, nB, cn, f, label)
#' @export
split_alt_states = function(state, grid = cn_grid(), min_df = 0.15) {
  st0 = as.character(state)[1]
  g = as.data.table(grid)[cn >= 1 & nA >= nB]
  s0 = g[g$state %in% st0]
  if (nrow(s0) != 1) return(g[0])
  g[g$state != s0$state & abs(cn - s0$cn) <= 1 & (abs(f - s0$f) >= min_df - 1e-9 | pmin(cn, s0$cn) <= 2)]
}

# contrast data restricted to a set of regions (columns of L, the phase blocks of those regions)
subset_contrast_regions = function(cd, regs) {
  R = match(regs, cd$regions); R = R[!is.na(R)]
  out = cd
  out$regions = cd$regions[R]
  out$L = cd$L[, R, drop = FALSE]
  out$agg = NULL
  if (!is.null(cd$trip)) {
    B = which(cd$blk_region %in% R)
    if (!length(B)) {
      out$trip = NULL; out$MA = NULL; out$MB = NULL; out$blk_region = NULL; out$N = NULL
      return(out)
    }
    tr = cd$trip[blk %in% B]
    tr[, blk := match(blk, B)]
    out$trip = tr
    out$MA = cd$MA[, B, drop = FALSE]; out$MB = cd$MB[, B, drop = FALSE]
    out$blk_region = match(cd$blk_region[B], R)
    if (!is.null(cd$N)) out$N = cd$N[, R, drop = FALSE]
  }
  out
}

#' Split clones whose cells disagree with the clone state on a segment
#'
#' Targeted split pass of [caller_tumour()] (`split_pass = "disagree"`). Allele-driven discovery can put cells that share
#' a haplotype imbalance but differ by one copy in one clone (e.g. a chr21 loss subclone inside a chr21 CN-LOH clone;
#' a chr21 3:0 / 4:0 LOH-gain subclone inside a 2:1 clone, fitted 3:1). For each clone with >= 2 `min_size`
#' members and each of its non-neutral segments (ploidy 2, >= `split_min_len`; the copy-neutral 2:0 included), every
#' member is scored on the discovery halves (`win_disc`, `votes_disc`) under the clone state and under each alternative of
#' [split_alt_states()]. Depth: the cell's segment log2 ratio, centred on the members' median (the clone state; clone
#' profiles are relative to the clone baseline), in units of its `sd_null` scaled to the half depth and inflated by the
#' members' own robust spread, t density at `c * log2(cn / cn_state)` (clone c from `par`). Alleles: the cell's votes
#' oriented by the clone's leave-one-out pseudobulk, beta-binomial at `0.5 + a_eff (f - 0.5)`, a_eff = the attenuation of
#' the members' pooled oriented fraction. The contradiction of a segment is the fraction of members whose best state is
#' not the clone state; its noise floor is the same fraction on data simulated under the clone state at each cell's SD
#' and vote count (5 replicates). The alternative is the state with the largest excess preference over its noise floor.
#' A segment is a candidate when the contradiction is >= 10 % and exceeds its noise floor by >= `split_min_excess`, and a
#' two-means split of the per-cell log-likelihood difference (alternative - clone state) has both parts >=
#' max(`min_size`, `small_clone`, 10 %) and a separation >= `split_sep`. Candidates are ranked by the excess
#' contradiction (not by separation: the 8q amplicon's depth spread is unimodal and wide) over all clones; up to
#' `split_cand` per clone are tested, at most one split per clone and `split_max` splits in all. A candidate is applied
#' when (1) the permutation split test (`split_test()`, budget `q`) on the test halves, restricted to the `cdata` regions
#' the segment overlaps (>= 80 %, else >= 50 % of the shorter of the two; depth + allele), says different, and (2) the two parts re-fitted on the test
#' halves ([clone_profiles()], `fit_states`, [segment_clone()]) differ on the segment (copy number or major fraction,
#' [clone_state_distance()] with `min_df = 0.01, min_dcn = 1` as the caller's same-state merge) over >= max(`same_dist`,
#' half the segment), with the re-fitted states at >= 1 copy, both within 2 and one within 1 copies (total) of the clone state (the clone
#' state of a mixture is an average: e.g. a fitted 4:1 / 3:1 clone that is a mix of 2:1 and 3:0 / 4:0 cells).
#' @param members data.table cell, clone (confident members)
#' @param segments segments table of the members' clones
#' @param par data.table clone, c, a ([fit_clone_states()] `par`)
#' @param inp scn_features list or `clone_input` (regions, atoms, cells, features for `sd_null`, tpar, rho)
#' @param win_disc,votes_disc discovery halves of the window depth and the votes ([split_win()], [split_votes()])
#' @param cdata [contrast_data()] of the test halves
#' @param win_test,votes_test test halves for the re-fit (default: the full `inp$win` / `inp$votes`)
#' @param ref_cells depth reference cells of the re-fit (used with >= `min_ref` cells)
#' @param fit_states optional replacement of [fit_clone_states()] (`function(profile, grid)`)
#' @param opts [clone_defaults()] (`min_size`, `split_sep`, `split_max`, `split_cand`, `split_min_len`, `split_min_df`,
#'   `split_min_excess`, `q`, `n_perm`, `max_cells`, `depth_delta`, `same_dist`, `small_clone`, `size_small`, `min_ref`, `p_switch`,
#'   `min_win`, `min_llr`)
#' @param params [caller_defaults()] (`sd_floor`, `eps`)
#' @param seed RNG seed
#' @param frac fraction of the counts in the discovery half
#' @return list(members = relabelled members (split clones get suffixes `a` / `b`, `b` = the part that prefers the
#'   alternative), tests = one row per candidate: clone, region, chrom, start, end, state, alt, n, contra (fraction of
#'   members whose best state differs from the clone's), contra_null (its noise floor), excess, contra_alt /
#'   contra_alt_null (the same for the alternative), a_eff, infl (depth SD inflation), sep, n_a, n_b, d_a, d_b (mean
#'   log-likelihood difference per part), outcome, rank, stat, thr, n_perm, p (permutation p of the test-half split
#'   test), differ, state_a, state_b (re-fitted states at the segment midpoint), dist_bp, split), cells_b = named list of the moved cells per applied split (names = the `b` clones),
#'   sets = per candidate (rows of `tests`) the proposed parts list(a, b)
#' @export
split_disagree = function(members, segments, par, inp, win_disc, votes_disc, cdata, win_test = NULL, votes_test = NULL,
                          ref_cells = character(), fit_states = NULL, opts = clone_defaults(), params = caller_defaults(),
                          seed = 1, frac = 0.5) {
  inp = clone_input(inp)
  opts = utils::modifyList(clone_defaults(), as.list(opts))
  cells = inp$cells$cell
  members = copy(as.data.table(members))[, .(cell, clone = as.character(clone))]
  sg = as.data.table(segments)
  par = as.data.table(par)
  grid = cn_grid()
  fit_states = fit_states %||% function(profile, grid) fit_clone_states(profile, grid)
  inp_d = inp; inp_d$win = win_disc %||% inp$win
  vd = votes_disc %||% inp$votes
  win_t = win_test %||% inp$win; votes_t = votes_test %||% inp$votes
  ms = max(opts$min_size, opts$small_clone, 1)
  cands = list(); cell_sets = list()
  for (k in sort(unique(members$clone))) {
    mk = members[clone == k, cell]
    pk = par[clone == k]
    if (length(mk) < 2 * opts$min_size || !nrow(pk)) next
    ev_all = sg[clone == k & ploidy == 2 & (end - start) >= opts$split_min_len]
    ev = ev_all[label != "neutral"]
    if (!nrow(ev)) next
    regs = seg_regions(ev_all[, .(chrom, start, end)], inp$regions, tag = k)
    d = region_features(inp_d, regs, cells = mk, sd_floor = params$sd_floor)
    d[, sdh := pmax(sd_null, params$sd_floor) / sqrt(frac)]
    d[, `:=`(oA = 0, oN = 0)]
    tv = trip_votes(vd, inp$atoms, regs, length(cells), inp$allelic_chroms)
    if (nrow(tv$trip)) {
      oc = oriented_counts(tv, match(mk, cells))
      oc[, `:=`(cell = cells[ci], region = regs$region[ri])]
      d[oc, `:=`(oA = i.oA, oN = i.oN), on = .(cell, region)]
    }
    # the clone's own noise, as medians over its segments (one mixed segment does not move them): depth = robust SD of
    # the members' standardised segment depth (all segments >= split_min_len, neutral ones included; >= 0.5: sd_null
    # scaled to the half depth can over- or underestimate the spread), alleles =
    # method-of-moments beta-binomial rho of the members' oriented votes on the imbalanced segments
    nz = d[is.finite(log2ratio), .(infl = stats::mad((log2ratio - stats::median(log2ratio)) / sdh)), by = region]
    infl_k = max(0.5, stats::median(nz$infl, na.rm = TRUE), na.rm = TRUE)
    ev_reg = seg_regions(ev[, .(chrom, start, end)], inp$regions, tag = k)$region
    imb_reg = ev_reg[nA_of(ev$state) != nB_of(ev$state)]
    rr = d[region %in% imb_reg & oN >= 2, {
      f_ = min(max(sum(oA) / sum(oN), 0.02), 0.98)
      .(num = sum((oA - oN * f_)^2 / (oN * f_ * (1 - f_)) - 1), den = sum(oN - 1))
    }, by = region][den > 0]
    rho_k = if (nrow(rr)) max(stats::median(pmax(rr$num / rr$den, 0)), 1e-4) else inp$rho
    for (i in seq_len(nrow(ev))) {
      s0 = ev$state[i]
      alts = split_alt_states(s0, grid, opts$split_min_df)
      x = d[region == ev_reg[i] & is.finite(log2ratio)]
      if (!nrow(alts) || nrow(x) < 2 * opts$min_size) next
      st = rbind(grid[state == s0], alts)
      tp = cl_match_tpar(x, inp$tpar)
      set.seed(seed + length(cands))
      ds = disagree_ll(x$log2ratio, x$sdh, x$oA, x$oN, st, pk$c[1], pk$a[1], rho_k, params$eps, if (is.null(tp)) 5 else tp$t_df,
                       infl = infl_k)
      ja = ds$alt
      dd = ds$L[, ja] - ds$L[, 1]
      tm = if (length(unique(dd)) > 1) two_means(dd) else list(g = rep(1L, length(dd)), m = c(0, 0), sep = 0)
      hi = if (tm$m[2] >= tm$m[1]) 2L else 1L
      cb = x$cell[tm$g == hi]; ca_ = setdiff(mk, cb)
      n_b = length(cb); n_a = sum(tm$g != hi)
      ok = min(n_a, n_b) >= max(ms, 0.1 * nrow(x))
      id = length(cands) + 1L
      cands[[id]] = data.table(clone = k, region = ev_reg[i], chrom = ev$chrom[i], start = ev$start[i], end = ev$end[i],
                               state = s0, alt = st$state[ja], n = nrow(x), contra = ds$contra, contra_null = ds$contra_null,
                               excess = ds$contra - ds$contra_null, contra_alt = ds$by_state[ja], contra_alt_null = ds$by_state_null[ja],
                               a_eff = ds$a_eff, infl = ds$infl, rho = rho_k, sep = tm$sep, n_a = n_a, n_b = n_b, d_a = mean(dd[tm$g != hi]),
                               d_b = mean(dd[tm$g == hi]),
                               outcome = if (ds$contra < 0.1 || ds$contra - ds$contra_null < opts$split_min_excess) "contradiction within the noise floor" else
                                 if (!ok) "parts below min_size / 10 %" else if (tm$sep < opts$split_sep) "separation below split_sep" else "candidate")
      cell_sets[[id]] = list(a = ca_, b = cb)
    }
  }
  if (!length(cands)) return(list(members = members, tests = NULL, cells_b = list(), sets = list()))
  tests = rbindlist(cands)
  tests[, `:=`(id = .I, rank = NA_integer_, stat = NA_real_, thr = NA_real_, n_perm = NA_integer_, p = NA_real_, differ = NA,
               state_a = NA_character_, state_b = NA_character_, dist_bp = NA_real_, split = FALSE)]
  ord = tests[outcome == "candidate"][order(-excess, -sep), id]
  tests[ord, rank := seq_along(ord)]
  crg = as.data.table(inp$regions)[region %in% cdata$regions, .(region, chrom, start, end)]
  refit = function(m2) {
    sz = m2[, .N, by = clone]
    size = stats::setNames(ifelse(sz$N < opts$small_clone, opts$size_small, 1e6), sz$clone)
    prof = clone_profiles(win_t, votes_t, inp$atoms, m2, ref_cells = if (length(ref_cells) >= opts$min_ref) ref_cells else NULL,
                          size = size, seed = seed, cells = cells, rho = inp$rho, allelic_chroms = inp$allelic_chroms,
                          haploid = inp$haploid, min_ref = opts$min_ref)
    segment_clone(fit_states(prof, grid), opts$p_switch, opts$min_win, grid, opts$min_llr, opts$min_excess, opts$bal_frac,
                  opts$min_bal_len, opts$max_contra)
  }
  done = character(); n_tested = list(); moved = list()
  for (ii in ord) {
    tk = tests[ii]
    k = tk$clone
    if (k %in% done) { tests[ii, outcome := "clone already split"]; next }
    if (sum(tests$split) >= opts$split_max) { tests[ii, outcome := "split_max reached"]; next }
    if ((n_tested[[k]] %||% 0L) >= opts$split_cand) { tests[ii, outcome := "not tested (split_cand)"]; next }
    n_tested[[k]] = (n_tested[[k]] %||% 0L) + 1L
    ov = crg[chrom == tk$chrom & end > tk$start & start < tk$end]
    ov[, fr := (pmin(end, tk$end) - pmax(start, tk$start)) / pmin(end - start, tk$end - tk$start)]
    treg = if (any(ov$fr >= 0.8)) ov[fr >= 0.8, region] else ov[fr >= 0.5, region]
    if (!length(treg)) { tests[ii, outcome := "no test region"]; next }
    cs = cell_sets[[ii]]
    ia = match(cs$a, cdata$cells); ib = match(cs$b, cdata$cells)
    tt = split_test(subset_contrast_regions(cdata, treg), ia[!is.na(ia)], ib[!is.na(ib)], inp$rho, opts$q, opts$n_perm,
                    opts$max_cells, seed = seed + ii, depth_delta = opts$depth_delta, with_p = TRUE)
    tests[ii, `:=`(stat = tt$stat, thr = tt$thr, n_perm = tt$n_perm, p = tt$p, differ = tt$differ)]
    if (!tt$differ) { tests[ii, outcome := "not different on the test halves"]; next }
    ka = paste0(k, "a"); kb = paste0(k, "b")
    sgr = refit(data.table(cell = c(cs$a, cs$b), clone = c(rep(ka, length(cs$a)), rep(kb, length(cs$b)))))
    mid = (tk$start + tk$end) / 2
    at = function(kk) { y = sgr[clone == kk & chrom == tk$chrom & start <= mid & end > mid]; if (nrow(y)) y$state[1] else NA_character_ }
    clip = sgr[chrom == tk$chrom & end > tk$start & start < tk$end]
    clip[, `:=`(start = pmax(start, tk$start), end = pmin(end, tk$end))]
    dm = if (uniqueN(clip$clone) == 2) clone_state_distance(clip, min_df = 0.01, min_dcn = 1) else matrix(0, 2, 2)
    dist = dm[1, 2]
    tests[ii, `:=`(state_a = at(ka), state_b = at(kb), dist_bp = dist)]
    if (dist < max(opts$same_dist, 0.5 * (tk$end - tk$start))) { tests[ii, outcome := "same states after re-fit"]; next }
    cn_ = function(s_) if (is.na(s_)) NA_integer_ else nA_of(s_) + nB_of(s_)
    dcn = abs(c(cn_(at(ka)), cn_(at(kb))) - cn_(tk$state))
    if (!isTRUE(all(dcn <= 2) && min(dcn) <= 1 && all(c(cn_(at(ka)), cn_(at(kb))) >= 1))) {
      tests[ii, outcome := "re-fitted states too far from the clone state"]; next
    }
    members[clone == k & cell %in% cs$b, clone := kb]
    members[clone == k, clone := ka]
    tests[ii, `:=`(split = TRUE, outcome = "split")]
    done = c(done, k)
    moved[[kb]] = cs$b
  }
  tests[, id := NULL]
  list(members = members, tests = tests, cells_b = moved, sets = cell_sets)
}
utils::globalVariables(c("sdh", "zh", "oA", "outcome", "contra", "contra_alt", "fr", "dist_bp", "state_a", "state_b", "excess"))

# per-cell log-likelihoods of one clone segment under the clone state (row 1 of `st`) and its alternatives, with a
# parametric noise floor. Depth: the cells' log2 ratio `l` centred on its median (= the clone state; the clone profiles
# are centred on the clone baseline, not on a reference) with the per-cell SD `sdh` inflated by `infl` (default the
# robust spread of these values, >= 1), a t(df) density at the state shifts c * log2(cn / cn_state). Alleles:
# beta-binomial (rho) at 0.5 + a_eff (f - 0.5), with a_eff from the clone's pooled oriented fraction (the attenuation the
# per-cell orientation actually has; `a_` when the clone state is balanced or there are no votes). The contradiction
# (fraction of cells whose best state is not the clone state) is compared with its value on `n_null` replicates of data
# simulated under the clone state at each cell's own SD and vote count (t noise, beta-binomial votes).
disagree_ll = function(l, sdh, oA, oN, st, c_, a_, rho, eps = 0.02, df = 5, n_null = 5, infl = NULL) {
  n = length(l)
  m = stats::median(l)
  z0 = (l - m) / sdh
  infl = infl %||% max(1, stats::mad(z0))
  s = sdh * infl
  nh = sum(oN)
  fo = if (nh > 0) sum(oA) / nh else NA_real_
  a_eff = if (abs(st$f[1] - 0.5) > 1e-9 && is.finite(fo)) min(max((fo - 0.5) / (st$f[1] - 0.5), 0.2), 1) else a_
  mu = c_ * (log2(pmax(st$cn, 0.2)) - log2(pmax(st$cn[1], 0.2)))
  fe = 0.5 + a_eff * (st$f - 0.5)
  ll = function(zz, xa) matrix(vapply(seq_len(nrow(st)), function(j)
    stats::dt(zz - mu[j] / s, df, log = TRUE) + ll_allelic(xa, oN, fe[j], rho, eps, two_sided = FALSE, per_row = TRUE), numeric(n)), n)
  L = ll(z0 / infl, oA)
  best = max.col(L, ties.method = "first")
  r_ = max(rho, 1e-6)
  bn = unlist(lapply(seq_len(n_null), function(r) {
    pb = stats::rbeta(n, fe[1] * (1 / r_ - 1), (1 - fe[1]) * (1 / r_ - 1))
    max.col(ll(stats::rt(n, df), stats::rbinom(n, oN, pb)), ties.method = "first")
  }))
  by_s = tabulate(best, nrow(st)) / n
  by_n = tabulate(bn, nrow(st)) / length(bn)
  ex = by_s - by_n
  ex[1] = -Inf
  list(L = L, alt = which.max(ex), contra = mean(best != 1L), contra_null = mean(bn != 1L), by_state = by_s, by_state_null = by_n,
       a_eff = a_eff, infl = infl)
}
