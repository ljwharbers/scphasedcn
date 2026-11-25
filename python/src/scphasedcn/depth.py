"""Per-cell x bin depth layers from fragment arrays (vectorised numpy accumulators)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .bam import ChromFragments, F_BLACKLIST
from .bins import Tiling

ASSIGN_MODES = ("start", "midpoint", "end5", "overlap")
FRAGLEN_EDGES = np.array([0, 200, 300, 500, 750, 1000, 2000, 5000, 10000, np.inf])
FRAGLEN_LABELS = ["lt200", "200_300", "300_500", "500_750", "750_1000", "1k_2k", "2k_5k", "5k_10k", "ge10k"]


def assigned_position(fr: ChromFragments, assign: str) -> np.ndarray:
    if assign == "start":
        return fr.start
    if assign == "midpoint":
        return (fr.start + fr.end) // 2
    if assign == "end5":
        return np.where(fr.strand == 0, fr.start, fr.end - 1)
    raise ValueError(assign)


DEPTH_LAYERS = ["counts", "span", "frag_short", "frag_long"]
HAP_LAYERS = ["hapA", "hapB", "hap_ambig"]


def accumulate(cell: np.ndarray, bin_idx: np.ndarray, n_cells: int, n_bins: int, weights=None) -> sp.csr_matrix:
    """cells x bins CSR of counts (int64) or summed weights (float64) for entries with bin_idx >= 0, built sparse
    (no dense cells x bins array): canonical (duplicates summed, sorted indices), entries summing to 0 not stored."""
    ok = bin_idx >= 0
    row = cell[ok].astype(np.int64); col = np.asarray(bin_idx[ok], dtype=np.int64)
    data = np.ones(len(row), dtype=np.int64) if weights is None else np.asarray(weights[ok], dtype=np.float64)
    m = sp.coo_matrix((data, (row, col)), shape=(n_cells, n_bins)).tocsr()
    m.sum_duplicates()
    m.eliminate_zeros()
    return m


def depth_layers(fr: ChromFragments, tiling: Tiling, n_cells: int, assign: str = "start",
                 use: np.ndarray | None = None, layers: list[str] | None = None) -> dict:
    """counts / span / frag_short / frag_long (DEPTH_LAYERS) and hapA / hapB / hap_ambig (HAP_LAYERS) layers, or only
    those named in `layers`; `use` selects fragments (default: not blacklisted)."""
    if use is None:
        use = (fr.flags & F_BLACKLIST) == 0
    want = set(DEPTH_LAYERS + HAP_LAYERS if layers is None else layers)
    nb = len(tiling)
    cell = fr.cell[use]; start = fr.start[use]; end = fr.end[use]; L = end - start
    out = {}
    if assign == "overlap":
        fi, bi = tiling.bins_overlapping(start, end)
        c = cell[fi]; Lc = L[fi]
        span_w = lambda: (np.minimum(end[fi], tiling.ends[bi])                      # bases in the bin
                          - np.maximum(start[fi], tiling.starts[bi])).astype(np.float64)
        hap = fr.hap[use][fi]
    else:
        bi = tiling.bin_of(assigned_position(fr, assign)[use])
        c = cell; Lc = L
        span_w = lambda: L.astype(np.float64)
        hap = fr.hap[use]
    if "counts" in want:
        out["counts"] = accumulate(c, bi, n_cells, nb)
    if "span" in want:
        out["span"] = accumulate(c, bi, n_cells, nb, weights=span_w())
    if "frag_short" in want:
        out["frag_short"] = accumulate(c[Lc < 300], bi[Lc < 300], n_cells, nb)
    if "frag_long" in want:
        out["frag_long"] = accumulate(c[Lc >= 1000], bi[Lc >= 1000], n_cells, nb)
    for code, name in ((1, "hapA"), (2, "hapB"), (3, "hap_ambig")):
        if name in want:
            sel = hap == code
            out[name] = accumulate(c[sel], bi[sel], n_cells, nb)
    for k in out:
        out[k] = out[k].astype(np.float64 if k == "span" else np.int32)
    return out


def fraglen_hist(fr: ChromFragments, n_cells: int, use: np.ndarray | None = None) -> np.ndarray:
    if use is None:
        use = (fr.flags & F_BLACKLIST) == 0
    L = (fr.end - fr.start)[use]
    b = np.searchsorted(FRAGLEN_EDGES, L, side="right") - 1
    b = np.clip(b, 0, len(FRAGLEN_LABELS) - 1)
    return np.bincount(fr.cell[use].astype(np.int64) * len(FRAGLEN_LABELS) + b,
                       minlength=n_cells * len(FRAGLEN_LABELS)).reshape(n_cells, len(FRAGLEN_LABELS))
