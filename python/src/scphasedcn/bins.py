"""Genomic tilings (regular bins), masks and per-bin annotation."""
from __future__ import annotations

import gzip
from dataclasses import dataclass

import numpy as np


@dataclass
class Tiling:
    """Non-overlapping, sorted tiles on one chromosome (0-based half-open)."""

    chrom: str
    starts: np.ndarray
    ends: np.ndarray
    regular: int = 0  # bin size when the tiling is a regular grid from 0

    @classmethod
    def regular_grid(cls, chrom: str, length: int, binsize: int) -> "Tiling":
        n = int(np.ceil(length / binsize))
        starts = np.arange(n, dtype=np.int64) * binsize
        ends = np.minimum(starts + binsize, length)
        return cls(chrom, starts, ends, regular=binsize)

    def __len__(self):
        return len(self.starts)

    def bin_of(self, pos: np.ndarray) -> np.ndarray:
        """Tile index of each position, -1 outside all tiles."""
        pos = np.asarray(pos, dtype=np.int64)
        if self.regular:
            idx = pos // self.regular
            bad = (pos < 0) | (idx >= len(self))
            idx[bad] = -1
            ok = ~bad
            idx[ok] = np.where(pos[ok] < self.ends[idx[ok]], idx[ok], -1)
            return idx
        idx = np.searchsorted(self.starts, pos, side="right") - 1
        ok = (idx >= 0) & (pos < self.ends[np.clip(idx, 0, max(len(self) - 1, 0))])
        idx[~ok] = -1
        return idx

    def bins_overlapping(self, start: np.ndarray, end: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(fragment index, tile index) for every tile a fragment [start, end) overlaps (>= 1 bp)."""
        start = np.asarray(start, dtype=np.int64); end = np.asarray(end, dtype=np.int64)
        lo = np.searchsorted(self.ends, start, side="right")     # first tile with end > start
        hi = np.searchsorted(self.starts, end, side="left")      # tiles with start < end are [0, hi)
        n = np.maximum(hi - lo, 0)
        fi = np.repeat(np.arange(len(start)), n)
        off = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
        return fi, np.repeat(lo, n) + off

    def frame(self):
        import pandas as pd
        return pd.DataFrame({"chrom": self.chrom, "start": self.starts, "end": self.ends})


def read_bed_intervals(path: str | None, chrom: str) -> tuple[np.ndarray, np.ndarray]:
    """Sorted, merged intervals of one chromosome from a BED (3+ columns); empty arrays when path is None."""
    if path is None:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    s, e = [], []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            s.append(int(f[1])); e.append(int(f[2]))
    return merge_intervals(np.asarray(s, dtype=np.int64), np.asarray(e, dtype=np.int64))


def merge_intervals(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(starts) == 0:
        return starts.astype(np.int64), ends.astype(np.int64)
    o = np.argsort(starts, kind="stable"); starts, ends = starts[o], ends[o]
    ms, me = [starts[0]], [ends[0]]
    for s, e in zip(starts[1:], ends[1:]):
        if s <= me[-1]:
            me[-1] = max(me[-1], e)
        else:
            ms.append(s); me.append(e)
    return np.asarray(ms, dtype=np.int64), np.asarray(me, dtype=np.int64)


def in_intervals(pos: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Boolean: position falls in a (merged, sorted) interval set."""
    pos = np.asarray(pos, dtype=np.int64)
    if len(starts) == 0:
        return np.zeros(len(pos), dtype=bool)
    i = np.searchsorted(starts, pos, side="right") - 1
    return (i >= 0) & (pos < ends[np.clip(i, 0, len(starts) - 1)])


def overlap_bp(tiling: Tiling, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Bases of each tile covered by a (merged) interval set."""
    out = np.zeros(len(tiling), dtype=np.int64)
    if len(starts) == 0:
        return out
    fi, bi = tiling.bins_overlapping(starts, ends)
    ov = np.minimum(ends[fi], tiling.ends[bi]) - np.maximum(starts[fi], tiling.starts[bi])
    np.add.at(out, bi, ov)
    return out


def read_centromere(path: str, chrom: str) -> tuple[int, int] | None:
    """Battenberg chromosome_coordinates file: chr start cen.left.base cen.right.base end (1-based)."""
    names = {chrom, {"chrX": "chr23", "chrY": "chr24"}.get(chrom, chrom)}   # Battenberg names chrX / chrY chr23 / chr24
    with open(path) as fh:
        header = fh.readline().split()
        for line in fh:
            f = line.split()
            if f and f[0] in names:
                d = dict(zip(header, f))
                return int(d["cen.left.base"]) - 1, int(d["cen.right.base"])
    return None


def gc_per_tile(fasta_path: str | None, tiling: Tiling) -> tuple[np.ndarray, np.ndarray]:
    """GC fraction (NaN when no ACGT) and number of N bases per tile."""
    gc = np.full(len(tiling), np.nan); nn = np.zeros(len(tiling), dtype=np.int64)
    if fasta_path is None:
        return gc, nn
    import pysam
    with pysam.FastaFile(fasta_path) as fa:
        if tiling.chrom not in fa.references:   # checked before the BAM pass (cli.check_chrom_reference)
            raise SystemExit(f"{tiling.chrom} is not in the FASTA {fasta_path}: no GC content for its bins")
        for i, (s, e) in enumerate(zip(tiling.starts, tiling.ends)):
            seq = fa.fetch(tiling.chrom, int(s), int(e)).upper()
            g = seq.count("G") + seq.count("C"); at = seq.count("A") + seq.count("T")
            nn[i] = seq.count("N")
            gc[i] = g / (g + at) if (g + at) else np.nan
    return gc, nn


def annotate_tiling(tiling: Tiling, ref: dict, peaks: str | None = None):
    """Per-tile annotation frame: gc, n_N, frac_blacklist, frac_trf, is_centromere, peak_bp, width."""
    import pandas as pd
    df = tiling.frame()
    df["width"] = df["end"] - df["start"]
    gc, nn = gc_per_tile(ref.get("fasta"), tiling)
    df["gc"] = gc; df["n_N"] = nn
    bs, be = read_bed_intervals(ref.get("blacklist"), tiling.chrom)
    df["frac_blacklist"] = overlap_bp(tiling, bs, be) / df["width"].values
    ts, te = read_bed_intervals(ref.get("trf"), tiling.chrom)
    df["frac_trf"] = overlap_bp(tiling, ts, te) / df["width"].values
    cen = read_centromere(ref["centromeres"], tiling.chrom) if ref.get("centromeres") else None
    if cen is None:
        df["is_centromere"] = False
    else:
        df["is_centromere"] = (df["start"].values < cen[1]) & (df["end"].values > cen[0])
    ps, pe = read_bed_intervals(peaks, tiling.chrom)
    df["peak_bp"] = overlap_bp(tiling, ps, pe)
    df.index = [f"{tiling.chrom}:{s}-{e}" for s, e in zip(df["start"], df["end"])]
    return df
