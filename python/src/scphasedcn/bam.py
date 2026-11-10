"""One pass over a chromosome of a SPLONGGET DNA BAM: record filter, fragment arrays, per-cell record counters.

Filter (see vignette("scphasedcn-inputs", package = "scphasedcn")): primary, not duplicate, MAPQ >= min_mapq, CB in the whitelist; "-" (unassigned)
is counted and dropped; secondaries are counted and skipped; supplementaries are counted and skipped but the primary's
SA tag sets `has_SA`. Blacklist / TRF / centromere / peak membership is flagged after the loop (bins.py). Fragments
whose start is in the blacklist are not "clean" (cli.py: `clean = (flags & F_BLACKLIST) == 0`): they are left out of the
depth layers, the fragment-length histogram and the per-bin hapA / hapB / hap_ambig layers, and kept in the fragments
table with the flag bit set. The allele pass runs inside the loop, before the flags are set, so their SNP observations
are in the allele counts (AD / DP) and their `hap` / `nsnp` / `ps` are in the fragments table (and in the per-cell
n_hapA / n_hapB / n_snp_obs QC) like those of any other fragment; the R caller's `votes_clean` option drops them from
the molecule votes.
"""
from __future__ import annotations

import time
from array import array
from dataclasses import dataclass, field

import numpy as np
import pysam

from .names import parse_name

# fragment flag bits
F_HAS_SA = 1
F_MAPQ60 = 2
F_PIECE_MULTI = 4
F_BLACKLIST = 8                  # start in the blacklist; "clean" = this bit unset (depth layers; R votes_clean)
F_TRF = 16
F_CENTROMERE = 32
F_PEAK = 64
FLAG_BITS = {"has_SA": F_HAS_SA, "mapq60": F_MAPQ60, "piece_multi": F_PIECE_MULTI, "blacklist": F_BLACKLIST,
             "trf": F_TRF, "centromere": F_CENTROMERE, "peak": F_PEAK}

MAPQ_BUCKETS = ["mapq0", "mapq1_29", "mapq30_59", "mapq60"]
LOOP_COUNTERS = ["n_records", "n_secondary", "n_supplementary", "n_duplicate", "n_primary", "n_mapq_fail", "n_pass"]


def mapq_bucket(mq: int) -> int:
    return 0 if mq == 0 else (1 if mq < 30 else (2 if mq < 60 else 3))


@dataclass
class ChromFragments:
    """Fragments that passed the record filter on one chromosome (arrays aligned by index)."""

    chrom: str
    length: int
    cell: np.ndarray            # int32 index into the whitelist
    start: np.ndarray           # int64, 0-based
    end: np.ndarray             # int64, exclusive
    strand: np.ndarray          # int8: 0 forward, 1 reverse
    mapq: np.ndarray            # uint8
    flags: np.ndarray           # uint16 bitfield (FLAG_BITS)
    piece: np.ndarray           # int16 piece index within the raw read
    npieces: np.ndarray         # int16
    read_id: np.ndarray         # uint64 (first 64 bits of the raw-read uuid)
    hap: np.ndarray             # int8: 0 none, 1 hapA, 2 hapB, 3 ambiguous
    nsnp: np.ndarray            # int16 SNP observations on the read
    ps: np.ndarray              # int64 phase set (-1)

    def __len__(self):
        return len(self.start)

    @property
    def length_bp(self):
        return self.end - self.start


@dataclass
class LoopCounters:
    n_cells: int
    per_cell: np.ndarray                       # n_cells x len(LOOP_COUNTERS)
    mapq_hist: np.ndarray                      # n_cells x 4 (primary, non-dup records)
    dash: dict = field(default_factory=lambda: {k: 0 for k in LOOP_COUNTERS})
    unlisted: dict = field(default_factory=lambda: {k: 0 for k in LOOP_COUNTERS})
    n_no_cb: int = 0
    n_unmapped: int = 0
    seconds: float = 0.0


def _count(d: dict, flag: int, mq: int, min_mapq: int):
    d["n_records"] += 1
    if flag & 0x100:
        d["n_secondary"] += 1
    elif flag & 0x800:
        d["n_supplementary"] += 1
    elif flag & 0x400:
        d["n_duplicate"] += 1
    else:
        d["n_primary"] += 1
        if mq < min_mapq:
            d["n_mapq_fail"] += 1
        else:
            d["n_pass"] += 1


def scan_chrom(bam_path: str, chrom: str, cells: list[str], min_mapq: int = 30, index_filename: str | None = None,
               region: tuple[int, int] | None = None, alleles=None,
               cb_tag: str = "CB") -> tuple[ChromFragments, LoopCounters]:
    """Scan one chromosome (or a region: fragments with start in [rstart, rend)) and return fragment arrays + counters.

    `alleles` is an AlleleAccumulator (alleles.py) or None; it is called for every passing fragment and returns
    (hap, nsnp, ps).
    """
    t0 = time.time()
    n_cells = len(cells)
    cidx = {c: i for i, c in enumerate(cells)}
    cnt = np.zeros((n_cells, len(LOOP_COUNTERS)), dtype=np.int64)
    mqh = np.zeros((n_cells, 4), dtype=np.int64)
    ctr = LoopCounters(n_cells, cnt, mqh)
    I_REC, I_SEC, I_SUP, I_DUP, I_PRI, I_MQF, I_PASS = range(7)

    a_cell, a_start, a_end, a_strand = array("i"), array("q"), array("q"), array("b")
    a_mapq, a_flags, a_piece, a_np = array("B"), array("H"), array("h"), array("h")
    a_rid, a_hap, a_nsnp, a_ps = array("Q"), array("b"), array("h"), array("q")

    with pysam.AlignmentFile(bam_path, "rb", index_filename=index_filename) as bam:
        length = bam.get_reference_length(chrom)
        rstart, rend = (region if region else (None, None))
        it = bam.fetch(chrom, rstart, rend) if region else bam.fetch(chrom)
        for r in it:
            flag = r.flag
            if flag & 0x4:
                ctr.n_unmapped += 1
                continue
            try:
                cb = r.get_tag(cb_tag)
            except KeyError:
                ctr.n_no_cb += 1
                continue
            mq = r.mapping_quality
            if cb == "-":
                _count(ctr.dash, flag, mq, min_mapq)
                continue
            ci = cidx.get(cb)
            if ci is None:
                _count(ctr.unlisted, flag, mq, min_mapq)
                continue
            row = cnt[ci]
            row[I_REC] += 1
            if flag & 0x100:
                row[I_SEC] += 1
                continue
            if flag & 0x800:
                row[I_SUP] += 1
                continue
            if flag & 0x400:
                row[I_DUP] += 1
                continue
            row[I_PRI] += 1
            mqh[ci, mapq_bucket(mq)] += 1
            if mq < min_mapq:
                row[I_MQF] += 1
                continue
            start = r.reference_start
            if region is not None and not (rstart <= start < rend):
                continue
            end = r.reference_end
            row[I_PASS] += 1
            _, rid, piece, npieces, _sign = parse_name(r.query_name)
            fl = 0
            if r.has_tag("SA"):
                fl |= F_HAS_SA
            if mq >= 60:
                fl |= F_MAPQ60
            if npieces > 1:
                fl |= F_PIECE_MULTI
            if alleles is not None:
                hap, nsnp, ps = alleles.process(r, ci, start, end)
            else:
                hap, nsnp, ps = 0, 0, -1
            a_cell.append(ci); a_start.append(start); a_end.append(end); a_strand.append(1 if flag & 0x10 else 0)
            a_mapq.append(mq if mq < 256 else 255); a_flags.append(fl); a_piece.append(piece); a_np.append(npieces)
            a_rid.append(rid); a_hap.append(hap); a_nsnp.append(nsnp); a_ps.append(ps)

    ctr.seconds = time.time() - t0
    fr = ChromFragments(
        chrom=chrom, length=length,
        cell=np.frombuffer(a_cell, dtype=np.int32).copy() if len(a_cell) else np.zeros(0, np.int32),
        start=np.frombuffer(a_start, dtype=np.int64).copy() if len(a_start) else np.zeros(0, np.int64),
        end=np.frombuffer(a_end, dtype=np.int64).copy() if len(a_end) else np.zeros(0, np.int64),
        strand=np.frombuffer(a_strand, dtype=np.int8).copy() if len(a_strand) else np.zeros(0, np.int8),
        mapq=np.frombuffer(a_mapq, dtype=np.uint8).copy() if len(a_mapq) else np.zeros(0, np.uint8),
        flags=np.frombuffer(a_flags, dtype=np.uint16).copy() if len(a_flags) else np.zeros(0, np.uint16),
        piece=np.frombuffer(a_piece, dtype=np.int16).copy() if len(a_piece) else np.zeros(0, np.int16),
        npieces=np.frombuffer(a_np, dtype=np.int16).copy() if len(a_np) else np.zeros(0, np.int16),
        read_id=np.frombuffer(a_rid, dtype=np.uint64).copy() if len(a_rid) else np.zeros(0, np.uint64),
        hap=np.frombuffer(a_hap, dtype=np.int8).copy() if len(a_hap) else np.zeros(0, np.int8),
        nsnp=np.frombuffer(a_nsnp, dtype=np.int16).copy() if len(a_nsnp) else np.zeros(0, np.int16),
        ps=np.frombuffer(a_ps, dtype=np.int64).copy() if len(a_ps) else np.zeros(0, np.int64),
    )
    return fr, ctr


def flag_fragments(fr: ChromFragments, blacklist, trf, centromere, peaks) -> None:
    """Set BLACKLIST (start in blacklist), TRF (start in TRF), CENTROMERE (overlap) and PEAK (overlap) flag bits."""
    from .bins import in_intervals
    if len(fr) == 0:
        return
    fl = fr.flags.astype(np.uint16)
    bs, be = blacklist
    fl |= (in_intervals(fr.start, bs, be) * F_BLACKLIST).astype(np.uint16)
    ts, te = trf
    fl |= (in_intervals(fr.start, ts, te) * F_TRF).astype(np.uint16)
    if centromere is not None:
        cs, ce = centromere
        fl |= (((fr.start < ce) & (fr.end > cs)) * F_CENTROMERE).astype(np.uint16)
    ps, pe = peaks
    if len(ps):
        lo = np.searchsorted(pe, fr.start, side="right")
        hi = np.searchsorted(ps, fr.end, side="left")
        fl |= ((hi > lo) * F_PEAK).astype(np.uint16)
    fr.flags = fl


def check_index(bam_path: str, index_filename: str | None = None) -> dict:
    """Report index presence and staleness (index older than the BAM)."""
    import os
    bai = index_filename or (bam_path + ".bai" if os.path.exists(bam_path + ".bai") else bam_path[:-4] + ".bai")
    info = {"bam": bam_path, "bai": bai, "exists": os.path.exists(bai)}
    if info["exists"]:
        info["stale"] = os.path.getmtime(bai) < os.path.getmtime(bam_path)
        info["bai_mtime"] = os.path.getmtime(bai)
    info["bam_mtime"] = os.path.getmtime(bam_path)
    return info


def build_index(bam_path: str, out_bai: str, threads: int = 8) -> str:
    """Build a .bai with samtools index -o (htslib >= 1.16) into out_bai (a separate path, never next to the input BAM)."""
    pysam.samtools.index("-@", str(threads), "-o", out_bai, bam_path)
    return out_bai
