"""Per-cell QC columns for one chromosome (merged by summation in merge.py)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .bam import ChromFragments, LoopCounters, LOOP_COUNTERS, MAPQ_BUCKETS, FLAG_BITS, F_BLACKLIST, F_PEAK

SUMMED = None  # all numeric columns are summed over chromosomes


def chrom_qc(fr: ChromFragments, ctr: LoopCounters, cells: list[str], hap_counts: bool = True) -> pd.DataFrame:
    n = len(cells)
    df = pd.DataFrame(ctr.per_cell, columns=LOOP_COUNTERS, index=cells)
    for k, name in enumerate(MAPQ_BUCKETS):
        df[name] = ctr.mapq_hist[:, k]
    c = fr.cell.astype(np.int64)
    bc = lambda sel: np.bincount(c[sel], minlength=n)
    clean = (fr.flags & F_BLACKLIST) == 0
    df["n_clean"] = bc(clean)
    for name, bit in FLAG_BITS.items():
        df[f"n_{name}"] = bc((fr.flags & bit) > 0)
    L = fr.end - fr.start
    df["n_short_lt300"] = bc(clean & (L < 300))
    df["n_long_ge1000"] = bc(clean & (L >= 1000))
    df["span_sum"] = np.bincount(c[clean], weights=L[clean], minlength=n).astype(np.int64)
    # FRiP numerator: peak fragments among the clean ones (n_peak counts all passing fragments); used by finalize_qc
    # and not written to the merged QC table
    df["n_peak_clean"] = bc(clean & ((fr.flags & F_PEAK) > 0))
    if hap_counts:
        for code, name in ((1, "n_hapA"), (2, "n_hapB"), (3, "n_hap_ambig")):
            df[name] = bc(fr.hap == code)
        df["n_snp_obs"] = np.bincount(c, weights=fr.nsnp.astype(np.float64), minlength=n).astype(np.int64)
    df.index.name = "cell"
    return df


def finalize_qc(qc: pd.DataFrame, per_chrom_pass: dict[str, np.ndarray], peaks_given: bool) -> pd.DataFrame:
    """Derived columns after summing over chromosomes."""
    qc = qc.copy()
    for ch in ("chrM", "chrX", "chrY"):
        qc[f"n_pass_{ch}"] = per_chrom_pass.get(ch, np.zeros(len(qc), dtype=np.int64))
    tot = qc["n_primary"] + qc["n_duplicate"]
    qc["frac_duplicate"] = np.where(tot > 0, qc["n_duplicate"] / np.maximum(tot, 1), np.nan)
    qc["frac_mapq60"] = np.where(qc["n_pass"] > 0, qc["n_mapq60"] / np.maximum(qc["n_pass"], 1), np.nan)
    qc["frac_has_SA"] = np.where(qc["n_pass"] > 0, qc["n_has_SA"] / np.maximum(qc["n_pass"], 1), np.nan)
    qc["frac_piece_multi"] = np.where(qc["n_pass"] > 0, qc["n_piece_multi"] / np.maximum(qc["n_pass"], 1), np.nan)
    qc["frac_blacklist"] = np.where(qc["n_pass"] > 0, qc["n_blacklist"] / np.maximum(qc["n_pass"], 1), np.nan)
    qc["mean_span"] = np.where(qc["n_clean"] > 0, qc["span_sum"] / np.maximum(qc["n_clean"], 1), np.nan)
    qc["frac_chrM"] = np.where(qc["n_pass"] > 0, qc["n_pass_chrM"] / np.maximum(qc["n_pass"], 1), np.nan)
    # FRiP = clean peak fragments / clean fragments (partials without n_peak_clean, written before it: n_peak)
    n_peak_clean = qc.pop("n_peak_clean") if "n_peak_clean" in qc else qc["n_peak"]
    if peaks_given:
        qc["frip"] = np.where(qc["n_clean"] > 0, n_peak_clean / np.maximum(qc["n_clean"], 1), np.nan)
    return qc
