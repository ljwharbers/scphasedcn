"""Partial (per-chromosome) outputs and bgzipped tables."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp

# files of each layer of a partial (partials/<chrom>/); cells.txt and params.json belong to every partial
LAYER_FILES = {"depth": ["layers.npz", "obsm.npz", "bins.tsv", "qc.tsv"],
               "alleles": ["alleles.npz", "snps.tsv", "alleles_qc.tsv"],
               "fragments": ["fragments.tsv.gz"]}
DONE_FORMAT = 1

FRAGMENT_COLUMNS = ["chrom", "start", "end", "cell", "strand", "mapq", "flags", "len", "piece", "npieces", "hap", "nsnp", "ps", "read_id"]


def save_csr_dict(path: str | os.PathLike, mats: dict[str, sp.csr_matrix]) -> None:
    d = {}
    for k, m in mats.items():
        m = sp.csr_matrix(m)
        d[f"{k}__data"] = m.data; d[f"{k}__indices"] = m.indices; d[f"{k}__indptr"] = m.indptr
        d[f"{k}__shape"] = np.asarray(m.shape, dtype=np.int64)
    np.savez_compressed(path, **d)


def load_csr_dict(path: str | os.PathLike) -> dict[str, sp.csr_matrix]:
    z = np.load(path)
    names = sorted({k.split("__")[0] for k in z.files})
    return {k: sp.csr_matrix((z[f"{k}__data"], z[f"{k}__indices"], z[f"{k}__indptr"]), shape=tuple(z[f"{k}__shape"])) for k in names}


def write_bgzf_table(path: str | os.PathLike, df: pd.DataFrame, header: bool = False, chunk: int = 1_000_000) -> None:
    """Partial tables carry no header; merge.concat writes one header block in front of the raw BGZF parts."""
    with pysam.BGZFile(str(path), "wb") as fh:
        if header:
            fh.write(("#" + "\t".join(df.columns) + "\n").encode())
        for i in range(0, max(len(df), 1), chunk):
            part = df.iloc[i:i + chunk]
            if len(part):
                fh.write(part.to_csv(sep="\t", header=False, index=False).encode())


def fragments_frame(fr, cells: list[str]) -> pd.DataFrame:
    bc = np.asarray(cells, dtype=object)
    return pd.DataFrame({
        "chrom": fr.chrom, "start": fr.start, "end": fr.end, "cell": bc[fr.cell],
        "strand": np.where(fr.strand == 0, "+", "-"), "mapq": fr.mapq, "flags": fr.flags, "len": fr.end - fr.start,
        "piece": fr.piece, "npieces": fr.npieces, "hap": fr.hap, "nsnp": fr.nsnp, "ps": fr.ps, "read_id": fr.read_id,
    }, columns=FRAGMENT_COLUMNS)


def concat_bgzf(parts: list[str], out: str, header: list[str] | None = None) -> None:
    """Concatenate header-less BGZF files (a concatenation of BGZF files is a valid BGZF file; the EOF marker is an
    empty block). `header` column names are written as one leading '#...' line in their own BGZF block."""
    with open(out, "wb") as oh:
        if header:
            with pysam.BGZFile(out + ".hdr", "wb") as hh:
                hh.write(("#" + "\t".join(header) + "\n").encode())
            with open(out + ".hdr", "rb") as ih:
                oh.write(ih.read())
            os.remove(out + ".hdr")
        for p in parts:
            with open(p, "rb") as ih:
                for blk in iter(lambda: ih.read(1 << 24), b""):
                    oh.write(blk)


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def file_fingerprint(path) -> dict | None:
    """Content identity of an input file that does not depend on its path (size + sha256 of the first and last
    64 kB), recorded per partial so that `merge` can tell partials made from different inputs apart."""
    if not path:
        return None
    p = str(path)
    n = os.path.getsize(p)
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        h.update(fh.read(1 << 16))
        if n > (1 << 17):
            fh.seek(n - (1 << 16))
            h.update(fh.read(1 << 16))
    return {"size": n, "sha256_ends": h.hexdigest()}


def cells_digest(cells: list[str]) -> str:
    return hashlib.sha256(("\n".join(cells) + "\n").encode()).hexdigest()


def write_done(pdir, meta: dict) -> None:
    """Mark a partial finished: DONE holds the partial's metadata (JSON), written atomically."""
    tmp = Path(pdir) / "DONE.tmp"
    write_json(tmp, meta)
    os.replace(tmp, Path(pdir) / "DONE")


def read_done(pdir) -> dict | None:
    """The metadata in a partial's DONE, or None when DONE is missing or has none (written by scphasedcn <= 0.1.1)."""
    f = Path(pdir) / "DONE"
    try:
        meta = json.loads(f.read_text())
    except (OSError, ValueError):
        return None
    return meta if isinstance(meta, dict) and "layers" in meta else None
