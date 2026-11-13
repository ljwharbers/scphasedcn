"""Merge per-chromosome partials into the per-sample outputs (layout: the scphasedcn-inputs vignette)."""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp

from .qc import finalize_qc
from .writers import load_csr_dict, concat_bgzf, write_json, read_done, LAYER_FILES, FRAGMENT_COLUMNS


def _read_cells(pdir: Path) -> list[str]:
    return [l.rstrip("\n") for l in open(pdir / "cells.txt") if l.strip()]


def _npz_shapes(path: Path) -> dict[str, tuple]:
    with np.load(path) as z:
        return {k[:-len("__shape")]: tuple(int(x) for x in z[k]) for k in z.files if k.endswith("__shape")}


# DONE inputs that describe the SNPs a layer was made from (fragments: haplotype votes; alleles: allele counts)
SNP_INPUTS = ("snps", "phased", "vcf_sample", "min_bq")


def _fmt(v) -> str:
    s = json.dumps(v, sort_keys=True, default=str)
    return s if len(s) <= 120 else s[:117] + "..."


def check_partials(dirs: list[Path], log) -> None:
    """Fail (SystemExit) on partials that cannot be merged into one sample: a --region partial, different cell
    lists, a layer present on some chromosomes only, layer files the DONE metadata does not list (stale), matrices whose
    shape does not fit the cells or the depth bins, or layers made from different inputs (bin size, assignment, MAPQ,
    BAM, SNP file, masks; recorded in DONE by extract), or, within a partial, a fragments table whose haplotype votes
    come from other SNP inputs than its alleles layer. Partials of scphasedcn <= 0.1.1 (no DONE metadata) are
    merged with a warning: only their cells, layer presence and shapes are checked."""
    first = dirs[0]
    cells0 = _read_cells(first)
    metas = {d.name: read_done(d) for d in dirs}
    for d in dirs:
        m = metas[d.name]
        region = m.get("region") if m else None
        if m is None and (d / "params.json").exists():
            region = json.load(open(d / "params.json")).get("region")
        if region:
            raise SystemExit(f"partial {d.name} is a --region extraction ({region}), not a whole chromosome; "
                             f"re-extract {d.name} without --region before merging")
        c = _read_cells(d)
        if c != cells0:
            raise SystemExit(f"partials {first.name} and {d.name} have different cell lists ({len(cells0)} vs {len(c)} "
                             f"cells): they were extracted with different whitelists; re-extract with one whitelist")
    legacy = [d.name for d in dirs if metas[d.name] is None]
    if legacy:
        log(f"warning: partials written by scphasedcn <= 0.1.1 (no DONE metadata): {', '.join(legacy)}; their inputs "
            f"cannot be checked (re-extract them to be sure they match)")
    # layer presence: every partial has the same layers
    has = {L: [d.name for d in dirs if (d / LAYER_FILES[L][0]).exists()] for L in LAYER_FILES}
    for L, names in has.items():
        if names and len(names) != len(dirs):
            lack = [d.name for d in dirs if d.name not in names]
            raise SystemExit(f"the {L} layer ({LAYER_FILES[L][0]}) is on {', '.join(names)} but not on "
                             f"{', '.join(lack)}: extract the same layers on every chromosome")
    for d in dirs:
        m = metas[d.name]
        if m is not None:
            for L, files in LAYER_FILES.items():
                stale = [f for f in files if (d / f).exists()]
                if stale and L not in m["layers"]:
                    raise SystemExit(f"partial {d.name}: {', '.join(stale)} not listed in its DONE (left from an earlier "
                                     f"run); delete them or re-extract {d.name} with the {L} layer")
        # shapes: rows = cells; hap layers have the depth bins
        ncol = None
        if (d / "layers.npz").exists():
            sh = _npz_shapes(d / "layers.npz")
            if any(v[0] != len(cells0) for v in sh.values()):
                raise SystemExit(f"partial {d.name}: layers.npz has a different number of cells than cells.txt "
                                 f"(stale depth layer); re-extract {d.name}")
            ncol = sh["counts"][1]
        if (d / "alleles.npz").exists():
            sh = _npz_shapes(d / "alleles.npz")
            if any(v[0] != len(cells0) for v in sh.values()):
                raise SystemExit(f"partial {d.name}: alleles.npz has a different number of cells than cells.txt "
                                 f"(stale allele layer); re-extract {d.name}")
            if ncol is not None and "hapA" in sh and sh["hapA"][1] != ncol:
                raise SystemExit(f"partial {d.name}: the hap layers of alleles.npz have {sh['hapA'][1]} bins, "
                                 f"layers.npz {ncol} (different bin sizes); re-extract {d.name}")
    # inputs per layer: the same on every partial, and within a partial the cells of the current cells.txt
    for L in LAYER_FILES:
        ref_name = ref = None
        unknown = []
        for d in dirs:
            m = metas[d.name]
            if m is None or L not in m["layers"]:
                continue
            inp = m["layers"][L].get("inputs")
            if inp is None:
                unknown.append(d.name)
                continue
            if inp.get("cells") != m.get("cells"):
                raise SystemExit(f"partial {d.name}: its {L} layer was extracted for another cell list than its "
                                 f"cells.txt (kept from an earlier run); re-extract {d.name} with the {L} layer")
            if ref is None:
                ref_name, ref = d.name, inp
                continue
            diff = sorted(k for k in set(ref) | set(inp) if ref.get(k) != inp.get(k))
            if diff:
                k = diff[0]
                raise SystemExit(f"partials disagree on the inputs of the {L} layer: {k} is {_fmt(ref.get(k))} on "
                                 f"{ref_name} but {_fmt(inp.get(k))} on {d.name}"
                                 f"{' (also ' + ', '.join(diff[1:]) + ')' if len(diff) > 1 else ''}; "
                                 f"re-extract the partials with the same inputs")
        if unknown:
            log(f"warning: the {L} layer of {', '.join(unknown)} was kept from an extraction without metadata; "
                f"its inputs cannot be checked")
    for d in dirs:   # within a partial, the allele hap layers use the depth bins
        m = metas[d.name]
        if m is None:
            continue
        b = {L: (m["layers"].get(L) or {}).get("inputs") for L in ("depth", "alleles")}
        if b["depth"] and b["alleles"] and b["depth"].get("binsize") != b["alleles"].get("binsize"):
            raise SystemExit(f"partial {d.name}: depth and alleles layers have different bin sizes "
                             f"({b['depth'].get('binsize')} vs {b['alleles'].get('binsize')}); re-extract {d.name}")
        # the caller takes haplotype votes from the fragments table and allele counts from the alleles layer: both
        # must come from the same SNPs (a layer kept from an extraction without metadata was warned about above)
        f, al = (m["layers"].get(L, {}).get("inputs") for L in ("fragments", "alleles"))
        if f and al:
            if f.get("snps") is None and al.get("snps") is not None:
                raise SystemExit(f"partial {d.name}: its fragments layer was made without SNPs (no haplotype votes) "
                                 f"but its alleles layer with SNPs; re-extract {d.name} with --layers alleles,fragments")
            diff = [k for k in SNP_INPUTS if f.get(k) != al.get(k)]
            if diff:
                raise SystemExit(f"partial {d.name}: its fragments layer (haplotype votes) and its alleles layer were "
                                 f"made from different SNP inputs ({', '.join(diff)}); re-extract {d.name} with "
                                 f"--layers alleles,fragments")


def _stack_layers(dirs: list[Path], fname: str) -> dict[str, sp.csr_matrix]:
    """hstack per-chromosome CSR layers; a layer missing on some chromosome is zero there."""
    per = [load_csr_dict(d / fname) for d in dirs]
    names = []
    for p in per:
        names += [k for k in p if k not in names]
    out = {}
    for k in names:
        mats = []
        for p in per:
            if k in p:
                mats.append(p[k])
            else:
                ref = next(iter(p.values()))
                mats.append(sp.csr_matrix(ref.shape, dtype=ref.dtype))
        out[k] = sp.hstack(mats, format="csr")
    return out


def _sum_obsm(dirs: list[Path], fname: str) -> dict[str, np.ndarray]:
    out = {}
    for d in dirs:
        z = np.load(d / fname)
        for k in z.files:
            out[k] = out.get(k, 0) + z[k]
    return out


def _bins(dirs: list[Path], fname: str) -> pd.DataFrame:
    return pd.concat([pd.read_csv(d / fname, sep="\t", index_col=0) for d in dirs])


def _qc(dirs: list[Path], fname: str) -> tuple[pd.DataFrame, dict]:
    per_chrom_pass = {}
    tot = None
    for d in dirs:
        q = pd.read_csv(d / fname, sep="\t", index_col=0)
        per_chrom_pass[d.name] = q["n_pass"].values if "n_pass" in q else None
        tot = q if tot is None else tot.add(q, fill_value=0)
    return tot, per_chrom_pass


def _snps(dirs: list[Path]) -> pd.DataFrame:
    df = pd.concat([pd.read_csv(d / "snps.tsv", sep="\t", index_col=0) for d in dirs])
    for c, t in (("pos", np.int64), ("hap_of_alt", np.int8), ("ps", np.int64), ("n_other", np.int64), ("n_lowbq", np.int64)):
        if c in df:
            df[c] = df[c].astype(t)
    for c in ("chrom", "ref", "alt"):
        df[c] = df[c].astype(str)
    df.index = df.index.astype(str)
    return df


def _uns_safe(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out[k] = _uns_safe(v)
        elif v is None:
            out[k] = "None"
        elif isinstance(v, (list, tuple)):
            out[k] = [str(x) if x is None else x for x in v]
        else:
            out[k] = v
    return out


def _write_via_tmp(dest: Path, write) -> None:
    """Write `dest` through a temporary file in the same directory (`write(path)`), then move it into place: an
    interrupted merge never leaves a cut-off output under its final name. The temporary file is removed on failure."""
    tmp = dest.with_name(f".tmp.{dest.name}")
    try:
        write(tmp)
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)


def merge_sample(out: str, chroms: list[str], tag: str = "100kb", keep_partials: bool = True,
                 allow_missing: bool = False, log=None) -> dict:
    """Merge the finished partials of `chroms`. A chromosome without a finished partial (partials/<chrom>/DONE) is an
    error unless `allow_missing` (then it is listed as `missing` in summary.json). Partials are checked with
    check_partials() first. Every output is written to a temporary file and moved into place; summary.json is removed
    at the start and written last, so it is present only after a complete merge."""
    t0 = time.time()
    log = log or (lambda *a: print(*a, file=sys.stderr, flush=True))
    out = Path(out); P = out / "partials"
    dirs = [P / c for c in chroms if (P / c / "DONE").exists()]
    missing = [c for c in chroms if not (P / c / "DONE").exists()]
    if not dirs:
        raise SystemExit(f"no finished partials under {P}")
    if missing and not allow_missing:
        raise SystemExit(f"merge: no finished partial (partials/<chrom>/DONE) for {', '.join(missing)}; extract them, "
                         f"or merge without them with --allow-missing (or leave them out of --chroms)")
    if missing:
        log(f"warning: merging without {', '.join(missing)} (--allow-missing)")
    check_partials(dirs, log)
    (out / "summary.json").unlink(missing_ok=True)
    first = dirs[0]
    cells = _read_cells(first)
    params = json.load(open(first / "params.json"))
    summary = {"chroms": [d.name for d in dirs], "missing": missing, "n_cells": len(cells), "params": params}
    # --- bins (depth layers) ---
    if (dirs[0] / "layers.npz").exists():
        ldirs = [d for d in dirs if (d / "layers.npz").exists()]
        layers = _stack_layers(ldirs, "layers.npz")
        if all((d / "alleles.npz").exists() for d in ldirs):
            al = _stack_layers(ldirs, "alleles.npz")
            for k in ("hapA", "hapB", "hap_ambig"):
                if k in al:
                    layers[k] = al[k]
        obsm = _sum_obsm(ldirs, "obsm.npz")
        bins = _bins(ldirs, "bins.tsv")
        qc, per_chrom_pass = _qc(ldirs, "qc.tsv")
        if all((d / "alleles_qc.tsv").exists() for d in ldirs):
            aq, _ = _qc(ldirs, "alleles_qc.tsv")
            for c in aq.columns:
                qc[c] = aq[c]
        qc = finalize_qc(qc, per_chrom_pass, peaks_given=bool(params.get("peaks")))
        qc = qc.loc[cells]
        X = layers.pop("counts")
        a = ad.AnnData(X=X, obs=qc, var=bins, layers=layers, obsm={k: v.astype(np.int32) for k, v in obsm.items()})
        a.uns["params"] = _uns_safe(params)
        a.uns["obsm_layout"] = _uns_safe(obsm_layout())
        _write_via_tmp(out / f"bins_{tag}.h5ad", lambda p: a.write_h5ad(p, compression="gzip"))
        _write_via_tmp(out / "cells_qc.tsv", lambda p: qc.to_csv(p, sep="\t"))
        summary["bins"] = {"n_bins": int(a.n_vars), "total_counts": int(X.sum()), "layers": list(layers.keys())}
        for k in ("n_pass", "n_clean"):
            if k in qc:
                summary[f"median_{k}"] = float(np.median(qc[k]))
    # --- alleles (DNA) ---
    adirs = [d for d in dirs if (d / "alleles.npz").exists()]
    if adirs:
        al = _stack_layers(adirs, "alleles.npz")
        snps = _snps(adirs)
        a = ad.AnnData(X=al["DP"], obs=pd.DataFrame(index=cells), var=snps, layers={"AD": al["AD"]})
        a.uns["params"] = _uns_safe(params)
        _write_via_tmp(out / "alleles_dna.h5ad", lambda p: a.write_h5ad(p, compression="gzip"))
        summary["alleles_dna"] = {"n_snps": int(a.n_vars), "n_phased": int((snps["hap_of_alt"] > 0).sum()),
                                  "total_dp": int(al["DP"].sum()), "total_ad": int(al["AD"].sum())}
    # --- fragments table ---
    parts = [str(d / "fragments.tsv.gz") for d in dirs if (d / "fragments.tsv.gz").exists()]
    if parts:
        fn = out / "fragments.tsv.gz"
        tmp = fn.with_name(f".tmp.{fn.name}")
        tmp_idx = Path(f"{tmp}.tbi")
        try:   # table and index are built under a temporary name, then moved: data first, then its own index
            concat_bgzf(parts, str(tmp), header=FRAGMENT_COLUMNS)
            pysam.tabix_index(str(tmp), seq_col=0, start_col=1, end_col=2, zerobased=True, force=True, meta_char="#")
            Path(f"{fn}.tbi").unlink(missing_ok=True)   # never the new table next to the old index
            os.replace(tmp, fn)
            os.replace(tmp_idx, f"{fn}.tbi")
        finally:
            for f in (tmp, tmp_idx, Path(f"{tmp}.hdr")):
                f.unlink(missing_ok=True)
        summary["fragments"] = str(fn)
    summary["merge_seconds"] = round(time.time() - t0, 1)
    _write_via_tmp(out / "summary.json", lambda p: write_json(p, summary))
    return summary


def obsm_layout() -> dict:
    from .depth import FRAGLEN_LABELS
    return {
        "fraglen_hist": FRAGLEN_LABELS,
        "mapq_hist": ["mapq0", "mapq1_29", "mapq30_59", "mapq60"],
    }
