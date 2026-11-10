"""Sample config (a sample yaml, see inst/extdata/sample_template.yaml) and resource paths."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import yaml

RESOURCES = Path(__file__).parent / "resources"
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def load_yaml(path: str | os.PathLike) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


@dataclass
class SampleConfig:
    name: str
    bam: str
    barcodes: str | None = None
    sex: str = "unknown"
    type: str = "unknown"
    chroms: list[str] = field(default_factory=lambda: list(MAIN_CHROMS))
    snps: str | None = None
    phased: bool = False
    vcf_sample: str | None = None
    peaks: str | None = None
    barcode_map: str | None = None
    out: str | None = None
    extra: dict = field(default_factory=dict)


def sample_from_config(cfg: dict, name: str, out_root: str | None = None) -> SampleConfig:
    """Resolve a sample block of the sample yaml into a SampleConfig (whitelist and VCF rules applied)."""
    s = cfg["samples"][name]
    chroms = list(s.get("chroms") or MAIN_CHROMS)
    if "chroms" in s and s["chroms"]:
        chroms = list(s["chroms"])
    root = out_root or cfg.get("outputs", {}).get("root", ".")
    sc = SampleConfig(
        name=name,
        bam=s["bam"],
        barcodes=s.get("barcodes") or s.get("whitelist"),
        sex=s.get("sex", "unknown"),
        type=s.get("type", "unknown"),
        chroms=chroms,
        peaks=s.get("peaks"),
        barcode_map=s.get("barcode_map"),
        vcf_sample=s.get("vcf_sample"),
        out=str(Path(root) / name / "extract"),
        extra={k: v for k, v in s.items()},
    )
    # SNP sources: phased longphase VCF preferred, then an unphased het list
    if s.get("phased_vcf"):
        sc.snps, sc.phased = s["phased_vcf"], True
    elif s.get("het_vcf"):
        sc.snps, sc.phased = s["het_vcf"], False
    return sc


def reference_paths(cfg: dict | None = None) -> dict:
    """Reference files: yaml `reference:` values win; chromosome sizes, blacklist and centromeres fall back to the
    hg38 files shipped in `resources/`. The FASTA (GC per bin) and the TRF bed have no fallback (None)."""
    ref = (cfg or {}).get("reference", {}) or {}
    return {
        "chrom_sizes": ref.get("chrom_sizes", str(RESOURCES / "hg38.chrom.sizes")),
        "fasta": ref.get("fasta"),
        "blacklist": ref.get("blacklist_v2", ref.get("blacklist", str(RESOURCES / "hg38-blacklist.v2.bed"))),
        "centromeres": ref.get("centromeres", str(RESOURCES / "centromeres_hg38.tsv")),
        "trf": ref.get("trf"),
    }


def check_reference(ref: dict, require_fasta: bool = True) -> list[str]:
    """Fail early (before the BAM pass) on reference paths that do not exist. Returns warnings."""
    warn = []
    for k, v in ref.items():
        if v and not Path(v).exists():
            raise SystemExit(f"reference file `{k}` not found: {v}")
    if ref.get("fasta"):
        if not Path(str(ref["fasta"]) + ".fai").exists():
            warn.append(f"no {ref['fasta']}.fai: pysam builds it next to the FASTA (needs write access)")
    elif require_fasta:
        raise SystemExit("need --fasta (or reference.fasta in the yaml): the caller needs the GC content per bin; "
                         "use --no-gc to extract without it (that output is not callable by the scphasedcn R caller)")
    return warn


def read_chrom_sizes(path: str | os.PathLike) -> dict[str, int]:
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            c, n = line.split()[:2]
            out[c] = int(n)
    return out


def read_barcodes(path: str | os.PathLike) -> list[str]:
    """One barcode per line (gz or plain); strips a trailing '-1' (the BAM lookup accepts CB with or without it);
    drops '-' and header-like lines."""
    import gzip
    op = gzip.open if str(path).endswith(".gz") else open
    bcs = []
    with op(path, "rt") as fh:
        for line in fh:
            t = line.strip().split("\t")[0]
            if not t or t == "-" or t.lower() in ("barcode", "cell", "cell_id"):
                continue
            if t.endswith("-1"):
                t = t[:-2]
            bcs.append(t)
    return sorted(set(bcs))


def translate_barcodes(bcs: list[str], barcode_map: str | None) -> tuple[list[str], str]:
    """A whitelist may be in RNA-barcode space (e.g. 10x Multiome GEX barcodes) while the BAM `CB` is the ATAC
    barcode. With a barcode map (a TSV with columns atac_barcode, rna_barcode; `--barcode-map`) whitelist entries are
    matched against both columns; when more are found in the RNA column they are translated to the ATAC column.
    Returns (barcodes, space) with space "atac" (used as is) or "rna->atac" (translated)."""
    if not barcode_map:
        return bcs, "atac"
    import pandas as pd
    m = pd.read_csv(barcode_map, sep="\t", usecols=["atac_barcode", "rna_barcode"])
    strip = lambda x: x[:-2] if x.endswith("-1") else x
    atac = {strip(x) for x in m["atac_barcode"]}
    rna2atac = {strip(r): strip(a) for a, r in zip(m["atac_barcode"], m["rna_barcode"])}
    n_atac = sum(b in atac for b in bcs); n_rna = sum(b in rna2atac for b in bcs)
    if n_atac >= n_rna:
        return bcs, "atac"
    out = sorted({rna2atac[b] for b in bcs if b in rna2atac})
    return out, "rna->atac"
