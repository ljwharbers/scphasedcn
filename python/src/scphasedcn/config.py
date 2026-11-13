"""Sample config (a sample yaml, see inst/extdata/sample_template.yaml) and resource paths."""
from __future__ import annotations

import os
import sys
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
    chroms_set: bool = False   # True when the yaml sample gives `chroms` (else `chroms` is the default list)
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
        chroms_set=bool(s.get("chroms")),
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


MIN_MAP_MATCH = 0.9   # warn when the chosen barcode space holds fewer of the whitelist barcodes than this


def translate_barcodes(bcs: list[str], barcode_map: str | None, log=None) -> tuple[list[str], str]:
    """A whitelist may be in RNA-barcode space (e.g. 10x Multiome GEX barcodes) while the BAM `CB` is the ATAC
    barcode. With a barcode map (a TSV with columns atac_barcode, rna_barcode; `--barcode-map`) whitelist entries are
    matched against both columns; when more are found in the RNA column they are translated to the ATAC column
    (whitelist barcodes not in the map are then dropped). Returns (barcodes, space) with space "atac" (used as is) or
    "rna->atac" (translated). Both match counts are logged; a map matching no whitelist barcode is an error, and a tie
    or a chosen column holding fewer than MIN_MAP_MATCH of the whitelist is warned about."""
    if not barcode_map:
        return bcs, "atac"
    log = log or (lambda *a: print(*a, file=sys.stderr, flush=True))
    import pandas as pd
    m = pd.read_csv(barcode_map, sep="\t", usecols=["atac_barcode", "rna_barcode"])
    strip = lambda x: x[:-2] if x.endswith("-1") else x
    atac = {strip(x) for x in m["atac_barcode"]}
    rna2atac = {strip(r): strip(a) for a, r in zip(m["atac_barcode"], m["rna_barcode"])}
    n_atac = sum(b in atac for b in bcs); n_rna = sum(b in rna2atac for b in bcs)
    n = len(bcs)
    if n_atac == 0 and n_rna == 0:
        raise SystemExit(f"--barcode-map {barcode_map}: none of the {n} whitelist barcodes is in its atac_barcode or "
                         f"rna_barcode column; is it the map of this sample's whitelist?")
    space = "atac" if n_atac >= n_rna else "rna->atac"
    n_hit = n_atac if space == "atac" else n_rna
    log(f"barcode map: {n_atac} of {n} whitelist barcodes in its atac_barcode column, {n_rna} of {n} in its "
        f"rna_barcode column: " + ("whitelist used as is (ATAC barcodes)" if space == "atac" else
                                   f"translated to ATAC barcodes; {n - n_rna} whitelist barcodes not in the map "
                                   f"are dropped"))
    if n_atac == n_rna:
        log(f"warning: barcode map: a tie ({n_atac} whitelist barcodes in each column); the whitelist is used as is "
            f"(ATAC barcodes): check that it is the map of this whitelist")
    if n_hit < MIN_MAP_MATCH * n:
        log(f"warning: barcode map: only {n_hit} of {n} whitelist barcodes ({100 * n_hit / n:.0f}%) are in its "
            f"{'atac' if space == 'atac' else 'rna'}_barcode column: check that it is the map of this whitelist")
    if space == "atac":
        return bcs, "atac"
    out = sorted({rna2atac[b] for b in bcs if b in rna2atac})
    return out, "rna->atac"
