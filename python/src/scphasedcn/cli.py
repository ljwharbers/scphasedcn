"""scphasedcn command line: run | extract | merge | select-hets | hybrid-phase | index-check.

`scphasedcn run` extracts one sample (all chromosomes, then merge) into the directory the scphasedcn R package reads; `extract` (one
chromosome) + `merge` do the same as separate steps, e.g. one cluster job per chromosome."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from . import __version__
from .config import (MAIN_CHROMS, load_yaml, sample_from_config, reference_paths, check_reference, read_chrom_sizes,
                     read_barcodes, translate_barcodes)

LAYER_CHOICES = ("depth", "alleles", "fragments", "qc")
DEFAULT_LAYERS = "depth,alleles,fragments,qc"


def _layers(s: str) -> set[str]:
    """--layers value -> set of layer names (unknown names are an error)."""
    on = {x for x in s.split(",") if x}
    bad = sorted(on - set(LAYER_CHOICES))
    if bad:
        raise SystemExit(f"--layers: unknown layer(s) {', '.join(bad)}; choose from {','.join(LAYER_CHOICES)}")
    return on


def _log(*a):
    print(time.strftime("[%Y-%m-%d %H:%M:%S]"), *a, file=sys.stderr, flush=True)


def _region(s: str | None):
    if not s:
        return None
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, (int(a), int(b))


# ----------------------------------------------------------------------------------------------- index-check
def cmd_index_check(a):
    from .bam import check_index, build_index
    info = check_index(a.bam, a.index)
    _log(json.dumps(info))
    if a.reindex and (not info["exists"] or info.get("stale") or a.force):
        out = a.reindex
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        build_index(a.bam, out, threads=a.threads)
        _log(f"index written to {out} in {time.time() - t0:.0f} s")
    elif a.reindex:
        _log("index is fresh; nothing to do (use --force to rebuild)")
    return 0


# ----------------------------------------------------------------------------------------------- extract
def resolve_inputs(a):
    """Merge yaml sample config and command-line overrides into one dict."""
    cfg = load_yaml(a.config) if a.config else {}
    ref = reference_paths(cfg)
    d = {"sample": a.sample, "bam": a.bam, "barcodes": a.barcodes, "snps": a.snps, "phased": a.phased,
         "peaks": a.peaks, "out": a.out, "chroms": None, "barcode_map": getattr(a, "barcode_map", None)}
    if a.config and a.sample:
        sc = sample_from_config(cfg, a.sample)
        for k in ("bam", "barcodes", "snps", "peaks", "out", "barcode_map"):
            if d[k] is None:
                d[k] = getattr(sc, k)
        if not a.phased and sc.snps and d["snps"] == sc.snps:
            d["phased"] = sc.phased
        d["chroms"] = sc.chroms
    if getattr(a, "no_snps", False):
        d["snps"] = None; d["phased"] = False
    for k in ("fasta", "blacklist", "trf", "centromeres", "chrom_sizes"):
        v = getattr(a, k, None)
        if v:
            ref[k] = v
    if getattr(a, "no_gc", False):
        ref["fasta"] = None
    if getattr(a, "no_trf", False):
        ref["trf"] = None
    if getattr(a, "no_masks", False):
        ref["blacklist"] = None; ref["trf"] = None; ref["centromeres"] = None
    for w in check_reference(ref, require_fasta=not getattr(a, "no_gc", False)):
        _log("warning:", w)
    if getattr(a, "no_gc", False):
        _log("warning: --no-gc: bins have no GC content; the scphasedcn R caller cannot call this extraction")
    d["ref"] = ref
    if d["bam"] is None or d["out"] is None:
        raise SystemExit("need --bam and --out (or --config with --sample)")
    if d["barcodes"] is None:
        raise SystemExit("need --barcodes (whitelist) or a sample with `barcodes` in the yaml")
    for k in ("bam", "barcodes", "snps", "peaks", "barcode_map"):
        if d[k] and not Path(d[k]).exists():
            raise SystemExit(f"--{k.replace('_', '-')} not found: {d[k]}")
    if d["peaks"] and not Path(d["peaks"]).is_file():
        raise SystemExit(f"--peaks must be a BED file, not a directory: {d['peaks']}")
    return d


def cmd_extract(a):
    from .bam import scan_chrom, flag_fragments, F_BLACKLIST
    from .bins import Tiling, annotate_tiling, read_bed_intervals, read_centromere
    from .depth import depth_layers, fraglen_hist
    from .qc import chrom_qc
    from .writers import save_csr_dict, write_bgzf_table, fragments_frame, write_json
    from .alleles import load_snps, AlleleAccumulator

    t0 = time.time()
    layers_on = _layers(a.layers)
    d = resolve_inputs(a)
    reg = _region(a.region)
    chrom = reg[0] if reg else a.chrom
    if chrom is None:
        raise SystemExit("need --chrom or --region")
    cells, space = translate_barcodes(read_barcodes(d["barcodes"]), d["barcode_map"])
    if space != "atac":
        _log(f"whitelist translated {space} via {d['barcode_map']}: {len(cells)} cells")
    sizes = read_chrom_sizes(d["ref"]["chrom_sizes"])
    if chrom not in sizes:
        import pysam
        with pysam.AlignmentFile(d["bam"], "rb", index_filename=a.index) as bf:
            sizes[chrom] = bf.get_reference_length(chrom)
    length = sizes[chrom]
    pdir = Path(d["out"]) / "partials" / chrom
    pdir.mkdir(parents=True, exist_ok=True)
    (pdir / "DONE").unlink(missing_ok=True)
    _log(f"scphasedcn {__version__} extract {d['sample']} {chrom} ({length} bp) cells={len(cells)} layers={sorted(layers_on)}")

    # SNPs
    acc = None
    if "alleles" in layers_on and d["snps"]:
        snps = load_snps(d["snps"], chrom, phased=d["phased"], region=reg[1] if reg else None)
        acc = AlleleAccumulator(snps, min_bq=a.min_bq)
        _log(f"SNPs on {chrom}: {len(snps)} ({int((snps.hap_of_alt > 0).sum())} phased) from {d['snps']}")

    # BAM pass
    fr, ctr = scan_chrom(d["bam"], chrom, cells, min_mapq=a.min_mapq, index_filename=a.index,
                         region=reg[1] if reg else None, alleles=acc)
    _log(f"scan: {int(ctr.per_cell[:, 0].sum())} records of listed cells, {len(fr)} passing fragments, "
         f"'-' records {ctr.dash['n_records']} (pass {ctr.dash['n_pass']}), unlisted {ctr.unlisted['n_records']}, "
         f"{ctr.seconds:.0f} s")

    # masks
    ref = d["ref"]
    blacklist = read_bed_intervals(ref.get("blacklist"), chrom)
    trf = read_bed_intervals(ref.get("trf"), chrom)
    cen = read_centromere(ref["centromeres"], chrom) if ref.get("centromeres") else None
    peaks = read_bed_intervals(d["peaks"], chrom)
    flag_fragments(fr, blacklist, trf, cen, peaks)
    clean = (fr.flags & F_BLACKLIST) == 0
    n_cells = len(cells)

    params = {"version": __version__, "sample": d["sample"], "bam": d["bam"], "barcodes": d["barcodes"],
              "barcode_space": space, "barcode_map": d["barcode_map"], "chrom": chrom,
              "region": a.region, "binsize": a.binsize, "assign": a.assign, "min_mapq": a.min_mapq, "min_bq": a.min_bq,
              "snps": d["snps"], "phased": d["phased"], "peaks": d["peaks"],
              "layers": sorted(layers_on), "ref": ref, "command": " ".join(sys.argv)}

    tiling = Tiling.regular_grid(chrom, length, a.binsize)
    with open(pdir / "cells.txt", "w") as fh:
        fh.write("\n".join(cells) + "\n")

    if "depth" in layers_on:
        bins_df = annotate_tiling(tiling, ref, d["peaks"])
        layers = depth_layers(fr, tiling, n_cells, assign=a.assign, use=clean)
        for k in ("hapA", "hapB", "hap_ambig"):
            layers.pop(k, None)
        obsm = {"fraglen_hist": fraglen_hist(fr, n_cells, use=clean), "mapq_hist": ctr.mapq_hist}
        _log(f"depth: {int(layers['counts'].sum())} counts in {len(tiling)} bins")
        save_csr_dict(pdir / "layers.npz", layers)
        np.savez_compressed(pdir / "obsm.npz", **obsm)
        bins_df.to_csv(pdir / "bins.tsv", sep="\t")
        qc = chrom_qc(fr, ctr, cells, hap_counts=False)
        qc.to_csv(pdir / "qc.tsv", sep="\t")

    if "alleles" in layers_on and acc is not None:
        al = acc.matrices(n_cells)
        hl = depth_layers(fr, tiling, n_cells, assign=a.assign if a.assign != "overlap" else "start", use=clean)
        for k in ("hapA", "hapB", "hap_ambig"):
            al[k] = hl[k]
        save_csr_dict(pdir / "alleles.npz", al)
        sdf = acc.snps.frame(); sdf["n_other"] = acc.n_other; sdf["n_lowbq"] = acc.n_lowbq
        sdf.to_csv(pdir / "snps.tsv", sep="\t")
        aq = chrom_qc(fr, ctr, cells, hap_counts=True)[["n_hapA", "n_hapB", "n_hap_ambig", "n_snp_obs"]]
        aq.to_csv(pdir / "alleles_qc.tsv", sep="\t")
        _log(f"alleles: DP {int(al['DP'].sum())}, AD {int(al['AD'].sum())}, hap-assigned fragments "
             f"{int((fr.hap == 1).sum() + (fr.hap == 2).sum())}, ambiguous {int((fr.hap == 3).sum())}")

    if "fragments" in layers_on:
        write_bgzf_table(pdir / "fragments.tsv.gz", fragments_frame(fr, cells))

    params["seconds"] = round(time.time() - t0, 1)
    params["n_fragments"] = int(len(fr)); params["n_clean"] = int(clean.sum())
    params["dash"] = ctr.dash; params["unlisted"] = ctr.unlisted; params["n_no_cb"] = ctr.n_no_cb
    write_json(pdir / "params.json", params)
    (pdir / "DONE").touch()
    _log(f"done {chrom} in {params['seconds']} s -> {pdir}")
    return 0


# ----------------------------------------------------------------------------------------------- merge
def cmd_merge(a):
    from .merge import merge_sample
    chroms = a.chroms.split(",") if a.chroms else None
    out = a.out
    if a.config and a.sample:
        cfg = load_yaml(a.config); sc = sample_from_config(cfg, a.sample)
        out = out or sc.out
        chroms = chroms or sc.chroms
    chroms = chroms or MAIN_CHROMS
    s = merge_sample(out, chroms, tag=a.tag)
    _log(json.dumps({k: v for k, v in s.items() if k != "params"}, default=str))
    return 0


# ----------------------------------------------------------------------------------------------- run
def _run_one(ns):
    """One chromosome of `run` (in a worker process); returns (chrom, exit code, message)."""
    try:
        return ns.chrom, cmd_extract(ns), ""
    except SystemExit as e:
        return ns.chrom, 1, str(e)
    except Exception as e:  # noqa: BLE001 - reported per chromosome, the run fails at the end
        return ns.chrom, 1, f"{type(e).__name__}: {e}"


def cmd_run(a):
    """extract on every chromosome (a process pool of --jobs workers), then merge: the portable alternative to one
    cluster job per chromosome."""
    import copy
    from concurrent.futures import ProcessPoolExecutor
    import pysam
    _layers(a.layers)
    d = resolve_inputs(a)
    if a.chroms:
        chroms = a.chroms.split(",")
    elif d["chroms"]:
        chroms = d["chroms"]
    else:
        with pysam.AlignmentFile(d["bam"], "rb", index_filename=a.index) as bf:
            present = set(bf.references)
        chroms = [c for c in MAIN_CHROMS if c in present]
    if not chroms:
        raise SystemExit("run: no chromosomes (none of chr1-22, X, Y, M in the BAM header; give --chroms)")
    jobs = max(1, min(a.jobs, len(chroms)))
    _log(f"run: {len(chroms)} chromosomes, {jobs} parallel jobs -> {d['out']}")
    nss = []
    for c in chroms:
        ns = copy.copy(a); ns.chrom = c; ns.region = None; ns.out = d["out"]
        nss.append(ns)
    fails = []
    if jobs == 1:
        res = [_run_one(ns) for ns in nss]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            res = list(ex.map(_run_one, nss))
    for c, code, msg in res:
        if code != 0:
            fails.append(c); _log(f"run: {c} failed {msg}")
    if fails:
        raise SystemExit(f"run: extraction failed on {', '.join(fails)}; partials in {d['out']}/partials")
    from .merge import merge_sample
    s = merge_sample(d["out"], chroms, tag=a.tag)
    _log(json.dumps({k: v for k, v in s.items() if k != "params"}, default=str))
    return 0


# ----------------------------------------------------------------------------------------------- select-hets
def cmd_select_hets(a):
    from .hets import select_hets
    split = lambda s: [c for c in s.split(",") if c]
    select_hets(a.vcf, a.out, loci=a.loci, contigs=split(a.contigs) if a.contigs else None,
                no_loci_contigs=split(a.no_loci_contigs), min_gq=a.min_gq, min_dp=a.min_dp, dp_quantile=a.dp_quantile,
                stats=a.stats, log=_log)
    return 0


# ----------------------------------------------------------------------------------------------- hybrid-phase
def cmd_hybrid_phase(a):
    from .hybrid import hybrid_phase
    hybrid_phase(a.blocks, a.population, a.out, contigs=[c for c in a.contigs.split(",") if c] if a.contigs else None,
                 sample=a.sample, population_sample=a.population_sample, imputed_flag=a.imputed_flag or None,
                 min_votes=a.min_votes, fill_unphased=a.fill_unphased, block_tsv=a.block_tsv, chrom_tsv=a.chrom_tsv,
                 log=_log)
    return 0


# ----------------------------------------------------------------------------------------------- parser
def build_parser():
    p = argparse.ArgumentParser(prog="scphasedcn", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("index-check", help="check the BAM index (and optionally rebuild it at a separate path, e.g. in the output directory)")
    s.add_argument("--bam", required=True); s.add_argument("--index", help="existing .bai to check")
    s.add_argument("--reindex", help="write a fresh .bai to this path when missing/stale (use with --index in extract / run)")
    s.add_argument("--force", action="store_true"); s.add_argument("--threads", type=int, default=8)
    s.set_defaults(func=cmd_index_check)

    s = sub.add_parser("extract", help="one chromosome: depth, alleles, fragments, qc partials (then `merge`)")
    s.add_argument("--chrom"); s.add_argument("--region", help="chrom:start-end (0-based half-open; fragments by start)")
    add_extract_args(s)
    s.set_defaults(func=cmd_extract)

    s = sub.add_parser("run", help="extract every chromosome in parallel, then merge (one command per sample)")
    s.add_argument("--chroms", help="comma list (default: the yaml sample's chroms, else chr1-22, X, Y, M in the BAM)")
    s.add_argument("--jobs", type=int, default=4, help="chromosomes extracted in parallel (memory: ~2-9 GB each)")
    s.add_argument("--tag", default="100kb", help="merge tag (bins_<tag>.h5ad; the caller reads bins_100kb.h5ad)")
    add_extract_args(s)
    s.set_defaults(func=cmd_run, chrom=None, region=None)

    s = sub.add_parser("select-hets", help="heterozygous SNVs of a phased VCF for the allele layer (see scphasedcn.hets)")
    s.add_argument("--vcf", required=True, help="phased VCF (bgzip + tabix), e.g. Clair3 + longphase")
    s.add_argument("--out", required=True, help="output .vcf.gz (a .tbi is written next to it)")
    s.add_argument("--loci", help="population sites VCF (bgzip + tabix; recommended, e.g. 1000 Genomes SNVs)")
    s.add_argument("--contigs", help="comma list (default chr1-22, chrX)")
    s.add_argument("--no-loci-contigs", default="chrX", help="contigs kept without the loci intersect ('' for none)")
    s.add_argument("--min-gq", type=int, default=10); s.add_argument("--min-dp", type=int, default=20)
    s.add_argument("--dp-quantile", type=float, default=0.99); s.add_argument("--stats", help="per-contig counts TSV")
    s.set_defaults(func=cmd_select_hets)

    s = sub.add_parser("hybrid-phase", help="orient read-backed phase blocks with a population-phased VCF: one phase set "
                                            "per contig (normal-tissue mode; see scphasedcn.hybrid)")
    s.add_argument("--blocks", required=True, help="read-backed phased VCF (GT | + PS), e.g. longphase or the select-hets "
                                                   "output")
    s.add_argument("--population", required=True, help="population-phased VCF of the same sample (e.g. Beagle 5, "
                                                       "SHAPEIT), or a pattern with {chrom} for one file per contig")
    s.add_argument("--out", required=True, help="output .vcf.gz (a .tbi is written next to it)")
    s.add_argument("--contigs", help="comma list of contigs to orient (default all; others are written unchanged)")
    s.add_argument("--sample", help="sample of --blocks (default: the only one)")
    s.add_argument("--population-sample", help="sample of --population (default: the --blocks sample, else the only one)")
    s.add_argument("--imputed-flag", default="IMP", help="INFO flag of imputed population markers, which never vote "
                                                         "('' for none)")
    s.add_argument("--min-votes", type=int, default=1, help="minimum voting sites to orient a block")
    s.add_argument("--fill-unphased", action="store_true", help="hets outside the blocks take the population phase")
    s.add_argument("--block-tsv", help="per-block report TSV"); s.add_argument("--chrom-tsv", help="per-contig report TSV")
    s.set_defaults(func=cmd_hybrid_phase)

    s = sub.add_parser("merge", help="merge chromosome partials into h5ad / tsv outputs")
    s.add_argument("--config"); s.add_argument("--sample"); s.add_argument("--out"); s.add_argument("--chroms")
    s.add_argument("--tag", default="100kb")
    s.set_defaults(func=cmd_merge)
    return p


def add_extract_args(s):
    """Options shared by `extract` (one chromosome) and `run` (all chromosomes)."""
    s.add_argument("--config", help="sample yaml (layout: inst/extdata/sample_template.yaml of the scphasedcn R package)")
    s.add_argument("--sample", help="sample key of the yaml")
    s.add_argument("--bam"); s.add_argument("--index", help=".bai to use instead of the one next to the BAM")
    s.add_argument("--barcodes", help="cell whitelist (one barcode per line)")
    s.add_argument("--barcode-map", help="TSV with columns atac_barcode, rna_barcode: a whitelist in RNA-barcode space "
                                         "(e.g. 10x Multiome GEX barcodes) is translated to the ATAC barcodes of the "
                                         "BAM's CB tag")
    s.add_argument("--out", help="sample output directory (partials/<chrom>/ is created)")
    s.add_argument("--binsize", type=int, default=100_000)
    s.add_argument("--assign", default="start", choices=["start", "midpoint", "end5", "overlap"],
                   help="fragment position counted in a bin: start, midpoint, 5' end, or every bin it overlaps "
                        "(default start)")
    s.add_argument("--snps", help="VCF(.gz, indexed) of het SNPs"); s.add_argument("--phased", action="store_true")
    s.add_argument("--no-snps", action="store_true", help="ignore the yaml VCF")
    s.add_argument("--peaks", help="peak BED file for FRiP / peak_bp (optional)")
    s.add_argument("--min-mapq", type=int, default=30); s.add_argument("--min-bq", type=int, default=10)
    s.add_argument("--layers", default=DEFAULT_LAYERS,
                   help=f"comma list of {','.join(LAYER_CHOICES)} (default {DEFAULT_LAYERS}: what the R caller reads)")
    s.add_argument("--fasta", help="reference FASTA (hg38, chr names; required: GC per bin)")
    s.add_argument("--blacklist", help="blacklist BED (default: the bundled ENCODE hg38 blacklist v2)")
    s.add_argument("--trf", help="tandem-repeat BED (optional)")
    s.add_argument("--centromeres", help="centromere table (default: the bundled hg38 table)")
    s.add_argument("--chrom-sizes", help="chromosome sizes (default: the bundled hg38 sizes)")
    s.add_argument("--no-gc", action="store_true", help="no FASTA: no GC per bin (the output is not callable)")
    s.add_argument("--no-trf", action="store_true", help="no tandem-repeat mask")
    s.add_argument("--no-masks", action="store_true", help="no blacklist, TRF or centromere masks")


def main(argv=None):
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
