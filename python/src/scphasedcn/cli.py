"""scphasedcn command line: run | extract | merge | select-hets | index-check.

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

LAYER_CHOICES = ("depth", "alleles", "fragments")
DEFAULT_LAYERS = "depth,alleles,fragments"


def _layers(s: str) -> set[str]:
    """--layers value -> set of layer names (unknown names are an error). `qc` (a layer name of older versions) is
    accepted and ignored: the per-cell QC tables are always written with depth (qc.tsv) and alleles (alleles_qc.tsv)."""
    on = {x for x in s.split(",") if x}
    if "qc" in on:
        on.discard("qc")
        _log("note: --layers qc is not a layer of its own (ignored): the per-cell QC tables are written with depth "
             "(qc.tsv) and alleles (alleles_qc.tsv)")
    bad = sorted(on - set(LAYER_CHOICES))
    if bad:
        raise SystemExit(f"--layers: unknown layer(s) {', '.join(bad)}; choose from {','.join(LAYER_CHOICES)}")
    if not on:
        raise SystemExit(f"--layers: nothing to extract; choose from {','.join(LAYER_CHOICES)}")
    return on


def _layers_written(on: set[str]) -> list[str]:
    """Layer names recorded in params.json: `qc` is listed when QC tables are written (with depth or alleles)."""
    return sorted(on | ({"qc"} if on & {"depth", "alleles"} else set()))


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
         "peaks": a.peaks, "out": a.out, "chroms": None, "barcode_map": getattr(a, "barcode_map", None),
         "vcf_sample": getattr(a, "vcf_sample", None)}
    if a.config and a.sample:
        sc = sample_from_config(cfg, a.sample)
        for k in ("bam", "barcodes", "snps", "peaks", "out", "barcode_map"):
            if d[k] is None:
                d[k] = getattr(sc, k)
        if not a.phased and sc.snps and d["snps"] == sc.snps:
            d["phased"] = sc.phased
        if d["vcf_sample"] is None and sc.snps and d["snps"] == sc.snps:
            d["vcf_sample"] = sc.vcf_sample
        d["chroms"] = sc.chroms
    if getattr(a, "no_snps", False):
        d["snps"] = None; d["phased"] = False; d["vcf_sample"] = None
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


def _bam_references(bam, index):
    import pysam
    with pysam.AlignmentFile(bam, "rb", index_filename=index) as bf:
        return dict(zip(bf.references, bf.lengths))


def cmd_extract(a):
    from .bam import scan_chrom, flag_fragments, check_barcode_match, F_BLACKLIST
    from .bins import Tiling, annotate_tiling, read_bed_intervals, read_centromere
    from .depth import depth_layers, fraglen_hist
    from .qc import chrom_qc
    from .writers import (save_csr_dict, write_bgzf_table, fragments_frame, write_json, write_done, read_done,
                          file_fingerprint, cells_digest, LAYER_FILES, DONE_FORMAT)
    from .alleles import check_vcf, load_snps, AlleleAccumulator, CONTIG_PS

    t0 = time.time()
    layers_on = _layers(a.layers)
    d = resolve_inputs(a)
    reg = _region(a.region)
    chrom = reg[0] if reg else a.chrom
    if chrom is None:
        raise SystemExit("need --chrom or --region")
    cells, space = translate_barcodes(read_barcodes(d["barcodes"]), d["barcode_map"], log=_log)
    if space != "atac":
        _log(f"whitelist translated {space} via {d['barcode_map']}: {len(cells)} cells")
    sizes = read_chrom_sizes(d["ref"]["chrom_sizes"])
    use_snps = "alleles" in layers_on and bool(d["snps"])
    bam_refs = _bam_references(d["bam"], a.index) if (chrom not in sizes or use_snps) else None
    if chrom not in sizes:
        if chrom not in bam_refs:
            raise SystemExit(f"{chrom} is neither in the chromosome sizes nor a reference of the BAM")
        sizes[chrom] = bam_refs[chrom]
    length = sizes[chrom]
    vcf_sample = None
    if use_snps:   # before the BAM pass: an unusable VCF fails here, not after hours with an empty allele layer
        vcf_sample, warns = check_vcf(d["snps"], [chrom], sample=d["vcf_sample"], reference_names=list(bam_refs))
        for w in warns:
            _log("warning:", w)
    # a --region extraction has its own partial, which `merge` refuses: it never replaces the chromosome's partial
    pdir = Path(d["out"]) / "partials" / (f"region_{chrom}_{reg[1][0]}-{reg[1][1]}" if reg else chrom)
    pdir.mkdir(parents=True, exist_ok=True)
    old_meta = read_done(pdir)
    (pdir / "DONE").unlink(missing_ok=True)
    produce = layers_on & {"depth", "fragments"} | ({"alleles"} if use_snps else set())
    for layer in sorted(layers_on - produce):   # requested but not produced (alleles without --snps): drop stale files
        for f in LAYER_FILES[layer]:
            (pdir / f).unlink(missing_ok=True)
    _log(f"scphasedcn {__version__} extract {d['sample']} {chrom} ({length} bp) cells={len(cells)} layers={sorted(layers_on)}")

    # SNPs
    acc = None
    if use_snps:
        snps = load_snps(d["snps"], chrom, phased=d["phased"], region=reg[1] if reg else None, sample=vcf_sample)
        acc = AlleleAccumulator(snps, min_bq=a.min_bq)
        _log(f"SNPs on {chrom}: {len(snps)} ({int((snps.hap_of_alt > 0).sum())} phased) from {d['snps']}")
        if snps.n_phased_gt:
            _log(f"warning: {snps.n_phased_gt} SNPs on {chrom} have a phased genotype but --phased was not given: "
                 f"their phase is ignored (no haplotype votes)")
        if snps.n_contig_ps:
            _log(f"warning: {snps.n_contig_ps} phased SNPs on {chrom} have no PS: one phase set per contig is "
                 f"assumed (PS {CONTIG_PS}, VCF specification), i.e. the whole chromosome is taken as correctly "
                 f"phased. Population-only (statistical) phasing of one sample carries switch errors every few "
                 f"megabases, which attenuate the allele signal over large regions; read-based phase sets (PS from "
                 f"e.g. longphase or WhatsHap) are preferable where available")

    # BAM pass
    fr, ctr = scan_chrom(d["bam"], chrom, cells, min_mapq=a.min_mapq, index_filename=a.index,
                         region=reg[1] if reg else None, alleles=acc)
    _log(f"scan: {int(ctr.per_cell[:, 0].sum())} records of listed cells, {len(fr)} passing fragments, "
         f"'-' records {ctr.dash['n_records']} (pass {ctr.dash['n_pass']}), unlisted {ctr.unlisted['n_records']}, "
         f"{ctr.seconds:.0f} s")
    bc_warn = check_barcode_match(ctr, cells, chrom)
    if bc_warn:
        _log("WARNING:", bc_warn)

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
              "layers": _layers_written(layers_on), "ref": ref, "command": " ".join(sys.argv)}

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

    hap_assign = a.assign if a.assign != "overlap" else "start"
    if "alleles" in layers_on and acc is not None:
        al = acc.matrices(n_cells)
        hl = depth_layers(fr, tiling, n_cells, assign=hap_assign, use=clean)
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

    # DONE carries what the partial holds and the inputs of each layer; `merge` checks they agree across partials
    cdig = cells_digest(cells)
    common = {"cells": cdig, "bam": file_fingerprint(d["bam"]), "min_mapq": a.min_mapq}
    masks = {k: file_fingerprint(ref.get(k)) for k in ("blacklist", "trf", "centromeres")}
    masks["peaks"] = file_fingerprint(d["peaks"])
    snp_in = {"snps": file_fingerprint(d["snps"]) if acc is not None else None, "phased": bool(d["phased"]) and acc is not None,
              "vcf_sample": vcf_sample, "min_bq": a.min_bq if acc is not None else None}
    inputs = {"depth": {**common, "binsize": a.binsize, "assign": a.assign, "fasta": file_fingerprint(ref.get("fasta")),
                        **masks},
              "alleles": {**common, "binsize": a.binsize, "assign": hap_assign, **snp_in},
              "fragments": {**common, **masks, **snp_in}}
    lmeta = {}
    for layer, files in LAYER_FILES.items():
        present = [f for f in files if (pdir / f).exists()]
        if layer in produce:
            lmeta[layer] = {"files": present, "inputs": inputs[layer]}
        elif present:   # not requested in this run: kept from an earlier run (inputs unknown if that run was older)
            lmeta[layer] = (old_meta["layers"][layer] if old_meta and layer in old_meta["layers"]
                            else {"files": present, "inputs": None})
    write_done(pdir, {"format": DONE_FORMAT, "version": __version__, "chrom": chrom, "region": a.region,
                      "length": int(length), "binsize": a.binsize, "cells": cdig, "n_cells": n_cells, "layers": lmeta})
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
    s = merge_sample(out, chroms, tag=a.tag, allow_missing=a.allow_missing, log=_log)
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
    from .alleles import check_vcf
    layers_on = _layers(a.layers)
    d = resolve_inputs(a)
    present = _bam_references(d["bam"], a.index)
    if a.chroms or d["chroms"]:   # requested explicitly: every one must be a reference of the BAM, checked up front
        chroms = a.chroms.split(",") if a.chroms else list(d["chroms"])
        absent = [c for c in chroms if c not in present]
        if absent:
            raise SystemExit(f"run: {', '.join(absent)} not in the BAM header (requested with "
                             f"{'--chroms' if a.chroms else 'the yaml chroms'}); leave them out")
    else:
        chroms = [c for c in MAIN_CHROMS if c in present]
    if not chroms:
        raise SystemExit("run: no chromosomes (none of chr1-22, X, Y, M in the BAM header; give --chroms)")
    if "alleles" in layers_on and d["snps"]:   # once, before any BAM pass (per-chromosome warnings come from extract)
        check_vcf(d["snps"], [], sample=d["vcf_sample"], reference_names=list(present))
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
    s.add_argument("--chrom"); s.add_argument("--region", help="chrom:start-end (0-based half-open; records by start); "
                                                              "written to partials/region_<chrom>_<start>-<end>, "
                                                              "which `merge` does not merge (testing / inspection)")
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

    s = sub.add_parser("merge", help="merge chromosome partials into h5ad / tsv outputs")
    s.add_argument("--config"); s.add_argument("--sample"); s.add_argument("--out"); s.add_argument("--chroms")
    s.add_argument("--tag", default="100kb")
    s.add_argument("--allow-missing", action="store_true",
                   help="merge although some chromosomes have no finished partial (default: fail and list them)")
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
    s.add_argument("--vcf-sample", help="sample column of --snps to use (required when the VCF has several samples; "
                                        "yaml: vcf_sample)")
    s.add_argument("--no-snps", action="store_true", help="ignore the yaml VCF")
    s.add_argument("--peaks", help="peak BED file for FRiP / peak_bp (optional)")
    s.add_argument("--min-mapq", type=int, default=30); s.add_argument("--min-bq", type=int, default=10)
    s.add_argument("--layers", default=DEFAULT_LAYERS,
                   help=f"comma list of {','.join(LAYER_CHOICES)} (default {DEFAULT_LAYERS}: what the R caller reads; "
                        f"the per-cell QC tables come with depth and alleles)")
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
