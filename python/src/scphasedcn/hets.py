"""Select heterozygous SNVs from a phased VCF (e.g. Clair3 + longphase) for the allele layer of `extract`.

Rule: FILTER PASS, biallelic SNV, heterozygous GT, GQ >= min_gq, DP in [min_dp, q(dp_quantile)]. The DP upper bound is
the quantile of DP over the candidate sites of the processed contigs (first pass); it drops collapsed repeats. There
is no allele-frequency window (it would delete LOH and amplified sites in tumours). All FORMAT fields are kept.

Optionally (recommended) the sites are intersected with a population sites VCF (`loci`, e.g. 1000 Genomes SNVs:
position + REF + ALT), which removes most systematic long-read genotyping errors. Contigs in `no_loci_contigs` (default
chrX: population panels are often autosomes only) are kept without the intersect.
"""
from __future__ import annotations

import sys
from collections import Counter

import numpy as np

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def is_candidate(rec, min_gq: int, min_dp: int) -> tuple[bool, str]:
    """PASS, biallelic SNV, het GT, GQ >= min_gq, DP >= min_dp. Returns (ok, reason)."""
    if "PASS" not in rec.filter:
        return False, "not_pass"
    if rec.alts is None or len(rec.alts) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
        return False, "not_biallelic_snv"
    s = rec.samples[0]
    gt = s.get("GT")
    if gt is None or len(gt) != 2 or None in gt or gt[0] == gt[1]:
        return False, "not_het"
    gq = s.get("GQ")
    dp = s.get("DP")
    if gq is None or gq < min_gq:
        return False, "low_gq"
    if dp is None or dp < min_dp:
        return False, "low_dp"
    return True, "ok"


def select_hets(vcf_path: str, out: str, loci: str | None = None, contigs: list[str] | None = None,
                no_loci_contigs: list[str] | None = None, min_gq: int = 10, min_dp: int = 20,
                dp_quantile: float = 0.99, stats: str | None = None, log=None) -> dict:
    """Write the selected sites to `out` (.vcf.gz + .tbi). Returns the total counts (candidate, selected, ...)."""
    import pysam
    log = log or (lambda *a: print(*a, file=sys.stderr, flush=True))
    contigs = list(contigs) if contigs else AUTOSOMES + ["chrX"]
    no_loci = set(["chrX"] if no_loci_contigs is None else no_loci_contigs)
    vcf = pysam.VariantFile(vcf_path)
    lf = pysam.VariantFile(loci) if loci else None
    loci_contigs = set(lf.header.contigs) if lf is not None else set()
    missing = [c for c in contigs if c not in set(vcf.header.contigs)]
    if missing:
        log(f"[warn] contigs absent from the VCF header, skipped: {missing}")
        contigs = [c for c in contigs if c not in missing]

    # pass 1: DP distribution of the candidates -> upper bound
    dps, reasons = [], Counter()
    for c in contigs:
        for rec in vcf.fetch(c):
            ok, why = is_candidate(rec, min_gq, min_dp)
            reasons[why] += 1
            if ok:
                dps.append(rec.samples[0]["DP"])
    if not dps:
        raise SystemExit("select-hets: no candidate sites")
    dps = np.asarray(dps)
    dp_max = int(np.quantile(dps, dp_quantile))
    log(f"[pass1] candidates={len(dps)} median_DP={int(np.median(dps))} DP_q{dp_quantile:.2f}={dp_max} reasons={dict(reasons)}")

    # pass 2: DP upper bound and the optional population-sites intersect (merge join of sorted streams)
    hdr = vcf.header
    rule = f"PASS, het, biallelic SNV, GQ>={min_gq}, DP in [{min_dp},{dp_max}] (q{dp_quantile:.2f} of candidates)"
    rule += f", loci intersect (pos+REF+ALT) except {sorted(no_loci)}; loci={loci}" if loci else ", no loci intersect"
    hdr.add_line(f'##select_hets="{rule}; no AF window"')
    ofh = pysam.VariantFile(out, "wz", header=hdr)
    rows, tot = [], Counter()
    for c in contigs:
        n = Counter()
        use_loci = lf is not None and c in loci_contigs and c not in no_loci
        if lf is not None and c not in loci_contigs and c not in no_loci:
            log(f"[warn] {c}: not in the loci file and not in no_loci_contigs -> all sites dropped")
        lit = lf.fetch(c) if use_loci else None
        lrec = next(lit, None) if lit is not None else None
        lbuf_pos, lbuf = None, []          # loci at the current position (kept for a second record at that position)
        for rec in vcf.fetch(c):
            ok, _ = is_candidate(rec, min_gq, min_dp)
            if not ok:
                continue
            n["candidate"] += 1
            if rec.samples[0]["DP"] > dp_max:
                n["dp_above_q"] += 1
                continue
            if use_loci:
                if rec.pos != lbuf_pos:
                    while lrec is not None and lrec.pos < rec.pos:
                        lrec = next(lit, None)
                    lbuf = []
                    while lrec is not None and lrec.pos == rec.pos:   # several loci can share a position
                        lbuf.append((lrec.ref, lrec.alts))
                        lrec = next(lit, None)
                    lbuf_pos = rec.pos
                pos_hit = bool(lbuf)
                hit = any(r == rec.ref and alts and rec.alts[0] in alts for r, alts in lbuf)
                if not hit:
                    n["pos_in_loci_allele_mismatch" if pos_hit else "not_in_loci"] += 1
                    continue
            elif lf is not None and c not in no_loci:
                continue
            n["selected"] += 1
            if rec.samples[0].phased:
                n["selected_phased"] += 1
            ofh.write(rec)
        tot.update(n)
        rows.append((c, n))
        log(f"[{c}] " + " ".join(f"{k}={v}" for k, v in sorted(n.items())))
    ofh.close()
    pysam.tabix_index(out, preset="vcf", force=True)
    log("[total] " + " ".join(f"{k}={v}" for k, v in sorted(tot.items())))
    if stats:
        keys = ["candidate", "dp_above_q", "not_in_loci", "pos_in_loci_allele_mismatch", "selected", "selected_phased"]
        with open(stats, "w") as fh:
            fh.write("contig\t" + "\t".join(keys) + "\n")
            for c, n in rows:
                fh.write(c + "\t" + "\t".join(str(n[k]) for k in keys) + "\n")
            fh.write("total\t" + "\t".join(str(tot[k]) for k in keys) + "\n")
            fh.write(f"# dp_max={dp_max} min_gq={min_gq} min_dp={min_dp} n_candidates_pass1={len(dps)}\n")
    return dict(tot) | {"dp_max": dp_max}
