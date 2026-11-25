"""Select heterozygous SNVs from a phased VCF (e.g. Clair3 + longphase) for the allele layer of `extract`.

Rule: FILTER PASS, biallelic SNV, heterozygous GT, GQ >= min_gq, DP in [min_dp, q(dp_quantile)]. The DP upper bound is
the quantile of DP over the candidate sites of the processed contigs (first pass); it drops collapsed repeats. There
is no allele-frequency window (it would delete LOH and amplified sites in tumours). All FORMAT fields are kept.

Optionally (recommended) the sites are intersected with a population sites VCF (`loci`, e.g. 1000 Genomes SNVs:
position + REF + ALT), which removes most systematic long-read genotyping errors. Contigs in `no_loci_contigs` (default
chrX: population panels are often autosomes only) are kept without the intersect. Loci contig names may differ from the
VCF's by the `chr` prefix (`1` vs `chr1`). Selecting no site at all is an error and leaves no output (nor a `stats`
file of an earlier run); selecting sites but none on an autosome while autosomes had records (or were named in
`contigs`) is a loud warning (e.g. a loci file that matches chrX only).
"""
from __future__ import annotations

import sys
from collections import Counter

import numpy as np

from .hybrid import _bare, _names, is_snv

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def is_candidate(rec, min_gq: int, min_dp: int) -> tuple[bool, str]:
    """PASS, biallelic SNV, het GT, GQ >= min_gq, DP >= min_dp. Returns (ok, reason)."""
    if "PASS" not in rec.filter:
        return False, "not_pass"
    if not is_snv(rec):                 # also REF / ALT `*`, N, IUPAC codes
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
    import os
    import pysam
    log = log or (lambda *a: print(*a, file=sys.stderr, flush=True))
    if stats:   # never a stats file of an earlier run next to this run's error
        try:
            os.unlink(stats)
        except FileNotFoundError:
            pass

    def fail(msg):   # no output that looks like a selection (nor an index of an earlier one)
        for f in (out, f"{out}.tbi"):
            if os.path.exists(f):
                os.unlink(f)
        raise SystemExit(msg)
    given = list(contigs) if contigs else []
    contigs = given or AUTOSOMES + ["chrX"]
    no_loci = set(["chrX"] if no_loci_contigs is None else no_loci_contigs)
    vcf = pysam.VariantFile(vcf_path)
    lf = pysam.VariantFile(loci) if loci else None
    loci_contigs = set(lf.header.contigs) if lf is not None else set()
    loci_name = lambda c: next((x for x in _names(c) if x in loci_contigs), None)   # `chr1` here may be `1` there
    missing = [c for c in contigs if c not in set(vcf.header.contigs)]
    if missing:
        log(f"[warn] contigs absent from the VCF header, skipped: {missing}")
        contigs = [c for c in contigs if c not in missing]

    # pass 1: DP distribution of the candidates -> upper bound
    dps, reasons, nrec = [], Counter(), Counter()
    for c in contigs:
        for rec in vcf.fetch(c):
            nrec[c] += 1
            ok, why = is_candidate(rec, min_gq, min_dp)
            reasons[why] += 1
            if ok:
                dps.append(rec.samples[0]["DP"])
    if not dps:
        fail("select-hets: no candidate sites")
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
        lc = loci_name(c) if lf is not None else None
        use_loci = lc is not None and c not in no_loci
        if lf is not None and lc is None and c not in no_loci:
            log(f"[warn] {c}: not in the loci file and not in no_loci_contigs -> all sites dropped")
        elif use_loci and lc != c:
            log(f"[{c}] loci contig {lc}")
        lit = lf.fetch(lc) if use_loci else None
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
    log("[total] " + " ".join(f"{k}={v}" for k, v in sorted(tot.items())))
    is_auto = lambda c: _bare(c) in {str(i) for i in range(1, 23)}
    auto = [c for c in contigs if is_auto(c)]
    hint = ""
    if lf is not None:
        lcs = sorted(loci_contigs)
        hint = (f"; the loci file has contigs {', '.join(lcs[:4]) or 'none'}{', ...' if len(lcs) > 4 else ''}: do "
                f"they match the VCF's ({', '.join(contigs[:4])}{', ...' if len(contigs) > 4 else ''})?")
    if tot["selected"] == 0:
        fail(f"select-hets: {'no autosomal site' if auto else 'no site'} selected ("
             + " ".join(f"{k}={v}" for k, v in sorted(tot.items())) + f"){hint}")
    with_rec = [c for c in auto if nrec[c]]
    named = [c for c in given if is_auto(c)]
    if sum(n["selected"] for c, n in rows if c in auto) == 0 and (with_rec or named):
        log(f"[WARN] select-hets: NO AUTOSOMAL SITE selected ({tot['selected']} sites on "
            f"{', '.join(c for c, n in rows if n['selected'])}), although "
            + (f"{', '.join(with_rec[:4])}{', ...' if len(with_rec) > 4 else ''} had records" if with_rec
               else f"--contigs names {', '.join(named[:4])}{', ...' if len(named) > 4 else ''}")
            + f": the allele layer will have no autosomal site{hint}")
    pysam.tabix_index(out, preset="vcf", force=True)
    if stats:
        keys = ["candidate", "dp_above_q", "not_in_loci", "pos_in_loci_allele_mismatch", "selected", "selected_phased"]
        with open(stats, "w") as fh:
            fh.write("contig\t" + "\t".join(keys) + "\n")
            for c, n in rows:
                fh.write(c + "\t" + "\t".join(str(n[k]) for k in keys) + "\n")
            fh.write("total\t" + "\t".join(str(tot[k]) for k in keys) + "\n")
            fh.write(f"# dp_max={dp_max} min_gq={min_gq} min_dp={min_dp} n_candidates_pass1={len(dps)}\n")
    return dict(tot) | {"dp_max": dp_max}
