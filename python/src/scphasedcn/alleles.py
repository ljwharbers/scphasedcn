"""Allele layer: per-read ref/alt observations at heterozygous SNPs and per-read haplotype votes.

SNPs come from a bgzipped, indexed VCF (bi-allelic SNVs; with a sample column only GT-het records are used, genotypes
with a missing allele are skipped; `phased=True` reads the `|` phase and the PS field: hapA = first allele of the GT;
phased genotypes without a PS value form one phase set per contig, as in the VCF specification, and get PS
`CONTIG_PS`). A VCF without sample columns (e.g. a het-site list) yields unphased sites; a VCF with several sample
columns needs the sample named (`check_vcf`, `--vcf-sample`). Per read, the CIGAR is converted to cumulative reference/query offsets with numpy and only the
SNPs inside the aligned span are looked up (I/D/S/H/N handled; deletions at a SNP count as `other`). Bases with
quality < min_bq are ignored. The per-read haplotype is a majority vote over phased SNPs (ambiguous on ties).
"""
from __future__ import annotations

from array import array
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

REF_CONSUME = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1, 0], dtype=np.int64)   # M I D N S H P = X B
Q_CONSUME = np.array([1, 1, 0, 0, 1, 0, 0, 1, 1, 0], dtype=np.int64)
MATCH_OPS = np.array([1, 0, 0, 0, 0, 0, 0, 1, 1, 0], dtype=bool)
CONTIG_PS = 1   # phase set of phased genotypes without a PS value (VCF spec: one phase set per contig); must be > 0


@dataclass
class SnpSet:
    chrom: str
    pos: np.ndarray          # int64, 0-based, sorted
    ref: list
    alt: list
    hap_of_alt: np.ndarray   # int8: 0 unphased, 1 alt on hapA, 2 alt on hapB
    ps: np.ndarray           # int64 phase set id, -1 when unphased
    ids: list
    n_phased_gt: int = 0     # phased het genotypes seen while phased=False (their phase is ignored)
    n_contig_ps: int = 0     # phased het genotypes without a PS value (given CONTIG_PS)

    def __len__(self):
        return len(self.pos)

    def frame(self):
        import pandas as pd
        df = pd.DataFrame({"chrom": self.chrom, "pos": self.pos + 1, "ref": self.ref, "alt": self.alt,
                           "hap_of_alt": self.hap_of_alt, "ps": self.ps})
        df.index = [f"{self.chrom}:{p}{r}>{a}" for p, r, a in zip(df["pos"], df["ref"], df["alt"])]
        return df


def _no_index(vcf_path) -> str:
    return (f"--snps {vcf_path}: no index (.tbi / .csi) found; the VCF must be bgzip-compressed and indexed, e.g. "
            f"`bgzip -c snps.vcf > snps.vcf.gz && tabix -p vcf snps.vcf.gz`")


def _alias(c: str) -> str:
    """The other contig naming convention of `c` (chr2 <-> 2, chrM <-> MT)."""
    if c in ("chrM", "MT"):
        return "MT" if c == "chrM" else "chrM"
    return c[3:] if c.startswith("chr") else "chr" + c


def check_vcf(vcf_path: str, chroms: list[str], sample: str | None = None,
              reference_names: list[str] | None = None) -> tuple[str | None, list[str]]:
    """Checks of the SNP VCF before the BAM pass. Returns (sample column to use, warnings).

    Fails (SystemExit) when the VCF has no index; when none of `reference_names` (default `chroms`; the BAM's
    reference names, so that one chromosome of `extract` is judged against the whole BAM) has records in the VCF
    (contig naming mismatch, e.g. `2` vs `chr2`); when the VCF has several sample columns and `sample` is None; or
    when `sample` is not a sample of the VCF. A chromosome of `chroms` without records in an otherwise matching VCF
    is a warning (no SNPs there)."""
    import pysam
    warn = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        if vf.index is None:
            raise SystemExit(_no_index(vcf_path))
        indexed = set(vf.index)                 # contigs with records
        refs = list(reference_names) if reference_names is not None else list(chroms)
        if refs and not any(c in indexed for c in refs):
            seen = sorted(indexed)
            alias = [c for c in refs if _alias(c) in indexed]
            hint = (f"; the names differ by the `chr` prefix ({alias[0]} here, {_alias(alias[0])} in the VCF): rename "
                    f"the VCF contigs, e.g. `bcftools annotate --rename-chrs`" if alias else "")
            raise SystemExit(f"--snps {vcf_path}: none of the chromosomes ({', '.join(refs[:4])}"
                             f"{', ...' if len(refs) > 4 else ''}) has records in the VCF (VCF contigs with records: "
                             f"{', '.join(seen[:4]) or 'none'}{', ...' if len(seen) > 4 else ''}){hint}")
        for c in chroms:
            if c not in indexed:
                warn.append(f"--snps {vcf_path} has no records on {c}: no SNPs there")
        samples = list(vf.header.samples)
    if sample is not None:
        if sample not in samples:
            raise SystemExit(f"--vcf-sample {sample} is not a sample of {vcf_path} (samples: {', '.join(samples) or 'none'})")
        return sample, warn
    if len(samples) > 1:
        raise SystemExit(f"--snps {vcf_path} has {len(samples)} samples ({', '.join(samples[:4])}"
                         f"{', ...' if len(samples) > 4 else ''}): name the one of this BAM with --vcf-sample")
    return (samples[0] if samples else None), warn


def load_snps(vcf_path: str, chrom: str, phased: bool = False, region: tuple[int, int] | None = None,
              sample: str | None = None) -> SnpSet:
    """Het SNVs of `chrom` (optionally `region`, 0-based half-open) from an indexed VCF. `sample` names the sample
    column (required when there are several). A contig without records gives an empty set."""
    import pysam
    pos, ref, alt, hoa, ps, ids = [], [], [], [], [], []
    n_phased_gt = n_contig_ps = 0
    with pysam.VariantFile(str(vcf_path)) as vf:
        if vf.index is None:
            raise SystemExit(_no_index(vcf_path))
        samples = list(vf.header.samples)
        si = None
        if samples:
            if sample is None and len(samples) > 1:
                raise SystemExit(f"--snps {vcf_path} has {len(samples)} samples: name one with --vcf-sample")
            si = sample if sample is not None else samples[0]
        if chrom in vf.index:
            it = vf.fetch(chrom, region[0], region[1]) if region else vf.fetch(chrom)
        else:
            it = iter(())
        for rec in it:
            if rec.alts is None or len(rec.alts) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue
            flt = list(rec.filter.keys())
            if flt and flt != ["PASS"]:
                continue
            h, p = 0, -1
            if si is not None:
                s = rec.samples[si]
                gt = s.get("GT")
                if gt is None or len(gt) != 2 or None in gt or sorted(gt) != [0, 1]:
                    continue
                if s.phased:
                    if not phased:
                        n_phased_gt += 1
                    else:
                        h = 1 if gt[0] == 1 else 2
                        pv = s.get("PS") if "PS" in s else None
                        if pv is None:
                            p = CONTIG_PS; n_contig_ps += 1
                        else:
                            p = int(pv)
            pos.append(rec.pos - 1); ref.append(rec.ref.upper()); alt.append(rec.alts[0].upper()); hoa.append(h); ps.append(p)
            ids.append(rec.id or ".")
    pos = np.asarray(pos, dtype=np.int64)
    o = np.argsort(pos, kind="stable")
    return SnpSet(chrom, pos[o], [ref[i] for i in o], [alt[i] for i in o], np.asarray(hoa, dtype=np.int8)[o],
                  np.asarray(ps, dtype=np.int64)[o], [ids[i] for i in o], n_phased_gt=n_phased_gt,
                  n_contig_ps=n_contig_ps)


class AlleleAccumulator:
    """Collects (cell, snp, allele) observations and per-read haplotype votes."""

    VOTE_PS = ("first", "majority")   # phase-set rule of the per-read vote (see vote(); cli --vote-ps)

    def __init__(self, snps: SnpSet, min_bq: int = 10, vote_ps: str = "first"):
        if vote_ps not in self.VOTE_PS:
            raise ValueError(f"vote_ps must be one of {', '.join(self.VOTE_PS)}, not {vote_ps!r}")
        self.snps = snps
        self.min_bq = min_bq
        self.vote_ps = vote_ps
        self.o_cell = array("i"); self.o_snp = array("i"); self.o_alt = array("b")
        self.n_other = np.zeros(len(snps), dtype=np.int64)
        self.n_lowbq = np.zeros(len(snps), dtype=np.int64)

    def observe(self, r, start: int, end: int) -> list[tuple[int, int]]:
        """[(snp_index, allele)] for SNPs on the read: allele 0 ref, 1 alt, 2 other/deleted; low-BQ bases dropped."""
        pos = self.snps.pos
        lo = int(np.searchsorted(pos, start)); hi = int(np.searchsorted(pos, end))
        if hi <= lo:
            return []
        cig = np.asarray(r.cigartuples, dtype=np.int64).reshape(-1, 2)
        ops, lens = cig[:, 0], cig[:, 1]
        rlen = lens * REF_CONSUME[ops]; qlen = lens * Q_CONSUME[ops]
        rend = start + np.cumsum(rlen); qend = np.cumsum(qlen)
        p = pos[lo:hi]
        k = np.searchsorted(rend, p, side="right")
        k = np.minimum(k, len(ops) - 1)
        ism = MATCH_OPS[ops[k]] & (p < rend[k])
        qpos = qend[k] - qlen[k] + (p - (rend[k] - rlen[k]))
        seq = r.query_sequence
        quals = r.query_qualities
        out = []
        for t in range(hi - lo):
            s = lo + t
            if not ism[t]:
                self.n_other[s] += 1
                out.append((s, 2))
                continue
            q = int(qpos[t])
            if quals is not None and quals[q] < self.min_bq:
                self.n_lowbq[s] += 1
                continue
            b = seq[q]
            if b == self.snps.ref[s]:
                out.append((s, 0))
            elif b == self.snps.alt[s]:
                out.append((s, 1))
            else:
                self.n_other[s] += 1
                out.append((s, 2))
        return out

    def vote(self, obs) -> tuple[int, int, int]:
        """(hap, n_snp_obs, ps) from observations: hap 0 none, 1 A, 2 B, 3 ambiguous.

        `n_snp_obs` counts the ref / alt observations of all SNPs (phased or not). vote_ps "first": the haplotype is
        the majority over all phased SNPs whatever their phase set, `ps` the phase set of the first phased SNP that
        has one. vote_ps "majority": only the phased SNPs of the phase set with the most phased observations vote and
        `ps` is that phase set; a tie between phase sets is ambiguous (hap 3, `ps` the first of the tied sets).
        """
        if self.vote_ps == "majority":
            return self._vote_majority(obs)
        na = nb = n = 0; ps = -1
        hoa = self.snps.hap_of_alt; pss = self.snps.ps
        for s, a in obs:
            if a == 2:
                continue
            n += 1
            h = hoa[s]
            if h == 0:
                continue
            if (a == 1) == (h == 1):
                na += 1
            else:
                nb += 1
            if ps == -1:
                ps = int(pss[s])
        hap = 0 if na + nb == 0 else (1 if na > nb else (2 if nb > na else 3))
        return hap, n, ps

    def _vote_majority(self, obs) -> tuple[int, int, int]:
        """vote() with vote_ps "majority": the vote of the phase set holding most phased observations."""
        n = 0
        tally = {}                      # phase set -> [n hapA, n hapB], in order of first observation
        hoa = self.snps.hap_of_alt; pss = self.snps.ps
        for s, a in obs:
            if a == 2:
                continue
            n += 1
            h = hoa[s]
            if h == 0:
                continue
            t = tally.setdefault(int(pss[s]), [0, 0])
            t[0 if (a == 1) == (h == 1) else 1] += 1
        if not tally:
            return 0, n, -1
        best = max(a + b for a, b in tally.values())
        top = [p for p, (a, b) in tally.items() if a + b == best]
        ps = top[0]
        if len(top) > 1:
            return 3, n, ps
        na, nb = tally[ps]
        return (1 if na > nb else (2 if nb > na else 3)), n, ps

    def process(self, r, ci: int, start: int, end: int) -> tuple[int, int, int]:
        obs = self.observe(r, start, end)
        if not obs:
            return 0, 0, -1
        for s, a in obs:
            if a < 2:
                self.o_cell.append(ci); self.o_snp.append(s); self.o_alt.append(a)
        return self.vote(obs)

    def matrices(self, n_cells: int) -> dict:
        """AD (alt) and DP (ref + alt) as cells x snps CSR."""
        c = np.frombuffer(self.o_cell, dtype=np.int32).astype(np.int64) if len(self.o_cell) else np.zeros(0, np.int64)
        s = np.frombuffer(self.o_snp, dtype=np.int32).astype(np.int64) if len(self.o_snp) else np.zeros(0, np.int64)
        a = np.frombuffer(self.o_alt, dtype=np.int8).astype(np.int64) if len(self.o_alt) else np.zeros(0, np.int64)
        shape = (n_cells, len(self.snps))
        dp = sp.coo_matrix((np.ones(len(c), dtype=np.int32), (c, s)), shape=shape).tocsr()
        ad = sp.coo_matrix((a.astype(np.int32), (c, s)), shape=shape).tocsr()
        dp.sum_duplicates(); ad.sum_duplicates(); ad.eliminate_zeros()
        return {"AD": ad, "DP": dp}
