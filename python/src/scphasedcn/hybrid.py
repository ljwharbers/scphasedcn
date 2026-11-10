"""Hybrid phasing: orient read-backed phase blocks with a population-phased VCF (one phase set per contig).

A read-backed phaser (e.g. longphase on the bulk long reads) phases heterozygous sites into blocks (GT `|` + PS) whose
orientation relative to each other is unknown. A population phaser (e.g. Beagle 5 or SHAPEIT with a reference panel)
phases the same sample's genotypes along the whole chromosome, with occasional switch errors. `hybrid_phase` puts every
read-backed block into the population haplotype frame by majority vote and merges the oriented blocks into one
chromosome-wide phase set per contig, which the normal-tissue mode of the scphasedcn caller needs.

Votes. A site votes when it is a heterozygous biallelic SNV in both inputs with the same REF and ALT, phased in the
population file and inside a block (phased GT with a PS) in the blocks file; population records flagged imputed (INFO
flag `imputed_flag`, default `IMP`, as Beagle writes) never vote. The vote is "agree" (A) when the block's first GT
allele equals the population's first allele, else "disagree" (D).

Orientation. Per block (contig + PS): with at least `min_votes` votes and a strict majority the block is oriented: GT
flipped when the majority disagrees, PS set to the contig's merged phase set (the smallest position of an oriented site,
moved up if a kept block already uses that value). Blocks without votes, with fewer than `min_votes` or with a tie keep
their own GT and PS. With `fill_unphased`, heterozygous SNVs outside any block (unphased, or phased without PS) take the
population phase and the merged phase set when the population file has them (and the contig has one). Every other
record and FORMAT field is written unchanged.

INFO tags of the output: PHASE_SRC (O = block oriented by the population phase, PS merged; B = block kept with its own
PS; P = population phase of a site outside the blocks) and POP_VOTE (A / D = the site's population vote relative to the
output phase, . = no vote).

Reports. Per block: sites, votes, agree / disagree, switches (consecutive votes that change inside the block),
oriented, flipped, output PS. Per contig: population markers used, blocks oriented / flipped, the switch-disagreement
rate (switches / consecutive vote pairs inside blocks) and the minority fraction (votes against the majority in oriented
blocks / votes in oriented blocks). Both measure how often the population phase disagrees with the reads: the
switch rate is roughly the population phaser's switch-error rate at the spacing of the voting sites.

Inputs. `blocks`: the read-backed phased VCF (plain or bgzip, indexed or not; sorted). `population`: one phased VCF, or a
pattern with `{chrom}` for one file per contig (the contig name is tried as in the blocks file, then with / without the
`chr` prefix). Contig names may differ by the `chr` prefix between the two inputs; population files may be plain or
gzip / bgzip, with or without an index (an unindexed single file is read once into memory).
"""
from __future__ import annotations

import contextlib
import gzip
import os
import shutil
import sys
import tempfile
from collections import Counter, OrderedDict

SRC_TAG = "PHASE_SRC"
VOTE_TAG = "POP_VOTE"


def _bare(c: str) -> str:
    """Contig name without a leading `chr` (case-insensitive)."""
    return c[3:] if c[:3].lower() == "chr" else c


def _names(c: str) -> list[str]:
    """The contig name as given, then the other `chr` form."""
    alt = _bare(c) if _bare(c) != c else "chr" + c
    return [c, alt]


def _snv(rec) -> bool:
    return (rec.alts is not None and len(rec.alts) == 1 and len(rec.ref) == 1 and len(rec.alts[0]) == 1
            and rec.ref.upper() in "ACGT" and rec.alts[0].upper() in "ACGT")


def _het(gt) -> bool:
    return gt is not None and len(gt) == 2 and None not in gt and sorted(gt) == [0, 1]


def _sample_index(vf, sample: str | None, what: str, prefer: str | None = None) -> str:
    """Sample name to read: `sample`, else `prefer` when present, else the only sample."""
    samples = list(vf.header.samples)
    if not samples:
        raise SystemExit(f"hybrid-phase: the {what} VCF has no sample column")
    if sample is not None:
        if sample not in samples:
            raise SystemExit(f"hybrid-phase: sample {sample!r} not in the {what} VCF (samples: {', '.join(samples[:5])}...)")
        return sample
    if prefer is not None and prefer in samples:
        return prefer
    if len(samples) > 1:
        raise SystemExit(f"hybrid-phase: the {what} VCF has {len(samples)} samples; name one "
                         f"(--{'population-sample' if what == 'population' else 'sample'})")
    return samples[0]


@contextlib.contextmanager
def _open_vcf(path: str):
    """pysam.VariantFile of a plain, BGZF or BCF file; a plain-gzip VCF (which htslib cannot open) is decompressed
    to a temporary file first."""
    import pysam
    tmp = None
    try:
        try:
            vf = pysam.VariantFile(path)
        except (NotImplementedError, OSError, ValueError):
            with open(path, "rb") as fh:
                if fh.read(2) != b"\x1f\x8b":
                    raise
            fd, tmp = tempfile.mkstemp(suffix=".vcf", prefix="hybrid_phase_")
            with os.fdopen(fd, "wb") as dst, gzip.open(path, "rb") as src:
                shutil.copyfileobj(src, dst)
            vf = pysam.VariantFile(tmp)
        with vf:
            yield vf
    finally:
        if tmp is not None:
            os.unlink(tmp)


class Population:
    """Phased heterozygous SNV markers of the population VCF(s), per contig of the blocks file."""

    def __init__(self, spec: str, sample: str | None = None, prefer_sample: str | None = None,
                 imputed_flag: str | None = "IMP", log=None):
        self.spec = spec
        self.sample = sample
        self.prefer = prefer_sample
        self.flag = imputed_flag or None
        self.log = log or (lambda *a: None)
        self.pattern = "{chrom}" in spec
        self._cache = None          # unindexed single file: bare contig -> (markers, stats)
        if not self.pattern and not os.path.exists(spec):
            raise SystemExit(f"hybrid-phase: --population not found: {spec}")

    def path_for(self, contig: str) -> str | None:
        if not self.pattern:
            return self.spec
        for n in _names(contig):
            p = self.spec.format(chrom=n)
            if os.path.exists(p):
                return p
        return None

    def _parse(self, rec, si, st: Counter, out: dict):
        st["markers"] += 1
        if not _snv(rec):
            return
        if self.flag is not None and self.flag in rec.info:
            st["imputed"] += 1
            return
        s = rec.samples[si]
        gt = s.get("GT")
        if not _het(gt):
            st["not_het"] += 1
            return
        if not s.phased:
            st["unphased"] += 1
            return
        out[(rec.pos, rec.ref.upper(), rec.alts[0].upper())] = gt[0] == 1
        st["het_used"] += 1

    def _read(self, path: str, contig: str | None):
        """Markers of one contig (None: all contigs, keyed by bare name) from `path`."""
        with _open_vcf(path) as vf:
            si = _sample_index(vf, self.sample, "population", self.prefer)
            if contig is not None and getattr(vf, "index", None) is not None:
                it = None
                for n in _names(contig):
                    try:
                        it = vf.fetch(n)
                        break
                    except ValueError:
                        continue
                st, out = Counter(), {}
                for rec in it or ():
                    self._parse(rec, si, st, out)
                return out, st
            res = {}
            want = _bare(contig) if contig is not None else None
            for rec in vf:
                b = _bare(rec.chrom)
                if want is not None and b != want:
                    continue
                if b not in res:
                    res[b] = ({}, Counter())
                out, st = res[b]
                self._parse(rec, si, st, out)
            if contig is not None:
                return res.get(want, ({}, Counter()))
            return res

    def markers(self, contig: str):
        """(dict (pos, ref, alt) -> first allele is ALT, stats Counter, file name or None)."""
        path = self.path_for(contig)
        if path is None:
            return {}, Counter(), None
        if self.pattern:
            out, st = self._read(path, contig)
            return out, st, os.path.basename(path)
        if self._cache is None:
            with _open_vcf(path) as vf:
                indexed = getattr(vf, "index", None) is not None
            if indexed:
                out, st = self._read(path, contig)
                return out, st, os.path.basename(path)
            self.log(f"[population] {path} has no index: reading it once into memory")
            self._cache = self._read(path, None)
        out, st = self._cache.get(_bare(contig), ({}, Counter()))
        return out, st, os.path.basename(path)


def _decide(contig: str, sites: list, pop: dict, min_votes: int):
    """Orientation of the blocks of one contig.

    sites: (pos, ref, alt, ps or None, first_is_alt) of the contig's het SNVs (ps None = outside the blocks).
    Returns dict with flip (ps -> bool, oriented blocks only), merged (PS or None), votes ((pos, ref, alt) -> A / D
    relative to the input GT), fill ((pos, ref, alt) -> first allele ALT), block rows and contig counts."""
    blocks = OrderedDict()          # ps -> ([positions], [(position, vote)])
    votes = {}
    fill_cand = {}
    for pos, ref, alt, ps, first_alt in sites:
        key = (pos, ref, alt)
        if ps is None:
            if key in pop:
                fill_cand[key] = pop[key]
            continue
        poss, vs = blocks.setdefault(ps, ([], []))
        poss.append(pos)
        if key in pop:
            v = "A" if first_alt == pop[key] else "D"
            votes[key] = v
            vs.append((pos, v))
    flip, rows = {}, []
    n = Counter()
    for ps, (poss, vs) in blocks.items():
        vs = sorted(vs)
        n_a = sum(1 for _, v in vs if v == "A")
        n_d = len(vs) - n_a
        sw = sum(1 for i in range(1, len(vs)) if vs[i][1] != vs[i - 1][1])
        n["sites_in_blocks"] += len(poss)
        n["sites_voting"] += len(vs)
        n["switches"] += sw
        n["vote_pairs"] += max(len(vs) - 1, 0)
        ok = len(vs) >= max(min_votes, 1) and n_a != n_d
        if ok:
            flip[ps] = n_d > n_a
            n["sites_voting_oriented"] += len(vs)
            n["sites_minority"] += min(n_a, n_d)
        rows.append([contig, ps, len(poss), min(poss), max(poss), len(vs), n_a, n_d, sw, int(ok), int(flip.get(ps, False))])
    merged = None
    if flip:
        kept = {ps for ps in blocks if ps not in flip}
        merged = min(min(blocks[ps][0]) for ps in flip)
        while merged in kept:
            merged += 1
    for r in rows:
        r.append(merged if r[1] in flip else r[1])
    n["blocks"] = len(blocks)
    n["blocks_oriented"] = len(flip)
    n["blocks_flipped"] = sum(flip.values())
    return dict(flip=flip, merged=merged, votes=votes, fill=fill_cand if merged is not None else {}, rows=rows, n=n)


def _site(rec, si):
    """(pos, ref, alt, ps or None, first allele ALT) for a het biallelic SNV, else None."""
    if not _snv(rec):
        return None
    s = rec.samples[si]
    gt = s.get("GT")
    if not _het(gt):
        return None
    ps = s.get("PS") if s.phased else None
    return rec.pos, rec.ref.upper(), rec.alts[0].upper(), ps, gt[0] == 1


CHROM_KEYS = ["contig", "population_file", "population_markers", "population_imputed", "population_het_used",
              "records", "het_snvs", "blocks", "blocks_oriented", "blocks_flipped", "sites_in_blocks", "sites_voting",
              "sites_minority", "minority_frac", "switches", "vote_pairs", "switch_rate", "filled", "merged_ps"]
BLOCK_KEYS = ["contig", "ps", "n_sites", "start", "end", "n_votes", "agree", "disagree", "switches", "oriented",
              "flipped", "out_ps"]


def hybrid_phase(blocks: str, population: str, out: str, contigs: list[str] | None = None, sample: str | None = None,
                 population_sample: str | None = None, imputed_flag: str | None = "IMP", min_votes: int = 1,
                 fill_unphased: bool = False, block_tsv: str | None = None, chrom_tsv: str | None = None,
                 log=None) -> dict:
    """Write the hybrid-phased VCF to `out` (.vcf.gz + .tbi) and the optional reports. Returns the total counts."""
    import pysam
    log = log or (lambda *a: print(*a, file=sys.stderr, flush=True))
    if not os.path.exists(blocks):
        raise SystemExit(f"hybrid-phase: --blocks not found: {blocks}")
    if not out.endswith(".gz"):
        raise SystemExit(f"hybrid-phase: --out must end in .gz (bgzip VCF with a tabix index): {out}")
    if min_votes < 1:
        raise SystemExit("hybrid-phase: --min-votes must be at least 1")
    want = set(contigs) if contigs else None

    # pass 1: the het SNVs of each contig -> orientation decisions (contig by contig; the VCF must be sorted)
    vf = pysam.VariantFile(blocks)
    si = _sample_index(vf, sample, "blocks")
    pop = Population(population, sample=population_sample, prefer_sample=si, imputed_flag=imputed_flag, log=log)
    dec, rows, chrom_rows, seen = {}, [], [], []
    n_rec = Counter()

    def finish(c, sites):
        mk, st, fname = pop.markers(c)
        if fname is None:
            log(f"[warn] {c}: no population file ({population}); blocks kept as they are")
        d = _decide(c, sites, mk, min_votes)
        dec[c] = d
        rows.extend(d["rows"])
        n = d["n"]
        chrom_rows.append(OrderedDict(
            contig=c, population_file=fname or "missing", population_markers=st["markers"],
            population_imputed=st["imputed"], population_het_used=st["het_used"], records=n_rec[c], het_snvs=len(sites),
            blocks=n["blocks"], blocks_oriented=n["blocks_oriented"], blocks_flipped=n["blocks_flipped"],
            sites_in_blocks=n["sites_in_blocks"], sites_voting=n["sites_voting"], sites_minority=n["sites_minority"],
            minority_frac=round(n["sites_minority"] / n["sites_voting_oriented"], 5) if n["sites_voting_oriented"] else 0,
            switches=n["switches"], vote_pairs=n["vote_pairs"],
            switch_rate=round(n["switches"] / n["vote_pairs"], 5) if n["vote_pairs"] else 0,
            filled=len(d["fill"]) if fill_unphased else 0, merged_ps=d["merged"] if d["merged"] is not None else "."))

    cur, sites = None, []
    for rec in vf:
        c = rec.chrom
        if want is not None and c not in want:
            continue
        if c != cur:
            if cur is not None:
                finish(cur, sites)
            if c in dec:
                raise SystemExit(f"hybrid-phase: --blocks is not sorted by contig ({c} appears twice)")
            cur, sites = c, []
            seen.append(c)
        n_rec[c] += 1
        x = _site(rec, si)
        if x is not None:
            sites.append(x)
    if cur is not None:
        finish(cur, sites)
    vf.close()
    if not seen:
        raise SystemExit("hybrid-phase: no records in --blocks" + (f" on {', '.join(sorted(want))}" if want else ""))
    if pop.pattern and all(r["population_file"] == "missing" for r in chrom_rows):
        raise SystemExit(f"hybrid-phase: no population file matches {population} for any contig")
    if sum(r["sites_voting"] for r in chrom_rows) == 0:
        raise SystemExit("hybrid-phase: no overlapping sites on any contig (phased het SNVs with the same position, REF "
                         "and ALT in both inputs, not imputed): check the sample, the genome build, the contig names and "
                         "that the population file is phased")

    # pass 2: write
    vf = pysam.VariantFile(blocks)
    hdr = vf.header
    if SRC_TAG not in hdr.info:
        hdr.add_line(f'##INFO=<ID={SRC_TAG},Number=1,Type=String,Description="Hybrid phase source: O=read-backed block '
                     f'oriented by the population phase (PS merged per contig), B=read-backed block kept with its own PS '
                     f'(no population majority), P=population phase of a het outside the blocks">')
    if VOTE_TAG not in hdr.info:
        hdr.add_line(f'##INFO=<ID={VOTE_TAG},Number=1,Type=String,Description="Population phase vote of the site '
                     f'relative to the output phase: A=agree, D=disagree, .=no vote">')
    if "PS" not in hdr.formats:
        hdr.add_line('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">')
    hdr.add_line(f'##scphasedcn_hybrid_phase="population={os.path.basename(population)} min_votes={min_votes} '
                 f'imputed_flag={imputed_flag or "none"} fill_unphased={fill_unphased}"')
    ofh = pysam.VariantFile(out, "wz", header=hdr)
    tot = Counter(records=0, sites_oriented=0, sites_kept=0, sites_filled=0)
    for rec in vf:
        d = dec.get(rec.chrom)
        x = _site(rec, si) if d is not None else None
        if x is not None:
            pos, ref, alt, ps, _ = x
            key = (pos, ref, alt)
            s = rec.samples[si]
            if ps is not None:
                v = d["votes"].get(key, ".")
                if ps in d["flip"]:
                    if d["flip"][ps]:
                        gt = s["GT"]
                        s["GT"] = (gt[1], gt[0])
                        v = {"A": "D", "D": "A"}.get(v, v)
                    s.phased = True
                    s["PS"] = d["merged"]
                    rec.info[SRC_TAG] = "O"
                    tot["sites_oriented"] += 1
                else:
                    rec.info[SRC_TAG] = "B"
                    tot["sites_kept"] += 1
                rec.info[VOTE_TAG] = v
            elif fill_unphased and key in d["fill"]:
                s["GT"] = (1, 0) if d["fill"][key] else (0, 1)
                s.phased = True
                s["PS"] = d["merged"]
                rec.info[SRC_TAG] = "P"
                rec.info[VOTE_TAG] = "."
                tot["sites_filled"] += 1
        ofh.write(rec)
        tot["records"] += 1
    ofh.close()
    vf.close()
    pysam.tabix_index(out, preset="vcf", force=True)

    for r in chrom_rows:
        log(f"[{r['contig']}] " + " ".join(f"{k}={v}" for k, v in r.items() if k != "contig"))
    sums = Counter()
    for r in chrom_rows:
        for k, v in r.items():
            if isinstance(v, int) and k not in ("merged_ps",):
                sums[k] += v
    nvo = sum(dec[c]["n"]["sites_voting_oriented"] for c in dec)
    sums["minority_frac"] = round(sums["sites_minority"] / nvo, 5) if nvo else 0
    sums["switch_rate"] = round(sums["switches"] / sums["vote_pairs"], 5) if sums["vote_pairs"] else 0
    log("[total] " + " ".join(f"{k}={sums[k]}" for k in CHROM_KEYS[2:-1]) +
        f" sites_oriented={tot['sites_oriented']} sites_kept={tot['sites_kept']} sites_filled={tot['sites_filled']}")
    if chrom_tsv:
        with open(chrom_tsv, "w") as fh:
            fh.write("\t".join(CHROM_KEYS) + "\n")
            for r in chrom_rows:
                fh.write("\t".join(str(r[k]) for k in CHROM_KEYS) + "\n")
            fh.write("\t".join(["total", "."] + [str(sums[k]) for k in CHROM_KEYS[2:-1]] + ["."]) + "\n")
    if block_tsv:
        with open(block_tsv, "w") as fh:
            fh.write("\t".join(BLOCK_KEYS) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
    return dict(sums) | dict(tot)
