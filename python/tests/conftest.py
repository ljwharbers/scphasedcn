"""Synthetic SPLONGGET-like BAM for the scphasedcn tests (built with pysam; see docstrings for the expected answers)."""
from __future__ import annotations

import hashlib
import os
import random
from pathlib import Path

import numpy as np
import pysam
import pytest

CHROMS = [("chr1", 1_000_000), ("chr2", 500_000), ("chrM", 16_569), ("chrX", 300_000)]
CELL_A = "AAAACCCCGGGGTTTT"
CELL_B = "CCCCGGGGTTTTAAAA"
CELL_C = "GGGGTTTTAAAACCCC"   # not whitelisted
UUIDS = {}


def uuid_for(key: str) -> str:
    """Deterministic uuid-looking string per key."""
    if key not in UUIDS:
        dig = hashlib.blake2b(f"u|{key}|{len(UUIDS)}".encode(), digest_size=8).digest()
        h = int.from_bytes(dig, "big") % (16 ** 12)
        base = f"{len(UUIDS):08x}-{(h >> 32) & 0xFFFF:04x}-4{(h >> 20) & 0xFFF:03x}-8{(h >> 8) & 0xFFF:03x}-{h & 0xFFFFFFFFFFFF:012x}"
        UUIDS[key] = base
    return UUIDS[key]


def cigar_ops(cigar: str):
    return [(int(n), c) for n, c in __import__("re").findall(r"(\d+)([MIDNSHP=X])", cigar)]


def cigar_lengths(cigar: str):
    ref = q = 0
    for n, c in cigar_ops(cigar):
        if c in "MDN=X":
            ref += n
        if c in "MIS=X":
            q += n
    return ref, q


def qpos_of(cigar: str, start: int, refpos: int) -> int | None:
    """Independent CIGAR walk: query index of a reference position (None if in a deletion)."""
    r, q = start, 0
    for n, c in cigar_ops(cigar):
        if c in "M=X":
            if r <= refpos < r + n:
                return q + (refpos - r)
            r += n; q += n
        elif c in "DN":
            if r <= refpos < r + n:
                return None
            r += n
        elif c in "IS":
            q += n
    return None


def make_record(header, name, chrom, start, cigar, strand="+", mapq=60, cb=None, flag_extra=0, bases=None, quals=None, tags=()):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = header.get_tid(chrom)
    a.reference_start = start
    a.cigarstring = cigar
    _, qlen = cigar_lengths(cigar)
    seq = list("A" * qlen)
    for pos, b in (bases or {}).items():
        seq[pos] = b
    a.query_sequence = "".join(seq)
    qq = np.full(qlen, 30, dtype=np.uint8)
    for pos, v in (quals or {}).items():
        qq[pos] = v
    a.query_qualities = pysam.qualitystring_to_array("".join(chr(x + 33) for x in qq))
    a.mapping_quality = mapq
    a.flag = (16 if strand == "-" else 0) | flag_extra
    a.set_tag("CB", cb if cb is not None else name.split("_")[0])
    for k, v in tags:
        a.set_tag(k, v)
    return a


def name(cb, key, piece=1, npieces=1, sign="+"):
    return f"{cb}_#{uuid_for(key)}_{sign}{piece}of{npieces}"


@pytest.fixture(scope="session")
def genome(tmp_path_factory):
    d = tmp_path_factory.mktemp("genome")
    rng = random.Random(1)
    fa = d / "genome.fa"
    with open(fa, "w") as fh:
        for c, n in CHROMS:
            fh.write(f">{c}\n")
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            for i in range(0, n, 80):
                fh.write(seq[i:i + 80] + "\n")
    pysam.faidx(str(fa))
    sizes = d / "chrom.sizes"
    sizes.write_text("".join(f"{c}\t{n}\n" for c, n in CHROMS))
    (d / "blacklist.bed").write_text("chr1\t900000\t910000\tHigh Signal Region\n")
    (d / "trf.bed").write_text("chr1\t700000\t700100\n")
    (d / "centromeres.tsv").write_text("chr\tstart\tcen.left.base\tcen.right.base\tend\nchr1\t1\t500001\t520000\t1000000\n")
    (d / "peaks.bed").write_text("chr1\t10000\t10500\tpeak1\nchr2\t1000\t1100\tpeak2\n")
    (d / "whitelist.tsv").write_text(f"{CELL_A}\n{CELL_B}-1\n")
    return d


@pytest.fixture(scope="session")
def bam(genome):
    """Synthetic BAM. Expected per-cell facts (whitelist = A, B):

    CELL_A chr1 (all +/- strands as listed, MAPQ 60 unless stated):
      bin 0      a=[10000,12000)+ b=[10500,13000)- c=[20000,22000)+ d=[20100,21000)+ e=[30000,32000)+
                 f=[30005,32500)- g=[40000,42000)+ h=[40200,42100)+ k=[60000,61000)+ l=[60960,62000)+
      same read  i=[50000,51000)+ (1of2) j=[50400,51500)+ (2of2) (two pieces of one raw read)
      bin edge   m=[199950,200150)- (start bin 1, midpoint bin 2, 5' end bin 2, overlaps bins 1 and 2)
      blacklist  n=[900100,901000)+ (dropped from depth, flagged)
      trf start  o=[700050,700300)+ ; MAPQ0 p=[300000,300500) ; dup q (flag 1024) ; secondary r (256) ;
      supplementary s (2048) of t=[400000,400800)+ which carries SA ; u=[120000,121000) mapq 45 (not 60)
    CELL_B chr1: v=[10000,12000)+ (same position as a, other cell), w=[800000,800500)-; chrM: 2 reads; chrX: 1 read.
    CELL_C (unlisted) and "-" barcodes: 2 reads each on chr1.
    chr2 SNP reads (VCF: 1001 A>G 0|1 PS1000; 1501 C>T 1|0 PS1000; 2001 G>A 0/1 unphased):
      R1 CELL_A [900,2105) 20S400M5D400M4I400M: SNP1 ref (hapA), SNP2 alt (hapA), SNP3 ref -> hap A, nsnp 3
      R2 CELL_B [950,1600) 650M: SNP1 alt (hapB), SNP2 ref (hapB) -> hap B, nsnp 2
      R3 CELL_A [1400,1600) 200M: SNP2 with BQ 5 -> ignored (nsnp 0, hap 0)
      R4 CELL_B [990,1100) 5M20D85M: SNP1 in a deletion -> other
      R5 CELL_A [1950,2100) 150M: SNP3 base T -> other
    """
    d = genome
    header = pysam.AlignmentHeader.from_references([c for c, _ in CHROMS], [n for _, n in CHROMS])
    R = []
    A, B, C = CELL_A, CELL_B, CELL_C
    mk = lambda *a, **k: R.append(make_record(header, *a, **k))
    # CELL_A chr1
    mk(name(A, "a"), "chr1", 10000, "2000M", "+")
    mk(name(A, "b"), "chr1", 10500, "2500M", "-")
    mk(name(A, "c"), "chr1", 20000, "2000M", "+")
    mk(name(A, "d"), "chr1", 20100, "900M", "+")
    mk(name(A, "e"), "chr1", 30000, "2000M", "+")
    mk(name(A, "f"), "chr1", 30005, "2495M", "-")
    mk(name(A, "g"), "chr1", 40000, "2000M", "+")
    mk(name(A, "h"), "chr1", 40200, "1900M", "+")
    mk(name(A, "ij", 1, 2, "+"), "chr1", 50000, "1000M", "+")
    mk(name(A, "ij", 2, 2, "-"), "chr1", 50400, "1100M", "+")
    mk(name(A, "k"), "chr1", 60000, "1000M", "+")
    mk(name(A, "l"), "chr1", 60960, "1040M", "+")
    mk(name(A, "m"), "chr1", 199950, "200M", "-")
    mk(name(A, "n"), "chr1", 900100, "900M", "+")
    mk(name(A, "o"), "chr1", 700050, "250M", "+")
    mk(name(A, "p"), "chr1", 300000, "500M", "+", mapq=0)
    mk(name(A, "q"), "chr1", 310000, "500M", "+", flag_extra=1024)
    mk(name(A, "r"), "chr1", 320000, "500M", "+", flag_extra=256)
    mk(name(A, "t"), "chr1", 400000, "800M", "+", tags=[("SA", "chr1,450000,+,800S200M,60,0;")])
    mk(name(A, "t"), "chr1", 450000, "800S200M", "+", flag_extra=2048, tags=[("SA", "chr1,400000,+,800M,60,0;")])
    mk(name(A, "u"), "chr1", 120000, "1000M", "+", mapq=45)
    # CELL_B
    mk(name(B, "v"), "chr1", 10000, "2000M", "+")
    mk(name(B, "w"), "chr1", 800000, "500M", "-")
    mk(name(B, "mt1"), "chrM", 100, "300M", "+")
    mk(name(B, "mt2"), "chrM", 5000, "300M", "-")
    mk(name(B, "x1"), "chrX", 100000, "1500M", "+")
    # unlisted and "-"
    mk(name(C, "c1"), "chr1", 15000, "500M", "+")
    mk(name(C, "c2"), "chr1", 16000, "500M", "+")
    mk(name("-", "d1"), "chr1", 17000, "500M", "+")
    mk(name("-", "d2"), "chr1", 18000, "500M", "+", flag_extra=1024)
    # SNP reads chr2
    cig1 = "20S400M5D400M4I400M"
    b1 = {qpos_of(cig1, 900, 1000): "A", qpos_of(cig1, 900, 1500): "T", qpos_of(cig1, 900, 2000): "G"}
    mk(name(A, "R1"), "chr2", 900, cig1, "+", bases=b1)
    mk(name(B, "R2"), "chr2", 950, "650M", "-", bases={qpos_of("650M", 950, 1000): "G", qpos_of("650M", 950, 1500): "C"})
    mk(name(A, "R3"), "chr2", 1400, "200M", "+", bases={100: "T"}, quals={100: 5})
    mk(name(B, "R4"), "chr2", 990, "5M20D85M", "+")
    mk(name(A, "R5"), "chr2", 1950, "150M", "+", bases={50: "T"})
    R.sort(key=lambda r: (r.reference_id, r.reference_start))
    path = d / "synthetic.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for r in R:
            fh.write(r)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="session")
def vcf(genome):
    p = genome / "snps.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        '##FILTER=<ID=LowQual,Description="low quality">\n'
        "##contig=<ID=chr1,length=1000000>\n##contig=<ID=chr2,length=500000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr2\t1001\t.\tA\tG\t50\tPASS\t.\tGT:PS\t0|1:1000\n"
        "chr2\t1501\t.\tC\tT\t50\tPASS\t.\tGT:PS\t1|0:1000\n"
        "chr2\t2001\t.\tG\tA\t50\tPASS\t.\tGT:PS\t0/1:.\n"
        "chr2\t3001\t.\tG\tA\t50\tLowQual\t.\tGT:PS\t0/1:.\n"
        "chr2\t4001\t.\tGT\tA\t50\tPASS\t.\tGT:PS\t0/1:.\n"
    )
    pysam.tabix_compress(str(p), str(p) + ".gz", force=True)
    pysam.tabix_index(str(p) + ".gz", preset="vcf", force=True)
    return Path(str(p) + ".gz")


def write_vcf(path, rows, samples=("SAMPLE",), contigs=(("chr1", 1_000_000), ("chr2", 500_000)), fmt="GT:PS",
              index=True):
    """Write a small VCF (bgzipped + tabix-indexed when `index`; plain text otherwise). `rows` are
    (chrom, pos, ref, alt, filter, [sample fields...]) with the sample fields formatted for `fmt`. Returns the path."""
    path = Path(path)
    head = ("##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
            '##FILTER=<ID=LowQual,Description="low quality">\n'
            + "".join(f"##contig=<ID={c},length={n}>\n" for c, n in contigs)
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + ("\tFORMAT\t" + "\t".join(samples) if samples else "")
            + "\n")
    body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t50\t{f}\t." + (f"\t{fmt}\t" + "\t".join(sf) if samples else "") + "\n"
                   for c, p, r, a, f, sf in rows)
    path.write_text(head + body)
    if not index:
        return path
    pysam.tabix_compress(str(path), str(path) + ".gz", force=True)
    pysam.tabix_index(str(path) + ".gz", preset="vcf", force=True)
    return Path(str(path) + ".gz")


def rewrite_bam(src, dst, edit):
    """Copy a BAM with `edit(record)` applied to every record (returns the record or None to drop it); indexed."""
    with pysam.AlignmentFile(str(src), "rb") as fi, pysam.AlignmentFile(str(dst), "wb", template=fi) as fo:
        for r in fi:
            r = edit(r)
            if r is not None:
                fo.write(r)
    pysam.index(str(dst))
    return Path(dst)


@pytest.fixture(scope="session")
def bam_dash1(bam, tmp_path_factory):
    """The synthetic BAM with a 10x-style "-1" suffix on every CB except "-" (the whitelist has it on CELL_B only)."""
    def edit(r):
        cb = r.get_tag("CB")
        if cb != "-":
            r.set_tag("CB", cb + "-1")
        return r
    return rewrite_bam(bam, tmp_path_factory.mktemp("dash1") / "dash1.bam", edit)


@pytest.fixture(scope="session")
def extracted(genome, bam, vcf, tmp_path_factory):
    """Run extract (default layers) on chr1, chr2, chrM, chrX and merge; return the output directory."""
    from scphasedcn.cli import main
    out = tmp_path_factory.mktemp("out")
    common = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(out), "--binsize", "100000",
              "--chrom-sizes", str(genome / "chrom.sizes"), "--fasta", str(genome / "genome.fa"),
              "--blacklist", str(genome / "blacklist.bed"), "--trf", str(genome / "trf.bed"),
              "--centromeres", str(genome / "centromeres.tsv"), "--peaks", str(genome / "peaks.bed"),
              "--snps", str(vcf), "--phased"]
    for c in ("chr1", "chr2", "chrM", "chrX"):
        assert main(["extract", "--chrom", c] + common) == 0
    assert main(["merge", "--out", str(out), "--chroms", "chr1,chr2,chrM,chrX"]) == 0
    return out
