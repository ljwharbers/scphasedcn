"""hybrid-phase (scphasedcn.hybrid): orientation of read-backed phase blocks by a population-phased VCF, on small
hand-made VCFs with known answers."""
from __future__ import annotations

import gzip

import pysam
import pytest

from scphasedcn.alleles import load_snps
from scphasedcn.cli import main
from scphasedcn.hybrid import hybrid_phase

HEAD = ("##fileformat=VCFv4.2\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n')
POP_HEAD = ("##fileformat=VCFv4.2\n"
            '##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

# Blocks (read-backed) on chr1:
#   PS 100  : 100 0|1, 150 A>G 0|1 (population has A>T: no vote), 200 1|0, 300 0|1; population says the opposite at
#             100, 200, 300 -> 3 disagree -> flipped
#   PS 1000 : 1000-1300 all 0|1; population 0|1, 0|1, 1|0, 1|0 -> tie -> kept (PS 1000)
#   PS 2000 : 2000, 2100 not in the population -> no vote -> kept (PS 2000)
#   PS 3000 : 3000 0|1 (population 0|1: agree), 3100 / 3200 0|1 (population 1|0 but IMP: ignored) -> oriented as is
#   4000 0/1 unphased (population 1|0: filled with --fill-unphased); 4500 indel; 4600 hom alt
# chr2: two blocks, both orientable (PS 10 agrees, PS 500 disagrees) -> one phase set (10)
BLOCK_ROWS = [
    ("chr1", 100, "A", "C", "0|1", 100), ("chr1", 150, "A", "G", "0|1", 100), ("chr1", 200, "G", "T", "1|0", 100),
    ("chr1", 300, "C", "T", "0|1", 100),
    ("chr1", 1000, "A", "C", "0|1", 1000), ("chr1", 1100, "A", "C", "0|1", 1000), ("chr1", 1200, "A", "C", "0|1", 1000),
    ("chr1", 1300, "A", "C", "0|1", 1000),
    ("chr1", 2000, "A", "C", "0|1", 2000), ("chr1", 2100, "A", "C", "1|0", 2000),
    ("chr1", 3000, "A", "C", "0|1", 3000), ("chr1", 3100, "A", "C", "0|1", 3000), ("chr1", 3200, "A", "C", "0|1", 3000),
    ("chr1", 4000, "A", "C", "0/1", None), ("chr1", 4500, "AT", "A", "0|1", 3000), ("chr1", 4600, "A", "C", "1/1", None),
    ("chr2", 10, "A", "C", "0|1", 10), ("chr2", 20, "A", "C", "1|0", 10),
    ("chr2", 500, "A", "C", "0|1", 500), ("chr2", 600, "A", "C", "0|1", 500),
]
# population: (bare contig, pos, ref, alt, gt, imputed)
POP_ROWS = [
    ("1", 100, "A", "C", "1|0", False), ("1", 150, "A", "T", "1|0", False), ("1", 200, "G", "T", "0|1", False),
    ("1", 300, "C", "T", "1|0", False),
    ("1", 1000, "A", "C", "0|1", False), ("1", 1100, "A", "C", "0|1", False), ("1", 1200, "A", "C", "1|0", False),
    ("1", 1300, "A", "C", "1|0", False),
    ("1", 3000, "A", "C", "0|1", False), ("1", 3100, "A", "C", "1|0", True), ("1", 3200, "A", "C", "1|0", True),
    ("1", 4000, "A", "C", "1|0", False), ("1", 5000, "A", "C", "0|1", False),
    ("2", 10, "A", "C", "0|1", False), ("2", 20, "A", "C", "1|0", False),
    ("2", 500, "A", "C", "1|0", False), ("2", 600, "A", "C", "1|0", False),
]


def write_blocks(path, rows=BLOCK_ROWS):
    lines = [HEAD, "##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n", COLS]
    for c, p, r, a, gt, ps in rows:
        fmt, val = ("GT:PS", f"{gt}:{ps}") if ps is not None else ("GT", gt)
        lines.append(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\t{fmt}\t{val}\n")
    txt = path.with_suffix("")
    txt.write_text("".join(lines))
    pysam.tabix_compress(str(txt), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def pop_text(rows, chr_prefix=False, contig_lines=True):
    cs = sorted({r[0] for r in rows}, key=str)
    nm = (lambda c: "chr" + c) if chr_prefix else (lambda c: c)
    head = POP_HEAD + ("".join(f"##contig=<ID={nm(c)},length=100000>\n" for c in cs) if contig_lines else "") + COLS
    body = "".join(f"{nm(c)}\t{p}\t.\t{r}\t{a}\t.\tPASS\t{'IMP' if imp else '.'}\tGT\t{gt}\n" for c, p, r, a, gt, imp in rows)
    return head + body


def write_pop(path, rows=POP_ROWS, chr_prefix=False, mode="bgzip"):
    """mode: bgzip (+ tabix), gzip (plain gzip, no index) or plain (.vcf)."""
    txt = pop_text(rows, chr_prefix)
    if mode == "plain":
        path.write_text(txt)
    elif mode == "gzip":
        with gzip.open(path, "wt") as fh:
            fh.write(txt)
    else:
        raw = path.with_suffix("")
        raw.write_text(txt)
        pysam.tabix_compress(str(raw), str(path), force=True)
        pysam.tabix_index(str(path), preset="vcf", force=True)
    return path


def read_out(path):
    out = {}
    with pysam.VariantFile(str(path)) as vf:
        for rec in vf:
            s = rec.samples[0]
            out[(rec.chrom, rec.pos)] = dict(gt=s["GT"], phased=s.phased, ps=s.get("PS"),
                                             src=rec.info.get("PHASE_SRC"), vote=rec.info.get("POP_VOTE"))
    return out


def read_tsv(path):
    lines = path.read_text().splitlines()
    keys = lines[0].split("\t")
    return [dict(zip(keys, ln.split("\t"))) for ln in lines[1:]]


@pytest.fixture()
def inputs(tmp_path):
    return write_blocks(tmp_path / "blocks.vcf.gz"), write_pop(tmp_path / "pop.vcf.gz")


def test_orientation_tie_no_overlap_imputed(inputs, tmp_path):
    blocks, pop = inputs
    out = tmp_path / "hybrid.vcf.gz"
    tot = hybrid_phase(str(blocks), str(pop), str(out), block_tsv=str(tmp_path / "b.tsv"), chrom_tsv=str(tmp_path / "c.tsv"),
                       log=lambda *a: None)
    assert (tmp_path / "hybrid.vcf.gz.tbi").exists()
    o = read_out(out)
    # PS 100 flipped into the population frame, merged PS = first oriented site = 100
    assert o[("chr1", 100)]["gt"] == (1, 0) and o[("chr1", 200)]["gt"] == (0, 1) and o[("chr1", 300)]["gt"] == (1, 0)
    assert o[("chr1", 150)]["gt"] == (1, 0)                     # flipped with its block, no vote (allele mismatch)
    assert o[("chr1", 150)]["vote"] == "." and o[("chr1", 100)]["vote"] == "A"
    assert all(o[("chr1", p)]["ps"] == 100 and o[("chr1", p)]["src"] == "O" for p in (100, 150, 200, 300))
    # tie: kept as is, votes relative to the block
    assert all(o[("chr1", p)]["ps"] == 1000 and o[("chr1", p)]["src"] == "B" for p in (1000, 1100, 1200, 1300))
    assert [o[("chr1", p)]["vote"] for p in (1000, 1100, 1200, 1300)] == ["A", "A", "D", "D"]
    assert o[("chr1", 1000)]["gt"] == (0, 1)
    # no overlap: kept
    assert o[("chr1", 2100)]["ps"] == 2000 and o[("chr1", 2100)]["gt"] == (1, 0) and o[("chr1", 2100)]["vote"] == "."
    # imputed markers never vote: one agreeing vote orients PS 3000 unflipped into the merged set
    assert all(o[("chr1", p)]["ps"] == 100 and o[("chr1", p)]["gt"] == (0, 1) for p in (3000, 3100, 3200))
    assert [o[("chr1", p)]["vote"] for p in (3000, 3100, 3200)] == ["A", ".", "."]
    # the indel of PS 3000 follows its oriented block into the merged phase set (unflipped, no vote)
    assert o[("chr1", 4500)]["ps"] == 100 and o[("chr1", 4500)]["gt"] == (0, 1) and o[("chr1", 4500)]["src"] == "O"
    # untouched records: unphased het without --fill-unphased, hom alt
    assert o[("chr1", 4000)]["phased"] is False and o[("chr1", 4000)]["src"] is None
    assert o[("chr1", 4600)]["gt"] == (1, 1)
    # chr2: PS 500 flipped, both blocks in one phase set
    assert {o[("chr2", p)]["ps"] for p in (10, 20, 500, 600)} == {10}
    assert o[("chr2", 500)]["gt"] == (1, 0) and o[("chr2", 10)]["gt"] == (0, 1)
    # reports
    b = {(r["contig"], r["ps"]): r for r in read_tsv(tmp_path / "b.tsv")}
    assert (b[("chr1", "100")]["oriented"], b[("chr1", "100")]["flipped"], b[("chr1", "100")]["disagree"]) == ("1", "1", "3")
    assert (b[("chr1", "1000")]["oriented"], b[("chr1", "1000")]["switches"], b[("chr1", "1000")]["out_ps"]) == ("0", "1", "1000")
    assert (b[("chr1", "2000")]["n_votes"], b[("chr1", "2000")]["oriented"]) == ("0", "0")
    assert (b[("chr1", "3000")]["n_votes"], b[("chr1", "3000")]["out_ps"]) == ("1", "100")
    c = {r["contig"]: r for r in read_tsv(tmp_path / "c.tsv")}
    assert set(c) == {"chr1", "chr2", "total"}
    assert (c["chr1"]["blocks"], c["chr1"]["blocks_oriented"], c["chr1"]["blocks_flipped"]) == ("4", "2", "1")
    assert c["chr1"]["population_imputed"] == "2" and c["chr1"]["merged_ps"] == "100"
    # chr1 votes: PS 100 D D D, PS 1000 A A D D (1 switch / 3 pairs), PS 3000 A -> switch rate 1 / 5
    assert c["chr1"]["sites_voting"] == "8" and float(c["chr1"]["switch_rate"]) == pytest.approx(0.2)
    assert float(c["chr1"]["minority_frac"]) == 0
    assert c["total"]["blocks_oriented"] == "4"
    assert tot["sites_filled"] == 0


def test_fill_unphased_and_min_votes(inputs, tmp_path):
    blocks, pop = inputs
    out = tmp_path / "fill.vcf.gz"
    hybrid_phase(str(blocks), str(pop), str(out), fill_unphased=True, min_votes=2, log=lambda *a: None)
    o = read_out(out)
    assert o[("chr1", 4000)]["gt"] == (1, 0) and o[("chr1", 4000)]["phased"] and o[("chr1", 4000)]["ps"] == 100
    assert o[("chr1", 4000)]["src"] == "P"
    # min_votes 2: PS 3000 (one vote) keeps its own PS
    assert o[("chr1", 3000)]["ps"] == 3000 and o[("chr1", 3000)]["src"] == "B"
    assert o[("chr1", 100)]["ps"] == 100 and o[("chr1", 100)]["src"] == "O"


@pytest.mark.parametrize("mode", ["plain", "gzip", "bgzip"])
def test_chr_prefix_mismatch_and_file_forms(tmp_path, mode):
    blocks = write_blocks(tmp_path / "blocks.vcf.gz")
    ext = {"plain": ".vcf", "gzip": ".vcf.gz", "bgzip": ".vcf.gz"}[mode]
    pop = write_pop(tmp_path / f"pop{ext}", mode=mode)                    # contigs "1", "2" vs "chr1", "chr2"
    out = tmp_path / "h.vcf.gz"
    hybrid_phase(str(blocks), str(pop), str(out), log=lambda *a: None)
    o = read_out(out)
    assert o[("chr1", 100)]["gt"] == (1, 0) and o[("chr2", 600)]["ps"] == 10


def test_population_pattern_per_contig(tmp_path):
    # blocks without "chr", per-contig files named with "chr" and records with "chr"
    rows = [(c.replace("chr", ""), *rest) for c, *rest in BLOCK_ROWS]
    lines = [HEAD, "##contig=<ID=1,length=100000>\n##contig=<ID=2,length=100000>\n", COLS]
    for c, p, r, a, gt, ps in rows:
        fmt, val = ("GT:PS", f"{gt}:{ps}") if ps is not None else ("GT", gt)
        lines.append(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\t{fmt}\t{val}\n")
    blocks = tmp_path / "b.vcf"
    blocks.write_text("".join(lines))                                       # plain, unindexed
    for c in ("1", "2"):
        raw = tmp_path / f"pop_chr{c}.vcf"
        raw.write_text(pop_text([r for r in POP_ROWS if r[0] == c], chr_prefix=True))
        pysam.tabix_compress(str(raw), str(raw) + ".gz", force=True)
        pysam.tabix_index(str(raw) + ".gz", preset="vcf", force=True)
    out = tmp_path / "h.vcf.gz"
    assert main(["hybrid-phase", "--blocks", str(blocks), "--population", str(tmp_path / "pop_{chrom}.vcf.gz"),
                 "--out", str(out), "--chrom-tsv", str(tmp_path / "c.tsv")]) == 0
    o = read_out(out)
    assert o[("1", 100)]["gt"] == (1, 0) and o[("2", 500)]["gt"] == (1, 0)
    assert {r["population_file"] for r in read_tsv(tmp_path / "c.tsv")} == {"pop_chr1.vcf.gz", "pop_chr2.vcf.gz", "."}


def test_imputed_flag_configurable(inputs, tmp_path):
    blocks, pop = inputs
    out = tmp_path / "noimp.vcf.gz"
    # with no imputed flag, the two imputed disagreeing markers outvote the one agreeing marker of PS 3000
    hybrid_phase(str(blocks), str(pop), str(out), imputed_flag=None, log=lambda *a: None)
    o = read_out(out)
    assert o[("chr1", 3000)]["gt"] == (1, 0) and o[("chr1", 3000)]["ps"] == 100


def test_errors(inputs, tmp_path):
    blocks, pop = inputs
    with pytest.raises(SystemExit, match="--blocks not found"):
        hybrid_phase(str(tmp_path / "nope.vcf.gz"), str(pop), str(tmp_path / "o.vcf.gz"))
    with pytest.raises(SystemExit, match="--population not found"):
        hybrid_phase(str(blocks), str(tmp_path / "nope.vcf.gz"), str(tmp_path / "o.vcf.gz"))
    with pytest.raises(SystemExit, match="no population file matches"):
        hybrid_phase(str(blocks), str(tmp_path / "nope_{chrom}.vcf.gz"), str(tmp_path / "o.vcf.gz"), log=lambda *a: None)
    with pytest.raises(SystemExit, match="must end in .gz"):
        hybrid_phase(str(blocks), str(pop), str(tmp_path / "o.vcf"))
    far = write_pop(tmp_path / "far.vcf.gz", rows=[("1", 99999, "A", "C", "0|1", False)])
    with pytest.raises(SystemExit, match="no overlapping sites on any contig"):
        hybrid_phase(str(blocks), str(far), str(tmp_path / "o.vcf.gz"), log=lambda *a: None)


def test_extractor_reads_one_phase_set_per_contig(inputs, tmp_path):
    blocks, pop = inputs
    out = tmp_path / "hybrid.vcf.gz"
    hybrid_phase(str(blocks), str(pop), str(out), log=lambda *a: None)
    s = load_snps(str(out), "chr2", phased=True)
    assert len(s) == 4 and set(s.ps.tolist()) == {10}
    # hapA = population haplotype 1: alt on hapA where the population GT is 1|0
    assert s.hap_of_alt.tolist() == [2, 1, 1, 1]
    s1 = load_snps(str(out), "chr1", phased=True)
    assert set(s1.ps[s1.hap_of_alt > 0].tolist()) == {100, 1000, 2000}       # tie and no-overlap blocks keep their PS
    before = load_snps(str(blocks), "chr2", phased=True)
    assert set(before.ps.tolist()) == {10, 500}


def test_non_snv_records_follow_their_oriented_block(tmp_path):
    # PS 100: SNVs at 100, 200, 300 disagree with the population -> flipped; the indel (150) and the multi-allelic
    # het (250) of the same block carry no vote but must be flipped with it and take the merged PS.
    # PS 1000: a tie -> kept; its indel (1150) stays as it was.
    rows = [("chr1", 100, "A", "C", "0|1", 100), ("chr1", 150, "AT", "A", "0|1", 100),
            ("chr1", 200, "G", "T", "1|0", 100), ("chr1", 250, "A", "C,G", "1|2", 100),
            ("chr1", 300, "C", "T", "0|1", 100),
            ("chr1", 1000, "A", "C", "0|1", 1000), ("chr1", 1100, "A", "C", "0|1", 1000),
            ("chr1", 1150, "AT", "A", "1|0", 1000)]
    blocks = write_blocks(tmp_path / "blocks.vcf.gz", rows=rows)
    pop = write_pop(tmp_path / "pop.vcf.gz", rows=[
        ("1", 100, "A", "C", "1|0", False), ("1", 200, "G", "T", "0|1", False), ("1", 300, "C", "T", "1|0", False),
        ("1", 1000, "A", "C", "0|1", False), ("1", 1100, "A", "C", "1|0", False)])
    out = tmp_path / "h.vcf.gz"
    hybrid_phase(str(blocks), str(pop), str(out), log=lambda *a: None)
    o = read_out(out)
    assert o[("chr1", 100)]["gt"] == (1, 0) and o[("chr1", 200)]["gt"] == (0, 1)
    # the indel and the multi-allelic het are flipped with their block and carry its merged PS
    assert o[("chr1", 150)]["gt"] == (1, 0) and o[("chr1", 150)]["ps"] == 100 and o[("chr1", 150)]["src"] == "O"
    assert o[("chr1", 250)]["gt"] == (2, 1) and o[("chr1", 250)]["ps"] == 100 and o[("chr1", 250)]["src"] == "O"
    assert o[("chr1", 150)]["vote"] == "." and o[("chr1", 250)]["vote"] == "."
    # the tie block keeps its own PS, and its indel is untouched
    assert o[("chr1", 1150)]["gt"] == (1, 0) and o[("chr1", 1150)]["ps"] == 1000 and o[("chr1", 1150)]["src"] is None


def test_contig_without_population_markers_warns(tmp_path):
    blocks = write_blocks(tmp_path / "blocks.vcf.gz")
    pop = write_pop(tmp_path / "pop1.vcf.gz", rows=[r for r in POP_ROWS if r[0] == "1"])   # no chr2 markers
    out = tmp_path / "h.vcf.gz"
    msgs = []
    tot = hybrid_phase(str(blocks), str(pop), str(out), log=lambda *a: msgs.append(" ".join(map(str, a))))
    warn = [m for m in msgs if m.startswith("[warn]")]
    assert len([m for m in warn if "chr2" in m]) == 1 and "no markers" in warn[0]
    assert not any("chr1" in m for m in warn)
    assert tot["contigs_without_markers"] == 1
    assert any(m.startswith("[total]") and "contigs_without_markers=1" in m for m in msgs)
    o = read_out(out)
    assert o[("chr2", 500)]["ps"] == 500 and o[("chr2", 500)]["src"] == "B"
