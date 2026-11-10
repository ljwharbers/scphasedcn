"""`scphasedcn run` (all chromosomes + merge), `select-hets`, and the early reference checks."""
import anndata as ad
import pysam
import pytest

from scphasedcn.cli import main


def _common(genome, bam, out):
    return ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(out),
            "--chrom-sizes", str(genome / "chrom.sizes"), "--fasta", str(genome / "genome.fa"),
            "--blacklist", str(genome / "blacklist.bed"), "--centromeres", str(genome / "centromeres.tsv"), "--no-trf"]


def test_run_equals_extract_plus_merge(genome, bam, vcf, tmp_path):
    a_out, b_out = tmp_path / "run", tmp_path / "steps"
    chroms = "chr1,chr2,chrX"
    assert main(["run", "--chroms", chroms, "--jobs", "2", "--snps", str(vcf), "--phased"] + _common(genome, bam, a_out)) == 0
    for c in chroms.split(","):
        assert main(["extract", "--chrom", c, "--snps", str(vcf), "--phased"] + _common(genome, bam, b_out)) == 0
    assert main(["merge", "--out", str(b_out), "--chroms", chroms]) == 0
    a = ad.read_h5ad(a_out / "bins_100kb.h5ad"); b = ad.read_h5ad(b_out / "bins_100kb.h5ad")
    assert (a.X != b.X).nnz == 0 and list(a.var_names) == list(b.var_names)
    assert a.var["gc"].notna().any()                      # the FASTA gives GC (the caller needs it)
    assert (a_out / "fragments.tsv.gz").exists() and (a_out / "fragments.tsv.gz.tbi").exists()
    assert (a_out / "cells_qc.tsv").exists()


def test_reference_checks_fail_before_the_bam_pass(genome, bam, tmp_path):
    base = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(tmp_path),
            "--chrom-sizes", str(genome / "chrom.sizes"), "--chrom", "chr1"]
    with pytest.raises(SystemExit, match="need --fasta"):
        main(["extract"] + base)
    with pytest.raises(SystemExit, match="not found"):
        main(["extract", "--fasta", str(tmp_path / "missing.fa")] + base)
    with pytest.raises(SystemExit, match="not found"):
        main(["extract", "--fasta", str(genome / "genome.fa"), "--trf", str(tmp_path / "missing.bed")] + base)
    assert not (tmp_path / "partials").exists()
    assert main(["extract", "--no-gc", "--no-masks"] + base) == 0


def test_select_hets(tmp_path):
    raw = tmp_path / "calls.vcf"
    rows = [  # pos, ref, alt, filter, gt, gq, dp
        (1001, "A", "G", "PASS", "0|1", 30, 40), (2001, "C", "T", "PASS", "1|0", 30, 50),
        (3001, "G", "A", "PASS", "0|1", 5, 40),            # low GQ
        (4001, "G", "A", "LowQual", "0|1", 30, 40),         # not PASS
        (5001, "GT", "A", "PASS", "0|1", 30, 40),           # indel
        (6001, "T", "C", "PASS", "1|1", 30, 40),            # hom
        (7001, "T", "C", "PASS", "0|1", 30, 10),            # low DP
        (8001, "T", "C", "PASS", "0|1", 30, 45),            # not in loci
    ]
    hdr = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##FILTER=<ID=LowQual,Description=\"x\">\n"
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"x\">\n"
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"x\">\n"
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    raw.write_text(hdr + "".join(f"chr1\t{p}\t.\t{r}\t{a}\t50\t{f}\t.\tGT:GQ:DP:PS\t{g}:{q}:{d}:1000\n"
                                 for p, r, a, f, g, q, d in rows))
    pysam.tabix_compress(str(raw), str(raw) + ".gz", force=True); pysam.tabix_index(str(raw) + ".gz", preset="vcf", force=True)
    loci = tmp_path / "loci.vcf"
    loci.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    + "".join(f"chr1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a, *_ in rows if p != 8001))
    pysam.tabix_compress(str(loci), str(loci) + ".gz", force=True); pysam.tabix_index(str(loci) + ".gz", preset="vcf", force=True)

    def selected(out):
        return [r.pos for r in pysam.VariantFile(str(out))]

    o1 = tmp_path / "het.vcf.gz"
    assert main(["select-hets", "--vcf", str(raw) + ".gz", "--out", str(o1), "--contigs", "chr1", "--dp-quantile", "1"]) == 0
    assert selected(o1) == [1001, 2001, 8001]               # no loci: all passing hets
    o2 = tmp_path / "het_loci.vcf.gz"
    assert main(["select-hets", "--vcf", str(raw) + ".gz", "--out", str(o2), "--contigs", "chr1", "--dp-quantile", "1",
                 "--loci", str(loci) + ".gz"]) == 0
    assert selected(o2) == [1001, 2001]
    assert (tmp_path / "het_loci.vcf.gz.tbi").exists()


def test_select_hets_two_records_at_one_position(tmp_path):
    hdr = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"x\">\n"
           "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    raw = tmp_path / "calls.vcf"
    raw.write_text(hdr + "chr1\t1001\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n"
                         "chr1\t1001\t.\tA\tT\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n"
                         "chr1\t2001\t.\tC\tT\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n")
    loci = tmp_path / "loci.vcf"
    loci.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t1001\t.\tA\tG,T\t.\t.\t.\nchr1\t2001\t.\tC\tT\t.\t.\t.\n")
    for p in (raw, loci):
        pysam.tabix_compress(str(p), str(p) + ".gz", force=True); pysam.tabix_index(str(p) + ".gz", preset="vcf", force=True)
    out = tmp_path / "het.vcf.gz"
    assert main(["select-hets", "--vcf", str(raw) + ".gz", "--out", str(out), "--contigs", "chr1", "--dp-quantile", "1",
                 "--loci", str(loci) + ".gz"]) == 0
    assert [(r.pos, r.alts[0]) for r in pysam.VariantFile(str(out))] == [(1001, "G"), (1001, "T"), (2001, "T")]
