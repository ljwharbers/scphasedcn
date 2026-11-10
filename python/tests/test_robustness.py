"""Inputs that used to give silent empty or wrong output: they now fail loudly or are handled (expected values of the
synthetic BAM in conftest.bam)."""
import json
import shutil

import anndata as ad
import pandas as pd
import pytest

from conftest import CELL_A, CELL_B, name, rewrite_bam, write_vcf
from scphasedcn.cli import main


def _common(genome, bam, out, *extra):
    return ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(out),
            "--chrom-sizes", str(genome / "chrom.sizes"), "--no-gc", "--no-masks", *extra]


def _frags(out):
    return pd.read_csv(out / "fragments.tsv.gz", sep="\t")


# --------------------------------------------------------------------------------------------- C2 VCF index / contigs
def test_unindexed_vcf_fails_before_the_bam_pass(genome, bam, tmp_path):
    plain = genome / "snps.vcf"                       # the fixture's uncompressed text VCF
    gz = tmp_path / "noindex.vcf.gz"
    shutil.copy(genome / "snps.vcf.gz", gz)           # bgzipped, no .tbi
    for v in (plain, gz):
        with pytest.raises(SystemExit, match="no index"):
            main(["extract", "--chrom", "chr2", "--snps", str(v), "--phased"] + _common(genome, bam, tmp_path / "o"))
        with pytest.raises(SystemExit, match="no index"):
            main(["run", "--chroms", "chr2", "--jobs", "1", "--snps", str(v)] + _common(genome, bam, tmp_path / "o"))
    assert not (tmp_path / "o" / "partials").exists()


def test_contig_name_mismatch_fails_and_absent_contig_warns(genome, bam, tmp_path, capsys):
    nochr = write_vcf(tmp_path / "nochr.vcf", [("2", 1001, "A", "G", "PASS", ["0|1:1000"])],
                      contigs=(("1", 1_000_000), ("2", 500_000)))
    with pytest.raises(SystemExit, match="`chr` prefix"):
        main(["extract", "--chrom", "chr2", "--snps", str(nochr), "--phased"] + _common(genome, bam, tmp_path / "o"))
    assert not (tmp_path / "o" / "partials").exists()
    # chrX has no records in an otherwise matching VCF: a warning, the run continues (no SNPs there)
    vcf = genome / "snps.vcf.gz"
    assert main(["extract", "--chrom", "chrX", "--snps", str(vcf), "--phased"] + _common(genome, bam, tmp_path / "o")) == 0
    assert "has no records on chrX" in capsys.readouterr().err
    assert len(pd.read_csv(tmp_path / "o" / "partials" / "chrX" / "snps.tsv", sep="\t")) == 0


# --------------------------------------------------------------------------------------------- C3 missing partials
def test_merge_fails_on_missing_chromosomes(genome, bam, tmp_path):
    assert main(["extract", "--chrom", "chr1"] + _common(genome, bam, tmp_path)) == 0
    with pytest.raises(SystemExit, match="chr2, chrX"):
        main(["merge", "--out", str(tmp_path), "--chroms", "chr1,chr2,chrX"])
    assert not (tmp_path / "bins_100kb.h5ad").exists()
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr1,chr2,chrX", "--allow-missing"]) == 0
    assert json.load(open(tmp_path / "summary.json"))["missing"] == ["chr2", "chrX"]


# --------------------------------------------------------------------------------------------- C4 stale partials
def test_rerun_without_snps_removes_the_stale_allele_layer(genome, bam, vcf, tmp_path):
    assert main(["extract", "--chrom", "chr2", "--snps", str(vcf), "--phased"] + _common(genome, bam, tmp_path)) == 0
    assert (tmp_path / "partials" / "chr2" / "alleles.npz").exists()
    assert main(["extract", "--chrom", "chr2", "--no-snps"] + _common(genome, bam, tmp_path)) == 0
    assert not (tmp_path / "partials" / "chr2" / "alleles.npz").exists()
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr2"]) == 0
    assert not (tmp_path / "alleles_dna.h5ad").exists()
    done = json.load(open(tmp_path / "partials" / "chr2" / "DONE"))
    assert sorted(done["layers"]) == ["depth", "fragments"] and done["binsize"] == 100_000


def test_merge_checks_the_inputs_of_the_partials(genome, bam, vcf, tmp_path):
    c = lambda out, *x: _common(genome, bam, out, *x)
    # bin size
    o = tmp_path / "binsize"
    assert main(["extract", "--chrom", "chr1"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--binsize", "50000"] + c(o)) == 0
    with pytest.raises(SystemExit, match="binsize is 100000 on chr1 but 50000 on chr2"):
        main(["merge", "--out", str(o), "--chroms", "chr1,chr2"])
    # SNP file: chr1 and chr2 from different VCFs
    o = tmp_path / "vcf"
    other = write_vcf(tmp_path / "other.vcf", [("chr1", 5001, "A", "G", "PASS", ["0|1:5000"]),
                                                ("chr2", 1001, "A", "G", "PASS", ["0|1:1000"])])
    assert main(["extract", "--chrom", "chr1", "--snps", str(other), "--phased"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--snps", str(vcf), "--phased"] + c(o)) == 0
    with pytest.raises(SystemExit, match="inputs of the alleles layer: snps"):
        main(["merge", "--out", str(o), "--chroms", "chr1,chr2"])
    # cell list: a whitelist with one cell on chr2
    o = tmp_path / "cells"
    wl = tmp_path / "one.tsv"; wl.write_text(f"{CELL_A}\n")
    assert main(["extract", "--chrom", "chr1"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2"] + c(o) + ["--barcodes", str(wl)]) == 0
    with pytest.raises(SystemExit, match="different cell lists"):
        main(["merge", "--out", str(o), "--chroms", "chr1,chr2"])
    # an allele-only rerun at another bin size than the kept depth layer
    o = tmp_path / "rerun"
    assert main(["extract", "--chrom", "chr2"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--layers", "alleles", "--binsize", "50000", "--snps", str(vcf),
                 "--phased"] + c(o)) == 0
    with pytest.raises(SystemExit, match="bins"):
        main(["merge", "--out", str(o), "--chroms", "chr2"])
    # partials of an older version (empty DONE) still merge, with a warning
    o = tmp_path / "legacy"
    for ch in ("chr1", "chr2"):
        assert main(["extract", "--chrom", ch] + c(o)) == 0
    (o / "partials" / "chr1" / "DONE").write_text("")
    from scphasedcn.merge import merge_sample
    msgs = []
    merge_sample(str(o), ["chr1", "chr2"], log=msgs.append)
    assert any("no DONE metadata): chr1" in m for m in msgs) and (o / "bins_100kb.h5ad").exists()


# --------------------------------------------------------------------------------------------- C5 phased without PS
def test_phased_genotypes_without_ps_get_one_phase_set_per_contig(genome, bam, tmp_path, capsys):
    from scphasedcn.alleles import CONTIG_PS, load_snps
    rows = [("chr2", 1001, "A", "G", "PASS", ["0|1"]), ("chr2", 1501, "C", "T", "PASS", ["1|0"]),
            ("chr2", 2001, "G", "A", "PASS", ["0/1"])]
    v = write_vcf(tmp_path / "beagle.vcf", rows, fmt="GT")
    s = load_snps(str(v), "chr2", phased=True)
    assert s.ps.tolist() == [CONTIG_PS, CONTIG_PS, -1] and s.hap_of_alt.tolist() == [2, 1, 0] and CONTIG_PS > 0
    out = tmp_path / "o"
    assert main(["extract", "--chrom", "chr2", "--snps", str(v), "--phased"] + _common(genome, bam, out)) == 0
    f = pd.read_csv(out / "partials" / "chr2" / "fragments.tsv.gz", sep="\t", header=None)
    hap, ps = f[10], f[12]
    assert ((hap == 1) | (hap == 2)).sum() == 2 and (ps[(hap == 1) | (hap == 2)] == CONTIG_PS).all()
    # phased genotypes in --snps without --phased: warned
    capsys.readouterr()
    assert main(["extract", "--chrom", "chr2", "--snps", str(v)] + _common(genome, bam, tmp_path / "u")) == 0
    assert "--phased was not given" in capsys.readouterr().err


# --------------------------------------------------------------------------------------------- C6 GT / samples
def test_missing_allele_genotypes_are_skipped(tmp_path):
    from scphasedcn.alleles import load_snps
    rows = [("chr2", 1001, "A", "G", "PASS", [".|1:1000"]), ("chr2", 1501, "C", "T", "PASS", ["./.:."]),
            ("chr2", 2001, "G", "A", "PASS", ["0|1:1000"])]
    s = load_snps(str(write_vcf(tmp_path / "miss.vcf", rows)), "chr2", phased=True)
    assert (s.pos + 1).tolist() == [2001]


def test_multi_sample_vcf_needs_vcf_sample(genome, bam, tmp_path):
    rows = [("chr2", 1001, "A", "G", "PASS", ["0|0:1000", "0|1:1000"]),
            ("chr2", 1501, "C", "T", "PASS", ["1|1:1000", "1|0:1000"])]
    v = write_vcf(tmp_path / "two.vcf", rows, samples=("S1", "S2"))
    base = ["extract", "--chrom", "chr2", "--snps", str(v), "--phased"]
    with pytest.raises(SystemExit, match="2 samples .* --vcf-sample"):
        main(base + _common(genome, bam, tmp_path / "o"))
    with pytest.raises(SystemExit, match="not a sample"):
        main(base + ["--vcf-sample", "S3"] + _common(genome, bam, tmp_path / "o"))
    assert main(base + ["--vcf-sample", "S2"] + _common(genome, bam, tmp_path / "o")) == 0
    snps = pd.read_csv(tmp_path / "o" / "partials" / "chr2" / "snps.tsv", sep="\t")
    assert snps["pos"].tolist() == [1001, 1501] and snps["hap_of_alt"].tolist() == [2, 1]   # S1 is homozygous


# --------------------------------------------------------------------------------------------- C8 CB "-1" suffix
def test_cb_with_dash1_suffix_counts_for_the_listed_cell(genome, bam, bam_dash1, vcf, tmp_path):
    for b, o in ((bam, tmp_path / "plain"), (bam_dash1, tmp_path / "dash1")):
        for ch in ("chr1", "chr2"):
            assert main(["extract", "--chrom", ch, "--snps", str(vcf), "--phased"] + _common(genome, b, o)) == 0
        assert main(["merge", "--out", str(o), "--chroms", "chr1,chr2"]) == 0
    assert _frags(tmp_path / "plain").equals(_frags(tmp_path / "dash1"))
    a, b = (ad.read_h5ad(tmp_path / x / "bins_100kb.h5ad") for x in ("plain", "dash1"))
    assert (a.X != b.X).nnz == 0 and a.X.sum() > 0 and a.obs.equals(b.obs)


def test_whitelist_matching_no_record_fails(genome, bam, tmp_path, monkeypatch, capsys):
    import scphasedcn.bam
    wl = tmp_path / "other.tsv"; wl.write_text("TTTTTTTTTTTTTTTT\nACGTACGTACGTACGT\n")
    args = ["extract", "--chrom", "chr1"] + _common(genome, bam, tmp_path / "o") + ["--barcodes", str(wl)]
    assert main(args) == 0                                      # few unlisted records: a loud warning only
    assert "whitelist does not match" in capsys.readouterr().err
    monkeypatch.setattr(scphasedcn.bam, "MIN_UNLISTED_MISMATCH", 10)
    with pytest.raises(SystemExit, match="whitelist does not match"):
        main(args)
    assert not (tmp_path / "o" / "partials" / "chr1" / "DONE").exists()


# --------------------------------------------------------------------------------------------- C9 --layers qc, --region
def test_layers_qc_is_not_a_layer(genome, bam, tmp_path):
    with pytest.raises(SystemExit, match="nothing to extract"):
        main(["extract", "--chrom", "chr1", "--layers", "qc"] + _common(genome, bam, tmp_path))
    assert main(["extract", "--chrom", "chr1", "--layers", "depth,qc"] + _common(genome, bam, tmp_path)) == 0
    p = json.load(open(tmp_path / "partials" / "chr1" / "params.json"))
    assert p["layers"] == ["depth", "qc"] and (tmp_path / "partials" / "chr1" / "qc.tsv").exists()


def test_region_counts_records_by_start_and_has_its_own_partial(genome, bam, tmp_path):
    assert main(["extract", "--chrom", "chr1"] + _common(genome, bam, tmp_path)) == 0
    full = (tmp_path / "partials" / "chr1" / "qc.tsv").read_bytes()
    assert main(["extract", "--region", "chr1:10600-50000"] + _common(genome, bam, tmp_path)) == 0
    rdir = tmp_path / "partials" / "region_chr1_10600-50000"
    q = pd.read_csv(rdir / "qc.tsv", sep="\t", index_col=0)
    assert q.loc[CELL_A, "n_records"] == 6 and q.loc[CELL_A, "n_pass"] == 6   # c d e f g h (not a, b: start before)
    assert q.loc[CELL_B, "n_records"] == 0
    p = json.load(open(rdir / "params.json"))
    assert p["unlisted"]["n_records"] == 2 and p["dash"]["n_records"] == 2
    # the chromosome's partial is untouched and still merges; the region partial is refused
    assert (tmp_path / "partials" / "chr1" / "qc.tsv").read_bytes() == full
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr1"]) == 0
    with pytest.raises(SystemExit, match="--region extraction"):
        main(["merge", "--out", str(tmp_path), "--chroms", "region_chr1_10600-50000"])


# --------------------------------------------------------------------------------------------- extras
def test_piece_numbers_beyond_int16(genome, bam, tmp_path):
    def edit(r):
        if r.reference_name == "chr1" and r.reference_start == 10000 and r.get_tag("CB") == CELL_A:
            r.query_name = name(CELL_A, "big", 40000, 40000)
        return r
    big = rewrite_bam(bam, tmp_path / "big.bam", edit)
    assert main(["extract", "--chrom", "chr1"] + _common(genome, big, tmp_path / "o")) == 0
    f = pd.read_csv(tmp_path / "o" / "partials" / "chr1" / "fragments.tsv.gz", sep="\t", header=None)
    row = f[(f[1] == 10000) & (f[3] == CELL_A)].iloc[0]
    assert row[8] == 40000 and row[9] == 40000
