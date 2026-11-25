"""Inputs that used to give silent empty or wrong output: they now fail loudly or are handled (expected values of the
synthetic BAM in conftest.bam)."""
import json
import shutil
from pathlib import Path

import anndata as ad
import pandas as pd
import pysam
import pytest

from conftest import CELL_A, CELL_B, name, rewrite_bam, write_vcf
from scphasedcn.cli import main, _run_one as _real_run_one   # before any monkeypatch


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



# --------------------------------------------------------------------------------------------- SNP inputs within a partial
def test_merge_refuses_fragments_and_alleles_from_different_snp_inputs(genome, bam, vcf, tmp_path):
    c = lambda out, *x: _common(genome, bam, out, *x)
    flip = write_vcf(tmp_path / "flip.vcf", [("chr2", 1001, "A", "G", "PASS", ["1|0:1000"]),
                                             ("chr2", 1501, "C", "T", "PASS", ["0|1:1000"])])
    # an allele-only rerun with another VCF: the fragments keep the phase of the first one
    o = tmp_path / "newvcf"
    assert main(["extract", "--chrom", "chr2", "--snps", str(vcf), "--phased"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--layers", "alleles", "--snps", str(flip), "--phased"] + c(o)) == 0
    with pytest.raises(SystemExit, match="chr2: .*snps.*re-extract chr2 with --layers alleles,fragments"):
        main(["merge", "--out", str(o), "--chroms", "chr2"])
    assert not (o / "summary.json").exists()
    # the documented rerun of both layers merges
    assert main(["extract", "--chrom", "chr2", "--layers", "alleles,fragments", "--snps", str(flip), "--phased"] + c(o)) == 0
    assert main(["merge", "--out", str(o), "--chroms", "chr2"]) == 0
    # fragments made without SNPs, alleles added later: every fragment would have hap 0
    o = tmp_path / "later"
    assert main(["extract", "--chrom", "chr2", "--layers", "depth,fragments"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--layers", "alleles", "--snps", str(vcf), "--phased"] + c(o)) == 0
    with pytest.raises(SystemExit, match="fragments layer was made without SNPs.*--layers alleles,fragments"):
        main(["merge", "--out", str(o), "--chroms", "chr2"])
    # a fragments layer kept from an extraction without metadata cannot be checked: warned, not refused
    done = json.load(open(o / "partials" / "chr2" / "DONE"))
    done["layers"]["fragments"]["inputs"] = None
    (o / "partials" / "chr2" / "DONE").write_text(json.dumps(done))
    from scphasedcn.merge import merge_sample
    msgs = []
    merge_sample(str(o), ["chr2"], log=msgs.append)
    assert any("fragments layer of chr2 was kept" in m for m in msgs)


def test_merge_checks_the_vote_ps_of_the_partials(genome, bam, vcf, tmp_path):
    c = lambda out, *x: _common(genome, bam, out, "--snps", str(vcf), "--phased", *x)
    # chromosomes extracted with and without --vote-ps majority
    o = tmp_path / "mix"
    assert main(["extract", "--chrom", "chr1"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--vote-ps", "majority"] + c(o)) == 0
    with pytest.raises(SystemExit, match='inputs of the alleles layer: vote_ps is "first" on chr1 but "majority" on chr2'):
        main(["merge", "--out", str(o), "--chroms", "chr1,chr2"])
    assert not (o / "summary.json").exists()
    # within a partial: an allele-only rerun with --vote-ps majority next to the kept default fragments
    o = tmp_path / "rerun"
    assert main(["extract", "--chrom", "chr2"] + c(o)) == 0
    assert main(["extract", "--chrom", "chr2", "--layers", "alleles", "--vote-ps", "majority"] + c(o)) == 0
    with pytest.raises(SystemExit, match="chr2: .*different SNP inputs \\(vote_ps\\)"):
        main(["merge", "--out", str(o), "--chroms", "chr2"])
    # partials written before vote_ps was recorded (no key) are the default: they merge with new default partials
    # and are refused next to --vote-ps majority
    o = tmp_path / "old"
    for ch in ("chr1", "chr2"):
        assert main(["extract", "--chrom", ch] + c(o)) == 0
    done = json.load(open(o / "partials" / "chr1" / "DONE"))
    for L in ("alleles", "fragments"):
        assert done["layers"][L]["inputs"].pop("vote_ps") == "first"
    (o / "partials" / "chr1" / "DONE").write_text(json.dumps(done))
    assert main(["merge", "--out", str(o), "--chroms", "chr1,chr2"]) == 0
    assert main(["extract", "--chrom", "chr2", "--vote-ps", "majority"] + c(o)) == 0
    with pytest.raises(SystemExit, match='vote_ps is "first" on chr1 but "majority" on chr2'):
        main(["merge", "--out", str(o), "--chroms", "chr1,chr2"])


def test_merge_after_moving_the_output_and_extracting_out_of_order(genome, bam, vcf, tmp_path):
    o = tmp_path / "o"
    for ch in ("chr2", "chr1"):
        assert main(["extract", "--chrom", ch, "--snps", str(vcf), "--phased"] + _common(genome, bam, o)) == 0
    shutil.move(str(o), str(tmp_path / "moved"))
    assert main(["merge", "--out", str(tmp_path / "moved"), "--chroms", "chr1,chr2"]) == 0



# --------------------------------------------------------------------------------------------- interrupted merge
def test_interrupted_merge_leaves_no_summary_and_no_mixed_fragments(genome, bam, vcf, tmp_path, monkeypatch):
    import scphasedcn.merge as M
    for ch in ("chr1", "chr2"):
        assert main(["extract", "--chrom", ch, "--snps", str(vcf), "--phased"] + _common(genome, bam, tmp_path)) == 0
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr1,chr2"]) == 0
    frag, tbi = (tmp_path / "fragments.tsv.gz").read_bytes(), (tmp_path / "fragments.tsv.gz.tbi").read_bytes()
    files = sorted(p.name for p in tmp_path.iterdir())

    def cut_concat(parts, out, header=None):        # killed while writing the fragments table
        with open(out, "wb") as fh:
            fh.write(open(parts[0], "rb").read()[:100])
        raise KeyboardInterrupt("killed")

    def failed_index(*a, **k):                       # killed while indexing it
        raise MemoryError("killed")

    for target, attr, fake, err in ((M, "concat_bgzf", cut_concat, KeyboardInterrupt),
                                    (M.pysam, "tabix_index", failed_index, MemoryError)):
        with monkeypatch.context() as mp:
            mp.setattr(target, attr, fake)
            with pytest.raises(err):
                M.merge_sample(str(tmp_path), ["chr1", "chr2"], log=lambda *a: None)
        assert not (tmp_path / "summary.json").exists()
        # the previous table and its index are intact; no temporary file is left
        assert (tmp_path / "fragments.tsv.gz").read_bytes() == frag
        assert (tmp_path / "fragments.tsv.gz.tbi").read_bytes() == tbi
        assert sorted(p.name for p in tmp_path.iterdir()) == [f for f in files if f != "summary.json"]
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr1,chr2"]) == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == files
    assert (tmp_path / "fragments.tsv.gz").read_bytes() == frag
    with pysam.TabixFile(str(tmp_path / "fragments.tsv.gz")) as t:
        assert len(list(t.fetch("chr2"))) > 0


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
    err = capsys.readouterr().err   # one phase set per contig: warned (long-range switch errors of population phasing)
    assert "warning: 2 phased SNPs on chr2 have no PS" in err and "switch errors" in err
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


# --------------------------------------------------------------------------------------------- reference consistency
def _fasta(path, seqs):
    with open(path, "w") as fh:
        for c, n in seqs:
            fh.write(f">{c}\n" + "".join("ACGT"[i % 4] for i in range(n)) + "\n")
    pysam.faidx(str(path))
    return path


def test_chromosome_missing_from_the_bam_fails_before_the_bam_pass(genome, bam, tmp_path):
    sizes = tmp_path / "sizes"
    sizes.write_text((genome / "chrom.sizes").read_text() + "chr5\t181538259\n")
    args = ["extract", "--chrom", "chr5"] + _common(genome, bam, tmp_path / "o") + ["--chrom-sizes", str(sizes)]
    with pytest.raises(SystemExit, match="chr5 is not a reference of the BAM"):
        main(args)
    assert not (tmp_path / "o" / "partials").exists()


def test_chromosome_missing_from_the_fasta_fails(genome, bam, tmp_path):
    fa = _fasta(tmp_path / "one.fa", [("1", 1_000_000)])           # other naming, chr2 absent
    base = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(tmp_path / "o"),
            "--chrom-sizes", str(genome / "chrom.sizes"), "--no-masks", "--fasta", str(fa)]
    with pytest.raises(SystemExit, match="chr2 is not in the FASTA"):
        main(["extract", "--chrom", "chr2"] + base)
    with pytest.raises(SystemExit, match="chr1 is not in the FASTA .*`chr` prefix"):
        main(["extract", "--chrom", "chr1"] + base)
    assert not (tmp_path / "o" / "partials").exists()


def test_chromosome_lengths_must_agree(genome, bam, tmp_path):
    sizes = tmp_path / "sizes"
    sizes.write_text((genome / "chrom.sizes").read_text().replace("chr1\t1000000", "chr1\t2000000"))
    base = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(tmp_path / "o"), "--no-masks"]
    with pytest.raises(SystemExit, match=r"chr1 .*2000000 bp in the chromosome sizes .* 1000000 bp in the BAM"):
        main(["extract", "--chrom", "chr1", "--no-gc", "--chrom-sizes", str(sizes)] + base)
    fa = _fasta(tmp_path / "short.fa", [("chr1", 900_000), ("chr2", 500_000)])
    with pytest.raises(SystemExit, match=r"chr1 .*900000 bp in the FASTA .* 1000000 bp in the BAM"):
        main(["extract", "--chrom", "chr1", "--fasta", str(fa), "--chrom-sizes", str(genome / "chrom.sizes")] + base)
    assert not (tmp_path / "o" / "partials").exists()
    # run: one error for all chromosomes, before any worker starts
    sizes.write_text("chr1\t2000000\nchr2\t600000\nchrM\t16569\nchrX\t300000\n")
    with pytest.raises(SystemExit, match=r"(?s)chr1 .*2000000.*chr2 .*600000"):
        main(["run", "--chroms", "chr1,chr2,chrX", "--jobs", "2", "--no-gc", "--chrom-sizes", str(sizes)] + base)
    assert not (tmp_path / "o" / "partials").exists()


def test_fasta_without_chrm_or_with_another_chrm(genome, bam, tmp_path, capsys):
    """chrM only: a FASTA without it, or with another chrM length, is a warning and chrM's bins have no GC (NaN, as
    the caller tolerates); the run goes on and merges. Any other chromosome keeps the hard error."""
    main_chroms = [("chr1", 1_000_000), ("chr2", 500_000), ("chrX", 300_000)]
    for x, extra in (("nochrM", []), ("otherchrM", [("chrM", 16_571)])):
        fa = _fasta(tmp_path / f"{x}.fa", main_chroms + extra)
        o = tmp_path / x
        capsys.readouterr()
        assert main(["run", "--jobs", "1", "--fasta", str(fa)] + _common(genome, bam, o)[:-2] + ["--no-masks"]) == 0
        err = capsys.readouterr().err
        assert ("chrM is not in the FASTA" if x == "nochrM" else "chrM has 16571 bp in the FASTA") in err
        a = ad.read_h5ad(o / "bins_100kb.h5ad")
        gc = a.var["gc"]
        m = (a.var["chrom"] == "chrM").values
        assert m.sum() == 1 and gc[m].isna().all() and gc[~m].notna().all()
        assert (a.var.loc[m, "n_N"] == 0).all()
        assert json.load(open(o / "summary.json"))["chroms"] == ["chr1", "chr2", "chrX", "chrM"]
    fa = _fasta(tmp_path / "nochr2.fa", [("chr1", 1_000_000), ("chrX", 300_000)])
    with pytest.raises(SystemExit, match="chr2 is not in the FASTA"):
        main(["run", "--jobs", "1", "--fasta", str(fa)] + _common(genome, bam, tmp_path / "o")[:-2] + ["--no-masks"])


# --------------------------------------------------------------------------------------------- select-hets loci naming
_HETS_HDR = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n"
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">\n"
             "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"x\">\n"
             "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n"
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")


def _bgzip(p):
    pysam.tabix_compress(str(p), str(p) + ".gz", force=True); pysam.tabix_index(str(p) + ".gz", preset="vcf", force=True)
    return str(p) + ".gz"


def test_select_hets_loci_with_other_contig_names(tmp_path):
    raw = tmp_path / "calls.vcf"
    raw.write_text(_HETS_HDR + "chr1\t1001\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n"
                               "chr1\t2001\t.\tC\tT\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n")

    def loci(name, contig):
        p = tmp_path / f"{name}.vcf"
        p.write_text(f"##fileformat=VCFv4.2\n##contig=<ID={contig},length=100000>\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                     f"{contig}\t1001\t.\tA\tG\t.\t.\t.\n{contig}\t2001\t.\tC\tT\t.\t.\t.\n")
        return _bgzip(p)
    base = ["select-hets", "--vcf", _bgzip(raw), "--contigs", "chr1", "--dp-quantile", "1"]
    out = tmp_path / "het.vcf.gz"
    assert main(base + ["--out", str(out), "--loci", loci("nochr", "1")]) == 0     # `1` in the loci file, `chr1` here
    assert [r.pos for r in pysam.VariantFile(str(out))] == [1001, 2001]
    out2 = tmp_path / "none.vcf.gz"                                                  # no loci contig matches at all
    (tmp_path / "none.vcf.gz.tbi").write_text("index of an earlier run")
    with pytest.raises(SystemExit, match="no autosomal site"):
        main(base + ["--out", str(out2), "--loci", loci("other", "chr9")])
    assert not out2.exists() and not (tmp_path / "none.vcf.gz.tbi").exists()


def test_select_hets_chrx_only_call_set(tmp_path, capsys):
    """A chrX-only call set whose header lists the autosomes too selects its chrX sites (no autosome had records);
    sites but none on an autosome that had records, or that --contigs names, is a loud warning, not an error. An
    error leaves no output and no --stats file of an earlier run."""
    hdr = _HETS_HDR.replace("##contig=<ID=chr2,length=100000>\n",
                            "##contig=<ID=chr2,length=100000>\n##contig=<ID=chrX,length=100000>\n")
    x = "".join(f"chrX\t{p}\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n" for p in (1001, 2001))
    raw = tmp_path / "x.vcf"
    raw.write_text(hdr + x)
    out, st = tmp_path / "het.vcf.gz", tmp_path / "stats.tsv"
    base = ["select-hets", "--vcf", _bgzip(raw), "--out", str(out), "--dp-quantile", "1", "--stats", str(st)]
    capsys.readouterr()
    assert main(base) == 0
    assert [(r.chrom, r.pos) for r in pysam.VariantFile(str(out))] == [("chrX", 1001), ("chrX", 2001)]
    assert "NO AUTOSOMAL SITE" not in capsys.readouterr().err and st.exists()
    assert main(base + ["--contigs", "chr1,chrX"]) == 0                                # an autosome named: warn
    assert "NO AUTOSOMAL SITE selected (2 sites on chrX), although --contigs names chr1" in capsys.readouterr().err
    raw2 = tmp_path / "x2.vcf"                                                         # chr1 records, none selected
    raw2.write_text(hdr + "chr1\t1001\t.\tA\tG\t50\tPASS\t.\tGT:GQ:DP\t0|1:5:40\n" + x)
    base2 = ["select-hets", "--vcf", _bgzip(raw2), "--out", str(out), "--dp-quantile", "1", "--stats", str(st)]
    assert main(base2) == 0
    assert "NO AUTOSOMAL SITE selected (2 sites on chrX), although chr1 had records" in capsys.readouterr().err
    assert st.exists() and out.exists()
    lo = tmp_path / "loci.vcf"                                                         # chrX only, none selected
    lo.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                  "chr1\t1001\t.\tA\tG\t.\t.\t.\n")
    with pytest.raises(SystemExit, match="no site selected"):
        main(base2 + ["--contigs", "chrX", "--loci", _bgzip(lo), "--no-loci-contigs", ""])
    assert not st.exists() and not out.exists() and not Path(f"{out}.tbi").exists()
    st.write_text("stats of an earlier run")
    with pytest.raises(SystemExit, match="no candidate sites"):
        main(base2 + ["--min-dp", "100"])
    assert not st.exists()


# --------------------------------------------------------------------------------------------- biallelic SNV test
def test_spanning_deletion_and_n_alleles_are_not_snvs(tmp_path):
    from scphasedcn.alleles import load_snps
    from scphasedcn.hets import is_candidate
    sites = [(1001, "A", "G"), (1101, "A", "*"), (1201, "N", "A"), (1301, "c", "t"), (1401, "G", "N")]
    v = write_vcf(tmp_path / "star.vcf", [("chr2", p, r, a, "PASS", ["0|1:1000"]) for p, r, a in sites])
    s = load_snps(str(v), "chr2", phased=True)
    assert (s.pos + 1).tolist() == [1001, 1301] and s.ref == ["A", "C"] and s.alt == ["G", "T"]
    raw = tmp_path / "calls.vcf"
    raw.write_text(_HETS_HDR + "".join(f"chr2\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\tGT:GQ:DP\t0|1:30:40\n" for p, r, a in sites))
    with pysam.VariantFile(str(raw)) as vf:
        why = [is_candidate(r, 10, 20)[1] for r in vf]
    assert why == ["ok", "not_biallelic_snv", "not_biallelic_snv", "ok", "not_biallelic_snv"]


# --------------------------------------------------------------------------------------------- a worker killed
def _die_on_chr2(ns):
    """Stand-in for cli._run_one in a worker: chr1 is extracted for real; chr2 dies like an OOM-killed process once
    chr1's partial is finished (a broken pool terminates the other workers, so it must not die earlier)."""
    import os
    import time
    if ns.chrom == "chr2":
        done = Path(ns.out) / "partials" / "chr1" / "DONE"
        t0 = time.time()
        while time.time() - t0 < 120:
            try:
                if "layers" in json.loads(done.read_text()):     # DONE is written atomically, last
                    break
            except (OSError, ValueError):
                pass
            time.sleep(0.05)
        os._exit(9)
    return _real_run_one(ns)


def test_run_reports_a_killed_worker(genome, bam, tmp_path, monkeypatch):
    import concurrent.futures
    import multiprocessing
    import scphasedcn.cli as cli
    out = tmp_path / "o"
    assert main(["extract", "--chrom", "chr2"] + _common(genome, bam, out)) == 0      # finished in an earlier run
    real = concurrent.futures.ProcessPoolExecutor

    class ForkPool(real):
        def __init__(self, max_workers=None):
            super().__init__(max_workers=max_workers, mp_context=multiprocessing.get_context("fork"))
    monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", ForkPool)
    monkeypatch.setattr(cli, "_run_one", _die_on_chr2)
    with pytest.raises(SystemExit, match=r"(?s)worker process died.*without a finished partial \(no DONE\): chr2\..*--jobs"):
        main(["run", "--chroms", "chr1,chr2", "--jobs", "2"] + _common(genome, bam, out))
    assert not (out / "partials" / "chr2" / "DONE").exists()                         # not the earlier run's DONE
    # the surviving chromosome was extracted for real: its partial is complete and merges like a separate extraction
    ref = tmp_path / "ref"
    assert main(["extract", "--chrom", "chr1"] + _common(genome, bam, ref)) == 0
    for x in (out, ref):
        assert main(["merge", "--out", str(x), "--chroms", "chr1"]) == 0
    a, b = (ad.read_h5ad(x / "bins_100kb.h5ad") for x in (out, ref))
    assert a.X.sum() > 0 and (a.X != b.X).nnz == 0 and all((a.layers[k] != b.layers[k]).nnz == 0 for k in b.layers)
    pd.testing.assert_frame_equal(a.obs, b.obs)
    pd.testing.assert_frame_equal(_frags(out), _frags(ref))


def test_run_removes_done_of_earlier_runs_before_extracting(genome, bam, tmp_path, monkeypatch):
    """A chromosome queued but never started when the pool breaks must not keep the DONE of an earlier run in the same
    output directory (it would be reported, and merged, as finished)."""
    import concurrent.futures
    from concurrent.futures.process import BrokenProcessPool
    out = tmp_path / "o"
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1"] + _common(genome, bam, out)) == 0   # earlier run

    class BrokenPool:                                  # breaks before any chromosome starts
        def __init__(self, max_workers=None):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def map(self, fn, it):
            raise BrokenProcessPool("a worker died")
    with monkeypatch.context() as mp:
        mp.setattr(concurrent.futures, "ProcessPoolExecutor", BrokenPool)
        with pytest.raises(SystemExit, match=r"without a finished partial \(no DONE\): chr1, chr2\."):
            main(["run", "--chroms", "chr1,chr2", "--jobs", "2"] + _common(genome, bam, out))
    assert not any((out / "partials" / c / "DONE").exists() for c in ("chr1", "chr2"))
    # a later run rewriting one layer keeps the inputs of the others (passed on from the DONE it removed)
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1"] + _common(genome, bam, out)) == 0
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1", "--layers", "fragments"]
                + _common(genome, bam, out)) == 0
    done = json.load(open(out / "partials" / "chr1" / "DONE"))
    assert done["layers"]["depth"]["inputs"]["binsize"] == 100_000


# --------------------------------------------------------------------------------------------- records without SEQ / CIGAR
def _strip_seq(r):
    """SEQ `*` on CELL_A's chr1 record c ([20000, 22000)) and on CELL_B's chr2 R2 (over SNPs 1001 and 1501)."""
    if (r.get_tag("CB"), r.reference_name, r.reference_start) in ((CELL_A, "chr1", 20000), (CELL_B, "chr2", 950)):
        r.query_sequence = None
    return r


def test_records_without_seq_count_in_depth_and_fragments(genome, bam, tmp_path):
    """Depth and fragments need only the CIGAR: a BAM with SEQ stripped gives the same depth layers, QC and fragments
    table as the full BAM (without SNPs, --layers depth,fragments)."""
    b = rewrite_bam(bam, tmp_path / "noseq.bam", _strip_seq)
    for x, bb in (("full", bam), ("noseq", b)):
        o = tmp_path / x
        assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1", "--layers", "depth,fragments"]
                    + _common(genome, bb, o)) == 0
    full, noseq = (ad.read_h5ad(tmp_path / x / "bins_100kb.h5ad") for x in ("full", "noseq"))
    assert full.X.sum() > 0 and (full.X != noseq.X).nnz == 0
    assert all((full.layers[k] != noseq.layers[k]).nnz == 0 for k in full.layers)
    pd.testing.assert_frame_equal(full.obs, noseq.obs)
    pd.testing.assert_frame_equal(_frags(tmp_path / "full"), _frags(tmp_path / "noseq"))
    s = json.load(open(tmp_path / "noseq" / "summary.json"))
    assert s["n_no_seq"] == 2 and s["n_no_cigar"] == 0


def test_records_without_seq_have_no_snp_observation(genome, bam, vcf, tmp_path):
    """In an allele run, a record without SEQ over SNPs is a fragment without SNP observation (hap 0, nsnp 0, ps -1):
    the allele counts lose its observations, nothing else changes."""
    b = rewrite_bam(bam, tmp_path / "noseq.bam", _strip_seq)
    for x, bb in (("full", bam), ("noseq", b)):
        assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1", "--snps", str(vcf), "--phased"]
                    + _common(genome, bb, tmp_path / x)) == 0
    f0, f1 = _frags(tmp_path / "full"), _frags(tmp_path / "noseq")
    r2 = (f1.cell == CELL_B) & (f1.start == 950)
    assert f1.loc[r2, ["hap", "nsnp", "ps"]].values.tolist() == [[0, 0, -1]]
    assert f0.loc[r2, ["hap", "nsnp", "ps"]].values.tolist() == [[2, 2, 1000]]
    pd.testing.assert_frame_equal(f0.loc[~r2], f1.loc[~r2])
    a0, a1 = (ad.read_h5ad(tmp_path / x / "alleles_dna.h5ad") for x in ("full", "noseq"))
    pd.testing.assert_frame_equal(a0.var, a1.var)                       # n_other / n_lowbq unchanged
    d = (a0.X - a1.X).toarray()
    snp = list(a0.var["pos"])
    assert d[1, snp.index(1001)] == 1 and d[1, snp.index(1501)] == 1 and abs(d).sum() == 2   # R2 only (cell B)
    assert ((a0.layers["AD"] - a1.layers["AD"]).toarray()[1, snp.index(1001)]) == 1           # R2 was alt at 1001
    q0, q1 = (pd.read_csv(tmp_path / x / "cells_qc.tsv", sep="\t", index_col=0) for x in ("full", "noseq"))
    assert q1.loc[CELL_B, "n_hapB"] == q0.loc[CELL_B, "n_hapB"] - 1 and q1.loc[CELL_B, "n_pass"] == q0.loc[CELL_B, "n_pass"]
    assert json.load(open(tmp_path / "noseq" / "summary.json"))["n_no_seq"] == 2


def test_records_without_cigar_are_skipped_and_counted(genome, bam, tmp_path):
    """A primary record without CIGAR has no reference span: skipped before the per-record counters, the same way for
    listed cells, "-" and unlisted barcodes, so n_records stays the sum of its categories; counted in n_no_cigar."""
    def edit(r):
        if (r.get_tag("CB"), r.reference_name, r.reference_start) in (
                (CELL_A, "chr1", 30000), ("-", "chr1", 17000), ("GGGGTTTTAAAACCCC", "chr1", 15000)):
            r.cigartuples = None
        return r
    b = rewrite_bam(bam, tmp_path / "nocigar.bam", edit)
    for x, bb in (("full", bam), ("o", b)):
        assert main(["extract", "--chrom", "chr1"] + _common(genome, bb, tmp_path / x)) == 0
        assert main(["merge", "--out", str(tmp_path / x), "--chroms", "chr1"]) == 0
    o = tmp_path / "o"
    p = json.load(open(o / "partials" / "chr1" / "params.json"))
    assert p["n_no_cigar"] == 3 and p["n_no_seq"] == 0
    assert json.load(open(o / "summary.json"))["n_no_cigar"] == 3
    cats = ["n_secondary", "n_supplementary", "n_duplicate", "n_primary", "n_qcfail"]
    for d in (p["dash"], p["unlisted"]):
        assert d["n_records"] == sum(d[k] for k in cats)
    assert p["dash"]["n_records"] == 1 and p["dash"]["n_pass"] == 0                    # d2 (duplicate) only
    assert p["unlisted"]["n_records"] == 1 and p["unlisted"]["n_pass"] == 1
    q = pd.read_csv(o / "cells_qc.tsv", sep="\t", index_col=0)
    assert (q["n_records"] == q[cats].sum(axis=1)).all()
    q0 = pd.read_csv(tmp_path / "full" / "cells_qc.tsv", sep="\t", index_col=0)
    for k in ("n_records", "n_primary", "n_pass"):                                    # one record less: e
        assert q.loc[CELL_A, k] == q0.loc[CELL_A, k] - 1 and q.loc[CELL_B, k] == q0.loc[CELL_B, k]
    assert q.loc[CELL_A, ["mapq0", "mapq1_29", "mapq30_59", "mapq60"]].sum() == q.loc[CELL_A, "n_primary"]
    assert not _frags(o)["start"].isin([30000]).any()


# --------------------------------------------------------------------------------------------- bins tag
def test_tag_follows_the_bin_size(genome, bam, tmp_path):
    from scphasedcn.cli import default_tag
    assert [default_tag(n) for n in (100_000, 50_000, 1_000_000, 500_000, 2_500, 12_345, 1_500_000)] == \
        ["100kb", "50kb", "1Mb", "500kb", "2500bp", "12345bp", "1500kb"]
    o = tmp_path / "run"
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1", "--binsize", "50000"] + _common(genome, bam, o)) == 0
    assert (o / "bins_50kb.h5ad").exists() and not (o / "bins_100kb.h5ad").exists()
    assert ad.read_h5ad(o / "bins_50kb.h5ad").n_vars == 30
    (o / "bins_50kb.h5ad").unlink()
    assert main(["merge", "--out", str(o), "--chroms", "chr1,chr2"]) == 0            # merge: the tag of the partials
    assert (o / "bins_50kb.h5ad").exists()
    assert main(["merge", "--out", str(o), "--chroms", "chr1,chr2", "--tag", "x"]) == 0
    assert sorted(f.name for f in o.glob("bins_*.h5ad")) == ["bins_x.h5ad"]           # one bins file per directory
    # a later run that rewrites the fragments only (default --binsize): the depth layer keeps its 50 kb bins
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1", "--layers", "fragments"] + _common(genome, bam, o)) == 0
    assert sorted(f.name for f in o.glob("bins_*.h5ad")) == ["bins_50kb.h5ad"]


def test_merge_removes_the_bins_file_of_another_bin_size(genome, bam, tmp_path, capsys):
    """Re-extracting at another bin size into the same directory: the R caller prefers bins_100kb.h5ad, so a stale
    one beside the new bins_50kb.h5ad would be paired with the new cells_qc.tsv and fragments.tsv.gz."""
    o = tmp_path / "o"
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1"] + _common(genome, bam, o)) == 0
    assert (o / "bins_100kb.h5ad").exists()
    (o / "bins_other.h5ad").write_text("")
    (o / "keep_bins_100kb.h5ad").write_text("")                                          # not a bins file name
    capsys.readouterr()
    assert main(["run", "--chroms", "chr1,chr2", "--jobs", "1", "--binsize", "50000"] + _common(genome, bam, o)) == 0
    assert sorted(f.name for f in o.glob("*.h5ad")) == ["bins_50kb.h5ad", "keep_bins_100kb.h5ad"]
    err = capsys.readouterr().err
    assert "removed bins_100kb.h5ad" in err and "removed bins_other.h5ad" in err


def test_tag_of_a_depth_layer_without_per_layer_inputs(genome, bam, tmp_path):
    """A depth layer kept from a partial of scphasedcn <= 0.1.1 (no per-layer inputs in DONE): the tag follows the bin
    size of the partial's own bins.tsv, not the top-level `binsize` of the latest extract run (which may have
    rewritten other layers only)."""
    o = tmp_path / "o"
    assert main(["run", "--chroms", "chrM,chr1,chr2", "--jobs", "1", "--binsize", "50000"]
                + _common(genome, bam, o)) == 0
    for c in ("chrM", "chr1", "chr2"):                  # a later fragments-only run at the default bin size, by hand
        P = o / "partials" / c
        done = json.load(open(P / "DONE"))
        done["binsize"] = 100_000; done["layers"]["depth"]["inputs"] = None
        (P / "DONE").write_text(json.dumps(done))
        prm = json.load(open(P / "params.json")); prm["binsize"] = 100_000
        (P / "params.json").write_text(json.dumps(prm))
    (o / "bins_50kb.h5ad").unlink()
    assert main(["merge", "--out", str(o), "--chroms", "chrM,chr1,chr2"]) == 0
    assert sorted(f.name for f in o.glob("bins_*.h5ad")) == ["bins_50kb.h5ad"]
    for c in ("chrM", "chr1", "chr2"):                  # partials of scphasedcn <= 0.1.1: empty DONE
        (o / "partials" / c / "DONE").write_text("")
    (o / "bins_50kb.h5ad").unlink()
    assert main(["merge", "--out", str(o), "--chroms", "chrM,chr1,chr2"]) == 0
    assert sorted(f.name for f in o.glob("bins_*.h5ad")) == ["bins_50kb.h5ad"]


# --------------------------------------------------------------------------------------------- FRiP
def test_frip_counts_clean_peak_fragments(genome, bam, tmp_path):
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t10000\t10500\tp1\nchr1\t900000\t900200\tp2\n")            # p2 covers blacklisted n
    o = tmp_path / "o"
    args = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(o), "--no-gc",
            "--chrom-sizes", str(genome / "chrom.sizes"), "--blacklist", str(genome / "blacklist.bed"),
            "--peaks", str(peaks)]
    assert main(["extract", "--chrom", "chr1"] + args) == 0
    assert main(["merge", "--out", str(o), "--chroms", "chr1"]) == 0
    q = pd.read_csv(o / "cells_qc.tsv", sep="\t", index_col=0)
    a = q.loc[CELL_A]
    assert a["n_peak"] == 2 and a["n_clean"] == 16                    # a and n overlap a peak; n is blacklisted
    assert a["frip"] == pytest.approx(1 / 16)
    assert "n_peak_clean" not in q.columns
    # a partial written before (no n_peak_clean column) keeps the old ratio
    pq = o / "partials" / "chr1" / "qc.tsv"
    pd.read_csv(pq, sep="\t", index_col=0).drop(columns="n_peak_clean").to_csv(pq, sep="\t")
    assert main(["merge", "--out", str(o), "--chroms", "chr1"]) == 0
    assert pd.read_csv(o / "cells_qc.tsv", sep="\t", index_col=0).loc[CELL_A, "frip"] == pytest.approx(2 / 16)


# --------------------------------------------------------------------------------------------- merge defaults
def test_merge_needs_an_output_directory():
    with pytest.raises(SystemExit, match="merge: need --out"):
        main(["merge"])


def test_merge_without_chroms_merges_the_chromosomes_of_the_bam(genome, bam, tmp_path, capsys):
    """Without --chroms, merge takes the set `run` extracts: the main chromosomes of the BAM header (the fixture BAM
    has chr1, chr2, chrM, chrX and no chrY). One without a finished partial (an extract job that died or never ran)
    is missing: an error naming it, or listed in summary.json with --allow-missing."""
    o = tmp_path / "o"
    for ch in ("chr2", "chr1"):
        assert main(["extract", "--chrom", ch] + _common(genome, bam, o)) == 0
    assert main(["extract", "--region", "chr1:0-50000"] + _common(genome, bam, o)) == 0   # never merged
    (o / "partials" / "chrX").mkdir()                                                    # died: no DONE
    with pytest.raises(SystemExit, match=r"no finished partial .* for chrX, chrM; .*--allow-missing"):
        main(["merge", "--out", str(o)])
    assert not (o / "summary.json").exists()
    capsys.readouterr()
    assert main(["merge", "--out", str(o), "--allow-missing"]) == 0
    s = json.load(open(o / "summary.json"))
    assert s["chroms"] == ["chr1", "chr2"] and s["missing"] == ["chrX", "chrM"]
    err = capsys.readouterr().err
    assert "main chromosomes of the BAM" in err and "chr1, chr2, chrX, chrM" in err and "chrY" not in err
    # every chromosome of the BAM extracted: merges without --chroms although the BAM has no chrY
    for ch in ("chrM", "chrX"):
        assert main(["extract", "--chrom", ch] + _common(genome, bam, o)) == 0
    assert main(["merge", "--out", str(o)]) == 0
    s = json.load(open(o / "summary.json"))
    assert s["chroms"] == ["chr1", "chr2", "chrX", "chrM"] and s["missing"] == []
    # an explicit list keeps the missing-partial error and --allow-missing; an unfinished partial outside it is reported
    (o / "partials" / "chrX" / "DONE").unlink()
    capsys.readouterr()
    assert main(["merge", "--out", str(o), "--chroms", "chr1,chr2"]) == 0
    assert "unfinished partials (no DONE; not merged): chrX" in capsys.readouterr().err
    with pytest.raises(SystemExit, match="chrX"):
        main(["merge", "--out", str(o), "--chroms", "chr1,chr2,chrX"])
    assert main(["merge", "--out", str(o), "--chroms", "chr1,chr2,chrX", "--allow-missing"]) == 0
    with pytest.raises(SystemExit, match="no finished partials"):
        main(["merge", "--out", str(tmp_path / "empty")])


def test_merge_without_chroms_and_without_the_bam(genome, bam, tmp_path, capsys):
    """When the recorded BAM cannot be read (moved, or another file at its path), the default set is the main
    chromosomes with a partial directory, finished or not."""
    b = tmp_path / "moved.bam"
    shutil.copy(bam, b); shutil.copy(str(bam) + ".bai", str(b) + ".bai")
    o = tmp_path / "o"
    for ch in ("chr1", "chr2"):
        assert main(["extract", "--chrom", ch] + _common(genome, b, o)) == 0
    (o / "partials" / "chrX").mkdir()                                                    # died: no DONE
    b.unlink()
    with pytest.raises(SystemExit, match=r"no finished partial .* for chrX;"):
        main(["merge", "--out", str(o)])
    capsys.readouterr()
    assert main(["merge", "--out", str(o), "--allow-missing"]) == 0
    assert json.load(open(o / "summary.json"))["missing"] == ["chrX"]
    assert "with a partial directory, finished or not" in capsys.readouterr().err
    b.write_bytes(b"not the BAM the partials were made from")                         # another file at its path
    capsys.readouterr()
    assert main(["merge", "--out", str(o), "--allow-missing"]) == 0
    assert json.load(open(o / "summary.json"))["missing"] == ["chrX"]
    assert "not the file the partials were extracted from" in capsys.readouterr().err
