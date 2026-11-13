"""End-to-end checks on the synthetic BAM (expected values documented in conftest.bam)."""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pysam
import pytest

from conftest import CELL_A, CELL_B


@pytest.fixture(scope="module")
def adata(extracted):
    return ad.read_h5ad(extracted / "bins_100kb.h5ad")


@pytest.fixture(scope="module")
def frags(extracted):
    return pd.read_csv(extracted / "fragments.tsv.gz", sep="\t")


def test_cells_and_bins(adata):
    assert list(adata.obs_names) == sorted([CELL_A, CELL_B])
    assert "-" not in adata.obs_names
    chroms = adata.var["chrom"].value_counts()
    assert chroms["chr1"] == 10 and chroms["chr2"] == 5 and chroms["chrM"] == 1 and chroms["chrX"] == 3
    assert adata.var["width"].iloc[0] == 100_000 and adata.var.loc["chrM:0-16569", "width"] == 16569
    assert 0.3 < adata.var["gc"].iloc[0] < 0.7
    assert adata.var.loc["chr1:900000-1000000", "frac_blacklist"] == pytest.approx(0.1)
    assert adata.var.loc["chr1:700000-800000", "frac_trf"] == pytest.approx(0.001)
    assert adata.var.loc["chr1:500000-600000", "is_centromere"] and not adata.var.loc["chr1:600000-700000", "is_centromere"]
    assert adata.var.loc["chr1:0-100000", "peak_bp"] == 500


def test_record_filter_and_qc(adata):
    q = adata.obs
    a = q.loc[CELL_A]; b = q.loc[CELL_B]
    # CELL_A chr1: 21 records listed (a..u + supplementary), chr2: 3 -> records = 24
    assert a["n_records"] == 24 and a["n_secondary"] == 1 and a["n_supplementary"] == 1 and a["n_duplicate"] == 1
    assert a["n_primary"] == 21 and a["n_mapq_fail"] == 1 and a["n_pass"] == 20
    assert a["mapq0"] == 1 and a["mapq30_59"] == 1 and a["mapq60"] == 19
    assert a["n_blacklist"] == 1 and a["n_clean"] == 19 and a["n_trf"] == 1 and a["n_has_SA"] == 1
    assert a["n_mapq60"] == 19 and a["n_piece_multi"] == 2
    assert a["n_peak"] == 2  # a and b overlap chr1 peak (10000-10500); R1 overlaps chr2 peak
    assert a["n_pass_chrM"] == 0 and b["n_pass_chrM"] == 2 and b["n_pass_chrX"] == 1
    assert b["n_pass"] == 7 and b["frac_chrM"] == pytest.approx(2 / 7)
    assert a["frac_duplicate"] == pytest.approx(1 / 22)
    assert a["n_short_lt300"] == 4  # m (200), o (250), R3 (200), R5 (150)
    assert a["n_hapA"] == 1 and a["n_hapB"] == 0 and b["n_hapB"] == 1 and a["n_snp_obs"] == 3 and b["n_snp_obs"] == 2
    assert a["frip"] == pytest.approx(2 / 19)


def test_depth_layers(adata):
    X = adata.X.toarray(); ia = list(adata.obs_names).index(CELL_A)
    v = list(adata.var_names)
    # start assignment: m=[199950,200150) counted in bin 1 (chr1:100000-200000)
    assert X[ia, v.index("chr1:100000-200000")] == 2   # u and m
    assert X[ia, v.index("chr1:200000-300000")] == 0
    assert X[ia, v.index("chr1:0-100000")] == 12       # a b c d e f g h i j k l
    assert X[ia, v.index("chr1:900000-1000000")] == 0  # blacklisted n dropped
    assert X[ia, v.index("chr1:300000-400000")] == 0   # p mapq0, q dup, r secondary
    assert X[ia, v.index("chr1:400000-500000")] == 1   # t (supplementary piece not counted)
    assert adata.layers["span"][ia, v.index("chr1:400000-500000")] == 800
    assert adata.layers["frag_long"][ia, v.index("chr1:0-100000")] == 11   # d is 900 bp
    assert adata.layers["frag_short"][ia, v.index("chr1:100000-200000")] == 1
    assert adata.layers["hapA"].sum() == 1 and adata.layers["hapB"].sum() == 1 and adata.layers["hap_ambig"].sum() == 0
    fl = adata.obsm["fraglen_hist"][ia]
    assert fl.sum() == 19


def test_fragments_table(extracted, frags):
    assert list(frags.columns) == ["#chrom", "start", "end", "cell", "strand", "mapq", "flags", "len", "piece", "npieces",
                                   "hap", "nsnp", "ps", "read_id"]
    assert set(frags["cell"]) == {CELL_A, CELL_B}
    a = frags[frags.cell == CELL_A]
    assert len(a) == 20 and (a["#chrom"] == "chr1").sum() == 17 and (a["#chrom"] == "chr2").sum() == 3
    n = a[a.start == 900100].iloc[0]; assert n["flags"] & 8          # blacklist flag kept in the table
    t = a[a.start == 400000].iloc[0]; assert t["flags"] & 1 and t["flags"] & 2
    o = a[a.start == 700050].iloc[0]; assert o["flags"] & 16
    u = a[a.start == 120000].iloc[0]; assert u["mapq"] == 45 and not (u["flags"] & 2)
    ij = a[a.start.isin([50000, 50400])]
    assert ij["read_id"].nunique() == 1 and ij["npieces"].tolist() == [2, 2] and ij["piece"].tolist() == [1, 2]
    r1 = a[(a["#chrom"] == "chr2") & (a.start == 900)].iloc[0]
    assert r1["hap"] == 1 and r1["nsnp"] == 3 and r1["ps"] == 1000
    r2 = frags[(frags.cell == CELL_B) & (frags.start == 950)].iloc[0]
    assert r2["hap"] == 2 and r2["nsnp"] == 2
    # within-chromosome sort order and a working tabix index
    for c, g in frags.groupby("#chrom"):
        assert g["start"].is_monotonic_increasing
    tb = pysam.TabixFile(str(extracted / "fragments.tsv.gz"))
    assert len(list(tb.fetch("chr1", 10000, 10001))) == 2     # a and v


def test_alleles(extracted):
    a = ad.read_h5ad(extracted / "alleles_dna.h5ad")
    assert a.n_vars == 3 and list(a.var["pos"]) == [1001, 1501, 2001]
    assert list(a.var["hap_of_alt"]) == [2, 1, 0] and list(a.var["ps"]) == [1000, 1000, -1]
    DP = a.X.toarray(); AD = a.layers["AD"].toarray()
    ia = list(a.obs_names).index(CELL_A); ib = list(a.obs_names).index(CELL_B)
    assert DP[ia].tolist() == [1, 1, 1] and AD[ia].tolist() == [0, 1, 0]
    assert DP[ib].tolist() == [1, 1, 0] and AD[ib].tolist() == [1, 0, 0]
    assert list(a.var["n_other"]) == [1, 0, 1] and list(a.var["n_lowbq"]) == [0, 1, 0]


def test_summary_and_params(extracted):
    s = json.load(open(extracted / "summary.json"))
    assert s["n_cells"] == 2 and s["chroms"] == ["chr1", "chr2", "chrM", "chrX"]
    p = json.load(open(extracted / "partials" / "chr1" / "params.json"))
    assert p["dash"]["n_records"] == 2 and p["dash"]["n_pass"] == 1 and p["unlisted"]["n_pass"] == 2
    assert p["assign"] == "start"


def test_region_and_assign(genome, bam, tmp_path):
    from scphasedcn.cli import main
    from scphasedcn.writers import load_csr_dict
    common = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(tmp_path),
              "--chrom-sizes", str(genome / "chrom.sizes"), "--no-gc", "--no-masks"]
    assert main(["extract", "--region", "chr1:150000-250000", "--assign", "midpoint"] + common) == 0
    rdir = tmp_path / "partials" / "region_chr1_150000-250000"   # a region has its own partial (merge refuses it)
    ia = sorted([CELL_A, CELL_B]).index(CELL_A)
    X = load_csr_dict(rdir / "layers.npz")["counts"].toarray()
    assert X[ia].sum() == 1 and X[ia, 2] == 1           # only m, midpoint 200050 -> bin 2
    assert main(["extract", "--chrom", "chr1", "--assign", "end5"] + common) == 0
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr1"]) == 0
    a = ad.read_h5ad(tmp_path / "bins_100kb.h5ad"); X = a.X.toarray()
    assert X[ia, 2] == 1 and X[ia, 1] == 1              # m 5' end (reverse) 200149 -> bin 2; u in bin 1
    assert main(["extract", "--chrom", "chr1", "--assign", "overlap"] + common) == 0
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr1"]) == 0
    a = ad.read_h5ad(tmp_path / "bins_100kb.h5ad"); X = a.X.toarray()
    assert X[ia, 1] == 2 and X[ia, 2] == 1 and a.layers["span"][ia, 1] == 1050 and a.layers["span"][ia, 2] == 150


def test_alleles_only_rerun(genome, bam, vcf, tmp_path):
    """The allele layer can be recomputed alone (new VCF later) without touching the depth partials."""
    from scphasedcn.cli import main
    common = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(tmp_path),
              "--chrom-sizes", str(genome / "chrom.sizes"), "--no-gc", "--no-masks"]
    assert main(["extract", "--chrom", "chr2"] + common) == 0
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr2"]) == 0
    assert not (tmp_path / "alleles_dna.h5ad").exists()
    assert main(["extract", "--chrom", "chr2", "--layers", "alleles,fragments", "--snps", str(vcf), "--phased"] + common) == 0
    assert main(["merge", "--out", str(tmp_path), "--chroms", "chr2"]) == 0
    a = ad.read_h5ad(tmp_path / "alleles_dna.h5ad"); assert a.n_vars == 3
    b = ad.read_h5ad(tmp_path / "bins_100kb.h5ad"); assert b.layers["hapA"].sum() == 1 and "n_hapA" in b.obs
    f = pd.read_csv(tmp_path / "fragments.tsv.gz", sep="\t"); assert (f["hap"] > 0).sum() == 2


def test_index_check(bam, tmp_path):
    from scphasedcn.cli import main
    from scphasedcn.bam import check_index
    info = check_index(str(bam))
    assert info["exists"] and not info["stale"]
    out = tmp_path / "idx" / "synthetic.bam.bai"
    assert main(["index-check", "--bam", str(bam), "--reindex", str(out), "--force", "--threads", "1"]) == 0
    assert out.exists()
    with pysam.AlignmentFile(str(bam), "rb", index_filename=str(out)) as f:
        assert len(list(f.fetch("chr2"))) == 5


def test_outputs_and_layers(extracted, adata):
    files = sorted(p.name for p in extracted.iterdir() if p.is_file())
    assert files == ["alleles_dna.h5ad", "bins_100kb.h5ad", "cells_qc.tsv", "fragments.tsv.gz", "fragments.tsv.gz.tbi",
                     "summary.json"]
    assert sorted(k for k in adata.layers.keys() if k is not None) == ["frag_long", "frag_short", "hapA", "hapB", "hap_ambig", "span"]
    assert sorted(adata.obsm) == ["fraglen_hist", "mapq_hist"]
    assert adata.X[:, list(adata.var_names).index("chrM:0-16569")].sum() == 2


def test_bad_layers_and_peaks_dir(genome, bam, tmp_path):
    from scphasedcn.cli import main
    common = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--out", str(tmp_path / "o"),
              "--chrom-sizes", str(genome / "chrom.sizes"), "--no-gc", "--no-masks", "--chrom", "chr1"]
    with pytest.raises(SystemExit, match="unknown layer"):
        main(["extract", "--layers", "depth,overlap"] + common)
    with pytest.raises(SystemExit, match="must be a BED file"):
        main(["extract", "--peaks", str(genome)] + common)
    assert not (tmp_path / "o").exists()


