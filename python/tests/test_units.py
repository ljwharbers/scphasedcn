import os
import subprocess
import sys

import numpy as np
import pytest

from scphasedcn.names import parse_name
from scphasedcn.bins import Tiling, merge_intervals, in_intervals, overlap_bp


def test_parse_name():
    cb, rid, p, n, s = parse_name("TTGTTCCTCCCGCAAT_#97e7e0ae-4945-4cf6-bbf4-be87b7e3f675_+1of1")
    assert cb == "TTGTTCCTCCCGCAAT" and (p, n, s) == (1, 1, 1)
    assert rid == int("97e7e0ae49454cf6", 16)
    cb2, rid2, p2, n2, s2 = parse_name("-_#97e7e0ae-4945-4cf6-bbf4-be87b7e3f675_-3of3")
    assert cb2 == "-" and rid2 == rid and (p2, n2, s2) == (3, 3, -1)
    cb3, rid3, p3, n3, s3 = parse_name("plainread")
    assert (p3, n3, s3) == (1, 1, 0)
    assert rid3 == 11265827583416625683   # blake2b-64 of the name, the same in every run


def test_parse_name_fallback_stable_across_processes():
    """The read id of a name outside the convention must not depend on PYTHONHASHSEED (salted hash())."""
    code = "from scphasedcn.names import parse_name; print(parse_name('plainread')[1])"
    out = set()
    for seed in ("1", "2"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        out.add(subprocess.run([sys.executable, "-c", code], env=env, check=True, capture_output=True,
                               text=True).stdout.strip())
    assert len(out) == 1


def test_regular_tiling():
    t = Tiling.regular_grid("chr1", 250_000, 100_000)
    assert len(t) == 3 and t.ends[-1] == 250_000
    assert t.bin_of(np.array([0, 99_999, 100_000, 249_999, 250_000, -1])).tolist() == [0, 0, 1, 2, -1, -1]
    fi, bi = t.bins_overlapping(np.array([99_950, 10]), np.array([200_150, 20]))
    assert fi.tolist() == [0, 0, 0, 1] and bi.tolist() == [0, 1, 2, 0]


def test_bundled_centromeres_cover_chrx():
    from scphasedcn.bins import read_centromere
    from scphasedcn.config import RESOURCES
    p = str(RESOURCES / "centromeres_hg38.tsv")
    assert read_centromere(p, "chrX") == (58605579, 62412542)   # Battenberg names chrX chr23
    assert read_centromere(p, "chr1") == (122026459, 125184587)


def test_intervals():
    s, e = merge_intervals(np.array([10, 0, 5, 50]), np.array([20, 6, 12, 60]))
    assert s.tolist() == [0, 50] and e.tolist() == [20, 60]
    assert in_intervals(np.array([0, 19, 20, 55, 60]), s, e).tolist() == [True, True, False, True, False]
    t = Tiling.regular_grid("chr1", 100, 25)
    assert overlap_bp(t, s, e).tolist() == [20, 0, 10, 0]


def _acc(hoa, ps, vote_ps):
    """AlleleAccumulator on hand-made SNPs (all with ref A, alt G)."""
    from scphasedcn.alleles import SnpSet, AlleleAccumulator
    n = len(hoa)
    snps = SnpSet("chr1", np.arange(n, dtype=np.int64) * 100, ["A"] * n, ["G"] * n, np.asarray(hoa, dtype=np.int8),
                  np.asarray(ps, dtype=np.int64), ["."] * n)
    return AlleleAccumulator(snps, vote_ps=vote_ps)


def test_vote_ps_junction():
    """A fragment across a phase-block junction: 1 SNP votes A in PS 50, 2 vote B in PS 4000."""
    hoa, ps = [1, 1, 1, 0], [50, 4000, 4000, -1]          # alt on hapA; the last SNP unphased
    obs = [(0, 1), (1, 0), (2, 0), (3, 1)]                 # alt (A), ref (B), ref (B), unphased alt
    assert _acc(hoa, ps, "first").vote(obs) == (2, 4, 50)        # B filed under the first phase set
    assert _acc(hoa, ps, "majority").vote(obs) == (2, 4, 4000)   # only PS 4000 votes
    # the majority phase set decides the haplotype: PS 50 votes A, A, B (3 obs), PS 4000 votes B, B
    hoa2, ps2 = [1] * 5, [50, 50, 50, 4000, 4000]
    obs2 = [(0, 1), (1, 1), (2, 0), (3, 0), (4, 0)]
    assert _acc(hoa2, ps2, "first").vote(obs2) == (2, 5, 50)
    assert _acc(hoa2, ps2, "majority").vote(obs2) == (1, 5, 50)


def test_vote_ps_single_phase_set_identical():
    hoa, ps = [1, 2, 1, 0, 2], [7, 7, 7, -1, 7]
    for obs in ([(0, 1), (1, 1), (2, 0), (3, 0)], [(0, 1), (1, 0)], [(0, 1), (1, 1)], [(3, 1)], [(0, 2), (1, 1)], []):
        assert _acc(hoa, ps, "first").vote(obs) == _acc(hoa, ps, "majority").vote(obs)


def test_vote_ps_tie_between_phase_sets():
    hoa, ps = [1, 1, 1, 1], [50, 50, 4000, 4000]
    obs = [(0, 1), (1, 1), (2, 1), (3, 0)]                 # PS 50: A, A; PS 4000: A, B
    assert _acc(hoa, ps, "first").vote(obs) == (1, 4, 50)
    assert _acc(hoa, ps, "majority").vote(obs) == (3, 4, 50)   # 2 vs 2 phased observations: ambiguous
    import pytest
    with pytest.raises(ValueError, match="vote_ps"):
        _acc(hoa, ps, "last")


def test_vote_ps_cli(genome, bam, vcf, tmp_path):
    """--vote-ps is recorded in params.json; on single-phase-set reads both settings write the same fragments."""
    import json
    from scphasedcn.cli import main
    common = ["--bam", str(bam), "--barcodes", str(genome / "whitelist.tsv"), "--chrom-sizes", str(genome / "chrom.sizes"),
              "--no-gc", "--no-masks", "--chrom", "chr2", "--snps", str(vcf), "--phased"]
    assert main(["extract", "--out", str(tmp_path / "f")] + common) == 0
    assert main(["extract", "--out", str(tmp_path / "m"), "--vote-ps", "majority"] + common) == 0
    pf, pm = (tmp_path / k / "partials" / "chr2" for k in ("f", "m"))
    assert json.load(open(pf / "params.json"))["vote_ps"] == "first"
    assert json.load(open(pm / "params.json"))["vote_ps"] == "majority"
    assert (pf / "fragments.tsv.gz").read_bytes() == (pm / "fragments.tsv.gz").read_bytes()
    from scphasedcn.writers import load_csr_dict
    lf, lm = load_csr_dict(pf / "alleles.npz"), load_csr_dict(pm / "alleles.npz")   # npz zips carry timestamps
    assert sorted(lf) == sorted(lm) and all((lf[k] != lm[k]).nnz == 0 for k in lf)


def test_vote_ps_cli_across_a_phase_set_junction(genome, tmp_path):
    """End to end: fragments across a phase-set junction get another haplotype (read 1) or phase set (read 2) with
    --vote-ps majority than with first. SNPs C>G, alt on hapA: PS 10000 at 10001, 10101, 10201; PS 20000 at 10301,
    10401. Read 1 (10000-10500): alt, alt, ref | ref, ref -> first: B (2 A vs 3 B) under PS 10000; majority: PS 10000
    (3 observations) votes A. Read 2 (10150-10500): ref | ref, ref -> first: B under PS 10000; majority: B under 20000."""
    import json
    import pandas as pd
    from conftest import CELL_A, name, write_bam, write_vcf
    from scphasedcn.cli import main
    from scphasedcn.writers import load_csr_dict
    snps = [(10001, 10000), (10101, 10000), (10201, 10000), (10301, 20000), (10401, 20000)]
    v = write_vcf(tmp_path / "junction.vcf", [("chr2", p, "C", "G", "PASS", [f"1|0:{ps}"]) for p, ps in snps])
    r1 = {p - 1 - 10000: b for (p, _), b in zip(snps, "GGCCC")}
    r2 = {p - 1 - 10150: b for (p, _), b in zip(snps[2:], "CCC")}
    b = write_bam(tmp_path / "junction.bam", [((name(CELL_A, "j1"), "chr2", 10000, "500M", "+"), {"bases": r1}),
                                              ((name(CELL_A, "j2"), "chr2", 10150, "350M", "+"), {"bases": r2})])
    common = ["--bam", str(b), "--barcodes", str(genome / "whitelist.tsv"), "--chrom-sizes", str(genome / "chrom.sizes"),
              "--no-gc", "--no-masks", "--chrom", "chr2", "--snps", str(v), "--phased"]
    out = {}
    for mode in ("first", "majority"):
        assert main(["extract", "--out", str(tmp_path / mode), "--vote-ps", mode] + common) == 0
        p = tmp_path / mode / "partials" / "chr2"
        f = pd.read_csv(p / "fragments.tsv.gz", sep="\t", header=None)
        out[mode] = (f[1].tolist(), f[10].tolist(), f[11].tolist(), f[12].tolist(), load_csr_dict(p / "alleles.npz"))
        assert json.load(open(p / "DONE"))["layers"]["fragments"]["inputs"]["vote_ps"] == mode
    assert out["first"][:4] == ([10000, 10150], [2, 2], [5, 3], [10000, 10000])
    assert out["majority"][:4] == ([10000, 10150], [1, 2], [5, 3], [10000, 20000])
    hf, hm = out["first"][4], out["majority"][4]
    assert (hf["hapA"].sum(), hf["hapB"].sum()) == (0, 2) and (hm["hapA"].sum(), hm["hapB"].sum()) == (1, 1)
    assert (hf["DP"] != hm["DP"]).nnz == 0 and (hf["AD"] != hm["AD"]).nnz == 0     # allele counts do not depend on it


def _accumulate_dense(cell, bin_idx, n_cells, n_bins, weights=None):
    """The dense reference implementation of depth.accumulate (before it was made sparse)."""
    import scipy.sparse as sp
    ok = bin_idx >= 0
    key = cell[ok].astype(np.int64) * n_bins + bin_idx[ok]
    w = None if weights is None else weights[ok]
    dense = np.bincount(key, weights=w, minlength=n_cells * n_bins).reshape(n_cells, n_bins)
    m = sp.csr_matrix(dense.astype(np.int64 if w is None else np.float64))
    m.eliminate_zeros()
    return m


def _same_csr(a, b):
    return (a.shape == b.shape and a.dtype == b.dtype and a.indices.dtype == b.indices.dtype
            and a.indptr.dtype == b.indptr.dtype and np.array_equal(a.data, b.data)
            and np.array_equal(a.indices, b.indices) and np.array_equal(a.indptr, b.indptr) and a.has_sorted_indices)


def test_accumulate_matches_the_dense_reference():
    from scphasedcn.depth import accumulate
    rng = np.random.default_rng(3)
    for n in (0, 1, 50, 500):
        cell = rng.integers(0, 7, n).astype(np.int32)
        bi = rng.integers(-1, 13, n).astype(np.int64)
        w = rng.integers(0, 3, n).astype(np.float64)         # zeros: a bin whose weights sum to 0 is not stored
        for weights in (None, w):
            assert _same_csr(accumulate(cell, bi, 7, 13, weights=weights), _accumulate_dense(cell, bi, 7, 13, weights))


def test_accumulate_is_sparse():
    """No dense cells x bins array: 10^7 x 10^7 would need 800 TB."""
    from scphasedcn.depth import accumulate
    m = accumulate(np.array([0, 9_999_999, 0], np.int32), np.array([5, 9_999_998, 5]), 10_000_000, 10_000_000)
    assert m.shape == (10_000_000, 10_000_000) and m.nnz == 2 and m[0, 5] == 2 and m.dtype == np.int64


def _depth_layers_old(fr, tiling, n_cells, assign="start", use=None):
    """Frozen copy of depth.depth_layers before the layers were built sparse and could be selected (dense
    accumulator); the reference the current depth_layers is compared with."""
    from scphasedcn.bam import F_BLACKLIST
    from scphasedcn.depth import assigned_position
    accumulate = _accumulate_dense
    if use is None:
        use = (fr.flags & F_BLACKLIST) == 0
    nb = len(tiling)
    cell = fr.cell[use]; start = fr.start[use]; end = fr.end[use]; L = end - start
    out = {}
    if assign == "overlap":
        fi, bi = tiling.bins_overlapping(start, end)
        c = cell[fi]
        bases = np.minimum(end[fi], tiling.ends[bi]) - np.maximum(start[fi], tiling.starts[bi])
        out["counts"] = accumulate(c, bi, n_cells, nb)
        out["span"] = accumulate(c, bi, n_cells, nb, weights=bases.astype(np.float64))
        out["frag_short"] = accumulate(c[L[fi] < 300], bi[L[fi] < 300], n_cells, nb)
        out["frag_long"] = accumulate(c[L[fi] >= 1000], bi[L[fi] >= 1000], n_cells, nb)
        hap = fr.hap[use][fi]; bi_h = bi
    else:
        pos = assigned_position(fr, assign)[use]
        bi = tiling.bin_of(pos)
        out["counts"] = accumulate(cell, bi, n_cells, nb)
        out["span"] = accumulate(cell, bi, n_cells, nb, weights=L.astype(np.float64))
        out["frag_short"] = accumulate(cell[L < 300], bi[L < 300], n_cells, nb)
        out["frag_long"] = accumulate(cell[L >= 1000], bi[L >= 1000], n_cells, nb)
        hap = fr.hap[use]; c = cell; bi_h = bi
    for code, name in ((1, "hapA"), (2, "hapB"), (3, "hap_ambig")):
        sel = hap == code
        out[name] = accumulate(c[sel], bi_h[sel], n_cells, nb)
    for k in out:
        out[k] = out[k].astype(np.float64 if k == "span" else np.int32)
    return out


def test_depth_layers_match_the_dense_reference(genome, bam, vcf):
    """The current depth_layers equals the frozen pre-sparse one, for every assignment and bin size, whole and for a
    subset of layers (also with an explicit `use`)."""
    import scphasedcn.depth as D
    from scphasedcn.alleles import AlleleAccumulator, load_snps
    from scphasedcn.bam import scan_chrom, flag_fragments
    from scphasedcn.bins import Tiling, read_bed_intervals
    from conftest import CELL_A, CELL_B
    cells = [CELL_A, CELL_B]
    for ch, n in (("chr1", 1_000_000), ("chr2", 500_000)):
        acc = AlleleAccumulator(load_snps(str(vcf), ch, phased=True))
        fr, _ = scan_chrom(str(bam), ch, cells, alleles=acc)
        flag_fragments(fr, read_bed_intervals(str(genome / "blacklist.bed"), ch), read_bed_intervals(None, ch), None,
                       read_bed_intervals(None, ch))
        assert len(fr) > 0 and ((fr.flags & 8) > 0).any() == (ch == "chr1")     # chr1 has a blacklisted fragment
        for binsize in (100_000, 7_000):
            t = Tiling.regular_grid(ch, n, binsize)
            for assign in D.ASSIGN_MODES:
                for use in (None, fr.mapq >= 60):
                    new = D.depth_layers(fr, t, 2, assign=assign, use=use)
                    ref = _depth_layers_old(fr, t, 2, assign=assign, use=use)
                    assert sorted(new) == sorted(ref) == sorted(D.DEPTH_LAYERS + D.HAP_LAYERS)
                    assert all(_same_csr(new[k], ref[k]) for k in ref)
                    for sub in (D.HAP_LAYERS, D.DEPTH_LAYERS, ["span", "hapB"]):
                        part = D.depth_layers(fr, t, 2, assign=assign, use=use, layers=sub)
                        assert sorted(part) == sorted(sub) and all(_same_csr(part[k], ref[k]) for k in part)


def test_translate_barcodes(tmp_path):
    from scphasedcn.config import translate_barcodes
    m = tmp_path / "map.tsv"
    m.write_text("index\tatac_barcode\trna_barcode\n1\tAAAA\tCCCC\n2\tGGGG\tTTTT\n")
    assert translate_barcodes(["CCCC", "TTTT", "NNNN"], str(m)) == (["AAAA", "GGGG"], "rna->atac")
    assert translate_barcodes(["AAAA", "GGGG"], str(m)) == (["AAAA", "GGGG"], "atac")
    assert translate_barcodes(["AAAA"], None) == (["AAAA"], "atac")


def test_sample_config_chroms_default_and_set():
    from scphasedcn.config import MAIN_CHROMS, sample_from_config
    cfg = {"samples": {"s1": {"bam": "x.bam"}, "s2": {"bam": "x.bam", "chroms": ["chr1", "chr2"]}}}
    s1, s2 = sample_from_config(cfg, "s1"), sample_from_config(cfg, "s2")
    assert s1.chroms == MAIN_CHROMS and s1.chroms is not MAIN_CHROMS and not s1.chroms_set
    assert s2.chroms == ["chr1", "chr2"] and s2.chroms_set


def test_translate_barcodes_reports_and_checks_the_match(tmp_path, capsys):
    from scphasedcn.config import translate_barcodes
    m = tmp_path / "map.tsv"
    m.write_text("atac_barcode\trna_barcode\n" + "".join(f"A{i}-1\tR{i}-1\n" for i in range(20)))
    # no whitelist barcode in either column: a wrong map
    with pytest.raises(SystemExit, match="none of the 2 whitelist barcodes"):
        translate_barcodes(["X1", "X2"], str(m))
    # a tie: warned, used as is (ATAC space) as before
    capsys.readouterr()
    assert translate_barcodes(["A1", "R2"], str(m)) == (["A1", "R2"], "atac")
    err = capsys.readouterr().err
    assert "warning" in err and "tie" in err
    # RNA whitelist, all in the map: translated, both counts logged, no warning
    rna = [f"R{i}" for i in range(20)]
    assert translate_barcodes(rna, str(m)) == (sorted(f"A{i}" for i in range(20)), "rna->atac")
    err = capsys.readouterr().err
    assert "0 of 20" in err and "20 of 20" in err and "warning" not in err
    # RNA whitelist, 15 of 20 in the map: 5 dropped, warned
    out, space = translate_barcodes(rna[:15] + [f"Z{i}" for i in range(5)], str(m))
    err = capsys.readouterr().err
    assert space == "rna->atac" and len(out) == 15
    assert "5 whitelist barcodes not in the map are dropped" in err and "warning" in err and "75%" in err
