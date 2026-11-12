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
