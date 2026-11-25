#!/usr/bin/env python
"""tiny_extraction.py: write the extraction directory of the Python test suite's synthetic BAM (python/tests/conftest.py,
fixture `extracted`: two cells, chr1 / chr2 / chrM / chrX, phased SNPs on chr2) to a directory, so that other tools
(the R readers, tools/check_extraction_io.R) can read real extractor output.

The fixture code runs through pytest (one test of test_extract.py, which uses the fixture); its output directory is
copied to OUT. A second extraction of the same sample at a 50 kb bin size (`scphasedcn run --binsize 50000`, so its
bins file is bins_50kb.h5ad) is written to OUT_50kb, next to OUT. scphasedcn must be importable (installed, or
PYTHONPATH=python/src).
Usage: python tools/tiny_extraction.py OUT
"""
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent.parent / "python" / "tests"


BINSIZE2 = 50_000


class Grab:
    """pytest plugin: copy the `extracted` fixture's directory to `out` when the test that uses it runs, and extract the
    same sample at BINSIZE2 into `out2` (while the fixture's BAM and reference files exist)."""

    def __init__(self, out, out2):
        self.out, self.out2 = out, out2
        self.path = None
        self.rc2 = None

    @pytest.hookimpl(tryfirst=True)
    def pytest_runtest_call(self, item):
        fa = getattr(item, "funcargs", {})
        if "extracted" in fa:
            self.path = Path(fa["extracted"])
            shutil.copytree(self.path, self.out, dirs_exist_ok=False)
            from conftest import sample_args
            from scphasedcn.cli import main
            self.rc2 = main(["run", "--chroms", "chr1,chr2,chrM,chrX", "--jobs", "1"]
                            + sample_args(fa["genome"], fa["bam"], fa["vcf"], self.out2, binsize=BINSIZE2))


def main(argv):
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    out = Path(argv[0]).resolve()
    out2 = out.with_name(out.name + "_50kb")
    for o in (out, out2):
        if o.exists():
            print(f"tiny_extraction: {o} exists", file=sys.stderr)
            return 2
    grab = Grab(out, out2)
    with tempfile.TemporaryDirectory(prefix="tiny_extraction_") as tmp:
        rc = pytest.main(["-q", "-p", "no:cacheprovider", f"--basetemp={tmp}/pytest",
                          f"{TESTS}/test_extract.py::test_summary_and_params"], plugins=[grab])
    if rc != 0 or grab.path is None or not out.is_dir() or grab.rc2 != 0 or not out2.is_dir():
        print(f"tiny_extraction: the extraction failed (pytest exit {rc}, {BINSIZE2} bp extraction {grab.rc2})",
              file=sys.stderr)
        return 1
    for o in (out, out2):
        print(f"extraction directory: {o}")
        for p in sorted(o.iterdir()):
            print(f"  {p.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
