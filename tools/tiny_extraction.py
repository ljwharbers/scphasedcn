#!/usr/bin/env python
"""tiny_extraction.py: write the extraction directory of the Python test suite's synthetic BAM (python/tests/conftest.py,
fixture `extracted`: two cells, chr1 / chr2 / chrM / chrX, phased SNPs on chr2) to a directory, so that other tools
(the R readers, tools/check_extraction_io.R) can read real extractor output.

The fixture code runs through pytest (one test of test_extract.py, which uses the fixture); its output directory is
copied to OUT. scphasedcn must be importable (installed, or PYTHONPATH=python/src).
Usage: python tools/tiny_extraction.py OUT
"""
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent.parent / "python" / "tests"


class Grab:
    """pytest plugin: copy the `extracted` fixture's directory to `out` when the test that uses it runs."""

    def __init__(self, out):
        self.out = out
        self.path = None

    @pytest.hookimpl(tryfirst=True)
    def pytest_runtest_call(self, item):
        if "extracted" in getattr(item, "funcargs", {}):
            self.path = Path(item.funcargs["extracted"])
            shutil.copytree(self.path, self.out, dirs_exist_ok=False)


def main(argv):
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    out = Path(argv[0]).resolve()
    if out.exists():
        print(f"tiny_extraction: {out} exists", file=sys.stderr)
        return 2
    grab = Grab(out)
    with tempfile.TemporaryDirectory(prefix="tiny_extraction_") as tmp:
        rc = pytest.main(["-q", "-p", "no:cacheprovider", f"--basetemp={tmp}/pytest",
                          f"{TESTS}/test_extract.py::test_summary_and_params"], plugins=[grab])
    if rc != 0 or grab.path is None or not out.is_dir():
        print(f"tiny_extraction: the extraction failed (pytest exit {rc})", file=sys.stderr)
        return 1
    print(f"extraction directory: {out}")
    for p in sorted(out.iterdir()):
        print(f"  {p.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
