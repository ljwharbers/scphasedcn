#!/usr/bin/env bash
# check_versions.sh: the R and Python versions agree with the git ref being built.
#   On a release tag vX.Y.Z: DESCRIPTION Version and scphasedcn.__version__ must both be X.Y.Z.
#   Anywhere else: both must carry a development marker (R: four-part X.Y.Z.9000; Python: X.Y.Z.devN), so that output
#   written by an unreleased build never claims a released version.
# Usage: tools/check_versions.sh [ref]   (default: $GITHUB_REF; e.g. refs/tags/v0.1.2, refs/heads/main, v0.1.2)
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
REF=${1:-${GITHUB_REF:-}}
R_VER=$(sed -n 's/^Version:[[:space:]]*//p' "$ROOT/DESCRIPTION" | tr -d '[:space:]')
PY_VER=$(sed -n 's/^__version__[[:space:]]*=[[:space:]]*["'\'']\([^"'\'']*\)["'\''].*/\1/p' "$ROOT/python/src/scphasedcn/__init__.py")
[ -n "$R_VER" ] || { echo "check_versions: no Version in DESCRIPTION" >&2; exit 1; }
[ -n "$PY_VER" ] || { echo "check_versions: no __version__ in python/src/scphasedcn/__init__.py" >&2; exit 1; }
echo "R package: $R_VER; Python package: $PY_VER; ref: ${REF:-<none>}"
TAG=${REF#refs/tags/}
fail=0
if [[ "$TAG" =~ ^v([0-9]+\.[0-9]+\.[0-9]+)$ ]]; then
  want=${BASH_REMATCH[1]}
  [ "$R_VER" = "$want" ] || { echo "check_versions: tag $TAG but DESCRIPTION Version is $R_VER" >&2; fail=1; }
  [ "$PY_VER" = "$want" ] || { echo "check_versions: tag $TAG but scphasedcn.__version__ is $PY_VER" >&2; fail=1; }
else
  [[ "$R_VER" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.9[0-9]{3}$ ]] ||
    { echo "check_versions: not a release tag, but the R version $R_VER has no development part (.9000)" >&2; fail=1; }
  [[ "$PY_VER" =~ ^[0-9]+\.[0-9]+\.[0-9]+\.dev[0-9]+$ ]] ||
    { echo "check_versions: not a release tag, but the Python version $PY_VER has no development part (.devN)" >&2; fail=1; }
fi
[ "$fail" = 0 ] && echo "versions ok"
exit "$fail"
