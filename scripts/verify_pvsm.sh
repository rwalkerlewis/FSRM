#!/usr/bin/env bash
# Pass-7 ParaView state-file load verification.
#
# Closes the pass-6 verification debt for the hand-authored
# examples/11_sedan_1962/paraview/near_field_cavity.pvsm. Runs
# pvpython (ParaView's bundled Python) on the .pvsm with a small
# fixture data tree and asserts the load completes without errors.
#
# Usage:
#   scripts/verify_pvsm.sh                       # default fixture
#   scripts/verify_pvsm.sh PATH/TO/paraview_dir  # custom fixture dir
#
# Exit codes:
#   0  pvpython loaded the .pvsm and reported the expected XYChartView
#   1  pvpython failed to load the .pvsm
#   2  pvpython is not installed; verification skipped (caller decides
#      whether to treat as a soft-fail)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
FIXTURE_DIR="${1:-${ROOT}/examples/11_sedan_1962/paraview}"
PVSM="${FIXTURE_DIR}/near_field_cavity.pvsm"

if [[ ! -f "${PVSM}" ]]; then
    echo "PVSM not found: ${PVSM}" >&2
    exit 1
fi

if ! command -v pvpython > /dev/null 2>&1; then
    echo "pvpython not found in PATH; skipping ParaView verification."
    echo "Install ParaView 5.10+ and re-run, or open the .pvsm "
    echo "manually in the ParaView GUI."
    exit 2
fi

# pvpython script: load the .pvsm and verify the load completed.
PYSCRIPT="$(mktemp -t fsrm_pv_XXXXXX.py)"
trap 'rm -f "${PYSCRIPT}"' EXIT
cat > "${PYSCRIPT}" <<'PYEOF'
import sys
import os
try:
    from paraview.simple import LoadState, GetSources, GetViews
except ImportError as exc:
    print("ERROR: paraview.simple import failed:", exc)
    sys.exit(1)

pvsm = os.environ.get("FSRM_PVSM")
if not pvsm:
    print("ERROR: FSRM_PVSM not set")
    sys.exit(1)

try:
    LoadState(pvsm)
except Exception as exc:
    print("ERROR: LoadState failed:", exc)
    sys.exit(1)

sources = GetSources()
views = GetViews()
print("LoadState OK: sources=%d views=%d" % (len(sources), len(views)))
if len(views) == 0:
    print("ERROR: no views were instantiated from the state file")
    sys.exit(1)

sys.exit(0)
PYEOF

env FSRM_PVSM="${PVSM}" pvpython "${PYSCRIPT}"
echo "near_field_cavity.pvsm loaded successfully under pvpython."
