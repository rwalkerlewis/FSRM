#!/bin/bash
# Example 34: DPRK 2009
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/dprk_2009.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 34: DPRK 2009 ==="
echo "  ~4 kt (mb 4.5), ~550 m depth, Mt. Mantap granite, Punggye-ri NK"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/dprk_2009
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/dprk_2009/ 2>/dev/null || echo "No output files generated."
