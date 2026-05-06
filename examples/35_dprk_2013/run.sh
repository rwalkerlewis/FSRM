#!/bin/bash
# Example 35: DPRK 2013
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/dprk_2013.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 35: DPRK 2013 ==="
echo "  ~10 kt (mb 5.1), ~600 m depth, Mt. Mantap granite, Punggye-ri NK"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/dprk_2013
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/dprk_2013/ 2>/dev/null || echo "No output files generated."
