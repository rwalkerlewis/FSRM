#!/bin/bash
# Example 36: DPRK 2016a (January)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/dprk_2016a.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 36: DPRK 2016a (January) ==="
echo "  ~10 kt (mb 5.1), ~600 m depth, Mt. Mantap granite, Punggye-ri NK"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/dprk_2016a
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/dprk_2016a/ 2>/dev/null || echo "No output files generated."
