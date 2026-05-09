#!/bin/bash
# Example 37: DPRK 2016b (September)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/dprk_2016b.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 37: DPRK 2016b (September) ==="
echo "  ~25 kt (mb 5.3), ~700 m depth, Mt. Mantap granite, Punggye-ri NK"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/dprk_2016b
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/dprk_2016b/ 2>/dev/null || echo "No output files generated."
