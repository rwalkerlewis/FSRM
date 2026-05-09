#!/bin/bash
# Example 24: Cannikin (1971)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/cannikin_1971.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 24: Cannikin (1971) ==="
echo "  ~5 Mt, 1860 m depth, andesite, Amchitka Island AK (deepest US UGT)"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/cannikin_1971
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/cannikin_1971/ 2>/dev/null || echo "No output files generated."
