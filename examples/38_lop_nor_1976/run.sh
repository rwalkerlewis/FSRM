#!/bin/bash
# Example 38: Lop Nor 1976-10-17
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/lop_nor_1976.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 38: Lop Nor 1976-10-17 ==="
echo "  ~4 kt, ~250 m depth, Beishan-equivalent granite, Lop Nor XJ CN"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/lop_nor_1976
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/lop_nor_1976/ 2>/dev/null || echo "No output files generated."
