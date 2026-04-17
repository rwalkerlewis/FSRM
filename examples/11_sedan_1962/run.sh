#!/bin/bash
# Example 11: Sedan Crater (1962)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/sedan_1962.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 11: Sedan Crater (1962) ==="
echo "  104 kt, 194m depth, alluvium, NTS Yucca Flat"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/sedan_1962
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/sedan_1962/ 2>/dev/null || echo "No output files generated."
