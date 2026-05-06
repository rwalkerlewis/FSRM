#!/bin/bash
# Example 32: Pokhran-I (Smiling Buddha, 1974)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/pokhran_i_1974.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 32: Pokhran-I (Smiling Buddha, 1974) ==="
echo "  8 kt, 107 m depth, Vindhyan sandstone, Pokhran Rajasthan IN"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/pokhran_i_1974
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/pokhran_i_1974/ 2>/dev/null || echo "No output files generated."
