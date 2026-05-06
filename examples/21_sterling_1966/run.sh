#!/bin/bash
# Example 21: Sterling (1966)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/sterling_1966.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 21: Sterling (1966) ==="
echo "  0.38 kt, 828 m depth, decoupled in Salmon cavity, Tatum Salt Dome"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/sterling_1966
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/sterling_1966/ 2>/dev/null || echo "No output files generated."
