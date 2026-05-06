#!/bin/bash
# Example 29: Rio Blanco (1973)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/rio_blanco_1973.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 29: Rio Blanco (1973) ==="
echo "  3 x 33 kt -> single equiv 99 kt, 1903 m depth, Mesaverde, Piceance Basin CO"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/rio_blanco_1973
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/rio_blanco_1973/ 2>/dev/null || echo "No output files generated."
