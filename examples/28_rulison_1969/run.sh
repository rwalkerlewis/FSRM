#!/bin/bash
# Example 28: Rulison (1969)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/rulison_1969.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 28: Rulison (1969) ==="
echo "  40 kt, 2568 m depth, Williams Fork (Mesaverde) Sandstone, Piceance Basin CO"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/rulison_1969
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/rulison_1969/ 2>/dev/null || echo "No output files generated."
