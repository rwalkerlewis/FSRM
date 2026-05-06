#!/bin/bash
# Example 31: Azgir A-1 (1966)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/azgir_a1_1966.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 31: Azgir A-1 (1966) ==="
echo "  1.1 kt, 165 m depth, halite (Caspian salt dome), Astrakhan RU"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/azgir_a1_1966
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/azgir_a1_1966/ 2>/dev/null || echo "No output files generated."
