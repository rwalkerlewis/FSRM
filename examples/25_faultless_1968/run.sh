#!/bin/bash
# Example 25: Faultless (1968)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/faultless_1968.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 25: Faultless (1968) ==="
echo "  ~1 Mt, 975 m depth, alluvium/Tertiary volcanics, Hot Creek Valley NV"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/faultless_1968
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/faultless_1968/ 2>/dev/null || echo "No output files generated."
