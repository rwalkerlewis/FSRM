#!/bin/bash
# Example 26: Baneberry (1970)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/baneberry_1970.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 26: Baneberry (1970) ==="
echo "  10 kt, 278 m depth, saturated tuff, NTS Yucca Flat U-8d"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/baneberry_1970
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/baneberry_1970/ 2>/dev/null || echo "No output files generated."
