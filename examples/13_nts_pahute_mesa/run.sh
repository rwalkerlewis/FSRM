#!/bin/bash
# Example 13: NTS Pahute Mesa, Nevada
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/nts_pahute_mesa.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 13: NTS Pahute Mesa, Nevada ==="
echo "  150 kt, 600m depth, zeolitized tuff, Nevada Test Site"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/nts_pahute_mesa
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/nts_pahute_mesa/ 2>/dev/null || echo "No output files generated."
