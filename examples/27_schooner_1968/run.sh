#!/bin/bash
# Example 27: Schooner (1968)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/schooner_1968.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 27: Schooner (1968) ==="
echo "  30 kt, 111 m depth, welded tuff/basalt, NTS Pahute Mesa (cratering)"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/schooner_1968
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/schooner_1968/ 2>/dev/null || echo "No output files generated."
