#!/bin/bash
# Example 12: Degelen Mountain, Semipalatinsk (Kazakhstan)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/degelen_mountain.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 12: Degelen Mountain, Semipalatinsk ==="
echo "  50 kt, 300m depth, granite, Semipalatinsk Test Site"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/degelen_mountain
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/degelen_mountain/ 2>/dev/null || echo "No output files generated."
