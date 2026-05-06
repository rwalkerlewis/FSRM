#!/bin/bash
# Example 19: Rainier (1957)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/rainier_1957.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 19: Rainier (1957) ==="
echo "  1.7 kt, 274 m depth, bedded tuff, NTS Area 12 (Rainier Mesa)"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/rainier_1957
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/rainier_1957/ 2>/dev/null || echo "No output files generated."
