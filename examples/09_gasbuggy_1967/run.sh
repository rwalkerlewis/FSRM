#!/bin/bash
# Example 09: Project Gasbuggy (1967)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/gasbuggy_1967.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 09: Project Gasbuggy (1967) ==="
echo "  29 kt, 1280m depth, Lewis Shale, San Juan Basin NM"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/gasbuggy_1967
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/gasbuggy_1967/ 2>/dev/null || echo "No output files generated."
