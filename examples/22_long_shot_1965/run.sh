#!/bin/bash
# Example 22: Long Shot (1965)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/long_shot_1965.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 22: Long Shot (1965) ==="
echo "  80 kt, 716 m depth, andesite/breccia, Amchitka Island AK"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/long_shot_1965
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/long_shot_1965/ 2>/dev/null || echo "No output files generated."
