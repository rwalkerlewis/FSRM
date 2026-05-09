#!/bin/bash
# Example 20: Salmon (1964)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/salmon_1964.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 20: Salmon (1964) ==="
echo "  5.3 kt, 828 m depth, Tatum Salt Dome, MS"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/salmon_1964
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/salmon_1964/ 2>/dev/null || echo "No output files generated."
