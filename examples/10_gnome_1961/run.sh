#!/bin/bash
# Example 10: Project Gnome (1961)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/gnome_1961.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 10: Project Gnome (1961) ==="
echo "  3.1 kt, 361m depth, Salado Salt, Carlsbad NM"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/gnome_1961
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/gnome_1961/ 2>/dev/null || echo "No output files generated."
