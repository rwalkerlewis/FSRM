#!/bin/bash
# Example 30: Chagan (1965)
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/chagan_1965.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 30: Chagan ('Atomic Lake', 1965) ==="
echo "  140 kt, 178 m depth, sandstone/siltstone with permafrost, Semipalatinsk KZ"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/chagan_1965
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/chagan_1965/ 2>/dev/null || echo "No output files generated."
