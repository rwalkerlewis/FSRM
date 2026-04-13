#!/bin/bash
# Example 04: Locked Fault
# Cohesive cell insertion with locked-mode constraint
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/locked_fault_compression.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

echo "=== Example 04: Locked Fault ==="
echo "Config: ${CONFIG}"
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/solution.h5 2>/dev/null && echo "  solution.h5 -- HDF5 displacement field"
echo ""
echo "To visualize: Open output/solution.h5 in ParaView"
