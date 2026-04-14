#!/bin/bash
# Example 05: Punggye-ri Nuclear Test
# Layered geology with explosion source and seismograms
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
# Use the quick version (8x8x6 mesh, 1.0s) for fast demonstration.
# For full resolution, use punggye_ri_layered.config (15x15x10, 5.0s).
CONFIG="${REPO_DIR}/config/examples/punggye_ri_layered_quick.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

echo "=== Example 05: Punggye-ri Nuclear Test ==="
echo "Config: ${CONFIG}"
echo "Note: This example may take 10-60 seconds depending on hardware."
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/solution.h5 2>/dev/null && echo "  solution.h5 -- HDF5 displacement snapshots"
ls output/punggye_ri/*.sac 2>/dev/null && echo "  SAC seismogram files"
echo ""
echo "To visualize seismograms:"
echo "  python3 ${REPO_DIR}/scripts/plot_seismograms.py ${BUILD_DIR}/output/punggye_ri/"
