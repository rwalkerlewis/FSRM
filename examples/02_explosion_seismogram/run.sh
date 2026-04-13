#!/bin/bash
# Example 02: Explosion Seismogram
# Elastodynamic wave propagation with seismometer recording
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/explosion_seismogram.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first. Run from repository root:"
    echo "  mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON && make -j\$(nproc)"
    exit 1
fi

echo "=== Example 02: Explosion Seismogram ==="
echo "Config: ${CONFIG}"
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/solution.h5 2>/dev/null && echo "  solution.h5 -- HDF5 displacement snapshots"
ls output/seismograms/*.sac 2>/dev/null && echo "  SAC seismogram files for each station"
echo ""
echo "To visualize seismograms:"
echo "  python3 ${REPO_DIR}/scripts/plot_seismograms.py ${BUILD_DIR}/output/seismograms/"
