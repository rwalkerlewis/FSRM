#!/bin/bash
# Example 01: Uniaxial Compression
# Quasi-static linear elastostatics
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/uniaxial_compression.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first. Run from repository root:"
    echo "  mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON && make -j\$(nproc)"
    exit 1
fi

echo "=== Example 01: Uniaxial Compression ==="
echo "Config: ${CONFIG}"
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/solution.h5 2>/dev/null && echo "  solution.h5 -- HDF5 displacement field"
ls -lh output_SUMMARY.txt 2>/dev/null && echo "  output_SUMMARY.txt -- simulation summary"
echo ""
echo "To visualize: Open output/solution.h5 in ParaView or use scripts/plot_wavefield.py"
