#!/bin/bash
# Example 39: DPRK Sixth Underground Nuclear Test (2017-09-03)
# Punggye-ri, Mt. Mantap, ~250 kt in granite under tuff overburden
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${SCRIPT_DIR}/config.config"
OUT_DIR="${SCRIPT_DIR}/output"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    echo "  mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON && make -j\$(nproc)"
    exit 1
fi

source "${REPO_DIR}/scripts/run_with_mpi.sh"

mkdir -p "${OUT_DIR}"
cd "${SCRIPT_DIR}"

echo "=== Example 39: DPRK 2017 (Punggye-ri, Mt. Mantap) ==="
echo "  ~250 kt, 600 m depth, granite under tuff overburden"
echo "Config:  ${CONFIG}"
echo "Output:  ${OUT_DIR}"
echo "Ranks:   ${MPI_RANKS:-4}"
echo ""

run_with_mpi "${BUILD_DIR}/fsrm" -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh "${OUT_DIR}" 2>/dev/null || echo "No output files generated."
echo ""
echo "Cached IRIS waveforms (if populated):"
echo "  tools/waveform_vv/cache/dprk_2017/"
echo "Refresh with: python scripts/fetch_dprk_2017_waveforms.py"
