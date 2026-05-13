#!/bin/bash
# Example 26: Baneberry (1970)
# 
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

echo "=== Example 26: Baneberry (1970) ==="
echo "Config:  ${CONFIG}"
echo "Output:  ${OUT_DIR}"
echo "Ranks:   ${MPI_RANKS:-4}"
echo ""

run_with_mpi "${BUILD_DIR}/fsrm" -c "${CONFIG}"

echo ""

echo ""
echo "=== Post-process: build wavefield XDMF + render figures ==="
bash "${REPO_DIR}/tools/postprocess_example.sh" "${SCRIPT_DIR}" "${CONFIG}" "output" "26_baneberry_1970" || true
echo "=== Output Files ==="
ls -lh "${OUT_DIR}" 2>/dev/null || echo "No output files generated."
