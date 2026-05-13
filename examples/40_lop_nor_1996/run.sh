#!/bin/bash
# Example 40: Lop Nor 1996 (Final Chinese Underground Nuclear Test)
# 1996-07-29, Xinjiang, ~5 kt mb 5.0
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${SCRIPT_DIR}/config.config"
OUT_DIR="${SCRIPT_DIR}/output"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

source "${REPO_DIR}/scripts/run_with_mpi.sh"

mkdir -p "${OUT_DIR}"
cd "${SCRIPT_DIR}"

echo "=== Example 40: Lop Nor 1996 (Final Chinese Underground Test) ==="
echo "  ~5 kt, 800 m depth, weathered granite under Tertiary sediment"
echo "Config:  ${CONFIG}"
echo "Output:  ${OUT_DIR}"
echo "Ranks:   ${MPI_RANKS:-4}"
echo ""

run_with_mpi "${BUILD_DIR}/fsrm" -c "${CONFIG}"

echo ""

echo ""
echo "=== Post-process: build wavefield XDMF + render figures ==="
bash "${REPO_DIR}/tools/postprocess_example.sh" "${SCRIPT_DIR}" "${CONFIG}" "output" "40_lop_nor_1996" || true
echo "=== Output Files ==="
ls -lh "${OUT_DIR}" 2>/dev/null || echo "No output files generated."
