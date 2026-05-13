#!/bin/bash
# Example 41: Pokhran II Shakti-I (1998-05-11)
# Indian thermonuclear test in granitic gneiss; Marwar craton
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

echo "=== Example 41: Pokhran II Shakti-I (1998) ==="
echo "  ~20 kt central, 210 m depth, granitic gneiss"
echo "Config:  ${CONFIG}"
echo "Output:  ${OUT_DIR}"
echo "Ranks:   ${MPI_RANKS:-4}"
echo ""

run_with_mpi "${BUILD_DIR}/fsrm" -c "${CONFIG}"

echo ""

echo ""
echo "=== Post-process: build wavefield XDMF + render figures ==="
bash "${REPO_DIR}/tools/postprocess_example.sh" "${SCRIPT_DIR}" "${CONFIG}" "output" "41_pokhran_ii_1998" || true
echo "=== Output Files ==="
ls -lh "${OUT_DIR}" 2>/dev/null || echo "No output files generated."
