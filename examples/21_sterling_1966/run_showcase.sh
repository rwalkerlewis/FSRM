#!/bin/bash
# examples/21_sterling_1966/run_showcase.sh
#
# Sterling 1966 showcase: 0.38 kt detonated in the Salmon-generated
# cavity. Anchors the decoupling-factor recovery story alongside
# Salmon 1964 and the low end of the yield-scaling pair with
# Cannikin 1971.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

source "${REPO_DIR}/scripts/run_with_mpi.sh"
export MPI_RANKS="${MPI_RANKS:-8}"

echo "=== Sterling 1966 showcase: ${MPI_RANKS} ranks ==="

mkdir -p "${SCRIPT_DIR}/output"
cd "${SCRIPT_DIR}"
run_with_mpi "${BUILD_DIR}/fsrm" -c "${SCRIPT_DIR}/config.config"

cd "${SCRIPT_DIR}"
"${SCRIPT_DIR}/figures/regenerate.sh"

echo ""
echo "=== Showcase complete ==="
echo "Outputs: ${SCRIPT_DIR}/output"
echo "Figures: ${SCRIPT_DIR}/figures/*.png"
