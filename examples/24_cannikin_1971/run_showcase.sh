#!/bin/bash
# examples/24_cannikin_1971/run_showcase.sh
#
# Cannikin 1971 showcase: largest US underground nuclear test (~5 Mt).
# Anchors the high end of the yield-scaling story across five orders
# of magnitude when paired with Sterling 1966.
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

echo "=== Cannikin 1971 showcase: ${MPI_RANKS} ranks ==="

mkdir -p "${SCRIPT_DIR}/output"
cd "${SCRIPT_DIR}"
run_with_mpi "${BUILD_DIR}/fsrm" -c "${SCRIPT_DIR}/config.config"

cd "${SCRIPT_DIR}"
"${SCRIPT_DIR}/figures/regenerate.sh"

echo ""
echo "=== Showcase complete ==="
echo "Outputs: ${SCRIPT_DIR}/output"
echo "Figures: ${SCRIPT_DIR}/figures/*.png"
