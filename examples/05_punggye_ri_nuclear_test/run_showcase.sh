#!/bin/bash
# examples/05_punggye_ri_nuclear_test/run_showcase.sh
#
# Punggye-ri showcase: motivates axis-2 topography. Anchored on the
# layered velocity model from this example; the production-resolution
# DPRK 2017 event lives at examples/39_dprk_2017/.
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

echo "=== Punggye-ri showcase: ${MPI_RANKS} ranks ==="

run_tier() {
    local cfg="$1"
    local outdir="$2"
    echo ""
    echo "--- ${cfg} -> ${outdir} ---"
    mkdir -p "${SCRIPT_DIR}/${outdir}"
    cd "${SCRIPT_DIR}"
    OUTDIR="${outdir}" run_with_mpi "${BUILD_DIR}/fsrm" \
        -c "${SCRIPT_DIR}/${cfg}"
}

run_tier "config.config" "output"
[ -f "${SCRIPT_DIR}/config_full.config" ] && \
    run_tier "config_full.config" "output_full"

cd "${SCRIPT_DIR}"
"${SCRIPT_DIR}/figures/regenerate.sh"

echo ""
echo "=== Showcase complete ==="
