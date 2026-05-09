#!/bin/bash
# examples/20_salmon_1964/run_showcase.sh
#
# End-to-end showcase regeneration. Runs the simulation at the LOW,
# MED, and HIGHEST fidelity tiers (all three Salmon configs ship in
# this directory), then regenerates the figure pack.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

source "${REPO_DIR}/scripts/run_with_mpi.sh"

# Showcase ranks default 8 (override with MPI_RANKS).
export MPI_RANKS="${MPI_RANKS:-8}"

echo "=== Salmon 1964 showcase: ${MPI_RANKS} ranks ==="

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
run_tier "config_marshak.config" "output_marshak"
[ -f "${SCRIPT_DIR}/config_highest.config" ] && \
    run_tier "config_highest.config" "output_highest"

cd "${SCRIPT_DIR}"
"${SCRIPT_DIR}/figures/regenerate.sh"

echo ""
echo "=== Showcase complete ==="
echo "Outputs:    ${SCRIPT_DIR}/output*"
echo "Figures:    ${SCRIPT_DIR}/figures/*.png"
