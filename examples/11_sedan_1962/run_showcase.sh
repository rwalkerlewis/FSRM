#!/bin/bash
# examples/11_sedan_1962/run_showcase.sh
#
# Sedan showcase: cratering shot with Mueller-Murphy + dynamic source.
# Runs the pass-5 dynamic-plastic config alongside the kinematic
# default config so the cavity-formation animation can be regenerated
# from per-snapshot HDF5.
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

echo "=== Sedan 1962 showcase: ${MPI_RANKS} ranks ==="

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
[ -f "${SCRIPT_DIR}/config_dynamic.config" ] && \
    run_tier "config_dynamic.config" "output_dynamic"

cd "${SCRIPT_DIR}"
"${SCRIPT_DIR}/figures/regenerate.sh"

echo ""
echo "=== Showcase complete ==="
