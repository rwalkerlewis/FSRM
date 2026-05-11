#!/bin/bash
# examples/24_cannikin_1971/run_radial.sh
#
# Cannikin 1971 -- DYNAMIC_PLASTIC + RADIAL_LAGRANGIAN variant.
# Runs the 1D radial Lagrangian shock solver at setup time (writes
# output/near_field_history.csv and output/near_field_profile.h5),
# then drives the FEM elastodynamic solve from the recorded Mdot(t)
# history with 4 MPI ranks.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${SCRIPT_DIR}/config_radial.config"
OUT_DIR="${SCRIPT_DIR}/output"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    echo "  mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON && make -j\$(nproc)"
    exit 1
fi

source "${REPO_DIR}/scripts/run_with_mpi.sh"
export MPI_RANKS="${MPI_RANKS:-4}"

echo "=== Example 24: Cannikin (1971) -- RADIAL_LAGRANGIAN ==="
echo "Config:  ${CONFIG}"
echo "Output:  ${OUT_DIR}"
echo "Ranks:   ${MPI_RANKS}"
echo ""

mkdir -p "${OUT_DIR}"
cd "${SCRIPT_DIR}"

run_with_mpi "${BUILD_DIR}/fsrm" -c "${CONFIG}" \
    -ksp_type gmres -pc_type bjacobi -sub_pc_type lu \
    -ksp_max_it 500 -ksp_gmres_restart 100 -ksp_rtol 1e-4 \
    -ksp_atol 1e-50 -ksp_pc_side right \
    -ts_max_snes_failures unlimited -ts_adapt_type none

# Ensure seismograms symlink is in place for figure scripts
if [ ! -e "${OUT_DIR}/seismograms" ]; then
    ln -sfn "${OUT_DIR}/cannikin_1971" "${OUT_DIR}/seismograms"
    echo "Linked output/seismograms -> output/cannikin_1971"
fi

echo ""
echo "=== Near-field outputs ==="
ls -lh "${OUT_DIR}/cannikin_1971/near_field_history.csv" 2>/dev/null || echo "  near_field_history.csv not found"
ls -lh "${OUT_DIR}/cannikin_1971/near_field_profile.h5"  2>/dev/null || echo "  near_field_profile.h5 not found"

echo ""
echo "=== Output Files ==="
ls -lh "${OUT_DIR}" 2>/dev/null || echo "No output files generated."
