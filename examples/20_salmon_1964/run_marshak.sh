#!/bin/bash
# Example 20 (Pass-8 Marshak): Salmon (1964) -- Marshak variant
# Runs the dynamic-plastic near-field source with the MARSHAK_GREY
# radiation-diffusion phase enabled (pass-8 MED tier) instead of the
# pass-7 Zel'dovich-Raizer end-state approximation.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${SCRIPT_DIR}/config_marshak.config"
OUT_DIR="${SCRIPT_DIR}/output_marshak"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

source "${REPO_DIR}/scripts/run_with_mpi.sh"

mkdir -p "${OUT_DIR}"
cd "${SCRIPT_DIR}"

echo "=== Example 20 (Pass-8 Marshak): Salmon (1964) ==="
echo "  5.3 kt, 828 m depth, Tatum Salt Dome, MS"
echo "  Radiation phase: MARSHAK_GREY (Pass-8 MED tier)"
echo "Config:  ${CONFIG}"
echo "Output:  ${OUT_DIR}"
echo "Ranks:   ${MPI_RANKS:-4}"
echo ""

run_with_mpi "${BUILD_DIR}/fsrm" -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh "${OUT_DIR}" 2>/dev/null || echo "No output files generated."
echo ""
echo "Pass-8 V&V: cached observed traces under tools/waveform_vv/cache/Salmon1964/"
echo "Run python tools/waveform_vv/refresh.py --event Salmon1964 to populate the cache."
