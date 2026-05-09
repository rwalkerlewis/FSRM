#!/bin/bash
# Example 11: Sedan Crater (1962) -- pass-5 DYNAMIC_PLASTIC variant
# Runs the same Sedan event with [NEAR_FIELD_SOURCE] mode =
# DYNAMIC_PLASTIC, which records a near_field_history.csv alongside
# the SAC seismograms.
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/sedan_1962_dynamic.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 11: Sedan Crater (1962) -- DYNAMIC_PLASTIC ==="
echo "  104 kt, 194m depth, alluvium, NTS Yucca Flat"
echo "  [NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/sedan_1962_dynamic
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/sedan_1962_dynamic/ 2>/dev/null || \
    echo "No output files generated."
echo ""
echo "Near-field history CSV (pass-5):"
echo "  output/sedan_1962_dynamic/near_field_history.csv"
echo ""
echo "ParaView state files:"
echo "  ${SCRIPT_DIR}/paraview/near_field_cavity.pvsm"
echo "  ${SCRIPT_DIR}/paraview/far_field_propagation.pvsm"
echo "  ${SCRIPT_DIR}/paraview/combined.pvsm"
