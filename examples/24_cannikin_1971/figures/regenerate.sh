#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS="${SCRIPT_DIR}/scripts"
cd "${SCRIPT_DIR}"
echo "=== Regenerating Cannikin 1971 showcase figures ==="
for fig in 01_geology_cross_section 02_cavity_formation_radial_profile \
           03_moment_tensor_history 04_synthetic_seismograms_grid \
           05_synthetic_vs_observed 06_fidelity_ladder_comparison; do
    echo "  ${fig}.py ..."
    python3 "${SCRIPTS}/${fig}.py" || \
        echo "    (skipped: missing input data)"
done
echo ""
echo "Figures regenerated under ${SCRIPT_DIR}/"
