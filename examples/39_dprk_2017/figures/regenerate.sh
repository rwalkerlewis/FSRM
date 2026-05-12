#!/bin/bash
# examples/39_dprk_2017/figures/regenerate.sh
#
# Generate all six showcase figures for DPRK 2017 (Punggye-ri, 6th test).
# Run from the figures/ directory; requires the simulation to have been
# executed first (output/*.sac must exist).
#
# Scripts that need optional data will print WARN/INFO and continue
# rather than aborting, so a partial run still produces useful figures.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS="${SCRIPT_DIR}/scripts"

cd "${SCRIPT_DIR}"

echo "=== Regenerating DPRK 2017 showcase figures ==="
echo "    example: examples/39_dprk_2017/"
echo ""

for fig in 01_geology_cross_section \
           02_cavity_formation_radial_profile \
           03_moment_tensor_history \
           04_synthetic_seismograms_grid \
           05_synthetic_vs_observed \
           06_fidelity_ladder_comparison \
           07_synthetic_seismograms_by_station; do
    echo "  ${fig}.py ..."
    python3 "${SCRIPTS}/${fig}.py" || \
        echo "    (skipped: missing input data -- see WARN above)"
done

echo ""
echo "Figures written under ${SCRIPT_DIR}/"
ls "${SCRIPT_DIR}"/*.png 2>/dev/null || echo "  (no .png files yet)"
