#!/bin/bash
# Generic post-processing: build wavefield XDMF + render Sedan-style figure pack.
#
# Usage:
#   postprocess_example.sh <example_dir> <config_path> [<output_dir>] [<title>]
#
# Looks for <output_dir>/solution.h5 (default output_dir = "output") and emits:
#   <output_dir>/wavefield.xdmf
#   <output_dir>/wavefield_aux.h5
#   <output_dir>/figures/{wavefield_panels,surface_panels,seismograms}.png
#   <output_dir>/figures/{wavefield_xz_slice,wavefield_surface}.gif
set -e

EX_DIR="${1:?example dir required}"
CFG="${2:?config path required}"
OUT_DIR="${3:-output}"
TITLE="${4:-$(basename "$EX_DIR")}"

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS="${REPO_DIR}/tools"

cd "$EX_DIR"

if [ ! -f "${OUT_DIR}/solution.h5" ]; then
    echo "postprocess: no ${OUT_DIR}/solution.h5; skipping figure render"
    exit 0
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "postprocess: python3 not found; skipping figure render"
    exit 0
fi
if ! python3 -c 'import h5py, numpy, matplotlib' 2>/dev/null; then
    echo "postprocess: python3 missing h5py/numpy/matplotlib; skipping"
    exit 0
fi

echo ""
echo "=== Post-process: build wavefield XDMF ==="
python3 "${TOOLS}/build_wavefield_xdmf.py" "${OUT_DIR}"

echo ""
echo "=== Post-process: render figure pack ==="
python3 "${TOOLS}/render_wavefield.py" "${OUT_DIR}" "${CFG}" "${TITLE}"

echo ""
echo "Figures in: ${EX_DIR}/${OUT_DIR}/figures"
