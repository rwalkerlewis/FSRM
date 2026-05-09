#!/bin/bash
# Example 20: Salmon (1964) -- Pass-8 Marshak variant.
# Runs the pass-8 dynamic-plastic near-field source path with the
# Marshak grey radiation-diffusion phase enabled, instead of the
# pass-7 Zel'dovich-Raizer end-state approximation.
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/salmon_1964_marshak.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first (see build instructions in repository root)."
    exit 1
fi

echo "=== Example 20 (Pass-8 Marshak): Salmon (1964) ==="
echo "  5.3 kt, 828 m depth, Tatum Salt Dome, MS"
echo "  Radiation phase: MARSHAK_GREY (Pass-8 MED tier)"
echo "  Config: ${CONFIG}"
cd "${BUILD_DIR}"
mkdir -p output/salmon_1964_marshak
./fsrm -c "${CONFIG}"
echo ""
echo "=== Output Files ==="
ls -lh output/salmon_1964_marshak/ 2>/dev/null || echo "No output files generated."
echo ""
echo "Pass-8 V&V: cached observed traces under tools/waveform_vv/cache/Salmon1964/"
echo "Run python tools/waveform_vv/refresh.py --event Salmon1964 to populate the cache."
