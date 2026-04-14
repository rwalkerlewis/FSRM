#!/bin/bash
# Run the thermal expansion example
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Running thermal expansion simulation..."
cd ../../build
./fsrm -c ../examples/18_thermal_expansion/config.config \
    -ts_type beuler \
    -pc_type lu \
    -ksp_type preonly

echo "Done."
