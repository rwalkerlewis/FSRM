#!/bin/bash
# Run the velocity model example
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Generate velocity model if not present
if [ ! -f velocity_model.bin ]; then
    echo "Generating velocity model..."
    python3 ../../scripts/generate_velocity_model.py velocity_model.bin \
        --nx 10 --ny 10 --nz 10 --xmax 5000 --ymax 5000 --zmax 5000
fi

# Run simulation
echo "Running simulation..."
cd ../../build
./fsrm -c ../examples/17_velocity_model/config.config \
    -ts_type alpha2 \
    -pc_type lu \
    -ksp_type preonly

echo "Done."
