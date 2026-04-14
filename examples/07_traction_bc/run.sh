#!/bin/bash
# Example 07: Traction Boundary Condition
# Uniaxial stress with 1 MPa traction on top face

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/traction_bc.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first. Run: cd build && cmake .. && make -j\$(nproc)"
    exit 1
fi

echo "=== Example 07: Traction Boundary Condition ==="
echo "Config: ${CONFIG}"
echo ""
echo "Physics: Uniaxial stress with 1 MPa traction on top face"
echo "Analytical: u_z(100m) = -0.01 m"
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/solution*.h5 2>/dev/null || echo "  (no HDF5 output)"
ls -lh output_SUMMARY.txt 2>/dev/null || echo "  (no summary)"
