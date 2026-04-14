#!/bin/bash
# Example 06: Gmsh Multi-Material
# Gmsh mesh import with per-region material properties
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/nuclear_twin_gmsh.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first."
    exit 1
fi

if [ ! -f "${REPO_DIR}/meshes/test_two_material.msh" ]; then
    echo "Error: Mesh file meshes/test_two_material.msh not found."
    exit 1
fi

echo "=== Example 06: Gmsh Multi-Material ==="
echo "Config: ${CONFIG}"
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/nuclear_twin_gmsh/solution.h5 2>/dev/null && echo "  solution.h5 -- HDF5 displacement snapshots"
ls -lh output/solution.h5 2>/dev/null && echo "  solution.h5 -- HDF5 displacement snapshots"
echo ""
echo "To visualize: Open the HDF5 output in ParaView"
