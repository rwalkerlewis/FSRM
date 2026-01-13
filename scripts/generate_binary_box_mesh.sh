#!/usr/bin/env bash
set -euo pipefail

# Generates a small Gmsh *binary* .msh used by config/gmsh_binary_box_example.config.
#
# Requirements:
#   - gmsh must be installed and on PATH
#
# Usage:
#   ./scripts/generate_binary_box_mesh.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

GEO="${ROOT_DIR}/meshes/binary_box/binary_box.geo"
OUT_DIR="${ROOT_DIR}/meshes/binary_box"
OUT_MSH="${OUT_DIR}/binary_box_3d_bin.msh"

if ! command -v gmsh >/dev/null 2>&1; then
  echo "error: gmsh not found on PATH" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

gmsh -3 "${GEO}" -format msh4 -bin 1 -o "${OUT_MSH}"

echo "wrote: ${OUT_MSH}"
