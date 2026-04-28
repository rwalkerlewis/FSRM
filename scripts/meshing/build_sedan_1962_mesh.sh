#!/usr/bin/env bash
# =============================================================================
# Build the Sedan 1962 near-field mesh.
#
# Produces meshes/historical/sedan_1962_nearfield.msh in Gmsh MSH version 2.2
# (the format consumed by DMPlexCreateGmshFromFile in the FSRM Simulator).
#
# Usage:
#     bash scripts/meshing/build_sedan_1962_mesh.sh
#
# Requires gmsh (>= 4.10) on PATH.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

GEO_FILE="${SCRIPT_DIR}/sedan_1962_nearfield.geo"
OUT_DIR="${REPO_ROOT}/meshes/historical"
OUT_FILE="${OUT_DIR}/sedan_1962_nearfield.msh"

if ! command -v gmsh >/dev/null 2>&1; then
    echo "error: gmsh not found on PATH" >&2
    echo "       install gmsh (apt-get install gmsh) and retry" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

echo "[build_sedan_1962_mesh] gmsh -3 -format msh2"
echo "    geo : ${GEO_FILE}"
echo "    out : ${OUT_FILE}"

gmsh -3 -format msh2 -o "${OUT_FILE}" "${GEO_FILE}"

echo "[build_sedan_1962_mesh] done."
