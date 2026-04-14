#!/bin/bash
# Run viscoelastic attenuation example
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE="$(cd "$SCRIPT_DIR/../.." && pwd)"

docker run --rm -v "$WORKSPACE:/workspace" -w /workspace/build fsrm-ci:local \
  ./fsrm -c ../examples/15_viscoelastic_attenuation/config.config \
  -ts_type alpha2 -pc_type lu -ksp_type preonly \
  -snes_rtol 1e-6 -snes_atol 1e-8
