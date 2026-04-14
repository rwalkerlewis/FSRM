#!/bin/bash
# Run single-phase pressure diffusion example
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE="$(cd "$SCRIPT_DIR/../.." && pwd)"

docker run --rm -v "$WORKSPACE:/workspace" -w /workspace/build fsrm-ci:local \
  ./fsrm -c ../examples/14_single_phase_flow/config.config \
  -ts_type beuler -pc_type lu -ksp_type preonly \
  -snes_rtol 1e-8 -snes_atol 1e-10
