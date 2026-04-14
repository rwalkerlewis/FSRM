#!/bin/bash
# Run SCEC TPV5 dynamic rupture benchmark
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/../../build"

if [ ! -x "${BUILD_DIR}/fsrm" ]; then
    echo "Error: fsrm executable not found at ${BUILD_DIR}/fsrm"
    echo "Build first: cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j\$(nproc)"
    exit 1
fi

cd "${BUILD_DIR}"
./fsrm -c "${SCRIPT_DIR}/config.config" \
    -ts_type alpha2 \
    -snes_max_it 50 \
    -snes_rtol 1e-6 \
    -snes_atol 1e-6 \
    -pc_type lu \
    -ksp_type preonly \
    -snes_linesearch_type basic \
    -ts_max_snes_failures 100 \
    "$@"
