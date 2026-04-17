#!/bin/bash
# Example 08: Time-dependent Prescribed Slip
# Strike-slip fault with linear time ramp

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${REPO_DIR}/build"
CONFIG="${REPO_DIR}/config/examples/time_dependent_slip.config"

if [ ! -f "${BUILD_DIR}/fsrm" ]; then
    echo "Error: Build FSRM first. Run: cd build && cmake .. && make -j\$(nproc)"
    exit 1
fi

echo "=== Example 08: Time-dependent Prescribed Slip ==="
echo "Config: ${CONFIG}"
echo ""
echo "Physics: Strike-slip fault with linear ramp from t=0.1 to t=0.6"
echo "Maximum slip: 1.0 m right-lateral"
echo ""

cd "${BUILD_DIR}"
./fsrm -c "${CONFIG}"

echo ""
echo "=== Output Files ==="
ls -lh output/solution*.h5 2>/dev/null || echo "  (no HDF5 output)"
ls -lh output_SUMMARY.txt 2>/dev/null || echo "  (no summary)"
