#!/bin/bash
# Compare CPU vs GPU performance on the same problem
set -e

CONFIG="${1:-config/examples/punggye_ri_layered.config}"

echo "=== CPU run ==="
docker run --rm -v "$(pwd):/workspace" -w /workspace/build fsrm-ci:local \
  ./fsrm "../${CONFIG}" -log_view 2>&1 | grep "Time (sec):" | tee cpu_timing.txt

echo "=== GPU run ==="
docker run --gpus all --rm -v "$(pwd):/workspace" -w /workspace/build-cuda fsrm-cuda:local \
  ./fsrm "../${CONFIG}" -vec_type cuda -mat_type aijcusparse -log_view 2>&1 | grep "Time (sec):" | tee gpu_timing.txt

echo "=== Comparison ==="
paste cpu_timing.txt gpu_timing.txt
