#!/bin/bash
# Run FSRM with PETSc CUDA acceleration
# Requires: Dockerfile.cuda image, NVIDIA GPU, nvidia-container-toolkit

set -e

CONFIG="${1:?Usage: $0 <config_file>}"

echo "Building CUDA image..."
docker build -f Dockerfile.cuda -t fsrm-cuda:local .

echo "Building FSRM (CUDA)..."
docker run --gpus all --rm -v "$(pwd):/workspace" -w /workspace fsrm-cuda:local bash -c \
  'mkdir -p build-cuda && cd build-cuda && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF && make -j$(nproc)'

echo "Running with GPU acceleration..."
docker run --gpus all --rm -v "$(pwd):/workspace" -w /workspace/build-cuda fsrm-cuda:local \
  ./fsrm "../${CONFIG}" -vec_type cuda -mat_type aijcusparse \
  -log_view 2>&1 | tee gpu_run.log

echo "GPU execution complete. See gpu_run.log for PETSc performance summary."
