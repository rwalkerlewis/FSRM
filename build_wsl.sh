#!/bin/bash
# FSRM Build Script for WSL
# Run this script after rebooting and starting WSL: wsl bash build_wsl.sh

set -e  # Exit on error

echo "================================="
echo "FSRM Build Script for WSL"
echo "================================="

# Navigate to project directory
cd /mnt/c/Users/User/Projects/FSRM

echo ""
echo "Step 1: Installing dependencies..."
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    libopenmpi-dev \
    openmpi-bin \
    libhdf5-openmpi-dev \
    libpetsc-real-dev \
    liblapack-dev \
    libblas-dev \
    gfortran \
    gnuplot

echo ""
echo "Step 2: Configuring build environment..."
export PETSC_DIR=/usr/lib/petsc
export PETSC_ARCH=""

echo ""
echo "Step 3: Creating build directory..."
mkdir -p build
cd build

echo ""
echo "Step 4: Running CMake configuration..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_EXAMPLES=ON \
    -DENABLE_TESTING=ON

echo ""
echo "Step 5: Building FSRM (this may take several minutes)..."
make -j$(nproc)

echo ""
echo "Step 6: Running tests..."
ctest --output-on-failure

echo ""
echo "================================="
echo "Build Complete!"
echo "================================="
echo ""
echo "Example binaries are in: build/examples/"
echo ""
echo "To run a simple example:"
echo "  cd examples"
echo "  ./simulator -c ../../config/default.config"
echo ""
echo "To run with MPI (4 processes):"
echo "  mpirun -np 4 ./simulator -c ../../config/hydraulic_fracturing.config"
echo ""
