#!/bin/bash
# Script to build FSRM
# Assumes PETSc and dependencies are already installed

set -e

echo "=========================================="
echo "FSRM Build Script"
echo "=========================================="
echo ""

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check for PETSc
if [ -z "$PETSC_DIR" ]; then
    print_error "PETSC_DIR not set. Please install PETSc first."
    exit 1
fi

if [ ! -d "$PETSC_DIR" ]; then
    print_error "PETSC_DIR ($PETSC_DIR) does not exist."
    exit 1
fi

print_info "Using PETSc from: $PETSC_DIR"

# Determine source directory
if [ -d "/workspace" ]; then
    SOURCE_DIR="/workspace"
elif [ -d "$HOME/FSRM" ]; then
    SOURCE_DIR="$HOME/FSRM"
elif [ -d "$(pwd)" ] && [ -f "$(pwd)/CMakeLists.txt" ]; then
    SOURCE_DIR="$(pwd)"
else
    print_error "Cannot find FSRM source directory"
    exit 1
fi

print_info "Source directory: $SOURCE_DIR"

# Build type
BUILD_TYPE=${BUILD_TYPE:-Release}
print_info "Build type: $BUILD_TYPE"

# Number of cores for parallel build
NCORES=$(nproc)
print_info "Building with $NCORES cores"

# Create build directory
BUILD_DIR="$SOURCE_DIR/build"
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# Configure with CMake
print_info "Configuring with CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DENABLE_TESTING=ON \
    -DBUILD_EXAMPLES=ON

# Build
print_info "Building FSRM..."
make -j$NCORES

# Run tests
print_info "Running tests..."
if ctest --output-on-failure; then
    print_info "All tests passed!"
else
    print_warning "Some tests failed. Check the output above."
fi

# Create symlink in user's PATH
print_info "Creating symlink..."
sudo ln -sf $BUILD_DIR/fsrm /usr/local/bin/fsrm

print_info "Build complete!"
echo ""
echo "FSRM installed to: $BUILD_DIR/fsrm"
echo "You can now run: fsrm -help"
echo ""
echo "Examples are in: $BUILD_DIR/examples/"
echo "Configuration files are in: $SOURCE_DIR/config/"
