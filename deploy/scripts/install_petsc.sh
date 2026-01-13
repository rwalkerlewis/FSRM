#!/bin/bash
# Script to install PETSc on any Linux system
# This can be used on cloud instances or local machines

set -e

echo "=========================================="
echo "PETSc Installation Script"
echo "=========================================="
echo ""

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Configuration
PETSC_VERSION=${PETSC_VERSION:-3.20.0}
INSTALL_DIR=${INSTALL_DIR:-/opt}
PETSC_ARCH=${PETSC_ARCH:-arch-linux-c-opt}

print_info "Installing PETSc $PETSC_VERSION to $INSTALL_DIR"

# Check for required tools
if ! command -v mpicc &> /dev/null; then
    print_warning "MPI not found. Installing OpenMPI..."
    sudo apt-get update
    sudo apt-get install -y libopenmpi-dev openmpi-bin
fi

if ! command -v gfortran &> /dev/null; then
    print_warning "Fortran compiler not found. Installing gfortran..."
    sudo apt-get install -y gfortran
fi

# Install build dependencies
print_info "Installing build dependencies..."
sudo apt-get install -y \
    build-essential \
    python3 \
    wget \
    liblapack-dev \
    libblas-dev

# Download PETSc
cd $INSTALL_DIR
if [ ! -d "petsc-$PETSC_VERSION" ]; then
    print_info "Downloading PETSc $PETSC_VERSION..."
    wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$PETSC_VERSION.tar.gz
    tar -xzf petsc-$PETSC_VERSION.tar.gz
    rm petsc-$PETSC_VERSION.tar.gz
else
    print_info "PETSc directory already exists, skipping download"
fi

cd petsc-$PETSC_VERSION

# Configure PETSc
print_info "Configuring PETSc..."
./configure \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=mpif90 \
    --download-fblaslapack \
    --with-debugging=0 \
    COPTFLAGS='-O3 -march=native' \
    CXXOPTFLAGS='-O3 -march=native' \
    FOPTFLAGS='-O3 -march=native'

# Build PETSc
print_info "Building PETSc (this may take 10-20 minutes)..."
make PETSC_DIR=$INSTALL_DIR/petsc-$PETSC_VERSION PETSC_ARCH=$PETSC_ARCH all

# Test installation
print_info "Testing PETSc installation..."
make PETSC_DIR=$INSTALL_DIR/petsc-$PETSC_VERSION PETSC_ARCH=$PETSC_ARCH check

# Set up environment
print_info "Setting up environment variables..."

ENV_FILE="/etc/profile.d/petsc.sh"
sudo bash -c "cat > $ENV_FILE << EOF
export PETSC_DIR=$INSTALL_DIR/petsc-$PETSC_VERSION
export PETSC_ARCH=$PETSC_ARCH
export PKG_CONFIG_PATH=\\\$PETSC_DIR/\\\$PETSC_ARCH/lib/pkgconfig:\\\$PKG_CONFIG_PATH
EOF"

# Source for current session
export PETSC_DIR=$INSTALL_DIR/petsc-$PETSC_VERSION
export PETSC_ARCH=$PETSC_ARCH
export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig:$PKG_CONFIG_PATH

print_info "PETSc installation complete!"
echo ""
echo "Environment variables set:"
echo "  PETSC_DIR=$PETSC_DIR"
echo "  PETSC_ARCH=$PETSC_ARCH"
echo ""
echo "To use in a new shell, run: source $ENV_FILE"
