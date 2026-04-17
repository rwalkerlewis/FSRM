#!/bin/bash
# ============================================================================
# FSRM Development Container Post-Create Script
# ============================================================================
# This script runs after the container is created to set up the development
# environment. It configures git, sets up aliases, and pre-configures CMake.
# ============================================================================

set -e

echo "=========================================="
echo "FSRM Development Environment Setup"
echo "=========================================="

# Configure git if user identity is not set
if ! git config --global user.email > /dev/null 2>&1; then
    echo "Note: Git user.email not configured."
    echo "Run: git config --global user.email 'you@example.com'"
fi

if ! git config --global user.name > /dev/null 2>&1; then
    echo "Note: Git user.name not configured."
    echo "Run: git config --global user.name 'Your Name'"
fi

# Create useful aliases
cat >> ~/.bashrc << 'EOF'

# FSRM Development Aliases
alias build='cmake --build /workspaces/FSRM/build -j$(nproc)'
alias rebuild='rm -rf /workspaces/FSRM/build && mkdir /workspaces/FSRM/build && cd /workspaces/FSRM/build && cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_TESTING=ON -DBUILD_EXAMPLES=ON -DENABLE_CUDA=ON && build'
alias test='ctest --test-dir /workspaces/FSRM/build --output-on-failure'
alias testv='ctest --test-dir /workspaces/FSRM/build -V'
alias sim='./build/examples/simulator -c config/default.config'
alias mpisim='mpirun -np 4 ./build/examples/simulator -c config/default.config'

# Change to workspace
cd /workspaces/FSRM 2>/dev/null || true
EOF

# Create build directory if it doesn't exist
if [ ! -d "/workspaces/FSRM/build" ]; then
    echo "Creating build directory..."
    mkdir -p /workspaces/FSRM/build
fi

# Pre-configure CMake if not already configured
if [ ! -f "/workspaces/FSRM/build/CMakeCache.txt" ]; then
    echo "Pre-configuring CMake..."
    cd /workspaces/FSRM/build
    cmake .. \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DENABLE_TESTING=ON \
        -DBUILD_EXAMPLES=ON \
        -DENABLE_CUDA=ON \
        || echo "CMake configuration will be completed when you open the project."
fi

# Create output directory for simulation results
mkdir -p /workspaces/FSRM/output

# Verify CUDA installation
echo ""
echo "Verifying CUDA installation..."
if command -v nvcc &> /dev/null; then
    echo "  CUDA: $(nvcc --version | grep release | sed 's/.*release //' | sed 's/,.*//')"
    echo "  nvcc: $(which nvcc)"
else
    echo "  Warning: nvcc not found — CUDA toolkit may not be installed"
    echo "  GPU builds will be disabled"
fi

# Verify PETSc installation
echo ""
echo "Verifying PETSc installation..."
if [ -n "$PETSC_DIR" ] && [ -d "$PETSC_DIR" ]; then
    echo "  PETSC_DIR: $PETSC_DIR"
    echo "  PETSC_ARCH: $PETSC_ARCH"
else
    echo "  Warning: PETSC_DIR not set or directory doesn't exist"
fi

# Ensure libgmt.so symlink exists for PyGMT
if [ -f /usr/lib/x86_64-linux-gnu/libgmt.so.6 ] && [ ! -f /usr/lib/x86_64-linux-gnu/libgmt.so ]; then
    sudo ln -sf /usr/lib/x86_64-linux-gnu/libgmt.so.6 /usr/lib/x86_64-linux-gnu/libgmt.so
    sudo ldconfig
fi

# Show available tools
echo ""
echo "Development tools available:"
echo "  - CMake $(cmake --version | head -1 | cut -d' ' -f3)"
echo "  - GCC $(gcc --version | head -1 | cut -d' ' -f4)"
echo "  - GDB $(gdb --version | head -1 | cut -d' ' -f4)"
echo "  - Valgrind $(valgrind --version)"
echo "  - clang-format $(clang-format --version | cut -d' ' -f3)"
echo "  - GMT $(gmt --version 2>/dev/null || echo 'not found')"
echo ""
echo "Python analysis packages:"
python3 -c "import matplotlib; print(f'  - matplotlib {matplotlib.__version__}')" 2>/dev/null || echo "  - matplotlib: not installed"
python3 -c "import scipy; print(f'  - scipy {scipy.__version__}')" 2>/dev/null || echo "  - scipy: not installed"
python3 -c "import pygmt; print(f'  - pygmt {pygmt.__version__}')" 2>/dev/null || echo "  - pygmt: not installed"
python3 -c "import obspy; print(f'  - obspy {obspy.__version__}')" 2>/dev/null || echo "  - obspy: not installed"
python3 -c "import pandas; print(f'  - pandas {pandas.__version__}')" 2>/dev/null || echo "  - pandas: not installed"

echo ""
echo "=========================================="
echo "Quick Start Commands:"
echo "=========================================="
echo "  build     - Build the project"
echo "  rebuild   - Clean and rebuild from scratch"
echo "  test      - Run all tests"
echo "  testv     - Run tests with verbose output"
echo "  sim       - Run simulator with default config"
echo "  mpisim    - Run simulator with MPI (4 processes)"
echo ""
echo "Setup complete!"
