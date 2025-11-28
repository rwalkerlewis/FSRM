# FSRM Development Guide

This guide covers building, testing, and contributing to FSRM.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Building from Source](#building-from-source)
- [Testing](#testing)
- [Continuous Integration](#continuous-integration)
- [Development Workflow](#development-workflow)
- [Code Style](#code-style)
- [Contributing](#contributing)

## Prerequisites

### Required Dependencies

- **CMake** >= 3.15
- **C++17 compiler**: GCC 7+, Clang 6+, or MSVC 2019+
- **PETSc** >= 3.15 with MPI support
- **MPI**: OpenMPI, MPICH, or Intel MPI

### Optional Dependencies

- **HDF5** >= 1.10 (with MPI support) - For efficient I/O
- **GTest** >= 1.10 - For unit testing
- **gnuplot** - For automatic plotting
- **ParaView** >= 5.0 - For visualization
- **PROJ** >= 6.0 - For coordinate transformations
- **Gmsh** - For mesh generation

### GPU Support (Optional)

For NVIDIA GPUs:
- **CUDA Toolkit** >= 11.0
- Compatible GPU: Pascal (SM 6.0) or newer
- Recommended: CUDA 12.0+ for best performance

For AMD GPUs:
- **ROCm** >= 5.0
- Compatible GPU: MI50, MI100, MI200 series

## Building from Source

### Ubuntu/Debian

```bash
# Install dependencies
sudo apt update
sudo apt install -y \
    cmake \
    g++ \
    libpetsc-real-dev \
    libhdf5-mpi-dev \
    libgtest-dev \
    gnuplot \
    libproj-dev

# Set PETSc environment
export PETSC_DIR=/usr/lib/petsc
export PETSC_ARCH=""

# Clone and build
git clone https://github.com/your-repo/fsrm.git
cd fsrm
mkdir build && cd build

# Configure (CPU only)
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# Test
make test

# Install
sudo make install
```

### macOS

```bash
# Install dependencies with Homebrew
brew install cmake gcc open-mpi petsc hdf5-mpi googletest gnuplot proj

# Set environment
export PETSC_DIR=$(brew --prefix petsc)
export PETSC_ARCH=""

# Clone and build
git clone https://github.com/your-repo/fsrm.git
cd fsrm
mkdir build && cd build

# Configure
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_C_COMPILER=gcc-13 \
         -DCMAKE_CXX_COMPILER=g++-13

# Build
make -j$(sysctl -n hw.ncpu)

# Test
make test
```

### RHEL/CentOS

```bash
# Enable EPEL
sudo yum install -y epel-release

# Install dependencies
sudo yum install -y \
    cmake3 \
    gcc-c++ \
    petsc-openmpi-devel \
    hdf5-openmpi-devel \
    gtest-devel \
    gnuplot \
    proj-devel

# Load MPI
module load mpi/openmpi-x86_64

# Set PETSc environment
export PETSC_DIR=/usr/lib64/openmpi
export PETSC_ARCH=""

# Build as above
```

## CMake Build Options

### Basic Options

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release|Debug|RelWithDebInfo \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DENABLE_TESTING=ON|OFF \
  -DENABLE_EXAMPLES=ON|OFF
```

### GPU Options

```bash
# NVIDIA CUDA
cmake .. \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="80;86;89" \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda

# AMD ROCm
cmake .. \
  -DENABLE_HIP=ON \
  -DHIP_ROOT_DIR=/opt/rocm
```

### Advanced Options

```bash
cmake .. \
  -DENABLE_PROJ=ON \              # Coordinate transformation support
  -DENABLE_GMSH=ON \              # Gmsh mesh I/O
  -DENABLE_PROFILING=ON \         # Performance profiling
  -DBUILD_SHARED_LIBS=ON \        # Build shared libraries
  -DUSE_DOUBLE_PRECISION=ON \     # Use double precision (default)
  -DENABLE_OPENMP=ON              # OpenMP threading
```

## GPU Build

### NVIDIA CUDA

```bash
mkdir build && cd build

# Configure with CUDA
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES="80;86"  # Adjust for your GPU

# Detect GPU architecture automatically
nvidia-smi --query-gpu=compute_cap --format=csv,noheader

# Build
make -j$(nproc)

# Verify GPU support
./fsrm --gpu-info
```

**Supported Architectures**:
- `60`: Pascal (P100, GTX 10 series)
- `70`: Volta (V100)
- `75`: Turing (RTX 20 series, T4)
- `80`: Ampere (A100, RTX 30 series)
- `86`: Ampere (RTX 3090, A40)
- `89`: Ada Lovelace (RTX 40 series)
- `90`: Hopper (H100)

### AMD ROCm

```bash
# Install ROCm
wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | sudo apt-key add -
echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/debian/ ubuntu main' | sudo tee /etc/apt/sources.list.d/rocm.list
sudo apt update
sudo apt install -y rocm-dev hipcc

# Configure with HIP
mkdir build && cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_HIP=ON \
  -DCMAKE_CXX_COMPILER=hipcc

# Build
make -j$(nproc)

# Verify
./fsrm --gpu-info
```

## Testing

### Unit Tests

```bash
cd build

# Run all tests
make test
# or
ctest

# Run with verbose output
ctest -V

# Run specific test
ctest -R test_unit_system

# Run tests in parallel
ctest -j$(nproc)
```

### Manual Test Execution

```bash
cd build/tests

# Run individual test
./test_unit_system
./test_config_reader
./test_fluid_models
./test_fracture_models

# With GTest filters
./test_unit_system --gtest_filter=UnitSystemTest.*
./test_unit_system --gtest_filter=*Permeability*
```

### Integration Tests

```bash
# Run method of manufactured solutions tests
cd tests
./run_mms_tests.sh

# Run SCEC benchmark suite
cd benchmarks
./run_scec_suite.sh --all

# Compare with reference solutions
./run_scec_suite.sh --all --verify
```

### Memory Leak Detection

```bash
# With Valgrind
valgrind --leak-check=full --show-leak-kinds=all \
  ./test_unit_system

# With AddressSanitizer
cmake .. -DCMAKE_BUILD_TYPE=Debug \
         -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer"
make
./test_unit_system
```

### Coverage Analysis

```bash
# Configure with coverage
cmake .. -DCMAKE_BUILD_TYPE=Debug \
         -DENABLE_COVERAGE=ON

# Build and run tests
make
make test

# Generate coverage report
make coverage

# View report
firefox coverage/index.html
```

## Continuous Integration

### GitHub Actions

FSRM uses GitHub Actions for CI/CD. Configuration: `.github/workflows/ci.yml`

**CI Pipeline**:
1. **Build Testing**: Ubuntu, macOS, Windows
2. **Compiler Testing**: GCC, Clang, MSVC
3. **Unit Tests**: All test suites
4. **Benchmark Validation**: Key benchmarks
5. **GPU Tests**: CUDA builds (GPU runners)
6. **Documentation**: Build and deploy docs

### Running CI Locally

```bash
# Install act (GitHub Actions local runner)
brew install act  # macOS
# or download from https://github.com/nektos/act

# Run CI locally
act push

# Run specific job
act -j build-linux
```

### CI Environment Variables

Set these secrets in GitHub repository settings:
- `DOCKER_USERNAME`: Docker Hub username (for GPU images)
- `DOCKER_PASSWORD`: Docker Hub token
- `CODECOV_TOKEN`: Codecov token (for coverage)

## Development Workflow

### Branch Strategy

- `main`: Stable release branch
- `develop`: Development integration branch
- `feature/*`: New feature branches
- `bugfix/*`: Bug fix branches
- `hotfix/*`: Production hotfixes

### Creating a New Feature

```bash
# Create feature branch
git checkout -b feature/my-new-feature develop

# Make changes and commit
git add .
git commit -m "Add my new feature"

# Push to remote
git push -u origin feature/my-new-feature

# Create pull request on GitHub
```

### Running Pre-commit Checks

```bash
# Format code
clang-format -i src/*.cpp include/*.hpp

# Run linter
cppcheck --enable=all --inconclusive src/ include/

# Run tests
cd build && make test

# Check for memory leaks
valgrind --leak-check=full ./test_suite
```

## Code Style

### C++ Style Guide

FSRM follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) with modifications:

**Naming Conventions**:
- Classes: `PascalCase` (e.g., `FluidModel`, `FaultModel`)
- Functions: `camelCase` (e.g., `computeFlux`, `updateState`)
- Variables: `snake_case` (e.g., `permeability_x`, `total_stress`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `MAX_ITERATIONS`, `PI`)
- File names: `PascalCase.cpp` / `PascalCase.hpp`

**Formatting**:
- Indentation: 2 spaces (no tabs)
- Line length: 100 characters maximum
- Braces: Allman style (opening brace on new line)
- Comments: Use `//` for single line, `/* */` for multi-line

**Example**:
```cpp
class FluidModel
{
public:
  FluidModel(double density, double viscosity)
    : density_(density),
      viscosity_(viscosity)
  {
  }

  double computeVelocity(double pressure_gradient)
  {
    const double permeability = 1.0e-13;  // m^2
    return -permeability / viscosity_ * pressure_gradient;
  }

private:
  double density_;
  double viscosity_;
};
```

### Documentation

Use Doxygen-style comments:

```cpp
/**
 * @brief Compute Darcy flux from pressure gradient
 * 
 * @param pressure_gradient Pressure gradient in Pa/m
 * @param permeability Permeability in m^2
 * @param viscosity Dynamic viscosity in Pa·s
 * @return Darcy flux in m/s
 */
double computeDarcyFlux(
  double pressure_gradient,
  double permeability,
  double viscosity
);
```

### Configuration Files

- Use INI format with `[SECTIONS]`
- Include comments explaining parameters
- Specify units explicitly (e.g., `permeability = 150 mD`)
- Group related parameters in sections

### Git Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <subject>

<body>

<footer>
```

**Types**:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting)
- `refactor`: Code refactoring
- `test`: Adding or modifying tests
- `perf`: Performance improvements
- `build`: Build system changes
- `ci`: CI configuration changes

**Example**:
```
feat(gpu): add CUDA support for elastodynamics kernel

Implement GPU-accelerated elastodynamics kernel using CUDA.
Achieves 30x speedup on large problems (>1M cells).

Closes #123
```

## Contributing

### Reporting Bugs

1. Check if bug already reported in [Issues](https://github.com/your-repo/fsrm/issues)
2. Create new issue with:
   - Clear title
   - Steps to reproduce
   - Expected vs actual behavior
   - System information (OS, compiler, PETSc version)
   - Configuration file (if applicable)
   - Error messages / logs

### Suggesting Features

1. Check [Issues](https://github.com/your-repo/fsrm/issues) for existing requests
2. Create new issue with:
   - Clear title and description
   - Use case / motivation
   - Proposed implementation (if any)
   - References (papers, other software)

### Pull Requests

1. Fork the repository
2. Create feature branch from `develop`
3. Make changes with tests and documentation
4. Ensure all tests pass
5. Submit pull request with:
   - Clear description of changes
   - Link to related issues
   - Screenshots (if UI changes)

**PR Checklist**:
- [ ] Code follows style guide
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] All tests pass
- [ ] No compiler warnings
- [ ] Commits follow convention
- [ ] Branch up to date with `develop`

## Troubleshooting

### PETSc Not Found

```bash
# Set environment variables
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=linux-gnu-c-debug  # or your arch

# Or specify in CMake
cmake .. -DPETSC_DIR=/path/to/petsc
```

### MPI Errors

```bash
# Check MPI installation
which mpirun
mpirun --version

# Try single process first
./fsrm -c config.ini

# Use correct mpirun
/usr/lib64/openmpi/bin/mpirun -np 4 ./fsrm -c config.ini
```

### CUDA Compilation Errors

```bash
# Check CUDA version
nvcc --version
nvidia-smi

# Specify CUDA path
cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.0

# Try older CUDA if needed
module load cuda/11.8
```

### Test Failures

```bash
# Run single test with verbose output
ctest -V -R test_name

# Check for library conflicts
ldd ./test_unit_system

# Run with debugging
gdb ./test_unit_system
```

## Performance Profiling

### CPU Profiling

```bash
# With gprof
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo \
         -DCMAKE_CXX_FLAGS="-pg"
make
./fsrm -c config.ini
gprof ./fsrm gmon.out > profile.txt

# With perf
perf record -g ./fsrm -c config.ini
perf report
```

### GPU Profiling

```bash
# NVIDIA Nsight Compute
ncu --set full -o profile ./fsrm -c config.ini --use-gpu

# NVIDIA Nsight Systems
nsys profile -o profile ./fsrm -c config.ini --use-gpu

# View results
ncu-ui profile.ncu-rep
nsys-ui profile.qdrep
```

### PETSc Profiling

```bash
# Enable PETSc logging
./fsrm -c config.ini -log_view

# Detailed profiling
./fsrm -c config.ini -log_view :profile.txt:ascii_info_detail
```

## Building Documentation

```bash
# Install Doxygen
sudo apt install doxygen graphviz

# Generate docs
cd build
make docs

# View
firefox docs/html/index.html
```

## Docker Development

```bash
# Build Docker image
docker build -t fsrm:dev .

# Run container
docker run -it --rm \
  -v $(pwd):/workspace \
  fsrm:dev

# With GPU support
docker run --gpus all -it --rm \
  -v $(pwd):/workspace \
  fsrm:dev
```

## Resources

- [PETSc Documentation](https://petsc.org/release/docs/manual/)
- [CMake Documentation](https://cmake.org/documentation/)
- [CUDA Programming Guide](https://docs.nvidia.com/cuda/cuda-c-programming-guide/)
- [Google Test Documentation](https://google.github.io/googletest/)
- [FSRM GitHub](https://github.com/your-repo/fsrm)

## See Also

- [Quick Start](QUICK_START.md) - Get started quickly
- [User Guide](USER_GUIDE.md) - Running simulations
- [Benchmarks](BENCHMARKS.md) - Validation and testing
- [Configuration Reference](CONFIGURATION.md) - All options
- [Deployment Guide](DEPLOYMENT.md) - Cloud and HPC deployment
