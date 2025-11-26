# Build Requirements for FSRM

## Required System Dependencies

The FSRM project requires the following system packages to build successfully:

### Essential Build Tools
- **CMake** >= 3.15
- **C++17 compiler** (g++ or clang++)
- **GNU Make** or **Ninja**

### Required Libraries

1. **MPI (Message Passing Interface)**
   - OpenMPI or MPICH
   - Install commands:
     ```bash
     # Ubuntu/Debian
     sudo apt-get install libopenmpi-dev openmpi-bin
     
     # CentOS/RHEL
     sudo yum install openmpi-devel
     
     # macOS
     brew install open-mpi
     ```

2. **PETSc** >= 3.15 (Portable, Extensible Toolkit for Scientific Computation)
   - Required for parallel computing and scientific solvers
   - Install commands:
     ```bash
     # Ubuntu/Debian
     sudo apt-get install petsc-dev
     
     # Or build from source:
     git clone https://gitlab.com/petsc/petsc.git
     cd petsc
     ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran \
                 --download-fblaslapack --download-mpich
     make all test
     ```

3. **HDF5** (with MPI support)
   - Required for output
   - Install commands:
     ```bash
     # Ubuntu/Debian
     sudo apt-get install libhdf5-openmpi-dev
     
     # CentOS/RHEL
     sudo yum install hdf5-openmpi-devel
     ```

4. **pkg-config**
   - Required for finding libraries
   - Usually pre-installed, or:
     ```bash
     sudo apt-get install pkg-config
     ```

### Optional Dependencies

1. **CUDA Toolkit** >= 11.0 (for GPU support)
   - Only needed if building with `-DENABLE_CUDA=ON`
   
2. **ROCm/HIP** (for AMD GPU support)
   - Only needed if building with `-DENABLE_HIP=ON`

3. **Google Test** (for testing)
   - Only needed if building with `-DENABLE_TESTING=ON`
   - Install:
     ```bash
     sudo apt-get install libgtest-dev
     ```

4. **gnuplot** (for visualization)
   - Optional for plotting
   - Install:
     ```bash
     sudo apt-get install gnuplot
     ```

## Building the Project

### Standard Build (CPU only)

```bash
# Install dependencies first (Ubuntu example)
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    libopenmpi-dev \
    openmpi-bin \
    petsc-dev \
    libhdf5-openmpi-dev \
    pkg-config

# Configure
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_CUDA=OFF

# Build
make -j$(nproc)

# Run tests
make test

# Install
sudo make install
```

### Build with CUDA Support

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DENABLE_CUDA=ON \
         -DCMAKE_CUDA_ARCHITECTURES="80;86"

make -j$(nproc)
```

## Troubleshooting

### CMake cannot find MPI
- Ensure MPI is installed: `which mpicc mpicxx`
- Set MPI compilers explicitly:
  ```bash
  cmake .. -DMPI_C_COMPILER=mpicc -DMPI_CXX_COMPILER=mpicxx
  ```

### CMake cannot find PETSc
- Check PETSc pkg-config: `pkg-config --modversion PETSc`
- Set PETSc directory:
  ```bash
  export PETSC_DIR=/path/to/petsc
  export PETSC_ARCH=arch-linux-c-debug
  ```

### Linker errors with libstdc++
- Ensure C++ standard library is installed:
  ```bash
  sudo apt-get install libstdc++-12-dev
  ```
- Or use g++ instead of clang++:
  ```bash
  CC=gcc CXX=g++ cmake ..
  ```

## Docker Build

For a reproducible build environment, use Docker:

```bash
docker build -t fsrm:latest .
docker run -it fsrm:latest
```

## CI/CD Build

The project includes GitHub Actions workflows in `.github/workflows/` that demonstrate the complete build process in a clean environment.
