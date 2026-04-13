# Building FSRM on Windows

FSRM requires Linux-specific dependencies (PETSc, MPI, HDF5) that are not natively available on Windows. You have two options:

## Option 1: WSL (Windows Subsystem for Linux) - Recommended for Development

WSL allows you to run a full Linux environment directly on Windows.

### Prerequisites
- Windows 10/11 with WSL2 enabled
- **System reboot required** if WSL was just installed

### Setup Steps

1. **Reboot your computer** (required after WSL installation)

2. **Start Ubuntu in WSL:**
   ```powershell
   wsl
   ```

3. **Run the automated build script:**
   ```bash
   cd /mnt/c/Users/User/Projects/FSRM
   bash build_wsl.sh
   ```

   Or manually:
   ```bash
   # Install dependencies
   sudo apt update
   sudo apt install -y build-essential cmake libopenmpi-dev libhdf5-openmpi-dev libpetsc-real-dev

   # Configure environment
   export PETSC_DIR=/usr/lib/petsc
   export PETSC_ARCH=""

   # Build
   cd /mnt/c/Users/User/Projects/FSRM
   mkdir -p build && cd build
   cmake .. -DBUILD_EXAMPLES=ON -DENABLE_TESTING=ON
   make -j$(nproc)

   # Test
   ctest
   ```

4. **Run examples:**
   ```bash
   cd /mnt/c/Users/User/Projects/FSRM/build/examples
   ./simulator -c ../../config/default.config
   ```

### Advantages
- Native Linux performance
- Full access to all FSRM features
- Easy debugging and development
- Direct file access between Windows and Linux

### Disadvantages
- Requires reboot for initial setup
- Uses more system resources

---

## Option 2: Docker - Easiest for Quick Testing

Docker provides a containerized Linux environment.

### Prerequisites
- Docker Desktop for Windows
- Download from: https://www.docker.com/products/docker-desktop/

### Setup Steps

1. **Install and start Docker Desktop**

2. **Run the automated build script:**
   ```powershell
   cd c:\Users\User\Projects\FSRM
   .\build_docker.ps1
   ```

   Or manually:
   ```powershell
   # Build image
   docker build -t fsrm:latest .

   # Start container
   docker-compose up -d

   # Run tests
   docker exec fsrm bash -c "cd /app/build && ctest"

   # Enter container for interactive work
   docker exec -it fsrm bash
   ```

3. **Run examples inside container:**
   ```bash
   cd /app/build/examples
   ./simulator -c /app/config/default.config
   ```

### Advantages
- No reboot required
- Isolated environment
- Easy to share/deploy
- Works on any Docker-capable system

### Disadvantages
- Slower than native WSL
- Additional abstraction layer
- File I/O overhead for shared volumes

---

## Current System Status

**WSL Status:** Installed, requires reboot  
**Docker Status:** Not installed  
**Build Tools:** None available on Windows natively

---

## Recommended Next Steps

1. **For development work:** Reboot, then use WSL with `build_wsl.sh`
2. **For quick testing:** Install Docker Desktop, then use `build_docker.ps1`

---

## Running Examples After Build

### Available Examples
- **SPE Benchmarks:** `spe1` through `spe11` (reservoir simulation standards)
- **SCEC Benchmarks:** `scec_tpv*`, `scec_loh*` (earthquake dynamics)
- **Config-driven:** `simulator` (unified executable for all configs)

### Example Commands

```bash
# Simple single-phase flow
./simulator -c /path/to/config/default.config

# Hydraulic fracturing with MPI (4 processes)
mpirun -np 4 ./simulator -c /path/to/config/hydraulic_fracturing.config

# SCEC earthquake benchmark
./scec_tpv10

# SPE reservoir benchmark
./spe1
```

### Test Suite

```bash
cd build
ctest                          # Run all tests
ctest -L unit                  # Unit tests only
ctest -R scec_tpv10 -V        # Specific test with verbose output
```

---

## Troubleshooting

### WSL Issues
- **"Virtual Machine Platform" error:** Run as Administrator:
  ```powershell
  wsl.exe --install --no-distribution
  ```
  Then reboot.

- **Ubuntu won't start:** Enable virtualization in BIOS/UEFI

### Docker Issues
- **Docker daemon not running:** Start Docker Desktop
- **Build fails:** Ensure you have at least 4GB RAM and 20GB disk space
- **Permission errors:** Run Docker Desktop as Administrator

### Build Issues
- **PETSc not found:** Check `PETSC_DIR` and `PETSC_ARCH` environment variables
- **MPI errors:** Ensure `libopenmpi-dev` is installed
- **Out of memory:** Reduce parallel jobs: `make -j2` instead of `make -j$(nproc)`

---

## Performance Comparison

| Method | Build Time | Runtime Performance | Development Experience |
|--------|-----------|-------------------|----------------------|
| WSL    | ~10-15 min | Native (100%) | Excellent |
| Docker | ~15-20 min | Good (90-95%) | Good |

---

For more information, see:
- [DEVELOPMENT.md](docs/DEVELOPMENT.md) - Full development guide
- [README.md](README.md) - Project overview
- [examples/README.md](examples/README.md) - Example documentation
