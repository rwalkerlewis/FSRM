# CI/CD Build Fix - Detailed Instructions

## Root Cause Analysis

The CI/CD builds are failing because **PETSc pkg-config is not available** in Ubuntu's petsc-dev package, even though PETSc itself is installed. This is a known limitation of the Ubuntu PETSc packaging.

## Solution Implemented

### 1. **Robust CMakeLists.txt** (Updated)

The CMakeLists.txt now has a 3-tier fallback system for finding PETSc:

```cmake
# Tier 1: Try pkg-config (if available)
pkg_check_modules(PETSC PETSc>=3.15)

# Tier 2: Try environment variables (PETSC_DIR, PETSC_ARCH)
if(NOT PETSC_FOUND)
    if(DEFINED ENV{PETSC_DIR})
        set(PETSC_DIR $ENV{PETSC_DIR})
        # Use PETSC_DIR and PETSC_ARCH
    endif()
endif()

# Tier 3: Search common installation paths
if(NOT PETSC_FOUND)
    find_path(PETSC_INCLUDE_DIR petsc.h
        PATHS /usr/include/petsc /usr/lib/petscdir/petsc-*/...)
    find_library(PETSC_LIBRARY NAMES petsc ...)
endif()
```

### 2. **GitHub Actions Workflow** (Updated)

Added a new step to automatically find and set PETSC_DIR:

```yaml
- name: Set PETSc Environment
  run: |
    if [ -d "/usr/lib/petscdir" ]; then
      PETSC_FOUND=$(find /usr/lib/petscdir -maxdepth 2 -name "petsc-*" -type d | head -1)
      if [ -n "$PETSC_FOUND" ]; then
        echo "PETSC_DIR=$PETSC_FOUND" >> $GITHUB_ENV
        PETSC_ARCH_FOUND=$(ls -1 "$PETSC_FOUND" | grep "x86_64" | head -1)
        if [ -n "$PETSC_ARCH_FOUND" ]; then
          echo "PETSC_ARCH=$PETSC_ARCH_FOUND" >> $GITHUB_ENV
        fi
      fi
    fi
```

This automatically finds PETSc in Ubuntu's standard location: `/usr/lib/petscdir/petsc-X.XX/x86_64-linux-gnu-real/`

### 3. **Enhanced Debugging**

Added verbose output to help diagnose issues:
- Lists all installed compilers and versions
- Searches for petsc.h and libpetsc.* files
- Shows PETSC_DIR if found
- Enables CMAKE_VERBOSE_MAKEFILE for detailed build output

## Expected Build Flow

### Step 1: Install Dependencies
```bash
sudo apt-get install -y \
  build-essential cmake g++ gcc gfortran \
  mpich libmpich-dev libhdf5-mpich-dev \
  libgtest-dev petsc-dev libpetsc-real-dev \
  pkg-config
```

### Step 2: Find PETSc
```bash
# Ubuntu installs PETSc here:
/usr/lib/petscdir/petsc-3.XX/x86_64-linux-gnu-real/
```

The workflow automatically finds this and sets `PETSC_DIR` and `PETSC_ARCH`.

### Step 3: Configure CMake
```bash
CC=gcc CXX=g++ cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TESTING=ON \
  -DCMAKE_VERBOSE_MAKEFILE=ON
```

CMake will:
1. Try pkg-config for PETSc (fails, but that's OK)
2. Check PETSC_DIR environment variable (succeeds!)
3. Set appropriate include dirs and libraries

### Step 4: Build
```bash
make -j$(nproc)
```

Should complete successfully.

## Manual Testing

To test this locally:

```bash
# Install dependencies
sudo apt-get update
sudo apt-get install -y \
    build-essential cmake g++ gcc gfortran \
    mpich libmpich-dev libhdf5-mpich-dev \
    petsc-dev libpetsc-real-dev \
    pkg-config libgtest-dev

# Find PETSc installation
find /usr -name "petsc.h" 2>/dev/null
find /usr -name "libpetsc.so" 2>/dev/null

# Set PETSC_DIR (example)
export PETSC_DIR=/usr/lib/petscdir/petsc-3.18
export PETSC_ARCH=x86_64-linux-gnu-real

# Build
mkdir build && cd build
CC=gcc CXX=g++ cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_TESTING=ON \
    -DCMAKE_VERBOSE_MAKEFILE=ON

make -j$(nproc)
```

## What Changed

### Files Modified:

1. **CMakeLists.txt**
   - Added robust PETSc finding with 3 fallback methods
   - Added HDF5 fallback detection
   - Better error messages

2. **.github/workflows/tests.yml**
   - Added "Set PETSc Environment" step to all jobs
   - Added dependency verification with file searches
   - Added CMAKE_VERBOSE_MAKEFILE for debugging
   - Added || true to prevent failures during verification

## Common PETSc Locations by Platform

| Platform | Typical PETSc Location |
|----------|------------------------|
| **Ubuntu/Debian** | `/usr/lib/petscdir/petsc-X.XX/x86_64-linux-gnu-real/` |
| **Fedora/RHEL** | `/usr/lib64/petsc/` or `/usr/lib/petsc/` |
| **macOS (Homebrew)** | `/opt/homebrew/opt/petsc/` or `/usr/local/opt/petsc/` |
| **Custom Install** | Use `PETSC_DIR` and `PETSC_ARCH` environment variables |

## Troubleshooting

### If build still fails:

1. **Check PETSc installation in logs**:
   ```
   === Finding PETSc installation ===
   /usr/lib/petscdir/petsc-3.18/x86_64-linux-gnu-real/include/petsc.h
   ```

2. **Verify PETSC_DIR was set**:
   ```
   Found PETSc at: /usr/lib/petscdir/petsc-3.18
   ```

3. **Check CMake configuration output**:
   ```
   -- Found PETSc: /usr/lib/libpetsc.so
   ```

### If MPI not found:

```bash
# Reinstall MPI
sudo apt-get install --reinstall mpich libmpich-dev

# Verify
which mpicc mpicxx
mpicc --version
```

### If HDF5 not found:

HDF5 is now optional (warning only), so this shouldn't block the build.

## Success Criteria

Build is successful when you see:

```
-- Found PETSc: /usr/lib/...
-- Configuring done
-- Generating done
-- Build files written to: /workspace/build
[100%] Built target fsrm
[100%] Built target fsrmlib
```

## Next Steps After Fix

Once the build succeeds:

1. ✅ All 6 CI/CD jobs should complete
2. ✅ Unit tests will run
3. ✅ Integration tests will run  
4. ✅ DirectionalWell tests will run
5. ✅ Code coverage will be generated
6. ✅ Performance benchmarks will run

The PR will then be ready for merge!

---

**Summary**: The issue was Ubuntu's PETSc package not providing pkg-config support. The fix adds automatic PETSC_DIR detection and multiple fallback methods for finding PETSc.
