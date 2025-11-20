# FSRM Testing Quick Reference

## Quick Start

```bash
# Build and run all tests
cd /home/dockimble/Projects/FSRM
mkdir -p build && cd build
cmake .. -DENABLE_TESTING=ON
make -j$(nproc)
ctest

# Or use the convenient script
cd /home/dockimble/Projects/FSRM
./tests/run_tests.sh --quick
```

## Test Categories

| Category | Command | Duration | Description |
|----------|---------|----------|-------------|
| **Quick Unit Tests** | `./tests/run_tests.sh --quick` | ~1 min | Fast unit tests for development |
| **Convergence Tests** | `./tests/run_tests.sh --convergence` | ~5 min | MMS and convergence rate verification |
| **Integration Tests** | `./tests/run_tests.sh --integration` | ~10 min | Full simulation workflows |
| **All Tests** | `./tests/run_tests.sh` | ~15 min | Complete test suite |

## Common Test Commands

### Run Specific Physics Tests

```bash
# Single phase flow
ctest -R SinglePhaseFlow

# Poroelasticity
ctest -R Poroelasticity

# Elastodynamics (wave propagation)
ctest -R Elastodynamics

# Hydraulic fracturing
ctest -R HydraulicFracture

# Fracture mechanics
ctest -R FractureModel
```

### Run with MPI

```bash
# 4 processes
./tests/run_tests.sh --nprocs 4

# Quick tests with 2 processes
./tests/run_tests.sh --quick --nprocs 2
```

### Verbose Output

```bash
# See detailed test output
./tests/run_tests.sh --verbose

# Or with ctest
ctest -V
```

## Test Files Overview

| File | Purpose | Key Tests |
|------|---------|-----------|
| `test_main.cpp` | Basic framework | Unit test infrastructure |
| `test_physics_kernels.cpp` | Physics kernels | Residual/Jacobian evaluation, wave speeds |
| `test_fracture_model.cpp` | Fracture mechanics | LEFM, KGD model, network connectivity |
| `test_convergence_mms.cpp` | Convergence | Spatial O(h²), temporal O(dt²) convergence |
| `test_integration.cpp` | Full simulations | Workflows, benchmarks, regressions |

## What Gets Tested

### Physics Kernels ✓
- [x] Single phase flow (Darcy)
- [x] Black oil model
- [x] Poroelasticity (Biot)
- [x] Elastodynamics (wave equation)
- [x] Poroelastodynamics (Biot waves)
- [x] Thermal flow

### Fracture Mechanics ✓
- [x] LEFM stress intensity factors
- [x] Fracture propagation criteria
- [x] Hydraulic fracture growth (KGD)
- [x] Fracture width/conductivity
- [x] Network connectivity
- [x] Proppant transport
- [x] Leak-off (Carter model)

### Numerical Methods ✓
- [x] Spatial convergence (FEM)
- [x] Temporal convergence (time stepping)
- [x] MMS verification
- [x] Jacobian accuracy (vs finite differences)

### Workflows ✓
- [x] Simple domain simulation
- [x] Wells (injection/production)
- [x] Consolidation
- [x] Wave propagation
- [x] Induced seismicity
- [x] Restart/checkpoint

## Expected Results

### Convergence Rates
- **Spatial**: O(h²) for linear finite elements
- **Temporal**: 
  - O(dt) for Backward Euler
  - O(dt²) for BDF2/Crank-Nicolson

### Conservation Errors
- Mass balance: < 1%
- Energy conservation: < 0.1% (elastodynamics)

### Performance Targets
- Single phase flow: > 10,000 DOFs/sec
- Poroelasticity: > 5,000 DOFs/sec
- Elastodynamics: > 1,000 DOFs/sec

## Debugging Failed Tests

### View test output
```bash
cd build
cat Testing/Temporary/LastTest.log
```

### Run specific test with GDB
```bash
cd build
gdb --args ./tests/run_tests --gtest_filter=MyFailingTest*
```

### Check memory with Valgrind
```bash
cd build
valgrind --leak-check=full ./tests/run_tests --gtest_filter=MyTest*
```

### PETSc debugging
```bash
export PETSC_OPTIONS="-on_error_attach_debugger -log_view"
./tests/run_tests
```

## Adding Your Own Tests

### 1. Create test file
```bash
touch tests/test_my_feature.cpp
```

### 2. Add test
```cpp
#include <gtest/gtest.h>
#include "MyFeature.hpp"

TEST(MyFeatureTest, BasicTest) {
    MyFeature feature;
    EXPECT_EQ(feature.compute(), 42);
}
```

### 3. Rebuild and run
```bash
cd build
make run_tests
ctest -R MyFeature
```

## Continuous Integration

Tests run automatically on:
- Every commit to main/develop
- Every pull request
- Nightly builds (2 AM UTC)

View results: GitHub Actions tab

## Test Coverage

Generate coverage report:
```bash
cd build
cmake .. -DENABLE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug
make
ctest
lcov --capture --directory . -o coverage.info
genhtml coverage.info -o coverage_html
firefox coverage_html/index.html
```

## Benchmark Comparisons

Track performance over time:
```bash
# Run benchmarks
ctest -R Performance

# Compare with baseline
python3 scripts/compare_benchmarks.py \
    baseline_results.json \
    current_results.json
```

## Common Issues

### Test timeout
Increase timeout in `tests/CMakeLists.txt`:
```cmake
set_tests_properties(MyTest PROPERTIES TIMEOUT 300)
```

### MPI issues
```bash
# Check MPI installation
which mpirun
mpirun --version

# Test MPI hello world
mpirun -np 2 echo "Hello from MPI"
```

### Missing dependencies
```bash
# Ubuntu/Debian
sudo apt-get install libgtest-dev libhdf5-dev

# Build GoogleTest if needed
cd /usr/src/gtest
sudo cmake .
sudo make
sudo cp lib/*.a /usr/lib
```

## Resources

- Full test documentation: `tests/README.md`
- GoogleTest docs: https://google.github.io/googletest/
- PETSc testing: https://petsc.org/release/docs/manual/developers/
- MMS methodology: Roache (2002)

## Need Help?

- Check `tests/README.md` for detailed documentation
- View example tests in `test_*.cpp` files
- Open an issue on GitHub
- Contact the development team
