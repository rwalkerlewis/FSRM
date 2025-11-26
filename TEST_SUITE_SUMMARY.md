# FSRM Test Suite Summary

## Overview

The FSRM (Full-Scale Reservoir Simulator with Mechanics) test suite has been organized into five categories to ensure comprehensive code coverage and functionality verification.

## Test Categories

### 1. Unit Tests (`tests/unit/`)
Fast, isolated tests for individual components.

| Test File | Description | Tests |
|-----------|-------------|-------|
| `test_config_reader.cpp` | ConfigReader class parsing and retrieval | 12 tests |
| `test_physics_kernels.cpp` | Physics kernel creation, property setting, mathematical verification | 18 tests |
| `test_well_model.cpp` | Well types, completions, well index calculations, Peaceman formula | 16 tests |
| `test_fracture_model.cpp` | Fracture network, hydraulic fractures, faults, stress analysis | 25 tests |
| `test_eclipse_io.cpp` | Eclipse file I/O, parsing, grid config retrieval | 14 tests |

### 2. Functional Tests (`tests/functional/`)
Tests for specific features and module integration.

| Test File | Description | Tests |
|-----------|-------------|-------|
| `test_simulator_init.cpp` | Simulator configuration and initialization | 10 tests |
| `test_solver_convergence.cpp` | Numerical solver convergence (Newton-Raphson, Jacobi, etc.) | 15 tests |
| `test_boundary_conditions.cpp` | Boundary condition formulations | 10 tests |
| `test_well_operations.cpp` | Well calculations (Peaceman, PI, IPR, pressure drops) | 15 tests |

### 3. Physics/MMS Tests (`tests/physics/`)
Method of Manufactured Solutions tests for numerical verification.

| Test File | Description | Tests |
|-----------|-------------|-------|
| `test_mms_diffusion.cpp` | Diffusion equation MMS tests (spatial/temporal convergence) | 10 tests |
| `test_mms_elasticity.cpp` | Linear elasticity patch tests and convergence | 10 tests |
| `test_mms_wave_propagation.cpp` | Wave equation verification (P-waves, S-waves, energy) | 10 tests |
| `test_analytical_solutions.cpp` | Validation against Theis, Terzaghi, Buckley-Leverett, etc. | 10 tests |

### 4. Integration Tests (`tests/integration/`)
End-to-end workflow tests.

| Test File | Description | Tests |
|-----------|-------------|-------|
| `test_full_simulation.cpp` | Complete simulation pipeline tests | 10 tests |
| `test_restart.cpp` | Checkpoint/restart functionality | 8 tests |

### 5. Performance Tests (`tests/performance/`)
Benchmarking and scaling tests.

| Test File | Description | Tests |
|-----------|-------------|-------|
| `test_benchmarks.cpp` | Kernel evaluation timing, memory access patterns | 10 tests |
| `test_scaling.cpp` | MPI parallel scaling (strong/weak) | 8 tests |

## Running Tests

### Run All Tests
```bash
cd build
ctest --output-on-failure
```

### Run by Category
```bash
# Unit tests only (fast)
ctest -L "unit" --output-on-failure

# Functional tests
ctest -L "functional" --output-on-failure

# Physics/MMS tests
ctest -L "physics" --output-on-failure

# Integration tests
ctest -L "integration" --output-on-failure

# Performance tests
ctest -L "performance" --output-on-failure
```

### Run Quick Tests (unit + functional)
```bash
ctest -L "quick" --output-on-failure
```

### Run Specific Test
```bash
ctest -R "Unit.ConfigReader" --output-on-failure
```

## Test Labels

| Label | Description |
|-------|-------------|
| `unit` | Unit tests - isolated component tests |
| `functional` | Functional tests - feature tests |
| `physics` | Physics tests - MMS and analytical solutions |
| `mms` | Method of Manufactured Solutions tests |
| `integration` | Integration tests - end-to-end workflows |
| `performance` | Performance/benchmark tests |
| `quick` | Fast tests (unit + some functional) |
| `slow` | Slow tests (integration + performance) |

## CI/CD Integration

The GitHub Actions workflow (`.github/workflows/tests.yml`) runs tests in stages:

1. **Build** - Compiles the project and caches artifacts
2. **Unit Tests** - Runs first (fast feedback)
3. **Functional Tests** - Runs in parallel with unit tests
4. **Physics Tests** - Runs in parallel with other tests
5. **Integration Tests** - Runs after unit and functional tests pass
6. **Performance Tests** - Only on main branch pushes

## Test Results (Local Run)

```
Test project /workspace/build
      Start  1: Unit.ConfigReader ................   Passed
      Start  2: Unit.PhysicsKernels ..............   Passed
      Start  3: Unit.WellModel ...................   Passed
      Start  4: Unit.FractureModel ...............   Passed
      Start  5: Unit.EclipseIO ...................   Passed
      Start  6: Functional.SimulatorInit .........   Passed
      Start  7: Functional.SolverConvergence .....   Passed
      Start  8: Functional.BoundaryConditions ....   Passed
      Start  9: Functional.WellOperations ........   Passed
      Start 10: Physics.MMS.Diffusion ............   Passed
      Start 11: Physics.MMS.Elasticity ...........   Passed
      Start 12: Physics.MMS.WavePropagation ......   Passed
      Start 13: Physics.AnalyticalSolutions ......   Passed
      Start 14: Integration.FullSimulation .......   Passed
      Start 15: Integration.Restart ..............   Passed
      Start 16: Performance.Benchmarks ...........   Passed
      Start 17: Performance.Scaling ..............   Passed

100% tests passed, 0 tests failed out of 17
```

## Test Development Guidelines

### Adding New Tests

1. Place test files in the appropriate category directory
2. Include in the corresponding source list in `tests/CMakeLists.txt`
3. Use consistent naming: `test_<component>.cpp`
4. Register test cases with meaningful names

### Test Structure

```cpp
#include <gtest/gtest.h>
#include "YourHeader.hpp"

class YourTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test fixtures
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

TEST_F(YourTest, TestName) {
    // Test implementation
    EXPECT_TRUE(condition);
}
```

### Best Practices

1. **Unit Tests**: Test single classes/functions in isolation
2. **Functional Tests**: Test feature workflows without full simulations
3. **Physics Tests**: Verify mathematical correctness with known solutions
4. **Integration Tests**: Test complete pipelines with realistic scenarios
5. **Performance Tests**: Include timing assertions with reasonable thresholds

## Dependencies

- **Google Test**: Testing framework
- **PETSc**: Parallel scientific computing (required for MPI initialization)
- **MPI**: Message passing for parallel tests
