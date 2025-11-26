# FSRM Test Suite

## Directory Structure

```
tests/
├── CMakeLists.txt              # Test build configuration
├── README.md                   # This file
├── test_main_gtest.cpp         # Main entry point for Google Test
├── data/                       # Test data files
│   └── test_simple.config
├── unit/                       # Unit tests
│   ├── test_config_reader.cpp
│   ├── test_physics_kernels.cpp
│   ├── test_well_model.cpp
│   ├── test_fracture_model.cpp
│   └── test_eclipse_io.cpp
├── functional/                 # Functional tests
│   ├── test_simulator_init.cpp
│   ├── test_solver_convergence.cpp
│   ├── test_boundary_conditions.cpp
│   └── test_well_operations.cpp
├── physics/                    # Physics/MMS tests
│   ├── test_mms_diffusion.cpp
│   ├── test_mms_elasticity.cpp
│   ├── test_mms_wave_propagation.cpp
│   └── test_analytical_solutions.cpp
├── integration/                # Integration tests
│   ├── test_full_simulation.cpp
│   └── test_restart.cpp
└── performance/                # Performance tests
    ├── test_benchmarks.cpp
    └── test_scaling.cpp
```

## Test Categories

### Unit Tests
Fast, isolated tests for individual classes and functions:
- **ConfigReader**: Configuration file parsing
- **PhysicsKernels**: Physics kernel creation and property setting
- **WellModel**: Well types, completions, calculations
- **FractureModel**: Fracture networks, faults, stress calculations
- **EclipseIO**: Eclipse format file I/O

### Functional Tests
Tests for specific features and module integration:
- **SimulatorInit**: Simulator configuration and initialization
- **SolverConvergence**: Numerical solver algorithms
- **BoundaryConditions**: Boundary condition handling
- **WellOperations**: Well index and productivity calculations

### Physics Tests (MMS)
Method of Manufactured Solutions for numerical verification:
- **Diffusion**: Pressure diffusion equation
- **Elasticity**: Linear elasticity and poroelasticity
- **WavePropagation**: Elastic wave equations
- **AnalyticalSolutions**: Theis, Terzaghi, Buckley-Leverett, etc.

### Integration Tests
End-to-end workflow tests:
- **FullSimulation**: Complete simulation pipeline
- **Restart**: Checkpoint and restart functionality

### Performance Tests
Benchmarking and parallel scaling:
- **Benchmarks**: Kernel timing and memory access
- **Scaling**: MPI parallel performance

## Building Tests

```bash
mkdir build && cd build
cmake .. -DENABLE_TESTING=ON
make -j$(nproc)
```

## Running Tests

### All Tests
```bash
ctest --output-on-failure
```

### By Category
```bash
ctest -L "unit"          # Unit tests
ctest -L "functional"    # Functional tests
ctest -L "physics"       # Physics/MMS tests
ctest -L "integration"   # Integration tests
ctest -L "performance"   # Performance tests
```

### Specific Test
```bash
ctest -R "Unit.ConfigReader"
```

### Verbose Output
```bash
ctest -V
```

## Writing New Tests

1. Create test file in appropriate category directory
2. Include in `CMakeLists.txt` source lists
3. Use Google Test framework:

```cpp
#include <gtest/gtest.h>
#include "YourComponent.hpp"

class YourTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank;
};

TEST_F(YourTest, TestName) {
    // Test code
    EXPECT_TRUE(condition);
}
```

## Test Dependencies

- **Google Test**: C++ testing framework
- **PETSc**: For MPI initialization
- **MPI**: Parallel execution support
