# FSRM Test Suite

Comprehensive test suite organized by test type for the FSRM simulator.

## Test Organization

Tests are organized into categories by purpose and execution time:

| Category | Purpose | Typical Runtime | CI Pipeline |
|----------|---------|-----------------|-------------|
| **Unit** | Test individual components in isolation | < 1 min | Always |
| **Functional** | Test feature integration | 1-5 min | Always |
| **Physics** | Verify numerical accuracy (MMS) | 2-10 min | Always |
| **Integration** | End-to-end workflows | 5-20 min | On merge |
| **Performance** | Benchmarking and scaling | 10-30 min | Main only |

## Directory Structure

```
tests/
├── CMakeLists.txt              # Test build configuration
├── README.md                   # This file
├── test_main_gtest.cpp         # Main entry point for Google Test
├── data/                       # Test data files
│   └── test_simple.config
├── unit/                       # Unit tests (fast, isolated)
│   ├── test_config_reader.cpp  # ConfigReader parsing
│   ├── test_physics_kernels.cpp # Physics kernel classes
│   ├── test_well_model.cpp     # Well model operations
│   ├── test_fracture_model.cpp # Fracture/fault models
│   ├── test_eclipse_io.cpp     # Eclipse I/O
│   ├── test_fluid_model.cpp    # Generic fluid models
│   ├── test_material_model.cpp # Generic material models
│   ├── test_high_fidelity_numerics.cpp    # SUPG, time integration, p-adapt
│   ├── test_high_fidelity_fluid_flow.cpp  # Non-Darcy, dual porosity
│   ├── test_high_fidelity_geomechanics.cpp # Finite strain, creep, damage
│   └── test_advanced_coupled_physics.cpp  # Biot, THM, THMC, multiscale
├── functional/                 # Functional tests
│   ├── test_simulator_init.cpp
│   ├── test_solver_convergence.cpp
│   ├── test_boundary_conditions.cpp
│   └── test_well_operations.cpp
├── physics/                    # Physics/MMS tests
│   ├── test_mms_diffusion.cpp
│   ├── test_mms_elasticity.cpp
│   ├── test_mms_wave_propagation.cpp
│   ├── test_analytical_solutions.cpp
│   ├── test_mms_high_fidelity.cpp     # MMS for high-fidelity modules
│   └── test_coupling_consistency.cpp  # Cross-module coupling tests
├── integration/                # Integration tests
│   ├── test_full_simulation.cpp
│   ├── test_restart.cpp
│   └── test_high_fidelity_integration.cpp # End-to-end high-fidelity
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
- **FluidModel**: Generic fluid models (SinglePhase, BlackOil, Compositional, Brine, CO2)
- **MaterialModel**: Generic material models (LinearElastic, Poroelastic, Viscoelastic, Elastoplastic, Anisotropic)

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

### High-Fidelity Module Tests
Comprehensive tests for advanced physics and numerical methods:

#### Unit Tests
- **HighFidelityNumerics**: SUPG stabilization, PETSc time integration, p-adaptivity
- **HighFidelityFluidFlow**: Non-Darcy, dual/triple porosity, dynamic rel perm, miscible flow
- **HighFidelityGeomechanics**: Finite strain, viscoplasticity, creep, hypoplasticity, gradient damage
- **AdvancedCoupledPhysics**: Full Biot dynamics, THM, THMC, unsaturated flow, multi-scale

#### MMS Tests
- Non-Darcy flow convergence
- Dual porosity transfer consistency
- THM coupling verification
- Biot poroelasticity (Terzaghi analytical)
- Viscoplasticity overstress formulation
- Gradient damage regularization

#### Coupling Consistency Tests
- Gassmann relation verification
- Skempton-Biot consistency
- THM reduces to Biot when isothermal
- THMC reduces to THM without chemistry
- Voigt-Reuss-Hashin-Shtrikman bounds

#### Integration Tests
- Reservoir simulation scenarios (non-Darcy, dual porosity, miscible)
- Geomechanics scenarios (finite strain, creep, damage)
- Coupled physics scenarios (Biot dynamics, THM, THMC)
- End-to-end scenarios (CO2 sequestration, geothermal EGS)

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

### High-Fidelity Tests
```bash
# All high-fidelity tests
make test_high_fidelity_all

# By category
make test_high_fidelity_unit
make test_high_fidelity_mms
make test_high_fidelity_integration

# With CTest
ctest -R "HighFidelity" --output-on-failure
ctest -R "THM" --output-on-failure
ctest -R "Biot" --output-on-failure
ctest -R "Coupling" --output-on-failure
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
