# FSRM Test Suite

This directory contains comprehensive tests for the Finite-volume Simulator for Reservoir Mechanics (FSRM).

## Test Organization

### Unit Tests
Located in `test_main.cpp`, `test_physics_kernels.cpp`, and `test_fracture_model.cpp`

**Purpose**: Test individual components in isolation
- `test_main.cpp`: Basic testing framework and simple unit tests
- `test_physics_kernels.cpp`: Physics kernel implementations
  - Single phase flow residual/Jacobian
  - Poroelasticity coupling terms
  - Elastodynamics wave equations
  - Thermal flow
- `test_fracture_model.cpp`: Fracture mechanics
  - LEFM stress intensity factors
  - Hydraulic fracture growth (KGD model)
  - Fracture network connectivity
  - Proppant transport

### Convergence Tests
Located in `test_convergence_mms.cpp`

**Purpose**: Verify numerical accuracy using Method of Manufactured Solutions (MMS)

**Test Categories**:
- Spatial convergence (expected: O(h²) for linear elements)
- Temporal convergence (O(dt) for BE, O(dt²) for BDF2)
- MMS tests with known analytical solutions

**Test Cases**:
1. **Single Phase Flow MMS**
   - Linear solution: `u = x + 2y + 3z + t`
   - Quadratic solution: `u = x² + y² + z²`
   - Trigonometric: `u = sin(πx)sin(πy)exp(-t)`

2. **Poroelasticity Analytical**
   - Mandel-Cryer problem
   - Terzaghi consolidation

3. **Elastodynamics MMS**
   - Plane wave propagation
   - Rayleigh surface waves

### Integration Tests
Located in `test_integration.cpp`

**Purpose**: Test complete simulation workflows

**Test Scenarios**:
1. Single phase flow in simple domain
2. Single phase flow with injection/production wells
3. Poroelasticity consolidation
4. Elastodynamic wave propagation
5. Hydraulic fracturing simulation
6. Induced seismicity
7. Performance/scalability benchmarks
8. Restart/checkpoint functionality
9. SPE1 regression test

## Building Tests

### Prerequisites
- CMake >= 3.15
- GoogleTest (optional but recommended)
- PETSc >= 3.15
- MPI

### Build Commands

```bash
# From project root
mkdir -p build
cd build
cmake .. -DENABLE_TESTING=ON
make

# Build tests only
make run_tests
```

## Running Tests

### Run All Tests
```bash
cd build
ctest
# or
make test
```

### Run Specific Test Categories

```bash
# Quick unit tests only (~1 minute)
make test_quick
# or
ctest -R "UnitTest|EclipseIO|PhysicsKernel|WellModel"

# Convergence tests (~5 minutes)
make test_convergence
# or
ctest -R "MMS|Convergence"

# Integration tests (~10 minutes)
make test_integration
# or
ctest -R "Integration"

# All tests
make test_all
```

### Run Individual Test Groups

```bash
# Physics kernel tests
ctest -R PhysicsKernelTestFixture

# Fracture model tests
ctest -R FractureModelTest

# Single phase flow tests (all)
ctest -R SinglePhaseFlow

# Poroelasticity tests
ctest -R Poroelasticity

# Elastodynamics tests
ctest -R Elastodynamics

# Hydraulic fracture tests
ctest -R HydraulicFracture
```

### Run with Verbose Output

```bash
ctest -V
# or
ctest --verbose
```

### Run in Parallel (MPI)

```bash
# Run tests with MPI
mpirun -np 4 ./run_tests

# Or use CTest's parallel execution
ctest -j 4
```

## Test Output

### Test Reports
Tests generate reports in the build directory:
- `test_report.txt`: Summary of all unit tests
- Convergence plots: `*_convergence.png`
- Comparison plots: `*_comparison.png`

### Interpreting Results

**Success Criteria**:
- All unit tests pass
- Convergence rates meet theoretical predictions:
  - Spatial: O(h²) for linear FEM
  - Temporal: O(dt) for BE, O(dt²) for BDF2
- Integration tests complete without errors
- Mass/energy conservation errors < 1%

**Common Failures**:
1. **Convergence rate below expected**: 
   - Check mesh quality
   - Verify boundary conditions
   - Review source term computation

2. **Integration test timeouts**:
   - Reduce problem size
   - Adjust solver tolerances
   - Check for infinite loops

3. **Mass balance errors**:
   - Verify well rates
   - Check boundary conditions
   - Review time step size

## Adding New Tests

### Unit Test Example

```cpp
TEST_F(MyTestFixture, NewFeatureTest) {
    // Setup
    MyClass obj;
    obj.setParameter(value);
    
    // Execute
    double result = obj.compute();
    
    // Verify
    EXPECT_NEAR(result, expected, tolerance);
    EXPECT_GT(result, 0.0);
}
```

### MMS Test Example

```cpp
TEST_F(MMSConvergenceTest, NewPhysicsMMS) {
    MMSTest mms_test("TestName");
    
    auto u_exact = [](double x, double y, double z, double t) {
        return x*x + y*y;  // Manufactured solution
    };
    
    mms_test.setManufacturedSolution(u_exact);
    
    std::vector<int> mesh_sizes = {8, 16, 32, 64};
    auto result = mms_test.runConvergenceTest(sim, mesh_sizes);
    
    EXPECT_TRUE(result.passed);
    EXPECT_GT(result.achieved_rate, 1.9);
}
```

### Integration Test Example

```cpp
TEST_F(IntegrationTest, NewScenario) {
    // Create config file
    std::ofstream config("test_config.config");
    config << "physics_type = FLUID_FLOW\n";
    config << "grid_nx = 10\n";
    // ... more parameters
    config.close();
    
    // Run simulation
    ConfigReader reader("test_config.config");
    Simulator sim(PETSC_COMM_WORLD);
    sim.initialize(reader);
    sim.solve();
    
    // Verify results
    EXPECT_TRUE(sim.converged());
}
```

## Continuous Integration

Tests are automatically run on:
- Every commit to main branch
- Every pull request
- Nightly builds

CI Configuration: `.github/workflows/tests.yml`

## Performance Benchmarks

Benchmark tests measure:
- DOFs/second throughput
- Memory usage
- Parallel scalability
- Strong/weak scaling

Results are tracked over time to detect performance regressions.

## Test Coverage

Current coverage (approximate):
- Physics kernels: 85%
- Fracture models: 75%
- Well models: 70%
- I/O operations: 60%
- Utilities: 80%

Generate coverage report:
```bash
cmake .. -DENABLE_COVERAGE=ON
make coverage
```

## Debugging Failed Tests

### GDB with GoogleTest

```bash
gdb --args ./run_tests --gtest_filter=MyFailingTest
(gdb) run
(gdb) backtrace
```

### Valgrind for Memory Errors

```bash
valgrind --leak-check=full ./run_tests --gtest_filter=MyTest
```

### PETSc Debugging

Set environment variables:
```bash
export PETSC_OPTIONS="-on_error_attach_debugger -log_view"
./run_tests
```

## References

1. **Method of Manufactured Solutions**
   - Roache, P.J. (2002). "Code Verification by the Method of Manufactured Solutions"

2. **Convergence Testing**
   - Oberkampf & Roy (2010). "Verification and Validation in Scientific Computing"

3. **Benchmark Problems**
   - SPE Comparative Solution Projects
   - Mandel-Cryer problem (poroelasticity)
   - Terzaghi consolidation

4. **GoogleTest Documentation**
   - https://google.github.io/googletest/

## Contact

For test-related questions, contact the development team or open an issue on GitHub.
