# High-Fidelity Modules Testing Guide

Comprehensive guide for testing the high-fidelity physics and numerical modules.

---

## Test Organization

### Test Categories

| Category | Purpose | Files |
|----------|---------|-------|
| **Unit Tests** | Test individual class methods | `tests/unit/test_high_fidelity_*.cpp` |
| **MMS Tests** | Verify numerical convergence | `tests/physics/test_mms_high_fidelity.cpp` |
| **Coupling Tests** | Verify physics consistency | `tests/physics/test_coupling_consistency.cpp` |
| **Integration Tests** | End-to-end scenarios | `tests/integration/test_high_fidelity_integration.cpp` |

### Test File Structure

```
tests/
├── unit/
│   ├── test_high_fidelity_numerics.cpp      # SUPG, time integration, p-adapt
│   ├── test_high_fidelity_fluid_flow.cpp    # Non-Darcy, dual porosity, etc.
│   ├── test_high_fidelity_geomechanics.cpp  # Finite strain, creep, damage
│   └── test_advanced_coupled_physics.cpp    # Biot, THM, THMC, multiscale
├── physics/
│   ├── test_mms_high_fidelity.cpp           # MMS verification tests
│   └── test_coupling_consistency.cpp        # Cross-module consistency
└── integration/
    └── test_high_fidelity_integration.cpp   # Full-scale scenarios
```

---

## Running Tests

### Quick Reference

```bash
# Build tests
cd build
cmake .. -DBUILD_TESTS=ON
make

# Run all high-fidelity tests
make test_high_fidelity_all

# Run by category
make test_high_fidelity_unit
make test_high_fidelity_mms
make test_high_fidelity_integration

# Run with CTest
ctest -R "HighFidelity" --output-on-failure
ctest -R "THM" --output-on-failure
ctest -R "Coupling" --output-on-failure

# Run specific test executable with filter
./run_unit_tests --gtest_filter="*NonDarcy*"
./run_physics_tests --gtest_filter="*BiotCoupling*"
```

### Test Labels

| Label | Description |
|-------|-------------|
| `unit` | Unit tests |
| `quick` | Fast tests (< 1 min) |
| `physics` | Physics verification |
| `mms` | Method of Manufactured Solutions |
| `integration` | End-to-end tests |
| `slow` | Long-running tests |

```bash
# Run by label
ctest -L "unit" --output-on-failure
ctest -L "mms" --output-on-failure
ctest -L "integration" --output-on-failure
```

---

## Unit Test Coverage

### HighFidelityNumerics Tests

| Test Suite | Tests | Description |
|------------|-------|-------------|
| `SUPGStabilizationTest` | 10 | SUPG τ calculation, Peclet number |
| `PETScTimeIntegrationTest` | 8 | Time scheme orders, stability |
| `PAdaptivityTest` | 6 | Order change decisions |
| `PreconditionerTest` | 2 | Configuration |
| `ConfigValidationTest` | 4 | Parameter validation |
| `FactoryTest` | 3 | Factory function creation |

### HighFidelityFluidFlow Tests

| Test Suite | Tests | Description |
|------------|-------|-------------|
| `NonDarcyFlowTest` | 9 | Forchheimer, Klinkenberg, Reynolds |
| `DualPorosityTest` | 8 | Transfer rates, storativity |
| `NonIsothermalFlowTest` | 7 | Thermal properties, viscosity |
| `DynamicRelPermTest` | 9 | Capillary number, desaturation |
| `MiscibleFlowTest` | 10 | MMP, mixing, dispersion |
| `FluidFlowConfigTest` | 4 | Configuration validation |

### HighFidelityGeomechanics Tests

| Test Suite | Tests | Description |
|------------|-------|-------------|
| `FiniteStrainTest` | 10 | Neo-Hookean, Mooney-Rivlin, energy |
| `FiniteStrainPlasticityTest` | 6 | Yield function, hardening |
| `ViscoplasticityTest` | 8 | Perzyna, Duvaut-Lions |
| `CreepModelTest` | 8 | Power law, Burgers |
| `HypoplasticityTest` | 7 | Stress rate, void ratio limits |
| `GradientDamageTest` | 9 | Damage evolution, regularization |

### AdvancedCoupledPhysics Tests

| Test Suite | Tests | Description |
|------------|-------|-------------|
| `FullBiotDynamicsTest` | 14 | Wave velocities, Skempton, Gassmann |
| `THMCouplingTest` | 9 | Thermal stress, pressurization |
| `THMCCouplingTest` | 8 | Chemical strain, reactions |
| `UnsaturatedFlowTest` | 10 | Retention, effective stress |
| `MultiScaleCouplingTest` | 8 | Bounds, homogenization |
| `CouplingConsistencyTest` | 5 | Cross-module checks |

---

## MMS Test Details

### Purpose

Method of Manufactured Solutions (MMS) tests verify that numerical implementations converge at the expected rate.

### MMS Methodology

1. **Choose analytical solution**: Pick a smooth function
2. **Compute source term**: Substitute into PDE to find required source
3. **Solve numerically**: With manufactured source
4. **Compare**: Numerical vs analytical solution
5. **Verify convergence**: Error decreases with mesh refinement

### MMS Tests Included

#### Non-Darcy Flow MMS
- Forchheimer correction verification
- Klinkenberg source term consistency

#### Dual Porosity MMS
- Transfer rate consistency
- Storativity ratio physics

#### THM Coupling MMS
- Thermal strain consistency
- Effective stress verification
- Thermal diffusivity consistency

#### Biot Poroelasticity MMS
- Terzaghi consolidation (analytical solution)
- Wave velocity bounds verification

#### Viscoplasticity MMS
- Overstress formulation verification
- Strain rate scaling

#### Gradient Damage MMS
- Damage evolution bounds
- Internal length effect

#### Finite Strain MMS
- Neo-Hookean uniaxial verification
- Energy consistency

---

## Coupling Consistency Tests

### Purpose

Verify that coupled physics modules give consistent results.

### Test Categories

#### Biot Coupling Consistency
- Gassmann relation: `K_u = K_d + α²M`
- Skempton-Biot: `B = αM/K_u`
- Effective stress principle
- Wave velocity relations

#### THM Coupling Consistency
- Thermal strain symmetry
- Stress-strain consistency
- Reduces to Biot when isothermal

#### THMC Coupling Consistency
- Chemical strain linearity
- Reaction heat conservation
- Arrhenius rate at reference T

#### Unsaturated Coupling Consistency
- Saturation bounds
- Reduces to saturated at zero suction
- Bishop parameter bounds

#### Multi-Scale Consistency
- Voigt-Reuss bounds
- Hashin-Shtrikman tighter than VR
- Mori-Tanaka within bounds

---

## Integration Test Scenarios

### Reservoir Simulation

| Test | Scenario | Key Checks |
|------|----------|------------|
| `NonDarcyHighVelocityFlow` | Wellbore region | Correction factors |
| `DualPorosityExchange` | NFR production | Transfer rates |
| `MiscibleCO2Injection` | CO2 EOR | Miscibility, viscosity |
| `ThermalRecoveryProcess` | Steam injection | Temperature effects |

### Geomechanics

| Test | Scenario | Key Checks |
|------|----------|------------|
| `FiniteStrainRubberLike` | Soft rock | Stress-strain, energy |
| `ViscoplasticDeformation` | Salt creep | Time-dependent strain |
| `CreepUnderConstantLoad` | Long-term | Strain accumulation |
| `GradientDamageLocalization` | Fracturing | Regularization |

### Coupled Physics

| Test | Scenario | Key Checks |
|------|----------|------------|
| `FullBiotDynamicResponse` | Seismic | Wave velocities |
| `THMGeothermalReservoir` | EGS | Thermal stress |
| `THMCReactiveTransport` | CO2 storage | Chemical effects |
| `UnsaturatedSoilMechanics` | Near-surface | Bishop stress |
| `MultiScaleHomogenization` | Heterogeneous | Bounds |

### End-to-End

| Test | Application | Modules Used |
|------|-------------|--------------|
| `CO2SequestrationScenario` | CCS | Miscible + Geomech + THMC |
| `GeothermalEGSScenario` | EGS | NonDarcy + THM + Biot |

---

## Adding New Tests

### Unit Test Template

```cpp
#include <gtest/gtest.h>
#include "HighFidelityModule.hpp"

class NewModelTest : public ::testing::Test {
protected:
    NewModel model;
    
    void SetUp() override {
        NewModel::Parameters params;
        params.parameter1 = value1;
        model.setParameters(params);
    }
};

TEST_F(NewModelTest, BasicFunctionality) {
    double result = model.compute(input);
    EXPECT_NEAR(result, expected, tolerance);
}

TEST_F(NewModelTest, EdgeCase) {
    // Test edge cases
    EXPECT_GT(model.compute(0.0), 0.0);
}

TEST_F(NewModelTest, Configure) {
    std::map<std::string, std::string> config;
    config["param"] = "value";
    
    NewModel model2;
    model2.configure(config);
    
    SUCCEED();
}
```

### MMS Test Template

```cpp
class MMSNewModelTest : public ::testing::Test {
protected:
    // Analytical solution
    double analyticalSolution(double x, double t) {
        return std::sin(M_PI * x) * std::exp(-t);
    }
    
    // Compute L2 error
    double computeError(int n_cells) {
        // Implementation
    }
};

TEST_F(MMSNewModelTest, ConvergenceRate) {
    std::vector<int> mesh_sizes = {10, 20, 40, 80};
    std::vector<double> errors;
    
    for (int n : mesh_sizes) {
        errors.push_back(computeError(n));
    }
    
    // Verify convergence rate
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_GT(rate, expected_order - 0.5);
    }
}
```

### Registering Tests in CMake

```cmake
# Add to appropriate source list
set(UNIT_TEST_SOURCES
    ...
    unit/test_new_model.cpp
)

# Register with CTest
add_test(NAME Unit.NewModel
    COMMAND run_unit_tests --gtest_filter=*NewModel*.*)

# Set properties
set_tests_properties(Unit.NewModel
    PROPERTIES
        TIMEOUT 60
        LABELS "unit;quick"
)
```

---

## Test Best Practices

### General Guidelines

1. **One concept per test**: Each test should verify one thing
2. **Clear naming**: `TestSuite.TestName` describes what's tested
3. **Independence**: Tests should not depend on order
4. **Fast**: Unit tests should complete in milliseconds
5. **Deterministic**: Same result every run

### Numerical Testing

1. **Use tolerances**: `EXPECT_NEAR` instead of `EXPECT_EQ`
2. **Check bounds**: Verify physical constraints
3. **Test edge cases**: Zero, negative, extreme values
4. **Verify conservation**: Mass, energy, momentum

### Coupling Tests

1. **Check limiting cases**: Module should reduce to simpler model
2. **Verify consistency**: Same physics from different modules
3. **Test both directions**: A→B and B→A coupling
4. **Conservation**: Coupled system conserves appropriately

---

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| Test timeout | Increase `TIMEOUT` property |
| Memory errors | Run with Valgrind |
| MPI issues | Check `PETSC_COMM_WORLD` initialization |
| Numerical precision | Adjust tolerance in `EXPECT_NEAR` |

### Debugging Tests

```bash
# Run with verbose output
./run_unit_tests --gtest_filter="*FailingTest*" --gtest_print_time=1

# Run under GDB
gdb --args ./run_unit_tests --gtest_filter="*FailingTest*"

# Run with Valgrind
valgrind --leak-check=full ./run_unit_tests --gtest_filter="*FailingTest*"
```

### Test Discovery

```bash
# List all tests
./run_unit_tests --gtest_list_tests

# Filter syntax
./run_unit_tests --gtest_filter="Suite.*"           # All in suite
./run_unit_tests --gtest_filter="*.TestName"        # Specific test
./run_unit_tests --gtest_filter="*Pattern*"         # Pattern match
./run_unit_tests --gtest_filter="-Slow*"            # Exclude pattern
```
