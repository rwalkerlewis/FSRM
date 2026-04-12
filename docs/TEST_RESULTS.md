# FSRM Test Results

Generated from `ctest --output-on-failure` on the `local_fix` branch.

**68/68 tests pass. 0 failures. 0 skips.**

## Test Summary by Label

| Label | Tests | Time |
|-------|-------|------|
| unit | 30 | 13.10 sec |
| physics_validation | 18 | 21.14 sec |
| integration | 12 | 16.49 sec |
| functional | 4 | 1.48 sec |
| performance | 3 | 1.89 sec |
| experimental | 1 | 3.03 sec |

## Complete Test Matrix

### Unit Tests (30 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 1 | Unit.ConfigReader | PASS | 2.18s |
| 2 | Unit.PhysicsModuleRegistry | PASS | 0.37s |
| 3 | Unit.SimulationConfig | PASS | 0.38s |
| 4 | Unit.PhysicsKernelBase | PASS | 0.37s |
| 5 | Unit.SinglePhaseFlow | PASS | 0.38s |
| 6 | Unit.Geomechanics | PASS | 0.37s |
| 7 | Unit.Thermal | PASS | 0.38s |
| 8 | Unit.Elastodynamics | PASS | 0.41s |
| 9 | Unit.Poroelastodynamics | PASS | 0.39s |
| 10 | Unit.BoundaryConditions | PASS | 0.37s |
| 11 | Unit.WellModel | PASS | 0.36s |
| 12 | Unit.FractureModel | PASS | 0.38s |
| 13 | Unit.MaterialModel | PASS | 0.37s |
| 14 | Unit.SeismicSource | PASS | 0.36s |
| 15 | Unit.ExplosionSource | PASS | 0.37s |
| 16 | Unit.EclipseIO | PASS | 0.37s |
| 17 | Unit.MLRegistry | PASS | 0.37s |
| 18 | Unit.UnitSystem | PASS | 0.40s |
| 19 | Unit.CoordinateSystem | PASS | 0.38s |
| 20 | Unit.CO2Properties | PASS | 0.36s |
| 21 | Unit.FlashCalculation | PASS | 0.36s |
| 22 | Unit.ResidualTrapping | PASS | 0.39s |
| 23 | Unit.CaprockIntegrity | PASS | 0.38s |
| 24 | Unit.MultiphaseFlow | PASS | 0.40s |
| 25 | Unit.FaultMechanics | PASS | 0.40s |
| 26 | Unit.FaultCohesiveDyn | PASS | 0.38s |
| 27 | Unit.CoulombStressTransfer | PASS | 0.38s |
| 28 | Unit.FaultMeshManager | PASS | 0.40s |
| 29 | Unit.CohesiveFaultKernel | PASS | 0.38s |
| 30 | Unit.ElasticityAux | PASS | 0.38s |

### Functional Tests (4 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 31 | Functional.SimulatorInit | PASS | 0.38s |
| 32 | Functional.ModuleLifecycle | PASS | 0.37s |
| 33 | Functional.SolverConvergence | PASS | 0.37s |
| 34 | Functional.WellOperations | PASS | 0.36s |

### Physics Validation Tests (18 tests)

| # | Test Name | Status | Time | Description |
|---|-----------|--------|------|-------------|
| 35 | Physics.MMS.Diffusion | PASS | 0.38s | Method of manufactured solutions for diffusion |
| 36 | Physics.MMS.Elasticity | PASS | 0.39s | MMS for elasticity |
| 37 | Physics.MMS.WavePropagation | PASS | 0.38s | MMS for wave propagation |
| 38 | Physics.AnalyticalSolutions | PASS | 0.36s | Analytical benchmark comparisons |
| 39 | Physics.Thermoporoelastic | PASS | 0.39s | Thermoporoelastic coupling |
| 40 | Physics.SCEC.TPV5 | PASS | 0.40s | SCEC TPV5 fault infrastructure (3 real tests, no skips) |
| 41 | Physics.ElastostaticsPatch | PASS | 0.41s | Elastostatics patch test |
| 42 | Physics.TerzaghiConsolidation | PASS | 4.54s | Terzaghi 1D consolidation vs analytical |
| 43 | Physics.AbsorbingBC | PASS | 5.61s | Clayton-Engquist absorbing boundary conditions |
| 44 | Physics.GravityLithostatic | PASS | 0.51s | Gravity and lithostatic stress |
| 45 | Physics.LithostaticStress | PASS | 0.97s | Lithostatic stress verification |
| 46 | Physics.LambsProblem | PASS | 1.55s | Lamb's problem (point force on halfspace) |
| 47 | Physics.GarvinsProblem | PASS | 1.58s | Garvin's problem (buried explosion) |
| 48 | Physics.MuellerMurphy | PASS | 0.38s | Mueller-Murphy seismic source (9 test cases) |
| 49 | Physics.AtmosphericExplosion | PASS | 0.38s | Sedov-Taylor, Brode, EMP, overpressure (8 test cases) |
| 50 | Physics.NearFieldExplosion | PASS | 0.38s | Cavity, damage, spall, solver (13 test cases) |
| 51 | Physics.DruckerPrager | PASS | 0.41s | DP/VM/MC yield + hydrostatic check (8 test cases) |
| 52 | Physics.MomentTensorSource | PASS | 2.13s | FEM source injection, solution norm (3 test cases) |

### Integration Tests (12 tests)

| # | Test Name | Status | Time | Description |
|---|-----------|--------|------|-------------|
| 53 | Integration.FullSimulation | PASS | 0.39s | Full simulation lifecycle |
| 54 | Integration.Restart | PASS | 3.03s | Checkpoint/restart |
| 55 | Integration.CoupledPhysics | PASS | 0.40s | Coupled physics modules |
| 56 | Integration.InjectionRuptureChain | PASS | 0.40s | Injection-to-rupture chain |
| 57 | Integration.LockedFault | PASS | 0.38s | Locked fault behavior |
| 58 | Integration.LayeredElastostatics | PASS | 0.40s | Layered material elastostatics |
| 59 | Integration.ExplosionSeismogram | PASS | 8.77s | Explosion source to seismogram pipeline |
| 60 | Integration.PrescribedSlip | PASS | 0.44s | Prescribed slip on cohesive fault |
| 61 | Integration.OutputFile | PASS | 0.63s | Output file generation |
| 62 | Integration.GmshImport | PASS | 0.64s | Gmsh mesh import |
| 63 | Integration.DerivedFields | PASS | 0.63s | Derived field computation (stress, strain, CFS) |
| 64 | Integration.DPRK2017Comparison | PASS | 0.38s | DPRK 2017 synthetic vs observed mb |

### Performance Tests (3 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 65 | Performance.Benchmarks | PASS | 0.41s |
| 66 | Performance.Scaling | PASS | 1.10s |
| 67 | Performance.Memory | PASS | 0.38s |

### Experimental Tests (1 test)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 68 | Experimental.NeuralStubs | PASS | 3.03s |

## Known Limitations Documented in Tests

1. **Drucker-Prager return mapping** (Physics.DruckerPrager): Yield function evaluation works correctly for Drucker-Prager, von Mises, and Mohr-Coulomb. Return mapping algorithm does not produce nonzero plastic strain. Tests document this gap. Not wired into PETSc FEM pipeline.

2. **SCEC TPV5 full benchmark** (Physics.SCEC.TPV5): Infrastructure (parameters, CohesiveFaultKernel, FaultMeshManager, friction parameter setup) verified with 3 real tests. Full dynamic rupture solve is a work in progress.

3. **Mohr-Coulomb principal stress computation** (Physics.DruckerPrager): The computePrincipalStresses function returns NaN for pure hydrostatic stress (q=0) due to division by zero in the Cardano formula. Tests use stress states with nonzero deviatoric component to avoid this.

## Test Environment

- **Docker image**: fsrm-ci:local (Dockerfile.ci)
- **PETSc**: 3.22.2 with --with-debugging=0
- **Compiler**: g++ (C++17)
- **Build**: CMake Release, ENABLE_TESTING=ON, ENABLE_CUDA=OFF, BUILD_EXAMPLES=ON
- **Branch**: local_fix
