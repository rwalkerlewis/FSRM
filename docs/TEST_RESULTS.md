# FSRM Test Results

Generated from `ctest --output-on-failure` on the `feature/hydrofrac-fem` branch.

**84/84 tests pass. 0 failures. 0 skips.**

Recent additions in this session:
- `Physics.PressurizedFracture` (Phase 1 hydrofracture validation)
- `Integration.DynamicRuptureBasic.AbsorbingCoexist` (fault + absorbing callback coexistence)
- `Integration.FaultAbsorbingCoexist` (region-specific PetscDS setup)
- `Integration.FractureFlow` (Phase 2 lubrication callback validation)
- `Integration.CoupledHydrofrac` (Phase 3 coupling utility checks)
- `Integration.FracturePropagation` (Phase 4 propagation criterion checks)
- `Integration.StressShadowing` (Phase 5 multi-cluster stress shadow)
- `Integration.InducedSeismicity` (Phase 7 moment tensor and Mw)
- `Integration.ProppantTransport` (Phase 6 settling, bridging, mass balance)
- `Integration.LeakoffCoupling` (Phase 8 Carter leak-off)
- `Integration.ProductionForecast` (Phase 9 Arps decline curves)

## Test Summary by Label

| Label | Tests | Time |
|-------|-------|------|
| unit | 30 | ~14 sec |
| physics_validation | 20 | ~22 sec |
| integration | 26 | ~25 sec |
| functional | 4 | ~1.6 sec |
| performance | 3 | ~1.9 sec |
| experimental | 1 | ~0.4 sec |

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

### Physics Validation Tests (20 tests)

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
| 52 | Physics.Elastoplasticity | PASS | 0.37s | Elastoplasticity return mapping |
| 53 | Physics.MomentTensorSource | PASS | 2.13s | FEM source injection, solution norm (3 test cases) |
| 54 | Physics.PressurizedFracture | PASS | 0.36s | Sneddon aperture validation (Phase 1) |

### Integration Tests (26 tests)

| # | Test Name | Status | Time | Description |
|---|-----------|--------|------|-------------|
| 55 | Integration.FullSimulation | PASS | 0.39s | Full simulation lifecycle |
| 56 | Integration.Restart | PASS | 3.03s | Checkpoint/restart |
| 57 | Integration.CoupledPhysics | PASS | 0.40s | Coupled physics modules |
| 58 | Integration.InjectionPressure | PASS | 0.79s | Injection pressure buildup |
| 59 | Integration.InjectionRuptureChain | PASS | 0.40s | Injection-to-rupture chain |
| 60 | Integration.LockedFault | PASS | 0.38s | Locked fault behavior |
| 61 | Integration.LayeredElastostatics | PASS | 0.40s | Layered material elastostatics |
| 62 | Integration.ExplosionSeismogram | PASS | 8.77s | Explosion source to seismogram pipeline |
| 63 | Integration.PrescribedSlip | PASS | 0.44s | Prescribed slip on cohesive fault |
| 64 | Integration.OutputFile | PASS | 0.63s | Output file generation |
| 65 | Integration.GmshImport | PASS | 0.64s | Gmsh mesh import |
| 66 | Integration.DerivedFields | PASS | 0.63s | Derived field computation (stress, strain, CFS) |
| 67 | Integration.DPRK2017Comparison | PASS | 0.38s | DPRK 2017 synthetic vs observed mb |
| 68 | Integration.DynamicRuptureBasic.MeshSplitting | PASS | 0.38s | Cohesive cell mesh splitting |
| 69 | Integration.DynamicRuptureBasic.LockedFault | PASS | 0.45s | Locked fault setup |
| 70 | Integration.DynamicRuptureBasic.SlippingFault | PASS | 3.04s | Slipping fault setup |
| 71 | Integration.DynamicRuptureBasic.AbsorbingCoexist | PASS | 0.39s | Fault + absorbing BC coexistence |
| 72 | Integration.FaultAbsorbingCoexist | PASS | 0.39s | Region-specific PetscDS setup |
| 73 | Integration.FractureFlow | PASS | 0.37s | Poiseuille lubrication validation (Phase 2) |
| 74 | Integration.CoupledHydrofrac | PASS | 0.36s | PKN width scaling (Phase 3) |
| 75 | Integration.FracturePropagation | PASS | 0.35s | Cohesive strength + opening criterion (Phase 4) |
| 76 | Integration.StressShadowing | PASS | 2.10s | Multi-cluster stress shadow (Phase 5, 6 tests) |
| 77 | Integration.InducedSeismicity | PASS | 0.39s | Moment tensor + Mw (Phase 7, 6 tests) |
| 78 | Integration.ProppantTransport | PASS | 0.38s | Stokes settling + bridging (Phase 6, 8 tests) |
| 79 | Integration.LeakoffCoupling | PASS | 0.35s | Carter leak-off sqrt(t) scaling (Phase 8, 7 tests) |
| 80 | Integration.ProductionForecast | PASS | 0.38s | Arps decline + PI (Phase 9, 9 tests) |

### Performance Tests (3 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 81 | Performance.Benchmarks | PASS | 0.41s |
| 82 | Performance.Scaling | PASS | 1.10s |
| 83 | Performance.Memory | PASS | 0.38s |

### Experimental Tests (1 test)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 84 | Experimental.NeuralStubs | PASS | 0.38s |

## Known Limitations Documented in Tests

1. **Drucker-Prager return mapping** (Physics.DruckerPrager): Yield function evaluation works correctly for Drucker-Prager, von Mises, and Mohr-Coulomb. Return mapping algorithm does not produce nonzero plastic strain. Tests document this gap. Not wired into PETSc FEM pipeline.

2. **SCEC TPV5 full benchmark** (Physics.SCEC.TPV5): Infrastructure (parameters, CohesiveFaultKernel, FaultMeshManager, friction parameter setup) verified with 3 real tests. Full dynamic rupture solve is a work in progress.

3. **Mohr-Coulomb principal stress computation** (Physics.DruckerPrager): The computePrincipalStresses function returns NaN for pure hydrostatic stress (q=0) due to division by zero in the Cardano formula. Tests use stress states with nonzero deviatoric component to avoid this.

## Test Environment

- **Docker image**: fsrm-ci:local (Dockerfile.ci)
- **PETSc**: 3.22.2 with --with-debugging=0
- **Compiler**: g++ (C++17)
- **Build**: CMake Release, ENABLE_TESTING=ON, ENABLE_CUDA=OFF, BUILD_EXAMPLES=ON
- **Branch**: feature/hydrofrac-fem
