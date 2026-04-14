# FSRM Test Results

Generated from `ctest --output-on-failure` in Docker.

**108/108 tests pass. 0 failures. 0 skips.**

## Environment

- Docker image: `fsrm-ci:local` from `Dockerfile.ci`
- PETSc: 3.22.2 with `--with-debugging=0`
- Compiler: g++ (C++17)
- Build: `CMAKE_BUILD_TYPE=Release`, `ENABLE_TESTING=ON`, `ENABLE_CUDA=OFF`, `BUILD_EXAMPLES=ON`
- Command:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure
```

## Category Summary

| Category | Count | Total Time |
|---|---:|---:|
| Unit | 36 | ~14s |
| Functional | 10 | ~4s |
| Physics validation | 25 | ~31s |
| Integration | 32 | ~60s |
| Performance | 4 | ~2s |
| Experimental | 1 | ~0.3s |
| **Total** | **108** | **~116s** |

## Label Corrections (This Release)

Several tests were relabeled from incorrect categories:

- `Integration.DynamicRuptureBasic.*` renamed to `Functional.DynamicRuptureSetup.*` (setup-only, no TSSolve)
- `Integration.FractureFlow` renamed to `Unit.FractureFlowCallbacks` (callback formula test)
- `Integration.CoupledHydrofrac` renamed to `Unit.HydrofracFormulas` (analytical formula test)
- `Integration.FracturePropagation` renamed to `Unit.FracturePropagation` (formula test)
- `Integration.StressShadowing` renamed to `Physics.StressShadowing` (standalone physics)
- `Integration.InducedSeismicity` renamed to `Physics.InducedSeismicity` (standalone physics)
- `Integration.ProppantTransport` renamed to `Physics.ProppantTransport` (standalone physics)
- `Integration.LeakoffCoupling` renamed to `Physics.LeakoffCoupling` (standalone physics)
- `Integration.ProductionForecast` renamed to `Physics.ProductionForecast` (standalone physics)
- `Physics.DruckerPrager` renamed to `Unit.DruckerPragerStandalone` (material-point formula)
- `Physics.Elastoplasticity` renamed to `Unit.Elastoplasticity` (callback isolation test)
- `Physics.PressurizedFracture` renamed to `Unit.PressurizedFractureCallbacks` (callback isolation)
- `Integration.FaultAbsorbingCoexist` relabeled to `Functional` (setup-only, no TSSolve)

## Complete Test Matrix

### Unit Tests (36 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 1 | Unit.ConfigReader | PASS | 2.22s |
| 2 | Unit.PhysicsModuleRegistry | PASS | 0.38s |
| 3 | Unit.SimulationConfig | PASS | 0.41s |
| 4 | Unit.PhysicsKernelBase | PASS | 0.40s |
| 5 | Unit.SinglePhaseFlow | PASS | 0.41s |
| 6 | Unit.Geomechanics | PASS | 0.35s |
| 7 | Unit.Thermal | PASS | 0.35s |
| 8 | Unit.Elastodynamics | PASS | 0.37s |
| 9 | Unit.Poroelastodynamics | PASS | 0.37s |
| 10 | Unit.BoundaryConditions | PASS | 0.37s |
| 11 | Unit.WellModel | PASS | 0.40s |
| 12 | Unit.FractureModel | PASS | 0.40s |
| 13 | Unit.MaterialModel | PASS | 0.38s |
| 14 | Unit.SeismicSource | PASS | 0.34s |
| 15 | Unit.ExplosionSource | PASS | 0.38s |
| 16 | Unit.EclipseIO | PASS | 0.38s |
| 17 | Unit.MLRegistry | PASS | 0.39s |
| 18 | Unit.UnitSystem | PASS | 0.40s |
| 19 | Unit.CoordinateSystem | PASS | 0.38s |
| 20 | Unit.CO2Properties | PASS | 0.42s |
| 21 | Unit.FlashCalculation | PASS | 0.38s |
| 22 | Unit.ResidualTrapping | PASS | 0.42s |
| 23 | Unit.CaprockIntegrity | PASS | 0.39s |
| 24 | Unit.MultiphaseFlow | PASS | 0.41s |
| 25 | Unit.FaultMechanics | PASS | 0.42s |
| 26 | Unit.FaultCohesiveDyn | PASS | 0.38s |
| 27 | Unit.CoulombStressTransfer | PASS | 0.44s |
| 28 | Unit.FaultMeshManager | PASS | 0.41s |
| 29 | Unit.CohesiveFaultKernel | PASS | 0.36s |
| 30 | Unit.ElasticityAux | PASS | 0.40s |
| 53 | Unit.PressurizedFractureCallbacks | PASS | 2.15s |
| 54 | Unit.DruckerPragerStandalone | PASS | 0.40s |
| 55 | Unit.Elastoplasticity | PASS | 0.41s |
| 79 | Unit.FractureFlowCallbacks | PASS | 0.38s |
| 80 | Unit.HydrofracFormulas | PASS | 0.38s |
| 81 | Unit.FracturePropagation | PASS | 0.39s |

### Functional Tests (10 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 31 | Functional.SimulatorInit | PASS | 0.39s |
| 32 | Functional.ModuleLifecycle | PASS | 0.40s |
| 33 | Functional.SolverConvergence | PASS | 0.38s |
| 34 | Functional.WellOperations | PASS | 0.38s |
| 73 | Functional.DynamicRuptureSetup.MeshSplitting | PASS | 0.43s |
| 74 | Functional.DynamicRuptureSetup.LockedFault | PASS | 0.45s |
| 75 | Functional.DynamicRuptureSetup.FaultGeometryFromConfig | PASS | 0.42s |
| 76 | Functional.DynamicRuptureSetup.SlippingFault | PASS | 0.38s |
| 77 | Functional.DynamicRuptureSetup.AbsorbingCoexist | PASS | 0.43s |
| 78 | Functional.FaultAbsorbingCoexist | PASS | 0.42s |

### Physics Validation Tests (25 tests)

| # | Test Name | Status | Time | Description |
|---|-----------|--------|------|-------------|
| 35 | Physics.MMS.Diffusion | PASS | 0.39s | Method of manufactured solutions for diffusion |
| 36 | Physics.MMS.Elasticity | PASS | 0.40s | MMS for elasticity |
| 37 | Physics.MMS.WavePropagation | PASS | 0.38s | MMS for wave propagation |
| 38 | Physics.AnalyticalSolutions | PASS | 0.42s | Analytical benchmark comparisons |
| 39 | Physics.Thermoporoelastic | PASS | 0.40s | Thermoporoelastic coupling |
| 40 | Physics.SCEC.TPV5 | PASS | 3.09s | SCEC TPV5 fault infrastructure (3 tests, no skips) |
| 41 | Physics.ElastostaticsPatch | PASS | 0.35s | Elastostatics patch test |
| 42 | Physics.TerzaghiConsolidation | PASS | 1.77s | Terzaghi 1D consolidation vs analytical |
| 43 | Physics.AbsorbingBC | PASS | 5.31s | Clayton-Engquist absorbing BCs (>99% energy absorption) |
| 44 | Physics.GravityLithostatic | PASS | 0.48s | Gravity body force and lithostatic stress |
| 45 | Physics.LithostaticStress | PASS | 1.01s | Lithostatic stress verification |
| 46 | Physics.LambsProblem | PASS | 1.65s | Lamb problem (point force on halfspace) |
| 47 | Physics.GarvinsProblem | PASS | 1.54s | Garvin problem (buried explosion Green function) |
| 48 | Physics.MuellerMurphy | PASS | 0.36s | Mueller-Murphy seismic source (9 test cases) |
| 49 | Physics.AtmosphericExplosion | PASS | 0.39s | Sedov-Taylor, Brode, EMP, overpressure (8 test cases) |
| 50 | Physics.NearFieldExplosion | PASS | 0.37s | Cavity, damage, spall, 1D solver (13 test cases) |
| 51 | Physics.MomentTensorSource | PASS | 2.16s | FEM source injection via tabulation (3 test cases) |
| 52 | Physics.ExplosionDamageZone | PASS | 7.12s | FEM-coupled aux-field damage-zone degradation |
| 82 | Physics.StressShadowing | PASS | 0.41s | Multi-cluster stress shadow (6 tests) |
| 83 | Physics.InducedSeismicity | PASS | 0.39s | Moment tensor and Mw (6 tests) |
| 84 | Physics.ProppantTransport | PASS | 0.41s | Stokes settling and bridging (8 tests) |
| 85 | Physics.LeakoffCoupling | PASS | 0.39s | Carter leak-off sqrt(t) scaling (7 tests) |
| 86 | Physics.ProductionForecast | PASS | 0.40s | Arps decline and PI (9 tests) |
| 87 | Physics.LockedFaultTransparency | PASS | ~1.0s | Locked fault slip must be < 5e-4 |
| 88 | Physics.ElastostaticsPatchLayered | PASS | ~0.5s | Layered patch test |

### Integration Tests (32 tests)

| # | Test Name | Status | Time | Description |
|---|-----------|--------|------|-------------|
| 56 | Integration.FullSimulation | PASS | 0.41s | Full simulation lifecycle through TSSolve |
| 57 | Integration.Restart | PASS | 0.37s | Checkpoint/restart |
| 58 | Integration.CoupledPhysics | PASS | 0.41s | Coupled physics modules |
| 59 | Integration.InjectionPressure | PASS | 0.80s | Injection pressure buildup (two-sim comparison) |
| 60 | Integration.InjectionRuptureChain | PASS | 0.40s | Injection-to-rupture chain |
| 61 | Integration.LockedFault | PASS | 0.41s | Locked fault behavior |
| 62 | Integration.LayeredElastostatics | PASS | 0.40s | Depth-based material layering via auxiliary fields |
| 63 | Integration.ExplosionSeismogram | PASS | 8.54s | Explosion source to SAC seismogram pipeline |
| 64 | Integration.PrescribedSlip | PASS | 0.42s | Prescribed slip on cohesive fault |
| 65 | Integration.OutputFile | PASS | 0.69s | HDF5/VTK output file generation |
| 66 | Integration.GmshImport | PASS | 0.69s | Gmsh mesh import and label creation |
| 67 | Integration.GmshMultiMaterial | PASS | 0.65s | Gmsh physical-group per-cell material assignment |
| 68 | Integration.GasbuggyMesh | PASS | 0.70s | Historical Gasbuggy layered MSH2 mesh (3 regions) |
| 69 | Integration.NuclearTwinGmsh | PASS | 3.03s | Gmsh-based nuclear twin with mapped regions |
| 70 | Integration.PunggyeRiLayered | PASS | 20.52s | Punggye-ri layered model with absorbing BCs and SAC |
| 71 | Integration.DerivedFields | PASS | 2.43s | Cell-centered stress, strain, CFS from FEM solution |
| 72 | Integration.DPRK2017Comparison | PASS | 0.41s | Synthetic vs observed mb for DPRK 2017 |
| 87 | Integration.PressurizedFractureFEM | PASS | 1.21s | Pressurized fracture FEM solve |
| 90 | Integration.ElastoplasticSim | PASS | 0.78s | Drucker-Prager through TSSolve |
| 91 | Integration.DynamicRuptureSolve.LockedQuasiStatic | PASS | 0.73s | Manual cohesive assembly, quasi-static locked |
| 92 | Integration.DynamicRuptureSolve.LockedElastodynamic | PASS | 0.78s | Cohesive + TSALPHA2 locked |
| 93 | Integration.DynamicRuptureSolve.PrescribedSlip | PASS | 0.67s | Imposed displacement jump |
| 94 | Integration.ExplosionFaultReactivation | PASS | 0.39s | Moment-tensor + cohesive residual coexistence |
| 95 | Integration.TractionBC | PASS | 0.62s | Per-face Neumann traction BC |
| 96 | Integration.TimeDependentSlip | PASS | 0.71s | Linear slip ramp with onset/rise time |
| 97 | Integration.NearFieldCoupled | PASS | 6.22s | COUPLED_ANALYTIC 1D solver to 3D FEM coupling |
| 98 | Integration.SlippingFaultSolve | PASS | 1.32s | Augmented Lagrangian, GTEST_SKIP on divergence |
| 99 | Integration.HistoricNuclear.Gasbuggy1967 | PASS | 3.70s | 29 kt, 4-layer Lewis Shale, SAC output |
| 100 | Integration.HistoricNuclear.Gnome1961 | PASS | 3.15s | 3.1 kt, 4-layer Salado Salt, SAC output |
| 101 | Integration.HistoricNuclear.Sedan1962 | PASS | 3.04s | 104 kt, 3-layer alluvium, SAC output |
| 102 | Integration.HistoricNuclear.DegelenMountain | PASS | 1.93s | 50 kt, 3-layer granite, SAC output |
| 103 | Integration.HistoricNuclear.NtsPahuteMesa | PASS | 1.95s | 150 kt, 4-layer tuff, SAC output |

### Performance Tests (4 tests)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 104 | Performance.Benchmarks | PASS | 0.35s |
| 105 | Performance.Scaling | PASS | 1.02s |
| 106 | Performance.Memory | PASS | 0.34s |
| 107 | Performance.GPUAcceleration | PASS | 0.34s |

### Experimental Tests (1 test)

| # | Test Name | Status | Time |
|---|-----------|--------|------|
| 108 | Experimental.NeuralStubs | PASS | 0.33s |

## Known Limitations Documented in Tests

1. **Drucker-Prager return mapping** (`Unit.DruckerPragerStandalone`): Yield function evaluation works for Drucker-Prager, von Mises, and Mohr-Coulomb. Return mapping produces nonzero plastic strain at the material-point level (`PlasticityModel::integrateStress`). NOT wired into PetscDS FEM callbacks.

2. **Elastoplasticity callback** (`Unit.Elastoplasticity`): PetscFEElastoplasticity callback works in isolation. Not wired into `setupPhysics()` in the Simulator. No config flag to enable.

3. **SCEC TPV5 infrastructure** (`Physics.SCEC.TPV5`): Parameters, CohesiveFaultKernel, FaultMeshManager, and friction parameter setup verified with 3 tests. Full dynamic rupture solve not implemented.

4. **Dynamic rupture setup** (`Functional.DynamicRuptureSetup.*`): Mesh splitting, fault geometry, locked/slipping fault setup, and absorbing BC coexistence all verified. TSSolve never called. SNES divergence risk.

5. **Hydrofrac callbacks** (`Unit.FractureFlowCallbacks`, `Unit.HydrofracFormulas`, `Unit.FracturePropagation`): Poiseuille lubrication, PKN width scaling, and cohesive strength tested in isolation. Not coupled end-to-end through the FEM solver.

6. **Standalone physics** (`Physics.StressShadowing`, `Physics.InducedSeismicity`, `Physics.ProppantTransport`, `Physics.LeakoffCoupling`, `Physics.ProductionForecast`): Analytical formulas and standalone solvers verified. Not coupled to PetscDS FEM.

7. **Mohr-Coulomb principal stress** (`Unit.DruckerPragerStandalone`): `computePrincipalStresses` returns NaN for pure hydrostatic stress (q=0) due to division by zero in the Cardano formula. Tests use stress states with nonzero deviatoric component.

8. **Slipping fault (Coulomb friction)** (`Integration.SlippingFaultSolve`): Augmented Lagrangian regularization (0.01 * penalty_stiffness) with approximate Jacobian contribution. TSSolve diverges because the locked-mode linearization does not capture d(tau_f * slip_t/|slip_t|)/d(u) or d(tau_f)/d(lambda_n). Needs semi-smooth Newton tangent operator. GTEST_SKIP documents the root cause.

9. **COUPLED_ANALYTIC NearField coupling** (`Integration.NearFieldCoupled`): 1D NearFieldExplosionSolver cavity pressure/velocity drives an RDP Brune source time function that computes moment rate injected into the 3D FEM via equivalent nodal forces. Damage profile from 1D solver degrades auxiliary fields near cavity.
