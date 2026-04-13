# CLAUDE.md -- FSRM Persistent Context

This file provides all context needed to work on FSRM. Read it fully before making any changes.

## Project Summary

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed.

Repository: github.com/rwalkerlewis/FSRM
Branch: main

### Codebase Size

- 100 source files (.cpp): 86,658 lines
- 117 header files (.hpp): 53,174 lines
- Total: ~140,000 lines
- **Live code** (referenced by Simulator and tested): ~35,000 lines in ~32 source files
- **Tested standalone** (not Simulator-integrated): ~8,000 lines in ~8 source files
- **Dead code** (never referenced, never tested): ~44,000 lines in ~60 source files

### Test Suite

95 registered tests. All pass or GTEST_SKIP. Zero failures.

| Category | Tests | Description |
|----------|------:|-------------|
| Unit | 36 | Standalone formula, callback, and component tests |
| Functional | 10 | Setup pipeline verification (no TSSolve) |
| Physics | 24 | Analytical solutions, FEM-coupled benchmarks, standalone physics |
| Integration | 21 | Full Simulator pipeline through TSSolve |
| Performance | 3 | Benchmarks, scaling, memory |
| Experimental | 1 | Neural operator stubs |

Note: 4 integration tests (DynamicRuptureSolve, ExplosionFaultReactivation) use GTEST_SKIP when TSSolve diverges due to PetscDSSetBdResidual on cohesive cells in PETSc 3.22.

## Build Environment

Docker-based. PETSc 3.22.2 with ctetgen.

```bash
# Build
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'

# Test
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

# Interactive shell
docker run --rm -it -v $(pwd):/workspace -w /workspace fsrm-ci:local bash
```

## MANDATORY: Check PETSc 3.22.2 API Before Use

Before calling ANY PETSc function, verify its signature exists in the installed headers:

```bash
grep -rn "FunctionName" /opt/petsc-3.22.2/include/
```

Do not assume PETSc API signatures from memory. They change between versions.

## Architecture

### Pipeline (main.cpp -> Simulator)

```
initializeFromConfigFile()  -> parse config, materials, fluids, explosion, injection, hydrofrac
setupDM()                   -> create mesh (simplex if faults enabled)
setupFaultNetwork()         -> PyLith workflow: submesh + DMPlexConstructCohesiveCells
labelBoundaries()           -> label 6 faces of bounding box by coordinate
setupFields()               -> create PetscFE, DMCreateDS, setupBoundaryConditions, rebuild section
setupPhysics()              -> register PetscDS callbacks + unified constants
setupTimeStepper()          -> create TS with IFunction/IJacobian
setupSolvers()              -> configure SNES/KSP
setInitialConditions()      -> zero solution, locate injection/explosion cells
run()                       -> TSSolve
```

FormFunction uses DMPlexTSComputeIFunctionFEM for volume residual, plus addExplosionSourceToResidual (displacement DOFs via moment tensor equivalent nodal forces) and addInjectionToResidual (pressure DOFs) for point sources, plus addFaultPressureToResidual for pressurized fractures.

### CRITICAL: DS/BC Ordering in setupFields()

```cpp
DMCreateDS(dm);                          // 1. PETSc 3.22 requires DS before DMAddBoundary
setupBoundaryConditions();               // 2. DMAddBoundary calls (needs DS)
DMSetLocalSection(dm, nullptr);          // 3. Clear cached section
DMSetUp(dm);                             // 4. Rebuild section WITH BC constraints
DMGetDS(dm, &prob);                      // 5. Get DS for later use
```

**NEVER change this ordering. It was debugged over multiple sessions.**

### PDE Assembly

PETSc DMPlex unstructured FEM with PetscDS pointwise callbacks (f0, f1, g0, g3).

### Unified Constants Array (36 elements)

Set once in `setupPhysics()`:

```
[0]  lambda           [1]  mu              [2]  rho_s
[3]  phi              [4]  kx              [5]  ky             [6]  kz
[7]  cw               [8]  co              [9]  cg
[10] mu_w             [11] mu_o            [12] mu_g
[13] Swr              [14] Sor             [15] Sgr
[16] nw               [17] no              [18] ng
[19] krw0             [20] kro0            [21] krg0
[22] biot_alpha        [23] 1/M             [24] rho_f
[25] cohesive_mode    [26] cohesive_mu_f   [27-31] reserved
[32] ep_cohesion      [33] ep_friction_angle
[34] ep_dilation_angle [35] ep_hardening_modulus
```

When auxiliary fields are used, callbacks read material properties from `a[]`/`a_x[]` via `aOff[]` instead of from `constants[]`. The constants array remains for non-material parameters.

## File Locations

| Component | Source | Header |
|---|---|---|
| Simulator | src/core/Simulator.cpp | include/core/Simulator.hpp |
| Config parsing | src/core/ConfigReader.cpp | include/core/ConfigReader.hpp |
| Main | src/core/main.cpp | -- |
| Derived fields | src/core/DerivedFieldComputer.cpp | include/core/DerivedFieldComputer.hpp |
| Elasticity callbacks | src/numerics/PetscFEElasticity.cpp | include/numerics/PetscFEElasticity.hpp |
| Elasticity aux callbacks | src/numerics/PetscFEElasticityAux.cpp | include/numerics/PetscFEElasticityAux.hpp |
| Poroelasticity callbacks | src/numerics/PetscFEPoroelasticity.cpp | include/numerics/PetscFEPoroelasticity.hpp |
| Fluid flow callbacks | src/numerics/PetscFEFluidFlow.cpp | include/numerics/PetscFEFluidFlow.hpp |
| Elastoplasticity callback | src/numerics/PetscFEElastoplasticity.cpp | include/numerics/PetscFEElastoplasticity.hpp |
| Hydrofrac callbacks | src/numerics/PetscFEHydrofrac.cpp | include/numerics/PetscFEHydrofrac.hpp |
| Absorbing BC | src/numerics/AbsorbingBC.cpp | include/numerics/AbsorbingBC.hpp |
| Gravity body force | src/numerics/PetscFEElasticityGravity.cpp | include/numerics/PetscFEElasticityGravity.hpp |
| Boundary conditions | src/numerics/BoundaryConditions.cpp | include/numerics/BoundaryConditions.hpp |
| Fault mesh | src/numerics/FaultMeshManager.cpp | include/numerics/FaultMeshManager.hpp |
| Cohesive kernel | src/physics/CohesiveFaultKernel.cpp | include/physics/CohesiveFaultKernel.hpp |
| Explosion physics | src/domain/explosion/ExplosionImpactPhysics.cpp | include/domain/explosion/ExplosionImpactPhysics.hpp |
| Seismometers | src/domain/seismic/SeismometerNetwork.cpp | include/domain/seismic/SeismometerNetwork.hpp |
| Fracture model | src/domain/geomechanics/FractureModel.cpp | include/domain/geomechanics/FractureModel.hpp |
| Gmsh I/O | src/io/GmshIO.cpp | include/io/GmshIO.hpp |
| Unit tests | tests/unit/ | -- |
| Integration tests | tests/integration/ | -- |
| Physics validation | tests/physics_validation/ | -- |
| Example configs | config/examples/ | -- |
| Aspirational configs | config/aspirational/ | DO NOT USE |

## What Runs Through TSSolve (Integration-Tested)

Every feature below has at least one integration test that calls `sim.run()` (TSSolve) and checks quantitative results.

1. **Elastostatics**: PetscFEElasticity f0/f1/g3 callbacks. Uniaxial compression, patch test, lithostatic stress. Tests: `Physics.ElastostaticsPatch`, `Physics.LithostaticStress`, `Integration.FullSimulation`.
2. **Poroelasticity**: PetscFEPoroelasticity f0/f1 for pressure+displacement, all 4 Jacobian blocks. Tests: `Physics.TerzaghiConsolidation` (analytical comparison).
3. **Elastodynamics**: Implicit time stepping with TSALPHA2. Tests: `Physics.LambsProblem` (point force on halfspace), `Physics.GarvinsProblem` (buried explosion Green function).
4. **Explosion seismogram pipeline**: Moment tensor source injection -> elastodynamic solve -> seismometer sampling -> SAC output. Tests: `Integration.ExplosionSeismogram`.
5. **Absorbing BCs**: Clayton-Engquist first-order Lysmer-Kuhlmeyer. Energy absorption >99% at normal incidence. Tests: `Physics.AbsorbingBC`.
6. **Gravity body force**: Density-scaled gravity callback. Lithostatic stress column with K0 ratio within 5% tolerance. Tests: `Physics.GravityLithostatic`.
7. **Moment tensor source injection**: Proper FEM equivalent nodal forces via PetscFECreateTabulation. Tests: `Physics.MomentTensorSource`.
8. **Injection pressure buildup**: Poroelastic Simulator with point injection. Two-simulation comparison. Tests: `Integration.InjectionPressure`.
9. **Derived fields**: Cell-centered stress, strain, CFS from FEM solution. Tests: `Integration.DerivedFields`.
10. **Layered elastostatics**: Depth-based material layering via auxiliary fields. Tests: `Integration.LayeredElastostatics`.
11. **Simulation lifecycle and restart**: Tests: `Integration.FullSimulation`, `Integration.Restart`.
12. **DPRK 2017 synthetic mb**: Synthetic vs observed body-wave magnitude. Tests: `Integration.DPRK2017Comparison`.
13. **Boundary conditions**: 6-face labeling, DMAddBoundary registration. Tests: `Unit.BoundaryConditions`.
14. **Gmsh label-based material assignment**: Tests: `Integration.GmshMultiMaterial`, `Integration.GasbuggyMesh`.
15. **Gmsh mesh import**: Tests: `Integration.GmshImport`.
16. **Gmsh nuclear twin workflow**: Tests: `Integration.NuclearTwinGmsh`.
17. **Explosion damage-zone degradation**: Tests: `Physics.ExplosionDamageZone`.
18. **Punggye-ri layered model**: Tests: `Integration.PunggyeRiLayered`.
19. **Pressurized fracture FEM**: Tests: `Integration.PressurizedFractureFEM`.
20. **Prescribed slip**: Tests: `Integration.PrescribedSlip`.
21. **HDF5/VTK output**: Tests: `Integration.OutputFile`.
22. **Elastoplasticity (Drucker-Prager)**: PetscFEElastoplasticity f1/g3 callbacks wired into setupPhysics via `[PLASTICITY]` config. Unified constants indices 32-35. Tests: `Integration.ElastoplasticSim` (quasi-static beyond/below yield), `Physics.ElastoplasticBearingCapacity` (below-yield matches elastic, above-yield converges).

## What Has Correct Callbacks But is Not Fully Solver-Coupled

Setup pipeline completes, but TSSolve is either not called or not tested end-to-end.

1. **Cohesive fault mesh splitting**: PyLith workflow. Tests: `Functional.DynamicRuptureSetup.MeshSplitting`. TSSolve not called.
2. **Dynamic rupture setup**: Tests: `Functional.DynamicRuptureSetup.LockedFault`, `.SlippingFault`, `.FaultGeometryFromConfig`. TSSolve not called -- SNES diverges due to PetscDSSetBdResidual on cohesive cells in PETSc 3.22 (BdResidual support[0] totDim mismatch with bulk tets).
3. **Dynamic rupture TSSolve attempts**: Tests: `Integration.DynamicRuptureSolve` (locked quasi-static, locked elastodynamic, prescribed slip). All use GTEST_SKIP because SNES diverges (ierr=91). PetscReturnErrorHandler is used to prevent cascading fatal errors.
4. **Explosion-induced fault reactivation**: Tests: `Integration.ExplosionFaultReactivation`. Setup succeeds (explosion + cohesive callbacks coexist), TSSolve diverges. Uses GTEST_SKIP.
5. **Fault + absorbing BC coexistence**: Tests: `Functional.DynamicRuptureSetup.AbsorbingCoexist`, `Functional.FaultAbsorbingCoexist`. TSSolve not called.
6. **Fracture lubrication flow callbacks**: Tests: `Unit.FractureFlowCallbacks`. Callback tested in isolation.

## Standalone Verified Code (Not FEM-Coupled)

Standalone analytical formulas and solvers. Correct and tested, but not coupled to PETSc FEM.

1. **PlasticityModel return mapping**: integrateStress works and produces nonzero plastic strain. Tests: `Unit.DruckerPragerStandalone`. Also wired into PetscDS callbacks via PetscFEElastoplasticity (see item 22 above).
2. **NearFieldExplosionSolver**: 1D Lagrangian with Mie-Gruneisen EOS, damage, spall. Tests: `Physics.NearFieldExplosion`.
3. **NuclearAirburstEffects**: Sedov-Taylor, Brode, EMP. Tests: `Physics.AtmosphericExplosion`.
4. **Mueller-Murphy source**: Corner frequency, mb, RDP, moment rate. Tests: `Physics.MuellerMurphy`.
5. **Hydrofrac utilities**: Sneddon, PKN, Carter, stress shadowing, induced seismicity, proppant, Arps. Tests: `Unit.HydrofracFormulas`, `Unit.FracturePropagation`, `Unit.PressurizedFractureCallbacks`, `Physics.StressShadowing`, `Physics.InducedSeismicity`, `Physics.ProppantTransport`, `Physics.LeakoffCoupling`, `Physics.ProductionForecast`.
6. **Friction laws**: Tests: `Unit.FaultMechanics`.
7. **CoulombStressTransfer**: Tests: `Unit.CoulombStressTransfer`.
8. **SeismicSource**: Tests: `Unit.SeismicSource`.
9. **MaterialModel**: Tests: `Unit.MaterialModel`.
10. **SCEC TPV5 infrastructure**: Tests: `Physics.SCEC.TPV5`. Full solve not implemented.
11. **Fluid flow callbacks**: Tests: `Unit.SinglePhaseFlow`, `Unit.MultiphaseFlow`. NOT end-to-end verified.

## Dead Code (Never Referenced, Never Tested)

~44,000 lines across ~51 source files compile into libfsrmlib.so but are never called. Categories:

| Category | Lines | Examples |
|----------|------:|---------|
| Volcano modeling | ~4,200 | VolcanoModel.cpp, VolcanoModule.cpp |
| Ocean/tsunami | ~5,000 | HydrodynamicModel, OceanPhysics, TsunamiModel |
| ML/neural operators | ~9,200 | FNO, NeuralSurrogates, GraphNeuralOperator |
| DG/ADER | ~1,900 | DiscontinuousGalerkin.cpp |
| High-fidelity stubs | ~4,000 | HighFidelityNumerics, AdvancedCoupledPhysics |
| Material library | ~7,700 | MaterialLibrary.cpp |
| Radiation | ~1,650 | RadiationTransport.cpp |
| AMR | ~1,470 | AdaptiveMeshRefinement.cpp |
| Reservoir stubs | ~3,300 | Tracer, Aquifer, RelativePermeability |
| Other | ~5,200 | Module stubs, VFP, GnuplotViz, standalone solvers |

Do not reference these in documentation or claims. Do not try to use them.

## What Does NOT Work (Known Gaps)

1. **Hydrofrac not fully coupled**: PressurizedFractureFEM passes, but full lubrication+deformation coupled solve is callback-tested only.
2. **Dynamic rupture TSSolve diverges**: Setup completes, TSSolve diverges (ierr=91) due to PetscDSSetBdResidual on cohesive cells in PETSc 3.22. Tests use GTEST_SKIP. Pending upstream PETSc fix.
3. **Fault + absorbing coexistence never solved**: Setup succeeds, TSSolve never called.
4. **Dead code stubs**: DG/ADER, GPU, FNO/ML, volcano, tsunami, ocean, radiation, AMR. Not functional.

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.22.2 API signatures before calling any PETSc function.
3. All existing tests must continue to pass after every change.
4. NEVER change the DS/BC ordering in setupFields().
5. Do NOT modify callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp.
6. Do NOT modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS.
7. Plasticity return mapping works at the material-point level (PlasticityModel::integrateStress). PetscFEElastoplasticity callback exists but is not wired into setupPhysics(). DG, ADER, GPU, ML features are stubs. Volcano, ocean, tsunami, radiation code is dead (never referenced). Do not claim any of these are functional.
8. No Python. Everything in C++ within the Simulator.
9. Ignore everything in config/aspirational/. Those configs reference features that do not exist.
10. Working examples are in config/examples/ with `.config` extension.
11. The executable is `fsrm`, not `fsrm_simulator`. Run as `./fsrm ../config/examples/example.config`.
12. No em dashes or contractions in code comments or documentation.
13. Verification tests must have quantitative pass/fail criteria with numerical tolerances.
14. Update CLAUDE.md and README.md after every session to reflect actual code state. Every claim must be backed by a specific test name or code reference.

## Parallel Development and Execution

### Building

```bash
make -j$(nproc)
```

### Testing

```bash
ctest -j$(nproc) --output-on-failure
```

Or by category:
```bash
ctest -j$(nproc) -L "unit" --output-on-failure &
ctest -j$(nproc) -L "functional" --output-on-failure &
ctest -j$(nproc) -L "physics_validation" --output-on-failure &
ctest -j$(nproc) -L "integration" --output-on-failure &
wait
```

Tests within the same executable share a process and run sequentially. Parallelism is across executables and CTest registrations.

### Build + Test in One Docker Command

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'cd build && make -j$(nproc) && ctest -j$(nproc) --output-on-failure'
```

## Reference Implementations

- PyLith (github.com/geodynamics/pylith): DMPlex cohesive cells
- PETSc examples: ex17.c (elasticity), ex56.c (elasticity with BCs), ex62.c (cohesive)
- Auxiliary field pattern: PETSc ex17.c shows DMAux setup for heterogeneous material
