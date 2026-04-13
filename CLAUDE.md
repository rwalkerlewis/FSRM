# CLAUDE.md -- FSRM Persistent Context

This file provides all context needed to work on FSRM. Read it fully before making any changes.

## Project Summary

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed. ~134k lines.

Repository: github.com/rwalkerlewis/FSRM
Branch: main

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

`Simulator::FormFunction` converts global to local vectors, calls `DMPlexInsertBoundaryValues`, then `DMPlexTSComputeIFunctionFEM`. After FEM residual, it adds injection source (pressure DOFs) and explosion source (displacement DOFs) via `DMPlexVecGetClosure`/`SetClosure`.

### Unified Constants Array (27 elements)

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
[25] cohesive_mode    [26] cohesive_mu_f
```

When auxiliary fields are used (Phase 2+), callbacks read material properties from `a[]`/`a_x[]` via `aOff[]` instead of from `constants[]`. The constants array remains for non-material parameters.

## File Locations

| Component | Source | Header |
|---|---|---|
| Simulator | src/core/Simulator.cpp | include/core/Simulator.hpp |
| Config parsing | src/core/ConfigReader.cpp | include/core/ConfigReader.hpp |
| Main | src/core/main.cpp | -- |
| Elasticity callbacks | src/numerics/PetscFEElasticity.cpp | include/numerics/PetscFEElasticity.hpp |
| Poroelasticity callbacks | src/numerics/PetscFEPoroelasticity.cpp | include/numerics/PetscFEPoroelasticity.hpp |
| Fluid flow callbacks | src/numerics/PetscFEFluidFlow.cpp | include/numerics/PetscFEFluidFlow.hpp |
| Explosion physics | src/domain/explosion/ExplosionImpactPhysics.cpp | include/domain/explosion/ExplosionImpactPhysics.hpp |
| Fault mesh | src/numerics/FaultMeshManager.cpp | include/numerics/FaultMeshManager.hpp |
| Seismometers | src/domain/seismic/SeismometerNetwork.cpp | include/domain/seismic/SeismometerNetwork.hpp |
| Fracture model | src/domain/geomechanics/FractureModel.cpp | include/domain/geomechanics/FractureModel.hpp |
| LEFM | src/physics/PhysicsKernel_LEFM.cpp | include/physics/PhysicsKernel_LEFM.hpp |
| Cohesive kernel | src/physics/CohesiveFaultKernel.cpp | include/physics/CohesiveFaultKernel.hpp |
| Unit tests | tests/unit/ | -- |
| Integration tests | tests/integration/ | -- |
| Physics validation | tests/physics_validation/ | -- |
| Example configs | config/examples/ | -- |
| Aspirational configs | config/aspirational/ | DO NOT USE |

## What Works

### Verified End-to-End Through TSSolve (91 Tests Pass)

These features are run through the Simulator with TSSolve and produce verified results:

1. **Elastostatics**: PetscFEElasticity f0/f1/g3 callbacks. Uniaxial compression, patch test, lithostatic stress. Integration-tested.
2. **Poroelasticity**: PetscFEPoroelasticity f0/f1 for pressure+displacement, all 4 Jacobian blocks. Terzaghi consolidation passes with analytical comparison. Integration-tested.
3. **Elastodynamics**: Lamb's problem and Garvin's problem pass with quantitative error norms. Integration-tested.
4. **Explosion seismogram pipeline**: Moment tensor source injection -> elastodynamic solve -> seismometer sampling -> SAC output. Integration-tested.
5. **Absorbing BCs**: Clayton-Engquist first-order. Energy absorption >99% at normal incidence. Integration-tested.
6. **Gravity body force**: Density-scaled gravity callback. Lithostatic stress column with analytical comparison passes within 5% tolerance. K0 ratio verified. Integration-tested.
7. **Moment tensor source injection**: Proper FEM equivalent nodal forces via PetscFECreateTabulation. Produces nonzero displacement. Integration-tested.
8. **Injection pressure buildup**: Poroelastic Simulator with point injection source runs end-to-end. Two-simulation comparison confirms injection affects solution. Integration-tested.
9. **Derived fields**: Cell-centered stress, strain, CFS from FEM solution. Integration-tested.
10. **Layered elastostatics**: Depth-based material layering via auxiliary fields. Integration-tested.
11. **Simulation lifecycle and restart**: Checkpoint/restart verified. Integration-tested.
12. **DPRK 2017 synthetic mb**: Synthetic body-wave magnitude vs observed for 250 kt. 4 tests pass.
13. **Seismometer network SAC output**: Coupled to explosion seismogram pipeline. Integration-tested.
14. **Boundary conditions**: labelBoundaries() labels 6 box faces, setupBoundaryConditions() registers via DMAddBoundary. 5 unit tests pass.

Additional verified coverage in the current Docker test suite:

- **Gmsh label-based material assignment**: `populateAuxFieldsByMaterialLabel()` maps Gmsh physical-group labels to per-cell auxiliary `lambda`, `mu`, and `rho` values. Verified by `Integration.GmshMultiMaterial` and `Integration.GasbuggyMesh`.
- **Historical Gasbuggy mesh**: `meshes/historical/gasbuggy_layered_3d.msh` loads through PETSc DMPlex with three validated material regions.
- **Gmsh nuclear twin workflow**: A compact mapped Gmsh multi-material mesh with underground source and HDF5 output is verified by `Integration.NuclearTwinGmsh`.
- **Explosion damage-zone aux degradation**: `applyExplosionDamageToAuxFields()` degrades auxiliary properties near the source and changes the dynamic response. Verified by `Physics.ExplosionDamageZone`.
- **Punggye-ri layered workflow**: A compact three-layer underground nuclear test case with absorbing boundaries, HDF5 output, and SAC seismograms is verified by `Integration.PunggyeRiLayered`.

### Callback-Verified Only (Not Run Through TSSolve)

These have correct callbacks wired into setupPhysics, but no test calls sim.run():

15. **Cohesive fault mesh splitting**: PyLith workflow verified. 32 cohesive cells on 4x4x4 simplex mesh. Setup completes but TSSolve never called with cohesive cells.
16. **Dynamic rupture setup**: Simulator setup pipeline with cohesive cells completes without error. TSSolve never called.
17. **Fault plus absorbing BC coexistence**: DMSetRegionDS region split verified. setupPhysics and setupTimeStepper succeed. TSSolve never called.
18. **Pressurized hydrofracture callbacks**: f0_lagrange_pressure_balance registered on cohesive DS. Callback tested in isolation. Never solved.
19. **Fracture lubrication flow callbacks**: f0_fracture_pressure, f1_fracture_pressure registered. Callback tested in isolation. Never solved.
20. **Elastoplasticity callback**: f1_elastoplastic_aux with Drucker-Prager return mapping. Callback returns stress to yield surface. Not wired into setupPhysics.

### Unit-Tested Utility Functions (No FEM Context)

These are standalone analytical formulas, correct but not FEM-integrated:

21. **Fluid flow callbacks**: PetscFEFluidFlow single-phase and black-oil. Unit tested. NOT verified end-to-end.
22. **Friction laws**: Slip-weakening and rate-state (aging). Unit tested.
23. **CoulombStressTransfer**: Hooke stress, fault projection, delta_CFS. Unit tested.
24. **Mueller-Murphy source**: Corner frequency, mb-yield scaling, RDP spectrum. 9 physics tests pass.
25. **Atmospheric explosion effects**: Sedov-Taylor blast, Brode fireball, EMP E1. 8 tests pass.
26. **Near-field explosion phenomenology**: Cavity radius, damage zones, spall. 13 tests pass.
27. **SCEC TPV5 infrastructure**: Parameters, kernel construction, friction setup. 3 tests pass. Full solve WIP.
28. **Plasticity yield evaluation**: DP, VM, MC yield surface detection works. Return mapping works at material point level. 8 tests pass.
29. **Hydrofrac utility functions**: Sneddon aperture, PKN width scaling, Carter leak-off, stress shadowing, cluster efficiency, induced seismicity (M0/Mw), proppant transport, Arps decline curves, fracture productivity index. All correct. None assembled/solved via PetscDS.

## What Does NOT Work (Known Gaps)

1. **Elastoplasticity not wired into Simulator**: f1_elastoplastic_aux callback exists and works in isolation. PlasticityModel::integrateStress works. Neither is wired into setupPhysics(). No config flag to enable.
2. **Hydrofrac never solved through TSSolve**: All hydrofrac PetscDS callbacks are registered but never assembled and solved. Pressurized fracture, lubrication flow, coupled flow-deformation are all callback-tested only.
3. **Dynamic rupture never solved**: Cohesive cells inserted into mesh, callbacks registered, but TSSolve never called. Reason: "SNES may diverge at the first time step because the cohesive Jacobian blocks may not be fully consistent."
4. **Fault + absorbing coexistence never solved**: DMSetRegionDS wired, setup succeeds, TSSolve never called.
5. **Stubs only**: DG/ADER, GPU (CUDA/HIP), FNO/ML solvers, volcano, tsunami, ocean, infrasound, radiation transport, hypervelocity impacts. Headers/docs exist but no functional implementation.
6. **Dead hydrofrac code in MonitorFunction**: The hydrofrac_ member uses a standalone analytical Thiem equation estimate disconnected from the FEM solution. Redundant with PetscFEHydrofrac callbacks.

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.22.2 API signatures before calling any PETSc function.
3. All existing tests must continue to pass after every change.
4. NEVER change the DS/BC ordering in setupFields().
5. Do NOT modify callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp. New callbacks for auxiliary fields go in new files.
6. Do NOT modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS.
7. Do NOT use DG, ADER, GPU, or plasticity features. They are stubs. Plasticity note: DruckerPragerModel, VonMisesModel, MohrCoulombModel exist in PlasticityModel.hpp with yield function evaluation (works) and return mapping algorithms (broken -- does not produce nonzero plastic strain). Not wired into the PETSc FEM pipeline. See tests/physics_validation/test_drucker_prager.cpp for details.
8. No Python. Everything in C++ within the Simulator.
9. Ignore everything in config/aspirational/. Those configs reference features that do not exist.
10. Working examples are in config/examples/ with `.config` extension (INI-style format with `[SECTION]\nkey = value`): uniaxial_compression, terzaghi_consolidation, injection_pressure_buildup, hydraulic_fracture_pkn, cohesive_hydraulic_fracture, explosion_seismogram, dprk_2017_quick, underground_explosion_template.
13. The executable is `fsrm`, not `fsrm_simulator`. Run as `./fsrm ../config/examples/example.config`.
11. No em dashes or contractions in code comments or documentation.
12. Verification tests must have quantitative pass/fail criteria with numerical tolerances.
13. At the end of each prompt run review the current state of the code and update CLAUDE.md as well as README.md to reflect the current state of the code and what works and is verified.

## Reference Implementations

- PyLith (github.com/geodynamics/pylith): DMPlex cohesive cells, libsrc/pylith/faults/TopologyOps.cc
- PETSc examples: ex17.c (elasticity), ex56.c (elasticity with BCs), ex62.c (cohesive)
- Auxiliary field pattern: PETSc ex17.c shows DMAux setup for heterogeneous material