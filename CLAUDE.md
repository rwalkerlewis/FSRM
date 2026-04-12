# CLAUDE.md -- FSRM Persistent Context

This file provides all context needed to work on FSRM. Read it fully before making any changes.

## Project Summary

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed. ~134k lines.

Repository: github.com/rwalkerlewis/FSRM
Branch: local_fix

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

## What Works (Verified, All 67 Tests Pass)

1. **Elastostatics**: PetscFEElasticity f0/f1/g3 callbacks. Uniaxial compression converges with nonzero FNORM in 1 iteration. Unit tested.
2. **Poroelasticity callbacks**: PetscFEPoroelasticity f0/f1 for pressure+displacement, all 4 Jacobian blocks. Unit tested. Terzaghi consolidation passes with analytical comparison.
3. **Fluid flow callbacks**: PetscFEFluidFlow single-phase and black-oil. Unit tested. NOT verified end-to-end.
4. **Cohesive fault mesh splitting**: PyLith workflow (DMPlexCreateSubmesh -> subpoint map -> DMPlexLabelCohesiveComplete -> DMPlexConstructCohesiveCells). 32 cohesive cells, 96 fault vertices on 4x4x4 simplex mesh. All fault tests pass.
5. **Friction laws**: Slip-weakening and rate-state (aging). Tested.
6. **CoulombStressTransfer**: Hooke stress, fault projection, delta_CFS. Tested.
7. **Boundary conditions**: labelBoundaries() labels 6 box faces, setupBoundaryConditions() registers via DMAddBoundary, section rebuilt. Elastostatics BCs verified.
8. **Mueller-Murphy source**: Corner frequency, mb-yield scaling, RDP spectrum, cavity radius scaling, moment rate function. 6 physics tests pass.
9. **Absorbing BCs**: Clayton-Engquist first-order. Energy absorption >99% at normal incidence. Tested.
10. **Elastodynamics**: Lamb's problem and Garvin's problem pass with quantitative error norms.
11. **Explosion seismogram pipeline**: Source injection -> elastodynamic solve -> seismometer sampling -> SAC output. Integration-tested.
12. **DPRK 2017 synthetic mb**: Synthetic body-wave magnitude vs observed for 250 kt. 4 tests pass.
13. **Atmospheric explosion effects**: Sedov-Taylor blast, Brode fireball, EMP E1, overpressure. 6 tests pass.
14. **Near-field explosion phenomenology**: Cavity radius, damage zones, spall velocity/thickness. 6 tests pass.
15. **SCEC TPV5 infrastructure**: Parameters, CohesiveFaultKernel construction, FaultMeshManager. Verified. Full benchmark solve is WIP.
16. **Derived fields**: Cell-centered stress, strain, CFS from FEM solution. Integration-tested.
17. **Plasticity yield evaluation**: Drucker-Prager, von Mises yield surface detection works. Return mapping is broken (see gap #1).

## What Does NOT Work (Known Gaps)

1. **Plasticity return mapping broken**: DruckerPragerModel returnMapping does not produce nonzero plastic strain. Not wired into PETSc FEM pipeline. See tests/physics_validation/test_drucker_prager.cpp.
2. **Homogeneous material only**: Single constants array for the entire mesh. No per-cell material properties.
3. **No gravity / lithostatic pre-stress**: Solver starts from zero stress state.
4. **Unverified end-to-end**: Injection, hydraulic fracture, seismometer SAC output format have not been run to completion.
5. **No Gmsh mesh import tested**: DMPlexCreateGmsh exists, config stubs exist, but material region assignment from physical groups is untested.
6. **Stubs only**: DG/ADER, GPU (CUDA/HIP), FNO/ML solvers, volcano, tsunami, ocean, infrasound, radiation transport, hypervelocity impacts, ResFrac-equivalent fracturing. Headers/docs exist but no functional implementation.

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.22.2 API signatures before calling any PETSc function.
3. All existing 67 tests must continue to pass after every change.
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

## Reference Implementations

- PyLith (github.com/geodynamics/pylith): DMPlex cohesive cells, libsrc/pylith/faults/TopologyOps.cc
- PETSc examples: ex17.c (elasticity), ex56.c (elasticity with BCs), ex62.c (cohesive)
- Auxiliary field pattern: PETSc ex17.c shows DMAux setup for heterogeneous material