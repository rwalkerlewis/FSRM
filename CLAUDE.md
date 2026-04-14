# CLAUDE.md -- FSRM Persistent Context

This file provides all context needed to work on FSRM. Read it fully before making any changes.

## Project Summary

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed.

Repository: github.com/rwalkerlewis/FSRM
Branch: main

### Codebase Size

After cleanup, the live codebase is:
- 58 source files (.cpp): ~50,000 lines
- 69 header files (.hpp): ~32,000 lines
- Total live: ~82,000 lines

Dead code (~45,000 lines across ~60 files) has been moved to `archive/src/` and `archive/include/`.

### Test Suite

101 registered tests. All pass or GTEST_SKIP. Zero failures.

| Category | Tests | Description |
|----------|------:|-------------|
| Unit | 36 | Standalone formula, callback, and component tests |
| Functional | 10 | Setup pipeline verification (no TSSolve) |
| Physics | 25 | Analytical solutions, FEM-coupled benchmarks, standalone physics |
| Integration | 25 | Full Simulator pipeline through TSSolve, plus traction BC, time-dependent slip, and explosion-fault residual coexistence |
| Performance | 4 | Benchmarks, scaling, memory, GPU |
| Experimental | 1 | Neural operator stubs |

Note: the performance test `GPUAcceleration` uses `GTEST_SKIP` when CUDA is not available.
Note: some tests may fail when run in parallel (`ctest -j`) due to HDF5 output file conflicts. All pass when run individually.

## Build Environment

Docker-based. PETSc 3.22.2 with ctetgen.

```bash
# Build
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF && make -j$(nproc)'

# Test
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

# Interactive shell
docker run --rm -it -v $(pwd):/workspace -w /workspace fsrm-ci:local bash
```

### GPU Acceleration (PETSc CUDA)

FSRM supports GPU acceleration via PETSc native CUDA backend. No FSRM code changes needed. Build PETSc with `--with-cuda` (Dockerfile.cuda) and add runtime flags:

```bash
docker build -f Dockerfile.cuda -t fsrm-cuda:local .
docker run --gpus all --rm -v $(pwd):/workspace -w /workspace/build-cuda fsrm-cuda:local \
  ./fsrm -c ../config/examples/punggye_ri_layered.config \
  -vec_type cuda -mat_type aijcusparse -log_view
```

The PetscDS pointwise callbacks (f0, f1, g0, g3) run on CPU. PETSc handles all vector operations via cuBLAS, matrix operations via cuSPARSE, and KSP solves on GPU. Data transfer is automatic.

## MANDATORY: Check PETSc 3.22.2 API Before Use

Before calling ANY PETSc function, verify its signature exists in the installed headers:

```bash
grep -rn "FunctionName" /opt/petsc-3.22.2/include/
```

Do not assume PETSc API signatures from memory. They change between versions.

## Architecture

### Pipeline (main.cpp -> Simulator)

```
initializeFromConfigFile() -> parse config, materials, fluids, explosion, injection, hydrofrac
setupDM() -> create mesh (simplex if faults enabled)
setupFaultNetwork() -> PyLith workflow: submesh + DMPlexConstructCohesiveCells
labelBoundaries() -> label 6 faces of bounding box by coordinate
setupFields() -> create PetscFE, DMCreateDS, setupBoundaryConditions, rebuild section
setupPhysics() -> register PetscDS callbacks + unified constants
setupTimeStepper() -> create TS with IFunction/IJacobian
setupSolvers() -> configure SNES/KSP
setInitialConditions() -> zero solution, locate injection/explosion cells
run() -> TSSolve
```

FormFunction uses DMPlexTSComputeIFunctionFEM for volume residual, plus addExplosionSourceToResidual (displacement DOFs via moment tensor equivalent nodal forces) and addInjectionToResidual (pressure DOFs) for point sources, plus addFaultPressureToResidual for pressurized fractures.

### CRITICAL: DS/BC Ordering in setupFields()

```cpp
DMCreateDS(dm); // 1. PETSc 3.22 requires DS before DMAddBoundary
setupBoundaryConditions(); // 2. DMAddBoundary calls (needs DS)
DMSetLocalSection(dm, nullptr); // 3. Clear cached section
DMSetUp(dm); // 4. Rebuild section WITH BC constraints
DMGetDS(dm, &prob); // 5. Get DS for later use
```

**NEVER change this ordering. It was debugged over multiple sessions.**

### PDE Assembly

PETSc DMPlex unstructured FEM with PetscDS pointwise callbacks (f0, f1, g0, g3).

### Unified Constants Array (36 elements)

Set once in `setupPhysics()`:

```
[0] lambda [1] mu [2] rho_s
[3] phi [4] kx [5] ky [6] kz
[7] cw [8] co [9] cg
[10] mu_w [11] mu_o [12] mu_g
[13] Swr [14] Sor [15] Sgr
[16] nw [17] no [18] ng
[19] krw0 [20] kro0 [21] krg0
[22] biot_alpha [23] 1/M [24] rho_f
[25] cohesive_mode [26] cohesive_mu_f [27-31] reserved
[32] ep_cohesion [33] ep_friction_angle
[34] ep_dilation_angle [35] ep_hardening_modulus
```

When auxiliary fields are used, callbacks read material properties from `a[]`/`a_x[]` via `aOff[]` instead of from `constants[]`. The constants array remains for non-material parameters.

## Feature Inventory

### WORKS: Integration-tested through TSSolve (23+ tests)

| Feature | Test Name | What It Proves |
|---------|-----------|----------------|
| Elastostatics | Physics.ElastostaticsPatch, Physics.LithostaticStress | Hooke stress, patch test, K0 ratio |
| Elastodynamics | Physics.LambsProblem, Physics.GarvinsProblem | Wave propagation, analytical error norms |
| Poroelasticity | Physics.TerzaghiConsolidation | Biot coupling, analytical consolidation |
| Absorbing BCs | Physics.AbsorbingBC | Clayton-Engquist, >99% energy absorption |
| Gravity body force | Physics.GravityLithostatic | Lithostatic column, K0 within 5% |
| Moment tensor source | Physics.MomentTensorSource | FEM equivalent nodal forces |
| Explosion seismograms | Integration.ExplosionSeismogram | Source -> waves -> SAC output |
| Injection pressure | Integration.InjectionPressure | Poroelastic injection end-to-end |
| Depth-layered material | Integration.LayeredElastostatics | Aux field material assignment |
| Gmsh mesh import | Integration.GmshImport | MSH2 physical names, tet cells |
| Gmsh multi-material | Integration.GmshMultiMaterial, Integration.GasbuggyMesh | Per-label lambda/mu/rho |
| Gmsh nuclear twin | Integration.NuclearTwinGmsh | Mapped materials + explosion + HDF5 |
| Explosion damage zones | Physics.ExplosionDamageZone | Degraded aux near cavity |
| Punggye-ri layered | Integration.PunggyeRiLayered | 3-layer + absorbing + SAC + HDF5 |
| Pressurized fracture | Integration.PressurizedFractureFEM | Cohesive traction through TSSolve |
| Elastoplasticity | Integration.ElastoplasticSim | Drucker-Prager through TSSolve |
| Locked fault (quasi-static) | Integration.DynamicRuptureSolve.LockedQuasiStatic | Manual cohesive assembly |
| Locked fault (elastodynamic) | Integration.DynamicRuptureSolve.LockedElastodynamic | Cohesive + TSALPHA2 |
| Prescribed slip | Integration.DynamicRuptureSolve.PrescribedSlip | Imposed displacement jump |
| Locked fault transparency | Physics.LockedFaultTransparency | Fault slip < 5e-4 |
| Per-face traction BC | Integration.TractionBC | Manual assembly, uniaxial analytical |
| Time-dependent slip ramp | Integration.TimeDependentSlip | Linear slip ramp with onset/rise time |
| Derived fields | Integration.DerivedFields | Stress/strain/CFS from solution |
| HDF5/VTK output | Integration.OutputFile | PetscViewerHDF5, VTK |
| Restart | Integration.Restart | Checkpoint/restore lifecycle |
| DPRK 2017 mb | Integration.DPRK2017Comparison | Synthetic vs observed body-wave magnitude |
| Explosion+fault residual | Integration.ExplosionFaultReactivation | Coexistence of moment-tensor and cohesive residual |

### WORKS: Standalone (correct, tested, not FEM-coupled)

| Feature | Test Name |
|---------|-----------|
| Mueller-Murphy source | Physics.MuellerMurphy |
| Near-field explosion (1D) | Physics.NearFieldExplosion |
| Atmospheric explosion | Physics.AtmosphericExplosion |
| Plasticity return mapping | Unit.DruckerPragerStandalone |
| Friction laws | Unit.FaultMechanics |
| Coulomb stress transfer | Unit.CoulombStressTransfer |
| Hydrofrac formulas | Unit.HydrofracFormulas + 8 Physics tests |
| Fluid flow callbacks | Unit.SinglePhaseFlow, Unit.MultiphaseFlow |

### DOES NOT WORK: Code Exists but No TSSolve Test

| Feature | What is Missing |
|---------|-----------------|
| Multiphase flow end-to-end | Buckley-Leverett or waterflood through TSSolve |
| Hydraulic fracture coupled solve | Full lubrication + deformation coupling |
| Slipping fault (Coulomb friction) | TSSolve test with shear loading |
| Explosion + fault full TSSolve | Diverges; residual coexistence verified only |

### DEAD CODE (Archived to archive/src/)

~45,000 lines across ~60 files. Categories:

| Category | Examples |
|----------|---------|
| Volcano modeling | VolcanoModel.cpp, VolcanoModule.cpp |
| Ocean/tsunami | HydrodynamicModel, OceanPhysics, TsunamiModel |
| DG/ADER | DiscontinuousGalerkin.cpp |
| High-fidelity stubs | HighFidelityNumerics, AdvancedCoupledPhysics |
| Material library | MaterialLibrary.cpp |
| Radiation | RadiationPhysics.cpp |
| AMR | AdaptiveMeshRefinement.cpp (header kept for test dep) |
| Reservoir stubs | TracerModel, AquiferModel, TwoPhaseFlow |
| Other | Module stubs, VFPTables, GnuplotViz |

Do not reference these in documentation or claims. Do not try to use them.

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
| Runnable examples | examples/ | -- |
| Visualization scripts | scripts/plot_seismograms.py, scripts/plot_wavefield.py | -- |

## Runnable Examples

Eight examples in `examples/`, each with README.md and run.sh:

| # | Directory | Config | Physics |
|---|-----------|--------|---------|
| 01 | examples/01_uniaxial_compression | uniaxial_compression.config | Elastostatics, BCs |
| 02 | examples/02_explosion_seismogram | explosion_seismogram_quick.config | Elastodynamics, SAC |
| 03 | examples/03_elastoplastic_compression | elastoplastic_compression.config | Drucker-Prager |
| 04 | examples/04_locked_fault | locked_fault_compression.config | Cohesive cells |
| 05 | examples/05_punggye_ri_nuclear_test | punggye_ri_layered_quick.config | Layered, SAC, HDF5 |
| 06 | examples/06_gmsh_multimaterial | nuclear_twin_gmsh.config | Gmsh, per-region |
| 07 | examples/07_traction_bc | traction_bc.config | Per-face traction BC |
| 08 | examples/08_time_dependent_slip | time_dependent_slip.config | Slip ramp |

## Roadmap: Features to Implement

Each item requires: PetscDS callbacks integrated into setupPhysics(), integration tests
through TSSolve, example config, and visualization. Source code in archive/src/ may provide
a starting point but must be rewritten to use the PetscDS callback pattern.

1. Slipping fault TSSolve (Coulomb friction through manual assembly)
2. Multiphase flow end-to-end (Buckley-Leverett waterflood)
3. Full coupled hydraulic fracturing (lubrication + deformation)
4. Viscoelastic attenuation (Q-factor memory variables)
5. Thermal coupling (heat equation + THM Biot)
6. Radiation transport (advection-diffusion for fallout)
7. Per-cell material from velocity model files
8. SCEC TPV5 dynamic rupture benchmark
9. Multi-stage hydraulic fracturing with stress shadowing
10. Production forecasting through propped fracture

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.22.2 API signatures before calling any PETSc function.
3. All existing tests must continue to pass after every change.
4. NEVER change the DS/BC ordering in setupFields().
5. Do NOT modify callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp.
6. Do NOT modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS.
7. Dead code has been moved to archive/. Do not reference it in documentation or claims.
8. No Python. Everything in C++ within the Simulator. Python is ONLY for post-processing visualization.
9. Ignore everything in config/aspirational/ (archived). Those configs reference features that do not exist.
10. Working examples are in config/examples/ with `.config` extension.
11. The executable is `fsrm`, not `fsrm_simulator`. Run as `./fsrm -c ../config/examples/example.config`.
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
ctest -j$(nproc) -L "unit" --output-on-failure
ctest -j$(nproc) -L "integration" --output-on-failure
```

### Build + Test in One Docker Command

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'cd build && make -j$(nproc) && ctest -j$(nproc) --output-on-failure'
```

## Reference Implementations

- PyLith (github.com/geodynamics/pylith): DMPlex cohesive cells
- PETSc examples: ex17.c (elasticity), ex56.c (elasticity with BCs), ex62.c (cohesive)
- Auxiliary field pattern: PETSc ex17.c shows DMAux setup for heterogeneous material
