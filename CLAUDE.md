# CLAUDE.md -- FSRM Persistent Context

This file provides all context needed to work on FSRM. Read it fully before making any changes.

## Project Summary

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed.

Repository: github.com/rwalkerlewis/FSRM
Branch: main

### Codebase Size

After cleanup, the live codebase is:
- 63 source files (.cpp): ~54,000 lines
- 73 header files (.hpp): ~34,000 lines
- Total live: ~88,000 lines

Dead code (~45,000 lines across ~60 files) has been moved to `archive/src/` and `archive/include/`.

### Test Suite

116 registered tests. 115 pass, 1 fail honestly (no fake skips).

Known failure:
- Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic: prescribed slip produces
  near-zero displacement jump on coarse mesh. The penalty-scaled Lagrange diagonal in
  addCohesivePenaltyToJacobian slows lambda convergence for prescribed slip mode.

The quasi-static locked fault tests (LockedFaultTransparency, LockedQuasiStatic) were
fixed by replacing the old f=lambda (epsilon=1) volume regularization with weak
regularization f=epsilon*lambda, g=epsilon*I where epsilon=1e-4. The weak regularization
is negligible compared to the BdResidual constraint on cohesive faces. A penalty-scaled
Lagrange diagonal (penalty*coeff) is added manually at cohesive vertices in
addCohesivePenaltyToJacobian to keep the LU factorization well-conditioned.

Physics.SCEC.TPV5 was fixed by adding the displacement field for elastodynamics
(previously only added for geomechanics, causing error 63 when faults were enabled
without geomechanics).

| Category | Tests | Description |
|----------|------:|-------------|
| Unit | 36 | Standalone formula, callback, and component tests |
| Functional | 10 | Setup pipeline verification (no TSSolve) |
| Physics | 27 | Analytical solutions, FEM-coupled benchmarks, standalone physics, viscoelastic relaxation, cohesive BdResidual verification, SCEC TPV5 |
| Integration | 38 | Full Simulator pipeline through TSSolve, plus NearField coupling, slipping fault, slip-weakening fault, 5 historic nuclear tests, traction BC, time-dependent slip, explosion-fault residual coexistence, single-phase flow, viscoelastic wave, velocity model material, thermal diffusion, thermal expansion |
| Performance | 4 | Benchmarks, scaling, memory, GPU |
| Experimental | 1 | Neural operator stubs |

Note: `GTEST_SKIP` is used ONLY for hardware-dependent tests (GPU/CUDA detection) and one known segfault bug (explosion+fault TSSolve). A test that produces a zero solution is a FAILURE, not a skip. If SNES converges to zero, the physics setup is wrong -- fix it.
Note: some tests may fail when run in parallel (`ctest -j`) due to HDF5 output file conflicts. All pass when run individually.

## Build Environment

Docker-based. PETSc 3.25.0 with ctetgen.

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

## MANDATORY: Check PETSc 3.25.0 API Before Use

Before calling ANY PETSc function, verify its signature exists in the installed headers:

```bash
grep -rn "FunctionName" /opt/petsc-3.25.0/include/
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
applyInitialFaultStress() -> set Lagrange DOFs to initial traction (optional, for TPV5-type benchmarks)
run() -> TSSolve
```

FormFunction uses DMPlexTSComputeIFunctionFEM for volume residual (including BdResidual on cohesive cells for fault constraints), plus addExplosionSourceToResidual (displacement DOFs via moment tensor equivalent nodal forces) and addInjectionToResidual (pressure DOFs) for point sources, plus addFaultPressureToResidual for pressurized fractures. The Jacobian is assembled by DMPlexTSComputeIJacobianFEM for volume terms, plus addCohesivePenaltyToJacobian for cohesive fault terms (BdJacobian is not functional in PETSc 3.25).

### CRITICAL: DS/BC Ordering in setupFields()

```cpp
DMCreateDS(dm); // 1. PETSc 3.25 requires DS before DMAddBoundary
setupBoundaryConditions(); // 2. DMAddBoundary calls (needs DS)
DMSetLocalSection(dm, nullptr); // 3. Clear cached section
DMSetUp(dm); // 4. Rebuild section WITH BC constraints
DMGetDS(dm, &prob); // 5. Get DS for later use
```

**NEVER change this ordering. It was debugged over multiple sessions.**

### Lagrange Multiplier Field and Cohesive Cells

The Lagrange multiplier field is added via `DMAddField(dm, nullptr, fe_lagrange)` on
ALL cells. The displacement field is added for both `enable_geomechanics` and
`enable_elastodynamics`.

Weak volume regularization (f = epsilon * lambda, g = epsilon * I with epsilon = 1e-4)
keeps interior Lagrange DOFs non-singular for LU factorization. The regularization
is negligible compared to the BdResidual constraint on cohesive faces. A penalty-scaled
Lagrange diagonal (penalty * coeff ~ O(E/h)) is added manually at cohesive vertices
in addCohesivePenaltyToJacobian to match the displacement stiffness scale.

Restricting the Lagrange field to cohesive cells via `DMAddField(dm, label, ...)`
was investigated but PETSc 3.25 region DS does not support volume assembly via
DMPlexTSComputeIFunctionFEM. See docs/LAGRANGE_FIX_STATUS.md for details.

### PDE Assembly

PETSc DMPlex unstructured FEM with PetscDS pointwise callbacks (f0, f1, g0, g3).

### Unified Constants Array (up to 80 elements)

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
[36-53] traction BC (6 faces x 3 components)
[54-69] viscoelastic (N, tau, delta_mu, delta_kappa)
[70] friction_model (0=constant, 1=slip_weakening)
[71] mu_s (static friction) [72] mu_d (dynamic friction)
[73] D_c (critical slip distance)
[74] thermal_conductivity (W/(m*K))
[75] specific_heat (J/(kg*K))
[76] thermal_expansion_coeff (1/K)
[77] reference_temperature (K)
[78] thermal_field_index
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
| NearField-FEM coupling | Integration.NearFieldCoupled | COUPLED_ANALYTIC 1D solver to 3D FEM moment rate |
| Slipping fault (Coulomb) | Integration.SlippingFaultSolve | Full semi-smooth Newton Jacobian for Coulomb friction |
| Slip-weakening friction | Integration.SlipWeakeningFault | Linear slip-weakening (mu_s, mu_d, Dc) through TSSolve |
| SCEC TPV5 benchmark | Physics.SCEC.TPV5 | Dynamic rupture with initial stress and nucleation patch |
| Single-phase flow | Integration.SinglePhaseFlow | Pressure diffusion with Dirichlet pressure BCs |
| Viscoelastic attenuation | Integration.ViscoelasticWave | Generalized Maxwell body, memory variables in aux fields |
| Historic: Gasbuggy 1967 | Integration.HistoricNuclear.Gasbuggy1967 | 29 kt, 4-layer Lewis Shale, SAC output |
| Historic: Gnome 1961 | Integration.HistoricNuclear.Gnome1961 | 3.1 kt, 4-layer Salado Salt, SAC output |
| Historic: Sedan 1962 | Integration.HistoricNuclear.Sedan1962 | 104 kt, 3-layer alluvium, SAC output |
| Historic: Degelen Mountain | Integration.HistoricNuclear.DegelenMountain | 50 kt, 3-layer granite, SAC output |
| Historic: NTS Pahute Mesa | Integration.HistoricNuclear.NtsPahuteMesa | 150 kt, 4-layer tuff, SAC output |
| Per-cell velocity model | Integration.VelocityModelMaterial | Binary Vp/Vs/rho file, trilinear interp, aux fields |
| Thermal diffusion | Integration.ThermalDiffusion | Steady-state heat equation via PetscDS |
| Thermal expansion (THM) | Integration.ThermalExpansion | Thermoelastic stress, alpha_T*dT*Lz |

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
| Viscoelastic relaxation | Physics.ViscoelasticRelaxation |
| PetscDS BdResidual on cohesive | Physics.CohesiveBdResidual |

### DOES NOT WORK: Code Exists but No TSSolve Test

| Feature | What is Missing |
|---------|-----------------|
| Multiphase flow end-to-end | Buckley-Leverett or waterflood through TSSolve |
| Hydraulic fracture coupled solve | Full lubrication + deformation coupling |

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
| Viscoelastic callbacks | src/numerics/PetscFEViscoelastic.cpp | include/numerics/PetscFEViscoelastic.hpp |
| Absorbing BC | src/numerics/AbsorbingBC.cpp | include/numerics/AbsorbingBC.hpp |
| Gravity body force | src/numerics/PetscFEElasticityGravity.cpp | include/numerics/PetscFEElasticityGravity.hpp |
| Boundary conditions | src/numerics/BoundaryConditions.cpp | include/numerics/BoundaryConditions.hpp |
| Fault mesh | src/numerics/FaultMeshManager.cpp | include/numerics/FaultMeshManager.hpp |
| Cohesive kernel | src/physics/CohesiveFaultKernel.cpp | include/physics/CohesiveFaultKernel.hpp |
| Explosion physics | src/domain/explosion/ExplosionImpactPhysics.cpp | include/domain/explosion/ExplosionImpactPhysics.hpp |
| Seismometers | src/domain/seismic/SeismometerNetwork.cpp | include/domain/seismic/SeismometerNetwork.hpp |
| Fracture model | src/domain/geomechanics/FractureModel.cpp | include/domain/geomechanics/FractureModel.hpp |
| Gmsh I/O | src/io/GmshIO.cpp | include/io/GmshIO.hpp |
| Velocity model I/O | src/io/VelocityModelReader.cpp | include/io/VelocityModelReader.hpp |
| Thermal callbacks | src/numerics/PetscFEThermal.cpp | include/numerics/PetscFEThermal.hpp |
| Unit tests | tests/unit/ | -- |
| Integration tests | tests/integration/ | -- |
| Physics validation | tests/physics_validation/ | -- |
| Example configs | config/examples/ | -- |
| Runnable examples | examples/ | -- |
| Visualization scripts | scripts/plot_seismograms.py, scripts/plot_wavefield.py | -- |

## Runnable Examples

Eighteen examples in `examples/`, each with README.md and run.sh:

| # | Directory | Config | Physics |
|---|-----------|--------|--------|
| 01 | examples/01_uniaxial_compression | uniaxial_compression.config | Elastostatics, BCs |
| 02 | examples/02_explosion_seismogram | explosion_seismogram_quick.config | Elastodynamics, SAC |
| 03 | examples/03_elastoplastic_compression | elastoplastic_compression.config | Drucker-Prager |
| 04 | examples/04_locked_fault | locked_fault_compression.config | Cohesive cells |
| 05 | examples/05_punggye_ri_nuclear_test | punggye_ri_layered_quick.config | Layered, SAC, HDF5 |
| 06 | examples/06_gmsh_multimaterial | nuclear_twin_gmsh.config | Gmsh, per-region |
| 07 | examples/07_traction_bc | traction_bc.config | Per-face traction BC |
| 08 | examples/08_time_dependent_slip | time_dependent_slip.config | Slip ramp |
| 09 | examples/09_gasbuggy_1967 | gasbuggy_1967.config | 29 kt, 4-layer Lewis Shale |
| 10 | examples/10_gnome_1961 | gnome_1961.config | 3.1 kt, 4-layer Salado Salt |
| 11 | examples/11_sedan_1962 | sedan_1962.config | 104 kt, 3-layer alluvium |
| 12 | examples/12_degelen_mountain | degelen_mountain.config | 50 kt, 3-layer granite |
| 13 | examples/13_nts_pahute_mesa | nts_pahute_mesa.config | 150 kt, 4-layer tuff |
| 14 | examples/14_single_phase_flow | config.config | Darcy pressure diffusion |
| 15 | examples/15_viscoelastic_attenuation | config.config | GMB attenuation, seismograms |
| 16 | examples/16_scec_tpv5 | config.config | TPV5 dynamic rupture, slip-weakening |
| 17 | examples/17_velocity_model | config.config | Per-cell material from velocity file |
| 18 | examples/18_thermal_expansion | config.config | THM coupling, thermal stress |

## Roadmap: Features to Implement

Each item requires: PetscDS callbacks integrated into setupPhysics(), integration tests
through TSSolve, example config, and visualization. Source code in archive/src/ may provide
a starting point but must be rewritten to use the PetscDS callback pattern.

1. ~~Slipping fault convergence~~ DONE (full semi-smooth Newton Jacobian with off-diagonal slip direction derivatives and friction-normal coupling)
2. Multiphase flow end-to-end (Buckley-Leverett waterflood)
3. Full coupled hydraulic fracturing (lubrication + deformation)
4. ~~Viscoelastic attenuation~~ DONE (generalized Maxwell body, Q-factor memory variables)
5. ~~Thermal coupling~~ DONE (heat equation + THM Biot)
6. Radiation transport (advection-diffusion for fallout)
7. ~~Per-cell material from velocity model files~~ DONE
8. ~~SCEC TPV5 dynamic rupture benchmark~~ DONE (slip-weakening friction, initial fault stress, nucleation patch)
9. Multi-stage hydraulic fracturing with stress shadowing
10. Production forecasting through propped fracture

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.25.0 API signatures before calling any PETSc function.
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
15. GTEST_SKIP is ONLY for hardware-dependent tests (GPU, MPI rank count) or genuine crash bugs that would kill the test runner. A test that produces a zero solution is a FAILURE, not a skip. If SNES converges to zero, the physics setup is wrong -- fix it, do not hide it behind GTEST_SKIP.

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
