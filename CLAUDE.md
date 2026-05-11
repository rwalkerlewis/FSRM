# CLAUDE.md -- FSRM Persistent Context

This file provides all context needed to work on FSRM. Read it fully before making any changes.

## Project Summary

FSRM (Full Service Reservoir Model) is a C++17/PETSc/MPI coupled multiphysics simulator for nuclear explosion monitoring, seismic wave propagation, dynamic fault rupture, and coupled THM poroelasticity. MIT licensed.

Repository: github.com/rwalkerlewis/FSRM
Branch: main

## Documentation Layout

Read in this order when resuming work:

- `docs/SOLVER_STATE.md`: Current fault-solver state, pass/fail, matrix characterization, bottlenecks, change log.
- `docs/PYLITH_REFERENCE.md`: Verified PyLith architecture pins with file:line references.
- `docs/SESSION_RUNBOOK.md`: Env vars, commands, decision-branch template, report structure.
- `docs/SESSION_NN_REPORT.md`: Per-session changelog, numbered reports. Avoid reading these in sequence; prefer the three standing documents above and consult a specific report only when those direct you to one.

Secondary references:

- `docs/HISTORIC_NUCLEAR_FIDELITY.md`: standing truth about what each historic-nuclear pass shipped.
- `docs/HISTORIC_NUCLEAR_ROADMAP.md`: forward-looking roadmap with six fidelity axes. Pass-5 (in progress) marks axis 1 partial.
- `docs/archive/PYLITH_COMPATIBILITY.md`: Feature-parity wishlist (archived pass-12, no longer an active track).
- `docs/archive/LAGRANGE_FIX_STATUS.md`: Historical record of PETSc 3.25 architectural blockers (archived pass-12; superseded by `SOLVER_STATE.md`).
- `docs/archive/FAULT_TEST_REGRESSION_AUDIT.md`: Test inventory and history (archived pass-12; six fault tests now disabled by PR #119, see `SOLVER_STATE.md`).
- `docs/NUMERICAL_METHODS.md`, `docs/PHYSICS_MODELS.md`: physics reference.

### Codebase Size

After cleanup, the live codebase is:
- 63 source files (.cpp): ~54,000 lines
- 73 header files (.hpp): ~34,000 lines
- Total live: ~88,000 lines

Dead code (~45,000 lines across ~60 files) has been moved to `archive/src/` and `archive/include/`.

### Test Suite

116 registered tests. 110 pass, 6 fail honestly (no fake skips). Measured at Session 30 (default path, no env vars).

Known failures on default path (source of truth: `docs/SOLVER_STATE.md`):

- Physics.SCEC.TPV5: friction Jacobian port incomplete for slip-weakening + nucleation patch coupling.
- Physics.LockedFaultTransparency: PETSc 3.25 BdResidual on cohesive geometry; flaky between runs.
- Integration.PressurizedFractureFEM: hydrofrac rewire needed (separate scope from fault-solver bottleneck).
- Integration.DynamicRuptureSolve.PrescribedSlip: BdResidual on the Lagrange field does not fire on the cohesive geometry under PETSc 3.25. See `docs/SOLVER_STATE.md` and `docs/PYLITH_REFERENCE.md` for the current investigation state.
- Integration.SlippingFaultSolve: friction Jacobian port needed.
- Integration.SlipWeakeningFault: friction Jacobian port needed.

Three additional tests regress under the experimental fieldsplit path (`FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit`): `Integration.DynamicRuptureSolve.LockedQuasiStatic`, `Integration.DynamicRuptureSolve.LockedElastodynamic`, `Integration.TimeDependentSlip`. Under experimental the fault subset is 7 / 16. See `docs/SOLVER_STATE.md` "Extra failures under experimental path."

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
  ./fsrm -c ../examples/05_punggye_ri_nuclear_test/config_full.config \
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

addExplosionSourceToResidual has three branches:

1. Pre-pass-5 path: scalar moment rate from `RDPSeismicSource::psiDot` (COUPLED_ANALYTIC explosion_solve_mode) or `MuellerMurphySource::momentRate` (PROXY), applied as the trace of an isotropic moment tensor. This is the legacy KINEMATIC_RDP path.
2. Pass-4 distributed path: same scalar moment rate, but distributed over a multi-cell support ball (GAUSSIAN or UNIFORM_SPHERE) instead of a single cell. Selected via `[SOURCE_DISTRIBUTION] mode`.
3. Pass-5 dynamic-plastic path: under `[NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC`, a setup-time 1D solver records the full 6-component moment-rate tensor and cavity-radius history; the recorded history is interpolated and injected. Pass-6 added a `solver_kind` sub-key that selects between two backends: `CLOSED_FORM` (pass-5 RDP-driven byte-identical path; runs the existing closed-form `NearFieldExplosionSolver` and `RDPSeismicSource::momentRateTensor`) and `RADIAL_LAGRANGIAN` (the 1D radial Lagrangian elastoplastic shock solver in `src/domain/explosion/RadialLagrangian.cpp` with explicit CFL stepping, Wilkins AV, Drucker-Prager radial return, MG EOS, outgoing-characteristic outer BC, and surface-integral moment extraction at the elastic radius). Both backends write the pass-5 `near_field_history.csv`; `RADIAL_LAGRANGIAN` also writes `near_field_profile.h5` + `.xdmf` (HDF5 + ParaView wrapper) with per-snapshot radial state. Pass-7 (closing axis-1a) replaces the pass-6 ideal-gas inner-cavity placeholder with a Tillotson host-rock EOS (`include/domain/explosion/TillotsonEOS.hpp`, four parameter sets covering granite, tuff, salt, alluvium-placeholder), a first-principles Newton energy-partition solve for the cavity initial state at the radiation-to-hydrodynamic transition time (Zel'dovich-Raizer 1967 end-state approximation; vapor density at the solid-rock density), and Wilkins (1980) literature AV coefficients (`c_l = 0.06`, `c_q = 1.5`). The Sedan 1962 anchor lands at a 2.24x amplitude ratio relative to the closed-form RDP estimate, well inside the factor-5 envelope from the pass-7 spec. Pass-7 promotes `solver_kind = RADIAL_LAGRANGIAN` to the default for `DYNAMIC_PLASTIC`; the legacy `Sedan1962_Dynamic` historic-nuclear test fixture and `config/examples/sedan_1962_dynamic.config` are pinned to `solver_kind = CLOSED_FORM` so their pass-5 byte-identical assertions are preserved. Pass-8 (axis 1) makes the radiation-phase fidelity ladder explicit and selectable via a new `[NEAR_FIELD_SOURCE] radiation_phase` sub-key: `ZELDOVICH_RAIZER` (LOW, default, pass-7 closed-form end-state preserved byte-for-byte), `MARSHAK_GREY` (MED, new explicit grey radiation-diffusion + Newton T^4 coupling solver in `include/domain/explosion/MarshakRadiationDiffusion.hpp`), `MARSHAK_MULTIGROUP` (HIGH scaffold; throws clear runtime_error), `SN_TRANSPORT` (HIGHEST, named only). Pass-8 also ships an IRIS waveform V&V infrastructure (production SACReader at `include/io/SACReader.hpp`, comparison metrics at `include/diagnostics/WaveformComparison.hpp`, ObsPy refresh tool at `tools/waveform_vv/refresh.py`, cached waveform layout at `tools/waveform_vv/cache/`) anchored on Salmon 1964 with Chagan 1965 and Pokhran I 1974 cross-validation; the new CTest label `iris_validation` runs five integration gates that GTEST_SKIP cleanly when the cache is empty. See `docs/HISTORIC_NUCLEAR_FIDELITY.md` "Closed in pass 7" and "Closed in pass 8", `docs/HISTORIC_NUCLEAR_ROADMAP.md` axis 1 and the fidelity ladder section, `docs/EXPLOSION_IMPACT_PHYSICS.md` "Pass-7 inner-cavity physics" and "Pass-8 radiation transport phase", and `docs/WAVEFORM_VV.md`. Pass-9 (axis 1) advances the EOS and opacity ladders one rung each: `cavity_eos = TILLOTSON_TABULATED_PATCH` blends Tillotson with a Z-R partial-ionization plasma table (sin^2 in pressure across [5e10, 6e10] Pa) sourced from `tools/tabulated_data/tables/eos/<medium>_aneos.h5`; `opacity_model = TABULATED_PATCHED` blends the Z-R power law with Mihalas-Mihalas-corrected Rosseland/Planck means (sin^2 in temperature across [1e5, 1.26e5] K) sourced from `tools/tabulated_data/tables/opacity/<medium>_{rosseland,planck}.h5`. Both options ship with on-disk readers under `include/io/TabulatedData/TabulatedDataReader.hpp` (HDF5 with required `source_citation` metadata; out-of-table queries fall back to the analytic with one-time stderr warnings). Pass-9 also adds `operator_splitting = STRANG` (second-order Strang 1968 split between hydro and Marshak radiation; LIE remains the pass-8 byte-identical default) and a `radial_outer_radius_m` config field that exposes a direct outer-domain override (used by the pass-9 free-field gate triage to extend Salmon to 700 m). Five of eight pass-8 V&V gates were tightened to the spec target; three retain a documented residual (Salmon cavity factor 3 vs target 2; Marshak self-similar front factor 8 vs target 2; Pokhran mb ±0.4 vs target ±0.3). The Salmon Marshak example config moves to the pass-9 HIGH tier defaults; existing 27 historic tests continue to pass under their pinned configs (cavity_eos=TILLOTSON, opacity_model=POWER_LAW_ZR, operator_splitting=LIE). See `tools/tabulated_data/README.md` and the `Physics.TabulatedEOS.*`, `Physics.TabulatedOpacity.*`, `Physics.Marshak.OperatorSplittingConvergence`, `Unit.TabulatedDataReader` test families. Pass-10 (axis-1a closeout) promotes the three HIGHEST-tier scaffolds on the axis-1a ladders to working implementations: `MARSHAK_MULTIGROUP` is a working frequency-discretized counterpart to the pass-8 grey solver (G coupled tridiagonal solves per Newton iteration, default 16 log-spaced groups from 1e14 to 1e18 Hz, per-group analytic Mihalas-Mihalas + Kramers + Thomson opacity model in `include/domain/explosion/MultigroupOpacity.hpp` and the solver in `include/domain/explosion/MultigroupRadiationDiffusion.hpp`); `cavity_eos = TABULATED_FULL` is pure-tabulated EOS evaluation everywhere with a Tillotson safety net for out-of-table queries; `opacity_model = TABULATED_FULL` is pure-tabulated grey opacity (no Z-R fallback) with the dispatch matrix `(opacity_model, radiation_phase)` enforced at config time. Pass-10 also adds higher-order explicit time integrators `time_integrator = TVD_RK2` (Heun's method) and `RK3_SSP` (Shu-Osher 1988 SSP3) implemented as convex blends of the pass-9 explicit-Euler hydro operator. The Strang convergence-order test is instrumented at the substep level via the new `operator_splitting_convergence_diagnostic` flag. After pass-10 every cell on the axis-1a ladders has a working implementation; LOW/MED/HIGH defaults reproduce pass-N byte-for-byte under pinned configs. The 11 new physics-validation gates are in `tests/physics_validation/test_multigroup_radiation.cpp` (six gates) and `tests/physics_validation/test_pass10_integrators.cpp` (five gates). See `docs/AXIS_1A_FIDELITY_REPORT.md` for the canonical cross-pass gate-by-gate result and the named follow-up axes (axis-1b 3D source ball, axis-1c implicit time stepping for the diffusion solve, axis-1d 3D far-field FEM coupling). Pass-11 (axis 1c implicit diffusion + 549 m BC + ANEOS extension + 1b scaffold) closes the Marshak self-similar gate to factor 2 and the radiation-energy conservation gate to 2% by introducing the new `DiffusionTimeIntegrator` strategy (BackwardEuler / CrankNicolson / BDF2; gated by `time_integrator_diffusion` and the new HIGHEST-tier default BDF2). The 549 m free-field BC triage (sweep at radial_outer_radius_m = [700, 1000, 1500, 2000] m) shows the impedance BC at 700 m contaminates the gauge by factor ~3.4 at the closest two radii; pass-11 ships the Israeli-Orszag 1981 graded-damping sponge layer (gated by `sponge_layer_enabled`) as the BC fix. Extended ANEOS Hugoniot match gates (Marsh 1980 granite, McQueen 1970 salt) ship at 30% envelope; the 5% spec target requires a Tillotson parameter refit (axis-4). Axis-1b 3D source ball scaffolding lands: `Source3DBall.hpp` interface, `cavity_geometry` config knob (THREE_DIMENSIONAL throws pass-12 work), `docs/AXIS_1B_DESIGN.md` design stub. Pass-12 implements the 3D source ball on this scaffold. After pass-11 axis-1a has eight gates closed at spec target, three CavityRadius gates pinned awaiting axis-1b, three FarFieldBodyWaveMagnitude gates pinned awaiting axis-3, and the Hugoniot match gates pinned awaiting axis-4. After pass-11 axis-1a is closed for now. Pass-12 was housekeeping (documentation hygiene, per-event config relocation into `examples/N_<event>/config.config`, three new historic events DPRK 2017 / Lop Nor 1996 / Pokhran II 1998, showcase figure infrastructure for Sedan / Salmon / Punggye-ri, MPI standardization across all example run.sh scripts via `scripts/run_with_mpi.sh`); pass-13 implements axis-1b 3D source ball on the pass-11 scaffold. See `docs/HISTORIC_NUCLEAR_FIDELITY.md` "4k. Closed in pass 12" for the pass-12 closeout. Pass-13 is split into three slices because the full scope (TetGen integration, 3D DMPlex, 3D Drucker-Prager, asymmetric overburden, 3D radiation diffusion, surface-integral moment extraction, host delegation, MPI=4 Salmon end-to-end) is several months of senior-engineer work. Pass-13a (foundation, merged) ships the leaf data structures: TetGen pre-process tool at `tools/mesh_generation/build_source_ball_mesh.py`, 3D DMPlex mesh load + distribute via `include/domain/explosion/Source3DBallMesh.hpp`, and `Source3DBallImpl` skeleton (`include/domain/explosion/Source3DBallImpl.hpp`) replacing the pass-11 factory throw. Foundation `Source3DBallImpl::initialize` loads the mesh and validates geometry; `step()` throws with a clear pass-13b message; `getMomentTensor` returns zeros. Pass-13b (physics, this pass) lands 3D Drucker-Prager radial return (`include/domain/explosion/DruckerPrager3D.hpp`, Simo & Hughes 1998 sec 3.6 explicit projection with per-medium parameter sets), asymmetric overburden initial state (`include/domain/explosion/SourceBallOverburden.hpp`, Hoek & Brown 1980 K_0 default 0.5; Patton 1991 cavity-asymmetry rationale), cell-centred FV grey radiation diffusion (`include/domain/explosion/SourceBallRadiation3D.hpp` + `src/domain/explosion/SourceBallRadiation3D.cpp`, backward-Euler implicit assembly via PETSc Mat / Vec / KSP under the pass-12 parallel KSP convention, Mihalas & Mihalas 1984 secs 81-82), host-side `RadialLagrangianSolver -> Source3DBallImpl` delegation (`cavity_geometry = THREE_DIMENSIONAL` no longer throws), ConfigReader plumbing for six new `[NEAR_FIELD_SOURCE]` 3D sub-keys (`cavity_geometry`, `cavity_radius_m`, `outer_radius_m`, `mesh_path`, `overburden_K0`, `source_ball_radiation_discretization`), and the headline `Physics.Source3DBall.SphericalCellLevelEquivalenceVs1D` gate at near machine precision under spherical-symmetric loading. `getMomentTensor` and `getMomentRateTensor` continue to return zeros intentionally; surface-integral extraction is pass-13c. The `Source3DBallImpl::name()` tag bumps from `Source3DBallImpl_pass-13a_foundation_skeleton` to `Source3DBallImpl_v1_pass13b_physics`; `RadialLagrangianSolver::name()` reports `RadialLagrangianSolver+Source3DBallImpl_v1_pass13b_physics` when 3D delegation is active. Pass-13c (axis-1b validation slice + wavefield infrastructure + new examples) lands the surface-integral moment-tensor extraction (`Source3DBallImpl::recomputeMomentTensorFromSurface`, Day & McLaughlin 1991), the new `[OUTPUT]` config block (`wavefield_format`, `wavefield_cadence_steps`, `wavefield_fields`, `wavefield_output_directory`, `wavefield_basename`, `source_ball_3d_output_format`, `source_ball_3d_output_cadence_steps`; all default to NONE so the 32 historic-event tests stay byte-identical), the wavefield writer wired in `Simulator::MonitorFunction` (VTU + HDF5_XDMF with ParaView temporal-collection wrapper, Henderson 2007), the source-ball 3D per-cell snapshot writer, and two new academic verification anchors (`examples/44_lambs_problem/` with Eringen-Suhubi 1975 reference; `examples/45_layered_halfspace_explosion/` with Haskell-Thomson 1953/1950 reference). The headline pass-13c gate `Physics.Source3DBall.CLVDContentWithOverburden` measures CLVD content > 0.01 of isotropic with the CLVD axis aligned to vertical within 10 deg; on the cube_5tet fixture the measured ratio is 62.9 with cos(angle) = 0.999999. The pass-13b `name()` tag bumps from `Source3DBallImpl_v1_pass13b_physics` to `Source3DBallImpl_v1_pass13c_validation`. The end-to-end MPI=4 Salmon-with-overburden integration test (running to simulation completion) is deferred to a follow-up PR pending the Salmon-specific TetGen mesh. Cavity aspect ratio is reported via `Source3DBallImpl::cavityRadiusExtremes` (no gate; ~1.00 within solver noise in pass-13c, will populate asymmetric in pass-14). See `docs/HISTORIC_NUCLEAR_FIDELITY.md` "4n. Closed in pass 13c" for the canonical pass-13c closeout and `docs/AXIS_1B_DESIGN.md` for the slicing rationale.

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
DMPlexTSComputeIFunctionFEM. See docs/archive/LAGRANGE_FIX_STATUS.md for the historical decision trail.

Section B Option B from `docs/archive/PYLITH_COMPATIBILITY.md` is partially landed:
`addInteriorLagrangeResidual` (`src/core/Simulator.cpp:4180-4271`) zeros the Lagrange
residual at non-cohesive DOFs after `DMPlexTSComputeIFunctionFEM`, and the disjoint
loop inside `addCohesivePenaltyToJacobian` (`src/core/Simulator.cpp:4677-4720`) stamps
a canonical penalty-scaled diagonal on every interior (non-cohesive) Lagrange DOF.
The strict spec also calls for removing the PetscDS volume callback, but the volume
callback remains registered: empirically PETSc 3.25 needs at least one volume residual
on the Lagrange field for the BdResidual on cohesive cells to fire reliably, and
removing it regresses `Integration.DynamicRuptureSolve.LockedFaultElastodynamic`. The
manual residual zeroing overrides the small epsilon contribution at interior DOFs, so
the net interior behaviour matches the strict Section B form.

The prescribed-slip Cartesian vector is now copied from `cohesive_kernel_` into the
unified constants array in `setupPhysics` (`src/core/Simulator.cpp:2234-2249`). The
original code allocated `COHESIVE_CONST_PRESCRIBED_SLIP_X..Z` slots but never wrote
them, so `f0_prescribed_slip` always read a zero target jump. The push is gated on
`fault_mode_ == "prescribed_slip"` and uses the new
`CohesiveFaultKernel::getPrescribedSlip` accessor in
`include/physics/CohesiveFaultKernel.hpp`.

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
| Schema-anchor configs | config/ (default, complete_template, test_*) | -- |
| Per-physics templates | config/templates/ | -- |
| Per-event configs | examples/N_<event>/config.config | -- |
| Runnable examples | examples/ (41 directories pass-12) | -- |
| Showcase figure infra | tools/figures/figure_style.py | -- |
| MPI launcher | scripts/run_with_mpi.sh | -- |
| Visualization scripts | scripts/plot_seismograms.py, scripts/plot_wavefield.py | -- |

## Runnable Examples

Pass-12 relocated every per-event config into its example directory.
Each `examples/N_<event>/` is self-contained: `config.config`,
`README.md`, and `run.sh` live together; running `./run.sh` builds
output under `output/`. The runner sources
`scripts/run_with_mpi.sh`; override the rank count with the
`MPI_RANKS` env var (default 4).

| # | Directory | Physics |
|---|---|---|
| 01 | examples/01_uniaxial_compression | Elastostatics, BCs |
| 02 | examples/02_explosion_seismogram | Elastodynamics, SAC |
| 03 | examples/03_elastoplastic_compression | Drucker-Prager |
| 04 | examples/04_locked_fault | Cohesive cells |
| 05 | examples/05_punggye_ri_nuclear_test | Layered, SAC, HDF5 (Punggye-ri showcase) |
| 06 | examples/06_gmsh_multimaterial | Gmsh, per-region |
| 07 | examples/07_traction_bc | Per-face traction BC |
| 08 | examples/08_time_dependent_slip | Slip ramp |
| 09 | examples/09_gasbuggy_1967 | 29 kt, 4-layer Lewis Shale |
| 10 | examples/10_gnome_1961 | 3.1 kt, 4-layer Salado Salt |
| 11 | examples/11_sedan_1962 | 104 kt, 3-layer alluvium (Sedan showcase) |
| 12 | examples/12_degelen_mountain | 50 kt, 3-layer granite |
| 13 | examples/13_nts_pahute_mesa | 150 kt, 4-layer tuff |
| 14 | examples/14_single_phase_flow | Darcy pressure diffusion |
| 15 | examples/15_viscoelastic_attenuation | GMB attenuation, seismograms |
| 16 | examples/16_scec_tpv5 | TPV5 dynamic rupture, slip-weakening |
| 17 | examples/17_velocity_model | Per-cell material from velocity file |
| 18 | examples/18_thermal_expansion | THM coupling, thermal stress |
| 19-38 | examples/19_rainier_1957 through examples/38_lop_nor_1976 | 20 historic-nuclear events spanning 1957-2016 |
| 39 | examples/39_dprk_2017 | Pass-12: DPRK 2017 (~250 kt, granite, Mt. Mantap) |
| 40 | examples/40_lop_nor_1996 | Pass-12: Lop Nor 1996 (~5 kt, Chinese final test) |
| 41 | examples/41_pokhran_ii_1998 | Pass-12: Pokhran II Shakti-I (~20 kt thermonuclear) |

The Sedan, Salmon, and Punggye-ri showcase events also ship a
`figures/` directory with six per-event Python scripts that produce a
presentation-quality figure pack from simulation output (run
`./run_showcase.sh` in the example dir). See
`docs/FIDELITY_LADDER_GUIDE.md` for tier selection and
`tools/figures/README.md` for the shared style infrastructure.

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
10. Working examples are in `examples/N_<event>/` with `config.config`. Per-physics templates live in `config/templates/`. The top-level `config/` retains only schema-anchor configs (default, complete_template, test_*).
11. The executable is `fsrm`, not `fsrm_simulator`. Each example ships a `run.sh` that sources `scripts/run_with_mpi.sh` and runs from the example dir; override rank count with `MPI_RANKS`. Manual invocation: `./fsrm -c ../examples/N_<event>/config.config`.
12. No em dashes or contractions in code comments or documentation.
13. Verification tests must have quantitative pass/fail criteria with numerical tolerances.
14. Update CLAUDE.md and README.md after every session to reflect actual code state. Every claim must be backed by a specific test name or code reference.
15. GTEST_SKIP is ONLY for hardware-dependent tests (GPU, MPI rank count) or genuine crash bugs that would kill the test runner. A test that produces a zero solution is a FAILURE, not a skip. If SNES converges to zero, the physics setup is wrong -- fix it, do not hide it behind GTEST_SKIP.
16. Never truncate terminal output. Do not pipe `make`, `ctest`, `docker run`, `cmake`, or any other build / test / run command through `head`, `tail`, `grep`, `awk`, `sed`, or `cut`. Capture the full stdout / stderr stream to a log file (for example `... 2>&1 | tee /tmp/full.log`) and inspect the saved log directly with the `Read` tool when results are needed. Truncating output in-line hides build warnings, test names, SNES iterations, and stack traces that are diagnostic for failures.
17. No option tuning without diagnostic data. Before changing any PETSc option in the solver path (preconditioner type, factorization type, tolerances, damping, smoother, restart length, etc.), the session report must document a specific measurement (condition number, spectrum gap, KSP residual trajectory, FD vs hand-coded Jacobian difference) that justifies the specific option being changed. Speculative tuning loops without new measurement data are prohibited.
18. No Jacobian-adjacent modifications in a diagnostic session. A session whose prompt is labeled "diagnostic" must not modify `src/physics/`, `src/numerics/`, or the cohesive / bulk Jacobian registration code paths in `src/core/Simulator.cpp`. Diagnostic prints are allowed when guarded by `FSRM_S<N>_<TAG>=1` env-var gates.
19. Every session that changes solver behavior must append a one-line entry to the change log in `docs/SOLVER_STATE.md` and update the "Current test state" table if pass/fail counts changed.
20. Before starting a solver-related session, read `docs/SOLVER_STATE.md`. Before re-deriving any PyLith architectural fact, consult `docs/PYLITH_REFERENCE.md` and update it there rather than in a session report.
21. Every example under `examples/N_*/` must run end-to-end at `MPI_RANKS=4` (the production default) and produce at least one non-empty file under `examples/N/output/`. The `examples_runtime` CTest label asserts this on every shipped example. The full `MPI={1,2,4}` matrix runs when `FSRM_EXAMPLES_RUNTIME_FULL=ON` is passed at CMake configure time. New examples added under `examples/N_*/` must extend the gate, not bypass it. This rule ties to rule 15 (no fake skips): adding `GTEST_SKIP` or wrapping the smoke runner in a conditional that hides a real regression is prohibited.
22. Strict configuration validation runs by default in `Simulator::initializeFromConfigFile`. Adding a new `[SECTION]` or key to any config under `examples/` requires extending the registry in `src/io/ConfigValidator.cpp`. Documentation-grade templates in `config/` and `config/templates/` carry an explicit `[META] strict_validation = false` opt-out (always logged on stderr) so the schema-anchor templates can enumerate aspirational keys without rejecting at runtime. The `Unit.ConfigStrictValidation.AllShippedExampleConfigsValidate` gate ensures every shipped example validates. See `docs/CONFIGURATION_VALIDATION.md`.

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
