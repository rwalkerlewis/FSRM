# FSRM Session Prompt: CTBTO Abstract Readiness

## READ FIRST

Read `CLAUDE.md` in the repository root completely before starting. All rules there are binding.

Hard rules (violating any of these wastes the session):
- NEVER change DS/BC ordering in setupFields()
- NEVER modify existing callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp
- NEVER modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS
- Build and test in Docker using Dockerfile.ci. Always.
- Check PETSc 3.22.2 API signatures before calling any PETSc function: `grep -rn "FunctionName" /opt/petsc-3.22.2/include/`
- All existing tests must continue to pass after every change.

## Objective

Get the `local_fix` branch to a state where the following abstract is defensible, then PR to main.

> FSRM demonstrates AI-assisted development of scientific software for nuclear explosion monitoring. The code implements Mueller-Murphy seismic source models for underground tests, capturing cavity formation, spall, chimney collapse, and seismic wave generation essential for treaty verification. Atmospheric detonation models include Brode fireball evolution, Sedov-Taylor blast wave propagation, and EMP effects. Built on PETSc for parallel computing, the simulator handles coupling between elastodynamics, plasticity, fracture mechanics, and radiation transport. The tool validates against SCEC earthquake benchmarks and generates synthetic seismograms comparable to real IMS station data.

## Code Review Findings (Do Not Rediscover These)

An external review of the local_fix branch has been completed. These findings are ground truth. Do not re-derive them.

### Verified Correct by Inspection

1. **FormFunction / FormJacobian** -- Proper PETSc FEM assembly. Global-to-local scatter, DMPlexInsertBoundaryValues, DMPlexTSComputeIFunctionFEM/IJacobianFEM, injection/explosion source addition, local-to-global accumulation. This is correct.

2. **PetscFEElasticity callbacks** -- f1_elastostatics computes sigma = lambda*tr(eps)*I + 2*mu*eps from u_x (displacement gradient). Mathematically correct Hooke's law.

3. **Elastodynamics + TSALPHA2** -- f0_elastodynamics = rho * u_t (mass term), f1_elastodynamics delegates to f1_elastostatics (elastic stress). g0 = rho * u_tShift * I. With TSALPHA2, this correctly solves the second-order wave equation M*a + K*u = 0. The TSSetIFunction form is correct for TSALPHA2 in PETSc 3.22.

4. **AbsorbingBC** (src/numerics/AbsorbingBC.cpp) -- Lysmer-Kuhlemeyer absorbing BC. f0 decomposes velocity into normal (P-wave impedance rho*cp) and tangential (S-wave impedance rho*cs) components. g0 Jacobian is correct. Registered as DM_BC_NATURAL per-face in setupBoundaryConditions. This is a correct first-order absorbing BC.

5. **PetscFEPoroelasticity** -- All 4 Biot coupling blocks registered: f0/f1 pressure, f0/f1 displacement, g0_pp, g3_pp, g1_pu, g2_up, g3_uu. Correct for fully coupled Biot poroelasticity.

6. **Mueller-Murphy source** -- Psi function with 3 yield ranges. Cavity radius Rc = 12 * W^(1/3). Scalar moment M0 = 4*pi*rho*vp^2*Rc^3. mb = 4.45 + 0.75*log10(W). RDP spectral shape omega^-2. Reasonable approximation of Mueller-Murphy (1971), not exact but defensible.

7. **SeismometerNetwork** -- DMInterpolationEvaluate at station coordinates, time-displacement trace storage, velocity/acceleration from finite differences, SAC file output with proper header. Real end-to-end pipeline.

8. **SphericalCavitySource stress tensor** -- Radial/tangential cavity stresses mapped to Cartesian via sigma = sigma_tt*I + (sigma_rr - sigma_tt)*n*n^T. Correct for isotropic source.

9. **Atmospheric explosion** -- NuclearAirburstEffects.cpp has Sedov-Taylor blast radius R = 1.15*(E/rho)^0.2*t^0.4, Brode fireball scaling, EMP E1/E2/E3 with configurable parameters.

### Known Bug: Explosion Source Injection

**Location:** `Simulator::addExplosionSourceToResidual()` (src/core/Simulator.cpp, line ~2491)

**Problem:** The current code distributes M[0], M[1], M[2] as direct body forces to nodal displacement DOFs:
```cpp
closure[n * dim + d] += -M[d] / vol;
```
This is NOT a proper moment tensor source. A correct FEM moment tensor source requires contracting M_ij with shape function gradients: F_i^node = integral(dN/dx_j * M_ij dV). What exists is a uniform pressure force at every node of the source cell.

**Impact:** For an isotropic explosion (M[0]=M[1]=M[2]=scale, M[3-5]=0), this produces outward displacement in all three Cartesian directions at every node equally. This is qualitatively correct (radial expansion) but quantitatively wrong: amplitude, radiation pattern, and frequency content will be incorrect compared to the proper formulation. For conference demo purposes, this may be acceptable if the goal is "nonzero seismogram with correct travel time." For publication, it is not.

**Fix (if time permits):** Replace the body force injection with a proper equivalent nodal forces computation:
1. Tabulate shape function gradients dN_a/dx_j at the source element centroid (use DMPlexComputeCellGeometryFEM to get the inverse Jacobian, then evaluate basis gradients).
2. For each node a of the source cell: F_i^a = sum_j (dN_a/dx_j * M_ij) * vol
3. Insert these forces via DMPlexVecSetClosure.

This is a self-contained fix in addExplosionSourceToResidual. It does not affect any other code.

**Skip condition:** If the existing body-force injection produces nonzero seismograms at stations with correct P-wave travel time and amplitude decay with distance, it is sufficient for the CTBTO workshop abstract. Fix the moment tensor formulation later.

### Architectural Note: BC Conflicts

When `config.side_bc == "roller"` AND `config.absorbing_bc_enabled == true`, the same boundary faces may get both DM_BC_ESSENTIAL (roller) and DM_BC_NATURAL (absorbing) registrations. PETSc will apply the essential BC, masking the absorbing BC. For explosion seismograms, the correct config is:

```ini
[SIMULATION]
side_bc = free
bottom_bc = free
top_bc = free
absorbing_bc_enabled = true
absorbing_bc_x_min = true
absorbing_bc_x_max = true
absorbing_bc_y_min = true
absorbing_bc_y_max = true
absorbing_bc_z_min = true
absorbing_bc_z_max = false   # free surface on top
```

Verify that the explosion seismogram config and test configs use this pattern. If they use roller + absorbing, the absorbing BCs are silently ignored and waves reflect.

## Development Plan

### Phase 0: Build, Baseline, Fix Compilation

Build the project in Docker. Document which tests pass and which fail.

```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc) 2>&1 | tee /workspace/build.log'
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure 2>&1 | tee /workspace/test_baseline.txt
```

Fix any compilation errors. Fix any test infrastructure failures (missing files, bad paths, timeouts). DO NOT fix physics failures yet. Commit: "Phase 0: baseline build, N/M tests passing"

After this phase, you know the ground truth. Report the full test results before proceeding.

---

### Phase 1: Explosion Seismogram Pipeline End-to-End

This is the highest-priority deliverable. Everything else is secondary.

**Goal:** `Integration.ExplosionSeismogram` passes with these quantitative criteria:
1. SAC files are produced for each station and component
2. SAC data contains nonzero samples
3. Closer station has larger peak amplitude than farther station
4. P-wave first arrival time at each station is within 20% of (distance / vp)

**Likely failure modes and fixes:**

A) **Absorbing BCs not active** -- Check the test config. If it uses `side_bc = roller` alongside `absorbing_bc_enabled = true`, the absorbing BCs are masked. Change to `side_bc = free`.

B) **use_fem_time_residual_ is false** -- Means no PetscDS callbacks were registered. Check that config has `enable_geomechanics = true` and `enable_elastodynamics = true` and `solid_model = ELASTIC` and `fluid_model = NONE`.

C) **TSALPHA2 diverges** -- Check SNES convergence. May need `rtol = 1e-4` and `atol = 1e-6` with larger `max_nonlinear_iterations`. Can also try smaller dt_initial.

D) **Explosion source produces zero displacement** -- Check that `explosion_cell_ >= 0` (source point inside mesh). Check that `ExplosionCoupling::cavity.equivalentMoment()` returns nonzero. Check that the moment rate function `(M0/tau) * exp(-t/tau)` is nonzero during the simulation window.

E) **Seismometers not sampling** -- Check that `seismometers_->initialized_` is true. DMInterpolation setup requires station coordinates to lie inside the mesh domain. If stations are outside the bounding box, DMInterpolationEvaluate will silently skip them.

F) **SAC files not written** -- Check `seismometers_->writeAllTraces()` is called after `sim.run()`. If it is only called in the destructor, ensure the Simulator goes out of scope before test assertions.

**Files to examine/modify:**
- `tests/integration/test_explosion_seismogram.cpp` (the test)
- `config/examples/explosion_seismogram.config` (or inline config in test)
- `src/core/Simulator.cpp` (addExplosionSourceToResidual, MonitorFunction seismometer sampling)
- `src/domain/seismic/SeismometerNetwork.cpp` (setup, sample, writeSAC)

**Verification:** `ctest -R ExplosionSeismogram --output-on-failure` passes.

---

### Phase 2: Absorbing BC Energy Test

**Goal:** Verify absorbing BCs actually absorb waves (not just compile).

Create or verify `tests/physics_validation/test_absorbing_bc.cpp` includes an energy decay test:

1. 3D box, granite properties (lambda=30e9, mu=25e9, rho=2650).
2. Initial Gaussian displacement pulse at center.
3. Absorbing BCs on all 6 faces.
4. Run for time = 2 * (L/2) / cp (enough for wave to cross the domain twice).
5. Compute total kinetic energy: E = 0.5 * integral(rho * |v|^2 dV). Use VecNorm on velocity estimate or VecDot(U_t, M*U_t).
6. Assert: E(t_final) < 0.1 * E(t_peak). If energy is not absorbed, absorbing BCs are broken.

If the existing absorbing BC test only checks callback values (not full pipeline energy absorption), extend it.

**Verification:** `ctest -R AbsorbingBC --output-on-failure` passes.

---

### Phase 3: Mueller-Murphy Verification

**Goal:** Confirm the source physics produce correct numbers in isolation (no PDE solve needed).

Create `tests/physics_validation/test_mueller_murphy.cpp` if it does not exist:

1. **Corner frequency vs yield:** 1 kt granite 300m depth. Compute fc. Check against Mueller-Murphy 1971: fc should be in range 1-10 Hz. The psi function gives psi = 16.2/W^0.33 = 16.2 for 1 kt, and fc = psi * (rho_ref/rho)^(-1/3) / depth. For 300m depth in granite: fc ~ 16.2 / 300 ~ 0.054 Hz? That seems low. Actually the formula in the code is `corner_frequency = psi * pow(rho/2650, -1/3) / depth` which for 1 kt gives psi=16.2, so fc = 16.2 * 1.0 / 300 = 0.054 Hz. This seems too low by a factor of ~100. The actual Mueller-Murphy corner frequencies for 1 kt are 1-5 Hz. **This may be a units or formula bug.** Investigate.

2. **mb-yield relation:** Compute mb for W = 1, 10, 100, 250, 1000 kt. The code uses mb = 4.45 + 0.75*log10(W). Check:
   - 1 kt: mb = 4.45 (reasonable, observed ~4.0-4.5 for small tests)
   - 250 kt: mb = 4.45 + 0.75*2.398 = 6.25 (observed DPRK 2017: ~6.3, good)
   - 1000 kt: mb = 4.45 + 2.25 = 6.70 (reasonable)

3. **Scalar moment:** For 250 kt granite, Rc = 12 * 250^(1/3) = 12 * 6.3 = 75.6 m. M0 = 4*pi*2700*5500^2*75.6^3 = 4*pi*2700*3.025e7*4.32e5 = ~4.4e16 N*m. This is the right order of magnitude for a 250 kt explosion.

**Key investigation:** The corner frequency formula may have a bug. Mueller-Murphy (1971) corner frequencies for 1 kt hard rock are 2-5 Hz, not 0.05 Hz. If fc is 100x too low, the source time function is 100x too slow, and seismograms will have wrong frequency content. Check whether `depth` in the formula should be cavity radius (Rc ~ 12m for 1 kt) instead of depth of burial (300m). In the actual M-M paper, the corner frequency scales with Rc, not burial depth.

Register test as `Physics.MuellerMurphy` in CMakeLists.txt.

**Verification:** `ctest -R MuellerMurphy --output-on-failure` passes.

---

### Phase 4: Atmospheric Explosion Verification

**Goal:** Verify NuclearAirburstEffects produces correct Sedov-Taylor, Brode, and EMP values.

Create `tests/physics_validation/test_atmospheric_explosion.cpp`:

1. **Sedov-Taylor blast radius:** 20 kt (E = 20e3 * 4.184e12 J = 8.368e16 J) in standard atmosphere (rho_0 = 1.225 kg/m^3). R(0.01s) = 1.15 * (E/rho)^0.2 * t^0.4. Compute and verify within 5%.

2. **Brode fireball max radius:** R_max ~ 90 * W^0.4 m. For 20 kt: R_max ~ 90 * 20^0.4 ~ 90 * 3.31 ~ 298 m. Verify the code returns something in range 200-400 m.

3. **EMP E1 peak field:** Verify empE1Peak() returns a value in range 10-100 kV/m for 100 kt at 100 km. This is parametric, depends on configurable e0_vm.

4. **Overpressure:** Verify peak overpressure at known distances follows Glasstone-Dolan scaling within an order of magnitude.

Register as `Physics.AtmosphericExplosion`.

**Verification:** `ctest -R AtmosphericExplosion --output-on-failure` passes.

---

### Phase 5: Near-Field Explosion Phenomenology

**Goal:** Verify NearFieldExplosion and NuclearSourceParameters produce physically reasonable cavity/damage/spall numbers.

Create `tests/physics_validation/test_near_field_explosion.cpp`:

1. **Cavity radius:** 1 kt granite at 300m. Code gives Rc = 12 * 1^(1/3) = 12 m. Verify this is in the empirical range 8-15 m/kt^(1/3) for hard rock.

2. **Crushed zone radius:** Code gives 3 * Rc = 36 m. Verify > Rc.

3. **Fractured zone radius:** Code gives 10 * Rc = 120 m. Verify > crushed zone.

4. **SphericalCavitySource displacement:** At r = 100 m, t = 0.1 s, verify nonzero radial displacement with correct sign (outward).

Register as `Physics.NearFieldExplosion`.

**Verification:** `ctest -R NearFieldExplosion --output-on-failure` passes.

---

### Phase 6: Verify Existing Physics Validation Tests

Run all existing physics validation tests and fix any failures:

- `Physics.ElastostaticsPatch` -- elastostatics with known patch test
- `Physics.TerzaghiConsolidation` -- Biot poroelastic column
- `Physics.LambsProblem` -- point force on halfspace
- `Physics.GarvinsProblem` -- buried explosion Green's function
- `Physics.GravityLithostatic` -- gravity body force
- `Physics.LithostaticStress` -- lithostatic initial stress
- `Physics.SCEC.TPV5` -- dynamic rupture benchmark

For each failing test, determine if the failure is:
- **Config issue** (wrong BC combination, wrong field indices) -- fix config
- **Physics bug** (wrong callback, wrong constants) -- fix in a NEW file per CLAUDE.md rules
- **Missing feature** (needs code that does not exist) -- document and skip

**Verification:** Maximum number of physics tests passing. Document any that remain failing with a one-line explanation.

---

### Phase 7: Plasticity Assessment

CLAUDE.md Rule #7 says plasticity is a stub. The abstract says "plasticity." Determine the truth.

1. Read `src/domain/geomechanics/PlasticityModel.cpp` (1008 lines). Determine if it:
   a. Has a yield surface evaluation (Drucker-Prager, Mohr-Coulomb, von Mises)
   b. Has a return-mapping algorithm
   c. Is wired into a PetscDS callback (f1 that modifies the stress tensor)
   d. Or is standalone code that never gets called during a PDE solve

2. If it is standalone (never called from FormFunction), the abstract claim is "implementation exists but is not coupled to the FEM solver." This is defensible for a workshop abstract if stated honestly.

3. If it IS wired in, write a simple uniaxial compression test exceeding the yield stress and verify plastic strain > 0.

Report findings. Do not spend more than 30 minutes on this phase.

---

### Phase 8: README Cleanup

The main branch README claims volcano modeling, tsunami modeling, ocean physics, SeisSol-compatible DG, ResFrac-equivalent hydraulic fracturing, GPU acceleration with specific speedup numbers, and FNO neural operator solvers. None of these have passing tests.

1. Create a new section "Verified Capabilities (All Tests Pass)" listing ONLY features backed by green tests.
2. Move everything else to "Implemented (Not Yet Verified)" or "Planned."
3. Remove specific performance numbers (GPU speedups, scaling curves) unless backed by benchmark data in the repo.
4. Remove "ResFrac-Equivalent" section header. Replace with "Hydraulic Fracturing (Prototype)."
5. Remove volcano, tsunami, ocean physics from the main feature list. They can stay in docs/ as aspirational.

---

### Phase 9: PR to Main

After all phases complete:

1. Run full test suite: `ctest --output-on-failure`. Report results.
2. Create `docs/TEST_RESULTS.md` with every test name and PASS/FAIL.
3. Update CLAUDE.md "What Works" and "What Does NOT Work" to match reality.
4. Commit: "CTBTO abstract readiness: verified explosion monitoring pipeline"
5. Push local_fix and create PR:

```bash
gh pr create \
  --base main \
  --head local_fix \
  --title "Verified multi-physics explosion monitoring pipeline" \
  --body "Verified capabilities with passing tests:
- Elastodynamic wave propagation (TSALPHA2 generalized-alpha)
- Mueller-Murphy seismic source model
- Absorbing boundary conditions (Lysmer-Kuhlemeyer)
- End-to-end explosion seismogram pipeline (source -> FEM solve -> DMInterpolation -> SAC output)
- Biot poroelastic coupling (Terzaghi verification)
- Atmospheric explosion physics (Sedov-Taylor, Brode, EMP)
- Near-field explosion phenomenology (cavity, damage, spall)
- Cohesive fault mechanics with friction laws

Known limitations documented in CLAUDE.md and README.

Test results: docs/TEST_RESULTS.md"
```

Wait for CI. Fix any failures. Do not merge until green.

## Priority If Time-Limited

If you cannot complete all phases, do them in this order:
1. Phase 0 (build baseline) -- mandatory
2. Phase 1 (explosion seismogram) -- the core deliverable
3. Phase 3 (Mueller-Murphy, especially the corner frequency bug investigation)
4. Phase 2 (absorbing BC energy test)
5. Phase 6 (existing test fixes)
6. Everything else

## General Rules

1. Every new .cpp test file needs an entry in tests/CMakeLists.txt and a CTest registration.
2. Every new test must have quantitative pass/fail criteria with numerical tolerances.
3. New PetscDS callbacks go in NEW files. Do not modify existing callback files.
4. Build and test in Docker. Always.
5. Commit after each passing phase.
6. If a PETSc function does not exist in 3.22.2 headers, find the equivalent. Do not guess.
7. Do not fake passing tests. An honest failure with documented reason is better.
8. The explosion seismogram pipeline (Phase 1) is worth more than all other phases combined. Prioritize accordingly.