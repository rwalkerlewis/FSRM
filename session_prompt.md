# FSRM Session Prompt: local_fix -- Fix Everything, Test Everything, PR to Main

## READ FIRST

Read `CLAUDE.md` in the repository root completely before starting. All rules there are binding.

Hard rules:
- NEVER change DS/BC ordering in setupFields()
- NEVER modify existing callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp
- NEVER modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS
- Build and test in Docker using Dockerfile.ci. Always.
- Check PETSc 3.22.2 API signatures: `grep -rn "FunctionName" /opt/petsc-3.22.2/include/`
- All existing tests must continue to pass after every change.
- No GTEST_SKIP for tests that should work. If it should work, make it work. If it genuinely cannot work yet, delete the test and document why.

## What This Prompt Demands

Fix every known bug. Write tests for every feature the abstract claims. Wire every test to CI. Make the full suite pass. PR to main. No hedging. No "if time permits." No "skip condition." Do the work.

## External Code Review Findings

Two full code reviews have been completed. These findings are verified. Do not rediscover them.

### What Is Correct (Do Not Touch)

1. FormFunction / FormJacobian -- correct PETSc FEM assembly
2. PetscFEElasticity f0/f1/g3 -- correct Hooke's law
3. Elastodynamics + TSALPHA2 -- correct second-order wave equation
4. AbsorbingBC -- correct Lysmer-Kuhlemeyer
5. PetscFEPoroelasticity -- correct Biot coupling (4 blocks)
6. PetscFEElasticityAux -- correct aux-field callbacks, correct gravity sign
7. setupAuxiliaryDM / populateAuxFieldsByDepth -- correct aux field setup
8. SeismometerNetwork -- correct DMInterpolation, SAC output, field lookup
9. CohesiveFaultKernel prescribed slip -- correct constraint enforcement
10. Sedov-Taylor blast radius formula -- correct
11. SphericalCavitySource stress tensor -- correct radial/tangential to Cartesian mapping

### What Is Real Code But Not FEM-Coupled

1. **PlasticityModel** (src/domain/geomechanics/PlasticityModel.cpp, 1008 lines) -- Real Drucker-Prager, von Mises, Mohr-Coulomb, Cap model return-mapping. Never referenced from Simulator or any PetscDS callback. Standalone only.

2. **NearFieldExplosion** (src/domain/explosion/NearFieldExplosion.cpp, 822 lines) -- Real 1D Lagrangian solver with Mie-Gruneisen EOS, damage evolution, spall detection. Not coupled to PETSc FEM solver. Standalone only.

### Tests That Are Fake

These "tests" currently pass CI but test nothing:

1. **test_scec_tpv5.cpp** -- Two GTEST_SKIP calls. Zero assertions. Pure placeholder.
2. **test_boundary_conditions.cpp** -- Instantiates objects, calls SUCCEED(). Zero assertions.
3. **test_mms_elasticity.cpp** -- Key test skipped because GeomechanicsKernel uses "strong-form placeholder."
4. **test_mms_wave_propagation.cpp** -- Key test skipped because ElastodynamicsKernel uses "strong-form placeholder."
5. **test_analytical_solutions.cpp** -- 1 assertion, 1 skip.
6. **test_thermoporoelastic.cpp** -- 1 assertion.

## Bugs to Fix

### BUG 1: Explosion Seismogram Test Timing

**File:** tests/integration/test_explosion_seismogram.cpp, SetUp()

**Problem:** Simulation runs 0.05 seconds. P-wave velocity = 5213 m/s. NEAR station is 500m from source. P-wave reaches NEAR at 0.096s. The wave never arrives. SAC data will be zeros. Tests will fail.

**Fix:**
- Set end_time = 0.5
- Set max_timesteps = 200
- Set dt_initial = 0.0025
- Move NEAR to (2100, 2000, 2000) = 100m from source
- Move FAR to (2500, 2000, 2000) = 500m from source
- P-wave arrives at NEAR at 0.019s, at FAR at 0.096s. Both within 0.5s window.

### BUG 2: No Absorbing BCs in Explosion Test

**File:** tests/integration/test_explosion_seismogram.cpp, SetUp()

**Problem:** Config uses `sides = roller` with no `[ABSORBING_BC]` section. Waves reflect off all faces.

**Fix:** Replace the boundary section:
```ini
[BOUNDARY_CONDITIONS]
bottom = free
sides = free
top = free

[ABSORBING_BC]
enabled = true
x_min = true
x_max = true
y_min = true
y_max = true
z_min = true
z_max = false
```

When `sides = roller` AND `absorbing_bc_enabled = true`, the Dirichlet BC wins and absorbing is silently ignored. You must use `sides = free` for absorbing BCs to fire.

### BUG 3: Mueller-Murphy Corner Frequency Divides by Burial Depth

**File:** src/domain/explosion/ExplosionImpactPhysics.cpp, line ~189

**Problem:** `corner_frequency = psi / depth_of_burial`. Gives fc = 0.054 Hz for 1 kt at 300m. Should be fc = psi / cavity_radius = 16.2 / 12 = 1.35 Hz.

**Fix:**
```cpp
double Rc = source_params.cavity_radius(density);
corner_frequency = psi * std::pow(density / 2650.0, -1.0 / 3.0) / std::max(1.0, Rc);
```

### BUG 4: FEM Source Uses Wrong Scalar Moment

**File:** src/core/Simulator.cpp, addExplosionSourceToResidual()

**Problem:** Uses `SphericalCavitySource::equivalentMoment()` = (4/3)*pi*Rc^3*P0. For 250 kt this gives M0 = 7.2e17 N*m. The Mueller-Murphy M0 = 4*pi*rho*vp^2*Rc^3 = 4.4e16 N*m. Source is 16x too strong.

**Fix:** 
1. Add `double yield_kt` and `double depth_of_burial` as stored members in `ExplosionCoupling`.
2. Store them in `configureUndergroundNuclear`.
3. In `addExplosionSourceToResidual`, compute M0 from Mueller-Murphy formula:
```cpp
NuclearSourceParameters nsp;
nsp.yield_kt = explosion_->yield_kt;
nsp.depth_of_burial = explosion_->depth_of_burial;
MuellerMurphySource mm;
mm.setMediumProperties(explosion_->rho, explosion_->vp, explosion_->vs);
mm.setParameters(nsp);
double mr = mm.momentRate(elapsed);
```

### BUG 5: Explosion Source Injection Is Not a Proper Moment Tensor

**File:** src/core/Simulator.cpp, addExplosionSourceToResidual(), line ~2555

**Problem:** Distributes M[d]/vol as uniform body forces. A proper FEM moment tensor source requires: F_i^node = sum_j(dN_a/dx_j * M_ij) * vol.

**Fix:**
1. Use `DMPlexComputeCellGeometryFEM(dm, explosion_cell_, NULL, NULL, J, invJ, &detJ)` to get the inverse Jacobian.
2. Get the PetscFE from the DM. Tabulate basis function gradients at the cell centroid using `PetscFEGetTabulation`.
3. For each basis function a and each component i: `F_i^a = sum_j(dN_a/dx_j * M_ij) * detJ * weight`.
4. Insert via DMPlexVecSetClosure.

Check PETSc 3.22.2 API for exact signatures. The key functions are `PetscFECreateTabulation` or `PetscFEGetCellTabulation`.

### BUG 6: Constants Array Conflict with Gravity

**File:** src/core/Simulator.cpp, setupPhysics(), line ~1686

**Problem:** `unified_constants[0] = config.gravity` overwrites lambda when aux callbacks are enabled. AbsorbingBC fallback reads constants[0] as lambda.

**Fix:** When absorbing BCs and gravity are both enabled, store cp and cs explicitly as additional constants:
```cpp
if (config.absorbing_bc_enabled && !material_props.empty()) {
    double E = material_props[0].youngs_modulus;
    double nu = material_props[0].poisson_ratio;
    double rho = material_props[0].density;
    double lam = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));
    // Store wave speeds where AbsorbingBC can find them if aux fields fail
    // Use slots after cohesive constants
    unified_constants[ABSORBING_CP_SLOT] = std::sqrt((lam + 2 * mu) / rho);
    unified_constants[ABSORBING_CS_SLOT] = std::sqrt(mu / rho);
}
```

Or simpler: make AbsorbingBC always require aux fields when gravity is enabled. Add an assertion in setupPhysics that absorbing + gravity requires aux fields.

## Tests to Write

Every test listed below is mandatory. Each test must have quantitative pass/fail criteria with numerical tolerances. Each test must be registered in tests/CMakeLists.txt and wired to CTest with a label.

### NEW TEST 1: Mueller-Murphy Physics (tests/physics_validation/test_mueller_murphy.cpp)

Tests:
1. Corner frequency for 1 kt granite at 300m: 0.5 < fc < 10.0 Hz
2. Corner frequency for 250 kt granite at 800m: 0.01 < fc < 1.0 Hz
3. Corner frequency decreases with yield (fc(10kt) < fc(1kt))
4. mb for 1 kt: 4.0 < mb < 5.0
5. mb for 250 kt: 5.8 < mb < 6.8 (observed DPRK 2017: 6.3)
6. mb for 1000 kt: 6.0 < mb < 7.5
7. Scalar moment for 1 kt: 1e13 < M0 < 1e16 N*m
8. Scalar moment for 250 kt: 1e15 < M0 < 1e18 N*m
9. Moment rate at t=0 is nonzero
10. Moment rate at t=10*rise_time is < 1% of peak
11. RDP spectral amplitude at omega << omega_c equals M0 (flat low-frequency plateau)
12. RDP spectral amplitude at omega >> omega_c decays as omega^-2

Register as `Physics.MuellerMurphy` with label `physics_validation`.

### NEW TEST 2: Atmospheric Explosion Physics (tests/physics_validation/test_atmospheric_explosion.cpp)

Tests:
1. Sedov-Taylor blast radius at t=0.01s for 20 kt: check against R = 1.15*(E/rho)^0.2*t^0.4, tolerance 5%
2. Sedov-Taylor blast radius at t=0.1s for 20 kt: same formula, tolerance 5%
3. Blast radius increases with time (R(0.1) > R(0.01))
4. Blast radius increases with yield (R(100kt) > R(20kt) at same time)
5. Fireball max radius for 1 kt: 50 < R_max < 200 m
6. Fireball max radius for 20 kt: 150 < R_max < 500 m
7. Fireball max radius for 1000 kt: 500 < R_max < 3000 m
8. Peak overpressure at 1 km for 20 kt surface burst: 10 < P < 500 kPa
9. Peak overpressure decreases with distance
10. EMP E1 peak field for 100 kt: 1e3 < E < 1e5 V/m at 100 km
11. Seismic magnitude from airburst is nonzero

Register as `Physics.AtmosphericExplosion` with label `physics_validation`.

### NEW TEST 3: Near-Field Explosion Phenomenology (tests/physics_validation/test_near_field_explosion.cpp)

Tests:
1. Cavity radius for 1 kt granite: 8 < Rc < 20 m
2. Cavity radius for 250 kt granite: 50 < Rc < 120 m
3. Cavity radius scales as W^(1/3) within 10%
4. Crushed zone radius > cavity radius
5. Fractured zone radius > crushed zone radius
6. SphericalCavitySource displacement at r=100m, t=0.1s is nonzero and positive (outward)
7. SphericalCavitySource displacement decays with distance: u(200m) < u(100m)
8. SphericalCavitySource equivalent moment is nonzero and positive
9. SphericalCavitySource stress: sigma_rr < 0 (compressive), sigma_tt > 0 (tensile)
10. NearFieldExplosionSolver initializes without error
11. NearFieldExplosionSolver step produces nonzero stress state
12. Spall detection triggers for tensile failure (pressure < -tensile_strength)

Register as `Physics.NearFieldExplosion` with label `physics_validation`.

### NEW TEST 4: Plasticity Return Mapping (tests/physics_validation/test_plasticity_returnmap.cpp)

Tests:
1. DruckerPrager yield function returns negative for stress inside yield surface
2. DruckerPrager yield function returns positive for stress outside yield surface
3. DruckerPrager return mapping brings stress back to yield surface (yield function = 0 after return)
4. DruckerPrager return mapping produces nonzero plastic strain increment
5. VonMises yield function: sigma_eq < sigma_y returns negative
6. VonMises return mapping: verify equivalent stress = yield stress after return
7. MohrCoulomb return mapping: stress on yield surface after return
8. Uniaxial compression beyond yield: plastic strain > 0, total strain = elastic + plastic
9. Hydrostatic stress inside all yield surfaces (no yielding under pure pressure for DP and VM)
10. Document: PlasticityModel is NOT coupled to PetscDS callbacks. These are standalone tests only.

Register as `Physics.PlasticityReturnMapping` with label `physics_validation`.

### NEW TEST 5: Absorbing BC Energy Decay (tests/physics_validation/test_absorbing_bc_energy.cpp)

If test_absorbing_bc.cpp already has an energy test, verify it. If not, add one:

1. 3D box, granite properties (lambda=30e9, mu=25e9, rho=2650)
2. Initial Gaussian displacement pulse at center
3. Absorbing BCs on all 6 faces. Boundary config: sides=free, top=free, bottom=free, absorbing_bc_enabled=true, all faces true
4. Run for T = 2 * (L/2) / cp (enough for wave to cross domain and exit)
5. Measure total kinetic energy at t_final
6. Assert: E(t_final) < 0.1 * E(t_peak)
7. Compare with a control run using fixed BCs: E_fixed(t_final) should be >> E_absorbing(t_final)

Register as `Physics.AbsorbingBCEnergy` with label `physics_validation`.

### NEW TEST 6: Moment Tensor Source Verification (tests/physics_validation/test_moment_tensor_source.cpp)

After fixing BUG 5, verify the source produces correct physics:

1. Isotropic explosion source produces radially symmetric displacement pattern (|u_x| at (+r,0,0) approximately equals |u_y| at (0,+r,0))
2. Displacement amplitude decays as 1/r^2 in near field or 1/r in far field (geometric spreading)
3. P-wave arrival time at distance r matches r/vp within 20%
4. Source produces nonzero P-wave (radial displacement) and near-zero S-wave for isotropic source

This test runs the full Simulator pipeline with absorbing BCs, an explosion source, and seismometers at multiple distances and azimuths.

Register as `Physics.MomentTensorSource` with label `physics_validation`.

### REPLACE TEST 7: SCEC TPV5 (tests/physics_validation/test_scec_tpv5.cpp)

Delete the placeholder. Replace with a real test:

1. Set up a 3D box with a vertical fault at the center
2. Slip-weakening friction: mu_s = 0.677, mu_d = 0.525, D_c = 0.40 m
3. Background stress: tau_0 = 70 MPa shear, sigma_n = 120 MPa normal
4. Nucleation zone: reduced strength in a 3km x 3km patch
5. Run dynamic rupture through the Simulator (enable_faults=true, cohesive kernel)
6. Assert: rupture propagates (nonzero slip outside nucleation zone)
7. Assert: peak slip velocity > 0 at a station 5km along strike

If the cohesive fault dynamic rupture pipeline cannot run end-to-end, document exactly where it fails (mesh splitting? callback registration? SNES divergence?) and write the test to assert as far as the pipeline goes. Do NOT use GTEST_SKIP. Either the test passes or it fails with a documented reason.

Register as `Physics.SCEC.TPV5` with label `physics_validation`.

### FIX TEST 8: Boundary Conditions (tests/unit/physics/test_boundary_conditions.cpp)

Delete the current empty test. Replace with:

1. Test that bc_zero function returns zero for all components
2. Test that bc_compression function returns correct compression value
3. Test that bc_drained function returns zero pressure
4. Test that labelBoundaries correctly labels all 6 faces of a box mesh (verify label counts via DMGetLabelSize)
5. Test that setupBoundaryConditions registers the correct number of BCs for each configuration (elastic, poroelastic, elastodynamic)

Register as `Unit.BoundaryConditions` with label `unit`.

### EXISTING TESTS: Fix All Failures

After fixing bugs and writing new tests, run the full suite. Every test that should pass must pass. Specifically verify:

- `Physics.TerzaghiConsolidation` -- Poroelastic coupling end-to-end
- `Physics.GarvinsProblem` -- Buried explosion Green's function
- `Physics.LambsProblem` -- Point force on halfspace
- `Physics.AbsorbingBC` -- Callback correctness
- `Physics.GravityLithostatic` -- Gravity body force
- `Physics.LithostaticStress` -- Lithostatic initial stress
- `Physics.ElastostaticsPatch` -- Patch test
- `Integration.ExplosionSeismogram` -- End-to-end explosion seismogram
- `Integration.PrescribedSlip` -- Prescribed slip on cohesive fault
- `Integration.LayeredElastostatics` -- Layered material properties
- `Integration.FullSimulation` -- Basic simulation lifecycle
- `Integration.Restart` -- Checkpoint/restart

## CI Requirements

### Update .github/workflows/ci.yml

The CI must run ALL test categories and the PR must not merge unless all pass:

1. Verify the CI workflow has `required: true` or equivalent for all test jobs
2. Add the new tests to the appropriate labels in CMakeLists.txt
3. Verify timeout values are sufficient:
   - Unit tests: 120s
   - Functional tests: 180s
   - Physics validation: 600s (increase from 300s -- the energy decay and moment tensor tests need more time)
   - Integration tests: 600s
   - Performance tests: 900s

### Test Results Summary

After all tests pass, create `docs/TEST_RESULTS.md`:

```markdown
# FSRM Test Results

Generated: [date]
Branch: local_fix
PETSc: 3.22.2
Compiler: g++ [version]

## Summary
- Total tests: [N]
- Passed: [N]
- Failed: 0
- Skipped: 0

## Unit Tests
| Test | Status | Assertions |
|------|--------|------------|
| Unit.ConfigReader | PASS | 38 |
...

## Physics Validation Tests
...

## Integration Tests
...
```

Every test must be PASS. Zero FAIL. Zero SKIP (delete tests that cannot pass rather than skipping them).

## README Cleanup

1. Delete the entire "ResFrac-Equivalent Hydraulic Fracturing" section
2. Delete the volcano modeling section
3. Delete the tsunami modeling section
4. Delete the atmospheric infrasound section (keep the explosion infrasound mention)
5. Delete the ocean physics section
6. Delete specific GPU performance numbers (no benchmark data)
7. Delete the "Machine Learning Solvers" section (FNO is a stub)
8. Delete the "Earthquake Simulation (SeisSol-Compatible)" section (DG/ADER are stubs)
9. Delete the "Hypervelocity Impacts" subsection or mark it clearly as standalone analysis (not FEM-coupled)

Create a new top-level section "Verified Capabilities (All Tests Pass)" that lists ONLY features with corresponding green tests:
- Elastostatics and elastodynamics (TSALPHA2)
- Poroelastic coupling (Biot, Terzaghi verified)
- Mueller-Murphy seismic source (mb-yield, corner frequency, scalar moment)
- Explosion seismogram pipeline (source injection, wave propagation, absorbing BCs, seismometer sampling, SAC output)
- Atmospheric explosion physics (Sedov-Taylor, Brode, EMP)
- Near-field explosion phenomenology (cavity, damage, spall -- standalone)
- Cohesive fault mechanics (mesh splitting, prescribed slip, friction laws)
- Gravity and lithostatic prestress
- Plasticity models (standalone, not FEM-coupled)
- Lamb's and Garvin's problem verification

Move everything else to "Planned / In Development" with honest status.

## PR to Main

After everything above is done:

1. Run `ctest --output-on-failure` in Docker. Zero failures.
2. Update CLAUDE.md "What Works" and "What Does NOT Work" to match reality.
3. Commit with message: "Verified explosion monitoring pipeline: 6 bugs fixed, N new tests, all passing"
4. Push local_fix.
5. Create PR:

```bash
gh pr create \
  --base main \
  --head local_fix \
  --title "Verified multi-physics explosion monitoring pipeline" \
  --body "## Changes

### Bugs Fixed
1. Explosion seismogram test timing (P-wave never reached stations)
2. Absorbing BCs not activated in explosion test (roller masked them)
3. Mueller-Murphy corner frequency divided by burial depth instead of cavity radius
4. FEM source used wrong scalar moment (cavity M0 vs Mueller-Murphy M0, 16x discrepancy)
5. Explosion source injection: replaced body force with proper moment tensor equivalent nodal forces
6. Constants array conflict when gravity and absorbing BCs both enabled

### New Tests
- Physics.MuellerMurphy (12 assertions)
- Physics.AtmosphericExplosion (11 assertions)
- Physics.NearFieldExplosion (12 assertions)
- Physics.PlasticityReturnMapping (10 assertions)
- Physics.AbsorbingBCEnergy (energy decay verification)
- Physics.MomentTensorSource (radiation pattern, amplitude decay, travel time)
- Physics.SCEC.TPV5 (real dynamic rupture, not placeholder)
- Unit.BoundaryConditions (label verification, BC registration)

### Test Results
All [N] tests pass. Zero failures. Zero skips.
See docs/TEST_RESULTS.md for full matrix.

### README
Cleaned up to reflect only verified capabilities.
Removed unverified claims (volcano, tsunami, ocean, GPU perf, FNO, SeisSol DG).

### Known Limitations (documented in CLAUDE.md)
- PlasticityModel has real return-mapping algorithms but is not coupled to the PetscDS FEM solver
- NearFieldExplosion is a standalone 1D Lagrangian solver, not FEM-coupled
- Single homogeneous material per domain (no per-cell heterogeneity beyond depth layering)
- Gmsh mesh import exists but is not end-to-end verified
- DG, ADER, GPU acceleration, FNO are stubs"
```

6. Wait for CI. All jobs must be green. Fix anything that fails. Do not merge until every CI check passes.

## Execution Order

1. Phase 0: Build baseline, document current pass/fail
2. Fix BUG 3 (corner frequency) -- one-line fix, no risk
3. Fix BUG 4 (moment consistency) -- add members to ExplosionCoupling, use Mueller-Murphy momentRate
4. Fix BUG 5 (moment tensor source) -- proper equivalent nodal forces
5. Fix BUG 6 (constants conflict) -- add assertion or explicit wave speed storage
6. Fix BUG 1 + BUG 2 (explosion test config) -- timing and absorbing BCs
7. Write NEW TEST 1 (Mueller-Murphy)
8. Write NEW TEST 2 (atmospheric explosion)
9. Write NEW TEST 3 (near-field explosion)
10. Write NEW TEST 4 (plasticity return mapping)
11. Write NEW TEST 5 (absorbing BC energy)
12. Write NEW TEST 6 (moment tensor source verification)
13. Replace TEST 7 (SCEC TPV5)
14. Fix TEST 8 (boundary conditions)
15. Fix all remaining test failures
16. Update CI timeouts
17. README cleanup
18. CLAUDE.md update
19. Create docs/TEST_RESULTS.md
20. PR to main
21. CI green. Merge.

Do not skip steps. Do not hedge. Do not add "skip conditions." Every step produces a commit. Every commit passes the existing test suite.