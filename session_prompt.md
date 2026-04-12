# FSRM Development Prompt: CTBTO HPC Workshop Abstract Readiness

## READ FIRST

Read `CLAUDE.md` in the repository root completely before starting. All rules there are binding. In particular:

- NEVER change DS/BC ordering in setupFields()
- NEVER modify existing callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp
- NEVER modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS
- Build and test everything in Docker using Dockerfile.ci
- All existing tests must continue to pass after every change
- Check PETSc 3.22.2 API signatures before calling any PETSc function: `grep -rn "FunctionName" /opt/petsc-3.22.2/include/`

## Objective

Make the FSRM codebase on the `local_fix` branch truthfully represent this abstract:

> FSRM (Full Service Reservoir Model) demonstrates the power of AI-assisted development to create sophisticated scientific software for nuclear explosion monitoring. Using "vibecoding" -- iterative collaboration with large language models -- a comprehensive multi-physics simulator was developed that accurately reproduces historical nuclear test signatures for both underground and atmospheric detonations. The code implements Mueller-Murphy seismic source models for underground tests, capturing cavity formation, spall, chimney collapse, and seismic wave generation essential for treaty verification. Atmospheric detonation models include Brode fireball evolution, Sedov-Taylor blast wave propagation, and electromagnetic pulse (EMP) effects. Built on PETSc for parallel computing, the simulator handles complex coupling between elastodynamics, plasticity, fracture mechanics, and radiation transport. This "vibecoded" approach enabled rapid development of capabilities traditionally requiring years of specialized coding, producing a tool that validates against SCEC earthquake benchmarks and generates synthetic seismograms comparable to real IMS station data. The resulting open-source code provides nuclear monitoring agencies with a flexible platform for testing detection algorithms, analyzing historical events, and training analysts on explosion phenomenology. The project illustrates how AI-assisted development can democratize access to sophisticated geophysical simulation tools previously limited to national laboratories.

## Current State Assessment

### What WORKS (verified, tests exist and pass)

1. **Elastostatics PetscDS callbacks** (f0, f1, g3) with verified convergence
2. **Cohesive fault mesh splitting** (PyLith workflow, 32 cohesive cells on 4x4x4 simplex)
3. **Friction laws** (slip-weakening, rate-state aging)
4. **CoulombStressTransfer** (Hooke stress, fault projection, delta-CFS)
5. **Boundary condition labeling** on structured grids
6. **Mueller-Murphy source** (corner frequency, mb relation, moment computation) -- implementation exists, not end-to-end verified
7. **Unit tests** for config reader, physics kernels, domain models, etc.

### What EXISTS but is UNVERIFIED end-to-end

1. **Poroelasticity callbacks** (PetscFEPoroelasticity) -- unit tested, NOT verified through Simulator
2. **Absorbing boundary conditions** (AbsorbingBC.cpp) -- added on local_fix, test exists
3. **Gravity/lithostatic prestress** (PetscFEElasticityGravity.cpp) -- added on local_fix, test exists
4. **Terzaghi consolidation** -- test_terzaghi.cpp added on local_fix, includes Simulator pipeline test
5. **Explosion seismogram pipeline** -- test_explosion_seismogram.cpp added on local_fix
6. **Prescribed slip** -- test_prescribed_slip.cpp added on local_fix
7. **Lamb's problem** -- test_lambs_problem.cpp added on local_fix
8. **Garvin's problem** (buried explosion Green's function) -- test_garvins_problem.cpp added on local_fix
9. **Atmospheric explosion** (Brode fireball, Sedov-Taylor blast, EMP) -- NuclearAirburstEffects.cpp (707 lines), no verification test
10. **Near-field explosion** (cavity, spall, chimney) -- NearFieldExplosion.cpp (822 lines), no verification test
11. **Radiation transport** -- RadiationPhysics.cpp (1649 lines), no verification test
12. **Plasticity** -- PlasticityModel.cpp (1008 lines), CLAUDE.md says "stubs, do not use"
13. **Seismometer network** (DMInterpolation, SAC output) -- coded, wired to MonitorFunction
14. **SCEC TPV5 benchmark** -- test_scec_tpv5.cpp exists, unclear if it passes with real physics

### Known Gaps Between Abstract and Code

| Abstract Claim | Code State | Gap |
|---|---|---|
| "Accurately reproduces historical nuclear test signatures" | Mueller-Murphy coded, explosion seismogram test exists but unverified | No quantitative validation against real mb/Ms data |
| "Cavity formation, spall, chimney collapse" | NearFieldExplosion.cpp exists | No verification test |
| "Brode fireball evolution" | NuclearAirburstEffects.cpp exists | No verification test |
| "Sedov-Taylor blast wave propagation" | NuclearAirburstEffects.cpp exists | No verification test comparing R(t) to analytical |
| "EMP effects" | NuclearAirburstEffects.cpp has E1/E2/E3 | No verification test |
| "Coupling between elastodynamics, plasticity, fracture mechanics, and radiation transport" | Elastodynamics works; plasticity/radiation are untested | Plasticity marked as stub in CLAUDE.md |
| "Validates against SCEC earthquake benchmarks" | test_scec_tpv5.cpp exists | Unknown if test passes with quantitative comparison |
| "Generates synthetic seismograms comparable to real IMS station data" | AbsorbingBC + SeismometerNetwork + explosion source exist | No comparison against real IMS waveforms |

## Development Plan -- Execute in Order

Each phase produces a testable artifact. Do not start a phase until the previous phase's test passes. After each phase, run the full test suite (`ctest --output-on-failure`) to confirm no regressions.

### Phase 0: Establish Baseline -- Build and Run All Existing Tests

**Task:** Build the project in Docker and run ALL existing tests. Document which tests pass and which fail. This is the ground truth.

```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc) 2>&1 | tee build.log'
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure 2>&1 | tee test_results.txt
```

**Deliverable:** A file `test_baseline.txt` at the repo root listing every test name and PASS/FAIL. Fix any compilation errors first. If tests fail, categorize failures as:
- (A) Real bug -- fix it
- (B) Test expects physics that do not work yet -- mark as EXPECTED_FAIL and move on
- (C) Infrastructure issue (missing file, timeout) -- fix it

**Constraint:** Do NOT skip this phase. Everything downstream depends on knowing the actual baseline.

---

### Phase 1: Verify Absorbing Boundary Conditions

**Why:** Without absorbing BCs, seismic waves reflect off domain boundaries. CLAUDE.md explicitly says "seismograms are garbage after one transit time." This blocks every seismogram-related abstract claim.

**Task:** Ensure `tests/physics_validation/test_absorbing_bc.cpp` passes with quantitative criteria:

1. The Clayton-Engquist first-order absorbing BC callback `f0_absorbing` must produce traction `t_i = rho * c * v_i` for normal-incidence waves (P-wave: c = cp, S-wave: c = cs).
2. The Jacobian callback `g0_absorbing` must produce the correct diagonal: `u_tShift * rho * c`.
3. Pipeline test: A Simulator run with absorbing BCs enabled must complete without error.
4. **Energy test (new if not present):** A Gaussian pulse in a 3D box with absorbing BCs on all faces. Measure total kinetic energy at t=0 and after the pulse has crossed the domain. Energy at late time should be < 5% of initial energy (wave exits the domain instead of reflecting). If this test does not exist, create it.

**Verification:** `Physics.AbsorbingBC` passes in ctest.

**Files to check/modify:**
- `src/numerics/AbsorbingBC.cpp`
- `include/numerics/AbsorbingBC.hpp`
- `tests/physics_validation/test_absorbing_bc.cpp`
- `src/core/Simulator.cpp` (setupBoundaryConditions: verify absorbing BC registration path)

---

### Phase 2: Verify Terzaghi Consolidation (Poroelastic Coupling)

**Why:** The abstract claims "coupling between elastodynamics" which requires working Biot poroelasticity. Terzaghi is the canonical verification.

**Task:** Ensure `tests/physics_validation/test_terzaghi.cpp` passes:

1. **Callback-level test:** f0/f1 for poroelasticity produce correct values for known inputs matching analytical Terzaghi solution (pore pressure at mid-height at t = 0.1*H^2/cv).
2. **Pipeline test (TerzaghiPipelineTest):** The full Simulator pipeline (config -> setup -> solve) completes without error for a 1D poroelastic column.
3. **Quantitative validation:** If the pipeline test currently only checks "completes without error," extend it to extract the pressure solution at mid-height and compare against the analytical series solution. Tolerance: 10% relative error (coarse mesh).

**Verification:** `Physics.TerzaghiConsolidation` passes in ctest.

**Files to check/modify:**
- `tests/physics_validation/test_terzaghi.cpp`
- `src/numerics/PetscFEPoroelasticity.cpp`
- `src/core/Simulator.cpp` (poroelastic solve path)
- `config/examples/terzaghi_consolidation.config`

---

### Phase 3: Verify Lamb's and Garvin's Problems (Elastodynamic Wave Propagation)

**Why:** Lamb's problem (point force on halfspace) and Garvin's problem (buried explosion in halfspace) are the standard verification problems for elastodynamic codes. Garvin's problem directly validates the explosion seismogram pipeline.

**Task:**

1. Ensure `tests/physics_validation/test_lambs_problem.cpp` passes. At minimum: verify nonzero displacement at surface receivers, verify amplitude decay with distance, verify correct arrival time (distance/cp for P-wave).

2. Ensure `tests/physics_validation/test_garvins_problem.cpp` passes. This is the critical one: a buried isotropic expansion source (explosion) producing seismic waves in a halfspace. Quantitative checks:
   - P-wave first arrival at surface station at correct time
   - Radial pattern consistent with isotropic source
   - Amplitude at closer station > amplitude at farther station

**Verification:** `Physics.LambsProblem` and `Physics.GarvinsProblem` pass in ctest.

**Files to check/modify:**
- `tests/physics_validation/test_lambs_problem.cpp`
- `tests/physics_validation/test_garvins_problem.cpp`

---

### Phase 4: Verify End-to-End Explosion Seismogram Pipeline

**Why:** This is the core deliverable of the abstract: "generates synthetic seismograms."

**Task:** Ensure `tests/integration/test_explosion_seismogram.cpp` passes:

1. Explosion moment tensor is injected at the source cell
2. TSALPHA2 second-order time integration advances the solution
3. Absorbing BCs prevent spurious reflections
4. SeismometerNetwork samples displacement at station locations via DMInterpolation
5. SAC files are written with nonzero data
6. Closer station has larger amplitude than farther station

This test already exists. Make it pass. If it fails due to absorbing BCs not being wired in, wire them in. If it fails due to time integration issues (TSALPHA2 not converging), debug the IFunction/IJacobian callbacks for the elastodynamic case.

**Verification:** `Integration.ExplosionSeismogram` passes in ctest.

**Files to check/modify:**
- `tests/integration/test_explosion_seismogram.cpp`
- `src/core/Simulator.cpp` (FormFunction, MonitorFunction, explosion source injection)
- `src/domain/seismic/SeismometerNetwork.cpp`

---

### Phase 5: Mueller-Murphy Source Validation

**Why:** The abstract specifically claims "Mueller-Murphy seismic source models" and "accurately reproduces historical nuclear test signatures."

**Task:** Create `tests/physics_validation/test_mueller_murphy.cpp`:

1. **Corner frequency test:** For a 1 kt granite shot at 300m depth, verify corner frequency matches Mueller-Murphy (1971) Table 2 within 10%.
2. **mb-yield scaling test:** Compute body-wave magnitude for yields 1 kt, 10 kt, 100 kt, 1000 kt in granite. Verify mb follows the Mueller-Murphy scaling relation (approximately mb = 4.05 + 0.75*log10(W) for tamped shots in hard rock) within 0.3 magnitude units.
3. **Spectral shape test:** Verify the reduced displacement potential (RDP) spectral shape has the correct low-frequency plateau and high-frequency rolloff (omega^-2 or omega^-3 depending on source model variant).
4. **DPRK 2017 calibration:** For W=250 kt, depth=800m, granite medium (rho=2650, vp=5500, vs=3100): verify predicted mb is in the range 6.0-6.5 (observed mb was ~6.3).

**Verification:** `Physics.MuellerMurphy` passes in ctest.

**Files to create:**
- `tests/physics_validation/test_mueller_murphy.cpp`

**Files to modify:**
- `tests/CMakeLists.txt` (add new test)

---

### Phase 6: Atmospheric Explosion Physics Verification

**Why:** The abstract claims "Brode fireball evolution, Sedov-Taylor blast wave propagation, and electromagnetic pulse (EMP) effects."

**Task:** Create `tests/physics_validation/test_atmospheric_explosion.cpp`:

1. **Sedov-Taylor blast radius:** For a 20 kt airburst in standard atmosphere (rho_0 = 1.225 kg/m^3), verify blast radius R(t) = 1.15 * (E/rho_0)^0.2 * t^0.4 at t = 0.01s, 0.1s, 1.0s. Tolerance: 5%.

2. **Brode fireball max radius:** Verify fireball maximum radius scales as R_max ~ 90 * W^0.4 meters (W in kt). Test for 1 kt, 20 kt, 1000 kt. Tolerance: 20% (Brode model is approximate).

3. **EMP E1 peak field:** For a 100 kt burst, verify E1 peak field at 100 km ground range is order 25-50 kV/m. This is a parametric check, not exact, since the implementation uses a configurable e0_vm parameter.

4. **Overpressure at distance:** For a 20 kt surface burst, verify peak overpressure at 1 km is in the range 50-200 kPa (Glasstone-Dolan scaling).

**Verification:** `Physics.AtmosphericExplosion` passes in ctest.

**Files to create:**
- `tests/physics_validation/test_atmospheric_explosion.cpp`

**Files to modify:**
- `tests/CMakeLists.txt` (add new test)

---

### Phase 7: Near-Field Explosion Phenomenology

**Why:** The abstract claims "cavity formation, spall, chimney collapse."

**Task:** Create `tests/physics_validation/test_near_field_explosion.cpp`:

1. **Cavity radius:** For a 1 kt tamped shot in granite at 300m depth, verify computed cavity radius is in the range 10-20m (empirical: Rc ~ 8-15 * W^(1/3) meters for hard rock).

2. **Damage zone extent:** Verify the damage zone radius is 2-5x the cavity radius.

3. **Spall depth:** For a shallow explosion (depth < 5*Rc), verify spall computation produces a nonzero spall velocity at the free surface.

These are parametric/phenomenological checks, not exact solutions. The NearFieldExplosion class must produce physically reasonable numbers, not zeros or infinities.

**Verification:** `Physics.NearFieldExplosion` passes in ctest.

**Files to create:**
- `tests/physics_validation/test_near_field_explosion.cpp`

**Files to modify:**
- `tests/CMakeLists.txt`

---

### Phase 8: Plasticity -- Verify or Remove Claim

**Why:** CLAUDE.md Rule #7 says "Do NOT use DG, ADER, GPU, or plasticity features. They are stubs." But the abstract claims "plasticity."

**Decision point:**

**Option A -- Verify plasticity works:** Create `tests/physics_validation/test_drucker_prager.cpp`. Simple uniaxial compression test on a single element. Apply stress beyond yield surface. Verify plastic strain is nonzero and stress stays on/inside yield surface. If PlasticityModel.cpp actually works, update CLAUDE.md to remove "plasticity" from the stubs list.

**Option B -- Plasticity is truly broken:** If PlasticityModel.cpp does not produce correct results (stress exceeds yield, plastic strain is zero, solver diverges), then:
- Do NOT try to fix it. That is a multi-week effort.
- Instead, update CLAUDE.md to document exactly what is broken.
- The abstract claim about "plasticity" stands on the existence of the implementation and the Drucker-Prager/von Mises/Mohr-Coulomb code structure, even if end-to-end verification is incomplete. This is defensible for a workshop abstract (stating capability, not validated results) but be prepared to discuss limitations if asked.

**Deliverable:** Either a passing test or an honest note in CLAUDE.md documenting the gap.

---

### Phase 9: SCEC TPV5 Benchmark Validation

**Why:** The abstract claims "validates against SCEC earthquake benchmarks."

**Task:** Verify `tests/physics_validation/test_scec_tpv5.cpp` passes with quantitative comparison:

1. The test should set up a 2D strike-slip fault with linear slip-weakening friction (SCEC TPV5 parameters).
2. Rupture should propagate along the fault.
3. Slip velocity time history at on-fault stations should be compared against SCEC reference solutions.
4. If the existing test only checks "solver converges" or "slip is nonzero," extend it to compare against the Day et al. reference: peak slip velocity within 20%, rupture arrival time within 10%.

If the SCEC test is not passing because dynamic rupture through the CohesiveFaultKernel is broken, document the failure mode. At minimum, the test must demonstrate that:
- The fault mesh splits correctly
- Friction law is applied
- Slip initiates in the nucleation zone
- Rupture propagates (slip velocity > 0 at stations away from nucleation)

**Verification:** `Physics.SCEC.TPV5` passes in ctest.

---

### Phase 10: Synthetic vs. Real IMS Data Comparison

**Why:** The abstract claims seismograms "comparable to real IMS station data."

**Task:** Create `tests/integration/test_dprk_2017_comparison.cpp` or a Python validation script `scripts/validate_dprk_2017.py`:

1. Run FSRM with `config/examples/dprk_2017_quick.config` (or equivalent).
2. Extract synthetic seismogram at a station ~1000 km distance.
3. Compute synthetic body-wave magnitude (mb) from the peak displacement on the vertical component in the 0.5-5 Hz band.
4. Compare synthetic mb against the known observed mb of ~6.3 for the DPRK 2017 test.
5. Tolerance: synthetic mb should be in the range 5.5-7.0 (within 1 magnitude unit).

This does NOT require waveform-level match (that would require realistic 3D Earth structure). It requires that the Mueller-Murphy source, coupled through the elastodynamic solver with absorbing BCs, produces a seismogram whose amplitude is in the correct ballpark.

If the full Simulator pipeline cannot yet run this end-to-end, document exactly where it fails and create a standalone test that:
- Uses the Mueller-Murphy source to generate a moment rate function
- Convolves with a simple Green's function (far-field P-wave: 1/r * moment_rate(t - r/vp))
- Computes mb from the resulting synthetic
- Verifies against observed value

**Verification:** Synthetic mb for DPRK 2017 is in [5.5, 7.0].

---

### Phase 11: Clean Up README and Documentation

**Why:** The main branch README claims volcano modeling, tsunami modeling, ocean physics, FNO solvers, ResFrac-equivalent hydraulic fracturing, GPU acceleration, and dozens of other features that do not work. This is dishonest and undermines credibility.

**Task:**

1. Create a new README section "Verified Capabilities" that lists ONLY features with passing tests.
2. Move all unverified features to a "Planned / In Development" section with honest status labels.
3. Remove performance claims (GPU speedups, scaling numbers) that have no supporting benchmark data.
4. Remove the "ResFrac-Equivalent" section entirely unless hydraulic fracturing has a passing test.
5. Remove volcano, tsunami, ocean physics sections unless they have passing tests.
6. Update the config table to clearly mark which configs are tested and which are aspirational.
7. Keep `config/aspirational/` separate from `config/examples/` as CLAUDE.md already requires.

**Verification:** README accurately reflects the tested state of the code.

---

### Phase 12: Final Test Suite and PR

**Task:**

1. Run the complete test suite: `ctest --output-on-failure`
2. Every test must either PASS or be explicitly marked as EXPECTED_FAIL with a documented reason.
3. Zero compilation warnings (or document unavoidable ones).
4. Create a test summary document `docs/TEST_RESULTS.md` listing every test and its status.
5. Update CLAUDE.md "What Works" and "What Does NOT Work" sections to match reality.
6. Commit all changes with message: "CTBTO abstract readiness: verified explosion seismogram pipeline, Mueller-Murphy, absorbing BCs, Terzaghi, SCEC TPV5, atmospheric explosion physics"

Then create a PR from `local_fix` to `main`:

```bash
git checkout local_fix
git push origin local_fix
# Create PR via GitHub CLI or web interface:
gh pr create \
  --base main \
  --head local_fix \
  --title "CTBTO HPC Workshop: Verified multi-physics explosion monitoring pipeline" \
  --body "## Summary

This PR brings verified nuclear explosion monitoring capabilities to main.

### Verified (all tests pass):
- Mueller-Murphy seismic source model with mb-yield validation
- End-to-end explosion seismogram pipeline (source injection -> elastodynamic solve -> absorbing BCs -> seismometer sampling -> SAC output)
- Absorbing boundary conditions (Clayton-Engquist first-order)
- Terzaghi consolidation (Biot poroelastic coupling)
- Lamb's problem (point force on halfspace)
- Garvin's problem (buried explosion Green's function)
- Atmospheric explosion physics (Sedov-Taylor, Brode, EMP)
- Near-field explosion phenomenology (cavity, damage zone, spall)
- SCEC TPV5 dynamic rupture benchmark
- Prescribed slip on cohesive faults
- Lithostatic prestress with gravity

### Test results:
See docs/TEST_RESULTS.md for complete test matrix.

### Known limitations:
- Plasticity models exist but are not end-to-end verified
- Radiation transport physics exist but are not verified
- Synthetic seismograms are qualitatively correct (right order of magnitude) but not waveform-matched against real data
- Single material (homogeneous) only; no per-cell heterogeneity
- No Gmsh mesh import verified end-to-end

### Abstract supported by this PR:
CTBTO 3rd HPC Workshop for Nuclear Explosion Monitoring, Vienna, May 2026"
```

Wait for CI to pass on the PR. If CI fails, fix the failures and push again. Do not merge until all CI checks are green.

---

## General Rules for All Phases

1. Every new .cpp file needs a corresponding entry in `tests/CMakeLists.txt` (for tests) and a CTest registration.
2. Every new test must have quantitative pass/fail criteria with numerical tolerances. No "looks right" tests.
3. New PetscDS callbacks go in NEW files. Do not modify existing callback files.
4. Use the Docker build for all compilation and testing.
5. After each phase, run `ctest --output-on-failure` and confirm all tests (old + new) pass before proceeding.
6. Commit after each passing phase with a descriptive message referencing the phase number.
7. If a PETSc function does not exist in 3.22.2 headers, find the equivalent that does. Do not guess.
8. No Python in the C++ codebase. Python scripts are separate post-processing tools in `scripts/`.
9. When in doubt about PETSc DMPlex API, look at `ex17.c` (elasticity with aux fields), `ex56.c` (elasticity with BCs), and `ex62.c` (cohesive cells) in the PETSc source.
10. The explosion seismogram pipeline (Phase 4) is the most important deliverable. If time is limited, prioritize Phases 0-1-4-5 over everything else.
11. Do NOT lie about test results. If a test fails, document WHY it fails and what would be needed to fix it. An honest "this does not yet work" is infinitely better than a fake passing test.
12. Verification tests must have quantitative pass/fail criteria with numerical tolerances, not just "runs without crashing."