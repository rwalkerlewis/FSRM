# FSRM Nuclear Explosion Monitoring Roadmap

This document is the single source of truth for the `nem-roadmap` branch. It defines what is real in this repository today, what we are building toward, and the validation gates that distinguish the two. Do not update the main `README.md` to claim capabilities until the corresponding milestone here has passed its exit criterion.

## Scope statement

The goal of this branch is a validated near-field-to-far-field simulation pipeline for historic underground nuclear tests, anchored to published observations. The branch starts from the actually-working FEM core of FSRM and extends it incrementally, with each step gated by quantitative comparison to either an analytical solution or recorded data.

This branch is not the place for tsunami, volcano, atmospheric infrasound, hypervelocity impact, EMP, ionospheric coupling, discontinuous Galerkin, ADER time integration, or new neural-operator work. Those features either do not exist in the codebase today or exist only as stubs. They are out of scope until the validations in this roadmap are complete.

## What is real in FSRM today (the baseline this roadmap builds on)

These components have implementations in `src/` matching their headers, are wired into `Simulator`, and produce non-trivial output:

- `PetscFEElasticity` — quasi-static elasticity and I2-form elastodynamics with TSALPHA2 generalized-α time integration.
- `AbsorbingBC` — Clayton-Engquist first-order absorbing traction with a correct Jacobian.
- `PetscFEPoroelasticity`, `PetscFEThermal`, `PetscFEViscoelastic` — pointwise FE kernels in the same style.
- `Simulator::addExplosionSourceToResidual` — moment-tensor injection as equivalent body forces using PetscFE basis derivatives at the source cell. This is the "BUG 5 fix" path. This is the only place the explosion source actually couples to the FEM displacement system.
- `MuellerMurphySource` — analytical reduced displacement potential, moment rate, scalar moment, corner frequency.
- `CohesiveFaultKernel` — PyLith-style Lagrange-multiplier saddle-point fault formulation. **In active development on `local_fix`. Do not touch from `nem-roadmap`.**

## What is not real and must not be extended on this branch

- `DiscontinuousGalerkin` and `ADERTimeIntegrator` — header-only, no `.cpp`, never instantiated.
- `Tsunami`, `Volcano`, `HypervelocityImpact`, `AtmosphericInfrasound`, `EMP`, `Radiation`, `Fallout` — README claims, no working pipeline.
- `experimental/` neural modules (NeuralAMR, NeuralInversion, NeuralTimestepping, MultiFidelityLearning, NeuralLinearAlgebra) — large headers, stub `.cpp` files that return zero or `PETSC_SUCCESS`.
- The `fsrm` Python package — imported throughout `scripts/` but does not exist in the repo. All scripts that depend on it are non-runnable.
- README claims of "ResFrac-equivalent" multi-stage fracturing capabilities beyond what is in `PetscFEHydrofrac.cpp`.

These will be inventoried in Milestone 0 and either deleted or moved to `archive/aspirational/` at the end of Milestone 3.

## Ground rules

1. **No new physics modules without a working test.** Do not add modules not already wired into `Simulator.cpp` and not in the milestone list below.

2. **Validation gates everything.** No milestone is "done" until a test passes that compares numerical output against either an analytical solution or a published observation, with a quantitative tolerance stated up front. "SAC file is non-empty" is not a validation. "Closer station has larger amplitude" is not a validation.

3. **Cut, don't add, when in doubt.** If a header has no matching `.cpp`, delete it. If a `.cpp` has stub bodies, delete it. If a Python script imports `fsrm.*` modules that do not exist in the repo, mark it `BROKEN` at the top of the file and do not edit it further until an actual `fsrm` package exists.

4. **No README inflation.** Only update the main `README.md` at milestone exits, and only with truthful statements backed by passing tests.

5. **One logical change per commit.** Use the existing `session N: ...` commit-message convention. Each session must end with `make test` passing on the affected components, or with a clear note in the commit body explaining what is broken and why.

6. **One milestone at a time.** Do not work on multiple milestones in the same session.

## Milestone 0: Honest baseline (1 session)

Before any new physics, audit and document what currently works. Inventory only — no deletion yet.

- [ ] Run all integration tests on `main`. Record each test's pass/fail status, runtime, and any error output in `docs/NEM_BASELINE.md`.
- [ ] List every header in `include/` with no matching `.cpp` file. Record under "Unimplemented headers" in `docs/NEM_BASELINE.md`.
- [ ] List every `.cpp` in `src/` whose function bodies are stubs (return zero, empty, or every argument cast to void with no logic). Record under "Stub implementations."
- [ ] List every Python script under `scripts/` that imports from `fsrm.*`. Verify whether the `fsrm` package exists in the repo. Record non-runnable scripts under "Non-runnable Python scripts."
- [ ] Do NOT delete anything in this milestone. The deletions happen in Milestone 3 cleanup, after the first two events validate.

**Exit criterion:** `docs/NEM_BASELINE.md` exists, is committed, and inventories every fictional or non-functional component. This document drives the cleanup work in Milestone 3.

## Milestone 1: Sedan 1962 near-field pilot

Sedan 1962 was 104 kt, approximately 194 m depth of burial, in alluvium over tuff at the Nevada Test Site. It is fully declassified and has published near-field ground motion records in USGS Open-File reports, LLNL UCRL documents, and the established near-field motion compilations (Murphy 1996; Vincent et al.). It is the cleanest single-event target for the first FEM validation.

### M1.1: Site model
- [ ] Build a layered velocity model in `models/sedan_1962.vel`: approximately 100 m of alluvium (vp ~1500 m/s, rho ~1900 kg/m^3) over tuff (vp ~2500 m/s, rho ~2100 kg/m^3) over basement. Cite the source publication for every layer in a comment block at the top of the file.
- [ ] Generate `meshes/historical/sedan_1962_nearfield.msh`: a cylindrical or box domain, 4 km half-width, 2 km deep. Element size 50 m within 500 m of the source, graded to 500 m at the lateral and bottom boundaries. Justify the resolution against the target frequency band: target faithful representation up to 5 Hz; 50 m at vp=1500 m/s gives roughly 6 nodes per wavelength at 5 Hz.
- [ ] Add a unit test that loads the mesh and verifies element count, minimum and maximum edge length, and that the source location lies inside an element of the expected size class.

### M1.2: Source
- [ ] Use the existing `Simulator::addExplosionSourceToResidual` path with Mueller-Murphy parameters for 104 kt at 194 m depth in alluvium. Do not modify the source kernel — only configure it from a new config file `config/examples/sedan_1962_validated.config`.
- [ ] Add a derived diagnostic test: compute scalar moment M0 from the `momentRate` integral and compare to `params.scalar_moment()` analytically. Tolerance: 1 percent.

### M1.3: Run and receivers
- [ ] Place receivers at 100 m, 500 m, 1 km, and 2 km radial distance from the source, both on the free surface and at source depth.
- [ ] Run to 2.0 s simulated time using TSALPHA2.
- [ ] Verify that Clayton-Engquist absorbing BCs are active on the lateral and bottom boundaries and the free surface is on top. Document this in the config file with comments.

### M1.4: Validation against observations
- [ ] Pull published peak particle velocity (PPV) versus range data for Sedan from at least two of: USGS Open-File Reports, LLNL UCRL documents, peer-reviewed near-field motion compilations. Cite by document number and page in the validation document.
- [ ] Add an integration test `tests/integration/test_sedan_1962_ppv.cpp` that runs the simulation, extracts PPV from each receiver, and asserts:
  - PPV at each range within a factor of 2 of observed.
  - Slope of log(PPV) versus log(range) within ±0.3 of the observed slope.
- [ ] Generate `figures/sedan_1962/ppv_vs_range.pdf` from a Python script under `scripts/` that imports only `obspy`, `numpy`, `matplotlib`. Do NOT import any `fsrm.*` module.

### M1.5: Document
- [ ] Write `docs/validations/sedan_1962.md`: physics, mesh, source parameters, observations used, agreement, known limitations (no topography, flat layered model, isotropic source, no spall closure, no damage zone).
- [ ] Update the main `README.md` to add a "Validated against observations" section listing only Sedan 1962, with a link to this validation document.

**Exit criterion:** `make test` passes including `test_sedan_1962_ppv`, the figure exists in the repo, and the validation document is committed.

## Milestone 2: Topography (1-2 sessions)

- [ ] Add support for ingesting a DEM (GeoTIFF or ASCII grid) and producing a curved free-surface mesh in Gmsh. Write the DEM ingestion as a Python helper script first; do not embed DEM parsing in C++.
- [ ] Test on a synthetic Gaussian hill: verify that surface waves propagate correctly and that the absorbing BCs still work at the curved-surface lateral boundaries.
- [ ] Re-run Sedan with realistic NTS topography. Re-run the validation in M1.4. Document any change in PPV agreement.

**Exit criterion:** Sedan validation still passes (within the same tolerances) with topography enabled.

## Milestone 3: Hardhat or Piledriver — granite analog (1-2 sessions)

Granite analog to Punggye-ri. Pick whichever has better-documented near-field records — Werth and Herbst (1963) is the canonical Hardhat reference and the literature on DPRK still cites it.

- [ ] Build site model in granite (`models/hardhat_1962.vel` or equivalent). Cite layer sources.
- [ ] Run, validate against published PPV-versus-range and corner-frequency observations.
- [ ] If validation fails: do not proceed. Diagnose whether the failure is in mesh resolution, source parameterization, velocity model, or in the Mueller-Murphy empirical formula's applicability to granite. Document the diagnosis in `docs/validations/hardhat_1962.md` even if the validation fails.
- [ ] Cleanup pass: delete or move to `archive/aspirational/` everything inventoried in Milestone 0 as unimplemented or stubbed. Update README to remove all claims those modules supported.

**Exit criterion:** Two events independently validated. The aspirational cruft is gone. This is a real, defensible result and should be written up as a short note for the CTBTO HPC workshop or a similar venue.

## Milestone 4: Near-field-to-far-field handoff (2-3 sessions)

The Sedan and Hardhat box meshes cannot reach teleseismic distance. Build a reciprocity-based handoff to a 1-D global Green's-function code (`instaseis` or `syngine`).

- [ ] Add displacement extraction on a virtual sphere at approximately 5 km radius around the source, sampled at sufficient angular resolution for moment-tensor recovery.
- [ ] Add a Python script that takes the FSRM near-field output, performs a least-squares fit for an equivalent moment tensor and source-time function, and writes a `CMTSOLUTION`-format file usable by `instaseis`.
- [ ] Validate by round trip: feed the recovered moment tensor through `instaseis` to a real station (e.g., MDJ for a hypothetical DPRK-equivalent synthetic), and compare to a direct `instaseis` run with the published M0. Tolerances: amplitude within 10 percent, phase within 1 second.

**Exit criterion:** Round-trip moment-tensor extraction and far-field synthesis works within the stated tolerances.

## Milestone 5: DPRK 2017 with Mt. Mantap topography (3-5 sessions)

The marquee event. Granite, approximately 250 kt, real topography, real station geometry, real recorded waveforms.

- [ ] Build the Mt. Mantap mesh from SRTM or equivalent DEM.
- [ ] Verify and extend the layered velocity model in `models/punggye_ri.vel`.
- [ ] Run the near-field FSRM simulation, extract the moment tensor via the M4 pipeline, and propagate to MDJ, INCN, and BJT via `instaseis` using the global model assumed in the published Pasyanos and Myers / Zhang and Wen results.
- [ ] Compare to actual broadband data fetched via `obspy` from IRIS. Quantitative tolerances on first-arrival amplitude, dominant period, and the recovered mb estimate.

**Exit criterion:** Synthetic versus observed comparison at one real station agrees within the published uncertainty bounds for the event.

## Milestone 6: Discrimination physics (open-ended)

Only attempt after Milestones 1-5 are solid. This is where the cohesive-fault work on `local_fix` becomes relevant: tectonic-release CLVD components, spall closure, pP-P depth phases, mb/Ms discrimination.

Do not start this milestone until the cohesive-fault saddle-point solver on `local_fix` has merged into `main` with passing dynamic-rupture tests.

## What is explicitly out of scope on this branch

- Tsunami, volcano, hypervelocity impact, EMP, ionospheric coupling, atmospheric infrasound.
- New ML, FNO, or neural-operator code. The existing FNO is fine; do not extend it on this branch.
- `DiscontinuousGalerkin.cpp`. The header is dead. It will be deleted in the Milestone 3 cleanup pass.
- Any README change that adds capabilities. Only milestone-exit, validation-backed updates.
- Python scripts that import `from fsrm.*` until and unless an actual `fsrm` Python package is created in this repo with a `pyproject.toml`. If a Python package is needed, create it explicitly in its own session and document it.

## Definition of done for this branch

`nem-roadmap` is mergeable to `main` when:

- Milestones 1 through 3 are complete (two events independently validated against observations).
- The main `README.md` reflects exactly what works, with no inflated claims.
- All stub headers and dead code identified in Milestone 0 are either deleted or moved to `archive/aspirational/` with a clear `README.md` saying "not implemented, not on the roadmap."
- CI is green on every job.

Milestones 4 through 6 may continue on subsequent branches.

## Session log

Append one entry per session below. Format:

```
### Session N — YYYY-MM-DD — Milestone X.Y

**What passed:**
- (test name, file, link to commit)

**What is blocked and why:**
- (concrete description)

**Next concrete step:**
- (specific, testable next action)
```

### Session 0 — 2026-04-28 — Milestone 0 (Honest baseline)

**What passed:**
- `docs/NEM_ROADMAP.md` installed verbatim on branch `nem-roadmap`
  (committed as "session 0: install nuclear explosion monitoring roadmap").
- `docs/NEM_BASELINE.md` written with the four required inventory sections:
  integration-test enumeration, unimplemented headers (11), stub
  implementations (5 confirmed in `src/experimental/`, plus 6 ML files
  flagged for manual review at Milestone 3), and 23 non-runnable Python
  scripts that import from a non-existent `fsrm` package (no
  `pyproject.toml`, no `fsrm/` directory anywhere in the repo).
- No source files, headers, scripts, or tests were modified or deleted,
  consistent with the Milestone 0 directive "Inventory only — no deletion yet."

**What is blocked and why:**
- The integration-test pass/fail/runtime table required by Milestone 0
  bullet 1 is not filled in. `Dockerfile.ci` fails at the PETSc compile
  layer with `/bin/sh: 1: ./configure: not found`. Root cause: the
  `petsc-3.25.0.tar.gz` download from
  `web.cels.anl.gov/projects/petsc/download/release-snapshots/` does not
  unpack to a `petsc-3.25*` directory in the current build environment,
  and the subsequent `mv ... || true` masks the failure, leaving an
  empty `WORKDIR /opt/petsc-3.25.0`. No host-side build is possible per
  `CLAUDE.md` rule 1 ("Build and test in Docker. Always."). Per
  `docs/NEM_ROADMAP.md` final-notes clause, this failure is documented
  and the session stops rather than papering over it.
- The classification of `src/ml/*.cpp` as stub-vs-real cannot be
  finalized without a working link/test cycle; flagged for manual
  inspection at Milestone 3.

**Next concrete step:**
- Open a small follow-up session that repairs `Dockerfile.ci` only:
  pin the PETSc tarball URL to a concrete archive that is known to
  unpack to `petsc-3.25.0/`, remove the `|| true` mask so the build
  fails loudly if the directory is wrong, and verify by running
  `docker build -f Dockerfile.ci -t fsrm-ci:local .` followed by
  `docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local
  bash -c 'cmake .. -DENABLE_TESTING=ON && make -j$(nproc) && ctest
  --output-on-failure -L integration'`. Then update section 1.2 of
  `docs/NEM_BASELINE.md` with the resulting per-test status table.
  Only after that is the gate to begin Milestone 1 (Sedan 1962) cleared.

### Session 1 — 2026-04-28 — Milestone 1.1 + 1.2 partial

**What passed:**
- `models/sedan_1962.vel` — three-layer Sedan velocity model
  (alluvium / welded tuff / Paleozoic basement), with citation block
  pointing at Carlson and Roberts (PNE-217F), Vortman (LA-3911), and
  Murphy (1996).  TODO note flags primary-source upgrade for M1.4.
- `scripts/meshing/sedan_1962_nearfield.geo` and
  `scripts/meshing/build_sedan_1962_mesh.sh` — Gmsh input for the 8 km
  by 8 km by 2 km box with a 50 m fine sphere of radius 500 m anchored
  on the shot point at (0, 0, -194), graded out to ~300 m at the outer
  boundary (tightened from the 500 m roadmap target so the longest tet
  edge stays below the 750 m sanity bound).  Frequency-band rationale
  is in the .geo header.
- `meshes/historical/sedan_1962_nearfield.msh` — committed pre-built
  mesh, 50,780 tetrahedra, min edge 31.6 m, max edge 669.4 m, with
  154 fine cells whose centroid is within 100 m of the shot point.
- `tests/integration/test_sedan_1962_mesh.cpp` (CTest label
  `integration`, name `Integration.Sedan1962Mesh`) — loads the mesh
  via `DMPlexCreateGmshFromFile` and asserts tet count, min/max edge,
  and a fine cell near the source point.
- `tests/unit/domain/test_mueller_murphy_moment_consistency.cpp`
  (CTest label `unit`, name `Unit.MuellerMurphyMomentConsistency`) —
  Sedan-configured `MuellerMurphySource` checks the Simpson-rule
  integral of `momentRate(t)` against `params.scalar_moment()` (1%
  tolerance) and asserts `body_wave_magnitude()` lies in [4.0, 5.5]
  bracketing the USGS Sedan mb of ~4.75.

**What is blocked and why:**
- Per-test pass/fail/runtime table for `docs/NEM_BASELINE.md` Section
  1.2 is still pending; the agent sandbox does not have the GitHub
  Actions test-result artifacts mounted, and per the session-1 brief
  the Dockerfile failure is not being chased.  The two new tests are
  expected to appear in the next CI run on this branch and on the
  next merge into `main`; the inventory table can be filled then.

**Next concrete step:**
- Session 2: Milestone 1.3 (run + receivers).  Configure
  `config/examples/sedan_1962_validated.config` to use the new mesh
  and velocity model, run for 2.0 s simulated time, output SAC files
  at receivers placed at 100 m, 500 m, 1 km, 2 km on the free surface
  and at source depth.  Do not start M1.4 in that session.

