# Session 2 Report

Date: 2026-04-18
Branch: `local_fix`

## 1. PETSc commit selected

- Repo: gitlab.com/petsc/petsc
- Branch: `main`
- Commit SHA: `267b8824abdb82fedbefdec8ea8d0cb3c505ddfd`
- Commit date: 2026-04-17 11:09:30 -0500
- Commit subject: "Merge remote-tracking branch 'origin/release'"
- PETSc git describe: `v3.25.0-117-g267b8824abd` (PETSc main is now in the
  3.25.0 line, 117 commits ahead of the v3.25.0 tag).

How chosen: PyLith does not build against PETSc main; its `DEPENDENCIES`
file and `docker/pylith-testenv` Dockerfile track the `knepley/pylith`
fork branch, which has no specific SHA pinned. The user instructed to
use PETSc `main` HEAD instead of the PyLith fork. The HEAD of
`refs/heads/main` on gitlab.com/petsc/petsc was captured via
`git ls-remote` and pinned in every Dockerfile via the `PETSC_SHA` ARG.

## 2. Dockerfile changes

All four Dockerfiles now build the same PETSc commit from source. Each
file's distinctive purpose was preserved:

| File | Status | What changed |
|------|--------|--------------|
| `Dockerfile` | rewritten | Was PETSc 3.20.0 tarball. Now clones gitlab `main` at the pinned SHA, builds with `-O3 -march=native`, otherwise unchanged (multi-stage production build that bakes FSRM source in and ships with `obspy`/`pygmt`/`gmt` runtime). |
| `Dockerfile.ci` | rewritten | Was PETSc 3.25.0 tarball. Now clones gitlab `main` at the pinned SHA. CI build environment, no FSRM source baked in. |
| `Dockerfile.cuda` | rewritten | Was PETSc 3.22.2 tarball. Now clones gitlab `main` at the pinned SHA. Adds `--with-cuda=1 --with-cuda-dir=/usr/local/cuda`. |
| `Dockerfile.dev` | rewritten | Was apt `petsc-dev` (Ubuntu repository version). Now builds PETSc from source at the same SHA so the development environment matches CI/production. Also adds Google Test build, since `petsc-dev` was previously providing PETSc only and `gtest` was relied on from apt. |

Stale comment headers fixed in all four files (each file now accurately
describes which PETSc it builds and what it is for).

No Dockerfile was deleted. A new `DOCKER.md` at the repository root
documents the four-file structure (production / CI / CUDA / dev) and
explains how to bump the PETSc SHA. The user's prompt explicitly
permitted this layout: "A single canonical production Dockerfile plus
Dockerfile.ci plus Dockerfile.cuda plus Dockerfile.dev is acceptable".

Each PETSc clone records its checked-out SHA at
`/opt/petsc-main/COMMIT_SHA` for runtime introspection.

## 3. FSRM build outcome

Build was clean: exit 0, only the two pre-existing
`-Wunused-variable` warnings in
`tests/integration/test_viscoelastic_wave.cpp:143-144` (`auxDM`,
`auxVec`). No PETSc API incompatibilities surfaced when moving from
PETSc 3.25.0 release to PETSc main `v3.25.0-117-g267b8824abd`.

Logs:
- `/tmp/build.log` — Docker image build (PETSc compile from source, ~270 s)
- `/tmp/fsrm_build.log` — FSRM cmake + make
- `/tmp/full_test.log` — ctest output (2001 lines)

## 4. Test count

| Run | Pass | Fail | Total |
|-----|-----:|-----:|-----:|
| Baseline (per `docs/FAULT_TEST_REGRESSION_AUDIT.md`, PETSc 3.25.0 release) | 110 | 6 | 116 |
| This session, `ctest -j$(nproc)` | 108 | 8 | 116 |
| This session, sequential re-run of suspect failures | 109 | 7 | 116 |

The two extra failures under `ctest -j` (Physics.AbsorbingBC,
Integration.TractionBC) both pass in a sequential rerun. Their failure
mode is HDF5 file-handle contention when two tests write to the same
`output/solution.h5` concurrently — the documented parallel-test
limitation in CLAUDE.md ("some tests may fail when run in parallel
(`ctest -j`) due to HDF5 output file conflicts. All pass when run
individually"). They are not regressions from this session's changes.
Confirmed by re-running the three suspect tests serially:

```
1/3 Test #43: Physics.AbsorbingBC ..............   Passed    7.68 sec
2/3 Test #59: Physics.CohesiveBdResidual .......***Failed    0.44 sec
3/3 Test #97: Integration.TractionBC ...........   Passed    0.70 sec
```

(Captured at `/tmp/sequential_recheck.log`.)

### Newly-passing tests (vs. baseline)
- `Integration.TimeDependentSlip` — was in the baseline failure list
  (the audit lists it under "Known failures"). Now passes (1.43 sec).
  This is the clearest signal that the cohesive-complete label change
  is doing what the diagnosis predicted.

### Newly-failing tests (vs. baseline)
- `Physics.CohesiveBdResidual` — was in the baseline pass list
  (CLAUDE.md "Standalone (correct, tested, not FEM-coupled)" table,
  "PetscDS BdResidual on cohesive"). The new failure persists in the
  sequential rerun, so it is a real regression caused by switching
  the `fault_constraint` BC from `fault_label` (face-only) to
  `interfaces_label_` (cohesive cell + closure faces, edges, vertices).
  Assertion: `tests/physics_validation/test_cohesive_bd_residual.cpp:233`,
  `EXPECT_TRUE(std::isfinite(f_norm))` — `f_norm` evaluates to NaN/Inf.
  Likely cause: the BdResidual callback is being invoked on label
  points that are not facets (e.g. cohesive cells themselves, or
  vertices), feeding ungeometric quadrature data to the f0 callback
  and producing NaN.

## 5. Status of the six fault tests

All six remain failures, but two show qualitative changes from the
label fix:

| Test | Pass/Fail | Failed assertion |
|------|-----------|------------------|
| `Physics.LockedFaultTransparency` | FAIL | `tests/physics_validation/test_locked_fault_transparency.cpp:341`, `EXPECT_GT(fault_norm, 0.0)` (`fault_norm = 0`). Same mode as baseline. |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | FAIL | `tests/integration/test_dynamic_rupture_solve.cpp:389`, `EXPECT_GT(sol_norm, 0.0)` (`sol_norm = 0`). SNES now iterates 200 times decreasing slowly (4 reattempts each step, all `DIVERGED_MAX_IT`) instead of the baseline `DIVERGED_FUNCTION_NANORINF iterations 0` failure mode. Reaching `sol_norm = 0` is from TS step truncation, not SNES success. Different failure mode but still failing. |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | `tests/integration/test_dynamic_rupture_solve.cpp:433`, `EXPECT_GT(max_fault_slip, 5.0e-4)`. **`max_fault_slip = 2.15241e-23`** (see section 6). |
| `Integration.TimeDependentSlip` | **PASS** | n/a — this was in the baseline failure list and is now passing. |
| `Integration.SlippingFaultSolve` | FAIL | `tests/integration/test_slipping_fault_solve.cpp:333`, `EXPECT_EQ(ierr, 0)`. SNES diverges with `DIVERGED_FUNCTION_NANORINF` repeatedly. |
| `Integration.SlipWeakeningFault` | FAIL | `tests/integration/test_slip_weakening.cpp:156`, `EXPECT_EQ(ierr, 0)`. SNES diverges with `DIVERGED_FUNCTION_NANORINF` repeatedly. |

Note: the `TimeDependentSlip` improvement and the `LockedQuasiStatic`
change in failure mode (from immediate `NANORINF` to slow-decay
`MAX_IT`) are both consistent with the BdResidual now firing on
cohesive cells (which it did not under `fault_label`). The cohesive
constraint is producing finite values, but is either too weak or
miscalibrated to drive the displacement jump to the prescribed value.

## 6. PrescribedSlipQuasiStatic max_fault_slip

`max_fault_slip = 2.15241e-23`

(Test failure line: `tests/integration/test_dynamic_rupture_solve.cpp:433`,
`Expected: (max_fault_slip) > (5.0e-4), actual: 2.15241e-23 vs 0.0005`.)

Compared to the audit's reported baseline of `0`, the value is now
non-zero but ~20 orders of magnitude below the target. The SNES log
shows the solve converging at iteration 2 (`Nonlinear solve converged
due to CONVERGED_SNORM_RELATIVE`), so the BdResidual is not pushing
the Lagrange multiplier toward the prescribed jump — the solve is
satisfied with essentially zero slip. The label fix has put the
cohesive geometry into the BdResidual reach, but the residual
contribution from `f0_prescribed_slip` is not large enough to dominate
the manual penalty diagonal added in `addCohesivePenaltyToJacobian` at
the new interior-Lagrange points. This points back to the unresolved
Section B Option B work in `docs/PYLITH_COMPATIBILITY.md`.

## 7. Decision gate

Result: mixed (some pass, some fail). Per the gate:

> If some pass some fail: report which and stop. Do not begin a second
> fix in this session.

Stopping. No further code changes in this session beyond the label
plumbing (TASK B) and the Dockerfile unification (TASK A). No
Section B revisits. No CLAUDE.md test-count edit (the count moved by
zero net: `+TimeDependentSlip` cancels `-CohesiveBdResidual`).

`docs/PYLITH_COMPATIBILITY.md` Section B is **not** marked superseded —
the label fix did not, on its own, drive the Lagrange residual to the
prescribed-slip target.

## 8. Suggested follow-ups (not done in this session)

These are recommendations, not a Session 3 plan. The user will pick:

1. Investigate why `Physics.CohesiveBdResidual` produces NaN under the
   new label. Likely the BdResidual callback is being dispatched on
   label points that are cells or vertices rather than facets. The
   `getOrCreateInterfacesLabel` walks `iDim = 0..dim` and labels
   everything from `pMax` to `pEnd` at each height — including the
   cohesive cells themselves at `iDim = 0`. PyLith uses this label for
   `IntegratorInterface` (per-stratum dispatch), not `DMAddBoundary`.
   For `DMAddBoundary` we may need a label that contains *only* the
   cohesive faces, or to pre-filter the label before passing to PETSc.
2. Once the residual fires cleanly (finite f_norm), re-investigate
   why `max_fault_slip` is 2e-23 instead of 5e-4 — likely the
   penalty-scaled diagonal is over-stiff and the Lagrange
   prescribed-slip residual cannot pull against it.
3. The `LockedQuasiStatic` change from `NANORINF iter 0` to slow
   `MAX_IT` is positive and suggests the SNES Jacobian is now
   well-formed but ill-conditioned. Consider stronger preconditioning
   or `-snes_linesearch_type basic`.

## 9. Files changed in this session

- `Dockerfile` (rewritten — PETSc main from source)
- `Dockerfile.ci` (rewritten — PETSc main from source)
- `Dockerfile.cuda` (rewritten — PETSc main from source + CUDA)
- `Dockerfile.dev` (rewritten — PETSc main from source, drops apt petsc-dev)
- `DOCKER.md` (new — explains the four-file Docker layout)
- `include/core/Simulator.hpp` (added `DMLabel interfaces_label_ = nullptr;` member)
- `src/core/Simulator.cpp`:
  - line 3558: store `interfaces_label_` instead of discarding
  - lines 6470-6480: `fault_constraint` BC now uses
    `interfaces_label_` (with `fault_label` fallback if not yet
    populated)
- `docs/SESSION_2_REPORT.md` (this file)

`FaultMeshManager::splitMeshAlongFault` and
`CohesiveFaultKernel::registerWithDS` were not modified (Rule 6).
The DS/BC ordering in `setupFields()` was not modified (Rule 4).
The `f0_weak_lagrange / g0_weak_lagrange` PetscDS volume callback
from commit 8c1ebbf was left in place per the session brief.
