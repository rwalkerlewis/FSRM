# Session 11 Report

## Goal

Fix the Lagrange FE degree so that Lagrange DOFs are placed on vertices of
the cohesive surface (depth 0), not on the cohesive cells themselves
(depth 3). Session 10 showed the DOFs were landing on depth 3 because
`PetscFECreateDefault(..., PETSC_DETERMINE, ...)` skips the space-degree
setup block in PETSc 3.25 (gated on `if (degree >= 0)`), and FSRM had no
`-lagrange_petscspace_degree 1` option set to override the default.

## Code change

`src/core/Simulator.cpp` (around line 1705). Replaced
`PetscFECreateDefault(..., PETSC_DETERMINE, ...)` with
`PetscFECreateLagrange(..., 1, -1, ...)` so the polynomial degree is set
explicitly at construction time.

Before:

```cpp
PetscFE fe_lagrange = nullptr;
ierr = PetscFECreateDefault(PETSC_COMM_SELF, dim - 1, dim, isSimplex,
                            "lagrange_", PETSC_DETERMINE,
                            &fe_lagrange); CHKERRQ(ierr);
ierr = PetscObjectSetName((PetscObject)fe_lagrange,
                          "lagrange_multiplier_fault"); CHKERRQ(ierr);
```

After:

```cpp
PetscFE fe_lagrange = nullptr;
ierr = PetscFECreateLagrange(PETSC_COMM_SELF, dim - 1, dim, isSimplex,
                             1 /* degree */, -1 /* qorder */,
                             &fe_lagrange); CHKERRQ(ierr);
ierr = PetscObjectSetName((PetscObject)fe_lagrange,
                          "lagrange_multiplier_fault"); CHKERRQ(ierr);
```

This is the exact signature from
`/opt/petsc-main/src/dm/dt/fe/interface/fe.c:2181`.

## Verification results

### Build

Clean rebuild against PETSc 3.25 / `fsrm-ci:main` Docker image. No new
warnings or errors. Log at `/tmp/s11_build.log`.

### Depth distribution (Lagrange-bearing points)

| State | d0 | d1 | d2 | d3 | d4 | total points | DOFs (x3) |
|-------|---:|---:|---:|---:|---:|-------------:|----------:|
| Session 9 / Session 10 baseline (before fix) | 0 | 0 | 0 | 32 | 0 | 32 | 96 |
| Session 11 (after fix)                       | 0 | 25 | 0 | 0 | 0 | 25 | 75 |
| Prompt's expected target                     | 50 | 0 | 0 | 0 | 0 | 50 | 150 |

The fix moved all Lagrange DOFs off the cohesive cells (`d3=32 -> d3=0`),
which is progress, but the DOFs landed on edges (`d1=25`), not on
vertices (`d0`) as anticipated. The total count dropped from 96 to 75 and
the section width shrank accordingly (`total dof=525` vs the prior
section; field 1: `dof=75, free=75`).

No `PetscFEGetNumDof` was probed directly; the depth label walk in
`cacheFaultBoundaryLagrangeDofs` serves the same purpose and is what the
prompt cited as the diagnostic of record.

### Why not d0

Hypothesis: PETSc 3.25 `PetscFECreateLagrange` with `dim=dim-1`
(surface FE, `dim_topo=2`), `Nc=dim`, `k=1` on a cohesive cell inserts
the dual basis on the edges of the two-sided interface rather than on
vertex pairs. This does not match what the Session 11 prompt predicted
based on PyLith and on `$PETSC_DIR/src/dm/impls/plex/tests/ex5.c`; the
prompt assumed `d0>0, d3==0`. Diagnostic work to confirm or refute the
edge-placement hypothesis (by reading back `PetscFEGetNumDof` and
`PetscDualSpaceGetNumDof`) is Session 12 scope.

### PCSVD condition number

| State | condition number | smallest 5 singular values |
|-------|-----------------:|----------------------------|
| Session 9 baseline | `5.36e+17` | `3.72e-08 6.11e-08 7.02e-08 1.08e-07 1.26e-07` |
| Session 10 Step C (pin-at-cells) | `2.36e+18` | `8.33e-09 2.66e-08 6.02e-08 8.61e-08 9.52e-08` |
| Session 11 (after fix) | `2.54e+17` | `7.88e-08 9.55e-08 1.01e-07 1.35e-07 2.02e-07` |

Condition number dropped from `5.36e+17` to `2.53e+17` (about 2x).
Smallest singular value rose from `3.72e-08` to `7.88e-08`. Both shifts
are in the right direction but neither is the multi-order-of-magnitude
improvement the prompt hoped for. The matrix is still rank deficient in
the saddle-point sense.

### PrescribedSlipQuasiStatic SNES trace

The test (`DynamicRuptureSolveTest.PrescribedSlipQuasiStatic`) still
fails with `max_fault_slip = 0` and `sol_norm = 0`. First SNES step
residual is `2.12e+01` and SNES reports
`DIVERGED_MAX_IT iterations 1` for every attempted dt. The TS adaptivity
then retries with smaller dt, all of which diverge the same way. Full
SNES/SVD transcript is at `/tmp/s11_svd.log`.

```
0 TS dt 1. time 0.
    0 SNES Function norm 2.118962305111e+01
        SVD: condition number 2.533601461364e+17, 0 of 297 singular values are (nearly) zero
        SVD: smallest singular values: 7.88e-08 9.55e-08 1.01e-07 1.35e-07 2.02e-07
        SVD: largest singular values : 1.70e+10 1.76e+10 1.76e+10 2.00e+10 2.00e+10
    1 SNES Function norm 1.400732819250e+01
    Nonlinear solve did not converge due to DIVERGED_MAX_IT iterations 1
```

### Fault-test subset (sequential, 16 tests)

Comparison before and after the fix. Both columns are from sequential
ctest runs against the current `local_fix` tree. The "before" column is
what the same subset produced on this branch with the fix reverted (i.e.
with `PetscFECreateDefault` restored) -- captured in
`/tmp/s11_fault_baseline.log`. The "after" column is with
`PetscFECreateLagrange` applied -- captured in `/tmp/s11_fault.log`.

| Test | Before fix | After fix |
|------|:----------:|:---------:|
| Physics.SCEC.TPV5 | FAIL | FAIL |
| Physics.LockedFaultTransparency | FAIL | FAIL |
| Physics.CohesiveBdResidual | pass | pass |
| Functional.DynamicRuptureSetup.MeshSplitting | pass | pass |
| Functional.DynamicRuptureSetup.LockedFault | pass | pass |
| Functional.DynamicRuptureSetup.FaultGeometryFromConfig | pass | pass |
| Functional.DynamicRuptureSetup.SlippingFault | pass | pass |
| Functional.DynamicRuptureSetup.AbsorbingCoexist | pass | pass |
| Integration.PressurizedFractureFEM | FAIL | FAIL |
| Integration.DynamicRuptureSolve.LockedQuasiStatic | FAIL | FAIL |
| Integration.DynamicRuptureSolve.LockedElastodynamic | FAIL | FAIL |
| Integration.DynamicRuptureSolve.PrescribedSlip | FAIL | FAIL |
| Integration.ExplosionFaultReactivation | pass | pass |
| Integration.TimeDependentSlip | FAIL | FAIL |
| Integration.SlippingFaultSolve | FAIL | FAIL |
| Integration.SlipWeakeningFault | FAIL | FAIL |

Identical: 9 failing, 7 passing, before and after. `diff` over the
failure lists is empty. No regressions in the fault subset.

Important: the Session 10 CLAUDE.md baseline claimed only 6 failures
(missing `Physics.SCEC.TPV5`, `Integration.PressurizedFractureFEM`,
`Integration.DynamicRuptureSolve.LockedElastodynamic`). Those three were
already broken at the `local_fix` tip before Session 11 touched the
code; the CLAUDE.md list is stale. See the Session 10 report's own table
which lists the same count.

### Full suite (parallel, 116 tests)

`ctest -j$(nproc)` against the full suite after the fix: **16 failures
out of 116** (`/tmp/s11_test.log`):

```
40  Physics.SCEC.TPV5                                    (fault-related, baseline)
42  Physics.TerzaghiConsolidation                        (HDF5 parallel conflict)
43  Physics.AbsorbingBC                                  (HDF5 parallel conflict)
57  Physics.LockedFaultTransparency                      (fault-related, baseline)
67  Integration.ExplosionSeismogram                      (HDF5 parallel conflict)
91  Integration.PressurizedFractureFEM                   (fault-related, baseline)
93  Integration.DynamicRuptureSolve.LockedQuasiStatic    (fault-related, baseline)
94  Integration.DynamicRuptureSolve.LockedElastodynamic  (fault-related, baseline)
95  Integration.DynamicRuptureSolve.PrescribedSlip       (fault-related, baseline)
97  Integration.TractionBC                               (HDF5 parallel conflict)
98  Integration.TimeDependentSlip                        (fault-related, baseline)
99  Integration.NearFieldCoupled                         (HDF5 parallel conflict)
100 Integration.SlippingFaultSolve                       (fault-related, baseline)
101 Integration.SlipWeakeningFault                       (fault-related, baseline)
109 Integration.HistoricNuclear.Sedan1962                (HDF5 parallel conflict)
111 Integration.HistoricNuclear.NtsPahuteMesa            (HDF5 parallel conflict)
```

Per CLAUDE.md rule ("some tests may fail when run in parallel (`ctest
-j`) due to HDF5 output file conflicts; all pass when run
individually"), the 7 non-fault failures are expected parallel noise and
are not attributable to the FE change. The 9 fault-related failures
match the sequential fault-subset run exactly. No new failures.

## Decision-gate outcome

Re-reading the prompt's gates:

- `max_fault_slip > 5e-4` AND no regressions: NOT MET. Slip is still 0.
- `max_fault_slip == 0` but condition number dropped *significantly*:
  NOT cleanly met. Condition number dropped by ~2x (`5.36e+17` ->
  `2.54e+17`), smallest sigma doubled. These shifts are real but not the
  multi-order-of-magnitude improvement the gate language implied.
- Condition number essentially unchanged: partial match. The spectrum
  moved, but not enough to unlock a solve.
- Regressions: NOT MET. No regressions detected.

Closest match: gate 2 with the caveat that the improvement is modest.
The DOF placement is **not** on vertices (`d0=0`); it landed on edges
(`d1=25`). So Session 10's rim-pin BC is not yet directly applicable in
the form anticipated by the Session 11 prompt; a Session 12 rim-pin
would have to pin Lagrange DOFs at edges on the fault boundary, not at
vertices. The existing `cacheFaultBoundaryLagrangeDofs` diagnostic
already loops over Lagrange-bearing points of any depth via the closure
walk, so adapting to edges is mechanical.

## Session 12 scope proposal

1. Confirm the edge-placement hypothesis by calling
   `PetscFEGetNumDof(fe_lagrange, &numDof)` and
   `PetscDualSpaceGetNumDof(sp, &numDof)` on the newly constructed
   Lagrange FE and logging the array. If
   `numDof = [0, Nc, 0, 0]` (one set per edge), the hypothesis is
   confirmed and Session 10's rim-pin is viable against edges.
2. Decide between:
   a. Accept edge DOFs and port Session 10's rim-pin to edges
      (mechanical change, `cacheFaultBoundaryLagrangeDofs` already
      works on any depth via closure). Then re-run PCSVD and check
      for the 5+ order-of-magnitude drop.
   b. Force vertex DOFs by building a bespoke `PetscDualSpace` with
      `PetscDualSpaceSetNumDof` to pin DOFs at vertices only. More
      invasive but matches PyLith exactly.
3. If neither yields an invertible saddle point, move to the
   fieldsplit + Schur-complement preconditioner path that Session 10
   identified as a final fallback (`-pc_type fieldsplit
   -pc_fieldsplit_type schur`).

## Commit

Committed as: `session 11: use PetscFECreateLagrange with explicit
degree 1 for Lagrange field`. Push to `origin/local_fix` per prompt.

## Artefacts

- `/tmp/s11_build.log` -- clean rebuild against `fsrm-ci:main`.
- `/tmp/s11_svd.log` -- `PrescribedSlipQuasiStatic` with `-pc_svd_monitor`
  showing depth distribution, condition number, and SNES trace.
- `/tmp/s11_fault.log` -- sequential fault-subset ctest run with fix.
- `/tmp/s11_fault_baseline.log` -- sequential fault-subset ctest run
  with fix reverted (for the before/after table).
- `/tmp/s11_test.log` -- parallel full-suite ctest run with fix.
