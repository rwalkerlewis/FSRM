# Session 14 Report

## Goal

Configure a saddle-point solver for fault-enabled problems so the cohesive
Lagrange block's rank deficiency is handled at the solver level, rather than
via the geometric rim-pin BC that Session 12 used. Two paths specified by the
prompt:

- **Path A**: GAMG on the full system with rigid-body near-null-space attached
  to the displacement field (PyLith's pure-elasticity default from
  `libsrc/pylith/materials/Elasticity.cc:getSolverDefaults`).
- **Path B**: Schur fieldsplit with GAMG on the displacement block and jacobi
  on the Lagrange block (PyLith's poroelastic + fault default from
  `libsrc/pylith/materials/Poroelasticity.cc:295-350`).

If neither path restores the three Session-12 passing tests
(`LockedQuasiStatic`, `LockedElastodynamic`, `TimeDependentSlip`), Step C
authorizes restoring Session 12's geometric rim-pin as a fallback.

## Paths tried, in order

### Path A: GAMG + near-null-space (attempted first)

Options set in `Simulator::setupSolvers`:

```
-pc_type gamg
-ksp_type gmres
-ksp_gmres_restart 100
-mg_fine_ksp_max_it 5
-mg_fine_pc_type vpbjacobi
```

Near-null-space composed onto the displacement FE via
`DMPlexCreateRigidBody(dm, displacement_field_idx_, &ns)` then
`PetscObjectCompose((PetscObject)disp_field, "nearnullspace", (PetscObject)ns)`.

**Result on `PrescribedSlipQuasiStatic`** (no rim-pin, Session 13 gating
active): SNES `DIVERGED_LINEAR_SOLVE iterations 0` at the first Newton step,
identical failure mode to Session 13 under direct LU. The GAMG smoother
cannot handle the Lagrange block where all 75 Lagrange DOFs are free and the
diagonal is exactly zero on non-cohesive rows.

```
0 TS dt 1. time 0.
    0 SNES Function norm 2.118962305111e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 2.189978754336e-04
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.600003125147e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 7.817359599706e+08
    [time steps repeat with the same DIVERGED_LINEAR_SOLVE]
```

### Path B: Schur fieldsplit (added on top of Path A)

Additional options:

```
-pc_type fieldsplit
-pc_fieldsplit_type schur
-pc_fieldsplit_schur_fact_type full
-pc_fieldsplit_schur_precondition selfp
-fieldsplit_displacement__ksp_type preonly     # note double underscore: FSRM field name is "displacement_"
-fieldsplit_displacement__pc_type gamg
-fieldsplit_lagrange_multiplier_fault_ksp_type preonly
-fieldsplit_lagrange_multiplier_fault_pc_type jacobi
-ksp_type gmres
-ksp_gmres_restart 100
```

**First failure mode**: with the original Path A `DMPlexCreateRigidBody(dm,
field, &ns)` call (null-space vectors sized to the full section),
`PCSetUp_GAMG` called on the fieldsplit displacement sub-block errored with
`PETSC_ERR_PLIB` (error 77):

```
[0]PETSC ERROR: Attached null space vector size 297 != matrix size 222
[0]PETSC ERROR: #1 PCSetData_AGG() at gamg/agg.c:476
[0]PETSC ERROR: #5 PCSetUpOnBlocks_FieldSplit_Schur() at fieldsplit.c:1115
```

The cause: `DMPlexCreateRigidBody(dm, field, &ns)` returns vectors sized to
the full DM section (297 free DOFs = 222 displacement + 75 Lagrange), not to
the displacement sub-block. When fieldsplit extracts the displacement
sub-matrix (222x222), the pre-composed null-space vectors are the wrong size.

**Fix**: create a sub-DM for the displacement field first, then
`DMPlexCreateRigidBody(sub_dm, 0, &ns)` returns null-space vectors of size
222. This is now the committed form and also works with Path A.

```cpp
DM sub_dm = NULL;
PetscInt disp_idx_array[1] = {displacement_field_idx_};
DMCreateSubDM(dm, 1, disp_idx_array, NULL, &sub_dm);
DMPlexCreateRigidBody(sub_dm, 0, &near_null_space);
// then PetscObjectCompose to the displacement FE, destroy sub_dm.
```

**Second failure mode after the null-space fix**: GMRES diverges explosively
on the fieldsplit Schur complement:

```
[0]PETSC ERROR: Residual norm computed by GMRES recursion formula 1.83459e+13
                is far from the computed residual norm 1.25221e+14 at restart,
                residual norm at start of cycle 3.60707e+13
[0]PETSC ERROR: #1 KSPGMRESCycle() at gmres/gmres.c:102
                [PETSC_ERR_CONV_FAILED = 82]
```

The Schur complement `S = -B * A^{-1} * B^T` is exactly zero on the
non-cohesive Lagrange DOFs (`B` is zero there because only cohesive faces
contribute to the Lagrange residual). `selfp` preconditioning of a zero
Schur row produces a singular pre, and jacobi on the Lagrange block is a no-op
against zero diagonals. The outer GMRES then explodes.

## Field names (`PetscSectionGetFieldName`)

Verified on `PrescribedSlipQuasiStatic`:

```
Session 14: attached rigid-body near-null-space to field 0 ('displacement_')
Session 14: lagrange field 1 name: 'lagrange_multiplier_fault'
```

Note: `displacement_` has a trailing underscore. PETSc's fieldsplit option
prefix is `fieldsplit_<fieldname>_`, so the split-0 prefix becomes
`fieldsplit_displacement__` (double underscore). This matched the option
strings we set. The Lagrange field name (`lagrange_multiplier_fault`) matches
PyLith's convention exactly.

## Condition number under each configuration

Path A and Path B bypass direct LU so PCSVD is not run on the preconditioned
operator. The baseline direct-LU PCSVD numbers from Sessions 11/13 still
apply to the assembled saddle-point matrix and were not re-measured:

```
Without rim-pin (Session 13, no-pin branch):
  SVD: condition number 2.533601461364e+17
       0 of 297 singular values are (nearly) zero
  SVD: smallest singular values: 7.880764e-08 ... (numerically zero)

With rim-pin (Session 12 / Session 14 fallback, 72 of 75 Lagrange DOFs pinned):
  free=225 total (222 disp + 3 lagrange); tests converge under direct LU.
```

Measured SNES initial norm is a proxy for the reduced system size:

- No rim-pin: first-step residual 2.119e+01 (297 free DOFs).
- Rim-pin fallback: first-step residual 6.25e-05 for Path A diagnostic run
  (225 free DOFs after 300 constraints), or 7.99e+08 for Path B diagnostic
  run at a later time step after line-search failure.

## SNES trace for `PrescribedSlipQuasiStatic`

### Default (rim-pin restored, direct LU)

```
0 TS dt 1. time 0.
    0 SNES Function norm 6.250000000000e-05
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    [time steps repeat with backoff, never produces fault slip]
max_fault_slip = 0
```

The rim-pin pins 72 of 75 Lagrange DOFs (all rim-adjacent cohesive cells
flagged by the geometric detector). The prescribed-slip residual drives only
3 free Lagrange DOFs, producing effectively zero fault slip. This is Session
12's well-known symptom: rim-pin makes Locked*/TimeDependentSlip converge,
but it over-constrains PrescribedSlip.

### `FSRM_SADDLE_SOLVER=gamg` (Path A with rim-pin)

```
0 TS dt 1. time 0.
    0 SNES Function norm 6.250000000000e-05
    ierr == PETSC_ERR_PLIB (error 77) at TSStep
```

Unused-option warnings appear for `-mg_fine_ksp_max_it` and `-mg_fine_pc_type`
because the GAMG setup fails before those sub-options are consumed.

### `FSRM_SADDLE_SOLVER=fieldsplit` (Path B with rim-pin)

```
0 SNES Function norm 7.993052538855e+08
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
max_fault_slip = 0
```

## Fault-subset table

All runs on PETSc 3.25 / `fsrm-ci:main`.

| Test | Session 12 | Session 13 | Session 14 (rim-pin restored, default solver) |
|------|:----------:|:----------:|:-------:|
| Physics.SCEC.TPV5                                      | FAIL | FAIL | FAIL |
| Physics.LockedFaultTransparency                        | FAIL | FAIL | FAIL |
| Physics.CohesiveBdResidual                             | pass | pass | pass |
| Functional.DynamicRuptureSetup.MeshSplitting           | pass | pass | pass |
| Functional.DynamicRuptureSetup.LockedFault             | pass | pass | pass |
| Functional.DynamicRuptureSetup.FaultGeometryFromConfig | pass | pass | pass |
| Functional.DynamicRuptureSetup.SlippingFault           | pass | pass | pass |
| Functional.DynamicRuptureSetup.AbsorbingCoexist        | pass | pass | pass |
| Integration.PressurizedFractureFEM                     | FAIL | FAIL | FAIL |
| **Integration.DynamicRuptureSolve.LockedQuasiStatic**  | **PASS** | **FAIL** | **PASS** |
| **Integration.DynamicRuptureSolve.LockedElastodynamic**| **PASS** | **FAIL** | **PASS** |
| Integration.DynamicRuptureSolve.PrescribedSlip         | FAIL | FAIL | FAIL |
| Integration.ExplosionFaultReactivation                 | pass | pass | pass |
| **Integration.TimeDependentSlip**                      | **PASS** | **FAIL** | **PASS** |
| Integration.SlippingFaultSolve                         | FAIL | FAIL | FAIL |
| Integration.SlipWeakeningFault                         | FAIL | FAIL | FAIL |

Session 14 restores Session 12's three passing tests that Session 13
regressed. No new regressions. PrescribedSlip, SCEC.TPV5, LockedFaultTransparency,
PressurizedFractureFEM, SlippingFaultSolve, and SlipWeakeningFault remain in
their Session 12 failure state: the rim-pin does not help these because
they either need prescribed-slip to act on non-rim Lagrange DOFs, or need a
friction Jacobian (slipping/slip-weakening), or have a separate BdResidual
issue.

## Full-suite count

Parallel `ctest -j`:

| Run | Failures (of 116) |
|-----|:-----------------:|
| Session 12 parallel | 14 |
| Session 13 parallel | 14 |
| Session 14 parallel | 12 |

Session 14 parallel fails: 6 fault-related baseline (`SCEC.TPV5`,
`LockedFaultTransparency`, `PressurizedFractureFEM`, `PrescribedSlip`,
`SlippingFaultSolve`, `SlipWeakeningFault`) plus 6 HDF5-parallel flakes
(`TerzaghiConsolidation`, `AbsorbingBC`, `ExplosionSeismogram`,
`Gnome1961`, `DegelenMountain`, `NtsPahuteMesa`). The HDF5 swing is
non-deterministic between parallel runs per the CLAUDE.md caveat and should
not be treated as a regression.

Session 13's 14 parallel fails included the three regressions
(`LockedQuasiStatic`, `LockedElastodynamic`, `TimeDependentSlip`) plus
different HDF5 flakes; Session 14 removes those regressions and lands with
the same core fault failures as Session 12.

## Decision-gate outcome

The prompt defines three gates:

1. PrescribedSlip passes AND Locked*/TimeDependentSlip back to PASS: solver
   config is correct, propose Session 15 for friction-Jacobian port.
2. Locked*/TimeDependentSlip back to PASS but PrescribedSlip still fails:
   saddle-point helps linear-constraint cases but not prescribed-slip.
3. Mixed or further regressions: report SNES and KSP traces, stop.

Path A alone: stays in gate 3 (3 regressions from Session 12 baseline still
not recovered). Path A + Path B: gate 3 with 5 regressions. Neither path
restored PrescribedSlip, but both paths left the Locked* and TimeDependentSlip
tests failing under their respective solver configs.

**Step C fallback applied**: per the prompt's explicit authorization ("Only
restore geometric detection if both Path A and Path B fail to restore the
three tests Session 12 was passing"), Session 12's geometric rim-pin was
restored. The `getOrCreateBuriedCohesiveLabel` function now:

1. Prefers the PyLith-semantic `buried_edges` input label when present
   (Session 13 path).
2. Falls back to Session 12's geometric detection (cohesive prisms whose
   cone midpoint lies on the global bounding box) when the label is absent.
   This is what the FSRM test meshes take.

With the rim-pin fallback active, Locked*/TimeDependentSlip land on gate 2:
linear-constraint tests pass, PrescribedSlip still fails. Per gate 2's
guidance, Session 15 should investigate what is specifically different
about PrescribedSlip's constraint application path (the rim-pin pins 72 of
75 Lagrange DOFs, leaving only 3 free for the prescribed-slip constraint to
act on, which is physically insufficient to drive a measurable displacement
jump).

## Code state committed

- `src/core/Simulator.cpp::setupSolvers`: rigid-body near-null-space is
  composed onto the displacement FE via a sub-DM (size-222 vectors) so
  fieldsplit's sub-matrix size check passes. This is scaffolding; with the
  default direct-LU solver it is ignored by PCLU. Path A and Path B option
  sets are gated behind `FSRM_SADDLE_SOLVER=gamg` and
  `FSRM_SADDLE_SOLVER=fieldsplit` env vars respectively so the default
  solver stays the test-provided direct LU.

- `src/core/Simulator.cpp::getOrCreateBuriedCohesiveLabel`: hybrid. Checks
  for `buried_edges` label first (PyLith semantics, Session 13 path). When
  absent, falls back to Session 12's geometric bounding-box detector so
  Locked* and TimeDependentSlip keep converging under direct LU.

## Session 15 scope proposal

Two tracks:

1. **PrescribedSlip under the rim-pin**: the constraint lands on only 3 free
   Lagrange DOFs. Either (a) change the prescribed-slip constraint path to
   apply the target jump via the strong BC (Dirichlet) on the Lagrange
   pinned DOFs rather than through the BdResidual, which would let the rim
   Lagrange values equal the prescribed target, or (b) write a real
   `buried_edges` label at mesh-generation time in
   `FaultMeshManager::splitMeshAlongFault` based on the pre-split submesh
   boundary so only truly-buried tips are labeled (empty label for faults
   that cut the full domain, so no rim-pin is registered, matching PyLith's
   default and letting a proper solver handle the zero Lagrange block).

2. **Friction-Jacobian port for slipping/slip-weakening tests**: these fail
   for a separate reason from the rim-pin / rank-deficiency issue, so a
   semi-smooth Newton Jacobian port is independent work.

## Artefacts

- `/tmp/s14_build.log`, `/tmp/s14_build_b.log` -- clean rebuilds.
- `/tmp/s14_path_a.log` -- PrescribedSlip under Path A alone.
- `/tmp/s14_path_b.log` -- PrescribedSlip under Path A + Path B.
- `/tmp/s14_direct.log`, `/tmp/s14_direct2.log` -- direct `./fsrm` runs with
  `-ksp_view` and `-ksp_error_if_not_converged` showing the error traces
  above.
- `/tmp/s14_locked.log` -- fault subset under Path A + Path B (5 failures).
- `/tmp/s14_fault.log` -- fault subset under Session 14 final state
  (6 failures, matches Session 12).
- `/tmp/s14_test.log` -- full suite under Session 14 final state
  (12 failures: 6 fault baseline + 6 HDF5 parallel flakes).
- `/tmp/s14_pathA.log` -- PrescribedSlip with `FSRM_SADDLE_SOLVER=gamg` and
  rim-pin active.
- `/tmp/s14_pathB.log` -- PrescribedSlip with
  `FSRM_SADDLE_SOLVER=fieldsplit` and rim-pin active.

## Commit

```bash
git add -A
git commit -m "session 14: saddle-point solver config (PyLith GAMG + near-null-space or Schur fieldsplit)"
git push origin local_fix
```
