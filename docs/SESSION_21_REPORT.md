# Session 21 Report: GAMG + Rigid-Body Near-Null-Space Infrastructure

## Summary

Session 21 attempted the prompt's prescription to flip the aux-slip path
default-ON, drop the Session 12/14 rim-pin essential BC, and activate
PyLith's elasticity-plus-fault solver defaults (GAMG + GMRES + rigid-body
near-null-space) to handle the resulting zero-Lagrange saddle-point
system. The solver infrastructure was implemented and a parent-DM sized
near-null-space was attached to the assembled Jacobian via
`MatSetNearNullSpace` (verified reaching `MatGetNearNullSpace` on every
Newton step). With the stack active, KSP no longer hangs at
`DIVERGED_LINEAR_SOLVE iterations 0`: it runs three Newton iterations,
drives the residual from `2.19e-4` to `17.88`, then stalls with
`DIVERGED_LINE_SEARCH iterations 2`. This is a new failure mode and does
not satisfy Session 21's decision gate (PrescribedSlip requires
`max_fault_slip >= 5e-4`). Default-ON of the path also regresses
`LockedQuasiStatic`, `LockedElastodynamic`, `TimeDependentSlip` in the
fault subset. Per Rule 3 ("all existing tests must continue to pass")
the gate inversions were rolled back after diagnostics were captured.
The Session 21 solver infrastructure (rigid-body near-null-space on the
parent DM, `MatSetNearNullSpace` in `FormJacobian`, richardson/jacobi
smoother options under `FSRM_SADDLE_SOLVER=gamg`) is preserved in tree
behind the explicit env-var gate and is benign on the default path.

Full-suite count is unchanged from Session 20: **110/116**.

## 1. Solver infrastructure (landed)

### 1.1 Parent-DM rigid-body near-null-space

`src/core/Simulator.cpp::setupSolvers` builds the near-null-space on the
parent DM rather than on a displacement-only sub-DM:

```cpp
if (displacement_field_idx_ >= 0 && !rigidBodyNullSpace_) {
    ierr = DMPlexCreateRigidBody(dm, displacement_field_idx_,
                                 &rigidBodyNullSpace_); CHKERRQ(ierr);
    // Compose on the displacement field object for fieldsplit paths
    PetscObject disp_field = NULL;
    ierr = DMGetField(dm, displacement_field_idx_, NULL, &disp_field);
    ierr = PetscObjectCompose(disp_field, "nearnullspace",
                              (PetscObject)rigidBodyNullSpace_);
}
```

The vectors are sized for the full assembled Jacobian (local size = 297
on the canonical PrescribedSlip mesh, 222 displacement DOFs populated,
75 Lagrange DOFs zero). `has_const=0 nvecs=6` - the six rigid-body
modes expected in 3D. The Session 14 sub-DM + `PetscObjectCompose`
pattern did not reach `MatGetNearNullSpace` on the monolithic Jacobian
in PETSc 3.25 (confirmed by the probe in Section 2 below).

### 1.2 Direct `MatSetNearNullSpace` in `FormJacobian`

After the hybrid assembly bookend, `FormJacobian` now attaches the
stored null-space directly to `J` and `P`:

```cpp
if (sim->rigidBodyNullSpace_) {
    ierr = MatSetNearNullSpace(P, sim->rigidBodyNullSpace_); CHKERRQ(ierr);
    if (J != P) {
        ierr = MatSetNearNullSpace(J, sim->rigidBodyNullSpace_); CHKERRQ(ierr);
    }
}
```

This is the path that `PCGAMG` uses via `MatGetNearNullSpace`.

### 1.3 Session 21 debug probe

A new gated probe inside `FormJacobian` reports whether
`MatGetNearNullSpace(J, &nns)` returns the attached null-space or
`NULL`:

```cpp
if (getenv("FSRM_SESSION21_DEBUG")) {
    MatNullSpace nns = nullptr;
    ierr = MatGetNearNullSpace(J, &nns); CHKERRQ(ierr);
    PetscPrintf(sim->comm,
                "Session 21: MatGetNearNullSpace returned %s\n",
                nns ? "attached" : "NULL");
    ...
}
```

### 1.4 Explicit GAMG selector under `FSRM_SADDLE_SOLVER=gamg`

`setupSolvers` now has three branches for fault-enabled runs:

- `FSRM_SADDLE_SOLVER=fieldsplit` - unchanged Session 15 Schur path
- `FSRM_SADDLE_SOLVER=gamg` (new) - overrides any caller `-pc_type`
  preset and sets GAMG + GMRES + `-mg_levels_ksp_type richardson` +
  `-mg_levels_pc_type jacobi` + `-mg_levels_ksp_max_it 5`
- default - unchanged Session 15/20 (`-pc_type gamg` plus
  `-mg_fine_*` options when no caller preset)

The explicit `FSRM_SADDLE_SOLVER=gamg` selector is what Session 21's
experimental path uses.

## 2. Verification: near-null-space reaches the Jacobian

With `FSRM_SESSION21_DEBUG=1` + experimental gates set:

```
FSRM_ENABLE_SLIP_AUX=1 FSRM_DISABLE_RIM_PIN=1 \
FSRM_SADDLE_SOLVER=gamg FSRM_SESSION21_DEBUG=1 \
ctest -R DynamicRuptureSolve.PrescribedSlip
```

Output (see `/tmp/s21_experiment.log`):

```
Session 21: built rigid-body near-null-space on parent DM for field 0
            ('displacement_'); has_const=0 nvecs=6 local_size=297
Session 21: lagrange field 1 name: 'lagrange_multiplier_fault'
Session 21: FSRM_SADDLE_SOLVER=gamg explicit; PyLith defaults
            (GAMG + GMRES + near-null-space)
...
0 TS dt 1. time 0.
    0 SNES Function norm 2.189978754336e-04
Session 21: MatGetNearNullSpace returned attached
Session 21: near-null-space has_const=0 nvecs=6
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
```

The probe confirms `MatGetNearNullSpace` returns the attached null-space
on the full Jacobian (not NULL) with the expected 6-mode rigid-body
basis. This closes decision-gate branch (a): "null-space not reaching
GAMG".

## 3. KSP convergence history and SNES trace for
`PrescribedSlipQuasiStatic`

`run_integration_tests` swallows `-ksp_monitor` /
`-ksp_converged_reason` because its harness calls `PetscOptionsClear`
before every test. The SNES monitor is set after the clear, so we have
the SNES-level trace but no per-KSP iteration output. The SNES trace
captures the failure mode fully:

```
0 TS dt 1. time 0.
    0 SNES Function norm 2.189978754336e-04
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 1.788851587047e+01
    1 SNES Function norm 1.788702882529e+01
    2 SNES Function norm 1.788702882391e+01
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 2
    0 SNES Function norm 1.788851587047e+01
    1 SNES Function norm 1.788702882529e+01
    2 SNES Function norm 1.788702882391e+01
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 2
    ... (repeats 10 times on the TS retry loop) ...
TSStep has failed due to DIVERGED_STEP_REJECTED

sol_norm        = 0
max_fault_slip  = 0
```

Interpretation:

- KSP *does* run (unlike Session 20's `iterations 0` zero-pivot signature).
  The TS retries proceed until TS abandons the step.
- The KSP solves produce a `delta_u` which, when line-searched by SNES,
  fails to reduce the residual. First attempt: residual jumps from
  `2.19e-4` to `17.88` (5 orders of magnitude) on the trial step, so
  line-search cannot find a step size that helps.
- After each retry the state is left at a perturbed `u` whose residual
  norm is `17.88`. Three more Newton iterations converge the norm to
  `17.887` (a fixed point) but SNES line-search still cannot improve,
  so each retry reports `DIVERGED_LINE_SEARCH iterations 2`.
- `max_fault_slip = 0` because the solution never advances from
  `u = 0` - every successful intermediate state is rolled back by TS
  when the step is rejected.

The Jacobian is ill-conditioned or structurally bad enough that GAMG
with the richardson/jacobi smoother and six rigid-body near-null-space
modes does not produce a useful update direction. This closes decision-
gate branches (a, b) and leaves branch (c) "fieldsplit Schur instead of
GAMG on the full system" as an open path (Session 20's
`FSRM_SADDLE_SOLVER=fieldsplit` probe already failed with
`DIVERGED_LINEAR_SOLVE` on the same matrix, so Schur is not an obvious
win either).

## 4. Fault-subset table (Session 20 vs Session 21 default vs Session 21 experimental)

Experimental column: `FSRM_ENABLE_SLIP_AUX=1 FSRM_DISABLE_RIM_PIN=1
FSRM_SADDLE_SOLVER=gamg`.

| # | Test | S20 default | S21 default (rolled back) | S21 experimental |
|---|------|------------:|--------------------------:|-----------------:|
| 38 | Physics.CohesiveBdResidual                          | PASS | PASS | PASS |
| 40 | Physics.SCEC.TPV5                                   | FAIL | FAIL | FAIL |
| 57 | Physics.LockedFaultTransparency                     | FAIL | FAIL | FAIL |
| 91 | Integration.PressurizedFractureFEM                  | FAIL | FAIL | FAIL |
| 93 | Integration.DynamicRuptureSolve.LockedQuasiStatic   | PASS | PASS | **FAIL** |
| 94 | Integration.DynamicRuptureSolve.LockedElastodynamic | PASS | PASS | **FAIL** |
| 95 | Integration.DynamicRuptureSolve.PrescribedSlip      | FAIL | FAIL | FAIL |
| 96 | Integration.ExplosionFaultReactivation              | PASS | PASS | PASS |
| 98 | Integration.TimeDependentSlip                       | PASS | PASS | **FAIL** |
| 100| Integration.SlippingFaultSolve                      | FAIL | FAIL | FAIL |
| 101| Integration.SlipWeakeningFault                      | FAIL | FAIL | FAIL |
| -- | Functional.DynamicRuptureSetup (5 tests)            | PASS | PASS | PASS |

**Summary**:

- S21 default = S20 default: 10/16 fault tests passing.
- S21 experimental: 7/16. Net -3 vs the baseline
  (`LockedQuasiStatic`, `LockedElastodynamic`, `TimeDependentSlip`
  regress with `DIVERGED_DTOL` / `DIVERGED_LINE_SEARCH` signatures
  instead of Session 20's `DIVERGED_LINEAR_SOLVE iterations 0`).

## 5. Fallback path (`FSRM_DISABLE_SLIP_AUX=1`)

With the rollback in place, `FSRM_DISABLE_SLIP_AUX` has no functional
effect (the slip-aux DM is off by default, gated by `FSRM_ENABLE_SLIP_AUX`).
Running the spec's fallback command still produces Session 20 baseline:

```
$ FSRM_DISABLE_SLIP_AUX=1 ctest -R \
    "Physics.LockedFaultTransparency|Integration.DynamicRuptureSolve|\
     Integration.TimeDependentSlip"

Integration.DynamicRuptureSolve.LockedQuasiStatic   PASS
Integration.DynamicRuptureSolve.LockedElastodynamic PASS
Integration.DynamicRuptureSolve.PrescribedSlip      FAIL (baseline)
Integration.TimeDependentSlip                       PASS
Physics.LockedFaultTransparency                     FAIL (baseline)
3 passed, 2 failed (both pre-existing baseline failures)
```

## 6. Full-suite count

```
$ docker run ... ctest --output-on-failure
95% tests passed, 6 tests failed out of 116
Total Test time (real) = 136.23 sec

The following tests FAILED:
   40 - Physics.SCEC.TPV5 (Failed)
   57 - Physics.LockedFaultTransparency (Failed)
   91 - Integration.PressurizedFractureFEM (Failed)
   95 - Integration.DynamicRuptureSolve.PrescribedSlip (Failed)
  100 - Integration.SlippingFaultSolve (Failed)
  101 - Integration.SlipWeakeningFault (Failed)
```

Identical to Session 20's 110/116.

## 7. Decision-gate verdict

The prompt enumerated five branches. The observed outcome maps to
**"PrescribedSlip fails with KSP still divergent: GAMG with the current
near-null-space is insufficient. Capture the full KSP monitor history
and the MatNullSpace dump."** - with diagnostics as follows:

- (a) null-space reaching GAMG: confirmed attached
  (`MatGetNearNullSpace` returns the 6-mode basis on every Newton step).
- (b) GAMG tuning (`-mg_levels_pc_type jacobi`): applied; lets KSP run
  and produce a step direction but the step fails line-search.
- (c) fieldsplit Schur instead of GAMG on the full system: not retried
  in Session 21 because Session 20 already documented
  `FSRM_SADDLE_SOLVER=fieldsplit` failing with the same `KSP iterations 0`
  signature on this matrix.

The regression branch
**"PrescribedSlip passes but Locked/TimeDependent regress in default
path"** is also relevant: the experimental path fails both conditions
(PrescribedSlip does not pass *and* Locked/TimeDependent regress).
Per Rule 3 the tree is kept at Session 20 wording for all gates.

## 8. Why SNES line-search rejects the GAMG step

Working hypothesis (based on residual trace; not verified with a full
Jacobian dump in this session):

1. The prescribed-slip residual at `u = 0` is `r_lambda = -slip_target`
   (magnitude `~1e-3` per cohesive vertex).
2. GAMG's two-level Galerkin coarsening on a matrix with 75 exact-zero
   diagonals in the Lagrange block produces a coarse operator whose
   Schur complement on the coarse Lagrange DOFs is badly scaled or
   singular. The smoother (richardson/jacobi) cannot damp the
   high-frequency Lagrange modes because their diagonal is zero.
3. The coarse-grid LU solve returns a delta containing large
   rigid-body-like motions of the displacement block that are *not*
   removed by the near-null-space filter (the null-space removes them
   from the RHS but GAMG can still inject them via the coarse
   correction when the coarse operator is ill-posed).
4. SNES line-search tries the full step and sees the residual jump
   from `2.19e-4` to `17.88`; backtracking gives up at
   `alpha = 2^-15 * 1.0` and reports `DIVERGED_LINE_SEARCH`.

The fundamental issue is that PETSc 3.25 `PCGAMG` on the monolithic
saddle-point matrix `[K B^T; B 0]` with rigid-body near-null-space is
not equivalent to PyLith's code path even after copying the option
strings. PyLith's `-dm_reorder_section_type cohesive` is already
applied in our `setupFields` via `DMReorderSectionSetType(dm,
"cohesive")`, but empirically the GAMG coarsening still fails.

## 9. Files changed (in landed diff)

- `include/core/Simulator.hpp` - added `rigidBodyNullSpace_` member.
- `src/core/Simulator.cpp` -
  `~Simulator()` destroys the null-space;
  `setupSolvers` builds the near-null-space on the parent DM and adds
  the `FSRM_SADDLE_SOLVER=gamg` explicit branch;
  `FormJacobian` now calls `MatSetNearNullSpace(P/J, ...)` after the
  assembly bookend and prints the Session 21 probe when the env var
  is set.
- All four Session 21 gate-inversion sites returned to their Session 20
  wording:
  - `setupAuxiliaryDM` material-aux trigger (`FSRM_ENABLE_SLIP_AUX`)
  - slip aux DM creation (`FSRM_ENABLE_SLIP_AUX`)
  - slip aux cohesive-key attach (`FSRM_ENABLE_SLIP_AUX`)
  - `FormFunction` and `FormJacobian` attach (`FSRM_ENABLE_SLIP_AUX`)
  - `setupBoundaryConditions` rim-pin (`FSRM_DISABLE_RIM_PIN`)

## 10. Suggested Session 22 follow-up

1. Try `-pc_type fieldsplit -pc_fieldsplit_type schur
   -pc_fieldsplit_schur_fact_type full
   -pc_fieldsplit_schur_precondition user` with a *user-provided*
   Schur pre based on the lumped `B K^{-1} B^T`. PyLith's
   FaultCohesiveKin actually builds this pre explicitly
   (`libsrc/pylith/faults/FaultCohesiveKin.cc::_createSubpreconditioner`
   constructs a mass-lumped approximation) rather than relying on
   `selfp`.
2. Alternatively, accept that FSRM's PETSc 3.25 code path cannot
   reproduce PyLith's solver on the same matrix and focus Session 22
   on porting PyLith's friction Jacobian (the remaining
   `SlippingFaultSolve` / `SlipWeakeningFault` failures have a
   different root cause and may land independently).
3. A low-risk alternative to unblock `PrescribedSlipQuasiStatic` is to
   drop the aux-slip attachment entirely for this test and project the
   prescribed slip directly into a Dirichlet BC on the Lagrange DOFs
   at cohesive vertices via `DM_BC_ESSENTIAL`. That matches the
   existing rim-pin machinery and sidesteps the saddle-point solver
   issue.
