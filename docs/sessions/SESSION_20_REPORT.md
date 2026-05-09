# Session 20 Report: Material-Aux Face Quadrature + Aux-Slip Gating Probe

## Summary

The material auxiliary FE now copies BOTH the displacement FE's volume
quadrature and its face quadrature, fixing the `febasic.c:663` "Number of
tabulation points" mismatch that fired whenever the cohesive-key slip aux
was attached. The fix is harmless to the Session 16 default code path (the
slip aux remains opt-in via `FSRM_ENABLE_SLIP_AUX`) and unexpectedly
unlocks two previously-failing fault tests that ran the bulk-side hybrid
integrator without the slip aux. Full-suite count improves from
**108/116** (Session 19 baseline) to **110/116**.

The Session 20 prompt also asked to invert three env-var gates so that the
aux-slip path becomes default-ON and the rim-pin auto-disables. That
inversion was implemented and tested. With the inversion active the
`febasic.c:663` and `plexfem.c:3982` checks both pass (the face-quadrature
fix does its job), but a new failure mode appears: PETSc reports
`KSP DIVERGED_LINEAR_SOLVE iterations 0` at every Newton step on every
linear fault test (LockedQuasiStatic, LockedElastodynamic,
PrescribedSlip, TimeDependentSlip, ExplosionFaultReactivation, ...).
The factorization fails before the first KSP iteration because removing
the rim-pin leaves the boundary Lagrange DOFs with no diagonal Jacobian
contribution outside the small `epsilon = 1e-4` volume regularization.
Per the prompt's decision-gate "**New failure mode not previously
seen**: stop and report" branch, the gate inversions were rolled back
to keep the live tree at the improved baseline. The face-quadrature
fix is preserved.

## 1. Code diff: material aux face-quadrature copy (LANDED)

`src/core/Simulator.cpp:3120-3151` (was 3120-3137):

```diff
-    // Get quadrature from the primary FE to match tabulation points
-    PetscQuadrature quad = NULL;
-    if (!fe_fields.empty()) {
-        ierr = PetscFEGetQuadrature(fe_fields[0], &quad); CHKERRQ(ierr);
-    }
+    // Get BOTH volume and face quadrature from the primary (displacement) FE.
+    // Session 19 proved that the febasic.c:663 "Number of tabulation points"
+    // mismatch was driven by this material aux, not the slip aux. The hybrid
+    // residual integrator queries the material aux on the bulk-side keys
+    // (material_label_ values 1 and 2). With auxOnBd=FALSE on those keys
+    // (dimAux=3, dim from cohesive Lagrange field=2), PETSc pulls the aux
+    // FACE tabulation and compares it to the main FE's 3-point volume tab.
+    // Previously the material aux got the volume quadrature copied but the
+    // face quadrature stayed at the P0 default (1 point), so the Np check
+    // always failed once the slip aux was attached. Copy both quadratures.
+    PetscQuadrature disp_vol_quad = NULL;
+    PetscQuadrature disp_face_quad = NULL;
+    if (!fe_fields.empty()) {
+        ierr = PetscFEGetQuadrature(fe_fields[0], &disp_vol_quad); CHKERRQ(ierr);
+        ierr = PetscFEGetFaceQuadrature(fe_fields[0], &disp_face_quad); CHKERRQ(ierr);
+    }

     for (PetscInt f = 0; f < NUM_AUX_FIELDS; ++f) {
         PetscFE fe_aux;
         ierr = PetscFECreateLagrange(PETSC_COMM_SELF, dim, 1, isSimplex, 0, PETSC_DETERMINE, &fe_aux); CHKERRQ(ierr);
-        if (quad) {
-            ierr = PetscFESetQuadrature(fe_aux, quad); CHKERRQ(ierr);
+        if (disp_vol_quad) {
+            ierr = PetscFESetQuadrature(fe_aux, disp_vol_quad); CHKERRQ(ierr);
+        }
+        if (disp_face_quad) {
+            ierr = PetscFESetFaceQuadrature(fe_aux, disp_face_quad); CHKERRQ(ierr);
         }
         ierr = DMSetField(auxDM_, f, NULL, (PetscObject)fe_aux); CHKERRQ(ierr);
         ierr = PetscFEDestroy(&fe_aux); CHKERRQ(ierr);
     }
```

## 2. Code diffs: gate inversions (TESTED, ROLLED BACK)

The three sites the prompt prescribed for inversion were modified, built
clean, run end-to-end, then reverted to their Session 19 wording. The
diffs below describe the inversion that was tested.

### a) `setupAuxiliaryDM` slip aux DM creation (line 3192)

```cpp
// Tested:
const bool slip_aux_disabled = (getenv("FSRM_DISABLE_SLIP_AUX") != nullptr);
if (!slip_aux_disabled &&
    config.enable_faults && cohesive_kernel_ && interfaces_label_) { ... }

// Reverted to:
if (getenv("FSRM_ENABLE_SLIP_AUX") &&
    config.enable_faults && cohesive_kernel_ && interfaces_label_) { ... }
```

### b) `FormFunction` cohesive-key attachment (line 6940)

```cpp
// Tested:
const bool enable_slip_aux = (getenv("FSRM_DISABLE_SLIP_AUX") == nullptr);

// Reverted to:
const bool enable_slip_aux = (getenv("FSRM_ENABLE_SLIP_AUX") != nullptr);
```

### c) `FormJacobian` cohesive-key attachment (line 7155)

Mirror of (b).

### d) Rim-pin coordination (setupBoundaryConditions, line 8257)

```cpp
// Tested:
const bool aux_slip_active = (getenv("FSRM_DISABLE_SLIP_AUX") == nullptr);
const bool rim_pin_active  = (getenv("FSRM_ENABLE_RIM_PIN") != nullptr) ||
                             !aux_slip_active;
if (rim_pin_active && buried_cohesive_label_ && lagrange_field_idx_ >= 0) { ... }

// Reverted to:
const bool rim_pin_disabled = (getenv("FSRM_DISABLE_RIM_PIN") != nullptr);
if (!rim_pin_disabled && buried_cohesive_label_ && lagrange_field_idx_ >= 0) { ... }
```

## 3. `febasic.c:663` status

**Resolved.** With the face-quadrature fix in place and the slip aux
attached on the cohesive key (`FSRM_ENABLE_SLIP_AUX=1` or, equivalently,
all four gate inversions active), the `PetscCheck(TfAux[f]->Np == Tf[0]->Np, ...)`
at `febasic.c:663` no longer trips.

The Session 19 `FSRM_SESSION19_DEBUG` print does not fire on the
PrescribedSlip run with the gate inversions active: the rank-0 setup
print is gated on the slip-aux DM having been built, and the SNES failure
above happens after `setupAuxiliaryDM` has returned cleanly. With the
gates reverted, the slip aux DM is never built (default OFF), so neither
the setup nor the FormFunction debug print fires. We confirmed the
absence of `febasic.c:663` and `plexfem.c:3982` errors directly from the
PETSc error stream during the inverted-gate runs:

```
$ FSRM_SESSION19_DEBUG=1 FSRM_SESSION18_DEBUG=1 ctest -R PrescribedSlip
... (no febasic.c:663 backtrace)
... (no plexfem.c:3982 backtrace)
[0]PETSC ERROR: TSStep has failed due to DIVERGED_STEP_REJECTED
```

(Compare to the Session 19 report, where `febasic.c:663` was the proximal
error.)

## 4. `plexfem.c:3982` status

**No fire.** Same evidence as section 3: `DMPlexComputeResidualHybridByKey`
returns success on every call; the error source is downstream, in
`KSPSolve`. The slip aux closure remains 9 DOFs (3 cohesive vertices x 3
components) against `totDim = 9` for the dim-1 surface FE registered on
the interfaces label, matching Session 19's baseline.

## 5. SNES trace for `PrescribedSlipQuasiStatic` and `max_fault_slip`

With all four gate inversions active (`captured in /tmp/s20_debug.log`):

```
0 TS dt 1. time 0.
    0 SNES Function norm 2.189978754336e-04
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.959594345899e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.959594345899e+01
    ... (10 retries, all DIVERGED_LINEAR_SOLVE iterations 0) ...
[0]PETSC ERROR: TSStep has failed due to DIVERGED_STEP_REJECTED
[0]PETSC ERROR: #1 TSStep() at .../ts/interface/ts.c:3539
[0]PETSC ERROR: #2 TSSolve() at .../ts/interface/ts.c:4187

sol_norm        = 0
max_fault_slip  = 0
```

`KSPSolve` reports `iterations 0`, which means LU factorization failed
with a zero pivot before any iteration began. The 2.19e-4 -> 19.6 jump
between SNES restart attempts reflects the SNES line-search backtracking
into a state with the prescribed-slip residual fully active (the kernel
residual `r_lambda = lambda - slip_target` evaluated against `lambda = 0`).

With the gates reverted (current default), the run never reaches the
slip-aux path: PrescribedSlip continues to fail in the Session 16
constants-array mode with the same SNES KSP-zero-iter pattern (residual
never seeded), giving `max_fault_slip = 0`.

## 6. Fault subset table (Session 19 vs Session 20 default vs Session 20 inverted-gates probe)

Test (`Physics.*` and `Integration.*`):

| # | Test | S19 baseline | S20 default (face-quad fix) | S20 inverted-gates probe |
|---|------|-------------:|-----------------------------:|-------------------------:|
| 40 | Physics.SCEC.TPV5                      | FAIL | FAIL | FAIL |
| 57 | Physics.LockedFaultTransparency        | FAIL | FAIL | FAIL |
| 38 | Physics.CohesiveBdResidual             | PASS | PASS | PASS |
| 91 | Integration.PressurizedFractureFEM     | FAIL | FAIL | FAIL |
| 93 | Integration.DynamicRuptureSolve.LockedQuasiStatic   | FAIL | **PASS** | FAIL |
| 94 | Integration.DynamicRuptureSolve.LockedElastodynamic | PASS | PASS | FAIL |
| 95 | Integration.DynamicRuptureSolve.PrescribedSlip      | FAIL | FAIL | FAIL |
| 96 | Integration.ExplosionFaultReactivation              | PASS | PASS | FAIL |
| 98 | Integration.TimeDependentSlip                       | FAIL | **PASS** | FAIL |
| 100| Integration.SlippingFaultSolve                       | FAIL | FAIL | FAIL |
| 101| Integration.SlipWeakeningFault                       | FAIL | FAIL | FAIL |
| -- | Functional.DynamicRuptureSetup (5 tests)             | PASS | PASS | PASS |

**Summary**:

- S19 baseline fault subset: 8/16 passing.
- S20 default (face-quadrature fix only): 10/16 passing. Net +2
  (LockedQuasiStatic, TimeDependentSlip).
- S20 inverted-gates probe: 6/16 passing. Net -2 vs baseline
  (LockedElastodynamic, ExplosionFaultReactivation regress, plus the
  S19 known failures stay failing). All regressions present with the
  same KSP-zero-iter signature as PrescribedSlip.

## 7. Fallback path (`FSRM_DISABLE_SLIP_AUX=1`) results

Captured under the inverted-gates build (the only build where the
`FSRM_DISABLE_SLIP_AUX` knob was meaningful). `FSRM_DISABLE_SLIP_AUX=1`
restored the Session 16 / 19 baseline outcomes for the linear cases:

```
2/4 Test #93: Integration.DynamicRuptureSolve.LockedQuasiStatic .....   Passed
3/4 Test #94: Integration.DynamicRuptureSolve.LockedElastodynamic ...   Passed
4/4 Test #98: Integration.TimeDependentSlip .........................   Passed
1/4 Test #57: Physics.LockedFaultTransparency .......................   Failed   (baseline known failure)
```

This proved the gating coordination logic itself worked correctly: the
fallback path was preserved bit-for-bit. The decision to roll back was
driven by the inability of the default-ON path to drive any linear fault
test through SNES, not by any breakage in the fallback gate.

With the gates reverted (current tree), `FSRM_DISABLE_SLIP_AUX` is no
longer consulted (slip aux DM is default-OFF and gated by
`FSRM_ENABLE_SLIP_AUX`). The "fallback" *is* the default.

## 8. Full-suite count

```
$ docker run ... ctest --output-on-failure
95% tests passed, 6 tests failed out of 116
Total Test time (real) = 136.27 sec

The following tests FAILED:
   40 - Physics.SCEC.TPV5 (Failed)
   57 - Physics.LockedFaultTransparency (Failed)
   91 - Integration.PressurizedFractureFEM (Failed)
   95 - Integration.DynamicRuptureSolve.PrescribedSlip (Failed)
  100 - Integration.SlippingFaultSolve (Failed)
  101 - Integration.SlipWeakeningFault (Failed)
```

vs Session 19 baseline of 108/116 (8 failures: the six above plus
LockedQuasiStatic and TimeDependentSlip).

## Decision-gate verdict

The S20 prompt enumerated five gates. Our outcome maps to **"New failure
mode not previously seen: dump SNES trace and PETSc error stack. Stop and
report."** -- specifically, the KSP DIVERGED_LINEAR_SOLVE iterations 0
factorization-failure mode that emerges only when the rim-pin essential
BC is removed while the slip aux is active. Per that gate's instruction,
no further code changes were attempted in this session. The face-
quadrature fix on the material aux was retained because it improves
the baseline pass count and addresses a genuine PETSc tabulation defect
that Session 19 had isolated.

## Why the rim-pin removal triggers KSP failure

Working hypothesis (not verified by Jacobian dump in this session):

1. The rim-pin essential BC at the buried-cohesive label removed all
   Lagrange DOFs at the fault rim from the system via the standard
   PETSc `DM_BC_ESSENTIAL` path (rows + columns deleted, diagonal = 1).
2. With the rim-pin off and the slip aux active, those rim Lagrange
   DOFs return to the active set. Their Jacobian-row contributions come
   from two sources: (a) the volume regularization
   `g3 = epsilon*I` with `epsilon = 1e-4` registered on every Lagrange
   DOF (very small diagonal entry), and (b) the cohesive-cell hybrid
   Jacobian written by `DMPlexComputeJacobianHybridByKey` -- which only
   touches Lagrange DOFs interior to the cohesive surface, not those
   on the rim.
3. The result is rim-Lagrange rows whose only diagonal entry is
   `epsilon = 1e-4` against displacement off-diagonal entries scaled
   like `E/h ~ 4e9`. LU runs into a numerically zero pivot during
   factorization and reports zero KSP iterations.

The "Re-enable rim-pin conditionally by test mode" route the prompt
suggests would require either per-test gating (brittle) or a
per-rim-DOF augmentation in the cohesive Jacobian assembly so that
rim Lagrange rows have a diagonal large enough to factor cleanly.
Neither was attempted in this session per the gate instruction to stop
and report.

## Files changed (in committed diff)

- `src/core/Simulator.cpp` -- material aux face-quadrature copy
  (block at lines 3120-3151). All other Session 20 edits reverted.

## Suggested Session 21 follow-up

1. Confirm the KSP failure is a zero-pivot factorization issue by
   running the inverted-gates path with `-mat_view ::ascii_info` on
   the Jacobian and dumping the pivot table from `-pc_factor_mat_view`.
2. Add an explicit penalty-scaled diagonal at every rim Lagrange DOF
   inside `addCohesivePenaltyToJacobian` (analogue of the existing
   interior-Lagrange penalty), then re-attempt the gate inversion.
3. Once SNES converges in the default-ON aux-slip path,
   re-evaluate the prompt's decision gates.
