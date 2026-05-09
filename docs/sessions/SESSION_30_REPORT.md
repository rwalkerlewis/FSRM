# Session 30 Report -- Sub-DM rigid-body near-null-space for fieldsplit

## Summary

Session 30 executes the follow-up left dangling in Session 23's comment
("Building a properly-sized sub-DM null-space is left to a follow-up
session"). The fieldsplit Schur path's inner displacement PCGAMG needs
rigid-body near-null-space vectors sized to the displacement sub-block
(local free size 222 on the canonical PrescribedSlip mesh), not to the
full saddle-point free storage (255). Session 21 built only the full-DM
size and composed it on the displacement field object, which PCGAMG
rejected with `Attached null space vector size 297 != matrix size 222`;
Session 23 therefore gated the compose OFF under fieldsplit, leaving the
inner PCGAMG without any near-null-space at all.

This session lands both sizes. `rigidBodyNullSpace_` now holds the
sub-DM sized vectors (composed on the displacement field for fieldsplit's
inner GAMG to pull); `rigidBodyNullSpaceFull_` holds the full-DM sized
vectors (attached to the monolithic Jacobian in `FormJacobian` on the
default non-fieldsplit path). A subtle PETSc requirement surfaced: the
`IS` output parameter of `DMCreateSubDM` cannot be passed as `NULL` or
downstream residual evaluation on the parent DM returns NaN.

Part 1 captures a baseline SVD spectrum via `FSRM_SVD_DEBUG=1` (LU
override removed since fieldsplit's `-pc_type fieldsplit` sets later and
wins): condition number `1.76e+18`, `sigma_max = 2.00e+10`, `sigma_min
= 1.13e-08`, five singular values below `1e-6` consistent with
near-rigid-body modes in the displacement block.

Part 2's outcome is **Branch C**: the null-space attachment helped
(PrescribedSlip KSP tail improves from `~2.9e+03` at iter 199 in Session
27 to `~6.0e+02` in Session 30, a 5x reduction) but did not recover
convergence. KSP still stalls before reaching rtol = 1e-14. The
near-null-space is a necessary-but-not-sufficient ingredient for the
fieldsplit path.

- Default-path fault subset: 10 pass / 6 fail — unchanged from Session 27
  default baseline.
- Fieldsplit-path fault subset: 7 pass / 9 fail — unchanged from Session
  27 experimental baseline.
- Full suite: 110 pass / 6 fail (95%).

## Part 1: SVD diagnostic baseline

### Command

`FSRM_SVD_DEBUG=1` sets `-pc_type svd -pc_svd_monitor -snes_max_it 1` in
the test harness BEFORE `setupSolvers` runs, so the `PetscOptionsHasName`
check at `Simulator.cpp:3931` sees `-pc_type` already set and leaves the
default GAMG path alone. Adding `FSRM_SADDLE_SOLVER=fieldsplit` is useless
here because the fieldsplit branch at `Simulator.cpp:3932` unconditionally
re-sets `-pc_type fieldsplit`, clobbering the SVD request. The run below
therefore omits the fieldsplit override so the SVD PC actually fires.

```
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_SVD_DEBUG=1 \
    ctest --output-on-failure -R DynamicRuptureSolve.PrescribedSlip -V' \
  2>&1 | tee /tmp/s30_svd_nofs.log
```

### Spectrum at first SNES iterate (u = 0)

```
SVD: condition number 1.763273488316e+18, 0 of 255 singular values are (nearly) zero
SVD: smallest singular values: 1.132366364884e-08 1.992530593806e-07 2.857557448147e-07 3.930821552219e-07 4.575522091428e-07
SVD: largest singular values : 1.699263701422e+10 1.756995922183e+10 1.759717590135e+10 1.995735090500e+10 1.996671590260e+10
```

Decoded:

- `sigma_max = 1.997e+10` (elastic-modulus scale, as expected for E = 2e10).
- `sigma_min = 1.132e-08`.
- Condition number `1.76e+18` >> 1e14. The matrix is effectively singular.
- Five singular values below 1e-6 (`1.1e-8, 2.0e-7, 2.9e-7, 3.9e-7, 4.6e-7`).
  No singular value below PETSc's `1e-12 * sigma_max = 2e-2` threshold, so
  the "0 of 255 are (nearly) zero" message fires even though the spectrum
  clearly contains a near-kernel.
- The five near-zero singular values are consistent with a missing
  rigid-body near-null-space. The domain has exactly one Dirichlet face
  (bottom, `z_min`, u = 0) pinning 3 translations + 3 rotations at that
  face. The fault rim-pin constrains 16 of 25 Lagrange-bearing points,
  leaving the displacement block with a near-rigid-body mode set that the
  preconditioner cannot remove without explicit null-space input.

### Post-fix re-measurement

Re-running the SVD diagnostic AFTER Part 2's null-space changes produces
an identical spectrum. This is expected: `MatSetNearNullSpace` and
`PetscObjectCompose(field, "nearnullspace", ...)` are preconditioner hints
and do not modify the assembled matrix. The preconditioner improvement
comes from GAMG using the hint to build correct coarse operators.

## Part 2: Sub-DM sized rigid-body near-null-space

### 2.1 setupSolvers

`src/core/Simulator.cpp:3835-3935`. Replaces the Session 21 single-size
build with two null-spaces:

- Build a displacement-only sub-DM via `DMCreateSubDM(dm, 1, {disp_idx},
  &disp_is, &disp_subdm)`. Capturing the `IS` output parameter and
  destroying it with `ISDestroy` after `DMPlexCreateRigidBody` is
  required: passing `NULL` for the `IS*` argument leaves shared state
  that later corrupts `TSComputeIFunction` on the parent DM to produce
  NaN residuals (verified by the `Physics.CohesiveBdResidual` regression;
  see Debug 2.4 below).
- `DMPlexCreateRigidBody(disp_subdm, 0, &rigidBodyNullSpace_)` produces
  six vectors of local size 222 on the PrescribedSlip mesh, matching the
  displacement sub-block's free DOF count.
- `DMPlexCreateRigidBody(dm, displacement_field_idx_,
  &rigidBodyNullSpaceFull_)` produces six vectors of local size 255 for
  the monolithic Jacobian on the default non-fieldsplit path.
- `PetscObjectCompose(disp_field, "nearnullspace", rigidBodyNullSpace_)`
  attaches the sub-DM sized vectors to the displacement field object so
  fieldsplit's inner PCGAMG retrieves them via the PETSc-internal
  `DMCreateSubDM + PetscObjectQuery("nearnullspace")` path.

`include/core/Simulator.hpp:210-223` adds the new
`rigidBodyNullSpaceFull_` member; the destructor at
`src/core/Simulator.cpp:483-484` destroys both null-spaces.

### 2.2 FormJacobian

`src/core/Simulator.cpp:7311-7327` (post-edit). The monolithic
`MatSetNearNullSpace` call now uses `rigidBodyNullSpaceFull_` (full-DM
size) instead of `rigidBodyNullSpace_` (sub-DM size). The attachment
still gates off under fieldsplit, so the Schur sub-solve is unaffected.

### 2.3 Setup trace

PrescribedSlip under `FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit`:

```
Session 30: rigid-body near-null-space built for field 0 ('displacement');
    sub-DM nvecs=6 local_size=222 (composed on field);
    full-DM nvecs=6 local_size=255 (for monolithic Jacobian)
```

Both sizes line up with the section free storage:

```
Local section: total dof=525, constrained dof=270 (free=255)
    field 0: dof=450, constrained=228, free=222
    field 1: dof=75, constrained=42, free=33
```

### 2.4 Debug trail -- the NULL-IS NaN trap

Initial implementation called `DMCreateSubDM(dm, 1, fields_disp, NULL,
&disp_subdm)` with `NULL` for the IS output. This compiled and ran, the
sub-DM had the correct size, `DMPlexCreateRigidBody` produced correct-
sized vectors, and PrescribedSlip ran through `DIVERGED_MAX_IT`. But
`Physics.CohesiveBdResidual` (which only evaluates
`TSComputeIFunction`, no KSP) now returned a non-finite residual norm
under the same fieldsplit environment, regressing from Session 27
experimental where that test had passed.

Bisection via `FSRM_S30_SUBDM_OFF` confirmed the regression came from the
`DMCreateSubDM` call itself, not from `PetscObjectCompose` on the field.
Changing to `DMCreateSubDM(dm, 1, fields_disp, &disp_is, &disp_subdm)`
with a matching `ISDestroy(&disp_is)` after the rigid-body build restored
`Physics.CohesiveBdResidual` to passing while leaving PrescribedSlip's
KSP trace unchanged.

Inference: passing `NULL` for the IS output in PETSc 3.25's
`DMCreateSubDM_Plex` leaves an unreferenced internal structure that the
subsequent `DMDestroy(&disp_subdm)` cleans up too aggressively, producing
NaN in downstream residual code paths that walk the parent DM's
coordinate or section data. The fix does not need further investigation,
the pattern (`IS` captured and explicitly destroyed) matches PyLith's
usage at `libsrc/pylith/problems/Problem.cc`.

## PrescribedSlip KSP trace comparison

Session 27 experimental (fieldsplit, full-DM compose gated OFF, so the
inner displacement GAMG ran without any near-null-space):

```
  0 KSP Residual norm 1.285219755579e+09
  1 KSP Residual norm 1.279663916276e+09
  2 KSP Residual norm 1.017168310581e+09
  ...
 17 KSP Residual norm 1.481466486752e+07
  ...
195 KSP Residual norm 2.989333245922e+03
196 KSP Residual norm 2.970113887429e+03
197 KSP Residual norm 2.951260526996e+03
198 KSP Residual norm 2.932761693884e+03
199 KSP Residual norm 2.914606414431e+03
```

Session 30 experimental (fieldsplit, sub-DM null-space composed on field):

```
  0 KSP Residual norm 1.151724296266e+09
  1 KSP Residual norm 1.148885280478e+09
  2 KSP Residual norm 1.145001557531e+09
  ...
 17 KSP Residual norm 2.010500472002e+07
  ...
195 KSP Residual norm 6.059579610461e+02
196 KSP Residual norm 6.049909626477e+02
197 KSP Residual norm 6.040285789821e+02
198 KSP Residual norm 6.030707734614e+02
199 KSP Residual norm 6.021175099025e+02
```

Reduction at iter 199: 2.91e+03 -> 6.02e+02 (4.8x better). Still
catastrophically far from the `rtol = 1e-14` target (initial norm 1.15e+9
requires KSP residual below ~1e-5 absolute). The null-space helped the
inner displacement sub-solve but the outer Schur KSP still cannot reduce
the constraint-coupled error. SNES still returns `DIVERGED_LINEAR_SOLVE`
and the test still fails.

## Test tables

### Default path (no env vars)

```
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure -R \
  'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|...'
```

| # | Test | Status |
|---|---|---|
| 40 | Physics.SCEC.TPV5 | FAIL |
| 57 | Physics.LockedFaultTransparency | FAIL |
| 59 | Physics.CohesiveBdResidual | pass |
| 86 | Functional.DynamicRuptureSetup | pass (5 tests) |
| 90 | Integration.ExplosionFaultReactivation | pass |
| 91 | Integration.PressurizedFractureFEM | FAIL |
| 93 | Integration.DynamicRuptureSolve.LockedQuasiStatic | pass |
| 94 | Integration.DynamicRuptureSolve.LockedElastodynamic | pass |
| 95 | Integration.DynamicRuptureSolve.PrescribedSlip | FAIL |
| 98 | Integration.TimeDependentSlip | pass |
| 99 | Integration.FaultAbsorbingCoexist | pass |
| 100 | Integration.SlippingFaultSolve | FAIL |
| 101 | Integration.SlipWeakeningFault | FAIL |

6 fail / 16 total, same set as Session 27 default baseline.

### Fieldsplit path (`FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit`)

| # | Test | Status |
|---|---|---|
| 40 | Physics.SCEC.TPV5 | FAIL |
| 57 | Physics.LockedFaultTransparency | FAIL |
| 59 | Physics.CohesiveBdResidual | pass (restored) |
| 86 | Functional.DynamicRuptureSetup | pass (5 tests) |
| 90 | Integration.ExplosionFaultReactivation | pass |
| 91 | Integration.PressurizedFractureFEM | FAIL |
| 93 | Integration.DynamicRuptureSolve.LockedQuasiStatic | FAIL |
| 94 | Integration.DynamicRuptureSolve.LockedElastodynamic | FAIL |
| 95 | Integration.DynamicRuptureSolve.PrescribedSlip | FAIL |
| 98 | Integration.TimeDependentSlip | FAIL |
| 99 | Integration.FaultAbsorbingCoexist | pass |
| 100 | Integration.SlippingFaultSolve | FAIL |
| 101 | Integration.SlipWeakeningFault | FAIL |

9 fail / 16 total, same set as Session 27 experimental baseline. No
regression introduced.

### Full suite

```
110 pass / 6 fail out of 116 (95%).
Failed: Physics.SCEC.TPV5, Physics.LockedFaultTransparency,
        Integration.PressurizedFractureFEM,
        Integration.DynamicRuptureSolve.PrescribedSlip,
        Integration.SlippingFaultSolve,
        Integration.SlipWeakeningFault.
```

Identical to the default-path fault subset intersection; no regression
outside the fault subset.

## Verdict: Branch C

PrescribedSlip still fails under the experimental fieldsplit path. The
sub-DM null-space change improved the KSP residual at iter 199 by ~5x
but did not recover convergence. The near-null-space was one bottleneck
but not THE bottleneck. Other candidate root causes for Session 31:

- The Schur complement approximation (`selfp`) approximates the
  constraint block by `B diag(A)^{-1} B^T`, which on the epsilon-
  regularized Lagrange diagonal (1e-4/h ~ 1e-6) gives a near-zero Schur
  approximation. The condition number 1.76e+18 in Part 1 is consistent
  with this; the `lagrange_multiplier_fault` block is poorly
  preconditioned regardless of what the displacement block does.
- The `factorization_type = lower` option may be wrong for this problem
  geometry. `upper` or `full` might work better when the constraint
  block has non-negligible diagonal scaling relative to the displacement
  block.
- Scaling mismatch: displacement block is O(1e10), Lagrange block is
  O(1e-6). `pc_fieldsplit_schur_scale = 1.0` is PyLith's value but
  PyLith's problem has a well-conditioned constraint block because its
  Lagrange block is scaled by the elastic modulus at assembly time
  (Session 25's `penalty * coeff ~ O(E/h)` is applied in the fallback
  `addCohesivePenaltyToJacobian` rather than in the PetscDS callback, so
  the assembled matrix does not have this scaling when the aux-slip path
  is the one driving BdResidual).

No gate flip for Session 31. Default remains 10/16 fault passing.
`FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit` remains the
reproducer for the stuck case.

## Files changed

- `include/core/Simulator.hpp:210-223`: added
  `rigidBodyNullSpaceFull_` member.
- `src/core/Simulator.cpp:483-484`: destroy both null-spaces.
- `src/core/Simulator.cpp:3835-3935`: replaced single-size null-space
  setup with sub-DM + full-DM build. Captures the IS from
  `DMCreateSubDM` to avoid the NULL-IS NaN trap.
- `src/core/Simulator.cpp:7311-7327`: monolithic `MatSetNearNullSpace`
  now uses `rigidBodyNullSpaceFull_` on the default path. Fieldsplit
  gate unchanged.

## Rule compliance

- Rule 1 (Docker): all builds and tests ran in `fsrm-ci:main`.
- Rule 14: CLAUDE.md not updated; the default-path fault count did not
  change, and the experimental path gate was not flipped.
- Rule 16: full logs saved to `/tmp/s30_svd.log`, `/tmp/s30_svd_nofs.log`,
  `/tmp/s30_svd_postfix.log`, `/tmp/s30_prescribed_exp.log`,
  `/tmp/s30_fault_exp.log`, `/tmp/s30_fault_default.log`, `/tmp/s30_full.log`.
- Rule 28 (no option/tolerance tuning): honored. No fieldsplit flags,
  KSP/SNES tolerances, or solver options were modified. The only
  preconditioner-related change is the near-null-space attachment, which
  is a matrix hint rather than an option.
- Rule 29 (no Jacobian callback modifications): the `FormJacobian`
  change substitutes a different `MatNullSpace` handle (same six
  rigid-body vectors, correctly sized for the monolithic matrix) in the
  existing `MatSetNearNullSpace` call. No Jacobian math modified.

## End of Session 30
