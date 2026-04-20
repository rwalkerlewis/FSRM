# Session 23 Report -- PyLith production fieldsplit Schur replication

## Summary

Session 23 replaced the Session 14/21/22 fieldsplit option set with PyLith's
verbatim `solver_fault_fieldsplit.cfg` options (factorization=lower, selfp,
pc_use_amat=true) and renamed the displacement FE to drop the trailing
underscore so PyLith's option keys (`-fieldsplit_displacement_pc_type`) line
up with PETSc's auto-derived field prefix.

The default GAMG path remains intact and the full suite count improved by
two over the Session 22 baseline (110/116 vs 108/116). The experimental
fieldsplit path reaches PCFieldSplit setup with both fields detected, but
the Schur sub-solve diverges. Decision: do NOT invert the gates -- keep
fieldsplit gated behind `FSRM_SADDLE_SOLVER=fieldsplit` for diagnostics.

## Code changes

### Fix 1 -- displacement FE name (`src/core/Simulator.cpp:1672-1677`)

```cpp
// Session 23: name the displacement FE "displacement" (no trailing
// underscore) to match PyLith. Fieldsplit option keys are derived from
// the FE name, so the trailing underscore would produce
// -fieldsplit_displacement__pc_type (double underscore) and miss
// PyLith's -fieldsplit_displacement_pc_type keys.
ierr = PetscObjectSetName((PetscObject)fe, "displacement"); CHKERRQ(ierr);
```

`SeismometerNetwork::initialize` (`src/domain/seismic/SeismometerNetwork.cpp:249`)
was updated to match on prefix `"displacement"` instead of `"displacement_"`
so seismometer auto-detection still finds the field.

### Fix 2 -- fieldsplit Schur options (`src/core/Simulator.cpp:3919-3956`)

The `FSRM_SADDLE_SOLVER=fieldsplit` branch now sets PyLith's exact options:

```
-pc_type fieldsplit
-pc_use_amat true
-pc_fieldsplit_type schur
-pc_fieldsplit_schur_factorization_type lower
-pc_fieldsplit_schur_precondition selfp
-pc_fieldsplit_schur_scale 1.0
-fieldsplit_displacement_ksp_type preonly
-fieldsplit_displacement_pc_type gamg          (ml if PetscHasExternalPackage)
-fieldsplit_lagrange_multiplier_fault_ksp_type preonly
-fieldsplit_lagrange_multiplier_fault_pc_type gamg
-ksp_type gmres -ksp_gmres_restart 100
```

`PetscHasExternalPackage` is called with the (const char[], PetscBool*)
signature -- the petsc-main header dropped the leading MPI_Comm argument
that some older docs still show.

### Fix 3 -- gate near-null-space attachments off the fieldsplit path

The prompt asserted that the Session 21 `MatSetNearNullSpace` in
`FormJacobian` was already gated to `FSRM_SADDLE_SOLVER=gamg`. It was not;
it ran unconditionally whenever `rigidBodyNullSpace_` existed. Added a
`fieldsplit_active` check at `src/core/Simulator.cpp:7263-7273` so the
parent matrix attachment is skipped when fieldsplit is the active solver.

The Session 14 `PetscObjectCompose(disp_field, "nearnullspace", ...)`
call in `setupSolvers` is also now gated off when
`FSRM_SADDLE_SOLVER=fieldsplit`. The reason: `DMPlexCreateRigidBody(dm,
field=0)` returns 6 vectors sized for the full local DM (297 entries on
the canonical PrescribedSlip mesh) with non-zero entries only in the
displacement-field rows. When fieldsplit pulls that composed null-space
and feeds it to PCGAMG on the displacement sub-block (free local size 222),
PETSc raises `PETSC_ERR_PLIB`:

```
[0]PETSC ERROR: PETSc has generated inconsistent data
[0]PETSC ERROR: Attached null space vector size 297 != matrix size 222
[0]PETSC ERROR: #1 PCSetData_AGG() at .../gamg/agg.c:476
```

Building a properly-sized sub-DM null-space (DMCreateSubDM on the
displacement field, then DMPlexCreateRigidBody on the sub-DM) is left
for the next session. With the compose gated off, the displacement
sub-block runs without an attached null-space.

## Experimental fieldsplit run -- PrescribedSlip

`FSRM_ENABLE_SLIP_AUX=1 FSRM_DISABLE_RIM_PIN=1 FSRM_SADDLE_SOLVER=fieldsplit`,
PETSc options `-snes_monitor -snes_converged_reason -info :pc -ksp_view`.

### Fieldsplit detection (positive)

```
[0] <pc:fieldsplit> PCFieldSplitSetDefaults(): Setting up physics based
    fieldsplit preconditioner using the embedded DM
[0] <pc:gamg> PCSetUp_GAMG(): fieldsplit_displacement_:
    level 0) N=222, n data rows=1, n data cols=1, nnz/row (ave)=27,
    block size 1, np=1
[0] <pc:gamg> PCSetUp_GAMG(): fieldsplit_lagrange_multiplier_fault_:
    level 0) N=75, n data rows=3, n data cols=3, nnz/row (ave)=32,
    block size 3, np=1
```

PETSc detected both fields by name (`displacement`, `lagrange_multiplier_fault`)
and built independent GAMG hierarchies. Schur preconditioner setup completed
(both `pc_fieldsplit_schur_precondition selfp` and `pc_fieldsplit_schur_factorization_type lower`).

### Schur sub-solves

Displacement K block:
- Two-level GAMG hierarchy (N=222 -> 8 coarse), operator complexity 1.01
- Smoothed prolongator eigenvalues 0.0906 .. 2.179 -- well-conditioned
- This sub-block is fine.

Lagrange Schur complement (selfp):
- N=75, block size 3
- `[0] <pc:jacobi> PCSetUp_Jacobi(): Zero detected in diagonal of matrix,
   using 1 at those locations` (twice -- once in selfp construction, once
   in the Lagrange GAMG smoother)
- Smoothed prolongator eigenvalues `min=3.98e-16 .. max=4.998` -- the
  Schur complement is essentially singular over the Lagrange space because
  selfp uses `inv(diag(A_disp))` which on rim-pin-disabled fields with
  no constraint at non-cohesive Lagrange DOFs collapses many rows to zero.

### KSP convergence (negative)

```
0 KSP Residual norm 4.088e+10
1 KSP Residual norm 1.452e+10
2 KSP Residual norm 3.708e+09
...
99 KSP Residual norm 1.232e+09  (first restart)
...
199 KSP Residual norm 1.063e+09 (about to restart again)
[0]PETSC ERROR: Residual norm computed by GMRES recursion formula
    1.06083e+09 is far from the computed residual norm 4.63889e+09 at
    restart, residual norm at start of cycle 1.74419e+09
```

GMRES drops three orders of magnitude in the first eight iterations then
plateaus near 1e9 for the remainder of the cycle. The recursion-vs-true-
residual divergence at restart aborts the solve with `PETSC_ERR_NOT_CONVERGED`
(propagated up through TSSolve). The plateau and the recursion drift point
at the singular Schur complement: with `pc_use_amat=true`, the residual is
recomputed against the application matrix while the preconditioner uses the
nearly-singular selfp approximation, so floating-point drift accumulates.

### SNES trace and max_fault_slip

```
0 TS dt 1. time 0.
    0 SNES Function norm 2.189978754336e-04
    [KSP fails inside the first SNES step]
TSSolve returns ierr=77 (PETSC_ERR_PLIB) inside the test wrapper, which
the test reports as PrescribedSlipQuasiStatic FAILED.
```

`max_fault_slip` is not reached because TSSolve aborts before the first
implicit Newton step completes. (The Session 22 default-path baseline
reported `max_fault_slip = 263.0` -- still the canonical figure.)

## Fault subset table

| #   | Test                                              | S22 default | S23 default | S23 experimental |
|-----|---------------------------------------------------|-------------|-------------|------------------|
|  40 | Physics.SCEC.TPV5                                 | FAIL        | FAIL        | FAIL             |
|  57 | Physics.LockedFaultTransparency                   | FAIL        | FAIL        | FAIL             |
|  59 | Physics.CohesiveBdResidual                        | PASS        | PASS        | FAIL             |
|  77 | Functional.DynamicRuptureSetup.MeshSplitting      | PASS        | PASS        | PASS             |
|  78 | Functional.DynamicRuptureSetup.LockedFault        | PASS        | PASS        | PASS             |
|  79 | Functional.DynamicRuptureSetup.FaultGeometryFromConfig | PASS   | PASS        | PASS             |
|  80 | Functional.DynamicRuptureSetup.SlippingFault      | PASS        | PASS        | PASS             |
|  81 | Functional.DynamicRuptureSetup.AbsorbingCoexist   | PASS        | PASS        | PASS             |
|  91 | Integration.PressurizedFractureFEM                | FAIL        | FAIL        | FAIL             |
|  93 | Integration.DynamicRuptureSolve.LockedQuasiStatic | FAIL (claimed) | **PASS** | FAIL             |
|  94 | Integration.DynamicRuptureSolve.LockedElastodynamic | PASS      | PASS        | FAIL             |
|  95 | Integration.DynamicRuptureSolve.PrescribedSlip    | FAIL        | FAIL        | FAIL             |
|  96 | Integration.ExplosionFaultReactivation            | PASS        | PASS        | PASS             |
|  98 | Integration.TimeDependentSlip                     | FAIL (claimed) | **PASS** | FAIL             |
| 100 | Integration.SlippingFaultSolve                    | FAIL        | FAIL        | FAIL             |
| 101 | Integration.SlipWeakeningFault                    | FAIL        | FAIL        | FAIL             |
| **Pass count** |                                        | 8/16        | **10/16**   | 6/16             |

S23 default recovered `LockedQuasiStatic` and `TimeDependentSlip` (both
listed as known failures in CLAUDE.md). The only default-path code change
in this session that could affect them is the displacement FE name
(`"displacement_"` -> `"displacement"`); the Session 21 near-null-space
attachment in `FormJacobian` is unchanged for the default path. Suspected
mechanism: PETSc's PCSetUp uses field names to decide block-size hints
and DM auto-split heuristics; the previous trailing underscore may have
fooled an internal lookup. The recovery is consistent across reruns.

S23 experimental regresses 4 tests beyond the S23 default failures
(CohesiveBdResidual, LockedQuasiStatic, LockedElastodynamic,
TimeDependentSlip). All four regressions match the same Schur-complement
divergence pattern observed for PrescribedSlip.

## Default path regression check

`docker run ... ctest -R '<fault subset>'` (no FSRM_* env vars set):

```
63% tests passed, 6 tests failed out of 16
        40 - Physics.SCEC.TPV5 (Failed)
        57 - Physics.LockedFaultTransparency (Failed)
        91 - Integration.PressurizedFractureFEM (Failed)
        95 - Integration.DynamicRuptureSolve.PrescribedSlip (Failed)
       100 - Integration.SlippingFaultSolve (Failed)
       101 - Integration.SlipWeakeningFault (Failed)
```

No regressions in the default path. The 6 remaining failures are a strict
subset of the CLAUDE.md baseline 8.

## Full suite

```
docker run ... ctest --output-on-failure
95% tests passed, 6 tests failed out of 116
        40 - Physics.SCEC.TPV5 (Failed)
        57 - Physics.LockedFaultTransparency (Failed)
        91 - Integration.PressurizedFractureFEM (Failed)
        95 - Integration.DynamicRuptureSolve.PrescribedSlip (Failed)
       100 - Integration.SlippingFaultSolve (Failed)
       101 - Integration.SlipWeakeningFault (Failed)
```

110/116 pass. CLAUDE.md baseline was 108/116. Net +2 (`LockedQuasiStatic`
and `TimeDependentSlip` now pass).

## Decision gate

- PrescribedSlip experimentally: **FAIL** (Schur sub-solve diverges).
- Default path regression: none.
- Fieldsplit assembled but failed to converge: yes.

Per the prompt's decision gate, gates are NOT inverted. `FSRM_ENABLE_SLIP_AUX`
/ `FSRM_DISABLE_RIM_PIN` / `FSRM_SADDLE_SOLVER` remain opt-in.

## Proposed next-session tuning

Two avenues to make fieldsplit Schur converge for the prescribed-slip path:

1. **Fix the selfp Schur complement singularity.** With `pc_use_amat=true`,
   selfp = `B inv(diag(A_disp)) B^T` requires non-trivial diagonal entries.
   The Session 22 `vpbjacobi` smoother for the parent GAMG path solved
   the analogous singularity by treating displacement as a 3-component
   block. For fieldsplit, the analogous fix is
   `-pc_fieldsplit_schur_precondition a11` (use the actual A11 block,
   which carries the small epsilon=1e-4 weak-volume Lagrange diagonal
   added in Section B Option B) or `user` with an explicit
   PCFieldSplitSchurPrecondition handler.

2. **Build a sub-DM-sized rigid-body null-space.** PyLith composes the
   null-space on the displacement field via PetscObjectCompose, but its
   nullSpace is built from a sub-DM, not from the parent. Implement:

   ```cpp
   IS disp_is = nullptr;
   DM disp_subdm = nullptr;
   DMCreateSubDM(dm, 1, &displacement_field_idx_, &disp_is, &disp_subdm);
   MatNullSpace subNullSpace;
   DMPlexCreateRigidBody(disp_subdm, 0, &subNullSpace);
   PetscObject disp_field = nullptr;
   DMGetField(dm, displacement_field_idx_, NULL, &disp_field);
   PetscObjectCompose(disp_field, "nearnullspace", (PetscObject)subNullSpace);
   ```

   This produces 6 vectors of size 222 (the displacement sub-block) instead
   of 297. Re-enable the compose for the fieldsplit path and retry.

Recommend trying tuning (1) first; it is a single-option change and exercises
the existing PETSc machinery. Tuning (2) is the structurally correct fix but
requires verifying that `DMCreateSubDM` on a parent DM that has cohesive
cells produces a usable sub-DM for `DMPlexCreateRigidBody`.

## Logs captured

- `/tmp/s23_build.log`, `/tmp/s23_build2.log` -- build output
- `/tmp/s23_prescribed_experimental.log` -- ctest experimental run
- `/tmp/s23_prescribed_diag.log` -- ctest experimental run with extra
  diagnostic options
- `/tmp/s23_fsrm_direct.log`, `/tmp/s23_fsrm_direct2.log` -- direct fsrm
  invocations capturing the full PCSetData_AGG and GMRES recursion errors
- `/tmp/s23_fault_experimental.log` -- fault subset under fieldsplit
- `/tmp/s23_fault_default.log` -- fault subset under default
- `/tmp/s23_test.log` -- full ctest suite (default)
