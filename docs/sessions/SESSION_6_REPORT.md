# Session 6 Report

Branch: `local_fix`.
Scope: remove the manual `addCohesivePenaltyToJacobian` augmentation so the
hybrid driver `DMPlexComputeJacobianHybridByKey` is the sole source of
displacement-Lagrange coupling on cohesive cells, and verify whether the
hybrid path produces a solvable Jacobian on `PrescribedSlipQuasiStatic`.

**Decision gate outcome: STOP. Some pass, some still
DIVERGED_LINEAR_SOLVE iterations 0. Per the gate this triggers the
"matrix dump and stop" branch, not a CLAUDE.md update or
`PYLITH_COMPATIBILITY.md` Section B retirement. The Session 5 hypothesis
("KSP fails because the manual penalty double-writes the
displacement-Lagrange entries") is **falsified** by the matrix dump in
Section 4 below: the lambda-u block is correct after the call is
removed, but the displacement-displacement block is empty. The actual
root cause is upstream in the volume Jacobian assembly path, not in the
hybrid kernels.**

---

## 1. Code changes landed

`src/core/Simulator.cpp` (FormJacobian):

- Removed both call sites that invoked `addCohesivePenaltyToJacobian` on
  `P` and `J` (formerly `Simulator.cpp:5816` and `:5820`). The block
  is replaced with an explanatory comment; no other logic touched.
- Added an explicit `MatAssemblyBegin/End(MAT_FINAL_ASSEMBLY)` after
  `DMPlexComputeJacobianHybridByKey` on `P` (and on `J` when the two
  matrices are distinct). The hybrid driver writes via
  `DMPlexMatSetClosure` but does not bookend with assembly; the previous
  `addCohesivePenaltyToJacobian` body did the assembly as a side effect.
  Without this fix, the next consumer (LU factor) trips
  `PETSC_ERR_ARG_WRONGSTATE` (73) on an unassembled matrix. PETSc
  `src/dm/impls/plex/tests/ex5.c:1228-1231` does the same MatAssembly
  bookend after the hybrid Jacobian call.
- Also added a matching `MatSetOption(MAT_NEW_NONZERO_ALLOCATION_ERR,
  PETSC_FALSE)` on `J` when `J != P`. The previous code only set this
  on `P`.

`src/core/Simulator.cpp` (`addCohesivePenaltyToJacobian`):

- Function body untouched, but the routine now emits a one-time stderr
  warning on entry. No production code path calls it; the warning is a
  trip wire so any test that still routes through it will surface dead
  code that needs removing in Session 7.

`include/core/Simulator.hpp`:

- Marked `addCohesivePenaltyToJacobian(...)` `[[deprecated]]` with a
  message pointing at the Session 5 `DIVERGED_LINEAR_SOLVE iterations 0`
  failure mode that the manual augmentation produced.

`tests/integration/test_dynamic_rupture_solve.cpp`:

- Added an opt-in `-ksp_view_mat` Petsc option gated on the
  `FSRM_DUMP_JAC` environment variable, used to capture the diagnostic
  Jacobian dump in Section 4. The instrumentation is dormant unless the
  env var is set, so it does not affect regular CI runs.

No changes to:

- `FaultMeshManager::splitMeshAlongFault`
- `CohesiveFaultKernel::registerWithDS` / `registerHybridWithDS`
- `CohesiveFaultKernel::g0_hybrid_*` callbacks
- The `(dim-1)` Lagrange FE registration landed in Session 5
- DS / boundary ordering in `setupFields`

Per the task, Rule 6 was not invoked.

## 2. Verification of the hybrid g0 callbacks against our F sign

Session 5's residual sign convention (`src/physics/CohesiveFaultKernel.cpp:729-883`):

| Callback                | Residual                                        |
|-------------------------|--------------------------------------------------|
| `f0_hybrid_u_neg[c]`    | `-u[uOff[1]+c]`  =  `-lambda`                    |
| `f0_hybrid_u_pos[c]`    | `+u[uOff[1]+c]`  =  `+lambda`                    |
| `f0_hybrid_lambda[c]`   | `u[Nc+c] - u[c] - delta`  =  `u_pos - u_neg - delta` (locked: `delta = 0`) |

Required Jacobian blocks (`J[i,j] = dR[i]/du[j]`):

| Block                          | Expected        | g0 callback                     | Code (line)                                              |
|--------------------------------|-----------------|----------------------------------|----------------------------------------------------------|
| `d(R_u_neg)/d(lambda)`         | `-I`            | `g0_hybrid_u_lambda_neg`         | `g0[c*dim+c]=-1.0` (`CohesiveFaultKernel.cpp:885-908`)   |
| `d(R_u_pos)/d(lambda)`         | `+I`            | `g0_hybrid_u_lambda_pos`         | `g0[c*dim+c]=+1.0` (`CohesiveFaultKernel.cpp:910-933`)   |
| `d(R_lambda)/d(u_neg)`         | `-I`            | `g0_hybrid_lambda_u` (1st block) | `g0[c*Nc+c]=-1.0` (`CohesiveFaultKernel.cpp:962-964`)    |
| `d(R_lambda)/d(u_pos)`         | `+I`            | `g0_hybrid_lambda_u` (2nd block) | `g0[Nc*Nc+c*Nc+c]=+1.0` (`CohesiveFaultKernel.cpp:962-964`) |
| `d(R_lambda)/d(lambda)`        | `0` (locked / prescribed) | `g0_hybrid_lambda_lambda` | early return at `mode<0.5 || mode>1.5` (`CohesiveFaultKernel.cpp:1076-1079`) |

PETSc's hybrid integration layout was cross-checked at
`src/dm/dt/fe/impls/basic/febasic.c:1171-1187` and
`src/dm/dt/fe/interface/fe.c:2712-2770`: for a cohesive test field
against a non-cohesive trial field on a cohesive cell, the integrator
calls `PetscFEUpdateElementMat_Hybrid_Internal` twice with `(s,t) =
(0,0)` then `(1,1)`; the first call reads `g0[0..NcI*NcJ-1]` (trial side
0 = negative) and the second reads `g0[NcI*NcJ..2*NcI*NcJ-1]` (trial
side 1 = positive). DMPlex sets cone[0] as the negative side and
cone[1] as the positive side. Our `g0_hybrid_lambda_u` layout therefore
matches PETSc's expectation, and the assembled lambda-u block in the
matrix dump (Section 4) confirms this empirically.

PETSc ex5's `g0_bd_lu` writes `[neg=-1, pos=+1]` against an
`f0_bd_l = delta + u^- - u^+` whose direct linearization is `[neg=+1,
pos=-1]`. Either ex5 carries a sign convention (or an internal sign flip
in the integrator) we cannot reproduce from source alone, or ex5 itself
has a sign ambiguity that its CI tests do not surface. We did not copy
ex5's signs blindly. Our F is consistent with our J, the residual
norm matches the prescribed-slip scale (Section 3), and the dump in
Section 4 shows the expected lambda-u stamp. The remaining failure is
elsewhere.

## 3. Exact `PrescribedSlipQuasiStatic` SNES trace

Captured at `/tmp/s6_pres2.log` from
`./tests/run_integration_tests --gtest_filter=DynamicRuptureSolveTest.PrescribedSlipQuasiStatic`
in Docker, full output, no truncation:

```
Step 0, Time = 0.
HDF5 output written at step 0 to output/solution.h5
0 TS dt 1. time 0.
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.000000015625e+00
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.000000015625e+00
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.251587298439e-02
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
[gtest assertions]
Expected: (sol_norm) > (0.0), actual: 0 vs 0
Expected: (max_fault_slip) > (5.0e-4), actual: 0 vs 0.0005
[FAILED] DynamicRuptureSolveTest.PrescribedSlipQuasiStatic (534 ms)
```

Specifically:

- `DIVERGED_LINEAR_SOLVE iterations 0` is **still present**, and
  identical to the Session 5 trace.
- Initial SNES function norm: **`1.767766952966e-04`**, identical to
  Session 5. Removing the manual penalty did not change the residual
  (expected -- the residual was never the question).
- TS adaptive cutbacks produce the same cycle (`1.0` → `6.25e-2`) seen
  in Session 5.
- Final `max_fault_slip` (test assertion): **`0`**, identical to
  Session 5. SNES never produces an update.

The Session 5 prediction "remove the double-write and the linear solve
will succeed" did not hold. The matrix dump in Section 4 explains why.

## 4. Jacobian dump (`-ksp_view_mat`)

Captured via `FSRM_DUMP_JAC=ascii:/workspace/build/s6_jac.txt` on the
2x2x2 prescribed-slip configuration. Stored as `/tmp/s6_jac.txt`
(1.1 MB, 3520 lines). Representative rows:

Row 0 (a displacement DOF, `u_x` at vertex 0):

```
row 0: (0, 0.) (1, 0.) (2, 0.) (120, -0.0104167) (121, 0.) (122, 0.)
       (131, 0.) (132, 0.) (133, 0.) (294, 0.0104167) (295, 0.) (296, 0.)
```

Row 117 (a Lagrange DOF, `lambda_x` at a fault vertex):

```
row 117: (6, -0.0104167) (7, 0.) (8, 0.) (27, -0.0104167) (28, 0.) (29, 0.)
         (57, -0.0104167) (58, 0.) (59, 0.) (75, -0.0104167) (76, 0.) (77, 0.)
         (81, -0.0104167) (82, 0.) (83, 0.) (90, -0.0104167) (91, 0.) (92, 0.)
         (104..173, all 0.)  (291..305, all 0.)
```

Row 316 (a Lagrange DOF on the positive-side trial space):

```
row 316: (15, 0.) (16, 0.0104167) (17, 0.) (21, 0.) (22, 0.0104167) (23, 0.)
         ...
```

Structural observations:

1. **Lambda-u block is correct.** Lagrange rows have `-coeff` stamps
   against negative-side displacement DOFs and `+coeff` stamps against
   positive-side displacement DOFs (`coeff = face_area / n_pairs =
   1/96 ~= 0.01042` for the 2x2x2 mesh). This matches the expected
   `g0_hybrid_lambda_u` layout `[d(R_lambda)/d(u_neg) = -I,
   d(R_lambda)/d(u_pos) = +I]`. Removing the manual penalty did not
   corrupt these entries; on the contrary the resulting block is exactly
   what we wanted.
2. **U-lambda block is correct.** Row 0 has `-0.01` at column 120
   (lambda DOF on the same fault vertex, negative side) and `+0.01` at
   column 294 (lambda DOF on the partner positive side). This matches
   `g0_hybrid_u_lambda_neg = -I` and `g0_hybrid_u_lambda_pos = +I`.
3. **Lambda-lambda block is zero.** All rows in `[117, 122]` and
   subsequent Lagrange ranges show `(117..122, 0.)`, etc., consistent
   with the locked / prescribed branch of `g0_hybrid_lambda_lambda`
   that returns zero. Correct.
4. **Displacement-displacement block is empty.** Every displacement row
   inspected (rows 0, 1, 2, 37, 77, 97, 197) shows zero values on the
   diagonal and zero values on every off-diagonal that should carry
   elastic stiffness. Example: row 97 has 21 column entries
   `(96..162)`, all displacement-side, **every value `0.`**. This is
   the entire 3D linear elasticity stiffness for the
   `DMPlexTSComputeIJacobianFEM` pass, and it is missing.

The `MatAssembly` bookend in Section 1 ensures the matrix is in
assembled state, so the values shown are the genuine assembled values
-- not a viewing artifact. The displacement diagonal is zero, the
elastic coupling is zero, and `LU` of this matrix is singular at the
first row, hence `DIVERGED_LINEAR_SOLVE iterations 0`.

The Session 5 manual penalty was **masking** this gap: its
`penalty * coeff` augmentation on
`(neg_disp, neg_disp), (pos_disp, pos_disp)` and the
`-penalty * coeff` cross terms were the only reason any of the locked /
prescribed test variants ever produced a nonsingular displacement
diagonal. Even with the masking the linear solve still diverged in
Session 5 because the penalty was not a substitute for true elastic
stiffness; it just kept the matrix from being trivially zero on the
diagonal.

The actual question Session 7 needs to answer is **why
`DMPlexTSComputeIJacobianFEM` is producing zero displacement entries
on this configuration**. Plausible candidates, ranked by suspicion:

1. The PetscDS volume Jacobian for elasticity (`g3_uu_3d`) is registered
   on the **default DS** (which only carries displacement after Session
   5's region restriction) but the cell loop in
   `DMPlexTSComputeIJacobianFEM` is iterating cells whose
   `DMGetCellDS` returns the cohesive-region DS (which has displacement
   + Lagrange but does not have `g3_uu_3d` registered). `setupPhysics`
   (`src/core/Simulator.cpp:2300+`) registers elastic Jacobian on the
   default DS only.
2. The Lagrange field's `(dim-1)` surface FE inadvertently changes the
   PetscFE quadrature on cohesive cells in a way that makes the
   displacement basis derivative integration return zero.
3. `DMPlexTSComputeIJacobianFEM` is being short-circuited at runtime by
   one of the new region keys (the cohesive cells live in the cohesive
   DS region but the displacement field is global, so the iteration
   plan may be confused).

Recommended Session 7 first probe: in `setupPhysics`, after the volume
Jacobian registration, walk all DSes via `DMGetNumDS` /
`DMGetRegionNumDS` and confirm `g3_uu_3d` is attached to every DS that
contains the displacement field (mirroring the Session 5 fix that
walked DSes for `PetscDSSetCohesive`).

Per the gate, this report does **not** attempt that fix. The dump and
diagnosis are the deliverable; Session 7 owns the implementation.

## 5. Status table of all thirteen fault-related tests

"Session 5" column from `docs/SESSION_5_REPORT.md` Section 4
(`/tmp/s5_fault.log`). "Session 6" column from `/tmp/s6_test.log` and
`/tmp/s6_fault.log` (full ctest run plus focused fault subset).

| Test                                                            | Session 5  | Session 6  | Notes |
|-----------------------------------------------------------------|------------|------------|-------|
| Physics.CohesiveBdResidual                                      | PASS       | PASS       | unchanged |
| Physics.LockedFaultTransparency                                 | FAIL       | FAIL       | DIVERGED_LINEAR_SOLVE, same as S5 |
| Physics.SCEC.TPV5                                               | FAIL       | FAIL       | unchanged |
| Integration.PressurizedFractureFEM                              | FAIL       | FAIL       | unchanged |
| Integration.DynamicRuptureSolve.LockedQuasiStatic               | FAIL       | FAIL       | unchanged |
| Integration.DynamicRuptureSolve.LockedElastodynamic             | FAIL       | FAIL       | unchanged |
| Integration.DynamicRuptureSolve.PrescribedSlip                  | FAIL       | FAIL       | DIVERGED_LINEAR_SOLVE iters 0; matrix dump in Section 4 |
| Integration.TimeDependentSlip                                   | FAIL       | FAIL       | unchanged |
| Integration.SlippingFaultSolve                                  | FAIL       | FAIL       | unchanged |
| Integration.SlipWeakeningFault                                  | FAIL       | FAIL       | unchanged |
| Integration.ExplosionFaultReactivation                          | PASS       | PASS       | unchanged |
| Functional.FaultAbsorbingCoexist                                | PASS       | PASS       | unchanged |
| Functional.DynamicRuptureSetup.LockedFault                      | PASS       | PASS       | unchanged |
| Functional.DynamicRuptureSetup.FaultGeometryFromConfig          | PASS       | PASS       | unchanged |
| Functional.DynamicRuptureSetup.SlippingFault                    | PASS       | PASS       | unchanged |
| Functional.DynamicRuptureSetup.MeshSplitting                    | PASS       | PASS       | unchanged |

**No fault-test regressions vs Session 5.** The same nine tests fail
with the same root cause; the same seven tests still pass.

Full-suite summary, Session 6 (`ctest -j`, `/tmp/s6_test.log`):
**13 failed of 116** (vs Session 5: 15 failed of 116). The two
differences are `Physics.AbsorbingBC` and
`Integration.HistoricNuclear.DegelenMountain`, both flagged in
CLAUDE.md as parallel HDF5 artifacts. Session 6 also picked up
`Integration.NearFieldCoupled` and `Physics.LithostaticStress` as
parallel artifacts (these were also flagged as parallel-flake in
Session 5). Net: no fault regressions, same elastic baseline; the
parallel-HDF5 artifact set rotated by the usual one or two tests.

## 6. Decision gate evaluation

Applying the Session 6 task gate block by block:

- **"LockedFaultTransparency AND LockedQuasiStatic AND PrescribedSlip
  AND LockedElastodynamic AND CohesiveBdResidual all pass AND no new
  regressions: declare core fault architecture working"** -- NOT MET.
  Only `CohesiveBdResidual` passes. The other four still fail with the
  same `DIVERGED_LINEAR_SOLVE iterations 0`.
- **"Some pass, some still DIVERGED_LINEAR_SOLVE: report and stop with
  the matrix dump from item 4"** -- **MATCHES**. This report and
  `/tmp/s6_jac.txt` are the deliverables.
- **"Regressions vs Session 5: revert and report"** -- not applicable;
  no fault-test regressions.

Per the gate, **CLAUDE.md and `PYLITH_COMPATIBILITY.md` are not
updated this session.** The hybrid Jacobian path is now provably
self-contained (no double-write), but the Jacobian itself is still
singular because of an upstream gap in volume assembly that is unrelated
to the cohesive kernels. Section B of `PYLITH_COMPATIBILITY.md`
(interior-Lagrange penalty diagonal) is still dead in the source -- the
disjoint loop was removed in Session 5 -- but cannot be formally retired
in the doc until the Section 7 follow-up confirms the hybrid Jacobian is
the sole driver of all cohesive entries with no new fallback added.

## 7. Followup work for Session 7

Strict order; do not skip:

1. Diagnose why `DMPlexTSComputeIJacobianFEM` produces no
   displacement-displacement entries on the prescribed-slip mesh. Start
   with the DS-walk hypothesis in Section 4 item 1 (register the
   elastic `g3_uu_3d` callback on every DS that contains the
   displacement field, mirroring the Session 5 PetscDSSetCohesive fix).
2. Once the displacement diagonal is non-zero, rerun
   `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic` and
   confirm `max_fault_slip > 5e-4`. The hybrid lambda-u and u-lambda
   blocks are already correct (Section 4); the remaining unknown is
   only the volume stiffness.
3. With prescribed slip working, revisit the locked elastodynamic
   regression (`LockedElastodynamic`,
   `Integration.PressurizedFractureFEM`, `Physics.SCEC.TPV5`). Each
   exercises a path the dim-1 Lagrange field has not been validated on.
4. Once all linear-constraint variants pass, the slipping-fault tests
   (`SlippingFaultSolve`, `SlipWeakeningFault`) become the next
   target. The slipping branches of `g0_hybrid_lambda_u` and
   `g0_hybrid_lambda_lambda` (`CohesiveFaultKernel.cpp:966-1041`,
   `1081-1136`) port the pre-Session-4 friction tangent math; they
   may need their own diagnostic pass once the underlying linear-solve
   issue is gone.
5. After (4), delete `addCohesivePenaltyToJacobian` entirely (the
   `[[deprecated]]` and stderr trip wire are in place this session) and
   formally retire `PYLITH_COMPATIBILITY.md` Section B.

The Session 5-era hypothesis "Session 6 sign-check on
`f0_hybrid_lambda` against PETSc ex5.c:1098" turned out to be the
**wrong** root cause; it was promoted to the Session 5 followup based
on the residual norm being on-scale, but the matrix dump in Section 4
shows the Lagrange block is in fact correct. The DIVERGED_LINEAR_SOLVE
is driven entirely by the missing displacement-displacement stiffness.

## 8. Log artifacts

- `/tmp/s6_build.log` -- clean build, only pre-existing `-Wunused-*`
  warnings.
- `/tmp/s6_build2.log` -- incremental rebuild after the MatAssembly
  fix.
- `/tmp/s6_test.log` -- full `ctest -j` run, 13/116 failures.
- `/tmp/s6_fault.log` -- focused fault-subset run, 9/16 failures
  matching the Session 5 set.
- `/tmp/s6_pres.log` -- isolated PrescribedSlip run before the
  MatAssembly fix (failure mode `ierr=73 PETSC_ERR_ARG_WRONGSTATE`,
  TSSolve aborts after one SNES print -- this is the regression that
  the MatAssembly fix resolves).
- `/tmp/s6_pres2.log` -- isolated PrescribedSlip run after the
  MatAssembly fix (full DIVERGED_LINEAR_SOLVE trace shown in
  Section 3).
- `/tmp/s6_pres_jac.log` -- PrescribedSlip run with
  `FSRM_DUMP_JAC=ascii:/workspace/build/s6_jac.txt`.
- `/tmp/s6_jac.txt` -- assembled Jacobian dump (1.1 MB, 3520 lines)
  for offline analysis. Source for the structural observations in
  Section 4.
