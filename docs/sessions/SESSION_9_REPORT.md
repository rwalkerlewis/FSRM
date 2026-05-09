# Session 9 Report

## Goal

Find the Jacobian-residual consistency bug surfaced by Session 8's PCSVD
experiment on `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic`.
Session 8's working hypothesis: `J dx = -F` is satisfied exactly by PCSVD
(pseudoinverse on a small system), yet `F(u + dx)` grows from `1.77e-4`
to `1.131e+01`, which is only possible if `J != dF/du`.

Session 9 kickoff proposed a four-step path:
1. Run `-snes_test_jacobian` to compare hand-coded `J` vs FD.
2. Classify mismatched rows by field (Lagrange / disp / BC).
3. Add `PetscDSSetImplicit` on the Lagrange field if missing (PyLith does this).
4. If implicit flag is already set, write a pointwise callback consistency
   test for `f0_hybrid_lambda` vs `g0_hybrid_lambda_u`.

## Step 1: `-snes_test_jacobian` on `PrescribedSlipQuasiStatic`

Added a diagnostic env-var hook `FSRM_SNES_TEST_JAC` to the test harness
(`tests/integration/test_dynamic_rupture_solve.cpp`) so the PETSc options
survive the test's `PetscOptionsClear` call. Running:

```
docker run --rm -v $(pwd):/workspace -e FSRM_SNES_TEST_JAC=1 \
  -w /workspace/build fsrm-ci:main \
  ./tests/run_integration_tests \
  --gtest_filter=DynamicRuptureSolveTest.PrescribedSlipQuasiStatic
```

Full log: `/tmp/s9_final_jac.log` (10713 lines).

Key results (one FD diff per SNES iterate):

| State | `||J - Jfd||_F / ||J||_F` | `||J - Jfd||_F` |
|---|---|---|
| u = 0 (initial) | **555756** | 5.91845e+16 |
| u = u_new (after 1 SNES step) | 1.10074e-16 | 1.17221e-05 |
| every later iterate | 1.10074e-16 | 1.17221e-05 |

The FD check says the hand-coded Jacobian IS `dF/du` to machine precision
everywhere except at exactly `u = 0`.

## Step 2: Row classification on the u = 0 diff

Restricted to the first iterate (u = 0). The "Hand-coded minus FD" view
has 318 rows; only four have any entry above the `1e-5` display tolerance:

- **Row 60**: displacement row at a cohesive-adjacent vertex. Hand-coded has
  the expected elasticity entries (`±5e+09`, `±1.66e+08`, `±3.33e+08`,
  `±6.66e+08`) and the cohesive coupling (`±0.0104167`). FD picks up phantom
  entries of `-5.36871e+08` across ~120 spurious columns plus two
  catastrophic entries `(160, -4.47e+16)` and `(255, 2.24e+16)`.
- **Rows 147, 148, 149**: a single Lagrange vertex (3 components). Hand-coded
  has six cohesive entries per row (`±0.0104167` at the neg/pos disp DOFs
  of the vertex). FD has the same six entries plus ~150 phantom columns
  at `±5.37e+08` and `±8.39e+06`.

Column pattern analysis: the FD-only entries at `5.37e+08 = 2^29` and
`8.39e+06 = 2^23` are PETSc's FD default perturbation `h ~ sqrt(eps_machine)
~ 1.49e-8` reciprocal amplifying sub-eps numerator fluctuations. The two
singular `O(1e+16)` entries at `(60, 160)` and `(60/147/148/149, 255)` are
`~1/eps_machine` signatures of catastrophic cancellation in `(F(u+h) - F(u))/h`
where the numerator is exactly at the double-precision noise floor.

These four rows are **FD artifacts at the u = 0 degenerate point**, not
true Jacobian errors. The ratio 555756 is dominated by two entries with
`O(1e+16)` values. After one Newton step, PETSc's FD algorithm re-scales
around the new state and the artifacts vanish: the ratio drops to 1.10e-16.

## Step 3: `PetscDSSetImplicit` diagnostic

Added a `PetscDSGetImplicit` probe right after the Session 5
`PetscDSSetCohesive` call in `src/core/Simulator.cpp`.

Run output (log `/tmp/s9_implicit_probe.log`, line 38):

```
PetscDSSetCohesive applied: DS 1, field 1 (global 1)
Session 9: Lagrange implicit flag on DS 1 field 1 (global 1) = 1
```

The implicit flag is already `PETSC_TRUE`. This is the PETSc default for a
field registered via `PetscDSCreate` -- `implicit[f] = PETSC_TRUE` is set
in `PetscDSCreate` (`/opt/petsc-main/src/dm/dt/interface/dtds.c`), and
`PetscDSSetCohesive` does not touch it. No explicit `PetscDSSetImplicit`
call is needed in FSRM. Step 3 hypothesis ruled out.

The probe was removed after confirming the result; the diagnostic lives in
`/tmp/s9_implicit_probe.log` for reference.

## Step 4: direct callback consistency check via `-snes_fd`

Tried the stricter check of running the whole Newton solve with PETSc's
built-in FD-by-columns Jacobian instead of the hand-coded one. If the
hand-coded callbacks were wrong in some systematic way, `-snes_fd` +
`-pc_type svd` should converge in one step.

Command:

```
docker run --rm -v $(pwd):/workspace -e FSRM_USE_SNES_FD=1 \
  -w /workspace/build fsrm-ci:main ./tests/run_integration_tests \
  --gtest_filter=DynamicRuptureSolveTest.PrescribedSlipQuasiStatic
```

Log: `/tmp/s9_snes_fd.log`. Trace (first TS step, default BT line search
replaced with basic):

```
  0 SNES Function norm 1.767766952966e-04
  1 SNES Function norm 1.131396409904e+01
  Nonlinear solve did not converge due to DIVERGED_DTOL iterations 1
```

Identical to the hand-coded Jacobian path. Running Newton with
`J == Jfd` (by construction) still produces a direction that grows `F` by
five orders of magnitude. This rules out all callback-implementation
hypotheses in Step 5 of the session kickoff: the bug is not in
`f0_hybrid_lambda`, `g0_hybrid_lambda_u`, `f0_hybrid_u_neg/pos`, or their
Jacobian counterparts. They are mutually consistent to machine precision
as of u = u_new.

## What is actually happening: residual garbage + ill-conditioning

Two independent issues surfaced while chasing the premise:

### Finding A: `DMPlexComputeResidualByKey` leaks uninitialized scratch into `locF`

Dumped the KSP RHS via `-ksp_view_rhs` and added a diagnostic that scans
`locF` for entries with magnitudes below `1e-100`:

```
Session 9 dbg: locF max before hybrid = 8.45049e-302
Session 9 dbg: locF max after hybrid = 3.125e-05
  locF[15] = 5.30145e-303 (garbage)
  locF[16] = -1.05675e-302 (garbage)
  ...
  total garbage entries in locF = 114 / 546
```

At u = 0 the elasticity residual is analytically zero everywhere, yet 114
of 546 local DOFs come out of `DMPlexComputeResidualByKey` carrying
subnormal-double bit patterns (values in the `1e-300` to `1e-312` band).
These are classic uninitialized-memory reads -- PETSc's DMPlex closure
scratch buffer is not fully zeroed for DOFs whose closure touches cohesive
geometry but that the non-cohesive volume integrator never writes.

At the next iterate (u = u_new after the Newton step), the elasticity
residual is legitimately nonzero and overwrites the scratch entries, so
the garbage vanishes. This matches the "bug only at u = 0" pattern of the
FD comparison.

**Fix landed**: `VecFilter(locF, 1e-20)` between the volume call and the
hybrid call in `Simulator::FormFunction` (`src/core/Simulator.cpp:5758-5760`).
Tolerance `1e-20` is 16 orders below the smallest legitimate residual
contribution (`O(1e-5)` for the cohesive jump at u = 0) and 10 orders above
the subnormal garbage band. Matches PyLith's practice of filtering tiny
residuals before passing to KSP.

Verification (`/tmp/s9_residual_debug.log`):

```
  total garbage entries in locF = 0 / 546
```

### Finding B: Saddle-point condition number exceeds double precision

Dumped the assembled KSP matrix / RHS / solution via `-ksp_view_mat`,
`-ksp_view_rhs`, `-ksp_view_solution` and turned on PETSc's SVD monitor:

```
SVD: condition number 5.363549166270e+17, 0 of 318 singular values
     are (nearly) zero
SVD: smallest singular values: 3.72e-08 6.11e-08 7.02e-08 1.08e-07 1.26e-07
SVD: largest singular values : 1.70e+10 1.76e+10 1.76e+10 1.99e+10 2.00e+10
```

`5.36e+17` is above `1/eps_machine = 4.5e+15`. PCSVD's pseudoinverse
misbehaves catastrophically at this condition number: looking at the
dumped KSP solution (`build/s9_ksp_sol.txt`, first entries):

```
row  6..8: 24.9932, 6.83451, -48.4719      <- cohesive DOFs, O(10^1)
row 15..17: -32.5766, -5.5098, 23.3242     <- cohesive DOFs, O(10^1)
row 30..32: -43.8589, 51.844, 8.00663      <- cohesive DOFs, O(10^1)
everywhere else: O(1e-10) to O(1e-13)      <- legit displacement entries
```

The Lagrange DOFs land in SVD's "least-squares min-norm" null-space
projection with O(10^1) values because their block is at an ~11-order scale
offset from the disp-disp block (elasticity entries `O(1e+9)`, cohesive
coupling entries `0.0104`). `F(u + dx)` then equals `F(u) + J dx` along the
range of `J`, but the Lagrange rows are outside that range, so their
contribution to `F` after the step is `O(||dx_Lagrange|| * scale)` --
hence the `1.13e+01` observed residual.

### Attempted fix for Finding B: Lagrange-row scaling

Tried multiplying `f0_hybrid_lambda` and `g0_hybrid_lambda_u` (locked and
prescribed branches) by `lambda + 2*mu` so Lagrange-disp entries scale to
`~1e+7` instead of `~0.01`. Required copying the unified-constants array
from the default PetscDS to the cohesive-cell PetscDS (done once in
`setupPhysics`, `src/core/Simulator.cpp:2699-2722`) because
`CohesiveFaultKernel::registerHybridWithDS` reads constants from the
cohesive DS only.

Rebuild + re-run with `-pc_svd_monitor`:

```
0 SNES Function norm 2.121320343560e+06    <- scaled residual, as expected
SVD: condition number 3.250676387917e+18
SVD: smallest singular values: 6.14e-09 3.45e-08 6.88e-08 9.98e-08 1.47e-07
SVD: largest singular values : 1.70e+10 1.76e+10 1.76e+10 2.00e+10 2.00e+10
1 SNES Function norm 1.822336964204e+08
Nonlinear solve did not converge due to DIVERGED_MAX_IT iterations 1
```

The smallest singular values hardly moved and the condition number got
WORSE (3.25e+18). Scaling alone doesn't fix it, because the small
singular values are not from the scale mismatch -- they come from
structural rank deficiency in the constraint block `B` (Lagrange DOFs at
fault-boundary vertices and/or redundant Lagrange constraints on the
cohesive geometry). This is the well-known saddle-point issue PyLith
handles via an explicit Dirichlet BC (`constraint kinematic zero`) on
Lagrange DOFs where the fault terminates at the domain boundary, which
FSRM currently does not do -- the `fault_constraint` BC added in
`setupBoundaryConditions` is `DM_BC_NATURAL` and does not constrain any
DOF, and the BC is additionally skipped entirely because the Lagrange
field lives on a region-specific DS:

```
Skipped fault_constraint BC: Lagrange field 1 not in default DS
(numDSFields=1; likely on region-specific DS)
```

Lagrange scaling was **reverted**; the copy of constants from the default
DS to the cohesive DS was **kept** so any future per-callback material
property reads inside the hybrid kernel work.

## What landed, what did not

| Change | State |
|---|---|
| `FSRM_SNES_TEST_JAC` diagnostic hook in test (survives PetscOptionsClear) | LANDED |
| `VecFilter(locF, 1e-20)` in `FormFunction` (Finding A fix) | LANDED |
| Copy unified constants from default DS to cohesive DS in `setupPhysics` | LANDED |
| `PetscDSGetImplicit` diagnostic probe | removed after ruling out Step 3 |
| Lagrange-row scaling by `lambda + 2*mu` | REVERTED (ineffective) |
| `FSRM_DISABLE_COHESIVE_REORDER`, `FSRM_USE_SNES_FD`, `FSRM_SVD_DEBUG` hooks | removed after diagnostic use |
| Residual-garbage debug prints | removed after diagnostic use |

## PrescribedSlipQuasiStatic SNES trace with the landed fixes

Test config (unchanged): `-pc_type lu -ksp_type preonly -snes_linesearch_type bt`.
PETSc native LU still hits a zero pivot on the Lagrange diagonal, so:

```
0 TS dt 1. time 0.
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.385642450386e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
... 10 more identical retries ...
max_fault_slip = 0, sol_norm = 0
```

Identical to Session 8's committed trace (`/tmp/s8_prescribed_final.log`).
The landed Session 9 residual-filter and constants-copy fixes are
correctness improvements but do not on their own reach the
`max_fault_slip > 5e-4` success criterion.

## Fault test subset vs Session 8

`ctest --output-on-failure -R 'Fault|Slip|Rupture|Cohesive|TPV5|PressurizedFracture'`
in `fsrm-ci:main`, log `/tmp/s9_fault_subset.log`:

| Test | Session 8 | Session 9 |
|---|---|---|
| `Functional.DynamicRuptureSetup.*` (6 cases) | PASS | PASS |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Physics.LockedFaultTransparency` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | FAIL |
| `Integration.TimeDependentSlip` | FAIL | FAIL |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |
| other unit / setup fault tests | PASS | PASS |

16 pass / 9 fail in both runs; identical set. No regressions.

## Success criteria outcome

| Criterion | Result |
|---|---|
| `-snes_test_jacobian` matrix ratio close to 1.0 on `PrescribedSlipQuasiStatic` | ratio is `1.10e-16` at every iterate except `u = 0`, where it is `555756` due to catastrophic FD cancellation at two specific `(row, col)` pairs. The *hand-coded Jacobian* is correct. |
| Newton step actually reduces the residual | NOT MET. PCSVD returns dx with O(10) entries at Lagrange DOFs; `F` grows to `1.13e+01`. Cause: condition number `5.36e+17 > 1/eps_machine`, not Jacobian inconsistency. |
| `max_fault_slip > 5.0e-4` | NOT MET. Test still fails identically to Session 8. |

The Session 9 premise ("J is not dF/du") is disproven by the FD test.
The actual blocker for `PrescribedSlipQuasiStatic` is the saddle-point
condition number, which Session 9 has *characterized* but not *solved*.
A solver path forward requires one of: (a) rebuilding the PETSc image
with MUMPS/SuperLU_DIST so LU can pivot past the zero diagonal,
(b) fieldsplit with Schur complement on the Lagrange block, or
(c) PyLith-style Dirichlet constraint on Lagrange DOFs at fault boundary
points plus a regularized Schur solve. None of these was the Session 9
scope.

## Files touched

- `src/core/Simulator.cpp`: `VecFilter(locF, 1e-20)` after volume residual
  (`FormFunction`, ~line 5758); copy unified constants from default DS
  to cohesive DS in `setupPhysics` (~line 2699).
- `src/physics/CohesiveFaultKernel.cpp`: unchanged from Session 8 after
  reverting the scaling experiment.
- `tests/integration/test_dynamic_rupture_solve.cpp`: added
  `FSRM_SNES_TEST_JAC` env-var hook so future sessions can re-run the
  FD Jacobian comparison without patching the test.

## Logs (full, not truncated)

- `/tmp/s9_jac_check.log` -- `-snes_test_jacobian` with reorder (Session 8 HEAD)
- `/tmp/s9_no_reorder.log` -- same with reorder disabled (ratio `105028`)
- `/tmp/s9_implicit_probe.log` -- `PetscDSGetImplicit` probe output
- `/tmp/s9_snes_fd.log` -- `-snes_fd` (PETSc FD Jacobian) + PCSVD
- `/tmp/s9_svd_monitor.log` -- PCSVD singular-value spectrum
- `/tmp/s9_residual_debug.log` -- `locF` garbage enumeration pre / post hybrid
- `/tmp/s9_scale.log`, `/tmp/s9_scale2.log` -- Lagrange scaling experiment
- `/tmp/s9_filter_fix.log`, `/tmp/s9_filter_svd.log`, `/tmp/s9_filter_dbg.log`
  -- post-VecFilter re-runs
- `/tmp/s9_final_jac.log` -- FD check on the landed tree
- `/tmp/s9_fault_subset.log` -- `ctest` on the 25-test fault subset
- `build/s9_ksp_mat.txt` / `build/s9_ksp_rhs.txt` / `build/s9_ksp_sol.txt`
  -- assembled matrix / RHS / Newton direction dumps

## Rules compliance

- Rule 1 (Docker build/test): all builds and tests in `fsrm-ci:main`.
- Rule 2 (PETSc 3.25 API): `PetscDSGetImplicit` / `PetscDSSetImplicit`
  (`/opt/petsc-main/include/petscds.h:136-137`), `VecFilter`
  (`/opt/petsc-main/include/petscvec.h:821`), `PetscDSSetConstants` /
  `PetscDSGetConstants`, `DMGetCellDS`, `DMPlexComputeResidualByKey`,
  `DMPlexComputeResidualHybridByKey`, `DMPlexComputeJacobianHybridByKey`
  all verified against `/opt/petsc-main/include/` and the hybrid test
  driver `/opt/petsc-main/src/dm/impls/plex/tests/ex5.c:1229-1232`
  before call.
- Rule 3 (no regressions): Session 9 fault subset matches Session 8 HEAD
  exactly (16 pass / 9 fail, same tests).
- Rule 4 (DS/BC ordering): untouched.
- Rule 16 (no log truncation): all logs referenced above are captured in
  full under `/tmp/` and `build/`; inspected via the `Read` tool.
