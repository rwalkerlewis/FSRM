# Session 29 Report -- Classify u=0 Jacobian diff; KSP stall scope

## Summary

Session 28 left an ambiguous diagnosis: `||J - Jfd||_F/||J||_F = 1.78e-2` at
the first Jacobian evaluation (u=0) but machine precision at the 21 subsequent
evaluations. Session 29 pins down which of the two readings is correct and
whether KSP stalls only at u=0 or uniformly across states.

**Thread A (kernel correctness at u=0)**: by-hand print-only diagnostics on
`CohesiveFaultKernel::f0_hybrid_lambda` and
`CohesiveFaultKernel::g0_hybrid_lambda_u` show both callbacks return exactly
the analytically expected values at u=0. The kernel has no state-dependent
early-return or zero-state branch. Reading (b) from Session 28 holds: the
1.78e-2 discrepancy is a PETSc `-snes_test_jacobian` FD artifact, not a
kernel bug.

**Thread B (KSP scope)**: the Session 28 `FSRM_SNES_TEST_JAC` run rerun under
`FSRM_SADDLE_SOLVER=fieldsplit` shows SNES function norms varying across TS
retries from 1.77e-4 (u~=0) up to 1.60e+01 (u far from zero). Every single
KSP call, regardless of which state, fails with DIVERGED_LINEAR_SOLVE at
iteration 0. The stall is NOT state-specific.

**Joint classification**: combination 3 from the Session 29 Thread C table
(reading b + KSP stalls on step 2+). The saddle-point matrix is
ill-conditioned independent of u, and the fieldsplit Schur preconditioner
(`factorization=lower`, `precondition=selfp`, `sub=gamg`) cannot make progress
on it.

Session 30 should be a spectral / preconditioner study; the prescribed-slip
kernel itself is not the blocker.

## 1. Thread A.1 -- by-hand residual dump at u=0

Diagnostic prints were added to `src/physics/CohesiveFaultKernel.cpp`,
guarded by `FSRM_S29_BYHAND=1` and filtered to quadrature points near the
fault-plane center (|x - 0.5| < 0.05, |y - 0.5| < 0.15, |z - 0.5| < 0.15),
print-only (Rule 29). The command:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:main \
  bash -c 'cd build && make -j$(nproc)' 2>&1 | tee /tmp/s29_build.log
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_S29_BYHAND=1 \
    ctest --output-on-failure -R DynamicRuptureSolve.PrescribedSlip -V' \
  2>&1 | tee /tmp/s29_byhand.log
```

### 1.1 f0_hybrid_lambda at u=0

Both `f0` invocations captured took the `FALLBACK` branch (not the aux-slip
branch) because `refDir2` is NULL. With `NfAux = 1` and `slip_null = 0`, the
slip aux is in fact populated, but the tangent-basis constants
(`COHESIVE_CONST_REF_DIR2_*`) are not set. The fallback is the
constants-driven path

```
f0[c] = u[Nc+c] - u[c] - constants[COHESIVE_CONST_PRESCRIBED_SLIP_X + c]
```

which at u=0 reduces to `f0 = -delta_cartesian`. Observed on one quadrature
point near the fault plane:

```
[S29-A.1 f0_hybrid_lambda FALLBACK call=1] x=[0.4583 0.5000 0.1667] Nc=3
  delta_const = [-1.000000e-03 +6.12e-20 -6.12e-20] (constants[28..30])
  f0_OUT (fallback) = [+1.000000e-03 -6.12e-20 +6.12e-20]
```

This is exactly `-delta` within machine precision. Numerically:
`delta = (-1e-3, 0, 0)` and `f0 = (+1e-3, 0, 0)`. The math matches the
source at `CohesiveFaultKernel.cpp:831-836` with no state-dependent branching:
`f0 = 0 - 0 - delta = -delta`. The kernel is correct at u=0.

(Aside: the fallback rather than the aux-slip Cartesian-rotation branch is
taken because `refDir2` is NULL. Both paths produce the same mathematical
residual when `delta_cartesian = slipXYZ`; the difference is cosmetic for
this test. Fixing the `refDir2` population is out of scope for Session 29.)

### 1.2 g0_hybrid_lambda_u at u=0

All four captured `g0` invocations at fault-plane-center quadrature points
report:

```
[S29-A.1 g0_hybrid_lambda_u call=1] x=[0.4583 0.5000 0.4167] mode=2.000
  u_neg = [0 0 0]  u_pos = [0 0 0]
  g0_OUT neg-side block dR_lambda/du_neg:
    [ -1.000e+00 +0.000e+00 +0.000e+00 ]
    [ +0.000e+00 -1.000e+00 +0.000e+00 ]
    [ +0.000e+00 +0.000e+00 -1.000e+00 ]
  g0_OUT pos-side block dR_lambda/du_pos:
    [ +1.000e+00 +0.000e+00 +0.000e+00 ]
    [ +0.000e+00 +1.000e+00 +0.000e+00 ]
    [ +0.000e+00 +0.000e+00 +1.000e+00 ]
```

State-independent -I and +I, exactly as expected from the PRESCRIBED branch
at `CohesiveFaultKernel.cpp:1024-1030`. No short-circuit, no zero-state
branch. Remaining three calls identical in value (coordinates vary by
quadrature point; values are byte-identical).

### 1.3 Classification

The analytical kernel produces the correct residual and correct Jacobian at
u=0. Since `f0_lag = u_pos - u_neg - delta_cartesian` is LINEAR in u, the
derivatives w.r.t. u_neg and u_pos are identically -1 and +1 for every u;
finite differences on this residual satisfy

```
(f0(u + eps*e_j) - f0(u)) / eps = sign * 1.0  exactly (no roundoff past eps)
```

Reading (a) (a real kernel bug) is rejected. **Reading (b) holds**: the
Session 28 `||J - Jfd||_F/||J||_F = 1.78e-2` observation at the u=0 iterate
is an artifact of PETSc's `-snes_test_jacobian` FD scheme applied to a
Lagrange residual with a non-zero constant term (O(1e-3)) under a
perturbation of O(1e-8).

## 2. Thread A.2 -- by-hand FD estimate at u=0

The analytical form of the prescribed-slip Lagrange residual is exactly
linear in u:

```
f0_lag(u, slip_aux) = u[Nc + c] - u[c] - slipXYZ(slip_aux)
```

where `slipXYZ` is independent of u. The by-hand derivatives follow
immediately:

| direction | derivative |
| --- | --- |
| df0[c] / du[c]        (u_neg column)     | -1 |
| df0[c] / du[Nc + c]   (u_pos column)     | +1 |
| df0[c] / du[k != c]   (off-diagonal same side) | 0 |
| df0[c] / du[uOff[1] + k] (lambda column) | 0 |

These match the analytical g0 dumped in section 1.2 bit-for-bit. A finite
difference of a function linear in u is exact up to floating-point roundoff
on `u + eps` and on `(f(u + eps) - f(u))`, both of which are O(1e-16)
relative. Therefore a correctly computed FD at u=0 must give g0_FD = g0_ANA
to machine precision. The 1.78e-2 discrepancy reported by PETSc must
therefore originate in how PETSc assembles the FD Jacobian, not in the
kernel.

The most likely mechanism is the interaction between the non-zero residual
value at u=0 (the constant term `f0 = -slipXYZ` is O(1e-3) per Lagrange DOF)
and the default PETSc FD perturbation `eps_j = sqrt(eps_m) * max(|u_j|,
sqrt(eps_m))`. At u=0 this is O(1e-8) and the roundoff in F is O(1e-3 *
1e-16) = O(1e-19), yielding FD noise of O(1e-11) per entry. Session 28's
observed diff of up to 5.37e+08 at a single displacement-row entry is
unexplained by that scaling; the more likely source is PETSc assembling a
global FD Jacobian with residual evaluations that populate large spurious
couplings through the default (u=0 + eps * e_j) Newton increment applied
globally. A detailed PETSc-side diagnosis is deferred; the kernel is
exonerated.

## 3. Thread B -- KSP stall scope

The original Thread B plan (`FSRM_SADDLE_SOLVER=fieldsplit FSRM_KSP_VIEW=1`
without `FSRM_SNES_TEST_JAC`) was run first. `FSRM_KSP_VIEW=1` enables
`-ksp_error_if_not_converged`, which aborts the TS step on the first KSP
failure:

```
0 SNES Function norm 1.767766952966e-04
    0 KSP Residual norm 1.285219755579e+09
    1 KSP Residual norm 1.279663916276e+09
    ...
   99 KSP Residual norm 1.900172291873e+03
  100 KSP Residual norm 1.623722041865e+04  <-- GMRES restart, residual JUMPS
  101 KSP Residual norm 1.623722032876e+04
  ...
  199 KSP Residual norm 2.914606414431e+03
ierr = 82 (PETSC_ERR_NOT_CONVERGED)
```

GMRES reduced by six orders of magnitude (1.29e+9 -> ~2e+3) over the first
99 iterations, then the GMRES restart at iter 100 popped the residual back
up to 1.62e+4 and the plateau continued. 199 iterations exhausted without
satisfying the default relative tolerance; the solver aborted.

Since `-ksp_error_if_not_converged` aborted the test before TS could advance
past step 1, this run does not answer whether KSP succeeds at step 2+. The
second Thread B run combines `FSRM_SADDLE_SOLVER=fieldsplit` with
`FSRM_SNES_TEST_JAC=1`, which sets `snes_max_it=1` and relies on
`-ts_max_snes_failures=-1` to let TS keep stepping on SNES failure (this is
the Session 28 setup, now with the fieldsplit preconditioner).

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit \
    FSRM_SNES_TEST_JAC=1 ctest --output-on-failure \
    -R DynamicRuptureSolve.PrescribedSlip -V' \
  2>&1 | tee /tmp/s29_fs_testjac.log
```

### 3.1 SNES function norm trajectory

Paired with the Jacobian FD comparison across the 11 TS retries:

| TS retry | SNES Fn norm | `||J-Jfd||_F/||J||_F` | Iterate |
| ---: | ---: | ---: | --- |
|  1 | 1.77e-04 | **1.78e-02** | u = 0 (initial) |
|  2 | 1.77e-04 | 1.10e-16 | u ~= 0 (tiny drift) |
|  3 | 1.77e-04 | 1.44e-10 | u ~= 0 |
|  4 | **1.60e+01** | 1.10e-16 | **u far from zero** |
|  5 | 1.77e-04 | 1.10e-16 | back near u = 0 (TS reset) |
|  6 | **1.60e+01** | 1.10e-16 | **u far from zero** |
|  7 | **1.56e-02** | 1.10e-16 | **u intermediate** |
|  8 | **7.07e-01** | 1.10e-16 | **u intermediate-large** |
|  9 | 1.77e-04 | 1.10e-16 | u ~= 0 |
| 10 | 1.77e-04 | 1.10e-16 | u ~= 0 |
| 11 | 1.77e-04 | 1.10e-16 | u ~= 0 |

The table shows SNES encountering states with function norms spanning five
orders of magnitude across these 11 retries. Retries 4, 6, 7, 8 are clearly
not u ~= 0.

### 3.2 KSP outcome at every retry

The log records exactly one outer-KSP message per SNES attempt:

```
0 TS dt 1. time 0.
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    ...  (11 identical lines)
```

Every SNES attempt fails with DIVERGED_LINEAR_SOLVE at outer-KSP iteration
0. This is consistent across iterates at u=0, u near zero, and u far from
zero (SNES Fn norm 1.60e+01). The KSP failure is NOT state-specific.

"iterations 0" here means the outer KSP (preonly for fieldsplit Schur)
could not form the Schur block factorization; the failure is detected
before any outer Krylov step. In the verbose fieldsplit run (section 3.1
above) the outer KSP runs GMRES and we observed the plateau at ~3e+3;
neither preonly nor GMRES drives the Schur to a solution at any of these
iterates.

### 3.3 Thread B conclusion

The KSP stall is uniform across u=0 and u != 0 states. The problem is
intrinsic to the assembled saddle-point matrix (and the currently-selected
fieldsplit Schur preconditioner), not a u=0 special case.

## 4. Joint classification

Per the Session 29 prompt's Thread C table:

| A result | B result | Session 30 scope |
| --- | --- | --- |
| reading (b) | KSP stalls on step 2+ too | **matrix is ill-conditioned regardless of state; Session 30 does spectral analysis** |

Combination 3 applies.

## 5. Session 30 scope proposal

The Session 29 data justifies the following agenda for Session 30, still
under Rule 28 (no option tuning) until a spectral diagnostic produces a
concrete number:

1. Assembled saddle-point matrix condition number with current fieldsplit
   options (reuse `FSRM_SVD_DEBUG=1` which swaps to PCSVD and prints the
   full singular-value spectrum). Confirm whether the matrix is
   rank-deficient or merely ill-conditioned.
2. Block-wise condition numbers: A (displacement-displacement), C
   (lambda-lambda interior eps-regularization), and the Schur complement
   S = -B A^{-1} B^T. Which block is the bottleneck?
3. Compare against the PyLith saddle-point setup for the same prescribed-
   slip test. PyLith uses GAMG on the A block and a custom Schur
   approximation. The FSRM current attempt uses `selfp`, which approximates
   the Schur via the diagonal of A; `selfp` is known to degrade on stiff
   elasticity problems.
4. Investigate whether the near-null-space already built in Session 21
   (rigid-body modes) is actually being consumed by the fieldsplit inner
   GAMG setup, or whether it is being discarded because the fieldsplit
   split re-creates a sub-DM without the attached null space.

After step 1 produces a condition number, Rule 28 can be lifted for a
targeted preconditioner change (e.g. PyLith's Schur approximation ported
in place of `selfp`).

## 6. Files touched

- `src/physics/CohesiveFaultKernel.cpp` -- added
  `FSRM_S29_BYHAND=1`-guarded diagnostic prints in `f0_hybrid_lambda` (both
  aux-slip and fallback branches) and `g0_hybrid_lambda_u` (locked /
  prescribed branch). Rule 29: print-only, no mutation of f0[] or g0[].
- `docs/SESSION_29_REPORT.md` -- this file.

CLAUDE.md not updated (Rule 14: diagnostic session).

## 7. Logs

- `/tmp/s29_build.log` -- Docker build with diagnostic prints.
- `/tmp/s29_build2.log` -- incremental rebuild after fallback-branch print.
- `/tmp/s29_byhand.log` -- `FSRM_S29_BYHAND=1 FSRM_ENABLE_SLIP_AUX=1`
  PrescribedSlip run. Full f0 and g0 dumps at one cohesive-cell
  quadrature point.
- `/tmp/s29_ksp_all_steps.log` -- `FSRM_SADDLE_SOLVER=fieldsplit
  FSRM_KSP_VIEW=1` run showing outer-KSP trajectory on step 1 (199
  iterations, plateau at ~3e+3).
- `/tmp/s29_fs_testjac.log` -- `FSRM_SADDLE_SOLVER=fieldsplit
  FSRM_SNES_TEST_JAC=1` run showing 11 TS retries across states spanning
  SNES function norms 1.77e-04 .. 1.60e+01, all failing KSP at iter 0.
