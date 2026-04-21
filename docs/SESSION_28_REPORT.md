# Session 28 Report -- Jacobian Audit: FD check

## Summary

Step 1 (FD vs hand-coded Jacobian) ran to completion and produced a
surprising, non-uniform result:

| TS step | `||J_hand - J_fd||_F / ||J_hand||_F` | `||J_hand - J_fd||_F` |
| ------: | --- | --- |
| 1 (u=0) | **1.78e-2** | 1.90e+9 |
| 2 | 1.10e-16 | 1.17e-5 |
| 3 | 1.44e-10 | 15.3 |
| 4..22 | ~1.10e-16 | ~1.17e-5 |

The decision gate in the Session 28 prompt says "`>= 1e-8` => a kernel
bug, stop for a targeted kernel fix." Per that gate, Step 1 fails:
1.78e-2 is six orders of magnitude above threshold.

However, the failure is confined to the FIRST time step, before any
Newton update is applied. At every subsequent state (u != 0) the hand-
coded Jacobian agrees with FD to machine precision (1e-16). This splits
the diagnosis:

- Reading (a): a real kernel bug whose manifestation requires u = 0. The
  kernel is correctly computing derivatives from a subsequent iterate
  but misses something intrinsic to u = 0 (for example, branch on
  `|slip| > eps` that defaults wrong at zero slip).
- Reading (b): PETSc's `-snes_test_jacobian` FD epsilon scaling degrades
  at u = 0. The default perturbation is `eps_j = 1e-8 * max(|u_j|, 1e-8)`
  so at u = 0 the absolute perturbation is 1e-8 and the residual noise
  floor dominates for any F whose non-constraint magnitude is small.
  The apparent 1.78e-2 then is an FD artifact, not a Jacobian bug.

Which reading is correct matters for Session 29's scope: reading (a)
wants a kernel fix; reading (b) reclassifies Session 27's KSP
convergence problem as pure ill-conditioning, since KSP is applied to
iterates 2..22 where the Jacobian matches FD at machine precision.

Per Session 28 Rule 28 (no option tuning, no speculative fixes), this
session stops after Step 1. Steps 2 (structure) and 3 (condition
number) are deferred so Session 29 can proceed from a clear
classification of the u = 0 discrepancy first.

## 1. Command and environment

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  bash -c 'FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit \
    FSRM_SNES_TEST_JAC=1 \
    ctest --output-on-failure -R DynamicRuptureSolve.PrescribedSlip -V' \
  2>&1 | tee /tmp/s28_jac_test.log
```

Test: `DynamicRuptureSolveTest.PrescribedSlipQuasiStatic` (ctest id 95).
`FSRM_SNES_TEST_JAC=1` at
`tests/integration/test_dynamic_rupture_solve.cpp:325-328` sets
`-snes_test_jacobian`, `-snes_test_jacobian_view`, `-snes_max_it 1`.

Mesh (from the test log header): 384 bulk cells, 416 total after
`DMPlexConstructCohesiveCells` (32 prisms inserted), 96 fault vertices.
The assembled Jacobian viewed here is a `seqaij` matrix of rank 255
(indices 0..254 in the log).

## 2. Full Frobenius-norm sequence

Extracted with `grep -E "\\|\\|J - Jfd" /tmp/s28_jac_test.log`:

```
||J - Jfd||_F/||J||_F = 0.0178018,   ||J - Jfd||_F = 1.89578e+09
||J - Jfd||_F/||J||_F = 1.10021e-16, ||J - Jfd||_F = 1.17165e-05
||J - Jfd||_F/||J||_F = 1.43853e-10, ||J - Jfd||_F = 15.3194
||J - Jfd||_F/||J||_F = 1.10321e-16, ||J - Jfd||_F = 1.17485e-05
||J - Jfd||_F/||J||_F = 1.10021e-16, ||J - Jfd||_F = 1.17165e-05
||J - Jfd||_F/||J||_F = 1.10021e-16, ||J - Jfd||_F = 1.17165e-05
||J - Jfd||_F/||J||_F = 1.10233e-16, ||J - Jfd||_F = 1.17391e-05
... (remaining 14 entries all at ~1.10e-16)
```

22 Jacobian evaluations in total. Only the first (TS step at t=0 with
the zero initial guess) shows any meaningful disagreement. TS step 3
spikes to 1.44e-10; that is still well above machine precision but
many orders below the Step-1 decision threshold of 1e-8.

## 3. Which rows are bad at u = 0

From the Step-1 "Hand-coded minus finite-difference Jacobian with
tolerance 1e-05" block (log lines 609-870). Rows with non-zero diff
entries (256 rows total, 11 had diffs):

| Row | Diff nnz | Min col | Max col | Diff value range | Likely DOF |
| --: | ------: | ------: | ------: | --- | --- |
| 29 | 1 | 163 | 163 | +5.37e+08 (single entry) | displacement |
| 30 | 254 | 1 | 254 | [-3.46e+07, -3.36e+07] | Lagrange |
| 31 | 253 | 1 | 254 | -1.05e+06 (uniform) | Lagrange |
| 32 | 254 | 1 | 254 | [-1.05e+06, +5.37e+08] | Lagrange |
| 57 | 254 | 1 | 254 | [-5.70e+08, -3.46e+07] | Lagrange |
| 60 | 254 | 1 | 254 | [-5.37e+08, -1.05e+06] | Lagrange |
| 61 | 254 | 1 | 254 | [-3.46e+07, +5.03e+08] | Lagrange |
| 62 | 254 | 1 | 254 | [-3.46e+07, -3.36e+07] | Lagrange |
| 84 | 254 | 1 | 254 | [-5.70e+08, -3.56e+07] | Lagrange |
| 85 | 254 | 1 | 254 | [-5.37e+08, -1.05e+06] | Lagrange |
| 86 | 253 | 1 | 254 | -1.05e+06 (uniform) | Lagrange |

Assignment "likely DOF" is inferred from the hand-coded matrix row
structure (section 4 below):
- row 29 has large elasticity-scale entries (1.5e9, 5e9, 6.67e8) and a
  few small (~0.03) Lagrange-coupling entries -- displacement-dominant.
- rows 30, 31, 32, 60-62, 84-86 contain only small-magnitude (1/32 and
  1/192 = ~0.031 and ~0.0052) mass-matrix-like entries with no
  elasticity-scale couplings -- Lagrange-dominant.

Interpretation of the 254-entry dense diffs: these rows in the hand-
coded Jacobian carry sparse mass-matrix entries only, but FD perceives
a near-constant coupling across almost every column in the mesh. A
uniform per-row constant (e.g., the entire row 30 showing -3.46e+07 at
every column from 1..254, with two small exceptions at the true mass-
matrix positions) is not a physically plausible Jacobian structure and
is strong evidence for an FD artifact when the residual has non-trivial
values at u = 0 (reading b).

Interpretation of row 29's single missing entry at column 163: column
163 is a low-connectivity displacement row itself (9 non-zeros, diag
2e9, versus the bulk displacement diagonals of 1e10 at 30+ non-zeros),
consistent with a near-fault vertex whose support is cut by the
cohesive prism insertion. A single +5.37e+08 entry missing at (29, 163)
would be a cross-fault elastic coupling (positive side to negative
side through the cohesive kinematic constraint) that the hand-coded
hybrid Jacobian may be failing to stamp when u = 0. This is a
plausible real bug candidate.

## 4. Hand-coded rows at u=0 (ground truth for comparison)

From log lines 96-351 (the "Hand-coded Jacobian" block at TS step 1):

| Row | Hand nnz | Hand diag | Max abs | Characterization |
| --: | ------: | --- | --- | --- |
| 29 | 33 | (row 29 diag ~= 5e9 typical) | 1.5e+09 | displacement, elasticity-scale |
| 30 | 14 | 0.03125 (= 1/32) | 0.031 | Lagrange, eps-regularization only |
| 31 | 14 | 0.03125 | 0.031 | Lagrange |
| 32 | 14 | 0.03125 | 0.031 | Lagrange |
| 60 | 10 | 0.03125 | 0.031 | Lagrange |
| 61 | 10 | 0.03125 | 0.031 | Lagrange |
| 62 | 10 | 0.03125 | 0.031 | Lagrange |
| 84 | 10 | 0.03125 | 0.031 | Lagrange |
| 85 | 10 | 0.03125 | 0.031 | Lagrange |
| 86 | 10 | 0.03125 | 0.031 | Lagrange |
| 163 | 9 | 2e+09 | 2e+09 | low-connectivity displacement (near fault) |
| 144 | 33 | 1e+10 | 1e+10 | typical bulk displacement |
| 160 | 29 | 1e+10 | 1e+10 | typical bulk displacement |
| 200 | 32 | 1e+10 | 1e+10 | typical bulk displacement |

The Lagrange rows 30-32, 60-62, 84-86 carry only the weak `eps=1e-4`
mass-matrix regularization (values 1/32 and 1/192 are the reference-
element integrals `eps * M` on the cohesive prism). They do NOT carry
any displacement coupling in the hand-coded Jacobian at u=0, because
the `DMPlexComputeJacobianHybridByKey` call with `FSRM_ENABLE_SLIP_AUX=1`
produces the prescribed-slip kinematic constraint as `u_pos - u_neg -
slip_aux`; its derivative w.r.t. the jump is ±I, which is local (only
at the same fault vertex), and w.r.t. `slip_aux` is identically zero
because `slip_aux` is read from a DM aux vec and not a column of u.

The FD perceiving a uniform -3.46e+07 across every column of row 30 is
therefore NOT what the prescribed-slip physics predicts. It is either
(a) a bug in the hybrid Jacobian setup that only appears when
`|slip_aux| + |u_jump|` is zero, or (b) PETSc's FD using a too-small
eps at u=0 and picking up the Lagrange residual F_lag = -slip_aux ~=
(0.001, 0, 0) as noise divided by eps = 1e-8.

## 5. Which kernel is at fault (if reading a holds)

The 9 Lagrange rows with dense FD-only diffs and the 1 displacement
row with a single missing cross-fault entry both point to
`DMPlexComputeJacobianHybridByKey`'s cohesive keys:

```
// src/core/Simulator.cpp:7264-7276
keys[0]: material_label_ / value=1 / field=displacement (neg side bulk)
keys[1]: material_label_ / value=2 / field=displacement (pos side bulk)
keys[2]: interfaces_label_ / value=1 / field=lagrange
```

The `g0_hybrid_*` kernels in
`src/physics/CohesiveFaultKernel.cpp` are consumed by key [2] to produce
the Lagrange-vs-displacement coupling. If these kernels short-circuit
at u=0 (for instance, `if (slip_magnitude < tol) return;` style guards
inherited from the slipping-friction branch), they would stamp zero
coupling rows for interior Lagrange DOFs only when the state is zero,
and correct ±I couplings otherwise -- matching the observed pattern.
Session 29 should grep `CohesiveFaultKernel::g0_hybrid_*` for any
state-dependent early return or zero-state branch.

## 6. What this does NOT say about Session 27's KSP failure

Session 27 reported that KSP cannot reduce the residual by 14 orders
of magnitude within 500 iterations. That was measured on the first
TS step (t=0, u=0), which per this session's Step 1 is the one state
where the hand-coded Jacobian disagrees with FD. So the Session 27
KSP failure could be caused by the u=0 Jacobian mismatch.

HOWEVER: in this test, with `snes_max_it=1`, PETSc advanced TS 21
additional steps despite the SNES failure on step 1, and on those 21
subsequent steps the hand-coded Jacobian matches FD to 1e-16. If the
experimental `FSRM_SADDLE_SOLVER=fieldsplit` path also fails KSP on
those subsequent steps (where the Jacobian is provably correct), then
the convergence problem is NOT a Jacobian bug -- it is preconditioning
or matrix conditioning, as the Session 27 prompt's options (1) or (b)
suggested.

Session 28 did not run the fieldsplit + KSP-verbose path a second time
to check whether the subsequent-step KSP also fails, because Rule 28
forbids option tuning and the Step-1 decision gate said to stop. That
check is on the Session 29 agenda and requires no option change (same
env vars as Session 27's run).

## 7. Session 29 recommendation

Based on the above, Session 29 has two threads that should be run in
order:

1. Classify the u=0 discrepancy. Build a reduced 2-element smoke test
   and print the hybrid-Jacobian closure at one cohesive cell at u=0
   by hand (no PETSc FD machinery), then compare against the
   hand-coded closure. If they agree, reading (b) is confirmed and the
   1.78e-2 is an FD artifact. If they differ, reading (a) is confirmed
   and the culprit is a `g0_hybrid_*` state-dependent branch.

2. Then, regardless of reading, check whether the fieldsplit KSP also
   stalls on TS steps 2..22 where the Jacobian is provably correct.
   If yes: the problem is preconditioning/conditioning. If no: the
   u=0 issue is the entire story.

Option tuning (Rule 28) remains banned until one of those two produces
a concrete number (e.g., "the hybrid closure at row 30 col 200
evaluates to X by hand and Y in the assembled matrix, difference Z").

## 8. Files touched

- `docs/SESSION_28_REPORT.md` -- this file.

No source changes. CLAUDE.md not updated (Rule 14: diagnostic session
only).

## Logs

- `/tmp/s28_build.log` -- Docker build (no-op, image already built).
- `/tmp/s28_jac_test.log` -- full `-snes_test_jacobian` run; 7.1 MB,
  22 Jacobian dumps, one per TS step.
