# Session 24 Report -- fieldsplit Schur selfp + rim-pin ON

## Summary

Session 24 ran the never-before-tested combination predicted by the prompt:
**rim-pin ON + aux-slip ON + fieldsplit Schur ON**. With the rim-pin active,
the Lagrange block collapses from 75 DOF to 3 free DOF (the rim-pin
essential-constrains 24 of 25 interior Lagrange vertices' 3 components, leaving
1 interior vertex). PCFieldSplit detects two fields by name and Schur-selfp
preconditions a 3x3 Lagrange Schur complement.

The crucial linear-solver result: **KSP now converges in 3 GMRES iterations
every Newton step** (residual drops 1.37e+08 → 4.27e-04, CONVERGED_RTOL).
This is the structural breakthrough the prompt predicted. Session 23
experimental KSP plateaued at 1e+9 with rim-pin OFF; with rim-pin ON the
selfp Sp is non-singular and GMRES converges cleanly.

However, the SNES Newton step is rejected at the line search:
`DIVERGED_LINE_SEARCH iterations 0 → 1`. The Newton update produced by
the converged KSP is too large (function norm jumps from 6.25e-05 to 5.75e+00),
so the line search backtracks to no-progress. `max_fault_slip = 0`, expected
> 5.0e-4. The test still fails.

Per the decision gate ("KSP converges but slip off-tolerance: tuning issue,
propose Session 25"), the three aux-slip gates are NOT inverted. The
fieldsplit + aux-slip + rim-pin combination is now reachable behind the
existing env vars for Session 25 tuning work.

## Code changes

None. All required infrastructure was already in place from Sessions 21-23:

- Rim-pin: ON by default; opt-out via `FSRM_DISABLE_RIM_PIN=1`
  (`src/core/Simulator.cpp:8391`)
- Aux-slip: opt-in via `FSRM_ENABLE_SLIP_AUX=1`
  (`src/core/Simulator.cpp:2978, 3204, 3352, 7027, 7241`)
- Fieldsplit Schur with PyLith options: opt-in via
  `FSRM_SADDLE_SOLVER=fieldsplit` (`src/core/Simulator.cpp:3922-3961`)
- Displacement FE named `"displacement"` so fieldsplit option keys match
  PyLith's `-fieldsplit_displacement_*` keys (`src/core/Simulator.cpp:1677`)

## Verification 1 -- rim-pin is ON in the experimental run

`buried_cohesive_label_` populated by geometric fallback:

```
Session 14 (kept in 15): no 'buried_edges' label on DM; falling back to
    geometric rim-prism detection (Session 12 pattern).
Session 14 (kept in 15): buried_cohesive label (geometric fallback):
    56 points labeled (candidates 81: seg_prism=56, point_prism=25)
```

DM section after BC insertion confirms the `lagrange_fault_rim` essential BC
is active:

```
Local section: total dof=525, constrained dof=300 (free=225)
    field 0: dof=450, constrained=228, free=222
    field 1: dof=75, constrained=72, free=3
```

field 1 (lagrange) has 75 = 25 × 3 DOFs total, of which 72 = 24 × 3 are
essential-constrained by the rim-pin. Only 1 interior vertex's 3 Lagrange
DOFs remain free in the assembled system. Total assembled matrix is 225 × 225
(222 displacement + 3 lagrange), confirmed in the KSP view at
`/tmp/s24_prescribed.log:4881`.

Rim-pin is ON. The `FSRM_DISABLE_RIM_PIN` env var was not set on the
experimental run.

## Verification 2 -- fieldsplit detected two fields by name

```
Session 23: fieldsplit Schur with PyLith options
    (factorization=lower, precondition=selfp, sub=gamg)
```

KSP view (`/tmp/s24_prescribed.log:93-99`):

```
type: fieldsplit
    FieldSplit with Schur preconditioner, factorization LOWER
    using Amat (not Pmat) as operator for blocks
    Preconditioner for the Schur complement formed from Sp, an assembled
        approximation to S, which uses A00's diagonal's inverse
    Split info:
    Split number 0 Defined by IS
    Split number 1 Defined by IS
    KSP solver for A00 block
      KSP Object: (fieldsplit_displacement_) 1 MPI process
        type: preonly
      PC Object: (fieldsplit_displacement_) 1 MPI process
        type: gamg
```

A00 is the 222 × 222 displacement block, A11 is the 3 × 3 Lagrange block.
Both PyLith-style sub-PCs are GAMG (no Trilinos ml available in the
fsrm-ci:main image). Subkeys `fieldsplit_displacement_*` and
`fieldsplit_lagrange_multiplier_fault_*` are correctly derived from the
field names.

## Verification 3 -- KSP convergence

Every Newton step's GMRES converges in 3 iterations with `CONVERGED_RTOL`:

```
0 SNES Function norm 6.250000000000e-05
  0 KSP Residual norm 1.371428571429e+08
  1 KSP Residual norm 1.411589936782e+07
  2 KSP Residual norm 5.680894673220e+05
  3 KSP Residual norm 4.271702224756e-04
  Linear solve converged due to CONVERGED_RTOL iterations 3
```

This is the qualitative breakthrough vs Session 23. Session 23 experimental
GMRES plateaued at 1e+9 with recursion drift because selfp's
`B inv(diag(A_disp)) B^T` was singular over the 25 unconstrained interior
Lagrange vertices. With rim-pin ON, only 1 Lagrange vertex (3 DOFs) survives
and selfp is non-singular over that tiny block.

Displacement sub-block GAMG (`/tmp/s24_prescribed.log:107-198`):

```
type: gamg
    type is MULTIPLICATIVE, levels=2 cycles=v
    Complexity:    grid = 1.03604    operator = 1.00983
        #equations  | nnz/row op | nnz/row int
              8     |     8      |     0
            222     |    27      |     3
    eigenvalue targets used: min 0.21787, max 2.39657
    eigenvalues provided (min 0.0906, max 2.179) with transform: [0. 0.1; 0. 1.1]
```

Two-level GAMG, smoothed prolongator eigenvalues 0.09 .. 2.18 -- well-conditioned.
The displacement block is now genuinely well-solved.

Lagrange sub-block (3 × 3 Schur complement) is preconditioned by GAMG over
selfp; setup completes without "Zero detected in diagonal" warnings (Session
23 had two such warnings).

## Verification 4 -- SNES trace and max_fault_slip

```
0 TS dt 1. time 0.
    0 SNES Function norm 6.250000000000e-05
        [KSP CONVERGED_RTOL 3 iters]
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    (line search rejected the step)

    0 SNES Function norm 6.250000000000e-05
        [KSP CONVERGED_RTOL 3 iters]
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0

    0 SNES Function norm 1.385637037625e+01    (after timestep cut/retry)
        [KSP CONVERGED_RTOL 3 iters]
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0

    0 SNES Function norm 1.385644254504e+01
        [KSP CONVERGED_RTOL 3 iters]
    1 SNES Function norm 5.750669116518e+00    (line search took partial step)

    0 SNES Function norm 8.000000000244e+00    (subsequent retries)
    1 SNES Function norm 7.944943148909e+00    (line search took partial step)
    [DIVERGED_LINE_SEARCH iterations 1; repeats 7 times for 7 cuts]
```

After several TS step cuts, SNES residual locks near 8.0 with line search
making essentially no progress per Newton step. TSAdaptCheckStage gives up
after the maximum number of step rejections.

Test result:
```
/workspace/tests/integration/test_dynamic_rupture_solve.cpp:455: Failure
    Expected: (sol_norm) > (0.0), actual: 0 vs 0
    Solution must be nonzero
/workspace/tests/integration/test_dynamic_rupture_solve.cpp:457: Failure
    Expected: (max_fault_slip) > (5.0e-4), actual: 0 vs 0.0005
    Prescribed slip must create a measurable displacement jump
```

`max_fault_slip = 0.0`. Failure mode: SNES line search rejects, no solution
update lands.

The structural problem: KSP returns a Newton step that is correct in
direction but enormous in magnitude (drives the function norm from 6.25e-05
to ~5.75e+00 -- five orders of magnitude overshoot). The line search
backtracks but cannot find a step that decreases the residual, because the
gradient direction is dominated by the 1 surviving interior Lagrange DOF
trying to enforce the prescribed slip target on a single vertex while the
neighboring 24 rim-pinned vertices stay at zero -- the resulting jump
discontinuity is non-physical and the bulk residual blows up.

## Fault subset table

| #   | Test                                              | S22 default | S23 default | S24 default | S24 experimental |
|-----|---------------------------------------------------|-------------|-------------|-------------|------------------|
|  40 | Physics.SCEC.TPV5                                 | FAIL        | FAIL        | FAIL        | FAIL             |
|  57 | Physics.LockedFaultTransparency                   | FAIL        | FAIL        | FAIL        | FAIL             |
|  59 | Physics.CohesiveBdResidual                        | PASS        | PASS        | PASS        | FAIL             |
|  77 | Functional.DynamicRuptureSetup.MeshSplitting      | PASS        | PASS        | PASS        | PASS             |
|  78 | Functional.DynamicRuptureSetup.LockedFault        | PASS        | PASS        | PASS        | PASS             |
|  79 | Functional.DynamicRuptureSetup.FaultGeometryFromConfig | PASS   | PASS        | PASS        | PASS             |
|  80 | Functional.DynamicRuptureSetup.SlippingFault      | PASS        | PASS        | PASS        | PASS             |
|  81 | Functional.DynamicRuptureSetup.AbsorbingCoexist   | PASS        | PASS        | PASS        | PASS             |
|  91 | Integration.PressurizedFractureFEM                | FAIL        | FAIL        | FAIL        | FAIL             |
|  93 | Integration.DynamicRuptureSolve.LockedQuasiStatic | FAIL        | PASS        | PASS        | FAIL             |
|  94 | Integration.DynamicRuptureSolve.LockedElastodynamic | PASS      | PASS        | PASS        | FAIL             |
|  95 | Integration.DynamicRuptureSolve.PrescribedSlip    | FAIL        | FAIL        | FAIL        | FAIL             |
|  96 | Integration.ExplosionFaultReactivation            | PASS        | PASS        | PASS        | FAIL             |
|  98 | Integration.TimeDependentSlip                     | FAIL        | PASS        | PASS        | FAIL             |
| 100 | Integration.SlippingFaultSolve                    | FAIL        | FAIL        | FAIL        | FAIL             |
| 101 | Integration.SlipWeakeningFault                    | FAIL        | FAIL        | FAIL        | FAIL             |
| **Pass count** |                                        | 8/16        | **10/16**   | **10/16**   | 5/16             |

S24 default holds at 10/16, identical to S23. S24 experimental regresses by 5
beyond the S24 default failures (CohesiveBdResidual, LockedQuasiStatic,
LockedElastodynamic, ExplosionFaultReactivation, TimeDependentSlip). All five
new failures share the same SNES `DIVERGED_LINE_SEARCH` pattern: KSP converges,
the Newton step magnitude is too large, line search rejects.

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

10/16 pass. No regression vs Session 23.

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

110/116 pass. CLAUDE.md baseline was 108/116 (`LockedQuasiStatic` and
`TimeDependentSlip` recovered in S23 default and held in S24 default).
Total wall time 136.93 sec.

## Decision gate

- PrescribedSlip experimentally: **FAIL** (SNES `DIVERGED_LINE_SEARCH`,
  max_fault_slip = 0).
- KSP convergence: **YES** (CONVERGED_RTOL in 3 iterations every Newton step).
- Default path regression: **none** (10/16 fault subset, 110/116 full suite).

This matches the prompt's decision-gate row: "PrescribedSlip fails but KSP
now converges (max_fault_slip off-tolerance): tuning issue, propose Session 25."

Gates are NOT inverted. `FSRM_ENABLE_SLIP_AUX` /
`FSRM_DISABLE_RIM_PIN` / `FSRM_SADDLE_SOLVER` remain opt-in.

## Diagnosis -- why the Newton step overshoots

The rim-pin essential-constrains 24 of 25 interior Lagrange vertices' 3
components. The 1 surviving interior Lagrange vertex carries the entire
prescribed-slip constraint for the whole fault. KSP solves the resulting
3-equation Schur system to high accuracy, but the step it returns is the
amount of slip that vertex would need to take to satisfy the (-0.001, 0, 0)
target -- which is the full slip in one Newton iteration starting from u=0.
This Newton step is also incompatible with the displacement field on the
adjacent rim-pinned vertices (which are required to satisfy lambda = 0
there), creating a spurious jump discontinuity at the rim. The bulk residual
explodes by 5 orders of magnitude when the step is taken, so the line
search rejects.

This is a structural mismatch: an aggressive interior-pin that keeps only
1 free Lagrange vertex is geometrically incompatible with a smooth
prescribed-slip field across the whole 4 × 4 fault face. PyLith gets away
with this on its `twoblocks` test because the fault has only a few interior
nodes anyway, and because PyLith uses a non-zero initial guess that is
physically close to the prescribed slip rather than starting from u=0.

## Proposed Session 25 tuning

Three avenues, in increasing order of structural change:

1. **Damp the line search.** Switch SNES line search to `bt` with
   `damping=0.1` and increase `snes_max_funcs`/`snes_max_it`. The aim is
   to take 10 small steps of 10% magnitude instead of one full step.
   Single-flag change:
   ```
   -snes_linesearch_type bt
   -snes_linesearch_damping 0.1
   -snes_max_it 50
   ```
   Likely to converge but slowly.

2. **Initialize the prescribed-slip target into the Lagrange field.** Set
   the 1 free interior Lagrange DOF to the target slip in
   `setInitialConditions` (or equivalently in `applyInitialFaultStress`).
   The first Newton residual then starts close to zero on that DOF and the
   step magnitude is small. This requires a small extension to
   `setInitialConditions` to write Lagrange components based on the
   prescribed-slip Cartesian vector for the 1 free vertex, which can be
   identified by intersecting `cohesive_interface` with NOT `buried_cohesive`.

3. **Coarsen the rim-pin to leave more interior Lagrange DOFs free.**
   The current geometric fallback labels all `seg_prism` candidates,
   which over-constrains the interior. PyLith's actual behavior in the
   `twoblocks` test labels only the rim *vertices* (point_prism), not the
   rim *edges* (seg_prism). Switch the geometric fallback to label only
   point_prism candidates and leave seg_prism free. This keeps more
   interior Lagrange DOFs available so the prescribed slip can be enforced
   smoothly across the fault.

Recommend tuning (3) as the structurally correct fix. Tunings (1) and (2)
are diagnostic complements: if (3) regresses the existing 10/16 fault
subset (e.g., LockedQuasiStatic re-fails because the relaxed rim-pin
re-introduces the original singularity), fall back to (1) + (2).

## Logs captured

- `/tmp/s24_build.log` -- build (incremental, no rebuild needed)
- `/tmp/s24_prescribed.log` -- experimental PrescribedSlip ctest with
  FSRM_KSP_VIEW (4913 lines, full GMRES traces and KSP views per Newton step)
- `/tmp/s24_fault.log` -- fault subset under FSRM_ENABLE_SLIP_AUX=1
  FSRM_SADDLE_SOLVER=fieldsplit (5/16 pass)
- `/tmp/s24_fault_default.log` -- fault subset, default env (10/16 pass)
- `/tmp/s24_test.log` -- full ctest suite, default env (110/116 pass)
