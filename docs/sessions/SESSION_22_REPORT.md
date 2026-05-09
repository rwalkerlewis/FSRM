# Session 22 Report

## Goal

Replace Session 21's incorrect GAMG smoother (`-mg_levels_pc_type jacobi`, scalar)
with PyLith's exact set (`-mg_fine_pc_type vpbjacobi`, variable point-block) and
re-run the aux-slip + rim-pin-off experiment on `DynamicRuptureSolve.PrescribedSlip`.

## Change in tree

`src/core/Simulator.cpp:3913-3928` (branch `FSRM_SADDLE_SOLVER=gamg`):

- REMOVED: `-mg_levels_ksp_type richardson`, `-mg_levels_pc_type jacobi`,
  `-mg_levels_ksp_max_it 5` (Session 21 scalar-jacobi smoother on every MG
  level).
- ADDED: `-mg_fine_ksp_max_it 5`, `-mg_fine_pc_type vpbjacobi` (PyLith
  `libsrc/pylith/materials/Elasticity.cc:309-316` set), leaving the coarse
  levels on the GAMG default (chebyshev / sor).

The explicit experimental path now matches the default (`!pc_already_set`)
branch byte-for-byte. The printed banner is updated to
`Session 22: FSRM_SADDLE_SOLVER=gamg with PyLith options (vpbjacobi on fine
level, GAMG defaults on coarse)`.

Test harness `runFullPipeline` (`tests/integration/test_dynamic_rupture_solve.cpp:295-315`)
calls `PetscOptionsClear` then installs its `-pc_type=lu -ksp_type=preonly`
defaults before invoking `Simulator::initializeFromConfigFile` and the rest of
the pipeline. `setupSolvers` fires later and the `FSRM_SADDLE_SOLVER=gamg`
branch overrides `-pc_type` back to `gamg`, confirmed in both experimental
logs by the banner `Session 22: FSRM_SADDLE_SOLVER=gamg ...` appearing before
the first SNES Function norm. Ordering is correct.

## Experimental PrescribedSlip run

Command:

```bash
FSRM_ENABLE_SLIP_AUX=1 FSRM_DISABLE_RIM_PIN=1 FSRM_SADDLE_SOLVER=gamg \
  FSRM_KSP_VIEW=1 \
  ctest -R DynamicRuptureSolve.PrescribedSlip
```

Log: `/tmp/s22_prescribed_experimental.log`, `/tmp/s22_gamg_verbose.log`.

### SNES residual history (experimental)

```
0 TS dt 1. time 0.
    0 SNES Function norm 2.189978754336e-04   <-- same signature as Session 21
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 2.189978754336e-04   (TS step rejected, retry)
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 6.666666666667e+08   (residual explodes on further retry)
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 4.082482904638e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    ... (oscillating 3.7e+08, 7.9e+08, 1.15e+09) ...
```

Final: `sol_norm = 0`, `max_fault_slip = 0`. TSSolve returned `ierr=1`; the
gtest assertions at lines 455 and 457 fired on zero displacement and zero
slip.

### KSP convergence (experimental)

KSP runs **zero iterations** then returns `DIVERGED_LINEAR_SOLVE`. There is no
per-iteration `-ksp_monitor` output because the PC setup succeeds but the
first residual reduction fails outright. `-info :pc` output (`/tmp/s22_gamg_verbose.log`):

```
<pc:gamg> PCSetUp_GAMG: level 0) N=297, n data cols=6, nnz/row=32
<pc:gamg> PCGAMGOptimizeProlongator_AGG: Smooth P0 level 0: emax=2.158573e+00 emin=1.602106e-12 PC=jacobi
<pc:jacobi> PCSetUp_Jacobi: Zero detected in diagonal of matrix, using 1 at those locations
<pc:gamg> PCSetUp_GAMG: 1) N=60, n data cols=6, nnz/row=49
<pc:gamg> PCGAMGOptimizeProlongator_AGG: Smooth P0 level 1: emax=2.688970e+11 emin=1.036947e-05
<pc:gamg> PCSetUp_GAMG: 2) N=12, 3 levels, operator complexity=1.32926
<pc:gamg> PCSetUp_MG: Using outer operators to define finest grid operator
  because PCMGGetSmoother(pc,nlevels-1,&ksp);KSPSetOperators(ksp,...); was not called.
<pc:jacobi> PCSetUp_Jacobi: Zero detected in diagonal of matrix, using 1 at those locations
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
```

Key observations:

1. `emax/emin = 2.69e+11 / 1.04e-05 ~ 2.6e+16` on the level-1 (coarse) matrix.
   The coarse operator is catastrophically ill-conditioned.
2. `PCSetUp_Jacobi: Zero detected in diagonal` fires both for the prolongator
   smoother and on the finest level, confirming null rows/columns in the
   coarse Lagrange block.
3. Removing the rim pin (Session 10 block) pushes the zero-diagonal pattern
   into the coarse operator, where neither `vpbjacobi` nor the GAMG default
   Chebyshev+SOR can condition it. The linear solve terminates before any
   iteration count is reported.

### Decision-gate signature match

Session 21 experimental: SNES residual went `2.19e-4 -> 17.88`
(DIVERGED_LINE_SEARCH at iteration 0 after one KSP step, line-search backed
the update out to origin).

Session 22 experimental: SNES residual goes
`2.19e-4 -> 2.19e-4 -> 6.67e+08`, where each intermediate is a TS step
retry and KSP does zero iterations before `DIVERGED_LINEAR_SOLVE`. The
failure mode shifted from "KSP returns garbage step, line-search rejects"
to "KSP refuses to iterate at all", but the initial residual (2.19e-4) and
the final verdict (SNES unable to reduce the residual, test-level
`max_fault_slip = 0`) are unchanged. This falls under the "PrescribedSlip
still fails with same 2.19e-4 signature -> structural issue" decision
gate in the Session 22 prompt.

## Fault subset: Session 21 default vs Session 22 default vs Session 22 experimental

All 16 tests in the fault-focused regex
`Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup`.

| Test | S21 default | S22 default | S22 experimental |
|------|:-----------:|:-----------:|:----------------:|
| Physics.SCEC.TPV5 | FAIL | FAIL | FAIL |
| Physics.LockedFaultTransparency | FAIL | FAIL | FAIL |
| Physics.CohesiveBdResidual | PASS | PASS | FAIL |
| Integration.PressurizedFractureFEM | FAIL | FAIL | FAIL |
| Integration.DynamicRuptureSolve.LockedQuasiStatic | PASS | PASS | FAIL |
| Integration.DynamicRuptureSolve.LockedElastodynamic | PASS | PASS | FAIL |
| Integration.DynamicRuptureSolve.PrescribedSlip | FAIL | FAIL | FAIL |
| Integration.ExplosionFaultReactivation | PASS | PASS | FAIL |
| Integration.TimeDependentSlip | PASS | PASS | FAIL |
| Integration.SlippingFaultSolve | FAIL | FAIL | FAIL |
| Integration.SlipWeakeningFault | FAIL | FAIL | FAIL |
| Functional.DynamicRuptureSetup (5 configs) | PASS x5 | PASS x5 | PASS x5 |
| Pass total | 10/16 | 10/16 | 5/16 |

Notes:

- S21 default column reflects the baseline after Session 21's GAMG +
  near-null-space activation, not the pre-Session-21 CLAUDE.md header
  (which is now partially stale; the baseline list of 8 known failures
  is already down to 6 as recorded in `/tmp/s22_test.log`).
- S22 default column is identical to S21 default. The Session 22 change
  does not touch the `!pc_already_set` branch.
- S22 experimental is worse (5 passes) because disabling the rim pin
  removes the only regularization that kept the coarse Lagrange block
  non-singular; the coarser levels GAMG constructs for LockedQuasiStatic,
  LockedElastodynamic, ExplosionFaultReactivation, TimeDependentSlip, and
  CohesiveBdResidual all suffer the same zero-diagonal / ill-conditioned
  coarse operator.

## Default-path regressions vs Session 21

None. `/tmp/s22_fault_default.log` reports `10/16` fault-subset passes,
identical to Session 21's post-activation baseline.

## Full-suite count (default path, Session 22)

`/tmp/s22_test.log`: **110/116 tests passed** (95%). Six failures, all in the
fault-related subset listed above:

```
40  - Physics.SCEC.TPV5
57  - Physics.LockedFaultTransparency
91  - Integration.PressurizedFractureFEM
95  - Integration.DynamicRuptureSolve.PrescribedSlip
100 - Integration.SlippingFaultSolve
101 - Integration.SlipWeakeningFault
```

This matches Session 21's post-activation baseline. The pre-Session-21 header
in `CLAUDE.md` listing 8 known failures is stale; `LockedQuasiStatic` and
`TimeDependentSlip` have been passing since Session 21.

## Decision gate

Per the Session 22 prompt the outcome is: **"PrescribedSlip still fails
with same 2.19e-4 signature -> structural issue; stop and propose Session
23."**

Gate inversion is therefore **NOT** performed. The five call sites Session
21 identified (`setupAuxiliaryDM` material-aux trigger, slip aux DM
creation, slip aux attach, `FormFunction`/`FormJacobian` attach,
`setupBoundaryConditions` rim-pin) plus the `saddle` selector in
`setupSolvers` remain on the Session 21 gating scheme.

## Proposed Session 23

The smoother option is now PyLith's exact set, and the failure mode has
shifted from `DIVERGED_LINE_SEARCH` (Session 21) to
`DIVERGED_LINEAR_SOLVE at 0 iters` (Session 22), but the underlying
coarse-operator ill-conditioning (`emax/emin ~ 2.6e+16` on level 1,
`Zero detected in diagonal`) is structural, not smoother-level. The
problem is that FSRM's DM state / section ordering / residual assembly
does not reproduce exactly what PyLith feeds into the same solver
options. Two branches for Session 23, as the prompt anticipated:

### 23a: Comparative KSP trace against PyLith

Build a minimal PyLith prescribed-slip case on the same 4x4x4 cube
mesh with one fault plane and run PyLith with `-ksp_monitor
-ksp_view -info :pc`. Diff the GAMG coarsening output (level sizes,
eigenvalue estimates, prolongator statistics) against the FSRM
experimental trace captured in `/tmp/s22_gamg_verbose.log`. The
discrepancy will point to one of: (a) PyLith applies
`-dm_reorder_section cohesive` and FSRM does not, (b) PyLith's
cohesive residual zeroes the interior Lagrange rows before assembly
and FSRM does not, or (c) PyLith attaches a different near-null-space
(one that includes the Lagrange block).

### 23b: Skip to friction Jacobian port (Session 22b)

Accept that the solver path is blocked on something subtle about the
saddle-point matrix structure and port PyLith's semi-smooth Newton
friction Jacobian (off-diagonal slip-direction derivatives,
friction-normal coupling) directly. The current
`addCohesivePenaltyToJacobian` has a simplified Coulomb block that is
already known to be insufficient for SlippingFaultSolve /
SlipWeakeningFault convergence. Porting PyLith's block would
simultaneously address the two slipping-fault failures AND remove
FSRM's dependence on the exact penalty-scaled diagonal we currently
rely on, which may change the conditioning enough to unblock
PrescribedSlip.

Recommendation: 23a first (a one-day diagnostic with clear pass/fail -
either the PyLith trace differs structurally or it matches and we have
direct evidence that the problem is outside the smoother), then 23b
regardless of 23a's outcome because the friction Jacobian is wanted for
its own sake.

## Summary

- Smoother change in tree, confirmed by runtime banner.
- Experimental PrescribedSlip: same 2.19e-4 structural signature, now
  `DIVERGED_LINEAR_SOLVE at 0 iters` with catastrophic coarse
  conditioning (`emax/emin ~ 2.6e+16`).
- Default path unchanged: `10/16` fault subset, `110/116` full suite.
- Do not invert gates. Propose Session 23 per 23a / 23b above.
