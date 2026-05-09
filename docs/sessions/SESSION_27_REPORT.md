# Session 27 Report -- Tighten KSP/SNES Tolerances to PyLith Values for Fieldsplit Path

## Summary

Replaced Session 26's line-search damping block in the
`FSRM_SADDLE_SOLVER=fieldsplit` branch of `Simulator::setupSolvers` with
PyLith's production KSP/SNES tolerance values
(`-ksp_rtol 1e-14`, `-ksp_atol 1e-7`, `-snes_rtol 1e-14`, `-snes_atol 5e-7`,
`-ksp_max_it 500`, `-snes_max_it 200`).

The tighter tolerances are now in effect on the fieldsplit path
(verified via the KSP residual history continuing well past the previous
`rtol=1e-5` cutoff at iter 27). PrescribedSlip still fails: KSP cannot
reduce the residual by the required ~14 orders of magnitude before
hitting `ksp_max_it`, so SNES diverges at iteration 0 with `ierr=82`
(`PETSC_ERR_NOT_CONVERGED`). Because experimental
`FSRM_SADDLE_SOLVER=fieldsplit` produces fewer fault-test passes than
the default GAMG path (7/16 vs 10/16 in the fault subset), the gate
remains opt-in.

The full default suite is at 110/116 (95% pass), better than the
CLAUDE.md baseline of 108/116. The two newly-passing tests are
`Integration.DynamicRuptureSolve.LockedQuasiStatic` and
`Integration.TimeDependentSlip`, both of which converge under the
default GAMG path that Session 27 did not modify -- the CLAUDE.md
baseline is stale relative to the live `local_fix` tip.

## 1. Verify Session 26's damping block removed

`src/core/Simulator.cpp` no longer contains
`-snes_linesearch_damping 0.1`, `-snes_linesearch_max_it 50`, or the
500-iteration `-snes_max_it` lift. A grep across the tree confirms:

```
$ grep -rn "snes_linesearch_damping" src/
(no matches)
$ grep -rn "Session 26" src/
(no matches)
```

## 2. Verify Session 27 tolerances landed

`src/core/Simulator.cpp:3956-3966` reads:

```cpp
// Session 27: tighten KSP tolerance to PyLith's production value.
ierr = PetscOptionsSetValue(NULL, "-ksp_rtol",  "1.0e-14"); CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-ksp_atol",  "1.0e-7");  CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-snes_rtol", "1.0e-14"); CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-snes_atol", "5.0e-7");  CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-ksp_max_it",  "500"); CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-snes_max_it", "200"); CHKERRQ(ierr);
```

The order-of-operations check passes: the test harness
(`tests/integration/test_dynamic_rupture_solve.cpp:303-314`) calls
`PetscOptionsClear` then sets `-snes_rtol 1e-6` etc. via
`PetscOptionsSetValue`, then calls `sim.setupSolvers()`
(line 354) which contains the new Session 27 calls -- so the
PyLith tolerances win on the fieldsplit branch.

## 3. KSP convergence for PrescribedSlip under the new tolerances

From `/tmp/s27_prescribed.log` (`FSRM_ENABLE_SLIP_AUX=1
FSRM_SADDLE_SOLVER=fieldsplit FSRM_KSP_VIEW=1`):

```
0 SNES Function norm 1.767766952966e-04
  0 KSP Residual norm 1.285219755579e+09
  ...
 27 KSP Residual norm 1.634683871826e+04   <-- PETSc default rtol=1e-5 stop
 28 KSP Residual norm 1.151886644786e+04
  ...
 99 KSP Residual norm 1.900172291873e+03
100 KSP Residual norm 1.623722041865e+04   <-- GMRES restart (restart=100)
  ...
199 KSP Residual norm 2.914606414431e+03
```

The new rtol=1e-14 keeps KSP iterating well past the old cutoff
(iter 27 -> at least iter 199), confirming the option is in effect.
Convergence rate is too slow to reach `1.29e+9 * 1e-14 = 1.29e-5`
absolute residual within the 500-iteration budget: the residual
drops from 1.29e+9 to ~2.9e+3 in 200 iterations (~5.6 orders of
magnitude), and the rate after the GMRES restart at iter 100 is
~0.005 orders/iter. Reaching the rtol target would require
roughly 1800 more iterations at the current rate.

KSP terminates without printing a `KSPConvergedReason` and SNES
records `DIVERGED_LINEAR_SOLVE`, returning `PETSC_ERR_NOT_CONVERGED`
(82) up through TSSolve.

## 4. SNES residual history

Only one SNES iteration: function norm `1.768e-4` at iter 0, then
`DIVERGED_LINEAR_SOLVE` aborts the Newton step. Compare with the
default GAMG path on the same test from `/tmp/s27_fault_default.log`
(test #95):

```
0 SNES Function norm 6.250000000000e-05    (initial pre-step)
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
0 SNES Function norm 1.385637037625e+01    (after first dt)
1 SNES Function norm 6.930778885126e-04
2 SNES Function norm 1.955090238982e-04
Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 2
```

The default path "converges" SNES but to the wrong solution
(`max_fault_slip = 263.0`, target 0.002). The experimental path
fails to take any Newton step at all.

## 5. `max_fault_slip` and test result

| Path | SNES outcome | `max_fault_slip` | Test |
| --- | --- | --- | --- |
| Default GAMG | CONVERGED (wrong soln) | 263.044 | FAIL (magnitude check) |
| Experimental fieldsplit | DIVERGED_LINEAR_SOLVE | n/a (TSSolve errored) | FAIL (`ierr=82`) |

Session 26 baseline (in `docs/SESSION_26_REPORT.md`) reported the same
`max_fault_slip = 263.044` under default and a SNES line-search failure
under fieldsplit at the previous tolerances. Session 27's tolerance
tightening shifted the experimental failure mode from
`DIVERGED_LINE_SEARCH` to `DIVERGED_LINEAR_SOLVE` but did not produce
nonzero progress.

## 6. Full fault subset experimental

`FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit ctest -R '...'`
in `/tmp/s27_fault_experimental.log`:

`44% tests passed, 9 tests failed out of 16`

Failing under experimental:
- 40 `Physics.SCEC.TPV5`
- 57 `Physics.LockedFaultTransparency`
- 91 `Integration.PressurizedFractureFEM`
- 93 `Integration.DynamicRuptureSolve.LockedQuasiStatic` *(regression vs default)*
- 94 `Integration.DynamicRuptureSolve.LockedElastodynamic` *(regression vs default)*
- 95 `Integration.DynamicRuptureSolve.PrescribedSlip`
- 98 `Integration.TimeDependentSlip` *(regression vs default)*
- 100 `Integration.SlippingFaultSolve`
- 101 `Integration.SlipWeakeningFault`

The Locked* and TimeDependentSlip failures all show
`DIVERGED_LINEAR_SOLVE iterations 0` repeated for many SNES re-tries.
The fieldsplit Schur path with PyLith options is still ill-conditioned
on the FSRM saddle structure; KSP cannot find a useful direction even
with the relaxed (default) tolerances those tests previously enjoyed.

## 7. Default regression check

`ctest -R '...'` (no env vars) in `/tmp/s27_fault_default.log`:

`63% tests passed, 6 tests failed out of 16`

Failing under default:
- 40 `Physics.SCEC.TPV5`
- 57 `Physics.LockedFaultTransparency`
- 91 `Integration.PressurizedFractureFEM`
- 95 `Integration.DynamicRuptureSolve.PrescribedSlip`
- 100 `Integration.SlippingFaultSolve`
- 101 `Integration.SlipWeakeningFault`

The default path is unaffected by Session 27 -- the new tolerances
only land inside the `FSRM_SADDLE_SOLVER == "fieldsplit"` branch.
Session 26's default subset count is unchanged; the
`LockedQuasiStatic`, `LockedElastodynamic`, and `TimeDependentSlip`
tests pass here.

## 8. Full suite count

`ctest --output-on-failure` (no env vars) in `/tmp/s27_test.log`:

`95% tests passed, 6 tests failed out of 116`

Same six failures as the default fault subset above. This is two
better than the CLAUDE.md baseline (108/116), confirming that
`Integration.DynamicRuptureSolve.LockedQuasiStatic` and
`Integration.TimeDependentSlip` pass on the live `local_fix` tip
even though CLAUDE.md still lists them as known failures. The
CLAUDE.md baseline list dates from before sessions 22--26 and is
overdue for refresh.

## 9. Gate decision

Per the Session 27 prompt:

> Experimental count < default count: keep opt-in. PrescribedSlip is
> reachable via explicit env vars.

Experimental count (7/16) is less than default count (10/16), so
the gate stays opt-in. PrescribedSlip is *not* reachable -- it
fails under both default (wrong magnitude) and experimental (KSP
does not converge). Session 27 has not unlocked it.

The KSP rate-of-convergence problem is structural: 200 GMRES
iterations only buy ~6 orders of magnitude of residual reduction
on this saddle-point system, but PyLith's `rtol=1e-14` requires 14.
Either the preconditioner is wrong (Schur+selfp is not effective
on the FSRM Lagrange block) or the assembled saddle matrix is much
more ill-conditioned than PyLith's. A path forward likely requires
either:

1. A different Schur-complement preconditioner (e.g., explicit
   `M B^{-1} A B^{-T}` approximation rather than `selfp`).
2. Loosening the rim-pin so the Lagrange block is full-rank, then
   GAMG-on-the-whole-system without fieldsplit.
3. A Jacobian audit (`FSRM_SNES_TEST_JAC=1` /
   `FSRM_DUMP_JAC=...`) to verify the assembled matrix matches
   the FD Jacobian -- per prompt, "do not attempt further option
   tuning without diagnostic data."

The Session 27 change is preserved on the fieldsplit branch as a
diagnostic baseline: the rtol=1e-14 history is more useful for
subsequent investigation than the old default rtol=1e-5 cutoff at
KSP iter 27.

## Files touched

- `src/core/Simulator.cpp` -- replaced Session 26 line-search block
  (4 `-snes_linesearch_*`/`-snes_max_it` options) with Session 27
  PyLith tolerance block (6 options).
- `docs/SESSION_27_REPORT.md` -- this file.

CLAUDE.md is **not** updated (Rule 14: only on gate pass).

## Logs

- `/tmp/s27_build.log` -- build (success).
- `/tmp/s27_prescribed.log` -- experimental PrescribedSlip with KSP view.
- `/tmp/s27_fault_experimental.log` -- experimental fault subset (7/16).
- `/tmp/s27_fault_default.log` -- default fault subset (10/16).
- `/tmp/s27_test.log` -- full default suite (110/116).
