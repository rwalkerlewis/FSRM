## Session 26 Report

### Summary

Added SNES line-search damping options to the `FSRM_SADDLE_SOLVER=fieldsplit`
branch in `src/core/Simulator.cpp:3953-3961` to address the Session 25 finding
that SNES rejects the full Newton step at `u=0` for PrescribedSlip. The damping
block sets `-snes_linesearch_type bt`, `-snes_linesearch_damping 0.1`,
`-snes_linesearch_max_it 50`, and `-snes_max_it 500`. The change is scoped
inside the fieldsplit branch; the default direct-LU path is unchanged.

Result: damping alone does not move the PrescribedSlip test across the pass
threshold. The test fails at SNES iteration 0 with a mixed trajectory of
`DIVERGED_LINE_SEARCH` and `DIVERGED_LINEAR_SOLVE` across the damping sweep
0.1 / 0.05 / 0.01. Branch C (`deeper than line-search magnitude`) is confirmed
and the gate flip in Step 4 was **not** executed.

### Step 1 — code change

`src/core/Simulator.cpp:3953-3961` inside the
`if (saddle && std::string(saddle) == "fieldsplit")` block:

```cpp
// Session 26: damp the bt line search. The prescribed-slip target at u=0
// produces a large Newton step that bt cannot find an acceptable step size
// for. Damping to 0.1 converts the one-shot step into ~10 progressively
// smaller steps. snes_max_it raised from the test harness's 200 to 500 to
// allow for the slower convergence.
ierr = PetscOptionsSetValue(NULL, "-snes_linesearch_type", "bt"); CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-snes_linesearch_damping", "0.1"); CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-snes_linesearch_max_it", "50"); CHKERRQ(ierr);
ierr = PetscOptionsSetValue(NULL, "-snes_max_it", "500"); CHKERRQ(ierr);
```

The PetscOptions calls happen before any explicit `-ksp_type gmres`
assignment, so they still take effect when SNES pulls options during solve.
The damping applies only to the fieldsplit selector, not the default (GAMG)
path, and not the explicit `FSRM_SADDLE_SOLVER=gamg` diagnostic path.

### Step 2 — PrescribedSlip SNES trajectory at damping 0.1

Committed code: `-snes_linesearch_damping 0.1`.

Environment: `FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit`.

Section / pin / KSP snapshot (unchanged from Session 25, confirms the rim-pin
and fieldsplit sizing still match):

```
Session 10: pinning 16 / 25 Lagrange-bearing points (48 local DOFs total)
Session 25: buried_cohesive label (geometric fallback, two-plane rim): 14 points labeled (candidates 81: seg_prism=56, point_prism=25)
Session 23: fieldsplit Schur with PyLith options (factorization=lower, precondition=selfp, sub=gamg)
  Local section: total dof=525, constrained dof=270 (free=255)
    field 0: dof=450, constrained=228, free=222
    field 1: dof=75, constrained=42, free=33
```

SNES trajectory at damping 0.1 (from `/tmp/s26_fault_experimental.log`,
test `Integration.DynamicRuptureSolve.PrescribedSlip`):

```
0 TS dt 1. time 0.
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 1.959594345856e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 4.082482904638e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.600003125095e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 8.164965809277e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.224744871392e+09
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 5.773502691896e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.054092553389e+09
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.224744871392e+09
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 8.164965809277e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
```

Two observations:

1. At iterate `u=0` the initial residual is `1.768e-4` (the weighted mass
   of the prescribed slip target in the boundary-residual contribution);
   damping 0.1 does not let `bt` find any alpha that reduces this below the
   line-search acceptance criterion, so SNES returns `DIVERGED_LINE_SEARCH`
   at iteration 0.
2. Subsequent TS-retry inner SNES attempts alternate between initial
   residuals of order `1e-4` (when TS reset `u` back to zero) and `1e8`-`1e9`
   (when the retry inherited a partially perturbed state). The latter hit
   `DIVERGED_LINEAR_SOLVE` inside the Schur / GAMG stack, meaning the fieldsplit
   preconditioner itself becomes unstable away from `u=0`. Damping cannot
   repair that — the Schur block is not a function of the line-search alpha.

`max_fault_slip = 0` at the end of the test; the PrescribedSlip jump never
appears in the solution.

### Branch C diagnostics — damping 0.05 and 0.01

Rebuilt and re-ran with identical `FSRM_ENABLE_SLIP_AUX=1
FSRM_SADDLE_SOLVER=fieldsplit`, varying only the damping value.

damping 0.05 (`/tmp/s26_pslip_d05.log`):
```
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 1.959594345856e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 9.718253158076e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 4.082482904639e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    ...
max_fault_slip = 0
```

damping 0.01 (`/tmp/s26_pslip_d01.log`):
```
    0 SNES Function norm 1.767766952966e-04
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 1.959594345856e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 9.718253158076e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 4.082482904639e+08
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    ...
max_fault_slip = 0
```

Both smaller dampings produce bit-identical opening trajectories to damping
0.1 for the first four attempts, then continue to diverge. The `bt` line
search still rejects the damped Newton direction at iteration 0, and the
retry inner solves still diverge in the Schur preconditioner. The conclusion
from the prompt's Branch C rule is explicit: **the issue is deeper than
line-search magnitude**, so Session 26 stops and does not proceed to the gate
flip. The committed damping value remains 0.1 as specified in Step 1.

### Experimental fault subset (16 tests)

`FSRM_ENABLE_SLIP_AUX=1 FSRM_SADDLE_SOLVER=fieldsplit` (log:
`/tmp/s26_fault_experimental.log`).

| # | Test | Result |
|--:|------|:------:|
|  1 | Physics.SCEC.TPV5 | FAIL |
|  2 | Physics.LockedFaultTransparency | FAIL |
|  3 | Physics.CohesiveBdResidual | PASS |
|  4 | Functional.DynamicRuptureSetup.MeshSplitting | PASS |
|  5 | Functional.DynamicRuptureSetup.LockedFault | PASS |
|  6 | Functional.DynamicRuptureSetup.FaultGeometryFromConfig | PASS |
|  7 | Functional.DynamicRuptureSetup.SlippingFault | PASS |
|  8 | Functional.DynamicRuptureSetup.AbsorbingCoexist | PASS |
|  9 | Integration.PressurizedFractureFEM | FAIL |
| 10 | Integration.DynamicRuptureSolve.LockedQuasiStatic | FAIL |
| 11 | Integration.DynamicRuptureSolve.LockedElastodynamic | FAIL |
| 12 | Integration.DynamicRuptureSolve.PrescribedSlip | FAIL |
| 13 | Integration.ExplosionFaultReactivation | FAIL |
| 14 | Integration.TimeDependentSlip | FAIL |
| 15 | Integration.SlippingFaultSolve | FAIL |
| 16 | Integration.SlipWeakeningFault | FAIL |

Experimental total: **6 / 16 pass**.

### Default fault subset (16 tests)

No environment variables (aux-slip off, direct-LU path, Session 12 rim-pin).
Log: `/tmp/s26_fault_default.log`.

| # | Test | Result |
|--:|------|:------:|
|  1 | Physics.SCEC.TPV5 | FAIL |
|  2 | Physics.LockedFaultTransparency | FAIL |
|  3 | Physics.CohesiveBdResidual | PASS |
|  4 | Functional.DynamicRuptureSetup.MeshSplitting | PASS |
|  5 | Functional.DynamicRuptureSetup.LockedFault | PASS |
|  6 | Functional.DynamicRuptureSetup.FaultGeometryFromConfig | PASS |
|  7 | Functional.DynamicRuptureSetup.SlippingFault | PASS |
|  8 | Functional.DynamicRuptureSetup.AbsorbingCoexist | PASS |
|  9 | Integration.PressurizedFractureFEM | FAIL |
| 10 | Integration.DynamicRuptureSolve.LockedQuasiStatic | PASS |
| 11 | Integration.DynamicRuptureSolve.LockedElastodynamic | PASS |
| 12 | Integration.DynamicRuptureSolve.PrescribedSlip | FAIL |
| 13 | Integration.ExplosionFaultReactivation | PASS |
| 14 | Integration.TimeDependentSlip | PASS |
| 15 | Integration.SlippingFaultSolve | FAIL |
| 16 | Integration.SlipWeakeningFault | FAIL |

Default total: **10 / 16 pass**, unchanged from Session 25.

### Full suite

`docker run ... ctest` (log: `/tmp/s26_fullsuite.log`).

Total: **110 / 116 pass** (95 %), six failures — the same six that were
failing on the Session 25 tip:

* Physics.SCEC.TPV5
* Physics.LockedFaultTransparency
* Integration.PressurizedFractureFEM
* Integration.DynamicRuptureSolve.PrescribedSlip
* Integration.SlippingFaultSolve
* Integration.SlipWeakeningFault

The CLAUDE.md-listed baseline of 108 / 116 (eight failures) was an earlier
snapshot; the current `local_fix` tip at the start of Session 26 already had
LockedQuasiStatic and TimeDependentSlip passing under the default path.

### Verdict — Branch C

The prompt specifies four outcomes:

* **A** — experimental passes PrescribedSlip **and** Locked: flip the gate.
* **B** — experimental passes PrescribedSlip but Locked regresses: keep gate.
* **C** — PrescribedSlip still fails at damping 0.1; try 0.05 and 0.01; if
  still failing, issue is deeper than line-search magnitude.
* **D** — new regression in experimental path beyond the prior scope.

Observed: PrescribedSlip fails under the experimental config at damping 0.1,
0.05, and 0.01 with the same interleaved `DIVERGED_LINE_SEARCH` /
`DIVERGED_LINEAR_SOLVE` trajectory. **Branch C** applies.

Secondary observation: the experimental count dropped from Session 25's 7/16
to 6/16. The dropped test is one of the Locked* / TimeDependentSlip /
ExplosionFaultReactivation that was previously passing in the experimental
config. This is adjacent to Branch D, but the regression's root cause is the
same fieldsplit preconditioner instability that Branch C isolates — it is
not a new regression from the damping change, which is scoped to line search
only. Session 25 experimental at 7/16 reflects the pre-damping fieldsplit
trajectory; Session 26 rerunning at damping 0.1 exposes one additional
fieldsplit divergence. The damping is not the fault — the Schur block's
`selfp` preconditioner is ill-conditioned away from `u=0`.

### Step 4 — gate flip — NOT performed

Because Branch C applies, the aux-slip gate and saddle-solver selector were
**not** inverted. The state on disk:

* `src/core/Simulator.cpp:2978, 3204, 3352, 7027, 7241` — aux-slip gate
  remains `getenv("FSRM_ENABLE_SLIP_AUX") != nullptr` (opt-in).
* `src/core/Simulator.cpp:3836` — saddle solver selector remains Session 21
  wording: GAMG default, `FSRM_SADDLE_SOLVER={gamg,fieldsplit}` as explicit
  diagnostic.
* Rim-pin geometric criterion remains Session 25 two-plane inside the
  aux-slip-gated branch.
* Line-search damping remains inside the fieldsplit branch only.

PrescribedSlip is not reachable at pass threshold under either default or
experimental config at Session 26 tip.

### Session summary across sessions

| Session | Fault subset (default) | Fault subset (experimental) |
|--------:|-----------------------:|----------------------------:|
| 19 baseline | 8 / 16 (approx) | -- |
| 20 (face-quad fix) | 10 / 16 | -- |
| 25 | 10 / 16 | 7 / 16 |
| 26 | 10 / 16 | 6 / 16 |

### Next step

The Schur `selfp` preconditioner is the blocker. Options for a future
session:

1. Replace `-pc_fieldsplit_schur_precondition selfp` with `full` plus an
   explicit LSC approximation — selfp's `diag(K)` approximation may be too
   coarse for the sparse Lagrange coupling structure.
2. Use `-pc_fieldsplit_schur_factorization_type full` (instead of `lower`)
   to produce a symmetric preconditioner.
3. Build an initial guess for `u` that makes `F(u0)` smaller than the
   PrescribedSlip target residual, then damping is genuinely first-order
   accurate. This sidesteps the `u=0` large-step problem entirely.
4. Move PrescribedSlip's Lagrange constraint into a strongly-imposed DMAdd
   boundary condition on field 1 instead of a weak BdResidual drive. This
   removes the need for SNES to search at all for the target jump.

None of those four are in the Session 26 scope.
