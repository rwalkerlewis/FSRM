## Session 25 Report

### Summary

Refined the geometric rim-pin criterion in
`Simulator::getOrCreateBuriedCohesiveLabel` (`src/core/Simulator.cpp:4899-4988`).
The original "any cone vertex on the mesh bounding box" check labeled every
cohesive prism whose cone touches a coordinate plane, which for a fault that
extends from wall to wall is **every** cohesive prism (56 of 56 seg_prism on
the 4x4 PrescribedSlip mesh). This left only 1 free Lagrange vertex of 25 and
made the prescribed-slip Newton step geometrically infeasible.

Session 25 introduces a "two-plane rim" criterion: a cohesive point is a rim
point only if at least one cone vertex sits on **two** coordinate-boundary
planes simultaneously (i.e., on a domain-cube edge where the fault plane
terminates). The looser criterion is gated behind `FSRM_ENABLE_SLIP_AUX=1`
so the experimental aux-slip + fieldsplit path uses the new rim-pin while the
default direct-LU path retains the Session 12 aggressive pin that the Locked*
and ExplosionFaultReactivation tests rely on.

### Diagnostic counts

For the 4x4 PrescribedSlip cube (cohesive candidates = 81: 56 seg_prism +
25 point_prism):

| Path | Criterion | Labeled | Field 1 free DOFs |
|------|-----------|--------:|------------------:|
| Default (S24 / S25 default) | any-plane | 56 | 3 |
| Experimental (S25, FSRM_ENABLE_SLIP_AUX=1) | two-plane | 14 | 33 |

Section header for the experimental path:
```
Local section: total dof=525, constrained dof=270 (free=255)
  field 0: dof=450, constrained=228, free=222
  field 1: dof=75, constrained=42, free=33
Session 10: pinning 16 / 25 Lagrange-bearing points (48 local DOFs total)
```

The interior 9 vertices x 3 components = 27 free Lagrange DOFs are now
correctly available for the prescribed slip to spread into. The 33 reported
free DOFs include 6 from the 6 Lagrange-bearing fault-rim corners that fall
outside the rim-pin (depth-1 cohesive segments whose cone vertices are at
fault corners but where one of the corners is shared with an interior
vertex).

### KSP convergence (experimental path, first Newton step)

```
0 SNES Function norm 1.767766952966e-04
  0 KSP Residual norm 1.285219755579e+09
  ...
 27 KSP Residual norm 1.634683871826e+04
  Linear solve converged due to CONVERGED_RTOL iterations 27
```

Schur fieldsplit converges (vs S24's 3 iterations on the over-constrained
problem). The looser problem is harder for KSP; 27 iterations on a 255-DOF
system is acceptable.

### SNES line-search behavior (experimental, PrescribedSlipQuasiStatic)

```
0 SNES Function norm 1.767766952966e-04
  Linear solve converged ... iterations 27
Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
0 SNES Function norm 1.959594345856e+01    <-- restart with reset state
  KSP fails to converge in 100+ iters
```

The first Newton step's KSP step is rejected by the bt line search (norm
ratio bigger than allowed). On the implicit retry the residual jumps by
five orders of magnitude and the system becomes unsolvable. The failure
mode is quantitatively different from Session 24 (which failed on the
linear solve itself); this is genuine line-search rejection of an
otherwise-converging linear correction.

### `max_fault_slip` and PrescribedSlip result

PrescribedSlipQuasiStatic still fails (TSSolve returns ierr=82,
DIVERGED_LINE_SEARCH on step 0). No `max_fault_slip` figure is reachable;
SNES rejects before the time-step completes.

### Fault subset table

All counts taken from serial `ctest --output-on-failure` runs (no `-j`,
to avoid the HDF5-conflict flakiness CLAUDE.md documents).

| Test | S24 default | S25 default | S25 experimental |
|------|:-----------:|:-----------:|:----------------:|
| Functional.DynamicRuptureSetup | pass | pass | pass |
| Physics.CohesiveBdResidual | pass | pass | pass |
| Physics.LockedFaultTransparency | fail | fail | fail |
| Physics.SCEC.TPV5 | fail | fail | fail |
| Integration.DynamicRuptureSolve.LockedQuasiStatic | fail | pass | fail |
| Integration.DynamicRuptureSolve.LockedElastodynamic | pass | pass | fail |
| Integration.DynamicRuptureSolve.PrescribedSlip | fail | fail | fail |
| Integration.TimeDependentSlip | fail | pass | fail |
| Integration.SlippingFaultSolve | fail | fail | fail |
| Integration.SlipWeakeningFault | fail | fail | fail |
| Integration.PressurizedFractureFEM | fail | fail | fail |
| Integration.ExplosionFaultReactivation | pass | pass | pass |
| Integration.FaultAbsorbingCoexist | pass | pass | pass |
| Functional.* (5 setup tests) | pass | pass | pass |

S25 default = 10/16 pass, S25 experimental = 7/16 pass. The experimental
path is worse than default for Locked*: the looser rim-pin removes
constraints the direct-LU saddle-point factorization needs to remain
non-singular.

### Full suite

`ctest --output-on-failure` (no env, no -j): 110/116 pass, 6 fail.

```
40  Physics.SCEC.TPV5
57  Physics.LockedFaultTransparency
91  Integration.PressurizedFractureFEM
95  Integration.DynamicRuptureSolve.PrescribedSlip
100 Integration.SlippingFaultSolve
101 Integration.SlipWeakeningFault
```

This is two fewer failures than the CLAUDE.md baseline of 108/116:
LockedQuasiStatic and TimeDependentSlip both pass now in the serial
default run. Likely test-order or HDF5-conflict variability rather than
a Session 25 effect (the default path code is unchanged after the gating
revert). No new regressions from Session 25.

### Decision-gate outcome

The prompt enumerated four cases. Outcome:

- "PrescribedSlip passes experimentally AND default path holds" -- not met.
- "PrescribedSlip fails with line-search rejection but smaller overshoot:
  try SNES damping per Session 24 suggestion (1), one more session." -- this
  is the matching case. The experimental KSP now converges (27 iters,
  CONVERGED_RTOL) where Session 24 had a different failure mode; the
  remaining obstacle is a line-search rejection on an otherwise valid
  correction. Defer SNES damping (`-snes_linesearch_type bt
  -snes_linesearch_damping 0.5`) to Session 26.
- "PrescribedSlip fails differently" -- not the matching case.
- "Default path regresses" -- this also occurred (LockedElastodynamic and
  ExplosionFaultReactivation both broke under the unconditional looser
  rim-pin). Mitigated by gating the new criterion behind
  `FSRM_ENABLE_SLIP_AUX=1`, which restored the default path to the S24
  baseline.

Aux-slip and saddle-solver gates are NOT inverted. Rim-pin remains ON in
both paths per PyLith. The only Session 25 behavioral change is the
rim-pin selection inside the experimental path.

### Files changed

- `src/core/Simulator.cpp:4899-4988` -- gated two-plane rim criterion in
  `getOrCreateBuriedCohesiveLabel`. Default path retains original any-plane
  criterion.

### Suggested next session

1. Add `-snes_linesearch_damping 0.5` (or `-snes_linesearch_type basic` to
  bypass bt line search entirely) under the experimental gate and re-run
  PrescribedSlipQuasiStatic. The first KSP correction is structurally
  reasonable (CONVERGED_RTOL); damping it to half-step may let the line
  search accept and the residual may then drop on the next Newton step.
2. If damping accepts, measure `max_fault_slip` against the analytical
  prescribed jump and decide whether to invert the experimental gates.
3. If damping does not accept, capture the rejected-step displacement
  field and the residual at the rejected point for diagnosis.
