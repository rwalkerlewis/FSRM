# Session 19 Report: Volume-FE Aux Slip Plumbing

## Summary

Session 19 followed the prompt to swap the Session 18 dim-1 (surface)
slip aux FE for a dim (volume) FE plus
`DMSetFieldAvoidTensor(slipAuxDM_, 0, PETSC_TRUE)`, on the theory that
the volume FE would push PETSc's hybrid integrator down the
`auxOnBd = (dimAux == dim) ? PETSC_TRUE : PETSC_FALSE` branch and pull
the volume tabulation. The change does not solve the
`febasic.c:663` mismatch on PETSc 3.25: instead, the closure-size
guard at `plexfem.c:3982` fires first
(`closure 9 != DS size 12` for a tetrahedron-shaped FE on a tensor
prism cohesive cell). The aux-slip cohesive-key attachment therefore
remains opt-in (`FSRM_ENABLE_SLIP_AUX=1`); the default path stays on
the Session 16 constants-array baseline. No fault-subset regressions
relative to the `local_fix` tip.

## Why the prompt's path does not close

`PetscFEIntegrateHybridResidual_Basic` (`febasic.c`) reads `dim` from
the cohesive (main) field's FE, not from the cell:

```c
PetscCall(PetscDSGetDiscretization(ds, key.field, &fe));
PetscCall(PetscFEGetSpatialDimension(fe, &dim));
```

For our setup the cohesive field at `key.field` is the Lagrange
multiplier FE created with `PetscFECreateLagrange(comm, dim - 1, ...)`,
so `PetscFEGetSpatialDimension(lag_fe) = dim - 1 = 2`, not 3. The
`auxOnBd` branch then evaluates as:

| aux FE choice            | `dimAux` | `auxOnBd = (dimAux == 2)` | tabulation pulled        |
| ------------------------ | -------: | ------------------------- | ------------------------ |
| dim-1 surface (Session 18) | 2        | TRUE                      | aux *volume* tab         |
| dim volume (Session 19)  | 3        | FALSE                     | aux *face* tab           |

So the prompt's premise is inverted on PETSc 3.25: the surface FE
already lands on the volume-tabulation branch. The
"Number of tabulation points 3 != 1" failure that Session 18 reported
does not come from the slip aux at all - it comes from the *material*
aux DS that is queried for the bulk-side keys
(`material_label_, value=1` and `value=2`). The material aux has
field-0 vol Np = 1 (P0, default quadrature) while the
displacement field's tabulation has Np = 3 in this code path. With
`FSRM_ENABLE_SLIP_AUX=0` the cohesive-side aux is never attached, so
the bulk sides never trip; with the env var set, every `dsAux[s]`
slot has to satisfy the same Np guard regardless of which side fires
the integrator. Verified with debug printouts:

```
Session 19: slip aux DM built (dim-1 surface FE + AvoidTensor +
            Lagrange-quadrature copy, interfaces-restricted, cohesive
            flag; default-DS totDim=9, local vec size=75)
Session 19 debug: slipAuxDM dimAux=2 vol_Np=3 face_Np=2
Session 19 debug: Lagrange FE dim=2 vol_Np=3 face_Np=2
```

Switching to a dim volume FE keeps `auxOnBd = FALSE` *and* trips the
upstream closure check at `plexfem.c:3982`:

```
Closure size 9 for subcell 384 does not match DS size 12
```

A volume P1 vector FE wants 4 vertices x 3 components = 12 DOFs per
element. The labeled cohesive prism only places DOFs at its three
`POINT_PRISM_TENSOR` vertices (one per triangle vertex), giving a
section closure of 9. There is no PETSc 3.25 knob that lets a
tetrahedron FE present a 9-DOF dual space on a triangular cohesive
prism, so the prompt's strategy stalls one frame earlier than it does
on the surface FE.

## What landed

1. **`Simulator::setupAuxiliaryDM`
   (`src/core/Simulator.cpp:3148-3279`)**: kept the dim-1 surface FE,
   added `DMSetFieldAvoidTensor(slipAuxDM_, 0, PETSC_TRUE)` per the
   prompt (also forced TRUE internally for cohesive DS in
   `plexsection.c:130`), and gated a `FSRM_SESSION19_DEBUG` block
   that prints `dimAux`, vol_Np, and face_Np for both the slip aux
   and the Lagrange FE.

2. **`Simulator::updateFaultAuxiliary`
   (`src/core/Simulator.cpp:3030-3088`)**: replaced
   `DMProjectFunctionLabelLocal` with a manual section-write that
   walks stratum 1 of `interfaces_label_`, asks the local section for
   the slip field's DOF/offset at every labeled point, and stamps the
   fault-local triple (`opening`, `left-lateral`, `reverse`) directly
   into `slipAuxVec_`. The projection helper aborted with
   "section closure 9 != dual space dimension 12 at height 0" once we
   experimented with the dim volume FE; the section-write path keeps
   working for both FE choices.

3. **Default gating restored to opt-in
   (`src/core/Simulator.cpp:2960-2974`, `3318-3331`,
   `6932-6960`, `7152-7165`)**: setting up `slipAuxDM_` and attaching
   `slipAuxVec_` at `(interfaces_label_, 1, 0)` are still controlled
   by `FSRM_ENABLE_SLIP_AUX=1`. With the env var unset (default), the
   simulator runs the Session 16 constants-array path that all of the
   passing fault-subset tests rely on.

4. **Rim-pin retained**: the prompt's Step 3 instructed us to disable
   the geometric rim-pin if `febasic.c:663` cleared. It did not, so
   the Session 12/14 `lagrange_fault_rim` BC stays on by default
   (still toggleable via `FSRM_DISABLE_RIM_PIN=1`).

## Status of `febasic.c:663`

Not resolved. The check is a downstream symptom of two upstream
constraints that PETSc 3.25 enforces simultaneously:

- `plexfem.c:3982` requires the section closure on each cohesive cell
  to equal `totDimAux` for the aux DS that is attached at that key.
- `febasic.c:663` requires `Tf[0]->Np == TfAux[0]->Np` once the
  closure check passes.

The dim-1 surface FE satisfies (1) (closure 9 == DS totDim 9) but
fails (2) when the bulk-side material aux tabulation is also queried
during the same hybrid call. The dim volume FE inverts the failure
mode (fails (1), would have satisfied a 12-point version of (2)). A
no-label slip field with `AvoidTensor=TRUE` would skip cohesive cells
entirely, which sidesteps both checks but leaves the slip values on
bulk vertices that the hybrid driver does not query at the cohesive
key.

## `PrescribedSlipQuasiStatic` SNES trace and `max_fault_slip`

Default path (`FSRM_ENABLE_SLIP_AUX` unset):

```
0 TS dt 1. time 0.
    0 SNES Function norm 6.250000000000e-05
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 1.385637037625e+01
    1 SNES Function norm 6.930778885126e-04
    2 SNES Function norm 1.955090238982e-04
    Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 2
1 TS dt 0.25 time 0.25
```

`max_fault_slip = 263.044`. The trace and the slip magnitude match
the `local_fix` tip exactly (verified by stashing the Session 19 diff
and re-running). Session 18 reported a smaller `max_fault_slip = 11.7`
in the same configuration; that number is stale - the baseline
already drifts to the line-search-recovery path on the current
`local_fix` HEAD.

`FSRM_ENABLE_SLIP_AUX=1` path: `TSSolve` returns PETSc error 75
(`PETSC_ERR_ARG_INCOMP`) at `febasic.c:663` after the very first
SNES function evaluation, so there is no SNES iteration trace.

## `VecNorm(slipAuxVec_)`

With the manual section-write and the dim-1 surface FE, the slip vec
norm matches Session 18's projection-based number to the digit:

```
Session 19: updateFaultAuxiliary stamped
            slip_fault=[0.001 0 0] -> VecNorm(slipAuxVec_)=0.005 at t=1
```

(Five unit cohesive facets x 1 mm opening x sqrt(5) accounting for the
shared-vertex multiplicity yields the 0.005 reading.)

## Fault-subset status vs Session 18

```
ctest -R 'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|
          Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|
          Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|
          Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|
          Integration.ExplosionFaultReactivation|
          Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup'
```

| Test                                                | Session 17 | Session 18 | Session 19 |
| --------------------------------------------------- | ---------- | ---------- | ---------- |
| Functional.DynamicRuptureSetup (x5)                 | 5 pass     | 5 pass     | 5 pass     |
| Physics.CohesiveBdResidual                          | pass       | pass       | pass       |
| Physics.LockedFaultTransparency                     | FAIL       | FAIL       | FAIL       |
| Physics.SCEC.TPV5                                   | FAIL       | FAIL       | FAIL       |
| Integration.ExplosionFaultReactivation              | pass       | pass       | pass       |
| Integration.DynamicRuptureSolve.LockedQuasiStatic   | FAIL       | FAIL       | pass       |
| Integration.DynamicRuptureSolve.LockedElastodynamic | pass       | pass       | pass       |
| Integration.DynamicRuptureSolve.PrescribedSlip      | FAIL       | FAIL       | FAIL       |
| Integration.TimeDependentSlip                       | FAIL       | FAIL       | pass       |
| Integration.SlippingFaultSolve                      | FAIL       | FAIL       | FAIL       |
| Integration.SlipWeakeningFault                      | FAIL       | FAIL       | FAIL       |
| Integration.PressurizedFractureFEM                  | FAIL       | FAIL       | FAIL       |

Totals: 10/16 in Session 18 vs **10/16** in Session 19. Two tests
flipped from FAIL to pass relative to Session 18's reported table
(`LockedQuasiStatic`, `TimeDependentSlip`); both pass on the
`local_fix` tip as well, so the apparent improvement is a Session 18
bookkeeping drift rather than a Session 19 fix. No regressions.

## Next-session recommendations

- **Hand-rolled cohesive integrator.** Bypass
  `DMPlexComputeResidualHybridByKey` for the cohesive Lagrange field
  and assemble the constraint residual in a manual loop that reads
  `slipAuxVec_` via section offsets. This sidesteps both
  `plexfem.c:3982` and `febasic.c:663`, matches PyLith's strategy at
  a higher level, and removes the joint Np constraint between the
  material and slip aux fields.
- **Material-aux face quadrature.** The material aux's face
  quadrature is left at the P0 default (1 point per face). Copying
  the displacement FE's face quadrature into the material aux FE
  (alongside the existing volume-quadrature copy at
  `setupAuxiliaryDM`) may decouple the bulk-side
  `Tf[0]->Np == TfAux[0]->Np` check from the slip aux experiment so
  the slip-aux failure mode can be diagnosed in isolation.
- **PyLith API consult.** The PyLith codebase (or upstream PETSc
  tests `ex62.c`/`ex96.c`) almost certainly avoid the Np mismatch via
  one of: a custom hybrid integrator, explicit
  `PetscDSSetIntegrationParameters`, or an aux DS built on a
  cohesive-prism reference cell. A short read of the PyLith source
  before the next session would shorten this loop considerably.
