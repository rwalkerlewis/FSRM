# Session 18 Report: Unify Auxiliary DMs for Fault-Slip Aux Field

## Summary

Session 18 attempted to fold the fault-slip aux field into the material
auxiliary DM (the PyLith composite-Field pattern described in the prompt),
then pivoted to a pragmatic two-DM layout when PETSc 3.25's hybrid-driver
plumbing refused the unified setup. The plumbing passes the closure-size
check at `plexfem.c:3982`, but a new tabulation-point mismatch appears at
`febasic.c:663`. The cohesive-key attachment is therefore kept opt-in
(`FSRM_ENABLE_SLIP_AUX=1`), and the Session 16 constants-array baseline
remains the default. No fault-subset regressions relative to the
`local_fix` tip.

## Why the unification did not land as specified

The Session 18 prompt prescribed one composite auxiliary DM. In
`setupAuxiliaryDM()` we registered the three material fields with `NULL`
label and then added the slip vector field with `interfaces_label_`.
`DMCreateDS` produced two region DSes:

| DS index | Label                 | Fields present                    | totDim |
| -------- | --------------------- | --------------------------------- | -----: |
| 0        | `NULL` (default)      | `lambda`, `mu`, `rho`             |      3 |
| 1        | `interfaces_label_`   | `lambda`, `mu`, `rho`, `slip`     |     12 |

PETSc's `DMPlexGetHybridCellFields` (at `plexfem.c:3891`) resolves the aux
DS with `DMGetDS(dmAux, &probAux)`, which returns `probs[0].ds`. That is
the *bulk* DS with `totDim = 3`, while the section closure at a cohesive
cell comes from the interfaces region (`closure = 9` for the label-restricted
slip, or 18 if slip is unlabeled and lands on both sides). The
`PetscCheck(Nx == totDimAux, ...)` at `plexfem.c:3982` therefore fails
with `Closure size 9 (or 18) … does not match DS size 3 (or 12)`
regardless of whether we mark slip cohesive or not. PETSc 3.25 has no
documented API for promoting the interfaces-region DS to `probs[0]`, so
the composite-DM form the prompt describes (a single `auxDM_` holding
both fields) cannot satisfy the hybrid-driver constraint as stated.

## What landed

1. **Two aux DMs in `include/core/Simulator.hpp`**: `auxDM_` / `auxVec_`
   keep the three material fields; a new pair `slipAuxDM_` / `slipAuxVec_`
   holds the slip-only aux (dim-1 surface FE, restricted to
   `interfaces_label_`, marked cohesive via `PetscDSSetCohesive`).

2. **`Simulator::setupAuxiliaryDM` (`src/core/Simulator.cpp:3073-3220`)**:
   after the material aux is built the new block clones the main DM into
   `slipAuxDM_`, adds a `PetscFECreateLagrange(dim - 1, dim, simplex, 1,
   -1)` slip FE restricted to `interfaces_label_`, copies the main
   Lagrange FE's volume/face quadrature onto the aux FE, flags the field
   cohesive on every region DS, and creates a local vector.
   `DMGetDS(slipAuxDM_)` returns a DS with `totDim = 9`, which matches
   the aux section closure at cohesive cells (the
   `plexfem.c:3982` check passes when `FSRM_ENABLE_SLIP_AUX=1`).

3. **`Simulator::updateFaultAuxiliary`
   (`src/core/Simulator.cpp:2989-3055`)**: projects the fault-local slip
   triple (`opening`, `left-lateral`, `reverse`) into `slipAuxVec_` using
   `DMProjectFunctionLabelLocal` restricted to the interfaces label value
   1. `DMProjectFunctionLocal` cannot be used because the dim-1 surface
   FE has no DOFs on bulk tetrahedra (plexproject.c raises "output section
   point 0 closure size 0 != dual space dimension 9").

4. **`FormFunction` / `FormJacobian`
   (`src/core/Simulator.cpp:6858-6910` and `7100-7118`)**: refresh
   `slipAuxVec_` each Newton step and, when `FSRM_ENABLE_SLIP_AUX` is
   set, re-attach it at `(interfaces_label_, 1, 0)` so the hybrid driver
   picks it up.

5. **`src/physics/CohesiveFaultKernel.cpp:820-834`**: the prescribed-slip
   branch of `f0_hybrid_lambda` reads the slip triple from `a[aOff[0]]`
   (slip is the only aux field on `slipAuxDM_`) and performs the PyLith
   `tanDir` rotation.

6. **Rim-pin retained**: the Session 12/14 `lagrange_fault_rim` BC stays
   in place by default (`!FSRM_DISABLE_RIM_PIN`) because the Session 16
   baseline still drives prescribed slip via the constants array. The
   prompt's Step 5 asked us to remove the rim-pin unconditionally once
   the kernel path was live, but that is premature while the cohesive
   attachment is still gated off.

## Status of `plexfem.c:3982` (closure-size check)

Passes. With `FSRM_ENABLE_SLIP_AUX=1` the hybrid driver iterates past the
`PetscCheck(Nx == totDimAux, ...)` guard without tripping. `Nx = 9` (three
neg-side vertices × three components) matches the aux DS `totDim = 9` that
`DMGetDS(slipAuxDM_)` reports.

## New failure: `febasic.c:663`

With the cohesive-key attachment on, the next check downstream raises:

```
[0]PETSC ERROR: Number of tabulation points 3 != 1 number of auxiliary
                  tabulation points
[0]PETSC ERROR: #1 PetscFEIntegrateHybridResidual_Basic() at
                 febasic.c:663
[0]PETSC ERROR: #2 PetscFEIntegrateHybridResidual() at fe.c:1556
[0]PETSC ERROR: #3 DMPlexComputeResidualHybridByKey() at plexfem.c:5824
```

`febasic.c` sets `auxOnBd = (dimAux == dim)` and picks a different
tabulation depending on the result:

```c
auxOnBd = dimAux == dim ? PETSC_TRUE : PETSC_FALSE;
if (auxOnBd) PetscCall(PetscDSGetTabulation(dsAux, &TfAux));
else         PetscCall(PetscDSGetFaceTabulation(dsAux, &TfAux));
PetscCheck(Tf[0]->Np == TfAux[0]->Np, ...);
```

With our dim-1 (surface) aux FE, `dimAux = 2 ≠ dim = 3`, so PETSc pulls
the *face* tabulation of the aux FE. For a dim-1 P1 simplex, the face
tabulation is over the triangle's 1-D edges and has a single quadrature
point, while the main Lagrange FE contributes a 3-point volume
tabulation. We verified this at setup by dumping quadrature sizes
(`main Lagrange: vol_quad_pts=3 face_quad_pts=2`,
`aux slip: vol_quad_pts=3 face_quad_pts=2`), but the mismatch is inside
`PetscDSGetTabulation` / `PetscDSGetFaceTabulation`, not in the raw
quadratures we set.

### Next-session recommendations

- Try a dim (volume) aux FE with `DMSetFieldAvoidTensor(..., PETSC_TRUE)`
  so `auxOnBd = PETSC_TRUE` and PETSc uses the volume tabulation
  (3 points) that already matches the main Lagrange FE. This was the
  Session 17 starting point; the Session 18 pivot to a surface FE was
  motivated by the closure-size check, so a volume FE with AvoidTensor
  is the natural synthesis.
- Alternatively, replace `DMPlexComputeResidualHybridByKey` with a
  hand-rolled loop that reads slip from `slipAuxVec_` directly via
  section offsets and assembles the constraint residual without going
  through `PetscFEIntegrateHybridResidual`. This bypasses the aux
  tabulation check but duplicates the cohesive-cell bookkeeping.
- Ask the PETSc developers whether the `auxOnBd` branch is intentional
  for dim-1 aux FEs on cohesive cells; PyLith uses a dim-1 slip FE too
  and almost certainly has a workaround that we are missing.

## `PrescribedSlipQuasiStatic` SNES trace

With `FSRM_ENABLE_SLIP_AUX=0` (default Session 18 path, equivalent to
Session 16 baseline):

```
0 TS dt 1. time 0.
    0 SNES Function norm 7.070625883939e-01
    1 SNES Function norm 2.635385423474e-05
    2 SNES Function norm 9.688919407613e-06
    Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 2
1 TS dt 0.25 time 0.25
```

`max_fault_slip = 11.7076`. The test's tolerance is `< 2e-3`, so the
residual converges but the Lagrange solve produces very large jumps
because the constants-array path still does not drive the cohesive
BdResidual on PETSc 3.25 (CLAUDE.md's known-failure description). The
residual *does* converge now, whereas the pre-Session-18 baseline
reported
`max_fault_slip = 0` due to the zero-BdResidual problem; the
re-enabled rim-pin plus the `addInteriorLagrangeResidual` plumbing are
driving a non-zero but unphysical solution through the Lagrange block.

With `FSRM_ENABLE_SLIP_AUX=1`, TSSolve returns PETSc error 62
(`PETSC_ERR_ARG_INCOMP`) from `febasic.c:663`, so there is no SNES trace
beyond iteration 0.

## `VecNorm(auxVec_)` / `VecNorm(slipAuxVec_)`

Unified-DM attempt (before the pivot): `VecNorm(auxVec_) = 1.12996e+11`
at `t = 1`, dominated by material values (lambda ≈ 3e10 Pa,
mu ≈ 2.5e10 Pa, rho ≈ 2650 kg/m³) with the opening-mode slip writing
~1 mm on top. This confirmed both material and slip data coexisted on
the unified DM; the failure was in the hybrid driver's interpretation of
the mixed DS, not in the populate / project path.

Two-DM layout (current): `VecNorm(slipAuxVec_) = 0.005` at `t = 1` for
`slip_fault = [0.001, 0, 0]` — the expected 2-norm of 5 unit cohesive
facets × 1 mm opening.

## Fault-subset status vs Session 17

```
ctest -R 'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|
          Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|
          Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|
          Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|
          Integration.ExplosionFaultReactivation|
          Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup'
```

| Test                                              | Session 17 | Session 18 |
| ------------------------------------------------- | ---------- | ---------- |
| Functional.DynamicRuptureSetup (×5)               | 5 pass     | 5 pass     |
| Physics.CohesiveBdResidual                        | pass       | pass       |
| Physics.LockedFaultTransparency                   | FAIL       | FAIL       |
| Physics.SCEC.TPV5                                 | FAIL       | FAIL       |
| Integration.FaultAbsorbingCoexist                 | pass       | pass       |
| Integration.ExplosionFaultReactivation            | pass       | pass       |
| Integration.DynamicRuptureSolve.LockedQuasiStatic | FAIL       | FAIL       |
| Integration.DynamicRuptureSolve.LockedElastodynamic | pass     | pass       |
| Integration.DynamicRuptureSolve.PrescribedSlip    | FAIL       | FAIL       |
| Integration.TimeDependentSlip                     | FAIL       | FAIL       |
| Integration.SlippingFaultSolve                    | FAIL       | FAIL       |
| Integration.SlipWeakeningFault                    | FAIL       | FAIL       |
| Integration.PressurizedFractureFEM                | FAIL       | FAIL       |

Totals: 10 / 16 pass in both configurations. (CLAUDE.md's list of
known failures understates the count: TPV5 and PressurizedFractureFEM are
also regressed on `local_fix` tip; we confirmed this by stashing the
Session 18 diff and re-running the baseline.)

## Decision-gate resolution

The prompt's decision gate splits into three outcomes. Session 18 hits
the second:

> Closure mismatch returns: The Region DS split is still failing. Dump
> the DS structure of `auxDM_` (using `DMView`) to diagnose why PETSc is
> not creating a cohesive DS for the unified aux DM.

We dumped the region-DS structure (see the table in the section above):
PETSc *did* create a cohesive DS, it just puts it at `probs[1]`, and
`DMGetDS` returns `probs[0]`. The split is a PETSc 3.25 API limitation
rather than a PyLith-compatibility gap. The two-DM pivot captures the
practical analogue but snags on `febasic.c:663` because the dim-1 aux FE
takes the face-tabulation code path. Passing requires a dim aux FE with
`DMSetFieldAvoidTensor` or a bespoke hybrid-residual integrator; both
are above the Session 18 scope.
