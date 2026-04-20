# Session 17 report

## Goal

Session 16 rolled back the aux-field prescribed-slip plumbing because every
aux-FE shape it tried tripped a PETSc 3.25 hybrid-driver check. The Session 17
prompt diagnosed the root cause as mismatched quadrature order between the
main Lagrange FE and the aux slip FE, and asked this session to replicate
PyLith's `FieldOps::createFE` verbatim so the tabulation-point count agrees
inside `PetscFEIntegrateHybridResidual_Basic`.

Specifically, land:

1. PyLith-style FE construction for the aux slip, with explicit
   `PetscSpace` + `PetscDualSpace` + `PetscDTCreateDefaultQuadrature` calls
   that mirror `libsrc/pylith/topology/FieldOps.cc:39-125`, using the same
   `basisOrder` and `quadOrder` as the main Lagrange FE.
2. `updateFaultAuxiliary(t)` rewritten to drive `DMProjectFunctionLocal` on
   the new aux DM.
3. Cohesive-key aux attachment via `DMSetAuxiliaryVec(dm, interfaces_label_,
   1, 0, faultAuxVec_)` in `FormFunction` / `FormJacobian`.
4. `f0_hybrid_lambda` prescribed branch rewritten to read the fault-local
   slip triple from `a[aOff[0]..+2]` and rotate into Cartesian via
   `refDir1 / refDir2` (slots 80..85, already populated in Session 16).
5. Rim-pin BC retired (no longer needed once the aux path drives the
   displacement jump at every Lagrange DOF).

Unlock `DynamicRuptureSolve.PrescribedSlip`, flip the fault subset to
12 / 16 passing, and commit the net delta.

## Outcome

All five code items from the prompt landed in tree. The PyLith FE
construction is exactly the recipe from `FieldOps.cc` and produces an aux
FE whose main quadrature (tet order-1, 4 points) and face quadrature
(triangle order-1, 3 points) are both reachable via the PyLith code path.
With the aux vec attached at `(interfaces_label_, 1, 0)`, however, PETSc
3.25's hybrid driver still rejects the plumbing: the closure-vs-DS check
at `DMPlexGetHybridCellFields` (`plexfem.c:3982`) trips
`Closure size 18 for subcell 384 does not match DS size 12` on the
`DMClone`-based aux DM. Removing the `DMClone` (i.e. reusing the main DM's
field structure the PyLith `Field` way) is a larger restructure than this
session could complete.

To preserve the Session 16 test baseline, the new aux-field path is gated
behind `FSRM_ENABLE_SLIP_AUX=1`. With the default (env-var off), the aux
DM is built exactly as in Session 16 (volume P1 via `PetscFECreateLagrange`
+ `DMAddField(NULL)`), the cohesive-key aux attachment is skipped,
`f0_hybrid_lambda` falls back to the constants-array prescribed slip, and
the rim-pin BC stays active. All five Session 16 fault failures remain;
no tests regress.

Fault-subset delta vs. Session 16 baseline (same Docker image, serial
`ctest`): **0 tests moved**.

| Test | Session 16 | Session 17 |
|---|---|---|
| `Functional.DynamicRuptureSetup` (5 cases) | 5 PASS | 5 PASS |
| `Physics.LockedFaultTransparency` | PASS | PASS |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | PASS | PASS |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | PASS | PASS |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | FAIL |
| `Integration.TimeDependentSlip` | PASS | PASS |
| `Integration.ExplosionFaultReactivation` | PASS | PASS |
| `Integration.FaultAbsorbingCoexist` | PASS | PASS |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| **Fault subset total** | **11 / 16 pass** | **11 / 16 pass** |

Full serial `ctest --output-on-failure`: **111 / 116 pass** with the same
five fault failures as Session 16. No other regressions.

## Diagnostic output

### Aux slip FE / main Lagrange FE parameters (with `FSRM_ENABLE_SLIP_AUX=1`)

PyLith-style FE construction with `main_basis_order = 1`, `main_quad_order
= 1`:

```
Session 17: fault aux DM built (path=pylith, local vec size=450)
Session 17: aux slip FE face quadrature nqp=3
Session 17: aux slip FE basisDeg=1 dualOrder=1 nqp=4
Session 17: main Lagrange FE basisDeg=1 dualOrder=1 nqp=3
```

Both FEs have `basisDeg=1`, `dualOrder=1`. The aux FE main quadrature has
4 points (tet volume, Stroud-conical order 1), its face quadrature has
3 points (triangle order 1). The main Lagrange FE (which is a surface FE
in 3D) has 3 main quadrature points (triangle order 1). The aux's face
quadrature matches the main Lagrange's main quadrature — which is what
the prompt predicted would satisfy
`PetscFEIntegrateHybridResidual_Basic`'s tabulation-point count check at
`febasic.c:660`.

### `VecNorm(faultAuxVec_)` before / after `updateFaultAuxiliary`

```
Session 17: updateFaultAuxiliary projected slip_fault=[0.001 0. 0.] ->
            VecNorm(faultAuxVec_)=0.0122474 at t=1.
Session 17: VecNorm(faultAuxVec_) = 0.0122474 at t=1.
```

`DMProjectFunctionLocal` writes the constant fault-local triple
`(0.001, 0, 0)` at every section point on the cloned DM. The resulting
VecNorm of 0.01225 matches sqrt(150 * 0.001^2) which is consistent with
150 dof-bearing section points each holding `slip_fault[0] = 0.001`.

### Residual-call failure with aux attached

The kernel-side debug print
(`Session 17: kernel sees slip=[...] n=[...] refDir2=[...]`) does NOT
fire because the hybrid driver rejects the closure before reaching
the kernel:

```
[0]PETSC ERROR: Arguments are incompatible
[0]PETSC ERROR: Closure size 18 for subcell 384 does not match DS size 12
[0]PETSC ERROR: #1 DMPlexGetHybridCellFields() at plexfem.c:3982
[0]PETSC ERROR: #2 DMPlexComputeResidualHybridByKey() at plexfem.c:5749
[0]PETSC ERROR: #3 FormFunction() at Simulator.cpp:6919
Session 17: DMPlexComputeResidualHybridByKey returned 75 (Arguments are incompatible)
```

This is the same `plexfem.c:3982` failure Session 16 saw for Shape 1 (volume
P1, NULL label). The Session 17 fix targeted the `febasic.c:660` tabulation
check (the one Session 16 reported at `febasic.c:663`), which IS now
satisfied by matching `basisOrder`/`quadOrder`. But it uncovers a prior
check: `DMPlexGetHybridCellFields`'s closure-vs-DS-totDim comparison still
fails on the `DMClone`-based aux DM because the aux DS built by
`DMCreateDS(faultAuxDM_)` has `totDim = 12` (tet P1 x Nc=3) while the
section closure on a cohesive prism cell is 18 (6 vertices x 3). PyLith
avoids this by reusing the main DM's fields through `pylith::topology::Field`
rather than cloning the DM; that restructure is out of scope for Session 17.

### SNES trace and `max_fault_slip` for `PrescribedSlipQuasiStatic`

With the default env (aux path off), the solve follows the Session 16
constants-array path:

```
0 TS dt 1. time 0.
    0 SNES Function norm 6.250000000000e-05
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    0 SNES Function norm 1.385637037625e+01
    1 SNES Function norm 6.930778885126e-04
    2 SNES Function norm 1.955090238982e-04
    Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 2
1 TS dt 0.25 time 0.25
EXPECT_LT failed: max_fault_slip = 263.044, expected < 0.002
```

Same failure mode as Session 16.

With `FSRM_ENABLE_SLIP_AUX=1`, the SNES never iterates because the first
`FormFunction` call returns `PETSC_ERR_ARG_INCOMP` from the hybrid residual
call.

## Decision-gate disposition

Per the prompt's gate matrix:

> Closure-size mismatch fires: `DMSetFieldAvoidTensor` did not propagate.
> Inspect whether section on cohesive cells has zero or nonzero slip DOFs;
> if nonzero, AvoidTensor isn't working in this code path.

Confirmed. With `DMSetField(faultAuxDM_, 0, NULL, fe_slip)` +
`DMSetFieldAvoidTensor(faultAuxDM_, 0, PETSC_TRUE)` on the cloned DM, the
section still produces closure size 18 on cohesive prism cells
(`dof = 3` per vertex, 6 vertices per prism) rather than the size-9 closure
(`dof = 3` per unique vertex, 3 vertices per side) AvoidTensor would give
for a `pylith::Field`-style construction on the main DM. Walking the aux
DM's region DSes confirmed there is only one DS, marked cohesive or not,
and `PetscDSSetCohesive` applied to it does not change `totDim`.

A second diagnostic path tried a surface FE restricted to
`interfaces_label_` (Session 16 Shape 2 equivalent) with the PyLith
`main_quad_order = 1` recipe. That shape produced
`totDim = 9` on the labeled DS (closure check would pass) and the aux's
main quadrature was 3 points matching the main Lagrange's 3. But
`PetscFEIntegrateHybridResidual_Basic` still reported aux tabulation =
1 point, i.e. the hybrid driver queries the FACE tabulation of the aux FE
in this configuration, which at `quadOrder=1` is an edge 1-point rule.
Setting `main_quad_order = 2` (per the gate matrix's common-fix
suggestion) moved the counts but did not reconcile them. The surface-FE
Shape 2 also triggers `DMProjectFunctionLocal` to assert
`closure size 0 != dual space dimension 9` on non-cohesive cells (the
aux has no dof there), so it additionally requires
`DMProjectFunctionLabelLocal`. Both of these are documented in-code as
reference material for the follow-up.

## What changed in tree

### `src/core/Simulator.cpp`

- `Simulator::createFieldsFromConfig` (`src/core/Simulator.cpp:1814 ff.`):
  the fault aux DM now clones the main DM and builds the slip FE via
  either Session 16's `PetscFECreateLagrange` path (default) or a
  PyLith-style `PetscSpace` + `PetscDualSpace` +
  `PetscDTCreateDefaultQuadrature` path (`FSRM_ENABLE_SLIP_AUX=1`). In the
  PyLith path, `DMSetField(faultAuxDM_, 0, NULL, fe_slip)` + `DMSetField
  AvoidTensor(faultAuxDM_, 0, PETSC_TRUE)` mirror
  `libsrc/pylith/topology/Field.cc:556-562` (`isFaultOnly = false`). In the
  default path the FE is added via `DMAddField(faultAuxDM_, nullptr, ...)`
  as Session 16.
- `Simulator::updateFaultAuxiliary` (`src/core/Simulator.cpp:3124 ff.`):
  rewritten to use `DMProjectFunctionLocal` with a static
  `prescribed_slip_constant_fn` callback that returns the fault-local slip
  triple `(opening, left-lateral, reverse) = (fault_slip_opening_,
  fault_slip_strike_, fault_slip_dip_)`. Gated diagnostic print behind
  `FSRM_SESSION17_DEBUG=1`.
- `Simulator::FormFunction` / `FormJacobian`
  (`src/core/Simulator.cpp:6919`, `7060`): when
  `FSRM_ENABLE_SLIP_AUX=1` is set, call
  `DMSetAuxiliaryVec(dm, interfaces_label_, 1, 0, faultAuxVec_)` before the
  hybrid residual / Jacobian call. Wrapped with a `PetscTraceBackErrorHandler`
  push/pop under `FSRM_SESSION17_DEBUG=1` so the hybrid driver's failure
  trace is visible under the test harness's own return-error handler.
- `Simulator::setupFields` `setupAuxiliaryDM` call
  (`src/core/Simulator.cpp:3017`): when `FSRM_ENABLE_SLIP_AUX=1` AND faults
  are enabled, force `setupAuxiliaryDM()` even if `use_aux_callbacks` is
  false. This satisfies `plexfem.c:5667`'s PetscCheck that demands a
  material aux at `(material_label_, 1, 0)` / `(material_label_, 2, 0)`
  whenever any cohesive-key aux is attached.
- `Simulator::setupBoundaryConditions` (`src/core/Simulator.cpp:8140`):
  the Session 12 / 14 `lagrange_fault_rim` BC is gated. It stays on by
  default so Session 16's locked-test convergence is preserved. Setting
  `FSRM_DISABLE_RIM_PIN=1` opts out (only meaningful when
  `FSRM_ENABLE_SLIP_AUX=1` is also set and the hybrid-driver issue has
  been resolved).
- `Simulator::run` (`src/core/Simulator.cpp:6775`): diagnostic
  `PetscTraceBackErrorHandler` push/pop around `TSSolve` under
  `FSRM_SESSION17_DEBUG=1` to surface the full stack trace that the test
  harness's `PetscReturnErrorHandler` otherwise hides.

### `src/physics/CohesiveFaultKernel.cpp`

- `f0_hybrid_lambda` (`src/physics/CohesiveFaultKernel.cpp:800-880`):
  prescribed-slip branch rewritten to read the fault-local triple from
  `a[aOff[0]..+2]` and rotate into Cartesian using the PyLith tangent
  basis
  `tanDir1 = normalize(refDir2 x n), tanDir2 = normalize(n x tanDir1)`.
  When `NfAux = 0` (aux not attached, default gating), falls back to the
  constants-array Cartesian slip from Session 15 / 16.
- Gated kernel-side print behind `FSRM_SESSION17_DEBUG=1`
  (`kernel sees slip=[...]`).

### No header changes

`include/physics/CohesiveFaultKernel.hpp` is unchanged; `COHESIVE_CONST_REF_DIR1_*`
and `COHESIVE_CONST_REF_DIR2_*` slots were already reserved in Session 15
and populated in Session 16.

## What did not land, why

- The unconditional `DMSetAuxiliaryVec(dm, interfaces_label_, 1, 0,
  faultAuxVec_)` that the Session 17 prompt specified.
- The unconditional rim-pin retirement.
- The unconditional `DMSetFieldAvoidTensor` on the default aux DM path.

All three are in-tree but gated behind `FSRM_ENABLE_SLIP_AUX=1` because
the cohesive-key attach still trips a PETSc 3.25 hybrid-driver check on
the cloned aux DM. The matching-quadOrder prescription from
`FieldOps::createFE` does pass `febasic.c:660` (tabulation-point count
check) with `basisOrder = 1, quadOrder = 1`, but a prior closure-size
check at `DMPlexGetHybridCellFields` (`plexfem.c:3982`) fires first and
reports `Closure size 18 != DS size 12`. That check compares the aux DM's
section closure on a cohesive prism (6 vertices x 3 components = 18) to
the aux DS's `totDim` (4 vertices x 3 = 12 for a tet P1 DS). Making the
aux DS cohesive (`PetscDSSetCohesive`) was attempted and did not change
`totDim` because the cloned DM exposes only one DS, not a separate
cohesive-region DS the way the main DM does.

The blocker is the same "PyLith uses a shared `Field` on the main DM, not
a `DMClone`" one that Session 16 documented. Landing the aux-field slip
path needs FSRM's aux infrastructure restructured so the slip field is a
subfield of the main DM's solution (or a parallel `pylith::Field`-style
view into it), not a separate cloned DM. That restructure is the
follow-up.

## Files touched

- `src/core/Simulator.cpp` (aux DM setup gating, `updateFaultAuxiliary`
  rewrite, `FormFunction` / `FormJacobian` cohesive-key attach gating,
  `setupAuxiliaryDM` force-call gating, rim-pin gating, diagnostic hooks).
- `src/physics/CohesiveFaultKernel.cpp` (`f0_hybrid_lambda` prescribed
  branch rewritten to read from `a[]` with constants-array fallback; gated
  kernel print).
- `docs/SESSION_17_REPORT.md` (this report).

## Reproducibility

Build:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:main \
  bash -c 'cd build && make -j$(nproc)' 2>&1 | tee /tmp/s17_build.log
```

Default-path fault subset (matches Session 16 baseline, 11 / 16 pass):

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure -R \
    'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup' \
  2>&1 | tee /tmp/s17_fault.log
```

Default-path full serial suite (111 / 116 pass):

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure 2>&1 | tee /tmp/s17_serial.log
```

Diagnostic run exercising the gated aux-field path
(`FSRM_ENABLE_SLIP_AUX=1 FSRM_SESSION17_DEBUG=1`):

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build \
  -e FSRM_ENABLE_SLIP_AUX=1 -e FSRM_SESSION17_DEBUG=1 fsrm-ci:main \
  ./tests/run_integration_tests \
  --gtest_filter=DynamicRuptureSolveTest.PrescribedSlipQuasiStatic \
  2>&1 | tee /tmp/s17_slipaux.log
```

The last command will fail with `PETSC_ERR_ARG_INCOMP` (75) and exercise
the full aux-field plumbing: PyLith FE construction,
`DMProjectFunctionLocal` writing fault-local slip into `faultAuxVec_`
(`VecNorm = 0.0122474`), `DMSetAuxiliaryVec` on the cohesive key, and the
hybrid driver's subsequent closure-size rejection at `plexfem.c:3982`.

Logs `/tmp/s17_build.log`, `/tmp/s17_fault3.log`, `/tmp/s17_serial.log`,
and `/tmp/s17_slipaux_run.log` were retained in the host workspace for
this session.
