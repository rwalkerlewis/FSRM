# Session 4 Report

Branch: `local_fix`.
Scope: replace the cohesive residual / Jacobian path with the PETSc hybrid
driver (`DMPlexComputeResidualHybridByKey` / `DMPlexComputeJacobianHybridByKey`)
instead of `DMPlexTSComputeIFunctionFEM` plus `PetscDSSetBdResidual` on
cohesive cells.

**Decision gate outcome: STOP AND REPORT — new regressions introduced.**

---

## 1. Material label

`Simulator::createMaterialLabel` (`src/core/Simulator.cpp` line ~3578) walks
every cohesive prism cell, takes `cone[0]` / `cone[1]` to reach the
negative-side and positive-side faces, finds the non-cohesive support cell
of each face, and writes label value 1 (negative side) / 2 (positive side)
onto that regular cell.

Diagnostic output (identical across every fault test run this session):

```
'material' label created: 32 cohesive cells -> 32 neg (value=1), 32 pos (value=2)
```

Every test produced non-zero counts on both sides. The fault is not on the
domain boundary and orientation logic is healthy.

---

## 2. Build outcome

Clean. `fsrm` executable and all test binaries link. Full log captured to
`/tmp/s4_build.log` (no errors, only the pre-existing `-Wunused-*` warnings
also present on `main`).

PETSc API signatures verified against the installed headers at
`/opt/petsc-main/include/` (PETSc `v3.25.0-117-g267b8824abd`, the version
referenced by the Session 4 ground-truth block):

- `PETSC_EXTERN PetscErrorCode DMPlexComputeResidualHybridByKey(DM, PetscFormKey[], IS, PetscReal, Vec, Vec, PetscReal, Vec, void *);`
- `PETSC_EXTERN PetscErrorCode DMPlexComputeJacobianHybridByKey(DM, PetscFormKey[], IS, PetscReal, PetscReal, Vec, Vec, Mat, Mat, void *);`
- `PETSC_EXTERN PetscErrorCode PetscWeakFormSetIndexBdResidual(PetscWeakForm, DMLabel, PetscInt, PetscInt, PetscInt, PetscInt, void (*)(...), PetscInt, void (*)(...));`
- `PETSC_EXTERN PetscErrorCode PetscWeakFormSetIndexBdJacobian(PetscWeakForm, DMLabel, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, ..., PetscInt, ..., PetscInt, ..., PetscInt, ...);`
- `PETSC_EXTERN PetscErrorCode DMPlexGetSimplexOrBoxCells(DM, PetscInt, PetscInt *, PetscInt *);`

No PETSc API mismatch during compile.

---

## 3. Fault test status (before / after)

Baseline values come from the pre-Session-4 tip of `local_fix` as documented
in `CLAUDE.md`. "After" values come from the Session 4 run captured to
`/tmp/s4_fault.log` and the full-suite run `/tmp/s4_test.log`.

| Test                                                         | Baseline | Session 4  |
|--------------------------------------------------------------|----------|------------|
| Physics.CohesiveBdResidual                                   | PASS     | FAIL (new) |
| Physics.LockedFaultTransparency                              | FAIL     | FAIL       |
| Physics.SCEC.TPV5                                            | PASS     | FAIL (new) |
| Integration.DynamicRuptureSolve.LockedQuasiStatic            | FAIL     | FAIL       |
| Integration.DynamicRuptureSolve.LockedElastodynamic          | PASS     | FAIL (new) |
| Integration.DynamicRuptureSolve.PrescribedSlip               | FAIL     | FAIL       |
| Integration.TimeDependentSlip                                | FAIL     | FAIL       |
| Integration.SlippingFaultSolve                               | FAIL     | FAIL       |
| Integration.SlipWeakeningFault                               | FAIL     | FAIL       |

Additional new regressions observed in the full-suite run that the eight-test
filter does not cover:

- Functional.DynamicRuptureSetup.LockedFault
- Functional.DynamicRuptureSetup.FaultGeometryFromConfig
- Functional.DynamicRuptureSetup.SlippingFault
- Integration.PressurizedFractureFEM
- Integration.ExplosionFaultReactivation

Full-suite summary (`ctest -j`): **19 failed of 116**. Of those 19, 6 match
the `CLAUDE.md` baseline fault failures, 5 are known parallel-HDF5 artifacts
(`Physics.AbsorbingBC`, `Integration.ExplosionSeismogram`,
`Integration.HistoricNuclear.{Gasbuggy1967, Gnome1961, NtsPahuteMesa}` — the
`CLAUDE.md` rule "some tests may fail when run in parallel (`ctest -j`) due
to HDF5 output file conflicts"), and 8 are **new fault-related regressions
introduced by the Session 4 change**.

PrescribedSlip `max_fault_slip`: not measurable. The test aborts during
`DMPlexInsertBoundaryValues` before `TSSolve` runs, and the written
`output/solution.h5` only contains step 0 (initial conditions).
LockedQuasiStatic `sol_norm`: zero (same reason — the test exits before the
first Newton iteration).

---

## 4. Root cause of the new regressions

Progression during the session:

1. **First run (hybrid driver only, no cohesive-field marking):**
   every fault test SEGVs with `Caught signal number 11 SEGV`. The fault
   stratum trace always completes (the material-label diagnostic always
   prints), so the crash is inside
   `DMPlexComputeResidualHybridByKey` / `...JacobianHybridByKey` themselves.

2. **Second run (after adding `PetscDSSetCohesive(prob, lagrange_field, PETSC_TRUE)` in `setupFields()`):**
   the SEGV is replaced with a clean PETSc error:

   ```
   [0]PETSC ERROR: Nonconforming object sizes
   [0]PETSC ERROR: The output section point (7) closure size 24 != dual space dimension 36 at height 0 in [0, 0]
   [0]PETSC ERROR: #1 DMProjectLocal_Generic_Plex() at plexproject.c:973
   [0]PETSC ERROR: #2 DMProjectFunctionLabelLocal_Plex() at plexproject.c:1123
   [0]PETSC ERROR: #4 DMPlexInsertBoundaryValuesEssential() at plexfem.c:926
   [0]PETSC ERROR: #7 setInitialConditions() at Simulator.cpp:5260
   ```

   The failure is now in the *projection* of Dirichlet BC values, not in the
   hybrid residual integration. Marking the Lagrange field cohesive switches
   the DS from volume-FE offset logic (`PetscDSGetFieldOffset`) to cohesive
   offset logic (`PetscDSGetFieldOffsetCohesive`), which changes the expected
   closure layout at every tet cell from 24 (4 tet vertices × (3 disp + 3
   Lagrange)) to 36 (dim-1 cohesive FE dual-space dimension).

The mismatch is structural:

- FSRM creates the Lagrange field as a **3-D volume** FE:
  `PetscFECreateLagrange(comm, 3, 3, isSimplex, 1, -1, &fe_lagrange)` and adds
  it with a `nullptr` label, so DOFs sit on **every** cell.
- PETSc's hybrid cohesive machinery (see `src/dm/impls/plex/tests/ex5.c`
  lines 766-795, and the `PetscDSGetCohesive` gate inside
  `DMPlexGetHybridCellFields` at `plexfem.c:3864`) expects the cohesive
  field to be a **(dim-1)-dimensional surface FE** added with a label that
  restricts it to cohesive cells, e.g.
  `PetscFECreateDefault(comm, dim-1, dim, …, "faulttraction_", …)` +
  `DMAddField(dm, fault_label, fe)`.

Until the Lagrange field is recreated with that shape, the hybrid driver
cannot extract DOFs into the expected `u[]` layout. The registered
`f0_hybrid_lambda` / `g0_hybrid_lambda_*` kernels never get invoked with
valid data — control aborts upstream.

The Session 4 exemption granted Rule 6 only for `registerWithDS`. Creating
the Lagrange FE with a different (dim, label) pair is outside that
exemption, so the work stops here per the decision gate.

---

## 5. Code landed on `local_fix`

Uncommitted changes on the working tree (to be committed as session 4):

- `include/core/Simulator.hpp`:
  - New members `DMLabel material_label_`, `PetscInt displacement_field_idx_`,
    `PetscInt lagrange_field_idx_`.
  - Declaration of `PetscErrorCode createMaterialLabel()`.
- `src/core/Simulator.cpp`:
  - `setupFields()` caches `displacement_field_idx_`, `lagrange_field_idx_`.
  - After `DMCreateDS`: `PetscDSSetCohesive(ds, lagrange_field_idx_, PETSC_TRUE)`.
  - `setupFaultNetwork()` calls `createMaterialLabel()`.
  - `createMaterialLabel()` populates the `"material"` DMLabel (1 = neg,
    2 = pos) by walking `cone(cohesive_cell)[0..1]` and their supports.
  - `setupPhysics()`: the non-`use_region_ds_split` cohesive block no
    longer calls `PetscDSSetBdResidual` on the displacement / Lagrange
    fields for non-hydrofrac fault modes. Instead it calls
    `cohesive_kernel_->registerHybridWithDS(...)`. The old BdResidual
    calls are retained as block comments per task spec.
  - The `use_region_ds_split` cohesive block comments out its cohesive
    BdResidual calls (absorbing BdResidual path preserved).
  - `FormFunction` runs `DMPlexComputeResidualHybridByKey` on the cohesive
    cell IS (built from `DMPlexGetSimplexOrBoxCells` vs.
    `DMPlexGetHeightStratum`) after the existing
    `DMPlexTSComputeIFunctionFEM` call and before the manual
    injection/explosion/traction/fault-pressure passes.
  - `FormJacobian` runs `DMPlexComputeJacobianHybridByKey` in the same
    place relative to `DMPlexTSComputeIJacobianFEM`.
- `include/physics/CohesiveFaultKernel.hpp`:
  - New `registerHybridWithDS(prob, wf, material_label, cohesive_label, disp_field, lagrange_field)`.
  - New static callbacks: `f0_hybrid_u_neg`, `f0_hybrid_u_pos`,
    `f0_hybrid_lambda`, `g0_hybrid_u_lambda_neg`, `g0_hybrid_u_lambda_pos`,
    `g0_hybrid_lambda_u`, `g0_hybrid_lambda_lambda`.
- `src/physics/CohesiveFaultKernel.cpp`:
  - `registerWithDS` body stubbed out (fprintf-warns and returns zero).
  - `registerHybridWithDS` registers the three weak forms against
    `(material, 1/2, disp_field)` and `(cohesive_label, 1, lagrange_field)`
    via `PetscWeakFormSetIndexBdResidual` / `...BdJacobian`.
  - Hybrid callbacks written per the Session 4 layout:
    `u[c]` = negative-side displacement, `u[Nc+c]` = positive-side
    displacement, `u[uOff[1]+c]` = Lagrange multiplier;
    `Nc = uOff[2] - uOff[1]`.
  - Slipping-mode Coulomb-friction semi-smooth Newton residual / Jacobian
    ported from the pre-Session-4 `f0_lagrange_constraint` and
    `g0_lagrange_displacement` / `g0_lagrange_lagrange`, adjusted for the
    hybrid `slip = u[Nc+d] - u[d]` layout and with the disp-neg / disp-pos
    block split in the `g0_hybrid_lambda_u` output.

What was **not** changed, per the Session 4 scope rules:

- `FaultMeshManager::splitMeshAlongFault` — untouched.
- Lagrange FE spatial dimension, components, degree, and
  `DMAddField` label — untouched (this is the blocking factor
  documented in Section 4 above).
- `addInteriorLagrangeResidual`, `addCohesivePenaltyToJacobian`, the
  `f0_weak_lagrange` / `g0_weak_lagrange` volume regularization — all
  left in place per task instruction E ("Leave in place for now. Do NOT
  remove in this session even though they will become redundant").

---

## 6. Decision-gate verdict

Applying the gate block-by-block:

- **"PrescribedSlip max_fault_slip > 5e-4 AND CohesiveBdResidual finite AND
  LockedFaultTransparency < 5e-4 AND no regression of currently-passing
  tests"** — **NOT MET**. Eight currently-passing tests regress (see
  Section 3).
- **"PrescribedSlip still 0 or near-zero: report the SNES log in full"** —
  this is the closest matching branch. `TSSolve` never executes; the
  process terminates in `DMPlexInsertBoundaryValues` during
  `setInitialConditions`. There is no SNES log to capture. The abort
  reason captured above.
- **"Build fails with PETSc API mismatch"** — N/A; build is clean.
- **"Any other regression: stop and report"** — **MATCHES**.

Outcome: **stop**. Do not declare architectural fix done. Do not clean up
per task instruction E. Leave `CLAUDE.md` / `README.md` unchanged per
Rule 14.

---

## 7. Recommended Session 5 entry point

Session 5 should request an explicit extension of the Rule 6 / scope
exemption to cover the Lagrange field FE construction. With that in hand,
the minimal change to unlock the hybrid driver is:

1. In `createFieldsFromConfig()`, replace

   ```cpp
   PetscFECreateLagrange(comm, 3, 3, isSimplex, 1, -1, &fe_lagrange);
   DMAddField(dm, nullptr, (PetscObject)fe_lagrange);
   ```

   with a dim-1 surface FE attached to the fault label, mirroring PETSc
   `tests/ex5.c:777-778`:

   ```cpp
   PetscInt dim;
   DMGetDimension(dm, &dim);
   PetscFECreateDefault(PETSC_COMM_SELF, dim - 1, dim, isSimplex,
                        "lagrange_", PETSC_DETERMINE, &fe_lagrange);
   DMAddField(dm, interfaces_label_, (PetscObject)fe_lagrange);
   ```

2. Remove the volume regularization (`f0_weak_lagrange` / `g0_weak_lagrange`
   `PetscDSSetResidual` / `PetscDSSetJacobian` calls in `setupPhysics`).
3. Remove `addInteriorLagrangeResidual(locF)` and the interior-DOF loop in
   `addCohesivePenaltyToJacobian` — with the Lagrange field restricted to
   cohesive geometry there are no interior Lagrange DOFs to zero out.
4. Verify `PetscDSSetCohesive(ds, lagrange_field_idx_, PETSC_TRUE)` (already
   added this session) still matches the new FE layout.
5. Re-run the eight fault tests; expect `DMPlexComputeResidualHybridByKey`
   to actually integrate the registered hybrid weak forms.

Steps 1 and 2 may also require revisiting the `use_region_ds_split` path
for fault + absorbing coexistence, which still carries legacy commented-out
BdResidual calls.

---

## 8. Log artifacts

- `/tmp/s4_build.log` — clean build.
- `/tmp/s4_fault.log` — initial fault-test run (SEGV in hybrid driver
  before `PetscDSSetCohesive` fix).
- `/tmp/s4_fault2.log` — retest after `PetscDSSetCohesive` fix,
  showing the "closure size 24 != dual space dimension 36" error.
- `/tmp/s4_test.log` — full `ctest -j` suite run after the fix
  (19 failed of 116).
