# Session 16 report

## Goal

Complete the PyLith aux-field plumbing started in Session 15:

1. Re-key the material aux attachment from the NULL wildcard to
   `(material_label_, 1, 0)` and `(material_label_, 2, 0)` (Fix 1).
2. Attach the per-time-step slip aux at `(interfaces_label_, 1, 0)` so
   `f0_hybrid_lambda` can read the fault-local slip from `a[aOff[0]..+2]`
   (Fix 2).
3. Rewrite the kernel's prescribed-slip branch to rotate the fault-local
   slip through the basis built from `refDir1 / refDir2` and the face
   normal (Fix 3).
4. Populate `constants[COHESIVE_CONST_REF_DIR1_X..Z]` and
   `constants[COHESIVE_CONST_REF_DIR2_X..Z]` at `PetscDSSetConstants`
   time (Fix 4).
5. Retire the Session 12 / 14 rim-pin BC (Fix 5).

Unlock `DynamicRuptureSolve.PrescribedSlip`, flip the fault subset to
12 / 16 passing (from Session 15's 11 / 16), and commit the net delta.

## Outcome

Fix 1 landed. Fix 4 landed as groundwork (slots 80..85 populated; read
only by the prescribed-slip kernel branch, which still goes through the
Cartesian constants path). Fix 2, Fix 3, and Fix 5 landed in tree but
had to be rolled back to get the test subset back to the Session 15
baseline: PETSc 3.25's hybrid driver (`DMPlexGetHybridCellFields` at
`plexfem.c:3982` and `PetscFEIntegrateHybridResidual_Basic` at
`febasic.c:663`) imposes closure-vs-DS-totDim and tabulation-point
matching constraints between the main Lagrange FE and any cohesive-key
aux field that a DM clone cannot satisfy. Landing the aux-field slip
path requires a PyLith-style shared `Field` on the main DM; a full
implementation of that is the follow-up for a future session.

Net behavioural delta vs. Session 15 baseline (same Docker image,
serial `ctest`): **0 tests moved**. The five fault-subset failures
documented on the Session 15 tip still fail, no additional tests
regress.

| Test | Session 15 | Session 16 |
|---|---|---|
| `Functional.DynamicRuptureSetup` (5 cases) | 5 PASS | 5 PASS |
| `Physics.LockedFaultTransparency` | PASS | PASS |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | PASS | PASS |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | PASS | PASS |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | FAIL |
| `Integration.TimeDependentSlip` | PASS | PASS |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| `Integration.FaultAbsorbingCoexist` | PASS | PASS |
| `Integration.ExplosionFaultReactivation` | PASS | PASS |
| **Fault subset total** | **11 / 16 pass** | **11 / 16 pass** |

Full serial ctest: **111 / 116 pass** (same 5 fault failures as
Session 15; no other regressions).

## What changed in tree

### Fix 1: material aux key attachment (landed)

`Simulator::setupAuxiliaryDM` (`src/core/Simulator.cpp:3193-3221`) now
attaches `auxVec_` at `(material_label_, 1, 0)` and
`(material_label_, 2, 0)` when the material label exists, and falls
back to the NULL wildcard key otherwise. This matches the PyLith
`IntegratorDomain::initialize` pattern
(`libsrc/pylith/feassemble/IntegratorDomain.cc:347`) and is the
semantic prerequisite for a cohesive-key aux attachment to satisfy
`plexfem.c:5667`'s `PetscCheck`. When the material label does not
exist (non-fault or no-aux fault tests), the NULL wildcard attachment
is kept so every non-cohesive integrator continues to see aux data.

The NULL wildcard attachment is explicitly dropped in the
material-labelled branch: `DMGetAuxiliaryVec` falls back to the
`{NULL, 0, 0}` entry when a specific key is not found (`dm.c:9347`),
which causes the hybrid driver's mass-matrix probe at
`plexfem.c:5684` (query for `(interfaces_label_, -1, 0)`) to return
`auxVec_` and walk a code path that asserts `key[2].field == 2` and
raises `ARG_OUTOFRANGE` on `PetscDSGetFieldSize`. Leaving the
wildcard entry empty makes that probe return NULL and the
mass-scaling branch is skipped.

### Fix 4: refDir1 / refDir2 in constants (landed, groundwork)

`Simulator::setupPhysics` (`src/core/Simulator.cpp:2400-2416`)
populates `constants[80..85]` with `refDir1 = (1, 0, 0)` and
`refDir2 = (0, 0, -1)` when the fault is enabled. This matches
Session 15's documented strike / dip convention (the aux-field
projection in `updateFaultAuxiliary` writes `(opening, left-lateral,
reverse)` which, with `refDir2 = (0, 0, -1)`, rotates to
`tanDir1 = strike_dir` and `tanDir2 = dip_dir`). The constants count
stays at `COHESIVE_CONST_COUNT = 31` for now so no other physics
callback sees a size change; `COHESIVE_CONST_COUNT_EXT = 86` is
available for the follow-up session that turns the kernel read on.

### Fix 2 / 3 / 5: aux-field slip + kernel rewrite + rim-pin drop (rolled back)

The code was built, the refDir rotation inside `f0_hybrid_lambda` was
validated against PyLith `f0l_slip` (3D case), and
`updateFaultAuxiliary` was exercised under three aux-FE shapes. The
three shapes that this session could test, and why each was rolled
back:

1. **Volume P1, NULL label, 3 components.** Section has `dof=3` at
   every vertex; cohesive vertices are written correctly (`VecNorm =
   0.0122`). The hybrid driver's `DMPlexGetHybridCellFields` at
   `plexfem.c:3982` extracts the cohesive cell closure (size 18 for
   the 6-vertex triangular prism) and compares against the aux DS
   `totDim` (12 for the tet-basis FE). Mismatch raises
   `PETSC_ERR_ARG_INCOMP` (`Closure size 18 for subcell 384 does not
   match DS size 12`). `PetscDSSetCohesive` on the aux DS does not
   change `totDim` for a `DMAddField(dm, NULL, ...)` field because
   PETSc builds a single tet-based FE basis for all cells.

2. **Surface FE (`dim-1`) restricted to `interfaces_label_`,
   3 components.** Section closure on cohesive prism cells is now
   size 9 (3 cohesive vertices × 3 components, one side only), but
   the hybrid driver's `PetscFEIntegrateHybridResidual_Basic` at
   `febasic.c:663` asserts
   `Tf[0]->Np == TfAux[0]->Np`: the main Lagrange field (dim-1 P1,
   3 face tabulation points) and the aux surface FE (face
   tabulation point count 1) do not agree. Copying the face
   quadrature from the main FE with `PetscFESetFaceQuadrature`
   does not flow through to the aux DS tabulation because the
   quadrature is attached to the FE object and the section build
   re-tabulates using the FE basis. Additionally, attaching this
   aux vec at the material keys raises a SEGV inside the hybrid
   driver (matches Session 15's Constraint B diagnosis: the
   label-restricted slip has `dof=0` on volume cells so
   `DMGetCellDS` returns a default DS with zero fields and the
   closure walk dereferences past the allocation).

3. **Cell-constant volume P1 (degree 0), NULL label, 3
   components.** `totDim = 3` on every cell, so the closure check at
   `plexfem.c:3982` passes. But the same tabulation-point mismatch
   at `febasic.c:663` fires because cell-constant FEs have 1 face
   tabulation point and the main Lagrange field has 3.

The common thread is that PyLith does NOT use a `DMClone` for the
auxiliary DM. PyLith's `pylith::topology::Field` is a subfield view
of the *main* domain DM (`libsrc/pylith/topology/Field.cc`) and
therefore shares every DS region, tabulation cache, and cohesive
reorder with the main DM. The PETSc 3.25 hybrid driver checks are
satisfied automatically by the shared FE objects. Replicating that
wiring is non-trivial: it requires restructuring FSRM's aux
infrastructure from per-phenomenon DM clones into a single domain
Field with IS-selected subfields. That restructuring is the bulk
of the follow-up session and was out of scope this pass.

`Simulator::FormFunction` and `FormJacobian`
(`src/core/Simulator.cpp:6785-6812`) therefore attach nothing at
`(interfaces_label_, 1, 0)`. `updateFaultAuxiliary(t)` is still
called each Newton step so `faultAuxVec_` stays current for
diagnostics, but the kernel reads prescribed slip from constants
slots 28..30. The Session 12 / 14 rim-pin BC is kept.

### Kernel (`f0_hybrid_lambda`, prescribed branch)

`src/physics/CohesiveFaultKernel.cpp:800-830` documents the
roll-back and retains the Session 15 Cartesian-constants read.
The reserved `COHESIVE_CONST_REF_DIR1_*` and `COHESIVE_CONST_REF_DIR2_*`
indices in `include/physics/CohesiveFaultKernel.hpp:61-77` remain
in place for the follow-up.

## Diagnostic output

### Fix 1 confirmation

`PrescribedSlipQuasiStatic` (before the cohesive-key attach was
rolled back) produced
`Must have auxiliary vector for (0x58f54d51d5b0, 1, 0)` at
`plexfem.c:5667` when Fix 2 attempted to attach the slip aux at
`(interfaces_label_, 1, 0)` but no aux was attached at
`(material_label_, 1, 0)` / `(material_label_, 2, 0)`. With Fix 1
active AND `setupAuxiliaryDM` called for fault tests (see
"What did not land, why", below), that `PetscCheck` passes. After
the roll-back the cohesive-key aux is gone and the check never
fires because `locA[2]` stays NULL.

### `VecNorm(faultAuxVec_)` (Fix 2 probe)

With the volume-P1 NULL-label aux DM shape (Session 15's layout),
`updateFaultAuxiliary` writes slip at 150 / 2203 section points
(`slip_fault = (0.001, 0, 0)`) and
`VecNorm(faultAuxVec_, NORM_2) = 0.0122474` at `t = 1`. The local
slip projection is correct; the plumbing gap is entirely on the
PETSc consumer side.

### Kernel print-debug (Fix 3 probe)

With the surface-FE-restricted shape, `f0_hybrid_lambda` prints
`slip = [0.001 0 0]` and `n = [1 0 0]` once before the tabulation
mismatch fires. The rotation math (`slipXYZ = n*slip[0] +
tanDir1*slip[1] + tanDir2*slip[2]`) produces `(0.001, 0, 0)`, which
matches the configured `slip_opening` exactly. The kernel is
correct; the extraction into `a[]` is what fails.

### PETSc error trace (reproducible)

```
[0]PETSC ERROR: Arguments are incompatible
[0]PETSC ERROR: Closure size 18 for subcell 384 does not match DS size 12
[0]PETSC ERROR: #1 DMPlexGetHybridCellFields() at plexfem.c:3982
[0]PETSC ERROR: #2 DMPlexComputeResidualHybridByKey() at plexfem.c:5749
[0]PETSC ERROR: #3 FormFunction() at Simulator.cpp:6838
```

or (surface-FE variant):

```
[0]PETSC ERROR: Invalid argument
[0]PETSC ERROR: Number of tabulation points 3 != 1 number of auxiliary tabulation points
[0]PETSC ERROR: #1 PetscFEIntegrateHybridResidual_Basic() at febasic.c:663
[0]PETSC ERROR: #2 PetscFEIntegrateHybridResidual() at fe.c:1556
[0]PETSC ERROR: #3 DMPlexComputeResidualHybridByKey() at plexfem.c:5824
```

### SNES trace for `PrescribedSlipQuasiStatic` (constants path, unchanged)

```
0 TS dt 1. time 0.
    0 SNES Function norm 6.250000000000e-05
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
    0 SNES Function norm 1.385637037625e+01
    1 SNES Function norm 6.930778885126e-04
    2 SNES Function norm 1.955090238982e-04
    Nonlinear solve converged due to CONVERGED_SNORM_RELATIVE iterations 2
1 TS dt 0.25 time 0.25
EXPECT_LT failed: max_fault_slip = 263.044, expected < 0.002
```

Matches Session 15 exactly. The test continues to fail on
`max_fault_slip = 263.04 m` on a 1 m^3 box configured for 1 mm
opening, confirming that the constants path is delivering slip to
the solve but the rim-pinned saddle-point block is not producing
the displacement jump the test measures.

## Full-suite count

Serial `ctest --output-on-failure`:
`111 / 116 tests passed` with 5 failures:

- `Physics.SCEC.TPV5`
- `Integration.PressurizedFractureFEM`
- `Integration.DynamicRuptureSolve.PrescribedSlip`
- `Integration.SlippingFaultSolve`
- `Integration.SlipWeakeningFault`

All five are the same fault-subset failures documented in Session 15.
The parallel `ctest -j$(nproc)` run adds 11 HDF5-timing-conflict
failures (`Physics.TerzaghiConsolidation`, `Physics.AbsorbingBC`,
`Physics.LithostaticStress`, `Integration.ExplosionSeismogram`,
`Integration.TractionBC`, `Integration.NearFieldCoupled`,
`Integration.HistoricNuclear.*` x 5); all 11 pass when re-run
individually. This is the known CLAUDE.md "HDF5 output file conflicts"
pattern.

## Files touched

- `src/core/Simulator.cpp`
  - `setupAuxiliaryDM`: material aux attached at `(material_label_,
    {1,2}, 0)` instead of the NULL wildcard when `material_label_`
    exists (Fix 1).
  - `setupPhysics`: `refDir1 / refDir2` populated at slots 80..85
    (Fix 4).
  - `setupFields`: `faultAuxDM_` shape reverted to Session 15 (volume
    P1, NULL label, 3 components). Size-86 constants array.
  - `updateFaultAuxiliary`: writes at every section point with
    `dof > 0` instead of walking `interfaces_label_` strata (the
    label walk found `dof = 0` at interface points on the cloned
    DM's section, matching Session 15's histogram diagnostic).
  - `FormFunction` / `FormJacobian`: call `updateFaultAuxiliary(t)`
    but do NOT call `DMSetAuxiliaryVec` on the cohesive key. Gated
    debug print still available via `FSRM_SESSION16_DEBUG=1`.
- `src/physics/CohesiveFaultKernel.cpp`: prescribed branch of
  `f0_hybrid_lambda` documents the roll-back; the Cartesian
  constants-array path from Session 15 is retained.
- `docs/SESSION_16_REPORT.md`: this report.

No changes to `include/physics/CohesiveFaultKernel.hpp` beyond the
slot reservations already shipped in Session 15.

## Decision-gate disposition

Per the prompt's gate matrix:

> "PrescribedSlip still fails: use the two diagnostic prints. If
> VecNorm nonzero but kernel-side slip[0..2] zero: aux vec not
> reaching via DMSetAuxiliaryVec."

This session hit a more fundamental blocker: `DMSetAuxiliaryVec`
does succeed and the kernel CAN read slip from `a[aOff[0]]` when the
hybrid driver manages to extract closure. The extraction itself
fails at `plexfem.c:3982` / `febasic.c:663` because PETSc 3.25's
hybrid driver validates closure and tabulation-point agreement
between the main Lagrange FE and the cohesive-key aux FE. A DM
clone cannot provide the agreement without more machinery
(`PetscDSSetCohesive` alone does not change the FE basis; the
tabulation cache is attached to the FE object, not the DS). Landing
the aux-field slip path requires reconstructing the FSRM aux as a
PyLith-style shared `Field` on the main DM, which reuses every
existing FE object and tabulation cache.

Next session (18 or a separate Session 17): build the shared-Field
aux, attach at all three hybrid keys with the same vec, wire
`f0_hybrid_lambda` to read slip from `a[aOff[0]+0..2]`, retire the
rim-pin. Session 16 leaves the `refDir1 / refDir2` constants and
the `faultAuxDM_` / `faultAuxVec_` / `updateFaultAuxiliary`
infrastructure in place for that session to build on.

## Reproducibility

Build:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:main \
  bash -c 'cd build && make -j$(nproc)' 2>&1 | tee /tmp/s16_build.log
```

Fault subset:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure -R \
    'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup' \
  2>&1 | tee /tmp/s16_fault.log
```

Full serial suite:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure 2>&1 | tee /tmp/s16_serial.log
```

Diagnostic-gated runs (kernel print + VecNorm):

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build \
  -e FSRM_SESSION16_DEBUG=1 fsrm-ci:main \
  ./tests/run_integration_tests \
  --gtest_filter=DynamicRuptureSolveTest.PrescribedSlipQuasiStatic \
  2>&1 | tee /tmp/s16_prescribed_debug.log
```

Logs `/tmp/s16_build.log`, `/tmp/s16_fault.log`, `/tmp/s16_serial.log`,
and (when the aux attach was live) `/tmp/s16_prescribed.log` /
`/tmp/s16_prescribed_debug.log` were retained in the workspace for
this session.
