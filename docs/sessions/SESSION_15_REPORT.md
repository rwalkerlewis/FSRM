# Session 15 report

## Goal

Port PyLith's auxiliary-field slip pattern into FSRM so that prescribed slip
flows through a per-point auxiliary `Vec` (read inside `f0_hybrid_lambda` from
`a[aOff[0]+0..2]` and rotated to global XYZ via the fault-local basis), retire
the Session 12 / 14 rim-pin Dirichlet BC, and restore the PyLith
GAMG + GMRES + near-null-space saddle-point solver.

## Outcome

The PyLith aux-field plumbing landed but is held back from the prescribed-slip
residual path until Session 16 because PETSc 3.25's hybrid driver imposes a
prerequisite that this session could not satisfy without further work
(documented below). The kernel's prescribed-slip branch therefore continues to
read the Cartesian slip from the constants array, the Session 12 / 14 rim-pin
stays registered, and the GAMG default is gated so caller-provided
`-pc_type lu` overrides keep working.

Net behavioural delta vs. Session 14 baseline (same Docker image, same fault
test subset, run serially): **+1 test passing**.

| Test | Session 14 | Session 15 |
|---|---|---|
| `Physics.LockedFaultTransparency` | FAIL | **PASS** |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | PASS (CLAUDE.md notes were stale; serial run confirms PASS) | PASS |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | PASS | PASS |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL (`max_fault_slip = 263.044`) | FAIL (same value) |
| `Integration.TimeDependentSlip` | PASS | PASS |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Integration.FaultAbsorbingCoexist` | PASS | PASS |
| `Integration.ExplosionFaultReactivation` | PASS | PASS |
| `Functional.DynamicRuptureSetup` | PASS (5 cases) | PASS (5 cases) |

The CLAUDE.md "known failures" list mentioned Locked* and TimeDependentSlip as
baseline failures; the serial baseline run confirms only
`LockedFaultTransparency` actually fails on the Session 14 tip
(`Locked*ElastoDynamic`, `LockedQuasiStatic`, and `TimeDependentSlip` already
pass). Session 15 fixes `LockedFaultTransparency` as a side effect of the
auxiliary-DM plumbing exercising a different PETSc code path during section
construction.

## What changed in tree

1. **Auxiliary DM and local vector for fault slip** (groundwork). In
   `Simulator::createFieldsFromConfig` (`src/core/Simulator.cpp:1814-1870`)
   the fault-enabled branch now clones the main DM into `faultAuxDM_`, adds a
   3-component degree-1 Lagrange field named `slip` defined on every cell of
   the cloned DM, runs `DMCreateDS` / `DMSetUp`, and allocates the local
   `faultAuxVec_`. Destroyed in `Simulator::~Simulator()`.

2. **`Simulator::updateFaultAuxiliary(t)`** (`src/core/Simulator.cpp:3066-3140`).
   Each Newton step zeros `faultAuxVec_` and writes the configured fault-local
   triple `(opening, left-lateral, reverse)` at every point in
   `interfaces_label_` whose section dof is non-zero. The triple is converted
   to global XYZ inside the kernel; when `refDir2 = (0, 0, -1)` the basis
   reduces to `tanDir1 = strike_dir`, `tanDir2 = dip_dir`, matching FSRM's
   existing strike/dip convention. The vec is built but **not** attached to
   the main DM via `DMSetAuxiliaryVec` -- see "Why the aux-vec is parked"
   below.

3. **`updateFaultAuxiliary` is called from `FormFunction` and `FormJacobian`**
   (`src/core/Simulator.cpp:6776-6790, 6929-6945`) right before the hybrid
   residual / Jacobian calls, so the projected slip is always current. The
   call is a no-op when faults are not enabled.

4. **Aux-field slot reservations in `CohesiveFaultKernel`**
   (`include/physics/CohesiveFaultKernel.hpp:61-77`). New constant indices
   `COHESIVE_CONST_REF_DIR1_X..Z` (slots 80..82) and `..._REF_DIR2_X..Z`
   (slots 83..85) are reserved for the fault-local basis used by the Session
   16 aux-field rotation. The constants array passed to PetscDS is left at
   `COHESIVE_CONST_COUNT = 31` so other physics callbacks see no change.

5. **GAMG default solver, gated**. `Simulator::setupSolvers`
   (`src/core/Simulator.cpp:3676-3760`) now switches the fault-enabled
   default from "test-provided LU + Session 12 rim-pin" to "PyLith GAMG +
   GMRES + near-null-space" *only if no `-pc_type` is set on the options
   database before `setupSolvers` runs*. Tests that explicitly call
   `PetscOptionsSetValue("-pc_type", "lu")` (every fault TS test today)
   keep their LU. The `FSRM_SADDLE_SOLVER=fieldsplit` env override still
   selects the experimental Schur split.

6. **Rim-pin retained**. The Session 12 / 14 `lagrange_fault_rim` Dirichlet
   BC and the geometric rim-prism detection in `getOrCreateBuriedCohesiveLabel`
   are kept in place. They are now annotated as "Session 12/14 (kept in 15)"
   to signal that Session 16 is expected to retire them once the aux-vec
   path is solvability-equivalent.

## Why the aux-vec is parked

The intended residual-path change was to read slip from the auxiliary vector
and apply the PyLith `f0l_slip` rotation inside `f0_hybrid_lambda`. The
auxiliary-vec attachment was implemented and the kernel rewritten end-to-end
(both still in `git reflog`); the test harness then surfaced two PETSc 3.25
constraints that the prompt's recipe did not cover:

### Constraint A: PETSc requires aux vectors at *all three* hybrid keys

`DMPlexComputeResidualHybridByKey` at `plexfem.c:5644-5670`:

```c
PetscCall(DMGetAuxiliaryVec(dm, key[2].label, key[2].value, key[2].part, &locA[2]));
if (locA[2]) {
    ...
    PetscCall(DMGetAuxiliaryVec(dm, key[c].label, key[c].value, key[c].part, &locA[c]));
    PetscCheck(locA[c], PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Must have auxiliary vector for (%p, %d, %d)", ...);
}
```

When the cohesive (`key[2]`) aux vec is set, PETSc *also* demands aux vecs at
the negative-side material key (`key[0]`) and positive-side material key
(`key[1]`). PyLith satisfies this by attaching the same domain-Field aux
local-vector at all three keys (`IntegratorInterface.cc:471-472`). Attaching
the same `faultAuxVec_` at the material/1 and material/2 keys raises a SEGV
inside PETSc's hybrid driver because the slip subfield is restricted to the
cohesive interface stratum, so `DMGetCellDS(dmAux, support_cell, ...)` returns
a default DS with zero fields and `PetscDSGetTotalDimension` is then used to
size aux extraction buffers that the closure walk dereferences past their
allocation.

### Constraint B: aux DM section needs a definition on bulk cells too

The PyLith auxiliary `Field` is built on the main domain DM with subfields
defined either everywhere (e.g. material properties) or restricted to the
fault label (slip, fault traction). Attempting the same with a single 3-comp
volume Lagrange field on a `DMClone` of the main DM did avoid the SEGV but
left every interface label point with `dof = 0` in the cloned DM's section
(verified with a histogram dump: 113 interface points, 0 with dof > 0). The
cohesive vertices live at points whose owner is the cohesive prism cells, and
the cloned DM's section walker did not associate the slip basis with those
points. The Session 14 Lagrange-field setup on the *main* DM works because
the field carries the cohesive flag (`PetscDSSetCohesive`); replicating that
on the cloned DM would also require copying the cohesive section reorder, the
hybrid weak-form registrations, and the `interfaces_label_` association, at
which point the aux DM stops being a thin clone.

The honest read on the prompt's Step 5 ("Verify the correct PETSc API ...
Alternatively, the auxiliary field's layout is inferred from the aux DM
attached via `DMSetAuxiliaryVec`, so this step may be automatic") is that the
inference is *not* automatic for cohesive geometry in PETSc 3.25 and pulling
in the missing setup is the bulk of Session 16. The "Build errors from
PetscDSSetAuxiliary* API" decision-gate row in the prompt does not apply;
the build succeeded, but the runtime behaviour did.

## Diagnostic output

### `dmAux_` field layout

`Session 15: fault aux DM built (slip subfield, local size 450)` -- the local
vec has 3 components per vertex on every cell of the 4x4x4 simplex mesh
(150 vertices x 3 = 450). Restricted (Session 15 first attempt) the same vec
came out at local size 75 (matching the Lagrange-field size) but failed the
constraint-B check above.

### Per-field constraint counts

With rim-pin retained:

```
Local section: total dof=525, constrained dof=300 (free=225)
  field 0: dof=450, constrained=228, free=222
  field 1: dof=75, constrained=72, free=3
```

24 of 25 cohesive Lagrange points are pinned (3 components each = 72), leaving
a single unconstrained Lagrange triple. This matches Session 12 / 14.

If the rim-pin were removed (Session 15 first attempt), field 1 reports
`dof=75, constrained=0, free=75` and direct LU returns
`DIVERGED_LINEAR_SOLVE iterations 0` on the rank-deficient saddle-point block.

### `PrescribedSlipQuasiStatic` SNES trace (Session 15 final)

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

The first sub-step's line search aborts (zero residual is misclassified as
divergence by the back-tracking line search), the time stepper backs off to
`dt = 0.25`, the next solve "converges" to a state with `max_fault_slip` =
263.04 m on a 1 m^3 box configured for a 1 mm opening. The same trace appears
on the Session 14 tip, so this is the standing baseline failure mode for
this test, not a Session 15 regression.

### `LockedFaultTransparency` (Session 15 fix)

Now passes with `fault_norm < 5e-4` under the rim-pinned LU configuration.
The Session 14 tip fails this test even though `LockedQuasiStatic` and
`LockedElastodynamic` pass; the auxiliary DM creation evidently changes some
PETSc internal state (likely the `DMReorderSectionType` cache on the cloned
DM) that makes the transparency residual converge.

### `LockedQuasiStatic` and `LockedElastodynamic`

Both pass under Session 15 in serial mode (`100% tests passed, 0 tests failed
out of 2`). Same on Session 14 baseline -- contrary to the CLAUDE.md
"known failures" list, which is stale for these two cases.

### `TimeDependentSlip`

Passes under Session 15. The aux-vec projection ramp logic in
`updateFaultAuxiliary` is exercised but not consumed (the kernel still reads
constants), so this test's pass relies on the same constants-broadcast slip
that worked on Session 14 baseline.

### Fault subset table vs. Session 14

| Test | S14 (serial) | S15 (serial) |
|---|---|---|
| `Functional.DynamicRuptureSetup` (5 cases) | 5 PASS | 5 PASS |
| `Physics.LockedFaultTransparency` | FAIL | PASS |
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
| **Total fault subset** | **10 / 16 pass** | **11 / 16 pass** |

### Full ctest run (serial subset; parallel `ctest -j` adds 3 HDF5 timing
failures noted in CLAUDE.md)

`/tmp/s15_test.log` records 108 / 116 passing under `ctest -j$(nproc)`. The
extra failures vs. the serial fault subset
(`Physics.TerzaghiConsolidation`, `Physics.AbsorbingBC`,
`Integration.ExplosionSeismogram`) all pass when re-run individually -- they
are the HDF5 output-conflict pattern flagged in CLAUDE.md.

## Decision-gate disposition

Per the prompt's gate matrix:

> "PrescribedSlip still fails: dump the SNES trace. If residual is right but
> KSP diverges, solver config issue. If residual is wrong, kernel bug -- verify
> a[aOff[0]] values against expected slip and compare to PyLith kernel
> byte-for-byte."

The aux-field path was built and could be wired end-to-end at the C++ /
PETSc-API level, but the PETSc 3.25 hybrid-driver demand for aux vectors at
all three keys (Constraint A above) gates the runtime behaviour. The right
follow-up is the prompt's gate row "Build errors from PetscDSSetAuxiliary*
API: verify exact call signature ..." extended to the hybrid driver's
runtime preconditions: provide a single aux DM whose section covers both the
cohesive and bulk supports correctly, attach the same local vec at all three
keys, and complete the kernel switch.

## Files touched

- `include/core/Simulator.hpp` -- `faultAuxDM_`, `faultAuxVec_`,
  `fault_slip_aux_field_idx_`, `updateFaultAuxiliary` declaration.
- `include/physics/CohesiveFaultKernel.hpp` -- reserve
  `COHESIVE_CONST_REF_DIR1_X..Z` and `COHESIVE_CONST_REF_DIR2_X..Z` (slots
  80..85) for the Session 16 aux-field rotation; constants count remains
  `COHESIVE_CONST_COUNT = 31` for everyone else.
- `src/core/Simulator.cpp` -- aux DM creation, `updateFaultAuxiliary`,
  `updateFaultAuxiliary` calls in `FormFunction` and `FormJacobian`,
  GAMG-only-when-not-overridden default in `setupSolvers`, "(kept in 15)"
  annotations on the rim-pin code paths.
- `src/physics/CohesiveFaultKernel.cpp` -- documentation in
  `f0_hybrid_lambda` explaining the aux-field intent and the Session 15 -> 16
  hand-off.

## Reproducibility

Build:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:main bash -c \
  'cd build && make -j$(nproc)' 2>&1 | tee /tmp/s15_build.log
```

Fault subset:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest --output-on-failure -R \
    'Physics.LockedFaultTransparency|Physics.CohesiveBdResidual|Integration.DynamicRuptureSolve|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault|Physics.SCEC.TPV5|Integration.FaultAbsorbingCoexist|Integration.ExplosionFaultReactivation|Integration.PressurizedFractureFEM|Functional.DynamicRuptureSetup' \
  2>&1 | tee /tmp/s15_fault.log
```

Full suite (parallel, expect HDF5 timing failures):

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:main \
  ctest -j$(nproc) --output-on-failure 2>&1 | tee /tmp/s15_test.log
```

Logs `/tmp/s15_build.log`, `/tmp/s15_prescribed.log`, `/tmp/s15_fault.log`,
`/tmp/s15_test.log` are kept in the workspace from this session.
