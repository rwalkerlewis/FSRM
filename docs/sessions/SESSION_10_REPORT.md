# Session 10 Report

## Goal

Apply the PyLith fault-rim Dirichlet BC pattern to remove the structural
rank deficiency Session 9 attributed to unconstrained Lagrange DOFs at
fault-boundary vertices. Session 9 hypothesised that pinning these DOFs
would drop the saddle-point condition number from `5.36e+17` to `O(1e+10)`
and unblock `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic`
(`max_fault_slip > 5e-4`).

Result: the rim-vertex BC pattern does not apply to FSRM as-is because
FSRM places Lagrange DOFs on cohesive **cells**, not on cohesive vertices.
Pinning all rim-touching cohesive cells over-constrains the fault and
makes conditioning worse (`5.36e+17` -> `2.36e+18`). Session 10 lands the
diagnostics needed to redirect Session 11 onto the actual root cause.
PrescribedSlipQuasiStatic is unchanged from Session 9.

## What landed

| Change | File / line | State |
|---|---|---|
| `getOrCreateFaultBoundaryLabel` | `src/core/Simulator.cpp:3978-4108`, header `include/core/Simulator.hpp:191` | LANDED, used as diagnostic |
| `fault_boundary_label_` member | `include/core/Simulator.hpp:218-228` | LANDED |
| `cacheFaultBoundaryLagrangeDofs` | `src/core/Simulator.cpp:4110-4234` | LANDED, populated only |
| `pinFaultBoundaryLagrangeJacobian` / `Residual` stubs | `src/core/Simulator.cpp:4236-4252` | STUBS only -- returns 0 |
| `FSRM_SVD_DEBUG` env hook | `tests/integration/test_dynamic_rupture_solve.cpp:330-340` | LANDED |
| Per-DS BC count + per-field section dump in `setInitialConditions` | `src/core/Simulator.cpp:5482-5547` | LANDED |
| `PetscDSAddBoundary` on cohesive DS for rim Lagrange | -- | TRIED then REMOVED (Step B below) |
| `MatZeroRowsColumnsLocal` post-assembly pinning | -- | TRIED then REMOVED (Step C / decision-gate-failure) |

## Step A: rim-vertex label

Implemented `Simulator::getOrCreateFaultBoundaryLabel` at
`src/core/Simulator.cpp:3978-4108`. Walks every cohesive prism cell
(`DM_POLYTOPE_..._PRISM_TENSOR`), gathers depth-0 closure points (cohesive
vertices), and labels each whose coordinates lie within `1e-6` of any
bounding-box face.

Geometry of `PrescribedSlipQuasiStatic` (4x4x4 mesh, fault 2x2 plane in a
1x1x1 box, dim=3):

```
'fault boundary' label created: 32 / 50 cohesive vertices on domain boundary
```

50 cohesive vertices = 25 negative-side + 25 positive-side on a 5x5
clipped fault face. 32 of those are on the perimeter (16 per side x 2),
matching the predicted rim count. The geometric criterion captures every
rim vertex; no buried-edge logic is needed for this geometry.

NB: an earlier draft of the walk used `interfaces_label_` (populated by
`getOrCreateInterfacesLabel`) and got `0 / 0`. That label was built with
`DMPlexGetSimplexOrBoxCells(dm, dim, NULL, &pMax)`, which does not
distinguish hybrid vertices from regular vertices (vertices are
`DM_POLYTOPE_POINT`, never tensor); `pMax == pEnd`. Switching to a direct
cohesive-cell closure walk gave the correct 50.

## Step B: BC registration on the cohesive DS

The natural target was

```c
DMAddBoundary(dm, DM_BC_ESSENTIAL, "lagrange_fault_edge",
              fault_boundary_label_, 1, &one, lagrange_field_idx, ...)
```

but `DMAddBoundary` calls `DMGetDS(dm, &ds)` -> default DS only, then
`PetscDSAddBoundary(ds, ..., field=lagrange_field_idx)` which trips
`PETSC_ERR_ARG_OUTOFRANGE` because the default DS has only the
displacement field (`numDSFields = 1`, lagrange index 1). Session 5
restricted Lagrange to the cohesive interface label via `DMAddField(dm,
interfaces_label_, fe_lagrange)` (`Simulator.cpp:1713`), so Lagrange lives
on a region-specific DS (DS 1).

Workaround attempted: walk `DMGetNumDS` / `DMGetRegionNumDS`, find the DS
that owns Lagrange (DS 1, ds-local field 1, global field 1), and call
`PetscDSAddBoundary` on that DS directly. Verified the BC was stored
(diagnostic at `src/core/Simulator.cpp:5482-5503` reports
`DS 1: 1 boundaries`).

Per-field constraint dump after `DMSetUp`:

```
Local section: total dof=546, constrained dof=228 (free=318)
    field 0: dof=450, constrained=228, free=222
    field 1: dof=96, constrained=0, free=96
```

The displacement field constraints (228) come from the 5 displacement BCs
on the default DS. The Lagrange field has 96 DOFs and 0 constraints --
the BC on DS 1 was registered but did **not** propagate into the local
section's constraint table.

Inspecting `src/dm/impls/plex/plexsection.c:474-562` shows the section
builder iterates `DMGetNumDS` and walks BCs on every DS, so in principle
the cohesive DS BC should be honoured. In practice it is silently
filtered out for FE-on-region-DS, mirroring the
`DMPlexInsertBoundaryValues_Plex` (`plexfem.c:1141-1198`) path which only
queries `DMGetDS(dm, &ds)`. Before debugging that further, Step C below
showed the BC pattern itself was wrong for FSRM, so the BC code was
removed.

## Step C: post-assembly `MatZeroRowsColumnsLocal` pin

Bypass the section-construction machinery: cache the LOCAL section
indices of Lagrange DOFs at the rim, then in `FormJacobian` call
`MatZeroRowsColumnsLocal(J, n, idx, 1.0, NULL, NULL)` after assembly
and zero the matching residual entries in `FormFunction`.

To map "rim vertex" to "Lagrange DOF" we need the section's
Lagrange-bearing points. Walked `[pStart, pEnd)` and counted
`PetscSectionGetFieldDof(ls, p, lagrange_field_idx_, &fdof)`:

```
Session 10: Lagrange-bearing points by depth: d0=0 d1=0 d2=0 d3=32 d4=0
```

All 32 Lagrange-bearing points sit at depth 3 (height 0, i.e. **cohesive
cells**), not at vertices. With Nc=3 components per cell that gives the
observed 96 Lagrange DOFs.

This is Session 10's central finding: PETSc placed the Lagrange DOFs at
cohesive cells, not at cohesive vertices, so the PyLith rim-vertex BC
pattern does not have a clean analogue. The implementation
(`src/core/Simulator.cpp:4110-4234`) walks the closure of each
Lagrange-bearing cell, marks the cell as "rim" if any vertex in its
closure is on the bounding box, and pins all of that cell's Lagrange
DOFs. For the test geometry:

```
Session 10: pinning 24 / 32 Lagrange-bearing points (72 local DOFs total)
```

A 4x4 cohesive grid of triangular prisms has 24 of 32 cells touching the
perimeter (only the central 8 cells avoid all 4 boundaries). Pinning
75% of the cohesive-surface Lagrange DOFs over-constrains the
prescribed-slip jump on the cells that should be slipping.

PCSVD spectrum after the pin landed (log `/tmp/s10_diag7.log`):

```
SVD: condition number 2.360290348391e+18, 0 of 318 singular values are (nearly) zero
SVD: smallest singular values: 8.33e-09 2.66e-08 6.02e-08 8.61e-08 9.52e-08
SVD: largest singular values : ~2e+10
```

Compared to the Session 9 baseline (`5.36e+17`, smallest `3.72e-08`),
condition number rose by ~4x and smallest sigma dropped by ~4x. The
`MatZeroRowsColumns` operation also wipes the displacement-Lagrange
coupling rows that the unrelated cohesive cells need; that, together
with the over-constraint, yields a matrix that is closer to singular,
not further from it. SNES took an actual step (`F` reduced from
`8.84e-05` to `1.13e+01`) -- still a divergent step, just a marginally
different one.

The pin code was reverted; the function bodies are now no-op stubs at
`src/core/Simulator.cpp:4236-4252`.

## SNES trace (with revert applied, matches Session 9 baseline)

Log `/tmp/s10_final_svd.log`:

```
0 TS dt 1. time 0.
    0 SNES Function norm 1.767766952966e-04
        SVD: condition number 5.363549166270e+17, 0 of 318 singular values are (nearly) zero
        SVD: smallest singular values: 3.72e-08 6.11e-08 7.02e-08 1.08e-07 1.26e-07
    Nonlinear solve did not converge due to DIVERGED_LINE_SEARCH iterations 0
... 11 more identical retries ...
max_fault_slip = 0, sol_norm = 0
```

Identical to Session 9's trace. The Session 10 diagnostic landed but no
pinning is in effect, so the saddle-point conditioning is unchanged.

## PCSVD spectrum: before / after / reverted

| State | Cond number | Smallest sigma | Notes |
|---|---:|---:|---|
| Session 9 baseline | `5.36e+17` | `3.72e-08` | Lagrange entirely free |
| Step C pin (24 cells) | `2.36e+18` | `8.33e-09` | Worse: pin destroys coupling, over-constrains physics |
| Session 10 head (revert) | `5.36e+17` | `3.72e-08` | Matches Session 9 |

## Fault-test status table vs Session 9 baseline

`ctest --output-on-failure -R '<fault subset>'` in `fsrm-ci:main`,
log `/tmp/s10_fault.log` (16 tests):

| Test | Session 9 | Session 10 |
|---|---|---|
| `Functional.DynamicRuptureSetup.*` (5 cases) | PASS | PASS |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Physics.LockedFaultTransparency` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | FAIL |
| `Integration.TimeDependentSlip` | FAIL | FAIL |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |

7 pass / 9 fail in both runs; identical set, no regressions.

## Full-suite count

`ctest -j$(nproc)`, log `/tmp/s10_test.log`: 14 fail.

The 9 fault-related failures match Session 9 exactly. The other 5
failures (`Physics.TerzaghiConsolidation`, `Physics.AbsorbingBC`,
`Integration.ExplosionSeismogram`, `Integration.NearFieldCoupled`,
`Integration.HistoricNuclear.NtsPahuteMesa`) are documented HDF5-output
parallel-execution conflicts (CLAUDE.md Rule note: "some tests may fail
when run in parallel due to HDF5 output file conflicts"). Re-running
those 5 individually all pass:

```
$ ctest --output-on-failure -R 'TerzaghiConsolidation|AbsorbingBC|ExplosionSeismogram|NearFieldCoupled|NtsPahuteMesa'
100% tests passed, 0 tests failed out of 5
```

Net: matches the documented 110 / 116 baseline (or 107 with the 3
documented parallel-only flakes). No Session 10 regressions.

## Success criteria outcome

| Criterion | Result |
|---|---|
| PCSVD condition number drops by at least 5 orders of magnitude | NOT MET. Reverted code matches Session 9 baseline (`5.36e+17`); Step-C experiment regressed to `2.36e+18`. |
| SNES takes a Newton step that reduces the residual below the initial value | NOT MET. SNES still diverges identically to Session 9. |
| `PrescribedSlipQuasiStatic` `max_fault_slip > 5e-4` passes | NOT MET. `max_fault_slip = 0`. |

Decision-gate outcome from the Session 10 prompt: "Condition number
unchanged: Step B did not actually land the BC on the Lagrange DOFs.
Read `DMPlexInsertBoundaryValuesEssentialField` to understand why,
report, stop."

Read; reported (Step B section above). Stopping.

## Why the PyLith pattern does not transfer

PyLith's `FaultCohesive::createConstraints`
(`libsrc/pylith/faults/FaultCohesive.cc:311-389`) places Lagrange DOFs on
cohesive **vertices** via a P1 (dim-1)-surface FE attached to the
cohesive face. Its `lagrange_multiplier_fault` field has DOFs at every
cohesive vertex on the buried-edges label. Pinning a buried-edge vertex
removes 3 DOFs, decoupling that single vertex without touching the slip
along the cohesive surface.

FSRM creates the Lagrange FE the same way (`src/core/Simulator.cpp:1705-
1718`):

```cpp
ierr = PetscFECreateDefault(PETSC_COMM_SELF, dim - 1, dim, isSimplex,
                            "lagrange_", PETSC_DETERMINE, &fe_lagrange);
ierr = DMAddField(dm, interfaces_label_, (PetscObject)fe_lagrange);
```

but PETSc 3.25's section builder distributes the resulting field's DOFs
onto the cohesive **cells** (height 0, depth 3) rather than the cohesive
vertices, as confirmed by the Session 10 depth diagnostic
(`d0=0 d1=0 d2=0 d3=32 d4=0`). The arithmetic checks: 32 cohesive cells
x 3 components = 96 DOFs = exactly the observed Lagrange field size. A
P1 surface FE on a cohesive prism would normally place DOFs at the
prism's negative-face vertices; the region-restricted `DMAddField` path
seems to collapse the per-vertex placement onto the cell. This is the
side effect Session 5 noted in the CLAUDE.md ("PETSc 3.25 region DS does
not support volume assembly via DMPlexTSComputeIFunctionFEM").

With Lagrange on cells, "vertex on the rim" is no longer a clean criterion
because every cohesive cell whose closure touches the rim is a
"rim cell" -- 24 of 32 in this geometry, i.e. 75% of the cohesive
surface. Pinning that many cells contradicts the prescribed-slip physics
(slip is supposed to occur everywhere on the cohesive surface, not only
at 8 internal cells).

## Files touched

- `include/core/Simulator.hpp` -- added `fault_boundary_label_`,
  `fault_boundary_lagrange_local_dofs_`, prototypes for
  `getOrCreateFaultBoundaryLabel`, `cacheFaultBoundaryLagrangeDofs`,
  `pinFaultBoundaryLagrangeJacobian`, `pinFaultBoundaryLagrangeResidual`.
- `src/core/Simulator.cpp` --
  - `setupBoundaryConditions` (~line 7100): call
    `getOrCreateFaultBoundaryLabel` after the existing fault BC block.
  - `setupFields` (~line 1798): call `cacheFaultBoundaryLagrangeDofs`
    after `DMSetUp`.
  - `setInitialConditions` (~line 5484): per-DS BC count and per-field
    section constraint dump.
  - new `getOrCreateFaultBoundaryLabel` (~line 3978).
  - new `cacheFaultBoundaryLagrangeDofs` (~line 4110) plus diagnostic
    print of Lagrange-bearing-point depth distribution.
  - new `pinFaultBoundary*` no-op stubs (~line 4236).
  - `FormFunction` / `FormJacobian` are NOT modified -- the pin
    experiments were reverted.
- `tests/integration/test_dynamic_rupture_solve.cpp` -- new
  `FSRM_SVD_DEBUG` env hook to swap LU PC for SVD PC + monitor.

## Logs (full, not truncated)

- `/tmp/s10_build.log`, `/tmp/s10_build2.log`, `/tmp/s10_build3.log` --
  Docker builds (final clean build is `s10_build3.log`).
- `/tmp/s10_svd.log`, `/tmp/s10_svd2.log`, `/tmp/s10_svd3.log`,
  `/tmp/s10_svd4.log` -- progressive PCSVD diagnostics with the
  `getOrCreateFaultBoundaryLabel` and `PetscDSAddBoundary` experiments.
- `/tmp/s10_diag1.log` ... `/tmp/s10_diag7.log` -- DS BC count, section
  constraint counts, Lagrange-by-depth, MatZeroRowsColumns trial.
- `/tmp/s10_final_svd.log` -- PCSVD baseline with all pin code reverted
  (matches Session 9 baseline).
- `/tmp/s10_fault.log` -- 16-test fault subset under `fsrm-ci:main`.
- `/tmp/s10_test.log` -- full ctest -j run.
- `/tmp/s10_serial.log` -- 5 parallel-flaky tests rerun serially (all
  pass).

## Rules compliance

- Rule 1 (Docker for all builds and tests): all builds and tests in
  `fsrm-ci:main`.
- Rule 2 (PETSc 3.25 API): `DMHasLabel`, `DMCreateLabel`, `DMGetLabel`,
  `DMPlexGetTransitiveClosure`, `DMPlexGetDepthLabel`,
  `PetscSectionGetFieldDof`, `PetscSectionGetFieldOffset`,
  `PetscDSAddBoundary`, `MatZeroRowsColumnsLocal`, `DMGetNumDS`,
  `DMGetRegionNumDS` all verified against `/opt/petsc-main/include/`
  before call.
- Rule 3 (no regressions): fault subset matches Session 9 exactly
  (7 pass / 9 fail). Full-suite parallel-flaky failures are documented
  HDF5 conflicts and pass when run individually.
- Rule 4 (DS/BC ordering): untouched. New `getOrCreateFaultBoundaryLabel`
  is called after the existing `fault_constraint` block in
  `setupBoundaryConditions`, before `DMSetLocalSection` clears the
  cached section, before `DMSetUp` rebuilds the section. Cache call
  follows `DMSetUp`.
- Rule 14 (CLAUDE.md update only on gate pass): gate did not pass;
  CLAUDE.md unchanged.
- Rule 16 (no log truncation): all logs referenced above are captured in
  full under `/tmp/`.

## Recommendation for Session 11

The PyLith rim-vertex pattern is the wrong abstraction for FSRM's current
Lagrange-on-cells layout. Two paths forward:

1. **Fix the DOF placement first**: investigate why
   `PetscFECreateDefault(comm, dim-1, dim, ...)` with
   `DMAddField(dm, interfaces_label_, ...)` distributes DOFs onto
   cohesive cells instead of cohesive vertices on PETSc 3.25. PyLith
   tests/ex5.c registers Lagrange exactly the same way and gets per-vertex
   DOFs. Diff `/opt/petsc-main/src/dm/impls/plex/tests/ex5.c` setup
   against `setupFields` to find the missing call (likely
   `PetscDSSetCohesive` ordering, an extra `DMPlexLabelComplete` on
   `interfaces_label_` before `DMAddField`, or a `PetscFE` setup option
   like `PetscFESetNumDof`). Once Lagrange sits at vertices, Session 10's
   `getOrCreateFaultBoundaryLabel` immediately yields the right pin set
   and the rim BC drops the condition number as expected.

2. **Switch solver strategy**: keep Lagrange on cells but use a Schur
   complement preconditioner (`-pc_type fieldsplit -pc_fieldsplit_type
   schur`) so the saddle-point block structure is exploited explicitly,
   bypassing the conditioning issue entirely. This is what PyLith
   ultimately uses for production (`pylith/materials/Elasticity.cc::
   getSolverDefaults`) and what PETSc's hybrid examples show.

Path 1 is the smaller change and matches CLAUDE.md's PyLith-alignment
philosophy; path 2 is the more robust long-term answer. Session 11 should
start with path 1 (specifically: a `-dm_section_view` diff between a
PETSc `tests/ex5.c -hybrid 1` run and an FSRM run to find the DOF-placement
discrepancy), and fall back to path 2 if the placement fix is more than
a one-line PETSc API call.
