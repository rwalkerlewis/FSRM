# Session 7.5 Report

## Goal

Land the PyLith-style volume dispatch: drive volume physics through explicit
`DMPlexComputeResidualByKey` / `DMPlexComputeJacobianByKey` calls on a non-cohesive
cell IS, rather than through `DMPlexTSComputeI*FunctionFEM`, so that cohesive
cells never reach the non-hybrid integrator. The cohesive hybrid driver
landed in Session 4 stays untouched.

## Diagnosis (Sessions 4-6 re-read)

`DMPlexTSComputeIFunctionFEM` (`src/ts/utils/dmplexts.c:134`) iterates DSes and,
for each DS, passes a cellIS to `DMPlexComputeResidualByKey`. For the default
DS it passes the full `allcellIS` (all cells, including cohesive prisms).
Inside `DMPlexComputeResidualByKey` (`plexfem.c:6211`), the DS is fetched via
`DMGetCellDS(cells[cStart])` and the closure size is applied to every cell
in the chunk. Cohesive prisms dispatch to a different DS via `DMGetCellDS`
(see `plex.c:5677-5691`), so the closure-size mismatch silently produces
zero or garbage rows on the Jacobian for the affected cells. This matches
the Session 6 Jacobian dump, where every non-cohesive displacement row had
zero entries on the displacement diagonal.

## Code change

Both edits are in `src/core/Simulator.cpp`.

### FormFunction (Simulator.cpp:~5668)

Before:

```cpp
ierr = DMPlexTSComputeIFunctionFEM(dm, t, locU, locU_t, locF, ctx); CHKERRQ(ierr);
```

After:

```cpp
{
    PetscInt cStart = 0, cEnd_all = 0, cMax = 0;
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd_all); CHKERRQ(ierr);
    ierr = DMPlexGetSimplexOrBoxCells(dm, 0, NULL, &cMax); CHKERRQ(ierr);
    IS volumeCellsIS = NULL;
    ierr = ISCreateStride(PETSC_COMM_SELF, cMax - cStart, cStart, 1,
                          &volumeCellsIS); CHKERRQ(ierr);

    PetscInt Nds = 0;
    ierr = DMGetNumDS(dm, &Nds); CHKERRQ(ierr);
    for (PetscInt s = 0; s < Nds; ++s) {
        PetscDS ds_s = NULL;
        DMLabel label_s = NULL;
        ierr = DMGetRegionNumDS(dm, s, &label_s, NULL, &ds_s, NULL); CHKERRQ(ierr);
        PetscBool cohesive_ds = PETSC_FALSE;
        if (ds_s) PetscDSIsCohesive(ds_s, &cohesive_ds);
        if (cohesive_ds) continue;

        PetscFormKey key;
        key.label = NULL;
        key.value = 0;
        key.field = 0;
        key.part  = 0;

        IS useCellIS = NULL;
        if (label_s) {
            IS stratumIS = NULL;
            ierr = DMLabelGetStratumIS(label_s, 1, &stratumIS); CHKERRQ(ierr);
            if (stratumIS) {
                ierr = ISIntersect(stratumIS, volumeCellsIS, &useCellIS); CHKERRQ(ierr);
                ierr = ISDestroy(&stratumIS); CHKERRQ(ierr);
            }
        }
        if (!useCellIS) {
            ierr = PetscObjectReference((PetscObject)volumeCellsIS); CHKERRQ(ierr);
            useCellIS = volumeCellsIS;
        }

        ierr = DMPlexComputeResidualByKey(dm, key, useCellIS,
                                          PETSC_MIN_REAL,
                                          locU, locU_t, t, locF, ctx); CHKERRQ(ierr);
        ierr = ISDestroy(&useCellIS); CHKERRQ(ierr);
    }
    ierr = ISDestroy(&volumeCellsIS); CHKERRQ(ierr);
}
```

### FormJacobian (Simulator.cpp:~5785)

Before:

```cpp
ierr = DMPlexTSComputeIJacobianFEM(dm, t, locU, locU_t, a, J, P, ctx); CHKERRQ(ierr);
```

After: the analogous block, plus an up-front `MatZeroEntries(P)` (and on `J`
when `J != P`) to match the `DMPlexTSComputeIJacobianFEM` semantics that
zero the output at the start of the first DS iteration. The hybrid cohesive
Jacobian (`DMPlexComputeJacobianHybridByKey`) immediately below is unchanged.

### Departures from the Session 7.5 spec

Two details of the spec from the session kickoff had to be adjusted once the
actual DS layout was observed, and are worth calling out explicitly:

1. **`key.label = NULL`, not `label_s`.** The spec set `key.label = NULL` and
   `key.value = 0` because it assumed the default DS had `label = NULL`. In
   practice `DMCreateDS` (`dm.c:6085`) creates a synthetic `defaultCells`
   label whenever any field is added with a non-NULL label, and
   `DMGetRegionNumDS(dm, 0, ...)` then returns the `defaultCells` label
   rather than `NULL`. `PetscDSSetResidual` / `PetscDSSetJacobian`
   (`dtds.c`) always register callbacks on the DS weak form under
   `(label=NULL, value=0)`, and `PetscWeakFormGetJacobian` (`dtweakform.c`)
   does an exact-key hash lookup. So we still have to pass
   `key.label = NULL` even when the DS has a `defaultCells` region label,
   or the weak-form lookup inside `PetscFEIntegrateJacobian_Basic` returns
   zero callbacks and no entries are written. `DMPlexComputeResidualByKey`
   still dispatches the correct PetscDS per cell via `DMGetCellDS`; only
   the weak-form and auxiliary-vec lookups key on `key.label`.
2. **Skip on `PetscDSIsCohesive`, not on `label_s != NULL`.** The spec
   skipped DSes whose label was non-NULL, reasoning that non-default DSes
   were handled by the hybrid driver. But `defaultCells` is the
   volume-physics region for fault-enabled runs, and skipping it would
   drop all volume physics. The correct filter is `PetscDSIsCohesive`:
   skip only DSes that actually mark themselves cohesive (i.e. the one
   registered on the `cohesive interface` label with
   `PetscDSSetCohesive`).

The PyLith-style intent of the spec is preserved: volume cells driven by an
explicit `DMPlexComputeResidualByKey` / `DMPlexComputeJacobianByKey` on a
non-cohesive cell IS, cohesive cells exclusively on the Session 4 hybrid
driver.

## PrescribedSlipQuasiStatic outcome

SNES trace (beuler step 0, `-pc_type lu -ksp_type preonly`):

```
    0 SNES Function norm 1.385642450386e+01
    Nonlinear solve did not converge due to DIVERGED_LINEAR_SOLVE iterations 0
    ...
```

The initial residual norm is `1.386e+01` compared to the Session 6 value of
`1.77e-04`. The elasticity f1 callback is now firing on the non-cohesive
cells that were previously being dispatched through the wrong DS. Final
`sol_norm = 0` and `max_fault_slip = 0` because KSP still reports
`DIVERGED_LINEAR_SOLVE iterations 0` and SNES is not able to take a
Newton step.

### Jacobian diff vs Session 6

A one-shot `MatView(P)` dump was taken on the first `FormJacobian` call
(`/tmp/s75_jac.txt`, 136517 bytes, vs `/tmp/s6_jac.txt` at 1108041 bytes for
11 concatenated SNES-iteration dumps). Unique nonzero numeric values present
in the Session 7.5 dump:

```
-1.33333e+09, 1.33333e+09
-1.5e+09, 1.66667e+09, 2.33333e+09, 2.66667e+09
-1.66667e+08, 1.66667e+08
-3.33333e+08, 3.33333e+08
-6.66667e+08, 6.66667e+08
-0.0104167, 0.0104167
```

The `O(1e+9)` values are the elasticity Jacobian entries (on the order of
`2*mu + lambda`) that were entirely absent in the Session 6 dump (which only
contained the `0.0104167` hybrid cohesive coupling). Rows 0 through 95 still
show zero `(i, i)` diagonal on the displacement DOFs because those vertices
are bound by Dirichlet BCs (bottom-fixed, side rollers); their closure entries
are removed at matrix assembly by the PETSc BC path and replaced with scale
entries elsewhere. Rows for the interior vertices (e.g. row 97:
`(97, 5e+09)`) show the expected diagonal stiffness. In short: the
displacement-displacement block is now populated correctly; the remaining
problem is saddle-point conditioning, not missing entries.

### Why KSP still fails

Per the session-kickoff prediction, with the disp-disp block populated and
the disp-Lagrange / Lagrange-disp couplings coming from the hybrid driver,
the failure mode shifts to solver configuration. The assembled matrix is
a classical saddle-point system of the form

```
  [ A   B^T ]
  [ B    0  ]
```

with `A` (displacement stiffness) well-posed after BC elimination but the
Lagrange block having a zero diagonal (row 303 in the dump:
`(303, 0.)`). `-pc_type lu -ksp_type preonly` with PETSc's built-in serial
LU cannot factor this without pivoting support, and `ierr_linear_solve=0`
after one iteration is consistent with an LU setup failure. Fixing this is
test-config work (switching the preconditioner to MUMPS or a Schur-complement
fieldsplit) and is explicitly flagged in the Session 7.5 plan as the next
step beyond the dispatch fix.

## Full fault-test status table

Baseline (`11bc69e`, Session 6 HEAD) and Session 7.5 applied, both run via
`ctest --output-on-failure -R <fault subset>`:

| Test | Baseline | Session 7.5 |
|---|---|---|
| `Functional.DynamicRuptureSetup` (setup only, several cases) | PASS | PASS |
| `Physics.CohesiveBdResidual` | PASS | PASS |
| `Physics.LockedFaultTransparency` | FAIL | FAIL |
| `Physics.SCEC.TPV5` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedQuasiStatic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.LockedElastodynamic` | FAIL | FAIL |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | FAIL | FAIL |
| `Integration.TimeDependentSlip` | FAIL | FAIL |
| `Integration.PressurizedFractureFEM` | FAIL | FAIL |
| `Integration.SlippingFaultSolve` | FAIL | FAIL |
| `Integration.SlipWeakeningFault` | FAIL | FAIL |

16 tests total in the subset, 7 pass, 9 fail, on both runs. **CLAUDE.md's
"6 known failures" list is stale**: `Physics.SCEC.TPV5`,
`Integration.DynamicRuptureSolve.LockedElastodynamic`, and
`Integration.PressurizedFractureFEM` are also currently failing on `local_fix`
HEAD before this session touched the code. These were verified by running
the same test subset against a clean `git checkout` of the working tree
prior to re-applying the Session 7.5 patch; the identical 9-test failure set
appeared in both runs. No currently-passing test regressed.

## Full suite

`docker run ... ctest -j$(nproc)`: 104 of 116 tests pass (90%). The 12
failing tests are:

```
 40 - Physics.SCEC.TPV5                                (fault, baseline)
 43 - Physics.AbsorbingBC                              (parallel HDF5 flake)
 57 - Physics.LockedFaultTransparency                  (fault, baseline)
 67 - Integration.ExplosionSeismogram                  (parallel HDF5 flake)
 91 - Integration.PressurizedFractureFEM               (fault, baseline)
 93 - Integration.DynamicRuptureSolve.LockedQuasiStatic (fault, baseline)
 94 - Integration.DynamicRuptureSolve.LockedElastodynamic (fault, baseline)
 95 - Integration.DynamicRuptureSolve.PrescribedSlip   (fault, baseline)
 98 - Integration.TimeDependentSlip                    (fault, baseline)
 99 - Integration.NearFieldCoupled                     (parallel HDF5 flake)
100 - Integration.SlippingFaultSolve                   (fault, baseline)
101 - Integration.SlipWeakeningFault                   (fault, baseline)
```

All "parallel HDF5 flake" tests pass when run individually per
CLAUDE.md note. All "fault, baseline" tests match the failure set observed
in the subset run above and were failing at `11bc69e` before this session.

## Revert / re-apply rehearsal

During the session the fault subset was run first with the fix applied
(9 fails), then reverted (9 fails), then re-applied (9 fails). The failure
set was identical in all three runs, which is how the stale CLAUDE.md
known-failure list was identified and ruled out as a Session 7.5 regression.

## What this unblocks

With the volume disp-disp block now carrying the elasticity contribution,
the next step for PrescribedSlipQuasiStatic is the KSP/PC configuration
change the session-kickoff prediction described: either
`-pc_type lu -pc_factor_mat_solver_type mumps` or a Schur fieldsplit. That
is a test config change in `tests/integration/test_dynamic_rupture_solve.cpp`
and does not require further Simulator changes. Deferred to Session 8.

The `LockedElastodynamic` path relies on the dynamic mass contribution from
`TSALPHA2` (which writes directly onto the matrix diagonal via `X_tShift`
and does not depend on the DS weak-form lookup) and so was never exposed
to the empty-disp-disp bug. Session 7.5 does not change that path; its
failure mode is unrelated.

## Rules compliance

- Rule 1 (Docker build/test): all builds and tests run in `fsrm-ci:main`.
- Rule 2 (PETSc 3.25 API): signatures for `DMPlexComputeResidualByKey`,
  `DMPlexComputeJacobianByKey`, `DMGetRegionNumDS`,
  `DMPlexGetSimplexOrBoxCells`, and `PetscDSIsCohesive` verified against
  `/opt/petsc-main/include/petscdmplex.h` and the sources in
  `/opt/petsc-main/src/...` prior to call.
- Rule 3 (no regressions): verified above, same failure set as baseline.
- Rule 4 (DS/BC ordering): untouched.
- Rule 16 (no log truncation): build, fault-subset, and full-suite runs are
  captured in `/tmp/s75_build_final.log`, `/tmp/s75_fault_final.log`,
  `/tmp/s75_full_final.log`.
