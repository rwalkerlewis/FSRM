# Fault Test Regression Audit

Session 1.5 diagnostic and baseline audit. No code changes. Produced from a
linear walk back along `git log --first-parent local_fix` of 10 commits plus
the baseline tip, each rebuilt in Docker image `fsrm-ci:local` and tested
through `ctest --output-on-failure -R` against the six fault tests.

Commands used:

```
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'cd build && make -j$(nproc) && ctest -j$(nproc) --output-on-failure' \
  > /tmp/baseline_tip.log 2>&1

# Per-commit bisect, capturing /tmp/bisect_<sha>.log each iteration
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  "cd build && make -j$(nproc) && ctest --output-on-failure -R \
   'Physics.LockedFaultTransparency|Integration.DynamicRuptureSolve.LockedQuasiStatic|Integration.DynamicRuptureSolve.PrescribedSlip|Integration.TimeDependentSlip|Integration.SlippingFaultSolve|Integration.SlipWeakeningFault'"
```

Log files retained: `/tmp/baseline_tip.log`, `/tmp/bisect_<sha>.log` for each
commit, and `/tmp/bisect_driver.log` for the driver output. None were
truncated per rule 16 of CLAUDE.md.

## Executive summary

All six fault tests FAIL at every commit walked back in this audit. There is
no last-known-good commit within the 10-commit window. The regression that
the commit messages of `31ffb29c` and `7bc858ea` refer to with
"115 pass, 1 fail" cannot be reproduced in the current environment.

The single environmental variable most likely to explain the difference is
the PETSc version. The current `fsrm-ci:local` image carries PETSc 3.22.2
(verified inside the container and in `build/CMakeCache.txt`), while
`docs/CLAUDE.md` and `Dockerfile.ci` target PETSc 3.25.0.

## 1. Baseline at tip (commit 8c1ebbf)

Build + full `ctest -j$(nproc)`, captured at `/tmp/baseline_tip.log`:

```
89% tests passed, 13 tests failed out of 116

The following tests FAILED:
   42 - Physics.TerzaghiConsolidation (Failed)
   43 - Physics.AbsorbingBC (Failed)
   57 - Physics.LockedFaultTransparency (Failed)
   67 - Integration.ExplosionSeismogram (Failed)
   93 - Integration.DynamicRuptureSolve.LockedQuasiStatic (Failed)
   95 - Integration.DynamicRuptureSolve.PrescribedSlip (Failed)
   97 - Integration.TractionBC (Failed)
   98 - Integration.TimeDependentSlip (Failed)
   99 - Integration.NearFieldCoupled (Failed)
  100 - Integration.SlippingFaultSolve (Failed)
  101 - Integration.SlipWeakeningFault (Failed)
  110 - Integration.HistoricNuclear.DegelenMountain (Failed)
  111 - Integration.HistoricNuclear.NtsPahuteMesa (Failed)
```

All six nominated fault tests (57, 93, 95, 98, 100, 101) are in the failure
set. The remaining seven failures are not in the scope of this audit.
CLAUDE.md documents several of them (HDF5 parallel conflicts, absorbing BC,
TerzaghiConsolidation, ExplosionSeismogram, TractionBC, NearFieldCoupled,
HistoricNuclear) as flaky under `ctest -j` or as known.

## 2. Per-commit bisect table

Ten first-parent commits walked back from the tip. Each was checked out,
rebuilt in Docker, and the six fault tests were run via `ctest -R`. Column
legend: `F-slow` indicates the test ran to completion or iteration limit and
returned a zero-slip / zero-displacement solution (ctest times tens of
seconds to over a minute per test); `F-fast` indicates the test crashed or
failed its nonzero-solution check in a fraction of a second. `F-fast`
corresponds to a `PETSC ERROR: Field 1 is not in [0, 1)` error reported in
`/tmp/bisect_be3a2b32.log` for the experimental commit.

| # | SHA | Date (UTC) | Author | Commit subject | Total test time | LT | LQS | PS | TDS | SFS | SWF |
|---|-----|------------|--------|----------------|-----------------|----|----|----|-----|-----|-----|
| 0 | 8c1ebbf | 2026-04-19 02:27 | Robert Walker | Section B Option B partial: interior Lagrange manual assembly | baseline full suite | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 1 | 71064adf | 2026-04-19 00:22 | Robert Walker | add check md (docs only) | 129.68 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 2 | 85f872c1 | 2026-04-18 23:45 | Robert Walker | Added fault gap analysis (docs only) | 137.10 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 3 | 6abdde5f | 2026-04-17 05:08 | Robert Walker | Remove Windows-reserved filename archive/deploy/deploy/gcp/nul (no code) | 128.96 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 4 | 7313da14 | 2026-04-16 19:49 | Cursor Agent | Revert "experiment: register Lagrange field with cohesive interface label" | 130.92 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 5 | 4526b494 | 2026-04-16 19:49 | Cursor Agent | defense: skip Lagrange BC when field is absent from default DS | 7.27 s | F-fast | F-fast | F-fast | F-fast | F-fast | F-fast |
| 6 | be3a2b32 | 2026-04-16 19:04 | Cursor Agent | experiment: register Lagrange field with cohesive interface label | 4.91 s | F-fast | F-fast | F-fast | F-fast | F-fast | F-fast |
| 7 | 8c8ae541 | 2026-04-16 17:45 | Cursor Agent | Call getOrCreateInterfacesLabel at end of setupFaultNetwork | 137.10 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 8 | 33aff547 | 2026-04-16 14:19 | Cursor Agent | feat: add getOrCreateInterfacesLabel method (unused, PyLith TopologyOps.cc port) | 141.31 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 9 | 31ffb29c | 2026-04-16 08:26 | Cursor Agent | docs: update CLAUDE.md and TEST_RESULTS.md for 115 pass, 1 fail | 140.89 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |
| 10 | d271f13f | 2026-04-16 08:18 | Cursor Agent | test: replace explosion+fault GTEST_SKIP with actual test body + segfault guard | 133.09 s | F-slow | F-slow | F-slow | F-slow | F-slow | F-slow |

Column keys: LT = `Physics.LockedFaultTransparency`, LQS =
`Integration.DynamicRuptureSolve.LockedQuasiStatic`, PS =
`Integration.DynamicRuptureSolve.PrescribedSlip`, TDS =
`Integration.TimeDependentSlip`, SFS = `Integration.SlippingFaultSolve`,
SWF = `Integration.SlipWeakeningFault`.

Note: commit `d271f13f` only touched `tests/integration/test_explosion_fault_reactivation.cpp`
(not one of the six audited tests) and did not modify `src/` or
`include/`. Its fault-test outcome is structurally consistent with its
parent on these six tests, confirmed by the 133.09 s F-slow result.

## 3. Last-known-good

There is no last-known-good commit within the 10-commit window.

All 10 commits (and the tip) fail all six fault tests in the current
`fsrm-ci:local` Docker environment. The commit-message claims of
"115 pass, 1 fail" at `31ffb29c` and `7bc858ea` are not reproducible here.
The most parsimonious explanation is environmental: the installed PETSc
is version 3.22.2, not the 3.25.0 that CLAUDE.md and Dockerfile.ci target.

```
$ docker run --rm fsrm-ci:local ls /opt/
petsc-3.22.2

$ grep petsc build/CMakeCache.txt | head
pkgcfg_lib_PETSC_petsc:FILEPATH=/opt/petsc-3.22.2/arch-linux-c-opt/lib/libpetsc.so
```

Two regression-mode transitions are visible within the window. Both transitions
change the failure mode but not the pass/fail verdict for the six tests.
Because all commits fail, neither transition identifies a regression
introduction point.

1. `be3a2b32` and `4526b494` exhibit `F-fast` failures (`PETSC ERROR:
   Field 1 is not in [0, 1)` at the `DMAddBoundary` call for the Lagrange
   field, surfaced by the experimental `DMAddField(dm, interfacesLabel,
   fe_lagrange)` that placed the Lagrange field on a region-specific DS
   absent from the default DS). These two commits straddle the revert at
   `7313da14`.

2. All other commits (`8c8ae541` through tip, and `33aff547`, `31ffb29c`)
   exhibit `F-slow` failures: the SNES solve completes but returns a
   zero-slip solution, tripping the test-side nonzero-solution assertion.

## 4. Environment discrepancy

`docs/CLAUDE.md` section "Build Environment" states "PETSc 3.25.0 with
ctetgen". `Dockerfile.ci` at tip also targets 3.25.0 (`wget
petsc-3.25.0.tar.gz`). The locally-tagged image `fsrm-ci:local` currently
contains a build of PETSc 3.22.2. The CMake cache generated by that image
points all PETSc paths at `/opt/petsc-3.22.2`. This image was presumably
built from an earlier version of `Dockerfile.ci` and was not rebuilt after
the upgrade.

Rebuild path that is cheapest and most likely to recover the "115 pass"
claim:

```
docker build -f Dockerfile.ci -t fsrm-ci:local .
rm -rf build && mkdir build
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON \
   -DENABLE_CUDA=OFF && make -j$(nproc)'
```

Re-running the bisect after this rebuild is the minimum action required to
identify a real last-known-good commit (if one exists on this branch).
`docs/PYLITH_COMPATIBILITY.md` Section B "Option A" and
`docs/LAGRANGE_FIX_STATUS.md` also list a PETSc main-branch upgrade as the
independently-motivated path forward for the region DS assembly fix that
`PrescribedSlipQuasiStatic` requires.

## PrescribedSlipQuasiStatic Label Topology Diagnosis

This section re-diagnoses `Integration.DynamicRuptureSolve.PrescribedSlip`
(specifically `DynamicRuptureSolveTest.PrescribedSlipQuasiStatic`) against
the finding in `docs/LAGRANGE_FIX_STATUS.md` under "Outstanding":

> The fault label registered for the Lagrange-field `DM_BC_NATURAL`
> boundary contains only depth less than dim points, so
> `DMPlexComputeBdResidual` has no cohesive cells to integrate over for
> that field.

### 1. Fault label creation and DMAddBoundary in FSRM

`FaultMeshManager::createPlanarFaultLabel` at
`src/numerics/FaultMeshManager.cpp:52-145` creates the label named `"fault"`:

- Line 81-82: `DMCreateLabel(dm, "fault"); DMGetLabel(dm, "fault", fault_label);`
- Line 85-86: iterates over the `height = 1` stratum only
  (`DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd)`).
- Line 121: `DMLabelSetValue(*fault_label, f, 1)` marks each face that lies
  on the fault plane.

This is a face-only label. No vertex, edge, or cell closure is added. There
is no call to `DMPlexLabelComplete` on this label anywhere in FSRM (a
temporary `DMLabelDuplicate` + `DMPlexLabelComplete` is used inside
`splitMeshAlongFault` at `src/numerics/FaultMeshManager.cpp:176-177`, but
that duplicate is destroyed immediately at line 182-209 and is not the label
passed to `DMAddBoundary` later).

The `DMAddBoundary` call that registers the Lagrange-field `DM_BC_NATURAL`
boundary in FSRM is at `src/core/Simulator.cpp:6468-6481`:

```cpp
ierr = DMAddBoundary(dm, DM_BC_NATURAL, "fault_traction",
    fault_label, 1, &label_value, displacement_field, 0, NULL,
    NULL, NULL, NULL, NULL); CHKERRQ(ierr);

if (lagrange_field_idx < numDSFields) {
    ierr = DMAddBoundary(dm, DM_BC_NATURAL, "fault_constraint",
        fault_label, 1, &label_value, lagrange_field_idx, 0, NULL,
        NULL, NULL, NULL, NULL); CHKERRQ(ierr);
} else if (rank == 0) {
    PetscPrintf(comm,
        "  Skipped fault_constraint BC: Lagrange field %d not in default DS "
        "(numDSFields=%d; likely on region-specific DS)\n",
        (int)lagrange_field_idx, (int)numDSFields);
}
```

`fault_label` here is the same face-only label from `createPlanarFaultLabel`.
Directly after the cohesive split in `setupFaultNetwork`, the label's height
stratum of 1 refers to points that on the PRE-SPLIT DM were faces
(`dim - 1 = 2` for a 3D simplex mesh). After `DMPlexConstructCohesiveCells`
the point numbering is remapped by PETSc, and the label tracks the new
numbering as the DM travels through the split. What it does NOT track are
the cohesive prism cells, which live at height 0 (not height 1), and the
cohesive vertices at height `dim`.

Diagnostic block at `src/core/Simulator.cpp:6483-6536` confirms this. It
iterates the fault-label stratum and classifies each point by its depth:

```
cells=N_cells, faces=N_faces, edges=N_edges, vertices=N_vertices
```

with the printed result under the runtime log showing only `faces > 0` in
the `label_value == 1` stratum. That matches what the LAGRANGE_FIX_STATUS
outstanding item claims.

### 2. PyLith reference: per-stratum label via `DMPlexGetSimplexOrBoxCells`

`/home/dockimble/Projects/pylith/libsrc/pylith/faults/TopologyOps.cc:438-467`:

```cpp
PetscDMLabel pylith::faults::TopologyOps::getInterfacesLabel(PetscDM dm) {
    PYLITH_METHOD_BEGIN;
    PetscErrorCode err;
    PetscDMLabel interfacesLabel = NULL;

    const char* interfacesLabelName = TopologyOps::getInterfacesLabelName();
    PetscBool hasInterfacesLabel = PETSC_FALSE;
    if (DMHasLabel(dm, interfacesLabelName, &hasInterfacesLabel)) {
        err = DMGetLabel(dm, interfacesLabelName, &interfacesLabel); ...
    } else {
        PetscInt dim = 0;
        PetscInt pStart = 0, pEnd = 0, pMax = 0;
        err = DMGetDimension(dm, &dim); ...
        err = DMCreateLabel(dm, interfacesLabelName); ...
        err = DMGetLabel(dm, interfacesLabelName, &interfacesLabel); ...
        for (PylithInt iDim = 0; iDim <= dim; ++iDim) {
            err = DMPlexGetHeightStratum(dm, iDim, &pStart, &pEnd); ...
            err = DMPlexGetSimplexOrBoxCells(dm, iDim, NULL, &pMax); ...
            for (PylithInt p = pMax; p < pEnd; ++p) {
                err = DMLabelSetValue(interfacesLabel, p, 1); ...
            }
        }
    }
    assert(interfacesLabel);
    PYLITH_METHOD_RETURN(interfacesLabel);
}
```

The pattern walks every topological dimension `iDim` from 0 to `dim`. At
each height stratum it calls `DMPlexGetSimplexOrBoxCells` to obtain `pMax`,
which is the first point at that stratum that belongs to a hybrid
(cohesive) cell. Every point from `pMax` to `pEnd` is marked with value 1.
The resulting label covers hybrid cohesive cells at depth `dim`, cohesive
faces at depth `dim - 1`, cohesive edges at depth 1, and cohesive vertices
at depth 0. It is topologically complete over the cohesive interface and
is exactly the label that PETSc's `DMPlexComputeBdResidual` expects for a
`DM_BC_NATURAL` boundary on a field restricted to the cohesive region.

### 3. Does FSRM use the per-stratum pattern for the Lagrange `DMAddBoundary`?

No. The per-stratum pattern is present in FSRM, but its output is not used.

- `Simulator::getOrCreateInterfacesLabel` at
  `src/core/Simulator.cpp:3625-3664` is a direct port of PyLith's
  `TopologyOps::getInterfacesLabel`. It creates a label named
  `"cohesive interface"` via the identical `DMPlexGetSimplexOrBoxCells`
  per-stratum loop (`src/core/Simulator.cpp:3643-3650`).

- `Simulator::setupFaultNetwork` calls it at
  `src/core/Simulator.cpp:3558-3560`:

  ```cpp
  DMLabel interfacesLabel = NULL;
  ierr = getOrCreateInterfacesLabel(&interfacesLabel); CHKERRQ(ierr);
  (void)interfacesLabel;
  ```

  The return value is immediately discarded with `(void)interfacesLabel;`.
  The `"cohesive interface"` label is created on the DM but it is never
  passed to `DMAddBoundary`, to `DMAddField`, or to any other DM operation
  that would consume it.

- The only `DMAddBoundary` that targets the Lagrange field is the one at
  `src/core/Simulator.cpp:6472-6475`, and it passes `fault_label`, which
  is the face-only label from `FaultMeshManager::createPlanarFaultLabel`.

So FSRM uses a face-only label for the Lagrange-field `DM_BC_NATURAL`
boundary, not a per-stratum cohesive-complete label. The LAGRANGE_FIX_STATUS
outstanding item is factually correct: the label registered for that
boundary contains only lower-depth points (faces pre-split, remapped
post-split but still not the cohesive cells themselves), and
`DMPlexComputeBdResidual` has no cohesive cells at height 0 to integrate
over for the Lagrange field.

History leading to the current state: commit `33aff547` added
`getOrCreateInterfacesLabel` unused; commit `8c8ae541` called it but
discarded the return value; commit `be3a2b32` wired it to
`DMAddField(dm, interfacesLabel, fe_lagrange)`, which crashed with
`PETSC ERROR: Field 1 is not in [0, 1)` at `DMAddBoundary` because the
Lagrange field was then on a region-specific DS and absent from the default
DS; commit `4526b494` added a defensive skip of the `DMAddBoundary` when
the Lagrange field was absent, avoiding the crash but still not driving
the cohesive BdResidual; commit `7313da14` reverted the experimental
`DMAddField` call, putting the Lagrange field back on the whole mesh but
leaving the defensive skip and the unused `getOrCreateInterfacesLabel`
in place.

### 4. Does Section B of docs/PYLITH_COMPATIBILITY.md need revision?

Yes.

Section B currently frames the prescribed-slip failure as "an augmented
Lagrangian scaling problem, not a physics problem", attributing it to the
penalty-damped diagonal in `addCohesivePenaltyToJacobian` and concluding
that Option B (manual Lagrange assembly bypassing PetscDS entirely) is the
correct fix because it removes the penalty damping and interior volume
regularization competition.

The label-topology diagnosis above identifies a distinct, upstream reason
the PetscDS BdResidual cannot drive the Lagrange multiplier to a nonzero
traction for prescribed slip: the `DMAddBoundary` that registers the
Lagrange-field natural BC is attached to a face-only `fault_label` that
contains no height-0 cohesive cells, so `DMPlexComputeBdResidual`
integrates over zero cohesive geometry. Whatever BdResidual callback is
registered on the Lagrange field is never invoked on a cohesive cell.

A revision of Section B would say, in addition to the existing content:

> There is a second, upstream failure mode for `PrescribedSlipQuasiStatic`
> that is independent of the penalty-damping argument. The label passed to
> `DMAddBoundary(dm, DM_BC_NATURAL, "fault_constraint", fault_label, ...,
> lagrange_field, ...)` at `src/core/Simulator.cpp:6472-6475` is the
> face-only `"fault"` label produced by
> `FaultMeshManager::createPlanarFaultLabel` at
> `src/numerics/FaultMeshManager.cpp:52-145`. This label contains no
> cohesive cells (height 0) and no cohesive vertices (height `dim`). PETSc
> `DMPlexComputeBdResidual` iterates over cohesive cells that appear in the
> label. Because no cohesive cells appear, the BdResidual kernel registered
> via `PetscDSSetBdResidual` is not invoked on any cohesive cell. Option B
> (full manual Lagrange assembly) sidesteps this by not using BdResidual at
> all, but a simpler intermediate fix is to replace `fault_label` in that
> `DMAddBoundary` call with the label produced by
> `Simulator::getOrCreateInterfacesLabel` at
> `src/core/Simulator.cpp:3625-3664`, which is already populated using the
> PyLith per-stratum `DMPlexGetSimplexOrBoxCells` pattern and is currently
> discarded unused at `src/core/Simulator.cpp:3558-3560`. This label
> includes cohesive cells at height 0. If the revised call succeeds in
> PETSc 3.25, the BdResidual driver path is restored without the
> architectural change Option B requires.

This revision is not written into `docs/PYLITH_COMPATIBILITY.md` in this
session. Session 1.5 is diagnostic only.
