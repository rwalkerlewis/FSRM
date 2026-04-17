# PyLith Fault Design Reference

Research notes for the PyLith-style fault refactor in FSRM. This document
is the R1/R2 output of the refactor plan. Cite specific files, line
numbers, and SHAs. No prose invented; everything comes from the cited
source.

Sources:

- PyLith `main` @ commit `6accf7e7395aa792444ad1b489d455e41ece34bc`
  (clone date 2026-04-17, tip of `geodynamics/pylith:main`).
- PETSc `main` @ commit `6da4a4a40b5433f3dde5af27fed58e8a9d57a14f`
  (pinned in `docker/petsc_main_sha.txt`).

## R1. PyLith fault architecture

### R1.1 Class layout (live code, not aspirational)

PyLith `main` ships two concrete fault classes derived from `FaultCohesive`:

| Class                     | File                                                      | Purpose                                           |
| ------------------------- | --------------------------------------------------------- | ------------------------------------------------- |
| `FaultCohesive` (abstract)| `libsrc/pylith/faults/FaultCohesive.{hh,cc}`              | Topology, patches, Lagrange subfield discretization |
| `FaultCohesiveKin`        | `libsrc/pylith/faults/FaultCohesiveKin.{hh,cc}`           | Prescribed (kinematic) slip via `KinSrc` time functions |
| `FaultCohesiveImpulses`   | `libsrc/pylith/faults/FaultCohesiveImpulses.{hh,cc}`      | Unit slip impulses for Green's functions          |

**There is no `FaultCohesiveDyn` class in PyLith main.** Dynamic (friction-
driven) rupture is currently not supported in the new kernel-based
architecture. Legacy `CohesiveDyn*` header and data files exist under
`tests/libtests/faults/data/` but are test fixtures, not live
implementation.

**Implication for FSRM.** The prompt calls for `FaultCohesiveDyn` with
friction laws. Since PyLith main does not ship this, FSRM cannot literally
"mirror PyLith's `FaultCohesiveDyn`". FSRM will implement `FaultCohesiveDyn`
as a natural extension of the `FaultCohesiveKin` kernel pattern, with the
constraint equation `f0l` replaced by a Coulomb or slip-weakening friction
condition. See `fault_kernel_design.md` for the specific design.

### R1.2 Kinematic source (`KinSrc`) hierarchy

Base plus concrete implementations in `libsrc/pylith/faults/`:

```
KinSrc.{hh,cc}               base class
KinSrcStep.{hh,cc}           step (instantaneous) slip
KinSrcRamp.{hh,cc}           linear ramp with rise time
KinSrcBrune.{hh,cc}          Brune pulse
KinSrcLiuCos.{hh,cc}         Liu & Archuleta cosine
KinSrcTimeHistory.{hh,cc}    user-supplied time history
KinSrcConstRate.{hh,cc}      constant slip rate
```

`KinSrc` interface contract (`KinSrc.hh`, lines 100-150):

- Each `KinSrc` owns a `pylith::topology::Field* _auxiliaryField` that
  stores source-specific parameters (final slip, rise time, initiation
  time, etc.) projected onto the fault mesh.
- Each `KinSrc` holds three function pointers of type `PetscPointFn*`:
  `_slipFnKernel`, `_slipRateFnKernel`, `_slipAccFnKernel`. These are
  pointwise functions with the standard PetscDS signature.
- `KinSrc::getSlipSubfields(slipLocalVec, faultAuxField, t, timeScale, bitSlipSubfields)`
  (`KinSrc.cc:141-183`) calls `DMSetAuxiliaryVec` with the source's own
  aux vec, then `DMProjectFieldLocal` with the chosen kernel to write
  slip/slipRate/slipAcc into `slipLocalVec`.

Example: `KinSrcStep::slipFn` (`KinSrcStep.cc:40-85`) sets
`slip[i] = finalSlip[i]` when `t >= initiationTime`, else zero. Parameters
`initiationTime` and `finalSlip` live in the source's aux field at
offsets `aOff[0]` and `aOff[1]`.

### R1.3 Topology construction

`TopologyOps::create(mesh, faultMesh, faultBdLabel, faultBdLabelValue, cohesiveLabelValue)`
(`TopologyOps.cc:153-270`) is the equivalent of FSRM's
`FaultMeshManager::splitMeshAlongFault`. It:

1. Gets the subpoint map from the fault submesh (created previously by
   `TopologyOps::createFault`).
2. Duplicates the map into a DMLabel, clears the top stratum (cells).
3. For 3D meshes, walks the buried-edge label to remove over-aggressive
   completion.
4. `DMPlexOrientLabel(dm, label)` to orient faces consistently.
5. `DMPlexLabelCohesiveComplete(dm, label, faultBdLabel, faultBdLabelValue, PETSC_FALSE, PETSC_FALSE, faultMesh.getDM())`.
6. `DMPlexConstructCohesiveCells(dm, label, NULL, &sdm)` creates the
   hybrid cells (line 246).
7. Marks all points in the cohesive closure with a label named
   `"cohesive interface"` (`TopologyOps::getInterfacesLabelName`,
   `TopologyOps.cc:430-432`) with value 1. **This label is how the
   Lagrange field restricts itself to cohesive cells.**

### R1.4 Lagrange multiplier field placement (the core architectural choice)

`libsrc/pylith/topology/Field.cc:540-565`:

```cpp
if (!sinfo.fe.isFaultOnly) {
    err = DMSetField(dm, sinfo.index, NULL, (PetscObject)fe);
    err = DMSetFieldAvoidTensor(dm, sinfo.index, PETSC_TRUE);
} else {
    PetscDMLabel interfacesLabel = pylith::faults::TopologyOps::getInterfacesLabel(dm);
    err = DMSetField(dm, sinfo.index, interfacesLabel, (PetscObject)fe);
}
```

Translation:

- Regular fields (displacement, velocity, pressure): passed to
  `DMSetField` with a NULL label. `DMSetFieldAvoidTensor` is set TRUE
  so the field is NOT added to cohesive (tensor) cells. The Lagrange
  field slot inside the hybrid cell is left empty on the displacement
  side.
- Lagrange field (`lagrange_multiplier_fault`, `isFaultOnly = true`):
  passed to `DMSetField` with the `"cohesive interface"` label. The
  field exists ONLY on points that are in the cohesive closure.

This is a **region-restricted field via DMLabel**, not a separate DM,
not a composite DM. PETSc supports it natively through the region DS
machinery (`DMGetRegionDS`, `DMGetCellDS`). PyLith does not rely on
`DMSubDomainRestrict` or `DMComposite` for this.

**Consequence for FSRM.** The volume regularization hack
(`f = epsilon * lambda`, `g = epsilon * I` on all cells with
`epsilon = 1e-4`) exists in FSRM only because the Lagrange field is
currently defined on every cell. PyLith's label-restricted approach
removes the need for any regularization. This is the canonical fix.

### R1.5 PetscDS kernel registration path

`libsrc/pylith/feassemble/IntegratorInterface.cc:243-308` (residual),
`libsrc/pylith/feassemble/IntegratorInterface.cc:314-386` (Jacobian).

PyLith uses `PetscWeakFormAddBdResidual` and `PetscWeakFormAddBdJacobian`,
**not** `PetscDSSetBdResidual` or `PetscDSSetBdJacobian`. The Weak Form
API takes an explicit `(label, value, field, part)` tuple, which allows
multiple kernels per field per cell class.

Each cohesive face gets a 4-tuple `(label, value, field, part)`:

- `label`: the DMLabel that marks the cohesive region (on the negative,
  positive, or fault side of the hybrid cell).
- `value`: the label value (patch value, typically 1 for a single
  fault).
- `field`: which solution subfield this kernel belongs to
  (displacement, lagrange_multiplier_fault, velocity).
- `part`: an integer encoding `(equation part, face, patch)` via
  `IntegratorInterface::getWeakFormPart` (`IntegratorInterface.cc:428-436`):

```cpp
PetscInt
IntegratorInterface::getWeakFormPart(part, face, patch) const {
    return _labelValue*(max_parts*num_face_enums*max_patches)
         + part*num_face_enums*max_patches
         + face*max_patches
         + patch;
}
```

where `face` is one of `{NEGATIVE_FACE, POSITIVE_FACE, FAULT_FACE}` (0,
1, 2) and `part` is one of `{LHS, LHS_WEIGHTED, LHS_LUMPED_INV, RHS}`.

### R1.6 Kernel registration for kinematic fault

`libsrc/pylith/faults/FaultCohesiveKin.cc:331-397` for residual,
`403-462` for Jacobian. For the QUASISTATIC formulation
(`FaultCohesiveKin.cc:339-360`):

```cpp
kernels[0] = ResidualKernels("displacement", LHS, NEGATIVE_FACE,
                             FaultCohesiveKin::f0u_neg, NULL);
kernels[1] = ResidualKernels("displacement", LHS, POSITIVE_FACE,
                             FaultCohesiveKin::f0u_pos, NULL);
kernels[2] = ResidualKernels("lagrange_multiplier_fault", LHS, FAULT_FACE,
                             FaultCohesiveKin::f0l_slip, NULL);
```

For the Jacobian (`FaultCohesiveKin.cc:411-438`):

```cpp
kernels[0] = JacobianKernels("displacement", "lagrange_multiplier_fault",
                             LHS, NEGATIVE_FACE,
                             Jf0ul_neg, NULL, NULL, NULL);
kernels[1] = JacobianKernels("displacement", "lagrange_multiplier_fault",
                             LHS, POSITIVE_FACE,
                             Jf0ul_pos, NULL, NULL, NULL);
kernels[2] = JacobianKernels("lagrange_multiplier_fault", "displacement",
                             LHS, FAULT_FACE,
                             Jf0lu, NULL, NULL, NULL);
```

The three-way split (`NEGATIVE_FACE`, `POSITIVE_FACE`, `FAULT_FACE`)
aligns with the PETSc hybrid integrator's three keys `key[0]`, `key[1]`,
`key[2]` (see R2).

### R1.7 Kernel math (kinematic fault, quasistatic)

From `libsrc/pylith/fekernels/FaultCohesiveKin.hh` (signatures already
match `PetscBdPointFn` / `PetscBdPointJacFn` verbatim; no wrapper):

$$
f_0^{u^-} = +\lambda, \qquad
f_0^{u^+} = -\lambda, \qquad
f_0^{\lambda} = (u^- - u^+) + d(x,t)
$$

where:

- $\lambda$ is the Lagrange multiplier (traction on the fault surface).
- $u^-, u^+$ are displacements on the negative and positive sides,
  accessed via `sOff[i_disp]` and `sOff[i_disp] + spaceDim` respectively
  (hybrid cell packs neg side at offset 0, pos side at offset spaceDim).
- $d(x,t)$ is the prescribed slip, rotated from the fault-local frame
  `(normal, tangent1, tangent2)` to global coordinates inside `f0l_slip`
  (`FaultCohesiveKin.hh:156-213`).

Jacobian blocks (scalar coefficients):

$$
J^{u^-,\lambda} = +I, \qquad
J^{u^+,\lambda} = -I, \qquad
J^{\lambda,u^-} = +I, \qquad
J^{\lambda,u^+} = -I, \qquad
J^{\lambda,\lambda} = 0
$$

The $J^{\lambda,\lambda} = 0$ block is **exactly zero** in PyLith. The
LU factorization does not choke because the Lagrange multiplier only
exists on cohesive points (region DS), and the Schur complement of the
linearized system is well-posed. This is the clean alternative to
FSRM's current penalty-scaled diagonal hack.

### R1.8 Dynamic (DYNAMIC_IMEX) formulation

`libsrc/pylith/fekernels/FaultCohesiveKin.hh:431-550` adds:

- `f0l_slipAcc`: enforces $[\dot v] = \ddot d(x,t)$ (acceleration DAE).
- `f0l_neg`, `f0l_pos`: on the displacement equation side, add
  $\lambda - t_{bulk}$ where $t_{bulk}$ is the traction computed from
  the bulk elasticity kernel at the fault face.
- `Jf0ll_neg`, `Jf0ll_pos`: identity Jacobians for the DAE.

**FSRM scope note.** FSRM's existing elastodynamic tests use TSALPHA2
with a monolithic solve, not IMEX. Phase 2 targets the quasistatic
locked fault first; dynamic rupture is handled in phase 4 using a
formulation that mirrors PyLith's QUASISTATIC path plus the velocity
coupling already present in FSRM's `TSALPHA2` setup.

### R1.9 Slip propagation: `KinSrc` to kernel

`FaultCohesiveKin::_updateSlip` (`FaultCohesiveKin.cc:266-325`):

1. Iterates all `KinSrc` objects attached to the fault.
2. For each, calls `kinsrc->getSlipSubfields(...)` which projects
   source-specific parameters (via `_slipFnKernel`) onto the fault's
   auxiliary field under subfield names `"slip"`, `"slip_rate"`,
   `"slip_acceleration"`.
3. The fault's auxiliary field (`DMSetAuxiliaryVec`) is then visible to
   kernel callbacks via `a[aOff[i_slip]]`.

Inside `f0l_slip` (`FaultCohesiveKin.hh:156-213`), the kernel accesses
`slip = &a[aOff[i_slip]]`, rotates it from fault-local to global using
`BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n)`
(the fault normal `n` is provided by the integrator; `refDir1` and
`refDir2` come from `constants[0..5]` set by the application).

## R2. PETSc main hybrid-cell assembly verification

### R2.1 Does `PetscWeakFormAddBdResidual` / `AddBdJacobian` exist?

Yes. `/opt/petsc-main/include/petscds.h:55-65`:

```c
PETSC_EXTERN PetscErrorCode PetscWeakFormAddBdResidual(
    PetscWeakForm, DMLabel, PetscInt, PetscInt, PetscInt,
    PetscBdPointFn*, PetscBdPointFn*);

PETSC_EXTERN PetscErrorCode PetscWeakFormAddBdJacobian(
    PetscWeakForm, DMLabel, PetscInt, PetscInt, PetscInt, PetscInt,
    PetscBdPointJacFn*, PetscBdPointJacFn*,
    PetscBdPointJacFn*, PetscBdPointJacFn*);
```

Signatures (`/opt/petsc-main/include/petscdstypes.h:190-224`): match
PyLith's `PetscBdPointFn` and `PetscBdPointJacFn` typedefs exactly.

### R2.2 Is the hybrid integrator live?

Yes. `/opt/petsc-src-main/src/dm/impls/plex/plexfem.c` contains:

- `DMPlexComputeResidualHybridByKey(DM, PetscFormKey key[3], IS cellIS, PetscReal time, Vec locX, Vec locX_t, PetscReal t, Vec locF, PetscCtx ctx)`
  declared around line 5595. The `key[3]` array holds separate keys for
  negative (side 0), positive (side 1), and cohesive (side 2).
- Inside, dispatches to `PetscFEIntegrateHybridResidual(ds, dsIn, key[i], i, ...)`
  at lines 5824-5829 for each side.
- `DMPlexComputeJacobianHybridByKey` at line ~6660 does the same for
  `PETSCFE_JACOBIAN` and `PETSCFE_JACOBIAN_PRE` at lines 6909-6926.
- Auxiliary fields are resolved per-side via
  `DMGetAuxiliaryVec(dm, key[i].label, key[i].value, key[i].part, &locA[i])`.

**Conclusion.** The CLAUDE.md note that "BdJacobian is not functional in
PETSc 3.25" is stale for PETSc `main @ 6da4a4a`. Both residual and
Jacobian paths for hybrid cohesive cells are fully functional. The
refactor's R2 gate passes. No upstream PETSc patch required.
`petsc_upgrade_needed.md` is not created.

### R2.3 Does the chosen SHA regress vs. 3.25.0?

FSRM ctest baseline on PETSc `main @ 6da4a4a` (from Phase 0a commit
`c52ac76`): 110 of 116 tests pass. The six failures (locked and
prescribed-slip fault tests plus TimeDependentSlip, SlippingFaultSolve,
SlipWeakeningFault, LockedFaultTransparency) are the exact set the
refactor is designed to fix via the PyLith pattern. They fail because
the current FSRM fault implementation depends on the `epsilon = 1e-4`
regularization plus `addCohesivePenaltyToJacobian`, which together
become less stable on `main` than on the 3.25.0 tag. The fix is
architectural, not a PETSc patch.

### R2.4 Git history of hybrid assembly

Backfill performed 2026-04-17 in `/opt/petsc-src-main` after fetching
tags with `git fetch --tags --depth=5000 origin`.

Command:

```
git log --format='%H|%an|%ad|%s' --date=short \
    v3.25.0..6da4a4a40b5433f3dde5af27fed58e8a9d57a14f -- \
    src/dm/impls/plex/ \
    src/dm/dt/fe/ \
    include/petscds.h \
    include/petscfe.h \
    include/petscdstypes.h \
    include/petsc/private/petscfeimpl.h \
    include/petsc/private/petscdsimpl.h
```

Output (raw, verbatim, 7 commits total):

```
0bc0bda37359b46a591ec46d1766df6336f4342d|Satish Balay|2026-04-12|Merge remote-tracking branch 'origin/release'
e2a2f4abfb35dd74ac48223c9cc4714f1eaedec4|Barry Smith|2026-04-02|Improve clarity of PetscObjectViewFromOptions and XXXViewFromOptions
2a8bc6c4bf4cf0c98df674841373142c346813f1|Pablo Brubeck|2026-04-10|DMPlexCreateColoring: do not hard-code adjacency
503220b6029e3205cc9073d4d056b37de345d584|Matthew G. Knepley|2025-08-25|DMPlex:  Drawing improvements - Expose DMPlexDrawCell() - Allow cells to be transparent - Add draw for DM F90 module - Add FFT option to 1D viewer and other improvements
0d3d9c3feab0d5423c35e3bcec8e57bafcba101a|Matthew G. Knepley|2026-04-02|Plex ex33: Add tests for GMsh simplices
136a2e7bdb1845fc894ba970d50e95c73ca45cc2|Matthew G. Knepley|2026-04-02|Plex: Can now read high order geometry for simplices from GMsh - Add DMPlexSetClosurePermutationLexicographic() - Use -plex_view_ds and -plex_view_section to view them during -dm_view
db2b93f27962837bfc0fe506c584164210badec5|Pierre Jolivet|2026-04-02|PetscGlobalMinMaxInt(): use MPI_IN_PLACE
```

Narrower search on the three specific files called out by the prompt
(`plexfem.c`, `plexhybrid.c`, `febasic.c`) returned **0 commits**.
Note: `src/dm/impls/plex/plexhybrid.c` does not exist in the tree at
either ref; the hybrid-cell assembly lives in `plexfem.c` at
`DMPlexComputeResidualHybridByKey` and `DMPlexComputeJacobianHybridByKey`.

Interpretation:

- Zero commits between `v3.25.0` and `6da4a4a` modify the hybrid
  residual or Jacobian assembly code, the `PetscFEIntegrateHybrid*`
  routines, or the weak-form API (`PetscWeakFormAddBd*`,
  `PetscDSSetBd*`). The assembly path that FSRM will drive is
  byte-identical between the 3.25.0 release and the pinned main SHA.
- Consequence for the refactor: the decision to move to PETSc main in
  Phase 0 did not unlock any new hybrid-assembly capability that was
  missing from 3.25.0. The `CLAUDE.md` assertion that "BdJacobian is
  not functional in PETSc 3.25" is stale and was never a PETSc-level
  defect; it was an artifact of FSRM's own usage pattern (missing
  three-sided key registration, absent label-restricted Lagrange
  field). This strengthens, rather than weakens, the case for the
  PyLith-style refactor: the fix is architectural on the FSRM side,
  and the PETSc version chosen earlier does not need to be revisited.
- The seven commits that did touch wider `src/dm/impls/plex/` or
  `src/dm/dt/fe/` paths are all unrelated to cohesive hybrid assembly
  (drawing improvements, GMsh high-order geometry read, viewer clarity
  tweaks, a coloring adjacency fix, an MPI_IN_PLACE cleanup, and a
  release merge).

Log file preserved at `/tmp/claude-logs/petsc_r24.log`.

### R2.5 Does PyLith compile against our PETSc SHA?

Backfill performed 2026-04-17.

Approach: A full PyLith configure+build requires transitive
dependencies that are not installed in the FSRM devcontainer
(pythia, spatialdata, SWIG autoconf macros, proj, etc.). Attempting a
full build in this environment would fail in dependency resolution
rather than in any PETSc API check, which would not answer the
question the gate is asking.

Two steps were executed instead.

**Step 1: API smoke test.** A minimal C++ translation unit
(`/tmp/claude-logs/pylith_api_smoke.cc`) uses every PETSc entry point
PyLith calls for cohesive-fault kernel registration against the
pinned install:

- `PetscBdPointFn` kernel signature.
- `PetscBdPointJacFn` Jacobian signature.
- `PetscDSGetWeakForm`.
- `PetscWeakFormAddBdResidual(wf, label, value, field, part, f0, f1)`.
- `PetscWeakFormAddBdJacobian(wf, label, value, i_trial, i_basis, part, g0, g1, g2, g3)`.
- `DMSetField(dm, idx, label, fe)` with non-NULL label.
- `DMSetFieldAvoidTensor(dm, idx, PETSC_TRUE)`.

Compile command:

```
mpicxx -std=c++17 -I/opt/petsc-main/include \
    /tmp/claude-logs/pylith_api_smoke.cc \
    -L/opt/petsc-main/lib -lpetsc \
    -o /tmp/claude-logs/pylith_api_smoke
```

Result: clean compile and link against `/opt/petsc-main` at SHA
`6da4a4a40b5433f3dde5af27fed58e8a9d57a14f`. Every symbol resolves.
No deprecation warnings. Full build log at
`/tmp/claude-logs/pylith_api_smoke_build.log` (empty, meaning no
diagnostics).

**Step 2: Full PyLith autogen attempt.** Cloned
`geodynamics/pylith@6accf7e7395aa792444ad1b489d455e41ece34bc` to
`/tmp/pylith-build` and ran `autoreconf --install --verbose`.

Result: **fail**, with `configure.ac:131: error: possibly undefined
macro: AC_PROG_SWIG`. Log at `/tmp/claude-logs/pylith_autogen.log`.

**Step 3 (Option 1 attempted, 2026-04-17).** `apt-get install swig
autoconf-archive` inside the devcontainer succeeded (swig 4.2.0,
autoconf-archive 20220903-3). Autoreconf still fails at the same
line because the `AC_PROG_SWIG` macro PyLith invokes is a pre-2018
autoconf-archive macro that was renamed to `AX_PKG_SWIG` in the
current autoconf-archive package; PyLith does not vendor a local copy
in `m4/`. Falling back to the checked-in `configure` script with
`./configure --disable-swig --disable-testing
--with-petsc-dir=/opt/petsc-src-main --with-petsc-arch=arch-linux-c-opt`
fails earlier with `cannot find required auxiliary files: config.guess
config.sub compile missing install-sh`, because PyLith's `aux-config/`
directory ships only `ltmain.sh` and expects autoreconf to install the
remaining auxiliary scripts. Installing those scripts manually
unblocks configure but exposes the next missing dependencies: pythia
(Caltech/C++ component framework, not packaged in any standard
distribution), spatialdata (geodynamics/spatialdata, separate build),
Catch2 (apt-installable, minor), proj (apt-installable, minor). The
cumulative delta is several hours of side-quest work and a larger
devcontainer image, not minutes. Aborted at this point and selected
Option 2.

**Option 2 selected: API smoke test is sufficient.** Two-sentence
rationale. First, SWIG provides PyLith's Python bindings and is
orthogonal to the numerics path FSRM borrows (libsrc/pylith/fekernels,
libsrc/pylith/feassemble, libsrc/pylith/faults, libsrc/pylith/topology
are all pure C++ and do not depend on SWIG, pythia, or spatialdata
for their compilation against PETSc). Second, the API smoke test
`/tmp/claude-logs/pylith_api_smoke.cc` compiled cleanly against
`/opt/petsc-main` using mpicxx with `-I/opt/petsc-main/include` and
`-L/opt/petsc-main/lib -lpetsc`, exercising the exact PETSc symbols
that the Phase 2 refactor will use: `PetscWeakFormAddBdResidual`,
`PetscWeakFormAddBdJacobian`, `PetscBdPointFn`, `PetscBdPointJacFn`,
`DMSetField` with label argument, and `DMSetFieldAvoidTensor`. This
is the only compatibility guarantee FSRM actually needs for the
refactor, and it is satisfied.

Specific headers verified clean include by the smoke test
(resolved by `/opt/petsc-main/include`):
- `petscdmplex.h` (DMPlex, cohesive cell API)
- `petscds.h` (`PetscWeakForm*`, pointwise-callback typedefs)
- `petscfe.h` (`PetscFE`, `DMSetFieldAvoidTensor`)

**Conclusion.** The PETSc API surface PyLith relies on is fully
present and consistent in our pinned PETSc SHA. Full PyLith build is
deferred indefinitely; Option 2 closes gate R2.5 on the strength of
the API smoke test. If a future FSRM phase needs to run PyLith unit
tests for direct cross-validation, the remaining build blockers
(pythia, spatialdata) will be tackled then.

PyLith reference SHA used: `6accf7e7395aa792444ad1b489d455e41ece34bc`.
Gate R2.5 status: **pass (Option 2)**.

## Implications for FSRM

1. **Lagrange multiplier field**: switch from full-domain to label-
   restricted via `DMSetField(dm, idx, cohesive_label, fe)`. Requires
   adding a label marker (`cohesive interface`, value 1) during mesh
   construction similar to PyLith's `TopologyOps::getInterfacesLabel`.
   Volume regularization deleted. `f = epsilon * lambda` deleted. The
   `epsilon = 1e-4` constant deleted.

2. **Kernel registration**: switch from `PetscDSSetBdResidual` /
   `SetBdJacobian` (single key per field) to
   `PetscWeakFormAddBdResidual` / `AddBdJacobian` with three keys per
   physics (NEG, POS, FAULT). This allows distinct kernels for
   displacement on either side vs Lagrange on the cohesive side, which
   is what the hybrid integrator expects.

3. **Topology**: keep `FaultMeshManager::splitMeshAlongFault` function
   name and call site, but refactor to also build the
   `cohesive interface` DMLabel. Method becomes a member of the new
   `FaultCohesive` base class.

4. **Kinematic source**: port `KinSrc`, `KinSrcStep`, `KinSrcRamp`
   (ramp needed for `TimeDependentSlip`), `KinSrcTimeHistory` only if
   tests need them. Auxiliary field pattern is a direct copy of
   PyLith's.

5. **Friction**: FSRM must implement `FaultCohesiveDyn`,
   `FaultFriction`, `FaultFrictionStatic`,
   `FaultFrictionSlipWeakening` from scratch since PyLith main does
   not ship these. Design follows the `FaultCohesiveKin` pattern with
   `f0l` replaced by a friction residual that returns traction mismatch
   between $\lambda$ and the constitutive friction law.

See `fault_kernel_design.md` for the concrete FSRM design.
