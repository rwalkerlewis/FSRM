# PyLith Reference (verified)

Every fact here is pinned to a specific file:line in the PyLith source tree
recorded below. When a future FSRM session re-derives a fact about PyLith,
the answer goes here so the next session does not re-derive it.

For a feature-parity wishlist (what FSRM aspires to port) see
`docs/PYLITH_COMPATIBILITY.md`. This document is a verified reference, not a
wishlist.

## PyLith source location and pin

- Path: `/home/dockimble/pylith-dev/src/pylith_official`
- Commit SHA: `52b3012a0f87992544059917a3cfedee12e7d02e`
- Branch: (as checked out during Session 31)

All file paths below are relative to `libsrc/pylith/` within that tree unless
otherwise specified. Line numbers are valid at the commit above; they drift
when PyLith is updated.

## 1. Fault-solver architecture

### 1.1 Weak form key encoding

File: `feassemble/IntegratorInterface.cc:419-428`
Function: `IntegratorInterface::getWeakFormPart(part, face, patch) const`

```cpp
return _labelValue*(max_parts*max_face_enums*max_patches)
     + part*num_face_enums*max_patches
     + face*max_patches
     + patch;
```

Constants (same file, lines 144-156):

- `max_parts      = 10`
- `max_face_enums = 10`
- `num_face_enums = 3`   (NEGATIVE_FACE=0, POSITIVE_FACE=1, FAULT_FACE=2)
- `max_patches    = 10`

Parameters:

- `_labelValue`: the fault mesh label value (fault identifier).
- `part`: equation part (RHS=0, LHS=1, LHS_WEIGHTED=2, LHS_LUMPED_INV=3).
- `face`: face enum (NEGATIVE_FACE=0, POSITIVE_FACE=1, FAULT_FACE=2).
- `patch`: integration patch value (0 to `max_patches - 1`).

### 1.2 Auxiliary vector attachment -- fault integrator

File: `feassemble/IntegratorInterface.cc:450-466`

```cpp
for (keysmap_t::const_iterator iter = keysmap.begin(); iter != keysmap.end(); ++iter) {
    const PetscInt patchValue = iter->second.cohesive.getValue();
    const size_t numParts = 3;
    const EquationPart equationParts[numParts] = {RHS, LHS, LHS_WEIGHTED};
    for (size_t i = 0; i < numParts; ++i) {
        PetscFormKey key = iter->second.cohesive.getPetscKey(solution, equationParts[i]);
        PetscInt part = getWeakFormPart(key.part, IntegratorInterface::FAULT_FACE, patchValue);
        err = DMSetAuxiliaryVec(dmSoln, key.label, key.value, part,
                                auxiliaryField->getLocalVector());
    }
}
```

- Label name: fault patch label (from `cohesive.getName()`).
- Values: `patchValue` (patch label value).
- Encoded key: `getWeakFormPart(equationPart, FAULT_FACE, patchValue)`.
- Equation parts attached: RHS, LHS, LHS_WEIGHTED (not LHS_LUMPED_INV for the
  fault integrator).

### 1.3 Auxiliary vector attachment -- material (domain) integrator

File: `feassemble/IntegratorDomain.cc:264-270`

```cpp
err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, LHS,            _auxiliaryField->getLocalVector());
err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, RHS,            _auxiliaryField->getLocalVector());
err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, LHS_LUMPED_INV, _auxiliaryField->getLocalVector());
```

For materials adjacent to a fault, the domain integrator also attaches
auxiliary vectors on the cohesive-side encoded keys (same file, 354-360):

```cpp
for (size_t iFace = 0; iFace < faceCount; ++iFace) {
    for (PetscInt iPart = 0; iPart < numParts; ++iPart) {
        const PetscInt part = integrator->getWeakFormPart(parts[iPart], faultFaces[iFace], patchValue);
        err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, part,
                                _auxiliaryField->getLocalVector());
    }
}
```

- Label name: material label (`dmLabel`).
- Value: `_labelValue` (material label value).
- Simple (non-fault) keys: plain `LHS` / `RHS` / `LHS_LUMPED_INV`.
- Fault-adjacent keys: `getWeakFormPart(equationPart, faultFace, patchValue)`
  for each `faultFace in {NEGATIVE_FACE, POSITIVE_FACE}`.

### 1.4 Integration patches

File: `feassemble/InterfacePatches.cc:86-153`
Function: `InterfacePatches::createMaterialPairs()`

Patches label name convention (line 92):

```cpp
const std::string& patchLabelName =
    fault->getSurfaceLabelName() + std::string("-integration-patches");
```

So for fault surface label `"fault"`, the patches label is
`"fault-integration-patches"`.

Patch assignment walks cohesive cells and groups them by
`(neg_material_value, pos_material_value)` pair (lines 109-123):

```cpp
for (PylithInt iCohesive = 0; iCohesive < numCohesiveCells; ++iCohesive) {
    const PetscInt cohesiveCell = cohesiveCells[iCohesive];
    PetscInt adjacentCellNegative = -1;
    PetscInt adjacentCellPositive = -1;
    TopologyOps::getAdjacentCells(&adjacentCellNegative, &adjacentCellPositive,
                                  dmSoln, cohesiveCell);
    std::pair<int,int> matPair;
    err = DMGetLabelValue(dmSoln, cellsLabelName, adjacentCellNegative, &matPair.first);
    err = DMGetLabelValue(dmSoln, cellsLabelName, adjacentCellPositive, &matPair.second);
    if (0 == integrationPatches.count(matPair)) {
        integrationPatches[matPair] = ++patchLabelValue;
    }
}
```

Single-material problems produce one patch; N distinct material-pair
boundaries produce N patches.

### 1.5 Kernel registration -- residual

File: `faults/FaultCohesiveKin.cc:349-416` (FSRM inline tree:
`/home/dockimble/pylith-dev/src/pylith-inline/libsrc/pylith/faults/...` same
lines in official tree).

QUASISTATIC residual registration (lines 359-380):

```cpp
const PetscBdPointFunc f0u_neg = pylith::fekernels::FaultCohesiveKin::f0u_neg;
const PetscBdPointFunc f0u_pos = pylith::fekernels::FaultCohesiveKin::f0u_pos;
const PetscBdPointFunc f0l     = pylith::fekernels::FaultCohesiveKin::f0l_slip;

kernels[0] = ResidualKernels("displacement",
                             integrator_t::LHS, integrator_t::NEGATIVE_FACE,
                             f0u_neg, NULL);
kernels[1] = ResidualKernels("displacement",
                             integrator_t::LHS, integrator_t::POSITIVE_FACE,
                             f0u_pos, NULL);
kernels[2] = ResidualKernels("lagrange_multiplier_fault",
                             integrator_t::LHS, integrator_t::FAULT_FACE,
                             f0l, NULL);
```

Kernel semantics (from `fekernels/FaultCohesiveKin.hh:80-160`):

- `f0u_neg`: elasticity residual contribution on negative face = `+lambda`
- `f0u_pos`: elasticity residual contribution on positive face = `-lambda`
- `f0l_slip`: constraint residual = `d - (u_pos - u_neg)`
  (prescribed slip minus jump).

### 1.6 Kernel registration -- Jacobian

File: `faults/FaultCohesiveKin.cc:422-481`

QUASISTATIC Jacobian blocks registered (lines 432-455):

```cpp
kernels[0] = JacobianKernels("displacement", "lagrange_multiplier_fault",
                             integrator_t::LHS, integrator_t::NEGATIVE_FACE,
                             Jf0ul_neg, NULL, NULL, NULL);
kernels[1] = JacobianKernels("displacement", "lagrange_multiplier_fault",
                             integrator_t::LHS, integrator_t::POSITIVE_FACE,
                             Jf0ul_pos, NULL, NULL, NULL);
kernels[2] = JacobianKernels("lagrange_multiplier_fault", "displacement",
                             integrator_t::LHS, integrator_t::FAULT_FACE,
                             Jf0lu, NULL, NULL, NULL);
```

Only three blocks are active: `Jf0ul_neg`, `Jf0ul_pos`, `Jf0lu`. All Jf1/Jf2/Jf3
slots are `NULL`. The Lagrange-Lagrange block (`Jf0ll` / `Jf1ll` / ...) is
**not registered at all**: PyLith's saddle-point matrix has an identically
zero Lagrange-Lagrange block and relies entirely on fieldsplit Schur with
`selfp` approximation to handle the saddle-point.

FSRM deviates here: `addCohesivePenaltyToJacobian`
(`src/core/Simulator.cpp:4677-4720`) stamps an epsilon-regularized
Lagrange-Lagrange diagonal so the direct-LU default path is non-singular.
This is the root source of scaling mismatch compared to PyLith (see
`docs/SOLVER_STATE.md` known bottleneck 3).

### 1.7 Residual hybrid-compute call

File: `feassemble/IntegratorInterface.cc:640-661`

```cpp
PetscFormKey weakFormKeys[3];
weakFormKeys[0] = iter->second.negative.getPetscKey(*solution, equationPart);
weakFormKeys[0].part = integrator->getWeakFormPart(
    equationPart, IntegratorInterface::NEGATIVE_FACE, patchValue);
weakFormKeys[1] = iter->second.positive.getPetscKey(*solution, equationPart);
weakFormKeys[1].part = integrator->getWeakFormPart(
    equationPart, IntegratorInterface::POSITIVE_FACE, patchValue);
weakFormKeys[2] = iter->second.cohesive.getPetscKey(*solution, equationPart);
weakFormKeys[2].part = integrator->getWeakFormPart(
    equationPart, IntegratorInterface::FAULT_FACE, patchValue);

err = DMPlexComputeResidual_Hybrid_Internal(dmSoln, weakFormKeys, patchCellsIS,
    t, solution->getLocalVector(), solutionDotVec, t,
    residual->getLocalVector(), NULL);
```

Each of the three keys carries the same equation part and patch but a
different face enum, so PETSc dispatches NEGATIVE_FACE, POSITIVE_FACE, and
FAULT_FACE kernels independently on the hybrid cells.

### 1.8 Rim-pin / buried edges constraint

File: `faults/FaultCohesive.cc:313-391`
Function: `FaultCohesive::createConstraints()`

Label name convention (line 338):

```cpp
buriedLabelName = getBuriedEdgesLabelName() + "_cohesive";
```

If the user supplies a `"buried_edges"` label, PyLith builds
`"buried_edges_cohesive"` for the cohesive supports and installs:

```cpp
pylith::feassemble::ConstraintSimple *constraint =
    new pylith::feassemble::ConstraintSimple(this);
constraint->setLabelName(buriedLabelName.c_str());
constraint->setConstrainedDOF(&constrainedDOF[0], constrainedDOF.size());
constraint->setSubfieldName(lagrangeName);   // "lagrange_multiplier_fault"
constraint->setUserFn(_zero);                 // lambda = 0 on buried edges
```

All components of `lambda` are pinned to zero on cohesive supports of buried
edges. FSRM's `getOrCreateFaultBoundaryLabel` + geometric rim-pin fallback
(Session 12, refined in Session 25) is the functional equivalent when no
user `buried_edges` label is supplied.

### 1.9 Null-space construction

File: `problems/Problem.cc:627-654`
Function: `Problem::createNullSpace(solution, subfieldName)`

```cpp
const int spaceDim = solution->getSpaceDim();
const PetscInt m = (spaceDim * (spaceDim + 1)) / 2;
const PetscDM dmSoln = solution->getDM();
const pylith::topology::Field::SubfieldInfo info =
    solution->getSubfieldInfo(subfieldName);
MatNullSpace nullSpace = NULL;
err = DMPlexCreateRigidBody(dmSoln, info.index, &nullSpace); PYLITH_CHECK_ERROR(err);

PetscObject field = NULL;
err = DMGetField(dmSoln, info.index, NULL, &field); PYLITH_CHECK_ERROR(err);
err = PetscObjectCompose(field, "nearnullspace", (PetscObject) nullSpace); PYLITH_CHECK_ERROR(err);
err = MatNullSpaceDestroy(&nullSpace); PYLITH_CHECK_ERROR(err);
```

PyLith builds the rigid-body null-space on the full DM (not a sub-DM),
composes it on the displacement field object, and releases its local
reference. The vectors produced are sized to the full solution section's
free storage -- this is the Session 30 "full-DM" sizing path.

FSRM carries an additional **sub-DM sized** null-space (Session 30,
`rigidBodyNullSpace_`) because FSRM uses `pc_use_amat = true` with
`pc_fieldsplit_schur_factorization_type lower` + `selfp`, and the
fieldsplit inner PCGAMG on the displacement sub-block rejects the full-DM
vectors (size 297) against the sub-block matrix (size 222). PyLith's ML
backend (see 2.3 below) handles the size mismatch internally via its own
coarsening; GAMG does not. The FSRM sub-DM sizing is a GAMG-specific
adaptation, not a deviation from PyLith's design.

## 2. Solver configuration

### 2.1 Solver defaults -- elasticity (no fault)

File: `materials/Elasticity.cc:245-301`
Function: `Elasticity::getSolverDefaults(isParallel, hasFault)`

Non-fault case (lines 256-261):

```cpp
if (!isParallel) {
    options->add("-pc_type", "lu");
} else {
    options->add("-pc_type", "gamg");
}
options->add("-ts_type", "beuler");
```

Serial: direct LU. Parallel: GAMG.

### 2.2 Solver defaults -- elasticity with faults

Same file, lines 262-290, programmatic options added when `hasFault == true`:

```
-pc_type                              fieldsplit
-pc_use_amat                          (flag)
-pc_fieldsplit_type                   schur
-pc_fieldsplit_schur_factorization_type lower
-pc_fieldsplit_schur_precondition     selfp
-pc_fieldsplit_schur_scale            1.0
-fieldsplit_displacement_ksp_type             preonly
-fieldsplit_lagrange_multiplier_fault_ksp_type preonly
-fieldsplit_displacement_pc_type              lu (serial) / ml (parallel)
-fieldsplit_lagrange_multiplier_fault_pc_type lu (serial) / ml (parallel)
```

Note: the inner PC is `lu` in serial, not `ml` or `gamg`. FSRM substitutes
`gamg` for `ml` because PETSc without Trilinos has no ML backend; this is
the closest native equivalent. FSRM has not tried the serial-LU inner PC
(doing so would defeat the scalability motivation for fieldsplit in the
first place).

### 2.3 Manual fault test config (ML-based fieldsplit)

File: `tests/manual/3d/cyclicfriction/fieldsplit.cfg:1-15`

```
[pylithapp.timedependent.formulation]
split_fields = True
matrix_type = aij
use_custom_constraint_pc = True

[pylithapp.petsc]
ksp_gmres_restart = 100
fs_pc_type = fieldsplit
fs_pc_use_amat = true
fs_pc_fieldsplit_type = multiplicative
fs_fieldsplit_displacement_pc_type = ml
fs_fieldsplit_lagrange_multiplier_pc_type = jacobi
fs_fieldsplit_displacement_ksp_type = preonly
fs_fieldsplit_lagrange_multiplier_ksp_type = preonly
```

This is NOT the Schur approach -- it is `pc_fieldsplit_type = multiplicative`
(Gauss-Seidel-like block iteration) with ML on displacement and Jacobi on
Lagrange. Present in PyLith source as a tested alternative for cyclic
friction problems.

### 2.4 Tolerances

File: `utils/PetscOptions.cc:317-330`
Function: `_PetscOptions::addSolverTolerances()`

```cpp
options->add("-ksp_rtol",  "1.0e-12");
options->add("-ksp_atol",  "1.0e-12");
options->add("-ksp_error_if_not_converged");
options->add("-snes_rtol", "1.0e-12");
options->add("-snes_atol", "1.0e-9");
options->add("-snes_error_if_not_converged");
```

**Historical discrepancy**: the Session 27 prompt and FSRM's current
fieldsplit branch use `ksp_rtol = 1e-14`, `ksp_atol = 1e-7`,
`snes_rtol = 1e-14`, `snes_atol = 5e-7`. These are TIGHTER than the current
PyLith source. Either a prior PyLith commit used these values, or the
Session 27 values came from a different PyLith config (e.g., a benchmark
cfg rather than the programmatic defaults). Session 27 landed them
nonetheless; KSP is the bottleneck, not tolerance, so loosening to the
actual PyLith 1e-12 is unlikely to unblock the PrescribedSlip case.

### 2.5 Initial guess (POD)

File: `utils/PetscOptions.cc:334-342`
Function: `_PetscOptions::addInitialGuess()`

```cpp
options->add("-ksp_guess_type", "pod");
options->add("-ksp_guess_pod_size", "8");
```

POD (Proper Orthogonal Decomposition) initial guess. Requires prior solves;
inactive on the first KSP call. FSRM has not enabled this; would require
N-step POD history cache.

## 3. Architecture differences between PyLith and FSRM

Non-exhaustive summary of where FSRM deliberately deviates:

| Area | PyLith | FSRM |
|---|---|---|
| Inner PC (serial with fault) | `lu` | `gamg` (no ML backend in PETSc without Trilinos) |
| Lambda-lambda block | identically zero | eps-regularized (`1e-4 * lambda`) for direct-LU default path |
| Null-space sizing | full-DM only | full-DM (default path) + sub-DM (fieldsplit path) |
| Rim-pin label | `buried_edges` user-supplied | geometric fallback when label absent |
| POD initial guess | default | not used |
| Tolerances | `1e-12` KSP / `1e-9` SNES | `1e-14` KSP / `5e-7` SNES (fieldsplit branch only) |

## Revision history

- Session 31: document created. PyLith commit `52b3012a0`.
