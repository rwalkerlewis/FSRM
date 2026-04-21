# Lagrange Field Fix -- Status Report

**HISTORICAL / SUPERSEDED** (Session 31).

This document is the decision trail for Sessions 10-14's initial exploration
of the PETSc 3.25 region-DS limitations that motivated FSRM's current
rim-pin + manual-interior-assembly strategy. It is retained as background
context for why the code evolved the way it did. For the current solver
state, failing tests, and open investigation bottlenecks see
`docs/SOLVER_STATE.md`. For verified PyLith architecture references see
`docs/PYLITH_REFERENCE.md`. New work should update those two documents, not
this one.

Phase 1 (cohesive-cell label) and Phase 2 (field restriction attempt) are
fully landed or reverted as described below. Phase 3 (epsilon sweep) is
historical; Sessions 10-30 have since established the working regime
(`epsilon = 1e-4`) and the outer Schur KSP is the current bottleneck, not
the epsilon choice. Phase 4 (`dm_reorder_section`) is landed and still
applicable (see `CLAUDE.md` DS/BC ordering rule 4).

Options A (PETSc upgrade), C (penalty dominance), D (fieldsplit LU) are
still orthogonally viable directions but each is subsumed into the current
`docs/SOLVER_STATE.md` bottleneck list (Schur approximation, factorization
type, Lagrange block scaling). Do not restart these investigations from
scratch; measure first, then pick from the remaining candidates.

## Phase: Partial resolution via Option B (interior assembly only)

A scoped subset of `docs/PYLITH_COMPATIBILITY.md` Section B Option B has landed.
Interior (non-cohesive) Lagrange degrees of freedom are now managed entirely by
manual assembly: `addInteriorLagrangeResidual` in `src/core/Simulator.cpp:4180-4271`
zeros the residual at every non-cohesive Lagrange DOF after the FEM volume residual
runs, and the disjoint loop inside `addCohesivePenaltyToJacobian`
(`src/core/Simulator.cpp:4677-4720`) stamps a canonical penalty-scaled diagonal on
those DOFs so the LU factorization stays non-singular. The PetscDS weak volume
callback is RETAINED (the strict spec called for removal but PETSc 3.25 in this
configuration regresses `LockedFaultElastodynamic` when the Lagrange field has no
volume residual; the manual zeroing overrides the small epsilon = 1e-4 contribution
at interior DOFs). The cohesive vertex Jacobian and PetscDS BdResidual paths for
locked, prescribed slip, and slipping modes are unchanged.

Separately, the prescribed-slip Cartesian vector is now copied from
`cohesive_kernel_->prescribed_slip_` into the unified constants array in
`setupPhysics` (`src/core/Simulator.cpp:2234-2249`). The previous code allocated the
COHESIVE_CONST_PRESCRIBED_SLIP_* slots in CohesiveFaultKernel::COHESIVE_CONST_COUNT
but never wrote them, so `f0_prescribed_slip` always read a zero target jump.

## Outstanding

Options A, C, and D as described below remain rejected.

`Integration.DynamicRuptureSolve.PrescribedSlip` (the
`DynamicRuptureSolveTest.PrescribedSlipQuasiStatic` test) is still failing in this
state. The PetscDS BdResidual on the Lagrange field does not produce a non-zero
constraint contribution on the 4x4x4 quasi-static configuration even after the
prescribed-slip constants are pushed correctly, indicating the failure is not in
the constants path. The fault label registered for the Lagrange-field
`DM_BC_NATURAL` boundary contains only depth less than dim points, so
`DMPlexComputeBdResidual` has no cohesive cells to integrate over for that field.
A full Option B implementation (manual cohesive Lagrange residual that bypasses
PetscDS BdResidual entirely for prescribed slip, paired with `penalty_scale = 0`
in the cohesive penalty so the Newton step is not damped) was prototyped but
introduced separate convergence failures in PETSc 3.25 LU and was reverted. A
follow-up session should either reformulate the prescribed-slip residual against
the displacement equation directly or pursue Option A (PETSc upgrade).

## Phase: Debugging (between Phase 5 and Phase 6)

## What was attempted

### Phase 1: Cohesive cells label (DONE)
- Enabled `createCohesiveCellLabel()` using PyLith's `DMPlexGetSimplexOrBoxCells`
  pattern (from `libsrc/pylith/faults/TopologyOps.cc` line ~438)
- Original `DMPlexLabelComplete` approach from the user prompt causes
  section/closure size mismatches in `DMPlexInsertBoundaryValues` because
  it labels shared vertices between regular and cohesive cells
- The `DMPlexGetSimplexOrBoxCells` approach labels ONLY hybrid points at each
  topological dimension, avoiding shared vertices

### Phase 2: Label restriction for DMAddField (REVERTED)
- `DMAddField(dm, cohesive_label, fe_lagrange)` creates a region-specific DS
- The cohesive region DS has both displacement and lagrange fields
- The default DS has only displacement
- PetscDS volume callbacks (`PetscDSSetResidual`, `PetscDSSetJacobian`) must
  be registered on the correct DS
- **Problem:** Even after copying displacement callbacks to the cohesive region
  DS and registering Lagrange callbacks there, `DMPlexTSComputeIFunctionFEM`
  produces zero solutions. The region DS mechanism in PETSc 3.25 does not
  properly integrate volume terms on labeled regions.
- **Reverted:** Field back on all cells (nullptr label)

### Phase 3: Regularization (INVESTIGATED multiple approaches)
- **f=lambda, g=I (original):** SNES converges to zero (regularization overwhelms
  BdResidual on coarse meshes). This is the pre-existing bug.
- **f=0, g=I:** First Newton step diverges to NaN. Identity Jacobian for the
  Lagrange block (magnitude 1) produces huge Newton updates when BdResidual
  (magnitude ~1e7) dominates the right-hand side.
- **f=0, g=0:** Singular matrix, LU fails immediately (zero diagonal rows).
- **f=eps*lambda, g=eps*I (eps=1e-10):** Too small for LU factorization.
- **f=eps*lambda, g=eps*I (eps=1e-8):** Condition number ~1e18, LU fails.
- **f=eps*lambda, g=eps*I (eps=1):** Same as original, converges to zero.
- **DMSetFieldAvoidTensor for Lagrange:** Prevents DOF allocation on cohesive
  cells, breaking BdResidual evaluation.

### Phase 4: dm_reorder_section options (DONE, but caused secondary issue)
- Added `-dm_reorder_section true` and `-dm_reorder_section_type cohesive`
  to all fault test files
- This caused `DMCreateDS` to create the section eagerly, blocking
  `setupBoundaryConditions` from calling `DMAddBoundary`
- **Fixed:** Added `DMSetLocalSection(dm, nullptr)` between `DMCreateDS` and
  `setupBoundaryConditions`

## Root cause analysis

The 3 failing tests (`LockedFaultTransparency`, `LockedQuasiStatic`,
`PrescribedSlipQuasiStatic`) produce zero solutions because:

1. The Lagrange multiplier field exists on ALL cells in the mesh
2. The volume regularization `f = lambda` on interior (non-cohesive) cells
   drives lambda toward zero
3. On cohesive cells, BdResidual provides the constraint equation, but the
   volume `f = lambda` competes with it
4. On coarse meshes (4x4x4 cells), interior cells vastly outnumber cohesive
   cells (~384 vs ~32), so the regularization overwhelms the constraint

The fix requires one of:
- **True label restriction:** Lagrange DOFs only on cohesive cell closure
  (eliminates interior Lagrange DOFs entirely)
- **Or:** Making the volume callback zero/negligible on ALL cells while keeping
  the system non-singular for LU

Both approaches hit PETSc 3.25 limitations:
- Label restriction: region DS does not support volume assembly
- Zero residual: creates singular or ill-conditioned Jacobian

## PyLith files read

1. `libsrc/pylith/topology/MeshOps.cc` line ~247: `DMPlexLabelComplete` for submesh
2. `libsrc/pylith/faults/TopologyOps.cc` line ~438: `getInterfacesLabel` using
   `DMPlexGetSimplexOrBoxCells` per height stratum
3. `libsrc/pylith/topology/Field.cc` line ~556-561: `DMSetField` with interfaces label
   for `isFaultOnly` subfields
4. `libsrc/pylith/materials/Elasticity.cc` line ~314: `dm_reorder_section` options
5. `share/settings/solver_elasticity_fault.cfg`: production PETSc options
6. `libsrc/pylith/problems/Problem.cc` line ~651: region DS setup for fault fields

## What would fix this

PyLith's approach works because PyLith's version of PETSc (recent main branch)
has fixes for region DS volume assembly that PETSc 3.25.0 (release) lacks. The
`DMSetField(dm, index, interfacesLabel, fe)` call in PyLith works end-to-end
because their PETSc properly handles:
1. Region DS creation with labeled fields
2. Volume residual assembly on labeled regions via `DMPlexTSComputeIFunctionFEM`
3. Boundary value insertion with mixed-region sections

### Option A: Upgrade to PETSc main branch
Building PETSc from the main branch (post-3.25.0) would gain the region DS fixes
that PyLith relies on. This is the most reliable path.

### Option B: Manual assembly for ALL Lagrange equations
Instead of relying on PetscDS for the Lagrange field, skip PetscDS registration
entirely and handle all Lagrange assembly manually (like `addCohesivePenaltyToJacobian`
does for the Jacobian). This would require a manual residual assembly function
that evaluates BdResidual-like equations on cohesive cells and zero on interior cells.

### Option C: Augmented Lagrangian with penalty dominance
Use a penalty approach where the cohesive penalty in `addCohesivePenaltyToJacobian`
is made large enough to dominate the volume regularization. This would require
scaling the penalty by 1/epsilon relative to the regularization.

## Additional investigation (epsilon sweep)

After the initial investigation, a systematic sweep of epsilon values for
the weak regularization `f = epsilon * lambda, g = epsilon * I` was performed:

| epsilon | LU factorization | Newton convergence | Result |
|---------|------------------|--------------------|--------|
| 1.0     | OK               | Converges to zero  | lambda overwhelms BdResidual |
| 1e-2    | OK               | NaN at iteration 1 | Jacobian mismatch on cohesive cells |
| 1e-3    | OK               | NaN at iteration 4 | Same, slower to diverge |
| 1e-6    | Fails (singular) | N/A                | Condition number ~1e16 |
| 1e-8    | Fails (singular) | N/A                | Condition number ~1e18 |
| 1e-10   | Fails (singular) | N/A                | Below machine epsilon ratio |
| 0       | Fails (singular) | N/A                | Zero diagonal rows |

The NaN divergence at epsilon < 1 is caused by the PetscDS Jacobian
(epsilon * I from volume integral) being evaluated on cohesive cells where it
conflicts with the manual Jacobian from addCohesivePenaltyToJacobian. The
manual Jacobian adds penalty entries O(E/h) ~ O(4e10), while the PetscDS adds
epsilon * I. The combined Jacobian on cohesive cells is inconsistent with the
combined residual (epsilon * lambda from PetscDS + BdResidual from PetscDS),
causing Newton to diverge.

Also attempted `subtractLagrangeRegularizationOnCohesive()` -- a manual
correction that subtracts the f=lambda contribution from cohesive vertices
after assembly. This did not work because cohesive vertices are also in the
closure of interior cells, and the subtraction only removes contributions from
cohesive cell quadrature points, not from the surrounding interior cells.

## What to try next

1. **Manual Lagrange assembly (Option B):** Most feasible without changing PETSc.
   Bypass PetscDS entirely for the Lagrange field: set no PetscDS callbacks,
   and add both residual AND Jacobian for the Lagrange field manually in
   FormFunction/FormJacobian. The manual residual would evaluate the BdResidual
   kernel on cohesive faces. The manual Jacobian already exists. Interior
   Lagrange DOFs would get a scaled identity diagonal (~E) in the manual
   Jacobian and zero residual, keeping the system non-singular.

2. **PETSc upgrade (Option A):** Build PETSc from main branch (post-3.25.0)
   and use DMAddField with label for true field restriction.

3. **Separate LU factors (Option D):** Use PETSc FieldSplit with separate LU
   on displacement and Jacobi on Lagrange. This avoids the conditioning issue
   entirely by not mixing the two field scales in one LU factorization.
