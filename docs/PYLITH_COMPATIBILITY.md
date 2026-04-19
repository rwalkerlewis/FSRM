# PyLith Compatibility Features

This document describes FSRM features that align with PyLith workflows and conventions,
enabling users familiar with PyLith to set up equivalent simulations in FSRM.

## Overview

The short summary table below is obsolete, superseded by Section A. It is retained for
continuity but must not be cited as authoritative. Section A is the canonical inventory.

| Feature | FSRM Status | PyLith Equivalent |
|---------|-------------|-------------------|
| Per-face Neumann traction BCs | Working (obsolete, superseded by section A) | Neumann BC on faces |
| Prescribed kinematic fault | Working (obsolete, superseded by section A) | FaultCohesiveKin |
| Time-dependent slip ramp | Working (obsolete, superseded by section A) | TimeHistorySlipFn |
| Locked fault (transparent) | Working (obsolete, superseded by section A) | FaultCohesiveKin (zero slip) |
| Roller (zero normal) BCs | Working (obsolete, superseded by section A) | Dirichlet BC (single component) |
| Cohesive cell mesh splitting | Working (obsolete, superseded by section A) | DMPlexConstructCohesiveCells |
| Lagrange multiplier fault traction | Working (obsolete, superseded by section A) | FaultCohesiveKin constraint |

## Per-face Boundary Conditions

FSRM supports per-face boundary condition specification, similar to PyLith. Each face
of the bounding box can be independently configured.

### Configuration

```ini
[BOUNDARY_CONDITIONS]
x_min = roller
x_max = roller
y_min = roller
y_max = roller
z_min = fixed
z_max = traction
z_max_traction = 0.0, 0.0, -1.0e6
```

### Supported BC Types

| Type | Description | PyLith Equivalent |
|------|-------------|-------------------|
| `fixed` | All displacement components = 0 | Dirichlet (all components) |
| `roller` | Normal component = 0 | Dirichlet (single component) |
| `free` | Traction-free (natural BC) | Default (no BC) |
| `traction` | Applied traction vector | Neumann BC |
| `compression` | Normal compressive stress | Neumann BC (normal) |
| `absorbing` | Clayton-Engquist absorbing | AbsorbingDampers |

### Legacy Format

The per-face system also supports the legacy shorthand format:

```ini
[BOUNDARY_CONDITIONS]
bottom = fixed
sides = roller
top = traction
z_max_traction = 0.0, 0.0, -1.0e6
```

## Traction Boundary Conditions

Neumann (traction) BCs apply a specified stress vector to a boundary face. The traction
vector is specified in Cartesian coordinates (Pa).

### Implementation

FSRM uses manual FEM assembly (DMPlex point closure) rather than PetscDS boundary
residual callbacks. This is the same proven pattern used for explosion sources, fault
pressure, and cohesive constraints.

For each traction face:
1. Get boundary label IS (e.g., `boundary_z_max`)
2. Filter for depth `dim-1` faces with `support_size == 1` (boundary faces only)
3. Compute face area via `DMPlexComputeCellGeometryFVM`
4. Get face closure to find associated nodes
5. Distribute force equally: `F_node = traction * face_area / n_nodes`

### Analytical Verification

The traction BC is verified against the analytical uniaxial stress solution:

- Column height L = 100 m, E = 10 GPa, sigma = -1 MPa
- Expected: u_z(L) = sigma * L / E = -0.01 m
- Test: `Integration.TractionBC` (max displacement > 1e-4 m, < 1 m)

## Time-dependent Prescribed Slip

FSRM supports kinematic fault slip with a linear time ramp, similar to PyLith
`TimeHistorySlipFn`.

### Configuration

```ini
[FAULT]
mode = prescribed_slip
strike = 90.0
dip = 90.0
center_x = 0.5
center_y = 0.5
center_z = 0.5
length = 2.0
width = 2.0
slip_strike = 0.0
slip_dip = 0.0
slip_opening = 0.001
friction_coefficient = 0.6
slip_onset_time = 0.1
slip_rise_time = 0.5
```

### Slip Time Function

```
slip(t) = 0                                          for t < onset_time
slip(t) = slip_max * (t - onset_time) / rise_time    for onset_time <= t < onset_time + rise_time
slip(t) = slip_max                                   for t >= onset_time + rise_time
```

When `rise_time = 0`, the slip is applied instantaneously at `t >= onset_time`.

### Test Coverage

| Test | What It Verifies |
|------|-----------------|
| `Integration.TimeDependentSlip` | Slip ramp completes, nonzero solution |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | Instant slip, displacement jump |
| `Physics.LockedFaultTransparency` | Zero slip gives slip < 5e-4 |

## Cohesive Cell Fault Representation

FSRM uses the same cohesive cell approach as PyLith:

1. `FaultMeshManager::markFaultFaces()` identifies interior faces on the fault plane
2. `DMPlexConstructCohesiveCells()` splits the mesh along the fault
3. Lagrange multiplier field enforces the fault constraint
4. PetscDS BdResidual callbacks on cohesive cells enforce the kinematic condition (PETSc 3.25+). Jacobian assembled manually via addCohesivePenaltyToJacobian.

### Fault Setup Pipeline

```
setupDM()            -- create hex/tet mesh (simplex if faults enabled)
setupFaultNetwork()  -- mark faces, split mesh, insert cohesive cells
labelBoundaries()    -- label boundary faces (after split)
setupFields()        -- add displacement + Lagrange multiplier fields
```

This ordering must not be changed. The fault mesh split must happen before boundary
labeling because it changes the mesh topology.

## Examples

| # | Directory | Config | Feature |
|---|-----------|--------|---------|
| 07 | `examples/07_traction_bc/` | `traction_bc.config` | Per-face traction BC |
| 08 | `examples/08_time_dependent_slip/` | `time_dependent_slip.config` | Slip ramp |

## Differences from PyLith

| Aspect | FSRM | PyLith |
|--------|------|--------|
| Configuration | INI-style `.config` files | `.cfg` files with Pyre framework |
| Mesh | Built-in hex + Gmsh import | Gmsh, CUBIT, LaGriT |
| Time function | Linear ramp only (obsolete, superseded by section A) | Step, constant rate, Liu cosine, time history |
| Fault friction | Coulomb (manual assembly) (obsolete, superseded by section A: FSRM also implements linear slip-weakening in src/physics/CohesiveFaultKernel.cpp:256-273) | Rate-and-state, slip weakening |
| Output | HDF5, VTK, SAC | HDF5/Xdmf |
| Parallelism | MPI + optional GPU (PETSc CUDA) | MPI only |

## A. Current State

This section replaces the obsolete Overview table. Every row cites a specific test or a
file:line reference that proves or locates the feature. Classes whose declarations exist
but are not wired into the PetscDS pipeline are marked "class exists, not wired" with a
statement of the missing wiring.

### A.1 Mesh and Field Infrastructure

| Feature | Status | Proof or Location |
|---------|--------|-------------------|
| Cohesive cell insertion via `DMPlexConstructCohesiveCells` | Working | `src/numerics/FaultMeshManager.cpp:147-229` (splitMeshAlongFault). Follows PyLith workflow documented at `/home/dockimble/Projects/pylith/libsrc/pylith/topology/MeshOps.cc` line 247 and `/home/dockimble/Projects/pylith/libsrc/pylith/faults/FaultOps.hh` (referenced). |
| Planar fault labeling | Working | `src/numerics/FaultMeshManager.cpp:52-145` (createPlanarFaultLabel). Skips boundary faces with support_size < 2. |
| Cohesive topology extraction (neg/pos vertex pairing) | Working | `src/numerics/FaultMeshManager.cpp:231-358` (extractCohesiveTopology). Pairs neg_verts[i] and pos_verts[i] by PETSc cone ordering. |
| Fault geometry (normal, strike, dip, area) | Working | `src/numerics/FaultMeshManager.cpp:360-591` (computeFaultGeometry). |
| Lagrange multiplier field on cohesive cells | Working (with caveats) | `src/core/Simulator.cpp:1680-1692` adds `fe_lagrange` via `DMAddField(dm, nullptr, ...)` on ALL cells (not label-restricted). See Section B. |
| Weak volume regularization on Lagrange field | Working | `src/core/Simulator.cpp:1722-1776` (f0_weak_lagrange, g0_weak_lagrange) with epsilon=1e-4. |
| Manual cohesive penalty assembly in Jacobian | Working | `src/core/Simulator.cpp:4182-4591` (addCohesivePenaltyToJacobian). Provides penalty-scaled diagonal at cohesive vertices (locked/prescribed modes) and semi-smooth Newton Jacobian for friction (slipping mode). |

### A.2 Fault Constraint Modes

| Mode | Status | Proof or Location |
|------|--------|-------------------|
| Locked fault (u+ = u- via Lagrange) | Working | `Physics.LockedFaultTransparency` (`tests/physics_validation/test_locked_fault_transparency.cpp:328-345`) and `Integration.DynamicRuptureSolve.LockedFaultQuasiStatic` (`tests/integration/test_dynamic_rupture_solve.cpp:374-392`). Max slip < 5e-4 m tolerance. Residual is `f0_lagrange_constraint` (`src/physics/CohesiveFaultKernel.cpp:183-215`). |
| Locked fault with elastodynamics (TSBEULER) | Working | `Integration.DynamicRuptureSolve.LockedFaultElastodynamic` (`tests/integration/test_dynamic_rupture_solve.cpp:395-411`). |
| Prescribed slip (u+ - u- = delta) | Partially working | Callback verified in isolation by `Integration.PrescribedSlip.CallbackEnforcesPrescribedJump` (`tests/integration/test_prescribed_slip.cpp:33-72`) and `.CallbackNonzeroResidualForMismatch` (lines 74-111). Full TSSolve integration `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic` (lines 418-435) is the known failure documented in CLAUDE.md. Callback at `src/physics/CohesiveFaultKernel.cpp:629-657`. |
| Time-dependent slip ramp (linear only) | Working | `Integration.TimeDependentSlip.SlipRampCompletesSuccessfully` (`tests/integration/test_time_dependent_slip.cpp:82-...`). Slip amplitude is scaled in `src/core/Simulator.cpp:setupFaultNetwork` via fixed Cartesian vector; no time-history callback wired. |
| Slipping fault with Coulomb friction (semi-smooth Newton) | Working | `Integration.SlippingFaultSolve` (`tests/integration/test_slipping_fault_solve.cpp`). Jacobian at `src/core/Simulator.cpp:4393-4548`. Residual branch at `src/physics/CohesiveFaultKernel.cpp:217-309`. |
| Slipping fault with linear slip-weakening | Working | `Physics.SCEC.TPV5` (`tests/physics_validation/test_scec_tpv5.cpp:113-304`). Friction model read from constants slot 70 (`COHESIVE_CONST_FRICTION_MODEL`). Callback branches at `src/physics/CohesiveFaultKernel.cpp:257-273` (residual), `src/physics/CohesiveFaultKernel.cpp:411-427` (J_lu), `src/physics/CohesiveFaultKernel.cpp:543-557` (J_ll). Manual Jacobian has matching branch at `src/core/Simulator.cpp:4443-4448` and the slip-weakening softening derivative at `src/core/Simulator.cpp:4515-4519`. |
| Tensile failure (opening when lambda_n > T_s) | Working (callback-level) | `src/physics/CohesiveFaultKernel.cpp:191-209` (tensile-failure branch in f0_lagrange_constraint). No dedicated integration test through TSSolve. |
| Kinematic Step / Brune / Liu / TimeHistory slip time functions | Class exists, not wired | Declared at `include/domain/geomechanics/PyLithFault.hpp:147-321` (`SlipTimeFnStep`, `SlipTimeFnRamp`, `SlipTimeFnBrune`, `SlipTimeFnLiu`, `SlipTimeFnConstantRate`, `SlipTimeFnTimeHistory`). `SlipTimeFnRamp` is logically referenced by the config keys `slip_onset_time` and `slip_rise_time`, but no call to `createSlipTimeFn` is emitted from `src/core/Simulator.cpp`; the prescribed slip vector is instead pushed straight into constants slots 28 through 30. The `SlipTimeFn` abstraction is not consulted by `f0_prescribed_slip`. |

### A.3 Friction Models

| Model | Status | Proof or Location |
|-------|--------|-------------------|
| Constant Coulomb friction | Working in Lagrange constraint kernel | `src/physics/CohesiveFaultKernel.cpp:270-273` and semi-smooth Newton Jacobian at `src/core/Simulator.cpp:4446-4448`. Integration-tested by `Integration.SlippingFaultSolve`. |
| Linear slip-weakening (mu_s, mu_d, Dc) | Working in Lagrange constraint kernel | `src/physics/CohesiveFaultKernel.cpp:260-269`, `src/physics/CohesiveFaultKernel.cpp:415-422`, `src/physics/CohesiveFaultKernel.cpp:546-553`. Integration-tested by `Physics.SCEC.TPV5`. The Coulomb-only row in "Differences from PyLith" is obsolete and already annotated above. |
| Static friction class | Class exists, partially wired | `include/domain/geomechanics/PyLithFault.hpp:479-497` (`StaticFriction`). Used only by the standalone `FaultCohesiveDyn` object path, not by `CohesiveFaultKernel`. |
| SlipWeakeningFriction C++ class | Class exists, used for standalone path only | `include/domain/geomechanics/PyLithFault.hpp:502-529`. Instantiated by `src/core/Simulator.cpp:3458-3462` and passed to `cohesive_fault_->setFrictionModel(...)` but the PetscDS kernels read their parameters directly from constants slots 70 through 73, not from this object. The object is effectively cosmetic for the PetscDS path. |
| RateStateFrictionAging C++ class | Class exists, not wired | `include/domain/geomechanics/PyLithFault.hpp:534-565`. No `CohesiveFaultKernel` constant slots are allocated for a, b, Dc, f0, V0. No auxiliary field for theta. No callback branch reads state_variable. Unit tests exercise the object only in isolation at `tests/unit/domain/test_fault_cohesive_dyn.cpp:47-75` and `tests/unit/domain/test_fault_mechanics.cpp:42-57`. |
| RateStateFrictionSlip C++ class | Class exists, not wired | `include/domain/geomechanics/PyLithFault.hpp:570-591`. Same missing wiring as the aging variant. |
| FaultCohesiveTHM (poroelastic / thermal fault coupling) | Class exists, not wired | `include/domain/geomechanics/PyLithFault.hpp:1432-1541`. The header declares `addPoroelasticResidual`, `addThermalResidual`, `addFlowResidual`, but no instance is constructed by `src/core/Simulator.cpp:setupFaultNetwork` and no PetscDS callback is registered for a fault-pressure Lagrange multiplier. See Section E. |

### A.4 PyLith Workflow Alignment (Kept)

These remain correctly described in the original text (sections above) and the row status
is unchanged under the new inventory:

- Per-face Neumann traction BCs: Working, `Integration.TractionBC` (`tests/integration/test_traction_bc.cpp`).
- Cohesive cell mesh splitting via `DMPlexConstructCohesiveCells`: Working (see A.1).
- Lagrange multiplier fault traction: Working for locked and Coulomb modes (see A.2).

## B. Lagrange Architecture Decision

The options investigated in `docs/LAGRANGE_FIX_STATUS.md` were:

- Option A. Upgrade PETSc to the main branch (post 3.25.0) to gain the region DS fixes
  that PyLith relies on. Referenced in `docs/LAGRANGE_FIX_STATUS.md` lines 92-95.
- Option B. Manual assembly for all Lagrange equations, bypassing PetscDS entirely for
  the Lagrange field. Referenced in `docs/LAGRANGE_FIX_STATUS.md` lines 97-101 and
  lines 139-145.
- Option C. Augmented Lagrangian with penalty dominance, scaling the cohesive penalty
  above the volume regularization. Referenced in `docs/LAGRANGE_FIX_STATUS.md`
  lines 103-106.
- Option D. FieldSplit with separate LU on displacement and Jacobi on Lagrange.
  Referenced in `docs/LAGRANGE_FIX_STATUS.md` lines 150-152.

### Current implementation (per CLAUDE.md and source)

The repository currently implements a hybrid of Options B and C. The concrete choices are:

1. The Lagrange field is added on the whole mesh via
   `DMAddField(dm, nullptr, fe_lagrange)` at `src/core/Simulator.cpp:1687`. No label
   restriction is used.
2. A weak volume regularization `f = epsilon * lambda`, `g = epsilon * I` with
   `epsilon = 1e-4` is registered on the Lagrange field via `PetscDSSetResidual` and
   `PetscDSSetJacobian` at `src/core/Simulator.cpp:2528-2532`, using the static
   `f0_weak_lagrange` and `g0_weak_lagrange` defined at
   `src/core/Simulator.cpp:1739-1776`.
3. `PetscDSSetBdResidual` registers the cohesive constraint callback on cohesive cells
   (locked, slipping, or prescribed) at `src/core/Simulator.cpp:2540-2548`. This works
   in PETSc 3.25.0, verified by `Physics.CohesiveBdResidual`
   (`tests/physics_validation/test_cohesive_bd_residual.cpp:103-243`).
4. `PetscDSSetBdJacobian` is not functional in PETSc 3.25.0. The cohesive Jacobian
   contribution is therefore assembled manually in `addCohesivePenaltyToJacobian`
   at `src/core/Simulator.cpp:4182-4591`, which adds a penalty-scaled diagonal
   `penalty * coeff` with `penalty = penalty_scale * 10.0 * youngs_modulus / h_char_jac`
   at cohesive vertices (`src/core/Simulator.cpp:4341-4342`).

CLAUDE.md (FSRM root) documents this as the current state under "Lagrange Multiplier
Field and Cohesive Cells" and "PDE Assembly".

### Why `PrescribedSlipQuasiStatic` fails under this choice

The failure is documented in CLAUDE.md under "Test Suite -> Known failure". The
PetscDS constraint residual for prescribed slip is
`f_lambda = (u+ - u-) - delta_prescribed` at `src/physics/CohesiveFaultKernel.cpp:651-656`.
Locked mode reaches the correct solution because the residual drives the jump to zero
against a small penalty-damped Lagrange update. For prescribed slip, the Lagrange
multiplier must take a nonzero traction that exactly supports the imposed jump. The
penalty-scaled diagonal `penalty * coeff` at
`src/core/Simulator.cpp:4573` and the matching displacement-displacement penalty block
at `src/core/Simulator.cpp:4574-4577` damp each Newton update to the Lagrange
multiplier. On a 4x4x4 coarse mesh with one quasi-static step, the damping is large
enough that the Lagrange update in one Newton iteration does not deliver the target
traction, so the jump response stays near zero. This is an augmented Lagrangian scaling
problem, not a physics problem. CLAUDE.md states this explicitly: "the penalty-scaled
Lagrange diagonal in addCohesivePenaltyToJacobian slows lambda convergence for
prescribed slip mode".

### Recommendation before adding rate-state

Commit to Option B (full manual Lagrange assembly) before adding rate-and-state. The
justification:

1. Option A requires a PETSc fork and breaks the Docker build constraint documented
   in CLAUDE.md rule 1 ("Build and test in Docker. Always."). The upgrade is also
   orthogonal to the planned work.
2. Option C alone is what we have, and it fails `PrescribedSlipQuasiStatic`. Tuning
   the penalty is forbidden by the spirit of rule 15 and by the project rule against
   ad-hoc tuning, and per `docs/LAGRANGE_FIX_STATUS.md` (epsilon sweep at lines 108-129)
   no epsilon satisfies both LU stability and Newton convergence simultaneously for
   a regularization-dominated system.
3. Option D (FieldSplit) does not address the fundamental residual assembly problem;
   the regularization still competes with BdResidual regardless of the linear solver
   partitioning.
4. Option B gives a clean separation: the PetscDS volume callback on the Lagrange field
   is removed entirely, interior Lagrange DOFs get a scaled identity in the manual
   Jacobian (keeping LU non-singular), and the BdResidual kernel is the sole driver
   of the Lagrange residual. This is the only option that makes the prescribed-slip
   constraint directly achievable by one Newton step.

Rate-and-state requires a well-behaved Lagrange system because the friction law reads
lambda_n to form the effective normal stress. The current augmented-Lagrangian damping
would corrupt that read-back and produce unphysical state evolution. Option B is the
load-bearing precondition for Section C.

## C. Rate-and-State Wiring Plan

Only actionable after Section B.

### C.1 Constant slot allocation

The existing constant slots that must not be collided with are documented in CLAUDE.md
(Unified Constants Array) and confirmed in the source:

- Slots 0-24 are fluid / elasticity / poroelasticity properties
  (`src/core/Simulator.cpp:2120-2163`).
- Slots 25-30 are cohesive kernel mode, friction coefficient, tensile strength, and
  prescribed slip components (`include/physics/CohesiveFaultKernel.hpp:47-53`).
- Slots 32-35 are elastoplasticity (CLAUDE.md, Unified Constants).
- Slots 36-53 are traction BC values (`src/core/Simulator.cpp:2232-2242`).
- Slots 54-69 are viscoelastic (CLAUDE.md, Unified Constants).
- Slots 70-73 are slip-weakening friction model parameters
  (`include/physics/CohesiveFaultKernel.hpp:56-59`).
- Slots 74-78 are thermal (CLAUDE.md, Unified Constants).
- MAX_UNIFIED_CONSTANTS is 80 (`src/core/Simulator.cpp:2106`).

Propose allocating slots 79 for rate-state model selector and extending
`MAX_UNIFIED_CONSTANTS` to at least 96 for the remaining rate-state scalar parameters:

- Slot 79: `COHESIVE_CONST_RS_EVOLUTION` (0 = aging law, 1 = slip law).
- Slot 80: `COHESIVE_CONST_RS_A` (direct effect, dimensionless).
- Slot 81: `COHESIVE_CONST_RS_B` (evolution effect, dimensionless).
- Slot 82: `COHESIVE_CONST_RS_DC` (critical slip distance, meters).
- Slot 83: `COHESIVE_CONST_RS_F0` (reference friction).
- Slot 84: `COHESIVE_CONST_RS_V0` (reference slip rate, m/s).
- Slot 85: `COHESIVE_CONST_RS_LIN_THRESH` (regularization slip rate).
- Slots 86-91 are reserved for thermal-pressurization coupling (Section E).
- Slots 92-95 are free scratch for nucleation perturbation.

The new block strictly follows the existing "disjoint contiguous band" pattern so no
existing callback or constant writer collides.

### C.2 Auxiliary field for theta on cohesive cells

Rate-state requires a per-vertex state variable theta that persists across time steps
and evolves nonlocally. Constants are not appropriate because theta is spatially
variable. PyLith uses an auxiliary subfield on the cohesive interface for fault state
(see `/home/dockimble/Projects/pylith/libsrc/pylith/faults/FaultCohesiveKin.cc:164-228`
for the analogous pattern used for slip, plus
`/home/dockimble/Projects/pylith/libsrc/pylith/faults/AuxiliaryFieldFactory.cc` which is
referenced for subfield registration).

Proposed FSRM approach that stays within rule 4 and rule 6:

1. Create a new auxiliary DM (or extend `auxDM_` referenced at
   `src/core/Simulator.cpp:2766`) with a scalar `theta_` subfield defined on the
   cohesive cell closure.
2. Register the aux vec via `DMSetAuxiliaryVec(dm, fault_label, 0, 0, theta_vec)`
   using the fault label already created at
   `src/numerics/FaultMeshManager.cpp:81-82`. The fault label is attached to the DM,
   survives mesh split, and defines exactly the cohesive region where theta is needed.
3. The callback reads theta via `a[aOff[theta_field]]` in the existing
   PetscDS kernel signature. This requires no new callback signatures; the aOff
   pathway is already present and currently unused at
   `src/physics/CohesiveFaultKernel.cpp:128-135` (the parameters `aOff`, `a` are
   cast to `(void)` in f0_displacement_cohesive and f0_lagrange_constraint).

### C.3 Callback signatures

No new signatures are required. The existing PetscDS boundary pointwise function
signature
(`PetscInt dim, PetscInt Nf, PetscInt NfAux, const PetscInt uOff[], const PetscInt uOff_x[], ...`)
is reused. The callback reads:

- `u[uOff[0]+d]` and `u[uOff[0]+dim+d]` for the negative and positive side displacement
  as already done at `src/physics/CohesiveFaultKernel.cpp:207-208` and :230-231.
- `u[uOff[1]+d]` for the Lagrange traction as at
  `src/physics/CohesiveFaultKernel.cpp:156` and :252-254.
- `a[aOff[theta_field]]` for the state variable. This is new.
- `u_t[uOff[0]+d]` and `u_t[uOff[0]+dim+d]` for the slip rate if a DYNAMIC_IMEX style
  formulation is used. Currently `u_t` is unused
  (`src/physics/CohesiveFaultKernel.cpp:134` marks it void).

### C.4 State-update strategy

Three candidates, with the selected choice last:

1. Aux-field time derivative (PetscDS native). Treats theta as a true solution
   variable with `du/dt = F(u, theta)` embedded in the implicit system. Advantages:
   monolithic convergence; no operator split error. Disadvantages: requires adding a
   theta field to the solution vector, not to aux, which would expand the Jacobian
   stencil and collide with the architecture decision in Section B (which reduces
   the Lagrange system to manual assembly). Rejected for this reason.
2. Fully implicit theta in a Schur-complement-style inner Newton. Advantages:
   physically accurate at long time steps. Disadvantages: the aging law derivative
   `dtheta/dt = 1 - V * theta / Dc` is stiff when V is small and the inner Newton
   cost is high. Feasible only after Option B of Section B is in place. Rejected
   for the first version because it compounds nonlinearity.
3. TSPostStep state update. At the end of each accepted time step, walk cohesive
   vertices, read V from the accepted displacement (slip rate estimated by finite
   difference from `solution` and `solution_old` already maintained at
   `src/core/Simulator.cpp:1575-1576`), and advance theta by an exponential
   integrator of the aging or slip evolution law. This is the choice.

Justification for the TSPostStep choice:

- Consistent with the existing asynchronous state-update pattern used by
  `Integration.SlippingFaultSolve` for cumulative slip tracking.
- PyLith uses a similar post-step hook for fault friction state
  (`/home/dockimble/Projects/pylith/libsrc/pylith/faults/FaultCohesiveKin.cc:234-260`
  `updateAuxiliaryField` is called per step with the accepted solution, analogous in
  role but invoked differently because PyLith treats slip as a pre-step auxiliary
  update rather than a post-step state integration).
- Keeps the PetscDS boundary callback pure with respect to theta: the callback
  reads theta as constant within a step, avoiding a nonlinear bootstrap problem on
  the first Newton iteration.
- Does not require modifying any file prohibited by rule 5 or rule 6.

## D. Kinematic Slip Time Function Wiring Plan

The kinematic slip time function classes declared at
`include/domain/geomechanics/PyLithFault.hpp:147-321` (`SlipTimeFnStep`,
`SlipTimeFnRamp`, `SlipTimeFnBrune`, `SlipTimeFnLiu`, `SlipTimeFnTimeHistory`) are
not consulted by `f0_prescribed_slip` at `src/physics/CohesiveFaultKernel.cpp:629-657`.
The current callback reads a single Cartesian slip vector from constants slots 28-30.
The slip magnitude ramp is applied externally by multiplying the constants vector
by a scalar time factor at prestep time, but that computation is not in the repository
for any time function other than Ramp, and Example 04 uses `mode = locked` so is
unaffected by the kinematic wiring.

Proposal to extend without breaking Example 04 or `test_prescribed_slip`:

1. Add three auxiliary subfields per cohesive vertex (scalar):
   `slip_onset`, `slip_rise`, `slip_amplitude`. Populate via the aux DM machinery
   described in C.2. Use the same fault label for restriction.
2. At TSPreStep, walk the `SlipTimeFn` instance owned by the simulator to compute
   the current slip amplitude at each vertex and write it into a fourth aux subfield
   `slip_current[3]` (Cartesian).
3. Change the residual kernel read path from
   `constants[COHESIVE_CONST_PRESCRIBED_SLIP_X + d]` (current, at
   `src/physics/CohesiveFaultKernel.cpp:653-654`) to
   `a[aOff[slip_current_field] + d]`.
4. Gate the behavior on the presence of the aux field: when
   `NfAux == 0`, fall back to reading from constants slots 28-30. This preserves
   `Integration.PrescribedSlip.CallbackEnforcesPrescribedJump`
   (`tests/integration/test_prescribed_slip.cpp:33-72`), which constructs a
   synthetic `uOff`, `constants`, `f` triple with no aux data, and Example 04 which
   is `mode = locked` (no prescribed slip read path exercised).

The `SlipTimeFn::computeSlip(t, slip_time, amplitude)` API defined at
`include/domain/geomechanics/PyLithFault.hpp:119-128` is already sufficient and
does not require new abstractions. Step / Ramp / Brune / Liu implementations exist.
Only the TimeHistory implementation needs populating from a spatial database or
ASCII file; the declaration at `include/domain/geomechanics/PyLithFault.hpp:301-321`
already anticipates this.

## E. Poroelastic Fault Coupling Plan

Highest architectural risk because it compounds with Section B. The `FaultCohesiveTHM`
class at `include/domain/geomechanics/PyLithFault.hpp:1432-1541` declares the methods
`addPoroelasticResidual`, `addThermalResidual`, `addFlowResidual` but has no
implementation file instantiated in the PetscDS pipeline. The pressure field on the
fault must satisfy its own constraint across the cohesive interface.

### E.1 Field layout options

Option F1. Second Lagrange field (one vector for traction, one scalar for pressure
discontinuity). Advantages: clean decomposition; each Lagrange multiplier has a
physically distinct conjugate constraint (displacement jump versus pressure jump or
fault-normal fluid flux balance). Disadvantages: adds a new field to the solution
vector, which extends the Jacobian stencil and interacts with the cohesive penalty
assembly at `src/core/Simulator.cpp:4182-4591`. Requires a new `registerWithDS` entry
point.

Option F2. Second component of the existing Lagrange field (extend the Lagrange
vector from dim to dim+1). Advantages: no new field declaration; the existing
penalty assembly naturally extends. Disadvantages: conflates a vector traction with
a scalar pressure constraint, making the Lagrange field heterogeneous in physical
meaning; complicates the displacement-side coupling because only the first `dim`
components of the Lagrange vector are tractions.

Option F1 is recommended despite the larger stencil impact, because the cleaner
decomposition aligns with the plan in Section B (manual assembly per Lagrange
component) and avoids special-casing the residual kernels by component index.

### E.2 Constraint equations

For Option F1, with `lambda_mech` the existing vector Lagrange multiplier
(the fault traction) and `lambda_flow` a new scalar Lagrange multiplier:

- Mechanical constraint (unchanged):
  `C_mech = (u+ - u-) - delta` for prescribed slip, or the Coulomb or slip-weakening
  or rate-state constraint for slipping modes.
- Flow constraint (new): one of the following, depending on the fault zone model:
  - Impermeable: `C_flow = p+ - p-` (continuous pressure across a sealed fault).
  - Leaky: `C_flow = n . (q+ - q-) - L * (p+ - p-)` where `L` is a leakage
    conductance proportional to fault transmissibility `k_normal /
    (mu * fault_thickness)` from `FaultHydraulicProperties`
    (`include/domain/geomechanics/PyLithFault.hpp:1050-1082`).

### E.3 Leak-off term

The leak-off term from `FaultPoroelastic::diffusivePressure`
(`include/domain/geomechanics/PyLithFault.hpp:1582-1585`) or an algebraic form
`q_leak = k_n / mu * (p_frac - p_matrix) / half_thickness` is added to the pressure
residual on both sides of the cohesive cell, with opposite signs. The resulting
assembly pattern matches the existing cohesive displacement assembly in
`addCohesivePenaltyToJacobian`, so the code path can be replicated with field
indices remapped from `disp_field` to the pressure field
(fluid model SINGLE_COMPONENT has the pressure at field 0, see
`src/core/Simulator.cpp:4215-4226`).

### E.4 Why this compounds with Section B

The same penalty damping that currently misbehaves for prescribed slip (Section B,
"Why PrescribedSlipQuasiStatic fails") will misbehave for the flow Lagrange
multiplier. Adding Option F1 before committing to Section B Option B produces a
second failing test with the same root cause and no independent progress. The
correct ordering is: Section B Option B, then rate-state (Section C), then
poroelastic coupling (this section).

## F. Known Failure: PrescribedSlipQuasiStatic

The failing test is `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic`
(`tests/integration/test_dynamic_rupture_solve.cpp:418-435`). The expected result is
`max_fault_slip > 5e-4` m with the configured 1 mm opening. The observed result per
CLAUDE.md is near-zero slip.

Diagnosis against the four options in `docs/LAGRANGE_FIX_STATUS.md`:

- Option A (upgrade PETSc to main branch, lines 92-95). Would resolve the failure
  indirectly by enabling region-restricted Lagrange fields, eliminating the volume
  regularization and its interaction with the penalty. Prohibited by the Docker
  build constraint in CLAUDE.md rule 1 without a companion Dockerfile update.
- Option B (manual Lagrange assembly, lines 97-101). Would resolve the failure
  directly by removing the regularization residual on interior cells and by
  allowing the cohesive BdResidual to drive the Lagrange multiplier to the correct
  prescribed-slip traction without competing against the penalty-damped diagonal.
  This is the recommended path.
- Option C (penalty dominance, lines 103-106). Cannot resolve the failure without
  ad-hoc penalty tuning, which is forbidden. The epsilon sweep at
  `docs/LAGRANGE_FIX_STATUS.md:108-129` confirms no scalar multiplier achieves both
  LU stability and Newton convergence.
- Option D (FieldSplit, lines 150-152). Does not resolve the failure because the
  residual imbalance exists before the linear solve; splitting the factorization
  does not change which system is being factored.

Option B resolves the failure. No other option is viable without violating a
project rule or requiring ad-hoc tuning.

## G. Rules Compliance

Checked against CLAUDE.md rules 4, 5, 6, 12, 14, 15. Each section above respects
each rule as follows:

- Rule 4 (never change DS/BC ordering in setupFields). Sections B, C, D, E do not
  reorder the five-step pattern `DMCreateDS -> setupBoundaryConditions ->
  DMSetLocalSection(nullptr) -> DMSetUp -> DMGetDS`
  (`src/core/Simulator.cpp:1703-1713`). All new PetscDS state writes (new constant
  slots, new aux field registrations) happen in `setupPhysics` and the aux DM
  construction path, both of which run after `DMSetUp`.
- Rule 5 (do not modify callback math in PetscFEElasticity.cpp,
  PetscFEPoroelasticity.cpp, PetscFEFluidFlow.cpp). Section C adds new reads from
  the aux vector inside `CohesiveFaultKernel` callbacks only. Section D adds reads
  inside `f0_prescribed_slip` only. Section E adds a new `FaultCohesiveTHM`-style
  callback file rather than editing the three prohibited files.
- Rule 6 (do not modify FaultMeshManager::splitMeshAlongFault or
  CohesiveFaultKernel::registerWithDS). Section C.1 adds new constant slots via a
  new entry point `registerRateStateWithDS`, not by editing `registerWithDS` at
  `src/physics/CohesiveFaultKernel.cpp:54-99`. Section C.2 uses the existing
  fault label without changing `splitMeshAlongFault` at
  `src/numerics/FaultMeshManager.cpp:147-229`. Section E adds a new registration
  function for the flow Lagrange multiplier rather than extending the mechanical
  one.
- Rule 12 (no em dashes, no contractions). This document uses hyphens and double
  hyphens, and "does not" and "cannot" in place of contractions.
- Rule 14 (update CLAUDE.md and README.md to reflect actual code state after every
  session). Session 0 only updates this file (`docs/PYLITH_COMPATIBILITY.md`)
  because the rule is scoped to actual code changes. Sessions that implement
  Section B, C, D, or E must update CLAUDE.md and README.md in the same session,
  with every claim cited to a test name or code reference.
- Rule 15 (GTEST_SKIP only for hardware-dependent tests or genuine crash bugs).
  Section F diagnoses `PrescribedSlipQuasiStatic` without proposing GTEST_SKIP.
  The failure is a physics setup issue traceable to the Lagrange architecture
  choice in Section B, and must be fixed, not hidden.
