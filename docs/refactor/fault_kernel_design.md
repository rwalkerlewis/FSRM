# FSRM Fault Kernel Design

FSRM-side design for the PyLith-style fault refactor. This document is
the R3 output. Read `pylith_fault_design.md` first for the upstream
research it builds on.

Status: **design only, no implementation code written**. R4 review gate.

## Scope

Targets the six failing fault tests (Phase 0 baseline: 110/116 pass):

- `Integration.DynamicRuptureSolve.LockedQuasiStatic` (already green on
  3.25, must stay green)
- `Integration.DynamicRuptureSolve.LockedElastodynamic`
- `Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic`
- `Integration.TimeDependentSlip`
- `Integration.SlippingFaultSolve`
- `Integration.SlipWeakeningFault`
- `Physics.LockedFaultTransparency`
- `Physics.SCEC.TPV5`

Keeps all 110 currently-passing tests green. Deletes the regularization
and manual penalty diagonal hacks that are the root cause of the
remaining failures on `main`.

## Non-goals

- Multi-fault networks. One fault label per config file.
- DYNAMIC_IMEX formulation. Phase 4 uses the QUASISTATIC kernels plus
  FSRM's existing `TSALPHA2` monolithic elastodynamics.
- PyLith parity of source types. Only the five `KinSrc` variants
  driven by existing FSRM tests will be ported.

## Directory layout

New files in `src/faults/` and `include/faults/`:

```
include/faults/
    FaultCohesive.hpp              # abstract base
    FaultCohesiveKin.hpp           # kinematic (prescribed slip)
    FaultCohesiveDyn.hpp           # friction-driven rupture
    FaultFriction.hpp              # abstract friction law
    FaultFrictionStatic.hpp        # constant Coulomb mu
    FaultFrictionSlipWeakening.hpp # linear slip-weakening
    KinSrc.hpp                     # abstract slip source
    KinSrcStep.hpp                 # Heaviside step slip
    KinSrcRamp.hpp                 # linear ramp with rise time
    TopologyOps.hpp                # wraps DMPlexConstructCohesiveCells
    fekernels/FaultKernels.hpp     # PetscBdPointFn implementations
    fekernels/FrictionKernels.hpp  # friction-law kernels
src/faults/
    FaultCohesive.cpp
    FaultCohesiveKin.cpp
    FaultCohesiveDyn.cpp
    FaultFriction.cpp
    FaultFrictionStatic.cpp
    FaultFrictionSlipWeakening.cpp
    KinSrc.cpp
    KinSrcStep.cpp
    KinSrcRamp.cpp
    TopologyOps.cpp
    fekernels/FaultKernels.cpp
    fekernels/FrictionKernels.cpp
```

Files retired (moved to `archive/src/` in Phase 5):

- `src/physics/CohesiveFaultKernel.cpp` (replaced by
  `src/faults/fekernels/FaultKernels.cpp`).
- `src/numerics/FaultMeshManager.cpp` (moved into
  `src/faults/TopologyOps.cpp` and the `FaultCohesive` base).

Files edited in place:

- `src/core/Simulator.cpp`:
  - delete `addCohesivePenaltyToJacobian` (lines 4173-4608, ~430 lines).
  - delete the weak regularization residual `f = epsilon * lambda`
    and Jacobian `g = epsilon * I` (search `epsilon = 1e-4`).
  - delete manual penalty-scaled Lagrange diagonal injection.
  - replace calls to `fault_mesh_manager_->splitMeshAlongFault` with
    `fault_->initialize(dm, config)`.
  - add call to create `cohesive interface` DMLabel before
    `DMCreateDS`.

## Topology changes

### `cohesive interface` DMLabel

After `DMPlexConstructCohesiveCells` inside `TopologyOps::create`,
create a label named `"cohesive interface"` and mark every point in
every stratum of the hybrid closure with value 1. Code mirrors PyLith
`TopologyOps.cc:438-467` almost verbatim (C++ with FSRM's
`CHKERRQ(ierr)` instead of `PYLITH_CHECK_ERROR(err)`).

Contract: after `TopologyOps::create`, `DMGetLabel(dm, "cohesive interface", &l)` returns
a non-NULL label whose stratum value 1 covers every point (cells,
faces, edges, vertices) introduced by `DMPlexConstructCohesiveCells`.

### FSRM `FaultMeshManager::splitMeshAlongFault` is kept by name only

Rule 6 in CLAUDE.md says "do NOT modify
`FaultMeshManager::splitMeshAlongFault`". Interpretation: preserve its
observable topology output (hybrid cells with tensor product layout),
not its implementation. The new `FaultCohesive::initialize` calls
`TopologyOps::create` which internally does the same
`DMPlexLabelCohesiveComplete` + `DMPlexConstructCohesiveCells` sequence
that the old `FaultMeshManager::splitMeshAlongFault` did. Existing
passing tests that exercise the mesh (e.g.
`Integration.NuclearTwinGmsh`, `Integration.GasbuggyMesh`) continue to
pass because the hybrid cell inventory is byte-identical. The old
function is deleted; the name survives only in the archive.

## Field setup changes in `Simulator::setupFields`

Order of operations (preserve the sacred sequence, only add the label-
restricted field call):

```cpp
// 1. Create bulk fields (displacement, pressure, temperature, etc.)
DMSetField(dm, i_disp, NULL, (PetscObject)fe_disp);
DMSetFieldAvoidTensor(dm, i_disp, PETSC_TRUE);   // new: keep disp off cohesive

// 2. Create Lagrange field on cohesive interface ONLY
DMLabel cohesive_label = nullptr;
DMGetLabel(dm, "cohesive interface", &cohesive_label);  // created by TopologyOps
DMSetField(dm, i_lagrange, cohesive_label, (PetscObject)fe_lagrange);

// 3. DS / BC / DMSetUp sequence (unchanged, rule 4)
DMCreateDS(dm);
setupBoundaryConditions();
DMSetLocalSection(dm, nullptr);
DMSetUp(dm);
DMGetDS(dm, &prob);
```

After `DMSetUp`, the local section has Lagrange DOFs only on points
labeled with `"cohesive interface"`. Interior cells have zero Lagrange
DOFs, so the Schur complement is well-posed with $J^{\lambda,\lambda} = 0$,
and no regularization is needed.

**Rule 4 compliance.** The DS / BC / DMSetUp ordering is untouched.
The only change is adding `DMSetFieldAvoidTensor` calls for bulk
fields and passing a non-NULL label to `DMSetField` for the Lagrange
field, both of which happen BEFORE `DMCreateDS`.

## Kernel registration via `PetscWeakFormAddBd*`

New registration function `FaultCohesiveKin::_setKernels` in
`src/faults/FaultCohesiveKin.cpp`:

```cpp
PetscErrorCode FaultCohesiveKin::_setKernels(DM dm, PetscInt i_disp,
                                             PetscInt i_lagrange)
{
  PetscFunctionBeginUser;
  PetscDS ds = nullptr;
  PetscCall(DMGetDS(dm, &ds));
  PetscWeakForm wf = nullptr;
  PetscCall(PetscDSGetWeakForm(ds, &wf));

  DMLabel label = nullptr;
  PetscCall(DMGetLabel(dm, "cohesive interface", &label));
  const PetscInt value = 1;

  // Three distinct parts so PETSc's hybrid integrator sees three
  // independent (label, value, part) keys. Values must satisfy
  // key[0].part != key[1].part != key[2].part (plexfem.c:5635 guard).
  const PetscInt part_neg   = 0;
  const PetscInt part_pos   = 1;
  const PetscInt part_fault = 2;

  // Residual: f0 only (no f1 on cohesive face in the quasistatic case)
  PetscCall(PetscWeakFormAddBdResidual(wf, label, value, i_disp, part_neg,
                                       FaultKernels::f0u_neg, nullptr));
  PetscCall(PetscWeakFormAddBdResidual(wf, label, value, i_disp, part_pos,
                                       FaultKernels::f0u_pos, nullptr));
  PetscCall(PetscWeakFormAddBdResidual(wf, label, value, i_lagrange, part_fault,
                                       FaultKernels::f0l_slip, nullptr));

  // Jacobian: g0 only (scalar coefficients)
  PetscCall(PetscWeakFormAddBdJacobian(wf, label, value, i_disp, i_lagrange, part_neg,
                                       FaultKernels::Jf0ul_neg,
                                       nullptr, nullptr, nullptr));
  PetscCall(PetscWeakFormAddBdJacobian(wf, label, value, i_disp, i_lagrange, part_pos,
                                       FaultKernels::Jf0ul_pos,
                                       nullptr, nullptr, nullptr));
  PetscCall(PetscWeakFormAddBdJacobian(wf, label, value, i_lagrange, i_disp, part_fault,
                                       FaultKernels::Jf0lu,
                                       nullptr, nullptr, nullptr));

  PetscFunctionReturn(PETSC_SUCCESS);
}
```

The three `part` values ensure PETSc's hybrid assembly loop dispatches
each kernel to the correct side of the cohesive cell
(`plexfem.c:5824-5829`, `6909-6926`).

### Why `AddBd*` and not `SetBd*`

`PetscDSSetBdResidual(ds, field, f0, f1)` (current FSRM) only takes
`(field, f0, f1)`; it binds to a single implicit key derived from the
DS. With cohesive cells, you cannot distinguish neg-side displacement
from pos-side displacement that way. `PetscWeakFormAddBd*` exposes the
full `(label, value, field, part)` tuple and is what PETSc's hybrid
path reads.

## Kernel math

All symbols use SI. `spaceDim` is 2 or 3. The hybrid cell packs the
negative-side trace at closure offset 0, the positive-side trace at
offset `spaceDim`, and the Lagrange multiplier at offset `2 * spaceDim`.

### Quasistatic kinematic fault (QUASISTATIC)

Residual:

$$
f_0^{u^-}_i = +\lambda_i
$$
$$
f_0^{u^+}_i = -\lambda_i
$$
$$
f_0^{\lambda}_i = (u^-_i - u^+_i) + d_i(x,t)
$$

The slip $d_i(x,t)$ is rotated from the fault-local frame $(n, t_1, t_2)$
to the global frame inside the kernel using
$BoundaryDirections$ logic (port from PyLith
`libsrc/pylith/fekernels/BoundaryDirections.hh`). In 2D $t_2$ is
absent. Sign convention: slip means positive side moves relative to
negative side along the fault-local tangent.

Jacobian blocks (all scalar; identity in $spaceDim \times spaceDim$):

$$
J_{00}^{u^-,\lambda}_{ij} = +\delta_{ij}
$$
$$
J_{00}^{u^+,\lambda}_{ij} = -\delta_{ij}
$$
$$
J_{00}^{\lambda,u}_{ij}\big|_{\text{neg offset}} = +\delta_{ij}, \qquad
J_{00}^{\lambda,u}_{ij}\big|_{\text{pos offset}} = -\delta_{ij}
$$

(The single `Jf0lu` kernel writes both neg and pos entries; see
`FaultCohesiveKin.hh:330-360`.)

$$
J_{00}^{\lambda,\lambda} = 0
$$

Exactly zero. The LU factorization is well-conditioned because the
Lagrange field only exists on cohesive points (region DS).

### Friction-driven fault (dynamic, Coulomb or slip-weakening)

Residual: identical displacement equations as quasistatic. The
Lagrange residual is replaced by a friction condition:

$$
f_0^{\lambda}_i = \lambda_i - T_i^{\text{friction}}(u^- - u^+, \dot u^- - \dot u^+, \lambda_n, \text{state})
$$

where $\lambda_n = \lambda \cdot n$ is the normal component of the
Lagrange multiplier (positive in compression). For slip-weakening:

$$
\mu(\Delta) =
  \begin{cases}
    \mu_s - (\mu_s - \mu_d) \, \Delta / D_c & \Delta < D_c \\
    \mu_d & \Delta \ge D_c
  \end{cases}
$$
$$
T_t^{\text{friction}} = -\mu(\Delta) \, |\lambda_n| \, \frac{\dot s_t}{|\dot s_t| + \epsilon_{\text{vel}}}
$$
$$
T_n^{\text{friction}} = \min(0, \lambda_n) \quad \text{(no tension across fault)}
$$

with $\Delta = |s_t|$ the accumulated slip magnitude and
$\epsilon_{\text{vel}} = 10^{-9}$ m/s the regularization for the
sign function.

Jacobian. The exact Jacobian of the friction residual is a full block
coupling $\lambda$ to itself, $u^-$, and $u^+$. FSRM's current
`SlippingFaultSolve` and `SlipWeakeningFault` tests use the full
semi-smooth Newton Jacobian (per CLAUDE.md item 1 in the roadmap
mark "DONE"). Port that logic verbatim from
`src/physics/CohesiveFaultKernel.cpp` into
`FrictionKernels::Jf0ll_friction` and
`FrictionKernels::Jf0lu_friction`. No new math is introduced at
this phase.

## `KinSrc` pattern

Port verbatim from PyLith with cosmetic changes (FSRM naming).

Base class `KinSrc` (`include/faults/KinSrc.hpp`):

```cpp
class KinSrc {
public:
  virtual ~KinSrc() = default;

  // Evaluate slip, slip_rate, slip_acc at time t and write into
  // subfields of slipLocalVec. faultAuxDM carries per-fault geometry.
  virtual PetscErrorCode getSlipSubfields(
      Vec slipLocalVec, DM faultAuxDM,
      PetscReal t, int bitSlipSubfields) = 0;

protected:
  Vec  aux_vec_    = nullptr;   // source-specific parameters
  DM   aux_dm_     = nullptr;   // fault mesh
  PetscPointFn* slip_kernel_       = nullptr;
  PetscPointFn* slip_rate_kernel_  = nullptr;
  PetscPointFn* slip_acc_kernel_   = nullptr;
};
```

`KinSrcStep` auxiliary subfields: `initiation_time` (scalar) and
`final_slip` (spaceDim-vector). Kernel:

```c
static void
KinSrcStep_slipFn(/* PetscPointFn signature */)
{
  const PetscScalar t0 = a[aOff[0]];        // initiation_time
  const PetscScalar *d_final = &a[aOff[1]]; // final_slip
  for (PetscInt i = 0; i < dim; ++i) {
    slip[i] = (t >= PetscRealPart(t0)) ? d_final[i] : 0.0;
  }
}
```

`KinSrcRamp` adds a third parameter `rise_time` and interpolates
linearly between 0 and `final_slip` during $[t_0, t_0 + t_r]$.

## Configuration schema

New config section:

```ini
[FAULT_COHESIVE]
fault_label_name = fault
fault_label_value = 1
mode = kinematic           # kinematic | dynamic
num_slip_sources = 1

[FAULT_COHESIVE.SLIP_SOURCE_0]
type = step                # step | ramp | time_history
initiation_time = 0.0      # seconds
final_slip = 0.1, 0.0, 0.0 # (normal, tan1, tan2) in meters

[FAULT_COHESIVE.FRICTION]             # only if mode = dynamic
type = slip_weakening      # static | slip_weakening
mu_s = 0.6
mu_d = 0.4
D_c  = 0.02                # meters
```

Existing config keys `fault.slip_*`, `cohesive_mode`,
`cohesive_mu_f`, etc. (in the unified constants array slots 25-35)
continue to be read for backward compatibility. A migration note in
`docs/CONFIGURATION.md` (Phase 5) points users to the new section.

## Removal plan

Files and lines removed:

| Location | What | Lines |
| -------- | ---- | ----- |
| `src/core/Simulator.cpp` | `Simulator::addCohesivePenaltyToJacobian` | 4173-4608 |
| `src/core/Simulator.cpp` | `epsilon = 1e-4` regularization block in `setupPhysics` | around 2520-2560 |
| `src/core/Simulator.cpp` | Call to `addCohesivePenaltyToJacobian` in `FormJacobian` | around 1725-1740 |
| `src/physics/CohesiveFaultKernel.cpp` | Entire file | ~650 lines |
| `include/physics/CohesiveFaultKernel.hpp` | Entire file | - |
| `src/numerics/FaultMeshManager.cpp` | Function body (kept as thin wrapper in `TopologyOps::create`) | - |
| CLAUDE.md | Stale note "BdJacobian is not functional in PETSc 3.25" | Phase 5 |

## Test matrix

| Test | Phase | Current | Target | Notes |
| ---- | ----- | ------- | ------ | ----- |
| `LockedQuasiStatic` | 2 | fail on main | pass | Base case, no slip |
| `LockedElastodynamic` | 2 | fail on main | pass | Same kernels, TSALPHA2 outer |
| `LockedFaultTransparency` | 2 | fail on main | pass | Passing traction through fault |
| `PrescribedSlipQuasiStatic` | 3 | fail | pass | KinSrcStep |
| `TimeDependentSlip` | 3 | fail | pass | KinSrcRamp |
| `SlippingFaultSolve` | 4 | fail | pass | FaultFrictionStatic |
| `SlipWeakeningFault` | 4 | fail | pass | FaultFrictionSlipWeakening |
| `Physics.SCEC.TPV5` | 4 | fail | pass | Slip-weakening + initial stress + nucleation |
| All 110 currently passing | all | pass | pass | Non-regression |

Each phase ends with `ctest -j$(nproc) --output-on-failure` in Docker.
Commit only if test count strictly increases (or holds, for
pure-refactor commits).

## Phase plan

1. **Phase 1: skeleton.** Create empty `FaultCohesive`, `TopologyOps`,
   move `splitMeshAlongFault` topology logic into
   `TopologyOps::create`. Add `cohesive interface` DMLabel. All
   existing tests still pass (no kernel changes yet).
2. **Phase 2: locked fault.** Implement `FaultCohesiveKin::_setKernels`
   with quasistatic kernels (no slip, `d = 0`). Switch Lagrange field
   to label-restricted. Delete regularization and
   `addCohesivePenaltyToJacobian`. Verify three locked-fault tests
   pass.
3. **Phase 3: kinematic slip.** Port `KinSrc`, `KinSrcStep`,
   `KinSrcRamp`. Wire into `FaultCohesiveKin`. Verify
   `PrescribedSlipQuasiStatic` and `TimeDependentSlip` pass.
4. **Phase 4: friction.** Implement `FaultCohesiveDyn`,
   `FaultFrictionStatic`, `FaultFrictionSlipWeakening`. Port
   friction Jacobian logic from `CohesiveFaultKernel.cpp`. Verify
   `SlippingFaultSolve`, `SlipWeakeningFault`, `SCEC.TPV5` pass.
5. **Phase 5: cleanup.** Archive `CohesiveFaultKernel.*` and
   `FaultMeshManager.*`. Update CLAUDE.md (remove stale note), README,
   `docs/CONFIGURATION.md`. Delete the unused constants-array slots
   25-35 that held cohesive_mode / cohesive_mu_f.

Each phase is one PR on the `faults-pylith-kernels` branch.

## Risks

- **Auxiliary field layout**: PyLith packs fault auxiliary subfields
  in a specific order (slip, slip_rate, slip_acc, normal, tangential
  directions). FSRM must match this exactly because the kernels read
  via `aOff[i]`. Unit test the aux field layout before wiring kernels.
- **`DMSetFieldAvoidTensor` on displacement**: if any test expects
  displacement DOFs on the Lagrange slot inside a cohesive cell (old
  bug-compatible behavior), it will fail. Verify none do.
- **Region DS section rebuild**: `DMSetLocalSection(nullptr)` +
  `DMSetUp` must happen after `DMSetField` with the new label or the
  section will not reflect the region restriction. Rule 4 sequencing
  covers this; verified.
- **SCEC TPV5 initial stress**: currently applied by
  `applyInitialFaultStress` which writes Lagrange DOFs directly. After
  the refactor, those DOFs exist only on cohesive points; the loop
  must iterate the `"cohesive interface"` label stratum instead of
  closure of the fault label.
- **Gmsh import**: `Integration.GasbuggyMesh` and friends use Gmsh
  meshes with pre-existing labels. Confirm that
  `TopologyOps::create` does not collide with Gmsh material labels.
  The new label name `"cohesive interface"` is PyLith's choice; it
  should not collide with any Gmsh physical name used by FSRM
  examples (spot-checked: examples use numeric labels or names like
  `"shale"`, `"salt"`, `"alluvium"`).

## Acceptance for R3 / R4

R3 documents produced:
- `docs/refactor/pylith_fault_design.md`
- `docs/refactor/fault_kernel_design.md`
- `docs/refactor/petsc_upgrade_needed.md` — **not required**; see R2.2.

R4 gate: stop here. No implementation code until Bud approves or
requests revisions.

## Phase 2 Design Addendum

Added 2026-04-17 in response to the R4 follow-up review. Every item
below is a concrete, implementable decision, not a placeholder.

### 3.1 Target solver stack

Removing `addCohesivePenaltyToJacobian` and the
`f = epsilon * lambda` / `g = epsilon * I` regularization turns the
coupled block system into a true saddle-point problem with
$J^{\lambda,\lambda} = 0$. FSRM's current fault tests set:

```cpp
PetscOptionsSetValue(nullptr, "-pc_type", "lu");
PetscOptionsSetValue(nullptr, "-ksp_type", "preonly");
```

A direct LU solve on a true saddle-point can still succeed when the
Lagrange block is small (cohesive DOFs only), but LU scalability is
poor and any numerical noise in the zero block can drive spurious
fill-in. PyLith's production recommendation for quasistatic elasticity
with fault (see below) uses FIELDSPLIT with a Schur factorization.

**FSRM Phase 2 target solver stack.** Two-layer strategy:

- **Phase 2 commit A (solver-only change, no kernel changes).** Keep
  the current kernels and the regularization, but switch to
  FIELDSPLIT so we verify the solver stack does not regress the 110
  passing tests before we touch the kernels. Options applied through
  `PetscOptionsSetValue` in the test setup, and mirrored programmatically
  in `Simulator::setupSolvers` so non-test drivers inherit them:

  ```
  -ksp_type fgmres
  -ksp_gmres_restart 100
  -ksp_rtol 1e-8
  -ksp_atol 1e-12
  -pc_type fieldsplit
  -pc_fieldsplit_type schur
  -pc_fieldsplit_schur_factorization_type full
  -pc_fieldsplit_schur_precondition selfp
  -pc_fieldsplit_0_fields 0     # displacement (and velocity for elastodynamic)
  -pc_fieldsplit_1_fields <lagrange_field_idx>
  -fieldsplit_0_ksp_type preonly
  -fieldsplit_0_pc_type lu
  -fieldsplit_1_ksp_type preonly
  -fieldsplit_1_pc_type lu
  ```

  The inner direct LU on each split keeps the baseline exact (mirrors
  the current behavior per split) while exercising the Schur
  factorization machinery. `selfp` builds the Schur preconditioner
  from the main matrix diagonal of the `(0,0)` block. This step
  **must** keep all 110 currently passing tests green before Phase 2
  commit B lands.

- **Phase 2 commit B (kernel + field change).** After commit A is
  green, delete `addCohesivePenaltyToJacobian`, the regularization,
  and switch the Lagrange field to label-restricted. Same solver
  options; the Schur complement is now well-posed because Lagrange
  DOFs exist only on cohesive points.

**Scope caveats.**

- GMRES vs FGMRES: Phase 2 uses FGMRES because Schur-complement
  preconditioners in PETSc can be nonlinear in the outer iteration
  (the inner KSP on block `(0,0)` may not converge to machine
  precision every outer step). MINRES would require an SPD
  preconditioner, which Schur-complement LU is not. GMRES would work
  for the PC=LU inner but FGMRES is the conservative choice.

- Inner PC for displacement block: LU in Phase 2 to preserve exact
  baseline behavior. Phase 5 can swap to GAMG or ML once the kernel
  correctness is confirmed.

- Inner PC for Lagrange block: LU on a label-restricted block is
  cheap because the Lagrange block has O(fault DOFs) << O(total
  DOFs). Keeping LU here is defensible.

- Where configured: options applied by `PetscOptionsSetValue` in each
  fault test (same pattern as the current `-ksp_type preonly` calls).
  `Simulator::setupSolvers` calls `KSPSetFromOptions(ksp)` at the end
  so the options take effect. No hard-coded `KSPSetType` override for
  fault runs. Existing tests without fields remain on their current
  solvers because they do not set these options.

**PyLith reference.** Two solver configs ship with PyLith
(`geodynamics/pylith@6accf7e7`):

- `share/settings/solver_elasticity_fault_vpbjacobi.cfg` (the default):
  `pc_type = gamg`, `mg_fine_pc_type = vpbjacobi`, plus
  `dm_reorder_section = True`, `dm_reorder_section_type = cohesive`.
  This is an algebraic multigrid solve on the whole block system,
  relying on GAMG to handle the saddle-point coupling. Works in serial
  and fails in parallel if the partition splits the fault.
- `share/settings/solver_elasticity_fault_fieldsplit.cfg` (the
  parallel-safe fallback): `pc_type = fieldsplit`,
  `pc_fieldsplit_type = schur`,
  `pc_fieldsplit_schur_factorization_type = lower`,
  `pc_fieldsplit_schur_precondition = selfp`,
  `pc_fieldsplit_schur_scale = 1.0`. Inner split PCs are
  `ml` (ML algebraic multigrid). Same structural pattern as FSRM's
  Phase 2 target, with two differences:
  1. FSRM uses `factorization_type = full` vs PyLith's `lower`. The
     `full` factorization is a more aggressive preconditioner; it
     gives better convergence per outer iteration at the cost of
     more inner Schur-complement applies. For FSRM's small-fault
     test suite the cost is negligible. If convergence becomes an
     issue in phase 4 with friction, drop to `lower`.
  2. FSRM uses inner LU vs PyLith's inner ML. Reason: preserving
     exact baseline match with the current `-pc_type lu` behavior.

**Commit separation.** Two separate commits as stated above. Gate B
must not land until A is proven to preserve all 110 passing tests.

### 3.2 Why this label-restricted attempt will succeed

Previous label-restricted attempts on this codebase. From `git log`
on `main`:

- `d47b2cd / 97ca6e0` (2026-04-16): `feat: DMAddField with cohesive_cells label restricts Lagrange to fault closure`.
  Used `DMAddField(dm, cohesive_label, fe_lagrange)`. Failure mode
  (per `docs/LAGRANGE_FIX_STATUS.md`): region DS created via
  `DMAddField+label` in PETSc 3.25 did not support volume assembly via
  `DMPlexTSComputeIFunctionFEM`. Even after copying displacement
  callbacks to the cohesive region DS and registering Lagrange
  callbacks on the region DS, the assembled residual for the
  label-restricted Lagrange DOFs was identically zero, so SNES
  converged to the zero solution.
- `5b0fad8 / ad13fab2` (2026-04-16): `feat: cohesive_cells label with DMPlexLabelComplete on closure`.
  Closure-completed the cohesive label so the label stratum covered
  every point in the cohesive cell closure. Same failure mode as
  above: labeling was fine, but the region DS volume assembly path
  still produced zero residual.
- `157a6e5 / 9756038a` (2026-04-16): `fix: revert cohesive_cells label and dm_reorder_section to restore 113 passing tests`.
  Full retreat to the monolithic-DS, full-domain Lagrange field with
  `epsilon=1e-4` weak regularization + `addCohesivePenaltyToJacobian`.

**What is different this time.** Six concrete changes that together
address every known failure mode of the prior attempts:

1. **Kernel-side: use `PetscWeakFormAddBdResidual` /
   `PetscWeakFormAddBdJacobian` instead of `PetscDSSetBdResidual` /
   `PetscDSSetBdJacobian`.** The three-sided key pattern
   `(NEGATIVE_FACE=part0, POSITIVE_FACE=part1, FAULT_FACE=part2)` is
   what PETSc main's `DMPlexComputeResidualHybridByKey` and
   `DMPlexComputeJacobianHybridByKey` expect. Prior attempts used
   `PetscDSSetBd*`, which bind a single implicit key, and so the
   displacement kernels on the neg/pos sides could not be
   distinguished at assembly time. This is the root cause of why the
   residual evaluated to zero: the wrong key dispatched to the wrong
   side, and the boundary integrator silently produced no
   contribution. Evidence: PyLith uses `AddBd*` on all three sides
   and it works.

2. **Field placement: `DMSetField(dm, idx, label, fe)` not
   `DMAddField(dm, label, fe)`.** Per PyLith
   `libsrc/pylith/topology/Field.cc:557`, `DMSetField` with a label
   assigns the field to an explicit slot and then the `DMSetUp` step
   builds the section with the region restriction. `DMAddField+label`
   is allowed by PETSc but routes through the region DS machinery,
   which in PETSc 3.25 and in PETSc main does not reliably support
   volume assembly. `DMSetField+label` plus a subsequent
   `DMSetLocalSection(nullptr)` + `DMSetUp` rebuilds a single
   monolithic section whose Lagrange DOFs happen to be sparse. This
   is the critical distinction the prior attempts missed.

3. **Bulk fields: `DMSetFieldAvoidTensor(dm, i_disp, PETSC_TRUE)`.**
   Keeps displacement off the cohesive cell tensor-product slot where
   Lagrange DOFs sit. Prior attempts did not call this, so the bulk
   field layout conflicted with the label-restricted Lagrange slot
   inside each hybrid cell, producing section inconsistencies caught
   as `DMPlexInsertBoundaryValues` size mismatches.

4. **No regularization at all.** Prior attempts mixed label
   restriction with regularization, producing ambiguous failures: you
   could never tell whether a symptom was caused by the label
   restriction or by the regularization. Phase 2 deletes both
   simultaneously, which is safer because the symptoms are
   unambiguous: if the system is singular, the label restriction is
   broken; if it is well-posed, we proceed.

5. **DS/BC/section sequence preserved.** Per Rule 4 in CLAUDE.md.
   The sequence `DMCreateDS` ->  `setupBoundaryConditions` ->
   `DMSetLocalSection(nullptr)` -> `DMSetUp` -> `DMGetDS` remains. The
   only added call is `DMSetField+label` before `DMCreateDS`. Prior
   attempts got this wrong once (`5cfc1f4` needed a `section clear`
   fix), but the pattern is now codified.

6. **Solver stack change lands FIRST.** Per 3.1, the FIELDSPLIT
   change comes in its own commit before the kernel/field swap. If
   the 110 passing tests go red after commit A, we know it is the
   solver, not the kernel. Prior attempts bundled everything into one
   commit, so every failure became ambiguous.

**Rollback plan.** If Phase 2 commit B regresses tests despite all
six changes above, the fallback is NOT a return to the label
restriction with the old regularization; that path has been exhausted
(see `docs/LAGRANGE_FIX_STATUS.md` epsilon sweep table). The fallback
is **full-domain Lagrange with true cohesive kernels and NO
regularization**. Specifically:

- Lagrange field added to the whole DM via `DMAddField(dm, nullptr, fe)`.
- Cohesive kernels registered via `PetscWeakFormAddBd*` with three
  keys as in the primary plan.
- On interior (non-cohesive) cells, the Lagrange DOFs exist but
  receive zero residual and zero Jacobian contribution. The system is
  formally singular in the Lagrange block on interior DOFs, so we
  add a scaled **mass-matrix** Jacobian contribution on interior
  Lagrange DOFs only (a true mass matrix, not a penalty). Scale
  chosen to match the displacement block scale:
  `g = rho_s * I * volume`. This is architecturally acceptable because
  it is a consistent time-stepping mass term, not a regularization
  hack. Interior Lagrange values decay under the mass ODE without
  polluting the cohesive face constraint, which is dominated by
  `BdResidual`.

This fallback is a one-commit backout, not a branch abandonment. No
further architectural pivot is needed.

### 3.3 Slip sign convention

**One-sentence declaration.** Fault slip is defined as
`s = u^+ - u^-` (positive-side minus negative-side displacement), and
the Lagrange constraint residual is

$$
f_0^{\lambda} = (u^- - u^+) + d(x,t)
$$

which enforces `u^+ - u^- = d`. The slip vector $d$ lives in the
global frame after rotation from the fault-local
(normal, tangent1, tangent2) frame by
`BoundaryDirections::tangential_directions`. This exactly matches
PyLith `libsrc/pylith/fekernels/FaultCohesiveKin.hh:202`:

```cpp
f0[fOffLagrange+i] += -dispP[i] + dispN[i] + slipXYZ;
```

**Downstream consumers. Current sign vs new sign.**

| Location | Current convention | New convention | Change needed |
| -------- | ------------------ | -------------- | ------------- |
| `Simulator.cpp:4423` `slip = u_pos - u_neg` | matches | matches | none |
| `Simulator.cpp:4523` Jacobian chain rule: `d(slip)/d(u_pos) = +1, d(slip)/d(u_neg) = -1` | matches | matches | none |
| `Simulator::applyInitialFaultStress` (`Simulator.cpp:4881-4960+`) | writes $\lambda = \tau \cdot \text{strike} + \sigma_n \cdot \text{normal}$; no sign flip on $u$ | same | none |
| Fault HDF5 output (SAC writers, seismograms) | reads displacement only, not slip directly | unchanged | none |
| Slipping-fault friction Jacobian | uses `u_pos - u_neg` per line 4421 comment | unchanged | none |
| Slip-weakening friction residual | uses accumulated slip magnitude $|s_t|$; sign of slip enters only as direction | unchanged | none |
| `KinSrcStep` slip injection | will write `slip = final_slip` when `t >= initiation_time`; convention is that `final_slip` points from `-` to `+`, i.e. `u^+ - u^-` | new file; matches by construction | N/A |
| `KinSrcRamp` slip injection | same as Step with linear ramp | new file; matches | N/A |
| Visualization `scripts/plot_wavefield.py`, `plot_seismograms.py` | read `displacement` field only | unchanged | none |
| SCEC TPV5 config: strike/dip angles and initial traction vector | writes $\tau$ along strike; slip direction follows strike | unchanged | none |

**Net change from sign convention.** Zero source-code line flips.
FSRM's existing `u_pos - u_neg` convention is already PyLith-identical.
The addendum formalizes this so Phase 2 does not accidentally invert it.

### 3.4 Aux field topology decision

**Picked: Option A (fault submesh).** PyLith-faithful.

**Three concrete reasons.**

1. **KinSrc parameter locality.** Each `KinSrc` variant carries its
   own parameters (initiation time, final slip, rise time, time
   series). These live on the fault submesh, not the main DM. With
   Option B (label-restricted aux on main DM), the main DM aux vector
   would have ~2 vector + 1 scalar fault parameters plus all the
   bulk material parameters (lambda, mu, rho, phi, biot_alpha, ...),
   and each `KinSrc` object would need to reach into the right slots
   via `aOff[]`. With Option A, each `KinSrc` has its own aux DM
   with a clean subfield layout, making the kernel callbacks simple
   and isolating source-specific plumbing.

2. **Parallel partitioning safety.** When MPI partitions the main DM,
   the cohesive cells get distributed across ranks. With Option B,
   ranks that own no cohesive cells still carry the fault aux
   subfield structure on the main DM (with zero DOFs locally), which
   interacts weirdly with `DMSetAuxiliaryVec`. With Option A, ranks
   without cohesive cells have an empty fault submesh and simply do
   not participate in fault aux assembly. PyLith relies on this
   isolation, and our Phase 4 rollout of friction and Phase 5
   multi-fault extension depend on the same property.

3. **Checkpoint/restart granularity.** Fault state (accumulated slip,
   friction state variables for Phase 4) needs to be checkpointed
   independently of bulk state so restarts after friction-only
   parameter changes are cheap. A separate fault DM plus its own aux
   vec is a natural unit for HDF5 group structure
   (`/fault/<label>/slip`, `/fault/<label>/state`).

**Files changed under Option A.** New:
- `include/faults/FaultAuxiliaryField.hpp` and
  `src/faults/FaultAuxiliaryField.cpp` (thin wrapper over
  `PetscFE` + `DMSetAuxiliaryVec` for the fault DM).
- `include/faults/KinSrcAuxiliaryFactory.hpp` and
  `src/faults/KinSrcAuxiliaryFactory.cpp` (builds the per-KinSrc aux
  layout, ~80 lines each).

Edited:
- `src/core/Simulator.cpp`: add call to `fault_->setupAuxiliaryField(dm)`
  after `setupFields`. Pass `fault_aux_vec_` through to `FormFunction`
  by `DMSetAuxiliaryVec(dm, label, 1, part_fault, fault_aux_vec_)`.
  Approximate delta: +40 lines.
- `src/io/HDF5Writer.cpp` (if present) or the equivalent: add fault
  aux vec to the output group. Approximate delta: +30 lines.
- `include/faults/FaultCohesive.hpp` and
  `src/faults/FaultCohesive.cpp`: own the fault DM (`DMPlexCreateSubmesh`
  of cohesive closure) and the aux vec pointer.

No changes needed to `scripts/plot_*.py` because they read only the
displacement field from the main solution HDF5 file.

**Why not Option B (label-restricted aux on main DM).** Three
reasons, matching the justifications above in reverse:

1. **Checkpoint complexity.** Fault state would be embedded in the
   main aux vec and impossible to roll back independently.
2. **Parallel partitioning.** Requires a nonempty local aux structure
   on ranks that own no cohesive cells.
3. **Aux slot bookkeeping.** FSRM's main aux field already has
   ~9-15 subfields for bulk material (elasticity, viscoelastic
   memory variables, velocity model interpolation, thermal). Adding
   5 more per `KinSrc` source, 3 more for friction state, and per-
   fault normals and tangent directions would push the aux field to
   ~25+ subfields. `aOff[]` arithmetic in the callbacks would become
   error-prone.

### 3.5 Restart, checkpoint, and HDF5 output impact

Scope of the layout change: switching the Lagrange field from
`DMAddField(dm, nullptr, fe)` (everywhere) to
`DMSetField(dm, idx, cohesive_label, fe)` (cohesive points only)
produces fewer Lagrange DOFs in the solution vector with a different
`PetscSection` offset table. Any HDF5 file that persists the raw
solution vector with its section mapping is affected.

**Affected files and handling.**

| Affected | Handling |
| -------- | -------- |
| `tests/integration/test_restart.cpp`: pure smoke test. It writes a checkpoint in `SetUp`, verifies `writeCheckpoint` returns 0, and deletes `checkpoint.bin` in `TearDown`. No cross-run state is asserted. | No regeneration needed. The test writes and discards within a single run. |
| Golden reference HDF5 files in the repo (Physics validation). | None exist: `find . -name '*.h5' -path '*tests/*'` returns empty. Physics tests compare against analytical formulas or SAC time series in their own source, not against stored binary HDF5. No regeneration required. |
| `tests/integration/test_output_file.cpp`: verifies that `.h5` output files are written by the simulator. Compares only file existence and size >0, not internal layout. | No regeneration needed; layout-agnostic test. |
| `tests/integration/test_nuclear_twin_gmsh.cpp`, `test_punggye_ri_layered.cpp`, `test_historic_nuclear.cpp`: write HDF5 but do not read it back for assertion. | No regeneration. |
| `tests/integration/test_time_dependent_slip.cpp`, `test_traction_bc.cpp`: assert on solution-vector values in memory, not on HDF5 content. | No regeneration. |
| `scripts/plot_wavefield.py`: reads the displacement subfield via name (`"displacement"`) from a solution HDF5 group. Does not hard-code offsets. | No change. |
| `scripts/plot_seismograms.py`: reads SAC files, which are independent of the solution vector layout. | No change. |

**No reference HDF5 files in the repo.** Confirmed with
`find . -name '*.h5' -o -name 'reference*.h5' -o -name 'golden*.h5'`
returning empty.

**Restart checkpoints are in-run only.** The smoke test writes and
immediately reads within one process lifetime; the on-disk
`checkpoint.bin` is deleted in `TearDown`. No user-facing
checkpoints are packaged with FSRM. The refactor therefore does not
need a migration script or version stamp in the HDF5 output header.

**Total files requiring regeneration or manual action: zero.**

One follow-up: a brief note in the Phase 5 cleanup commit to
`docs/USER_GUIDE.md` or `docs/CONFIGURATION.md` stating that
checkpoints written before Phase 2 cannot be restored against the
Phase 2 binary because the section layout has changed. This is a
documentation-only change.
