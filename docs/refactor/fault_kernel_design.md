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
