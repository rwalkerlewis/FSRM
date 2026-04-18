# Fault Extension Gap Analysis

Session 0, read-only audit. Baseline for the fault-extension work (rate-and-state
friction, kinematic slip time functions, poroelastic fault coupling, and a fix
for Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic).

All file and line references are to the FSRM tree unless prefixed with
`pylith:`. The PyLith tree at /home/dockimble/Projects/pylith is read-only
reference. No build or test was run as part of this audit.

## 1. PetscDS Constant Slot Map

The unified PetscDS constants array is populated in
src/core/Simulator.cpp near setupPhysics (see src/core/Simulator.cpp:2078 for
the documented header and src/core/Simulator.cpp:2216 for the cohesive writes).
CLAUDE.md advertises "up to 80 elements" but there is no hard cap in the code;
the effective upper bound is the `MAX_UNIFIED_CONSTANTS` buffer and whatever
PetscDSSetConstants accepts.

### 1.1 Occupied slots (as of branch local_fix)

| Slot  | Symbol / role                              | Owner                      | Source anchor |
|-------|--------------------------------------------|----------------------------|---------------|
| 0     | lambda (Lame, first)                       | elasticity, poro, aux gate | src/core/Simulator.cpp:2120 |
| 1     | mu (shear modulus)                         | elasticity, poro           | src/core/Simulator.cpp:2121 |
| 2     | rho_s (solid density)                      | elasticity                 | src/core/Simulator.cpp:2122 |
| 3     | phi (porosity)                             | flow, poro                 | src/core/Simulator.cpp:2138 |
| 4-6   | kx, ky, kz (permeability, m^2)             | flow, poro                 | src/core/Simulator.cpp:2139 |
| 7-9   | cw, co, cg (fluid compressibilities)       | flow                       | src/core/Simulator.cpp:2142 |
| 10-12 | mu_w, mu_o, mu_g (viscosities)             | flow, poro                 | src/core/Simulator.cpp:2145 |
| 13-15 | Swr, Sor, Sgr (residual saturations)       | flow                       | src/core/Simulator.cpp:2148 |
| 16-18 | nw, no, ng (Corey exponents)               | flow                       | src/core/Simulator.cpp:2151 |
| 19-21 | krw0, kro0, krg0 (kr maxima)               | flow                       | src/core/Simulator.cpp:2154 |
| 22    | biot_alpha                                 | poro                       | src/core/Simulator.cpp:2123 |
| 23    | 1/M (Biot storage inverse)                 | poro                       | src/core/Simulator.cpp:2161 |
| 24    | rho_f (fluid density)                      | flow, poro                 | src/core/Simulator.cpp:2157 |
| 25    | COHESIVE_CONST_MODE (0 locked, 1 slipping, 2 prescribed) | cohesive kernel | include/physics/CohesiveFaultKernel.hpp:47 |
| 26    | COHESIVE_CONST_MU_F (Coulomb friction)     | cohesive kernel           | include/physics/CohesiveFaultKernel.hpp:48 |
| 27    | COHESIVE_CONST_TENSILE_STRENGTH            | cohesive kernel           | include/physics/CohesiveFaultKernel.hpp:49 |
| 28-30 | COHESIVE_CONST_PRESCRIBED_SLIP_{X,Y,Z}     | cohesive kernel           | include/physics/CohesiveFaultKernel.hpp:50 |
| 31    | HYDROFRAC_CONST_PRESSURE                   | hydrofrac                  | src/core/Simulator.cpp:2189 |
| 32-35 | EP_CONST_{COHESION,FRICTION,DILATION,HARDENING} | elastoplasticity     | src/core/Simulator.cpp:2194 |
| 36-53 | TRACTION_CONST_BASE + 6 faces * 3 components | traction BC              | src/core/Simulator.cpp:2237 |
| 54    | VISCO_CONST_N (number of mechanisms)       | viscoelastic              | src/core/Simulator.cpp:2207 |
| 55-59 | VISCO_CONST_TAU_BASE + m, m in [0,5)       | viscoelastic              | src/core/Simulator.cpp:2209 |
| 60-64 | VISCO_CONST_DMU_BASE + m                   | viscoelastic              | src/core/Simulator.cpp:2210 |
| 65-69 | VISCO_CONST_DK_BASE + m                    | viscoelastic              | src/core/Simulator.cpp:2212 |
| 70    | COHESIVE_CONST_FRICTION_MODEL (0 constant, 1 slip-weakening) | cohesive kernel | include/physics/CohesiveFaultKernel.hpp:56 |
| 71    | COHESIVE_CONST_MU_S (static)               | cohesive kernel           | include/physics/CohesiveFaultKernel.hpp:57 |
| 72    | COHESIVE_CONST_MU_D (dynamic)              | cohesive kernel           | include/physics/CohesiveFaultKernel.hpp:58 |
| 73    | COHESIVE_CONST_DC (critical slip distance) | cohesive kernel           | include/physics/CohesiveFaultKernel.hpp:59 |
| 74    | thermal_conductivity                       | thermal                    | CLAUDE.md |
| 75    | specific_heat                              | thermal                    | CLAUDE.md |
| 76    | thermal_expansion_coeff                    | thermal                    | CLAUDE.md |
| 77    | reference_temperature                      | thermal                    | CLAUDE.md |
| 78    | thermal_field_index                        | thermal                    | CLAUDE.md |

### 1.2 Free ranges

- Slot 79 is the first guaranteed-free slot.
- Gravity, when enabled with aux callbacks, repurposes slots 0-2 for gravity
  magnitude and zeros (src/core/Simulator.cpp:2183). Rate-and-state extensions
  must not collide with that use case, so they must sit at or above slot 79.

### 1.3 Proposed extensions for rate-and-state friction

Place rate-and-state parameters in a contiguous block at slot 79 to leave a
small gap for future thermal expansion without pushing over 90:

| Slot | Proposed symbol                 | Meaning |
|------|---------------------------------|---------|
| 79   | COHESIVE_CONST_RS_LAW           | 0 aging, 1 slip, 2 regularized aging |
| 80   | COHESIVE_CONST_RS_A             | direct-effect parameter a |
| 81   | COHESIVE_CONST_RS_B             | evolution-effect parameter b |
| 82   | COHESIVE_CONST_RS_L             | state evolution length L (m) |
| 83   | COHESIVE_CONST_RS_F0            | reference friction |
| 84   | COHESIVE_CONST_RS_V0            | reference slip rate (m/s) |
| 85   | COHESIVE_CONST_RS_V_LIN         | regularization threshold (m/s) |
| 86   | COHESIVE_CONST_RS_THETA_INIT    | fallback initial state (s) |

Spatially varying (a, b, L, f0, v0) is out of scope here. When that is needed,
move these parameters to an auxiliary field keyed on fault vertices; see
section 3.

`COHESIVE_CONST_COUNT` stays at 31 so that the base cohesive kernel continues
to work without rate-and-state. Introduce `COHESIVE_CONST_RS_COUNT = 87` and
have Simulator grow `nconst_unified` to RS_COUNT only when the rate-and-state
model is active (parallel to the slip-weakening pattern at
src/core/Simulator.cpp:2250).

## 2. Auxiliary Field Layout

### 2.1 Volume cells (primary aux DM)

- `auxDM_` is created by cloning the primary DM (src/core/Simulator.cpp:2749)
  and given `NUM_AUX_FIELDS = 3` scalar fields (include/numerics/PetscFEElasticityAux.hpp:15):
  field 0 lambda, field 1 mu, field 2 rho.
- Each field is degree-0 Lagrange (piecewise constant per cell), quadrature
  matched to the primary FE (src/core/Simulator.cpp:2761).
- The aux vector is created, zero-initialised, populated by depth,
  gmsh-label, or velocity-model readers, and attached with
  `DMSetAuxiliaryVec(dm, NULL, 0, 0, auxVec_)` at src/core/Simulator.cpp:2795.
- Elastoplasticity adds 5 aux fields (EP_NUM_AUX_FIELDS = 5) that overlay the
  same aux DM machinery but are only consumed by the plasticity callbacks.
- Viscoelastic memory variables are explicitly NOT in the aux DM; they are
  updated by a TSPostStep callback against a separate storage mechanism
  (src/core/Simulator.cpp:2769 and src/numerics/PetscFEViscoelastic.cpp:197).

### 2.2 Cohesive cells

- No auxiliary field is currently attached to cohesive cells.
  `CohesiveFaultKernel::updateAuxiliaryState` is a stub
  (src/physics/CohesiveFaultKernel.cpp:101).
- Cohesive-cell friction state (theta, slip-rate history, cumulative slip) is
  held only in `FaultCohesiveDyn::dynamic_states` as a `std::vector` in host
  memory, and is not visible to PetscDS pointwise callbacks.
- Prescribed slip is passed through `constants[28..30]` rather than an aux
  field. Kinematic parameters do not yet vary per vertex.

### 2.3 What the current layout implies for extensions

- Adding a cohesive-only aux field requires either a label-restricted
  `DMAddField` (the approach reverted in docs/LAGRANGE_FIX_STATUS.md because
  PETSc 3.25 region DS does not support volume assembly via
  DMPlexTSComputeIFunctionFEM) or a separate aux DM attached by label via
  `DMSetAuxiliaryVec(dm, label, value, 0, faultAuxVec)`.
- PyLith handles this cleanly: `FaultCohesiveKin::createAuxiliaryField`
  (pylith:libsrc/pylith/faults/FaultCohesiveKin.cc:163) builds a dedicated
  auxiliary Field whose subfields (slip, slip_rate, slip_acc) are discretised
  on the fault submesh and attached via DMSetAuxiliaryVec. We should follow
  that pattern rather than forcing theta into the primary `auxDM_`.

## 3. Rate-and-State Wiring Plan

### 3.1 Constant slots

Use slots 79-86 as listed in section 1.3. Spatially uniform parameters ship
through `unified_constants` in Simulator::setupPhysics; spatially varying
parameters move to the fault aux field in section 3.2.

### 3.2 Auxiliary field for theta

Add a fault-local aux DM with one scalar subfield `theta` (and optionally
`slip_rate` so the callback does not have to numerically differentiate). This
mirrors PyLith's pattern: a separate topology (the fault submesh already
produced by FaultMeshManager) with a dedicated local Vec attached via
`DMSetAuxiliaryVec(dm, cohesive_label, 0, 0, faultAuxVec_)`. The cohesive label
is already created at src/core/Simulator.cpp:3549 (createCohesiveCellLabel),
though it is currently disabled at src/core/Simulator.cpp:3531. Re-enable it
for rate-and-state.

### 3.3 Callback signatures

Register a new boundary residual `f0_lagrange_rs_friction` and matching
Jacobian blocks via PetscDSSetBdResidual / PetscDSSetBdJacobian, parallel to
`f0_lagrange_constraint` (src/physics/CohesiveFaultKernel.cpp:160). The
signature matches the existing PetscBdPointFn prototype, and the callback
reads:

- `a[aOff[i_theta]]` from the fault aux vec,
- rate-and-state parameters from `constants[79..86]`,
- displacement jump and Lagrange traction from `u[uOff[...]]`.

Residual:
  tau_strength = a * sigma_n_eff * arcsinh( V / (2 V0) * exp((f0 + b ln(V0 theta / L)) / a) )
  f_lambda = lambda_tangential - tau_strength * slip_direction + max(-slip_n, 0) * n

For the first cut, keep the linear regularization used in the slip-weakening
branch and only switch in the full regularized form once the integration test
converges.

### 3.4 State update strategy

State evolution is stiff if integrated inside Newton. Evolve theta in
TSPostStep (after the nonlinear solve converges), not inside the residual:

1. In TSPostStep, pull the converged slip jump and Lagrange traction from the
   solution vector.
2. Compute V_{n+1} = |(u+ - u-)_{n+1} - (u+ - u-)_n| / dt per fault vertex.
3. Integrate the aging or slip law over [t_n, t_{n+1}] with V held fixed at
   V_{n+1}. Both laws have closed-form solutions for fixed V:
     aging: theta(t+dt) = L/V + (theta - L/V) exp(-V dt / L)
     slip:  theta(t+dt) = (L/V) (V theta / L)^{exp(-V dt / L)}
4. Write theta back into faultAuxVec_ so the next Newton solve sees it.

This matches the typical quasi-dynamic cycle-simulator pattern and avoids
making the Newton solve doubly nonlinear. The tradeoff is that theta lags the
displacement by one step; for rupture-scale problems (short dt) the lag is
acceptable. Document this explicitly in the test.

## 4. Kinematic Slip Time Function Wiring

### 4.1 Current behaviour

- `CohesiveFaultKernel::setPrescribedSlip` stores a Cartesian triple and
  `registerPrescribedSlipWithDS` writes it to constants[28..30]
  (src/physics/CohesiveFaultKernel.cpp:36 and :595).
- `f0_prescribed_slip` does `f[d] = (u+ - u-)[d] - constants[28+d]`, which is
  time-independent (src/physics/CohesiveFaultKernel.cpp:629).
- Simulator::setupFaultNetwork applies a strike / dip / opening rotation once
  at setup (src/core/Simulator.cpp:3474-3511) and never updates again. Example
  04 and tests/integration/test_prescribed_slip.cpp both rely on a static
  jump.

### 4.2 Backward-compatible extension

Add a slip-time-function selector slot and parameters without changing the
meaning of the existing constants:

| Slot  | Symbol                      | Meaning |
|-------|-----------------------------|---------|
| 87    | COHESIVE_CONST_SLIP_FN      | 0 step (default), 1 ramp, 2 brune, 3 liu, 4 constant_rate |
| 88    | COHESIVE_CONST_SLIP_TIME    | rupture initiation time (s) |
| 89    | COHESIVE_CONST_SLIP_TAU     | rise time or characteristic time (s) |

`f0_prescribed_slip` becomes:

  amp[d] = constants[COHESIVE_CONST_PRESCRIBED_SLIP_X + d]
  fn     = constants[COHESIVE_CONST_SLIP_FN]   // 0 if missing
  t0     = constants[COHESIVE_CONST_SLIP_TIME] // 0 if missing
  tau    = constants[COHESIVE_CONST_SLIP_TAU]  // 1 if missing
  s(t)   = slip_shape(fn, t - t0, tau)         // defaults to step in amp units
  f[d]   = (u+ - u-)[d] - s(t) * amp_hat[d] * |amp|

Because the default `fn = 0` returns `s(t) = 1 for t >= t0`, and `t0 = 0`, the
callback reduces to the current constant-jump behaviour. Example 04 and
tests/integration/test_prescribed_slip.cpp both pass parameters only for
slots 28-30, so they keep working as long as newly-allocated slots default
to zero (they do; Simulator zero-initialises unified_constants at
src/core/Simulator.cpp:2107).

Implementation notes:

- The shape functions already exist in include/domain/geomechanics/PyLithFault.hpp
  (SlipTimeFnRamp, SlipTimeFnBrune, SlipTimeFnLiu, etc.). Keep those as
  host-side configurators that write parameters into the constants array.
  The callback stays a small pure-C function for PetscDS compatibility.
- When spatially varying parameters are needed later (per-vertex t0 or rise
  time, Liu normalisation tables), promote them to the fault aux field from
  section 3.2 rather than padding the constants array.

## 5. Poroelastic Fault Coupling

### 5.1 Field layout options

Option A. Second Lagrange field for pressure (scalar).

- Add a 1-component Lagrange multiplier field `lagrange_p` alongside the
  existing 3-component `lagrange` field.
- Register a second pair of BdResidual and Jacobian callbacks for the
  pressure Lagrange equation.
- Matches PyLith separation of concerns. Easier to audit and to toggle on /
  off per test.

Option B. Add a 4th component to the existing Lagrange field.

- Extend fe_lagrange from 3 to 4 components (tx, ty, tz, p_lag) at
  src/core/Simulator.cpp:1684.
- Fewer field registrations, but blurs the meaning of the field and forces
  every existing BdResidual to loop over the wrong component count.
- Rejected as the default. Only revisit if PETSc region DS limitations in 3.25
  make option A unworkable.

### 5.2 Constraint equations

- Continuous pressure (impermeable fault):
    f_p_lag = p+ - p-
  Analogous to the locked displacement constraint and reuses the same
  Jacobian pattern.
- Resistive fault with leak-off:
    q = T_hyd (p+ - p-)          (Darcy across a thin fault zone)
    f_p_lag = q - q_prescribed
  The pressure Lagrange multiplier represents the mass flux across the
  interface.

### 5.3 Leak-off on the bulk pressure equation

PyLith and the existing hydrofrac code both drive the bulk pressure equation
from a boundary residual on cohesive cells: see
src/numerics/PetscFEHydrofrac::f0_lagrange_pressure_balance, invoked from
setupPhysics at src/core/Simulator.cpp:2537 and 2623. Reuse this machinery:

- Register `PetscDSSetBdResidual(prob, pressure_field, f0_p_leakoff, nullptr)`
  that injects `q` as a source in the bulk pressure equation.
- Register `PetscDSSetBdResidual(prob, pressure_lagrange_field, f0_p_constraint, nullptr)`
  to enforce either continuity or the resistive constraint.

### 5.4 Which PETSc DS machinery assembles the pressure-side residual

The bulk pressure residual on cohesive cells is assembled by
DMPlexTSComputeIFunctionFEM through PetscDSSetBdResidual on the pressure
field, the same path that src/core/Simulator.cpp:2547 uses for displacement.
The manual Jacobian path (addCohesivePenaltyToJacobian at
src/core/Simulator.cpp:4182) must grow a pressure-Lagrange block because
PETSc 3.25 BdJacobian is not functional (CLAUDE.md rule section). Plan to
extend the loop at src/core/Simulator.cpp:4292 to also write
`d(f_p_lag)/d(p)` and `d(f_p)/d(p_lag)` entries.

Keep this work behind a `config.enable_fault_poroelastic` flag so that every
existing fault test continues to use the 3-component Lagrange field.

## 6. Known Failure Diagnosis: PrescribedSlipQuasiStatic

### 6.1 CLAUDE.md hypothesis

CLAUDE.md states: "The penalty-scaled Lagrange diagonal (penalty*coeff) in
addCohesivePenaltyToJacobian slows lambda convergence for prescribed slip
mode."

### 6.2 What the code actually does

In src/core/Simulator.cpp:4190:

    const bool is_slipping = (fault_mode_ == "slipping");
    const PetscReal penalty_scale = is_slipping ? 0.0 : 1.0;

Prescribed slip therefore takes the full penalty branch. The penalty is
`10.0 * youngs_modulus / h_char_jac` (src/core/Simulator.cpp:4342), which for
E = 10 GPa and h ~ 0.25 m is on the order of 4e11. That penalty is added to
`d(f_lambda)/d(lambda)` at every cohesive vertex.

For the prescribed-slip constraint `f_lambda = (u+ - u-) - delta`, the
analytical `d(f_lambda)/d(lambda)` is zero. The manual penalty therefore
injects a large spurious compliance into an equation that should be a pure
kinematic constraint. Each Newton step resolves
`penalty * dLambda ~= residual`, which pushes lambda toward the reaction very
slowly when `penalty` is enormous. Residual in u is roughly constant each
iteration, so Newton stalls before the 1 mm jump is reached.

### 6.3 Minimal test of the hypothesis

Change src/core/Simulator.cpp:4190 to:

    const bool skip_penalty = (fault_mode_ == "slipping" ||
                               fault_mode_ == "prescribed_slip");
    const PetscReal penalty_scale = skip_penalty ? 0.0 : 1.0;

Re-run Integration.DynamicRuptureSolve.PrescribedSlipQuasiStatic in Docker.
Pass criterion is already encoded in the test
(tests/integration/test_dynamic_rupture_solve.cpp:433): tangential jump in
(5e-4, 2e-3) m.

Confounders to track before declaring the hypothesis confirmed:

- The locked fault test relies on the penalty to keep lambda bounded on
  coarse meshes. Run Integration.DynamicRuptureSolve.LockedFaultQuasiStatic
  in the same build to confirm the change does not regress the locked case.
- If Newton still stalls with penalty_scale = 0, the residual of
  `f0_prescribed_slip` may be starved by the weak volume regularisation
  (`f = epsilon * lambda` at src/core/Simulator.cpp:1739, epsilon = 1e-4).
  On a 4^3 mesh the cohesive BdResidual area-weighting can be smaller than
  the volume regularisation times domain size. If that turns out to be the
  next bottleneck, gate f0_weak_lagrange off on cohesive cells rather than
  raising epsilon.
- Optional second check: register a direct L1 norm of the Lagrange equation
  at the end of SNES (`-snes_monitor` plus a small custom monitor) and watch
  whether it decays by orders of magnitude once penalty_scale is zero.

## 7. Rules Compliance

Every phase of the plan maps back to specific CLAUDE.md rules. The phases are:

- Phase A. Diagnose PrescribedSlipQuasiStatic (section 6).
- Phase B. Kinematic slip time function extension (section 4).
- Phase C. Poroelastic fault coupling (section 5).
- Phase D. Rate-and-state friction (section 3).

### Rule 1. Build and test in Docker. Always.

Every phase builds and runs its tests exclusively through
`docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local`. No host
installs. Phase D touches TSPostStep state update; the one-step integration
test and the full TSSolve test both run in Docker.

### Rule 2. Check PETSc 3.25.0 API before use.

Phase A touches only existing addCohesivePenaltyToJacobian code that already
uses verified PETSc calls. Phase B adds no new PETSc calls. Phase C will
introduce new PetscDSSetBdResidual, PetscDSSetBdJacobian, and
DMSetAuxiliaryVec calls with a label; all three have been verified against
grep over /opt/petsc-3.25.0/include/ in prior sessions, but each phase will
re-grep before committing. Phase D uses PetscDSSetBdResidual with the
existing pattern.

### Rule 3. All existing tests must continue to pass after every change.

Phase A runs the full ctest suite in Docker. Specifically, the locked-fault,
slipping-fault, slip-weakening, TPV5, and prescribed-slip tests must stay
green. Phase B passes because the new constants default to zero, restoring
step behaviour. Phase C and D are behind feature flags
(`enable_fault_poroelastic`, `rs_friction_enabled`) so existing tests never
see the new code paths.

### Rule 4. NEVER change the DS/BC ordering in setupFields().

No phase touches the DMCreateDS / setupBoundaryConditions /
DMSetLocalSection(nullptr) / DMSetUp / DMGetDS sequence at
src/core/Simulator.cpp:1703-1713. Phase C adds the pressure Lagrange field
before DMCreateDS, preserving the existing ordering.

### Rule 5. Do NOT modify callback math in PetscFEElasticity.cpp,
PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp.

All new kernels go in src/physics/CohesiveFaultKernel.cpp or a new
src/physics/RateStateFrictionKernel.cpp. Existing volume kernels are
read-only.

### Rule 6. Do NOT modify FaultMeshManager::splitMeshAlongFault or
CohesiveFaultKernel::registerWithDS.

Phase D uses a new `registerRateStateWithDS` entry point rather than editing
the existing one. Phase C does the same for the pressure-side coupling.
`splitMeshAlongFault` is not touched.

### Rule 7. Do not reference archived dead code.

docs/design/fault_gap_analysis.md cites only live sources under src/,
include/, tests/, or pylith:libsrc. archive/ is not referenced.

### Rule 8. No Python in the Simulator. Python is ONLY for post-processing.

All phases stay in C++17. Any theta trajectory plotting is a Python
post-processing script under scripts/, kept out of the simulator binary.

### Rule 9, 10, 11. Config conventions (ignore aspirational, .config
extension, executable is `fsrm`).

New config keys go under [FAULT] and [FAULT_POROELASTIC] with a .config
example under config/examples/ and a matching runnable script under
examples/NN_rate_state_friction/run.sh invoking `./fsrm -c ...`.

### Rule 12. No em dashes or contractions in code comments or documentation.

This document uses only hyphens and written-out forms. New kernel comments
follow the same convention.

### Rule 13. Verification tests must have quantitative pass/fail criteria.

- Phase A pass criterion: tests/integration/test_dynamic_rupture_solve.cpp:433
  continues to require jump in (5e-4, 2e-3) m after the penalty fix.
- Phase B pass criterion: a new ramped-slip test requires that at t = 0.5 s of
  a rise time of 1 s, the jump is within 5 percent of 0.5 * final_slip.
- Phase C pass criterion: a Terzaghi-like fault-transverse pressure step
  decays to within 2 percent of the analytical diffusive profile.
- Phase D pass criterion: a rate-and-state spring-slider reproduces the
  steady-state friction `f = f0 + (a - b) ln(V / V0)` within 1 percent over
  three decades of V.

### Rule 14. Update CLAUDE.md and README.md after every session.

Each phase ends with a CLAUDE.md patch that records: which slots were
allocated, which tests were added, which tests were updated, and which
feature flags were introduced. README.md gets a one-line entry in the
feature inventory table with the matching test name.

### Rule 15. GTEST_SKIP is ONLY for hardware-dependent tests or genuine
crash bugs.

Any new test either passes with a quantitative criterion or fails honestly.
No new GTEST_SKIP calls. If a rate-and-state test cannot reach steady state
on the CI-mesh budget, reduce the parameter sweep rather than skipping.
