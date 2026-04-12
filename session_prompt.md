# FSRM Session Prompt: FEM-Coupled Hydraulic Fracturing

## Branch Setup

```bash
git checkout main
git pull origin main
git checkout -b feature/hydrofrac-fem
```

All work on this branch. PR to main when all tests pass.

## READ FIRST

Read `CLAUDE.md` completely. All rules binding. 82 tests currently registered. Do not break any.

## Current State After Review

Main has 82 registered tests, zero GTEST_SKIP calls. Key capabilities verified: elastostatics, elastodynamics (TSALPHA2), poroelasticity (Biot), Mueller-Murphy source, explosion seismograms, absorbing BCs, cohesive fault mesh splitting, elastoplasticity (PetscFEElastoplasticity with Drucker-Prager return mapping), prescribed slip, gravity/lithostatic stress.

### Remaining Architecture Issues

1. **BdResidual overwrite:** PetscDSSetBdResidual stores one callback per field. Cohesive fault and absorbing BC both register on displacement field. Second overwrites first. Current guard: SETERRQ rejects enable_faults + absorbing_bc_enabled. Proper fix (DMSetRegionDS) not done.

2. **Hydraulic fracturing is analytical, not FEM-coupled:** PKN/KGD/P3D models in FractureModel.cpp compute width and length analytically from net pressure. They are called from MonitorFunction using an estimated wellbore pressure, not the FEM solution. contributeToResidual() and contributeToJacobian() are empty stubs. MultiStageFracturing.cpp has detailed engineering data structures (perforations, pump schedules, stages, erosion models) but no FEM integration.

3. **Dynamic rupture SNES divergence:** test_dynamic_rupture_basic.cpp documents that SNES may diverge at first timestep with cohesive cells enabled.

## What "Beating ResFrac" Means

ResFrac uses planar fracture geometry on structured grids with DDM (Displacement Discontinuity Method) for stress interactions. FSRM's advantage is full 3D unstructured FEM with cohesive zone mechanics. The architecture below exploits this advantage.

ResFrac capabilities to match or exceed:
1. Pressurized fracture with aperture from geomechanics (not assumed)
2. Lubrication flow in fracture (Poiseuille)
3. Leak-off into formation (Carter model coupled to pore pressure)
4. Propagation criterion from FEM stress field (not analytical SIF)
5. Multi-cluster stress shadowing through the FEM stress field (better than DDM)
6. Proppant transport (settling + advection)
7. Induced seismicity monitoring (FSRM has Mueller-Murphy + seismometers -- ResFrac does not)

## Architecture: FEM-Coupled Hydraulic Fracturing

### Overview

The fracture is represented as a set of cohesive cells in the DMPlex mesh. Fluid flow in the fracture is governed by the lubrication equation, coupled to the poroelastic response of the surrounding rock. Fracture propagation occurs when cohesive cells at the tip transition from locked to open.

Three coupled fields:
- **Displacement** u (vector, dim components) -- rock deformation
- **Pore pressure** p (scalar) -- formation fluid pressure
- **Fracture pressure** p_f (scalar, on fracture faces only) -- fluid pressure inside the fracture

The fracture aperture w is computed from the displacement discontinuity: w = (u+ - u-) dot n, where n is the fracture normal.

### Field Layout

For a poroelastic simulation with a pressurized fracture:
- Field 0: pressure p (1 component, on all cells)
- Field 1: displacement u (3 components, on all cells)
- Field 2: Lagrange multiplier lambda (3 components, on cohesive cells only)

The Lagrange multiplier enforces the traction condition on the fracture face. For a pressurized fracture: lambda = -p_f * n (the Lagrange multiplier IS the fluid pressure traction).

### PetscDS Callbacks

**Volume cells (existing, verified):**
- Residual: PetscFEPoroelasticity f0/f1 for pressure and displacement
- Jacobian: PetscFEPoroelasticity g0/g1/g2/g3 (4 Biot coupling blocks)

**Cohesive cells (NEW):**
Create `src/numerics/PetscFEHydrofrac.cpp` and `include/numerics/PetscFEHydrofrac.hpp`:

```
f0_fracture_pressure:
  Lubrication equation for fracture pressure:
  R_pf = dw/dt + div_s(Q) - q_inj + q_leak = 0
  where:
    w = aperture = (u+ - u-) dot n
    Q = -w^3/(12*mu_f) * grad_s(p_f)  (Poiseuille flow)
    q_inj = injection source term
    q_leak = C_L / sqrt(t - t_open) * (p_f - p_formation)  (Carter leak-off)

f0_fracture_displacement:
  Traction on fracture faces from fluid pressure:
  R_u = lambda + p_f * n = 0  (traction balance)

f0_fracture_lagrange:
  Constraint equation:
  R_lambda = (u+ - u-) - w_min * n = 0  (minimum aperture constraint when open)
  OR
  R_lambda = lambda - t_cohesive(delta)  (cohesive traction-separation law at tips)
```

**Propagation criterion:**
At cohesive cells near the fracture tip, the traction from the FEM solution is compared against the cohesive strength. When the normal traction exceeds tensile strength T_s, the cell transitions from "locked" to "open" mode. This is tracked via a per-cell state flag in the auxiliary fields.

### Coupling Strategy

The system is monolithically coupled: all fields (p, u, lambda, p_f) are solved simultaneously by PETSc's SNES. The Jacobian includes all cross-coupling blocks. This gives quadratic convergence and handles the strong coupling between fracture aperture and fluid pressure naturally.

For the Jacobian:
- d(R_pf)/d(p_f): lubrication stiffness (diffusion in fracture plane)
- d(R_pf)/d(u): aperture sensitivity (dw/du couples fracture flow to rock deformation)
- d(R_u)/d(p_f): pressure loading on fracture faces
- d(R_u)/d(u): elastic stiffness (existing)
- d(R_lambda)/d(u): constraint sensitivity
- d(R_p)/d(u): Biot coupling (existing)
- d(R_u)/d(p): Biot coupling (existing)

## Development Plan

### Phase 1: Pressurized Fracture (Static, No Flow)

Simplest possible coupled problem: a pre-existing fracture in a 3D elastic domain, pressurized with a uniform known pressure. No fluid flow, no propagation. Verify that the aperture from the FEM solution matches the analytical Sneddon solution.

**Sneddon solution:** For a penny-shaped crack of radius a under internal pressure P in an infinite elastic body:
  w_max = 8*(1-nu^2)*P*a / (pi*E)

**Implementation:**
1. Create a 3D box mesh with a planar fracture at center (reuse FaultMeshManager for mesh splitting).
2. Set up elastostatics with cohesive cells in prescribed-slip mode, BUT instead of prescribed displacement jump, prescribe a normal traction p_f * n on the fracture faces.
3. This requires modifying the cohesive callback: instead of f0_lambda = (u+ - u-) - delta, use f0_lambda = lambda + p_f * n (traction balance). For this phase, p_f is a known constant stored in the constants array.
4. Solve the elastostatic problem.
5. Extract the aperture w = (u+ - u-) dot n at the fracture center.
6. Compare against Sneddon.

**Test:** `tests/physics_validation/test_pressurized_fracture.cpp`
- Sneddon aperture within 20% on coarse mesh (4x4x4 to 8x8x8)
- Aperture is positive (fracture opens)
- Aperture decreases toward fracture tips
- Stress is compressive (negative) away from fracture in the normal direction

Register as `Physics.PressurizedFracture`.

**Files to create:**
- `src/numerics/PetscFEHydrofrac.cpp`
- `include/numerics/PetscFEHydrofrac.hpp`
- `tests/physics_validation/test_pressurized_fracture.cpp`

**Files to modify:**
- `src/core/Simulator.cpp` (add hydrofrac physics registration path in setupPhysics)
- `tests/CMakeLists.txt`

---

### Phase 2: Lubrication Flow in Fracture

Add fluid flow in the fracture using the Reynolds lubrication equation.

**Physics:**
  dp_f/dt * w * c_f + dw/dt + div_s(-w^3/(12*mu_f) * grad_s(p_f)) = q_inj - q_leak

This is a nonlinear diffusion equation on the fracture surface, coupled to the rock deformation through the aperture w.

**Implementation:**
1. Add a fracture pressure field (scalar, defined only on cohesive cell faces). In PETSc DMPlex, this can be a field on the hybrid cells.
2. Implement f0_fracture_pressure (lubrication equation) and f1_fracture_pressure (Poiseuille flux: w^3/(12*mu) * grad(p_f)).
3. The aperture w is computed from u+ - u- (available in the closure of cohesive cells).
4. Register these callbacks via PetscDSSetResidual on the cohesive cell DS.
5. Injection is a source term at one cohesive cell (the wellbore intersection).

**Test:** `tests/integration/test_fracture_flow.cpp`
- Constant aperture, injection at center, no-flow boundaries at tips.
- Verify pressure diffuses outward from injection point.
- Verify steady-state pressure profile matches analytical Poiseuille solution for a slot.
- Verify mass conservation (integral of q_inj * dt = integral of w * p_f * c_f * dA + leaked volume).

Register as `Integration.FractureFlow`.

---

### Phase 3: Coupled Fracture Flow + Deformation

Couple the lubrication equation to the elastic/poroelastic deformation. The aperture is now computed from the FEM solution, and the pressure in the fracture loads the fracture faces.

**Implementation:**
1. Solve the monolithic system: [p, u, lambda, p_f] simultaneously.
2. The Jacobian must include cross-coupling terms: d(R_pf)/d(u) (aperture sensitivity) and d(R_u)/d(p_f) (pressure loading).
3. For the d(R_pf)/d(u) term: dw/du = n (the normal vector). This enters the lubrication equation through the dw/dt term and the w^3 nonlinearity.
4. For the d(R_u)/d(p_f) term: d(p_f * n)/d(p_f) = n.

**Test:** `tests/integration/test_coupled_hydrofrac.cpp`
- Inject fluid at constant rate into a pre-existing fracture in a poroelastic medium.
- Verify fracture aperture increases with injection.
- Verify pressure in the fracture is higher than pore pressure in the formation.
- Verify leak-off occurs (formation pore pressure increases near the fracture).
- Compare fracture width vs. time against PKN analytical solution for a viscosity-dominated regime: w ~ (mu*Q*t/E')^(1/4).

Register as `Integration.CoupledHydrofrac`.

---

### Phase 4: Fracture Propagation via Cohesive Zone Failure

Add propagation: cohesive cells at the fracture tip transition from locked to open when the normal traction exceeds the tensile strength.

**Implementation:**
1. Add a per-cell state variable: `fracture_state` (0=locked, 1=open). Store in auxiliary fields.
2. In the cohesive callback for locked cells: enforce zero displacement jump (locked constraint).
3. In the cohesive callback for open cells: enforce traction balance with fracture pressure.
4. After each converged time step (post-solve callback): check the traction on locked cells adjacent to open cells. If normal traction > T_s, transition to open.
5. When a cell transitions to open, the fracture pressure field is extended to include the new cell, with initial pressure from the formation pore pressure.
6. Re-solve the coupled system with the updated fracture geometry.

**Propagation criterion:** For LEFM consistency, the cohesive strength T_s should be related to the fracture toughness K_Ic and the element size h:
  T_s = K_Ic / sqrt(pi * h / 2)

This ensures the cohesive zone model converges to the LEFM solution as the mesh is refined.

**Test:** `tests/integration/test_fracture_propagation.cpp`
- Inject into a short initial fracture (2-3 cohesive cells).
- Verify fracture extends (more cells transition to open).
- Verify fracture length vs. time follows KGD or PKN scaling at early/late times.
- Verify pressure at wellbore decreases as fracture grows (larger compliant volume).

Register as `Integration.FracturePropagation`.

---

### Phase 5: Multi-Cluster Stress Shadowing

The key advantage over ResFrac: stress shadowing computed from the full 3D FEM stress field, not from DDM.

**Implementation:**
1. Create a mesh with multiple parallel fractures (N clusters) at specified spacing.
2. Each fracture is a set of cohesive cells, independently pressurized.
3. Inject into all fractures simultaneously (or sequentially, configurable).
4. The FEM stress field naturally captures the stress shadow: opening one fracture increases the normal stress on adjacent fractures, requiring higher pressure to open them.

**Test:** `tests/integration/test_stress_shadowing.cpp`
- 3 parallel fractures at 10m spacing in a 100m x 100m x 50m domain.
- Inject into the center fracture only.
- Verify: the outer fractures experience increased compressive stress (stress shadow).
- Verify: if all three are injected, the center one has the widest aperture (least stress shadow).
- Compare stress shadow magnitude against DDM solution for 3 parallel cracks (analytical).

Register as `Integration.StressShadowing`.

---

### Phase 6: Proppant Transport

**Implementation:**
1. Add a proppant concentration field c on cohesive cells (scalar, advection-diffusion in fracture plane).
2. Transport equation: d(w*c)/dt + div_s(c*Q) - div_s(D*w*grad_s(c)) + v_settle * dc/dz = 0
   where v_settle = Stokes settling velocity = d_p^2*(rho_p - rho_f)*g / (18*mu_f)
3. Proppant bridging: when c > c_max (pack concentration), set w_min = d_p (proppant holds fracture open).
4. After shut-in, fracture closes to w_min where proppant is present; closes fully where proppant is absent.

**Test:** `tests/integration/test_proppant_transport.cpp`
- Inject proppant-laden fluid into a fracture with known aperture profile.
- Verify proppant concentration decreases with distance from wellbore.
- Verify proppant settles (higher concentration at bottom of fracture).
- Verify proppant prevents full closure after shut-in.

Register as `Integration.ProppantTransport`.

---

### Phase 7: Induced Seismicity from Hydraulic Fracturing

The killer feature that ResFrac cannot match: combine the hydraulic fracturing simulation with the explosion monitoring pipeline to generate synthetic microseismic catalogs.

**Implementation:**
1. When a cohesive cell transitions from locked to open (fracture propagation event), compute the moment tensor from the displacement discontinuity and the area of the cell.
2. Feed the moment tensor to the SeismometerNetwork to generate synthetic seismograms.
3. Compute event magnitude from M0 = mu * delta_u * A (shear modulus * slip * area).
4. Record the event in a seismic catalog (time, location, magnitude, moment tensor).

**Test:** `tests/integration/test_induced_seismicity.cpp`
- Hydraulic fracturing simulation with fracture propagation.
- Verify at least one microseismic event is recorded when a cohesive cell opens.
- Verify event magnitude is in the range -3 < Mw < 1 (typical microseismic).
- Verify event location matches the fracture tip position.

Register as `Integration.InducedSeismicity`.

---

### Phase 8: Carter Leak-Off Coupling

**Implementation:**
1. On open fracture faces, add a leak-off source term to the formation pressure equation: q_leak = C_L * A_face / sqrt(t - t_open) * sign(p_f - p_formation)
2. The corresponding sink term appears in the fracture lubrication equation.
3. This couples the fracture and formation pressure fields bidirectionally.

**Test:** `tests/integration/test_leakoff_coupling.cpp`
- Inject into a fracture in a permeable formation.
- Verify formation pore pressure increases near the fracture.
- Verify fracture pressure is higher than formation pressure.
- Verify leak-off rate decreases with sqrt(t) (Carter scaling).
- Verify total leaked volume = integral of leak-off rate.

Register as `Integration.LeakoffCoupling`.

---

### Phase 9: Production Forecasting (Post-Frac)

After fracturing and proppant placement, simulate production by reversing the pressure boundary condition: formation pressure at the wellbore is lowered below reservoir pressure.

**Implementation:**
1. After shut-in, the fracture pressure equilibrates with formation pressure.
2. Set wellbore pressure below reservoir pressure (drawdown).
3. Fluid flows from formation through propped fracture to wellbore.
4. Production rate = sum of fracture flux at the wellbore intersection.
5. Decline curve analysis: fit production vs. time to Arps hyperbolic decline.

**Test:** `tests/integration/test_production_forecast.cpp`
- Complete hydraulic fracturing sequence: inject, propagate, place proppant, shut-in, produce.
- Verify initial production rate is nonzero.
- Verify production declines with time.
- Verify cumulative production approaches EUR (estimated ultimate recovery).

Register as `Integration.ProductionForecast`.

---

### Phase 10: BdResidual Architecture Fix

Before all the fracture phases above can work with absorbing BCs (needed for seismogram output during induced seismicity monitoring), fix the BdResidual overwrite.

**Implementation:**
Use PETSc's DMSetRegionDS to assign separate DS instances to:
- External boundary faces (absorbing BC callbacks)
- Cohesive interface cells (fracture/fault callbacks)

Check PETSc 3.22.2 API:
```bash
grep -rn "DMSetRegionDS\|DMGetRegionDS\|DMSetRegionNumDS" /opt/petsc-3.22.2/include/
```

If DMSetRegionDS is not available, use DMGetCellDS per-cell to assign different DS instances based on cell type (regular vs. hybrid).

**Test:** `tests/integration/test_fault_absorbing_coexist.cpp`
- Mesh with cohesive fault AND absorbing BCs.
- Verify cohesive traction is applied on fault faces.
- Verify absorbing BC traction is applied on external faces.
- Verify no callback overwrite (both callbacks fire).

Register as `Integration.FaultAbsorbingCoexist`.

---

### Phase 11: Integration Tests and Cleanup

1. Run full test suite: all 82+ existing tests plus all new tests.
2. Update docs/TEST_RESULTS.md with complete test matrix.
3. Update CLAUDE.md.
4. Update README.md: add "Verified Capabilities" entry for hydraulic fracturing.
5. Clean up dead code: remove the standalone PKN/KGD/P3D propagation from MonitorFunction (it is replaced by the FEM-coupled propagation).

---

### Phase 12: PR to Main

```bash
git push origin feature/hydrofrac-fem
gh pr create \
  --base main \
  --head feature/hydrofrac-fem \
  --title "FEM-coupled hydraulic fracturing with induced seismicity" \
  --body "## FEM-Coupled Hydraulic Fracturing

Replaces the standalone PKN/KGD/P3D analytical models with a fully FEM-coupled
cohesive zone fracture propagation system.

### Capabilities (all tested)
- Pressurized fracture with Sneddon verification
- Lubrication flow in fracture (Reynolds equation)
- Monolithic coupling: fracture pressure + rock deformation + pore pressure
- Fracture propagation via cohesive zone failure
- Multi-cluster stress shadowing from 3D FEM stress field
- Proppant transport with settling and bridging
- Carter leak-off coupling to formation
- Production forecasting with decline curves
- Induced seismicity monitoring (moment tensor + seismograms)

### Advantages over ResFrac
- Full 3D unstructured FEM (vs structured grid)
- Stress shadowing from FEM (vs DDM approximation)
- Elastoplastic damage at fracture tips
- Integrated seismic monitoring (Mueller-Murphy + SAC output)
- Open source (MIT license)

### Test Results
[N] tests pass. Zero failures. Zero skips."
```

Wait for CI green. Fix failures. Merge.

## Execution Order

Phases 1-4 are the critical path. Without them, nothing else works.

1. Phase 1 (pressurized fracture) -- foundation
2. Phase 2 (lubrication flow) -- fracture physics
3. Phase 3 (coupled flow + deformation) -- the real coupling
4. Phase 4 (propagation) -- the breakthrough feature
5. Phase 10 (BdResidual fix) -- unblocks Phase 7
6. Phase 5 (stress shadowing) -- key differentiator vs. ResFrac
7. Phase 7 (induced seismicity) -- killer feature
8. Phase 6 (proppant) -- completions engineering
9. Phase 8 (leak-off) -- formation coupling
10. Phase 9 (production) -- commercial application
11. Phase 11-12 (cleanup and PR)

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.22.2 API signatures before calling any function.
3. All existing 82 tests must continue to pass.
4. New PetscDS callbacks go in NEW files (PetscFEHydrofrac.cpp).
5. Every test has quantitative pass/fail with numerical tolerances.
6. Commit after each passing phase.
7. The cohesive cell infrastructure (FaultMeshManager, CohesiveFaultKernel) is the foundation for all fracture work. Understand it before writing new code. Read the PyLith cohesive cell workflow in CLAUDE.md.
8. The lubrication equation (Phase 2) is the hardest part. The w^3 nonlinearity in the Poiseuille flux makes the system stiff. Use continuation or adaptive time stepping.
9. Do not modify FaultMeshManager::splitMeshAlongFault. Create new cohesive cells by extending the fault label before calling split, not by modifying the split algorithm.
10. PetscFEHydrofrac callbacks must handle both 2D (surface integral on fracture face) and 3D (volume integral on cohesive cell) integration contexts. Check PETSc 3.22 documentation for how DMPlex handles hybrid cell integration.