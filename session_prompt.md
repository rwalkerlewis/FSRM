# FSRM Session Prompt: Phase 2 -- Dynamic Rupture, Plasticity, and Architecture Fixes

## READ FIRST

Read `CLAUDE.md` in the repository root completely before starting. All rules there are binding. 68 tests currently pass. Do not break any of them.

## Current State

The local_fix merge brought the explosion monitoring pipeline to a verified state. 68 tests pass covering elastostatics, elastodynamics, poroelasticity, Mueller-Murphy source, explosion seismograms, atmospheric/near-field explosion physics, absorbing BCs, gravity, fault infrastructure, and plasticity yield evaluation.

## Remaining Gaps

Ranked by impact:

1. **Cohesive fault forces are overwritten by absorbing BCs** -- Critical architecture bug
2. **Dynamic rupture never runs end-to-end** -- SCEC TPV5 tests infrastructure only
3. **Plasticity return mapping broken** -- Does not produce nonzero plastic strain
4. **Plasticity not FEM-coupled** -- Not wired to PetscDS callbacks
5. **20 GTEST_SKIP calls remain** -- Conditional skips in fault/MMS tests
6. **PhysicsKernel strong-form placeholder** -- MMS tests skip because standalone kernel residual methods are wrong (PetscDS callbacks are correct)

## Architecture Bug: BdResidual Overwrite

**This is the most important finding in this review. Read it carefully.**

In setupPhysics(), the cohesive fault kernel and absorbing BC both register callbacks via PetscDSSetBdResidual on the SAME displacement field index.

```
Line 1868: cohesive_kernel_->registerWithDS(prob, disp_field, lagrange_field)
  -> PetscDSSetBdResidual(prob, displacement_field, f0_displacement_cohesive, ...)

Line 1893: PetscDSSetBdResidual(prob, displacement_field_idx, AbsorbingBC::f0_absorbing, ...)
```

PetscDSSetBdResidual stores ONE callback per field. The second call OVERWRITES the first. When both faults and absorbing BCs are enabled, the cohesive displacement traction is silently replaced by the absorbing BC traction. The Lagrange multiplier BdResidual survives (absorbing BCs do not touch it), but the system becomes inconsistent: lambda enforces a constraint but the corresponding traction on displacement is missing.

**When absorbing BCs are disabled**, the cohesive BdResidual for displacement should be intact. Dynamic rupture might work without absorbing BCs. This should be tested first.

**The proper fix** requires one of:
A. Use PETSc's label-based DS system (DMSetRegionDS / DMGetRegionDS) to assign different DS instances to external boundary faces vs cohesive interface cells. This is the correct PETSc 3.22 approach.
B. Write a combined callback that dispatches to either absorbing or cohesive based on which label the current face belongs to.
C. Never enable both faults and absorbing BCs simultaneously (document as limitation).

Option A is correct but complex. Option C is pragmatic. Start with C to unblock dynamic rupture testing, then implement A.

## Development Plan

### Phase 1: Dynamic Rupture Without Absorbing BCs

**Goal:** Get a cohesive fault simulation to run through TSSolve and produce nonzero slip.

1. Create `tests/integration/test_dynamic_rupture_basic.cpp`:
   - 3D box, 4x4x4 simplex mesh (enable_faults=true forces simplex)
   - Vertical fault at center via createPlanarFaultLabel
   - FaultMeshManager splits mesh (32 cohesive cells verified by existing test)
   - Elastodynamics enabled, TSALPHA2
   - NO absorbing BCs (sides=fixed or roller on bottom, free on top and sides)
   - Cohesive kernel in slipping mode with slip-weakening friction
   - Initial stress state that exceeds yield on the fault (tau_0 > mu_s * sigma_n)
   - Run for 10 timesteps
   - Assert: simulation completes without SNES divergence
   - Assert: displacement is nonzero
   - Assert: displacement is discontinuous across the fault (u_plus != u_minus on at least one cohesive face)

2. If TSSolve diverges:
   - Check that the Lagrange multiplier field is created with correct number of components (dim=3)
   - Check that PetscDSSetBdResidual for the Lagrange field is not overwritten
   - Check that DMPlexTSComputeIFunctionFEM actually loops over cohesive cells (add a print inside f0_lagrange_cohesive to verify it gets called)
   - Try with locked fault first (mode=0, no slip) to verify the constraint alone is assembled
   - Then switch to slipping mode

3. If TSSolve converges but slip is zero:
   - The fault may be locked because the initial stress is below yield
   - Apply a large enough traction to exceed mu_s * sigma_n
   - Check that f0_lagrange_cohesive is computing the correct Coulomb criterion

Register as `Integration.DynamicRuptureBasic` with label `integration`.

---

### Phase 2: Fix Plasticity Return Mapping

**Problem:** DruckerPragerModel::returnMapping does not produce nonzero plastic strain. The test_drucker_prager.cpp documents this.

**Diagnosis steps:**

1. Read DruckerPragerModel::returnMapping in src/domain/geomechanics/PlasticityModel.cpp carefully.
2. Set up a simple test case: uniaxial compression sigma_zz = -200 MPa, all other components zero. With cohesion = 10 MPa, friction angle = 30 deg, the Drucker-Prager yield surface is:
   F = sqrt(J2) + alpha*I1 - k = 0
   where alpha = 2*sin(phi)/(sqrt(3)*(3-sin(phi))), k = 6*c*cos(phi)/(sqrt(3)*(3-sin(phi)))
   For phi=30, c=10 MPa: alpha = 0.236, k = 11.55 MPa
   I1 = -200 MPa, J2 = (1/3)*200^2 = 13333 MPa^2, sqrt(J2) = 115.5 MPa
   F = 115.5 + 0.236*(-200) - 11.55 = 115.5 - 47.2 - 11.55 = 56.7 MPa > 0 (outside yield surface)
3. Call returnMapping with this stress state. The returned stress should be ON the yield surface (F=0 after return) and accumulated_plastic_strain should be > 0.
4. If plastic strain is still zero, trace through the return mapping algorithm step by step. Common bugs:
   - Wrong sign convention (tension positive vs compression positive)
   - Division by zero in the return direction
   - Convergence tolerance too tight (Newton iterations hit max without converging)
   - Plastic strain increment computed but not accumulated

Fix the return mapping. Verify with the test.

**Files to modify:**
- `src/domain/geomechanics/PlasticityModel.cpp` (returnMapping methods)
- `tests/physics_validation/test_drucker_prager.cpp` (update assertions to require nonzero plastic strain)

---

### Phase 3: Couple Plasticity to PetscDS Callbacks

**Problem:** PlasticityModel is standalone code. It is not called during FEM assembly.

**Architecture:** The elasticity f1 callback computes sigma = C : eps (linear elastic). To add plasticity, the callback must:
1. Compute trial stress: sigma_trial = C : eps_total
2. Check yield function: F(sigma_trial)
3. If F > 0: apply return mapping to get sigma_corrected and plastic strain increment
4. Return sigma_corrected as f1

**Implementation:**

Create `src/numerics/PetscFEElastoplasticity.cpp` and `include/numerics/PetscFEElastoplasticity.hpp`:

```cpp
namespace PetscFEElastoplasticity {
    // f1_elastoplastic_aux: same as f1_elastostatics_aux but with return mapping
    void f1_elastoplastic_aux(PetscInt dim, ...);
    // g3: tangent stiffness (consistent tangent for plasticity, or elastic stiffness as approximation)
    void g3_elastoplastic_aux(PetscInt dim, ...);
}
```

The callback reads lambda, mu, rho from aux fields (same as PetscFEElasticityAux). It also needs yield parameters (cohesion, friction angle). Store these as additional aux fields (AUX_COHESION=3, AUX_FRICTION_ANGLE=4) or as additional constants.

**Key challenge:** Plasticity requires state history (accumulated plastic strain, back stress for kinematic hardening). PetscDS pointwise callbacks are stateless -- they see only the current field values and aux fields. To track plastic strain:
- Add a "plastic strain" field as an auxiliary field (6 components for symmetric tensor)
- Update it after each converged time step (not during Newton iterations)
- This is the standard approach in PETSc FEM plasticity (see PETSc ex17.c with plasticity)

This is a significant implementation effort. Verify the approach works with a simple 1-element test before scaling up.

Register a new test `Physics.Elastoplasticity` with a uniaxial compression exceeding yield. Assert: nonzero plastic strain in the aux field after solve.

---

### Phase 4: Eliminate GTEST_SKIP Calls

20 GTEST_SKIP calls remain across 6 files. Handle each:

**Fault tests (test_locked_fault.cpp, test_prescribed_slip.cpp, test_injection_rupture_chain.cpp, test_fault_mesh_manager.cpp):**
These skip when DMPlexCreateBoxMesh, createPlanarFaultLabel, or splitMeshAlongFault fail. In the Docker CI environment these should NOT fail (PETSc 3.22 with ctetgen). If they fail in CI, the test is broken and needs fixing, not skipping. Replace GTEST_SKIP with ASSERT_EQ(ierr, 0) and fix the underlying issue.

If any of these ACTUALLY fail in Docker CI: debug why. The fault mesh splitting is verified by Unit.FaultMeshManager (32 cohesive cells). If it works there but not in integration tests, the setup is different and needs alignment.

**MMS tests (test_mms_elasticity.cpp, test_mms_wave_propagation.cpp):**
These skip because GeomechanicsKernel::residual() and ElastodynamicsKernel::residual() use "strong-form placeholder" math. These PhysicsKernel residual() methods are NEVER used in the actual FEM assembly (which uses PetscDS callbacks). Two options:
A. Fix the PhysicsKernel::residual() methods to compute correct weak-form residuals. This is wasted effort since these methods are not used.
B. Delete the MMS tests that test standalone kernels. Replace with MMS tests that test the actual PetscDS callbacks by running the Simulator with a manufactured solution. This is the correct approach.

Option B is better. Create MMS tests that:
1. Define an exact solution u_exact(x) for elasticity
2. Compute the body force f = -div(C:grad(u_exact)) analytically
3. Set f as the f0 body force (create a new callback or use constants)
4. Solve the FEM problem
5. Compare the FEM solution against u_exact at nodes
6. Verify convergence: refine the mesh and check that the error decreases at the expected rate (O(h^k) for polynomial degree k)

**test_analytical_solutions.cpp:**
Skips because SinglePhaseFlowKernel::residual() uses "pointwise flux sum." Same issue as MMS tests. Delete the skip, fix the kernel residual method, or replace with a Simulator-based test.

---

### Phase 5: BdResidual Architecture Fix

After Phase 1 confirms dynamic rupture works without absorbing BCs, implement the proper fix for the overwrite issue.

**Approach: Per-region DS via DMSetRegionDS**

In PETSc 3.22, DMSetRegionDS allows assigning different PetscDS objects to different labeled regions of the mesh. This lets cohesive cells and external boundary faces have separate callback registrations.

1. Create a separate PetscDS for external boundary faces:
   ```
   PetscDS bdDS;
   PetscDSCreate(comm, &bdDS);
   // Register absorbing BC callbacks on bdDS
   PetscDSSetBdResidual(bdDS, disp_field, AbsorbingBC::f0_absorbing, NULL);
   PetscDSSetBdJacobian(bdDS, disp_field, disp_field, AbsorbingBC::g0_absorbing, ...);
   // Attach bdDS to the boundary labels
   DMSetRegionDS(dm, boundary_label, NULL, bdDS, NULL);
   ```

2. Keep the cohesive callbacks on the main DS (prob):
   ```
   PetscDSSetBdResidual(prob, disp_field, f0_displacement_cohesive, NULL);
   PetscDSSetBdResidual(prob, lagrange_field, f0_lagrange_cohesive, NULL);
   ```

3. DMPlexTSComputeIFunctionFEM will then use the correct DS for each cell/face.

**Before implementing:** Check PETSc 3.22.2 headers to verify DMSetRegionDS exists and has the expected signature. Look at PETSc examples ex17.c or ex62.c for usage patterns.

Create a test that enables BOTH faults and absorbing BCs. Verify that:
- Absorbing BCs fire on external boundary faces (energy leaves the domain)
- Cohesive forces fire on fault interfaces (displacement is discontinuous)

---

### Phase 6: End-to-End Injection Test

The injection source (addInjectionToResidual) exists but has never been end-to-end tested.

Create `tests/integration/test_injection_pressure.cpp`:
1. Single-phase fluid flow (fluid_model = SINGLE_COMPONENT)
2. Injection at center of domain
3. Assert: pressure at injection cell > initial pressure
4. Assert: pressure decays with distance from injection
5. Assert: pressure field is smooth (no NaN, no negative)

This verifies the PetscFEFluidFlow callbacks work end-to-end through the Simulator.

---

### Phase 7: Gmsh Import Verification

The Gmsh import path exists but is untested end-to-end.

1. Create a simple Gmsh .geo file that defines a box with two labeled material regions and one labeled fault surface.
2. Generate the .msh file.
3. Create `tests/integration/test_gmsh_import.cpp`:
   - Load the .msh file via setupGmshGrid
   - Verify physical group labels are mapped to DM labels
   - Verify material properties differ between regions
   - Run a simple elastostatic solve
   - Assert: solution is nonzero and converged

Note: An integration test `Integration.GmshImport` already exists (per TEST_RESULTS.md, test #62). Check if it does real verification or just checks "import succeeds without crash." If the latter, extend it.

---

## Execution Order

1. Phase 1 (dynamic rupture basic) -- highest physics impact
2. Phase 2 (plasticity return mapping fix) -- prerequisite for Phase 3
3. Phase 4 (eliminate GTEST_SKIP) -- test hygiene
4. Phase 3 (plasticity FEM coupling) -- new physics capability
5. Phase 5 (BdResidual architecture fix) -- enables faults + absorbing BCs simultaneously
6. Phase 6 (injection test) -- extends verified physics
7. Phase 7 (Gmsh import) -- extends verified mesh support

## Files Quick Reference

| What | Location |
|---|---|
| Cohesive kernel registration | src/physics/CohesiveFaultKernel.cpp:388 (registerWithDS) |
| Absorbing BC registration | src/core/Simulator.cpp:1893 (overwrites cohesive for same field) |
| Plasticity return mapping | src/domain/geomechanics/PlasticityModel.cpp:256 (DruckerPrager) |
| Plasticity yield function | src/domain/geomechanics/PlasticityModel.cpp:186 |
| FormFunction | src/core/Simulator.cpp:2935 |
| setupPhysics | src/core/Simulator.cpp:1605 |
| Injection source | src/core/Simulator.cpp:2444 (addInjectionToResidual) |
| Gmsh import | src/core/Simulator.cpp:1116 (setupGmshGrid) |
| PhysicsKernel residual (strong-form placeholder) | src/physics/PhysicsKernel.cpp:741 |
| Elasticity PetscDS callbacks (correct) | src/numerics/PetscFEElasticity.cpp |

## Rules

1. Build and test in Docker. Always.
2. Check PETSc 3.22.2 API signatures before calling any function.
3. 68 existing tests must continue to pass.
4. New PetscDS callbacks go in NEW files.
5. Every test has quantitative pass/fail with tolerances.
6. Commit after each passing phase.
7. Do not add GTEST_SKIP. Make tests pass or delete them with documented reason.
8. The BdResidual overwrite (Architecture Bug) is the root cause of cohesive fault failure. Fix it properly, do not paper over it.