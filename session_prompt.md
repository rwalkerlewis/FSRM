# FSRM Session Prompt: Coulomb Stress Transfer for Paleoseismicity-Gold Proxy Modeling

## READ FIRST

Read `CLAUDE.md` in the repository root completely before starting. All rules there are binding. In particular:

- NEVER change DS/BC ordering in setupFields()
- NEVER modify existing callback math in PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, or PetscFEFluidFlow.cpp
- NEVER modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS
- Build and test everything in Docker using Dockerfile.ci
- All existing 50 tests must continue to pass after every change
- Check PETSc 3.22.2 API signatures before calling any PETSc function: `grep -rn "FunctionName" /opt/petsc-3.22.2/include/`

## Project Context

FSRM is being adapted to model Coulomb stress transfer from prescribed fault slip events onto receiver faults in 3D, for comparison with gold deposit locations in the Yilgarn Craton, Western Australia. The scientific question: do gold deposit clusters spatially correlate with modeled zones of persistent positive delta-CFS (aftershock damage zones) at fault jogs and stepovers?

This is NOT dynamic rupture modeling. This is quasi-static forward modeling: prescribe slip on a source fault, compute the resulting stress perturbation throughout the domain, evaluate delta-CFS on receiver fault surfaces, and output the results for GIS overlay with deposit locations.

## Current State Assessment

**Working and tested:**
- Elastostatics PetscDS callbacks (f0, f1, g3) with verified convergence
- Cohesive fault mesh splitting (PyLith workflow: submesh, subpoint map, DMPlexConstructCohesiveCells)
- CoulombStressTransfer class: Hooke stress computation, fault projection, delta-CFS (unit tested)
- Auxiliary field DM for heterogeneous materials (depth-based layering only)
- Gmsh import code path (DMPlexCreateGmshFromFile + physical name mapping -- coded but UNTESTED)
- Rate-state and slip-weakening friction laws (unit tested)
- Boundary condition labeling on structured grids

**Broken or placeholder:**
- CoulombStressTransfer::sampleStressAtFaults() computes strain using displacement jumps divided by a hardcoded ref_len=1.0. This is garbage. It needs proper FEM strain recovery from the displacement gradient at quadrature points or cell centroids.
- No gravity or lithostatic prestress. Solver starts from zero stress state.
- Gmsh mesh import has never been tested end-to-end.
- Material assignment from Gmsh physical groups has never been tested.
- Only single-fault geometry supported (one planar fault via createPlanarFaultLabel). No multi-fault from Gmsh labels.
- No VTK output of volumetric stress/strain fields for GIS comparison.
- No prescribed-slip boundary condition on fault surfaces (only locked or friction-governed slipping).

**Architectural constraints you must respect:**
- The unified constants array is 27 elements [0-26]. Do NOT change its layout.
- PetscDS callback signatures are fixed. New physics go in NEW files.
- The DS/BC ordering in setupFields() (DMCreateDS -> DMAddBoundary -> DMSetLocalSection null -> DMSetUp) was debugged over multiple sessions. Do not touch it.
- FaultMeshManager::splitMeshAlongFault uses the PyLith workflow and is correct. Do not modify.

## Development Plan -- Execute in Order

Each phase produces a testable artifact. Do not start a phase until the previous phase's test passes. After each phase, run the full test suite (`ctest --output-on-failure`) to confirm no regressions.

### Phase 1: Fix CoulombStressTransfer Strain Recovery

**Problem:** `sampleStressAtFaults()` in `src/domain/geomechanics/CoulombStressTransfer.cpp` does not compute strain from the FEM solution. It uses `(u_pos - u_neg) / 1.0` as a strain proxy, which is dimensionally wrong and physically meaningless for anything except displacement-jump detection.

**Task:** Replace the strain computation with proper FEM gradient recovery. For each fault vertex:
1. Find an adjacent volume cell (walk the DAG: vertex -> edges -> faces -> cell, or use DMPlexGetTransitiveClosure on the vertex and filter for height-0 cells).
2. Use `DMPlexComputeCellGeometryFEM()` to get the Jacobian and inverse Jacobian at the cell centroid.
3. Use `DMPlexVecGetClosure()` to get the displacement DOFs for that cell.
4. Compute the displacement gradient using basis function gradients at the centroid (tabulate the FE basis, apply inverse Jacobian, contract with DOF values).
5. Symmetrize to get the strain tensor: eps_ij = 0.5 * (du_i/dx_j + du_j/dx_i).
6. The rest of the pipeline (Hooke's law, fault projection, delta-CFS) is already correct.

**Alternative simpler approach:** Use `DMPlexComputeGradientFVM()` if the mesh supports it (check PETSc 3.22.2 API), which gives cell-averaged gradients directly. This avoids basis function tabulation.

**Verification:** Modify `tests/unit/domain/test_coulomb_stress_transfer.cpp` to set up a small 3D mesh with a known linear displacement field (e.g., u_x = a*x, u_y = 0, u_z = 0), run sampleStressAtFaults(), and verify that the recovered stress matches Hooke's law analytically (sigma_xx = (lambda + 2*mu)*a, sigma_yy = sigma_zz = lambda*a, shear = 0). Tolerance: 1e-10 relative error.

**Files to modify:**
- `src/domain/geomechanics/CoulombStressTransfer.cpp` (sampleStressAtFaults only)
- `tests/unit/domain/test_coulomb_stress_transfer.cpp` (add gradient recovery test)

**Files NOT to modify:** Everything else.

---

### Phase 2: Prescribed Slip on Fault Surfaces

**Problem:** The CohesiveFaultKernel only supports two modes: locked (zero displacement jump) and slipping (friction-governed). There is no mode for prescribed displacement jump (kinematic slip), which is what we need to impose a known earthquake slip distribution on a source fault and compute the resulting stress field.

**Task:** Add a PRESCRIBED_SLIP mode to CohesiveFaultKernel.

1. Add a new constant slot. The constants array currently uses slots 0-27. Add slot 28 for the prescribed slip magnitude. Update COHESIVE_CONST_COUNT to 29. Add a `setPrescribedSlip(double strike_slip, double dip_slip, double opening)` method.

2. Add new static PetscDS callbacks `f0_prescribed_slip` and `g0_prescribed_slip` that enforce the constraint: `u_plus - u_minus = prescribed_delta`. The residual is: `R_lambda = (u+ - u-) - delta_prescribed`. The Jacobian block g0 for (lambda, u+) = I, for (lambda, u-) = -I. These callbacks go in a NEW section of CohesiveFaultKernel.cpp (or a new file `PrescribedSlipKernel.cpp`).

3. Wire into setupPhysics(): when mode is PRESCRIBED_SLIP, register the prescribed-slip callbacks instead of the locked or friction callbacks.

4. Add a config option: `[FAULT] mode = prescribed_slip` with `slip_strike = 1.0`, `slip_dip = 0.0`, `slip_opening = 0.0`.

**Verification:** Create `config/examples/prescribed_slip_test.config`: a 10x10x100m box with a vertical fault at x=5m, prescribe 1m of strike-slip (along-fault displacement jump = 1m in y). Verify:
- Displacement is discontinuous across the fault by exactly 1m in y.
- Stress field shows the expected elastic response (concentrated shear stress at fault tips).
- Far-field displacement approaches zero (or matches the elastic solution for a fault in a half-space).

Add `tests/integration/test_prescribed_slip.cpp` with quantitative checks.

**Files to create:**
- `tests/integration/test_prescribed_slip.cpp`
- `config/examples/prescribed_slip_test.config`

**Files to modify:**
- `include/physics/CohesiveFaultKernel.hpp` (add PRESCRIBED_SLIP mode, new constant slots)
- `src/physics/CohesiveFaultKernel.cpp` (add prescribed_slip callbacks)
- `src/core/Simulator.cpp` (setupFaultNetwork: read prescribed slip from config)
- `src/core/ConfigReader.cpp` (parse prescribed_slip config)
- `include/core/ConfigReader.hpp` (add prescribed_slip fields)

---

### Phase 3: Lithostatic Prestress (Gravity Body Force)

**Problem:** All simulations currently start from zero stress. For meaningful CFS computation on faults at depth, we need a lithostatic initial stress state: sigma_zz = -rho*g*z, sigma_xx = sigma_yy = K0 * sigma_zz, where K0 is the lateral earth pressure coefficient.

**Task:** 

1. The elastostatics f0 callback needs a gravity body force term: f0_i += rho * g_i. The auxiliary-field-aware callback `PetscFEElasticityAux::f0_elastostatics_aux` already reads rho from the auxiliary vector (slot AUX_RHO = 2). Add gravity as an additional constant (use an unused slot in the constants array, or add slot 28/29). The f0 body becomes: `f0[dim-1] += rho * g` where g = -9.81 m/s^2 (or configurable).

   IMPORTANT: Do NOT modify `PetscFEElasticity::f0_elastostatics`. Create a NEW callback `PetscFEElasticityAux::f0_elastostatics_gravity_aux` in a new file or in `PetscFEElasticityAux.cpp`.

2. Add a two-step initialization procedure:
   a. First, solve the elastostatic problem with gravity to obtain the lithostatic stress state.
   b. Store this as the reference/initial stress.
   c. Then apply tectonic loading or prescribed slip and compute delta-CFS relative to the lithostatic reference.

3. Add config options: `[SIMULATION] enable_gravity = true`, `[SIMULATION] gravity = 9.81`, `[SIMULATION] K0 = 0.5`.

**Verification:** Create `tests/physics_validation/test_lithostatic_stress.cpp`:
- 1D column (1x1x1000m), rho=2650, g=9.81, nu=0.25 (so K0 = nu/(1-nu) = 1/3).
- Expected: sigma_zz(z) = -rho*g*(H-z), sigma_xx = sigma_yy = K0*sigma_zz.
- Verify at 10 depth points. Tolerance: 1% relative error (FEM approximation).

**Files to create:**
- `src/numerics/PetscFEElasticityGravity.cpp`
- `include/numerics/PetscFEElasticityGravity.hpp`
- `tests/physics_validation/test_lithostatic_stress.cpp`
- `config/examples/lithostatic_column.config`

---

### Phase 4: Gmsh Multi-Fault Import (End-to-End Verification)

**Problem:** The Gmsh import path (`setupGmshGrid`) exists but has never been tested. For the WA gold problem, we need to import a mesh with multiple named fault surfaces and material domains from Gmsh physical groups.

**Task:**

1. Create a Gmsh .geo script (`meshes/two_fault_box.geo`) that defines:
   - A 20km x 20km x 15km box.
   - Two non-intersecting vertical fault surfaces as 2D physical surfaces ("fault_source", "fault_receiver") with a 5km stepover geometry.
   - Two material domains: "upper_crust" (0-10km depth) and "lower_crust" (10-15km depth).
   - Boundary surfaces labeled for BCs.
   - Generate the .msh file and commit it to `meshes/`.

2. Verify the full pipeline:
   a. `setupGmshGrid()` loads the mesh and maps physical names to DM labels.
   b. `setupFaultNetwork()` identifies fault surfaces by label name (not by planar geometry), calls `identifyFaultSurface()` with "fault_source", then splits the mesh.
   c. Material properties are assigned by region label (extend `populateAuxFieldsByDepth` or create `populateAuxFieldsByLabel` that reads from Gmsh physical group labels on cells).

3. Modify `setupFaultNetwork()` to support Gmsh-labeled faults:
   - If `grid_type = GMSH` and `[FAULT] label = fault_source`, use `identifyFaultSurface(dm, "fault_source", &label)` instead of `createPlanarFaultLabel()`.
   - Support multiple `[FAULT_N]` sections (N=1,2,...) for multi-fault problems. Each fault gets its own FaultMeshManager split. NOTE: Multiple sequential splits require care -- the DM topology changes after each split, so labels from earlier faults may shift. Test this carefully. If sequential splitting fails, fall back to a single-fault-at-a-time approach and document the limitation.

**Verification:** `tests/integration/test_gmsh_two_fault.cpp`:
- Load two_fault_box.msh.
- Verify both fault labels are found.
- Verify mesh splitting produces cohesive cells on the source fault.
- Verify material properties differ between upper and lower crust cells.
- Run a prescribed-slip problem on the source fault and verify non-zero delta-CFS on the receiver fault.

**Files to create:**
- `meshes/two_fault_box.geo`
- `meshes/two_fault_box.msh` (generated)
- `tests/integration/test_gmsh_two_fault.cpp`

**Files to modify:**
- `src/core/Simulator.cpp` (setupFaultNetwork: Gmsh label path; setupDM: material-by-label)
- `src/core/ConfigReader.cpp` (parse multiple [FAULT_N] sections)

---

### Phase 5: Volumetric Stress/CFS Output for GIS Overlay

**Problem:** Currently there is no way to output the 3D stress field or a CFS field on a grid for comparison with gold deposit coordinates.

**Task:**

1. After solving, compute the full stress tensor at every cell centroid (using the same gradient recovery from Phase 1, but over all cells, not just fault vertices). Store as a cell-based field.

2. Compute delta-CFS on a user-specified receiver fault orientation at every cell (not just on fault surfaces). This gives a volumetric CFS map: at each point in the domain, "if a fault with this orientation existed here, would it be brought closer to failure?" This is the standard approach in earthquake triggering studies (King et al. 1994).

3. Output to VTK (VTKViewer already has stubs in FSRM) or HDF5:
   - Cell-centered stress tensor (6 components in Voigt).
   - Cell-centered delta-CFS for the specified receiver orientation.
   - Cell centroid coordinates (for GIS projection).

4. Add config:
   ```
   [OUTPUT]
   output_stress_field = true
   output_coulomb_field = true
   receiver_strike = 0.0      # degrees
   receiver_dip = 90.0        # degrees
   receiver_friction = 0.4
   output_format = VTK
   ```

**Verification:** Run the two-fault stepover problem from Phase 4 with 1m of prescribed strike-slip on the source fault. Export the CFS field. Verify:
- Positive CFS lobes appear at the extensional quadrant of the fault tips (consistent with King et al. 1994 Figure 1).
- Negative CFS lobes appear at the compressional quadrant.
- The receiver fault, if located in the extensional stepover, shows positive CFS.

Write a Python script `scripts/plot_coulomb_field.py` that reads the VTK/HDF5 output and produces a map-view plot at a specified depth slice, overlaid with fault traces. This is the deliverable that gets compared with GSWA gold deposit coordinates.

**Files to create:**
- `src/io/StressFieldOutput.cpp`
- `include/io/StressFieldOutput.hpp`
- `scripts/plot_coulomb_field.py`

**Files to modify:**
- `src/core/Simulator.cpp` (call stress output after solve)
- `src/core/ConfigReader.cpp` (parse output options)

---

### Phase 6: Integration Example -- Yilgarn Fault Stepover

**Problem:** Need a complete working example that demonstrates the full pipeline from mesh to CFS map.

**Task:** Create `config/examples/yilgarn_stepover.config` and accompanying files:

1. A simplified 2D cross-section (or 3D slab) representing a portion of the Boulder-Lefroy fault zone with:
   - A main fault (source) with a compressional jog/stepover.
   - A secondary fault (receiver) in the stepover region.
   - Depth-dependent material properties (upper crust: E=70 GPa, nu=0.25, rho=2700; lower crust: E=100 GPa, nu=0.28, rho=2900).
   - Gravity-loaded lithostatic initial stress.

2. Prescribe 1m of reverse slip on the main fault (consistent with the Yilgarn's compressional regime).

3. Compute and output:
   - Delta-CFS on the receiver fault.
   - Volumetric CFS field at 5km depth slice.
   - Stress tensor field.

4. Python post-processing script that:
   - Plots the CFS map at 5km depth.
   - Overlays fault traces.
   - Marks hypothetical gold deposit locations in the positive-CFS lobes.
   - Produces a publication-quality figure.

This example serves as the template for the actual GSWA data analysis.

**Files to create:**
- `config/examples/yilgarn_stepover.config`
- `meshes/yilgarn_stepover.geo`
- `scripts/plot_yilgarn_stepover.py`

---

## General Rules for All Phases

1. Every new .cpp file needs a corresponding entry in `CMakeLists.txt` (either `tests/CMakeLists.txt` for tests or the main `CMakeLists.txt` for source).
2. Every new test must have quantitative pass/fail criteria with numerical tolerances. No "looks right" tests.
3. New PetscDS callbacks go in NEW files. Do not modify existing callback files.
4. Use the Docker build for all compilation and testing.
5. After each phase, run `ctest --output-on-failure` and confirm all tests (old + new) pass before proceeding.
6. Commit after each passing phase with a descriptive message referencing the phase number.
7. If a PETSc function does not exist in 3.22.2 headers, find the equivalent that does. Do not guess.
8. No Python in the C++ codebase. Python scripts are separate post-processing tools in `scripts/`.
9. When in doubt about PETSc DMPlex API, look at `ex17.c` (elasticity with aux fields), `ex56.c` (elasticity with BCs), and `ex62.c` (cohesive cells) in the PETSc source.
10. The strain recovery in Phase 1 is the foundation for everything else. Get it right. If the gradient at a cell centroid does not match a known linear field to machine precision, stop and debug before proceeding.