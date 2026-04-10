# FSRM Improvement Session Prompt

Read `CLAUDE.md` in the repo root before doing anything. All rules there are binding. Build and test everything in Docker using Dockerfile.ci. All existing 50 tests must pass after every change. Check PETSc 3.22.2 API signatures before calling any function: `grep -rn "FunctionName" /opt/petsc-3.22.2/include/`

## Pre-Flight: Confirm HDF5 Availability

Before starting any phase, verify HDF5 support inside the Docker container:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  # 1. Confirm PETSc was built with HDF5
  grep -q "PETSC_HAVE_HDF5" /opt/petsc-3.22.2/arch-linux-c-opt/include/petscconf.h && echo "HDF5: YES" || echo "HDF5: MISSING"
  
  # 2. Confirm the viewer API exists
  grep -rn "PetscViewerHDF5Open" /opt/petsc-3.22.2/include/ | head -1
  grep -rn "PetscViewerHDF5PushTimestepping" /opt/petsc-3.22.2/include/ | head -1
  grep -rn "PetscViewerHDF5SetTimestep" /opt/petsc-3.22.2/include/ | head -1
  
  # 3. Confirm the shared library links
  ls /opt/petsc-3.22.2/arch-linux-c-opt/lib/libhdf5* 2>/dev/null || \
    ls /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5* 2>/dev/null | head -3
'
```

**If `PETSC_HAVE_HDF5` is missing**, PETSc's configure did not find HDF5. Rebuild the Docker image after fixing Dockerfile.ci. The current Dockerfile.ci already has `libhdf5-openmpi-dev` installed (line 23) and passes `--with-hdf5-include=/usr/include/hdf5/openmpi --with-hdf5-lib="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5"` to PETSc configure (lines 53-54). If this still does not work, the most likely cause is an architecture mismatch on ARM (the lib path is x86_64-specific). Fix by replacing the hardcoded path with `$(dpkg -L libhdf5-openmpi-dev | grep libhdf5.so | head -1 | xargs dirname)`.

**If the `PetscViewerHDF5*` functions are not found in headers**, PETSc compiled without HDF5 despite the configure flags. Check `/opt/petsc-3.22.2/arch-linux-c-opt/lib/petsc/conf/petscvariables` for `PETSC_EXTERNAL_LIB_BASIC` and confirm it contains `-lhdf5`. If not, add `--download-hdf5` to the PETSc configure line in Dockerfile.ci (this downloads and builds HDF5 from source inside PETSc, which always works regardless of system packages):

```dockerfile
RUN ./configure \
        --with-cc=mpicc \
        --with-cxx=mpicxx \
        --with-fc=mpif90 \
        --download-fblaslapack \
        --download-ctetgen \
        --download-hdf5 \
        --with-debugging=0 \
        COPTFLAGS='-O3' \
        CXXOPTFLAGS='-O3' \
        FOPTFLAGS='-O3' \
    && make PETSC_DIR=/opt/petsc-3.22.2 PETSC_ARCH=arch-linux-c-opt all
```

Also add HDF5 to the CMake build. In the top-level `CMakeLists.txt`, after finding PETSc, add:

```cmake
# HDF5 is required for output (provided by PETSc's build or system)
# PETSc's pkg-config already includes HDF5 link flags, so no separate find needed,
# but confirm the header is available:
find_path(HDF5_INCLUDE_DIR hdf5.h HINTS ${PETSC_DIR}/${PETSC_ARCH}/include /usr/include/hdf5/openmpi)
if(HDF5_INCLUDE_DIR)
  message(STATUS "HDF5 headers found: ${HDF5_INCLUDE_DIR}")
  include_directories(${HDF5_INCLUDE_DIR})
else()
  message(WARNING "HDF5 headers not found -- HDF5 output will fail at link time")
endif()
```

Do not proceed to Phase 1 until `PetscViewerHDF5Open` is confirmed available in the Docker image.

## Current State (Honest Assessment)

FSRM compiles, 50 tests pass, but the simulator cannot actually produce usable output. The key gaps, in order of severity:

1. **writeOutput() is a no-op.** It prints "Output at step N would be written here" and returns. There is zero visualization output from any simulation. This means no run produces inspectable results.

2. **No end-to-end physics verification exists.** The Terzaghi pipeline test (`test_terzaghi.cpp` TerzaghiPipelineTest) runs the simulator but does NOT extract the solution or compare it to the analytical answer. The "physics validation" tests are unit tests that call individual callbacks with hand-crafted inputs. Nobody has verified that TSSolve produces a correct pressure field.

3. **CoulombStressTransfer::sampleStressAtFaults() computes garbage strain.** It divides displacement jumps by a hardcoded `ref_len=1.0` instead of recovering the FEM displacement gradient. Lines 150-170 of `src/domain/geomechanics/CoulombStressTransfer.cpp`. Everything downstream (stress, CFS) is wrong.

4. **Poroelasticity has never been verified end-to-end.** The callbacks are unit-tested. The pipeline test runs without crashing. Nobody has extracted the pressure solution and compared it to the Terzaghi series solution.

5. **No prescribed-slip fault mode.** CohesiveFaultKernel supports locked and friction-governed slipping. There is no mode to impose a known displacement jump for kinematic source modeling.

6. **Gmsh import is untested.** `setupGmshGrid()` calls `DMPlexCreateGmshFromFile` and has physical-name mapping code, but no test exercises this path.

7. **Gravity body force exists in callbacks but has no initialization procedure.** `PetscFEElasticityAux::f0_elastostatics_aux` reads gravity from `constants[0]` and rho from aux fields. But there is no two-step procedure to first solve for lithostatic equilibrium, store it as the reference state, then apply perturbations relative to that reference.

## Development Plan

Execute phases in order. Each phase has a verification gate. Do not start a phase until the previous gate passes.

---

### Phase 1: Implement VTK Output via DMView

**Problem:** `writeOutput()` at line 3002 of `Simulator.cpp` is empty. No simulation produces viewable results.

**Task:** Implement HDF5+Xdmf output using PETSc's built-in DMPlex viewers. This is the standard PyLith-compatible output format.

1. In `writeOutput(int step)`, use `PetscViewerHDF5Open` to write:
   - The DM topology (once, at step 0): `DMView(dm, viewer)`
   - The solution vector at each output step: `VecView(solution, viewer)` with proper timestep labeling via `PetscViewerHDF5PushTimestepping(viewer)` and `PetscViewerHDF5SetTimestep(viewer, step)`

2. Before writing, check that `PetscViewerHDF5Open`, `PetscViewerHDF5PushTimestepping`, and `PetscViewerHDF5SetTimestep` exist in PETSc 3.22.2 headers. If not, fall back to `PetscViewerVTKOpen` with `.vtu` extension (DMPlex supports VTK output natively via `DMView`).

3. Add config options:
   ```ini
   [OUTPUT]
   output_format = HDF5    # or VTK
   output_directory = output
   ```

4. Create the output directory in `run()` before `TSSolve`.

**Verification:** Run `config/examples/uniaxial_compression.config`. Confirm:
- An output file is created (`.h5` or `.vtu`).
- Open it in ParaView. Displacement field should be visible.
- For uniaxial compression, displacement should be linear in z.

Add `tests/integration/test_output_file.cpp`: run uniaxial compression, verify output file exists and is nonzero size.

**Files to modify:**
- `src/core/Simulator.cpp` (writeOutput)
- `src/core/ConfigReader.cpp` (parse output_format, output_directory)
- `include/core/ConfigReader.hpp`

**Files to create:**
- `tests/integration/test_output_file.cpp`

---

### Phase 2: End-to-End Terzaghi Verification

**Problem:** `TerzaghiPipelineTest::PoroelasticPipelineCompletes` runs TSSolve but never checks the answer. The test passes if the simulator doesn't crash. That proves nothing about physics.

**Task:** After `sim.run()`, extract the pressure field from the solution vector and compare against the analytical series solution.

1. After `sim.run()` returns, get the solution vector via `sim.getSolution()` (add a public getter if needed).
2. Get the DM and PetscSection to walk DOFs.
3. For each vertex at z = H/2 (mid-height), extract the pressure DOF (field 0 for poroelasticity).
4. Average the pressure values at mid-height.
5. Compare against the analytical Terzaghi series solution (already computed in the test file lines 67-73).
6. Tolerance: 10% relative error (FEM on a coarse 2x2x10 mesh will not be precise, but must be in the right ballpark).

Also fix the test config to include boundary conditions. The current `TerzaghiPipelineTest::SetUp()` writes a config file with NO `[BOUNDARY_CONDITIONS]` section. The Terzaghi problem requires:
- Bottom (z=0): fixed displacement, impermeable (no-flow)
- Top (z=Lz): applied traction (Neumann on displacement), drained (Dirichlet on pressure)
- Laterals: roller BCs (normal displacement fixed, tangential free), impermeable

Without these BCs the simulation is not a Terzaghi problem. Copy the boundary conditions from `config/test_poroelasticity.config` which has them defined correctly.

**Verification gate:** The test extracts pressure at mid-height after TSSolve and it falls between 0.2*P0 and 0.8*P0 (partial consolidation range). If the pressure is exactly zero or exactly P0, the Biot coupling is not working.

**Files to modify:**
- `tests/physics_validation/test_terzaghi.cpp` (TerzaghiPipelineTest: add BCs to config, extract and check solution)
- `include/core/Simulator.hpp` (add `Vec getSolution() const { return solution; }` if not present)

---

### Phase 3: Fix CoulombStressTransfer Strain Recovery

**Problem:** `sampleStressAtFaults()` at lines 110-173 of `CoulombStressTransfer.cpp` computes strain as `(u_pos - u_neg) / 1.0`. This is not FEM strain recovery. It is a displacement jump divided by an arbitrary constant.

**Task:** Replace with proper cell-averaged strain computation.

For each fault vertex:
1. Walk the DMPlex DAG from the vertex to find an adjacent volume cell (height-0 point). Use `DMPlexGetTransitiveClosure(dm, vertex, PETSC_FALSE, &closureSize, &closure)` and filter for cells in the height-0 stratum.
2. Once you have a cell, use `DMPlexVecGetClosure(dm, section, localVec, cell, &closureSize, &closureVals)` to get all DOFs for that cell.
3. Use `DMPlexComputeCellGeometryFEM(dm, cell, NULL, NULL, J, invJ, &detJ)` to get the inverse Jacobian.
4. Tabulate the FE basis function gradients at the cell centroid using `PetscFEGetTabulation` or `PetscFECreateTabulation` (check 3.22.2 API).
5. Compute `du_i/dx_j = sum_k (dN_k/dx_j * u_k_i)` where `dN_k/dx_j` comes from the reference gradient times invJ.
6. Symmetrize: `eps_ij = 0.5 * (du_i/dx_j + du_j/dx_i)`.

**Simpler alternative if basis tabulation is too complex:** Use `DMPlexComputeGradientFVM(dm, cell, localVec, gradient, dm2)` if available in 3.22.2. If not, use a least-squares fit: for each cell, get the closure vertices, get their coordinates and displacement values, fit a linear function `u_i(x) = A_ij * x_j + b_i` by least squares, and `A_ij` gives the displacement gradient directly. This is simpler and avoids FE tabulation.

**Verification:** In `tests/unit/domain/test_coulomb_stress_transfer.cpp`, add a test:
- Create a small 3D simplex mesh (4x4x4).
- Insert a fault, split mesh.
- Set the solution to a known linear displacement field: `u_x = epsilon * x`, `u_y = 0`, `u_z = 0` (uniform extensional strain `epsilon` in x).
- Call `sampleStressAtFaults()`.
- Verify recovered stress: `sigma_xx = (lambda + 2*mu) * epsilon`, `sigma_yy = sigma_zz = lambda * epsilon`, `sigma_xy = sigma_xz = sigma_yz = 0`.
- Tolerance: 1e-6 relative error (FEM should recover linear fields exactly on interior cells, not at boundaries).

**Files to modify:**
- `src/domain/geomechanics/CoulombStressTransfer.cpp` (sampleStressAtFaults only)
- `tests/unit/domain/test_coulomb_stress_transfer.cpp` (add linear field recovery test)

---

### Phase 4: Prescribed Slip on Cohesive Faults

**Problem:** No way to impose a known displacement jump on a fault surface.

**Task:** Add a PRESCRIBED_SLIP mode to CohesiveFaultKernel.

1. Add constant slot 28 for prescribed slip magnitude. Update `COHESIVE_CONST_COUNT` to 29.

2. Write new PetscDS callbacks `f0_cohesive_prescribed` that enforce the Lagrange multiplier constraint `u+ - u- = delta_prescribed` instead of `u+ - u- = 0` (locked) or friction-governed traction. The residual on the Lagrange multiplier field becomes:
   ```
   R_lambda[d] = (u_plus[d] - u_minus[d]) - delta[d]
   ```
   where `delta[d]` is the prescribed slip vector read from constants.

3. The Jacobian blocks are the same as locked mode (identity coupling between lambda and u+/u-). Only the residual changes (constant offset).

4. Add config support:
   ```ini
   [FAULT]
   mode = prescribed_slip
   slip_strike = 1.0
   slip_dip = 0.0
   slip_opening = 0.0
   ```

5. Wire into `setupFaultNetwork()` and `setupPhysics()`.

**Verification:** Create `tests/integration/test_prescribed_slip.cpp`:
- 10x10x100m box, vertical fault at x=5, elastic material.
- Prescribe 0.1m strike-slip.
- After solve, extract displacement on both sides of fault at the center.
- Verify displacement jump = 0.1m in the along-strike direction.
- Verify displacement jump = 0.0 in normal and dip directions.
- Tolerance: 1e-4 m.

**Files to create:**
- `tests/integration/test_prescribed_slip.cpp`
- `config/examples/prescribed_slip_test.config`

**Files to modify:**
- `include/physics/CohesiveFaultKernel.hpp` (add mode, constants)
- `src/physics/CohesiveFaultKernel.cpp` (add callbacks)
- `src/core/Simulator.cpp` (setupFaultNetwork: parse mode)
- `src/core/ConfigReader.cpp` (parse fault mode)

---

### Phase 5: Gmsh Import End-to-End Test

**Problem:** `setupGmshGrid()` and `applyGmshNameMappingsToDM()` exist but have never been exercised.

**Task:**

1. Create a minimal Gmsh .geo file (`meshes/test_box_gmsh.geo`):
   - 10x10x10m box.
   - Two material volumes: "material_upper" (z > 5) and "material_lower" (z <= 5).
   - Six boundary surfaces labeled "xmin", "xmax", "ymin", "ymax", "zmin", "zmax".
   - Mesh it: `gmsh -3 test_box_gmsh.geo -o test_box_gmsh.msh -format msh2` (MSH2 format for PETSc compatibility -- check if 3.22.2 supports MSH4).

2. Install gmsh in the Docker image or pre-generate the .msh and commit it.

3. Create `config/examples/gmsh_box.config`:
   ```ini
   [GRID]
   grid_type = GMSH
   mesh_file = ../meshes/test_box_gmsh.msh
   
   [GMSH_MAPPING]
   material_domains = material_upper:upper, material_lower:lower
   boundaries = zmin:bc_bottom, zmax:bc_top
   ```

4. Create `tests/integration/test_gmsh_import.cpp`:
   - Load the Gmsh mesh.
   - Verify `DMGetDimension` returns 3.
   - Verify cell count is nonzero.
   - Verify boundary labels "zmin" and "zmax" exist (using `DMGetLabel`).
   - Run a simple elastostatic solve (uniaxial compression via Dirichlet on zmin/zmax).
   - Verify nonzero displacement.

**Files to create:**
- `meshes/test_box_gmsh.geo`
- `meshes/test_box_gmsh.msh`
- `config/examples/gmsh_box.config`
- `tests/integration/test_gmsh_import.cpp`

**Files to modify:**
- `Dockerfile.ci` (add `apt-get install gmsh` if generating mesh in CI, or just commit the .msh)

---

### Phase 6: Derived Field Output (Stress, Strain, CFS)

**Problem:** Even with Phase 1 output working, only the raw solution vector (displacement + pressure) gets written. For any analysis, users need derived fields: stress tensor, strain tensor, CFS on a specified receiver orientation.

**Task:**

1. After TSSolve completes, compute cell-centered derived fields:
   - For each cell, get the closure DOFs, compute the displacement gradient (same approach as Phase 3 but over all cells, not just fault vertices), compute strain and stress via Hooke's law.
   - Store as a PETSc Vec on the auxiliary DM (or a new DM with 6-component cell-centered fields for Voigt stress).

2. Compute volumetric CFS for a user-specified receiver fault orientation:
   ```ini
   [OUTPUT]
   output_stress = true
   output_cfs = true
   cfs_receiver_strike = 0.0
   cfs_receiver_dip = 90.0
   cfs_friction = 0.4
   ```

3. Write these derived fields to the same HDF5/VTK output from Phase 1.

4. Create `scripts/plot_stress_field.py`: reads the output, plots a depth-slice of delta-CFS in map view using matplotlib. This is the primary analysis deliverable.

**Verification:** Run prescribed slip (Phase 4) through the full pipeline with stress output enabled. Verify:
- Stress file contains 6-component Voigt tensor.
- CFS file contains scalar values.
- Positive CFS lobes appear at extensional quadrants of fault tips.
- Python script produces a readable figure.

**Files to create:**
- `src/core/DerivedFieldComputer.cpp`
- `include/core/DerivedFieldComputer.hpp`
- `scripts/plot_stress_field.py`

**Files to modify:**
- `src/core/Simulator.cpp` (call derived field computation after run(), write to output)
- `src/core/ConfigReader.cpp` (parse output options)

---

### Phase 7: Lithostatic Equilibrium Initialization

**Problem:** All simulations start from zero stress. For problems with gravity, you need the lithostatic state as the starting point.

**Task:** Implement a two-step initialization:

1. When `config.gravity > 0.0`:
   a. Before the main TSSolve, solve a quasi-static elastostatic problem with the gravity body force (one Newton iteration, dt=1.0, single step). This uses the existing `PetscFEElasticityAux::f0_elastostatics_aux` callback which already supports gravity.
   b. Store the resulting solution as `solution_lithostatic`.
   c. Store the resulting stress field as `stress_reference` (using the derived field computer from Phase 6).
   d. Reset the solution to `solution_lithostatic` as the initial condition for the main simulation.
   e. When computing delta-CFS in Phase 6, subtract `stress_reference` to get the change from lithostatic.

2. Add config:
   ```ini
   [SIMULATION]
   gravity = 9.81
   K0 = 0.5    # lateral earth pressure coefficient (for initial horizontal stress)
   ```

3. K0 boundary conditions for the lithostatic step: lateral boundaries get roller BCs (zero normal displacement), bottom gets fixed, top gets free. This naturally produces sigma_zz = -rho*g*(H-z) and sigma_xx = sigma_yy = K0*sigma_zz if K0 is imposed through initial horizontal traction. Alternatively, let the Poisson effect determine the horizontal stress (K0 = nu/(1-nu)) which is what the current elastic solver naturally produces with roller BCs.

**Verification:** `tests/physics_validation/test_lithostatic_equilibrium.cpp`:
- 1x1x1000m column, rho=2650, g=9.81, nu=0.25.
- After lithostatic solve, extract sigma_zz at z=500m.
- Expected: sigma_zz = -rho*g*500 = -13.0 MPa.
- Tolerance: 5% (FEM discretization error on coarse mesh).

**Files to create:**
- `tests/physics_validation/test_lithostatic_equilibrium.cpp`
- `config/examples/lithostatic_column.config`

**Files to modify:**
- `src/core/Simulator.cpp` (add lithostatic pre-solve in run())

---

## General Rules

1. Every new .cpp needs a CMakeLists.txt entry.
2. Every test has quantitative pass/fail with tolerances. No "looks right."
3. New PetscDS callbacks go in NEW files. Never modify PetscFEElasticity.cpp, PetscFEPoroelasticity.cpp, PetscFEFluidFlow.cpp.
4. Never modify FaultMeshManager::splitMeshAlongFault or CohesiveFaultKernel::registerWithDS.
5. Never change the DS/BC ordering in setupFields().
6. Never change the unified constants array layout [0-26].
7. `ctest --output-on-failure` after every change. Zero regressions.
8. If a PETSc function does not exist in 3.22.2, find the one that does. `grep -rn` the headers.
9. No Python in the C++ codebase. Python scripts are separate tools in `scripts/`.
10. Commit after each passing phase.