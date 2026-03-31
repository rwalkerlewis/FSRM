# FSRM: Integrate Faults into Simulator Pipeline

Read CLAUDE.md first. Build and test everything in Docker.

## Goal

Wire the existing fault components (FaultMeshManager, CohesiveFaultKernel, FaultCohesiveDyn)
into the Simulator so that `./fsrm -c config/test_locked_fault.config` runs an elastostatic
solve with a locked cohesive fault embedded in the mesh.

## Bug Fix First: Constants Slot Collision

CohesiveFaultKernel::COHESIVE_CONST_MODE is currently 24, which collides with
rho_fluid at unified_constants[24] in Simulator::setupPhysics().

In `include/physics/CohesiveFaultKernel.hpp`, change:
```cpp
static constexpr PetscInt COHESIVE_CONST_MODE  = 25;  // was 24
static constexpr PetscInt COHESIVE_CONST_MU_F  = 26;  // was 25
static constexpr PetscInt COHESIVE_CONST_COUNT = 27;  // was 26
```

In `src/core/Simulator.cpp`, expand the unified constants array:
```cpp
PetscScalar unified_constants[27] = {0};  // was 25
```
Change the PetscDSSetConstants call to pass 27 instead of 25.
rho_fluid stays at [24]. Nothing else moves.

Do this FIRST before any other changes. Build and verify existing tests still pass.


## Subagent A: Wire FaultMeshManager into Simulator::setupFaultNetwork()

File: `src/core/Simulator.cpp`

The current setupFaultNetwork() creates an empty FaultNetwork and returns.
Replace it with real fault setup. The call order matters:

```cpp
PetscErrorCode Simulator::setupFaultNetwork() {
    PetscFunctionBeginUser;
    if (!config.enable_faults) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    if (rank == 0) PetscPrintf(comm, "Setting up fault network...\n");

    // 1. Create fault mesh manager
    fault_mesh_manager_ = std::make_unique<FaultMeshManager>(comm);

    // 2. Read fault geometry from config
    //    ConfigReader should have parsed [FAULT] section with:
    //    strike, dip, center_x/y/z, length, width, tolerance
    //    For now, use defaults if not in config
    double strike = 0.0;        // radians
    double dip = M_PI / 2.0;    // vertical fault
    double center[3] = {0.0, 0.0, 0.0};
    double length = 1e10;       // large enough to cut the whole mesh
    double width = 1e10;
    double tol = 0.15;

    // TODO: Read from config once FaultConfig parsing is implemented
    // For now, hardcode a vertical fault at x=Lx/2
    center[0] = grid_config.Lx / 2.0;
    center[1] = grid_config.Ly / 2.0;
    center[2] = grid_config.Lz / 2.0;
    length = 2.0 * std::max({grid_config.Lx, grid_config.Ly, grid_config.Lz});
    width = length;

    // 3. Label fault faces on the DM
    DMLabel fault_label = nullptr;
    ierr = fault_mesh_manager_->createPlanarFaultLabel(
        dm, &fault_label, strike, dip, center, length, width, tol);
    if (ierr != 0) {
        if (rank == 0) PetscPrintf(comm, "Warning: fault labeling failed, skipping fault setup\n");
        fault_mesh_manager_.reset();
        PetscFunctionReturn(0);
    }

    // 4. Split mesh along fault (inserts cohesive cells)
    //    THIS CHANGES THE DM TOPOLOGY
    ierr = fault_mesh_manager_->splitMeshAlongFault(&dm, "fault");
    if (ierr != 0) {
        if (rank == 0) PetscPrintf(comm, "Warning: mesh splitting failed, skipping fault setup\n");
        fault_mesh_manager_.reset();
        PetscFunctionReturn(0);
    }

    // 5. Create FaultCohesiveDyn and extract topology from cohesive cells
    cohesive_fault_ = std::make_unique<FaultCohesiveDyn>();
    ierr = fault_mesh_manager_->extractCohesiveTopology(dm, cohesive_fault_.get());
    CHKERRQ(ierr);

    // 6. Set friction model
    //    TODO: Read from config. Default: slip-weakening
    auto friction = std::make_unique<SlipWeakeningFriction>();
    friction->setStaticCoefficient(0.6);
    friction->setDynamicCoefficient(0.4);
    friction->setCriticalSlipDistance(0.4);
    cohesive_fault_->setFrictionModel(std::move(friction));

    // 7. Set initial traction (pre-stress on fault)
    //    For locked fault test: doesn't matter, but set something physical
    cohesive_fault_->setUniformInitialTraction(25e6, 0.0, -50e6);

    // 8. Initialize fault state
    cohesive_fault_->initialize();

    // 9. Create cohesive kernel
    cohesive_kernel_ = std::make_unique<CohesiveFaultKernel>();
    cohesive_kernel_->setMode(true);  // locked
    cohesive_kernel_->setFrictionCoefficient(0.6);

    if (rank == 0) {
        PetscPrintf(comm, "Fault network setup complete: %zu fault vertices\n",
                    cohesive_fault_->numVertices());
    }

    PetscFunctionReturn(0);
}
```

CRITICAL: This must be called AFTER setupDM() but BEFORE setupFields(), because
splitMeshAlongFault changes the DM topology. Check the call order in
Simulator::initializeFromConfigFile() or Simulator::run() and adjust if needed.

Also add the necessary includes at the top of Simulator.cpp if not already present:
```cpp
#include "domain/geomechanics/PyLithFault.hpp"  // FaultCohesiveDyn, SlipWeakeningFriction
#include "numerics/FaultMeshManager.hpp"
#include "physics/CohesiveFaultKernel.hpp"
```

These may already be included via Simulator.hpp -- check before adding.


## Subagent B: Add Lagrange Multiplier Field and Register Cohesive Callbacks

### In Simulator::setupFields():

When faults are enabled and cohesive_kernel_ exists, add a Lagrange multiplier field
for the traction on cohesive cells:

```cpp
// After adding the displacement field...
if (config.enable_faults && cohesive_kernel_) {
    PetscFE fe_lagrange;
    // Lagrange multiplier: 3-component vector (traction) on cohesive cells
    ierr = PetscFECreateDefault(comm, 3, 3, PETSC_FALSE, "lagrange_", -1, &fe_lagrange);
    CHKERRQ(ierr);
    ierr = DMAddField(dm, nullptr, (PetscObject)fe_lagrange); CHKERRQ(ierr);
    fe_fields.push_back(fe_lagrange);
    if (rank == 0) PetscPrintf(comm, "Added Lagrange multiplier field for fault traction\n");
}
```

NOTE: PETSc's DMPlex cohesive cell handling may require the Lagrange multiplier to be
set up differently -- it might need DMSetRegionDS or a separate DS for the cohesive
cells. Check PETSc examples that use DMPlexConstructCohesiveCells + PetscFE.
If the standard DMAddField approach doesn't work (PETSc errors about hybrid cells),
try using DMPlexGetHybridBounds or DMGetLabel("fault") to restrict the field.

### In Simulator::setupPhysics():

After registering elasticity callbacks, if faults are enabled:

```cpp
if (config.enable_faults && cohesive_kernel_) {
    // Determine field indices
    PetscInt disp_field = displacement_field_idx;  // computed earlier
    PetscInt lagrange_field = disp_field + 1;      // Lagrange multiplier is next field

    // Register cohesive callbacks via CohesiveFaultKernel::registerWithDS
    // This expands the constants array and registers PetscDSSetBdResidual/BdJacobian
    ierr = cohesive_kernel_->registerWithDS(prob, disp_field, lagrange_field);
    CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Registered cohesive fault callbacks on fields %d (disp), %d (lagrange)\n",
                    disp_field, lagrange_field);
    }
}
```


## Subagent C: Locked Fault Integration Test

Create `tests/integration/test_locked_fault.cpp`:

```cpp
/**
 * @file test_locked_fault.cpp
 * @brief Integration test: elastostatic solve with locked cohesive fault
 *
 * Verifies that inserting a locked fault (u+ = u-) into the mesh does not
 * change the elastostatic solution compared to a no-fault baseline.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <petscsys.h>
#include <petscdmplex.h>

#include "numerics/FaultMeshManager.hpp"
#include "numerics/PetscFEElasticity.hpp"
#include "physics/CohesiveFaultKernel.hpp"
#include "domain/geomechanics/PyLithFault.hpp"

using namespace FSRM;

class LockedFaultTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(LockedFaultTest, LockedFaultPreservesElastostatics) {
    // Create a 3D box mesh
    DM dm = nullptr;
    PetscInt faces[3] = {4, 4, 4};
    PetscReal lower[3] = {0.0, 0.0, 0.0};
    PetscReal upper[3] = {1.0, 1.0, 1.0};
    PetscErrorCode ierr = DMPlexCreateBoxMesh(
        PETSC_COMM_WORLD, 3, PETSC_FALSE, faces, lower, upper,
        nullptr, PETSC_TRUE, &dm);
    if (ierr != 0 || !dm) {
        GTEST_SKIP() << "DMPlexCreateBoxMesh failed";
    }

    // Insert a horizontal fault at z=0.5
    FaultMeshManager mgr(PETSC_COMM_WORLD);
    DMLabel fault_label = nullptr;
    const double center[3] = {0.5, 0.5, 0.5};
    ierr = mgr.createPlanarFaultLabel(dm, &fault_label, 0.0, 0.5*M_PI,
                                       center, 2.0, 2.0, 0.15);
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "createPlanarFaultLabel failed";
    }

    ierr = mgr.splitMeshAlongFault(&dm, "fault");
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "splitMeshAlongFault failed";
    }

    // Extract cohesive topology
    FaultCohesiveDyn fault;
    ierr = mgr.extractCohesiveTopology(dm, &fault);
    if (ierr != 0) {
        DMDestroy(&dm);
        GTEST_SKIP() << "extractCohesiveTopology failed";
    }

    fault.setFrictionModel(std::make_unique<SlipWeakeningFriction>());
    fault.initialize();

    // Create cohesive kernel in locked mode
    CohesiveFaultKernel kernel;
    kernel.setMode(true);  // locked: u+ - u- = 0
    kernel.setFrictionCoefficient(0.6);

    // Verify mesh was modified
    PetscInt num_cells;
    ierr = DMPlexGetHeightStratum(dm, 0, nullptr, &num_cells); CHKERRQ(ierr);
    EXPECT_GT(num_cells, 64) << "Mesh should have more cells after cohesive insertion";

    // If we got here, the mesh splitting and cohesive topology extraction worked
    // Full solve test would require setting up PetscFE + PetscDS + SNES,
    // which is complex. For now, verify the components wire together.
    EXPECT_GE(fault.numVertices(), 0u);

    DMDestroy(&dm);
}
```

Register this test in `tests/CMakeLists.txt`:
- Add `integration/test_locked_fault.cpp` to INTEGRATION_TEST_SOURCES
- Add CTest registration: `add_test(NAME Integration.LockedFault COMMAND run_integration_tests --gtest_filter=LockedFaultTest.*)`
- Add to set_tests_properties for integration label and timeout

Also create `config/test_locked_fault.config`:
```ini
[SIMULATION]
name = locked_fault_test
start_time = 0.0
end_time = 1.0
dt_initial = 1.0
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true
enable_faults = true
rtol = 1.0e-8
atol = 1.0e-10
max_nonlinear_iterations = 20

[GRID]
nx = 4
ny = 4
nz = 4
Lx = 1.0
Ly = 1.0
Lz = 1.0

[ROCK]
density = 2700.0
youngs_modulus = 50.0e9
poissons_ratio = 0.25

[FAULT]
enabled = true
strike = 0.0
dip = 90.0
center_x = 0.5
center_y = 0.5
center_z = 0.5
length = 2.0
width = 2.0
friction_static = 0.6
friction_dynamic = 0.4
critical_slip_distance = 0.4
```


## Main Agent: After Subagents Complete

### 1. Fix call order in Simulator

Check `Simulator::initializeFromConfigFile()` or wherever the setup methods are called.
The order MUST be:
```
initializeFromConfigFile()
  -> parse config
  -> setupDM()              // create base mesh
  -> setupFaultNetwork()    // split mesh (changes DM)
  -> setupFields()          // create PetscFE on final mesh
  -> setupPhysics()         // register PetscDS callbacks
  -> setupTimeStepper()     // create TS
  -> setupSolvers()         // configure SNES/KSP
```

If setupFaultNetwork() is currently called after setupFields(), move it before.

### 2. Build and test

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON &&
  make -j$(nproc) &&
  echo "=== ALL TESTS ===" &&
  ctest --output-on-failure 2>&1 | tail -40 &&
  echo "=== LOCKED FAULT CONFIG ===" &&
  ./fsrm -c ../config/test_locked_fault.config 2>&1 | tail -30
'
```

### 3. Debug common PETSc cohesive cell issues

If PETSc errors on DMAddField for the Lagrange multiplier on a hybrid mesh:
- Try DMPlexGetHybridBounds to check if PETSc recognizes the cohesive cells
- May need DMSetRegionDS to set a separate DS for cohesive vs bulk cells
- Check PETSc fault examples (ex52, ex62) for the correct DMPlex + cohesive setup

If registerWithDS fails because PetscDSSetBdResidual doesn't work with the field layout:
- Verify the Lagrange multiplier field index is correct
- Check that PetscDS has the right number of fields after DMAddField

### 4. If full solve works, verify locked fault is transparent

Run the same elastostatic problem with and without enable_faults:
- With enable_faults=false: baseline displacement
- With enable_faults=true, locked mode: should produce identical displacement
- If they differ, the cohesive constraint (u+ - u- = 0) isn't being applied correctly


## Do NOT
- Modify PetscFEElasticity, PetscFEPoroelasticity, or PetscFEFluidFlow
- Modify friction law implementations
- Modify CohesiveFaultKernel::registerWithDS internals
- Change IMEX transition logic
- Delete or gut existing tests