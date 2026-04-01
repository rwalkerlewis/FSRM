# FSRM: Verify BCs, Poroelastic Consolidation, Injection Source

Read CLAUDE.md first. Build and test everything in Docker.

## Phase 1: Verify Uniaxial Compression Works

Build and run:
```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON &&
  make -j$(nproc) &&
  echo "=== TESTS ===" &&
  ctest --output-on-failure 2>&1 | tail -10 &&
  echo "=== UNIAXIAL COMPRESSION ===" &&
  ./fsrm -c ../config/examples/uniaxial_compression.config \
    -snes_monitor -ksp_monitor_short -snes_converged_reason \
    -pc_type lu 2>&1 | tail -40
'
```

Expected:
- "Solution norm after BC insertion" should be NONZERO (the compression BC sets u_z=-0.001 on top)
- SNES should show nonzero initial FNORM
- SNES should converge in 1 iteration (linear problem with direct solver)
- NOT expected: FNORM = 0 or "converged due to FNORM_ABS"

If the solution norm is still zero after BC insertion, the constrained section
is not being built correctly. Debug by adding before the solve:
```cpp
PetscSection section;
DMGetLocalSection(dm, &section);
PetscSectionView(section, PETSC_VIEWER_STDOUT_WORLD);
```
Check that some DOFs show as constrained (cdof > 0).

STOP HERE if this doesn't work. Fix it before proceeding to Phase 2.


## Phase 2: Terzaghi Consolidation (Poroelasticity End-to-End)

This validates the Biot coupling callbacks. A 1D column under sudden load,
drained at the top. Pore pressure starts at the applied load value and
dissipates over time as fluid drains out the top.

### 2a: Add configurable BCs for the poroelastic case

The current setupBoundaryConditions() hardcodes elastostatics BCs. The
Terzaghi problem needs DIFFERENT BCs:

For poroelasticity (field 0 = pressure, field 1 = displacement):
- z_min (bottom): u = 0 (fixed), impermeable (no BC on pressure = natural zero-flux)
- z_max (top): p = 0 (drained), u_z = -0.001 (or applied load)
- x/y boundaries: roller (constrain normal displacement component only)

Modify setupBoundaryConditions() to handle the POROELASTIC case differently:

```cpp
if (config.solid_model == SolidModelType::POROELASTIC &&
    config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
    // Poroelasticity: field 0 = pressure (1 comp), field 1 = displacement (3 comp)

    // Bottom: fixed displacement (all components)
    if (label_zmin) {
        PetscInt comps_all[3] = {0, 1, 2};
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_bottom",
                             label_zmin, 1, &label_value,
                             1, 3, comps_all,  // field 1 = displacement
                             (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }

    // Top: drained (p = 0)
    if (label_zmax) {
        PetscInt comp = 0;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "drained_top",
                             label_zmax, 1, &label_value,
                             0, 1, &comp,      // field 0 = pressure
                             (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }

    // Top: applied compression displacement
    if (label_zmax) {
        PetscInt comps_all[3] = {0, 1, 2};
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "compression_top",
                             label_zmax, 1, &label_value,
                             1, 3, comps_all,  // field 1 = displacement
                             (void (*)(void))bc_compression, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }

    // Lateral rollers: constrain normal displacement only
    // x_min: u_x = 0 (component 0 of displacement field 1)
    if (label_xmin) {
        PetscInt comp = 0;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmin",
                             label_xmin, 1, &label_value,
                             1, 1, &comp,
                             (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }
    // x_max: u_x = 0
    if (label_xmax) {
        PetscInt comp = 0;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmax",
                             label_xmax, 1, &label_value,
                             1, 1, &comp,
                             (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }
    // y_min: u_y = 0 (component 1 of displacement field 1)
    if (label_ymin) {
        PetscInt comp = 1;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymin",
                             label_ymin, 1, &label_value,
                             1, 1, &comp,
                             (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }
    // y_max: u_y = 0
    if (label_ymax) {
        PetscInt comp = 1;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymax",
                             label_ymax, 1, &label_value,
                             1, 1, &comp,
                             (void (*)(void))bc_zero, nullptr, nullptr, nullptr);
        CHKERRQ(ierr);
    }
}
```

### 2b: Set Terzaghi initial conditions

For Terzaghi: at t=0, the applied load generates undrained pore pressure
p0 = B * sigma_applied, where B is Skempton's coefficient. For fully
coupled (alpha=1, incompressible grains): p0 approximately equals the
applied load.

Add a Terzaghi-specific initial condition. In setInitialConditions(),
after VecZeroEntries:

```cpp
// For poroelasticity: set initial pore pressure = applied load (undrained)
if (config.solid_model == SolidModelType::POROELASTIC) {
    // Get local vector, set pressure DOFs to initial pressure
    Vec localVec;
    ierr = DMGetLocalVector(dm, &localVec); CHKERRQ(ierr);
    ierr = VecZeroEntries(localVec); CHKERRQ(ierr);

    // Set pressure field (field 0) to initial value everywhere
    // Using DMProjectFunctionLocal for field-specific projection
    PetscErrorCode (*pressureIC)(PetscInt, PetscReal, const PetscReal[],
                                  PetscInt, PetscScalar*, void*) = nullptr;
    // For Terzaghi: p0 = applied_load (1 MPa)
    pressureIC = [](PetscInt dim, PetscReal time, const PetscReal x[],
                     PetscInt Nc, PetscScalar *u, void *ctx) -> PetscErrorCode {
        (void)dim; (void)time; (void)x; (void)Nc; (void)ctx;
        u[0] = 1.0e6;  // 1 MPa initial pore pressure
        return PETSC_SUCCESS;
    };

    // Project pressure IC onto field 0 only
    // NOTE: Check PETSc 3.22 API for DMProjectFieldLocal or DMProjectFunctionLocal
    // Alternative: manually iterate over pressure DOFs and set them

    ierr = DMRestoreLocalVector(dm, &localVec); CHKERRQ(ierr);
}
```

If DMProjectFunctionLocal is too complex, the simpler approach is to
iterate over DOFs using the PetscSection and set pressure DOFs manually.

### 2c: Create Terzaghi config

Create `config/examples/terzaghi_consolidation.config`:
```ini
# Terzaghi 1D Consolidation
# Column loaded on top, drained at top, impermeable sides/bottom
# Pressure dissipates over time as fluid drains upward
# Validates: PetscFEPoroelasticity Biot coupling

[SIMULATION]
name = terzaghi_consolidation
start_time = 0.0
end_time = 10.0
dt_initial = 0.1
dt_min = 0.001
dt_max = 1.0
max_timesteps = 100
output_frequency = 10
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = false
rtol = 1.0e-8
atol = 1.0e-10
max_nonlinear_iterations = 20

[GRID]
nx = 1
ny = 1
nz = 20
Lx = 1.0
Ly = 1.0
Lz = 10.0

[ROCK]
density = 2500.0
youngs_modulus = 1.0e9
poissons_ratio = 0.25
porosity = 0.3
permeability_x = 1.0
permeability_y = 1.0
permeability_z = 1.0
biot_coefficient = 1.0

[FLUID]
type = SINGLE_PHASE
density = 1000.0
viscosity = 0.001
compressibility = 4.5e-10
reference_pressure = 0.0
```

Run and verify:
```bash
./fsrm -c ../config/examples/terzaghi_consolidation.config \
  -snes_monitor -ts_monitor -pc_type lu 2>&1 | tail -40
```

Expected: pressure should decrease over time as fluid drains. The TS monitor
should show multiple timesteps progressing. SNES should converge at each step.


## Phase 3: Injection Point Source

Add a wellbore injection term to the pressure equation. This is needed for
hydraulic fracturing (pressurize a point in the formation).

### Approach: Add source term to FormFunction

Similar to the explosion source pattern. After DMPlexTSComputeIFunctionFEM
assembles the bulk residual, add a point source contribution.

Add to Simulator.hpp (private):
```cpp
PetscErrorCode addInjectionSourceToResidual(PetscReal t, Vec F);
```

Implementation: find the cell containing the injection point (DMLocatePoints),
get closure DOFs, add injection rate to the pressure DOF residual.

For the pressure equation, the injection source is:
  residual -= Q / V_cell  (mass injection rate per unit volume)

where Q is the volumetric injection rate (m^3/s) and V_cell is the cell volume.

This is simpler than the explosion moment tensor injection because it only
touches one scalar DOF (pressure) per node in the source cell.

Parse injection parameters from config:
```ini
[INJECTION]
enabled = true
x = 0.5
y = 0.5
z = 5.0
rate = 1.0e-4    # m^3/s
start_time = 0.0
end_time = 100.0
```

### Create injection example

Create `config/examples/injection_pressure_buildup.config`:
- Poroelastic domain, 10x10x10
- Single injection point at center
- All boundaries: roller + impermeable (closed system)
- Watch pressure build up around the injection point

This is the precursor to hydraulic fracturing: inject fluid, watch pressure
rise, and when it exceeds fracture gradient, trigger fracture propagation.


## Phase 4: Connect HydraulicFractureModel

In MonitorFunction (or a new post-timestep callback):

1. Get pressure at the injection point from the solution vector
2. Evaluate K_I using LEFM::computeStressIntensityFactor(pressure, crack_length, min_stress)
3. If K_I > K_Ic, call HydraulicFractureModel::propagate(pressure, dt)
4. Update fracture geometry (length, width, height)
5. Modify permeability in the fracture zone:
   - Identify cells along the fracture path
   - Increase their permeability (via auxiliary field or modified constants)
   - Fracture permeability from cubic law: k_f = w^2/12

The permeability modification is the hard part. Options:
a) Use PETSc auxiliary fields (DMAux) to store per-cell permeability
b) Modify the PetscDS constants at each timestep (only works for uniform)
c) Use a callback that reads fracture geometry and returns permeability

Option (a) is the PETSc way. Options (b) and (c) are simpler for a first pass.

Create `config/examples/hydraulic_fracture_pkn.config`:
- Poroelastic domain
- Injection at center
- PKN fracture model
- Watch fracture grow and pressure plateau

### Subagent breakdown for this session

Subagent A: Verify Phase 1 (build, run uniaxial compression, confirm nonzero)
Subagent B: Implement Phase 2 (poroelastic BCs, Terzaghi IC, config)
Subagent C: Implement Phase 3 (injection source, config)

Main agent: Phase 4 (HydraulicFractureModel coupling) after subagents complete.


## Do NOT
- Modify PetscFE callback files
- Modify FaultMeshManager or CohesiveFaultKernel
- Delete existing tests
- Skip Phase 1 verification