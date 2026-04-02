# FSRM: Poroelasticity, Hydraulic Fracturing, Cohesive Fracture Propagation, and Nuclear Explosion Modeling

Read CLAUDE.md first. Build and test everything in Docker. Execute phases sequentially.
Each phase must build and pass all tests before proceeding to the next.
Use subagents for parallel independent file creation within each phase.

## ============================================================================
## PHASE 1: Terzaghi Consolidation (Poroelasticity End-to-End)
## ============================================================================

The poroelastic callbacks (PetscFEPoroelasticity) are unit tested but have never
run end-to-end in the Simulator. The BCs are hardcoded for elastostatics only.

### 1a: Add poroelastic boundary conditions

In src/core/Simulator.cpp, modify setupBoundaryConditions() to handle the
POROELASTIC + SINGLE_COMPONENT case separately from pure elastostatics.

For Terzaghi consolidation (field 0 = pressure, field 1 = displacement):
- z_min (bottom): displacement fixed (u=0, all 3 components of field 1), impermeable (natural BC on pressure, no DMAddBoundary needed)
- z_max (top): drained (p=0, field 0 component 0), displacement fixed to applied compression (field 1, all 3 components)
- x_min/x_max: roller (constrain u_x only, component 0 of field 1)
- y_min/y_max: roller (constrain u_y only, component 1 of field 1)

The existing code already has a `displacement_field >= 0` block for elastostatics
and a separate poroelastic drained-top block. Restructure to:

```cpp
if (config.solid_model == SolidModelType::POROELASTIC &&
    config.fluid_model == FluidModelType::SINGLE_COMPONENT) {
    // POROELASTIC case: field 0 = pressure, field 1 = displacement
    PetscInt pfield = 0;
    PetscInt ufield = 1;
    PetscInt label_value = 1;

    // Bottom: fix all displacement components
    if (label_zmin) {
        PetscInt comps[3] = {0, 1, 2};
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "fixed_bottom",
            label_zmin, 1, &label_value, ufield, 3, comps,
            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
    // Top: drained (p = 0)
    if (label_zmax) {
        PetscInt comp = 0;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "drained_top",
            label_zmax, 1, &label_value, pfield, 1, &comp,
            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
    // Top: applied compression (u_z = -0.001)
    if (label_zmax) {
        PetscInt comps[3] = {0, 1, 2};
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "compression_top",
            label_zmax, 1, &label_value, ufield, 3, comps,
            (void (*)(void))bc_compression, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
    // Lateral rollers
    if (label_xmin) {
        PetscInt comp = 0;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmin",
            label_xmin, 1, &label_value, ufield, 1, &comp,
            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
    if (label_xmax) {
        PetscInt comp = 0;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_xmax",
            label_xmax, 1, &label_value, ufield, 1, &comp,
            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
    if (label_ymin) {
        PetscInt comp = 1;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymin",
            label_ymin, 1, &label_value, ufield, 1, &comp,
            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
    if (label_ymax) {
        PetscInt comp = 1;
        ierr = DMAddBoundary(dm, DM_BC_ESSENTIAL, "roller_ymax",
            label_ymax, 1, &label_value, ufield, 1, &comp,
            (void (*)(void))bc_zero, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    }
} else if (displacement_field >= 0) {
    // Pure elastostatics: existing code (field 0 = displacement)
    // ... keep existing elastostatics BC block unchanged ...
}
```

### 1b: Set poroelastic initial conditions

In setInitialConditions(), after VecZeroEntries, add initial pore pressure
for the Terzaghi problem. The undrained response to sudden loading produces
initial pore pressure p0 = applied_load (for alpha = 1).

The simplest approach: iterate over pressure DOFs using the PetscSection and
set them to the initial pressure value. Check how many DOFs field 0 has at
each point, and set them.

Alternatively, use DMProjectFunction with separate functions per field:
```cpp
if (config.solid_model == SolidModelType::POROELASTIC) {
    // Set initial pressure = 1 MPa everywhere (undrained response)
    // DMProjectFunction sets all fields; we only want pressure (field 0)
    // Use DMProjectFieldLocal for field-specific projection if available
    // Otherwise, iterate over section DOFs manually
}
```

For now, a simple approach: create a callback that returns 1e6 for pressure
and 0 for displacement, and use DMProjectFunction:
```cpp
static PetscErrorCode terzaghi_ic(PetscInt dim, PetscReal time, const PetscReal x[],
                                   PetscInt Nc, PetscScalar *u, void *ctx) {
    (void)dim; (void)time; (void)x; (void)ctx;
    // Nc=1 for pressure field, Nc=3 for displacement field
    if (Nc == 1) {
        u[0] = 1.0e6;  // 1 MPa initial pore pressure
    } else {
        for (PetscInt c = 0; c < Nc; c++) u[c] = 0.0;
    }
    return PETSC_SUCCESS;
}
```

### 1c: Create Terzaghi config

Create config/examples/terzaghi_consolidation.config:
```ini
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
nx = 2
ny = 2
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

### 1d: Verify

```bash
./fsrm -c ../config/examples/terzaghi_consolidation.config \
  -ts_monitor -snes_monitor -snes_converged_reason -pc_type lu 2>&1 | tail -50
```

Expected: TS advances through timesteps, SNES converges at each step.
Pressure should dissipate over time (the pore pressure at early times
is high, and decreases as fluid drains through the top).

If SNES doesn't converge:
- Check that constants are set correctly for poroelastic case (biot_alpha at [22], 1/M at [23])
- Check field ordering: pressure must be field 0, displacement field 1
- Try -pc_type lu -ksp_type preonly for direct solve
- Add -snes_fd to use finite difference Jacobian (bypasses analytical Jacobian bugs)


## ============================================================================
## PHASE 2: Injection Point Source
## ============================================================================

Add a volumetric injection source to the pressure equation so we can pump
fluid into the formation. This is the wellbore for hydraulic fracturing.

### 2a: Parse injection config

In initializeFromConfigFile(), parse an [INJECTION] section:
```cpp
if (reader.hasSection("INJECTION")) {
    injection_enabled_ = reader.getBool("INJECTION", "enabled", true);
    injection_x_ = reader.getDouble("INJECTION", "x", 0.0);
    injection_y_ = reader.getDouble("INJECTION", "y", 0.0);
    injection_z_ = reader.getDouble("INJECTION", "z", 0.0);
    injection_rate_ = reader.getDouble("INJECTION", "rate", 1e-4);  // m^3/s
    injection_start_ = reader.getDouble("INJECTION", "start_time", 0.0);
    injection_end_ = reader.getDouble("INJECTION", "end_time", 1e30);
}
```

Add these as members in Simulator.hpp:
```cpp
bool injection_enabled_ = false;
double injection_x_ = 0, injection_y_ = 0, injection_z_ = 0;
double injection_rate_ = 0;
double injection_start_ = 0, injection_end_ = 1e30;
PetscInt injection_cell_ = -1;  // cached cell index
```

### 2b: Locate injection cell

After setupDM(), find the cell containing the injection point:
```cpp
PetscErrorCode Simulator::locateInjectionCell() {
    if (!injection_enabled_) return 0;
    PetscErrorCode ierr;

    PetscReal point[3] = {injection_x_, injection_y_, injection_z_};
    Vec pointVec;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 3, &pointVec); CHKERRQ(ierr);
    PetscScalar *pv;
    ierr = VecGetArray(pointVec, &pv); CHKERRQ(ierr);
    pv[0] = point[0]; pv[1] = point[1]; pv[2] = point[2];
    ierr = VecRestoreArray(pointVec, &pv); CHKERRQ(ierr);

    IS cellsIS;
    ierr = DMLocatePoints(dm, pointVec, DM_POINTLOCATION_NONE, &cellsIS); CHKERRQ(ierr);
    const PetscInt *cells;
    ierr = ISGetIndices(cellsIS, &cells); CHKERRQ(ierr);
    injection_cell_ = cells[0];
    ierr = ISRestoreIndices(cellsIS, &cells); CHKERRQ(ierr);
    ierr = ISDestroy(&cellsIS); CHKERRQ(ierr);
    ierr = VecDestroy(&pointVec); CHKERRQ(ierr);

    if (rank == 0) {
        PetscPrintf(comm, "Injection point located in cell %d\n", injection_cell_);
    }
    return 0;
}
```

Call this after setupDM(). Check the PETSc 3.22 DMLocatePoints signature first.

### 2c: Add source to FormFunction

In FormFunction, after DMPlexTSComputeIFunctionFEM, add the injection:
```cpp
if (sim->injection_enabled_ && sim->injection_cell_ >= 0) {
    PetscReal t_local = static_cast<double>(t);
    if (t_local >= sim->injection_start_ && t_local <= sim->injection_end_) {
        ierr = sim->addInjectionToResidual(t, F); CHKERRQ(ierr);
    }
}
```

The injection modifies the pressure equation residual. For the pressure DOFs
in the injection cell, subtract Q/V_cell (positive injection = source):

```cpp
PetscErrorCode Simulator::addInjectionToResidual(PetscReal t, Vec F) {
    PetscFunctionBeginUser;
    if (injection_cell_ < 0) PetscFunctionReturn(0);
    PetscErrorCode ierr;

    // Compute cell volume
    PetscReal vol;
    ierr = DMPlexComputeCellGeometryFVM(dm, injection_cell_, &vol, nullptr, nullptr);
    CHKERRQ(ierr);

    // Get closure indices for the injection cell
    PetscSection section;
    ierr = DMGetLocalSection(dm, &section); CHKERRQ(ierr);

    // Pressure field is field 0 for poroelastic, field 0 for single-phase
    PetscInt pressure_field = 0;
    PetscInt pStart, pEnd;
    ierr = DMPlexGetClosureIndices(dm, section, section, injection_cell_,
                                    PETSC_TRUE, &pStart, nullptr, nullptr, nullptr);
    // ... this API is complex. Alternative: use DMPlexVecGetClosure/SetClosure

    // Simpler approach: add source via VecSetValue on the global residual
    // Find the global DOF for pressure at the injection cell
    PetscInt offset;
    ierr = PetscSectionGetOffset(section, injection_cell_, &offset); CHKERRQ(ierr);
    // offset gives the start of DOFs for this cell in the local vector

    // For now, use a simpler approach: modify the LOCAL residual in FormFunction
    // before DMLocalToGlobal
    PetscScalar source = -injection_rate_ / vol;  // negative = source in IFunction
    // Add to pressure DOF (first DOF in the cell for poroelastic)
    // ... need to figure out the right DOF index ...

    PetscFunctionReturn(0);
}
```

NOTE: Getting the right DOF index for a specific field at a specific cell is
nontrivial with PetscSection. The correct approach for adding a point source
to the FEM residual is to use the closure:

```cpp
// Get closure DOFs for the cell
PetscScalar *closure = nullptr;
PetscInt closureSize;
ierr = DMPlexVecGetClosure(dm, section, locF, injection_cell_, &closureSize, &closure);
CHKERRQ(ierr);

// The closure array has DOFs ordered by field then by basis function
// For poroelastic: pressure DOFs come first, then displacement DOFs
// Modify the pressure DOFs
PetscInt nPressureDofs = ...; // number of pressure basis functions in this cell
for (PetscInt i = 0; i < nPressureDofs; i++) {
    closure[i] += source / nPressureDofs;
}

ierr = DMPlexVecSetClosure(dm, section, locF, injection_cell_, closure, ADD_VALUES);
CHKERRQ(ierr);
ierr = DMPlexVecRestoreClosure(dm, section, locF, injection_cell_, &closureSize, &closure);
CHKERRQ(ierr);
```

Check the PETSc 3.22 API for DMPlexVecGetClosure/SetClosure.

### 2d: Create injection config and verify

Create config/examples/injection_pressure_buildup.config:
```ini
[SIMULATION]
name = injection_pressure_buildup
start_time = 0.0
end_time = 100.0
dt_initial = 1.0
dt_min = 0.01
dt_max = 10.0
max_timesteps = 100
output_frequency = 10
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_faults = false
rtol = 1.0e-6
atol = 1.0e-8
max_nonlinear_iterations = 20

[GRID]
nx = 10
ny = 10
nz = 10
Lx = 1000.0
Ly = 1000.0
Lz = 1000.0

[ROCK]
density = 2500.0
youngs_modulus = 10.0e9
poissons_ratio = 0.25
porosity = 0.2
permeability_x = 10.0
permeability_y = 10.0
permeability_z = 10.0
biot_coefficient = 0.8

[FLUID]
type = SINGLE_PHASE
density = 1000.0
viscosity = 0.001
compressibility = 4.5e-10
reference_pressure = 30.0e6

[INJECTION]
enabled = true
x = 500.0
y = 500.0
z = 500.0
rate = 1.0e-3
start_time = 0.0
end_time = 100.0
```


## ============================================================================
## PHASE 3: Hydraulic Fracture with PKN/KGD Model
## ============================================================================

Connect the existing HydraulicFractureModel to the solver. At each timestep,
evaluate pressure at the injection point, check fracture initiation, and
propagate the fracture if K_I > K_Ic.

### 3a: Wire fracture model to Simulator

In Simulator.hpp, add:
```cpp
std::unique_ptr<HydraulicFractureModel> hydrofrac_;
bool hydrofrac_initiated_ = false;
```

In initializeFromConfigFile(), parse [HYDRAULIC_FRACTURE] section:
```cpp
if (reader.hasSection("HYDRAULIC_FRACTURE")) {
    hydrofrac_ = std::make_unique<HydraulicFractureModel>();
    hydrofrac_->setFormationProperties(
        reader.getDouble("HYDRAULIC_FRACTURE", "youngs_modulus", 10e9),
        reader.getDouble("HYDRAULIC_FRACTURE", "poissons_ratio", 0.25));
    hydrofrac_->setFractureToughness(
        reader.getDouble("HYDRAULIC_FRACTURE", "toughness", 1.5e6));  // Pa*m^0.5
    hydrofrac_->setStressState(
        reader.getDouble("HYDRAULIC_FRACTURE", "min_horizontal_stress", 20e6),
        reader.getDouble("HYDRAULIC_FRACTURE", "max_horizontal_stress", 25e6),
        reader.getDouble("HYDRAULIC_FRACTURE", "vertical_stress", 30e6));
    hydrofrac_->setHeight(
        reader.getDouble("HYDRAULIC_FRACTURE", "height", 50.0));
    hydrofrac_->setFluidProperties(
        reader.getDouble("HYDRAULIC_FRACTURE", "fluid_density", 1000.0),
        reader.getDouble("HYDRAULIC_FRACTURE", "fluid_viscosity", 0.001));
    std::string model = reader.getString("HYDRAULIC_FRACTURE", "model", "PKN");
    if (model == "KGD") hydrofrac_->setModel("KGD");
    else if (model == "P3D") hydrofrac_->setModel("P3D");
    else hydrofrac_->setModel("PKN");
}
```

### 3b: Update fracture in MonitorFunction

In MonitorFunction, after seismometer sampling:
```cpp
if (sim->hydrofrac_ && sim->injection_enabled_) {
    // Get pressure at injection point from solution
    double wellbore_pressure = 0.0;
    // Use DMInterpolationEvaluate or closure to extract pressure at injection point
    // ... (similar to seismometer sampling) ...

    // Check initiation
    double K_I = FSRM::LEFM::computeStressIntensityFactor(
        wellbore_pressure, sim->hydrofrac_->getLength(),
        sim->hydrofrac_->getMinHorizontalStress());

    if (K_I > sim->hydrofrac_->getFractureToughness()) {
        if (!sim->hydrofrac_initiated_) {
            sim->hydrofrac_initiated_ = true;
            if (sim->rank == 0) {
                PetscPrintf(sim->comm,
                    "FRACTURE INITIATED at t=%.2f, K_I=%.2e > K_Ic=%.2e\n",
                    (double)t, K_I, sim->hydrofrac_->getFractureToughness());
            }
        }
        // Propagate
        sim->hydrofrac_->propagate(wellbore_pressure, sim->dt);

        // Modify permeability in fractured cells
        // The fracture extends from the injection point in the x-direction
        // (perpendicular to min horizontal stress)
        // For each cell within the fracture length, increase permeability
        // using cubic law: k_f = w^2 / 12
        double frac_perm = sim->hydrofrac_->getWidth() * sim->hydrofrac_->getWidth() / 12.0;
        // TODO: Update PetscDS constants or auxiliary field for affected cells
        // For a first pass, print the fracture state:
        if (sim->rank == 0 && ((int)step % 10 == 0)) {
            PetscPrintf(sim->comm,
                "FRACTURE: L=%.2f m, W=%.4e m, H=%.1f m, k_f=%.2e m^2\n",
                sim->hydrofrac_->getLength(), sim->hydrofrac_->getWidth(),
                sim->hydrofrac_->getHeight(), frac_perm);
        }
    }
}
```

### 3c: Create hydraulic fracture config

Create config/examples/hydraulic_fracture_pkn.config with the injection
section plus [HYDRAULIC_FRACTURE] section. Run and verify that fracture
initiates and grows, printing length/width at each step.


## ============================================================================
## PHASE 4: Cohesive Fracture Propagation
## ============================================================================

Use the working PyLith-style cohesive cell machinery to model fractures as
actual mesh discontinuities. The fracture plane is pre-defined (like a fault)
but starts locked and unlocks when the tensile stress exceeds strength.

### 4a: Tensile failure criterion in CohesiveFaultKernel

The CohesiveFaultKernel currently supports locked/slipping modes for shear
faults. Add a TENSILE mode:

In include/physics/CohesiveFaultKernel.hpp, add:
```cpp
static constexpr PetscInt COHESIVE_CONST_TENSILE_STRENGTH = 27;  // expand array to 28
```

Modify f0_lagrange_constraint to check for tensile failure:
When the normal traction on the cohesive surface exceeds the tensile strength,
switch from locked to opening mode. The Lagrange multiplier constraint becomes:
- Locked: [u+] - [u-] = 0 (displacement jump = 0)
- Tensile opening: lambda_n = 0 (normal traction = 0, free to open)

This models a pre-existing weakness plane that opens when pressurized by
injection. The fluid pressure on the fracture face drives the opening.

### 4b: Pre-define fracture plane

In setupFaultNetwork(), when the config specifies a fracture (not a tectonic
fault), create the cohesive cells along the planned fracture path. The fracture
starts locked and opens dynamically based on pressure.

Add a config option:
```ini
[FRACTURE_PLANE]
enabled = true
strike = 0.0
dip = 90.0
center_x = 500.0
center_y = 500.0
center_z = 500.0
length = 200.0
width = 100.0
tensile_strength = 5.0e6
```

### 4c: Create cohesive fracture example

Create config/examples/cohesive_hydraulic_fracture.config that uses the
injection source + pre-defined fracture plane. The fracture plane starts
locked. As injection raises pressure, the fracture opens and fluid flows
into the newly created aperture.


## ============================================================================
## PHASE 5: Explosion Source Physics and Seismograms
## ============================================================================

Fix the Mueller-Murphy source physics and inject the moment tensor into the
FEM residual to generate synthetic seismograms from underground explosions.

### 5a: Fix Mueller-Murphy corner frequency

In src/domain/explosion/ExplosionImpactPhysics.cpp,
MuellerMurphySource::computeDerivedQuantities():

Replace:
```cpp
corner_frequency = std::max(0.1, s_velocity / (3.0 * depth));
```

With the actual Mueller-Murphy (1971) psi function:
```cpp
const double W = std::max(1e-6, source_params.yield_kt);
const double depth = std::max(1.0, std::abs(source_params.depth_of_burial));
double psi;
if (W <= 10.0) {
    psi = 16.2 / std::pow(W, 0.33);
} else if (W <= 1000.0) {
    psi = 10.4 / std::pow(W, 0.25);
} else {
    psi = 5.85 / std::pow(W, 0.17);
}
corner_frequency = psi * std::pow(density / 2650.0, -1.0/3.0) / depth;
```

### 5b: Fix mb-yield relation

In NuclearSourceParameters::body_wave_magnitude():
Replace `4.5 + 0.80 * log10(W)` with `4.45 + 0.75 * log10(W)`.

For DPRK 2017 (250 kt): mb = 4.45 + 0.75*2.398 = 6.25 (observed: 6.3).

### 5c: Fix scalar moment from cavity mechanics

Replace the rough scaling with the cavity mechanics formula:
```cpp
double NuclearSourceParameters::scalar_moment() const {
    const double W = std::max(1e-12, yield_kt);
    double Rc = cavity_radius(2700.0);  // granite
    double rho = 2700.0;
    double vp = 5500.0;
    // M0 = 4*pi*rho*vp^2*Rc^3
    return 4.0 * M_PI * rho * vp * vp * Rc * Rc * Rc;
}
```

### 5d: Inject moment tensor into FEM residual

Add Simulator::addExplosionSourceToResidual() that injects the moment tensor
as equivalent nodal forces. The explosion source produces an isotropic moment
tensor: M_ij = M0(t)/3 * delta_ij.

The source injection follows the same pattern as the injection source (Phase 2)
but acts on displacement DOFs instead of pressure DOFs. Use DMLocatePoints to
find the source cell, then add forces via DMPlexVecGetClosure/SetClosure.

In FormFunction, after the injection source:
```cpp
if (sim->explosion_ && sim->use_fem_time_residual_) {
    ierr = sim->addExplosionSourceToResidual(t, locF); CHKERRQ(ierr);
}
```

### 5e: Create explosion + seismogram config

Create config/examples/explosion_seismogram.config:
- Elastodynamic domain (fluid_model = NONE, enable_elastodynamics = true)
- 10 kt explosion at 300m depth
- 3 surface seismometers at 1km, 3km, 5km
- TSCN time integration, 2000 timesteps, dt = 0.001
- Seismometer output in SAC format

IMPORTANT: Use ONLY working features. No DG, no ADER, no GPU, no plasticity.

### 5f: Add DPRK 2017 validation test

Add to tests/unit/domain/test_explosion_source.cpp:
```cpp
TEST_F(ExplosionSourceKernelTest, DPRK2017MagnitudeCheck) {
    NuclearSourceParameters params;
    params.yield_kt = 250.0;
    params.depth_of_burial = 800.0;
    double mb = params.body_wave_magnitude();
    EXPECT_NEAR(mb, 6.25, 0.15);
}
```


## ============================================================================
## PHASE 6: Nuclear Monitoring Validation Suite
## ============================================================================

Create a set of configs for historical nuclear tests that can be compared
against published seismic data.

### 6a: DPRK 2017 Punggye-ri quick config

Create config/examples/dprk_2017_quick.config:
- 250 kt, 800m depth in granite
- 20x20x20 mesh, 20km domain
- Elastodynamics (NOT DG/ADER -- use PetscFE with TSCN)
- 3 seismometers: near-field (5km), regional (50km placeholder), teleseismic (placeholder)
- End time: 5s (near-field only)

### 6b: Generic underground test template

Create config/examples/underground_explosion_template.config with comprehensive
comments explaining each parameter and how to customize for different tests.

### 6c: Update examples README

Update config/examples/README.md listing all working examples:
- uniaxial_compression.config (elastostatics)
- terzaghi_consolidation.config (poroelasticity)
- injection_pressure_buildup.config (injection)
- hydraulic_fracture_pkn.config (PKN fracture)
- cohesive_hydraulic_fracture.config (cohesive fracture)
- explosion_seismogram.config (explosion + seismograms)
- dprk_2017_quick.config (nuclear monitoring)


## ============================================================================
## VERIFICATION CHECKLIST
## ============================================================================

After all phases, run:
```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  cd build && make -j$(nproc) &&
  echo "=== TESTS ===" &&
  ctest --output-on-failure 2>&1 | tail -20 &&
  echo "=== UNIAXIAL ===" &&
  ./fsrm -c ../config/examples/uniaxial_compression.config -snes_monitor -pc_type lu 2>&1 | tail -10 &&
  echo "=== TERZAGHI ===" &&
  ./fsrm -c ../config/examples/terzaghi_consolidation.config -ts_monitor -snes_monitor -pc_type lu 2>&1 | tail -20 &&
  echo "=== INJECTION ===" &&
  ./fsrm -c ../config/examples/injection_pressure_buildup.config -ts_monitor -snes_monitor -pc_type lu 2>&1 | tail -20 &&
  echo "=== EXPLOSION ===" &&
  ./fsrm -c ../config/examples/explosion_seismogram.config -ts_monitor -snes_monitor -pc_type lu 2>&1 | tail -20
'
```

Each example must:
1. Not segfault
2. Show "Simulation completed successfully"
3. SNES converges at each timestep
4. Solution norms are nonzero

## Do NOT
- Modify PetscFE callback math
- Modify FaultMeshManager::splitMeshAlongFault
- Modify CohesiveFaultKernel::registerWithDS internals
- Modify friction law implementations
- Use DG, ADER, GPU, or plasticity features (they are stubs)
- Delete existing tests
- Change the unified constants layout [0]-[26]
- CHANGE THE ORDERING OF DMCreateDS / setupBoundaryConditions / DMSetLocalSection / DMSetUp
  IN setupFields(). The current order is:
    1. DMCreateDS(dm)                          -- PETSc 3.22 requires DS before DMAddBoundary
    2. setupBoundaryConditions()               -- calls DMAddBoundary (needs DS to exist)
    3. DMSetLocalSection(dm, nullptr)           -- clear cached section
    4. DMSetUp(dm)                             -- rebuild section WITH BC constraints
    5. DMGetDS(dm, &prob)                      -- get DS for later use
  This was debugged over multiple sessions. Reversing it causes either crashes
  (DMAddBoundary before DMCreateDS) or silent BC failure (no section rebuild).
  If you propose an edit that moves DMCreateDS, removes the section rebuild,
  or reorders these calls, you are introducing a known bug. Do not do it.