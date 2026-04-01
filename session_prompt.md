# FSRM: Boundary Conditions + Working Examples

Read CLAUDE.md first. Build and test in Docker.

## The Problem

Simulator::setupBoundaryConditions() is empty. No Dirichlet or Neumann BCs are
applied to the PetscDS. Without BCs, elastostatic systems are singular (rigid
body modes) and no config produces a physical result. setInitialConditions()
blindly sets VecSet(solution, 1e7) regardless of what physics are active.

## Goal

Implement boundary conditions via PetscDS, fix initial conditions, then create
3 working examples that produce verifiable physical output.

## Subagent A: Implement Boundary Labeling and Dirichlet BCs

### Step 1: Label boundary faces by coordinate

Add a new method Simulator::labelBoundaries() that creates DMLabels for each
face of the bounding box. Call it after setupDM()/setupFaultNetwork() but
BEFORE setupFields().

For a box mesh on [0,Lx] x [0,Ly] x [0,Lz], iterate over all faces (height 1
stratum), check support size == 1 (boundary), compute centroid, and assign to
the appropriate label based on which coordinate is at min/max.

Label names: "boundary_x_min", "boundary_x_max", "boundary_y_min",
"boundary_y_max", "boundary_z_min", "boundary_z_max"

After labeling faces, call DMPlexLabelComplete on each label to add vertices
and edges in the closure. PetscDSAddBoundary needs the complete closure.

Add the labelBoundaries() call to main.cpp between setupFaultNetwork() and
setupFields().

### Step 2: Register Dirichlet BCs in setupBoundaryConditions()

Replace the empty stub. The key PETSc function is PetscDSAddBoundary (or
DMAddBoundary in some versions).

IMPORTANT: Check the PETSc 3.22.2 API first:
```bash
grep -n "PetscDSAddBoundary\|DMAddBoundary" /opt/petsc-3.22.2/include/petscds.h /opt/petsc-3.22.2/include/petscdm.h
```

Also look at PETSc examples for the correct calling convention:
```bash
grep -rn "PetscDSAddBoundary\|DMAddBoundary" /opt/petsc-3.22.2/src/snes/tutorials/ex56.c /opt/petsc-3.22.2/src/snes/tutorials/ex17.c 2>/dev/null | head -20
```

The BC callback function signature for essential (Dirichlet) BCs:
```c
void bc_zero(PetscInt dim, PetscReal time, const PetscReal x[],
             PetscInt Nc, PetscScalar *u, void *ctx) {
    for (PetscInt c = 0; c < Nc; c++) u[c] = 0.0;
}
```

For elastostatics/elastodynamics, register:
- Fixed bottom (z_min): all displacement components = 0
- The top boundary forcing comes from the initial stress state or from
  a Neumann BC. For the first pass, use a nonzero Dirichlet BC on top
  (applied displacement) since it's simpler than Neumann.

For the top boundary (applied compression), create:
```c
void bc_compression(PetscInt dim, PetscReal time, const PetscReal x[],
                     PetscInt Nc, PetscScalar *u, void *ctx) {
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = -0.001;  // 1mm downward displacement
}
```

For poroelasticity, also register:
- Drained top (z_max): pressure = 0

For fluid flow only (no geomechanics):
- Fixed pressure on x_min: p = reference_pressure
- No-flow (natural BC) on all other faces

### Step 3: Fix initial conditions

Replace VecSet(solution, 1e7) with physics-appropriate initialization:

```cpp
PetscErrorCode Simulator::setInitialConditions() {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    // Start from zero (displacement at rest, zero pressure perturbation)
    ierr = VecZeroEntries(solution); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
```

Zero is appropriate because:
- Displacement: body starts at rest
- Pressure: for diffusion problems, IC is set via config or as perturbation
- The Dirichlet BCs on the boundary will drive the solution


## Subagent B: Example 1 - Uniaxial Compression (Elastostatics)

Create `config/examples/uniaxial_compression.config`:

```ini
# Uniaxial Compression: Applied displacement on top, fixed bottom
# Expected: uniform strain eps_zz, linear displacement u_z(z)
# Validates: PetscFEElasticity callbacks, Dirichlet BCs, SNES convergence

[SIMULATION]
name = uniaxial_compression
start_time = 0.0
end_time = 1.0
dt_initial = 1.0
output_frequency = 1
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true
enable_faults = false
rtol = 1.0e-10
atol = 1.0e-12
max_nonlinear_iterations = 1

[GRID]
nx = 4
ny = 4
nz = 8
Lx = 10.0
Ly = 10.0
Lz = 100.0

[ROCK]
density = 2650.0
youngs_modulus = 10.0e9
poissons_ratio = 0.25
```

This problem has an exact solution: u_z(z) = -0.001 * z/Lz (linear),
eps_zz = -1e-5, sigma_zz = (lambda + 2*mu) * eps_zz.

The SNES should converge in 1 iteration (linear problem).


## Subagent C: Example 2 - Locked Fault Under Compression

Create `config/examples/fault_compression.config`:

```ini
# Locked fault under uniaxial compression
# Expected: continuous displacement across fault (locked = no slip)
# Validates: Cohesive cell insertion, Lagrange multiplier constraint

[SIMULATION]
name = fault_compression
start_time = 0.0
end_time = 1.0
dt_initial = 1.0
output_frequency = 1
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true
enable_faults = true
rtol = 1.0e-10
atol = 1.0e-12
max_nonlinear_iterations = 10

[GRID]
nx = 4
ny = 4
nz = 8
Lx = 10.0
Ly = 10.0
Lz = 100.0

[ROCK]
density = 2650.0
youngs_modulus = 10.0e9
poissons_ratio = 0.25
```

Note: This uses simplex meshes (forced when enable_faults=true). The fault
plane is created by setupFaultNetwork() at the domain center.


## Main Agent: After Subagents Complete

### 1. Wire labelBoundaries() into the pipeline

In main.cpp, after setupFaultNetwork() and before setupFields():
```cpp
ierr = sim.labelBoundaries(); CHKERRQ(ierr);
```

Add the declaration to Simulator.hpp (public).

### 2. Wire setupBoundaryConditions() 

It's already called in the pipeline (check main.cpp). Just make sure it's
called AFTER setupPhysics() (needs the PetscDS prob to be created).

### 3. Clean up config directory

```bash
mkdir -p config/examples
mkdir -p config/aspirational
# Move new working examples
# Keep test configs in config/
# Move everything else to aspirational
git mv config/atmospheric_nuclear_test*.config config/aspirational/
git mv config/buckley_leverett*.config config/aspirational/
git mv config/cascadia*.config config/aspirational/
# ... etc for all non-test, non-example configs
```

### 4. Create config/examples/README.md

```markdown
# FSRM Working Examples

These examples exercise verified code paths and produce physical results.

## uniaxial_compression.config
Quasi-static uniaxial compression of an elastic cube.
Fixed bottom, applied 1mm displacement on top.
Validates elastostatic PetscFE callbacks and Dirichlet BCs.

## fault_compression.config  
Same as uniaxial compression but with a locked cohesive fault
at the domain center. Validates fault mesh splitting and
Lagrange multiplier constraints.

## Running
```
cd build
./fsrm -c ../config/examples/uniaxial_compression.config
./fsrm -c ../config/examples/fault_compression.config
```
```

### 5. Build, run, verify

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON &&
  make -j$(nproc) &&
  echo "=== ALL TESTS ===" &&
  ctest --output-on-failure 2>&1 | tail -20 &&
  echo "=== UNIAXIAL COMPRESSION ===" &&
  ./fsrm -c ../config/examples/uniaxial_compression.config -snes_monitor -snes_converged_reason -ksp_converged_reason 2>&1 | tail -30 &&
  echo "=== FAULT COMPRESSION ===" &&
  ./fsrm -c ../config/examples/fault_compression.config -snes_monitor -snes_converged_reason -ksp_converged_reason 2>&1 | tail -30
'
```

Expected output for uniaxial compression:
```
SNES Function norm 1.234e-05
  Linear solve converged due to CONVERGED_RTOL
SNES converged due to CONVERGED_FNORM_RELATIVE
Simulation completed successfully!
```

If SNES doesn't converge:
- Add -pc_type lu for direct solve (eliminates preconditioning issues)
- Check that boundary labels have nonzero point counts
- Print DMLabelView for boundary labels to verify
- Check that PetscDSAddBoundary is called with correct field indices


## Do NOT
- Modify PetscFE callback files
- Modify FaultMeshManager or CohesiveFaultKernel
- Delete existing tests
- Create examples for features that don't work (DG, ADER, GPU, plasticity, thermal)