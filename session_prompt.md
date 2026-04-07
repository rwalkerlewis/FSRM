# session_prompt.md -- FSRM Nuclear Test Modeling Development Session

Read CLAUDE.md in the repo root before doing anything. It contains build commands, architecture, constraints, and rules. Violating those rules will break the codebase.

## Pre-Session Setup

Before starting any phase, confirm the baseline:

```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure
```

**Expected**: All 50 tests pass. If any test fails, stop and fix it before proceeding. Do not proceed with failures.

---

## Phase 1: Verify All Existing Example Configs

**Goal**: Run every example config in config/examples/ and confirm SNES convergence. Identify and fix anything that crashes.

### Step 1.1: List all example configs

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'find config/examples/ -name "*.config" | sort'
```

Record the list. Each config must be run.

### Step 1.2: Run each example config

For each config file found in Step 1.1, run:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local bash -c \
  './fsrm ../config/examples/<CONFIG_FILE> 2>&1'
```

Replace `<CONFIG_FILE>` with the actual filename.

### Step 1.3: Triage results

For each config, record one of:
- **PASS**: SNES converges, nonzero FNORM, no crash
- **EXPECTED_UNVERIFIED**: Config runs but produces physically questionable output (this is acceptable for now)
- **CRASH**: Segfault, assertion failure, or PETSc error. Must be fixed.
- **CONFIG_ERROR**: Config references missing features or bad parameters. Note and skip.

### Step 1.4: Fix crashes

For each CRASH result:
1. Read the PETSc error message carefully. It usually names the function and error code.
2. Check if the crash is in a callback (PetscFEElasticity, PetscFEPoroelasticity, PetscFEFluidFlow). If so, the issue is likely a mismatched field count or incorrect aOff index.
3. Check if the crash is in mesh setup (DMPlex functions). If so, the issue is likely label or cone size mismatches.
4. Fix the minimal amount of code needed. Do not refactor.
5. Rebuild and rerun all 50 tests after each fix.

### Verification Gate 1

- All 50 unit/integration tests pass
- Every example config in config/examples/ either runs to completion or is documented as CONFIG_ERROR with explanation
- No regressions

**Do not proceed to Phase 2 until this gate passes.**

---

## Phase 2: Heterogeneous Material Properties via Auxiliary Fields

**Goal**: Replace the single-constants-array material property system with PETSc auxiliary fields that support per-cell material properties. This is the most important architectural change.

### Background

Currently, PetscDS callbacks read lambda, mu, rho from `constants[0]`, `constants[1]`, `constants[2]`. To support layered geology, callbacks must instead read material properties from auxiliary fields that are defined per-cell.

In PETSc, auxiliary fields are attached to a separate DM (the "auxiliary DM") via `DMSetAuxiliaryVec` or `PetscDSSetConstants` combined with `DMGetLocalVector` on the aux DM. The callback signature already supports this: the `a[]` array contains auxiliary field values and `aOff[]` contains offsets into the auxiliary fields.

Reference: PETSc ex17.c shows this pattern. Also see PyLith's MaterialNew::createAuxiliaryField().

### Step 2.1: Create new callback files for auxiliary-field-aware material

Create two new files:

```
src/numerics/PetscFEElasticityAux.cpp
include/numerics/PetscFEElasticityAux.hpp
```

These files contain new callback functions that read material properties from auxiliary fields instead of from the constants array. The callback signatures are identical to the existing ones, but the body reads from `a[aOff[AUX_LAMBDA]]`, `a[aOff[AUX_MU]]`, `a[aOff[AUX_RHO]]`.

Define auxiliary field indices:

```cpp
// In PetscFEElasticityAux.hpp
enum AuxFieldIndex {
    AUX_LAMBDA = 0,
    AUX_MU     = 1,
    AUX_RHO    = 2,
    NUM_AUX_FIELDS = 3
};
```

Write the following callbacks in PetscFEElasticityAux.cpp:

```cpp
// Residual: body force (gravity) using aux rho
static void f0_elasticity_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f0[])
{
    const PetscScalar rho = a[aOff[AUX_RHO]];
    for (PetscInt d = 0; d < dim; ++d) f0[d] = 0.0;
    // Gravity in negative z-direction (last coordinate)
    // constants[0] = gravity magnitude (positive value, e.g. 9.81)
    // Set to 0.0 if no gravity is configured
    f0[dim - 1] = -rho * constants[0];
}

// Residual: stress divergence using aux lambda, mu
static void f1_elasticity_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar f1[])
{
    const PetscScalar lambda = a[aOff[AUX_LAMBDA]];
    const PetscScalar mu     = a[aOff[AUX_MU]];
    // NOTE: These aux callbacks are written for the elastostatics case where
    // displacement is field 0 (fluid_model = NONE). For poroelastic cases
    // where displacement is field 1, separate aux callbacks with the
    // appropriate field offset would be needed.
    const PetscInt off = uOff_x[0];  // gradient offset for displacement field 0
    // Compute strain from u_x, then stress = lambda*tr(eps)*I + 2*mu*eps
    // Store in f1 as sigma_{ij}
    PetscScalar trace = 0.0;
    for (PetscInt d = 0; d < dim; ++d) trace += u_x[off + d * dim + d];
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            f1[i * dim + j] = mu * (u_x[off + i * dim + j] + u_x[off + j * dim + i]);
            if (i == j) f1[i * dim + j] += lambda * trace;
        }
    }
}

// Jacobian: g3 block using aux lambda, mu
static void g3_elasticity_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants,
    const PetscScalar constants[], PetscScalar g3[])
{
    const PetscScalar lambda = a[aOff[AUX_LAMBDA]];
    const PetscScalar mu     = a[aOff[AUX_MU]];
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            for (PetscInt k = 0; k < dim; ++k) {
                for (PetscInt l = 0; l < dim; ++l) {
                    g3[((i*dim+j)*dim+k)*dim+l] =
                        lambda * (i == j ? 1 : 0) * (k == l ? 1 : 0)
                        + mu * ((i == k ? 1 : 0) * (j == l ? 1 : 0)
                              + (i == l ? 1 : 0) * (j == k ? 1 : 0));
                }
            }
        }
    }
}
```

**IMPORTANT**: Before writing these callbacks, verify the exact callback typedef signature in PETSc 3.22.2:

```bash
grep -n "PetscPointFunc\|PetscPointJac" /opt/petsc-3.22.2/include/petscds.h
```

Match the signature exactly.

### Step 2.2: Add auxiliary DM setup to Simulator

Add a new method to `Simulator`:

```cpp
// In Simulator.hpp, add:
PetscErrorCode setupAuxiliaryDM();
DM             auxDM_ = nullptr;
Vec            auxVec_ = nullptr;
```

In `Simulator.cpp`, implement `setupAuxiliaryDM()`:

1. Clone the primary DM topology: `DMClone(dm_, &auxDM_)`
2. Create PetscFE fields for lambda, mu, rho on the auxiliary DM. Each is a scalar field (1 component).

   First, check which FE creation function exists in PETSc 3.22.2:
   ```bash
   grep -n "PetscFECreateLagrange\|PetscFECreateDefault" /opt/petsc-3.22.2/include/petscfe.h
   ```
   The existing FSRM code uses `PetscFECreateLagrange`. Use the same function for the auxiliary DM fields.

   ```cpp
   PetscFE fe_aux;
   // Use whichever FE creation function the grep above confirms exists.
   // The existing FSRM code uses PetscFECreateLagrange.
   PetscFECreateLagrange(PETSC_COMM_SELF, dim, 1, isSimplex, 1, PETSC_DETERMINE, &fe_aux);
   DMSetField(auxDM_, AUX_LAMBDA, NULL, (PetscObject)fe_aux);
   // Repeat for AUX_MU, AUX_RHO (create separate FE objects or reuse)
   DMCreateDS(auxDM_);
   ```
3. Create a local vector on the auxiliary DM: `DMCreateLocalVector(auxDM_, &auxVec_)`
4. Populate auxVec_ with material property values.

**Check PETSc 3.22.2 API**:

```bash
grep -n "DMClone\|DMSetField\|DMCreateDS\|DMCreateLocalVector\|DMSetAuxiliaryVec\|PetscDSSetObjective" /opt/petsc-3.22.2/include/petscdm*.h
```

### Step 2.3: Populate auxiliary fields by depth

Create a function that assigns material properties based on cell centroid z-coordinate:

```cpp
PetscErrorCode Simulator::populateAuxFieldsByDepth()
{
    // For each cell in the mesh:
    //   1. Compute centroid via DMPlexComputeCellGeometryFVM
    //   2. Look up z-coordinate
    //   3. Assign lambda, mu, rho based on depth intervals
    //   4. Set values in auxVec_ via DMPlexVecSetClosure on auxDM_

    // Example layered model (Punggye-ri-like):
    //   z > -500m:   alluvium     (lambda=3.0e9,  mu=1.5e9,  rho=1800)
    //   z > -2000m:  weathered granite (lambda=15e9, mu=12e9, rho=2400)
    //   z <= -2000m: competent granite  (lambda=30e9, mu=25e9, rho=2650)
}
```

Verify DMPlexComputeCellGeometryFVM signature:

```bash
grep -n "DMPlexComputeCellGeometryFVM" /opt/petsc-3.22.2/include/petscdmplex.h
```

### Step 2.4: Attach auxiliary data to primary DM

After populating auxVec_:

```cpp
DMSetAuxiliaryVec(dm_, NULL, 0, 0, auxVec_);
```

Or the PETSc 3.22-appropriate variant. Check:

```bash
grep -n "DMSetAuxiliaryVec\|PetscDSSetDiscretization" /opt/petsc-3.22.2/include/petscdm*.h
```

### Step 2.5: Register auxiliary-aware callbacks

In `setupPhysics()`, when `use_heterogeneous_material` config flag is true, register the new callbacks from PetscFEElasticityAux instead of PetscFEElasticity:

```cpp
if (config_.use_heterogeneous_material) {
    PetscDSSetResidual(prob, 0, f0_elasticity_aux, f1_elasticity_aux);
    PetscDSSetJacobian(prob, 0, 0, NULL, NULL, NULL, g3_elasticity_aux);
} else {
    // existing callback registration
}
```

### Step 2.6: Add config option

In ConfigReader, add parsing for:

```ini
[MATERIAL]
heterogeneous = true

[LAYER_1]
z_top = 0.0
z_bottom = -500.0
lambda = 3.0e9
mu = 1.5e9
rho = 1800.0

[LAYER_2]
z_top = -500.0
z_bottom = -2000.0
lambda = 15.0e9
mu = 12.0e9
rho = 2400.0

[LAYER_3]
z_top = -2000.0
z_bottom = -10000.0
lambda = 30.0e9
mu = 25.0e9
rho = 2650.0
```

### Step 2.7: Update CMakeLists.txt

Add PetscFEElasticityAux.cpp to the build. Confirm it compiles.

### Step 2.8: Write integration test for layered elastostatics

Create `tests/integration/test_layered_elastostatics.cpp`:

- 10x10x10 box mesh, z from 0 to -5000m
- 3-layer model (alluvium / weathered granite / granite)
- Dirichlet BC: fixed bottom (z = -5000), prescribed displacement on top (z = 0)
- Apply uniform downward displacement on top surface
- Solve elastostatics
- Verify: stress at a cell in each layer differs (because moduli differ)
- Quantitative check: for uniaxial strain eps_zz, sigma_zz = (lambda + 2*mu) * eps_zz. Compute expected stress for each layer and verify within 1% tolerance.

Create a config file: `config/examples/layered_elastostatics.config`

### Step 2.9: Write unit test for auxiliary callbacks

Create `tests/unit/test_petscfe_elasticity_aux.cpp`:

- Call f1_elasticity_aux with known input (lambda=30e9, mu=25e9, strain=0.001)
- Verify output stress matches analytical: sigma_zz = (30e9 + 2*25e9)*0.001 = 80e6
- Call g3_elasticity_aux, verify nonzero entries match elasticity tensor

### Verification Gate 2

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc) && ctest --output-on-failure'
```

- All previous 50 tests still pass
- New unit test for aux callbacks passes
- Layered elastostatics integration test passes with stress values matching analytical within 1%
- Homogeneous configs still work (fallback to constants array)

**Do not proceed to Phase 3 until this gate passes.**

### Failure Diagnosis (Phase 2)

- **DMSetAuxiliaryVec not found**: API may differ in 3.22. Check for `PetscObjectCompose` pattern instead. In older PETSc, aux data was attached via `PetscObjectCompose(dm, "A", auxVec)`.
- **Callback gets NULL a[] pointer**: Auxiliary vector not attached before TS solve. Ensure setupAuxiliaryDM() is called before setupTimeStepper().
- **Wrong values in a[]**: Check that auxiliary FE field ordering matches AUX_LAMBDA=0, AUX_MU=1, AUX_RHO=2. Print values inside callback with PetscPrintf for debugging.
- **SNES diverges with heterogeneous material**: Check that Jacobian g3 is consistent with residual f1. A common bug is using constants[] in one and a[] in the other.

---

## Phase 3: Absorbing Boundary Conditions

**Goal**: Implement first-order Clayton-Engquist absorbing boundary conditions so seismic waves exit the domain without reflection.

### Background

On each boundary face, the absorbing condition adds a damping traction:

```
t_i = -rho * c_p * (n_i * n_j) * v_j - rho * c_s * (delta_ij - n_i * n_j) * v_j
```

where:
- `v_j` is particle velocity (time derivative of displacement)
- `n_i` is the outward normal
- `c_p` is P-wave speed: sqrt((lambda + 2*mu) / rho)
- `c_s` is S-wave speed: sqrt(mu / rho)

This is a boundary integral term. In PETSc FEM, it is registered via `PetscDSSetBdResidual` on the boundary labels.

### Step 3.1: Create absorbing BC callback file

Create:

```
src/numerics/AbsorbingBC.cpp
include/numerics/AbsorbingBC.hpp
```

The boundary residual callback has the `PetscBdPointFunc` signature. Check:

```bash
grep -n "PetscBdPointFunc\|PetscDSSetBdResidual\|PetscDSSetBdJacobian" /opt/petsc-3.22.2/include/petscds.h
```

Implement:

```cpp
static void f0_absorbing_bc(PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[],
    const PetscScalar u_t[], const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[],
    const PetscScalar a_t[], const PetscScalar a_x[],
    PetscReal t, const PetscReal x[], const PetscReal n[],
    PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
    // Read material properties from aux fields if heterogeneous,
    // or from constants if homogeneous
    PetscScalar lambda, mu, rho;
    if (NfAux >= 3) {
        lambda = a[aOff[AUX_LAMBDA]];
        mu     = a[aOff[AUX_MU]];
        rho    = a[aOff[AUX_RHO]];
    } else {
        lambda = constants[0];
        mu     = constants[1];
        rho    = constants[2];
    }

    const PetscScalar cp = PetscSqrtReal((lambda + 2.0 * mu) / rho);
    const PetscScalar cs = PetscSqrtReal(mu / rho);

    // Absorbing traction: -rho*cp*(n.v)*n - rho*cs*(v - (n.v)*n)
    if (u_t == NULL) { for (PetscInt d = 0; d < dim; ++d) f0[d] = 0.0; return; }

    PetscScalar ndotv = 0.0;
    for (PetscInt d = 0; d < dim; ++d) ndotv += n[d] * u_t[d];

    for (PetscInt d = 0; d < dim; ++d) {
        f0[d] = rho * cp * ndotv * n[d]
              + rho * cs * (u_t[d] - ndotv * n[d]);
    }
}
```

Also implement the corresponding Jacobian contribution `g0_absorbing_bc` for implicit time stepping. The Jacobian of the absorbing term with respect to u_t is:

```
g0_{ij} = shift * (rho * cp * n_i * n_j + rho * cs * (delta_ij - n_i * n_j))
```

where `shift` is `u_tShift` from PETSc (the coefficient relating du_t/du).

### Step 3.2: Register absorbing BCs on boundary labels

In `setupPhysics()` or a new method `setupAbsorbingBCs()`:

1. Determine which boundary labels correspond to the 5 non-surface faces (or all 6 faces, depending on configuration).
2. For each absorbing boundary label, call:
   ```cpp
   PetscDSSetBdResidual(prob, field_displacement, bd_index, f0_absorbing_bc, NULL);
   PetscDSSetBdJacobian(prob, field_displacement, field_displacement, bd_index,
                         g0_absorbing_bc, NULL, NULL, NULL);
   ```

Check the exact API for `PetscDSSetBdResidual` and `PetscDSSetBdJacobian`:

```bash
grep -n "PetscDSSetBdResidual\|PetscDSSetBdJacobian\|PetscDSGetBoundary\|PetscDSGetNumBoundary" /opt/petsc-3.22.2/include/petscds.h
```

Note: PETSc may use a different mechanism for boundary residuals. If `PetscDSSetBdResidual` does not exist in 3.22, the alternative is to add the absorbing term as a surface integral via a custom FEM integration routine or via `DMAddBoundary` with a NATURAL (Neumann) type and a custom evaluation function.

### Step 3.3: Add config option

```ini
[ABSORBING_BC]
enabled = true
faces = x_min, x_max, y_min, y_max, z_min
# surface (z_max) stays free
```

### Step 3.4: Write verification test

Create `tests/physics_validation/test_absorbing_bc.cpp`:

**Test setup**:
- 20km x 20km x 20km box, homogeneous granite (vp=5500, vs=3175, rho=2650)
- Point force impulse at center (10km, 10km, 10km)
- Two seismometers: one at (15km, 10km, 10km) and one at (5km, 10km, 10km)
- Run for 8 seconds (enough for P-wave to reach boundary and reflect back)
- Transit time to boundary: 10000/5500 = 1.82s. Reflection arrives at ~3.64s.

**Without absorbing BCs**: Run and record seismogram. After 3.64s, reflected arrivals should be visible as secondary peaks.

**With absorbing BCs**: Run same setup. After 3.64s, amplitude should remain below 5% of direct P-wave peak amplitude.

**Quantitative pass/fail**: Let A_direct = max absolute amplitude of direct P arrival. Let A_reflected = max absolute amplitude in the window [3.5s, 5.0s]. Test passes if A_reflected / A_direct < 0.05.

If seismometer SAC output is not yet verified (depends on Phase 1 results), use VTK output at specific timesteps to check field values at receiver locations.

### Step 3.5: Update CMakeLists.txt

Add AbsorbingBC.cpp to the build.

### Verification Gate 3

- All previous tests still pass (50 original + Phase 2 additions)
- Absorbing BC test passes: reflection amplitude < 5% of direct arrival
- Configs without absorbing BCs still work identically (no regressions)

**Do not proceed to Phase 4 until this gate passes.**

### Failure Diagnosis (Phase 3)

- **PetscDSSetBdResidual not found**: Check if the function is named differently in 3.22 or if boundary residuals require a different registration path. Search: `grep -rn "BdResidual\|BdPoint\|BdIntegral" /opt/petsc-3.22.2/include/`.
- **No effect on seismograms**: Verify the boundary labels are correct. Print the label values used in registration and confirm they match labelBoundaries() output.
- **SNES diverges**: The absorbing BC Jacobian may be wrong. Try running with `-snes_fd` to use finite-difference Jacobian as a diagnostic.
- **Reflections still large**: Check the normal vector direction. PETSc may provide inward normals. Flip the sign if needed.

---

## Phase 4: Gravity and Lithostatic Pre-stress

**Goal**: Add gravity body force to the elasticity solver and compute a lithostatic pre-stress initial condition.

### Step 4.1: Add gravity to f0 callback

If using auxiliary-field callbacks (Phase 2), this is already included in `f0_elasticity_aux`. The gravity magnitude is passed via `constants[0]`.

If using the original homogeneous callbacks, create a wrapper or modify the config to pass gravity magnitude. Do NOT modify PetscFEElasticity.cpp. Instead, register `f0_elasticity_aux` even for homogeneous cases when gravity is enabled, with uniform auxiliary field values.

### Step 4.2: Static equilibrium solve

Create a config and example for computing lithostatic state:

```ini
[SIMULATION]
name = lithostatic_equilibrium
type = static_equilibrium
gravity = 9.81

[MESH]
type = box
x_extent = 10000.0
y_extent = 10000.0
z_extent = 5000.0
nx = 10
ny = 10
nz = 50

[MATERIAL]
heterogeneous = true

[LAYER_1]
z_top = 0.0
z_bottom = -500.0
lambda = 3.0e9
mu = 1.5e9
rho = 1800.0

[LAYER_2]
z_top = -500.0
z_bottom = -2000.0
lambda = 15.0e9
mu = 12.0e9
rho = 2400.0

[LAYER_3]
z_top = -2000.0
z_bottom = -5000.0
lambda = 30.0e9
mu = 25.0e9
rho = 2650.0

[BOUNDARY_CONDITIONS]
bottom = fixed
sides = roller
top = free
```

Roller BCs (zero normal displacement on side faces) already exist in `setupBoundaryConditions()` for the poroelastic case (roller_xmin, roller_xmax, roller_ymin, roller_ymax constraining individual displacement components). For the elastostatics-with-gravity case, add the same roller pattern using displacement field 0 instead of field 1.

### Step 4.3: Verify lithostatic stress

At depth z below the surface, the vertical stress should be:

```
sigma_zz(z) = integral from 0 to z of rho(z') * g * dz'
```

For the 3-layer model:
- At z = -250m (mid-alluvium): sigma_zz = 1800 * 9.81 * 250 = 4.4145 MPa
- At z = -1250m (mid-weathered): sigma_zz = 1800*9.81*500 + 2400*9.81*750 = 8.829 + 17.658 = 26.487 MPa
- At z = -3500m (mid-granite): sigma_zz = 8.829 + 2400*9.81*1500 + 2650*9.81*1500 = 8.829 + 35.316 + 38.99 = 83.135 MPa

**Quantitative pass/fail**: Extract sigma_zz at sample points via DMInterpolation or VTK output. Verify within 2% of analytical values.

### Step 4.4: Store lithostatic state for dynamic restart

Save the lithostatic solution vector to a PETSc binary file:

```cpp
PetscViewerBinaryOpen(comm, "lithostatic_state.bin", FILE_MODE_WRITE, &viewer);
VecView(solution, viewer);
PetscViewerDestroy(&viewer);
```

Add config option to load initial condition from file:

```ini
[INITIAL_CONDITIONS]
type = from_file
path = lithostatic_state.bin
```

This allows the dynamic explosion simulation to start from lithostatic equilibrium.

### Verification Gate 4

- All previous tests still pass
- Lithostatic stress at sample depths matches analytical values within 2%
- Save/load of solution vector works round-trip (save lithostatic state, load it, verify values unchanged)

**Do not proceed to Phase 5 until this gate passes.**

### Failure Diagnosis (Phase 4)

- **Roller BCs not working**: Roller BCs already exist in `setupBoundaryConditions()` for the poroelastic case (roller_xmin, roller_xmax, roller_ymin, roller_ymax constraining individual displacement components). For the elastostatics-with-gravity case, add the same roller pattern using displacement field 0 instead of field 1.
- **Stress incorrect at layer boundaries**: Check that material property assignment does not have off-by-one in depth lookup. Print centroid z-coordinates and assigned properties for boundary cells.
- **SNES diverges with gravity**: Gravity adds a body force that makes the problem non-trivial. May need better initial guess. Try ramping gravity from 0 to 9.81 over several load steps.

---

## Phase 5: Explosion Source End-to-End with Seismograms

**Goal**: Run the underground nuclear test simulation from source to synthetic seismograms and verify physical P-wave arrivals.

### Step 5.1: Run DPRK 2017 quick config

```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local bash -c \
  './fsrm ../config/examples/dprk_2017_quick.config 2>&1'
```

Check for:
- SNES convergence at each time step
- SAC files produced in the output directory
- No crashes or PETSc errors

### Step 5.2: If DPRK config does not run, debug

Read the config file to understand what it expects:

```bash
cat config/examples/dprk_2017_quick.config
```

Common issues:
- Config references heterogeneous material (needs Phase 2)
- Config references absorbing BCs (needs Phase 3)
- Config references features not yet implemented
- Explosion source parameters are incorrect

If the config requires Phase 2/3/4 features, create a simplified version that uses homogeneous material and no absorbing BCs, just to verify the explosion source and seismometer pipeline.

### Step 5.3: Verify P-wave arrival times

For a homogeneous medium with vp = 5500 m/s, the P-wave arrival time at distance d is:

```
t_P = d / vp
```

Place seismometers at known distances from the source. Read the SAC files (or VTK output at specific times). The first nonzero signal should arrive at t_P within dt_tolerance = 2 * dt (two time steps).

### Step 5.4: Verify mb magnitude

From Mueller-Murphy:
- mb = 4.45 + 0.75 * log10(W_kt)
- For W = 250 kt: mb = 4.45 + 0.75 * log10(250) = 4.45 + 0.75 * 2.398 = 4.45 + 1.799 = 6.249

Compute mb from synthetic seismograms:
1. Read vertical component velocity seismogram at a reference distance (e.g., 100 km equivalent or scaled)
2. Measure peak P-wave amplitude A
3. Apply standard mb formula: mb = log10(A/T) + Q(delta, h)
4. Compare against 6.25

**Quantitative pass/fail**: |mb_synthetic - 6.25| < 0.2

If SAC reading infrastructure does not exist in the test framework, use VTK output and extract peak displacement at a receiver location.

### Step 5.5: Create standalone explosion example

If the DPRK config is too complex, create a minimal config:

```ini
[SIMULATION]
name = minimal_explosion
type = dynamic
dt = 0.001
total_time = 5.0

[MESH]
type = box
x_extent = 20000.0
y_extent = 20000.0
z_extent = 20000.0
nx = 20
ny = 20
nz = 20

[MATERIAL]
heterogeneous = false
lambda = 30.0e9
mu = 25.0e9
rho = 2650.0

[EXPLOSION]
enabled = true
yield_kt = 250.0
depth_m = 800.0
x = 10000.0
y = 10000.0
z = 19200.0

[SEISMOMETER_1]
name = STA1
x = 15000.0
y = 10000.0
z = 20000.0

[SEISMOMETER_2]
name = STA2
x = 20000.0
y = 10000.0
z = 20000.0
```

### Verification Gate 5

- Explosion simulation runs to completion
- SAC or VTK files produced for each seismometer
- P-wave arrival at STA1 (5km) within 2*dt of expected 5000/5500 = 0.909s
- P-wave arrival at STA2 (10km) within 2*dt of expected 10000/5500 = 1.818s
- mb within 0.2 of 6.25 (if full Mueller-Murphy source is working)

**Do not proceed to Phase 6 until this gate passes.**

---

## Phase 6: Realistic Geology Examples

**Goal**: Create three geologically realistic example configs for well-known nuclear test sites.

### Step 6.1: Punggye-ri (DPRK)

Layered model based on published geology:

| Layer | Depth Range | Vp (m/s) | Vs (m/s) | Density (kg/m3) | Lambda (Pa) | Mu (Pa) |
|---|---|---|---|---|---|---|
| Rhyolite/tuff cap | 0 to -200m | 3500 | 2020 | 2200 | 8.53e9 | 8.98e9 |
| Fractured granite | -200 to -1000m | 4500 | 2598 | 2500 | 1.43e10 | 1.69e10 |
| Competent granite | -1000 to -5000m | 5800 | 3349 | 2650 | 2.73e10 | 2.97e10 |

Compute lambda = rho * (vp^2 - 2*vs^2), mu = rho * vs^2.

Create: `config/examples/punggye_ri_layered.config`

Explosion at 800m depth (within fractured granite layer), yield 250 kt.

### Step 6.2: Nevada Test Site (Pahute Mesa)

| Layer | Depth Range | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|---|---|---|---|---|
| Alluvium | 0 to -300m | 1800 | 1039 | 1600 |
| Welded tuff | -300 to -800m | 3200 | 1848 | 2100 |
| Zeolitized tuff | -800 to -1500m | 4000 | 2309 | 2300 |
| Rhyolite | -1500 to -5000m | 5000 | 2887 | 2500 |

Create: `config/examples/nts_pahute_mesa.config`

Explosion at 600m depth (within welded tuff), yield 150 kt.

### Step 6.3: Degelen Mountain (Semipalatinsk)

| Layer | Depth Range | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|---|---|---|---|---|
| Weathered granite | 0 to -100m | 3000 | 1732 | 2200 |
| Fractured granite | -100 to -500m | 4800 | 2771 | 2550 |
| Competent granite/gneiss | -500 to -5000m | 6000 | 3464 | 2700 |

Create: `config/examples/degelen_mountain.config`

Explosion at 300m depth (within fractured granite), yield 50 kt (typical tunnel test).

### Step 6.4: Place seismometers

For each example, place seismometers at:
- 5 km (local)
- 20 km (regional, if domain is large enough)
- Surface directly above explosion (spall check)

### Verification Gate 6

- All three example configs run to completion
- SNES converges at every time step
- Seismometer output files produced
- P-wave arrival times match expected values for each velocity model
- No test regressions

---

## Phase 7: Verification and Validation Tests

**Goal**: Implement rigorous V&V tests with analytical solutions.

### Step 7.1: Lamb's Problem

**Setup**: Elastic half-space, point vertical force on surface.

| Parameter | Value |
|---|---|
| Domain | 20km x 20km x 10km (z negative down) |
| Material | Homogeneous granite: vp=5500, vs=3175, rho=2650 |
| Source | Vertical point force at (10km, 10km, 0), Ricker wavelet, f0=2 Hz |
| Receivers | Surface line from (10km, 10km, 0) to (20km, 10km, 0), 20 receivers |
| Duration | 5 seconds |
| Absorbing BCs | On all faces except top (free surface) |

**Verification criteria**:
- P-wave arrival at each receiver: t_P = distance / 5500 (within 2*dt)
- S-wave arrival: t_S = distance / 3175 (within 2*dt)
- Rayleigh wave arrival: t_R = distance / (0.92 * 3175) (within 3*dt)
- Amplitude decay: P-wave amplitude should decrease as ~1/r for body waves on the surface

Create: `tests/physics_validation/test_lambs_problem.cpp`
Config: `config/examples/lambs_problem.config`

### Step 7.2: Garvin's Problem

**Setup**: Buried explosion in elastic half-space, verify pP phase.

| Parameter | Value |
|---|---|
| Domain | 20km x 20km x 10km |
| Material | Homogeneous: vp=5500, vs=3175, rho=2650 |
| Source | Isotropic explosion at (10km, 10km, -2000m), 1 kt |
| Receiver | Surface at (15km, 10km, 0) |
| Duration | 8 seconds |
| Absorbing BCs | All faces except top |

**Verification criteria**:
- Direct P path = sqrt(5000^2 + 2000^2) = 5385m, t_direct = 5385/5500 = 0.979s
- pP path: source to surface reflection point, then to receiver. For the vertical incidence approximation, pP delay relative to direct P is approximately 2*h*cos(i)/vp where h=2000m and i is the incidence angle. The exact pP time depends on Snell's law ray tracing, but for this test the important checks are:

Verify:
- Direct P arrives at t ~ 0.98s (within 0.05s)
- A second arrival (pP) exists between 1.0s and 1.5s
- pP has opposite polarity to direct P (free surface reflection reverses polarity)

**Quantitative pass/fail**: 
- Direct P arrival time within 0.05s of analytical
- pP arrival exists and has opposite polarity

Create: `tests/physics_validation/test_garvins_problem.cpp`

### Step 7.3: Terzaghi Consolidation

**Setup**: 1D consolidation column.

| Parameter | Value |
|---|---|
| Domain | 1m x 1m x 10m column |
| Material | lambda=1e9, mu=1e9, phi=0.3, k=1e-15 m^2, mu_f=1e-3 Pa.s, K_f=2.2e9 |
| Top BC | Drained (p=0), traction = -1 MPa |
| Bottom BC | Fixed displacement, undrained (dp/dn=0) |
| Sides | Roller + undrained |

**Analytical solution**: Pressure at depth z and time t:

```
p(z,t) = (4*p0/pi) * sum_{n=0}^{inf} [1/(2n+1)] * sin((2n+1)*pi*z/(2H)) * exp(-(2n+1)^2 * pi^2 * cv * t / (4*H^2))
```

where cv = k/(mu_f * (1/M + alpha^2/(lambda+2*mu))) and p0 is the initial undrained pressure response.

**Quantitative pass/fail**: At t = 0.1 * H^2 / cv, compare numerical pressure profile against first 10 terms of the Fourier series. L2 error < 5% of peak pressure.

Create: `tests/physics_validation/test_terzaghi.cpp`

### Step 7.4: Absorbing BC Quantification

This is the test from Phase 3, Step 3.4, formalized.

**Quantitative pass/fail**: Reflection coefficient < 5% for normal-incidence P-wave.

### Step 7.5: DPRK 2017 mb Check

Use the Punggye-ri config from Phase 6.

**Quantitative pass/fail**: |mb_synthetic - 6.25| < 0.2

### Verification Gate 7

- Lamb's problem: P, S, Rayleigh arrivals at correct times
- Garvin's problem: direct P and pP arrivals at correct times, pP polarity reversed
- Terzaghi: pressure profile L2 error < 5%
- Absorbing BC: reflection coefficient < 5%
- DPRK mb: within 0.2 of 6.25
- All previous tests still pass

---

## Summary of Files Created or Modified

### New files (Phase 2)
- src/numerics/PetscFEElasticityAux.cpp
- include/numerics/PetscFEElasticityAux.hpp
- tests/unit/test_petscfe_elasticity_aux.cpp
- tests/integration/test_layered_elastostatics.cpp
- config/examples/layered_elastostatics.config

### New files (Phase 3)
- src/numerics/AbsorbingBC.cpp
- include/numerics/AbsorbingBC.hpp
- tests/physics_validation/test_absorbing_bc.cpp
- config/examples/absorbing_bc_test.config

### New files (Phase 4)
- config/examples/lithostatic_equilibrium.config
- tests/physics_validation/test_lithostatic.cpp

### New files (Phase 5)
- config/examples/minimal_explosion.config (if DPRK config does not work directly)

### New files (Phase 6)
- config/examples/punggye_ri_layered.config
- config/examples/nts_pahute_mesa.config
- config/examples/degelen_mountain.config

### New files (Phase 7)
- tests/physics_validation/test_lambs_problem.cpp
- tests/physics_validation/test_garvins_problem.cpp
- tests/physics_validation/test_terzaghi.cpp
- config/examples/lambs_problem.config

### Modified files
- src/core/Simulator.cpp (aux DM setup, absorbing BC registration, gravity, initial condition load)
- include/core/Simulator.hpp (new members and methods)
- src/core/ConfigReader.cpp (new config options: heterogeneous material, absorbing BC, gravity, layers)
- CMakeLists.txt (new source files)

### Files NOT modified (per rules)
- src/numerics/PetscFEElasticity.cpp
- src/numerics/PetscFEPoroelasticity.cpp
- src/numerics/PetscFEFluidFlow.cpp
- src/numerics/FaultMeshManager.cpp
- src/physics/CohesiveFaultKernel.cpp

---

## Notes for Claude Code

1. After each phase, always run the full test suite. No exceptions.
2. If a PETSc function does not exist in 3.22.2, search the headers. Do not guess.
3. If SNES diverges, try `-snes_fd` first to rule out Jacobian bugs.
4. Print intermediate values inside callbacks during debugging, then remove prints before committing.
5. If a verification test fails marginally (e.g., 6% instead of 5% tolerance), report the actual value and suggest whether the tolerance should be loosened or the implementation fixed.
6. All git commits should have descriptive messages indicating which phase and step was completed.
7. Do not create any Python scripts. All verification is in C++ within the test framework.
8. When in doubt about PETSc patterns, check PETSc examples in /opt/petsc-3.22.2/src/dm/impls/plex/tutorials/ and /opt/petsc-3.22.2/src/snes/tutorials/.