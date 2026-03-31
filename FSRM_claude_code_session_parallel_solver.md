# FSRM: Solver Pipeline Completion

Read CLAUDE.md first. Build and test everything in Docker.

## Current State

PetscFEElasticity.hpp/.cpp exists with correct Hooke's law callbacks (f0, f1, g3) for
elastostatics and elastodynamics, plus I2-form callbacks for second-order dynamics.
Simulator::setupPhysics() registers these and sets use_fem_time_residual_ = true.
Config files test_elastostatics.config and test_elastodynamics.config exist.

## What's Broken

### Bug 1: PetscDS constants collision
In Simulator::setupPhysics(), fluid flow calls PetscDSSetConstants(prob, 19, constants)
at line ~1377, then elasticity calls PetscDSSetConstants(prob, 3, elast_constants) at
line ~1421. The second call OVERWRITES the first. PetscDS has ONE constants array.

Fix: create a unified constants layout. Define it in a header so all callback files agree:

```
// Constants layout (shared across all PetscDS callbacks):
// [0]  = lambda (first Lame parameter)
// [1]  = mu (shear modulus)
// [2]  = rho_solid
// [3]  = porosity
// [4]  = permeability (isotropic, m^2)
// [5]  = fluid_viscosity
// [6]  = fluid_compressibility
// [7]  = biot_coefficient
// [8]  = biot_modulus (1/M)
// [9]  = rho_fluid
// [10-18] = reserved for black oil / compositional (same layout as current PetscFEFluidFlow)
// Total: 19 constants (or more if needed)
```

Update PetscFEElasticity callbacks to read lambda from constants[0], mu from constants[1],
rho from constants[2]. They already do this. Then update PetscFEFluidFlow callbacks to read
from the same array (fluid properties start at index 3+). Call PetscDSSetConstants ONCE
with the full array.

### Bug 2: Dynamics stuck on TSBEULER
I2-form callbacks (f0_elastodynamics_I2, etc.) exist in PetscFEElasticity.cpp but are
dead code. Simulator uses TSBEULER with first-order TSSetIFunction for everything.

For wave propagation, TSBEULER is first-order accurate and massively dissipative.
It will smear out any wave within a few timesteps. Need TSALPHA2 or at minimum TSCN
(Crank-Nicolson, second-order).

Fix: Add FormFunction2 and FormJacobian2 to Simulator for the I2 interface:

```cpp
// In Simulator.hpp, add to private:
static PetscErrorCode FormFunction2(TS ts, PetscReal t, Vec U, Vec U_t,
                                     Vec U_tt, Vec F, void *ctx);
static PetscErrorCode FormJacobian2(TS ts, PetscReal t, Vec U, Vec U_t,
                                     Vec U_tt, PetscReal v, PetscReal a,
                                     Mat J, Mat P, void *ctx);
```

Implement them delegating to DMPlexTSComputeI2FunctionFEM / DMPlexTSComputeI2JacobianFEM.
IMPORTANT: Check PETSc 3.20 headers for exact signatures. They may differ from 3.21+.

In setupTimeStepper(), when elastodynamics is enabled:
```cpp
ierr = TSSetType(ts, TSALPHA2); CHKERRQ(ierr);
ierr = TSSetI2Function(ts, nullptr, FormFunction2, this); CHKERRQ(ierr);
ierr = TSSetI2Jacobian(ts, jacobian, jacobian, FormJacobian2, this); CHKERRQ(ierr);
```

If TSALPHA2 / I2Function is not available in PETSc 3.20 or the DMPlex I2 compute
functions don't exist, fall back to Newmark-beta implemented manually in the first-order
form, or use TSCN as a second-order alternative that uses the same IFunction interface.
Check what's actually available before writing code.

### Bug 3: No validation
Nothing proves the solver produces correct displacement. Need a test that:
1. Sets up a simple BVP (uniaxial tension or patch test)
2. Solves with SNES
3. Checks displacement against analytical solution

## Parallel Subagent Tasks

### Subagent A: PetscFEPoroelasticity callbacks

Create `include/numerics/PetscFEPoroelasticity.hpp` and `src/numerics/PetscFEPoroelasticity.cpp`.

Two-field system: field 0 = pressure (1 component), field 1 = displacement (3 components).

Uses the SAME unified constants layout as PetscFEElasticity (lambda at [0], mu at [1],
rho_solid at [2], porosity at [3], permeability at [4], fluid_viscosity at [5],
fluid_compressibility at [6], biot_coefficient at [7], biot_modulus at [8], rho_fluid at [9]).

```cpp
namespace FSRM {
namespace PetscFEPoroelasticity {

// Pressure equation: (1/M)*dp/dt + alpha*div(du/dt) + div((k/mu)*grad(p)) = q
void f0_pressure(...);   // (1/M)*p_t + alpha*div(u_t)
void f1_pressure(...);   // (k/mu)*grad(p)

// Displacement equation: -div(sigma - alpha*p*I) = f
void f0_displacement(...);  // 0 (quasi-static) or rho*u_tt (dynamic)
void f1_displacement(...);  // sigma_{cd} - alpha*p*delta_{cd}

// Jacobians (4 blocks):
void g0_pp(...);   // (1/M)*shift
void g3_pp(...);   // (k/mu)*delta_{dd}
void g1_pu(...);   // alpha*shift*delta_{cd}  (pressure eq coupling to displacement)
void g2_up(...);   // -alpha*delta_{cd}  (displacement eq coupling to pressure)
void g3_uu(...);   // C_{cdef} elasticity stiffness tensor

} // namespace
} // namespace
```

Reference: PoroelasticSolver.cpp lines 600-900. The math is identical, the index arithmetic
differs because PoroelasticSolver uses per-component scalar fields (P, Sw, ux, uy, uz = fields
0-4) while PetscFEPoroelasticity uses (P=field 0, u=field 1 with 3 components).

For a 3-component vector field 1:
- u_x[uOff_x[1] + c*dim + d] = d(u_c)/dx_d
- f1[c*dim + d] = stress contribution for component c, direction d

For scalar field 0:
- u[uOff[0]] = pressure
- u_x[uOff_x[0] + d] = dp/dx_d
- u_t[uOff[0]] = dp/dt


### Subagent B: Dockerfile.cuda

Create `Dockerfile.cuda` based on Dockerfile.ci but with CUDA support.

```dockerfile
FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

# Same build deps as Dockerfile.ci plus CUDA toolkit
# Build PETSc with: --with-cuda --with-cuda-dir=/usr/local/cuda
# Add: --download-kokkos --download-kokkos-kernels for GPU backend
# Set COPTFLAGS/CXXOPTFLAGS to -O3

# PETSc configure line:
# ./configure \
#     --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
#     --download-fblaslapack \
#     --with-cuda --with-cuda-dir=/usr/local/cuda \
#     --with-debugging=0 \
#     COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'
```

Also create `.github/workflows/cuda-ci.yml` that builds with this Dockerfile
but only runs on push to main (not every branch -- GPU runners are expensive).
Use `runs-on: ubuntu-latest` with Docker build only (no GPU test execution),
just verify compilation with -DENABLE_CUDA=ON.

In CMakeLists.txt, when ENABLE_CUDA=ON, add VecType/MatType selection:
```cmake
if(ENABLE_CUDA)
    add_definitions(-DUSE_CUDA)
    # PETSc handles GPU Vec/Mat via runtime options:
    # -vec_type cuda -mat_type aijcusparse
endif()
```

The PetscDS pointwise callbacks (PetscFEElasticity, PetscFEPoroelasticity, etc.)
run on GPU automatically when PETSc is built with CUDA and Vec/Mat types are set
to CUDA variants. No custom CUDA kernels needed.


### Subagent C: Validation tests

Create `tests/physics_validation/test_elastostatics_patch.cpp`:
```cpp
// Patch test: uniform strain should produce exact stress
// - Unit cube, Q1 hex elements
// - Prescribed displacement: u_x = eps_0 * x, u_y = u_z = 0
// - Check: sigma_xx = (lambda + 2*mu)*eps_0, sigma_yy = sigma_zz = lambda*eps_0
// - This tests that the stiffness tensor g3 is correct
```

Create `tests/physics_validation/test_terzaghi.cpp`:
```cpp
// Terzaghi consolidation: 1D poroelastic column
// - Column height H, load P0 on top, drained at top, impermeable bottom
// - Analytical: p(z,t) = sum_n (4*P0/(pi*(2n+1))) * sin((2n+1)*pi*z/(2*H)) * exp(-(2n+1)^2*pi^2*cv*t/(4*H^2))
// - cv = k/(mu_f * (1/M + alpha^2/(lambda+2*mu)))
// - Check pressure at mid-height at t = 0.1*H^2/cv against series solution
// This validates the Biot coupling Jacobian blocks
```

Register both in tests/CMakeLists.txt under PHYSICS_VALIDATION_TEST_SOURCES and add
CTest registrations.

Create `config/test_poroelasticity.config`:
```ini
[SIMULATION]
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
max_timesteps = 50
dt_initial = 0.001
end_time = 0.05

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
permeability_x = 1e-12
permeability_y = 1e-12
permeability_z = 1e-12
biot_coefficient = 1.0

[FLUID]
density = 1000.0
viscosity = 0.001
compressibility = 4.5e-10
```


## Main Agent Tasks (after subagents complete)

### 1. Unify PetscDS constants

In Simulator::setupPhysics(), build ONE constants array:

```cpp
PetscScalar unified_constants[20] = {0};
// Always set elastic properties
unified_constants[0] = lambda;
unified_constants[1] = mu;
unified_constants[2] = mat.density;
// Set fluid properties if fluid flow is enabled
if (config.fluid_model != FluidModelType::NONE) {
    unified_constants[3] = mat.porosity;
    unified_constants[4] = mat.permeability_x * 9.869233e-16; // mD -> m^2
    unified_constants[5] = flu.viscosity;
    unified_constants[6] = flu.water_compressibility;
    unified_constants[7] = mat.biot_coefficient;
    unified_constants[8] = 1.0 / (mat.porosity * flu.water_compressibility); // 1/M approx
    unified_constants[9] = flu.density;
    // [10-18] = black oil / compositional specific (existing layout)
}
ierr = PetscDSSetConstants(prob, 20, unified_constants); CHKERRQ(ierr);
```

Call it ONCE. Remove the separate PetscDSSetConstants calls from both the fluid and
elasticity blocks.

Update PetscFEFluidFlow callbacks to read from the unified layout. The pressure
callbacks currently read porosity from constants[0], permeability from constants[1], etc.
They need to read from the new offsets (constants[3], constants[4], etc.).

### 2. Wire PetscFEPoroelasticity into Simulator

In setupPhysics(), when solid_model == POROELASTIC and fluid_model != NONE:

```cpp
// Field 0: pressure (from fluid flow setup)
// Field 1: displacement (from geomechanics setup)
ierr = PetscDSSetResidual(prob, 0, PetscFEPoroelasticity::f0_pressure,
                           PetscFEPoroelasticity::f1_pressure); CHKERRQ(ierr);
ierr = PetscDSSetResidual(prob, 1, PetscFEPoroelasticity::f0_displacement,
                           PetscFEPoroelasticity::f1_displacement); CHKERRQ(ierr);

// All 4 Jacobian blocks
ierr = PetscDSSetJacobian(prob, 0, 0, PetscFEPoroelasticity::g0_pp,
                           nullptr, nullptr, PetscFEPoroelasticity::g3_pp); CHKERRQ(ierr);
ierr = PetscDSSetJacobian(prob, 0, 1, nullptr, PetscFEPoroelasticity::g1_pu,
                           nullptr, nullptr); CHKERRQ(ierr);
ierr = PetscDSSetJacobian(prob, 1, 0, nullptr, nullptr,
                           PetscFEPoroelasticity::g2_up, nullptr); CHKERRQ(ierr);
ierr = PetscDSSetJacobian(prob, 1, 1, nullptr, nullptr, nullptr,
                           PetscFEPoroelasticity::g3_uu); CHKERRQ(ierr);
use_fem_time_residual_ = true;
```

Make sure setupFields() creates pressure (1 comp) as field 0 and displacement (3 comp)
as field 1 when poroelastic mode is selected.

### 3. Wire TSALPHA2 or fix dynamics

Check PETSc 3.20 for TSSetI2Function availability:
```bash
docker run --rm fsrm-ci:local bash -c 'grep -r "TSSetI2Function" /opt/petsc-3.20.0/include/'
```

If it exists: add FormFunction2/FormJacobian2, wire TSALPHA2.
If not: use TSCN (Crank-Nicolson) which is second-order and uses the same IFunction
interface. Change TSSetType from TSBEULER to TSCN for dynamics.

### 4. Build, test, verify

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON &&
  make -j$(nproc) &&
  ctest --output-on-failure &&
  echo "=== Running elastostatics ===" &&
  ./fsrm -c ../config/test_elastostatics.config &&
  echo "=== Running elastodynamics ===" &&
  ./fsrm -c ../config/test_elastodynamics.config
'
```

Elastostatics should converge in 1 SNES iteration (linear problem).
Elastodynamics should run multiple timesteps without diverging.
All existing tests must still pass.

## Do NOT
- Modify PoroelasticSolver.hpp/.cpp
- Modify PetscFEFluidFlow if you can avoid it (but constants layout change may require it)
- Modify friction law implementations
- Delete or gut any tests
- Change CMakeLists.txt structure (GLOB_RECURSE handles new .cpp files automatically)