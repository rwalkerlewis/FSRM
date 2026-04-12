# FSRM: Wire the PETSc FEM Solver Pipeline

## The Problem

Simulator::setupPhysics() registers PetscDS weak-form callbacks (f0, f1, g0, g3) for 
fluid flow via PetscFEFluidFlow, sets `use_fem_time_residual_ = true`, and 
FormFunction correctly delegates to DMPlexTSComputeIFunctionFEM. This path works.

For elastostatics, elastodynamics, and poroelastodynamics, setupPhysics() creates 
PhysicsKernel objects but NEVER registers PetscDS callbacks. So `use_fem_time_residual_` 
stays false, FormFunction falls back to `F = U_t`, and the solver produces zero 
displacement. The code runs without error but computes nothing.

## What Already Exists

PoroelasticSolver (src/numerics/PoroelasticSolver.cpp, include/numerics/PoroelasticSolver.hpp) 
has COMPLETE, CORRECT PetscDS callbacks:

- f0_displacement, f1_displacement_{x,y,z} -- Hooke's law stress with Biot coupling
- f0_pressure, f1_pressure -- Darcy flow with storage
- f0_saturation, f1_saturation -- transport equation  
- g3_{uxux, uyuy, uzuz, uxuy, ...} -- full elasticity stiffness tensor Jacobian
- g2_{ux_p, uy_p, uz_p}, g1_{p_ux, p_uy, p_uz} -- Biot coupling Jacobian
- Proper PetscDS constant layout: [0..16] = fluid/rock props, [17]=biot, [18]=lambda, [19]=mu

This uses per-component scalar fields (field 0=P, 1=Sw, 2=ux, 3=uy, 4=uz) which is 
one valid PETSc DMPlex FEM approach.

PetscFEFluidFlow (include/numerics/PetscFEFluidFlow.hpp, src/numerics/PetscFEFluidFlow.cpp) 
provides the callback pattern already used by the Simulator for fluid flow.

## What Needs to Happen

### Step 1: Create PetscFEElasticity callbacks

Create `include/numerics/PetscFEElasticity.hpp` and `src/numerics/PetscFEElasticity.cpp` 
following the exact same pattern as PetscFEFluidFlow. These are free-standing static 
functions in a namespace, NOT class methods.

For elastostatics (single vector displacement field, 3 components):

```
namespace FSRM {
namespace PetscFEElasticity {

// Residual: -div(sigma) = f_body
// f0: body force (gravity, source term). For now, zero.
// f1: sigma_{ij} = lambda * eps_kk * delta_ij + 2*mu * eps_ij
//     where eps_ij = 0.5*(du_i/dx_j + du_j/dx_i)
//     u is a 3-component vector field, so u_x has layout [comp*dim + d]

void f0_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);
void f1_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);

// Jacobian: g3[i*dim*dim*dim + j*dim*dim + k*dim + l] = C_{ijkl}
// C_{ijkl} = lambda*delta_ij*delta_kl + mu*(delta_ik*delta_jl + delta_il*delta_jk)
void g3_elastostatics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);

// For elastodynamics, add inertial f0 term:
// f0_i = rho * u_tt_i (when using TS implicit)
void f0_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);
void f1_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);
void g0_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);
void g3_elastodynamics(PetscInt dim, PetscInt Nf, PetscInt NfAux, ...);

} // namespace PetscFEElasticity
} // namespace FSRM
```

Constants layout for these callbacks:
- constants[0] = lambda (first Lame parameter)
- constants[1] = mu (shear modulus, second Lame parameter)
- constants[2] = rho (density, used by elastodynamics f0)

IMPORTANT: With a 3-component vector field (not per-component scalars), the gradient 
layout is different from PoroelasticSolver. For a vector field with nc components:
- u[uOff[f] + c] = component c of field f
- u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d for field f
- f1[c*dim + d] = flux component for test function component c, direction d

Look at PETSc examples ex56.c (elasticity) and ex77.c (Stokes) for the correct 
vector field gradient layout.

### Step 2: Wire into Simulator::setupPhysics()

In the elastostatics/elastodynamics section of setupPhysics() (around line 1280, 
inside the `if (config.enable_geomechanics)` and `if (config.enable_elastodynamics)` 
blocks), add PetscDS registration AFTER the kernel is added:

```cpp
// After adding the elastodynamics kernel...
{
    const MaterialProperties& mat = material_props.empty() ? MaterialProperties() : material_props[0];
    double E = mat.youngs_modulus;
    double nu = mat.poisson_ratio;
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    PetscScalar constants[3];
    constants[0] = lambda;
    constants[1] = mu;
    constants[2] = mat.density;
    
    PetscErrorCode ierr;
    ierr = PetscDSSetConstants(prob, 3, constants); CHKERRQ(ierr);
    
    // For elastostatics:
    ierr = PetscDSSetResidual(prob, field_index, 
                              PetscFEElasticity::f0_elastostatics, 
                              PetscFEElasticity::f1_elastostatics); CHKERRQ(ierr);
    ierr = PetscDSSetJacobian(prob, field_index, field_index,
                              nullptr, nullptr, nullptr,
                              PetscFEElasticity::g3_elastostatics); CHKERRQ(ierr);
    
    use_fem_time_residual_ = true;
}
```

The field_index depends on how setupFields() adds the displacement field. Check 
what index the displacement PetscFE gets added at (it should be 0 for standalone 
elastostatics, or after pressure for coupled problems).

### Step 3: Create test config and verify

Create `config/test_elastostatics.config` (or use the existing one if it has the 
right settings):

```ini
[SIMULATION]
solid_model = ELASTIC
enable_geomechanics = true
fluid_model = NONE
max_timesteps = 1
dt_initial = 1.0
end_time = 1.0

[GRID]
nx = 4
ny = 4
nz = 4
Lx = 1.0
Ly = 1.0
Lz = 1.0
mesh_type = CARTESIAN

[MATERIAL]
youngs_modulus = 1e10
poisson_ratio = 0.25
density = 2700.0
```

Then build and run:
```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON \
    -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc) && 
  ./fsrm -c ../config/test_elastostatics.config
'
```

It should:
- Not segfault
- Create PetscFE displacement field
- Register PetscDS callbacks  
- Run SNES solve
- Report convergence (even if to zero solution with zero BCs -- that's correct)

### Step 4: Add a real boundary-value problem test

The zero-solution test proves the plumbing works. Now add a test that verifies 
actual physics. Create a test in tests/physics_validation/ that sets up:

- Unit cube, E=1, nu=0.25
- Fixed displacement on z=0 face (Dirichlet)
- Applied traction on z=1 face (Neumann) 
- Verify that the computed displacement at z=1 matches u_z = F*L/(E*A) for uniaxial tension

This is the simplest possible validation that the elasticity stiffness tensor is 
correct.

### Step 5: Remove GTEST_SKIP from MMS tests

After the PetscDS callbacks are correct, the MMS tests in 
test_mms_elasticity.cpp and test_mms_wave_propagation.cpp should pass with 
real weak-form residuals. Remove the GTEST_SKIP guards that were added because 
the placeholder strong-form residuals were wrong.

If they still don't pass (because the test uses PhysicsKernel::residual() rather 
than PetscDS callbacks), then update those tests to exercise the PetscDS callbacks 
directly.

## What NOT to do

- Do not refactor PoroelasticSolver. It works. Leave it alone.
- Do not try to implement dynamics (Newmark, TSALPHA2) in this session. 
  Elastostatics first. Dynamics is session 2.
- Do not touch fluid flow, explosion, fault, or any other physics module.
- Do not reorganize the Simulator class or change the callback delegation pattern.
  Just add the missing PetscDSSetResidual/PetscDSSetJacobian calls.

## Reference

PETSc examples with the same pattern:
- src/snes/tutorials/ex56.c (linear elasticity)
- src/snes/tutorials/ex77.c (Stokes, vector + scalar field)
- src/ts/tutorials/ex76.c (elastodynamics)
- PETSc manual section 5.7 (Pointwise Functions)

PoroelasticSolver callbacks in this repo:
- src/numerics/PoroelasticSolver.cpp lines 680-900 (working Hooke's law implementation)