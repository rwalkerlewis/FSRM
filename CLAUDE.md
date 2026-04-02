# FSRM Development Context

## Build (always use Docker)
```bash
docker build -f Dockerfile.ci -t fsrm-ci:local .
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c \
  'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release \
   -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
  ctest --output-on-failure
```

## PETSc Version: 3.22.2
ALWAYS check API signatures before using any PETSc function:
```bash
docker run --rm fsrm-ci:local bash -c 'grep -n "FUNCTION_NAME" /opt/petsc-3.22.2/include/*.h'
```

## Architecture
- C++17, PETSc 3.22.2 with ctetgen, DMPlex unstructured FEM
- PDE assembly via PetscDS pointwise callbacks (f0, f1, g0, g3)
- Callbacks registered in Simulator::setupPhysics()
- FormFunction converts global to local vectors, calls DMPlexInsertBoundaryValues,
  then DMPlexTSComputeIFunctionFEM
- Unified PetscDS constants array (27 elements), set once in setupPhysics()

## CRITICAL: DMCreateDS and Boundary Condition Ordering
DO NOT CHANGE THE ORDERING OF THESE THREE CALLS IN setupFields():
```cpp
ierr = DMCreateDS(dm); CHKERRQ(ierr);                    // 1. Create DS (required first in PETSc 3.22)
ierr = setupBoundaryConditions(); CHKERRQ(ierr);          // 2. Add BCs via DMAddBoundary (needs DS to exist)
ierr = DMSetLocalSection(dm, nullptr); CHKERRQ(ierr);     // 3. Clear cached section
ierr = DMSetUp(dm); CHKERRQ(ierr);                        // 4. Rebuild section WITH BC constraints
ierr = DMGetDS(dm, &prob); CHKERRQ(ierr);                 // 5. Get DS for later use
```
In PETSc 3.22, DMAddBoundary internally calls DMGetDS(), so DMCreateDS must come
first. But DMCreateDS creates the PetscSection without knowing about BCs, so after
adding BCs we must clear and rebuild the section (steps 3-4) so the constrained
DOFs are properly marked. Without steps 3-4, DMPlexInsertBoundaryValues finds no
constrained DOFs and inserts nothing.

This ordering was debugged over multiple sessions. It produces SNES convergence
with nonzero FNORM. NEVER move DMCreateDS after setupBoundaryConditions (crashes
because DS doesn't exist). NEVER remove the DMSetLocalSection/DMSetUp rebuild
(BCs silently have no effect). NEVER reorder these calls for any reason.

## Constants Layout
```
[0]=lambda [1]=mu [2]=rho_s [3]=phi [4]=kx [5]=ky [6]=kz
[7]=cw [8]=co [9]=cg [10]=mu_w [11]=mu_o [12]=mu_g
[13]=Swr [14]=Sor [15]=Sgr [16]=nw [17]=no [18]=ng
[19]=krw0 [20]=kro0 [21]=krg0 [22]=biot_alpha [23]=1/M [24]=rho_f
[25]=cohesive_mode [26]=cohesive_mu_f
```

## Pipeline Call Order (main.cpp)
```
initializeFromConfigFile()  -> parse config, materials, fluids, explosion
setupDM()                   -> create mesh (simplex if faults enabled)
setupFaultNetwork()         -> PyLith workflow: submesh + cohesive cells
labelBoundaries()           -> label 6 faces of bounding box
setupFields()               -> create PetscFE, DMCreateDS, add BCs, rebuild section
setupPhysics()              -> register PetscDS callbacks + constants
setupTimeStepper()          -> create TS with IFunction/IJacobian
setupSolvers()              -> configure SNES/KSP
setInitialConditions()      -> zero solution vector
run()                       -> TSSolve
```

## What Works Right Now
1. Elastostatics: PetscFEElasticity f0/f1/g3. Uniaxial compression converges with
   nonzero FNORM in 1 iteration. Field 0 = displacement (3 components).
2. PetscFEPoroelasticity: f0/f1 for pressure and displacement, all 4 Jacobian blocks.
   Unit tested. NOT yet run end-to-end in the Simulator.
3. PetscFEFluidFlow: single-phase and black-oil callbacks. Unit tested.
4. Cohesive faults: PyLith submesh workflow in FaultMeshManager. 32 cohesive cells
   inserted, 96 fault vertices. CohesiveFaultKernel registers on PetscDS.
5. Friction laws: slip-weakening and rate-state (aging). Tested.
6. CoulombStressTransfer: Hooke stress, fault projection, delta_CFS. Tested.
7. Boundary conditions: labelBoundaries() + setupBoundaryConditions() + section rebuild.
   Hardcoded BCs for elastostatics (fixed bottom, compression top).
8. Seismometer network: DMInterpolation sampling, SAC output, velocity/acceleration.
9. MuellerMurphySource: moment rate, RDP, moment tensor (simplified physics).
10. ExplosionCoupling: SphericalCavitySource stress at a point (analytical, not FEM).

## What Does NOT Work
1. Poroelasticity end-to-end: BCs only handle elastostatics case. No poroelastic BCs.
2. No injection source term in pressure equation.
3. Explosion source not injected into FEM residual (only analytical stress to faults).
4. Mueller-Murphy corner frequency is wrong: uses vs/(3*depth) instead of psi(W) function.
5. mb-yield relation is wrong: 4.5+0.80*log10(W) instead of 4.45+0.75*log10(W).
6. HydraulicFractureModel (PKN/KGD/P3D) exists standalone, not connected to solver.
7. No dynamic fracture propagation via cohesive cells.
8. Config [BOUNDARY_CONDITIONS] section is parsed but ignored (BCs are hardcoded).

## Required Poroelastic BCs (must be added in setupBoundaryConditions)
For poroelasticity: field 0 = pressure (1 component), field 1 = displacement (3 components).
The Terzaghi consolidation problem requires:

| Boundary | Pressure (field 0) | Displacement (field 1) |
|----------|-------------------|----------------------|
| z_min (bottom) | No BC (impermeable = natural zero-flux) | u=0 all 3 components (fixed) |
| z_max (top) | p=0 (drained, Dirichlet) | u_z=-0.001 applied compression (Dirichlet) |
| x_min | No BC (impermeable) | u_x=0 (roller, component 0 only) |
| x_max | No BC (impermeable) | u_x=0 (roller, component 0 only) |
| y_min | No BC (impermeable) | u_y=0 (roller, component 1 only) |
| y_max | No BC (impermeable) | u_y=0 (roller, component 1 only) |

Impermeable = no DMAddBoundary call for pressure on that face (natural zero-flux BC).
Drained = DMAddBoundary with DM_BC_ESSENTIAL, bc_zero callback, pressure field, component 0.
Roller = DMAddBoundary with DM_BC_ESSENTIAL, bc_zero callback, displacement field, single component.

The setupBoundaryConditions() function must detect `config.solid_model == SolidModelType::POROELASTIC
&& config.fluid_model == FluidModelType::SINGLE_COMPONENT` and use field indices 0 (pressure)
and 1 (displacement) instead of the elastostatics field layout where displacement is field 0.

## Rules
- Build and test in Docker
- Do not modify PetscFE callback math (PetscFEElasticity, PetscFEPoroelasticity, PetscFEFluidFlow)
- Do not modify FaultMeshManager::splitMeshAlongFault (PyLith workflow is correct)
- Do not modify CohesiveFaultKernel::registerWithDS
- Do not modify friction law implementations
- All existing tests must still pass after every phase
- Check PETSc 3.22.2 API signatures before using any function
- Use subagents for parallel independent file creation
- Each phase must build and test before proceeding to next phase
- NEVER change the DMCreateDS/setupBoundaryConditions/DMSetLocalSection/DMSetUp ordering in setupFields(). See the CRITICAL section above. This is a hard rule with no exceptions.