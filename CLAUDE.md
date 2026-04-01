# FSRM Development Rules

## Build
- docker build -f Dockerfile.ci -t fsrm-ci:local .
- docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
- docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

## Architecture
- PETSc 3.22.2 with ctetgen, DMPlex FEM
- PetscDS callbacks registered in setupPhysics()
- BCs registered via DMAddBoundary in setupBoundaryConditions(), BEFORE DMCreateDS
- DMPlexInsertBoundaryValues injects BC values into solution in setInitialConditions()
- FormFunction calls DMPlexTSComputeIFunctionFEM on local vectors with BC insertion
- Unified constants array (27 elements)

## Working physics paths
- Elastostatics: PetscFEElasticity f0/f1/g3, field 0 = displacement (3 components)
- Poroelasticity: PetscFEPoroelasticity, field 0 = pressure (1 component), field 1 = displacement (3 components)
- Single-phase flow: PetscFEFluidFlow, field 0 = pressure (1 component)
- Cohesive faults: PyLith workflow, CohesiveFaultKernel on hybrid cells

## Current BCs (hardcoded in setupBoundaryConditions)
- Elastostatics: fixed bottom (u=0 at z_min), applied compression (u_z=-0.001 at z_max)
- Poroelasticity: drained top (p=0 at z_max) added on top of elastostatics BCs
- Single-phase flow: fixed pressure (p=0 at x_min)

## Known gaps
- BCs are hardcoded, not read from config [BOUNDARY_CONDITIONS] section
- No Neumann (traction/flux) BCs, only Dirichlet
- setInitialConditions always zeros the solution
- No injection source term in the pressure equation
- HydraulicFractureModel exists standalone but is not connected to the solver

## Rules
- Build and test in Docker
- Do not modify PetscFE callback files
- Do not modify FaultMeshManager or CohesiveFaultKernel
- All existing tests must pass
- Check PETSc 3.22.2 API signatures before using any function