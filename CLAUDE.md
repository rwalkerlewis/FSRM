# FSRM Development Rules

## Build
- Build Docker image: docker build -f Dockerfile.ci -t fsrm-ci:local .
- Compile: docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
- Test: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

## Architecture
- PDE assembly: PETSc DMPlex FEM with PetscDS pointwise callbacks
- FormFunction delegates to DMPlexTSComputeIFunctionFEM when use_fem_time_residual_ = true
- Unified PetscDS constants (27 elements), set once in setupPhysics
- PETSc 3.22.2 with ctetgen, simplex meshes when faults enabled
- Cohesive cells working via PyLith submesh workflow

## Critical gap
setupBoundaryConditions() is empty. No Dirichlet or Neumann BCs applied.
setInitialConditions() does VecSet(solution, 1e7) regardless of physics.
Without BCs, elastostatic systems are singular. No end-to-end solve produces physical results.

## What works
- PetscFEElasticity callbacks (unit tested)
- PetscFEPoroelasticity callbacks (unit tested)
- PetscFEFluidFlow callbacks (unit tested)
- Mesh creation (structured hex, simplex for faults)
- Fault mesh splitting (PyLith workflow, 50/50 tests pass)
- Config parsing for simulation, grid, rock, fluid sections
- Seismometer network with DMInterpolation sampling

## Pipeline call order
1. initializeFromConfigFile() - parse config
2. setupDM() - create mesh
3. setupFaultNetwork() - split mesh if faults enabled
4. labelBoundaries() - NEW: label boundary faces by coordinate
5. setupFields() - create PetscFE
6. setupPhysics() - register PetscDS callbacks + constants
7. setupBoundaryConditions() - NEW: register Dirichlet BCs via PetscDSAddBoundary
8. setupTimeStepper() - create TS
9. setupSolvers() - configure SNES/KSP
10. setInitialConditions() - set IC
11. run() - TSSolve

## Rules
- Build and test in Docker
- Do not modify PetscFE callback files
- Do not modify FaultMeshManager or CohesiveFaultKernel
- Check PETSc 3.22.2 API signatures before using any function
- All 50 existing tests must still pass