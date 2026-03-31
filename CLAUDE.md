# FSRM Development Rules

## Build
- Build Docker image: docker build -f Dockerfile.ci -t fsrm-ci:local .
- Compile: docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
- Test: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure
- Run: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ./fsrm -c ../config/test_elastostatics.config

## Architecture
- PDE assembly uses PETSc DMPlex FEM (PetscDS pointwise callbacks)
- Simulator::FormFunction delegates to DMPlexTSComputeIFunctionFEM when use_fem_time_residual_ = true
- Callbacks must be registered via PetscDSSetResidual / PetscDSSetJacobian in Simulator::setupPhysics()
- If use_fem_time_residual_ is false, the solver does nothing (F = U_t placeholder)
- Working callback examples: PetscFEFluidFlow (fluid), PoroelasticSolver (Hooke + Biot)
- PETSc vector field gradient layout: u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d

## Rules
- Do not delete or gut tests. Fix code or GTEST_SKIP.
- Do not refactor, rename classes, or reorganize files.
- Do not modify CMakeLists.txt structure.
- Do not change the CI job structure.
- Physics test values are constraints. Fix implementations, not expected values.
- PETSc patterns: PetscErrorCode returns, PetscFunctionBeginUser/PetscFunctionReturn(0).
- Everything is in namespace FSRM. Test fixtures are global (friend class ::TestName syntax).
- Always build and test in Docker (Dockerfile.ci) to match CI environment.
- PoroelasticSolver works. Do not refactor it.

## Current priority
Make the Simulator actually solve PDEs. The solver pipeline (DM -> PetscFE -> PetscDS -> TS -> SNES) 
has plumbing but PetscDS callbacks are only registered for fluid flow. Elastostatics and 
elastodynamics need PetscDS callbacks registered so DMPlexTSComputeIFunctionFEM computes 
real residuals instead of the F=U_t placeholder.