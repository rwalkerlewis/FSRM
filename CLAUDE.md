# FSRM Development Rules

## Build
- Build Docker image: docker build -f Dockerfile.ci -t fsrm-ci:local .
- Compile: docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
- Test: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure
- Run example: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ./fsrm -c ../config/test_elastostatics.config

## Architecture
- PDE assembly uses PETSc DMPlex FEM (PetscDS pointwise callbacks)
- Simulator::FormFunction delegates to DMPlexTSComputeIFunctionFEM when use_fem_time_residual_ = true
- Callbacks registered via PetscDSSetResidual / PetscDSSetJacobian in Simulator::setupPhysics()
- If use_fem_time_residual_ is false, solver does nothing (F = U_t placeholder)
- Working callback files: PetscFEFluidFlow (fluid), PetscFEElasticity (elastostatics/dynamics), PoroelasticSolver (standalone Biot)
- PETSc vector field gradient layout: u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d
- PetscDS constants are GLOBAL per DS. A single constants[] array is shared by ALL callbacks on the same DS. Never call PetscDSSetConstants twice with different sizes.

## What works now
- PetscFEElasticity: f0/f1/g3 for elastostatics and elastodynamics (first-order form). I2-form callbacks exist but are not wired.
- Simulator registers elasticity callbacks and sets use_fem_time_residual_ = true
- Fluid flow PetscDS callbacks (PetscFEFluidFlow) are registered and work
- PoroelasticSolver: standalone solver with correct Biot callbacks (not routed through Simulator)
- CI workflow has || true removed from ctest steps
- Figures directory removed from git

## Known bugs in current code
- PetscDS constants collision: fluid sets 19, elasticity overwrites with 3. Coupled runs will crash.
- Dynamics uses TSBEULER (first-order, massively dissipative). I2 callbacks and TSALPHA2 are not wired.
- No FormFunction2 / FormJacobian2 for second-order TS
- No validation test proves the elasticity solve produces correct nonzero displacement
- No PetscFEPoroelasticity -- Biot coupling through Simulator doesn't exist

## Subagent strategy
Use subagents for independent file creation. Parallelizable:
- PetscFEPoroelasticity.hpp/.cpp (new, no dependencies on other new files)
- Dockerfile.cuda (new, independent)
- New test files in tests/physics_validation/
- Config files in config/

Sequential (main agent only):
- Simulator.cpp wiring (depends on callback files)
- Constants array unification (touches fluid + elasticity + poroelasticity code paths)
- Integration testing

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
- Do not modify PetscFEFluidFlow (working reference).
- Do not modify friction law implementations (tested and correct).