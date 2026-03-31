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
- Working callback files: PetscFEElasticity, PetscFEPoroelasticity, PetscFEFluidFlow
- Unified PetscDS constants layout (25 elements, set once in setupPhysics):
  [0]=lambda [1]=mu [2]=rho_s [3]=phi [4]=kx [5]=ky [6]=kz
  [7]=cw [8]=co [9]=cg [10]=mu_w [11]=mu_o [12]=mu_g
  [13]=Swr [14]=Sor [15]=Sgr [16]=nw [17]=no [18]=ng
  [19]=krw0 [20]=kro0 [21]=krg0 [22]=biot_alpha [23]=1/M [24]=rho_f
- CohesiveFaultKernel uses PetscDSSetBdResidual for cohesive cell callbacks
- CohesiveFaultKernel reads mode from constants[COHESIVE_CONST_MODE] and friction from constants[COHESIVE_CONST_MU_F]
- PETSc vector field gradient layout: u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d

## What works now
- PetscFEElasticity: f0/f1/g3 for elastostatics and elastodynamics. Constants layout correct.
- PetscFEPoroelasticity: f0/f1 for pressure and displacement, all 4 Jacobian blocks. Constants layout correct.
- PetscFEFluidFlow: single-phase and black-oil callbacks. Constants layout correct.
- CohesiveFaultKernel: registerWithDS() reads existing constants, expands array, patches cohesive slots, registers Bd callbacks.
- FaultMeshManager: createPlanarFaultLabel, splitMeshAlongFault (DMPlexConstructCohesiveCells), extractCohesiveTopology.
- FaultCohesiveDyn: friction models (SlipWeakeningFriction, RateStateFrictionAging) tested and correct.
- CoulombStressTransfer: Hooke stress from strain, fault projection, delta_CFS. Tested.
- Validation tests call PetscFE callbacks directly with known inputs.

## Known bugs
- COHESIVE_CONST_MODE (24) collides with rho_fluid at unified_constants[24]. Must bump cohesive slots to [25]/[26] and expand array to 27.
- Simulator::setupFaultNetwork() creates empty FaultNetwork and returns. Not wired to FaultMeshManager or CohesiveFaultKernel.
- No Lagrange multiplier field created in setupFields() for fault problems.
- IMEX transition exists but is not connected to the new PetscDS callback path.

## Fault component dependency chain
1. setupDM() -- create base mesh
2. setupFaultNetwork() -- split mesh, insert cohesive cells (CHANGES DM TOPOLOGY)
3. setupFields() -- create PetscFE for displacement + Lagrange multiplier
4. setupPhysics() -- register elasticity + cohesive callbacks on PetscDS
5. setupTimeStepper() -- TS with IFunction
6. solve()

## Subagent strategy
Parallelizable:
- FaultMeshManager wiring in setupFaultNetwork() (Simulator.cpp)
- Locked fault integration test (new test file)
- Config file for fault test

Sequential (main agent):
- Constants collision fix (header + Simulator)
- Lagrange multiplier field in setupFields()
- registerWithDS call in setupPhysics()
- Build and integration testing

## Rules
- Do not delete or gut tests. Fix code or GTEST_SKIP.
- Do not refactor, rename classes, or reorganize files.
- Do not modify PetscFEElasticity, PetscFEPoroelasticity, or PetscFEFluidFlow.
- Do not modify friction law implementations (tested and correct).
- Do not modify CohesiveFaultKernel::registerWithDS internals (it works).
- Do not change the IMEX transition logic (future session).
- Always build and test in Docker (Dockerfile.ci) to match CI environment.