# FSRM Development Rules

## Build
- Build Docker image: docker build -f Dockerfile.ci -t fsrm-ci:local .
- Compile: docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON && make -j$(nproc)'
- Test: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

## Critical context for this session
DMPlex cohesive cell insertion is broken. Every test that touches FaultMeshManager::splitMeshAlongFault() GTEST_SKIPs because DMPlexConstructCohesiveCells fails or produces zero cohesive cells.

Root cause: DMPlexLabelCohesiveComplete() is never called. PETSc requires the fault label to include the full closure (vertices + edges + faces) before DMPlexConstructCohesiveCells can split the mesh. Currently only faces are labeled.

Reference PETSc examples that do this correctly:
- $PETSC_DIR/src/dm/impls/plex/tutorials/ex12.c
- $PETSC_DIR/src/snes/tutorials/ex62.c

## Rules
- Build and test in Docker (Dockerfile.ci)
- Do not modify PetscFE callback files
- Do not modify friction law implementations
- Do not delete tests