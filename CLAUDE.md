# FSRM Development Rules

## Build
- Build Docker image: docker build -f Dockerfile.ci -t fsrm-ci:local .
- Compile: docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=OFF && make -j$(nproc)'
- Test: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

## Critical context
Cohesive cell insertion fails because createPlanarFaultLabel marks boundary faces. DMPlexConstructCohesiveCells cannot split vertices shared between the fault and domain boundary. The fix is to skip faces with support size < 2 (boundary faces) in createPlanarFaultLabel. Also, the 2x2x2 test meshes may be too small for the fault plane to have any interior (non-boundary) vertices.

## Rules
- Build and test in Docker
- Do not modify PetscFE callback files
- Do not modify friction law implementations
- Do not delete tests