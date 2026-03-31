# FSRM Development Rules

## Build
- Build with: docker build -f Dockerfile.ci -t fsrm-ci:local .
- Compile with: docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c 'mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=OFF && make -j$(nproc)'
- Test with: docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local ctest --output-on-failure

## Rules
- Do not delete or gut tests. Fix code or GTEST_SKIP.
- Do not refactor, rename classes, or reorganize files.
- Do not modify CMakeLists.txt structure.
- Do not change the CI job structure.
- Stub unimplemented methods with // STUB comment, return 0/false/{} as appropriate.
- Physics test values are constraints. Fix implementations, not expected values.
- PETSc patterns: PetscErrorCode returns, PetscFunctionBeginUser/PetscFunctionReturn(0).
- Everything is in namespace FSRM. Test fixtures are global (friend class ::TestName syntax).
