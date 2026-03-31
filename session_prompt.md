# FSRM: Debug PETSc Cohesive Cells in Docker CI

Read CLAUDE.md first. Everything runs in Docker.

## Problem

Every test that touches DMPlex cohesive cell operations GTEST_SKIPs on failure.
We have no idea which PETSc functions work and which don't in the CI image.
The fault pipeline (FaultMeshManager -> CohesiveFaultKernel -> Simulator) is
wired but possibly entirely non-functional.

## Step 1: Write a standalone diagnostic program

Create `tests/debug_cohesive.cpp` -- NOT a GTest, just a standalone PETSc main
that tests each step and reports PASS/FAIL. This program should:

1. Create a 3D hex box mesh (DMPlexCreateBoxMesh, 4x4x4)
2. Print mesh stats (cells, faces, vertices)
3. Create a "fault" label (DMCreateLabel)
4. Mark faces near z=0.5 (DMPlexComputeCellGeometryFVM for centroids, DMLabelSetValue)
5. Print how many faces were marked
6. Call DMPlexConstructCohesiveCells
7. Print mesh stats after split (count bulk vs cohesive cells by DMPlexGetCellType)
8. Try PetscFECreateDefault for displacement (3 comp) and Lagrange multiplier (3 comp)
9. Try DMSetField for both fields on the split mesh
10. Try DMCreateDS
11. Try PetscDSSetResidual on field 0 (bulk)
12. Try PetscDSSetBdResidual on field 1 (cohesive)
13. Try DMCreateGlobalVector -- print size
14. Try DMCreateMatrix -- print size

Each step should print PASS or FAIL with the PETSc error code. If a step fails,
print the error and continue to the next step where possible (don't abort).

Use a macro like:
```cpp
#define TRY_STEP(name, call) do { \
    PetscErrorCode _e = (call); \
    if (_e) { PetscPrintf(PETSC_COMM_WORLD, "FAIL: %s (error %d)\n", name, (int)_e); } \
    else    { PetscPrintf(PETSC_COMM_WORLD, "PASS: %s\n", name); } \
} while(0)
```

If the lambda for a dummy callback doesn't compile (C function pointer issue), use
a static function instead:
```cpp
static void dummy_f0(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[], const PetscInt aOff[],
                     const PetscInt aOff_x[], const PetscScalar a[],
                     const PetscScalar a_x[], const PetscScalar a_t[],
                     PetscReal t, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[],
                     PetscScalar f0[]) {
    (void)dim;(void)Nf;(void)NfAux;(void)uOff;(void)uOff_x;
    (void)u;(void)u_t;(void)u_x;(void)aOff;(void)aOff_x;
    (void)a;(void)a_x;(void)a_t;(void)t;(void)x;
    (void)numConstants;(void)constants;
    f0[0] = f0[1] = f0[2] = 0.0;
}
```

Add a build target. Either add to CMakeLists.txt:
```cmake
add_executable(debug_cohesive tests/debug_cohesive.cpp)
target_link_libraries(debug_cohesive ${PETSC_LIBRARIES} ${MPI_CXX_LIBRARIES})
target_include_directories(debug_cohesive PRIVATE ${PETSC_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
```

Or compile manually in Docker:
```bash
cd build && mpicxx -o debug_cohesive ../tests/debug_cohesive.cpp \
  -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include \
  -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lm
```

## Step 2: Build and run

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON &&
  make -j$(nproc) debug_cohesive 2>&1 | tail -10 &&
  echo "" &&
  echo "============================================" &&
  echo "  RUNNING COHESIVE CELL DIAGNOSTIC" &&
  echo "============================================" &&
  echo "" &&
  ./debug_cohesive 2>&1
'
```

Read every line of output. Report the FULL output to me.

## Step 3: Check PETSc 3.20 APIs

Regardless of diagnostic results, also run:
```bash
docker run --rm fsrm-ci:local bash -c '
  echo "=== DMPlexConstructCohesiveCells ===" &&
  grep -rn "DMPlexConstructCohesiveCells" /opt/petsc-3.20.0/include/ &&
  echo "" &&
  echo "=== PetscDSSetBdResidual ===" &&
  grep -rn "PetscDSSetBdResidual" /opt/petsc-3.20.0/include/ &&
  echo "" &&
  echo "=== DMPlexGetCellType ===" &&
  grep -rn "DMPlexGetCellType" /opt/petsc-3.20.0/include/ &&
  echo "" &&
  echo "=== Cohesive examples ===" &&
  find /opt/petsc-3.20.0/src -name "*.c" | xargs grep -l "CohesiveCells" 2>/dev/null &&
  echo "" &&
  echo "=== DMSetRegionDS ===" &&
  grep -rn "DMSetRegionDS" /opt/petsc-3.20.0/include/ &&
  echo "" &&
  echo "=== TENSOR_PRISM types ===" &&
  grep -rn "PRISM_TENSOR" /opt/petsc-3.20.0/include/petscdmplextypes.h
'
```

This tells us what's available. Report the FULL output.

## Step 4: Based on results

If everything passes: great. Remove GTEST_SKIPs from fault tests, add real
assertions, verify Simulator with config/test_locked_fault.config.

If DMPlexConstructCohesiveCells fails: check the error message. Common issues:
- Label must mark faces (height-1 entities), not vertices
- Label values must be exactly 1
- The fault surface must separate the mesh into two connected components
- PETSc 3.20 may require specific mesh distribution state

If PetscFE/PetscDS fails on hybrid mesh: need DMSetRegionDS or separate DS
for bulk vs cohesive cells. Check PETSc examples for the pattern.

If PetscDSSetBdResidual means "boundary" not "cohesive": need to find the
correct API for cohesive cell residual registration.

## Also run while you're in Docker

Check if the existing fault tests pass or skip:
```bash
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local bash -c '
  ctest -R "FaultMesh\|LockedFault\|InjectionRupture" --output-on-failure -V 2>&1
'
```

This tells us the actual test output, not just pass/skip counts.

## Do NOT
- GTEST_SKIP over failures. Report them.
- Modify PetscFE callback files.
- Modify friction laws.
- Try to fix everything. Diagnose first.