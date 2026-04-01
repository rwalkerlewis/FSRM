# FSRM: Fix Cohesive Cell Insertion -- Skip Boundary Faces

Read CLAUDE.md. Build and test in Docker.

## The Problem

DMPlexConstructCohesiveCells fails with "Could not locate point X in split or
unsplit points" because the fault label includes faces on the mesh boundary.
PETSc cannot split a vertex that is shared between the fault and the domain
boundary because it is ambiguous which side of the fault the boundary belongs to.

## The Fix

In `src/numerics/FaultMeshManager.cpp`, function `createPlanarFaultLabel()`,
add a boundary check INSIDE the face loop. A boundary face has exactly 1
supporting cell. An interior face has exactly 2.

Find this code in the face labeling loop (around line 92):

```cpp
    for (PetscInt f = fStart; f < fEnd; ++f) {
        // Compute face centroid
        PetscReal centroid[3] = {0.0, 0.0, 0.0};
        PetscReal vol;
        ierr = DMPlexComputeCellGeometryFVM(dm, f, &vol, centroid, nullptr); CHKERRQ(ierr);

        double point[3] = {centroid[0], centroid[1], centroid[2]};

        // Check if face centroid lies on the fault plane
        if (isPointOnFaultPlane(point, normal, center, tolerance) &&
            isPointWithinFaultExtent(point, center, strike_dir, dip_dir, length, width)) {
            ierr = DMLabelSetValue(*fault_label, f, 1); CHKERRQ(ierr);
            num_marked++;
        }
    }
```

Add the boundary check as the FIRST thing inside the loop:

```cpp
    PetscInt num_skipped_boundary = 0;
    for (PetscInt f = fStart; f < fEnd; ++f) {
        // Skip boundary faces: they have only 1 supporting cell.
        // Interior faces have exactly 2. Labeling boundary faces causes
        // DMPlexConstructCohesiveCells to fail because it cannot determine
        // which side of the fault a boundary vertex belongs to.
        PetscInt support_size;
        ierr = DMPlexGetSupportSize(dm, f, &support_size); CHKERRQ(ierr);
        if (support_size < 2) {
            // Check if it would have been labeled (for diagnostic count)
            PetscReal centroid[3] = {0.0, 0.0, 0.0};
            PetscReal vol;
            ierr = DMPlexComputeCellGeometryFVM(dm, f, &vol, centroid, nullptr); CHKERRQ(ierr);
            double point[3] = {centroid[0], centroid[1], centroid[2]};
            if (isPointOnFaultPlane(point, normal, center, tolerance) &&
                isPointWithinFaultExtent(point, center, strike_dir, dip_dir, length, width)) {
                num_skipped_boundary++;
            }
            continue;
        }

        // Compute face centroid
        PetscReal centroid[3] = {0.0, 0.0, 0.0};
        PetscReal vol;
        ierr = DMPlexComputeCellGeometryFVM(dm, f, &vol, centroid, nullptr); CHKERRQ(ierr);

        double point[3] = {centroid[0], centroid[1], centroid[2]};

        // Check if face centroid lies on the fault plane
        if (isPointOnFaultPlane(point, normal, center, tolerance) &&
            isPointWithinFaultExtent(point, center, strike_dir, dip_dir, length, width)) {
            ierr = DMLabelSetValue(*fault_label, f, 1); CHKERRQ(ierr);
            num_marked++;
        }
    }
```

Update the diagnostic print at the end to include the boundary count:

```cpp
    if (rank_ == 0) {
        PetscPrintf(comm_, "FaultMeshManager: Marked %d interior faces on fault plane "
                   "(skipped %d boundary faces) "
                   "(strike=%.1f deg, dip=%.1f deg, center=[%.0f,%.0f,%.0f], "
                   "L=%.0f, W=%.0f, tol=%.3f)\n",
                   (int)num_marked, (int)num_skipped_boundary,
                   strike * 180.0 / M_PI, dip * 180.0 / M_PI,
                   center[0], center[1], center[2], length, width, tolerance);
    }

    if (num_marked == 0) {
        if (rank_ == 0) {
            PetscPrintf(comm_, "WARNING: No interior faces marked. The fault plane may not "
                       "intersect any interior faces of the mesh. Check fault geometry, "
                       "mesh resolution, and tolerance (%.3f).\n", tolerance);
        }
    }
```

That is the complete fix for createPlanarFaultLabel.

## Also: Fix the test mesh sizes

The tests use a 2x2x2 mesh with a fault at the center. On a 2x2x2 hex mesh,
the interior face at z=0.5 may have ALL its vertices on the boundary (the
fault plane passes through the exact center, and every vertex on that plane
is also on one of the x or y boundary faces).

There are two fixes for the tests:

1. Use a larger mesh (4x4x4 or 6x6x6) so the fault plane has interior vertices
   that are not on the boundary.

2. Or use a fault that does not extend to the mesh boundary (set length and
   width smaller than the mesh domain).

Apply fix 1 to:
- tests/unit/numerics/test_fault_mesh_manager.cpp: change faces from {2,2,2} to {4,4,4}
- tests/integration/test_locked_fault.cpp: already uses {4,4,4}, should be fine
- tests/integration/test_injection_rupture_chain.cpp: change faces from {2,2,2} to {4,4,4}

## Verify

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  mkdir -p build && cd build &&
  cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=OFF &&
  make -j$(nproc) &&
  echo "=== FAULT MESH MANAGER TEST ===" &&
  ./tests/run_unit_tests --gtest_filter=FaultMeshManagerTest.SplitMeshAlongFaultRequiresDMPlex 2>&1 | tail -20 &&
  echo "=== LOCKED FAULT TEST ===" &&
  ./tests/run_integration_tests --gtest_filter=LockedFaultTest.* 2>&1 | tail -20 &&
  echo "=== INJECTION CHAIN TEST ===" &&
  ./tests/run_integration_tests --gtest_filter=InjectionRuptureChainTest.* 2>&1 | tail -20 &&
  echo "=== ALL TESTS ===" &&
  ctest --output-on-failure 2>&1 | tail -20
'
```

The fault tests should now PASS instead of SKIP. Expected output from the mesh manager:
```
FaultMeshManager: Marked N interior faces on fault plane (skipped M boundary faces)
FaultMeshManager: Split complete: X cells after split (K cohesive cells inserted)
```

Where N > 0 and K > 0.

If N == 0 (no interior faces marked):
- The mesh is too coarse. Increase from 4x4x4 to 6x6x6.
- Or the tolerance is wrong. For a 4x4x4 mesh on [0,1]^3, cell size is 0.25,
  face centroids at z=0.5 exist, tolerance 0.15 should capture them.

If N > 0 but DMPlexConstructCohesiveCells still fails:
- Print the label contents before splitting: use DMLabelView(fault_label, PETSC_VIEWER_STDOUT_WORLD)
- Check if DMPlexLabelComplete is adding the correct closure points
- The issue might be that DMPlexLabelComplete marks boundary vertices that
  DMPlexLabelCohesiveComplete would have handled differently. Try using
  DMPlexLabelCohesiveComplete instead if it exists in PETSc 3.20.

## Do NOT
- Modify any PetscFE callback files
- Modify friction laws or CohesiveFaultKernel
- Add new features
- Delete tests