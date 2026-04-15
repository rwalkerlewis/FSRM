# Example 04: Locked Fault

## Physics

Quasi-static compression of a domain containing a locked cohesive fault.

The mesh is split along a vertical fault at x=0.5 using the PyLith cohesive
cell workflow (DMPlexCreateSubmesh, DMPlexLabelCohesiveComplete,
DMPlexConstructCohesiveCells). The locked fault mode enforces displacement
continuity (u+ = u-) across the fault interface using manual cohesive
constraint assembly.

Under uniaxial compression, the locked fault should be transparent -- the
stress and displacement fields should be identical to a domain without a fault.

Uses:
- FaultMeshManager for cohesive cell insertion
- CohesiveFaultKernel for constraint assembly
- PetscDS BdResidual callbacks for fault constraint (PETSc 3.25+)
- Analytical interface Jacobian for locked mode

## Config

Uses `config/examples/locked_fault_compression.config`.

## Expected Output

- `output/solution.h5` -- HDF5 solution with displacement field

The fault slip should be negligible (< 5e-4 relative to applied displacement).

## Running

```bash
./run.sh
```

Or manually:
```bash
cd /path/to/build
./fsrm -c ../config/examples/locked_fault_compression.config
```

## Verified By

- `Integration.DynamicRuptureSolve.LockedQuasiStatic` -- locked solve through TSSolve
- `Physics.LockedFaultTransparency` -- fault slip < 5e-4
- `Functional.DynamicRuptureSetup.MeshSplitting` -- cohesive cell insertion
