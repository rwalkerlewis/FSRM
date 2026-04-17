# PyLith Compatibility Features

This document describes FSRM features that align with PyLith workflows and conventions,
enabling users familiar with PyLith to set up equivalent simulations in FSRM.

## Overview

FSRM implements several features commonly used in PyLith for crustal deformation and
fault mechanics:

| Feature | FSRM Status | PyLith Equivalent |
|---------|-------------|-------------------|
| Per-face Neumann traction BCs | Working | Neumann BC on faces |
| Prescribed kinematic fault | Working | FaultCohesiveKin |
| Time-dependent slip ramp | Working | TimeHistorySlipFn |
| Locked fault (transparent) | Working | FaultCohesiveKin (zero slip) |
| Roller (zero normal) BCs | Working | Dirichlet BC (single component) |
| Cohesive cell mesh splitting | Working | DMPlexConstructCohesiveCells |
| Lagrange multiplier fault traction | Working | FaultCohesiveKin constraint |

## Per-face Boundary Conditions

FSRM supports per-face boundary condition specification, similar to PyLith. Each face
of the bounding box can be independently configured.

### Configuration

```ini
[BOUNDARY_CONDITIONS]
x_min = roller
x_max = roller
y_min = roller
y_max = roller
z_min = fixed
z_max = traction
z_max_traction = 0.0, 0.0, -1.0e6
```

### Supported BC Types

| Type | Description | PyLith Equivalent |
|------|-------------|-------------------|
| `fixed` | All displacement components = 0 | Dirichlet (all components) |
| `roller` | Normal component = 0 | Dirichlet (single component) |
| `free` | Traction-free (natural BC) | Default (no BC) |
| `traction` | Applied traction vector | Neumann BC |
| `compression` | Normal compressive stress | Neumann BC (normal) |
| `absorbing` | Clayton-Engquist absorbing | AbsorbingDampers |

### Legacy Format

The per-face system also supports the legacy shorthand format:

```ini
[BOUNDARY_CONDITIONS]
bottom = fixed
sides = roller
top = traction
z_max_traction = 0.0, 0.0, -1.0e6
```

## Traction Boundary Conditions

Neumann (traction) BCs apply a specified stress vector to a boundary face. The traction
vector is specified in Cartesian coordinates (Pa).

### Implementation

FSRM uses manual FEM assembly (DMPlex point closure) rather than PetscDS boundary
residual callbacks. This is the same proven pattern used for explosion sources, fault
pressure, and cohesive constraints.

For each traction face:
1. Get boundary label IS (e.g., `boundary_z_max`)
2. Filter for depth `dim-1` faces with `support_size == 1` (boundary faces only)
3. Compute face area via `DMPlexComputeCellGeometryFVM`
4. Get face closure to find associated nodes
5. Distribute force equally: `F_node = traction * face_area / n_nodes`

### Analytical Verification

The traction BC is verified against the analytical uniaxial stress solution:

- Column height L = 100 m, E = 10 GPa, sigma = -1 MPa
- Expected: u_z(L) = sigma * L / E = -0.01 m
- Test: `Integration.TractionBC` (max displacement > 1e-4 m, < 1 m)

## Time-dependent Prescribed Slip

FSRM supports kinematic fault slip with a linear time ramp, similar to PyLith
`TimeHistorySlipFn`.

### Configuration

```ini
[FAULT]
mode = prescribed_slip
strike = 90.0
dip = 90.0
center_x = 0.5
center_y = 0.5
center_z = 0.5
length = 2.0
width = 2.0
slip_strike = 0.0
slip_dip = 0.0
slip_opening = 0.001
friction_coefficient = 0.6
slip_onset_time = 0.1
slip_rise_time = 0.5
```

### Slip Time Function

```
slip(t) = 0                                          for t < onset_time
slip(t) = slip_max * (t - onset_time) / rise_time    for onset_time <= t < onset_time + rise_time
slip(t) = slip_max                                   for t >= onset_time + rise_time
```

When `rise_time = 0`, the slip is applied instantaneously at `t >= onset_time`.

### Test Coverage

| Test | What It Verifies |
|------|-----------------|
| `Integration.TimeDependentSlip` | Slip ramp completes, nonzero solution |
| `Integration.DynamicRuptureSolve.PrescribedSlip` | Instant slip, displacement jump |
| `Physics.LockedFaultTransparency` | Zero slip gives slip < 5e-4 |

## Cohesive Cell Fault Representation

FSRM uses the same cohesive cell approach as PyLith:

1. `FaultMeshManager::markFaultFaces()` identifies interior faces on the fault plane
2. `DMPlexConstructCohesiveCells()` splits the mesh along the fault
3. Lagrange multiplier field enforces the fault constraint
4. PetscDS BdResidual callbacks on cohesive cells enforce the kinematic condition (PETSc 3.25+). Jacobian assembled manually via addCohesivePenaltyToJacobian.

### Fault Setup Pipeline

```
setupDM()            -- create hex/tet mesh (simplex if faults enabled)
setupFaultNetwork()  -- mark faces, split mesh, insert cohesive cells
labelBoundaries()    -- label boundary faces (after split)
setupFields()        -- add displacement + Lagrange multiplier fields
```

This ordering must not be changed. The fault mesh split must happen before boundary
labeling because it changes the mesh topology.

## Examples

| # | Directory | Config | Feature |
|---|-----------|--------|---------|
| 07 | `examples/07_traction_bc/` | `traction_bc.config` | Per-face traction BC |
| 08 | `examples/08_time_dependent_slip/` | `time_dependent_slip.config` | Slip ramp |

## Differences from PyLith

| Aspect | FSRM | PyLith |
|--------|------|--------|
| Configuration | INI-style `.config` files | `.cfg` files with Pyre framework |
| Mesh | Built-in hex + Gmsh import | Gmsh, CUBIT, LaGriT |
| Time function | Linear ramp only | Step, constant rate, Liu cosine, time history |
| Fault friction | Coulomb (manual assembly) | Rate-and-state, slip weakening |
| Output | HDF5, VTK, SAC | HDF5/Xdmf |
| Parallelism | MPI + optional GPU (PETSc CUDA) | MPI only |
