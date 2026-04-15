# Example 17: Per-Cell Material from Velocity Model

## Overview

This example demonstrates reading a 3D structured binary velocity model file
and assigning per-cell material properties (lambda, mu, rho) to the unstructured
FEM mesh via trilinear interpolation.

## Physics

- Elastodynamics with heterogeneous material
- Explosion source at the center of the domain
- Two-layer velocity structure: soft sediment over hard basement

## Velocity Model Format

Binary file (little-endian):

```
Header (44 bytes):
  nx, ny, nz       (3 x int32)
  x_min, x_max     (2 x float32)
  y_min, y_max     (2 x float32)
  z_min, z_max     (2 x float32)
  padding           (2 x float32)
Data (nx * ny * nz * 3 * float32):
  For each (ix, iy, iz) in C order:
    Vp, Vs, rho    (3 x float32)
```

## Running

```bash
# Generate the velocity model
python3 ../../scripts/generate_velocity_model.py velocity_model.bin \
  --nx 10 --ny 10 --nz 10 --xmax 5000 --ymax 5000 --zmax 5000

# Run the simulation
cd ../../build
./fsrm -c ../examples/17_velocity_model/config.config
```

## Expected Results

The P-wave velocity in the soft layer (Vp=3000 m/s) is half that of the hard
layer (Vp=6000 m/s), so waves should travel slower through the upper half.
The impedance contrast at the layer boundary will produce reflected waves.

## Material Conversion

Velocities are converted to Lame parameters:
- mu = rho * Vs^2
- lambda = rho * (Vp^2 - 2 * Vs^2)
