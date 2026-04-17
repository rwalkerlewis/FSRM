# Example 07: Traction Boundary Condition

## Physics

Applies a uniform downward compressive traction of 1 MPa on the top face of a 3D box,
with fixed bottom and roller (zero normal displacement) sides. This is a uniaxial stress
test that exercises the per-face boundary condition system.

### Analytical Solution

For a homogeneous isotropic elastic box under uniaxial stress:
- Applied traction: sigma_zz = -1 MPa on z_max face
- Young modulus: E = 10 GPa
- Vertical displacement: u_z(z) = sigma_zz * z / E = -1.0e-4 * z
- At z = 100 m: u_z = -0.01 m
- Lateral displacement: u_x = u_y = 0 (constrained by rollers)

## Config

`config/examples/traction_bc.config`

## Boundary Conditions

| Face   | Type    | Description                    |
|--------|---------|--------------------------------|
| z_min  | fixed   | All displacements = 0          |
| z_max  | traction| t = (0, 0, -1 MPa)            |
| x_min  | roller  | u_x = 0                        |
| x_max  | roller  | u_x = 0                        |
| y_min  | roller  | u_y = 0                        |
| y_max  | roller  | u_y = 0                        |

## Running

```bash
./run.sh
```

Or manually:

```bash
cd build
./fsrm -c ../config/examples/traction_bc.config
```

## Expected Output

- `output/solution.h5`: Solution with displacement field
- `output_SUMMARY.txt`: Run summary

## Verified By

- `Integration.TractionBC` test
