# Example 01: Uniaxial Compression

## Physics

Quasi-static linear elastostatics using PetscFE pointwise callbacks (f0, f1, g3).

A 10x10x100 m column of rock (E=10 GPa, nu=0.25, rho=2650 kg/m^3) is compressed
by applying a displacement boundary condition on the top face while the bottom face
is held fixed. The SNES (nonlinear solver) should converge in 1 iteration since
the problem is linear.

## Config

Uses `config/examples/uniaxial_compression.config`.

## Expected Output

- `output/solution.h5` -- HDF5 solution file with displacement field
- `output_SUMMARY.txt` -- simulation summary

The displacement field should show linear variation in z (u_z varies from 0 at
bottom to the applied value at top), with uniform axial strain eps_zz.

## PETSc Options

- `-snes_monitor` to see SNES convergence (should be 1 iteration)
- `-ksp_monitor` to see KSP convergence
- `-log_view` for performance summary

## Running

```bash
./run.sh
```

Or manually:
```bash
cd /path/to/build
./fsrm -c ../config/examples/uniaxial_compression.config
```

## Verified By

- `Physics.ElastostaticsPatch` -- patch test convergence
- `Physics.LithostaticStress` -- lithostatic stress column
- `Integration.FullSimulation` -- full pipeline test
