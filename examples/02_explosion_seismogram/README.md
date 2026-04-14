# Example 02: Explosion Seismogram

## Physics

Elastodynamic wave propagation from an underground nuclear explosion source.

A 10 kt explosion at 300 m depth generates seismic waves that propagate through
a homogeneous elastic halfspace (E=50 GPa, nu=0.25, rho=2700 kg/m^3). Three
surface seismometers at 1, 3, and 5 km record velocity waveforms in SAC format.

Uses:
- Mueller-Murphy source model for moment tensor
- TSALPHA2 implicit time integration (generalized-alpha, second-order)
- Clayton-Engquist absorbing boundary conditions
- SeismometerNetwork for waveform sampling via DMInterpolation

## Config

Uses `config/examples/explosion_seismogram.config`.

## Expected Output

- `output/solution.h5` -- HDF5 solution file with displacement snapshots
- `output/seismograms/*.sac` -- SAC waveform files for each station/component

The seismograms should show:
- P-wave arrival times proportional to distance
- Amplitude decreasing with distance (geometric spreading)
- Clean waveforms with minimal boundary reflections

## PETSc Options

- `-ts_monitor` to see time stepping progress
- `-ts_max_steps N` to limit number of time steps
- `-log_view` for performance summary

## Running

```bash
./run.sh
```

Or manually:
```bash
cd /path/to/build
./fsrm -c ../config/examples/explosion_seismogram.config
```

## Visualization

```bash
python3 scripts/plot_seismograms.py output/seismograms/ --output figures/explosion_seismograms.png
```

## Verified By

- `Integration.ExplosionSeismogram` -- pipeline completion, amplitude vs distance
- `Physics.MomentTensorSource` -- FEM equivalent nodal forces
- `Physics.AbsorbingBC` -- >99% energy absorption
