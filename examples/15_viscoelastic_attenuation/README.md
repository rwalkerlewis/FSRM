# Example 15: Viscoelastic Attenuation

## Physics

This example demonstrates seismic wave attenuation using the generalized Maxwell body (GMB). The constitutive relation is:

    sigma(t) = C_relaxed : eps(t) + sum_i R_i(t)

where R_i are memory variables satisfying:

    dR_i/dt = -R_i / tau_i + delta_mu_i * deps/dt

The quality factor Q describes amplitude loss per wavelength: after traveling one wavelength, the amplitude decreases by `exp(-pi/Q)`.

## Setup

- Domain: 30 km x 30 km x 10 km
- Grid: 15 x 15 x 10 hexahedral cells
- Explosion source: 10 kt at center (15 km, 15 km, 5 km)
- Seismometers: at 2 km and 8 km from source (surface)
- Attenuation: Q_p = 200, Q_s = 100, 3 Maxwell mechanisms
- Frequency band: 0.1 -- 10 Hz (constant-Q approximation)

## Expected Results

The seismograms at both stations should show:

1. P-wave and S-wave arrivals
2. Amplitude at the far station reduced beyond pure geometric spreading
3. For f=1 Hz, delta_r=6 km, vs=3.35 km/s, Q_s=100:
   - Attenuation factor = exp(-pi * f * delta_r / (Q_s * vs))
   - = exp(-pi * 1 * 6000 / (100 * 3349))
   - = exp(-0.056)
   - = 0.945

So the far station should show approximately 5.5% additional amplitude loss compared to an elastic simulation.

## Running

```bash
docker run --rm -v $(pwd)/../..:/workspace -w /workspace/build fsrm-ci:local \
  ./fsrm -c ../examples/15_viscoelastic_attenuation/config.config \
  -ts_type alpha2 -pc_type lu -ksp_type preonly
```

## Visualization

After running, use the provided Python script to plot seismograms:

```bash
python3 ../../scripts/plot_seismograms.py output/
```
