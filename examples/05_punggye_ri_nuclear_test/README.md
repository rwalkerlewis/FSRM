# Example 05: Punggye-ri Nuclear Test

## Physics

Three-layer elastic model of the Punggye-ri nuclear test site (DPRK, 2017 test,
estimated 250 kt). Elastodynamic wave propagation with layered geology, explosion
source, absorbing boundaries, and surface seismometer recording.

Layer structure (based on published geology of Mt. Mantap):
- **Rhyolite/tuff cap** (0-200 m depth): vp=3500, vs=2020, rho=2200
- **Fractured granite** (200-1000 m depth): vp=4500, vs=2598, rho=2500
- **Competent granite** (1000-5000 m depth): vp=5800, vs=3349, rho=2650

The explosion is placed at 800 m depth within the fractured granite layer.
Three surface seismometers record displacement waveforms.

Uses:
- Depth-based material layering via auxiliary fields
- Mueller-Murphy source model (250 kt)
- TSALPHA2 implicit time integration
- Clayton-Engquist absorbing BCs on 5 faces
- SeismometerNetwork with SAC output

## Config

Uses `config/examples/punggye_ri_layered.config`.

## Expected Output

- `output/solution.h5` -- HDF5 solution snapshots
- `output/punggye_ri/*.sac` -- SAC seismogram files

The seismograms should show:
- Complex waveforms from layer reflections and refractions
- P-wave followed by slower surface waves
- Amplitude decreasing with distance

## Running

```bash
./run.sh
```

Or manually (this example takes longer due to wave propagation):
```bash
cd /path/to/build
./fsrm -c ../config/examples/punggye_ri_layered.config
```

## Visualization

```bash
python3 scripts/plot_seismograms.py output/punggye_ri/ --output figures/punggye_ri_seismograms.png
```

## Verified By

- `Integration.PunggyeRiLayered` -- full layered workflow
- `Integration.LayeredElastostatics` -- depth-based material assignment
- `Integration.DPRK2017Comparison` -- synthetic vs observed mb
