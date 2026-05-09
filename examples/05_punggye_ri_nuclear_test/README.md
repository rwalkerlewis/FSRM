# Example 05: Punggye-ri Nuclear Test -- showcase event

The Punggye-ri showcase. Modern, mountain-emplacement, motivates
axis-2 topography. Production-resolution counterpart for DPRK 2017
specifically lives at `examples/39_dprk_2017/`; this example uses a
generic Punggye-ri layered crustal model for the showcase
visualization track.

## Physics summary

| Property | Value |
|---|---|
| Geology context | Mt. Mantap, Punggye-ri, DPRK |
| Approximate yield | ~250 kt (matches DPRK 2017 estimate) |
| Approximate depth | 800 m below ground |
| Host rock | Fractured granite under rhyolite/tuff cap |
| Topography note | Mt. Mantap surface relief contaminates teleseismic surface waves |

## What the simulator does for Punggye-ri

This example exercises the layered material path with a generic Punggye-
ri crustal model. The "where this work goes next" presentation slide
anchors here: the residual axis-2 (topography) and axis-1d (3D far-
field) work is what would tighten the body-wave magnitude gate to the
+/-0.2 spec target. Pass-11 leaves these residuals open.

## Geology / velocity model

| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| Rhyolite / tuff cap | 0-200 | 3500 | 2020 | 2200 |
| Fractured granite | 200-1000 | 4500 | 2598 | 2500 |
| Competent granite | 1000-5000 | 5800 | 3349 | 2650 |

The explosion is placed at 800 m depth within the fractured granite
layer. Three surface seismometers record displacement waveforms.

## Configs

- `config.config`: pass-12 default (CI-quick variant from
  `punggye_ri_layered_quick.config`).
- `config_full.config`: full-resolution variant, longer simulation.
- For the actual DPRK 2017 historic event configuration with citation
  blocks (Voytan 2019, Wen 2018, Pabian 2018, Tian 2018), see
  `examples/39_dprk_2017/`.

## Showcase figures

```bash
cd examples/05_punggye_ri_nuclear_test
MPI_RANKS=8 ./run_showcase.sh
```

The pack covers: layered velocity model with source location, cavity-
formation radial profile, moment tensor history, all-station
seismogram grid, synthetic-vs-observed comparison (when DPRK 2017
cache is populated at `tools/waveform_vv/cache/dprk_2017/`), and
fidelity-tier comparison.

## V&V

- `Integration.PunggyeRiLayered`: full layered workflow.
- `Integration.LayeredElastostatics`: depth-based material assignment.
- `Integration.DPRK2017Comparison`: synthetic vs observed mb, pinned
  to DPRK 2017 yield estimate.
- `Integration.HistoricNuclear.DPRK2017`: pass-12 smoke test on
  example 39's production config.

## References

- Pabian, F. V. and Coblentz, D. D. (2018), "The 6th NK Nuclear
  Test", CNS Occasional Paper 38.
- Wen, L., et al. (2018), "Topographic and source backscattering
  effects from the 2017 NK test", GRL 45.
- Tian, X., et al. (2018), "Source parameters of the 2017 North
  Korean nuclear test", GJI 213.
- Voytan, D. P., et al. (2019), "Yield estimates for the six North
  Korean nuclear tests", GRL 46.
