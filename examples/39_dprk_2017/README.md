# Example 39: DPRK Sixth Underground Nuclear Test (2017-09-03)

Modern, well-instrumented underground nuclear test: the most-recorded
event of the post-CTBTO era and the natural V&V anchor for any source-
physics simulation that targets contemporary monitoring scenarios. Mt.
Mantap topography contaminates teleseismic surface-wave recordings,
which is itself a useful test of how an axis-2 topography pass would
change the result.

## Event summary

| Property | Value |
|---|---|
| Date | 2017-09-03 |
| Site | Punggye-ri, Mt. Mantap, Kilju, North Hamgyong, DPRK |
| Body-wave magnitude | mb 6.3 (USGS) |
| Yield estimate | ~100-370 kt (analyses vary); ~250 kt central |
| Announced | "Thermonuclear" device |
| Approximate depth | 600-800 m below Mt. Mantap summit |
| Host rock | Competent granite under volcanic tuff overburden |

## Geology / velocity model

| Depth below summit | Layer | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| 0-300 m | Volcanic tuff / weathered overburden | 2400 | 1300 | 2200 |
| 300-3000 m | Competent granite host rock | 5800 | 3300 | 2700 |
| > 3000 m | Pre-Cambrian metamorphic basement | 6300 | 3650 | 2800 |

References: Pabian & Coblentz (2018) for geology, Wen et al. (2018)
for crustal model, Tian et al. (2018) for source parameters.

## Run

```bash
./run.sh                 # default 4 ranks
MPI_RANKS=8 ./run.sh     # production
```

Output lands in `examples/39_dprk_2017/output/`. SAC seismograms at
five stations (epicentral, near, and far in the cardinal-grid
directions) plus the standard `near_field_history.csv` and HDF5
wavefield.

## Variants

- `config.config`: pass-12 default. Pass-7 closed-form near-field
  source physics (TILLOTSON cavity EOS, ZELDOVICH_RAIZER radiation),
  3-layer geology, 5 seismometers, 8 s simulation.
- `config_quick.config`: pre-pass-12 minimal smoke variant. Single-
  layer granite, 5 s simulation, kept as a regression anchor.

## V&V

The IRIS waveform-comparison cache for DPRK 2017 lives at
`tools/waveform_vv/cache/dprk_2017/` and is populated by
`scripts/fetch_dprk_2017_waveforms.py`. Pass-12 ships a smoke-level
integration test (`Integration.HistoricNuclear.DPRK2017`) that
verifies the simulation runs end-to-end and produces SAC output. Gate-
authoring (synthetic-vs-observed cross-correlation under cached
station data) is a future physics-pass task.

`Integration.MPI.DPRK2017SerialVsParallelEquivalence` checks the
serial vs 4-rank SAC output for numerical equivalence (peak
amplitude relative difference under 1 %, cross-correlation > 0.99 in
the 0.5-5 Hz band).

## Why this event matters for the program

- **Modern monitoring scenario**: contemporary IRIS coverage and
  open-data tooling. This is the best public-data anchor for any
  forward-looking V&V work.
- **Topography contamination**: Mt. Mantap surface relief is large
  enough to contaminate teleseismic surface-wave recordings. Until
  axis-2 topography fidelity ships, body-wave magnitude and near-
  field physics are the part of the response we expect to capture
  cleanly.
- **Yield uncertainty**: bracketed academic estimates 100-370 kt
  illustrate how cavity-radius / radiation-phase fidelity propagates
  to mb. The example documents the mid-range yield as the default
  but the config knob is a single line edit.

## References

- Voytan, D. P., et al. (2019), "Yield estimates for the six North
  Korean nuclear tests from teleseismic P-wave modeling and
  intercorrelation of P- and Rg-wave receiver functions",
  GRL 46, 4137-4146.
- Wen, L., et al. (2018), "Topographic and source backscattering
  effects from the 2017 NK test", GRL 45.
- Pabian, F. V. and Coblentz, D. D. (2018), "The 6th North Korean
  Nuclear Test", James Martin Center for Nonproliferation Studies,
  CNS Occasional Paper 38.
- Tian, X., et al. (2018), "Source parameters of the 2017 North
  Korean nuclear test from regional and teleseismic moment tensor
  analysis", GJI 213.
- USGS event page: https://earthquake.usgs.gov/earthquakes/eventpage/us2000aert
