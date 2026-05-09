# Example 43: Mururoa "Xouthos" 1996 (final French underground test)

The last French nuclear test, signed off shortly before France
ratified the CTBT. The 1995-96 final test campaign was internationally
controversial and was one of the best-monitored Pacific atoll test
series. With this example, the FSRM example catalog now documents all
five weapon-state programs (US, USSR / Russia, China, India, France).

## Event summary

| Property | Value |
|---|---|
| Date | 1996-01-27 |
| Site | Mururoa Atoll, French Polynesia (~21.84 S, 138.96 W) |
| Yield estimate | ~120 kt |
| Approximate depth | 1100 m below atoll surface |
| Host rock | Cretaceous basaltic basement under volcanic tuff and coral cap |

## Geology / velocity model

Atoll geology is the distinctive feature: thin coral cap over a
subsiding volcanic edifice. The cavity expands into competent basalt;
the overlying volcanic / coral structure transmits the seismic signal
to the seafloor and the surrounding ocean.

| Depth | Layer | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| 0-300 m | Coral reef cap | 2000 | 1100 | 2100 |
| 300-900 m | Volcanic tuff / breccia | 3500 | 1900 | 2400 |
| > 900 m | Cretaceous basaltic basement | 5800 | 3300 | 2800 |

Composite from Massin et al. (2003) active-seismic atoll crustal model.

## Limitations

This config is a 3-D solid-Earth approximation that **omits the ocean
halfspace** (no fluid coupling). Atoll seismic propagation depends
significantly on water-column scattering and seafloor coupling that
require fluid-solid interface conditions; coupling the ocean response
is named axis-5 future work.

The smoke-test config places the source at 650 m depth (within the
volcanic-tuff layer) rather than the actual 1100 m depth (within the
basaltic basement). The deeper, stiffer-rock placement produces a
larger elastic-impedance jump at the basement interface that
exceeds the SNES nonlinear-solver's stability margin in the implicit
TSALPHA2 scheme on this CI-tractable mesh. Production-resolution
runs that recover the deeper placement are a future axis-5 follow-up
alongside the ocean-coupling work.

## Run

```bash
./run.sh
```

## V&V

Cached IRIS waveforms are not pre-populated;
`tools/waveform_vv/cache/mururoa/` is a placeholder. The 1995-96
series falls inside the IRIS open-data window and a fetcher is
authorable as a future contribution.

## References

- IFRC / IAEA (1998), "The radiological situation at Mururoa and
  Fangataufa atolls", IAEA Tech-Doc 1066.
- Massin, F., et al. (2003), "Crustal structure of Mururoa Atoll
  from active seismic data", Pure Appl. Geophys. 160.
- Adushkin, V. V. and Spivak, A. A. (2015), "Underground Explosions"
  (English transl.), Mir Publishers.
- USGS event service: us1996.01.27.
