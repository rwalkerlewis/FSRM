# Example 41: Pokhran II Shakti-I (1998-05-11)

The instrumented thermonuclear shot of the Indian Pokhran II series.
Pokhran II was a multi-shot operation across 1998-05-11 and
1998-05-13 with five named devices (Shakti-I through Shakti-V); this
example simulates Shakti-I, the announced thermonuclear test. The
example documents the multi-shot context but configures Shakti-I as
a single source.

## Event summary

| Property | Value |
|---|---|
| Date | 1998-05-11 |
| Site | Pokhran Test Range, Rajasthan, India (~27.07 N, 71.72 E) |
| Body-wave magnitude | mb 5.2 (Wallace 1998) |
| Indian announced yield | 43 kt thermonuclear |
| Western seismic estimate | 12-25 kt (Sykes & Wallace 1998; Wallace 1998) |
| Approximate depth | 210 m below surface |
| Host rock | Granitic gneiss (Marwar craton) under thin alluvium |

The yield discrepancy between the Indian announced yield (43 kt) and
the Western seismic estimate (12-25 kt mb 5.2) is itself a well-
documented case study in mb-to-yield conversion under non-NTS host
rock and is a useful frame for the axis-3 (layered-medium) and
axis-4 (regional refit) work.

## Geology / velocity model

| Depth | Layer | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| 0-50 m | Alluvium / weathered | 1900 | 900 | 2000 |
| 50-1500 m | Granitic gneiss host rock | 5600 | 3200 | 2700 |
| > 1500 m | Pre-Cambrian basement | 6300 | 3650 | 2800 |

Composite from Rodgers et al. (2002) regional crustal model and the
Marwar craton lithology summary in the Geological Survey of India
publications.

## Multi-shot context

Shakti-II was detonated near-simultaneously ~1 km away. Shakti-III,
IV, and V are reported as smaller fission tests on 1998-05-13. For
gauging close-in V&V the multi-shot superposition contaminates many
recordings. The far-field cardinal-grid stations in this config are
chosen at distances where the near-simultaneous shot's contribution
is smaller, but this example does not attempt to deconvolve the
multi-shot recording.

## Run

```bash
./run.sh
```

## V&V

Cached IRIS waveforms are not pre-populated; the
`tools/waveform_vv/cache/pokhran_ii_1998/` directory is a placeholder
awaiting an authored fetcher (no `scripts/fetch_pokhran_ii_*.py`
ships yet). The example smoke-tests the run pipeline without the
waveform-comparison gate.

## References

- Wallace, T. C. (1998), "The May 1998 India and Pakistan nuclear
  tests", Seism. Res. Lett. 69, 386-393.
- Sykes, L. R. and Wallace, T. C. (1998), "Re-evaluating yields of
  the Indian and Pakistani 1998 nuclear tests", GRL 25.
- Sikka, S. K., et al. (2000), "Update on Pokhran-II", BARC News
  Bulletin.
- Rodgers, A., et al. (2002), "An update on the lithospheric
  structure of the Indian shield", BSSA 92.
- USGS event service: us1998.05.11.
