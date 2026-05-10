# Example 42: Joint Verification Experiment, Semipalatinsk leg (1988)

The best-documented non-NTS event of the late Cold War. JVE was a
paired bilateral US-Soviet experiment with on-site monitoring at both
shots. Bilateral access enabled accurate yield estimation by both
parties using independent techniques (CORRTEX downhole instruments
host-side, plus open seismic recording visiting-side).

This example simulates the Soviet leg, "Shagan" at Semipalatinsk on
1988-09-14. The US leg ("Kearsarge" at NTS Pahute Mesa, 1988-08-17,
~150 kt) is implicitly covered by examples 13 and adjacent NTS
examples.

## Event summary

| Property | Value |
|---|---|
| Date | 1988-09-14 (Soviet leg) |
| Site | Shagan River testing area, Semipalatinsk Test Site |
| Coordinates | ~50.05 N, 78.99 E |
| Body-wave magnitude | mb 6.1 |
| Yield (CORRTEX) | ~115 kt |
| Approximate depth | 650 m |
| Host rock | Pre-Cambrian granitic basement under hornfels |

## Geology / velocity model

| Depth | Layer | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| 0-300 m | Weathered surface / shale | 2400 | 1300 | 2200 |
| 300-1500 m | Devonian sediment / hornfels | 4800 | 2700 | 2600 |
| > 1500 m | Pre-Cambrian granitic basement | 6100 | 3500 | 2750 |

References: Bocharov et al. (1989) for the layered model; Sykes &
Ekstrom (1989) for source-parameter analogy with other Shagan events;
Vergino & Mensing (1990) for JVE-specific Lg-coda magnitudes.

## Run

```bash
./run.sh                 # default 4 ranks
MPI_RANKS=8 ./run.sh     # production
```

## V&V

`tools/waveform_vv/cache/jve_1988/` is a placeholder. No fetcher
script ships pass-12; existing IRIS coverage of Shagan 1988 is
regional and would need a custom fetcher for the specific date and
station list.

## References

- Vergino, E. S. and Mensing, R. W. (1990), "Yield estimation using
  regional Lg-coda magnitudes", BSSA 80, 656-674.
- Sykes, L. R. and Ekstrom, G. (1989), "Yields and source mechanisms
  of Soviet underground tests from 1978-1985", GJI 99.
- Bocharov, V. S., et al. (1989), "Estimating yields of Soviet
  underground tests at Semipalatinsk", Soviet Physics-Doklady.
- USGS event service: us1988.09.14.
