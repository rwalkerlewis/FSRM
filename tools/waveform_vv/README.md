# FSRM Waveform V&V (Pass-8)

This directory holds the offline tools and cached data for the
pass-8 IRIS waveform V&V infrastructure. It is consumed by the
C++ tests under the `iris_validation` CTest label.

## Layout

```
tools/waveform_vv/
  events.yaml         # Event manifest (origin time, networks, stations)
  refresh.py          # ObsPy-based IRIS DMC refresh tool (offline)
  cache/              # Cached SAC files committed to the repo
    Salmon1964/
      <STA>.<CH>.sac   # Pre-processed traces, displacement, m
      metadata.yaml    # Provenance: query, fetch timestamp
    Chagan1965/
    PokhranI_1974/
```

## Cache state

The cache ships empty by default. Tests under the
`iris_validation` CTest label `GTEST_SKIP` cleanly with an explicit
message when the per-event cache directory is missing or empty.

To populate the cache, run the refresh tool **outside** the FSRM
Docker image (the image does not ship ObsPy):

```bash
python -m venv .venv
source .venv/bin/activate
pip install obspy>=1.4 pyyaml
python tools/waveform_vv/refresh.py --event Salmon1964
python tools/waveform_vv/refresh.py --event Chagan1965
python tools/waveform_vv/refresh.py --event PokhranI_1974
# or, all at once:
python tools/waveform_vv/refresh.py --all
```

After refresh, the cache holds one SAC file per
(event, station, channel) plus a `metadata.yaml` recording the
IRIS query parameters and download timestamp. Commit the cache
and re-run the CI suite (`ctest -L iris_validation`); the
previously-skipped tests now exercise the comparison library.

## Why the cache is committed

The pass-8 V&V is reproducible: the same set of stations, the same
window around origin time, the same instrument-response correction
must be applied for every CI run. Re-fetching from IRIS each run
is fragile (network outages, IRIS DMC maintenance, station
re-equipment) and would make the gates non-deterministic. Caching
is the correct tradeoff for a CI-runnable validation gate.

## Adding a new anchor event

1. Append an entry to `events.yaml` with origin time, networks,
   stations, and channels. Cite the published reference for the
   yield, depth, medium, and any measured-value gates (cavity
   radius, mb, etc).
2. Run `refresh.py --event <NewEvent>` to populate the cache.
3. Add a corresponding `IRISValidationTest.<NewEvent>...` test in
   `tests/integration/test_iris_validation.cpp`.
4. Register the test under the `iris_validation` CTest label.

## Pre-FDSN era (1964-1976) limitations

Salmon 1964, Chagan 1965, and Pokhran I 1974 predate the FDSN
standard and the modern IRIS DMC archive. Many stations from those
eras have analog WWSSN records that are not digitized, or have
been re-equipped with no provenance link to the original
instrument. Expect partial coverage. The C++ tests handle
empty / partial caches by GTEST_SKIP'ing per-station gates with
explicit messages so a partial fetch does not silently pass.

## References

- Goldstein, P. and Snoke, A. (2005), "SAC Availability for the
  IRIS Community", DMS Electronic Newsletter VII(1).
- Springer, D. L., Healy, J. H., Mickey, W. V. (1968), "Seismic
  source mechanism for the Salmon nuclear explosion in salt",
  Geophysics 33(4), pp. 581-588.
- Healy, J. H. (1971), "Seismic source mechanism studies of the
  Salmon and Sterling events", USGS Professional Paper 750-D.
- Patton, H. J. (1991), "Seismic moment estimation and the scaling
  of the long-period source spectrum at the Salmon site",
  BSSA 81(4), pp. 1376-1404.
- Adushkin, V. V. and Spivak, A. A. (2003), "Underground Explosions
  and Seismic Activity", Springer.
- Sykes, L. R. (1996), "Decade of seismic networks for nuclear
  test monitoring", Reviews of Geophysics 34(4).
