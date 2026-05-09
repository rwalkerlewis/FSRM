# FSRM Waveform V&V (Pass-8)

This document describes the IRIS waveform validation infrastructure
introduced in pass-8 of the historic-nuclear track. The infrastructure
is anchored on Salmon 1964 (Project Dribble, Mississippi salt dome)
with two cross-validation events (Chagan 1965, Pokhran I 1974).

The standing scope and current state of pass-8 deliverables is in
`docs/HISTORIC_NUCLEAR_FIDELITY.md` section 4g; the forward-looking
view is in `docs/HISTORIC_NUCLEAR_ROADMAP.md`. This document is the
operational guide for using the V&V infrastructure.

## What pass-8 ships

Three layers, each individually testable:

1. **Production SAC reader** (`include/io/SACReader.hpp`,
   `src/io/SACReader.cpp`). Full canonical-keyword extraction,
   bilateral endian support via the NVHDR sentinel, pre-processing
   helpers (resample / window / taper / demean+detrend). Unit
   tested by `Unit.SACReader` (7 cases).

2. **Comparison metrics library** (`include/diagnostics/WaveformComparison.hpp`,
   `src/diagnostics/WaveformComparison.cpp`). Five free functions
   (peakAmplitudeRatio, spectralAmplitudeRatio, dominantFrequency,
   envelopeMisfit, crossCorrelation) plus an aggregate
   `compareWaveforms`. Hand-rolled radix-2 Cooley-Tukey FFT
   (no external dep). Unit tested by `Unit.WaveformComparison`
   (8 cases).

3. **iris_validation CTest label**
   (`tests/integration/test_iris_validation.cpp`). Five integration
   gates anchored on three historic events. Each gate is cache-aware:
   `GTEST_SKIP`s cleanly when the relevant
   `tools/waveform_vv/cache/<EventName>/` directory is empty, with
   an explicit message pointing the user at refresh.py.

## Cache layout

```
tools/waveform_vv/
  events.yaml        # Event manifest (origin time, networks, stations)
  refresh.py         # ObsPy-based IRIS DMC refresh tool (offline)
  README.md
  cache/
    Salmon1964/
      <STA>.<CH>.sac      # one SAC per (station, channel)
      metadata.yaml       # provenance: query parameters, fetch timestamp
    Chagan1965/
    PokhranI_1974/
```

The cache ships empty by default and is committed to git when
populated. Testing should follow:

```bash
# 1. Run outside the FSRM Docker image (image does not ship ObsPy):
python -m venv .venv
source .venv/bin/activate
pip install obspy>=1.4 pyyaml

# 2. Refresh per-event:
python tools/waveform_vv/refresh.py --event Salmon1964
python tools/waveform_vv/refresh.py --event Chagan1965
python tools/waveform_vv/refresh.py --event PokhranI_1974

# 3. Run the V&V gates inside the FSRM Docker image:
docker run --rm -v $(pwd):/workspace -w /workspace/build fsrm-ci:local \
    ctest -L iris_validation --output-on-failure
```

## Comparison metrics

### `peakAmplitudeRatio(obs, syn)`

Returns `|peak(syn)| / |peak(obs)|`. The peak is the maximum
absolute sample value in the trace; window the trace before
calling if a sub-region peak is wanted.

### `spectralAmplitudeRatio(obs, syn, fmin, fmax)`

Both traces are resampled onto the smaller dt of the two and
zero-padded to the next power of two. The unfiltered amplitude
spectrum is integrated over `[fmin, fmax]` and the ratio is
returned. Used for the pass-8 `metric_freq_band_hz` 0.5-5 Hz band
on regional / teleseismic arrivals.

### `dominantFrequency(trace)`

Argmax of the smoothed amplitude spectrum (5-point boxcar). For
pass-8 cross-validation, dominant-frequency match within factor 3
gates the source corner-frequency physics.

### `envelopeMisfit(obs, syn)`

L2 norm of the difference between Hilbert-transform analytic-signal
envelopes, normalized by the L2 of the observed envelope. Robust
to small phase shifts (the envelope is invariant under symmetric
phase shifts) so it is a useful complement to crossCorrelation.

### `crossCorrelation(obs, syn, max_lag, &out_lag)`

Normalized cross-correlation, maximized over lags in
`[-max_lag, +max_lag]`. Both traces are demeaned before the CC.
The lag at the maximum is returned via the out-parameter. Used for
arrival-time consistency checks; pass-8 does NOT gate on full-
waveform CC (that is axis-5 work in a future pass).

## Adding a new anchor event

1. Append an entry to `tools/waveform_vv/events.yaml` with origin
   time, event coordinates, yield, medium, published mb, measured
   cavity radius (where available), and per-station channel
   lists. Cite the published reference for every measured value.
2. Run `refresh.py --event <NewEvent>` outside the Docker image.
3. Add an `IRISValidationTest.<NewEvent>...` test in
   `tests/integration/test_iris_validation.cpp`. Check the cache
   directory; `GTEST_SKIP` cleanly if empty.
4. Register the test under the `iris_validation` CTest label and
   add it to the `LABELS "iris_validation;integration"` block in
   `tests/CMakeLists.txt`.

## Pre-FDSN era coverage

Salmon 1964, Chagan 1965, and Pokhran I 1974 predate the FDSN
standard and the modern IRIS DMC archive. Many stations from those
eras have analog WWSSN records that are not digitized, or were
re-equipped with no provenance link. The `notes:` field of each
manifest entry documents known gaps.

When the IRIS query returns empty for a (station, channel) pair,
`refresh.py` records the failure in `metadata.yaml` and proceeds.
The C++ test reads the metadata.yaml on cache load (future work)
and `GTEST_SKIP`s the per-station gate with an explicit message.

## Why the cache is committed

Re-fetching from IRIS at every CI run is fragile: network outages,
IRIS DMC maintenance windows, and station re-equipment break
reproducibility. The cache pins the exact set of bytes the V&V
gates run against. The pass-8 spec rule applies: do not adjust
opacity coefficients or any physics parameter to make Salmon pass.
The cache is the immutable observation.

## Pass-8 V&V gates summary

| Gate | Event | Target | Source |
|------|-------|--------|--------|
| Salmon1964CavityRadius | Salmon 1964 | 17.4 m within factor 5 | Springer 1968; Patton 1991 |
| Salmon1964FreeFieldPeakVelocity | Salmon 1964 | informational (skip) | Healy 1971 |
| Salmon1964FarFieldMb | Salmon 1964 | 4.9 +/- 0.4 mb | Murphy 1981; Stump 1994 |
| Chagan1965CrossValidation | Chagan 1965 | cavity 75 m factor 10; mb 6.0 +/- 0.5 | Adushkin & Spivak 2003 |
| PokhranI1974CrossValidation | Pokhran I 1974 | mb 4.9 +/- 0.4 | Sykes 1998 |

Pass-8 V&V does not gate on full-waveform cross-correlation (that
depends on far-field propagation fidelity, axis 5+) nor on
free-field velocity at the Healy 1971 ranges (that depends on 3D
source physics, axis-1b). Both are named pass-9+ work.
