# DPRK 2017 IRIS waveform cache

Populated by `scripts/fetch_dprk_2017_waveforms.py` (ObsPy required).

## Event metadata

- Event ID: USGS us2000aert
- Origin time (UTC): 2017-09-03 03:30:01
- Hypocenter: 41.343 N, 129.036 E (Mt. Mantap, Punggye-ri)
- Magnitude: mb 6.3 (USGS)

## Cached stations

The fetcher populates this directory with SAC files from a small set of
representative IRIS stations on the Korean Peninsula and surrounding
regions. The exact station list is defined in
`scripts/fetch_dprk_2017_waveforms.py` and includes regional
broadband and short-period instruments; teleseismic stations are not
the right comparison for source-physics V&V at this scope.

## Layout

```
dprk_2017/
    inventory.xml          # ObsPy StationXML for the cached stations
    <NET>.<STA>.<LOC>.<CHA>.sac
    ...
```

## Time window

Default fetch window is from origin -10 s to origin +180 s. Override
in the fetcher.

## V&V usage

Comparison gates are applied via the `[WAVEFORM_VV]` config block. See
`docs/WAVEFORM_VV.md` for the format and the gate specifications.
