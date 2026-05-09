# Example 23: Milrow (1969)

## Physics
Elastodynamic wave propagation through the Aleutian-arc andesite
stratigraphy at Amchitka Island, AK. Models the ~1 Mt Milrow event of
1969-10-02 at 1218 m depth -- the engineering-validation precursor to
Cannikin 1971.

Medium choice: GENERIC (no andesite coefficient in the cavity table;
see Long Shot 1965 README for the rationale).

Features used:
- Layered heterogeneous material (aux field assignment)
- Underground nuclear explosion source (Mueller-Murphy, GENERIC medium)
- Absorbing boundary conditions (Clayton-Engquist)
- Seismometer network with SAC output

## Geology
Same Aleutian stratigraphy as Long Shot 1965, extended to 5 km depth:

| Layer | Depth (m) | Vp (m/s) | Vs (m/s) | Density (kg/m3) |
|-------|-----------|----------|----------|-----------------|
| Tundra / weathered overburden | 0-50    | 1500 | 750  | 1800 |
| Weathered andesite           | 50-300  | 2800 | 1500 | 2300 |
| Andesite breccia             | 300-1500 | 4000 | 2200 | 2500 |
| Competent andesite           | 1500-5000 | 4800 | 2700 | 2700 |

Reference: King, Foster, & Bingham 1972 (Amchitka).

## Config
Uses `config/examples/milrow_1969.config`.

## Expected Output
- `output/milrow_1969/*.SAC` -- synthetic seismograms at 3 stations

## Running
```bash
./run.sh
```

## Verified By
- `Integration.HistoricNuclear.Milrow1969` -- megaton andesite + explosion + SAC output
