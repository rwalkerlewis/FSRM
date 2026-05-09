# Example 40: Lop Nor 1996 (Final Chinese Underground Nuclear Test)

The last Chinese nuclear test before signing the CTBT in September
1996. The basin-context emplacement gives less topography
contamination than the 1970s tunnel sites in the Tian Shan; the open
desert keeps surface-wave excitation cleaner, which makes this a
good far-field-amplitude V&V anchor.

## Event summary

| Property | Value |
|---|---|
| Date | 1996-07-29 |
| Site | Lop Nor Test Site, Xinjiang, PRC (~41.66 N, 88.41 E) |
| Body-wave magnitude | mb 5.0 (USGS / IDC bulletin) |
| Yield estimate | ~5 kt central (range 1-10 kt; Sykes 2002) |
| Approximate depth | 700-900 m below surface |
| Host rock | Mesozoic / weathered granite under Tertiary sediment |

## Geology / velocity model

| Depth | Layer | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---:|---:|---:|
| 0-200 m | Tertiary sediment cover | 2200 | 1200 | 2100 |
| 200-1500 m | Mesozoic / weathered granite | 5400 | 3100 | 2650 |
| > 1500 m | Pre-Cambrian basement | 6200 | 3550 | 2750 |

Composite based on Bao et al. (2011) regional crustal model and Wen
et al. (1999) source-region geology.

## Run

```bash
./run.sh
```

Output lands in `examples/40_lop_nor_1996/output/`. Four
seismometers cover the near (4500 m), far (5500 m), and epicentral
geometry.

## V&V

Cached IRIS waveforms live at `tools/waveform_vv/cache/lop_nor_1996/`,
populated by `scripts/fetch_lop_nor_1996_waveforms.py`. Pass-12
ships a smoke-level integration test in
`tests/integration/test_historic_nuclear.cpp::Lop_Nor_1996_Basin`.

## References

- Wen, L. and Wen, Y. (1999), "Seismic source parameters for the
  1996 Lop Nor underground nuclear explosion", GRL 26, 2005-2008.
- Sykes, L. R. (2002), "Yields of the Soviet & Chinese explosions",
  CTBTO Bulletin / Conference Proceedings.
- Bao, X., et al. (2011), "Crust and upper mantle structure
  beneath the Tian Shan and Tarim basin from receiver functions",
  J. Geophys. Res. 116.
- IDC bulletin event 1996.211.
- USGS event service: us1996.07.29.
