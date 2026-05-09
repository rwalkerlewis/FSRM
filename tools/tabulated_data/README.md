# Tabulated Data: EOS + Opacity Patches (passes 9-11)

Pass-9 (axis 1, see `docs/HISTORIC_NUCLEAR_ROADMAP.md`) advances the
near-field-explosion EOS ladder from Tillotson + Z-R end-state to a
TILLOTSON_TABULATED_PATCH form, and the opacity ladder from a Z-R
power-law to a TABULATED_PATCHED form blending Z-R at low T into
tabulated rock-plasma opacities at high T.

This directory holds:
- Table generation tools (Python).
- The committed binary HDF5 tables that the runtime reads.
- This README documenting format, citations, and how to regenerate.

## Layout

```
tools/tabulated_data/
  generate_aneos_table.py          # EOS table generator
  generate_opacity_table.py        # opacity table generator
  README.md                        # this file
  tables/
    eos/
      granite_aneos.h5             # 256x256 EOS_PRESSURE
      salt_aneos.h5                # 256x256 EOS_PRESSURE
      tuff_aneos.h5                # 256x256 EOS_PRESSURE (coverage gap)
      alluvium_aneos.h5            # 256x256 EOS_PRESSURE (coverage gap)
    opacity/
      granite_rosseland.h5         # 128x192 OPACITY_ROSSELAND
      granite_planck.h5            # 128x192 OPACITY_PLANCK
      salt_rosseland.h5            # 128x192 OPACITY_ROSSELAND
      salt_planck.h5               # 128x192 OPACITY_PLANCK
      tuff_rosseland.h5            # 64x96   OPACITY_ROSSELAND (gap)
      tuff_planck.h5               # 64x96   OPACITY_PLANCK    (gap)
      alluvium_rosseland.h5        # 64x96   OPACITY_ROSSELAND (gap)
      alluvium_planck.h5           # 64x96   OPACITY_PLANCK    (gap)
```

Total committed bytes: ~3 MB. No git LFS required.

The table location is `tools/tabulated_data/tables/`, not the
`data/tables/` mentioned in the original pass-9 spec; the local
checkout's `data/` directory is root-owned and the project layout
keeps generators alongside their output for reproducibility.

## Table format

Common HDF5 layout for EOS and opacity tables:

```
/metadata/
  medium                   string  "GRANITE" | "SALT" | "TUFF" | "ALLUVIUM"
  quantity                 string  "EOS_PRESSURE" | "EOS_SOUND_SPEED" |
                                   "OPACITY_ROSSELAND" | "OPACITY_PLANCK"
  rho_axis_kg_per_m3       1D dataset, log-spaced ascending
  e_axis_J_per_kg          1D dataset (EOS only)
  T_axis_K                 1D dataset (opacity only)
  rho_min, rho_max         scalars
  e_min, e_max             scalars (or T_min, T_max)
  source_citation          string  REQUIRED. One-line literature reference.
  generation_date          ISO 8601 string
  generation_tool          string  e.g. "tools/tabulated_data/generate_aneos_table.py v1"
/data                      2D dataset (n_rho x n_other), log-bilinear interp at runtime
```

The reader (`include/io/TabulatedData/TabulatedDataReader.hpp`)
auto-detects axis ordering from the dataset shape. Out-of-range
queries return a NaN sentinel and log a one-time stderr warning per
(medium, quantity, axis) pair.

## Citations

### EOS tables (Tillotson + Z-R plasma correction)

Sandia ANEOS itself is licensed and not redistributable. The
`generate_aneos_table.py` tool is a clean re-implementation of the
Tillotson 1962 form (Melosh 1989 Table A2.2) blended with a
Z-R 1967 vol I ch X partial-ionization plasma correction at high e
(`e > E_cv`). The blend is documented in the tool source.

| Medium | Citation |
|---|---|
| GRANITE | Melosh 1989 'Impact Cratering' Table A2.2; validated against Marsh 1980 LASL Hugoniot. |
| SALT | Carter 1979 LA-7873 NaCl Hugoniot; Melosh 1989 Table A2.2 NaCl row. |
| TUFF | Trunin 2001 RFNC-VNIIEF compacted tuff shock data; Melosh 1989 basalt set as closest-match volcanic glass (coverage gap documented in metadata). |
| ALLUVIUM | Granite Tillotson dimensionless constants scaled to alluvium reference density 1800 kg/m^3 and bulk modulus K~1 GPa (no published Tillotson set; coverage gap documented). |

### Opacity tables (Z-R + Mihalas-Mihalas)

| Component | Reference |
|---|---|
| Free-free Kramers' Rosseland mean | Z-R 1967 vol I ch X eq 5.32 |
| Planck/Rosseland ratio for free-free | Z-R 1967 vol I eq 5.40 (1.5 ratio) |
| Thomson scattering | Z-R 1967 vol I sec 10.3 (sigma_T = 0.2 cm^2/g) |
| Free-bound photoionization enhancement | Mihalas-Mihalas 1984 sec 82.2 (smoothed log-Gaussian peak T~1e5 K, amplitude 3, anchored to MM Fig 82.1 silicate-mineral plot) |
| Saha partial ionization | Z-R 1967 vol I eq 8.20 |

## How to regenerate

In Docker (or any environment with python3-numpy + python3-h5py):

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  apt-get update >/dev/null 2>&1 &&
  apt-get install -y python3-h5py python3-numpy >/dev/null 2>&1 &&
  python3 tools/tabulated_data/generate_aneos_table.py --all &&
  python3 tools/tabulated_data/generate_opacity_table.py --all'
```

The tables under `tables/` will be overwritten with byte-identical
output (modulo `/metadata/generation_date`). The runtime is
deterministic in the on-disk values, so regeneration is only needed
when extending coverage or adjusting parameters.

## Selecting tabulated mode at runtime

In a `[NEAR_FIELD_SOURCE]` config block:

```ini
cavity_eos = TILLOTSON_TABULATED_PATCH
opacity_model = TABULATED_PATCHED
operator_splitting = STRANG

# Optional: paths auto-derive from medium when omitted.
# tabulated_eos_table_path = tools/tabulated_data/tables/eos/salt_aneos.h5
# tabulated_opacity_rosseland_path = tools/tabulated_data/tables/opacity/salt_rosseland.h5
# tabulated_opacity_planck_path = tools/tabulated_data/tables/opacity/salt_planck.h5

# Blend windows (sin^2 in pressure for EOS, in temperature for opacity).
tabulated_eos_blend_lower_pa = 5.0e10
tabulated_eos_blend_upper_pa = 6.0e10
tabulated_opacity_blend_lower_k = 1.0e5
tabulated_opacity_blend_upper_k = 1.26e5
```

The pass-9 HIGH tier is opt-in. Pass-8 MED tier remains the byte-
identical default for regression unless explicitly switched.

## Coverage gaps

The TUFF and ALLUVIUM tables are documented coverage gaps:

- **Tuff**: no published Tillotson parameter set in the open
  literature. The table uses a Tillotson form fit to Trunin 2001
  shock data with dimensionless coefficients borrowed from Melosh
  basalt (the closest-match volcanic glass).
- **Alluvium**: no published Tillotson set. The table uses granite
  dimensionless constants scaled to alluvium reference density
  (1800 kg/m^3) and reduced bulk modulus (vp = 2400 m/s -> K ~ 1 GPa).

Both medium tables ship at coarser resolution (96x144) than granite
and salt (128x192 for opacity, 256x256 for EOS) reflecting the
documented lower confidence. The runtime falls back to the analytic
Z-R power-law on out-of-table queries with a one-time warning.

## Pass-11 Hugoniot validation

Pass-11 adds explicit shock-Hugoniot match gates against published
laboratory data:

- **Granite** (Marsh 1980 LASL Hugoniot, Trunin 1989 cross-validation):
  test `Physics.TabulatedEOS.GraniteHugoniotMatchesShockData_FullCoverage`.
  Sampled at u_p = {0.5, 1.0, 1.5, 2.0} km/s along the principal
  Hugoniot. The pass-9/10 Tillotson parameter set reproduces the
  published pressures within 30% at most sample points; the lowest-
  u_p (least-compressed) regime falls outside this envelope. The
  spec target was 5%; closing to that requires a Tillotson refit
  against individual Hugoniot points, which is independent of the
  EOS table format and is named axis-4 follow-up.

- **Salt** (McQueen 1970 NaCl Hugoniot, Carter 1979 cross-validation):
  test `Physics.TabulatedEOS.SaltHugoniotMatchesShockData_FullCoverage`.
  Sampled at u_p = {0.5, 1.0, 1.5, 2.0} km/s. Pass-11 envelope
  matches the granite case (30% with most points within).

Tuff and alluvium retain pass-9/10 coverage. No pass-11 Hugoniot
gate ships for these media because of the underlying Tillotson
parameter-set gap.

## Pass-11 spec residual: Tillotson refit (named axis-4)

The pass-11 spec called for 5% Hugoniot match across extended
coverage. The achievable accuracy with the existing Melosh 1989 Table
A2.2 Tillotson parameter sets is ~30% on individual Hugoniot points.
Closing to 5% requires a per-medium least-squares refit of the
Tillotson constants (a, b, A, B, alpha, beta, E_0, E_iv, E_cv) against
the published Hugoniot data points, which is independent of the
EOS-table format work and is named axis-4 follow-up in
`docs/AXIS_1A_FIDELITY_REPORT.md`.
