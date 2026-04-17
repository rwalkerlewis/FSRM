# US Divider 1992 Nuclear Test

## Event Summary

| Parameter | Value |
|-----------|-------|
| **Date/Time** | 1992-09-23 15:04:00 UTC |
| **Location** | 37.021°N, 116.058°W |
| **Site** | Nevada Test Site, Area 4 |
| **Depth** | ~800 m below surface |
| **Yield** | <20 kt (announced); estimated ~1 kt |
| **Magnitude** | mb ~4.0 |
| **Character** | Last US underground nuclear test before moratorium |

## Background

Operation Divider (also designated "Julin/Divider") was the final US underground nuclear test, conducted on September 23, 1992, just days before President George H.W. Bush signed the nuclear testing moratorium into law. It marked the end of nearly 50 years of US nuclear testing.

### Historical Significance

- **Last US Test**: Final underground nuclear detonation in US history
- **Cold War End**: Conducted during the post-Soviet transition period
- **CTBT Precursor**: Led directly to Comprehensive Test Ban Treaty negotiations
- **Verification Research**: Data used for treaty verification method development

### Operation Julin

Divider was part of Operation Julin, the last US nuclear test series (1991–1992). The operation included:
- Junction (Sept 26, 1992) — simultaneous with Divider
- Hunters Trophy (Sept 18, 1992)
- Galena/Yellow (June 25, 1992)
- Diamond Fortune (April 30, 1992)

## Geology

### Nevada Test Site Setting

The Nevada Test Site (now Nevada National Security Site) is located in the Basin and Range province of southern Nevada. The regional geology includes:

- **Alluvium**: Quaternary basin fill (0–600 m)
- **Volcanic Tuff**: Miocene ash-flow tuffs from the Southwest Nevada Volcanic Field
- **Paleozoic Carbonates**: Limestones and dolomites
- **Precambrian Basement**: Crystalline metamorphic rocks

### Test Emplacement Medium

Area 4 tests were typically conducted in:
- **Rainier Mesa**: Zeolitized tuff with moderate water saturation
- **Depth**: 600–1000 m below surface
- **Rock Properties**:
  - Density: 1800–2200 kg/m³ (tuff)
  - P-wave Velocity: 3.0–4.5 km/s (tuff), 5.0–6.0 km/s (carbonate)
  - Porosity: 10–40% (high for tuff)

### Crustal Structure

The Nevada Test Site has well-characterized velocity structure from decades of nuclear test monitoring:

| Layer | Depth Range | Vp (km/s) | Vs (km/s) |
|-------|-------------|-----------|-----------|
| Alluvium/Tuff | 0–1 km | 2.5–3.5 | 1.4–2.0 |
| Volcanic | 1–5 km | 4.0–5.0 | 2.3–2.9 |
| Paleozoic | 5–15 km | 5.5–6.2 | 3.2–3.6 |
| Lower Crust | 15–30 km | 6.5–7.0 | 3.7–4.0 |
| Mantle | >30 km | 7.8–8.0 | 4.4–4.6 |

## Source Physics

### Yield Estimation

The announced yield was "<20 kt," with most seismological estimates placing it at ~1 kt based on:
- Body-wave magnitude (mb ~4.0)
- Comparison with other NTS tests of known yield
- Corner frequency analysis

### Cavity Radius

For a ~1 kt test in porous tuff at 800 m depth:
$$R_c = 55 \times W^{0.295} \times \left(\frac{\rho}{2650}\right)^{-1/3.4}$$

- **Estimated Cavity Radius**: ~20–30 m
- **Crushed Zone**: ~60–90 m
- **Fractured Zone**: ~150–200 m

### Tuff vs. Granite

Tests in tuff versus granite exhibit important differences:

| Property | Tuff | Granite |
|----------|------|---------|
| Porosity | 10–40% | <5% |
| Water Content | Variable | Low |
| Cavity Radius | Larger | Smaller |
| Corner Frequency | Lower | Higher |
| P/S Ratio | Lower | Higher |
| High-Freq Content | Less | More |

This is explored in the figure `fig05_tuff_vs_granite.png`.

## Physics of Underground Nuclear Tests in Porous Media

### Volcanic Tuff as a Test Medium

The Nevada Test Site used volcanic tuff for many tests because of its unique properties:

#### Material Properties

| Property | Zeolitized Tuff | Granite (for comparison) |
|----------|-----------------|-------------------------|
| Density | 1800–2200 kg/m³ | 2650–2750 kg/m³ |
| Porosity | 15–40% | <3% |
| P-wave velocity | 2.5–4.0 km/s | 5.5–6.5 km/s |
| S-wave velocity | 1.4–2.3 km/s | 3.2–3.8 km/s |
| Water saturation | 20–80% | <10% |
| Compressive strength | 10–50 MPa | 100–250 MPa |

#### Effect on Explosion Phenomenology

Porous tuff significantly modifies explosion behavior:

1. **Larger Cavity Radius**: Lower rock strength allows greater expansion
   $$R_c^{\text{tuff}} \approx 1.3 \times R_c^{\text{granite}}$$

2. **Lower Corner Frequency**: Larger cavity means longer source duration
   $$f_c^{\text{tuff}} < f_c^{\text{granite}}$$

3. **Reduced High-Frequency Radiation**: Porous medium attenuates high frequencies

4. **Different mb-Yield Scaling**: Tuff tests appear smaller for same yield
   $$\Delta m_b \approx -0.3 \text{ to } -0.5$$ relative to hard rock

### Detonation Process in Tuff

#### Phase 1: Initial Shock (0–1 ms)

The nuclear detonation creates a shock wave that:
- Collapses pore space ahead of the shock front
- Compresses water in pores (if saturated)
- Lower initial peak pressure than hard rock (energy absorbed by pore collapse)

#### Phase 2: Cavity Formation (1 ms–100 ms)

Cavity growth differs from hard rock:

$$\frac{dR}{dt} = \sqrt{\frac{2(P_{cavity} - P_{litho})}{\rho}}$$

In tuff:
- Lower lithostatic pressure (lower density)
- Lower rock strength
- Result: larger cavity, longer expansion time

#### Phase 3: Containment and Chimney

Tuff's properties affect containment:
- Higher porosity absorbs shock energy → better containment
- Weaker roof → earlier collapse
- Higher probability of chimney closure before reaching surface

### The Mueller-Murphy Source Model

#### Reduced Displacement Potential (RDP)

The Mueller-Murphy (1971) model describes the seismic source as:

$$\psi(t) = \psi_\infty \cdot \left[1 - e^{-t/\tau_1}\left(1 + \frac{t}{\tau_1}\right)\right] \cdot e^{-t/\tau_2}$$

**Physical parameters**:

| Parameter | Physical Meaning | Typical Value (1 kt in tuff) |
|-----------|-----------------|-----------------------------|
| $\psi_\infty$ | Steady-state displacement potential | ~10⁵ m³ |
| $\tau_1$ | Rise time (cavity expansion) | 15–25 ms |
| $\tau_2$ | Overshoot decay time | 50–100 ms |

#### Patton Corner Frequency Scaling

FSRM uses the Patton (1988) relation:

$$f_c = 1.5 \cdot W^{-0.3} \cdot \left(\frac{\rho}{2650}\right)^{0.4}$$

For 1 kt in tuff (ρ = 2000 kg/m³):
- $f_c$ ≈ 1.2 Hz (vs ~1.5 Hz for granite)

#### Source Spectrum

The far-field displacement spectrum:

$$|\hat{u}(f)| = \frac{\Omega_0 \cdot (f/f_c)^2}{1 + (f/f_c)^4}$$

Key features:
- **Long-period level**: $\Omega_0 \propto$ cavity volume
- **Corner frequency**: $f_c \approx 1.2$ Hz
- **High-frequency decay**: $f^{-2}$ above corner

### Seismic Wave Propagation in Basin and Range

The Nevada Test Site lies in the Basin and Range extensional province, which profoundly affects wave propagation.

#### Crustal Structure

```
    0 ────────────────── Surface
    │   Alluvium / Basin fill
    │   Vp = 2.0-3.5 km/s
  5 ├── ─ ─ ─ ─ ─ ─ ─ ─  
    │   Volcanic Section
    │   Vp = 4.0-5.0 km/s
 15 ├── ─ ─ ─ ─ ─ ─ ─ ─
    │   Paleozoic Carbonates
    │   Vp = 5.5-6.2 km/s
 30 ├── ═══════════════ Moho
    │   Upper Mantle
    │   Vp = 7.8-8.0 km/s
```

#### Regional Phase Propagation

Phase velocities in the Basin and Range:

| Phase | Velocity (km/s) | Path |
|-------|-----------------|------|
| **Pg** | 5.5–6.0 | Direct crustal P |
| **Pn** | 7.8–8.0 | Moho head wave |
| **Lg** | 3.0–3.4 | Crustal S waveguide |
| **Rg** | 2.5–2.8 | Rayleigh fundamental |

#### Basin Effects

Intervening basins cause:
1. **Lg Attenuation**: Extensional basins strongly attenuate Lg
2. **Amplification**: Basin sediments amplify ground motion
3. **Multipathing**: Complex arrivals from basin edges
4. **Basin-Edge Diffraction**: Secondary arrivals from basin margins

### Yield Estimation

#### Body-Wave Magnitude Method

The mb-yield relation for NTS tests in tuff:

$$m_b = 4.05 + 0.75 \log_{10}(W)$$

(Note: 0.4 units lower than hard-rock tests)

For Divider with $m_b$ ≈ 4.0:
- Inferred yield: ~1 kt
- Consistent with "<20 kt" announcement

#### Path Corrections

Yield estimation requires corrections for:
1. **Emplacement medium**: Tuff vs. granite
2. **Depth**: Deeper tests couple more efficiently
3. **Water table**: Saturated rock couples differently
4. **Station corrections**: Site amplification at each station

### Explosion-Earthquake Discrimination

Key discriminants for identifying Divider as an explosion:

1. **Ms:mb Ratio**
   - Explosions: $M_s - m_b < 0$
   - Earthquakes: $M_s - m_b > 0$
   - Divider: Characteristic explosion signature

2. **P/S Spectral Ratio**
   - Explosions: High P/S (compressional source)
   - Earthquakes: Low P/S (shear source)

3. **First-Motion Polarity**
   - Explosions: Compressional everywhere
   - Earthquakes: Quadrantal pattern

4. **Short-Period/Long-Period Ratio**
   - Explosions: Enhanced short-period energy
   - Earthquakes: More long-period energy

## Recording Stations

The following CI TERRAscope and regional stations recorded the event:

| Station | Network | Location | Distance (km) | Azimuth |
|---------|---------|----------|---------------|---------|
| GSC | CI | Goldstone, CA | 202 | 200° |
| PAS | CI | Pasadena, CA | 372 | 212° |
| PFO | CI | Piñon Flat, CA | 380 | 186° |
| SBC | CI | Santa Barbara, CA | 437 | 230° |

**Note**: 1992 broadband data availability is limited. Many stations (ISA, SVD, SAO, ELK, MNV, ANMO) were either not yet deployed, had only short-period instruments, or had data gaps for this event.

## Generated Figures

All figures are located in `figures/us_divider_1992/`:

### Station Map (`station_map.png`)

![Station Map](../../figures/us_divider_1992/station_map.png)

PyGMT map showing:
- Event epicenter (red star) at NTS Area 4
- CI TERRAscope stations (inverted triangles)
- Regional station coverage in California and Nevada
- Topography of the Basin and Range province

### Waveform Gather (`observed_waveform_gather.png`)

Standard record section showing:
- All retrieved station vertical-component (BHZ) waveforms
- Sorted by epicentral distance (202–437 km)
- Bandpass filtered 0.5–8.0 Hz
- Instrument response removed (velocity output)
- Time window: 120 s before to 600 s after origin

### Distance-Scaled Waveform Gather (`observed_waveform_gather_distance_scaled.png`)

![Distance-Scaled Waveform Gather](../../figures/us_divider_1992/observed_waveform_gather_distance_scaled.png)

Publication-quality record section with:
- Y-axis scaled to true epicentral distance
- Theoretical phase arrival lines:
  - **Pn** (blue solid): 8.1 km/s mantle P-wave
  - **Pg** (green dashed): 6.1 km/s crustal P-wave
  - **Sn** (red solid): 4.6 km/s mantle S-wave
  - **Lg** (orange dashed): 3.5 km/s crustal S-wave
- Clear identification of crustal phases at regional distances
- Waveform amplitudes normalized for visual comparison

### Zoomed Waveform Gather (`observed_waveform_gather_zoomed.png`)

Focused view highlighting:
- First-arriving P-waves
- P-to-S conversion timing
- Relative amplitude variations with distance
- Signal quality at each station

### Long-Period Waveform Gather (`observed_waveform_gather_long_period.png`)

Record section with long-period filtering:
- Bandpass: 0.02–0.10 Hz (10–50 s period)
- Surface wave analysis band
- Moment tensor inversion preparation
- Note: Small events like Divider have weak long-period energy

### Spectrograms (`observed_spectrograms.png`)

Time-frequency analysis showing:
- Frequency content at each station
- Arrival time variation with frequency
- Dispersion characteristics
- Useful for phase identification

### Key Station Close-up (`observed_closeup_key_stations.png`)

Detailed view of highest-quality stations:
- GSC and PAS waveforms
- Phase windows annotated
- Signal-to-noise assessment

### Moment Tensor Results (`moment_tensor_inversion.png`, `moment_tensor_summary.png`)

Full moment tensor inversion:
- Predominantly isotropic (explosion) source
- Small non-isotropic components from:
  - Spall effects
  - Local stress release
  - Path effects

### Tuff vs. Granite Comparison (`fig05_tuff_vs_granite.png`)

Comparison plot showing:
- Synthetic seismograms for same yield in different media
- Amplitude and frequency content differences
- Implications for yield estimation
- Discrimination characteristics

## FSRM Configuration

The FSRM simulation uses the configuration file:

```
config/us_divider_1992_far_field_seismograms.config
```

### Key Parameters

```ini
[SIMULATION]
name = us_divider_1992_far_field_seismograms
end_time = 40.0                      # 40 s — captures all phases to 75 km

# COUPLED_ANALYTIC mode — Mueller-Murphy RDP as analytical source
explosion_solve_mode = COUPLED_ANALYTIC

# Pure linear elastodynamics — far-field propagation only
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = false
enable_elastodynamics = true

# Numerical parameters for regional wave propagation
use_discontinuous_galerkin = true
dg_order = 3
use_ader_time_integration = true

[EXPLOSION_SOURCE]
type = NUCLEAR_UNDERGROUND
yield_kt = 1.0                       # Estimated actual yield
depth_of_burial = 800.0              # m
emplacement_medium = TUFF            # Volcanic tuff
porosity = 0.25                      # 25% porosity
water_saturation = 0.6               # 60% water saturated
```

### Near-Field vs. Far-Field Configurations

FSRM provides two configuration approaches for Divider:

1. **Near-Field** (`us_divider_1992_complex_geology_gmsh.config`):
   - Detailed Gmsh unstructured mesh
   - Complex NTS geology (tuff units, faults, intrusions)
   - Nonlinear material response
   - Critical for first few km of propagation

2. **Far-Field** (`us_divider_1992_far_field_seismograms.config`):
   - Structured grid with 500 m spacing
   - 1-D layered velocity model
   - Linear elastic propagation
   - Efficient for regional distances (10–75 km)

### Running the Simulation

```bash
# Far-field wave propagation simulation
mpirun -np 8 ./build/fsrm -c config/us_divider_1992_far_field_seismograms.config

# Near-field simulation with complex geology (requires Gmsh mesh)
python3 scripts/generate_divider_binary_mesh.py
mpirun -np 16 ./build/fsrm -c config/us_divider_1992_complex_geology_gmsh.config
```

## Data Processing Pipeline

### 1. Waveform Download

```bash
python3 scripts/fetch_us_divider_1992_waveforms.py
```

This script:
- Connects to IRIS FDSN web services
- Attempts multiple channel codes (BHZ, LHZ, VHZ, SHZ) due to 1992 data limitations
- Downloads available station data
- Saves processed waveforms to `figures/us_divider_1992/`

**Note**: 1992 data availability is limited. The script includes robust error handling for missing data.

### 2. Station Map Generation

```bash
python3 scripts/plot_us_divider_1992_station_map.py
```

Uses PyGMT to create a map of the western US showing:
- NTS location and surrounding stations
- Basin and Range topography
- State boundaries

### 3. Waveform Processing

Processing pipeline:
1. **Detrend**: Remove mean and linear trend
2. **Taper**: Apply 5% cosine taper
3. **Response Removal**: Deconvolve to velocity (water level = 60 dB)
4. **Bandpass Filter**: 0.5–8.0 Hz (4-pole Butterworth, zero-phase)

### 4. Visualization and Analysis

```bash
python3 scripts/visualize_fsrm_us_divider.py
```

Generates comparison plots between FSRM synthetics and observed data.

## Scientific Interpretation

### Phase Identification

At regional distances (200–450 km), the dominant phases are:

1. **Pg Phase**: Direct crustal P-wave (~6.1 km/s)
   - First arrival at distances <300 km
   - Strong amplitude in tuff-to-carbonate path

2. **Pn Phase**: Mantle-refracted P-wave (~8.1 km/s)
   - Becomes first arrival at distances >200 km
   - Head wave along Moho

3. **Lg Phase**: Crustal-guided S-wave (~3.5 km/s)
   - Highest amplitude at regional distances
   - Complex waveform due to multiple crustal reflections

4. **Rg Phase**: Rayleigh surface wave (~2.5 km/s)
   - Visible at shorter distances
   - Dispersive character

### Characteristic Path Effects

The Basin and Range province introduces:
- Strong lateral velocity variations
- Basin-induced amplification
- High Lg attenuation in extensional basins
- Complex multipathing

### Tuff Emplacement Effects

Tests in porous tuff exhibit:
- Lower P/S amplitude ratios than granite tests
- Reduced high-frequency content
- Modified cavity collapse dynamics
- Different corner frequency scaling

## References

1. **Mueller & Murphy (1971)** — Seismic characteristics of underground nuclear detonations. Part I: Seismic source. *BSSA* 61(6), 1675-1692.

2. **Patton (1988)** — Source models of the Harzer explosion from regional observations. *BSSA* 78(5), 1133-1157.

3. **Werth & Herbst (1963)** — NTS seismic velocity structure. *JGR* 68(6), 1875-1893.

4. **Priestley et al. (1988)** — Crustal structure beneath the Nevada Test Site. *JGR* 93(B3), 2247-2259.

5. **Springer (1966)** — NTS near-surface velocity profiles. *BSSA* 56(5), 1051-1073.

6. **Murphy (1981)** — P-wave coupling of underground explosions in various geologic media. *BSSA* 71(6), 1789-1806.

7. **US DOE (1994)** — United States Nuclear Tests, July 1945 through September 1992. DOE/NV-209 REV 15.

## Historical Notes

### The Testing Moratorium

On October 2, 1992, just nine days after Divider, President George H.W. Bush signed the Energy and Water Development Appropriations Act, which included a provision establishing a nine-month moratorium on US nuclear testing. This moratorium has been extended by every subsequent administration.

### Stockpile Stewardship

Following the moratorium, the US transitioned to the Stockpile Stewardship Program, which maintains nuclear weapons reliability through:
- Subcritical experiments
- Advanced computer simulations
- Surveillance and component testing
- Archived test data analysis

FSRM supports stockpile stewardship research by enabling high-fidelity simulations of nuclear test phenomenology without actual testing.
