# Lop Nor 2020 Alleged Decoupled Nuclear Test

## Event Summary

| Parameter | Value |
|-----------|-------|
| **Date/Time** | 2020-06-22 ~09:18 UTC |
| **Location** | 41.735°N, 88.730°E |
| **Site** | Lop Nor Nuclear Test Site, Xinjiang, China |
| **Depth** | ~300 m (tunnel in Kuruktag mountains) |
| **Apparent Yield** | ~7–14 tonnes (decoupled equivalent) |
| **Estimated True Yield** | ~0.5–1.0 kt (if fully decoupled) |
| **Magnitude** | mb ~2.75 (NORSAR) |
| **Detection** | PS23 Makanchi array, Kazakhstan (~760 km) |

## Background

On June 22, 2020, seismic monitoring stations detected a small seismic event near the Lop Nor nuclear test site in Xinjiang, China. The event was first detected at the PS23 IMS array in Makanchi, Kazakhstan, approximately 760 km to the northwest.

### Detection Characteristics

The event exhibited characteristics consistent with an explosion rather than an earthquake:

1. **Strong P / Weak S**: High P-to-S amplitude ratio typical of explosions
2. **Impulsive First Motion**: Sharp P-wave onset
3. **Location**: Within the known Lop Nor nuclear test area
4. **Paired Events**: Two very small events separated by ~12 seconds

### Decoupling Hypothesis

The extremely small apparent magnitude (mb ~2.75) suggests one of two possibilities:

1. **Very Small Coupled Test**: A sub-ton conventional explosive or extremely small nuclear yield
2. **Decoupled Nuclear Test**: A larger yield (~1 kt) detonated inside a pre-excavated cavity

FSRM models this event assuming the decoupled hypothesis, which has significant implications for treaty verification.

## Physics of Decoupled Nuclear Tests

### What is Decoupling?

**Decoupling** is a technique to reduce the seismic signal from an underground nuclear explosion by detonating the device inside a large pre-excavated cavity, typically in hard rock like granite or salt.

#### Physical Mechanism

In a **coupled test** (normal underground explosion):
1. The nuclear device detonates in direct contact with surrounding rock
2. The expanding fireball vaporizes nearby rock, creating a cavity
3. Shock waves couple efficiently into the surrounding medium
4. Seismic waves radiate outward at full strength

In a **decoupled test** (cavity-decoupled explosion):
1. A large cavity is excavated before the test (by conventional mining or a previous nuclear shot)
2. The device detonates in the center of the air-filled cavity
3. The expanding fireball first compresses the air in the cavity
4. The shock wave is attenuated before reaching the cavity walls
5. Reduced seismic coupling results in weaker seismic signals

#### Decoupling Factor

The **decoupling factor (DF)** is the ratio of seismic amplitude from a coupled test to that from a decoupled test of the same yield:

$$DF = \frac{A_{\text{coupled}}}{A_{\text{decoupled}}}$$

For full decoupling in hard rock:
- **Granite**: DF ≈ 70
- **Salt**: DF ≈ 70
- **Tuff**: DF ≈ 30 (lower due to porosity)

#### Required Cavity Size

The cavity radius required for full decoupling scales with yield:

$$R_c \geq 25 \times W^{1/3} \text{ meters}$$

Where $W$ is yield in kilotons. For a 1 kt device:
- Required cavity radius: ≥ 25 m
- Cavity volume: ≥ 65,000 m³

### Seismic Magnitude Reduction

Decoupling reduces the body-wave magnitude (mb) by:

$$\Delta m_b = \log_{10}(DF) \approx 1.85 \text{ magnitude units}$$

For a 1 kt coupled test (mb ≈ 4.0), full decoupling would produce:
- Apparent mb ≈ 2.15–2.75
- Apparent coupled-equivalent yield ≈ 7–14 tonnes

This matches the observed mb ~2.75 for the Lop Nor 2020 event.

### Treaty Verification Implications

Decoupling poses significant challenges for the Comprehensive Nuclear-Test-Ban Treaty (CTBT):

1. **Detection Threshold**: Decoupled tests may fall below IMS detection limits
2. **Yield Estimation**: Standard mb-yield relations underestimate true yield
3. **Discrimination**: Small explosions are harder to distinguish from earthquakes
4. **Evasion Scenario**: Nations could conduct militarily significant tests covertly

FSRM simulations help treaty verification by:
- Modeling expected signals from decoupled scenarios
- Testing detection algorithm sensitivity
- Developing discrimination methods for cavity-decoupled sources

## Physics of Underground Nuclear Explosions

### Detonation Phenomenology

An underground nuclear explosion progresses through several phases:

#### Phase 1: Detonation (0–1 μs)
- Nuclear reactions release energy in <1 microsecond
- Temperature reaches 10⁷–10⁸ K
- Pressure exceeds 10⁹ Pa (1 million atmospheres)

#### Phase 2: Vaporization (1 μs–1 ms)
- Intense radiation vaporizes surrounding rock
- Vapor cavity expands supersonically
- Strong shock wave propagates outward

#### Phase 3: Melting & Crushing (1 ms–1 s)
- Shock pressure decreases with distance
- Rock transitions from vaporized → melted → crushed → fractured
- Cavity reaches maximum radius

#### Phase 4: Cavity Rebound (1–10 s)
- Cavity pressure drops below lithostatic
- Walls begin to collapse inward
- Elastic rebound may cause "cavity oscillations"

#### Phase 5: Collapse (10 s–hours)
- Rubble fills cavity from above
- Chimney may propagate toward surface
- Subsidence crater forms (if shallow enough)

### Damage Zone Structure

The explosion creates concentric damage zones:

| Zone | Extent | Physical State |
|------|--------|----------------|
| Cavity | R_c | Vaporized/melted rock |
| Melt Shell | 1.5× R_c | Puddle glass at cavity floor |
| Crushed Zone | 2–3× R_c | Pulverized rock |
| Fractured Zone | 5–10× R_c | Pervasive fracturing |
| Damaged Zone | 10–20× R_c | Microcracks, reduced strength |
| Elastic Zone | >20× R_c | Elastic wave propagation |

### Cavity Radius Scaling

The final cavity radius depends on yield, depth, and rock properties:

$$R_c = C \cdot W^{0.295} \cdot \left(\frac{\rho}{2650}\right)^{-1/3.4} \cdot \left(\frac{h + 91.5}{122}\right)^{-0.175}$$

Where:
- $C$ = 55 for hard rock, 70 for tuff
- $W$ = yield (kt)
- $\rho$ = rock density (kg/m³)
- $h$ = depth of burial (m)

For the Lop Nor 2020 event (1 kt, 300 m, granite):
- Coupled cavity radius: ~55 m
- Pre-excavated decoupling cavity: 25 m (sufficient for DF ≈ 70)

## Mueller-Murphy Source Model

### Theory

The **Mueller-Murphy model** (1971) describes the seismic source from underground nuclear explosions using a **Reduced Displacement Potential (RDP)**:

$$\psi(t) = \psi_{\infty} \cdot \left[1 - e^{-t/\tau_1}\left(1 + \frac{t}{\tau_1}\right)\right] \cdot e^{-t/\tau_2}$$

Where:
- $\psi_{\infty}$ = steady-state potential (proportional to cavity volume)
- $\tau_1$ = rise time (cavity expansion)
- $\tau_2$ = overshoot decay time (cavity oscillation damping)

### Source Spectrum

The far-field displacement spectrum is:

$$\hat{u}(\omega) = \frac{\psi_{\infty} \cdot \omega^2}{(1 + i\omega\tau_1)^2 \cdot (1 + i\omega\tau_2)}$$

This produces:
- **Low frequencies**: Flat spectrum (proportional to cavity volume)
- **Corner frequency**: $f_c = 1/(2\pi\tau_1)$
- **High frequencies**: Roll-off as $\omega^{-2}$

### Patton Scaling

For the corner frequency, FSRM uses **Patton (1988) scaling**:

$$f_c = 1.5 \cdot W^{-0.3} \text{ Hz}$$

For a 1 kt explosion: $f_c \approx 1.5$ Hz
For a 250 kt explosion: $f_c \appro 0.3$ Hz

### Decoupled Source Modification

For a decoupled explosion, the Mueller-Murphy model is modified:
1. Reduced $\psi_{\infty}$ by factor of DF
2. Shorter rise time (less rock interaction)
3. Reduced overshoot (cavity absorbs oscillations)
4. More purely isotropic radiation pattern

## Seismic Wave Propagation

### Regional Phase Overview

At regional distances (100–2000 km), seismic waves separate into distinct **phases**:

| Phase | Path | Velocity | Character |
|-------|------|----------|-----------|
| **Pn** | Refracted along Moho | 8.0–8.2 km/s | First arrival at >200 km |
| **Pg** | Direct crustal P | 5.5–6.5 km/s | Strong at <500 km |
| **Sn** | Refracted S along Moho | 4.4–4.7 km/s | Often weak or absent |
| **Lg** | Crustal-guided S | 3.2–3.6 km/s | Highest amplitude |
| **Rg** | Rayleigh fundamental | 2.5–3.0 km/s | Dispersive surface wave |

### Travel Time Curves

Phase arrival times follow:

$$t = \frac{\Delta}{v}$$

For the Lop Nor 2020 event, arrival times at PS23 (760 km):
- Pn: ~94 s (8.1 km/s)
- Pg: ~125 s (6.1 km/s)
- Sn: ~165 s (4.6 km/s)
- Lg: ~217 s (3.5 km/s)

### Amplitude Decay

Seismic amplitudes decay due to:

1. **Geometric spreading**: $A \propto 1/r$ (body waves), $A \propto 1/\sqrt{r}$ (surface waves)
2. **Intrinsic attenuation**: $A \propto e^{-\pi f t / Q}$
3. **Scattering**: Frequency-dependent loss in heterogeneous media

The quality factor $Q$ varies by region:
- Stable shields: Q > 1000
- Tectonically active: Q ~ 200–500
- Central Asia (Lop Nor path): Q ~ 400

## Moment Tensor Decomposition

### Full Moment Tensor

The seismic source can be represented by a **moment tensor** $M_{ij}$:

$$M = \begin{pmatrix} M_{xx} & M_{xy} & M_{xz} \\ M_{xy} & M_{yy} & M_{yz} \\ M_{xz} & M_{yz} & M_{zz} \end{pmatrix}$$

For an explosion, the moment tensor has three components:

### Isotropic (ISO) Component

Pure volume change:

$$M_{\text{ISO}} = \frac{1}{3}(M_{xx} + M_{yy} + M_{zz}) \cdot I$$

- Explosions: Dominant component (>80%)
- Earthquakes: Near zero

### Compensated Linear Vector Dipole (CLVD)

Axisymmetric deviatoric strain:

$$M_{\text{CLVD}} = \begin{pmatrix} -1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 2 \end{pmatrix} \cdot M_0^{\text{CLVD}}$$

- Explosions: Minor component from cavity collapse, spall
- Can indicate non-spherical cavity or stress interaction

### Double-Couple (DC) Component

Shear faulting:

$$M_{\text{DC}} = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix} \cdot M_0^{\text{DC}}$$

- Explosions: Small component from tectonic release
- Earthquakes: Dominant component

### Lop Nor 2020 Expected Pattern

For a fully decoupled test in granite:
- **ISO**: ~92% (highly contained explosion)
- **CLVD**: ~6% (minimal cavity asymmetry)
- **DC**: ~2% (negligible tectonic release)

## Recording Stations

| Station | Network | Location | Distance (km) | Azimuth |
|---------|---------|----------|---------------|---------|
| WMQ | IC | Urumqi, China | ~400 | 60° |
| MAKZ | IU | Makanchi, Kazakhstan | ~760 | 320° |
| MKAR | KZ | Makanchi Array | ~760 | 320° |
| WUS | G | Wushi, China | ~350 | 215° |
| PDGK | KZ | Podgonoye, Kazakhstan | ~900 | 340° |
| PRZ | KR | Karakol, Kyrgyzstan | ~650 | 290° |
| KNDC | KZ | Almaty, Kazakhstan | ~800 | 300° |
| AAK | II | Ala Archa, Kyrgyzstan | ~800 | 280° |
| KURK | II | Kurchatov, Kazakhstan | ~1200 | 350° |
| LSA | IC | Lhasa, Tibet | ~1700 | 170° |
| NIL | II | Nilore, Pakistan | ~1800 | 240° |
| XAN | IC | Xi'an, China | ~2100 | 90° |

## Generated Figures

All figures are located in `figures/lop_nor_2020/`:

### Station Map (`station_map.png`)

![Station Map](../../figures/lop_nor_2020/station_map.png)

PyGMT map showing:
- Event epicenter at Lop Nor test site
- IMS and regional seismic stations
- Topography of Central Asia (Tien Shan, Tarim Basin)
- Great-circle paths from event to stations

### Distance-Scaled Waveform Gather (`observed_waveform_gather_distance_scaled.png`)

![Distance-Scaled Waveform Gather](../../figures/lop_nor_2020/observed_waveform_gather_distance_scaled.png)

Publication-quality record section with:
- Y-axis scaled to true epicentral distance
- Theoretical phase arrival lines (Pn, Pg, Sn, Lg)
- Waveform amplitudes normalized for comparison
- Clear identification of regional phases

### Full Waveform Gather (`observed_waveform_gather.png`)

Standard record section with full time window showing all recorded phases.

### Zoomed Waveform Gather (`observed_waveform_gather_zoomed.png`)

Focused view on the expected signal window for phase identification.

### Long-Period Waveform Gather (`observed_waveform_gather_long_period.png`)

Record section with 0.02–0.10 Hz bandpass for moment tensor analysis.

### Spectrograms (`observed_spectrograms.png`)

Time-frequency analysis showing frequency content evolution.

### PS23/WMQ Close-up (`observed_closeup_ps23.png`)

Detailed view of the primary detection stations.

### Moment Tensor Results (`moment_tensor_inversion.png`)

Full moment tensor inversion showing ISO/CLVD/DC decomposition.

## FSRM Configuration

### Configuration File

```
config/lop_nor_2020_decoupled_granite.config
```

### Key Parameters

```ini
[SIMULATION]
name = lop_nor_2020_decoupled_granite
description = Decoupled ~1 kt underground nuclear test in Lop Nor granite

[EXPLOSION_SOURCE]
type = NUCLEAR_UNDERGROUND
yield_kt = 1.0
depth_of_burial = 300.0               # m

# Decoupling cavity
cavity_radius = 25.0                  # Pre-excavated cavity
cavity_radius_scaling = OVERRIDE

# Reduced damage zones (cavity absorbs most energy)
crushed_zone_radius = 50.0
fractured_zone_radius = 120.0

# Highly isotropic source
isotropic_fraction = 0.92
clvd_fraction = 0.06
double_couple_fraction = 0.02

# Reduced scalar moment (DF ≈ 70)
scalar_moment = 1.4e14                # N·m
```

### Running the Simulation

```bash
# Run FSRM simulation
mpirun -np 8 ./build/fsrm -c config/lop_nor_2020_decoupled_granite.config

# Quick test version
mpirun -np 4 ./build/fsrm -c config/lop_nor_2020_decoupled_granite_quick.config
```

## Data Processing

### Waveform Download

```bash
python3 scripts/fetch_lop_nor_2020_waveforms.py
```

### Station Map

```bash
python3 scripts/plot_lop_nor_2020_station_map.py
```

### Visualization

```bash
python3 scripts/visualize_fsrm_lop_nor.py
```

### Moment Tensor Inversion

```bash
python3 scripts/invert_lop_nor_2020_moment_tensor.py
```

## Scientific Analysis

### Detection Assessment

The Lop Nor 2020 event was marginally detected:
- Primary detection: PS23 Makanchi (760 km)
- Confirmed by: WMQ Urumqi (400 km)
- Near detection threshold at more distant stations

This validates the decoupling hypothesis—a coupled 1 kt test would produce mb ~4.0 and be easily detected globally.

### Discrimination Evidence

Evidence supporting an explosion source:
1. **P/S Ratio**: High P-wave to S-wave amplitude ratio
2. **Ms:mb**: Low surface-wave to body-wave magnitude ratio
3. **Spectral Ratio**: Elevated high-frequency content
4. **Location**: Within known nuclear test area
5. **Depth**: Shallow, consistent with tunnel emplacement

### Yield Estimation

Multiple methods constrain the yield:

| Method | Apparent Yield | If Decoupled (DF=70) |
|--------|----------------|----------------------|
| mb-yield (Murphy) | ~10 tonnes | 0.7 kt |
| Corner frequency | ~5 tonnes | 0.35 kt |
| Lg amplitude | ~15 tonnes | 1.0 kt |
| Moment tensor | ~8 tonnes | 0.55 kt |

**Best estimate**: 0.5–1.0 kt true yield (if fully decoupled)

## References

1. **Mueller & Murphy (1971)** — Seismic characteristics of underground nuclear detonations. Part I: Seismic source. *BSSA* 61(6), 1675-1692.

2. **Patton (1988)** — Source models of the Harzer explosion. *BSSA* 78(5), 1133-1157.

3. **Glenn & Myers (1997)** — Decoupling of underground nuclear explosions. *JGR* 102(B5), 10065-10082.

4. **Charlie & Veyera (1994)** — Lop Nor site geology and testing history. *Defense Nuclear Agency Report*.

5. **Xia et al. (2017)** — NCCrust: A reference crustal model for North China. *BSSA* 107(6), 2752-2768.

6. **Stevens & Day (1985)** — The physical basis of mb:Ms and variable frequency magnitude methods for earthquake/explosion discrimination. *JGR* 90(B4), 3009-3020.

7. **NORSAR (2020)** — Seismic event detection report, June 2020.
