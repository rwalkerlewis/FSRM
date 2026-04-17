# DPRK 2017 Punggye-ri Nuclear Test

## Event Summary

| Parameter | Value |
|-----------|-------|
| **Date/Time** | 2017-09-03 03:30:01 UTC |
| **Location** | 41.300°N, 129.076°E |
| **Site** | Punggye-ri Nuclear Test Site, Mt. Mantap, DPRK |
| **Depth** | ~760 m below surface |
| **Yield** | ~250 kt (fully coupled) |
| **Magnitude** | mb 6.1–6.3, Mw 5.8–6.1 |
| **Character** | Largest DPRK nuclear test; 6th and final test at Punggye-ri |

## Background

The September 3, 2017 nuclear test was North Korea's sixth and largest underground nuclear detonation. It occurred at the Punggye-ri test site under Mt. Mantap, a ~2200 m peak in the mountainous eastern region of the country. The test was conducted in a horizontal tunnel system at approximately 760 m depth.

This test is notable for several reasons:

1. **Largest Yield**: At ~250 kt, it was roughly 10× larger than the previous DPRK tests
2. **Thermonuclear Device**: North Korea claimed it was a hydrogen bomb
3. **Mountain Collapse**: InSAR observations revealed ~0.5 m of surface subsidence over the test cavity
4. **Global Detection**: The event was recorded at seismic stations worldwide, with clear signals at distances >4000 km

## Geology

### Host Rock

The Punggye-ri test site is located in Mesozoic granite and syenite of the Korean Peninsula:

- **Rock Type**: Competent crystalline rock (granite/syenite)
- **Density**: ~2700 kg/m³
- **P-wave Velocity**: 5.6–7.0 km/s (crustal average)
- **Porosity**: Low (<1%)

### Site Characteristics

- **Mt. Mantap**: The test was conducted under a prominent mountain peak
- **Tunnel System**: Horizontal tunnel access at ~760 m depth
- **Post-Test Effects**: Satellite imagery showed surface cracks and subsidence

## Source Physics

### Cavity Formation

For a fully coupled nuclear test in hard rock, the cavity radius can be estimated using the empirical formula:

$$R_c = 55 \times W^{0.295} \times \left(\frac{\rho}{2650}\right)^{-1/3.4}$$

For W = 250 kt and ρ = 2700 kg/m³:
- **Cavity Radius**: ~280 m
- **Crushed Zone**: ~700 m (2.5× cavity)
- **Fractured Zone**: ~1400 m (5× cavity)

### Moment Tensor

The source mechanism includes:

1. **Isotropic Component**: Dominant explosion source from cavity expansion
2. **CLVD Component**: Compensated linear vector dipole from cavity collapse and asymmetry
3. **DC Component**: Small double-couple from tectonic stress release

The collapse of Mt. Mantap's summit region contributed a significant non-isotropic component to the seismic radiation.

## Physics of Coupled Underground Nuclear Explosions

### Detonation Phenomenology

A fully coupled underground nuclear test (like the DPRK 2017 event) progresses through distinct physical phases:

#### Phase 1: Nuclear Detonation (0–1 μs)

The nuclear device undergoes fission/fusion reactions in less than a microsecond:
- **Temperature**: 10⁷–10⁸ K (comparable to the Sun's core)
- **Pressure**: >10⁹ Pa (millions of atmospheres)
- **Energy release**: For 250 kt, approximately 10¹⁵ J (1 petajoule)
- **Radiation**: Intense X-rays and neutrons deposit energy in surrounding rock

#### Phase 2: Shock Wave and Vaporization (1 μs–10 ms)

The energy release creates an expanding shock wave:

$$P_{peak} = \rho_0 \cdot U_s \cdot u_p$$

Where:
- $P_{peak}$ = peak shock pressure
- $\rho_0$ = initial rock density
- $U_s$ = shock velocity
- $u_p$ = particle velocity

For granite at cavity wall: $P_{peak}$ ≈ 10–100 GPa

The shock wave:
1. Vaporizes rock within ~1 cavity radius
2. Melts rock within ~1.5 cavity radii
3. Crushes rock within ~3 cavity radii
4. Fractures rock within ~10 cavity radii

#### Phase 3: Cavity Expansion (10 ms–1 s)

The expanding gas cavity grows until internal pressure equals lithostatic stress:

$$R_c = 55 \cdot W^{0.295} \cdot \left(\frac{\rho}{2650}\right)^{-1/3.4}$$

For the DPRK 2017 test:
- Yield: 250 kt
- Density: 2700 kg/m³
- **Final cavity radius**: ~280 m
- **Cavity volume**: ~92 million m³

#### Phase 4: Cavity Collapse (1 s–1 hour)

After maximum expansion:
1. Cavity pressure drops below lithostatic
2. Weakened roof rock begins to fall
3. Rubble chimney propagates upward
4. Surface subsidence occurs if chimney reaches surface

For Mt. Mantap:
- Chimney reached near-surface
- InSAR measured ~0.5 m subsidence
- Collapse visible in satellite imagery

### Seismic Source Radiation

#### The Mueller-Murphy Model

FSRM uses the Mueller-Murphy (1971) **Reduced Displacement Potential (RDP)** to model the seismic source:

$$\psi(t) = \psi_\infty \cdot \left[1 - e^{-t/\tau_1}\left(1 + \frac{t}{\tau_1}\right)\right] \cdot e^{-t/\tau_2}$$

**Physical interpretation**:
- $\psi_\infty$ = steady-state potential (proportional to final cavity volume)
- $\tau_1$ = rise time (cavity expansion duration)
- $\tau_2$ = overshoot decay time (cavity oscillation damping)

#### Source Spectrum

The far-field P-wave displacement spectrum:

$$\hat{u}(f) \propto \frac{f^2}{(1 + (f/f_c)^2)^2}$$

Key spectral features:
- **Low frequency**: Amplitude ∝ cavity volume ∝ yield
- **Corner frequency**: $f_c \approx 0.3$ Hz for 250 kt
- **High frequency**: Roll-off as $f^{-2}$

#### Moment Tensor Representation

The seismic moment tensor for the DPRK 2017 test:

$$M = M_0 \begin{pmatrix} 1 + \epsilon & 0 & 0 \\ 0 & 1 + \epsilon & 0 \\ 0 & 0 & 1 - 2\epsilon \end{pmatrix} + M_{DC}$$

Where:
- $M_0$ = isotropic moment ≈ 5×10¹⁷ N·m
- $\epsilon$ = CLVD parameter (~0.1 for this event)
- $M_{DC}$ = double-couple component from tectonic release

The relatively large CLVD component reflects the Mt. Mantap collapse.

### Seismic Wave Propagation

#### Crustal Phase Propagation

At regional distances (300–2000 km), waves separate into distinct phases:

| Phase | Description | Velocity | Physical Path |
|-------|-------------|----------|---------------|
| **Pn** | Mantle-refracted P | 8.0–8.2 km/s | Dives to Moho, refracts along mantle lid |
| **Pg** | Crustal P | 5.5–6.5 km/s | Direct path through crust |
| **Sn** | Mantle-refracted S | 4.4–4.7 km/s | S-wave along Moho |
| **Lg** | Crustal-guided S | 3.2–3.6 km/s | Multiple S reflections in crustal waveguide |

#### Travel Time Equations

For a layered Earth model:

**Pg (direct crustal P)**:
$$t_{Pg} = \frac{\Delta}{V_c}$$

**Pn (head wave along Moho)**:
$$t_{Pn} = \frac{2h\cos(i_c)}{V_c} + \frac{\Delta - 2h\tan(i_c)}{V_m}$$

Where:
- $\Delta$ = epicentral distance
- $V_c$ = crustal P velocity (~6.1 km/s)
- $V_m$ = mantle P velocity (~8.1 km/s)
- $h$ = crustal thickness (~35 km)
- $i_c$ = critical angle = $\arcsin(V_c/V_m)$

The **Pn crossover distance** (where Pn becomes first arrival) is ~200 km for typical continental crust.

#### Amplitude Decay

Wave amplitudes decay with distance due to:

1. **Geometric spreading**:
   - Body waves: $A \propto 1/r$
   - Surface waves: $A \propto 1/\sqrt{r}$

2. **Intrinsic attenuation** (anelastic absorption):
   $$A(r) = A_0 \cdot e^{-\pi f r / (Q \cdot v)}$$
   
   Where $Q$ is the quality factor (~300–500 for East Asian crust)

3. **Scattering** from crustal heterogeneities

### Yield Estimation Methods

#### Body-Wave Magnitude (mb)

The standard mb-yield relation (Murphy 1981):

$$m_b = 4.45 + 0.75 \log_{10}(W)$$

For the DPRK 2017 test:
- Observed: $m_b$ ≈ 6.1–6.3
- Inverted yield: 150–350 kt
- **Best estimate**: ~250 kt

#### Surface-Wave Magnitude (Ms)

$$M_s = 3.30 + 1.35 \log_{10}(W)$$

Ms is lower than expected from mb for explosions (the "Ms:mb discriminant").

#### Corner Frequency Method

From spectral analysis:
$$f_c = 1.5 \cdot W^{-0.3}$$

For $f_c$ ≈ 0.3 Hz: W ≈ 200–300 kt

### Explosion vs. Earthquake Discrimination

Key discriminants identifying this event as an explosion:

1. **Ms:mb Ratio**: Low Ms relative to mb (explosions have weak surface waves)
2. **P/S Ratio**: High P-wave to S-wave amplitude (explosions are compressional)
3. **First Motion**: Compressional at all azimuths (earthquakes show quadrantal pattern)
4. **Spectral Ratio**: Elevated high-frequency P-wave content
5. **Depth**: Shallow (~760 m), inconsistent with tectonic earthquakes
6. **Location**: At known nuclear test site

## Recording Stations

The following stations recorded the event with good signal-to-noise ratio:

| Station | Network | Location | Distance (km) | Azimuth |
|---------|---------|----------|---------------|---------|
| MDJ | IC | Mudanjiang, China | 371 | 6° |
| MAJO | IU | Matsushiro, Japan | 951 | 121° |
| BJT | IC | Baijiatuan, China | 1100 | 267° |
| HIA | IC | Hailar, China | 1148 | 324° |
| ERM | II | Erimo, Japan | 1174 | 82° |
| YSS | IU | Yuzhno, Russia | 1260 | 56° |
| SSE | IC | Shanghai, China | 1335 | 215° |
| ULN | IU | Ulaanbaatar, Mongolia | 1887 | 300° |
| ENH | IC | Enshi, China | 2144 | 241° |
| AAK | II | Ala Archa, Kyrgyzstan | 4446 | 291° |

## Generated Figures

All figures are located in `figures/dprk_2017/`:

### Station Map (`station_map.png`)

![Station Map](../../figures/dprk_2017/station_map.png)

PyGMT map showing:
- Event epicenter (red star) at Punggye-ri
- Recording stations (inverted triangles) with network codes
- Topography and coastlines
- Distance annotations

### Waveform Gather (`observed_waveform_gather.png`)

Standard record section showing:
- All station vertical-component (BHZ) waveforms
- Sorted by epicentral distance
- Bandpass filtered 0.5–8.0 Hz
- Instrument response removed (velocity output)
- Full time window: 60 s before to 600 s after origin

### Distance-Scaled Waveform Gather (`observed_waveform_gather_distance_scaled.png`)

![Distance-Scaled Waveform Gather](../../figures/dprk_2017/observed_waveform_gather_distance_scaled.png)

Publication-quality record section with:
- Y-axis scaled to true epicentral distance
- Theoretical phase arrival lines:
  - **Pn** (blue solid): 8.1 km/s mantle P-wave
  - **Pg** (green dashed): 6.1 km/s crustal P-wave
  - **Sn** (red solid): 4.6 km/s mantle S-wave
  - **Lg** (orange dashed): 3.5 km/s crustal S-wave
- Waveform amplitudes normalized and scaled for visual clarity
- Positive/negative wiggle fills for phase identification

### Zoomed Waveform Gather (`observed_waveform_gather_zoomed.png`)

Focused view of the signal window:
- Truncated time axis to emphasize expected arrivals
- Enhanced visibility of first-arriving phases
- Station labels with distance annotations

### Long-Period Waveform Gather (`observed_waveform_gather_long_period.png`)

Record section with long-period filtering:
- Bandpass: 0.02–0.10 Hz (10–50 s period)
- Standard band for moment tensor inversion
- Emphasizes surface waves and long-period body waves
- Used for determining source mechanism

### Spectrograms (`observed_spectrograms.png`)

Time-frequency analysis showing:
- Spectrogram for each station
- Frequency range 0.5–10 Hz
- Arrival time picking assistance
- Frequency content evolution along ray path

### Key Station Close-up (`observed_closeup_key_stations.png`)

Detailed view of nearest/highest-SNR stations:
- Full waveform with P, S, and surface wave windows marked
- Annotated arrival times
- Useful for manual phase picking

### Moment Tensor Results (`moment_tensor_inversion.png`, `moment_tensor_summary.png`)

Results of waveform inversion:
- Full moment tensor solution
- Decomposition into ISO, CLVD, DC components
- Focal mechanism (beachball) representation
- Explained variance and fit quality

## FSRM Configuration

The FSRM simulation uses the configuration file:

```
config/dprk_2017_punggye_ri_coupled_granite.config
```

### Key Parameters

```ini
[SIMULATION]
name = dprk_2017_punggye_ri_coupled_granite
end_time = 60.0                      # 60 s simulation
solid_model = ELASTOPLASTIC_DP       # Drucker-Prager plasticity
enable_elastodynamics = true
use_discontinuous_galerkin = true
dg_order = 3
use_ader_time_integration = true

[EXPLOSION_SOURCE]
type = NUCLEAR_UNDERGROUND
yield_kt = 250.0
depth_of_burial = 760.0              # m
cavity_radius = 280.0                # m (empirical scaling)
crushed_zone_radius = 700.0          # m
fractured_zone_radius = 1400.0       # m
```

### Running the Simulation

```bash
# Run the coupled near-field simulation
mpirun -np 16 ./build/fsrm -c config/dprk_2017_punggye_ri_coupled_granite.config

# Generate synthetic seismograms at regional distances
mpirun -np 8 ./build/fsrm -c config/dprk_2017_punggye_ri_far_field_seismograms.config
```

## Data Processing Pipeline

### 1. Waveform Download

```bash
python3 scripts/fetch_dprk_2017_waveforms.py
```

This script:
- Connects to IRIS FDSN web services
- Downloads BHZ channel data for each station
- Retrieves station metadata (coordinates, instrument response)
- Saves raw and processed waveforms

### 2. Station Map Generation

```bash
python3 scripts/plot_dprk_2017_station_map.py
```

Uses PyGMT to create a topographic map with:
- SRTM30 elevation data
- Station locations from IRIS metadata
- Great-circle paths to event

### 3. Waveform Processing

The processing pipeline includes:
1. **Detrend**: Remove mean and linear trend
2. **Taper**: Apply 5% cosine taper
3. **Response Removal**: Deconvolve instrument response to velocity
4. **Bandpass Filter**: 0.5–8.0 Hz (4-pole Butterworth, zero-phase)
5. **Normalization**: Optional trace-by-trace normalization for display

### 4. Moment Tensor Inversion

Long-period waveforms (0.02–0.10 Hz) are inverted for the full moment tensor using:
- Green's functions from 1-D velocity model
- Deviatoric + isotropic components
- Grid search over depth and mechanism

## Scientific Interpretation

### Phase Identification

The distance-scaled record section allows identification of:

1. **Pn Phase**: First arrival at all distances; propagates along the Moho at ~8.1 km/s
2. **Pg Phase**: Secondary P arrival through the crust at ~6.1 km/s; dominant at <500 km
3. **Sn Phase**: Mantle-refracted S-wave at ~4.6 km/s
4. **Lg Phase**: Crustal-guided S-wave at ~3.5 km/s; high amplitude, long duration

### Explosion Characteristics

Key observations supporting an explosion source:
- Strong Pn/Pg arrivals relative to S-waves
- Impulsive P-wave onset
- Elevated high-frequency content (>2 Hz)
- Positive first motion at all azimuths (compressional)

### Yield Estimation

The yield can be estimated from:
- Body-wave magnitude (mb): Uses Murphy (1981) scaling
- Surface-wave magnitude (Ms): Requires longer-period data
- Corner frequency: Inversely proportional to cavity radius

For this event, multiple studies converge on ~250 kt using various methods.

## References

1. **Mueller & Murphy (1971)** — Seismic characteristics of underground nuclear detonations. Part I: Seismic source. *BSSA* 61(6), 1675-1692.

2. **Zhao et al. (2017)** — Yield estimation of the 2017 DPRK nuclear test. *GRL* 44(20), 10-138.

3. **Kim & Richards (2007)** — Large earthquake locations and source parameters of underground nuclear explosions. *BSSA* 97(1B), 130-142.

4. **Wei (2017)** — InSAR observations of the 2017 Punggye-ri nuclear test. *GRL* 44(19), 9739-9748.

5. **Yao et al. (2018)** — Source characteristics of the 2017 DPRK nuclear test. *SRL* 89(6), 2078-2084.

6. **Patton (1988)** — Source models of the Harzer explosion from regional observations of fundamental-mode and higher-mode surface waves. *BSSA* 78(5), 1133-1157.
