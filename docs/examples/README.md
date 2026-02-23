# FSRM Examples: Historical Nuclear Test Simulations

This directory contains detailed documentation for FSRM's historical nuclear test simulation examples. These examples demonstrate the simulator's capabilities for modeling underground nuclear explosions and validating synthetic seismograms against real recorded data.

## Available Examples

| Example | Date | Location | Yield | Description |
|---------|------|----------|-------|-------------|
| [DPRK 2017](DPRK_2017_PUNGGYE_RI.md) | 2017-09-03 | Punggye-ri, DPRK | ~250 kt | Largest DPRK nuclear test (coupled) |
| [US Divider 1992](US_DIVIDER_1992.md) | 1992-09-23 | NTS Area 4, Nevada | ~1 kt | Last US underground test (tuff medium) |
| [Lop Nor 2020](LOP_NOR_2020_DECOUPLED.md) | 2020-06-22 | Lop Nor, China | ~1 kt | Alleged decoupled test |

### Comparison of Test Scenarios

These three examples span a range of testing scenarios with distinct physics:

| Aspect | DPRK 2017 | US Divider 1992 | Lop Nor 2020 |
|--------|-----------|-----------------|--------------|
| **Coupling** | Fully coupled | Coupled | Decoupled (cavity) |
| **Medium** | Granite | Volcanic tuff | Granite |
| **Yield** | ~250 kt | ~1 kt | ~1 kt (apparent: ~10 t) |
| **Depth** | 760 m | 800 m | 300 m |
| **mb** | 6.1–6.3 | ~4.0 | ~2.75 |
| **Key Physics** | Cavity collapse, CLVD | Porous media effects | Decoupling factor |

## Overview

Each example includes:

1. **Event Parameters** — Origin time, location, depth, yield estimates
2. **Geology & Site Characterization** — Host rock, velocity structure, burial conditions
3. **FSRM Configuration** — Simulation setup, source model, mesh specifications
4. **Observational Data** — Real seismic waveforms from IRIS/FDSN
5. **Generated Figures** — Station maps, waveform gathers, spectrograms, moment tensors

## Running the Examples

### Prerequisites

```bash
# Install Python dependencies
pip install obspy pygmt matplotlib numpy scipy pandas

# Ensure FSRM is built
cd /workspaces/FSRM
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Workflow

Each example follows a standard workflow:

1. **Fetch Observational Data**
   ```bash
   python3 scripts/fetch_dprk_2017_waveforms.py
   python3 scripts/fetch_us_divider_1992_waveforms.py
   python3 scripts/fetch_lop_nor_2020_waveforms.py
   ```

2. **Generate Station Maps**
   ```bash
   python3 scripts/plot_dprk_2017_station_map.py
   python3 scripts/plot_us_divider_1992_station_map.py
   python3 scripts/plot_lop_nor_2020_station_map.py
   ```

3. **Run FSRM Simulation**
   ```bash
   mpirun -np 8 ./build/fsrm -c config/dprk_2017_punggye_ri_coupled_granite.config
   mpirun -np 8 ./build/fsrm -c config/us_divider_1992_far_field_seismograms.config
   mpirun -np 8 ./build/fsrm -c config/lop_nor_2020_decoupled_granite.config
   ```

4. **Compare Synthetics to Observations**
   ```bash
   python3 scripts/visualize_fsrm_dprk_2017.py
   python3 scripts/visualize_fsrm_us_divider.py
   python3 scripts/visualize_fsrm_lop_nor.py
   ```

## Key Figures

All figures are generated in `figures/<example_name>/`:

| Figure | Description |
|--------|-------------|
| `station_map.png` | PyGMT map showing event epicenter and recording stations |
| `observed_waveform_gather.png` | Record section of all station waveforms |
| `observed_waveform_gather_distance_scaled.png` | Distance-scaled record section with theoretical phase arrivals |
| `observed_waveform_gather_zoomed.png` | Zoomed view focused on expected signal window |
| `observed_waveform_gather_long_period.png` | Long-period filtered (0.02–0.1 Hz) for moment tensor analysis |
| `observed_spectrograms.png` | Time-frequency spectrograms for each station |
| `observed_closeup_key_stations.png` | Detailed view of closest/highest-SNR stations |
| `moment_tensor_inversion.png` | Full moment tensor inversion results |
| `fig01_source_model.png` | Mueller-Murphy source model visualization |
| `fig02_velocity_model.png` | Crustal velocity structure |
| `fig03_synthetic_seismograms.png` | FSRM synthetic waveforms |
| `fig07_comparison.png` | Synthetic vs. observed waveform comparison |
| `fig08_discrimination.png` | Explosion/earthquake discrimination analysis |

## Scientific Background

### Underground Nuclear Explosions

Underground nuclear tests produce distinctive seismic signatures:

- **Isotropic Source**: Explosions are predominantly spherically symmetric, producing strong compressional (P) waves
- **Shallow Depth**: Most tests are at depths of 0.5–1.5 km, generating surface waves and spall
- **Short Duration**: The source time function is typically <1 second
- **High-Frequency Content**: Explosions generate more high-frequency energy than earthquakes of similar magnitude

### Detonation Phenomenology

An underground nuclear explosion progresses through distinct phases:

1. **Nuclear Reactions** (0–1 μs): Energy release at 10⁷–10⁸ K
2. **Shock Wave & Vaporization** (1 μs–10 ms): Rock destroyed, cavity expands
3. **Cavity Expansion** (10 ms–1 s): Pressure-driven growth to final radius
4. **Rebound & Collapse** (1 s–hours): Cavity stabilization or chimney formation

### Cavity Radius Scaling

The final cavity radius depends on yield and rock properties:

$$R_c = 55 \cdot W^{0.295} \cdot \left(\frac{\rho}{2650}\right)^{-1/3.4}$$

| Yield | Hard Rock (granite) | Soft Rock (tuff) |
|-------|---------------------|------------------|
| 1 kt | ~55 m | ~70 m |
| 10 kt | ~110 m | ~140 m |
| 100 kt | ~220 m | ~280 m |
| 250 kt | ~280 m | ~360 m |

### Mueller-Murphy Source Model

FSRM uses the Mueller-Murphy (1971) reduced displacement potential (RDP) to model the explosion source:

$$\psi(t) = \psi_\infty \cdot \left[1 - e^{-t/\tau_1}\left(1 + \frac{t}{\tau_1}\right)\right] \cdot e^{-t/\tau_2}$$

Where:
- $\psi_\infty$ = steady-state potential (proportional to cavity volume)
- $\tau_1$ = rise time (cavity expansion duration)
- $\tau_2$ = overshoot decay time

### Source Spectrum

The far-field displacement spectrum has characteristic shape:

- **Low frequencies**: Flat, amplitude ∝ cavity volume ∝ yield
- **Corner frequency**: $f_c = 1.5 \cdot W^{-0.3}$ Hz (Patton scaling)
- **High frequencies**: Roll-off as $f^{-2}$

### Coupled vs. Decoupled Tests

| Aspect | Coupled Test | Decoupled Test |
|--------|--------------|----------------|
| Geometry | Device in rock | Device in pre-excavated cavity |
| Seismic Coupling | Full (~100%) | Reduced (~1–3%) |
| Decoupling Factor | 1 | 30–70 |
| Magnitude Reduction | 0 | 1.5–1.9 mb units |
| Example | DPRK 2017 | Lop Nor 2020 |

### Emplacement Medium Effects

Rock type significantly affects explosion signals:

| Medium | Density | Vp | Corner Freq | mb Bias |
|--------|---------|-----|-------------|---------|
| Granite | 2700 kg/m³ | 6.0 km/s | Higher | Baseline |
| Tuff | 2000 kg/m³ | 4.0 km/s | Lower | -0.3 to -0.5 |
| Salt | 2200 kg/m³ | 4.5 km/s | Medium | -0.2 |

### Phase Velocities

Theoretical phase arrival times are calculated using regional crustal velocities:

| Phase | Description | Velocity | Path |
|-------|-------------|----------|------|
| Pn | Mantle-refracted P | 8.0–8.2 km/s | Along Moho |
| Pg | Crustal P | 5.5–6.5 km/s | Direct through crust |
| Sn | Mantle-refracted S | 4.4–4.7 km/s | S-wave along Moho |
| Lg | Crustal-guided S | 3.2–3.6 km/s | Multiple crustal reflections |
| Rg | Rayleigh wave | 2.5–3.0 km/s | Surface wave |

### Explosion-Earthquake Discrimination

Key observables that distinguish explosions from earthquakes:

| Discriminant | Explosions | Earthquakes |
|--------------|------------|-------------|
| Ms:mb ratio | Low (Ms < mb) | High (Ms > mb) |
| P/S amplitude | High | Low |
| First motion | All compressional | Quadrantal |
| High-freq content | Enhanced | Reduced |
| Depth | Shallow (<2 km) | Variable |
| Isotropic moment | >50% | ~0% |

### Yield Estimation Methods

Multiple techniques constrain explosion yield:

1. **Body-wave magnitude**: $m_b = 4.45 + 0.75 \log_{10}(W)$ (hard rock)
2. **Corner frequency**: $f_c = 1.5 \cdot W^{-0.3}$
3. **Lg amplitude**: Regional calibration
4. **Isotropic moment**: Direct measure of cavity volume

## References

- Mueller, R.A. and Murphy, J.R. (1971). Seismic characteristics of underground nuclear detonations. Part I: Seismic source. *Bulletin of the Seismological Society of America*, 61(6), 1675-1692.
- Patton, H.J. (1988). Source models of the Harzer explosion from regional observations of fundamental-mode and higher-mode surface waves. *Bulletin of the Seismological Society of America*, 78(5), 1133-1157.
- Glenn, L.A. and Myers, S.C. (1997). Decoupling of seismic sources by large underground cavities. *Journal of Geophysical Research*, 102(B5), 10065-10082.
- Murphy, J.R. (1981). P-wave coupling of underground explosions in various geologic media. *Bulletin of the Seismological Society of America*, 71(6), 1789-1806.
- Stevens, J.L. and Day, S.M. (1985). The physical basis of mb:Ms and variable frequency magnitude methods for earthquake/explosion discrimination. *Journal of Geophysical Research*, 90(B4), 3009-3020.
- Zhao, L.F., et al. (2017). Yield estimation of the 2017 DPRK nuclear test. *Geophysical Research Letters*, 44(20), 10-138.
