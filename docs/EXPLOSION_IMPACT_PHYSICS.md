# Explosion and Impact Physics Models

FSRM includes comprehensive physics models for simulating explosions, nuclear detonations, and hypervelocity impacts. These models capture near-field damage, seismic wave generation, radiation transport, and electromagnetic pulse effects.

---

## Table of Contents

1. [Overview](#overview)
2. [Underground Nuclear Tests](#underground-nuclear-tests)
3. [Atmospheric Nuclear Detonations](#atmospheric-nuclear-detonations)
4. [Impact Events](#impact-events)
5. [Radiation Transport](#radiation-transport)
6. [Electromagnetic Pulse (EMP)](#electromagnetic-pulse-emp)
7. [Configuration Reference](#configuration-reference)
8. [Example Simulations](#example-simulations)
9. [Validation and Verification](#validation-and-verification)

---

## Overview

The explosion and impact physics module provides integrated multi-physics simulation capabilities for:

| Event Type | Key Physics | Applications |
|------------|-------------|--------------|
| Underground Nuclear Test | Cavity formation, damage zones, seismic coupling | Test monitoring, CTBT verification |
| Atmospheric Nuclear Detonation | Fireball, blast wave, thermal radiation | Effects assessment, historical analysis |
| Surface Explosion | Crater formation, ground shock | Engineering safety, forensics |
| Hypervelocity Impact | Shock metamorphism, crater excavation | Planetary defense, astrogeology |

### Energy Scaling

The physics spans many orders of magnitude:

| Event | Energy (Joules) | TNT Equivalent | Example |
|-------|-----------------|----------------|---------|
| Small meteorite (1 m) | 10^9 | 0.2 tons | Common bolide |
| Large HE explosion | 10^12 | 1 kt | Beirut 2020 |
| Nuclear weapon | 10^13 - 10^17 | 1 kt - 50 Mt | Hiroshima to Tsar Bomba |
| Chelyabinsk meteor (20 m) | 10^15 | 500 kt | 2013 event |
| Tunguska (60 m) | 10^16 | 15 Mt | 1908 event |
| Meteor Crater (50 m iron) | 10^16 | 10-20 Mt | Arizona |
| 1 km asteroid | 10^20 | 70 Gt | Regional extinction |
| Chicxulub (10 km) | 10^24 | 100 Tt | K-Pg extinction |

---

## Underground Nuclear Tests

Underground nuclear explosions create a sequence of phenomena that can be modeled with FSRM.

### Phenomenology

1. **Detonation** (< 1 μs): Nuclear reactions release energy
2. **Cavity Formation** (1 μs - 100 ms): Vaporization and melting create cavity
3. **Shock Propagation** (100 ms - 10 s): Shock wave expands through rock
4. **Chimney Collapse** (seconds - minutes): Cavity roof collapses
5. **Seismic Wave Generation**: P, S, and surface waves radiate outward

### Cavity and Damage Zones

The explosion creates concentric zones of decreasing damage:

```
        ┌──────────────────────────────────────┐
        │           Undamaged Rock             │
        │  ┌────────────────────────────────┐  │
        │  │       Fractured Zone           │  │
        │  │  ┌──────────────────────────┐  │  │
        │  │  │     Crushed Zone        │  │  │
        │  │  │  ┌────────────────────┐  │  │  │
        │  │  │  │      Cavity      │  │  │  │
        │  │  │  │    (melt pool)   │  │  │  │
        │  │  │  └────────────────────┘  │  │  │
        │  │  └──────────────────────────┘  │  │
        │  └────────────────────────────────┘  │
        └──────────────────────────────────────┘
```

Scaling relations (in granite, W in kilotons):

| Zone | Radius (m) | Formula |
|------|------------|---------|
| Cavity | 55 × W^0.3 | Empirical |
| Crushed | ~2 × cavity | Pervasive fracturing |
| Fractured | ~5 × cavity | Radial fractures |
| Damaged | ~10 × cavity | Detectable damage |

### Seismic Source Model (Mueller-Murphy)

The seismic source is characterized by the Reduced Displacement Potential:

$$\psi(\omega) = \frac{K \cdot W^n}{1 + (\omega/\omega_c)^2}$$

Where:
- K = scaling constant (medium-dependent)
- W = yield (kt)
- ω_c = corner frequency
- n ≈ 0.75-0.85 (yield exponent)

The moment tensor is predominantly **isotropic** with some CLVD (cavity collapse) and DC (tectonic release):

$$M_{ij} = M_{iso} \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} + M_{CLVD} \begin{pmatrix} -1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 2 \end{pmatrix} + M_{DC} \begin{pmatrix} ... \end{pmatrix}$$

### Configuration Example

```ini
[EXPLOSION_SOURCE]
type = NUCLEAR_UNDERGROUND
yield_kt = 150.0
depth_of_burial = 1000.0
location_x = 25000.0
location_y = 25000.0
location_z = -1000.0

source_model = MUELLER_MURPHY
corner_frequency_scaling = PATTON
isotropic_fraction = 0.70
clvd_fraction = 0.25
double_couple_fraction = 0.05
```

---

## Atmospheric Nuclear Detonations

Atmospheric nuclear explosions involve coupled phenomena across multiple domains.

### Energy Partition

For a typical thermonuclear weapon in air:

| Energy Channel | Fraction | Physical Mechanism |
|----------------|----------|-------------------|
| Blast (kinetic) | 50% | Hydrodynamic shock wave |
| Thermal radiation | 35% | Blackbody emission from fireball |
| Initial nuclear radiation | 5% | Prompt gamma and neutrons |
| Residual radiation | 10% | Fallout (fission products) |

For high-altitude bursts, EMP becomes significant.

### Fireball Evolution

The fireball goes through distinct phases:

1. **Initial Flash** (< 1 ms): X-ray emission, absorbed by air
2. **First Minimum** (~10 ms): Shock overtakes optical front
3. **Second Maximum** (~0.5 s): Shock becomes transparent
4. **Cooling** (seconds): Gradual temperature decrease

Maximum fireball radius (for 1 kt at sea level):

$$R_{max} \approx 66 \times W^{0.4} \text{ meters}$$

### Blast Wave Model

The Sedov-Taylor solution describes the expanding shock:

$$R(t) = \left(\frac{E}{\rho_0}\right)^{1/5} t^{2/5}$$

Peak overpressure scaling (Glasstone-Dolan):

| Overpressure | Distance Scaling | Effects |
|--------------|------------------|---------|
| 35 psi (240 kPa) | W^1/3 × 0.4 km/kt^1/3 | Severe structural damage |
| 10 psi (70 kPa) | W^1/3 × 0.8 km/kt^1/3 | Most buildings destroyed |
| 5 psi (35 kPa) | W^1/3 × 1.2 km/kt^1/3 | Moderate damage |
| 1 psi (7 kPa) | W^1/3 × 4.0 km/kt^1/3 | Light damage, broken windows |

### Atmospheric Coupling

FSRM couples the atmospheric blast with ground response:

```
        Atmosphere              Ground
       ┌─────────────┐       ┌─────────────┐
       │ Compressible│       │Elastodynamic│
       │ Navier-     │◄─────►│  +          │
       │ Stokes      │       │ Plasticity  │
       └─────────────┘       └─────────────┘
              │                     │
              ▼                     ▼
       Shock pressure ──► Ground motion & cratering
```

### Configuration Example

```ini
[NUCLEAR_DEVICE]
type = THERMONUCLEAR
yield_kt = 500.0
burst_height = 500.0
fission_fraction = 0.50
location_x = 50000.0
location_y = 50000.0

[FIREBALL]
model = BRODE
thermal_fraction = 0.35

[BLAST_WAVE]
model = SEDOV_TAYLOR
enable_ground_reflection = true
ground_reflection_factor = 1.8
```

---

## Impact Events

Hypervelocity impacts follow distinct phases with scaling laws derived from experiments and planetary observations.

### Impact Cratering Stages

1. **Contact & Compression** (microseconds)
   - Impactor contacts target
   - Shock waves propagate through both bodies
   - Peak pressures of 100s of GPa

2. **Excavation** (seconds)
   - Material flows outward and upward
   - Transient crater forms
   - Ejecta launched ballistically

3. **Modification** (seconds to minutes)
   - Crater walls collapse
   - Central uplift forms (complex craters)
   - Final crater morphology established

### Pi-Group Scaling

Crater dimensions follow the Holsapple-Schmidt scaling laws:

$$\frac{V}{m} = K_1 \left(\frac{\rho_t}{\rho_i}\right)^{1-3\nu} \left(\frac{Y}{\rho_t U^2}\right)^{-\frac{3\mu}{2}} \left(\frac{gR}{U^2}\right)^{-3\mu/2}$$

Where:
- V = crater volume
- m = impactor mass
- ρ = densities (target, impactor)
- Y = target strength
- U = impact velocity
- g = surface gravity
- R = impactor radius
- K₁, μ, ν = scaling constants

### Simple vs. Complex Craters

| Parameter | Simple | Complex |
|-----------|--------|---------|
| Transition diameter (Earth) | ~4 km | > 4 km |
| Depth/Diameter | 0.2 | 0.1 |
| Central peak | No | Yes |
| Terraced walls | No | Yes |
| Formation mechanism | Excavation only | + Gravitational collapse |

### Shock Metamorphism

Shock waves produce diagnostic features at different pressure levels:

| Pressure (GPa) | Effects |
|----------------|---------|
| 2-5 | Shatter cones |
| 10-25 | Planar deformation features (PDFs) in quartz |
| 25-50 | Diaplectic glass, high-pressure polymorphs |
| 50-80 | Incipient melting |
| > 80 | Complete melting |
| > 150 | Vaporization |

### Seismic Waves from Impact

Impact generates seismic waves with:

- **Seismic efficiency**: η ≈ 10^-5 to 10^-4 (fraction of KE)
- **Source mechanism**: Predominantly isotropic (explosion-like)
- **Equivalent magnitude**: M_w ≈ (2/3) log₁₀(η × KE) - 6.07

For a 1 km asteroid (KE ≈ 3 × 10^20 J):
- Seismic energy: ~3 × 10^16 J
- Equivalent magnitude: M_w ≈ 7-8

### Configuration Example

```ini
[IMPACTOR]
type = ROCKY
diameter = 1000.0
density = 3500.0
impact_velocity = 18000.0
impact_angle = 45.0
impact_x = 100000.0
impact_y = 100000.0

[CRATER]
model = PI_SCALING
scaling_law = HOLSAPPLE_ROCK
enable_central_peak = true
enable_crater_collapse = true
```

---

## Radiation Transport

### Prompt Radiation

Nuclear detonations emit prompt gamma rays and neutrons:

**Gamma Radiation:**
$$D_\gamma(r) = \frac{E_\gamma}{4\pi r^2} \cdot e^{-\mu r} \cdot B(\mu r)$$

Where B is the buildup factor for atmospheric scattering.

**Neutron Radiation:**
$$D_n(r) = \frac{N_0}{4\pi r^2} \cdot e^{-\Sigma r}$$

Where Σ is the macroscopic cross-section.

### Fallout Model

Fission products follow the Way-Wigner decay law:

$$A(t) = A_1 \cdot t^{-1.2}$$

Key isotopes:
- **I-131** (t₁/₂ = 8 days): Short-term hazard
- **Cs-137** (t₁/₂ = 30 years): Long-term contamination
- **Sr-90** (t₁/₂ = 29 years): Bone seeker

The HOTSPOT model is used for fallout deposition:

```ini
[RADIATION_TRANSPORT.FALLOUT]
enabled = true
particle_model = LOGNORMAL
median_particle_size = 100.0
deposition_model = HOTSPOT
```

---

## Electromagnetic Pulse (EMP)

Nuclear detonations, especially at high altitude, generate electromagnetic pulses through three mechanisms.

### E1 (Fast Component)

**Mechanism:** Prompt gamma rays Compton-scatter electrons in the upper atmosphere. These electrons spiral in Earth's magnetic field, creating a time-varying current that radiates.

**Characteristics:**
- Rise time: ~2.5 ns
- Peak field: 20,000-50,000 V/m
- Duration: ~1 μs
- Dominant at high altitude

### E2 (Intermediate Component)

**Mechanism:** Scattered gamma rays and neutrons continue to ionize the atmosphere.

**Characteristics:**
- Rise time: ~1 μs
- Peak field: ~100 V/m
- Duration: ~1 ms
- Similar to lightning EMP

### E3 (MHD-EMP)

**Mechanism:** The expanding conducting fireball distorts Earth's magnetic field.

**Characteristics:**
- Rise time: ~1 s
- Peak field: ~1-10 V/m
- Duration: ~100 s
- Couples to long conductors (power lines)

### EMP Configuration

```ini
[ELECTROMAGNETIC_PULSE]
enabled = true

[EMP.E1]
enabled = true
rise_time = 2.5e-9
peak_field = 50000.0
geomagnetic_field = 50.0e-6
geomagnetic_dip = 60.0

[EMP.E2]
enabled = true
duration = 1.0e-3

[EMP.E3]
enabled = true
fireball_conductivity = 1.0e6
magnetic_field_distortion = true

[EMP.GROUND_EFFECTS]
ground_conductivity = 0.01
enable_line_coupling = true
```

---

## Configuration Reference

### [EXPLOSION_SOURCE] Section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| type | enum | - | NUCLEAR_UNDERGROUND, NUCLEAR_ATMOSPHERIC, CHEMICAL, VOLCANIC |
| yield_kt | float | 100.0 | Explosive yield in kilotons TNT |
| depth_of_burial | float | 0.0 | Depth below surface (m) |
| location_x, y, z | float | 0.0 | Source location |
| source_model | enum | MUELLER_MURPHY | Seismic source model |
| cavity_radius | float | auto | Override cavity radius (m) |
| isotropic_fraction | float | 0.7 | Moment tensor decomposition |

### [IMPACTOR] Section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| type | enum | STONY | STONY, IRON, COMETARY, ICY |
| diameter | float | 1000.0 | Impactor diameter (m) |
| density | float | 3000.0 | Impactor density (kg/m³) |
| impact_velocity | float | 20000.0 | Impact velocity (m/s) |
| impact_angle | float | 45.0 | Angle from horizontal (degrees) |

### [RADIATION_TRANSPORT] Section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| enabled | bool | true | Enable radiation physics |
| method | enum | MONTE_CARLO | Transport method |
| num_particles | int | 1e6 | MC particles (if applicable) |

### [ELECTROMAGNETIC_PULSE] Section

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| enabled | bool | true | Enable EMP physics |
| source_mechanism | enum | COMPTON | EMP source type |
| geomagnetic_field | float | 50e-6 | Earth's field (T) |
| ground_conductivity | float | 0.01 | Soil conductivity (S/m) |

---

## Example Simulations

### Example 1: Underground Nuclear Test (Nevada-style)

```bash
# 150 kt contained underground test
mpirun -np 16 fsrm -c config/underground_nuclear_test.config

# Outputs:
# - Seismic waveforms at receiver stations
# - Damage zone evolution
# - Cavity and chimney formation
# - mb/Ms discrimination ratios
```

### Example 2: Atmospheric Detonation

```bash
# 500 kt airburst with radiation and EMP
mpirun -np 64 fsrm -c config/atmospheric_nuclear_test.config

# Outputs:
# - Blast overpressure contours
# - Thermal fluence maps
# - Prompt radiation doses
# - EMP field time histories
# - Fallout patterns (if surface burst)
```

### Example 3: Asteroid Impact

```bash
# 1 km asteroid impact
mpirun -np 128 fsrm -c config/impact_event.config

# Outputs:
# - Crater evolution movie
# - Shock metamorphism zones
# - Ejecta distribution
# - Seismic waveforms (global propagation)
# - Atmospheric plume dynamics
```

---

## Validation and Verification

### Nuclear Explosion Validation

| Benchmark | Description | Reference |
|-----------|-------------|-----------|
| Gnome | 3.1 kt underground, salt | NV, 1961 |
| Gasbuggy | 29 kt underground | NM, 1967 |
| Sedan | 104 kt crater | NTS, 1962 |
| mb-yield scaling | Seismic magnitude vs yield | Murphy (1996) |

### Impact Validation

| Benchmark | Description | Reference |
|-----------|-------------|-----------|
| Pi-scaling | Crater scaling laws | Holsapple (1993) |
| Meteor Crater | 50 m iron impact | Arizona |
| Chicxulub | 10 km asteroid | K-Pg boundary |
| Laboratory impacts | Vertical gun tests | NASA Ames |

### Analytical Solutions

- **Sedov-Taylor** blast wave (atmospheric)
- **Garvin's problem** (buried explosion in halfspace)
- **Lamb's problem** (point force on surface)

---

## References

1. **Nuclear Effects:**
   - Glasstone, S. & Dolan, P.J. (1977). *The Effects of Nuclear Weapons*
   - Mueller, R.A. & Murphy, J.R. (1971). Seismic characteristics of underground nuclear detonations

2. **Impact Cratering:**
   - Melosh, H.J. (1989). *Impact Cratering: A Geologic Process*
   - Holsapple, K.A. (1993). The scaling of impact processes in planetary sciences

3. **EMP:**
   - Longmire, C.L. (1978). On the electromagnetic pulse produced by nuclear explosions
   - Savage, E. et al. (2010). The Early-Time (E1) High-Altitude EMP

4. **Radiation Transport:**
   - Bridgman, C.J. (2001). *Introduction to the Physics of Nuclear Weapons Effects*

---

## See Also

- [Wave Propagation](./WAVE_PROPAGATION.md) - Seismic wave physics
- [Plasticity Models](./PLASTICITY_MODELS.md) - Material failure
- [Configuration Reference](./CONFIGURATION.md) - Full config options
- [Benchmarks](./BENCHMARKS.md) - Validation suite
