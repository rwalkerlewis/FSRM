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

### Dynamic near-field source (pass-5)

The far-field FEM problem can be driven by either a kinematic moment
rate (the legacy path) or by the recorded history of a 1D
NearFieldExplosionSolver (pass-5). The selection lives in the
`[NEAR_FIELD_SOURCE]` config section:

```ini
[NEAR_FIELD_SOURCE]
mode = DYNAMIC_PLASTIC          # or KINEMATIC_RDP (default)
elastic_radius_factor = 3.0     # extraction surface = factor * Rc
near_field_dt = 1e-5            # solver sub-step
damage_model = DRUCKER_PRAGER
output_cadence_microseconds = 100
```

`KINEMATIC_RDP` is the default and preserves byte-identical output
for configs that omit the section (the
`Integration.NearFieldSource.KinematicRDPLegacyByteIdentical` test
guards this guarantee).

`DYNAMIC_PLASTIC` runs the 1D `NearFieldExplosionSolver` at setup
time with the configured sub-step, samples the full 6-component
moment-rate tensor (with iso + CLVD + double-couple split) at the
configured cadence over a spherical extraction surface at
`elastic_radius_factor * Rc`, and writes the recorded history to
`<seismometer output_dir>/near_field_history.csv`.
`addExplosionSourceToResidual` interpolates from this history and
injects the FULL `Mdot_ij` tensor (not just the isotropic trace) into
the far-field linear-elastic problem.

The fidelity gain over `KINEMATIC_RDP` is:

  (a) the full moment-rate tensor including the CLVD content drives
      the far field, instead of the trace-only injection;
  (b) the elastic-radius extraction surface is configurable and
      reported via the recorded `R_cavity` and `R_plastic` time
      series;
  (c) the history CSV exposes the cavity-expansion and plastic-radius
      time series to downstream visualisation (ParaView, the Sedan
      1962 anchor example committed under
      `examples/11_sedan_1962/paraview/`).

The underlying `M(t)` shape is still RDP-derived in this build (the
1D solver couples its strength model, damage model, and Mie-Gruneisen
EOS to a closed-form cavity-expansion kernel rather than a true
radial Lagrangian solve). Roadmap axis 1 in
`docs/HISTORIC_NUCLEAR_ROADMAP.md` tracks the follow-up to replace
the closed-form kernel with a true 1D radial Lagrangian elastoplastic
shock solver.

### Pass-8 radiation transport phase

Pass-8 (axis 1) replaces the pass-7 Zel'dovich-Raizer end-state
approximation with an explicit numerical solve of the radiation-
transport phase. The new
`MarshakRadiationDiffusionSolver`
(`include/domain/explosion/MarshakRadiationDiffusion.hpp`) runs from
t = 0 (yield deposition) to t = t_rh (radiation-to-hydrodynamic
transition) on the existing 1D radial Lagrangian mesh, coupled to
the Tillotson host-rock matter via emission-absorption.

#### Governing equations

Grey radiation diffusion in 1D radial spherical symmetry:

```
dE_r/dt = (1/r^2) d/dr [ r^2 (c / (3 kappa_R rho)) dE_r/dr ]
        + c kappa_P rho ( a T_m^4 - E_r )

rho cv dT_m/dt = c kappa_P rho ( E_r - a T_m^4 )
```

with `a = 4 sigma_SB / c` the radiation constant, `kappa_R` the
Rosseland mean opacity (controls diffusion), `kappa_P` the Planck
mean opacity (controls emission/absorption), `T_m(rho, e)` the
matter temperature from the Tillotson `temperature(rho, e)`
accessor.

#### Numerical scheme

Operator-split per global timestep:

  1. Hydro substep (pass-7 path: momentum, deviatoric stress, EOS).
  2. Radiation-matter coupling: implicit backward-Euler in `E_r`.
     Tridiagonal Thomas solve in 1D. Outer Newton iteration on the
     `T_m^4` closure (Larsen 1988); typical 3-5 iterations to
     converge below `radiation_newton_tolerance` (default 1e-6).
  3. Energy update: per-cell matter energy increment
     `de = cv (T_m_new - T_m_old)` is added to `e_int_`; the host
     re-evaluates the EOS so the matter pressure picks up the
     deposited radiation energy.

Boundary conditions: zero-flux symmetry at `r = 0`; Marshak
reservoir (E_r = a T_amb^4 with T_amb = 300 K) at the outer face.

#### Hand-off criterion

Per substep the solver returns the radiation-front index (outermost
cell with `E_r > 1.5 a T_amb^4`) and `t_diff = (dr)^2 / D` at the
front. The host computes `t_hydro = dr / max(|v|, c_s)` at the same
cell. The radiation phase ends when `t_hydro < t_diff` for
`radiation_handoff_debounce_steps` consecutive substeps. After
hand-off the Marshak path is bypassed; the existing Lagrangian
hydro continues from the post-hand-off state.

#### Opacity model

Power-law (Kramers'-type) parameterization
(`include/domain/explosion/OpacityModel.hpp`):

```
kappa_R(rho, T) = kappa_R_0 (rho / rho_0)^a_R (T / T_0)^b_R
kappa_P(rho, T) = kappa_P_0 (rho / rho_0)^a_P (T / T_0)^b_P
```

Exponents from Zel'dovich-Raizer 1967 vol I sec 10 (free-free
ionized regime: a ~ 1, b ~ -3.5). Hard-coded sets per medium
(granite, tuff, salt, alluvium); the `Physics.Marshak.OpacityRegimeCoverage`
gate exercises all four across rho 1e2-5e3 kg/m^3 and T 1e4-1e7 K.

#### Configuration

```ini
[NEAR_FIELD_SOURCE]
radiation_phase = MARSHAK_GREY        # default ZELDOVICH_RAIZER
opacity_model = POWER_LAW_ZR
radiation_max_newton_iter = 10
radiation_newton_tolerance = 1.0e-6
radiation_handoff_debounce_steps = 3
tillotson_extrapolation_warning_threshold_pa = 5e10
```

`radiation_phase = MARSHAK_MULTIGROUP` and `SN_TRANSPORT` are
named scaffolds; selecting them throws a clear runtime_error at
solver construction time. Tabulated opacities
(`opacity_model = TABULATED_TOPS`) likewise throws.

#### References

  - Pomraning, G. C. (1973), ch IV (Marshak self-similar wave).
  - Mihalas & Mihalas (1984), sec 96-97 (operator splitting,
    backward-Euler stability for stiff coupling).
  - Marshak, R. E. (1958), Phys. Fluids 1(1), pp. 24-29 (boundary
    conditions for radiation diffusion).
  - Zel'dovich & Raizer (1967), vol I ch V sec 10 (Kramers'
    opacity), vol II ch X (radiation-to-hydrodynamic transition).
  - Larsen, E. W. (1988), J. Comp. Phys 78, pp. 459-480
    (linearised Newton on T^4 closure).

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

Scaling relations (granite-like host rock, \(W\) in kilotons TNT equivalent) are implemented in `include/NearFieldExplosion.hpp` and covered by unit tests.

| Zone | Radius (m) | Formula |
|------|------------|---------|
| Cavity | \(55 \times W^{0.295} \times (\rho/2.65)^{-1/3.4}\) | Empirical (NTS-style) |
| Crushed | \(\approx 2.5 \times\) cavity | Pervasive comminution |
| Fractured | \(\approx 5 \times\) cavity | Radial fractures |
| Damaged | \(\approx 10 \times\) cavity | Detectable damage |

See: `tests/unit/test_near_field_explosion.cpp` for regression coverage of these scalings.

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
[SIMULATION]
# Controls how explosions are coupled into the solve.
# - PROXY: reduced-order / triggering-only paths
# - COUPLED_ANALYTIC: integrated analytic coupling (e.g., spherical cavity stress)
# - FULL_PDE: full coupled blast+shock+wave PDE everywhere (not implemented in this build)
explosion_solve_mode = COUPLED_ANALYTIC

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

#### Variants (additional nuclear test example configs)

```bash
# Low-yield, shallow burial (monitoring/discrimination workflow demo)
mpirun -np 4 fsrm -c config/underground_nuclear_test_1kt_shallow.config

# Decoupled / weakly coupled salt-cavity scenario (lower coupling for given yield)
mpirun -np 8 fsrm -c config/underground_nuclear_test_decoupled_salt_10kt.config
```

#### Historical-inspired geology + geography (Gmsh + CRS)

These configs emphasize **site geology** (layered materials) and **geographic CRS**
handling (WGS84 lon/lat → UTM/local meters), suitable as bases for historically
anchored scenarios:

```bash
# Gasbuggy-style layered site model (Gmsh physical volumes + CRS transform)
mpirun -np 4 fsrm -c config/historical_gasbuggy_1967_site_geology_gmsh.config

# Gnome-style salt host model (Gmsh physical volume + CRS transform)
mpirun -np 4 fsrm -c config/historical_gnome_1961_site_geology_salt_gmsh.config

# Sedan-style near-surface alluvium+tuff model (Gmsh + CRS)
mpirun -np 4 fsrm -c config/historical_sedan_1962_site_geology_cratering.config

# Starfish-Prime-style high-altitude EMP context (geographic CRS + EMP knobs)
mpirun -np 64 fsrm -c config/historical_starfish_prime_1962_high_altitude_emp.config

# Faultless-style nearby fault network + induced seismicity trigger (Gmsh faults + IMEX)
mpirun -np 4 fsrm -c config/historical_faultless_1968_fault_network_induced_seismicity.config
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

#### Variants (additional atmospheric nuclear example configs)

```bash
# Surface burst emphasizing fallout + deposition mapping
mpirun -np 32 fsrm -c config/atmospheric_nuclear_test_surface_fallout_50kt.config

# High-altitude burst emphasizing EMP + ionospheric effects
mpirun -np 64 fsrm -c config/atmospheric_nuclear_test_high_altitude_emp_1400kt.config
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

## Pass-6 spatial profile output (RADIAL_LAGRANGIAN)

Pass-6 (axis-1a, see `docs/HISTORIC_NUCLEAR_ROADMAP.md`) lands a 1D
radial Lagrangian elastoplastic shock solver behind a new
`solver_kind` sub-key under `[NEAR_FIELD_SOURCE]`. When
`solver_kind = RADIAL_LAGRANGIAN`, the Simulator records spatial
snapshots of the radial state at the configured cadence and emits
two files alongside the pass-5 history CSV:

- `near_field_profile.h5` -- HDF5 file with the spatial profile
  time series.
- `near_field_profile.xdmf` -- ParaView wrapper around the HDF5
  exposing the cell-centred datasets as cell-data on a polyline
  mesh.

HDF5 schema (frozen for pass-6; consumers depend on it):

```
/time                 (n_snap,)  double, simulation time [s]
/num_cells            scalar     int, cell count N
/profiles/<i>/r        (N+1,)    double, face radii [m]
/profiles/<i>/r_cell   (N,)      double, cell-centred radii [m]
/profiles/<i>/v_r      (N+1,)    double, face velocities [m/s]
/profiles/<i>/rho      (N,)      double, density [kg/m^3]
/profiles/<i>/p        (N,)      double, pressure [Pa]
/profiles/<i>/sigma_rr (N,)      double, total radial stress [Pa]
/profiles/<i>/sigma_tt (N,)      double, total hoop stress [Pa]
/profiles/<i>/eps_p    (N,)      double, equivalent plastic strain
/profiles/<i>/damage   (N,)      double, scalar damage [0,1]
/profiles/<i>/yield_indicator (N,) double, 1.0 if cell yielded
```

The pass-6 solver is a 1D radial finite-volume Lagrangian formulation
in spherical symmetry. Per cell: density, specific internal energy,
pressure (compression positive), total radial and hoop stresses.
Per face: position and radial velocity (staggered grid). The per-step
update sequence runs CFL-bounded explicit time stepping with Wilkins
linear + quadratic artificial viscosity for shock capture, an elastic
predictor for the deviatoric stress, Drucker-Prager radial return on
the existing `PressureDependentStrength` data, Mie-Gruneisen EOS for
solid cells with an ideal-gas EOS for the inner cavity, energy update
including plastic dissipation, damage evolution from the existing
`DamageEvolutionModel`, and an outgoing-characteristic outer
boundary condition. The moment-tensor extraction integrates traction
over a fixed Eulerian sphere at the configured elastic radius;
spherical symmetry collapses the tensor to a purely diagonal
isotropic form, so CLVD content requires the deferred axis-1b 3D
subdomain. Implementation lives in
`src/domain/explosion/RadialLagrangian.cpp` (~600 lines).

`solver_kind = CLOSED_FORM` (kept available; pass-6 default that
pass-7 demoted) preserves the pass-5 RDP-driven path byte-for-byte
and writes only the pass-5 CSV. The spatial-profile HDF5 + XDMF
pair is `RADIAL_LAGRANGIAN`-only because the closed-form kernel is
0D analytic and has no meaningful radial profile.

## Pass-7 inner-cavity physics

Pass-7 closes the pass-6 amplitude calibration gap on axis 1a. The
pass-6 inner-cavity initialization deposited the entire yield as
ideal-gas internal energy in a hand-tuned cavity volume, with
gamma = 1.4. That state was an order-of-magnitude wrong on the
relevant physics: a nuclear cavity at the radiation-to-hydrodynamic
transition time contains vaporized rock plasma at temperatures four
to six orders of magnitude above chemical-detonation conditions, not
chemical-detonation gas. Pass-7 replaces the placeholder with a
real EOS for the host rock and a first-principles solve for the
initial cavity state.

### Tillotson host-rock EOS

The Tillotson form (Tillotson 1962, Melosh 1989 eqs 5.4.7-9) is an
analytic EOS calibrated against shock-Hugoniot data that handles
four physical regimes: cold compressed, cold expanded, hot expanded,
and the mixed (partial vaporization) intermediate. In compressed
cells (rho >= rho_0 or e < E_iv):

```
p = (a + b / (1 + e / (E_0 eta^2))) rho e + A mu + B mu^2
```

with eta = rho / rho_0 and mu = eta - 1. In hot-expanded cells
(rho < rho_0 and e > E_cv) the cold pressure decays exponentially
and the thermal term tends to the ideal-gas form. The mixed regime
is a linear interpolation in e between the two.

Four parameter sets ship in `include/domain/explosion/TillotsonEOS.hpp`:
granite (Melosh Table A2.2), tuff (volcanic-glass scaling fit to
Trunin 2001 shock data), salt (Carter 1979 + Melosh A2.2), and
alluvium (placeholder set documented as a known gap pending
purpose-built fit; pass-8). Each set is tested for thermodynamic
self-consistency by `Physics.TillotsonEOS.*`.

The host-rock parameter set is selected from the existing
`[EXPLOSION_SOURCE] medium_type` unless overridden by
`[NEAR_FIELD_SOURCE] tillotson_parameter_set`.

### First-principles cavity initial state

Premise (Zel'dovich and Raizer 1967, vol II, ch X). At the radiation-
to-hydrodynamic transition time `t_rh` (when radiation transport
stops outpacing hydrodynamic expansion), the cavity contains
fully-vaporized host rock at approximately the solid density: the
radiation wave heats material in place faster than the cavity can
hydrodynamically expand. Pass-7 imposes this Marshak-end-state
approximation as the initial condition (`rho_v = rho_0_solid`).

Energy partition. The total deposited yield equals the sum of
(latent vaporization heat, thermal internal energy of the vapor,
gravitational potential energy of the displaced overburden,
residual kinetic energy zero by definition at `t_rh`):

```
E_yield = m_v * E_cv
        + m_v * (e_v - E_cv)
        + m_v * g * h_eff
        + 0
```

With `rho_v = rho_0_solid` fixed, the cavity mass `m_v = (4/3) pi
R_v^3 rho_v` is the only unknown. Pass-7 solves the resulting
single-equation Newton iteration with target `e_v = E_cv` (just-
vaporized state, where the Tillotson evaluation gives ~50 GPa for
granite, plenty to drive the surrounding shock). Convergence in
5-10 iterations.

The default radiation-transition time is `t_rh ~ 1e-7 * W_kt^(1/3)`
seconds (Z-R vol II eq. 24.18); the user can override via
`radiation_transition_time_s`.

### Wilkins (1980) AV defaults

Pass-6 used `c_l = 0.5, c_q = 2.0` (early-development robustness
choice that over-dissipated the leading shock). Pass-7 defaults to
the Wilkins 1980 production prescription `c_l = 0.06, c_q = 1.5`.
The user can still override via `art_visc_linear` and
`art_visc_quadratic`.

### Result

Sedan 1962 anchor (104 kt alluvium, 194 m depth) under pass-7
defaults lands at ratio ~ 2.24x relative to the closed-form RDP
estimate at the elastic-radius extraction surface (Sedan 1962
amplitude diagnostic; `scripts/pass7_amplitude_diagnostic.sh`).
Within the factor-5 envelope from the original pass-7 spec.
`solver_kind = RADIAL_LAGRANGIAN` is now the default for
`[NEAR_FIELD_SOURCE] mode = DYNAMIC_PLASTIC`. The strict pass-7
spec tolerances on the standalone gates (5 percent peak amplitude,
1 percent reflected energy, 10 percent Sedov prefactor, 20 percent
NTS cavity, 2 percent energy conservation, documented convergence
order) remain pass-8 follow-up: they require either a tabulated
plasma-regime EOS (ANEOS / SESAME / QEOS) replacing extrapolated
Tillotson, an explicit Marshak-wave radiation-transport phase
replacing the Z-R end-state approximation, or a higher-order
numerical scheme replacing the explicit Wilkins-AV finite-volume
update.

References for this section: Tillotson 1962, Melosh 1989 (Impact
Cratering: A Geologic Process), Zel'dovich and Raizer 1967 (Physics
of Shock Waves and High-Temperature Hydrodynamic Phenomena, vol II),
Wilkins 1980 (Computer Simulation of Dynamic Phenomena), Marsh 1980
(LASL Shock Hugoniot Data), Trunin 2001 (RFNC-VNIIEF shock data),
Carter 1979 (LA-7873 NaCl shock data).

## Pass-10 axis-1a closeout

Pass-10 promotes the three HIGHEST-tier scaffolds on the axis-1a
fidelity ladders to working implementations. After this pass every
cell on the radiation-phase, cavity-EOS, opacity, operator-splitting,
and time-integrator ladders has a working implementation; LOW / MED /
HIGH defaults reproduce the corresponding pass-N behaviour
byte-for-byte under pinned configs.

### Multigroup radiation transport

`radiation_phase = MARSHAK_MULTIGROUP` engages the new
`MultigroupRadiationDiffusionSolver`. The grey diffusion equation
generalises to G coupled equations indexed by frequency group g:

```
dE_r^g/dt = (1/r^2) d/dr [ r^2 (c / (3 kappa_R^g rho)) dE_r^g/dr ]
          + c kappa_P^g rho ( 4 pi B_g(T_m) / c - E_r^g )
```

Per global timestep we run a Newton outer iteration: refresh
per-cell, per-group opacities at the current `T_iter`, solve G
separate tridiagonal systems for `E_g^{n+1}`, update `T_m` from the
linearised matter-energy balance summed over all groups, repeat
until the relative residual on `E_r^g` drops below
`radiation_newton_tolerance`. The matter-temperature equation
linearises around `T_iter` using `dB_g/dT` (computed by finite
difference from the band-integrated Planck integrals).

Per-group opacity follows path A from the pass-10 spec: the analytic
Mihalas-Mihalas 1984 sec 82.2 smoothed-continuum bound-bound + Kramers'
free-free + Thomson scattering model, evaluated at the cell's
`(rho, T_m)` per timestep. Per-group means are computed by Simpson
quadrature in log-frequency space over the group's
[`nu_g`, `nu_{g+1}`] band. `MultigroupOpacityEvaluator` ships the
quadrature; `FrequencyGroupGrid` parameterises G, `nu_min`, `nu_max`,
and the number of Simpson sub-points per group.

The default group grid is 16 log-spaced groups from 1e14 Hz (~6
micron IR) to 1e18 Hz (~3 nm soft X-ray), with 17 Simpson sub-points
per group. This covers the rock-plasma emission spectrum across the
cavity-formation regime (`T_m` = 1e4 to 1e7 K). Configurable via
the new `[NEAR_FIELD_SOURCE]` sub-keys `radiation_n_groups` (1-256),
`radiation_freq_min_hz`, `radiation_freq_max_hz`,
`radiation_simpson_points`.

### TABULATED_FULL EOS

`cavity_eos = TABULATED_FULL` removes the pass-9 sin^2 patch window.
Every EOS query goes through the tabulated reader with a Tillotson
safety net for out-of-table coverage. The pass-9 entry's residual
(small numerical artifact at the patch transition) is gone in this
mode; the residual factor-3 envelope on Salmon CavityRadius after
pass-10 is the spherical-symmetry assumption itself, named axis-1b
for follow-up.

### TABULATED_FULL opacity and the dispatch matrix

`opacity_model = TABULATED_FULL` (grey) removes the Z-R fallback at
in-table queries. The pass-10 dispatch matrix
`(opacity_model, radiation_phase)` is documented and enforced at
config time:

| opacity_model       | radiation_phase     | behaviour                  |
|---------------------|---------------------|----------------------------|
| `TABULATED_PATCHED` | `MARSHAK_GREY`      | pass-9 patched 2D table with Z-R fallback |
| `TABULATED_FULL`    | `MARSHAK_GREY`      | pass-9 2D table only, no fallback                 |
| `TABULATED_FULL`    | `MARSHAK_MULTIGROUP`| per-group analytic model (path A)                  |
| `TABULATED_PATCHED` | `MARSHAK_MULTIGROUP`| not allowed; throws clear error directing at FULL  |

### Higher-order explicit time integrators

`time_integrator = EXPLICIT_EULER` (default) preserves the pass-9
byte-identical hydro substep. `time_integrator = TVD_RK2` (Heun's
method) and `time_integrator = RK3_SSP` (Shu-Osher 1988
strong-stability-preserving Runge-Kutta) are implemented as convex
blends of the same explicit-Euler hydro operator:

  TVD_RK2:
    y1     = y_n + dt L(y_n)
    y1'    = y1 + dt L(y1)
    y_n+1  = (1/2) y_n + (1/2) y1'

  RK3_SSP:
    y1     = y_n + dt L(y_n)
    y1'    = y1 + dt L(y1)
    y2     = (3/4) y_n + (1/4) y1'
    y2'    = y2 + dt L(y2)
    y_n+1  = (1/3) y_n + (2/3) y2'

The SSP property of these schemes is preserved provided the
underlying explicit-Euler operator is monotone under the per-step
CFL bound, which the pass-6 Wilkins-AV Lagrangian update is.

### Strang inner-substep convergence diagnostic

Setting `operator_splitting_convergence_diagnostic = true` exposes
the inner-CFL substep dt at the first call of each step via the
new `RadialLagrangianSolver::getDiagnosticInnerSubstepDt()` /
`getDiagnosticInnerSubstepCount()` accessors. The
`Physics.Marshak.OperatorSplittingConvergence` test is then run with
this diagnostic enabled to expose the second-order Strang behaviour
even at CI-achievable resolution.

### Pass-10 references

- Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
  Hydrodynamics", Oxford University Press, sec 80 (multigroup
  diffusion), sec 82.2 (line opacity smoothed continuum, gaunt factor).
- Pomraning, G. C. (1973), "The Equations of Radiation
  Hydrodynamics", Pergamon Press, ch IV (multigroup formulation).
- Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
  Waves and High-Temperature Hydrodynamic Phenomena", vol I,
  ch V (free-free Kramers' opacity, frequency dependence).
- Shu, C.-W. and Osher, S. (1988), "Efficient implementation of
  essentially non-oscillatory shock-capturing schemes",
  J. Comp. Phys 77, pp 439-471 (SSP3 / SSP2 RK schemes).
- Strang, G. (1968), "On the construction and comparison of
  difference schemes", SIAM J. Num. Anal. 5(3), pp 506-517.

## Pass-9 tabulated EOS / opacity and Strang splitting

Pass-9 (axis 1a) advanced two rungs on the cavity-EOS and opacity
ladders and added second-order time coupling. Ships:

- `cavity_eos = TILLOTSON_TABULATED_PATCH`: Tillotson blended with a
  Z-R partial-ionization plasma table (sin^2 in pressure across
  [5e10, 6e10] Pa) sourced from
  `tools/tabulated_data/tables/eos/<medium>_aneos.h5`.
- `opacity_model = TABULATED_PATCHED`: Z-R blended with Mihalas-
  Mihalas-corrected Rosseland/Planck means (sin^2 in temperature
  across [1e5, 1.26e5] K) sourced from
  `tools/tabulated_data/tables/opacity/<medium>_{rosseland,planck}.h5`.
- `operator_splitting = STRANG`: second-order Strang splitting between
  hydro and Marshak radiation. `LIE` remains the pass-8 byte-identical
  default.
- `radial_outer_radius_m`: direct outer-domain override (used in pass-9
  to extend Salmon to 700 m for the free-field gate triage).

Tabulated readers under `include/io/TabulatedData/TabulatedDataReader.hpp`
require an HDF5 file with `source_citation` metadata; out-of-table
queries fall back to the analytic with one-time stderr warnings.

Pass-9 closed the Salmon FreeFieldPeakVelocity_166m and _322m gates at
spec target (factor 2). Three residuals remained at end of pass-9:
Salmon CavityRadius factor 3 (target 5 %), Marshak self-similar front
factor 8 (target factor 2), Pokhran mb +/-0.4 (target +/-0.3).

### Pass-9 references

- Marsh, S. P. (1980), "LASL Shock Hugoniot Data", University of
  California Press (granite Hugoniot reference).
- Strang, G. (1968), "On the construction and comparison of difference
  schemes", SIAM J. Num. Anal. 5(3), pp 506-517.

## Pass-11 implicit diffusion, sponge BC, ANEOS extension

Pass-11 (axis 1c implicit diffusion + 549 m BC + ANEOS extension)
closed the Marshak self-similar gate to factor 2 and the radiation-
energy conservation gate to 2 % via the new `DiffusionTimeIntegrator`
strategy in `MarshakRadiationDiffusion.hpp`:

- `time_integrator_diffusion = BACKWARD_EULER` (LOW/MED, default).
- `time_integrator_diffusion = CRANK_NICOLSON` (HIGH, second order,
  A-stable).
- `time_integrator_diffusion = BDF2` (HIGHEST, second order, A-stable).
  Pass-11 default for HIGH+ tiers.

The 549 m free-field gate triage (sweep at radial_outer_radius_m =
[700, 1000, 1500, 2000] m) shows the impedance BC at 700 m contaminates
the gauge by factor ~3.4 at the closest two radii; pass-11 ships the
Israeli-Orszag 1981 graded-damping sponge layer
(`sponge_layer_enabled = true`) as the BC fix. The default remains the
characteristic BC for byte-identical pass-9/10 reproducibility.

Extended ANEOS Hugoniot match gates ship at 30 % envelope (Marsh 1980
granite, McQueen 1970 salt). The 5 % spec target requires a Tillotson
parameter refit, named explicitly as axis-4 in
[HISTORIC_NUCLEAR_ROADMAP.md](HISTORIC_NUCLEAR_ROADMAP.md).

The axis-1b 3D source ball scaffolding lands in pass-11:
`Source3DBall.hpp` interface, `cavity_geometry` config knob
(`THREE_DIMENSIONAL` throws until pass-13), design stub at
[AXIS_1B_DESIGN.md](AXIS_1B_DESIGN.md). Pass-13 implements the 3D
source ball on this scaffold.

### Pass-11 references

- Israeli, M. and Orszag, S. A. (1981), "Approximation of radiation
  boundary conditions", J. Comp. Phys. 41, pp 115-135 (graded-damping
  sponge layer).
- Hairer, E. and Wanner, G. (1996), "Solving Ordinary Differential
  Equations II: Stiff and Differential-Algebraic Problems", Springer
  (BDF2, Crank-Nicolson stability).
- McQueen, R. G., et al. (1970), "The Equation of State of Solids
  from Shock Wave Studies", in High-Velocity Impact Phenomena (salt
  Hugoniot reference).

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

- [Wave Propagation](../tutorials/09_WAVE_PROPAGATION.md) - Seismic wave physics
- [Physics Models](./PHYSICS_MODELS.md) - Material failure and constitutive models
- [Configuration Reference](./CONFIGURATION.md) - Full config options
- [Benchmarks](./BENCHMARKS.md) - Validation suite
