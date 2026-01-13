# Tutorial 09: Wave Propagation and Dynamics

Model elastic and poroelastic wave propagation, seismic sources, and dynamic rupture.

## Table of Contents

1. [Wave Physics Overview](#wave-physics-overview)
2. [Elastodynamics](#elastodynamics)
3. [Poroelastodynamics](#poroelastodynamics)
4. [Seismic Sources](#seismic-sources)
5. [Absorbing Boundary Conditions](#absorbing-boundary-conditions)
6. [Attenuation and Damping](#attenuation-and-damping)
7. [Dynamic Rupture](#dynamic-rupture)
8. [SCEC Benchmarks](#scec-benchmarks)

---

## Wave Physics Overview

FSRM supports multiple wave propagation models:

```
┌──────────────────────────────────────────────────────────────┐
│                 Wave Propagation Physics                     │
├──────────────────────────────────────────────────────────────┤
│                                                               │
│  ELASTODYNAMICS (Elastic waves)                              │
│  ├─ P-waves (compressional)                                 │
│  ├─ S-waves (shear)                                         │
│  └─ Surface waves (Rayleigh, Love)                          │
│                                                               │
│  POROELASTODYNAMICS (Biot's theory)                          │
│  ├─ Fast P-wave (in-phase solid-fluid)                      │
│  ├─ Slow P-wave (out-of-phase, Biot wave II)               │
│  ├─ S-wave (shear)                                          │
│  └─ Fluid-solid coupling                                    │
│                                                               │
│  ACOUSTICS (Simplified)                                      │
│  └─ Pressure waves only                                      │
│                                                               │
│  SEISMIC SOURCES                                             │
│  ├─ Moment tensor                                           │
│  ├─ Double couple (earthquake)                              │
│  ├─ Point force                                             │
│  └─ Explosive (isotropic)                                   │
│                                                               │
└──────────────────────────────────────────────────────────────┘
```

### Enabling Dynamics

```ini
[SIMULATION]
enable_elastodynamics = true       # Elastic waves only
# OR
enable_poroelastodynamics = true   # Coupled waves (Biot)

# For dynamic time stepping
adaptive_dt = true
cfl_based_dt = true                # Use CFL condition
```

---

## Elastodynamics

### Basic Configuration

```ini
[SIMULATION]
solid_model = ELASTIC
enable_elastodynamics = true
end_time = 10.0                    # 10 seconds of wave propagation

[DYNAMICS]
enable = true
wave_type = ELASTIC                # ELASTIC, POROELASTIC, ACOUSTIC

# Time stepping
cfl_number = 0.5                   # CFL < 1 for stability
dt_initial = 1.0e-4                # Initial timestep

# Grid resolution: at least 10 points per wavelength
# For 10 Hz waves in 3000 m/s medium: λ = 300m, need dx ≤ 30m
```

### Material Properties for Waves

```ini
[ROCK]
# Elastic properties
youngs_modulus = 20.0e9            # Pa
poisson_ratio = 0.25
density = 2500.0                   # kg/m³

# Derived wave velocities:
# Vp = √((K + 4G/3) / ρ) ≈ 3464 m/s
# Vs = √(G / ρ) ≈ 2000 m/s

# Or specify directly:
# p_wave_velocity = 3464.0         # m/s
# s_wave_velocity = 2000.0         # m/s
```

### Wave Velocity Reference

| Material | Vp (m/s) | Vs (m/s) | ρ (kg/m³) |
|----------|----------|----------|-----------|
| Sandstone | 2500-4500 | 1200-2500 | 2200-2700 |
| Limestone | 4000-6500 | 2500-3500 | 2400-2700 |
| Shale | 2000-4500 | 1000-2500 | 2000-2700 |
| Granite | 5000-6500 | 2800-3500 | 2600-2700 |
| Salt | 4500-5500 | 2500-3100 | 2100-2200 |

---

## Poroelastodynamics

### Biot's Wave Equations

```ini
[SIMULATION]
solid_model = POROELASTIC
enable_poroelastodynamics = true

[DYNAMICS]
wave_type = POROELASTIC

[ROCK]
constitutive_model = POROELASTIC

# Solid properties
youngs_modulus = 20.0e9
poisson_ratio = 0.25
density = 2650.0                   # Grain density
porosity = 0.20

# Poroelastic parameters
biot_coefficient = 0.8
biot_modulus = 1.0e10

# Permeability (affects slow wave)
permeability_x = 100.0             # mD

[FLUID]
density = 1000.0
viscosity = 0.001
bulk_modulus = 2.2e9               # Water ~ 2.2 GPa
```

### Biot Wave Velocities

Three wave types in poroelastic media:

1. **Fast P-wave (Biot wave I)**: Similar to elastic P-wave
2. **Slow P-wave (Biot wave II)**: Out-of-phase solid-fluid motion
3. **S-wave**: Shear wave (similar to elastic)

```ini
# FSRM computes these automatically, but you can specify:
[DYNAMICS]
# vp_fast = 3500.0               # Fast P-wave (m/s)
# vp_slow = 1000.0               # Slow P-wave (m/s)
# vs = 2000.0                    # S-wave (m/s)
```

### Dynamic Permeability Enhancement

Wave passage can temporarily enhance permeability:

```ini
[DYNAMICS]
dynamic_permeability = true

# Model: k/k0 = 1 + α × |ε_vol| + β × |σ'_eff|
permeability_strain_coeff = 100.0  # α
permeability_stress_coeff = 1.0e-7 # β

# Recovery after wave passage
permeability_recovery_time = 10.0  # seconds
```

---

## Seismic Sources

### Moment Tensor Source

General seismic source representation:

```ini
[SOURCE1]
type = MOMENT_TENSOR

# Location
location = 5000.0, 5000.0, 5000.0  # x, y, z (m)

# Moment tensor components (N·m)
Mxx = 1.0e18
Myy = -1.0e18
Mzz = 0.0
Mxy = 0.0
Mxz = 0.0
Myz = 0.0

# Time function
time_function = RICKER             # RICKER, GAUSSIAN, STEP, CUSTOM
peak_frequency = 5.0               # Hz
amplitude = 1.0
onset_time = 0.1                   # Start time (s)
```

### Double Couple (Earthquake)

Shear faulting source:

```ini
[SOURCE1]
type = DOUBLE_COUPLE

location = 5000.0, 5000.0, 5000.0

# Fault orientation
strike = 0.0                       # degrees from N
dip = 90.0                         # degrees from horizontal
rake = 0.0                         # 0 = left-lateral, 180 = right-lateral

# Magnitude
seismic_moment = 1.0e18            # N·m
# OR specify moment magnitude
# magnitude = 5.5

# Time function
time_function = RICKER
peak_frequency = 2.0
```

### Point Force

Single force application:

```ini
[SOURCE1]
type = POINT_FORCE

location = 5000.0, 5000.0, 0.0     # Surface source

# Force direction and magnitude
force_x = 0.0
force_y = 0.0
force_z = -1.0e12                  # Downward force (N)

time_function = GAUSSIAN
peak_frequency = 10.0
```

### Explosive Source

Isotropic expansion:

```ini
[SOURCE1]
type = EXPLOSIVE

location = 5000.0, 5000.0, 100.0

# Isotropic moment
seismic_moment = 1.0e15            # N·m

time_function = RICKER
peak_frequency = 20.0
```

### Source Time Functions

```ini
[SOURCE1]
# RICKER wavelet (common for seismic)
time_function = RICKER
peak_frequency = 5.0               # f_p
onset_time = 0.0
# f(t) = (1 - 2π²f_p²(t-t_s)²) × exp(-π²f_p²(t-t_s)²)

# GAUSSIAN derivative
time_function = GAUSSIAN
peak_frequency = 5.0
# First derivative of Gaussian

# STEP function
time_function = STEP
step_time = 0.1                    # Time of step
rise_time = 0.01                   # Ramp duration

# CUSTOM from file
time_function = CUSTOM
time_function_file = source_time.csv
# Format: time, amplitude
```

### Multiple Sources

```ini
[SOURCE1]
type = DOUBLE_COUPLE
location = 5000.0, 3000.0, 5000.0
strike = 0.0
dip = 90.0
rake = 0.0
seismic_moment = 5.0e17
onset_time = 0.0

[SOURCE2]
type = DOUBLE_COUPLE
location = 5000.0, 5000.0, 5000.0
strike = 0.0
dip = 90.0
rake = 0.0
seismic_moment = 1.0e18
onset_time = 0.5                   # Delayed start

[SOURCE3]
type = DOUBLE_COUPLE
location = 5000.0, 7000.0, 5000.0
strike = 0.0
dip = 90.0
rake = 0.0
seismic_moment = 3.0e17
onset_time = 1.2
```

---

## Absorbing Boundary Conditions

Prevent artificial reflections at domain boundaries:

### Perfectly Matched Layer (PML)

```ini
[DYNAMICS]
# PML absorbing boundaries
absorbing_bc = PML

absorbing_layers = 15              # Number of PML cells
absorbing_strength = 100.0         # Damping coefficient

# PML on specific boundaries
pml_xmin = true
pml_xmax = true
pml_ymin = true
pml_ymax = true
pml_zmin = false                   # Free surface
pml_zmax = true
```

### Stacey Absorbing BC

First-order absorbing condition:

```ini
[DYNAMICS]
absorbing_bc = STACEY

# Applied to boundaries
stacey_xmin = true
stacey_xmax = true
# ...
```

### Clayton-Engquist

Higher-order absorbing condition:

```ini
[DYNAMICS]
absorbing_bc = CLAYTON_ENGQUIST
ce_order = 2                       # Order of approximation
```

### Free Surface

For surface wave modeling:

```ini
[DYNAMICS]
# Free surface at top (ZMAX)
free_surface = ZMAX

# Traction-free condition: σ·n = 0
```

---

## Attenuation and Damping

### Rayleigh Damping

Proportional damping:

```ini
[DYNAMICS]
damping_type = RAYLEIGH

# C = α × M + β × K
damping_alpha = 0.1                # Mass proportional
damping_beta = 0.001               # Stiffness proportional

# Or specify for target frequency range:
# damping_f1 = 1.0                 # Lower frequency (Hz)
# damping_f2 = 10.0                # Upper frequency (Hz)
# damping_xi = 0.05                # Damping ratio (5%)
```

### Quality Factor (Q)

Intrinsic attenuation:

```ini
[DYNAMICS]
attenuation = true

# Constant Q model
quality_factor_p = 100.0           # Q_p for P-waves
quality_factor_s = 50.0            # Q_s for S-waves

# Frequency-dependent Q
# q_model = POWER_LAW
# q_exponent = 0.3                 # Q ∝ f^0.3
```

### Viscoelastic Attenuation

Using generalized Maxwell/SLS model:

```ini
[ROCK]
constitutive_model = VISCOELASTIC_SLS

# For wave simulation with attenuation
num_mechanisms = 3                 # Number of relaxation mechanisms

# Relaxation times (cover frequency band of interest)
relaxation_times = 0.001, 0.01, 0.1  # seconds

# Modulus ratios
modulus_ratios = 0.1, 0.2, 0.1     # Fraction of unrelaxed modulus
```

---

## Dynamic Rupture

Earthquake simulation with fault slip:

### Configuration

```ini
[SIMULATION]
solid_model = ELASTIC
enable_elastodynamics = true
enable_faults = true

[FAULT1]
name = RUPTURE_PLANE
x = 5000.0
y = 5000.0
z = 7500.0
strike = 0.0
dip = 90.0
length = 20000.0
width = 15000.0

# Friction for dynamic rupture
friction_law = SLIP_WEAKENING
static_friction = 0.677
dynamic_friction = 0.525
critical_slip = 0.4                # Dc (m)

use_split_nodes = true

# Nucleation
nucleation_enabled = true
nucleation_x = 0.0                 # Along-strike from center
nucleation_z = 7500.0              # Depth
nucleation_radius = 3000.0         # Nucleation patch radius
nucleation_stress = 81.6e6         # Overstress to initiate

# Initial stress
[IC_STRESS]
stress_type = TENSOR
sxx = -120.0e6
syy = 0.0
szz = -120.0e6
sxy = 70.0e6                       # Initial shear on fault

[DYNAMICS]
cfl_number = 0.3                   # Lower CFL for stability
```

### Output for Rupture

```ini
[OUTPUT]
frequency = 10
displacement = true
velocity = true
fault_slip = true
fault_slip_rate = true
fault_traction = true
```

---

## SCEC Benchmarks

FSRM includes implementations of SCEC TPV benchmarks:

### TPV3: Vertical Strike-Slip

```bash
# Run TPV3 benchmark
./fsrm -c config/scec_tpv3.config

# Compare with reference solutions
python scripts/compare_tpv3.py output/
```

### TPV5: Through-Going Rupture

```ini
# config/scec_tpv5.config excerpt
[FAULT1]
name = TPV5_FAULT
x = 0.0
y = 0.0
z = 7500.0
strike = 90.0
dip = 90.0
length = 30000.0
width = 15000.0

friction_law = SLIP_WEAKENING
static_friction = 0.677
dynamic_friction = 0.373           # TPV5 uses lower μ_d
critical_slip = 0.4
```

### LOH.1: Layer Over Halfspace

Point source in layered medium:

```ini
# config/scec_loh1.config excerpt
[LAYERS]
num_layers = 2

# Layer 1: soft top layer
layer_1_thickness = 1000.0
layer_1_vp = 4000.0
layer_1_vs = 2000.0
layer_1_density = 2600.0

# Layer 2: halfspace
layer_2_vp = 6000.0
layer_2_vs = 3464.0
layer_2_density = 2700.0

[SOURCE1]
type = DOUBLE_COUPLE
location = 0.0, 0.0, 2000.0
strike = 0.0
dip = 90.0
rake = 0.0
seismic_moment = 1.0e18
```

### Running Benchmarks

```bash
# TPV benchmarks (dynamic rupture)
./fsrm -c config/scec_tpv3.config
./fsrm -c config/scec_tpv5.config
./fsrm -c config/scec_tpv10.config  # Normal fault

# LOH benchmarks (wave propagation)
./fsrm -c config/scec_loh1.config
./fsrm -c config/scec_loh2.config   # Anelastic
./fsrm -c config/scec_loh3.config   # Multiple sources

# Verification
python scripts/verify_scec.py --benchmark tpv3 --output output/
```

---

## Practical Examples

### Earthquake Wave Propagation

```ini
[SIMULATION]
name = earthquake_waves
solid_model = ELASTIC
enable_elastodynamics = true
end_time = 60.0                    # 60 seconds

[GRID]
nx = 200
ny = 200
nz = 100
Lx = 100000.0                      # 100 km
Ly = 100000.0
Lz = 50000.0                       # 50 km depth

[ROCK]
youngs_modulus = 80.0e9            # ~3500 m/s Vp
poisson_ratio = 0.25
density = 2700.0

[SOURCE1]
type = DOUBLE_COUPLE
location = 50000.0, 50000.0, 10000.0  # 10 km depth
strike = 45.0
dip = 80.0
rake = -90.0                       # Normal fault
magnitude = 6.0
time_function = RICKER
peak_frequency = 1.0

[DYNAMICS]
cfl_number = 0.5
absorbing_bc = PML
absorbing_layers = 20
free_surface = ZMAX

# Attenuation
attenuation = true
quality_factor_p = 200.0
quality_factor_s = 100.0

[OUTPUT]
frequency = 50
displacement = true
velocity = true
```

### Reservoir Monitoring (Induced Events)

```ini
[SIMULATION]
name = microseismic_monitoring
solid_model = ELASTIC
enable_elastodynamics = true
end_time = 2.0                     # 2 seconds for microseismic

[GRID]
# High resolution for small events
nx = 100
ny = 100
nz = 50
Lx = 2000.0                        # 2 km
Ly = 2000.0
Lz = 1000.0

[ROCK]
youngs_modulus = 30.0e9
poisson_ratio = 0.25
density = 2500.0

# Small injection-induced event
[SOURCE1]
type = DOUBLE_COUPLE
location = 1000.0, 1000.0, 500.0
strike = 30.0
dip = 70.0
rake = 10.0
magnitude = 1.5                    # Microseismic event
time_function = RICKER
peak_frequency = 50.0              # Higher frequency for small events

[DYNAMICS]
cfl_number = 0.5
absorbing_bc = PML
absorbing_layers = 10

[OUTPUT]
# High-frequency output for small events
frequency = 1
displacement = true
velocity = true

# Receiver locations (for seismogram output)
[RECEIVERS]
num_receivers = 5
receiver_1 = 500.0, 1000.0, 0.0    # Surface receiver
receiver_2 = 1500.0, 1000.0, 0.0
receiver_3 = 1000.0, 500.0, 0.0
receiver_4 = 1000.0, 1500.0, 0.0
receiver_5 = 1000.0, 1000.0, 100.0 # Downhole receiver
```

---

## Output and Post-Processing

### Seismogram Output

```ini
[OUTPUT]
# Time series at specific points
seismogram_output = true
seismogram_format = CSV            # CSV, SAC, MSEED

[RECEIVERS]
# Define receiver locations
num_receivers = 10
receiver_file = receivers.csv      # x, y, z coordinates
```

### Synthetic Seismogram Analysis

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Load seismograms
data = np.loadtxt('output/receiver_001.csv', delimiter=',', skiprows=1)
time = data[:, 0]
vx = data[:, 1]
vy = data[:, 2]
vz = data[:, 3]

dt = time[1] - time[0]

fig, axes = plt.subplots(4, 1, figsize=(12, 12))

# Time series
axes[0].plot(time, vx * 1e6, 'r', label='Vx')
axes[0].plot(time, vy * 1e6, 'g', label='Vy')
axes[0].plot(time, vz * 1e6, 'b', label='Vz')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('Velocity (μm/s)')
axes[0].legend()
axes[0].set_title('Three-Component Seismogram')

# Spectrum
n = len(time)
freq = fftfreq(n, dt)[:n//2]
spectrum = np.abs(fft(vz))[:n//2]

axes[1].semilogy(freq, spectrum)
axes[1].set_xlabel('Frequency (Hz)')
axes[1].set_ylabel('Amplitude')
axes[1].set_title('Frequency Spectrum')
axes[1].set_xlim(0, 20)

# Spectrogram
from scipy.signal import spectrogram
f, t, Sxx = spectrogram(vz, 1/dt, nperseg=256)
axes[2].pcolormesh(t, f, 10*np.log10(Sxx), shading='auto')
axes[2].set_xlabel('Time (s)')
axes[2].set_ylabel('Frequency (Hz)')
axes[2].set_title('Spectrogram')
axes[2].set_ylim(0, 20)

# Particle motion
axes[3].plot(vx * 1e6, vz * 1e6)
axes[3].set_xlabel('Vx (μm/s)')
axes[3].set_ylabel('Vz (μm/s)')
axes[3].set_title('Particle Motion')
axes[3].axis('equal')

plt.tight_layout()
plt.savefig('seismogram_analysis.png', dpi=150)
```

---

**Previous**: [← Adaptive Mesh Refinement](08_ADAPTIVE_MESH_REFINEMENT.md) | **Next**: [Advanced Simulations →](10_ADVANCED_SIMULATIONS.md)
