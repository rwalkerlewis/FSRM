# Implicit-Explicit Time Integration Transition

## Overview

The Implicit-Explicit (IMEX) Time Integration Transition feature enables automatic switching between quasi-static (implicit) and dynamic (explicit) time integration schemes during a simulation. This is essential for modeling induced seismicity where:

1. **Slow processes** (pore pressure diffusion, fault loading) occur over months-years
2. **Fast processes** (dynamic rupture, seismic waves) occur over milliseconds-seconds

Without IMEX, you would need to either:
- Use tiny time steps throughout (computationally prohibitive)
- Miss the dynamic physics entirely (physically incomplete)

## Physical Motivation

### Induced Seismicity Scenario

```
Time Scale Analysis:
├── Injection period: months-years
├── Pore pressure diffusion: days-months  
│   → Use IMPLICIT with dt ~ hours-days
│
├── Fault loading: weeks-months
│   → Use IMPLICIT with dt ~ hours
│
├── Dynamic rupture: milliseconds
│   → Use EXPLICIT with dt ~ 0.1 ms (CFL limited)
│
├── Wave propagation: seconds
│   → Use EXPLICIT with dt ~ 0.1 ms
│
└── Post-seismic relaxation: hours-days
    → Return to IMPLICIT with dt ~ hours
```

### Governing Equations

**Implicit (Quasi-static Poroelasticity):**
```
∇·σ' - α∇p = 0                    (Equilibrium)
∂/∂t(α∇·u + p/M) + ∇·q = f        (Mass conservation)
q = -(k/μ)∇p                       (Darcy's law)
```

**Explicit (Dynamic Poroelasticity - Biot's equations):**
```
ρ ü + C·u̇ = ∇·σ' - α∇p + ρg       (Momentum with inertia)
∂/∂t(α∇·u + p/M) + ∇·q = f        (Mass conservation)
```

The key difference is the inertia term `ρ ü` which enables wave propagation.

## Configuration

### Basic Setup

Add an `[IMEX]` section to your configuration file:

```ini
[IMEX]
enabled = true
initial_mode = IMPLICIT

# Implicit mode (quasi-static)
implicit_dt_initial = 3600.0      # 1 hour
implicit_dt_max = 86400.0         # 1 day
implicit_method = BEULER          # Backward Euler

# Explicit mode (dynamic)
explicit_dt_initial = 1e-4        # 0.1 ms
explicit_dt_max = 1e-3            # 1 ms
cfl_factor = 0.5                  # CFL safety
explicit_method = EULER

# Trigger condition (when to switch to explicit)
trigger_type = COULOMB_FAILURE
coulomb_threshold = 0.0

# Settling condition (when to return to implicit)
settling_type = COMBINED
min_dynamic_duration = 0.5
max_dynamic_duration = 30.0
settling_velocity = 1e-6
settling_energy_ratio = 0.01
```

### Trigger Types

Controls when to switch from IMPLICIT to EXPLICIT mode:

| Type | Description | Threshold Parameter |
|------|-------------|---------------------|
| `STRESS_THRESHOLD` | Von Mises stress exceeds threshold | `stress_threshold` (Pa) |
| `COULOMB_FAILURE` | Coulomb Failure Function ≥ threshold | `coulomb_threshold` (Pa) |
| `SLIP_RATE` | Fault slip rate exceeds seismic threshold | `slip_rate_threshold` (m/s) |
| `VELOCITY` | Particle velocity exceeds threshold | `velocity_threshold` (m/s) |
| `ACCELERATION` | Acceleration exceeds threshold | `acceleration_threshold` (m/s²) |
| `ENERGY_RATE` | Rate of kinetic energy change | `energy_rate_threshold` (J/s) |

### Settling Types

Controls when to return from EXPLICIT to IMPLICIT mode:

| Type | Description | Parameters |
|------|-------------|------------|
| `TIME_ELAPSED` | Fixed duration in explicit mode | `min_dynamic_duration` |
| `VELOCITY_DECAY` | Velocity drops below threshold | `settling_velocity` |
| `ENERGY_DECAY` | Kinetic energy drops relative to peak | `settling_energy_ratio` |
| `SLIP_RATE_DECAY` | Slip rate returns to interseismic | `settling_slip_rate` |
| `COMBINED` | All criteria must be met | All above parameters |

## Usage

### From Configuration File

```bash
# Run with the IMEX-enabled config
mpirun -np 4 ./simulator -c config/induced_seismicity_imex.config
```

### Programmatic Setup

```cpp
#include "ImplicitExplicitTransition.hpp"

// In simulator setup:
ConfigReader::IMEXConfig imex_config;
reader.parseIMEXConfig(imex_config);

// Setup IMEX manager
simulator.setupIMEXTransition(imex_config);
```

## Output

The IMEX manager generates several outputs:

### Transition Log

File: `output/imex_transitions.log`

```
# Time, From Mode, To Mode, Reason
15778800.0 IMPLICIT -> EXPLICIT : Trigger condition met
15778815.3 EXPLICIT -> IMPLICIT : Settling condition met
31557600.0 IMPLICIT -> EXPLICIT : Trigger condition met
...
```

### Phase Statistics

For each phase (implicit or explicit segment):
- Start/end time
- Number of steps
- Time step range (min, max, avg)
- Peak velocity
- Peak slip rate
- Total slip accumulated
- Energy released

### Dynamic Events

For each seismic event captured in explicit mode:
- Trigger time
- Duration  
- Peak moment rate
- Total moment
- Magnitude (Mw)
- Peak slip rate
- Maximum slip

## Performance Considerations

### Time Step Selection

**Implicit mode:**
- Large time steps (hours-days) are stable
- Governed by accuracy, not stability
- Adaptive time stepping based on convergence

**Explicit mode:**
- CFL condition limits time step:
  ```
  dt ≤ CFL * h_min / v_max
  ```
  where `h_min` is minimum cell size, `v_max` is maximum wave speed

### Computational Cost

Typical cost distribution for a 2-year induced seismicity simulation with 10 events:

| Phase | Time Steps | Wall Time |
|-------|------------|-----------|
| Implicit (98% of simulation time) | ~17,000 | ~80% |
| Explicit (2% of simulation time) | ~300,000 | ~20% |

The IMEX approach is ~100x faster than pure explicit time stepping.

### Memory Usage

Explicit mode requires additional storage for:
- Velocity vectors (3 components)
- Acceleration vectors (3 components)
- Lumped mass matrix (diagonal)

Approximately 2x memory compared to implicit-only.

## Physical Interpretation

### Event Detection

An event is detected when:
1. Trigger condition is met (e.g., CFF ≥ 0)
2. Slip rate exceeds seismic threshold (~10⁻³ m/s)

### Event Properties

During explicit phase, the manager tracks:
- Slip rate history at each fault node
- Moment rate: Ṁ₀ = μ A V (shear modulus × area × slip rate)
- Total moment: M₀ = ∫ Ṁ₀ dt
- Magnitude: Mw = (2/3) log₁₀(M₀) - 6.07

### Wave Propagation

In explicit mode, the simulation captures:
- P-waves (compressional): v_p ≈ 5-6 km/s
- S-waves (shear): v_s ≈ 3-4 km/s
- Biot slow wave (fluid diffusion): v_slow ≈ 100 m/s

## Validation

The IMEX transition has been validated against:

1. **Analytical solutions** for 1D wave propagation
2. **SCEC benchmarks** (TPV5, TPV10) for dynamic rupture
3. **Laboratory experiments** for rate-state friction

## Troubleshooting

### Common Issues

**Premature triggering:**
- Increase `stress_threshold` or `coulomb_threshold`
- Check initial stress state is not too close to failure

**Slow settling:**
- Decrease `settling_velocity` threshold
- Increase artificial damping (`damping_alpha`)
- Check absorbing boundary conditions

**Numerical instability in explicit mode:**
- Reduce `cfl_factor` (try 0.3 or 0.25)
- Check for very small cells in mesh
- Verify material properties (wave speeds)

**Memory issues:**
- Use lumped mass (`use_lumped_mass = true`)
- Reduce output frequency in explicit mode

## References

1. Biot, M.A. (1962). Mechanics of deformation and acoustic propagation in porous media.
2. Day, S.M. et al. (2005). Comparison of finite difference and boundary integral solutions to three-dimensional spontaneous rupture.
3. Dieterich, J.H. (1979). Modeling of rock friction.
4. Segall, P. (2010). Earthquake and Volcano Deformation.
