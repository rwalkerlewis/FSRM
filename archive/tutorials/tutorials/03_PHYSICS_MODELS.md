# Tutorial 03: Physics Models

Understand and configure the physics models available in FSRM for flow, mechanics, thermal, and coupled simulations.

## Table of Contents

1. [Overview of Physics Models](#overview-of-physics-models)
2. [Flow Models](#flow-models)
3. [Solid Mechanics Models](#solid-mechanics-models)
4. [Thermal Models](#thermal-models)
5. [Coupled Physics](#coupled-physics)
6. [Physics Kernels Architecture](#physics-kernels-architecture)
7. [Selecting the Right Model](#selecting-the-right-model)

---

## Overview of Physics Models

FSRM provides a modular physics system where different models can be combined:

```
┌─────────────────────────────────────────────────────────────┐
│                    FSRM Physics System                       │
├─────────────────┬─────────────────┬─────────────────────────┤
│   FLOW          │   MECHANICS     │   THERMAL               │
├─────────────────┼─────────────────┼─────────────────────────┤
│ • Single-phase  │ • Linear elastic│ • Heat conduction       │
│ • Black oil     │ • Viscoelastic  │ • Convection            │
│ • Compositional │ • Poroelastic   │ • THM coupling          │
│ • CO2/Brine     │ • Elastoplastic │                         │
│                 │ • Anisotropic   │                         │
└─────────────────┴─────────────────┴─────────────────────────┘
                            │
                   ┌────────┴────────┐
                   │  COUPLING       │
                   │ • HM (Hydro-Mech)│
                   │ • TH (Thermo-Hydro)│
                   │ • THM (Full)     │
                   │ • Poroelastic   │
                   └─────────────────┘
```

### Enabling Physics

```ini
[SIMULATION]
# Select flow model
fluid_model = SINGLE_COMPONENT   # or BLACK_OIL, COMPOSITIONAL, etc.

# Select solid model  
solid_model = ELASTIC            # or POROELASTIC, VISCOELASTIC, etc.

# Enable coupling
enable_geomechanics = true       # Flow-mechanics coupling
enable_thermal = true            # Heat transfer
enable_fractures = true          # Fracture modeling
enable_faults = true             # Fault mechanics
```

---

## Flow Models

### Single-Phase Flow

The simplest flow model for single-component liquid or gas.

**Governing Equation:**
```
φ cₜ ∂P/∂t = ∇·(k/μ ∇P) + q
```

**Configuration:**

```ini
[SIMULATION]
fluid_model = SINGLE_COMPONENT

[FLUID]
type = SINGLE_PHASE
density = 1000.0           # kg/m³
viscosity = 0.001          # Pa·s (1 cP)
compressibility = 4.5e-10  # 1/Pa
reference_pressure = 1.0e5 # Pa

[ROCK]
porosity = 0.20
permeability_x = 100.0     # mD
compressibility = 1.0e-9   # Rock compressibility
```

**Use Cases:**
- Water injection/production
- Single-phase gas reservoirs
- Pressure transient analysis
- Aquifer modeling

### Black Oil Model

Three-phase model with dissolved gas in oil.

**Governing Equations:**
```
∂/∂t(φ ρₒ Sₒ) + ∇·(ρₒ vₒ) = qₒ
∂/∂t(φ ρw Sw) + ∇·(ρw vw) = qw
∂/∂t(φ(ρg Sg + Rₛ ρₒ Sₒ)) + ∇·(ρg vg + Rₛ ρₒ vₒ) = qg
```

**Configuration:**

```ini
[SIMULATION]
fluid_model = BLACK_OIL

[FLUID]
type = BLACK_OIL

# Standard condition densities
oil_density_std = 850.0    # kg/m³
gas_density_std = 0.9      # kg/m³  
water_density_std = 1000.0 # kg/m³

# Solution properties
solution_gor = 100.0       # sm³/sm³
bubble_point = 15.0e6      # Pa

# PVT correlation
pvt_correlation = STANDING # Options: STANDING, VASQUEZ_BEGGS, GLASO, AL_MARHOUN

# Relative permeability model
relperm_model = COREY      # COREY, BROOKS_COREY, LET, TABLE
Swc = 0.2                  # Connate water saturation
Sor = 0.2                  # Residual oil saturation
krw_max = 0.3              # Max water rel perm
kro_max = 1.0              # Max oil rel perm
nw = 2.0                   # Water Corey exponent
no = 2.0                   # Oil Corey exponent
```

**PVT Correlations:**

| Correlation | Best For |
|------------|----------|
| STANDING | Light oils, API > 30 |
| VASQUEZ_BEGGS | General purpose |
| GLASO | North Sea oils |
| AL_MARHOUN | Middle East oils |

### Compositional Model

Multi-component model with thermodynamic flash calculations.

**Configuration:**

```ini
[SIMULATION]
fluid_model = COMPOSITIONAL

[FLUID]
type = COMPOSITIONAL

# Number and names of components
num_components = 5
component_names = C1, C2, C3, C4-6, C7+

# Critical properties (one value per component)
critical_temperatures = 190.6, 305.4, 369.8, 460.0, 580.0  # K
critical_pressures = 4.6e6, 4.88e6, 4.25e6, 3.5e6, 2.5e6   # Pa
acentric_factors = 0.011, 0.099, 0.152, 0.22, 0.35
molar_weights = 16.04, 30.07, 44.1, 72.0, 140.0             # kg/kmol

# Binary interaction coefficients (optional)
# kij matrix as flattened array

# Equation of state
eos_type = PENG_ROBINSON   # PENG_ROBINSON, SRK, VAN_DER_WAALS, IDEAL

# Flash calculation settings
flash_tolerance = 1.0e-8
flash_max_iterations = 100
```

**Initial Composition:**

```ini
[IC1]
field = COMPOSITION
component = C1
distribution = UNIFORM
value = 0.50               # 50% methane

[IC2]
field = COMPOSITION
component = C7+
distribution = UNIFORM
value = 0.20               # 20% heavy fraction
```

### CO2/Brine System

For CO2 sequestration studies:

```ini
[SIMULATION]
fluid_model = CO2

[FLUID]
type = CO2

# CO2 property model
co2_model = SPAN_WAGNER    # High-accuracy EOS

# Brine properties
salinity = 100000          # ppm
salt_type = NaCl

# CO2 dissolution
enable_dissolution = true
henry_coefficient = 1.5e-5 # mol/(Pa·m³)

# Mutual solubility
enable_water_in_co2 = true
```

---

## Solid Mechanics Models

### Linear Elastic

Standard Hookean elasticity:

```ini
[ROCK]
constitutive_model = LINEAR_ELASTIC

youngs_modulus = 20.0e9    # Pa
poisson_ratio = 0.25
density = 2650.0           # kg/m³
```

**Stress-Strain Relation:**
```
σᵢⱼ = λ εₖₖ δᵢⱼ + 2μ εᵢⱼ
```

### Viscoelastic Models

For time-dependent deformation:

#### Maxwell Model
```ini
[ROCK]
constitutive_model = VISCOELASTIC_MAXWELL

youngs_modulus = 20.0e9
poisson_ratio = 0.25
viscosity = 1.0e18         # Pa·s
relaxation_time = 1.0e12   # s = η/E
```

#### Kelvin-Voigt Model
```ini
[ROCK]
constitutive_model = VISCOELASTIC_KELVIN

youngs_modulus = 20.0e9
poisson_ratio = 0.25
viscosity = 1.0e17         # Pa·s
```

#### Standard Linear Solid (SLS)
```ini
[ROCK]
constitutive_model = VISCOELASTIC_SLS

# Spring-dashpot parameters
E1 = 20.0e9                # Instantaneous modulus
E2 = 15.0e9                # Long-term modulus
eta = 5.0e17               # Viscosity
```

### Poroelastic Model

Biot's theory for fluid-saturated porous media:

```ini
[SIMULATION]
solid_model = POROELASTIC
enable_geomechanics = true

[ROCK]
constitutive_model = POROELASTIC

# Drained elastic properties
youngs_modulus = 20.0e9
poisson_ratio = 0.25

# Poroelastic parameters
biot_coefficient = 0.8     # α (0 < α ≤ 1)
biot_modulus = 1.0e10      # M (Pa)

# Alternative specification
# undrained_poisson = 0.35  # νᵤ
# skempton = 0.9           # B
```

**Governing Equations:**
```
Equilibrium:  ∇·σ + ρb = 0
Constitutive: σᵢⱼ = λ εₖₖ δᵢⱼ + 2μ εᵢⱼ - α p δᵢⱼ
Flow:         ∂ζ/∂t + ∇·q = Q,  where ζ = α εᵥ + p/M
```

### Elastoplastic Models

For rock failure and permanent deformation:

#### Mohr-Coulomb
```ini
[ROCK]
constitutive_model = ELASTOPLASTIC_MC

failure_criterion = MOHR_COULOMB

# Elastic properties
youngs_modulus = 20.0e9
poisson_ratio = 0.25

# Yield surface parameters
cohesion = 5.0e6           # Pa
friction_angle = 30.0      # degrees
dilation_angle = 10.0      # degrees (non-associated flow)
tensile_strength = 2.0e6   # Pa (tensile cutoff)

# Hardening (optional)
hardening_type = NONE      # NONE, LINEAR, EXPONENTIAL
```

#### Drucker-Prager
```ini
[ROCK]
constitutive_model = ELASTOPLASTIC_DP

failure_criterion = DRUCKER_PRAGER

youngs_modulus = 20.0e9
poisson_ratio = 0.25
cohesion = 5.0e6
friction_angle = 30.0

# DP matching to M-C
dp_matching = INNER_EDGE   # INNER_EDGE, OUTER_EDGE, PLANE_STRAIN
```

### Anisotropic Models

#### VTI (Vertical Transverse Isotropy)

For horizontally layered media:

```ini
[ROCK]
constitutive_model = VTI

# Full stiffness tensor (Voigt notation)
c11 = 50.0e9               # Pa
c12 = 15.0e9
c13 = 12.0e9
c33 = 40.0e9               # Vertical stiffness
c44 = 12.0e9               # Shear stiffness

# Alternative: Thomsen parameters
# vp0 = 4000.0             # Vertical P-wave velocity
# vs0 = 2500.0             # Vertical S-wave velocity
# thomsen_epsilon = 0.15   # P-wave anisotropy
# thomsen_delta = 0.05     # Near-vertical P anisotropy
# thomsen_gamma = 0.10     # S-wave anisotropy
```

#### HTI (Horizontal Transverse Isotropy)

For vertically fractured media:

```ini
[ROCK]
constitutive_model = HTI

c11 = 40.0e9
c12 = 12.0e9
c13 = 15.0e9
c33 = 50.0e9
c44 = 15.0e9

# Symmetry axis azimuth (from North)
symmetry_azimuth = 45.0    # degrees
```

---

## Thermal Models

### Heat Conduction

```ini
[SIMULATION]
enable_thermal = true

[ROCK]
thermal_conductivity = 2.5     # W/(m·K)
specific_heat = 800.0          # J/(kg·K)
density = 2600.0               # kg/m³

[FLUID]
thermal_conductivity = 0.6     # W/(m·K)
specific_heat = 4180.0         # J/(kg·K) (water)
```

### Convective Heat Transfer

```ini
[SIMULATION]
enable_thermal = true

[THERMAL]
enable_convection = true       # Flow carries heat
enable_conduction = true       # Diffusive heat transfer

# Effective thermal properties
thermal_model = PARALLEL       # PARALLEL, SERIES, GEOMETRIC
```

### Thermal Expansion

```ini
[ROCK]
thermal_expansion = 1.0e-5     # 1/K (volumetric)
# or
thermal_expansion_linear = 3.3e-6  # 1/K (linear)

[FLUID]
thermal_expansion = 2.0e-4     # 1/K
```

---

## Coupled Physics

### Hydro-Mechanical (HM) Coupling

Pressure changes cause deformation, deformation affects permeability:

```ini
[SIMULATION]
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true

[ROCK]
constitutive_model = POROELASTIC
biot_coefficient = 0.8
biot_modulus = 1.0e10

# Permeability coupling
permeability_model = STRESS_DEPENDENT
# k/k₀ = exp(-c × Δσ')
perm_stress_coefficient = 1.0e-8  # 1/Pa
```

### Thermo-Hydro (TH) Coupling

Temperature affects fluid properties:

```ini
[SIMULATION]
enable_thermal = true

[FLUID]
# Temperature-dependent viscosity
viscosity_model = ARRHENIUS
viscosity_ref = 0.001          # Pa·s at T_ref
temperature_ref = 300.0        # K
activation_energy = 15000.0    # J/mol
```

### Thermo-Hydro-Mechanical (THM) Coupling

Full coupling for geothermal, nuclear waste, etc.:

```ini
[SIMULATION]
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_thermal = true

[ROCK]
constitutive_model = POROELASTIC
biot_coefficient = 0.8
thermal_expansion = 1.0e-5

# Temperature-dependent properties
youngs_modulus_T0 = 20.0e9     # At reference T
dE_dT = -1.0e7                 # Pa/K (softening with T)
```

---

## Physics Kernels Architecture

FSRM uses modular physics kernels that implement residuals and Jacobians:

```cpp
// Base physics kernel interface
class PhysicsKernel {
public:
    // Residual: F(u) = 0
    virtual void residual(const PetscScalar u[], 
                         const PetscScalar u_t[],
                         const PetscScalar u_x[], 
                         PetscScalar f[]) = 0;
    
    // Jacobian: dF/du
    virtual void jacobian(const PetscScalar u[], 
                         const PetscScalar u_t[],
                         const PetscScalar u_x[], 
                         PetscScalar J[]) = 0;
};

// Available kernels
SinglePhaseFlowKernel     // Darcy flow
BlackOilKernel            // Three-phase
CompositionalKernel       // Multi-component EOS
GeomechanicsKernel        // Elasticity/plasticity
ThermalKernel             // Heat transfer
ParticleTransportKernel   // Proppant/tracer
FracturePropagationKernel // Cohesive zone
ElastodynamicsKernel      // Wave propagation
PoroelastodynamicsKernel  // Coupled waves
```

### Adding Custom Physics (Advanced)

Create a new kernel:

```cpp
class MyCustomKernel : public PhysicsKernel {
public:
    MyCustomKernel() : PhysicsKernel(PhysicsType::CUSTOM) {}
    
    void residual(...) override {
        // Implement your physics
        f[0] = /* your equation */;
    }
    
    void jacobian(...) override {
        // Implement Jacobian
        J[0] = /* derivative */;
    }
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 1; }
};
```

---

## Selecting the Right Model

### Decision Tree

```
START
  │
  ├─ Single fluid type? ─────────────► SINGLE_COMPONENT
  │     │
  │     └─ Multiple phases? ─────────► BLACK_OIL
  │           │
  │           └─ Compositional effects? ─► COMPOSITIONAL
  │
  ├─ CO2 storage? ───────────────────► CO2
  │
  └─ Need mechanics?
        │
        ├─ Time-dependent? ──────────► VISCOELASTIC
        │
        ├─ Fluid coupling? ──────────► POROELASTIC
        │
        ├─ Failure/yielding? ────────► ELASTOPLASTIC
        │
        ├─ Layered/fractured? ───────► VTI/HTI
        │
        └─ Simple deformation? ──────► ELASTIC
```

### Model Selection by Application

| Application | Flow Model | Solid Model | Coupling |
|------------|------------|-------------|----------|
| Water flooding | BLACK_OIL | ELASTIC | None |
| Shale gas | SINGLE_COMPONENT | POROELASTIC | HM |
| Geothermal | SINGLE_COMPONENT | POROELASTIC | THM |
| CO2 storage | CO2 | ELASTOPLASTIC | HM |
| EOR/IOR | BLACK_OIL | POROELASTIC | HM |
| Gas injection | COMPOSITIONAL | ELASTIC | None |
| Seismicity | SINGLE_COMPONENT | POROELASTIC | HM+Faults |
| Salt creep | SINGLE_COMPONENT | VISCOELASTIC | HM |

### Computational Cost

| Model | Relative Cost | Memory |
|-------|--------------|--------|
| SINGLE_COMPONENT | 1× | Low |
| BLACK_OIL | 3-5× | Medium |
| COMPOSITIONAL | 10-20× | High |
| + ELASTIC | +1.5× | +Medium |
| + POROELASTIC | +3× | +High |
| + THERMAL | +1.5× | +Medium |

---

## Example Configurations

### Shale Gas Production

```ini
[SIMULATION]
name = shale_gas
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_fractures = true

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.05
permeability_x = 0.001         # Ultra-low perm
youngs_modulus = 30.0e9
biot_coefficient = 0.9

[FLUID]
type = SINGLE_PHASE
density = 100.0                # Gas
viscosity = 2.0e-5             # Gas viscosity
compressibility = 1.0e-7
```

### Geothermal Doublet

```ini
[SIMULATION]
name = geothermal_doublet
fluid_model = SINGLE_COMPONENT
solid_model = POROELASTIC
enable_geomechanics = true
enable_thermal = true

[ROCK]
constitutive_model = POROELASTIC
porosity = 0.15
permeability_x = 500.0
thermal_conductivity = 2.8
biot_coefficient = 0.7

[FLUID]
type = BRINE
salinity = 30000
thermal_conductivity = 0.6
specific_heat = 4000.0

[IC1]
field = TEMPERATURE
distribution = GRADIENT
value = 300.0                  # Surface T
gradient = 0.0, 0.0, 0.035     # 35°C/km
```

---

**Previous**: [← Configuration System](02_CONFIGURATION.md) | **Next**: [Well Modeling →](04_WELL_MODELING.md)
