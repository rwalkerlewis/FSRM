# Wave Physics Theory and Implementation

## Overview

This document describes the implementation of elastodynamic and poroelastodynamic wave physics in the reservoir simulator, including static-to-dynamic triggering for stress-induced seismicity and dynamic permeability changes from transient waves.

## 1. Elastodynamics

### 1.1 Governing Equations

The elastodynamic equation with inertia is:

```
ρ ∂²u/∂t² = ∇·σ + f
```

where:
- `ρ` is the density
- `u` is the displacement vector
- `σ` is the Cauchy stress tensor
- `f` is the body force

For linear elasticity, the stress-strain relationship is:

```
σ = λ(∇·u)I + 2μ∇ˢu
```

where `λ` and `μ` are Lamé parameters, and `∇ˢu` is the symmetric gradient (strain tensor).

### 1.2 Wave Velocities

P-wave (compressional) velocity:
```
vₚ = √[(λ + 2μ)/ρ]
```

S-wave (shear) velocity:
```
vₛ = √[μ/ρ]
```

### 1.3 Rayleigh Damping

To model energy dissipation, we use Rayleigh damping:

```
C = α·M + β·K
```

where:
- `C` is the damping matrix
- `M` is the mass matrix
- `K` is the stiffness matrix
- `α` is the mass-proportional damping coefficient
- `β` is the stiffness-proportional damping coefficient

The damping ratio `ζ` at frequency `ω` is:

```
ζ(ω) = (α/2ω) + (βω/2)
```

### 1.4 Quality Factor

The seismic quality factor Q relates to damping:

```
Q = 1/(2ζ)
```

Typical values:
- Granite: Q ~ 200-1000
- Sandstone: Q ~ 50-100
- Sediments: Q ~ 10-50

## 2. Poroelastodynamics (Biot's Theory)

### 2.1 Governing Equations

Biot's poroelastodynamic equations couple solid deformation and fluid flow with inertial effects:

**Momentum balance:**
```
ρ ∂²u/∂t² + ρₐ ∂²w/∂t² = ∇·σ - α∇p + f
```

**Darcy's law (dynamic):**
```
ρₐ ∂²u/∂t² + (ρf/ϕ) ∂²w/∂t² = -∇p - (η/k)∂w/∂t
```

**Mass conservation:**
```
∂/∂t[α∇·u + p/M] + ∇·(∂w/∂t) = Q
```

where:
- `u` = solid displacement
- `w` = relative fluid displacement (w = ϕ(uᶠ - uˢ))
- `p` = pore pressure
- `ρ = (1-ϕ)ρₛ + ϕρf` = bulk density
- `ρₐ = ϕρf` = added mass density
- `α` = Biot coefficient
- `M` = Biot modulus
- `η` = fluid viscosity
- `k` = permeability
- `ϕ` = porosity

### 2.2 Biot Wave Types

Poroelastic media support three wave types:

1. **Fast P-wave (Biot wave I)**: Similar to elastic P-wave, fluid and solid move in phase
2. **Shear wave (S-wave)**: Unchanged from elastic case
3. **Slow P-wave (Biot wave II)**: Diffusive wave, fluid and solid move out of phase

### 2.3 Wave Velocities

Fast P-wave velocity (simplified):
```
vₚ₁ ≈ √[(K_d + 4μ/3)/ρ]
```

where `K_d` is the drained bulk modulus.

Slow P-wave is highly attenuated and velocity depends strongly on frequency and permeability.

### 2.4 Frequency-Dependent Behavior

Characteristic frequency (Biot frequency):
```
ω_c = ηϕ/(ρf·k)
```

- Low frequency (ω << ω_c): Quasi-static Biot theory
- High frequency (ω >> ω_c): Dynamic decoupling, waves behave elastically

## 3. Static-to-Dynamic Triggering

### 3.1 Concept

Static stress changes from injection, tectonic loading, or tidal forces can trigger dynamic rupture when:

```
τ + Δτ_static ≥ τ_threshold
```

where:
- `τ` is the background shear stress
- `Δτ_static` is the static stress perturbation
- `τ_threshold` is the failure threshold

### 3.2 Implementation

The simulator monitors stress state during quasi-static evolution. When threshold is exceeded:

1. Switch from quasi-static to fully dynamic mode
2. Activate inertial terms in governing equations
3. Reduce timestep to resolve wave propagation
4. Run dynamic simulation for specified duration
5. Optionally return to quasi-static mode

### 3.3 Coulomb Failure Criterion

For fault triggering:

```
CFS = τ + μ(σₙ + Δp) - S₀
```

where:
- CFS = Coulomb failure stress
- `τ` = shear stress on fault
- `μ` = friction coefficient
- `σₙ` = normal stress (compression positive)
- `Δp` = pore pressure change
- `S₀` = fault strength

Failure occurs when CFS > 0.

### 3.4 Rate-and-State Friction

For more realistic fault behavior, rate-and-state friction:

```
τ = (μ₀ + a·ln(V/V₀) + b·ln(θ/θ₀))σₙ
```

State evolution:
```
dθ/dt = 1 - (Vθ/Dc)
```

where:
- `V` = slip velocity
- `θ` = state variable
- `a`, `b` = empirical parameters
- `Dc` = critical slip distance

Unstable (velocity weakening) when b > a.

## 4. Dynamic Permeability Changes

### 4.1 Mechanisms

Permeability can change instantaneously during wave passage due to:

1. **Volumetric strain**: Pore opening/closing
2. **Deviatoric strain**: Shear dilation
3. **Effective stress changes**: Stress-dependent permeability
4. **Damage/fracturing**: Irreversible changes

### 4.2 Strain-Based Model

Exponential model:

```
k/k₀ = exp(β·εᵥₒₗ)
```

where:
- `β` ~ 100-1000 for rocks
- `εᵥₒₗ` = volumetric strain

Linear approximation for small strains:

```
k/k₀ ≈ 1 + β·εᵥₒₗ
```

### 4.3 Stress-Based Model

```
k/k₀ = exp(-γ·σ')
```

where:
- `γ` ~ 10⁻⁸ - 10⁻⁷ Pa⁻¹
- `σ'` = effective stress

### 4.4 Shear Dilation Model

Particularly important in fault zones:

```
k/k₀ = 1 + C_d·|γ|
```

where:
- `C_d` = dilation coefficient (~ 10-100)
- `γ` = shear strain magnitude

### 4.5 Time-Dependent Recovery

Permeability recovers toward initial value with characteristic time τ_r:

```
k(t) = k₀ + [k(t₀) - k₀]·exp[-(t-t₀)/τᵣ]
```

Typical recovery times:
- Granular materials: τᵣ ~ 10-100 s
- Fractured rock: τᵣ ~ 100-1000 s
- Fault zones: τᵣ ~ hours to days

### 4.6 Combined Model

The implementation uses a combined model:

```
k_new = max(k_strain, k_stress, k_shear)
```

with bounds:
```
k_min ≤ k_new ≤ k_max
```

### 4.7 Physical Constraints

Permeability changes are limited by:
- Minimum: Complete pore closure (k_min ~ 10⁻²¹ m²)
- Maximum: Free fracture flow (k_max ~ 10⁻¹⁰ m²)
- Typical range: 0.1-1000 mD for reservoir rocks

## 5. Numerical Implementation

### 5.1 Time Integration

**Quasi-static problems:**
- Backward Euler (implicit, unconditionally stable)

**Dynamic problems:**
- Generalized-α method (second-order accurate, unconditionally stable)
- Newmark-β method (optional)

### 5.2 Spatial Discretization

- Finite elements (continuous Galerkin)
- Quadratic elements for better wave resolution
- Minimum 10-20 nodes per wavelength

### 5.3 Timestep Selection

For wave propagation:

```
Δt ≤ h/(c·CFL)
```

where:
- `h` = element size
- `c` = wave velocity (use maximum)
- `CFL` = Courant number (typically 0.5-0.8)

### 5.4 Absorbing Boundaries

To prevent spurious reflections:
- Perfectly Matched Layers (PML)
- Paraxial boundaries
- Robin boundary conditions

## 6. Applications

### 6.1 Induced Seismicity

- Wastewater injection
- Hydraulic fracturing
- Geothermal operations
- CO₂ storage

### 6.2 Enhanced Recovery

- Seismic stimulation of oil/gas reservoirs
- Wave-enhanced permeability
- Vibroseis for improved production

### 6.3 Geothermal Systems

- Fracture network activation
- Permeability enhancement
- Heat extraction optimization

### 6.4 Natural Hazards

- Earthquake triggering
- Aftershock sequences
- Fault zone evolution
- Landslide initiation

## 7. Model Parameters

### 7.1 Typical Rock Properties

**Granite (crystalline basement):**
- E = 40-70 GPa
- ν = 0.22-0.27
- ρ = 2650-2750 kg/m³
- vₚ = 5000-6000 m/s
- vₛ = 2800-3500 m/s
- Q = 200-1000

**Sandstone (reservoir):**
- E = 10-25 GPa
- ν = 0.15-0.25
- ρ = 2200-2450 kg/m³
- vₚ = 2500-4000 m/s
- vₛ = 1500-2300 m/s
- Q = 50-150

**Limestone:**
- E = 30-60 GPa
- ν = 0.20-0.30
- ρ = 2500-2700 kg/m³
- vₚ = 3500-6000 m/s
- vₛ = 2000-3000 m/s
- Q = 100-500

### 7.2 Permeability Sensitivity

From laboratory and field observations:
- β (strain coeff): 10² - 10³
- γ (stress coeff): 10⁻⁸ - 10⁻⁷ Pa⁻¹
- Enhancement factor: 2-100x (up to 1000x for earthquakes)
- Recovery time: 10 s - 1 month

## 8. Validation

The implementation can be validated against:
- Analytical solutions (Lamb's problem, point source)
- Laboratory experiments (wave propagation in cores)
- Field observations (seismic records, permeability logs)
- Benchmark problems (NURETH, SPE comparative studies)

## References

1. Biot, M.A. (1956). "Theory of propagation of elastic waves in a fluid-saturated porous solid."
2. Detournay, E. & Cheng, A.H.-D. (1993). "Fundamentals of poroelasticity."
3. Manga, M. et al. (2012). "Changes in permeability caused by transient stresses."
4. Elkhoury, J.E. et al. (2006). "Seismic waves increase permeability."
5. Rice, J.R. & Cleary, M.P. (1976). "Some basic stress diffusion solutions for fluid-saturated elastic porous media."
