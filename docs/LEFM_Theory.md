# Linear Elastic Fracture Mechanics (LEFM) in ReservoirSim

## Overview

ReservoirSim implements comprehensive Linear Elastic Fracture Mechanics (LEFM) for hydraulic fracture propagation. This document explains the theory and implementation.

## Fundamental Concepts

### Stress Intensity Factor (K)

The stress intensity factor characterizes the stress field near a crack tip:

```
K_I = Y · σ · √(π·a)
```

where:
- `K_I` = Mode I (opening) stress intensity factor
- `Y` = geometry factor (depends on crack and body geometry)
- `σ` = applied stress
- `a` = crack half-length

### Critical Stress Intensity Factor (Fracture Toughness)

Material property representing resistance to crack propagation:

```
K_Ic = fracture toughness
```

Typical values:
- Soft shale: 0.5 - 1.0 MPa·m^0.5
- Medium shale: 1.0 - 2.0 MPa·m^0.5
- Hard rock: 2.0 - 4.0 MPa·m^0.5

### Initiation Criterion

Fracture initiates when:

```
K_I ≥ K_Ic
```

### Propagation Criterion

For propagation to continue:

```
K_I > K_Ic
```

Propagation velocity:

```
v = C · ((K_I - K_Ic) / K_Ic)^m
```

where C and m are material constants.

## Energy Approach

### Energy Release Rate (G)

Energy available for crack growth per unit area:

```
G = K_I² · (1 - ν²) / E
```

### Critical Energy Release Rate (G_c)

Material property (related to K_Ic):

```
G_c = K_Ic² · (1 - ν²) / E
```

### Irwin Criterion

Alternative propagation criterion:

```
G ≥ G_c
```

## Fracture Width Profile

### Sneddon's Solution

For pressurized crack in infinite elastic medium:

```
w(x) = (4·(1-ν²)/E) · √(a² - x²) · (P - σ_min)
```

Maximum width at center:

```
w_max = (4·(1-ν²)/E) · a · (P - σ_min)
```

where:
- `w(x)` = width at position x from center
- `E` = Young's modulus
- `ν` = Poisson's ratio
- `P` = internal fluid pressure
- `σ_min` = minimum confining stress
- `a` = crack half-length

## Hydraulic Fracture Models

### PKN Model (Perkins-Kern-Nordgren)

Assumptions:
- Height constrained
- Plane strain in vertical cross-sections
- Length >> Height

Width:
```
w = (4·(P - σ_h)·H) / E'
```

where E' = E/(1-ν²)

### KGD Model (Khristianovic-Geertsma-de Klerk)

Assumptions:
- Height >> Length
- Plane strain in horizontal plane

Width:
```
w = (4·(P - σ_h)·L) / E'
```

### Radial/Penny-Shaped Model

Assumptions:
- Circular crack
- Length ≈ Height

Width:
```
w_max = (2·R·K_I) / (E'·√π)
```

## Scaling Laws

### Toughness-Dominated Regime

When fracture energy dominates:

PKN:
```
L(t) ∝ t^(4/5)
w ∝ (Q·μ·E'·t)^(1/5)
```

KGD:
```
L(t) ∝ t^(2/3)
```

Radial:
```
R(t) ∝ t^(2/5)
```

### Viscosity-Dominated Regime

When fluid viscosity dominates:

```
L(t) ∝ t^(4/5)
w ∝ (Q·μ·t)^(1/4)
```

Different coefficients than toughness regime.

### Leak-off Dominated

When Carter leak-off is significant:

```
L(t) ∝ t^(2/3)
```

## Fluid Flow in Fracture

### Cubic Law (Poiseuille Flow)

Flow rate between parallel plates:

```
q = -(w³)/(12·μ) · dP/dx
```

Fracture permeability:
```
k_f = w²/12
```

### Reynolds Number

Flow regime indicator:

```
Re = (ρ·v·w) / μ
```

- Re < 2300: Laminar (cubic law valid)
- Re > 2300: Turbulent (corrections needed)

## Leak-off Models

### Carter's Model

Leak-off velocity:

```
v_L = C_L / √t
```

Total leaked volume:

```
V_L = 2·C_L·A·√t
```

where:
- `C_L` = leak-off coefficient (m/s^0.5)
- `A` = fracture area
- `t` = time since fracture arrival

Typical C_L values:
- Tight shale: 10^-8 m/s^0.5
- Medium permeability: 10^-7 m/s^0.5
- High permeability: 10^-6 m/s^0.5

## Coupled Physics

### Fluid-Solid Coupling

1. **Fluid pressure** → determines fracture width via elasticity
2. **Fracture width** → affects fluid flow via cubic law
3. **Fluid flow** → controls pressure evolution
4. **Pressure gradient** → drives fracture propagation

### Governing Equations

**Elasticity** (solid):
```
∇·σ = 0
w = f(P, E, ν, geometry)
```

**Fluid flow** (fracture):
```
∂(ρ·w)/∂t + ∇·(ρ·v) = -q_L
v = -(w²/12μ)·∇P
```

**Propagation** (LEFM):
```
v_tip = f(K_I, K_Ic)
K_I = Y·(P - σ)·√(π·a)
```

## Implementation in ReservoirSim

### Step-by-Step Algorithm

1. **Initialize** with small notch/perforation
2. **Solve elasticity** for stress field
3. **Compute K_I** at crack tip
4. **Check initiation**: if K_I ≥ K_Ic, activate fracture
5. **Propagation loop**:
   - Inject fluid (specified rate)
   - Update pressure field in fracture
   - Compute width from elasticity
   - Calculate flow distribution (cubic law)
   - Account for leak-off
   - Compute K_I at tip
   - Update crack length: da/dt = v(K_I)
   - Update stress field
6. **Output** fracture geometry, pressure, K_I

### Adaptive Time Stepping

Time step controlled by:
- Maximum K_I change
- Maximum pressure change
- Maximum propagation distance per step
- Courant condition for fluid flow

## Verification

### Analytical Solutions

Compare against:

1. **Sneddon crack** - width profile
2. **PKN model** - asymptotic solution
3. **KGD model** - asymptotic solution
4. **Radial fracture** - self-similar solution

### Convergence Tests

Check spatial and temporal convergence:
- Refine mesh near crack tip
- Reduce time step
- Verify K_I calculation accuracy

## Practical Considerations

### Stress Shadow Effects

Multiple fractures interact through stress field:
- Parallel fractures: compressive stress shadow
- Offset fractures: complex interaction
- Spacing affects propagation

### In-Situ Stress

Principal stresses determine:
- Fracture orientation (perpendicular to σ_min)
- Required pressure (breakdown ≈ 3·σ_min - σ_max - P_p)
- Propagation direction

### Temperature Effects

- Thermal stresses from cool fluid
- Changes in material properties (E, K_Ic)
- Thermoelastic effects

## Example Results Interpretation

### Length vs. Time Plot

- Initial rapid growth (high K_I)
- Slowing with time (decreasing K_I)
- Power law: L ∝ t^n
- n ≈ 0.8 (toughness) or 0.67 (leak-off)

### Width Profile

- Elliptical shape confirms LEFM
- Maximum at center
- Zero at tips
- Compare to Sneddon solution

### Net Pressure

Net pressure = P_frac - σ_min:
- 2-4 MPa typical for field
- Higher early (breakdown)
- Decreases with growth
- Function of rate, fluid, K_Ic

### K_I Evolution

- Increases during pressurization
- Jumps at initiation
- Decreases during propagation
- Oscillates if unstable

## References

1. **Sneddon, I.N. (1946)** - "The distribution of stress in the neighbourhood of a crack in an elastic solid", *Proc. R. Soc. Lond. A*

2. **Irwin, G.R. (1957)** - "Analysis of stresses and strains near the end of a crack traversing a plate", *J. Appl. Mech.*

3. **Perkins, T.K. & Kern, L.R. (1961)** - "Widths of hydraulic fractures", *J. Pet. Tech.*

4. **Geertsma, J. & de Klerk, F. (1969)** - "A rapid method of predicting width and extent of hydraulically induced fractures", *J. Pet. Tech.*

5. **Adachi, J. & Detournay, E. (2007)** - "Self-similar solution of a plane-strain fracture driven by a power-law fluid", *Int. J. Num. Anal. Meth. Geomech.*

6. **Bunger, A.P. et al. (2013)** - "Toughness-dominated hydraulic fracture with leak-off", *Int. J. Fract.*

## See Also

- `examples/ex_lefm_fracture_growth.cpp` - Implementation example
- `config/lefm_fracture_growth.config` - Configuration example
- `src/PhysicsKernel_LEFM.cpp` - LEFM calculations
