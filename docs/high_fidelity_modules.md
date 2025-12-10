# High-Fidelity Modules Documentation

This document provides comprehensive documentation for the high-fidelity physics and numerical modules in FSRM. These modules extend the baseline capabilities with advanced physics models, sophisticated numerical methods, and fully-coupled multi-physics capabilities.

## Table of Contents

1. [Overview](#overview)
2. [High-Fidelity Fluid Flow](#high-fidelity-fluid-flow)
3. [High-Fidelity Geomechanics](#high-fidelity-geomechanics)
4. [Advanced Coupled Physics](#advanced-coupled-physics)
5. [High-Fidelity Numerics](#high-fidelity-numerics)
6. [Configuration Guide](#configuration-guide)
7. [Testing](#testing)
8. [Best Practices](#best-practices)

---

## Overview

The high-fidelity modules are designed as optional extensions that complement existing simpler models. They provide:

- **Higher physical fidelity**: More accurate representation of complex phenomena
- **Better coupling**: Seamless integration between different physics domains
- **Advanced numerics**: Stabilization, adaptivity, and PETSc-based solvers
- **Comprehensive testing**: Unit, MMS, and integration tests

### Design Philosophy

1. **Non-invasive**: Can be enabled/disabled without affecting base code
2. **Modular**: Each component can be used independently or combined
3. **Validated**: Verified against analytical solutions and benchmarks
4. **Well-documented**: Extensive API documentation and usage examples

---

## High-Fidelity Fluid Flow

### Module: `HighFidelityFluidFlow.hpp/cpp`

Provides advanced fluid flow models beyond standard Darcy flow.

### Non-Darcy Flow (`NonDarcyFlow`)

Models high-velocity flow effects in porous media.

#### Supported Models

| Model | Equation | Use Case |
|-------|----------|----------|
| **Forchheimer** | `-∇p = (μ/K)v + βρv²` | Near-wellbore flow, fractured media |
| **Klinkenberg** | `K_app = K_abs(1 + b/p)` | Gas flow in tight formations |
| **Ergun** | Combined viscous + inertial | High-porosity media |

#### Parameters

```cpp
struct Parameters {
    NonDarcyModel model;       // FORCHHEIMER, KLINKENBERG, ERGUN
    double forchheimer_beta;   // Inertial coefficient [1/m]
    double klinkenberg_b;      // Slip factor [Pa]
    double ergun_a, ergun_b;   // Ergun equation constants
};
```

#### Usage Example

```cpp
NonDarcyFlow ndf;
NonDarcyFlow::Parameters params;
params.model = NonDarcyModel::FORCHHEIMER;
params.forchheimer_beta = 1e8;
ndf.setParameters(params);

double correction = ndf.calculateCorrection(velocity, density, viscosity, permeability);
double effective_velocity = ndf.effectiveVelocity(grad_p, rho, mu, K);
```

### Dual/Triple Porosity (`DualPorosityPermeability`)

Models naturally fractured reservoirs with distinct matrix and fracture domains.

#### Supported Models

- **Warren-Root**: Classical dual-porosity
- **Kazemi**: Modified transfer function
- **Triple Porosity**: Matrix + fracture + vugs

#### Key Methods

```cpp
double calculateTransferRate(double p_matrix, double p_fracture);
double effectivePorosity();
double effectivePermeability();
double storativityRatio(double c_matrix, double c_fracture);
double interporosityFlowCoefficient(double length_scale);
```

### Non-Isothermal Multiphase Flow (`NonIsothermalMultiphaseFlow`)

Thermal effects in multiphase flow systems.

#### Features

- Temperature-dependent fluid properties (viscosity, density)
- Effective thermal conductivity (solid + fluids)
- Viscous heating
- Joule-Thomson cooling/heating

#### Key Methods

```cpp
double effectiveThermalConductivity(phi, S_w, S_o, S_g);
double effectiveHeatCapacity(phi, S_w, S_o, S_g, rho_s, rho_w, rho_o, rho_g);
double viscosityAtTemperature(mu_ref, T_ref, T, E_activation);
double viscousHeatingRate(mu, grad_p, K);
double jouleThomsonTemperatureChange(JT_coeff, delta_p);
```

### Dynamic Relative Permeability (`DynamicRelativePermeability`)

Velocity and capillary-number dependent relative permeability.

#### Models

- **Capillary Number**: Desaturation at high capillary numbers
- **Fingering**: Viscous fingering corrections

#### Key Methods

```cpp
double capillaryNumber(mu, v, sigma);
double desaturationCoefficient(Nc);
double dynamicResidualSaturation(S_or_static, D);
double modifyRelativePermeability(kr_static, Sw, mu, v, sigma);
```

### Miscible Flow (`MiscibleFlow`)

CO2 flooding, solvent injection, and miscible displacement.

#### Models

- **Todd-Longstaff**: Mixing parameter approach
- **First Contact Miscibility**: Sharp MMP transition
- **Multi-Contact Miscibility**: Gradual miscibility development

#### Key Methods

```cpp
double miscibilityFraction(double pressure);
double effectiveViscosity(mu_oil, mu_solvent, S_solvent, f_misc);
double effectiveDensity(rho_oil, rho_solvent, S_solvent, f_misc);
void dispersionTensor(velocity, D_molecular, D_effective);
```

---

## High-Fidelity Geomechanics

### Module: `HighFidelityGeomechanics.hpp/cpp`

Advanced constitutive models for geomaterials.

### Finite Strain Mechanics (`FiniteStrainMechanics`)

Large deformation hyperelastic models.

#### Supported Models

| Model | Strain Energy | Material Type |
|-------|---------------|---------------|
| **Neo-Hookean** | `W = μ/2(I₁-3) - μ·ln(J) + λ/2(ln J)²` | Rubber-like, soft rock |
| **Mooney-Rivlin** | `W = C₁₀(I₁-3) + C₀₁(I₂-3)` | More general rubber |
| **Ogden** | `W = Σ μₙ/αₙ(λ₁^αₙ + λ₂^αₙ + λ₃^αₙ - 3)` | Highly nonlinear |

#### Key Methods

```cpp
void computeStress(const double F[3][3], double stress[3][3]);
double computeStrainEnergy(const double F[3][3]);
void computeMaterialTangent(const double F[3][3], double C[3][3][3][3]);
```

### Finite Strain Plasticity (`FiniteStrainPlasticity`)

Multiplicative plasticity for large deformations.

#### Features

- Multiplicative decomposition: `F = Fᵉ · Fᵖ`
- J2 (von Mises) plasticity
- Isotropic/kinematic hardening

### Viscoplasticity (`Viscoplasticity`)

Rate-dependent plastic deformation.

#### Models

- **Perzyna**: `ε̇ᵖ = γ · ⟨f/σ_y⟩ⁿ · ∂f/∂σ`
- **Duvaut-Lions**: `σ̇ = (σ_∞ - σ)/τ`

#### Key Methods

```cpp
double computeOverstress(const double stress[6], double eps_p_eq);
void computePlasticStrainRate(const double stress[6], double eps_p_eq, double eps_p_dot[6]);
double equivalentPlasticStrainRate(const double stress[6], double eps_p_eq);
```

### Creep Models (`CreepModel`)

Time-dependent deformation under sustained load.

#### Supported Models

| Model | Equation | Application |
|-------|----------|-------------|
| **Power Law** | `ε̇ = A·σⁿ·exp(-Q/RT)` | High-temperature rock |
| **Norton** | `ε̇ = A·σⁿ` | Simplified power law |
| **Burgers** | Maxwell + Kelvin | Transient + steady creep |

#### Key Methods

```cpp
double creepStrainRate(double stress);
double burgersCreepStrain(double stress, double time);
void burgersCreepComponents(stress, time, elastic, transient, steady);
```

### Hypoplasticity (`Hypoplasticity`)

Non-linear, rate-independent model for granular materials.

#### Features

- Critical state framework
- Barotropy and pyknotropy
- Void ratio evolution

#### Key Methods

```cpp
void computeStressRate(stress, strain_rate, void_ratio, stress_rate);
double criticalVoidRatio(double mean_stress);
double relativeDensity(double e, double p);
```

### Gradient-Enhanced Damage (`GradientEnhancedDamage`)

Regularized continuum damage with nonlocal effects.

#### Features

- Internal length scale for mesh-independent damage zone
- Various damage evolution laws (exponential, linear softening, Mazars)
- Helmholtz-type gradient regularization

#### Key Methods

```cpp
double computeDamage(double kappa);
double equivalentStrain(const double strain[6]);
void effectiveStress(const double stress[6], double damage, double stress_eff[6]);
double gradientParameter();  // Returns l² (internal length squared)
```

---

## Advanced Coupled Physics

### Module: `AdvancedCoupledPhysics.hpp/cpp`

Multi-physics coupling frameworks.

### Full Biot Dynamic Poroelasticity (`FullBiotDynamics`)

Complete Biot theory including inertia and slow wave.

#### Formulations

- **u-p**: Displacement-pressure (standard)
- **u-w**: Displacement-relative fluid displacement
- **u-p-w**: Three-field formulation

#### Key Features

- Fast and slow P-waves
- Frequency-dependent behavior
- Tortuosity effects

#### Physical Relations

```cpp
// Gassmann's relation
K_u = K_d + α²M

// Skempton's coefficient
B = αM / K_u

// Wave velocities
V_p_undrained, V_p_drained, V_s, V_p_slow
```

#### Key Methods

```cpp
double biotCoefficient();
double biotModulus();
double skemptonCoefficient();
double undrainedBulkModulus();
double pwaveVelocityDrained();
double pwaveVelocityUndrained();
double swaveVelocity();
double slowPwaveVelocity();
double characteristicFrequency();
bool isHighFrequencyRegime(double frequency);
void effectiveStress(total_stress, pore_pressure, effective_stress);
```

### THM Coupling (`THMCoupling`)

Thermo-Hydro-Mechanical coupling for reservoirs and geothermal systems.

#### Coupling Terms

- **Thermal → Mechanical**: Thermal expansion/contraction stresses
- **Thermal → Hydraulic**: Thermal pressurization, viscosity changes
- **Mechanical → Hydraulic**: Stress-dependent permeability
- **Hydraulic → Mechanical**: Pore pressure (Biot) coupling

#### Key Methods

```cpp
double thermalStress(double delta_T);
double thermalExpansionStrain(double delta_T);
double thermalPressurization(double delta_T, double bulk_modulus);
double effectiveThermalExpansion();
double thermalDiffusivity(double density);
double effectiveThermalConductivity(phi, k_solid, k_fluid);
double consolidationCoefficient(K, mu, K_bulk, phi);
```

### THMC Coupling (`THMCCoupling`)

Thermo-Hydro-Mechanical-Chemical coupling.

#### Chemical Coupling Effects

- Mineral dissolution/precipitation → porosity change
- Chemical swelling (clay hydration)
- Reaction heat (exothermic/endothermic)
- Chemical osmosis

#### Key Methods

```cpp
double chemicalStrain(const std::vector<double>& concentrations);
double chemicalStrainRate(const std::vector<double>& concentration_rates);
double reactionHeat(const std::vector<double>& rates, const std::vector<double>& enthalpies);
double porosityChangeRate(double diss_rate, double precip_rate, double molar_volume);
double effectiveDiffusivity(double D_free, double phi, double tau);
double temperatureDependentRate(double k0, double Ea, double T, double T_ref);
```

### Unsaturated Flow Coupling (`UnsaturatedFlowCoupling`)

Variably saturated soil mechanics.

#### Retention Models

- **Van Genuchten**: Standard model with α, n, m parameters
- **Brooks-Corey**: Power-law model with λ parameter

#### Bishop's Effective Stress

`σ' = σ - [χ·p_w + (1-χ)·p_a]`

#### Key Methods

```cpp
double waterSaturation(double suction);
double relativePermeability(double effective_saturation);
double bishopParameter(double saturation);
double effectiveStress(sigma_total, p_water, p_air, chi);
double moistureCapacity(double suction);
```

### Multi-Scale Coupling (`MultiScaleCoupling`)

Homogenization and computational multi-scale methods.

#### Methods

- **FE²**: Fully computational (RVE at each integration point)
- **HMM**: Heterogeneous Multiscale Method
- **Asymptotic Homogenization**: Classical approach

#### Homogenization Bounds

```cpp
double voigtAverage(moduli, fractions);  // Upper bound
double reussAverage(moduli, fractions);  // Lower bound
void hashinShtrikmanBounds(K1, mu1, K2, mu2, f1, K_lower, K_upper, mu_lower, mu_upper);
void moriTanaka(K_m, mu_m, K_i, mu_i, f_i, K_eff, mu_eff);
void selfConsistent(K1, mu1, K2, mu2, f1, K_eff, mu_eff);
```

---

## High-Fidelity Numerics

### Module: `HighFidelityNumerics.hpp/cpp`

Advanced numerical methods for stability and accuracy.

### SUPG Stabilization (`SUPGStabilization`)

Streamline Upwind Petrov-Galerkin for advection-dominated problems.

#### Methods

- **SUPG**: Classic streamline-upwind
- **GLS**: Galerkin Least-Squares
- **PSPG**: Pressure-Stabilized Petrov-Galerkin (for Stokes)
- **VMS**: Variational Multi-Scale

#### Tau Formulas

- **Shakib**: Full transient + advection + diffusion + reaction
- **Tezduyar**: Simplified with element-based scaling
- **Franca-Valentin**: Original formulation

#### Key Methods

```cpp
double calculateTau(velocity, diffusivity, reaction, h, dt);
double pecletNumber(velocity_mag, h, diffusivity);
double optimalDiffusivity(velocity_mag, h, diffusivity);
void addStabilization(velocity, test_gradient, residual, tau, stab_term);
double pspgTerm(grad_p, grad_q, tau);
double shockCapturingDiffusivity(grad_u, u, h, residual);
```

### PETSc Time Integration (`PETScTimeIntegration`)

Wrapper for PETSc TS (Time Stepping) module.

#### Supported Schemes

| Category | Schemes |
|----------|---------|
| **BDF** | BDF1, BDF2, BDF3, BDF4 |
| **Runge-Kutta** | RK4, SSPRK3 |
| **Theta** | Forward Euler, Backward Euler, Crank-Nicolson |
| **IMEX** | ARKIMEX (additive Runge-Kutta) |
| **Second-order** | Newmark-Beta, HHT-Alpha, Generalized-Alpha |

#### Key Methods

```cpp
void create(DM dm);
void setType(TimeScheme scheme);
void setAdaptivity(bool adaptive, double atol, double rtol);
void setRHSFunction(TSRHSFunction func, void* ctx);
void setIFunction(TSIFunction func, void* ctx);
void setIJacobian(TSIJacobian func, void* ctx);
void solve(Vec u, double t_final);
int getOrder();
bool isAStable();
bool isLStable();
```

### p-Adaptivity (`PAdaptivityManager`)

Polynomial order adaptation for hp-FEM.

#### Key Methods

```cpp
int decideOrderChange(double smoothness, double error_indicator, int current_order);
double smoothnessIndicator(Vec solution, DM dm);
void changeOrder(DM dm, Vec solution, int new_order, Vec interpolated);
```

### Physics-Based Preconditioning (`PETScPhysicsPreconditioner`)

Block and Schur complement preconditioners.

#### Types

- **Block Jacobi/Gauss-Seidel**: Field-split block methods
- **Schur Complement**: For saddle-point systems
- **GAMG**: Algebraic multigrid
- **Custom Field-Split**: User-defined block structure

#### Key Methods

```cpp
void create(KSP ksp);
void setType(PreconditionerType type);
void setSchurFactorization(SchurFactType fact, SchurPreType pre);
void setupFieldSplit(DM dm, const std::vector<std::string>& field_names);
void setupMultigrid(int levels, const std::string& smoother);
```

---

## Configuration Guide

### Configuration File Format

```ini
[HighFidelity]
# Enable high-fidelity extensions
enable_high_fidelity = true

# Fluid Flow Options
enable_non_darcy = true
non_darcy_model = forchheimer
forchheimer_beta = 1e8

enable_dual_porosity = true
dual_porosity_model = warren_root
matrix_porosity = 0.15
fracture_porosity = 0.02

# Geomechanics Options
enable_finite_strain = true
hyperelastic_model = neo_hookean
shear_modulus = 1e6
bulk_modulus = 1e9

enable_creep = true
creep_model = power_law
creep_A = 1e-20
creep_n = 3.0

# Coupled Physics Options
enable_full_biot = true
biot_formulation = u_p
biot_alpha = 0.8

enable_thm = true
thermal_expansion_solid = 1e-5

# Numerics Options
enable_supg = true
supg_method = supg
supg_tau_formula = shakib

enable_advanced_time = true
time_scheme = bdf2
time_adaptive = true
```

### Programmatic Configuration

```cpp
// Create configurations
HighFidelityFluidFlowConfig fluid_config;
fluid_config.enable_non_darcy = true;
fluid_config.non_darcy_params.model = NonDarcyModel::FORCHHEIMER;
fluid_config.non_darcy_params.forchheimer_beta = 1e8;

HighFidelityGeomechanicsConfig geomech_config;
geomech_config.enable_finite_strain = true;
geomech_config.finite_strain_params.model = HyperelasticModel::NEO_HOOKEAN;

AdvancedCoupledPhysicsConfig coupled_config;
coupled_config.enable_full_biot = true;
coupled_config.biot_params.formulation = BiotFormulation::U_P;

// Validate configurations
std::string error_msg;
if (!fluid_config.validate(error_msg)) {
    throw std::runtime_error("Invalid configuration: " + error_msg);
}
```

---

## Testing

### Test Categories

| Category | Purpose | Location |
|----------|---------|----------|
| **Unit Tests** | Individual component testing | `tests/unit/test_high_fidelity_*.cpp` |
| **MMS Tests** | Convergence verification | `tests/physics/test_mms_high_fidelity.cpp` |
| **Coupling Tests** | Cross-module consistency | `tests/physics/test_coupling_consistency.cpp` |
| **Integration Tests** | End-to-end scenarios | `tests/integration/test_high_fidelity_integration.cpp` |

### Running Tests

```bash
# All high-fidelity tests
make test_high_fidelity_all

# Unit tests only
make test_high_fidelity_unit

# MMS and coupling tests
make test_high_fidelity_mms

# Integration tests
make test_high_fidelity_integration

# Specific test filter
./run_all_tests --gtest_filter="*FiniteStrain*"
```

### Adding New Tests

1. Create test file in appropriate directory
2. Add to `CMakeLists.txt` in the correct source list
3. Register with CTest
4. Set appropriate timeout and labels

---

## Best Practices

### When to Use High-Fidelity Models

| Scenario | Recommended Model |
|----------|-------------------|
| High-velocity wellbore flow | Non-Darcy (Forchheimer) |
| Naturally fractured reservoirs | Dual Porosity |
| CO2 sequestration | THMC Coupling + Miscible Flow |
| Geothermal systems | THM Coupling |
| Soft rock deformation | Finite Strain |
| Salt cavern creep | Viscoplastic + Creep |
| Hydraulic fracturing | Gradient Damage + Biot Dynamics |
| Unsaturated near-surface | Unsaturated Flow Coupling |

### Performance Considerations

1. **Start simple**: Begin with base models, add complexity as needed
2. **Check convergence**: MMS tests verify spatial accuracy
3. **Monitor timestep**: Adaptive time stepping recommended for stiff systems
4. **Use appropriate preconditioners**: Field-split for coupled systems
5. **Profile**: Identify computational bottlenecks

### Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| Non-convergence with large deformation | Reduce load increment, use arc-length |
| Oscillations in advection | Enable SUPG stabilization |
| Stiff coupled systems | Use implicit time integration, field-split PC |
| Damage localization | Increase internal length scale |
| Dual-porosity instability | Check transfer coefficient magnitude |

---

## References

### Biot Poroelasticity
- Biot, M.A. (1941). General theory of three-dimensional consolidation
- Wang, H.F. (2000). Theory of Linear Poroelasticity

### THM Coupling
- Coussy, O. (2004). Poromechanics
- Rutqvist, J. et al. (2002). Coupled THMC modeling

### Hyperelasticity
- Holzapfel, G.A. (2000). Nonlinear Solid Mechanics
- Bonet, J. & Wood, R.D. (2008). Nonlinear Continuum Mechanics for FEA

### Stabilization Methods
- Brooks, A.N. & Hughes, T.J.R. (1982). SUPG formulation
- Franca, L.P. & Frey, S.L. (1992). Stabilized finite elements

### Damage Mechanics
- Peerlings, R.H.J. et al. (1996). Gradient-enhanced damage
- Bazant, Z.P. & Jirasek, M. (2002). Nonlocal damage models
