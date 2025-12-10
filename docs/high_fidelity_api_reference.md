# High-Fidelity Modules API Reference

Complete API reference for all high-fidelity classes and functions.

---

## Table of Contents

- [HighFidelityFluidFlow.hpp](#highfidelityfluidflowhpp)
- [HighFidelityGeomechanics.hpp](#highfidelitygeomechanicshpp)
- [AdvancedCoupledPhysics.hpp](#advancedcoupledphysicshpp)
- [HighFidelityNumerics.hpp](#highfidelitynumericshpp)

---

## HighFidelityFluidFlow.hpp

### Namespace: `FSRM::HighFidelity`

### Enumerations

#### `NonDarcyModel`
```cpp
enum class NonDarcyModel {
    FORCHHEIMER,   // Forchheimer equation
    KLINKENBERG,   // Gas slip effect
    ERGUN,         // Ergun equation
    BARREE_CONWAY  // Barree-Conway model
};
```

#### `DualPorosityModel`
```cpp
enum class DualPorosityModel {
    WARREN_ROOT,     // Warren-Root transfer
    KAZEMI,          // Kazemi transfer function
    TRIPLE_POROSITY, // Matrix + fracture + vugs
    MINC             // Multiple interacting continua
};
```

#### `DynamicRelPermModel`
```cpp
enum class DynamicRelPermModel {
    CAPILLARY_NUMBER, // Nc-dependent desaturation
    FINGERING,        // Viscous fingering
    RATE_DEPENDENT    // Rate-dependent hysteresis
};
```

#### `MiscibleFlowModel`
```cpp
enum class MiscibleFlowModel {
    TODD_LONGSTAFF,  // Mixing parameter approach
    FIRST_CONTACT,   // FCM (sharp transition)
    MULTI_CONTACT    // MCM (gradual miscibility)
};
```

---

### Class: `NonDarcyFlow`

Models non-Darcy flow effects for high-velocity porous media flow.

#### Nested Structures

```cpp
struct Parameters {
    NonDarcyModel model;
    double forchheimer_beta;  // Inertial coefficient [1/m]
    double klinkenberg_b;     // Slip factor [Pa]
    double ergun_a;           // Ergun viscous constant
    double ergun_b;           // Ergun inertial constant
};
```

#### Methods

| Method | Description |
|--------|-------------|
| `void setParameters(const Parameters& params)` | Set model parameters |
| `double calculateCorrection(v, rho, mu, K)` | Compute permeability correction factor (0 to 1) |
| `double klinkenbergPermeability(K_abs, p)` | Apparent permeability with Klinkenberg effect |
| `double reynoldsNumber(v, rho, mu, K)` | Pore-scale Reynolds number |
| `double effectiveVelocity(grad_p, rho, mu, K)` | Velocity accounting for non-Darcy effects |
| `void configure(const std::map<std::string,std::string>& config)` | Configure from key-value pairs |

---

### Class: `DualPorosityPermeability`

Models naturally fractured reservoirs with matrix-fracture exchange.

#### Nested Structures

```cpp
struct Parameters {
    DualPorosityModel model;
    double matrix_porosity;
    double fracture_porosity;
    double vug_porosity;            // For triple porosity
    double matrix_permeability;     // [m²]
    double fracture_permeability;   // [m²]
    double vug_permeability;        // [m²]
    double shape_factor;            // [1/m²]
    double transfer_coeff;          // [1/s]
    double matrix_fracture_transfer;
    double fracture_vug_transfer;
};
```

#### Methods

| Method | Description |
|--------|-------------|
| `double calculateTransferRate(p_m, p_f)` | Matrix-fracture transfer rate |
| `double effectivePorosity()` | Combined porosity |
| `double effectivePermeability()` | Effective permeability tensor (simplified) |
| `double storativityRatio(c_m, c_f)` | Warren-Root storativity ratio ω |
| `double interporosityFlowCoefficient(L)` | Dimensionless interporosity flow λ |

---

### Class: `NonIsothermalMultiphaseFlow`

Temperature effects on multiphase flow.

#### Methods

| Method | Description |
|--------|-------------|
| `double effectiveThermalConductivity(phi, Sw, So, Sg)` | Weighted average conductivity |
| `double effectiveHeatCapacity(phi, Sw, So, Sg, rho_s, rho_w, rho_o, rho_g)` | Volumetric heat capacity |
| `double viscosityAtTemperature(mu_ref, T_ref, T, E_a)` | Arrhenius viscosity model |
| `double viscousHeatingRate(mu, grad_p, K)` | Dissipation heat source |
| `double jouleThomsonTemperatureChange(JT, dp)` | Temperature change from JT effect |

---

### Class: `DynamicRelativePermeability`

Velocity and capillary-number dependent relative permeability.

#### Methods

| Method | Description |
|--------|-------------|
| `double capillaryNumber(mu, v, sigma)` | Compute Nc = μv/σ |
| `double desaturationCoefficient(Nc)` | Desaturation vs capillary number |
| `double dynamicResidualSaturation(Sor_static, D)` | Modified residual saturation |
| `double dynamicEndpoint(kr_max_static, D)` | Modified endpoint rel perm |
| `double modifyRelativePermeability(kr, Sw, mu, v, sigma)` | Apply dynamic modification |
| `double mobilityRatio(kr_w, mu_w, kr_o, mu_o)` | Compute mobility ratio M |
| `double fractionalFlowWithFingering(kr_w, mu_w, kr_o, mu_o)` | Fractional flow with fingering |

---

### Class: `MiscibleFlow`

Miscible displacement and mixing models.

#### Methods

| Method | Description |
|--------|-------------|
| `double miscibilityFraction(double pressure)` | Fraction miscible (0 to 1) |
| `double effectiveViscosity(mu_o, mu_s, S_s, f_misc)` | Todd-Longstaff mixed viscosity |
| `double effectiveDensity(rho_o, rho_s, S_s, f_misc)` | Mixed density |
| `void dispersionTensor(v, D_mol, D_eff[3][3])` | Mechanical dispersion tensor |

---

### Configuration Structure

```cpp
struct HighFidelityFluidFlowConfig {
    bool enable_non_darcy;
    bool enable_dual_porosity;
    bool enable_non_isothermal;
    bool enable_dynamic_relperm;
    bool enable_miscible;
    
    NonDarcyFlow::Parameters non_darcy_params;
    DualPorosityPermeability::Parameters dual_porosity_params;
    NonIsothermalMultiphaseFlow::Parameters non_isothermal_params;
    DynamicRelativePermeability::Parameters dynamic_relperm_params;
    MiscibleFlow::Parameters miscible_params;
    
    bool validate(std::string& error_msg) const;
    void parseConfig(const std::map<std::string,std::string>& config);
};
```

---

## HighFidelityGeomechanics.hpp

### Namespace: `FSRM::HighFidelity`

### Enumerations

#### `HyperelasticModel`
```cpp
enum class HyperelasticModel {
    NEO_HOOKEAN,    // W = μ/2(I₁-3) + κ/2(J-1)²
    MOONEY_RIVLIN,  // W = C₁₀(I₁-3) + C₀₁(I₂-3)
    OGDEN,          // Ogden strain energy
    ARRUDA_BOYCE,   // Eight-chain network model
    GENT            // Gent limiting chain extensibility
};
```

#### `FinitePlasticityModel`
```cpp
enum class FinitePlasticityModel {
    MULTIPLICATIVE_J2,     // F = Fᵉ·Fᵖ with J2 yield
    MULTIPLICATIVE_DP,     // Drucker-Prager yield
    LOGARITHMIC_STRAIN     // Logarithmic strain plasticity
};
```

#### `ViscoplasticModel`
```cpp
enum class ViscoplasticModel {
    PERZYNA,       // Overstress model
    DUVAUT_LIONS,  // Relaxation model
    RATE_POWER     // Power law rate sensitivity
};
```

#### `CreepType`
```cpp
enum class CreepType {
    POWER_LAW,  // ε̇ = A·σⁿ·exp(-Q/RT)
    NORTON,     // ε̇ = A·σⁿ
    BURGERS,    // Maxwell + Kelvin-Voigt
    LEMAITRE    // Lemaitre creep
};
```

#### `DamageEvolutionModel`
```cpp
enum class DamageEvolutionModel {
    EXPONENTIAL,       // d = α(1 - exp(-β(κ-κ₀)))
    LINEAR_SOFTENING,  // Linear post-peak
    MAZARS,            // Mazars damage model
    LEMAITRE_DAMAGE    // Lemaitre coupled damage
};
```

---

### Class: `FiniteStrainMechanics`

Hyperelastic constitutive models for large deformations.

#### Methods

| Method | Description |
|--------|-------------|
| `void computeStress(F[3][3], stress[3][3])` | Cauchy stress from deformation gradient |
| `double computeStrainEnergy(F[3][3])` | Strain energy density |
| `void computeMaterialTangent(F[3][3], C[3][3][3][3])` | Material tangent tensor |

---

### Class: `FiniteStrainPlasticity`

Multiplicative finite strain plasticity.

#### Nested Structures

```cpp
struct State {
    double plastic_deformation_gradient[3][3];
    double equivalent_plastic_strain;
    double back_stress[6];  // Kinematic hardening
};
```

#### Methods

| Method | Description |
|--------|-------------|
| `void computeStress(F[3][3], State& state, stress[3][3])` | Stress with return mapping |
| `double yieldFunction(stress[6], eps_p_eq)` | Yield surface value |
| `double currentYieldStress(eps_p_eq)` | Hardened yield stress |

---

### Class: `Viscoplasticity`

Rate-dependent plasticity models.

#### Methods

| Method | Description |
|--------|-------------|
| `double computeOverstress(stress[6], eps_p_eq)` | Overstress (stress - yield) |
| `void computePlasticStrainRate(stress[6], eps_p_eq, eps_dot[6])` | Plastic strain rate tensor |
| `double equivalentPlasticStrainRate(stress[6], eps_p_eq)` | Scalar plastic strain rate |

---

### Class: `CreepModel`

Time-dependent creep deformation.

#### Methods

| Method | Description |
|--------|-------------|
| `double creepStrainRate(double stress)` | Steady-state creep rate |
| `double burgersCreepStrain(stress, time)` | Total Burgers creep strain |
| `void burgersCreepComponents(stress, t, elastic, transient, steady)` | Decomposed components |

---

### Class: `Hypoplasticity`

Non-linear model for granular materials.

#### Methods

| Method | Description |
|--------|-------------|
| `void computeStressRate(stress[6], D[6], e, stress_rate[6])` | Objective stress rate |
| `double criticalVoidRatio(double p)` | Critical state void ratio |
| `double minVoidRatio(double p)` | Minimum void ratio |
| `double maxVoidRatio(double p)` | Maximum void ratio |
| `double relativeDensity(double e, double p)` | Relative density index |
| `double barotropyFactor(double p, double e)` | Barotropy function |
| `double pyknotopyFactor(double e, double p)` | Pyknotropy function |

---

### Class: `GradientEnhancedDamage`

Regularized continuum damage with internal length.

#### Methods

| Method | Description |
|--------|-------------|
| `double computeDamage(double kappa)` | Damage from history variable |
| `double equivalentStrain(strain[6])` | Equivalent strain measure |
| `void effectiveStress(stress[6], d, stress_eff[6])` | Stress / (1-d) |
| `double gradientParameter()` | Returns l² for Helmholtz equation |

---

## AdvancedCoupledPhysics.hpp

### Namespace: `FSRM::AdvancedCoupled`

### Enumerations

#### `BiotFormulation`
```cpp
enum class BiotFormulation {
    U_P,     // Displacement-pressure
    U_W,     // Displacement-relative fluid displacement
    U_P_W    // Three-field formulation
};
```

#### `RetentionModel`
```cpp
enum class RetentionModel {
    VAN_GENUCHTEN,  // Se = [1 + (α·ψ)ⁿ]^(-m)
    BROOKS_COREY,   // Se = (ψ_ae/ψ)^λ
    LINEAR          // Simple linear
};
```

#### `BishopChiModel`
```cpp
enum class BishopChiModel {
    SATURATION,     // χ = S_w
    KHALILI,        // χ = (S_e)^κ
    EFFECTIVE       // Based on effective stress
};
```

#### `MultiScaleMethod`
```cpp
enum class MultiScaleMethod {
    FE_SQUARED,              // FE² computational homogenization
    HMM,                     // Heterogeneous Multiscale Method
    ASYMPTOTIC_HOMOGENIZATION // Classical asymptotic
};
```

---

### Class: `FullBiotDynamics`

Complete Biot poroelasticity with inertia.

#### Methods

| Method | Description |
|--------|-------------|
| `double biotCoefficient()` | Biot-Willis coefficient α |
| `double biotModulus()` | Biot modulus M |
| `double skemptonCoefficient()` | Skempton's B coefficient |
| `double drainedBulkModulus()` | Drained K |
| `double undrainedBulkModulus()` | Undrained Kᵤ = K + α²M |
| `double pwaveVelocityDrained()` | Drained P-wave velocity |
| `double pwaveVelocityUndrained()` | Undrained P-wave velocity |
| `double swaveVelocity()` | S-wave velocity |
| `double slowPwaveVelocity()` | Biot slow wave velocity |
| `double characteristicFrequency()` | Biot characteristic frequency |
| `bool isHighFrequencyRegime(double f)` | Check frequency regime |
| `void effectiveStress(sigma[6], p, sigma_eff[6])` | Terzaghi effective stress |
| `double fluidContent(eps_vol, p)` | Fluid content ζ |
| `void massMatrixCoefficients(M_uu, M_ww, M_uw)` | Dynamic mass terms |
| `double viscousDampingCoefficient()` | Viscous damping |

---

### Class: `THMCoupling`

Thermo-Hydro-Mechanical coupling.

#### Methods

| Method | Description |
|--------|-------------|
| `double thermalStress(double dT)` | Thermal stress for constrained heating |
| `double thermalExpansionStrain(double dT)` | Thermal strain |
| `double thermalPressurization(dT, K)` | Pore pressure change from heating |
| `double effectiveThermalExpansion()` | Combined solid-fluid expansion |
| `double thermalDiffusivity(double rho)` | κ = k/(ρc) |
| `double effectiveThermalConductivity(phi, k_s, k_f)` | Effective k |
| `double consolidationCoefficient(K, mu, K_bulk, phi)` | cv |
| `double characteristicTime(L, cv)` | Diffusion time L²/cv |

---

### Class: `THMCCoupling`

Thermo-Hydro-Mechanical-Chemical coupling.

#### Methods

| Method | Description |
|--------|-------------|
| `double chemicalStrain(concentrations)` | Strain from chemical change |
| `double chemicalStrainRate(concentration_rates)` | Strain rate |
| `double reactionHeat(rates, enthalpies)` | Heat from reactions |
| `double porosityChangeRate(diss, precip, V_mol)` | dφ/dt |
| `double effectiveDiffusivity(D, phi, tau)` | D_eff = φD/τ |
| `double temperatureDependentRate(k0, Ea, T, T_ref)` | Arrhenius rate |
| `double chemicalOsmosisFlux(grad_c, L_D)` | Osmotic flux |

---

### Class: `UnsaturatedFlowCoupling`

Variably saturated porous media.

#### Methods

| Method | Description |
|--------|-------------|
| `double waterSaturation(double suction)` | S_w from retention curve |
| `double relativePermeability(double S_e)` | kr(S_e) |
| `double bishopParameter(double S_w)` | χ for effective stress |
| `double effectiveStress(sigma, p_w, p_a, chi)` | Bishop effective stress |
| `double moistureCapacity(double suction)` | dS_w/dψ |

---

### Class: `MultiScaleCoupling`

Homogenization and multi-scale methods.

#### Methods

| Method | Description |
|--------|-------------|
| `double voigtAverage(moduli, fractions)` | Upper bound |
| `double reussAverage(moduli, fractions)` | Lower bound |
| `void hashinShtrikmanBounds(...)` | HS bounds |
| `void moriTanaka(K_m, mu_m, K_i, mu_i, f_i, K, mu)` | MT estimate |
| `void selfConsistent(K1, mu1, K2, mu2, f1, K, mu)` | SC estimate |
| `void strainLocalization(macro[6], A[6][6], micro[6])` | Strain localization |

---

## HighFidelityNumerics.hpp

### Namespace: `FSRM::HighFidelity`

### Enumerations

#### `StabilizationMethod`
```cpp
enum class StabilizationMethod {
    SUPG,  // Streamline Upwind Petrov-Galerkin
    GLS,   // Galerkin Least Squares
    PSPG,  // Pressure Stabilization Petrov-Galerkin
    VMS    // Variational Multi-Scale
};
```

#### `TauFormula`
```cpp
enum class TauFormula {
    SHAKIB,          // Full Shakib formula
    TEZDUYAR,        // Simplified Tezduyar
    FRANCA_VALENTIN, // Original FV
    CODINA           // Codina formula
};
```

#### `TimeScheme`
```cpp
enum class TimeScheme {
    FORWARD_EULER,
    BACKWARD_EULER,
    CRANK_NICOLSON,
    BDF1, BDF2, BDF3, BDF4,
    RK4, SSPRK3,
    NEWMARK_BETA,
    HHT_ALPHA,
    GENERALIZED_ALPHA,
    ARKIMEX  // Additive Runge-Kutta IMEX
};
```

#### `PreconditionerType`
```cpp
enum class PreconditionerType {
    BLOCK_JACOBI,
    BLOCK_GAUSS_SEIDEL,
    FIELDSPLIT_ADDITIVE,
    FIELDSPLIT_MULTIPLICATIVE,
    FIELDSPLIT_SCHUR,
    GAMG,
    HYPRE_BOOMERAMG,
    CUSTOM
};
```

---

### Class: `SUPGStabilization`

SUPG/GLS/PSPG stabilization for advection-diffusion.

#### Nested Structures

```cpp
struct Parameters {
    StabilizationMethod method;
    TauFormula tau_formula;
    PetscReal tau_multiplier;
    bool include_transient;
    bool include_reaction;
    bool shock_capturing;
    PetscReal shock_capture_coeff;
};
```

#### Methods

| Method | Signature |
|--------|-----------|
| `calculateTau` | `(velocity[3], κ, reaction, h, dt) → τ` |
| `pecletNumber` | `(v_mag, h, κ) → Pe` |
| `optimalDiffusivity` | `(v_mag, h, κ) → κ_art` |
| `addStabilization` | `(v[3], ∇φ[3], R, τ, term*)` |
| `pspgTerm` | `(∇p[3], ∇q[3], τ) → scalar` |
| `shockCapturingDiffusivity` | `(∇u[3], u, h, R) → κ_dc` |

---

### Class: `PETScTimeIntegration`

Wrapper for PETSc TS module.

#### Nested Structures

```cpp
struct Parameters {
    TimeScheme scheme;
    bool adaptive;
    PetscReal atol, rtol;
    PetscReal max_dt, min_dt;
    PetscReal theta;        // For theta method
    PetscReal rho_infinity; // For generalized-alpha
    PetscReal beta, gamma;  // For Newmark
};
```

#### Methods

| Method | Description |
|--------|-------------|
| `void create(DM dm)` | Create TS from DM |
| `void setType(TimeScheme scheme)` | Set time stepping type |
| `void setAdaptivity(bool, atol, rtol)` | Configure adaptivity |
| `void setRHSFunction(func, ctx)` | Set F in u̇ = F(t,u) |
| `void setIFunction(func, ctx)` | Set F in F(t,u,u̇) = 0 |
| `void setIJacobian(func, ctx)` | Set Jacobian |
| `void solve(Vec u, t_final)` | Integrate to t_final |
| `int getOrder()` | Temporal order |
| `bool isAStable()` | A-stability check |
| `bool isLStable()` | L-stability check |

---

### Class: `PAdaptivityManager`

Polynomial order adaptation.

#### Methods

| Method | Description |
|--------|-------------|
| `int decideOrderChange(smoothness, error, p)` | -1, 0, or +1 |
| `double smoothnessIndicator(Vec u, DM dm)` | Estimate regularity |
| `void changeOrder(dm, u, p_new, u_interp)` | Interpolate to new order |

---

### Class: `PETScPhysicsPreconditioner`

Physics-based preconditioning wrapper.

#### Methods

| Method | Description |
|--------|-------------|
| `void create(KSP ksp)` | Create PC from KSP |
| `void setType(PreconditionerType type)` | Set PC type |
| `void setSchurFactorization(fact, pre)` | Configure Schur complement |
| `void setupFieldSplit(dm, field_names)` | Field-split from DM fields |
| `void setupMultigrid(levels, smoother)` | GAMG configuration |

---

## Factory Functions

All modules provide factory functions for convenient object creation:

```cpp
// Fluid Flow
std::unique_ptr<NonDarcyFlow> createNonDarcyFlow(const std::map<std::string,std::string>& config);
std::unique_ptr<DualPorosityPermeability> createDualPorosityModel(const std::map<std::string,std::string>& config);
std::unique_ptr<MiscibleFlow> createMiscibleFlow(const std::map<std::string,std::string>& config);

// Geomechanics
std::unique_ptr<FiniteStrainMechanics> createFiniteStrainMechanics(const std::map<std::string,std::string>& config);
std::unique_ptr<Viscoplasticity> createViscoplasticity(const std::map<std::string,std::string>& config);
std::unique_ptr<CreepModel> createCreepModel(const std::map<std::string,std::string>& config);
std::unique_ptr<GradientEnhancedDamage> createGradientDamage(const std::map<std::string,std::string>& config);

// Coupled Physics
std::unique_ptr<FullBiotDynamics> createFullBiotDynamics(const std::map<std::string,std::string>& config);
std::unique_ptr<THMCoupling> createTHMCoupling(const std::map<std::string,std::string>& config);
std::unique_ptr<THMCCoupling> createTHMCCoupling(const std::map<std::string,std::string>& config);
std::unique_ptr<UnsaturatedFlowCoupling> createUnsaturatedFlowCoupling(const std::map<std::string,std::string>& config);
std::unique_ptr<MultiScaleCoupling> createMultiScaleCoupling(const std::map<std::string,std::string>& config);

// Numerics
std::unique_ptr<SUPGStabilization> createSUPGStabilization(const std::map<std::string,std::string>& config);
std::unique_ptr<PETScTimeIntegration> createTimeIntegration(MPI_Comm comm, const std::map<std::string,std::string>& config);
std::unique_ptr<PETScPhysicsPreconditioner> createPhysicsPreconditioner(MPI_Comm comm, const std::map<std::string,std::string>& config);
```
