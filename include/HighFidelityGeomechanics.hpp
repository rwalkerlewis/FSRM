#ifndef HIGH_FIDELITY_GEOMECHANICS_HPP
#define HIGH_FIDELITY_GEOMECHANICS_HPP

/**
 * @file HighFidelityGeomechanics.hpp
 * @brief High-fidelity geomechanics models for advanced stress analysis
 * 
 * This file provides advanced geomechanics formulations that extend the basic
 * capabilities in MaterialModel.hpp and PlasticityModel.hpp. These are optional
 * high-fidelity extensions for complex geomechanical analysis.
 * 
 * Features:
 * - Finite strain / large deformation formulations
 * - Rate-dependent plasticity (viscoplasticity)
 * - Creep models (power law, Norton, Lemaitre)
 * - Hypoplasticity for granular materials
 * - Gradient-enhanced damage (non-local)
 * - Cosserat / micropolar elasticity
 * - Hypoelastic-plastic formulations
 * - Advanced yield surfaces (Lade, Matsuoka-Nakai)
 * 
 * All models integrate with existing MaterialModel and PlasticityModel
 * classes without replacing them.
 * 
 * @note These are high-fidelity options - basic simulations should use
 *       the simpler models in MaterialModel.hpp
 */

#include "MaterialModel.hpp"
#include "PlasticityModel.hpp"
#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <map>
#include <cmath>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Tensor Operations for Finite Strain
// =============================================================================

/**
 * @brief Deformation gradient tensor F = ∂x/∂X
 * 
 * Stores the 3x3 deformation gradient for finite strain analysis.
 * F = R·U (polar decomposition) where R is rotation and U is stretch.
 */
struct DeformationGradient {
    std::array<std::array<double, 3>, 3> F;
    
    DeformationGradient();
    DeformationGradient(const std::array<std::array<double, 3>, 3>& components);
    
    // Tensor operations
    double determinant() const;  // J = det(F)
    DeformationGradient inverse() const;
    DeformationGradient transpose() const;
    
    // Strain measures
    std::array<std::array<double, 3>, 3> rightCauchyGreen() const;   // C = F^T·F
    std::array<std::array<double, 3>, 3> leftCauchyGreen() const;    // B = F·F^T
    std::array<std::array<double, 3>, 3> greenLagrangeStrain() const; // E = 0.5(C - I)
    std::array<std::array<double, 3>, 3> eulerAlmansiStrain() const;  // e = 0.5(I - B^-1)
    std::array<std::array<double, 3>, 3> logarithmicStrain() const;   // ln(U)
    
    // Polar decomposition
    void polarDecomposition(std::array<std::array<double, 3>, 3>& R,
                           std::array<std::array<double, 3>, 3>& U) const;
    
    // Rate of deformation
    static std::array<std::array<double, 3>, 3> 
    velocityGradient(const DeformationGradient& F, const DeformationGradient& F_dot);
    
    static std::array<std::array<double, 3>, 3>
    rateOfDeformation(const std::array<std::array<double, 3>, 3>& L);  // D = sym(L)
    
    static std::array<std::array<double, 3>, 3>
    spinTensor(const std::array<std::array<double, 3>, 3>& L);  // W = skew(L)
};

/**
 * @brief Stress measures for finite strain
 */
struct FiniteStrainStress {
    std::array<std::array<double, 3>, 3> cauchy;     // True stress σ
    std::array<std::array<double, 3>, 3> kirchhoff;  // τ = J·σ
    std::array<std::array<double, 3>, 3> piola1;     // P = J·σ·F^-T (1st Piola-Kirchhoff)
    std::array<std::array<double, 3>, 3> piola2;     // S = F^-1·P (2nd Piola-Kirchhoff)
    
    // Convert between stress measures
    static void cauchyToKirchhoff(const std::array<std::array<double, 3>, 3>& sigma,
                                  double J,
                                  std::array<std::array<double, 3>, 3>& tau);
    
    static void cauchyToPiola1(const std::array<std::array<double, 3>, 3>& sigma,
                               const DeformationGradient& F,
                               std::array<std::array<double, 3>, 3>& P);
    
    static void cauchyToPiola2(const std::array<std::array<double, 3>, 3>& sigma,
                               const DeformationGradient& F,
                               std::array<std::array<double, 3>, 3>& S);
};


// =============================================================================
// Finite Strain Hyperelastic Models
// =============================================================================

/**
 * @brief Hyperelastic model type
 */
enum class HyperelasticModel {
    NEO_HOOKEAN,         ///< Neo-Hookean (simplest)
    MOONEY_RIVLIN,       ///< Mooney-Rivlin (rubber-like)
    OGDEN,               ///< Ogden (highly nonlinear)
    SAINT_VENANT,        ///< St. Venant-Kirchhoff (geometrically nonlinear)
    YEOH,                ///< Yeoh model
    ARRUDA_BOYCE,        ///< Arruda-Boyce (8-chain)
    HENCKY                ///< Hencky (logarithmic strain)
};

/**
 * @brief Finite strain hyperelastic material
 * 
 * Derives stress from a strain energy function W(F) or W(C).
 * Provides consistent tangent modulus for Newton iterations.
 * 
 * Key models:
 * - Neo-Hookean: W = μ/2·(I₁-3) + λ/2·(J-1)² 
 * - Mooney-Rivlin: W = C₁·(I₁-3) + C₂·(I₂-3)
 * - Ogden: W = Σ μₚ/αₚ·(λ₁^αₚ + λ₂^αₚ + λ₃^αₚ - 3)
 */
class HyperelasticMaterial {
public:
    struct Parameters {
        HyperelasticModel model = HyperelasticModel::NEO_HOOKEAN;
        
        // Lamé parameters (for Neo-Hookean)
        double lambda = 1e9;     ///< First Lamé parameter [Pa]
        double mu = 1e9;         ///< Shear modulus [Pa]
        
        // Mooney-Rivlin parameters
        double C1 = 0.5e9;
        double C2 = 0.1e9;
        
        // Ogden parameters (up to 3 terms)
        std::array<double, 3> mu_ogden = {0.0, 0.0, 0.0};
        std::array<double, 3> alpha_ogden = {2.0, -2.0, 1.0};
        
        // Bulk modulus for volumetric response
        double kappa = 2e9;      ///< Bulk modulus [Pa]
        
        // Compressibility options
        bool nearly_incompressible = false;
        double penalty = 1e10;   ///< Penalty for incompressibility
    };
    
    HyperelasticMaterial();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate strain energy density
     * 
     * @param F Deformation gradient
     * @return Strain energy density [J/m³]
     */
    double strainEnergy(const DeformationGradient& F) const;
    
    /**
     * @brief Calculate 2nd Piola-Kirchhoff stress
     * 
     * S = 2·∂W/∂C
     * 
     * @param F Deformation gradient
     * @return 2nd Piola-Kirchhoff stress tensor
     */
    std::array<std::array<double, 3>, 3> secondPiolaKirchhoff(const DeformationGradient& F) const;
    
    /**
     * @brief Calculate Cauchy stress
     * 
     * σ = (1/J)·F·S·F^T
     * 
     * @param F Deformation gradient
     * @return Cauchy stress tensor
     */
    std::array<std::array<double, 3>, 3> cauchyStress(const DeformationGradient& F) const;
    
    /**
     * @brief Calculate material tangent modulus (for 2nd PK stress)
     * 
     * C = 4·∂²W/∂C∂C
     * 
     * @param F Deformation gradient
     * @return 6x6 material tangent (Voigt notation)
     */
    std::array<std::array<double, 6>, 6> materialTangent(const DeformationGradient& F) const;
    
    /**
     * @brief Calculate spatial tangent modulus (for Cauchy stress)
     * 
     * Push-forward of material tangent to current configuration.
     * 
     * @param F Deformation gradient
     * @return 6x6 spatial tangent (Voigt notation)
     */
    std::array<std::array<double, 6>, 6> spatialTangent(const DeformationGradient& F) const;
    
    /**
     * @brief Check if deformation is valid (J > 0)
     */
    bool isValidDeformation(const DeformationGradient& F) const;
    
private:
    Parameters params_;
    
    // Model-specific implementations
    double W_NeoHookean(const DeformationGradient& F) const;
    double W_MooneyRivlin(const DeformationGradient& F) const;
    double W_Ogden(const DeformationGradient& F) const;
    double W_StVenant(const DeformationGradient& F) const;
    
    std::array<std::array<double, 3>, 3> S_NeoHookean(const DeformationGradient& F) const;
    std::array<std::array<double, 3>, 3> S_MooneyRivlin(const DeformationGradient& F) const;
    
    // Principal stretches and invariants
    void computePrincipalStretches(const DeformationGradient& F,
                                   std::array<double, 3>& lambda) const;
    void computeInvariants(const DeformationGradient& F,
                          double& I1, double& I2, double& I3) const;
};


// =============================================================================
// Finite Strain Elastoplasticity
// =============================================================================

/**
 * @brief Finite strain plasticity formulation type
 */
enum class FinitePlasticityType {
    MULTIPLICATIVE,      ///< F = Fe·Fp decomposition
    ADDITIVE_HENCKY,     ///< Additive in logarithmic strain
    HYPOELASTIC          ///< Objective stress rate formulation
};

/**
 * @brief Objective stress rate type for hypoelastic-plastic
 */
enum class ObjectiveRate {
    JAUMANN,             ///< Jaumann rate (simplest)
    GREEN_NAGHDI,        ///< Green-Naghdi rate
    TRUESDELL,           ///< Truesdell rate
    OLDROYD,             ///< Oldroyd rate
    LOG_SPIN             ///< Logarithmic rate (best for large rotations)
};

/**
 * @brief Finite strain elastoplastic material
 * 
 * Implements large deformation plasticity using multiplicative
 * decomposition of the deformation gradient:
 * 
 * F = Fe·Fp
 * 
 * where Fe is elastic and Fp is plastic deformation.
 * 
 * The elastic response is computed from Fe, and plastic flow
 * updates Fp based on yield condition and flow rule.
 */
class FiniteStrainPlasticity {
public:
    struct Parameters {
        FinitePlasticityType formulation = FinitePlasticityType::MULTIPLICATIVE;
        ObjectiveRate objective_rate = ObjectiveRate::JAUMANN;
        
        // Yield criterion (uses existing PlasticityModel types)
        YieldCriterion yield_criterion = YieldCriterion::VON_MISES;
        
        // Initial yield stress
        double yield_stress = 100e6;    ///< σ_y [Pa]
        
        // Hardening type
        HardeningLaw hardening = HardeningLaw::ISOTROPIC_LINEAR;
        double hardening_modulus = 1e9;  ///< H [Pa]
        double saturation_stress = 500e6; ///< σ_sat [Pa]
        
        // Kinematic hardening (Bauschinger effect)
        bool kinematic_hardening = false;
        double kinematic_modulus = 5e8;  ///< Back-stress modulus [Pa]
        double kinematic_saturation = 2e8; ///< Back-stress saturation [Pa]
        
        // For pressure-dependent (Drucker-Prager in finite strain)
        double friction_angle = 30.0;    ///< φ [degrees]
        double cohesion = 1e6;           ///< c [Pa]
    };
    
    struct State {
        DeformationGradient Fp;                    ///< Plastic deformation gradient
        double accumulated_plastic_strain = 0.0;   ///< ε̄ₚ
        std::array<std::array<double, 3>, 3> back_stress; ///< Kinematic hardening back stress
        bool is_plastic = false;
    };
    
    FiniteStrainPlasticity();
    
    void setParameters(const Parameters& params);
    void setHyperelasticModel(std::shared_ptr<HyperelasticMaterial> elastic);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Initialize plastic state
     */
    void initializeState(State& state) const;
    
    /**
     * @brief Compute stress and update plastic state
     * 
     * Main stress integration routine using return mapping
     * in principal stretch space (Simo-style).
     * 
     * @param F_new Current deformation gradient
     * @param[in,out] state Plastic state variables
     * @return Cauchy stress
     */
    std::array<std::array<double, 3>, 3> computeStress(const DeformationGradient& F_new,
                                                        State& state) const;
    
    /**
     * @brief Compute consistent tangent modulus
     * 
     * Includes algorithmic contribution from return mapping.
     */
    std::array<std::array<double, 6>, 6> consistentTangent(const DeformationGradient& F,
                                                           const State& state) const;
    
    /**
     * @brief Evaluate yield function
     * 
     * f = ||dev(τ)|| - √(2/3)·(σ_y + H·ε̄ₚ)
     * 
     * where τ is Kirchhoff stress.
     */
    double yieldFunction(const std::array<std::array<double, 3>, 3>& kirchhoff_stress,
                        double accumulated_strain) const;
    
    /**
     * @brief Check if state is yielding
     */
    bool isYielding(const std::array<std::array<double, 3>, 3>& stress, double eps_p) const;
    
private:
    Parameters params_;
    std::shared_ptr<HyperelasticMaterial> elastic_;
    
    // Return mapping in principal space
    void returnMapping(const std::array<double, 3>& trial_principal_tau,
                      std::array<double, 3>& principal_tau,
                      double& delta_gamma,
                      double& eps_p_new,
                      double eps_p_old) const;
    
    // Objective stress rate computation
    std::array<std::array<double, 3>, 3> objectiveRate(
        const std::array<std::array<double, 3>, 3>& sigma,
        const std::array<std::array<double, 3>, 3>& D,
        const std::array<std::array<double, 3>, 3>& W) const;
};


// =============================================================================
// Rate-Dependent Plasticity (Viscoplasticity)
// =============================================================================

/**
 * @brief Viscoplasticity model type
 */
enum class ViscoplasticModel {
    PERZYNA,             ///< Perzyna (overstress) model
    DUVAUT_LIONS,        ///< Duvaut-Lions regularization
    CONSISTENCY,         ///< Consistency model (Wang)
    UNIFIED_CHABOCHE     ///< Chaboche unified model
};

/**
 * @brief Overstress function type for Perzyna model
 */
enum class OverstressFunction {
    LINEAR,              ///< φ(f) = f
    POWER_LAW,           ///< φ(f) = (f/σ₀)^n
    EXPONENTIAL,         ///< φ(f) = exp(f/σ₀) - 1
    SINH                 ///< φ(f) = sinh(f/σ₀)
};

/**
 * @brief Rate-dependent plasticity (viscoplasticity) model
 * 
 * Allows stress to temporarily exceed the yield surface,
 * with viscous relaxation back to the yield surface.
 * 
 * Perzyna model:
 * ε̇ₚ = γ·〈f〉/η · ∂f/∂σ
 * 
 * where 〈f〉 = max(0, f) is the overstress (stress above yield),
 * γ is flow direction, and η is viscosity.
 * 
 * This is important for:
 * - High strain rate loading
 * - Impact and blast problems
 * - Creep-plasticity interaction
 * - Regularization of strain localization
 */
class Viscoplasticity {
public:
    struct Parameters {
        ViscoplasticModel model = ViscoplasticModel::PERZYNA;
        OverstressFunction overstress_function = OverstressFunction::POWER_LAW;
        
        // Perzyna parameters
        double fluidity = 1e-6;          ///< γ [1/(Pa·s)]
        double reference_stress = 100e6; ///< σ₀ [Pa] (for normalization)
        double rate_exponent = 3.0;      ///< n (power law exponent)
        double viscosity = 1e12;         ///< η [Pa·s]
        
        // Yield surface (from PlasticityModel)
        YieldCriterion yield_criterion = YieldCriterion::VON_MISES;
        double yield_stress = 100e6;
        double hardening_modulus = 1e9;
        
        // Thermal effects
        bool thermal_activation = false;
        double activation_energy = 3e5;  ///< Q [J/mol]
        double reference_temperature = 293.15; ///< T₀ [K]
    };
    
    struct State {
        std::array<double, 6> plastic_strain;       ///< εₚ (Voigt)
        double accumulated_plastic_strain = 0.0;    ///< ε̄ₚ
        std::array<double, 6> back_stress;          ///< α (Voigt)
        double temperature = 293.15;                ///< T [K]
    };
    
    Viscoplasticity();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Initialize state variables
     */
    void initializeState(State& state, double T0 = 293.15) const;
    
    /**
     * @brief Compute stress with rate effects
     * 
     * Performs implicit integration of viscoplastic equations.
     * 
     * @param strain_increment Strain increment (Voigt notation)
     * @param dt Time step [s]
     * @param[in,out] stress Current stress (Voigt), updated in place
     * @param[in,out] state Viscoplastic state variables
     * @param C Elastic stiffness matrix
     */
    void integrateStress(const std::array<double, 6>& strain_increment,
                        double dt,
                        std::array<double, 6>& stress,
                        State& state,
                        const std::array<std::array<double, 6>, 6>& C) const;
    
    /**
     * @brief Calculate plastic strain rate
     * 
     * ε̇ₚ = γ·φ(f)·∂g/∂σ
     */
    std::array<double, 6> plasticStrainRate(const std::array<double, 6>& stress,
                                            const State& state) const;
    
    /**
     * @brief Evaluate overstress function
     * 
     * φ(f) for different models:
     * - Linear: f/σ₀
     * - Power: (f/σ₀)^n
     * - Exponential: exp(f/σ₀) - 1
     */
    double overstressFunction(double f) const;
    double overstressFunctionDerivative(double f) const;
    
    /**
     * @brief Calculate rate-dependent yield stress
     * 
     * For rate-dependent yield: σ_y(ε̇) = σ_y0·(1 + (ε̇/ε̇₀)^m)
     */
    double rateDependentYield(double strain_rate) const;
    
    /**
     * @brief Calculate consistent tangent with viscous effects
     */
    std::array<std::array<double, 6>, 6> consistentTangent(const std::array<double, 6>& stress,
                                                           const State& state,
                                                           double dt,
                                                           const std::array<std::array<double, 6>, 6>& C) const;
    
private:
    Parameters params_;
    
    double yieldFunction(const std::array<double, 6>& stress, double hardening) const;
    std::array<double, 6> flowDirection(const std::array<double, 6>& stress) const;
    double thermalFactor(double T) const;  // Arrhenius-type
};


// =============================================================================
// Creep Models
// =============================================================================

/**
 * @brief Creep model type
 */
enum class CreepModel {
    NONE,
    POWER_LAW,           ///< Norton/Bailey power law
    EXPONENTIAL,         ///< Exponential creep
    SINH,                ///< Hyperbolic sine (Garofalo)
    NORTON_HOFF,         ///< Norton-Hoff unified
    LEMAITRE,            ///< Lemaitre coupled damage-creep
    THETA_PROJECTION,    ///< θ-projection for long-term
    ANDRADE,             ///< Andrade transient creep
    BURGERS              ///< Burgers model (Maxwell + Kelvin)
};

/**
 * @brief Creep stage
 */
enum class CreepStage {
    PRIMARY,             ///< Transient (decreasing rate)
    SECONDARY,           ///< Steady-state (constant rate)
    TERTIARY             ///< Accelerating (damage)
};

/**
 * @brief Creep model for time-dependent deformation
 * 
 * Models long-term deformation under constant stress:
 * 
 * Norton power law: ε̇ = A·σⁿ·exp(-Q/RT)
 * 
 * where:
 * - A is material constant
 * - n is stress exponent (typically 3-5 for rocks)
 * - Q is activation energy
 * - R is gas constant
 * - T is temperature
 * 
 * Important for:
 * - Salt creep in cavern storage
 * - Long-term subsidence
 * - Clay squeezing
 * - Tunnel closure
 */
class CreepModel_Impl {
public:
    struct Parameters {
        CreepModel model = CreepModel::POWER_LAW;
        
        // Power law (Norton) parameters
        double A = 1e-30;            ///< Creep coefficient [1/Pa^n/s]
        double n = 4.0;              ///< Stress exponent [-]
        double activation_energy = 54e3; ///< Q [J/mol]
        double gas_constant = 8.314; ///< R [J/(mol·K)]
        
        // Reference values
        double reference_stress = 10e6;   ///< σ₀ [Pa]
        double reference_temperature = 293.15; ///< T₀ [K]
        double reference_strain_rate = 1e-10; ///< ε̇₀ [1/s]
        
        // Primary creep (transient)
        bool include_primary = true;
        double primary_exponent = 0.33;   ///< m (time exponent)
        double primary_coefficient = 1e-5; ///< B
        
        // Tertiary creep (damage)
        bool include_tertiary = false;
        double damage_rate = 1e-12;       ///< Damage evolution rate
        double critical_damage = 0.9;     ///< D_c (failure)
        
        // Theta-projection parameters
        double theta1 = 0.1;
        double theta2 = 1e-8;
        double theta3 = 0.5;
        double theta4 = 1e-10;
    };
    
    struct State {
        double creep_strain = 0.0;              ///< Total creep strain
        std::array<double, 6> creep_strain_tensor; ///< Creep strain tensor
        double time = 0.0;                       ///< Total creep time
        double damage = 0.0;                     ///< Creep damage (0-1)
        CreepStage stage = CreepStage::PRIMARY;
        
        // For Burgers model
        double maxwell_strain = 0.0;
        double kelvin_strain = 0.0;
    };
    
    CreepModel_Impl();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Initialize creep state
     */
    void initializeState(State& state) const;
    
    /**
     * @brief Calculate creep strain rate
     * 
     * Returns deviatoric creep strain rate tensor.
     * 
     * @param stress Stress tensor (Voigt)
     * @param T Temperature [K]
     * @param state Current creep state
     * @return Creep strain rate tensor (Voigt)
     */
    std::array<double, 6> creepStrainRate(const std::array<double, 6>& stress,
                                          double T,
                                          const State& state) const;
    
    /**
     * @brief Update creep state over time step
     * 
     * @param stress Applied stress (Voigt)
     * @param T Temperature [K]
     * @param dt Time step [s]
     * @param[in,out] state Creep state
     */
    void updateState(const std::array<double, 6>& stress,
                    double T,
                    double dt,
                    State& state) const;
    
    /**
     * @brief Calculate effective creep strain rate (scalar)
     * 
     * ε̇_eff = A·σ_eq^n·exp(-Q/RT)
     */
    double effectiveCreepRate(double sigma_eq, double T) const;
    
    /**
     * @brief Calculate time to rupture (tertiary creep)
     * 
     * Monkman-Grant relation: t_r = C·ε̇_min^(-m)
     */
    double timeToRupture(double sigma_eq, double T) const;
    
    /**
     * @brief Check creep stage
     */
    CreepStage determineStage(const State& state, double sigma_eq) const;
    
    /**
     * @brief Get viscosity for viscous analogy
     * 
     * η_creep = σ_eq / (3·ε̇_eq)
     */
    double getCreepViscosity(double sigma_eq, double T) const;
    
private:
    Parameters params_;
    
    // Individual creep components
    double primaryCreepRate(double sigma_eq, double T, double time) const;
    double secondaryCreepRate(double sigma_eq, double T) const;
    double tertiaryCreepRate(double sigma_eq, double T, double damage) const;
    
    // Theta-projection
    double thetaProjection(double sigma_eq, double T, double time) const;
    
    // Burgers model
    void burgersUpdate(double sigma_eq, double dt, double& maxwell, double& kelvin) const;
};


// =============================================================================
// Hypoplasticity (Granular Materials)
// =============================================================================

/**
 * @brief Hypoplasticity model type
 */
enum class HypoplasticModel {
    VON_WOLFFERSDORFF,   ///< Von Wolffersdorff (most common)
    GUDEHUS_BAUER,       ///< Gudehus-Bauer basic
    HERLE_GUDEHUS,       ///< With intergranular strain
    NIEMUNIS,            ///< Extended with viscosity
    MASIN_CLAY           ///< For clays
};

/**
 * @brief Hypoplastic material model
 * 
 * Non-linear constitutive model without explicit yield surface.
 * Stress rate depends on both stress and strain rate:
 * 
 * σ̊ = L(σ,e):D + N(σ,e)·||D||
 * 
 * Particularly suited for:
 * - Granular materials (sand, gravel)
 * - Large strain cyclic loading
 * - Stress-dilatancy behavior
 * - Void ratio dependency
 * 
 * Key features:
 * - No yield surface (continuous plasticity)
 * - Void ratio (density) dependency
 * - Automatic stress-dilatancy coupling
 * - Critical state mechanics built-in
 */
class HypoplasticMaterial {
public:
    struct Parameters {
        HypoplasticModel model = HypoplasticModel::VON_WOLFFERSDORFF;
        
        // Critical state parameters
        double phi_c = 33.0;         ///< Critical state friction angle [deg]
        double e_c0 = 0.95;          ///< Critical void ratio at p=0
        double e_d0 = 0.55;          ///< Minimum void ratio at p=0
        double e_i0 = 1.15;          ///< Maximum void ratio at p=0
        double h_s = 2e9;            ///< Granulate hardness [Pa]
        double n_h = 0.25;           ///< Exponent for hardness
        double alpha_h = 0.13;       ///< Pyknotropy exponent
        double beta_h = 1.0;         ///< Barotropic exponent
        
        // Intergranular strain (small strain stiffness)
        bool use_intergranular = true;
        double R = 1e-4;             ///< Intergranular strain range
        double m_R = 5.0;            ///< Interpolation exponent
        double m_T = 2.0;            ///< Transition factor
        double beta_r = 0.2;         ///< Rate factor
        double chi = 6.0;            ///< Stiffness multiplier
        
        // Rate effects
        bool rate_dependent = false;
        double reference_rate = 1e-5; ///< ε̇₀ [1/s]
        double rate_exponent = 0.1;   ///< Rate sensitivity
    };
    
    struct State {
        double void_ratio = 0.7;                   ///< e
        std::array<double, 6> intergranular_strain; ///< δ (Voigt)
        double rho = 0.5;                          ///< ||δ||/R
        double temperature = 293.15;
    };
    
    HypoplasticMaterial();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Initialize state
     */
    void initializeState(State& state, double e0 = 0.7) const;
    
    /**
     * @brief Compute stress rate from strain rate
     * 
     * σ̊ = fs·fb·fe·(L:D + N·||D||)
     * 
     * @param stress Current stress (Voigt)
     * @param strain_rate Strain rate (Voigt)
     * @param state Hypoplastic state
     * @return Stress rate (Voigt)
     */
    std::array<double, 6> stressRate(const std::array<double, 6>& stress,
                                     const std::array<double, 6>& strain_rate,
                                     const State& state) const;
    
    /**
     * @brief Integrate stress over finite strain increment
     * 
     * Uses sub-stepping for accuracy.
     */
    void integrateStress(const std::array<double, 6>& strain_increment,
                        std::array<double, 6>& stress,
                        State& state,
                        double dt) const;
    
    /**
     * @brief Calculate tangent stiffness matrix
     * 
     * For Newton iterations, linearization of hypoplastic equation.
     */
    std::array<std::array<double, 6>, 6> tangentStiffness(
        const std::array<double, 6>& stress,
        const std::array<double, 6>& strain_rate,
        const State& state) const;
    
    /**
     * @brief Calculate limit void ratios at pressure p
     */
    void limitVoidRatios(double p, double& e_c, double& e_d, double& e_i) const;
    
    /**
     * @brief Calculate relative density
     * 
     * D_r = (e_max - e) / (e_max - e_min)
     */
    double relativeDensity(double e, double p) const;
    
    /**
     * @brief Calculate density factor (pyknotropy)
     * 
     * f_e = (e_c/e)^α
     */
    double densityFactor(double e, double p) const;
    
    /**
     * @brief Calculate pressure factor (barotropy)
     * 
     * f_b = (h_s/σ)^n · (1 + e_i) · e_i / e · (3p/h_s)^(1-n)
     */
    double pressureFactor(double p, double e) const;
    
private:
    Parameters params_;
    
    // L and N tensors
    void computeLN(const std::array<double, 6>& stress,
                   const State& state,
                   std::array<std::array<double, 6>, 6>& L,
                   std::array<double, 6>& N) const;
    
    // Intergranular strain update
    void updateIntergranularStrain(const std::array<double, 6>& strain_rate,
                                   State& state,
                                   double dt) const;
    
    // Scalar factors
    double fStiffness(const State& state) const;  // Small strain stiffness factor
};


// =============================================================================
// Gradient-Enhanced Damage (Non-local Regularization)
// =============================================================================

/**
 * @brief Gradient damage formulation type
 */
enum class GradientDamageType {
    INTEGRAL,            ///< Integral non-local averaging
    IMPLICIT_GRADIENT,   ///< Implicit gradient (Helmholtz equation)
    EXPLICIT_GRADIENT,   ///< Explicit gradient terms
    MICROMORPHIC         ///< Micromorphic continuum
};

/**
 * @brief Gradient-enhanced damage model
 * 
 * Regularizes strain localization by introducing a characteristic length
 * scale. Prevents mesh-dependent results for softening materials.
 * 
 * Implicit gradient approach:
 * ε̄ - c·∇²ε̄ = ε_local
 * 
 * where ε̄ is the non-local strain and c = l²/2 determines the length scale.
 * 
 * Damage evolution:
 * Ḋ = f(ε̄) for ε̄ > κ (threshold)
 * 
 * Stress:
 * σ = (1 - D)·C:ε
 */
class GradientDamage {
public:
    struct Parameters {
        GradientDamageType type = GradientDamageType::IMPLICIT_GRADIENT;
        
        // Length scale
        double characteristic_length = 0.01;  ///< l [m]
        double gradient_parameter = 0.0001;   ///< c = l²/2 [m²]
        
        // Damage initiation
        double damage_threshold = 1e-4;       ///< κ₀ (strain threshold)
        double tensile_strength = 5e6;        ///< f_t [Pa]
        
        // Damage evolution
        double fracture_energy = 100.0;       ///< G_f [J/m²]
        double softening_exponent = 1.0;      ///< Exponential softening rate
        double max_damage = 0.99;             ///< D_max (numerical stability)
        
        // Regularization type
        bool use_viscous_regularization = false;
        double viscosity = 1e6;               ///< η_d [Pa·s]
        
        // Anisotropic damage
        bool anisotropic = false;
    };
    
    struct State {
        double damage = 0.0;                  ///< D (isotropic)
        std::array<double, 6> damage_tensor;  ///< Ω (anisotropic)
        double kappa = 0.0;                   ///< History variable (max non-local strain)
        double nonlocal_strain = 0.0;         ///< ε̄
        std::array<double, 3> grad_nonlocal;  ///< ∇ε̄
    };
    
    GradientDamage();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Initialize damage state
     */
    void initializeState(State& state) const;
    
    /**
     * @brief Calculate equivalent strain for damage
     * 
     * Various definitions:
     * - Mazars: ε_eq = √(Σ〈ε_i〉₊²)
     * - Modified von Mises
     * - Principal strain-based
     */
    double equivalentStrain(const std::array<double, 6>& strain) const;
    
    /**
     * @brief Solve non-local equation for ε̄
     * 
     * Implicit gradient: ε̄ - c·∇²ε̄ = ε_eq
     * 
     * Returns residual for FE assembly:
     * R = ε̄·v + c·∇ε̄·∇v - ε_eq·v
     */
    double nonlocalResidual(double nonlocal, const std::array<double, 3>& grad_nonlocal,
                           double local_strain,
                           double test_function,
                           const std::array<double, 3>& grad_test) const;
    
    /**
     * @brief Calculate damage from non-local strain
     * 
     * D = 1 - κ₀/κ · exp(-β·(κ - κ₀))
     */
    double calculateDamage(double nonlocal_strain, State& state) const;
    
    /**
     * @brief Calculate damaged stress
     * 
     * σ = (1 - D)·σ̃
     */
    std::array<double, 6> damagedStress(const std::array<double, 6>& effective_stress,
                                        const State& state) const;
    
    /**
     * @brief Calculate damaged tangent
     */
    std::array<std::array<double, 6>, 6> damagedTangent(
        const std::array<std::array<double, 6>, 6>& C,
        const std::array<double, 6>& stress,
        const std::array<double, 6>& strain,
        const State& state) const;
    
    /**
     * @brief Check if damage is localizing
     */
    bool isLocalizing(const State& state) const;
    
    /**
     * @brief Calculate dissipated energy density
     */
    double dissipatedEnergy(const State& state,
                           const std::array<double, 6>& stress,
                           const std::array<double, 6>& strain) const;
    
private:
    Parameters params_;
    
    double damageEvolutionLaw(double kappa) const;
    double damageDrivingForce(const std::array<double, 6>& strain,
                              const std::array<std::array<double, 6>, 6>& C) const;
};


// =============================================================================
// High-Fidelity Geomechanics Configuration
// =============================================================================

/**
 * @brief Configuration for high-fidelity geomechanics
 */
struct HighFidelityGeomechConfig {
    // Master enable
    bool enable_high_fidelity = false;
    
    // Finite strain
    bool enable_finite_strain = false;
    HyperelasticMaterial::Parameters hyperelastic_params;
    FiniteStrainPlasticity::Parameters finite_plasticity_params;
    
    // Viscoplasticity
    bool enable_viscoplasticity = false;
    Viscoplasticity::Parameters viscoplastic_params;
    
    // Creep
    bool enable_creep = false;
    CreepModel_Impl::Parameters creep_params;
    
    // Hypoplasticity
    bool enable_hypoplasticity = false;
    HypoplasticMaterial::Parameters hypoplastic_params;
    
    // Gradient damage
    bool enable_gradient_damage = false;
    GradientDamage::Parameters damage_params;
    
    // Parse from config
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Validate
    bool validate(std::string& error_msg) const;
};

/**
 * @brief Factory functions for high-fidelity geomechanics models
 */
std::unique_ptr<HyperelasticMaterial> createHyperelasticMaterial(
    const std::map<std::string, std::string>& config);

std::unique_ptr<FiniteStrainPlasticity> createFiniteStrainPlasticity(
    const std::map<std::string, std::string>& config);

std::unique_ptr<Viscoplasticity> createViscoplasticity(
    const std::map<std::string, std::string>& config);

std::unique_ptr<CreepModel_Impl> createCreepModel(
    const std::map<std::string, std::string>& config);

std::unique_ptr<HypoplasticMaterial> createHypoplasticMaterial(
    const std::map<std::string, std::string>& config);

std::unique_ptr<GradientDamage> createGradientDamage(
    const std::map<std::string, std::string>& config);

} // namespace HighFidelity
} // namespace FSRM

#endif // HIGH_FIDELITY_GEOMECHANICS_HPP
