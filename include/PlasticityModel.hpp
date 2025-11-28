#ifndef PLASTICITY_MODEL_HPP
#define PLASTICITY_MODEL_HPP

#include "MaterialModel.hpp"
#include <array>
#include <vector>
#include <map>
#include <string>
#include <memory>

namespace FSRM {

/**
 * @brief Plasticity model for off-fault deformation
 * 
 * Plasticity allows permanent deformation when stress exceeds a yield criterion.
 * Critical for modeling:
 * - Earthquake damage zones
 * - Fault gouge formation
 * - Off-fault deformation
 * - Energy dissipation
 * 
 * Implementation follows return mapping algorithm (Simo & Hughes, 1998).
 */

/**
 * @brief Yield criterion type
 */
enum class YieldCriterion {
    VON_MISES,          // J2 plasticity (pressure-independent)
    TRESCA,             // Maximum shear stress
    DRUCKER_PRAGER,     // Pressure-dependent (most common for rocks)
    MOHR_COULOMB,       // Classical soil mechanics
    CAM_CLAY,           // Critical state model
    HOEK_BROWN,         // Empirical rock model
    CAP_MODEL           // With hardening cap
};

/**
 * @brief Flow rule (relates plastic strain rate to stress)
 */
enum class FlowRule {
    ASSOCIATIVE,        // Plastic potential = yield function
    NON_ASSOCIATIVE     // Different plastic potential (dilation angle)
};

/**
 * @brief Hardening/softening behavior
 */
enum class HardeningLaw {
    NONE,               // Perfect plasticity
    ISOTROPIC_LINEAR,   // Linear isotropic hardening
    ISOTROPIC_EXPONENTIAL, // Exponential hardening
    KINEMATIC,          // Kinematic hardening (Bauschinger effect)
    COMBINED,           // Isotropic + kinematic
    STRAIN_SOFTENING    // Softening with accumulated strain
};

/**
 * @brief Drucker-Prager yield criterion
 * 
 * F = √J2 + α I1 - κ ≤ 0
 * 
 * where:
 * - J2 = second invariant of deviatoric stress
 * - I1 = first invariant of stress (mean stress)
 * - α, κ = material parameters related to friction angle and cohesion
 * 
 * This is the most common model for rocks and geological materials.
 */
class DruckerPragerModel {
public:
    struct Parameters {
        double friction_angle;      // φ [degrees or radians]
        double cohesion;            // c [Pa]
        double dilation_angle;      // ψ [degrees or radians] - for non-associative flow
        double tensile_strength;    // σ_t [Pa] - cutoff for tension
        
        // Hardening parameters
        double hardening_modulus;   // H [Pa]
        double saturation_yield;    // Maximum yield strength [Pa]
        
        // Softening parameters
        double softening_modulus;   // H_soft [Pa] (negative)
        double residual_cohesion;   // c_res [Pa] - after full softening
        double residual_friction;   // φ_res - residual friction angle
        
        // Viscoplastic regularization (optional)
        double viscosity;           // η [Pa·s] - for rate dependence
        
        Parameters() :
            friction_angle(30.0 * M_PI / 180.0),  // 30 degrees
            cohesion(1e6),                         // 1 MPa
            dilation_angle(10.0 * M_PI / 180.0),   // 10 degrees
            tensile_strength(0.1e6),               // 0.1 MPa
            hardening_modulus(0.0),
            saturation_yield(100e6),
            softening_modulus(0.0),
            residual_cohesion(0.1e6),
            residual_friction(20.0 * M_PI / 180.0),
            viscosity(0.0) {}
    };
    
    DruckerPragerModel();
    
    void setParameters(const Parameters& params);
    const Parameters& getParameters() const { return params; }
    
    // Compute yield function value
    double yieldFunction(const std::array<double, 6>& stress,
                        double hardening_variable) const;
    
    // Compute plastic flow direction (gradient of plastic potential)
    void flowDirection(const std::array<double, 6>& stress,
                      std::array<double, 6>& flow_dir) const;
    
    // Compute hardening/softening
    double hardeningFunction(double accumulated_plastic_strain) const;
    
    // Return mapping algorithm (stress integration)
    void returnMapping(std::array<double, 6>& stress,
                      std::array<double, 6>& plastic_strain,
                      double& accumulated_plastic_strain,
                      const std::array<double, 6>& strain_increment,
                      const std::array<std::array<double, 6>, 6>& elastic_modulus,
                      bool& plastic_loading);
    
private:
    Parameters params;
    
    // Internal parameters (computed from friction and cohesion)
    double alpha;   // Pressure sensitivity
    double kappa;   // Yield strength
    
    void updateInternalParameters();
};

/**
 * @brief von Mises (J2) plasticity
 * 
 * Pressure-independent plasticity. Simpler than Drucker-Prager but
 * less realistic for rocks. Good for metals and as a simplified model.
 */
class VonMisesModel {
public:
    struct Parameters {
        double yield_stress;        // σ_y [Pa]
        double hardening_modulus;   // H [Pa]
        double saturation_stress;   // σ_sat [Pa]
        
        Parameters() :
            yield_stress(10e6),
            hardening_modulus(1e9),
            saturation_stress(100e6) {}
    };
    
    VonMisesModel();
    
    void setParameters(const Parameters& params);
    
    double yieldFunction(const std::array<double, 6>& stress,
                        double hardening_variable) const;
    
    void returnMapping(std::array<double, 6>& stress,
                      std::array<double, 6>& plastic_strain,
                      double& accumulated_plastic_strain,
                      const std::array<double, 6>& strain_increment,
                      const std::array<std::array<double, 6>, 6>& elastic_modulus,
                      bool& plastic_loading);
    
private:
    Parameters params;
};

/**
 * @brief Mohr-Coulomb yield criterion
 * 
 * More complex than Drucker-Prager with corners in principal stress space.
 * Requires special treatment at corners.
 */
class MohrCoulombModel {
public:
    struct Parameters {
        double friction_angle;      // φ [radians]
        double cohesion;            // c [Pa]
        double dilation_angle;      // ψ [radians]
        double tension_cutoff;      // σ_t [Pa]
        
        Parameters() :
            friction_angle(30.0 * M_PI / 180.0),
            cohesion(1e6),
            dilation_angle(10.0 * M_PI / 180.0),
            tension_cutoff(0.1e6) {}
    };
    
    MohrCoulombModel();
    
    void setParameters(const Parameters& params);
    
    double yieldFunction(const std::array<double, 6>& stress,
                        double hardening_variable) const;
    
    void returnMapping(std::array<double, 6>& stress,
                      std::array<double, 6>& plastic_strain,
                      double& accumulated_plastic_strain,
                      const std::array<double, 6>& strain_increment,
                      const std::array<std::array<double, 6>, 6>& elastic_modulus,
                      bool& plastic_loading);
    
private:
    Parameters params;
};

/**
 * @brief Cap model (for porous rocks under compression)
 * 
 * Combines shear failure (Drucker-Prager) with a compression cap
 * to limit volumetric compaction. Important for:
 * - Sandstone compaction
 * - Pore collapse
 * - Reservoir geomechanics
 */
class CapModel {
public:
    struct Parameters {
        // Shear yield surface (Drucker-Prager)
        double friction_angle;
        double cohesion;
        
        // Cap surface (elliptical)
        double cap_eccentricity;    // R - shape of ellipse
        double initial_cap_position; // X0 - initial position
        double cap_hardening_modulus; // D - hardening rate
        
        // Transition region
        double transition_alpha;     // α - smooth transition parameter
        
        Parameters() :
            friction_angle(30.0 * M_PI / 180.0),
            cohesion(1e6),
            cap_eccentricity(0.3),
            initial_cap_position(-50e6),
            cap_hardening_modulus(1e10),
            transition_alpha(0.05) {}
    };
    
    CapModel();
    
    void setParameters(const Parameters& params);
    
    // Yield functions (two surfaces)
    double shearYield(const std::array<double, 6>& stress, double kappa) const;
    double capYield(const std::array<double, 6>& stress, double X) const;
    
    void returnMapping(std::array<double, 6>& stress,
                      std::array<double, 6>& plastic_strain,
                      double& accumulated_plastic_strain,
                      double& cap_position,
                      const std::array<double, 6>& strain_increment,
                      const std::array<std::array<double, 6>, 6>& elastic_modulus,
                      bool& plastic_loading);
    
private:
    Parameters params;
};

/**
 * @brief Plasticity state variables
 */
struct PlasticityState {
    std::array<double, 6> plastic_strain;        // ε^p
    double accumulated_plastic_strain;            // ε_acc^p = ∫||dε^p||dt
    double hardening_variable;                    // κ (internal state)
    bool is_plastic;                              // Currently yielding?
    double yield_function_value;                  // F(σ, κ)
    
    // For cap model
    double cap_position;                          // X
    
    PlasticityState() :
        plastic_strain({0, 0, 0, 0, 0, 0}),
        accumulated_plastic_strain(0.0),
        hardening_variable(0.0),
        is_plastic(false),
        yield_function_value(0.0),
        cap_position(0.0) {}
};

/**
 * @brief Unified plasticity interface
 * 
 * Provides a unified interface to different plasticity models
 * for integration into the main simulation code.
 */
class PlasticityModel {
public:
    PlasticityModel(YieldCriterion criterion = YieldCriterion::DRUCKER_PRAGER);
    
    // Set which model to use
    void setYieldCriterion(YieldCriterion criterion);
    void setFlowRule(FlowRule rule);
    void setHardeningLaw(HardeningLaw law);
    
    // Configure model
    void configure(const std::map<std::string, std::string>& config);
    void setDruckerPragerParameters(const DruckerPragerModel::Parameters& params);
    void setVonMisesParameters(const VonMisesModel::Parameters& params);
    void setMohrCoulombParameters(const MohrCoulombModel::Parameters& params);
    void setCapModelParameters(const CapModel::Parameters& params);
    
    // Stress integration (main function)
    void integrateStress(std::array<double, 6>& stress,
                        const std::array<double, 6>& strain_increment,
                        const std::array<std::array<double, 6>, 6>& elastic_modulus,
                        PlasticityState& state);
    
    // Check if yielding
    bool isYielding(const std::array<double, 6>& stress,
                   const PlasticityState& state) const;
    
    // Get yield function value
    double getYieldFunctionValue(const std::array<double, 6>& stress,
                                 const PlasticityState& state) const;
    
    // Get consistent tangent modulus (for Newton iterations)
    void getConsistentTangent(const std::array<double, 6>& stress,
                             const PlasticityState& state,
                             const std::array<std::array<double, 6>, 6>& elastic_modulus,
                             std::array<std::array<double, 6>, 6>& tangent_modulus);
    
    // Enable/disable plasticity
    void enable(bool enable) { enabled = enable; }
    bool isEnabled() const { return enabled; }
    
private:
    bool enabled;
    YieldCriterion criterion;
    FlowRule flow_rule;
    HardeningLaw hardening_law;
    
    // Specific model implementations
    std::unique_ptr<DruckerPragerModel> drucker_prager;
    std::unique_ptr<VonMisesModel> von_mises;
    std::unique_ptr<MohrCoulombModel> mohr_coulomb;
    std::unique_ptr<CapModel> cap_model;
};

/**
 * @brief Damage mechanics (continuum damage model)
 * 
 * Damage reduces elastic stiffness without explicit plasticity.
 * Useful for modeling:
 * - Microcrack accumulation
 * - Brittle failure
 * - Progressive damage
 */
class DamageModel {
public:
    struct Parameters {
        double damage_threshold;        // ε_D0 - strain threshold for damage
        double damage_exponent;         // α - damage evolution rate
        double maximum_damage;          // D_max - saturation damage
        double healing_rate;            // β - damage healing (if any)
        
        Parameters() :
            damage_threshold(1e-4),
            damage_exponent(2.0),
            maximum_damage(0.99),
            healing_rate(0.0) {}
    };
    
    DamageModel();
    
    void setParameters(const Parameters& params);
    
    // Update damage variable
    void updateDamage(double& damage,
                     const std::array<double, 6>& stress,
                     const std::array<double, 6>& strain,
                     double dt);
    
    // Reduce elastic modulus by damage
    void applyDamage(std::array<std::array<double, 6>, 6>& elastic_modulus,
                    double damage) const;
    
private:
    Parameters params;
};

/**
 * @brief Configuration helper for plasticity
 */
struct PlasticityConfig {
    bool enabled = false;
    std::string yield_criterion = "drucker_prager";  // drucker_prager, von_mises, mohr_coulomb, cap
    std::string flow_rule = "non_associative";        // associative, non_associative
    std::string hardening_law = "none";               // none, linear, exponential, softening
    
    // Drucker-Prager parameters
    double friction_angle = 30.0;       // degrees
    double cohesion = 1.0e6;            // Pa
    double dilation_angle = 10.0;       // degrees
    double tensile_strength = 0.1e6;    // Pa
    
    // Hardening/softening
    double hardening_modulus = 0.0;
    double softening_modulus = 0.0;
    double residual_cohesion = 0.1e6;
    
    // Parse from config file
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Create plasticity model
    PlasticityModel createModel() const;
};

} // namespace FSRM

#endif // PLASTICITY_MODEL_HPP
