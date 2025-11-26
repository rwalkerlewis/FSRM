#ifndef MATERIAL_MODEL_HPP
#define MATERIAL_MODEL_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <array>

namespace FSRM {

/**
 * @brief Types of constitutive models for rock/formation behavior
 */
enum class ConstitutiveModel {
    LINEAR_ELASTIC,         // Hookean elasticity
    VISCOELASTIC_MAXWELL,   // Maxwell viscoelastic
    VISCOELASTIC_KELVIN,    // Kelvin-Voigt viscoelastic
    VISCOELASTIC_SLS,       // Standard Linear Solid
    POROELASTIC,            // Biot poroelasticity
    THERMOELASTIC,          // Temperature-dependent elasticity
    ELASTOPLASTIC_MC,       // Mohr-Coulomb plasticity
    ELASTOPLASTIC_DP,       // Drucker-Prager plasticity
    ELASTOPLASTIC_CAP,      // Cap model with pore collapse
    DAMAGE,                 // Continuum damage mechanics
    ANISOTROPIC_VTI,        // Vertically transverse isotropic
    ANISOTROPIC_HTI         // Horizontally transverse isotropic
};

/**
 * @brief Types of failure/yield criteria
 */
enum class FailureCriterion {
    VON_MISES,              // J2 plasticity
    MOHR_COULOMB,           // Friction-based
    DRUCKER_PRAGER,         // Smoothed Mohr-Coulomb
    HOEK_BROWN,             // Rock masses
    MODIFIED_CAM_CLAY,      // Soils
    TENSILE_CUTOFF          // Maximum tensile stress
};

/**
 * @brief Stress tensor in Voigt notation
 */
struct StressTensor {
    double xx, yy, zz;      // Normal stresses
    double xy, xz, yz;      // Shear stresses
    
    double mean() const { return (xx + yy + zz) / 3.0; }
    double vonMises() const;
    double maxPrincipal() const;
    double minPrincipal() const;
    std::array<double, 3> principalStresses() const;
    double maxShear() const;
};

/**
 * @brief Strain tensor in Voigt notation
 */
struct StrainTensor {
    double xx, yy, zz;      // Normal strains
    double xy, xz, yz;      // Engineering shear strains (γ = 2ε)
    
    double volumetric() const { return xx + yy + zz; }
    double deviatoric() const;
};

/**
 * @brief Base class for material models
 * 
 * All properties and behaviors configurable from text files.
 */
class MaterialModelBase {
public:
    MaterialModelBase(ConstitutiveModel type) : model_type(type) {}
    virtual ~MaterialModelBase() = default;
    
    // Compute stress from strain
    virtual StressTensor computeStress(const StrainTensor& strain,
                                       double P_pore = 0.0,
                                       double T = 293.15) const = 0;
    
    // Compute elastic modulus tensor (6x6 in Voigt notation)
    virtual std::array<double, 36> getStiffnessMatrix(double P = 0.0, 
                                                       double T = 293.15) const = 0;
    
    // Configuration from key-value pairs
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    ConstitutiveModel getType() const { return model_type; }
    
    // Common property accessors
    virtual double getDensity() const = 0;
    virtual double getPorosity() const = 0;
    virtual double getYoungsModulus() const = 0;
    virtual double getPoissonsRatio() const = 0;
    
protected:
    ConstitutiveModel model_type;
};

/**
 * @brief Linear elastic material model
 * 
 * Standard Hookean elasticity with optional pressure/temperature dependence.
 */
class LinearElasticMaterial : public MaterialModelBase {
public:
    LinearElasticMaterial();
    
    StressTensor computeStress(const StrainTensor& strain,
                               double P_pore = 0.0,
                               double T = 293.15) const override;
    
    std::array<double, 36> getStiffnessMatrix(double P = 0.0,
                                               double T = 293.15) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Property accessors
    double getDensity() const override { return density; }
    double getPorosity() const override { return porosity; }
    double getYoungsModulus() const override { return youngs_modulus; }
    double getPoissonsRatio() const override { return poisson_ratio; }
    
    // Setters
    void setElasticModuli(double E, double nu);
    void setLameParameters(double lambda, double mu);
    void setDensity(double rho);
    void setPorosity(double phi);
    
    // Derived properties
    double getBulkModulus() const;
    double getShearModulus() const;
    double getLambda() const;
    double getPWaveVelocity() const;
    double getSWaveVelocity() const;
    
private:
    double youngs_modulus;      // E (Pa)
    double poisson_ratio;       // ν
    double density;             // ρ (kg/m³)
    double porosity;            // φ
    
    // Pressure/temperature dependence coefficients
    double E_P_coeff;           // dE/dP
    double E_T_coeff;           // dE/dT
    double nu_P_coeff;          // dν/dP
    double nu_T_coeff;          // dν/dT
    double P_reference;
    double T_reference;
};

/**
 * @brief Viscoelastic material model (Maxwell, Kelvin-Voigt, SLS)
 * 
 * Time-dependent response with configurable relaxation behavior.
 */
class ViscoelasticMaterial : public MaterialModelBase {
public:
    enum class ViscoType { MAXWELL, KELVIN_VOIGT, STANDARD_LINEAR_SOLID };
    
    ViscoelasticMaterial(ViscoType type = ViscoType::MAXWELL);
    
    StressTensor computeStress(const StrainTensor& strain,
                               double P_pore = 0.0,
                               double T = 293.15) const override;
    
    // Time-dependent stress update
    StressTensor computeStressRate(const StrainTensor& strain,
                                   const StrainTensor& strain_rate,
                                   const StressTensor& current_stress,
                                   double dt) const;
    
    std::array<double, 36> getStiffnessMatrix(double P = 0.0,
                                               double T = 293.15) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    double getDensity() const override { return density; }
    double getPorosity() const override { return porosity; }
    double getYoungsModulus() const override { return E_instantaneous; }
    double getPoissonsRatio() const override { return poisson_ratio; }
    
    // Setters
    void setInstantaneousModuli(double E, double nu);
    void setRelaxationTime(double tau);
    void setViscosity(double eta);
    void setLongTermModulus(double E_inf);  // For SLS
    
    // Relaxation function G(t)
    double getRelaxationModulus(double t) const;
    double getCreepCompliance(double t) const;
    
private:
    ViscoType visco_type;
    double E_instantaneous;     // Instantaneous Young's modulus
    double E_long_term;         // Long-term modulus (SLS only)
    double poisson_ratio;
    double relaxation_time;     // τ (s)
    double viscosity;           // η (Pa·s)
    double density;
    double porosity;
    
    // Multiple relaxation times (Prony series)
    std::vector<double> prony_g;    // Shear modulus ratios
    std::vector<double> prony_tau;  // Relaxation times
};

/**
 * @brief Poroelastic material model (Biot theory)
 * 
 * Coupled fluid-solid mechanics with configurable coupling parameters.
 */
class PoroelasticMaterial : public MaterialModelBase {
public:
    PoroelasticMaterial();
    
    StressTensor computeStress(const StrainTensor& strain,
                               double P_pore = 0.0,
                               double T = 293.15) const override;
    
    // Total stress = effective stress - α*P*I
    StressTensor computeEffectiveStress(const StrainTensor& strain) const;
    StressTensor computeTotalStress(const StressTensor& effective_stress,
                                    double P_pore) const;
    
    std::array<double, 36> getStiffnessMatrix(double P = 0.0,
                                               double T = 293.15) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    double getDensity() const override { return density_solid; }
    double getPorosity() const override { return porosity; }
    double getYoungsModulus() const override { return youngs_modulus; }
    double getPoissonsRatio() const override { return poisson_ratio; }
    
    // Poroelastic parameters
    double getBiotCoefficient() const { return biot_coefficient; }
    double getBiotModulus() const { return biot_modulus; }
    double getSkemptonCoefficient() const;
    double getUndrainedBulkModulus() const;
    double getStorageCoefficient() const;
    
    // Wave velocities
    double getFastPWaveVelocity() const;
    double getSlowPWaveVelocity() const;  // Biot slow wave
    double getSWaveVelocity() const;
    
    // Setters
    void setDrainedModuli(double E, double nu);
    void setBiotParameters(double alpha, double M);
    void setFluidProperties(double Kf, double rho_f, double mu_f);
    void setPermeability(double k);
    void setPorosity(double phi);
    
private:
    // Drained elastic properties
    double youngs_modulus;
    double poisson_ratio;
    double bulk_modulus_drained;
    double shear_modulus;
    
    // Poroelastic coupling
    double biot_coefficient;    // α (0 to 1)
    double biot_modulus;        // M (Pa)
    
    // Solid and fluid properties
    double density_solid;
    double density_fluid;
    double bulk_modulus_solid;  // K_s (grain bulk modulus)
    double bulk_modulus_fluid;  // K_f
    double fluid_viscosity;
    double porosity;
    double permeability;        // k (m²)
};

/**
 * @brief Elastoplastic material with Mohr-Coulomb or Drucker-Prager yield
 */
class ElastoplasticMaterial : public MaterialModelBase {
public:
    ElastoplasticMaterial(FailureCriterion criterion = FailureCriterion::MOHR_COULOMB);
    
    StressTensor computeStress(const StrainTensor& strain,
                               double P_pore = 0.0,
                               double T = 293.15) const override;
    
    // Plastic stress update (returns corrector)
    StressTensor computePlasticCorrection(const StressTensor& trial_stress,
                                          double& plastic_strain_eq) const;
    
    std::array<double, 36> getStiffnessMatrix(double P = 0.0,
                                               double T = 293.15) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    double getDensity() const override { return density; }
    double getPorosity() const override { return porosity; }
    double getYoungsModulus() const override { return youngs_modulus; }
    double getPoissonsRatio() const override { return poisson_ratio; }
    
    // Yield surface evaluation (negative = elastic, positive = yielded)
    double evaluateYieldFunction(const StressTensor& stress) const;
    bool isYielding(const StressTensor& stress) const;
    
    // Failure properties
    void setMohrCoulombParameters(double cohesion, double friction_angle,
                                  double dilation_angle = 0.0);
    void setDruckerPragerParameters(double d, double k);  // Matched to MC
    void setTensileStrength(double sigma_t);
    void setHardeningParameters(double H, double sigma_y);
    
    double getCohesion() const { return cohesion; }
    double getFrictionAngle() const { return friction_angle; }
    
private:
    FailureCriterion criterion;
    
    // Elastic properties
    double youngs_modulus;
    double poisson_ratio;
    double density;
    double porosity;
    
    // Mohr-Coulomb parameters
    double cohesion;            // c (Pa)
    double friction_angle;      // φ (radians)
    double dilation_angle;      // ψ (radians)
    double tensile_strength;    // σ_t (Pa)
    
    // Drucker-Prager parameters (converted from MC)
    double dp_d;
    double dp_k;
    
    // Hardening
    double hardening_modulus;   // H
    double initial_yield_stress;
    
    // Plastic state (mutable for const methods)
    mutable double accumulated_plastic_strain;
};

/**
 * @brief Anisotropic (VTI or HTI) elastic material
 */
class AnisotropicMaterial : public MaterialModelBase {
public:
    AnisotropicMaterial(bool vertical_ti = true);
    
    StressTensor computeStress(const StrainTensor& strain,
                               double P_pore = 0.0,
                               double T = 293.15) const override;
    
    std::array<double, 36> getStiffnessMatrix(double P = 0.0,
                                               double T = 293.15) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    double getDensity() const override { return density; }
    double getPorosity() const override { return porosity; }
    double getYoungsModulus() const override { return E_vertical; }
    double getPoissonsRatio() const override { return nu_vh; }
    
    // Set Thomsen parameters (for seismic anisotropy)
    void setThomsenParameters(double epsilon, double delta, double gamma);
    
    // Set engineering constants
    void setVerticalModuli(double Ev, double Gv, double nuv);
    void setHorizontalModuli(double Eh, double Gh, double nuh);
    
    // Wave velocities
    double getVP0() const;      // Vertical P-wave
    double getVS0() const;      // Vertical S-wave  
    double getVPh() const;      // Horizontal P-wave
    double getVSh() const;      // Horizontal S-wave
    
private:
    bool is_VTI;  // true = VTI, false = HTI
    
    // Stiffness tensor components (VTI notation)
    double C11, C13, C33, C44, C66;
    
    // Engineering constants (alternative parameterization)
    double E_vertical;
    double E_horizontal;
    double G_vertical;
    double nu_vh;  // Poisson's ratio for vertical loading
    double nu_hh;  // Poisson's ratio in horizontal plane
    
    // Thomsen anisotropy parameters
    double epsilon;  // P-wave anisotropy
    double delta;    // Near-vertical P-wave anisotropy
    double gamma;    // S-wave anisotropy
    
    double density;
    double porosity;
};

/**
 * @brief Thermal expansion properties
 */
struct ThermalProperties {
    double thermal_conductivity;    // k (W/m·K)
    double specific_heat;           // c_p (J/kg·K)
    double thermal_expansion;       // α_T (1/K)
    double reference_temperature;   // T_0 (K)
    
    // Temperature-dependent coefficients
    double k_T_coeff;   // dk/dT
    double cp_T_coeff;  // dcp/dT
    
    // Default values
    ThermalProperties() :
        thermal_conductivity(2.5),
        specific_heat(900.0),
        thermal_expansion(1e-5),
        reference_temperature(293.15),
        k_T_coeff(0.0),
        cp_T_coeff(0.0) {}
    
    void configure(const std::map<std::string, std::string>& config);
    
    double getConductivity(double T) const;
    double getSpecificHeat(double T) const;
    StressTensor getThermalStress(double E, double nu, double dT) const;
};

/**
 * @brief Permeability model with stress/strain dependence
 */
class PermeabilityModel {
public:
    enum class PermeabilityType {
        CONSTANT,           // No stress dependence
        KOZENY_CARMAN,      // Porosity-based
        EXPONENTIAL,        // Exponential stress dependence
        POWER_LAW,          // Power-law stress dependence
        CUBIC_LAW,          // For fractures (k ∝ w³)
        DYNAMIC             // Wave-induced changes
    };
    
    PermeabilityModel(PermeabilityType type = PermeabilityType::CONSTANT);
    
    void configure(const std::map<std::string, std::string>& config);
    
    // Get permeability (m²) given current conditions
    double getPermeability(double effective_stress, double porosity,
                          double aperture = 0.0) const;
    
    // For dynamic permeability changes
    double getDynamicPermeability(double k0, double volumetric_strain,
                                  double shear_strain, double dt) const;
    
    // Setters
    void setInitialPermeability(double k0);
    void setStressSensitivity(double coeff);
    void setPorosityExponent(double n);
    void setRecoveryTime(double tau);
    void setBounds(double k_min, double k_max);
    
    double getInitialPermeability() const { return k_initial; }
    
private:
    PermeabilityType type;
    
    double k_initial;           // k_0 (m²)
    double stress_coefficient;  // For exponential: k = k0 * exp(-c*σ')
    double porosity_exponent;   // For Kozeny-Carman: k ∝ φ^n / (1-φ)²
    double grain_diameter;      // For KC model
    
    // Dynamic permeability
    double strain_sensitivity;
    double recovery_time;
    
    // Bounds
    double k_min;
    double k_max;
};

/**
 * @brief Fracture mechanics properties for LEFM calculations
 */
struct FractureProperties {
    double toughness_KIc;       // Mode I fracture toughness (Pa·m^0.5)
    double toughness_KIIc;      // Mode II (optional)
    double toughness_KIIIc;     // Mode III (optional)
    double fracture_energy;     // G_c (J/m²) = K_Ic² / E'
    double tensile_strength;    // σ_t (Pa)
    double process_zone_length; // Cohesive zone length
    
    // Cohesive zone model parameters
    double cohesive_strength;   // Peak traction (Pa)
    double critical_opening;    // δ_c (m)
    
    FractureProperties() :
        toughness_KIc(1e6),
        toughness_KIIc(1.5e6),
        toughness_KIIIc(1e6),
        fracture_energy(100.0),
        tensile_strength(5e6),
        process_zone_length(0.01),
        cohesive_strength(5e6),
        critical_opening(2e-5) {}
    
    void configure(const std::map<std::string, std::string>& config);
    
    // Calculate critical stress intensity from energy
    double computeKIc(double E, double nu) const;
};

/**
 * @brief Complete rock/formation property set
 * 
 * Aggregates all material properties for a rock type.
 * Fully configurable from text files.
 */
class RockProperties {
public:
    RockProperties();
    
    void configure(const std::map<std::string, std::string>& config);
    
    // Sub-models (owned)
    std::unique_ptr<MaterialModelBase> mechanical;
    std::unique_ptr<PermeabilityModel> permeability;
    ThermalProperties thermal;
    FractureProperties fracture;
    
    // Common accessors
    double getDensity() const;
    double getPorosity() const;
    double getYoungsModulus() const;
    double getPoissonsRatio() const;
    double getPermeability() const;
    double getThermalConductivity() const;
    
    // Name/ID
    std::string name;
    int region_id;
    
private:
    void createMechanicalModel(const std::string& type,
                               const std::map<std::string, std::string>& config);
};

/**
 * @brief Factory function to create material models from configuration
 */
std::unique_ptr<MaterialModelBase> createMaterialModel(
    const std::string& type_str,
    const std::map<std::string, std::string>& config);

/**
 * @brief Parse constitutive model type from string
 */
ConstitutiveModel parseConstitutiveModel(const std::string& type_str);

} // namespace FSRM

#endif // MATERIAL_MODEL_HPP
