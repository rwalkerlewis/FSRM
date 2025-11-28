#ifndef VISCOELASTIC_ATTENUATION_HPP
#define VISCOELASTIC_ATTENUATION_HPP

#include "MaterialModel.hpp"
#include <vector>
#include <array>
#include <map>
#include <string>
#include <complex>

namespace FSRM {

/**
 * @brief Multi-mechanism viscoelastic attenuation
 * 
 * Implements frequency-dependent attenuation using multiple Maxwell bodies
 * (Standard Linear Solid or Generalized Maxwell Body).
 * 
 * This matches SeisSol's attenuation implementation and provides:
 * - Frequency-independent quality factor Q over bandwidth
 * - Proper dispersion relations
 * - Causality-preserving attenuation
 * 
 * Theory:
 * Each Maxwell body has a relaxation frequency f_i and contributes to
 * the complex modulus:
 * 
 * M(ω) = M_∞ [1 + Σ_i (Y_i * (iωτ_i)) / (1 + iωτ_i)]
 * 
 * where τ_i = 1/(2πf_i) is the relaxation time.
 * 
 * References:
 * - Moczo et al. (2007) - The Finite-Difference Modelling of Earthquake Motions
 * - Kaser & Dumbser (2006) - ADER-DG method for viscoelasticity
 * - SeisSol implementation
 */
class ViscoelasticAttenuation {
public:
    /**
     * @brief Maxwell mechanism (single relaxation process)
     */
    struct MaxwellMechanism {
        double relaxation_frequency;    // f_i [Hz]
        double relaxation_time;         // τ_i = 1/(2πf_i) [s]
        double anelastic_coefficient_p; // Y_i^p for P-waves
        double anelastic_coefficient_s; // Y_i^s for S-waves
        
        MaxwellMechanism(double freq = 1.0, double Yp = 0.0, double Ys = 0.0) :
            relaxation_frequency(freq),
            relaxation_time(1.0 / (2.0 * M_PI * freq)),
            anelastic_coefficient_p(Yp),
            anelastic_coefficient_s(Ys) {}
    };
    
    /**
     * @brief Anelastic state variables (memory variables)
     * 
     * Each Maxwell mechanism requires additional state variables to
     * track the history-dependent stress relaxation.
     */
    struct AnelasticState {
        std::vector<std::array<double, 6>> R;  // Memory variables for each mechanism
        
        AnelasticState(int num_mechanisms = 3) {
            R.resize(num_mechanisms, {0, 0, 0, 0, 0, 0});
        }
        
        void reset() {
            for (auto& r : R) {
                r = {0, 0, 0, 0, 0, 0};
            }
        }
    };
    
    /**
     * @brief Configuration parameters
     */
    struct Parameters {
        int num_mechanisms;             // Number of Maxwell bodies (typically 3)
        double freq_central;            // Central frequency [Hz]
        double freq_ratio;              // Ratio f_max/f_min
        double Qp_desired;              // Desired P-wave quality factor
        double Qs_desired;              // Desired S-wave quality factor
        
        // Optional: directly specify mechanisms
        std::vector<MaxwellMechanism> mechanisms;
        
        Parameters() :
            num_mechanisms(3),
            freq_central(0.5),
            freq_ratio(100.0),
            Qp_desired(50.0),
            Qs_desired(25.0) {}
    };
    
    ViscoelasticAttenuation();
    
    // Configuration
    void setParameters(const Parameters& params);
    const Parameters& getParameters() const { return params; }
    
    // Setup mechanisms from Q values
    void setupMechanismsFromQ(int num_mechanisms, double freq_central,
                              double freq_ratio, double Qp, double Qs);
    
    // Add mechanism manually
    void addMechanism(const MaxwellMechanism& mech);
    void clearMechanisms();
    
    // Get mechanisms
    const std::vector<MaxwellMechanism>& getMechanisms() const { return mechanisms; }
    int getNumMechanisms() const { return mechanisms.size(); }
    
    // Compute anelastic coefficients from desired Q
    void computeAnelasticCoefficients();
    
    // Update anelastic state variables (time integration)
    void updateAnelasticState(AnelasticState& state,
                             const std::array<double, 6>& stress,
                             double dt);
    
    // Compute stress contribution from anelastic mechanisms
    void computeAnelasticStress(const AnelasticState& state,
                               std::array<double, 6>& stress_correction);
    
    // Get effective (relaxed) moduli
    double getRelaxedModulus(double M_unrelaxed, bool is_shear) const;
    double getRelaxedVelocity(double v_unrelaxed, bool is_shear) const;
    
    // Frequency-dependent properties
    std::complex<double> getComplexModulus(double frequency, double M_unrelaxed,
                                          bool is_shear) const;
    double getQuality(double frequency, bool is_shear) const;
    double getVelocity(double frequency, double v_unrelaxed, bool is_shear) const;
    
    // Dispersion curve (for plotting/analysis)
    struct DispersionPoint {
        double frequency;
        double velocity_p;
        double velocity_s;
        double Qp;
        double Qs;
    };
    std::vector<DispersionPoint> computeDispersionCurve(
        double vp_ref, double vs_ref,
        double f_min, double f_max, int num_points) const;
    
    // Verify Q approximation quality
    double getMaxQError(double f_min, double f_max, bool is_shear) const;
    
    // Enable/disable attenuation
    void enable(bool enable) { enabled = enable; }
    bool isEnabled() const { return enabled; }
    
private:
    bool enabled;
    Parameters params;
    std::vector<MaxwellMechanism> mechanisms;
    
    // Precomputed factors for efficiency
    std::vector<double> exp_factors;  // exp(-dt/τ_i) precomputed
    
    void precomputeFactors(double dt);
    
    // Optimization: compute all anelastic coefficients to match target Q
    void solveForAnelasticCoefficients(const std::vector<double>& frequencies,
                                       double Q_target,
                                       std::vector<double>& Y_coefficients);
};

/**
 * @brief Attenuation models for different rock types
 * 
 * Provides pre-configured attenuation parameters for common materials.
 */
class AttenuationDatabase {
public:
    enum class RockType {
        GRANITE,
        BASALT,
        SANDSTONE,
        LIMESTONE,
        SHALE,
        SEDIMENT_SOFT,
        SEDIMENT_HARD,
        FAULT_ZONE,
        DAMAGE_ZONE,
        CUSTOM
    };
    
    static ViscoelasticAttenuation::Parameters getParameters(RockType rock);
    
    // Empirical Q relations
    static double estimateQs(double vs_km_per_s);  // Q_s ≈ 40-50 * V_s
    static double estimateQp(double Qs);           // Q_p ≈ 2 * Q_s
    
    // Depth dependence
    static double depthCorrection(double Q_surface, double depth_km);
    
    // Damage correction
    static double damageCorrection(double Q_intact, double damage);
};

/**
 * @brief Spatial variation of attenuation properties
 * 
 * Allows Q to vary spatially based on:
 * - Material type
 * - Depth
 * - Velocity
 * - Damage/fractures
 */
class SpatialAttenuation {
public:
    SpatialAttenuation();
    
    // Set base attenuation model
    void setBaseModel(const ViscoelasticAttenuation::Parameters& params);
    
    // Enable spatial variation
    void enableVelocityScaling(bool enable, double scaling_factor = 50.0);
    void enableDepthScaling(bool enable, double gradient = 0.1);
    void enableDamageScaling(bool enable);
    
    // Get parameters at specific location
    ViscoelasticAttenuation::Parameters getParametersAt(
        double x, double y, double z,
        double vp, double vs,
        double damage = 0.0) const;
    
private:
    ViscoelasticAttenuation::Parameters base_params;
    bool velocity_scaling;
    double velocity_scaling_factor;
    bool depth_scaling;
    double depth_gradient;
    bool damage_scaling;
};

/**
 * @brief Time integration schemes for anelastic state
 */
class AnelasticIntegrator {
public:
    enum class Scheme {
        FORWARD_EULER,      // First-order
        RUNGE_KUTTA_2,      // Second-order
        RUNGE_KUTTA_4,      // Fourth-order
        ADER                // ADER (matches spatial order)
    };
    
    AnelasticIntegrator(Scheme scheme = Scheme::RUNGE_KUTTA_2);
    
    void setScheme(Scheme scheme) { this->scheme = scheme; }
    
    // Integrate anelastic state forward in time
    void integrate(ViscoelasticAttenuation& attenuation,
                  ViscoelasticAttenuation::AnelasticState& state,
                  const std::array<double, 6>& stress,
                  double dt);
    
private:
    Scheme scheme;
    
    void forwardEuler(ViscoelasticAttenuation& attenuation,
                     ViscoelasticAttenuation::AnelasticState& state,
                     const std::array<double, 6>& stress,
                     double dt);
    
    void rungeKutta2(ViscoelasticAttenuation& attenuation,
                    ViscoelasticAttenuation::AnelasticState& state,
                    const std::array<double, 6>& stress,
                    double dt);
    
    void rungeKutta4(ViscoelasticAttenuation& attenuation,
                    ViscoelasticAttenuation::AnelasticState& state,
                    const std::array<double, 6>& stress,
                    double dt);
};

/**
 * @brief Configuration helper for viscoelastic attenuation
 */
struct ViscoelasticConfig {
    bool enabled = false;
    
    // Basic parameters
    int num_mechanisms = 3;
    double freq_central = 0.5;      // Hz
    double freq_ratio = 100.0;      // f_max/f_min
    
    // Q values
    double Qp = 50.0;
    double Qs = 25.0;
    
    // Spatial variation
    bool velocity_dependent = false;
    double velocity_scaling = 50.0;  // Q_s ≈ 50 * V_s (km/s)
    
    bool depth_dependent = false;
    double depth_gradient = 0.1;     // dQ/dz per km
    
    // Rock type (for database)
    std::string rock_type = "granite";
    
    // Parse from config file
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Create attenuation model
    ViscoelasticAttenuation createModel() const;
};

/**
 * @brief Utility functions for attenuation analysis
 */
namespace AttenuationUtils {
    
    // Verify causality (Kramers-Kronig relations)
    bool checkCausality(const ViscoelasticAttenuation& model,
                       double f_min, double f_max, double tolerance = 1e-3);
    
    // Compute phase velocity from group velocity and Q
    double phaseVelocity(double group_velocity, double Q, double frequency);
    
    // Compute attenuation coefficient α (Nepers/m)
    double attenuationCoefficient(double frequency, double velocity, double Q);
    
    // Amplitude reduction over distance
    double amplitudeReduction(double distance, double frequency,
                             double velocity, double Q);
    
    // Quality factor from decay measurement
    double computeQ(double amplitude_ratio, double distance,
                   double wavelength);
    
    // Spectral ratio method for Q estimation
    std::pair<double, double> estimateQFromSpectralRatio(
        const std::vector<double>& spectrum1,
        const std::vector<double>& spectrum2,
        const std::vector<double>& frequencies,
        double distance_diff);
    
    // Plot dispersion curves (output to file)
    void plotDispersionCurves(const ViscoelasticAttenuation& model,
                             double vp_ref, double vs_ref,
                             const std::string& filename);
}

/**
 * @brief Combine viscoelasticity with other physics
 */
class ViscoelasticCoupling {
public:
    ViscoelasticCoupling();
    
    // Couple with poroelasticity
    void setCoupledPoroelastic(bool coupled) { poroelastic_coupling = coupled; }
    
    // Couple with plasticity (Q reduction in plastic zone)
    void setCoupledPlasticity(bool coupled) { plasticity_coupling = coupled; }
    void setPlasticQReduction(double factor) { plastic_q_reduction = factor; }
    
    // Couple with damage
    void setCoupledDamage(bool coupled) { damage_coupling = coupled; }
    void setDamageQRelation(double exponent) { damage_q_exponent = exponent; }
    
    // Modify Q based on coupled physics
    double getModifiedQ(double Q_base, double porosity, double damage,
                       bool is_plastic) const;
    
private:
    bool poroelastic_coupling;
    bool plasticity_coupling;
    double plastic_q_reduction;
    bool damage_coupling;
    double damage_q_exponent;
};

} // namespace FSRM

#endif // VISCOELASTIC_ATTENUATION_HPP
