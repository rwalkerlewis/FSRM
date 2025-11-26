#ifndef FAULT_MODEL_HPP
#define FAULT_MODEL_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <array>

namespace FSRM {

/**
 * @brief Types of friction laws for fault mechanics
 */
enum class FrictionLaw {
    COULOMB,                // Static/dynamic Coulomb friction
    RATE_STATE_AGING,       // Rate-and-state with aging law
    RATE_STATE_SLIP,        // Rate-and-state with slip law
    SLIP_WEAKENING,         // Linear slip-weakening
    FLASH_HEATING,          // Flash heating at high slip rates
    THERMAL_PRESSURIZATION  // Pore pressure effects
};

/**
 * @brief Stress state on a fault plane
 */
struct FaultStressState {
    double sigma_n;         // Normal stress (compression positive, Pa)
    double sigma_n_eff;     // Effective normal stress (Pa)
    double tau;             // Shear stress (Pa)
    double tau_strength;    // Frictional strength (Pa)
    double CFF;             // Coulomb Failure Function (Pa)
    
    // Slip state
    double slip;            // Cumulative slip (m)
    double slip_rate;       // Current slip rate (m/s)
    double state_variable;  // Rate-state θ
    
    // Event properties
    double moment;          // Seismic moment (N·m)
    double stress_drop;     // Stress drop (Pa)
    bool is_slipping;       // Currently slipping?
    bool is_seismic;        // Seismic (vs aseismic) slip?
    
    // Helper functions
    double getMagnitude() const;  // Moment magnitude Mw
    bool isStable() const { return CFF < 0.0; }
    bool isCritical() const { return CFF >= 0.0 && !is_slipping; }
};

/**
 * @brief Seismic event record
 */
struct SeismicEvent {
    double time;            // Event time (s)
    double magnitude;       // Moment magnitude
    double moment;          // Seismic moment (N·m)
    double slip;            // Coseismic slip (m)
    double stress_drop;     // Stress drop (Pa)
    double duration;        // Rupture duration (s)
    double area;            // Rupture area (m²)
    
    // Location (center of rupture)
    double x, y, z;
    
    // Hypocenter (nucleation point)
    double hypo_x, hypo_y, hypo_z;
    
    std::string type;       // "earthquake", "slow_slip", "tremor"
};

/**
 * @brief Coulomb friction parameters
 */
struct CoulombParameters {
    double static_friction;     // μ_s
    double dynamic_friction;    // μ_d
    double cohesion;            // c (Pa)
    double slip_weakening_Dc;   // Critical slip distance (m)
    
    CoulombParameters() :
        static_friction(0.6),
        dynamic_friction(0.4),
        cohesion(1e6),
        slip_weakening_Dc(0.01) {}
    
    void configure(const std::map<std::string, std::string>& config);
};

/**
 * @brief Rate-and-state friction parameters
 */
struct RateStateParameters {
    double a;               // Direct effect parameter
    double b;               // Evolution effect parameter
    double Dc;              // State evolution distance (m)
    double V0;              // Reference velocity (m/s)
    double f0;              // Reference friction coefficient
    double theta0;          // Initial state variable (s)
    
    // Stability: a-b < 0 is velocity weakening (unstable)
    bool isVelocityWeakening() const { return a < b; }
    double getStabilityParameter() const { return a - b; }
    
    RateStateParameters() :
        a(0.010),
        b(0.015),
        Dc(1e-4),
        V0(1e-6),
        f0(0.6),
        theta0(1e6) {}
    
    void configure(const std::map<std::string, std::string>& config);
};

/**
 * @brief Base class for fault friction models
 */
class FrictionModelBase {
public:
    FrictionModelBase(FrictionLaw law) : friction_law(law) {}
    virtual ~FrictionModelBase() = default;
    
    // Compute friction coefficient
    virtual double getFriction(double slip_rate, double state_var,
                               double sigma_n_eff) const = 0;
    
    // Compute state variable evolution rate
    virtual double getStateEvolutionRate(double slip_rate, 
                                         double state_var) const = 0;
    
    // Configuration
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    FrictionLaw getLaw() const { return friction_law; }
    
protected:
    FrictionLaw friction_law;
};

/**
 * @brief Coulomb friction model with slip weakening
 */
class CoulombFriction : public FrictionModelBase {
public:
    CoulombFriction();
    
    double getFriction(double slip_rate, double state_var,
                       double sigma_n_eff) const override;
    
    double getStateEvolutionRate(double slip_rate, double state_var) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    CoulombParameters params;
    
private:
    double current_slip;  // For slip weakening
};

/**
 * @brief Rate-and-state friction model
 */
class RateStateFriction : public FrictionModelBase {
public:
    enum class EvolutionLaw { AGING, SLIP };
    
    RateStateFriction(EvolutionLaw law = EvolutionLaw::AGING);
    
    double getFriction(double slip_rate, double state_var,
                       double sigma_n_eff) const override;
    
    double getStateEvolutionRate(double slip_rate, double state_var) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    RateStateParameters params;
    
    // Steady-state friction at given slip rate
    double getSteadyStateFriction(double slip_rate) const;
    
    // Critical stiffness for stability (for spring-slider)
    double getCriticalStiffness(double sigma_n_eff) const;
    
private:
    EvolutionLaw evolution_law;
};

/**
 * @brief Fault geometry definition
 */
struct FaultGeometry {
    // Location (center point)
    double x, y, z;
    
    // Orientation
    double strike;          // Strike angle (radians from N)
    double dip;             // Dip angle (radians)
    double rake;            // Rake angle for slip direction (radians)
    
    // Dimensions
    double length;          // Along-strike (m)
    double width;           // Down-dip (m)
    
    // Derived
    double area() const { return length * width; }
    
    // Unit vectors
    std::array<double, 3> getStrikeVector() const;
    std::array<double, 3> getDipVector() const;
    std::array<double, 3> getNormalVector() const;
    
    // Transform stress tensor to fault-local coordinates
    void resolveStress(double sxx, double syy, double szz,
                      double sxy, double sxz, double syz,
                      double& sigma_n, double& tau) const;
    
    FaultGeometry() :
        x(0), y(0), z(0),
        strike(0), dip(M_PI/3),  // 60° dip default
        rake(0),
        length(1000), width(500) {}
    
    void configure(const std::map<std::string, std::string>& config);
};

/**
 * @brief Complete fault model with mechanics and event tracking
 * 
 * Note: This class provides advanced fault mechanics with rate-state friction.
 * For simple fault representation in fracture networks, see FractureModel.hpp::FaultModel
 */
class SeismicFaultModel {
public:
    SeismicFaultModel();
    
    // Configuration from key-value pairs
    void configure(const std::map<std::string, std::string>& config);
    
    // Set fault geometry
    void setGeometry(const FaultGeometry& geom);
    void setGeometry(double x, double y, double z,
                    double strike, double dip, double length, double width);
    
    // Set friction model
    void setFrictionModel(std::unique_ptr<FrictionModelBase> model);
    void setFrictionLaw(FrictionLaw law);
    
    // Compute stress state from regional stress field
    FaultStressState computeStressState(double sxx, double syy, double szz,
                                        double sxy, double sxz, double syz,
                                        double pore_pressure) const;
    
    // Time integration for slip evolution
    FaultStressState updateSlip(const FaultStressState& current_state,
                               double shear_modulus, double dt);
    
    // Check for event nucleation
    bool checkNucleation(const FaultStressState& state) const;
    
    // Compute event properties if slip occurs
    SeismicEvent computeEvent(const FaultStressState& initial_state,
                             const FaultStressState& final_state,
                             double shear_modulus,
                             double current_time) const;
    
    // Coulomb stress transfer from slip
    void computeCoulombStressChange(double slip, double receiver_x,
                                    double receiver_y, double receiver_z,
                                    double shear_modulus,
                                    double& dCFF) const;
    
    // Accessors
    const FaultGeometry& getGeometry() const { return geometry; }
    const std::vector<SeismicEvent>& getEventCatalog() const { return events; }
    
    // Event tracking
    void addEvent(const SeismicEvent& event);
    void clearCatalog();
    
    // Seismicity rate estimation (Dieterich 1994)
    double getSeismicityRate(double stressing_rate, double sigma_n_eff,
                            double temperature = 300.0) const;
    
    // Magnitude-frequency (Gutenberg-Richter)
    double getBValue() const { return b_value; }
    void setBValue(double b) { b_value = b; }
    
    std::string name;
    
private:
    FaultGeometry geometry;
    std::unique_ptr<FrictionModelBase> friction;
    std::vector<SeismicEvent> events;
    
    // Current state
    double current_slip;
    double current_state_var;
    double current_slip_rate;
    
    // Seismicity parameters
    double b_value;                 // G-R b-value (~1.0)
    double nucleation_size;         // Critical nucleation patch size
    double seismic_slip_rate;       // Threshold for seismic vs aseismic
    
    // Helper for Okada dislocation
    void okadaGreen(double x, double y, double z,
                   double strike, double dip, double slip,
                   double L, double W, double depth,
                   double& dux, double& duy, double& duz) const;
};

/**
 * @brief Network of faults with interaction
 */
class FaultNetwork {
public:
    FaultNetwork();
    
    // Add faults
    void addFault(std::unique_ptr<SeismicFaultModel> fault);
    void addFault(const std::map<std::string, std::string>& config);
    
    // Configure all faults from config sections
    void configure(const std::vector<std::map<std::string, std::string>>& fault_configs);
    
    // Update all faults given regional stress
    void update(double sxx, double syy, double szz,
               double sxy, double sxz, double syz,
               double pore_pressure, double shear_modulus, double dt);
    
    // Cascade: include stress transfer from slip events
    void updateWithCascade(double sxx, double syy, double szz,
                          double sxy, double sxz, double syz,
                          double pore_pressure, double shear_modulus, double dt);
    
    // Access faults
    SeismicFaultModel* getFault(size_t index);
    SeismicFaultModel* getFault(const std::string& name);
    size_t numFaults() const { return faults.size(); }
    
    // Combined event catalog
    std::vector<SeismicEvent> getAllEvents() const;
    
    // Seismicity analysis
    double getTotalMomentRate() const;
    double getMagnitudeOfCompleteness() const;
    
private:
    std::vector<std::unique_ptr<SeismicFaultModel>> faults;
    
    // Stress interaction matrix (optional precomputation)
    std::vector<std::vector<double>> stress_transfer_matrix;
    bool use_stress_transfer;
};

/**
 * @brief Factory to create friction models from configuration
 */
std::unique_ptr<FrictionModelBase> createFrictionModel(
    FrictionLaw law,
    const std::map<std::string, std::string>& config = {});

/**
 * @brief Parse friction law from string
 */
FrictionLaw parseFrictionLaw(const std::string& str);

/**
 * @brief Seismicity analysis utilities
 */
namespace SeismicityAnalysis {
    
    // Gutenberg-Richter analysis
    double estimateBValue(const std::vector<SeismicEvent>& events,
                         double Mc = 0.0);  // Mc = magnitude of completeness
    
    // Magnitude of completeness (maximum curvature method)
    double estimateMc(const std::vector<SeismicEvent>& events);
    
    // Omori aftershock decay
    double omoriRate(double t, double K, double c, double p);
    
    // ETAS model parameters fitting
    struct ETASParams {
        double mu;      // Background rate
        double K;       // Aftershock productivity
        double alpha;   // Magnitude scaling
        double c;       // Omori c-value
        double p;       // Omori p-value
    };
    ETASParams fitETAS(const std::vector<SeismicEvent>& events);
    
    // Seismicity rate change (Dieterich 1994)
    double rateChange(double stress_change, double stressing_rate,
                     double Asigma, double temperature = 300.0);
    
    // Stress transfer between faults (Coulomb)
    double coulombStressChange(double sxx, double syy, double sxy,
                              double receiver_strike, double receiver_dip,
                              double receiver_rake, double mu);
}

} // namespace FSRM

#endif // FAULT_MODEL_HPP
