#ifndef FAULT_MODEL_HPP
#define FAULT_MODEL_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <array>
#include <functional>

namespace FSRM {

/**
 * @brief Types of friction laws for fault mechanics
 */
enum class FrictionLaw {
    COULOMB,                    // Static/dynamic Coulomb friction
    RATE_STATE_AGING,           // Rate-and-state with aging law
    RATE_STATE_SLIP,            // Rate-and-state with slip law
    SLIP_WEAKENING,             // Linear slip-weakening
    FLASH_HEATING,              // Flash heating at high slip rates
    THERMAL_PRESSURIZATION,     // Pore pressure effects
    STRONG_VELOCITY_WEAKENING   // Strong velocity weakening (SVW)
};

/**
 * @brief Split node methods for fault modeling
 * 
 * Split nodes duplicate nodes along fault surfaces to allow discontinuous
 * displacement fields across the fault. Different methods handle the
 * traction conditions differently.
 */
enum class SplitNodeMethod {
    LAGRANGE_MULTIPLIER,    // Enforce contact via Lagrange multipliers
    PENALTY,                // Penalty method for contact enforcement
    NITSCHE,                // Nitsche's method (consistent penalty)
    MORTAR,                 // Mortar method for non-conforming meshes
    AUGMENTED_LAGRANGIAN    // Augmented Lagrangian for robustness
};

/**
 * @brief Traction type applied to split nodes
 */
enum class TractionType {
    PRESCRIBED,             // User-specified constant traction
    FRICTION_DEPENDENT,     // Computed from friction law
    COHESIVE_ZONE,          // Cohesive zone model traction
    RATE_STATE,             // Rate-and-state dependent traction
    DYNAMIC                 // Time-varying traction (e.g., seismic source)
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
    
    /**
     * @brief Get hypocenter coordinates for dynamic event
     * 
     * If events have occurred, returns the hypocenter of the most recent event.
     * Otherwise, returns the center of the fault geometry as an approximation.
     * 
     * @param[out] hypo_x X coordinate of hypocenter
     * @param[out] hypo_y Y coordinate of hypocenter  
     * @param[out] hypo_z Z coordinate of hypocenter
     */
    void getHypocenterCoordinates(double& hypo_x, double& hypo_y, double& hypo_z) const {
        if (!events.empty()) {
            // Use hypocenter of most recent event
            const auto& last_event = events.back();
            hypo_x = last_event.hypo_x;
            hypo_y = last_event.hypo_y;
            hypo_z = last_event.hypo_z;
        } else {
            // Use center of fault geometry as approximation
            hypo_x = geometry.x;
            hypo_y = geometry.y;
            hypo_z = geometry.z;
        }
    }
    
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
 * @brief Split node pair along a fault
 * 
 * Represents a pair of coincident nodes on opposite sides of the fault.
 * The plus/minus convention refers to the two sides of the fault surface.
 */
struct SplitNodePair {
    int node_plus;           // Node index on + side of fault
    int node_minus;          // Node index on - side of fault
    double x, y, z;          // Original coordinates
    std::array<double, 3> normal;      // Fault normal at this location
    std::array<double, 3> tangent1;    // First tangent direction (strike)
    std::array<double, 3> tangent2;    // Second tangent direction (dip)
    double weight;           // Integration weight for mortar method
    
    SplitNodePair() :
        node_plus(-1), node_minus(-1),
        x(0), y(0), z(0),
        normal({0, 0, 1}),
        tangent1({1, 0, 0}),
        tangent2({0, 1, 0}),
        weight(1.0) {}
};

/**
 * @brief Traction state at a split node
 */
struct SplitNodeTraction {
    double t_n;              // Normal traction (tension positive, Pa)
    double t_s;              // Strike-direction shear traction (Pa)
    double t_d;              // Dip-direction shear traction (Pa)
    
    // Slip state
    double slip_n;           // Normal opening (m)
    double slip_s;           // Strike-slip (m)
    double slip_d;           // Dip-slip (m)
    double slip_rate;        // Total slip rate (m/s)
    
    // Contact state
    bool in_contact;         // Are nodes in contact?
    bool is_slipping;        // Is fault slipping?
    bool is_open;            // Is fault open (tensile)?
    
    SplitNodeTraction() :
        t_n(0), t_s(0), t_d(0),
        slip_n(0), slip_s(0), slip_d(0), slip_rate(0),
        in_contact(true), is_slipping(false), is_open(false) {}
    
    // Total shear traction magnitude
    double shearTraction() const { return std::sqrt(t_s * t_s + t_d * t_d); }
    
    // Total slip magnitude
    double totalSlip() const { 
        return std::sqrt(slip_s * slip_s + slip_d * slip_d); 
    }
};

/**
 * @brief Configuration for split node fault modeling
 */
struct SplitNodeConfig {
    // Method selection
    SplitNodeMethod method = SplitNodeMethod::PENALTY;
    TractionType traction_type = TractionType::FRICTION_DEPENDENT;
    
    // Penalty/Nitsche parameters
    double penalty_normal = 1e12;     // Normal penalty stiffness (Pa/m)
    double penalty_tangent = 1e10;    // Tangential penalty stiffness (Pa/m)
    double nitsche_gamma = 100.0;     // Nitsche stabilization parameter
    
    // Prescribed tractions (for TractionType::PRESCRIBED)
    double prescribed_t_n = 0.0;      // Normal traction (Pa)
    double prescribed_t_s = 0.0;      // Strike-slip traction (Pa)
    double prescribed_t_d = 0.0;      // Dip-slip traction (Pa)
    
    // Time-varying traction (for TractionType::DYNAMIC)
    std::function<void(double, double&, double&, double&)> traction_function;
    
    // Cohesive zone parameters (for TractionType::COHESIVE_ZONE)
    double cohesive_strength = 5e6;   // Peak cohesive traction (Pa)
    double critical_opening = 1e-4;   // Critical opening (m)
    double critical_slip = 1e-3;      // Critical slip (m)
    
    // Contact parameters
    double contact_tolerance = 1e-10; // Gap tolerance for contact (m)
    bool allow_separation = true;     // Allow tensile opening?
    bool symmetric_contact = true;    // Symmetric contact formulation?
    
    // Augmented Lagrangian parameters
    double aug_lag_r = 1e10;          // Augmented Lagrangian parameter
    int aug_lag_max_iter = 10;        // Max inner iterations
    double aug_lag_tol = 1e-8;        // Convergence tolerance
    
    void configure(const std::map<std::string, std::string>& config);
};

/**
 * @brief Split node fault model with traction boundary conditions
 * 
 * This class implements fault modeling using split (duplicated) nodes
 * along the fault surface. Traction conditions are applied to enforce
 * the constitutive behavior of the fault (friction, cohesion, etc.).
 * 
 * Key features:
 * - Supports multiple enforcement methods (penalty, Lagrange, Nitsche)
 * - Handles stick-slip friction with rate-state laws
 * - Tracks slip accumulation and seismic events
 * - Compatible with DMPlex mesh representation
 */
class SplitNodeFault {
public:
    SplitNodeFault();
    
    // Configuration
    void configure(const std::map<std::string, std::string>& config);
    void setGeometry(const FaultGeometry& geom);
    void setFrictionModel(std::unique_ptr<FrictionModelBase> model);
    void setSplitNodeConfig(const SplitNodeConfig& config);
    
    // Mesh operations
    void identifySplitNodes(const std::vector<double>& node_coords,
                           const std::vector<int>& fault_face_nodes);
    void createSplitNodePairs(const FaultGeometry& geom);
    void duplicateNodesAlongFault(std::vector<double>& node_coords,
                                  std::vector<int>& connectivity);
    
    // Assembly contributions
    void assembleContactResidual(const double* displacement,
                                double* residual,
                                double pore_pressure) const;
    
    void assembleContactJacobian(const double* displacement,
                                double* jacobian,
                                double pore_pressure) const;
    
    // Traction computation
    void computeTractionFromGap(const SplitNodePair& pair,
                               const double* displacement,
                               double pore_pressure,
                               SplitNodeTraction& traction) const;
    
    void applyPrescribedTraction(double time,
                                SplitNodeTraction& traction) const;
    
    void applyCohesiveZoneTraction(double opening, double slip,
                                   SplitNodeTraction& traction) const;
    
    void applyFrictionTraction(double sigma_n_eff, double slip_rate,
                              double state_var,
                              SplitNodeTraction& traction) const;
    
    // State update
    void updateSlipState(double dt);
    void updateStateVariables(double dt);
    void checkSlipEvents(double current_time);
    
    // Accessors
    const std::vector<SplitNodePair>& getSplitNodePairs() const { return node_pairs; }
    const std::vector<SplitNodeTraction>& getNodeTractions() const { return tractions; }
    size_t numSplitNodes() const { return node_pairs.size(); }
    
    // Get accumulated slip at a node
    double getSlipAtNode(size_t node_idx) const;
    double getSlipRateAtNode(size_t node_idx) const;
    
    // Get traction at a node
    void getTractionAtNode(size_t node_idx, double& t_n, double& t_s, double& t_d) const;
    
    // Event tracking
    const std::vector<SeismicEvent>& getSlipEvents() const { return slip_events; }
    void clearSlipEvents();
    
    // Enable/disable features
    void enableRateState(bool enable) { use_rate_state = enable; }
    void enableCohesion(bool enable) { use_cohesion = enable; }
    void enableDilation(bool enable) { use_dilation = enable; }
    
    std::string name;
    
private:
    // Geometry
    FaultGeometry geometry;
    
    // Node pairs and tractions
    std::vector<SplitNodePair> node_pairs;
    std::vector<SplitNodeTraction> tractions;
    std::vector<double> state_variables;  // Rate-state θ at each node
    
    // Configuration
    SplitNodeConfig config;
    std::unique_ptr<FrictionModelBase> friction;
    
    // Flags
    bool use_rate_state = false;
    bool use_cohesion = true;
    bool use_dilation = false;
    
    // Seismic events from slip
    std::vector<SeismicEvent> slip_events;
    
    // Helper methods
    void computeGapVector(const SplitNodePair& pair,
                         const double* displacement,
                         double& gap_n, double& gap_s, double& gap_d) const;
    
    double penaltyNormalForce(double gap_n) const;
    double penaltyTangentForce(double gap_t, double sigma_n) const;
    
    void nitscheTerms(const SplitNodePair& pair,
                     const double* displacement,
                     double& term_consistency,
                     double& term_penalty) const;
    
    void augmentedLagrangianUpdate(double& lambda, double gap,
                                   double penalty) const;
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
 * @brief Parse split node method from string
 */
SplitNodeMethod parseSplitNodeMethod(const std::string& str);

/**
 * @brief Parse traction type from string
 */
TractionType parseTractionType(const std::string& str);

/**
 * @brief Factory to create split node fault from configuration
 */
std::unique_ptr<SplitNodeFault> createSplitNodeFault(
    const std::map<std::string, std::string>& config);

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
