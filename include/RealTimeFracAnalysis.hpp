#ifndef REAL_TIME_FRAC_ANALYSIS_HPP
#define REAL_TIME_FRAC_ANALYSIS_HPP

/**
 * @file RealTimeFracAnalysis.hpp
 * @brief Real-time hydraulic fracturing treatment analysis
 * 
 * ResFrac-equivalent capabilities:
 * - ISIP analysis
 * - Closure pressure detection
 * - Net pressure analysis
 * - Friction analysis (pipe, perforation, near-wellbore)
 * - Fracture geometry estimation
 * - Treatment efficiency
 * - Real-time model calibration
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <array>
#include <deque>
#include <cmath>

namespace FSRM {

/**
 * @brief Real-time treatment data point
 */
struct TreatmentDataPoint {
    double time;                        ///< Time from stage start (s)
    double surface_pressure;            ///< Surface treating pressure (Pa)
    double slurry_rate;                 ///< Slurry pump rate (m³/s)
    double proppant_concentration;      ///< Proppant concentration (kg/m³)
    double cumulative_slurry;           ///< Cumulative slurry volume (m³)
    double cumulative_proppant;         ///< Cumulative proppant mass (kg)
    double fluid_density;               ///< Current fluid density (kg/m³)
    
    // Optional measurements
    double bottomhole_pressure;         ///< BHP if measured (Pa), -1 if not
    double annulus_pressure;            ///< Annulus pressure (Pa)
    double surface_temperature;         ///< Surface temperature (K)
    
    TreatmentDataPoint();
    
    bool hasBHP() const { return bottomhole_pressure > 0; }
};

/**
 * @brief Treatment event
 */
struct TreatmentEvent {
    enum class EventType {
        BREAKDOWN,                      ///< Formation breakdown
        ISIP,                           ///< Instantaneous shut-in
        CLOSURE,                        ///< Fracture closure
        REOPENING,                      ///< Fracture reopening
        SCREENOUT,                      ///< Screenout detected
        RATE_CHANGE,                    ///< Pump rate change
        FLUID_CHANGE,                   ///< Fluid system change
        PROPPANT_START,                 ///< Start of proppant
        FLUSH_START,                    ///< Start of flush
        BALL_DROP                       ///< Ball/plug drop
    };
    
    EventType type;
    double time;                        ///< Time of event (s)
    double pressure;                    ///< Pressure at event (Pa)
    double rate;                        ///< Rate at event (m³/s)
    std::string description;            ///< Event description
    
    TreatmentEvent();
};

/**
 * @brief Friction component breakdown
 */
struct FrictionBreakdown {
    double pipe_friction;               ///< Pipe friction (Pa)
    double perforation_friction;        ///< Perforation friction (Pa)
    double near_wellbore_friction;      ///< Near-wellbore tortuosity (Pa)
    double total_friction;              ///< Total friction (Pa)
    
    // Friction per unit length or per perf
    double pipe_friction_gradient;      ///< Pa/m
    double perf_friction_per_perf;      ///< Pa per perforation
    
    FrictionBreakdown();
};

/**
 * @brief Net pressure analysis results
 */
struct NetPressureAnalysis {
    double bhp;                         ///< Bottomhole pressure (Pa)
    double closure_pressure;            ///< Closure pressure (Pa)
    double net_pressure;                ///< Net pressure = BHP - Pc (Pa)
    double process_zone_stress;         ///< Stress in process zone (Pa)
    
    // Net pressure trend
    std::string trend;                  ///< "increasing", "constant", "decreasing"
    double trend_slope;                 ///< Rate of change (Pa/s)
    
    // Interpretation
    std::string interpretation;         ///< e.g., "confined height", "height growth"
    
    NetPressureAnalysis();
};

/**
 * @brief Estimated fracture geometry from treatment
 */
struct EstimatedFractureGeometry {
    double half_length;                 ///< Estimated half-length (m)
    double height;                      ///< Estimated height (m)
    double average_width;               ///< Estimated average width (m)
    double max_width;                   ///< Estimated max width (m)
    double propped_length;              ///< Estimated propped length (m)
    double propped_height;              ///< Estimated propped height (m)
    
    // Confidence intervals
    double half_length_uncertainty;
    double height_uncertainty;
    
    // Method used
    std::string estimation_method;      ///< "material balance", "net pressure", "PKN", etc.
    
    EstimatedFractureGeometry();
};

/**
 * @brief ISIP analyzer
 */
class ISIPAnalyzer {
public:
    ISIPAnalyzer();
    
    /**
     * @brief Analyze shut-in for ISIP
     * 
     * @param times Time points during shut-in (s)
     * @param pressures Pressure points (Pa)
     * @return Identified ISIP (Pa)
     */
    double analyzeISIP(const std::vector<double>& times,
                       const std::vector<double>& pressures);
    
    /**
     * @brief Calculate ISIP using Horner method
     */
    double hornerISIP(const std::vector<double>& times,
                      const std::vector<double>& pressures,
                      double pumping_time);
    
    /**
     * @brief Get pressure decline rate at ISIP
     * 
     * @return Decline rate (Pa/s)
     */
    double getDeclineRate() const { return decline_rate_; }
    
    /**
     * @brief Estimate fracture gradient from ISIP
     * 
     * @param tvd True vertical depth (m)
     * @return Fracture gradient (Pa/m)
     */
    double estimateFractureGradient(double tvd) const;
    
private:
    double isip_;
    double decline_rate_;
    double identification_time_;
    
    // Detection algorithms
    double detectISIPFromDerivative(const std::vector<double>& times,
                                     const std::vector<double>& pressures);
};

/**
 * @brief Closure pressure analyzer
 */
class ClosureAnalyzer {
public:
    ClosureAnalyzer();
    
    /**
     * @brief Analyze for closure pressure
     * 
     * @param times Time points (s)
     * @param pressures Pressure points (Pa)
     * @return Closure pressure (Pa)
     */
    double analyzeClosure(const std::vector<double>& times,
                          const std::vector<double>& pressures);
    
    /**
     * @brief G-function analysis
     * 
     * @param times Time points (s)
     * @param pressures Pressure points (Pa)
     * @param pumping_time Total pumping time (s)
     * @return Closure from G-function (Pa)
     */
    double gFunctionClosure(const std::vector<double>& times,
                            const std::vector<double>& pressures,
                            double pumping_time);
    
    /**
     * @brief Square root of time analysis
     */
    double sqrtTimeClosure(const std::vector<double>& times,
                           const std::vector<double>& pressures);
    
    /**
     * @brief Log-log analysis
     */
    double logLogClosure(const std::vector<double>& times,
                         const std::vector<double>& pressures);
    
    /**
     * @brief Get closure mechanism interpretation
     */
    std::string getClosureMechanism() const { return closure_mechanism_; }
    
    /**
     * @brief Get closure time
     */
    double getClosureTime() const { return closure_time_; }
    
private:
    double closure_pressure_;
    double closure_time_;
    std::string closure_mechanism_;
    
    // G-function calculation
    std::vector<double> calculateGFunction(const std::vector<double>& times,
                                            double pumping_time);
};

/**
 * @brief Real-time friction analyzer
 */
class FrictionAnalyzer {
public:
    FrictionAnalyzer();
    
    /**
     * @brief Set wellbore properties for friction calculation
     */
    void setWellboreProperties(double md, double tvd,
                                double pipe_id, double roughness);
    
    /**
     * @brief Set perforation properties
     */
    void setPerforationProperties(int num_perfs, double perf_diameter,
                                   double discharge_coeff);
    
    /**
     * @brief Calculate friction breakdown
     * 
     * @param surface_pressure Surface pressure (Pa)
     * @param rate Pump rate (m³/s)
     * @param fluid_density Fluid density (kg/m³)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @param bhp_measured BHP if measured (Pa), -1 if not
     * @return Friction breakdown
     */
    FrictionBreakdown calculateFriction(double surface_pressure,
                                         double rate,
                                         double fluid_density,
                                         double fluid_viscosity,
                                         double bhp_measured = -1);
    
    /**
     * @brief Calibrate friction from step-down test
     * 
     * @param rates Rates during step-down (m³/s)
     * @param pressures Surface pressures (Pa)
     * @return Calibrated friction parameters
     */
    void calibrateFromStepDown(const std::vector<double>& rates,
                                const std::vector<double>& pressures);
    
    /**
     * @brief Track perforation erosion
     * 
     * @param dt Time step (s)
     * @param rate Flow rate (m³/s)
     * @param proppant_conc Proppant concentration (kg/m³)
     */
    void updatePerfErosion(double dt, double rate, double proppant_conc);
    
    /**
     * @brief Get current perforation diameter
     */
    double getCurrentPerfDiameter() const { return current_perf_diameter_; }
    
private:
    // Wellbore
    double md_, tvd_;
    double pipe_id_, roughness_;
    
    // Perforations
    int num_perfs_;
    double initial_perf_diameter_;
    double current_perf_diameter_;
    double discharge_coeff_;
    
    // Calibration
    double pipe_friction_multiplier_;
    double perf_friction_multiplier_;
    double nwb_friction_;
};

/**
 * @brief Net pressure analyzer
 */
class NetPressureAnalyzer {
public:
    NetPressureAnalyzer();
    
    /**
     * @brief Set closure pressure
     */
    void setClosurePressure(double closure) { closure_pressure_ = closure; }
    
    /**
     * @brief Analyze net pressure
     * 
     * @param bhp Current BHP (Pa)
     * @param time Current time (s)
     * @return Net pressure analysis
     */
    NetPressureAnalysis analyze(double bhp, double time);
    
    /**
     * @brief Interpret net pressure trend
     * 
     * @param history Recent net pressure history
     * @return Interpretation string
     */
    std::string interpretTrend(const std::vector<std::pair<double, double>>& history);
    
    /**
     * @brief Add data point to history
     */
    void addDataPoint(double time, double net_pressure);
    
    /**
     * @brief Get net pressure history
     */
    const std::deque<std::pair<double, double>>& getHistory() const { 
        return history_; 
    }
    
    // Nolte-Smith analysis
    enum class NolteSmithMode {
        MODE_I,                         ///< Confined height growth
        MODE_II,                        ///< Height growth with tip screenout
        MODE_III,                       ///< Unconfined height growth
        MODE_IV                         ///< Tip screenout
    };
    
    NolteSmithMode identifyNolteSmithMode();
    
private:
    double closure_pressure_;
    std::deque<std::pair<double, double>> history_;  // {time, net_pressure}
    size_t max_history_size_;
    
    double calculateTrendSlope();
};

/**
 * @brief Real-time fracture geometry estimator
 */
class RealTimeGeometryEstimator {
public:
    RealTimeGeometryEstimator();
    
    /**
     * @brief Set rock and fracture properties
     */
    void setRockProperties(double youngs_modulus, double poisson_ratio,
                           double leakoff_coefficient);
    
    /**
     * @brief Set fluid properties
     */
    void setFluidProperties(double viscosity, double density);
    
    /**
     * @brief Set fracture height constraint
     */
    void setHeightConstraint(double height);
    
    /**
     * @brief Estimate geometry from current treatment data
     * 
     * @param net_pressure Current net pressure (Pa)
     * @param cumulative_volume Cumulative volume pumped (m³)
     * @param pump_time Pumping time (s)
     * @param rate Current rate (m³/s)
     * @return Estimated geometry
     */
    EstimatedFractureGeometry estimateGeometry(double net_pressure,
                                                 double cumulative_volume,
                                                 double pump_time,
                                                 double rate);
    
    /**
     * @brief Material balance estimate
     * 
     * V_frac = V_pumped - V_leakoff
     */
    EstimatedFractureGeometry materialBalanceEstimate(
        double cumulative_volume,
        double pump_time);
    
    /**
     * @brief PKN model estimate
     */
    EstimatedFractureGeometry pknEstimate(double net_pressure,
                                           double rate,
                                           double pump_time);
    
    /**
     * @brief KGD model estimate
     */
    EstimatedFractureGeometry kgdEstimate(double net_pressure,
                                           double rate,
                                           double pump_time);
    
    /**
     * @brief Update estimates with proppant placement
     * 
     * @param proppant_mass Total proppant pumped (kg)
     * @param proppant_concentration Average concentration (kg/m³)
     * @param geometry Current geometry estimate
     */
    void updateProppedGeometry(double proppant_mass,
                                double proppant_concentration,
                                EstimatedFractureGeometry& geometry);
    
private:
    // Rock properties
    double youngs_modulus_;
    double poisson_ratio_;
    double plane_strain_modulus_;
    double leakoff_coefficient_;
    
    // Fluid properties
    double viscosity_;
    double density_;
    
    // Constraints
    double height_constraint_;
    bool use_height_constraint_;
};

/**
 * @brief Comprehensive real-time analysis engine
 */
class RealTimeAnalysisEngine {
public:
    RealTimeAnalysisEngine();
    
    /**
     * @brief Initialize for new stage
     * 
     * @param stage_number Stage number
     * @param stage_start_time Stage start time (s)
     */
    void initializeStage(int stage_number, double stage_start_time);
    
    /**
     * @brief Set wellbore and completion properties
     */
    void setWellboreProperties(double md, double tvd, double pipe_id,
                                int num_perfs, double perf_diameter);
    
    /**
     * @brief Set rock and stress properties
     */
    void setRockProperties(double youngs_modulus, double poisson_ratio,
                           double closure_pressure);
    
    /**
     * @brief Process new data point
     * 
     * @param data New treatment data point
     * @return Analysis results
     */
    struct RealTimeResults {
        FrictionBreakdown friction;
        double estimated_bhp;
        double net_pressure;
        EstimatedFractureGeometry geometry;
        std::string status_message;
        std::vector<TreatmentEvent> new_events;
    };
    
    RealTimeResults processDataPoint(const TreatmentDataPoint& data);
    
    /**
     * @brief Detect events (ISIP, closure, screenout, etc.)
     */
    std::vector<TreatmentEvent> detectEvents();
    
    /**
     * @brief Get current stage summary
     */
    struct StageSummary {
        int stage_number;
        double duration;
        double total_volume;
        double total_proppant;
        double max_pressure;
        double avg_rate;
        double isip;
        double closure_pressure;
        EstimatedFractureGeometry final_geometry;
        double treatment_efficiency;
        int num_events;
    };
    
    StageSummary getStageSummary() const;
    
    /**
     * @brief Get all detected events
     */
    const std::vector<TreatmentEvent>& getEvents() const { return events_; }
    
    /**
     * @brief Get treatment history
     */
    const std::vector<TreatmentDataPoint>& getHistory() const { return history_; }
    
    /**
     * @brief Calculate treatment efficiency
     * 
     * η = V_fracture / V_pumped
     */
    double calculateTreatmentEfficiency() const;
    
    /**
     * @brief Export real-time analysis to file
     */
    void exportAnalysis(const std::string& filename) const;
    
private:
    // Stage info
    int stage_number_;
    double stage_start_time_;
    
    // Data history
    std::vector<TreatmentDataPoint> history_;
    std::vector<TreatmentEvent> events_;
    
    // Analyzers
    ISIPAnalyzer isip_analyzer_;
    ClosureAnalyzer closure_analyzer_;
    FrictionAnalyzer friction_analyzer_;
    NetPressureAnalyzer net_pressure_analyzer_;
    RealTimeGeometryEstimator geometry_estimator_;
    
    // Current state
    double current_closure_pressure_;
    bool is_pumping_;
    bool breakdown_detected_;
    
    // Event detection
    void detectBreakdown(const TreatmentDataPoint& data);
    void detectScreenout(const TreatmentDataPoint& data);
    bool checkForShutIn(const TreatmentDataPoint& prev,
                        const TreatmentDataPoint& curr);
};

} // namespace FSRM

#endif // REAL_TIME_FRAC_ANALYSIS_HPP
