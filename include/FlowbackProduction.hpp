#ifndef FLOWBACK_PRODUCTION_HPP
#define FLOWBACK_PRODUCTION_HPP

/**
 * @file FlowbackProduction.hpp
 * @brief Flowback and production modeling for fractured wells
 * 
 * ResFrac-equivalent capabilities:
 * - Post-frac flowback simulation
 * - Water recovery prediction
 * - Cleanup efficiency
 * - Production decline analysis
 * - Multi-phase fracture flow
 * - Well performance prediction
 * - EUR estimation
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <array>
#include <cmath>

namespace FSRM {

/**
 * @brief Flowback fluid composition
 */
struct FlowbackFluid {
    double oil_rate;                    ///< Oil rate (m³/s)
    double water_rate;                  ///< Water rate (m³/s)
    double gas_rate;                    ///< Gas rate at SC (m³/s)
    
    double oil_cumulative;              ///< Cumulative oil (m³)
    double water_cumulative;            ///< Cumulative water (m³)
    double gas_cumulative;              ///< Cumulative gas at SC (m³)
    
    // Water source analysis
    double frac_water_fraction;         ///< Fraction from frac fluid
    double formation_water_fraction;    ///< Fraction from formation
    double condensate_water_fraction;   ///< Fraction from condensation
    
    // Quality metrics
    double total_dissolved_solids;      ///< TDS (ppm)
    double chloride_concentration;      ///< Cl- (ppm)
    double salinity;                    ///< Salinity (%)
    
    FlowbackFluid();
};

/**
 * @brief Fracture state during flowback
 */
struct FractureFlowbackState {
    double pressure;                    ///< Fracture pressure (Pa)
    double oil_saturation;              ///< Oil saturation in fracture
    double water_saturation;            ///< Water saturation in fracture
    double gas_saturation;              ///< Gas saturation in fracture
    
    // Proppant pack state
    double effective_permeability;      ///< Current permeability (m²)
    double effective_width;             ///< Current width after closure (m)
    double conductivity;                ///< Current conductivity (mD·ft)
    
    // Fluid distribution
    std::vector<double> pressure_profile;  ///< Pressure along fracture
    std::vector<double> sat_w_profile;     ///< Water saturation profile
    
    FractureFlowbackState();
};

/**
 * @brief Flowback simulator
 */
class FlowbackSimulator {
public:
    FlowbackSimulator();
    
    /**
     * @brief Initialize from completed treatment
     * 
     * @param fracture_half_length Propped half-length (m)
     * @param fracture_height Propped height (m)
     * @param fracture_width Propped width (m)
     * @param num_clusters Number of perforation clusters
     * @param total_fluid_pumped Total fluid volume pumped (m³)
     * @param total_proppant Total proppant mass (kg)
     */
    void initializeFromTreatment(double fracture_half_length,
                                  double fracture_height,
                                  double fracture_width,
                                  int num_clusters,
                                  double total_fluid_pumped,
                                  double total_proppant);
    
    /**
     * @brief Set reservoir properties
     */
    void setReservoirProperties(double porosity, double permeability,
                                 double initial_pressure, double temperature);
    
    /**
     * @brief Set fluid properties
     */
    void setFluidProperties(double oil_viscosity, double water_viscosity,
                            double oil_density, double water_density,
                            double oil_fvf, double water_fvf);
    
    /**
     * @brief Set proppant properties
     */
    void setProppantProperties(double proppant_diameter,
                                double proppant_permeability,
                                double stress_sensitivity);
    
    /**
     * @brief Set well operating conditions
     * 
     * @param mode "BHP" or "RATE"
     * @param target BHP (Pa) or rate (m³/s)
     */
    void setWellControl(const std::string& mode, double target);
    
    /**
     * @brief Set choke schedule
     * 
     * @param times Time points (s)
     * @param choke_sizes Choke sizes (m² or fraction open)
     */
    void setChokeSchedule(const std::vector<double>& times,
                          const std::vector<double>& choke_sizes);
    
    /**
     * @brief Run flowback simulation
     * 
     * @param duration Total flowback duration (s)
     * @param dt Time step (s)
     * @return Flowback history
     */
    struct FlowbackHistory {
        std::vector<double> time;
        std::vector<double> oil_rate;
        std::vector<double> water_rate;
        std::vector<double> gas_rate;
        std::vector<double> bhp;
        std::vector<double> whp;
        std::vector<double> cumulative_oil;
        std::vector<double> cumulative_water;
        std::vector<double> cumulative_gas;
        std::vector<double> water_recovery;
        std::vector<double> gor;
        std::vector<double> watercut;
    };
    
    FlowbackHistory simulate(double duration, double dt);
    
    /**
     * @brief Get current flowback state
     */
    const FlowbackFluid& getCurrentFluid() const { return current_fluid_; }
    
    /**
     * @brief Get fracture state
     */
    const FractureFlowbackState& getFractureState() const { return frac_state_; }
    
    /**
     * @brief Estimate water recovery
     * 
     * @return Expected total water recovery fraction (0-1)
     */
    double estimateWaterRecovery() const;
    
    /**
     * @brief Get cleanup efficiency
     * 
     * @return Current cleanup efficiency (0-1)
     */
    double getCleanupEfficiency() const;
    
private:
    // Fracture geometry
    double half_length_;
    double height_;
    double width_;
    int num_clusters_;
    
    // Treatment data
    double total_fluid_pumped_;
    double total_proppant_;
    double fluid_in_fracture_;          ///< Remaining frac fluid (m³)
    
    // Reservoir properties
    double reservoir_porosity_;
    double reservoir_permeability_;
    double reservoir_pressure_;
    double reservoir_temperature_;
    
    // Fluid properties
    double oil_viscosity_, water_viscosity_;
    double oil_density_, water_density_;
    double oil_fvf_, water_fvf_;
    
    // Proppant properties
    double proppant_diameter_;
    double proppant_permeability_;
    double stress_sensitivity_;
    
    // Well control
    std::string control_mode_;
    double control_target_;
    std::vector<double> choke_times_;
    std::vector<double> choke_sizes_;
    
    // State
    FlowbackFluid current_fluid_;
    FractureFlowbackState frac_state_;
    double current_time_;
    
    // Solver methods
    void updateFractureState(double dt);
    void calculateFlowRates(double bhp, double& oil_rate, 
                           double& water_rate, double& gas_rate);
    double calculateBHP(double whp, double liquid_rate, double gas_rate);
};

/**
 * @brief Production decline analysis
 */
class DeclineAnalysis {
public:
    DeclineAnalysis();
    
    /**
     * @brief Decline curve types
     */
    enum class DeclineType {
        EXPONENTIAL,                    ///< q = qi * exp(-Di*t)
        HYPERBOLIC,                     ///< q = qi / (1 + b*Di*t)^(1/b)
        HARMONIC,                       ///< q = qi / (1 + Di*t)
        ARPS_HYPERBOLIC,               ///< General Arps hyperbolic
        DUONG,                          ///< Duong model for tight oil
        STRETCHED_EXPONENTIAL,         ///< Stretched exponential
        MODIFIED_HYPERBOLIC            ///< Hyperbolic to exponential
    };
    
    /**
     * @brief Fit decline curve to production data
     * 
     * @param times Time points (days)
     * @param rates Production rates (m³/day)
     * @param type Decline curve type
     * @return Fitted parameters {qi, Di, b (if applicable)}
     */
    std::vector<double> fitDeclineCurve(
        const std::vector<double>& times,
        const std::vector<double>& rates,
        DeclineType type);
    
    /**
     * @brief Calculate EUR (Estimated Ultimate Recovery)
     * 
     * @param qi Initial rate (m³/day)
     * @param Di Initial decline rate (1/day)
     * @param b Decline exponent (for hyperbolic)
     * @param type Decline type
     * @param economic_limit Economic limit rate (m³/day)
     * @return EUR (m³)
     */
    double calculateEUR(double qi, double Di, double b,
                        DeclineType type, double economic_limit);
    
    /**
     * @brief Calculate cumulative production
     * 
     * @param qi Initial rate
     * @param Di Initial decline rate
     * @param b Decline exponent
     * @param type Decline type
     * @param time Time (days)
     * @return Cumulative production (m³)
     */
    double calculateCumulative(double qi, double Di, double b,
                                DeclineType type, double time);
    
    /**
     * @brief Forecast production rate
     * 
     * @param qi Initial rate
     * @param Di Initial decline rate
     * @param b Decline exponent
     * @param type Decline type
     * @param time Time (days)
     * @return Rate at time (m³/day)
     */
    double forecastRate(double qi, double Di, double b,
                        DeclineType type, double time);
    
    /**
     * @brief Generate production forecast
     * 
     * @param qi Initial rate
     * @param Di Initial decline rate
     * @param b Decline exponent
     * @param type Decline type
     * @param forecast_years Years to forecast
     * @return Time series of rates
     */
    std::pair<std::vector<double>, std::vector<double>> generateForecast(
        double qi, double Di, double b,
        DeclineType type, double forecast_years);
    
    /**
     * @brief Calculate type curve match
     * 
     * @param times Actual time data
     * @param rates Actual rate data
     * @param type_curve Type curve data
     * @return Match quality R²
     */
    double calculateTypeCurveMatch(
        const std::vector<double>& times,
        const std::vector<double>& rates,
        const std::pair<std::vector<double>, std::vector<double>>& type_curve);
    
private:
    // Optimization for curve fitting
    double objectiveFunction(const std::vector<double>& params,
                              const std::vector<double>& times,
                              const std::vector<double>& rates,
                              DeclineType type);
    
    std::vector<double> optimizeParameters(
        const std::vector<double>& times,
        const std::vector<double>& rates,
        DeclineType type);
};

/**
 * @brief Multi-fractured horizontal well production model
 */
class MFHWProductionModel {
public:
    MFHWProductionModel();
    
    /**
     * @brief Set well geometry
     * 
     * @param lateral_length Total lateral length (m)
     * @param num_stages Number of fracture stages
     * @param fractures_per_stage Fractures per stage
     */
    void setWellGeometry(double lateral_length, int num_stages,
                          int fractures_per_stage);
    
    /**
     * @brief Set fracture properties
     * 
     * @param half_length Average half-length (m)
     * @param height Average height (m)
     * @param conductivity Average conductivity (mD·ft)
     * @param spacing Fracture spacing (m)
     */
    void setFractureProperties(double half_length, double height,
                                double conductivity, double spacing);
    
    /**
     * @brief Set reservoir properties
     */
    void setReservoirProperties(double porosity, double permeability,
                                 double thickness, double initial_pressure,
                                 double temperature);
    
    /**
     * @brief Set fluid properties
     */
    void setFluidProperties(double oil_viscosity, double oil_fvf,
                            double total_compressibility);
    
    /**
     * @brief Calculate initial production rate (PI method)
     * 
     * @param drawdown Initial drawdown (Pa)
     * @return Initial rate (m³/s)
     */
    double calculateInitialRate(double drawdown);
    
    /**
     * @brief Calculate productivity index
     * 
     * @return PI (m³/s/Pa)
     */
    double calculatePI();
    
    /**
     * @brief Generate rate-time forecast
     * 
     * @param flowing_bhp Flowing BHP (Pa)
     * @param duration Forecast duration (s)
     * @param dt Time step (s)
     * @return {times, rates}
     */
    std::pair<std::vector<double>, std::vector<double>> forecastProduction(
        double flowing_bhp, double duration, double dt);
    
    /**
     * @brief Calculate drainage area
     * 
     * @return Effective drainage area (m²)
     */
    double calculateDrainageArea();
    
    /**
     * @brief Calculate OOIP/OGIP in drainage volume
     * 
     * @return OOIP or OGIP (m³ at standard conditions)
     */
    double calculateOOIP();
    
    /**
     * @brief Calculate recovery factor
     * 
     * @param cumulative_production Cumulative production (m³)
     * @return Recovery factor (0-1)
     */
    double calculateRecoveryFactor(double cumulative_production);
    
private:
    // Well geometry
    double lateral_length_;
    int num_stages_;
    int fracs_per_stage_;
    int total_fracs_;
    
    // Fracture properties
    double frac_half_length_;
    double frac_height_;
    double frac_conductivity_;
    double frac_spacing_;
    
    // Reservoir properties
    double porosity_;
    double permeability_;
    double thickness_;
    double initial_pressure_;
    double temperature_;
    
    // Fluid properties
    double oil_viscosity_;
    double oil_fvf_;
    double total_compressibility_;
    
    // Dimensionless fracture conductivity
    double calculateFcd();
    
    // Linear flow coefficient
    double calculateLinearFlowCoeff();
    
    // Time to boundary effects
    double calculateTimeToBoundary();
};

/**
 * @brief Type curve generator
 */
class TypeCurveGenerator {
public:
    TypeCurveGenerator();
    
    /**
     * @brief Generate Fetkovich type curves
     * 
     * @param re_rw_ratios Radius ratios for family
     * @return Type curve family {tD, qD} for each ratio
     */
    std::vector<std::pair<std::vector<double>, std::vector<double>>>
    generateFetkovichCurves(const std::vector<double>& re_rw_ratios);
    
    /**
     * @brief Generate Blasingame type curves
     */
    std::vector<std::pair<std::vector<double>, std::vector<double>>>
    generateBlasingameCurves(const std::vector<double>& params);
    
    /**
     * @brief Generate rate-normalized type curves
     */
    std::vector<std::pair<std::vector<double>, std::vector<double>>>
    generateRateNormalizedCurves();
    
    /**
     * @brief Generate MFHW-specific type curves
     */
    std::vector<std::pair<std::vector<double>, std::vector<double>>>
    generateMFHWTypeCurves(const std::vector<double>& fcd_values);
    
private:
    int num_points_;
    double tD_min_, tD_max_;
};

} // namespace FSRM

#endif // FLOWBACK_PRODUCTION_HPP
