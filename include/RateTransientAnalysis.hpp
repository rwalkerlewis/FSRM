#ifndef RATE_TRANSIENT_ANALYSIS_HPP
#define RATE_TRANSIENT_ANALYSIS_HPP

/**
 * @file RateTransientAnalysis.hpp
 * @brief Rate Transient Analysis (RTA) for unconventional wells
 * 
 * ResFrac-equivalent capabilities:
 * - Diagnostic plots (log-log, specialized)
 * - Flow regime identification
 * - Pressure normalized rate analysis
 * - Square root of time analysis
 * - Blasingame analysis
 * - Flowing material balance
 * - Fracture/reservoir property estimation
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
 * @brief Flow regime types
 */
enum class FlowRegime {
    WELLBORE_STORAGE,                   ///< Early-time storage effects
    LINEAR_FRACTURE,                    ///< Linear flow in fracture
    BILINEAR,                           ///< Bilinear flow (finite conductivity)
    LINEAR_FORMATION,                   ///< Linear flow in formation
    ELLIPTICAL,                         ///< Elliptical/transitional
    RADIAL,                             ///< Radial flow
    PSEUDO_STEADY_STATE,               ///< Boundary-dominated
    UNKNOWN                             ///< Cannot identify
};

/**
 * @brief Production data point
 */
struct ProductionDataPoint {
    double time;                        ///< Time (s)
    double oil_rate;                    ///< Oil rate (m³/s)
    double water_rate;                  ///< Water rate (m³/s)
    double gas_rate;                    ///< Gas rate (m³/s at SC)
    double bhp;                         ///< Bottomhole pressure (Pa)
    double whp;                         ///< Wellhead pressure (Pa)
    double cumulative_oil;              ///< Cumulative oil (m³)
    double cumulative_water;            ///< Cumulative water (m³)
    double cumulative_gas;              ///< Cumulative gas (m³)
    
    ProductionDataPoint();
};

/**
 * @brief Diagnostic plot data
 */
struct DiagnosticPlot {
    std::string name;                   ///< Plot name
    std::string x_label;                ///< X-axis label
    std::string y_label;                ///< Y-axis label
    std::vector<double> x_data;         ///< X values
    std::vector<double> y_data;         ///< Y values
    std::vector<double> derivative;     ///< Derivative (if applicable)
    
    // Flow regime markers
    std::vector<std::pair<double, double>> regime_boundaries;
    std::vector<FlowRegime> regimes;
};

/**
 * @brief Rate Transient Analysis engine
 */
class RTAEngine {
public:
    RTAEngine();
    
    /**
     * @brief Load production data
     * 
     * @param data Vector of production data points
     */
    void loadData(const std::vector<ProductionDataPoint>& data);
    
    /**
     * @brief Set reservoir properties (for analysis)
     */
    void setReservoirProperties(double porosity, double permeability,
                                 double thickness, double initial_pressure,
                                 double total_compressibility);
    
    /**
     * @brief Set fluid properties
     */
    void setFluidProperties(double oil_viscosity, double oil_fvf);
    
    /**
     * @brief Set well properties
     */
    void setWellProperties(double wellbore_radius, int num_fractures,
                           double fracture_half_length);
    
    // ================== Diagnostic Plots ==================
    
    /**
     * @brief Generate log-log rate-time plot
     * 
     * @return Diagnostic plot data
     */
    DiagnosticPlot generateLogLogPlot();
    
    /**
     * @brief Generate pressure-normalized rate plot
     * 
     * Plots q/(pi-pwf) vs time
     */
    DiagnosticPlot generatePressureNormalizedPlot();
    
    /**
     * @brief Generate square root of time plot
     * 
     * For linear flow regime identification
     */
    DiagnosticPlot generateSqrtTimePlot();
    
    /**
     * @brief Generate fourth root of time plot
     * 
     * For bilinear flow regime identification
     */
    DiagnosticPlot generateFourthRootTimePlot();
    
    /**
     * @brief Generate Blasingame plot
     * 
     * Rate-normalized cumulative production
     */
    DiagnosticPlot generateBlasingamePlot();
    
    /**
     * @brief Generate Fetkovich plot
     */
    DiagnosticPlot generateFetkovichPlot();
    
    /**
     * @brief Generate Agarwal-Gardner plot
     */
    DiagnosticPlot generateAgarwalGardnerPlot();
    
    /**
     * @brief Generate NPI (Normalized Pressure Integral) plot
     */
    DiagnosticPlot generateNPIPlot();
    
    // ================== Flow Regime Identification ==================
    
    /**
     * @brief Identify current flow regime
     * 
     * @param time Current time (s)
     * @return Identified flow regime
     */
    FlowRegime identifyFlowRegime(double time);
    
    /**
     * @brief Get all flow regime transitions
     * 
     * @return Vector of {time, regime} pairs
     */
    std::vector<std::pair<double, FlowRegime>> getFlowRegimeHistory();
    
    /**
     * @brief Calculate derivative for regime identification
     * 
     * Uses Bourdet derivative algorithm
     * 
     * @param x_data X values (typically time)
     * @param y_data Y values
     * @param L Smoothing parameter
     * @return Derivative values
     */
    std::vector<double> calculateBourdetDerivative(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        double L = 0.2);
    
    // ================== Property Estimation ==================
    
    /**
     * @brief Estimate fracture half-length from linear flow
     * 
     * Uses slope of sqrt(t) plot
     * 
     * @return Estimated half-length (m)
     */
    double estimateFractureHalfLength();
    
    /**
     * @brief Estimate fracture conductivity from bilinear flow
     * 
     * @return Estimated conductivity (mD·ft)
     */
    double estimateFractureConductivity();
    
    /**
     * @brief Estimate permeability from radial flow
     * 
     * @return Estimated permeability (m²)
     */
    double estimatePermeability();
    
    /**
     * @brief Estimate skin factor
     * 
     * @return Estimated skin
     */
    double estimateSkin();
    
    /**
     * @brief Estimate drainage area
     * 
     * From pseudo-steady state analysis
     * 
     * @return Drainage area (m²)
     */
    double estimateDrainageArea();
    
    /**
     * @brief Estimate OOIP from flowing material balance
     * 
     * @return Estimated OOIP (m³)
     */
    double estimateOOIP();
    
    /**
     * @brief Full property estimation
     * 
     * @return Map of property name to estimated value
     */
    std::map<std::string, double> estimateAllProperties();
    
    // ================== Flowing Material Balance ==================
    
    /**
     * @brief Calculate material balance time
     * 
     * tc = Np / q (cumulative / rate)
     * 
     * @return Vector of material balance times
     */
    std::vector<double> calculateMaterialBalanceTime();
    
    /**
     * @brief Generate flowing material balance plot
     * 
     * p/z vs Gp for gas wells
     * p vs Np for oil wells
     */
    DiagnosticPlot generateFlowingMaterialBalancePlot();
    
    /**
     * @brief Calculate average reservoir pressure
     * 
     * From flowing material balance
     * 
     * @return Estimated average pressure (Pa)
     */
    double estimateAverageReservoirPressure();
    
    // ================== Type Curve Matching ==================
    
    /**
     * @brief Match to type curve
     * 
     * @param type_curve Reference type curve
     * @return Match parameters and quality
     */
    struct TypeCurveMatch {
        double tD_match;                ///< Dimensionless time match
        double qD_match;                ///< Dimensionless rate match
        double quality;                 ///< Match quality (R²)
        std::map<std::string, double> derived_properties;
    };
    
    TypeCurveMatch matchTypeCurve(
        const std::pair<std::vector<double>, std::vector<double>>& type_curve);
    
    // ================== Output ==================
    
    /**
     * @brief Generate analysis report
     */
    void generateReport(std::ostream& os) const;
    
    /**
     * @brief Export diagnostic plots to file
     */
    void exportPlots(const std::string& directory) const;
    
private:
    // Production data
    std::vector<ProductionDataPoint> data_;
    
    // Reservoir properties
    double porosity_;
    double permeability_;
    double thickness_;
    double initial_pressure_;
    double total_compressibility_;
    
    // Fluid properties
    double oil_viscosity_;
    double oil_fvf_;
    
    // Well properties
    double wellbore_radius_;
    int num_fractures_;
    double fracture_half_length_;
    
    // Derived data
    std::vector<double> pressure_normalized_rate_;
    std::vector<double> material_balance_time_;
    std::vector<double> sqrt_time_;
    std::vector<double> fourth_root_time_;
    
    // Analysis results
    std::map<std::string, double> estimated_properties_;
    std::vector<std::pair<double, FlowRegime>> regime_history_;
    
    // Helper functions
    void calculateDerivedData();
    double linearRegression(const std::vector<double>& x,
                            const std::vector<double>& y,
                            double& slope, double& intercept);
    FlowRegime classifySlope(double slope, const std::string& plot_type);
};

/**
 * @brief DFIT (Diagnostic Fracture Injection Test) analyzer
 */
class DFITAnalysis {
public:
    DFITAnalysis();
    
    /**
     * @brief Load DFIT pressure data
     * 
     * @param times Time points (s)
     * @param pressures Pressures (Pa)
     */
    void loadData(const std::vector<double>& times,
                  const std::vector<double>& pressures);
    
    /**
     * @brief Set injection parameters
     */
    void setInjectionParameters(double volume, double rate, double duration);
    
    /**
     * @brief Generate G-function plot
     */
    DiagnosticPlot generateGFunctionPlot();
    
    /**
     * @brief Generate square root of time plot
     */
    DiagnosticPlot generateSqrtTimePlot();
    
    /**
     * @brief Generate log-log diagnostic plot
     */
    DiagnosticPlot generateLogLogPlot();
    
    /**
     * @brief Identify closure pressure
     * 
     * Using G-function derivative
     * 
     * @return Closure pressure (Pa)
     */
    double identifyClosurePressure();
    
    /**
     * @brief Calculate ISIP
     * 
     * @return Instantaneous shut-in pressure (Pa)
     */
    double calculateISIP();
    
    /**
     * @brief Estimate process zone stress
     * 
     * @return Process zone stress (Pa)
     */
    double estimateProcessZoneStress();
    
    /**
     * @brief Estimate minimum horizontal stress
     * 
     * @return Shmin (Pa)
     */
    double estimateMinHorizontalStress();
    
    /**
     * @brief Estimate permeability from after-closure analysis
     * 
     * @return Permeability (m²)
     */
    double estimatePermeability();
    
    /**
     * @brief Estimate initial reservoir pressure
     * 
     * From extrapolation to infinite time
     * 
     * @return Initial pressure (Pa)
     */
    double estimateInitialPressure();
    
    /**
     * @brief Get all DFIT results
     */
    struct DFITResults {
        double isip;
        double closure_pressure;
        double process_zone_stress;
        double min_horizontal_stress;
        double permeability;
        double initial_pressure;
        double leakoff_coefficient;
        std::string closure_mechanism;  ///< "tip extension", "fracture compliance", etc.
    };
    
    DFITResults analyze();
    
private:
    std::vector<double> times_;
    std::vector<double> pressures_;
    
    double injection_volume_;
    double injection_rate_;
    double injection_duration_;
    double shut_in_time_;
    
    // G-function calculation
    std::vector<double> calculateGFunction();
    std::vector<double> calculateGFunctionDerivative();
    
    // Superposition time
    std::vector<double> calculateSuperpositionTime();
};

/**
 * @brief Pressure buildup/falloff analysis
 */
class PressureTransientAnalysis {
public:
    PressureTransientAnalysis();
    
    /**
     * @brief Load pressure data
     */
    void loadData(const std::vector<double>& times,
                  const std::vector<double>& pressures);
    
    /**
     * @brief Set production history before test
     */
    void setProductionHistory(const std::vector<double>& times,
                               const std::vector<double>& rates);
    
    /**
     * @brief Generate Horner plot
     */
    DiagnosticPlot generateHornerPlot();
    
    /**
     * @brief Generate MDH plot
     */
    DiagnosticPlot generateMDHPlot();
    
    /**
     * @brief Generate superposition plot
     */
    DiagnosticPlot generateSuperpositionPlot();
    
    /**
     * @brief Estimate permeability
     */
    double estimatePermeability();
    
    /**
     * @brief Estimate skin
     */
    double estimateSkin();
    
    /**
     * @brief Estimate initial pressure
     */
    double estimateInitialPressure();
    
    /**
     * @brief Estimate wellbore storage coefficient
     */
    double estimateWellboreStorage();
    
private:
    std::vector<double> test_times_;
    std::vector<double> test_pressures_;
    std::vector<double> prod_times_;
    std::vector<double> prod_rates_;
    
    double calculateHornerTime(double delta_t, double tp);
    double calculateSuperpositionTime(double delta_t);
};

} // namespace FSRM

#endif // RATE_TRANSIENT_ANALYSIS_HPP
