#ifndef ECONOMIC_ANALYSIS_HPP
#define ECONOMIC_ANALYSIS_HPP

/**
 * @file EconomicAnalysis.hpp
 * @brief Economic analysis and optimization for completions
 * 
 * ResFrac-equivalent capabilities:
 * - NPV calculation
 * - EUR estimation
 * - Completion optimization (stages, clusters, proppant)
 * - Well spacing optimization
 * - Sensitivity analysis
 * - Monte Carlo uncertainty
 * - Economic screening
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <array>
#include <random>
#include <cmath>

namespace FSRM {

/**
 * @brief Commodity price forecast
 */
struct PriceForecast {
    std::string commodity;              ///< "OIL", "GAS", "NGL"
    std::vector<double> times;          ///< Time points (years)
    std::vector<double> prices;         ///< Prices at each time ($/bbl, $/mcf, etc.)
    
    double escalation_rate;             ///< Annual escalation rate
    double volatility;                  ///< Price volatility for Monte Carlo
    
    PriceForecast();
    
    // Get price at time
    double getPrice(double time) const;
    
    // Get random price for Monte Carlo
    double getRandomPrice(double time, std::mt19937& rng) const;
};

/**
 * @brief Cost breakdown
 */
struct CostBreakdown {
    // Drilling costs
    double drilling_cost;               ///< Total drilling cost ($)
    double cost_per_foot;               ///< $/ft drilled
    
    // Completion costs
    double completion_cost;             ///< Total completion cost ($)
    double cost_per_stage;              ///< $/stage
    double cost_per_cluster;            ///< $/cluster
    double proppant_cost_per_lb;        ///< $/lb proppant
    double fluid_cost_per_bbl;          ///< $/bbl frac fluid
    
    // Operating costs
    double fixed_opex;                  ///< Fixed operating cost ($/month)
    double variable_opex_oil;           ///< Variable opex ($/bbl oil)
    double variable_opex_gas;           ///< Variable opex ($/mcf gas)
    double water_disposal_cost;         ///< $/bbl water disposed
    
    // Taxes and royalties
    double severance_tax_rate;          ///< Severance tax (fraction)
    double royalty_rate;                ///< Royalty burden (fraction)
    double ad_valorem_rate;             ///< Ad valorem tax rate
    
    CostBreakdown();
    
    // Calculate total D&C cost
    double calculateDCCost(double lateral_length, int num_stages,
                          int clusters_per_stage, double proppant_lbs,
                          double fluid_bbls) const;
};

/**
 * @brief Production forecast
 */
struct ProductionForecast {
    std::vector<double> times;          ///< Time points (days)
    std::vector<double> oil_rate;       ///< Oil rate (bbl/day)
    std::vector<double> gas_rate;       ///< Gas rate (mcf/day)
    std::vector<double> water_rate;     ///< Water rate (bbl/day)
    std::vector<double> ngl_rate;       ///< NGL rate (bbl/day)
    
    double eur_oil;                     ///< EUR oil (bbls)
    double eur_gas;                     ///< EUR gas (mcf)
    double eur_ngl;                     ///< EUR NGL (bbls)
    
    ProductionForecast();
    
    // Get cumulative production at time
    double getCumulativeOil(double time) const;
    double getCumulativeGas(double time) const;
};

/**
 * @brief NPV calculation engine
 */
class NPVCalculator {
public:
    NPVCalculator();
    
    /**
     * @brief Set discount rate
     * 
     * @param rate Annual discount rate (e.g., 0.10 for 10%)
     */
    void setDiscountRate(double rate) { discount_rate_ = rate; }
    
    /**
     * @brief Set working interest
     * 
     * @param wi Working interest fraction (0-1)
     */
    void setWorkingInterest(double wi) { working_interest_ = wi; }
    
    /**
     * @brief Set price forecasts
     */
    void setOilPrice(const PriceForecast& forecast) { oil_price_ = forecast; }
    void setGasPrice(const PriceForecast& forecast) { gas_price_ = forecast; }
    void setNGLPrice(const PriceForecast& forecast) { ngl_price_ = forecast; }
    
    /**
     * @brief Set cost structure
     */
    void setCosts(const CostBreakdown& costs) { costs_ = costs; }
    
    /**
     * @brief Calculate NPV
     * 
     * @param production Production forecast
     * @param evaluation_years Years to evaluate
     * @return NPV in dollars
     */
    double calculateNPV(const ProductionForecast& production,
                        double evaluation_years);
    
    /**
     * @brief Calculate IRR
     * 
     * @param production Production forecast
     * @param initial_investment D&C cost
     * @param evaluation_years Evaluation period
     * @return IRR (decimal, e.g., 0.20 for 20%)
     */
    double calculateIRR(const ProductionForecast& production,
                        double initial_investment,
                        double evaluation_years);
    
    /**
     * @brief Calculate payout time
     * 
     * @param production Production forecast
     * @param initial_investment D&C cost
     * @return Payout time (years)
     */
    double calculatePayout(const ProductionForecast& production,
                           double initial_investment);
    
    /**
     * @brief Calculate ROI
     */
    double calculateROI(const ProductionForecast& production,
                        double initial_investment,
                        double evaluation_years);
    
    /**
     * @brief Get monthly cash flow
     * 
     * @param production Production forecast
     * @param initial_investment D&C cost
     * @return Vector of monthly cash flows
     */
    std::vector<double> getMonthlyCashFlow(
        const ProductionForecast& production,
        double initial_investment);
    
    /**
     * @brief Calculate breakeven price
     * 
     * @param production Production forecast
     * @param initial_investment D&C cost
     * @param commodity "OIL" or "GAS"
     * @return Breakeven price ($/bbl or $/mcf)
     */
    double calculateBreakevenPrice(const ProductionForecast& production,
                                    double initial_investment,
                                    const std::string& commodity);
    
private:
    double discount_rate_;
    double working_interest_;
    
    PriceForecast oil_price_;
    PriceForecast gas_price_;
    PriceForecast ngl_price_;
    CostBreakdown costs_;
    
    // Calculate net revenue for a month
    double calculateNetRevenue(double oil_prod, double gas_prod,
                               double water_prod, double ngl_prod,
                               double time);
};

/**
 * @brief Completion optimization
 */
class CompletionOptimizer {
public:
    CompletionOptimizer();
    
    /**
     * @brief Set economic calculator
     */
    void setEconomicCalculator(std::shared_ptr<NPVCalculator> calculator) {
        npv_calc_ = calculator;
    }
    
    /**
     * @brief Set production model
     * 
     * Function that takes completion parameters and returns production forecast
     */
    using ProductionModel = std::function<ProductionForecast(
        int num_stages, int clusters_per_stage,
        double proppant_per_stage, double fluid_per_stage)>;
    
    void setProductionModel(ProductionModel model) { prod_model_ = model; }
    
    /**
     * @brief Set lateral length
     */
    void setLateralLength(double length) { lateral_length_ = length; }
    
    /**
     * @brief Set completion constraints
     */
    struct CompletionConstraints {
        int min_stages, max_stages;
        int min_clusters, max_clusters;
        double min_proppant, max_proppant;  ///< Per stage (lbs)
        double min_fluid, max_fluid;        ///< Per stage (bbls)
        double max_total_cost;              ///< Total budget
        
        CompletionConstraints();
    };
    
    void setConstraints(const CompletionConstraints& constraints) {
        constraints_ = constraints;
    }
    
    /**
     * @brief Optimize completion for maximum NPV
     * 
     * @return Optimal completion parameters
     */
    struct OptimizationResult {
        int optimal_stages;
        int optimal_clusters;
        double optimal_proppant;            ///< Per stage (lbs)
        double optimal_fluid;               ///< Per stage (bbls)
        double total_cost;
        double npv;
        double irr;
        double payout;
        ProductionForecast forecast;
    };
    
    OptimizationResult optimizeNPV();
    
    /**
     * @brief Optimize for specific target
     * 
     * @param target "NPV", "IRR", "EUR", "PAYOUT"
     */
    OptimizationResult optimize(const std::string& target);
    
    /**
     * @brief Generate response surface
     * 
     * @param param1 First parameter name
     * @param values1 Values for first parameter
     * @param param2 Second parameter name
     * @param values2 Values for second parameter
     * @return 2D NPV surface
     */
    std::vector<std::vector<double>> generateResponseSurface(
        const std::string& param1, const std::vector<double>& values1,
        const std::string& param2, const std::vector<double>& values2);
    
private:
    std::shared_ptr<NPVCalculator> npv_calc_;
    ProductionModel prod_model_;
    double lateral_length_;
    CompletionConstraints constraints_;
    
    // Optimization methods
    OptimizationResult gridSearch();
    OptimizationResult gradientDescent();
    OptimizationResult geneticAlgorithm();
};

/**
 * @brief Sensitivity analysis
 */
class SensitivityAnalysis {
public:
    SensitivityAnalysis();
    
    /**
     * @brief Set base case parameters
     */
    void setBaseCase(const std::map<std::string, double>& params);
    
    /**
     * @brief Set economic calculator
     */
    void setEconomicCalculator(std::shared_ptr<NPVCalculator> calculator);
    
    /**
     * @brief Set production model
     */
    using ProductionModel = std::function<ProductionForecast(
        const std::map<std::string, double>& params)>;
    
    void setProductionModel(ProductionModel model);
    
    /**
     * @brief Run one-way sensitivity
     * 
     * @param parameter Parameter to vary
     * @param low_value Low case value
     * @param high_value High case value
     * @param num_points Number of points to evaluate
     * @return {parameter values, NPV values}
     */
    std::pair<std::vector<double>, std::vector<double>> oneWaySensitivity(
        const std::string& parameter,
        double low_value, double high_value, int num_points = 10);
    
    /**
     * @brief Run tornado chart analysis
     * 
     * @param parameters Parameters to analyze
     * @param low_multipliers Low case multipliers (e.g., 0.8 for -20%)
     * @param high_multipliers High case multipliers (e.g., 1.2 for +20%)
     * @return Sorted impact by parameter
     */
    struct TornadoResult {
        std::string parameter;
        double base_npv;
        double low_npv;
        double high_npv;
        double npv_swing;
    };
    
    std::vector<TornadoResult> tornadoAnalysis(
        const std::vector<std::string>& parameters,
        double low_mult = 0.8, double high_mult = 1.2);
    
    /**
     * @brief Generate spider plot data
     * 
     * @param parameters Parameters to analyze
     * @param range Percentage range to vary (+/-)
     * @return Spider plot data
     */
    std::map<std::string, std::vector<std::pair<double, double>>> spiderPlot(
        const std::vector<std::string>& parameters,
        double range = 0.3);
    
private:
    std::map<std::string, double> base_case_;
    std::shared_ptr<NPVCalculator> npv_calc_;
    ProductionModel prod_model_;
};

/**
 * @brief Monte Carlo uncertainty analysis
 */
class MonteCarloAnalysis {
public:
    MonteCarloAnalysis();
    
    /**
     * @brief Distribution type
     */
    enum class Distribution {
        NORMAL,
        LOGNORMAL,
        TRIANGULAR,
        UNIFORM
    };
    
    /**
     * @brief Add uncertain variable
     * 
     * @param name Variable name
     * @param distribution Distribution type
     * @param params Distribution parameters (mean, std for normal, etc.)
     */
    void addVariable(const std::string& name, Distribution dist,
                     const std::vector<double>& params);
    
    /**
     * @brief Set economic calculator
     */
    void setEconomicCalculator(std::shared_ptr<NPVCalculator> calculator);
    
    /**
     * @brief Set production model
     */
    using ProductionModel = std::function<ProductionForecast(
        const std::map<std::string, double>& params)>;
    
    void setProductionModel(ProductionModel model);
    
    /**
     * @brief Run Monte Carlo simulation
     * 
     * @param num_iterations Number of iterations
     * @param seed Random seed
     * @return Distribution of NPVs
     */
    struct MonteCarloResults {
        std::vector<double> npv_values;
        double p10, p50, p90;
        double mean, std_dev;
        double probability_positive;
        std::map<std::string, double> sensitivity_indices;
    };
    
    MonteCarloResults run(int num_iterations, int seed = 42);
    
    /**
     * @brief Get NPV histogram data
     */
    std::pair<std::vector<double>, std::vector<int>> getNPVHistogram(
        const std::vector<double>& npv_values, int num_bins = 50);
    
    /**
     * @brief Calculate Value at Risk
     * 
     * @param npv_values NPV distribution
     * @param confidence Confidence level (e.g., 0.95)
     * @return VaR value
     */
    double calculateVaR(const std::vector<double>& npv_values,
                        double confidence);
    
private:
    struct UncertainVariable {
        std::string name;
        Distribution distribution;
        std::vector<double> params;
    };
    
    std::vector<UncertainVariable> variables_;
    std::shared_ptr<NPVCalculator> npv_calc_;
    ProductionModel prod_model_;
    
    // Sample from distribution
    double sampleVariable(const UncertainVariable& var, std::mt19937& rng);
};

/**
 * @brief Well spacing optimizer
 */
class WellSpacingOptimizer {
public:
    WellSpacingOptimizer();
    
    /**
     * @brief Set section geometry
     * 
     * @param section_width Width of drilling unit (ft)
     * @param section_length Length of drilling unit (ft)
     */
    void setSectionGeometry(double width, double length);
    
    /**
     * @brief Set interference model
     * 
     * Function that calculates EUR reduction for given spacing
     */
    using InterferenceModel = std::function<double(double spacing, 
                                                    double fracture_half_length)>;
    
    void setInterferenceModel(InterferenceModel model);
    
    /**
     * @brief Set economic calculator
     */
    void setEconomicCalculator(std::shared_ptr<NPVCalculator> calculator);
    
    /**
     * @brief Set single well EUR
     */
    void setSingleWellEUR(double oil_eur, double gas_eur);
    
    /**
     * @brief Set drilling and completion cost per well
     */
    void setCostPerWell(double dc_cost);
    
    /**
     * @brief Optimize well spacing
     * 
     * @param fracture_half_length Fracture half-length (ft)
     * @param min_spacing Minimum well spacing (ft)
     * @param max_spacing Maximum well spacing (ft)
     * @return Optimal spacing and economics
     */
    struct SpacingResult {
        double optimal_spacing;
        int optimal_num_wells;
        double total_eur;
        double total_npv;
        double npv_per_well;
    };
    
    SpacingResult optimizeSpacing(double fracture_half_length,
                                   double min_spacing = 500,
                                   double max_spacing = 2000);
    
    /**
     * @brief Generate spacing vs NPV curve
     */
    std::pair<std::vector<double>, std::vector<double>> spacingCurve(
        double fracture_half_length,
        double min_spacing, double max_spacing, int num_points = 20);
    
private:
    double section_width_;
    double section_length_;
    InterferenceModel interference_model_;
    std::shared_ptr<NPVCalculator> npv_calc_;
    double single_well_oil_eur_;
    double single_well_gas_eur_;
    double dc_cost_per_well_;
};

} // namespace FSRM

#endif // ECONOMIC_ANALYSIS_HPP
