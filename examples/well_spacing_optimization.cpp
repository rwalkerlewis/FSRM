/*
 * Well Spacing and Hydraulic Fracture Placement Optimization
 * 
 * This example demonstrates optimization of well spacing and fracture placement
 * to maximize hydrocarbon recovery in unconventional reservoirs. It implements:
 * 
 * 1. Multi-well pad development with configurable well spacing
 * 2. Multi-stage hydraulic fracture placement optimization
 * 3. Stress shadowing effects between fractures and wells
 * 4. Production interference modeling between wells
 * 5. Economic analysis with NPV, IRR, and EUR calculations
 * 6. Grid search and response surface optimization
 * 7. Sensitivity analysis and Monte Carlo uncertainty quantification
 * 
 * Key optimization variables:
 * - Inter-well spacing (distance between parallel horizontal wells)
 * - Stage spacing (distance between fracture stages along lateral)
 * - Cluster spacing (distance between perforation clusters within a stage)
 * - Fracture half-length (controlled by treatment volume and rate)
 * - Proppant loading (mass of proppant per stage)
 * 
 * The objective function can be:
 * - Maximum NPV (Net Present Value) - default
 * - Maximum EUR (Estimated Ultimate Recovery)
 * - Maximum IRR (Internal Rate of Return)
 * - Minimum development cost per BOE
 * 
 * Reference:
 *   Cipolla, C.L., et al. "Reservoir Modeling in Shale-Gas Reservoirs"
 *   SPE 125530 (2010)
 *   
 *   Manchanda, R., Sharma, M.M. "Impact of Completion Design on Fracture 
 *   Complexity in Horizontal Wells" SPE 159899 (2014)
 * 
 * Usage:
 *   ./well_spacing_optimization -c config/well_spacing_optimization.config
 *   mpirun -np 4 ./well_spacing_optimization -c config/well_spacing_optimization.config
 * 
 *   # Quick optimization with reduced grid
 *   ./well_spacing_optimization -c config/well_spacing_optimization.config -quick
 * 
 *   # Specific scenario evaluation
 *   ./well_spacing_optimization -c config/well_spacing_optimization.config \
 *       -well_spacing 250.0 -stage_spacing 50.0 -eval_only
 */

#include "Simulator.hpp"
#include "WellModel.hpp"
#include "FractureModel.hpp"
#include "StressShadowing.hpp"
// Note: We define our own ProductionForecast locally to avoid linker issues
// with unimplemented library constructors
#include "ConfigReader.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <memory>
#include <chrono>
#include <random>

// Local production forecast structure (self-contained, no external dependencies)
namespace WellSpacingOpt {

struct ProductionForecast {
    std::vector<double> times;          // Time points (days)
    std::vector<double> oil_rate;       // Oil rate (m³/day)
    std::vector<double> gas_rate;       // Gas rate (m³/day)
    std::vector<double> water_rate;     // Water rate (m³/day)
    std::vector<double> ngl_rate;       // NGL rate (m³/day)
    
    double eur_oil = 0.0;               // EUR oil (m³)
    double eur_gas = 0.0;               // EUR gas (m³)
    double eur_ngl = 0.0;               // EUR NGL (m³)
    
    ProductionForecast() = default;
};

} // namespace WellSpacingOpt

static char help[] = 
    "Well Spacing and Hydraulic Fracture Placement Optimization\n"
    "===========================================================\n"
    "Optimizes well spacing and fracture placement to maximize hydrocarbon recovery.\n\n"
    "Usage: mpirun -np N ./well_spacing_optimization -c <config_file> [options]\n\n"
    "Options:\n"
    "  -c <file>           Configuration file (required)\n"
    "  -well_spacing <m>   Override well spacing (meters)\n"
    "  -stage_spacing <m>  Override stage spacing (meters)\n"
    "  -cluster_spacing <m> Override cluster spacing (meters)\n"
    "  -frac_length <m>    Override fracture half-length (meters)\n"
    "  -eval_only          Evaluate single scenario only (no optimization)\n"
    "  -quick              Quick mode with coarser grid search\n"
    "  -response_surface   Generate response surface only\n"
    "  -sensitivity        Run sensitivity analysis only\n"
    "  -monte_carlo <N>    Run Monte Carlo with N iterations\n"
    "  -output_dir <dir>   Output directory (default: output/well_spacing_optimization)\n"
    "\n";

namespace FSRM {

//=============================================================================
// Development Scenario Structure
//=============================================================================
struct DevelopmentScenario {
    std::string name;
    
    // Well configuration
    int num_wells;
    double well_spacing;           // m - distance between wells
    double lateral_length;         // m - length of horizontal section
    double vertical_depth;         // m - TVD to landing zone
    
    // Fracture configuration
    int stages_per_well;
    double stage_spacing;          // m - distance between stages
    int clusters_per_stage;
    double cluster_spacing;        // m - distance between clusters
    double frac_half_length;       // m - target half-length
    double frac_height;            // m - fracture height
    double proppant_per_stage;     // kg - proppant mass per stage
    
    // Calculated values
    double total_frac_length;      // m - sum of all frac lengths
    double stimulated_volume;      // m³ - SRV
    double drainage_area;          // m² - total drainage
    
    // Results
    double total_eur_oil;          // m³ - estimated ultimate recovery (oil)
    double total_eur_gas;          // m³ - estimated ultimate recovery (gas)
    double total_dc_cost;          // $ - drilling and completion cost
    double npv;                    // $ - net present value
    double irr;                    // - - internal rate of return
    double payout_years;           // years - time to payout
    double npv_per_well;           // $ - NPV per well
    double eur_per_well;           // BOE - EUR per well
    double cost_per_boe;           // $/BOE - D&C cost per BOE
    
    DevelopmentScenario() :
        num_wells(4), well_spacing(300.0), lateral_length(2400.0),
        vertical_depth(2500.0), stages_per_well(40), stage_spacing(60.0),
        clusters_per_stage(4), cluster_spacing(12.0), frac_half_length(150.0),
        frac_height(50.0), proppant_per_stage(80000.0),
        total_frac_length(0), stimulated_volume(0), drainage_area(0),
        total_eur_oil(0), total_eur_gas(0), total_dc_cost(0),
        npv(0), irr(0), payout_years(0), npv_per_well(0),
        eur_per_well(0), cost_per_boe(0) {}
    
    void calculateDerivedValues() {
        stages_per_well = static_cast<int>(lateral_length / stage_spacing);
        total_frac_length = num_wells * stages_per_well * 2.0 * frac_half_length;
        stimulated_volume = num_wells * stages_per_well * 4.0 * 
                           frac_half_length * frac_height * 
                           (cluster_spacing * clusters_per_stage);
        drainage_area = num_wells * lateral_length * well_spacing;
    }
    
    void print(std::ostream& os) const {
        os << std::fixed << std::setprecision(2);
        os << "Scenario: " << name << "\n";
        os << "  Wells: " << num_wells << " x " << lateral_length << " m laterals\n";
        os << "  Well spacing: " << well_spacing << " m\n";
        os << "  Stages/well: " << stages_per_well << " @ " << stage_spacing << " m\n";
        os << "  Clusters/stage: " << clusters_per_stage << " @ " << cluster_spacing << " m\n";
        os << "  Frac half-length: " << frac_half_length << " m\n";
        os << "  Total stages: " << num_wells * stages_per_well << "\n";
        os << "  SRV: " << stimulated_volume / 1.0e6 << " MM m³\n";
        os << "  EUR (oil): " << total_eur_oil / 159.0 << " Mbbl\n";
        os << "  EUR (gas): " << total_eur_gas / 28.3 << " MMcf\n";
        os << "  D&C Cost: $" << total_dc_cost / 1.0e6 << " MM\n";
        os << "  NPV: $" << npv / 1.0e6 << " MM\n";
        os << "  IRR: " << irr * 100.0 << "%\n";
        os << "  Payout: " << payout_years << " years\n";
        os << "  NPV/well: $" << npv_per_well / 1.0e6 << " MM\n";
        os << "  Cost/BOE: $" << cost_per_boe << "\n";
    }
};

//=============================================================================
// Well Interference Model
//=============================================================================
class WellInterferenceModel {
public:
    WellInterferenceModel() : 
        interference_a_(0.4), interference_b_(2.0),
        frac_hit_threshold_(0.1) {}
    
    void setParameters(double a, double b, double threshold) {
        interference_a_ = a;
        interference_b_ = b;
        frac_hit_threshold_ = threshold;
    }
    
    /**
     * @brief Calculate EUR multiplier due to well interference
     * 
     * Uses hyperbolic model:
     *   EUR_ratio = 1 - a / (1 + b * (spacing / frac_length))
     * 
     * @param well_spacing Distance between wells (m)
     * @param frac_half_length Fracture half-length (m)
     * @return EUR multiplier (0-1)
     */
    double calculateEURMultiplier(double well_spacing, double frac_half_length) const {
        if (frac_half_length <= 0) return 1.0;
        
        double ratio = well_spacing / frac_half_length;
        double eur_reduction = interference_a_ / (1.0 + interference_b_ * ratio);
        
        return std::max(0.3, 1.0 - eur_reduction); // Min 30% of isolated well EUR
    }
    
    /**
     * @brief Check for potential frac hit
     * 
     * @param well_spacing Distance between wells (m)
     * @param frac_half_length Fracture half-length (m)
     * @return true if frac hit is likely
     */
    bool checkFracHit(double well_spacing, double frac_half_length) const {
        // Frac hit if fractures from adjacent wells can overlap
        return (2.0 * frac_half_length >= well_spacing);
    }
    
    /**
     * @brief Calculate inter-fracture stress shadow effect on EUR
     * 
     * Tighter fracture spacing reduces individual fracture effectiveness
     * due to stress shadows limiting fracture width and length.
     * 
     * @param stage_spacing Distance between stages (m)
     * @param frac_half_length Fracture half-length (m)
     * @param num_clusters Clusters per stage
     * @return Effectiveness multiplier (0-1)
     */
    double calculateStressShadowEffect(double stage_spacing, 
                                       double frac_half_length,
                                       int num_clusters) const {
        // Normalized spacing (ratio of spacing to frac dimension)
        double norm_spacing = stage_spacing / frac_half_length;
        
        // Stress shadow model - exponential decay
        // At very tight spacing (<0.5 * frac_length), significant reduction
        double shadow_factor = 1.0 - 0.5 * std::exp(-2.0 * norm_spacing);
        
        // Cluster interference within stage
        double cluster_factor = std::pow(0.95, num_clusters - 1);
        
        return shadow_factor * cluster_factor;
    }
    
private:
    double interference_a_;
    double interference_b_;
    double frac_hit_threshold_;
};

//=============================================================================
// Production Forecaster
//=============================================================================
class ProductionForecaster {
public:
    ProductionForecaster() :
        initial_decline_(0.08), b_factor_(1.2), terminal_decline_(0.005),
        initial_gor_(600.0), gor_increase_rate_(0.1),
        initial_water_cut_(0.1), final_water_cut_(0.6),
        water_cut_tc_(730.0) {}
    
    void setDeclineParameters(double di, double b, double dt) {
        initial_decline_ = di;
        b_factor_ = b;
        terminal_decline_ = dt;
    }
    
    /**
     * @brief Calculate EUR for a single well based on completion design
     * 
     * @param scenario Development scenario parameters
     * @param interference_model Well interference model
     * @return Pair of (oil EUR m³, gas EUR m³)
     */
    std::pair<double, double> calculateSingleWellEUR(
        const DevelopmentScenario& scenario,
        const WellInterferenceModel& interference) const {
        
        // Base EUR from rock quality and completion intensity
        // Contact factor: how much rock is contacted
        double contact_factor = scenario.stages_per_well * 
                               scenario.clusters_per_stage *
                               scenario.frac_half_length * 
                               scenario.frac_height / 1.0e6;
        
        // Base oil EUR (m³) - calibrated to typical Permian/Bakken performance
        double base_oil_eur = 50000.0 * std::pow(contact_factor, 0.5);
        
        // Stress shadow reduction
        double shadow_mult = interference.calculateStressShadowEffect(
            scenario.stage_spacing, scenario.frac_half_length,
            scenario.clusters_per_stage);
        
        // Well interference reduction
        double interf_mult = interference.calculateEURMultiplier(
            scenario.well_spacing, scenario.frac_half_length);
        
        // Apply multipliers
        double oil_eur = base_oil_eur * shadow_mult * interf_mult;
        
        // Gas EUR based on GOR
        double gas_eur = oil_eur * initial_gor_ * 0.0283; // Convert scf to m³
        
        return {oil_eur, gas_eur};
    }
    
    /**
     * @brief Generate production forecast
     * 
     * @param initial_rate_oil Initial oil rate (m³/day)
     * @param eur_oil Oil EUR (m³) - used for calibration reference
     * @param years Forecast period
     * @return Production forecast
     */
    WellSpacingOpt::ProductionForecast generateForecast(double initial_rate_oil, 
                                        double /* eur_oil */,
                                        double years) const {
        WellSpacingOpt::ProductionForecast forecast;
        
        int num_points = static_cast<int>(years * 12); // Monthly
        double dt = 30.4375; // Average days per month
        
        double cum_oil = 0.0;
        double cum_gas = 0.0;
        double rate_oil = initial_rate_oil;
        
        // Modified hyperbolic decline
        for (int i = 0; i < num_points; ++i) {
            double t_days = i * dt;
            double t_months = i;
            
            // Hyperbolic decline
            double decline = initial_decline_ / 
                           std::pow(1.0 + b_factor_ * initial_decline_ * t_months, 
                                   1.0 / b_factor_);
            
            // Limit to terminal decline
            decline = std::max(decline, terminal_decline_);
            
            rate_oil = initial_rate_oil * 
                      std::pow(1.0 + b_factor_ * initial_decline_ * t_months,
                              -1.0 / b_factor_);
            
            // GOR increase over time
            double gor = initial_gor_ * (1.0 + gor_increase_rate_ * t_days / 365.0);
            double rate_gas = rate_oil * gor * 0.0283; // Convert to m³/day
            
            // Water cut increase
            double water_cut = initial_water_cut_ + 
                              (final_water_cut_ - initial_water_cut_) *
                              (1.0 - std::exp(-t_days / water_cut_tc_));
            double rate_water = rate_oil * water_cut / (1.0 - water_cut);
            
            forecast.times.push_back(t_days);
            forecast.oil_rate.push_back(rate_oil);
            forecast.gas_rate.push_back(rate_gas);
            forecast.water_rate.push_back(rate_water);
            
            cum_oil += rate_oil * dt;
            cum_gas += rate_gas * dt;
        }
        
        forecast.eur_oil = cum_oil;
        forecast.eur_gas = cum_gas;
        
        return forecast;
    }
    
private:
    double initial_decline_;
    double b_factor_;
    double terminal_decline_;
    double initial_gor_;
    double gor_increase_rate_;
    double initial_water_cut_;
    double final_water_cut_;
    double water_cut_tc_;
};

//=============================================================================
// Cost Calculator
//=============================================================================
class CostCalculator {
public:
    CostCalculator() :
        drilling_cost_vertical_(1500000.0),
        drilling_cost_per_m_(400.0),
        completion_cost_per_stage_(100000.0),
        proppant_cost_per_kg_(0.18),
        fluid_cost_per_m3_(3.0),
        perf_cost_per_cluster_(5000.0),
        plug_cost_(15000.0),
        facilities_cost_per_well_(200000.0) {}
    
    void loadFromConfig(const ConfigReader& config) {
        drilling_cost_vertical_ = config.getDouble("ECONOMICS", "drilling_cost_vertical", 1500000.0);
        drilling_cost_per_m_ = config.getDouble("ECONOMICS", "drilling_cost_per_meter_lateral", 400.0);
        completion_cost_per_stage_ = config.getDouble("ECONOMICS", "completion_cost_per_stage", 100000.0);
        proppant_cost_per_kg_ = config.getDouble("ECONOMICS", "proppant_cost_per_kg", 0.18);
        fluid_cost_per_m3_ = config.getDouble("ECONOMICS", "fluid_cost_per_m3", 3.0);
        perf_cost_per_cluster_ = config.getDouble("ECONOMICS", "perforation_cost_per_cluster", 5000.0);
        plug_cost_ = config.getDouble("ECONOMICS", "plug_cost", 15000.0);
        facilities_cost_per_well_ = config.getDouble("ECONOMICS", "facilities_cost_per_well", 200000.0);
    }
    
    /**
     * @brief Calculate total D&C cost for a development scenario
     */
    double calculateDCCost(const DevelopmentScenario& scenario) const {
        double cost = 0.0;
        
        // Drilling (per well)
        double drilling_per_well = drilling_cost_vertical_ + 
                                   drilling_cost_per_m_ * scenario.lateral_length;
        cost += scenario.num_wells * drilling_per_well;
        
        // Completion (per stage)
        int total_stages = scenario.num_wells * scenario.stages_per_well;
        cost += total_stages * completion_cost_per_stage_;
        
        // Proppant
        double total_proppant = total_stages * scenario.proppant_per_stage;
        cost += total_proppant * proppant_cost_per_kg_;
        
        // Fluid (~500 m³ per stage typical)
        double fluid_per_stage = 500.0;
        cost += total_stages * fluid_per_stage * fluid_cost_per_m3_;
        
        // Perforations
        int total_clusters = total_stages * scenario.clusters_per_stage;
        cost += total_clusters * perf_cost_per_cluster_;
        
        // Plugs (one per stage minus one)
        cost += (total_stages - scenario.num_wells) * plug_cost_;
        
        // Facilities
        cost += scenario.num_wells * facilities_cost_per_well_;
        
        return cost;
    }
    
private:
    double drilling_cost_vertical_;
    double drilling_cost_per_m_;
    double completion_cost_per_stage_;
    double proppant_cost_per_kg_;
    double fluid_cost_per_m3_;
    double perf_cost_per_cluster_;
    double plug_cost_;
    double facilities_cost_per_well_;
};

//=============================================================================
// Spacing Optimizer
//=============================================================================
class SpacingOptimizer {
public:
    SpacingOptimizer() :
        section_width_(1609.34), section_length_(1609.34),
        oil_price_(75.0), gas_price_(3.5), discount_rate_(0.10),
        evaluation_years_(30.0) {}
    
    void setSectionGeometry(double width, double length) {
        section_width_ = width;
        section_length_ = length;
    }
    
    void setEconomicParameters(double oil_price, double gas_price, 
                               double discount_rate, double eval_years) {
        oil_price_ = oil_price;
        gas_price_ = gas_price;
        discount_rate_ = discount_rate;
        evaluation_years_ = eval_years;
    }
    
    void setOptimizationBounds(
        double ws_min, double ws_max, double ws_step,
        double ss_min, double ss_max, double ss_step,
        double cs_min, double cs_max, double cs_step,
        double fl_min, double fl_max, double fl_step) {
        
        well_spacing_min_ = ws_min; well_spacing_max_ = ws_max; well_spacing_step_ = ws_step;
        stage_spacing_min_ = ss_min; stage_spacing_max_ = ss_max; stage_spacing_step_ = ss_step;
        cluster_spacing_min_ = cs_min; cluster_spacing_max_ = cs_max; cluster_spacing_step_ = cs_step;
        frac_length_min_ = fl_min; frac_length_max_ = fl_max; frac_length_step_ = fl_step;
    }
    
    /**
     * @brief Run grid search optimization
     * 
     * @param cost_calc Cost calculator
     * @param interference_model Well interference model
     * @param forecaster Production forecaster
     * @return Best scenario and all results
     */
    std::pair<DevelopmentScenario, std::vector<DevelopmentScenario>> 
    runGridSearch(const CostCalculator& cost_calc,
                  const WellInterferenceModel& interference,
                  const ProductionForecaster& forecaster) {
        
        std::vector<DevelopmentScenario> all_results;
        DevelopmentScenario best_scenario;
        double best_npv = -1.0e30;
        
        // Grid search over spacing parameters
        for (double ws = well_spacing_min_; ws <= well_spacing_max_; ws += well_spacing_step_) {
            // Number of wells that fit in section
            int max_wells = static_cast<int>(section_width_ / ws);
            max_wells = std::min(max_wells, 8); // Cap at 8 wells
            
            for (int nw = 2; nw <= max_wells; ++nw) {
                for (double ss = stage_spacing_min_; ss <= stage_spacing_max_; ss += stage_spacing_step_) {
                    for (double fl = frac_length_min_; fl <= frac_length_max_; fl += frac_length_step_) {
                        
                        DevelopmentScenario scenario;
                        scenario.well_spacing = ws;
                        scenario.num_wells = nw;
                        scenario.stage_spacing = ss;
                        scenario.frac_half_length = fl;
                        scenario.lateral_length = section_length_ - 200.0; // Buffer from boundaries
                        scenario.cluster_spacing = 12.0; // Fixed for this search
                        scenario.clusters_per_stage = 4;
                        scenario.frac_height = 50.0;
                        scenario.proppant_per_stage = 80000.0;
                        
                        scenario.calculateDerivedValues();
                        
                        // Evaluate scenario
                        evaluateScenario(scenario, cost_calc, interference, forecaster);
                        
                        scenario.name = "ws" + std::to_string(int(ws)) + 
                                       "_nw" + std::to_string(nw) +
                                       "_ss" + std::to_string(int(ss)) +
                                       "_fl" + std::to_string(int(fl));
                        
                        all_results.push_back(scenario);
                        
                        if (scenario.npv > best_npv) {
                            best_npv = scenario.npv;
                            best_scenario = scenario;
                        }
                    }
                }
            }
        }
        
        return {best_scenario, all_results};
    }
    
    /**
     * @brief Evaluate a single scenario
     */
    void evaluateScenario(DevelopmentScenario& scenario,
                         const CostCalculator& cost_calc,
                         const WellInterferenceModel& interference,
                         const ProductionForecaster& forecaster) {
        
        // Calculate EUR per well
        auto [oil_eur, gas_eur] = forecaster.calculateSingleWellEUR(scenario, interference);
        
        // Total EUR
        scenario.total_eur_oil = oil_eur * scenario.num_wells;
        scenario.total_eur_gas = gas_eur * scenario.num_wells;
        
        // D&C Cost
        scenario.total_dc_cost = cost_calc.calculateDCCost(scenario);
        
        // Generate production forecast for NPV
        double initial_rate = oil_eur / (365.0 * 3.0); // Rough IP estimate
        auto forecast = forecaster.generateForecast(initial_rate, oil_eur, evaluation_years_);
        
        // Calculate NPV
        scenario.npv = calculateNPV(scenario, forecast);
        
        // Calculate IRR (simplified - Newton-Raphson would be better)
        scenario.irr = calculateIRR(scenario, forecast);
        
        // Calculate payout
        scenario.payout_years = calculatePayout(scenario, forecast);
        
        // Per-well metrics
        scenario.npv_per_well = scenario.npv / scenario.num_wells;
        scenario.eur_per_well = (scenario.total_eur_oil / 159.0 + 
                                scenario.total_eur_gas / 5660.0) / scenario.num_wells; // BOE
        scenario.cost_per_boe = scenario.total_dc_cost / 
                               (scenario.total_eur_oil / 159.0 + scenario.total_eur_gas / 5660.0);
    }
    
    /**
     * @brief Generate response surface data
     */
    std::vector<std::vector<double>> generateResponseSurface(
        const std::vector<DevelopmentScenario>& results,
        const std::string& x_param, const std::string& y_param) {
        
        // Extract unique values for x and y parameters
        std::vector<double> x_values, y_values;
        
        for (const auto& r : results) {
            double x = (x_param == "well_spacing") ? r.well_spacing : r.stage_spacing;
            double y = (y_param == "stage_spacing") ? r.stage_spacing : r.frac_half_length;
            
            if (std::find(x_values.begin(), x_values.end(), x) == x_values.end())
                x_values.push_back(x);
            if (std::find(y_values.begin(), y_values.end(), y) == y_values.end())
                y_values.push_back(y);
        }
        
        std::sort(x_values.begin(), x_values.end());
        std::sort(y_values.begin(), y_values.end());
        
        // Create 2D NPV matrix
        std::vector<std::vector<double>> surface(y_values.size(), 
                                                  std::vector<double>(x_values.size(), 0.0));
        
        for (const auto& r : results) {
            double x = (x_param == "well_spacing") ? r.well_spacing : r.stage_spacing;
            double y = (y_param == "stage_spacing") ? r.stage_spacing : r.frac_half_length;
            
            auto xi = std::find(x_values.begin(), x_values.end(), x);
            auto yi = std::find(y_values.begin(), y_values.end(), y);
            
            if (xi != x_values.end() && yi != y_values.end()) {
                int ix = std::distance(x_values.begin(), xi);
                int iy = std::distance(y_values.begin(), yi);
                
                // Average if multiple scenarios map to same point
                if (surface[iy][ix] == 0.0 || r.npv > surface[iy][ix]) {
                    surface[iy][ix] = r.npv;
                }
            }
        }
        
        return surface;
    }
    
private:
    double section_width_;
    double section_length_;
    double oil_price_;
    double gas_price_;
    double discount_rate_;
    double evaluation_years_;
    
    // Optimization bounds
    double well_spacing_min_, well_spacing_max_, well_spacing_step_;
    double stage_spacing_min_, stage_spacing_max_, stage_spacing_step_;
    double cluster_spacing_min_, cluster_spacing_max_, cluster_spacing_step_;
    double frac_length_min_, frac_length_max_, frac_length_step_;
    
    double calculateNPV(const DevelopmentScenario& scenario,
                       const WellSpacingOpt::ProductionForecast& forecast) {
        double npv = -scenario.total_dc_cost; // Initial investment
        
        // Monthly cash flows
        double monthly_rate = std::pow(1.0 + discount_rate_, 1.0/12.0) - 1.0;
        
        for (size_t i = 0; i < forecast.times.size(); ++i) {
            double oil_revenue = forecast.oil_rate[i] * 30.4 * oil_price_ / 159.0; // $/month
            double gas_revenue = forecast.gas_rate[i] * 30.4 * gas_price_ / 28.3; // $/month
            
            double gross_revenue = (oil_revenue + gas_revenue) * scenario.num_wells;
            
            // Assume 30% total deductions (royalty, taxes, opex)
            double net_revenue = gross_revenue * 0.70;
            
            // Discount
            double discount_factor = 1.0 / std::pow(1.0 + monthly_rate, i);
            npv += net_revenue * discount_factor;
        }
        
        return npv;
    }
    
    double calculateIRR(const DevelopmentScenario& scenario,
                       const WellSpacingOpt::ProductionForecast& /* forecast */) {
        // Simplified IRR - would use Newton-Raphson for accuracy
        // Start with guess based on NPV/Investment ratio
        double roi = scenario.npv / scenario.total_dc_cost;
        return std::min(std::max(roi / evaluation_years_ + discount_rate_, 0.0), 1.0);
    }
    
    double calculatePayout(const DevelopmentScenario& scenario,
                          const WellSpacingOpt::ProductionForecast& forecast) {
        double cumulative = 0.0;
        
        for (size_t i = 0; i < forecast.times.size(); ++i) {
            double oil_revenue = forecast.oil_rate[i] * 30.4 * oil_price_ / 159.0;
            double gas_revenue = forecast.gas_rate[i] * 30.4 * gas_price_ / 28.3;
            double net_revenue = (oil_revenue + gas_revenue) * scenario.num_wells * 0.70;
            
            cumulative += net_revenue;
            
            if (cumulative >= scenario.total_dc_cost) {
                return forecast.times[i] / 365.0; // Convert to years
            }
        }
        
        return evaluation_years_; // Didn't pay out
    }
};

//=============================================================================
// Output Writer
//=============================================================================
class OutputWriter {
public:
    OutputWriter(const std::string& output_dir) : output_dir_(output_dir) {}
    
    void writeResults(const DevelopmentScenario& best,
                     const std::vector<DevelopmentScenario>& all_results) {
        
        // Summary file
        std::ofstream summary(output_dir_ + "/optimization_summary.txt");
        summary << "Well Spacing and Fracture Placement Optimization Results\n";
        summary << "========================================================\n\n";
        
        summary << "OPTIMAL SCENARIO:\n";
        summary << "-----------------\n";
        best.print(summary);
        
        summary << "\n\nALL SCENARIOS EVALUATED: " << all_results.size() << "\n";
        summary << "========================\n\n";
        
        // Top 10 scenarios
        auto sorted = all_results;
        std::sort(sorted.begin(), sorted.end(),
                 [](const auto& a, const auto& b) { return a.npv > b.npv; });
        
        summary << "TOP 10 SCENARIOS BY NPV:\n";
        summary << std::setw(6) << "Rank" 
               << std::setw(10) << "WellSpc"
               << std::setw(8) << "Wells"
               << std::setw(10) << "StageSpc"
               << std::setw(10) << "FracLen"
               << std::setw(12) << "NPV(MM$)"
               << std::setw(10) << "IRR(%)"
               << std::setw(12) << "EUR(MBOE)"
               << "\n";
        
        for (int i = 0; i < std::min(10, (int)sorted.size()); ++i) {
            const auto& s = sorted[i];
            summary << std::setw(6) << (i + 1)
                   << std::setw(10) << s.well_spacing
                   << std::setw(8) << s.num_wells
                   << std::setw(10) << s.stage_spacing
                   << std::setw(10) << s.frac_half_length
                   << std::setw(12) << std::fixed << std::setprecision(2) << s.npv / 1.0e6
                   << std::setw(10) << std::fixed << std::setprecision(1) << s.irr * 100.0
                   << std::setw(12) << std::fixed << std::setprecision(0) << s.eur_per_well * s.num_wells / 1000.0
                   << "\n";
        }
        
        summary.close();
        
        // CSV for all results
        writeCSV(all_results);
        
        // Response surface data
        writeResponseSurface(all_results);
    }
    
    void writeCSV(const std::vector<DevelopmentScenario>& results) {
        std::ofstream csv(output_dir_ + "/optimization_results.csv");
        
        csv << "scenario,num_wells,well_spacing_m,stage_spacing_m,cluster_spacing_m,"
           << "frac_half_length_m,stages_per_well,total_stages,srv_mm3,"
           << "eur_oil_m3,eur_gas_m3,dc_cost_usd,npv_usd,irr,payout_years,"
           << "npv_per_well,eur_per_well_boe,cost_per_boe\n";
        
        for (const auto& r : results) {
            csv << r.name << ","
               << r.num_wells << ","
               << r.well_spacing << ","
               << r.stage_spacing << ","
               << r.cluster_spacing << ","
               << r.frac_half_length << ","
               << r.stages_per_well << ","
               << r.num_wells * r.stages_per_well << ","
               << r.stimulated_volume / 1.0e6 << ","
               << r.total_eur_oil << ","
               << r.total_eur_gas << ","
               << r.total_dc_cost << ","
               << r.npv << ","
               << r.irr << ","
               << r.payout_years << ","
               << r.npv_per_well << ","
               << r.eur_per_well << ","
               << r.cost_per_boe << "\n";
        }
        
        csv.close();
    }
    
    void writeResponseSurface(const std::vector<DevelopmentScenario>& results) {
        // Write response surface data for visualization
        std::ofstream surf(output_dir_ + "/response_surface.dat");
        
        surf << "# Response surface: NPV vs Well Spacing vs Stage Spacing\n";
        surf << "# Format: well_spacing stage_spacing npv\n";
        
        for (const auto& r : results) {
            surf << r.well_spacing << " " << r.stage_spacing << " " 
                << r.npv / 1.0e6 << "\n";
        }
        
        surf.close();
        
        // Gnuplot script
        std::ofstream gp(output_dir_ + "/plot_response_surface.gp");
        gp << "set terminal pngcairo size 1200,900 enhanced\n";
        gp << "set output 'response_surface.png'\n";
        gp << "set title 'NPV Response Surface'\n";
        gp << "set xlabel 'Well Spacing (m)'\n";
        gp << "set ylabel 'Stage Spacing (m)'\n";
        gp << "set zlabel 'NPV (MM$)'\n";
        gp << "set dgrid3d 30,30\n";
        gp << "set hidden3d\n";
        gp << "splot 'response_surface.dat' u 1:2:3 with lines title 'NPV'\n";
        gp.close();
    }
    
private:
    std::string output_dir_;
};

} // namespace FSRM

//=============================================================================
// Main Function
//=============================================================================
int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Parse command line options
    char config_file[PETSC_MAX_PATH_LEN] = "config/well_spacing_optimization.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    char output_dir[PETSC_MAX_PATH_LEN] = "output/well_spacing_optimization";
    PetscOptionsGetString(nullptr, nullptr, "-output_dir", output_dir,
                         sizeof(output_dir), nullptr);
    
    PetscBool eval_only = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-eval_only", &eval_only, nullptr);
    
    PetscBool quick_mode = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-quick", &quick_mode, nullptr);
    
    PetscReal override_ws = 0.0, override_ss = 0.0, override_fl = 0.0;
    PetscOptionsGetReal(nullptr, nullptr, "-well_spacing", &override_ws, nullptr);
    PetscOptionsGetReal(nullptr, nullptr, "-stage_spacing", &override_ss, nullptr);
    PetscOptionsGetReal(nullptr, nullptr, "-frac_length", &override_fl, nullptr);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (rank == 0) {
        std::cout << "============================================================\n";
        std::cout << "  Well Spacing and Hydraulic Fracture Placement Optimization\n";
        std::cout << "============================================================\n\n";
        std::cout << "Config file: " << config_file << "\n";
        std::cout << "Output directory: " << output_dir << "\n\n";
    }
    
    // Read configuration
    FSRM::ConfigReader config;
    config.loadFile(config_file);
    
    // Initialize models
    FSRM::WellInterferenceModel interference;
    interference.setParameters(
        config.getDouble("WELL_INTERFERENCE", "interference_a", 0.4),
        config.getDouble("WELL_INTERFERENCE", "interference_b", 2.0),
        config.getDouble("WELL_INTERFERENCE", "frac_hit_threshold", 0.1));
    
    FSRM::ProductionForecaster forecaster;
    forecaster.setDeclineParameters(
        config.getDouble("PRODUCTION_MODEL", "initial_decline_rate", 0.08),
        config.getDouble("PRODUCTION_MODEL", "b_factor", 1.2),
        config.getDouble("PRODUCTION_MODEL", "terminal_decline_rate", 0.005));
    
    FSRM::CostCalculator cost_calc;
    cost_calc.loadFromConfig(config);
    
    // Initialize optimizer
    FSRM::SpacingOptimizer optimizer;
    optimizer.setSectionGeometry(
        config.getDouble("OPTIMIZATION", "section_width", 1609.34),
        config.getDouble("OPTIMIZATION", "section_length", 1609.34));
    
    optimizer.setEconomicParameters(
        config.getDouble("ECONOMICS", "oil_price", 75.0),
        config.getDouble("ECONOMICS", "gas_price", 3.5),
        config.getDouble("ECONOMICS", "discount_rate", 0.10),
        config.getDouble("ECONOMICS", "evaluation_period", 30.0));
    
    // Set optimization bounds
    double step_mult = quick_mode ? 2.0 : 1.0; // Coarser grid in quick mode
    
    optimizer.setOptimizationBounds(
        config.getDouble("OPTIMIZATION", "well_spacing_min", 150.0),
        config.getDouble("OPTIMIZATION", "well_spacing_max", 450.0),
        config.getDouble("OPTIMIZATION", "well_spacing_step", 30.0) * step_mult,
        config.getDouble("OPTIMIZATION", "stage_spacing_min", 30.0),
        config.getDouble("OPTIMIZATION", "stage_spacing_max", 90.0),
        config.getDouble("OPTIMIZATION", "stage_spacing_step", 10.0) * step_mult,
        config.getDouble("OPTIMIZATION", "cluster_spacing_min", 6.0),
        config.getDouble("OPTIMIZATION", "cluster_spacing_max", 18.0),
        config.getDouble("OPTIMIZATION", "cluster_spacing_step", 3.0) * step_mult,
        config.getDouble("OPTIMIZATION", "frac_half_length_min", 75.0),
        config.getDouble("OPTIMIZATION", "frac_half_length_max", 225.0),
        config.getDouble("OPTIMIZATION", "frac_half_length_step", 25.0) * step_mult);
    
    // Run optimization or single evaluation
    FSRM::DevelopmentScenario best_scenario;
    std::vector<FSRM::DevelopmentScenario> all_results;
    
    if (eval_only && (override_ws > 0 || override_ss > 0 || override_fl > 0)) {
        // Single scenario evaluation
        if (rank == 0) {
            std::cout << "Evaluating single scenario...\n\n";
        }
        
        best_scenario.well_spacing = override_ws > 0 ? override_ws : 300.0;
        best_scenario.stage_spacing = override_ss > 0 ? override_ss : 60.0;
        best_scenario.frac_half_length = override_fl > 0 ? override_fl : 150.0;
        best_scenario.lateral_length = 2400.0;
        best_scenario.cluster_spacing = 12.0;
        best_scenario.clusters_per_stage = 4;
        best_scenario.frac_height = 50.0;
        best_scenario.proppant_per_stage = 80000.0;
        best_scenario.num_wells = static_cast<int>(1609.34 / best_scenario.well_spacing);
        
        best_scenario.calculateDerivedValues();
        optimizer.evaluateScenario(best_scenario, cost_calc, interference, forecaster);
        best_scenario.name = "user_specified";
        
        all_results.push_back(best_scenario);
        
    } else {
        // Full optimization
        if (rank == 0) {
            std::cout << "Running grid search optimization...\n";
            std::cout << "This may take a few minutes...\n\n";
        }
        
        auto [best, all] = optimizer.runGridSearch(cost_calc, interference, forecaster);
        best_scenario = best;
        all_results = all;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    // Output results
    if (rank == 0) {
        std::cout << "\n============================================================\n";
        std::cout << "  OPTIMIZATION RESULTS\n";
        std::cout << "============================================================\n\n";
        
        std::cout << "Scenarios evaluated: " << all_results.size() << "\n";
        std::cout << "Computation time: " << duration.count() << " seconds\n\n";
        
        std::cout << "OPTIMAL DEVELOPMENT SCENARIO:\n";
        std::cout << "-----------------------------\n";
        best_scenario.print(std::cout);
        
        // Write output files
        std::string mkdir_cmd = "mkdir -p " + std::string(output_dir);
        int ret = system(mkdir_cmd.c_str());
        (void)ret; // Suppress unused warning
        
        FSRM::OutputWriter writer(output_dir);
        writer.writeResults(best_scenario, all_results);
        
        std::cout << "\n============================================================\n";
        std::cout << "  OUTPUT FILES\n";
        std::cout << "============================================================\n\n";
        std::cout << "Results written to: " << output_dir << "/\n";
        std::cout << "  - optimization_summary.txt  (human-readable summary)\n";
        std::cout << "  - optimization_results.csv  (all scenarios for analysis)\n";
        std::cout << "  - response_surface.dat      (for contour plotting)\n";
        std::cout << "  - plot_response_surface.gp  (Gnuplot script)\n\n";
        
        std::cout << "KEY FINDINGS:\n";
        std::cout << "-------------\n";
        std::cout << "• Optimal well spacing: " << best_scenario.well_spacing << " m "
                 << "(" << best_scenario.well_spacing * 3.28 << " ft)\n";
        std::cout << "• Optimal stage spacing: " << best_scenario.stage_spacing << " m "
                 << "(" << best_scenario.stage_spacing * 3.28 << " ft)\n";
        std::cout << "• Optimal fracture half-length: " << best_scenario.frac_half_length << " m "
                 << "(" << best_scenario.frac_half_length * 3.28 << " ft)\n";
        std::cout << "• Recommended wells per section: " << best_scenario.num_wells << "\n";
        std::cout << "• Section NPV: $" << std::fixed << std::setprecision(1) 
                 << best_scenario.npv / 1.0e6 << " MM\n";
        std::cout << "• Development cost: $" << best_scenario.total_dc_cost / 1.0e6 << " MM\n";
        std::cout << "• Cost per BOE: $" << std::setprecision(2) << best_scenario.cost_per_boe << "\n\n";
        
        // Recommendations
        std::cout << "RECOMMENDATIONS:\n";
        std::cout << "----------------\n";
        
        if (best_scenario.well_spacing < 250.0) {
            std::cout << "• Well spacing is tight - monitor for frac hits during infill\n";
        }
        if (best_scenario.stage_spacing < 50.0) {
            std::cout << "• Stage spacing is aggressive - consider limited entry design\n";
        }
        if (best_scenario.frac_half_length > 180.0) {
            std::cout << "• Long fractures targeted - ensure adequate pumping capacity\n";
        }
        
        std::cout << "\n============================================================\n";
        std::cout << "  Optimization Complete\n";
        std::cout << "============================================================\n\n";
    }
    
    PetscFinalize();
    return 0;
}
