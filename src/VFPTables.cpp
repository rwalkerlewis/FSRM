#include "VFPTables.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <numeric>

namespace FSRM {

// ============================================================================
// VFPProdTable Implementation
// ============================================================================

VFPProdTable::VFPProdTable(int table_num)
    : table_number_(table_num),
      datum_depth_(0.0),
      flow_type_(VFPFlowType::OIL),
      wct_type_(VFPWaterCutType::WCT),
      glr_type_(VFPGasLiquidType::GOR),
      alq_type_(VFPArtificialLiftType::NONE),
      unit_system_(VFPUnitSystem::DEFAULT) {}

void VFPProdTable::setTHPValues(const std::vector<double>& thp) {
    thp_values_ = thp;
}

void VFPProdTable::setFlowRateValues(const std::vector<double>& rates) {
    rate_values_ = rates;
}

void VFPProdTable::setWaterCutValues(const std::vector<double>& wct) {
    wct_values_ = wct;
}

void VFPProdTable::setGORValues(const std::vector<double>& gor) {
    gor_values_ = gor;
}

void VFPProdTable::setALQValues(const std::vector<double>& alq) {
    alq_values_ = alq;
}

void VFPProdTable::setBHPData(const std::vector<double>& bhp) {
    bhp_data_ = bhp;
}

void VFPProdTable::addDataPoint(const VFPDataPoint& point) {
    data_points_.push_back(point);
}

void VFPProdTable::buildTableFromPoints() {
    if (data_points_.empty()) return;
    
    // Extract unique axis values from data points
    std::vector<double> thp_set, rate_set, wct_set, gor_set, alq_set;
    
    for (const auto& p : data_points_) {
        if (std::find(thp_set.begin(), thp_set.end(), p.thp) == thp_set.end())
            thp_set.push_back(p.thp);
        if (std::find(rate_set.begin(), rate_set.end(), p.flow_rate) == rate_set.end())
            rate_set.push_back(p.flow_rate);
        if (std::find(wct_set.begin(), wct_set.end(), p.wct) == wct_set.end())
            wct_set.push_back(p.wct);
        if (std::find(gor_set.begin(), gor_set.end(), p.gor) == gor_set.end())
            gor_set.push_back(p.gor);
        if (std::find(alq_set.begin(), alq_set.end(), p.alq) == alq_set.end())
            alq_set.push_back(p.alq);
    }
    
    // Sort axis values
    std::sort(thp_set.begin(), thp_set.end());
    std::sort(rate_set.begin(), rate_set.end());
    std::sort(wct_set.begin(), wct_set.end());
    std::sort(gor_set.begin(), gor_set.end());
    std::sort(alq_set.begin(), alq_set.end());
    
    thp_values_ = thp_set;
    rate_values_ = rate_set;
    wct_values_ = wct_set;
    gor_values_ = gor_set;
    alq_values_ = alq_set;
    
    // Allocate and fill BHP data
    size_t n_total = thp_values_.size() * rate_values_.size() * 
                     wct_values_.size() * gor_values_.size() * alq_values_.size();
    bhp_data_.resize(n_total, 0.0);
    
    // Fill from data points
    for (const auto& p : data_points_) {
        auto it_thp = std::find(thp_values_.begin(), thp_values_.end(), p.thp);
        auto it_rate = std::find(rate_values_.begin(), rate_values_.end(), p.flow_rate);
        auto it_wct = std::find(wct_values_.begin(), wct_values_.end(), p.wct);
        auto it_gor = std::find(gor_values_.begin(), gor_values_.end(), p.gor);
        auto it_alq = std::find(alq_values_.begin(), alq_values_.end(), p.alq);
        
        size_t i_thp = it_thp - thp_values_.begin();
        size_t i_rate = it_rate - rate_values_.begin();
        size_t i_wct = it_wct - wct_values_.begin();
        size_t i_gor = it_gor - gor_values_.begin();
        size_t i_alq = it_alq - alq_values_.begin();
        
        bhp_data_[getIndex(i_thp, i_rate, i_wct, i_gor, i_alq)] = p.bhp;
    }
}

size_t VFPProdTable::getIndex(size_t i_thp, size_t i_rate, size_t i_wct,
                               size_t i_gor, size_t i_alq) const {
    size_t n_rate = rate_values_.size();
    size_t n_wct = wct_values_.size();
    size_t n_gor = gor_values_.size();
    size_t n_alq = alq_values_.size();
    
    return i_thp * (n_rate * n_wct * n_gor * n_alq) +
           i_rate * (n_wct * n_gor * n_alq) +
           i_wct * (n_gor * n_alq) +
           i_gor * n_alq +
           i_alq;
}

std::pair<size_t, double> VFPProdTable::findBracket(double value,
                                                     const std::vector<double>& axis) const {
    if (axis.empty()) return {0, 0.0};
    
    if (value <= axis.front()) return {0, 0.0};
    if (value >= axis.back()) return {axis.size() - 2, 1.0};
    
    for (size_t i = 0; i < axis.size() - 1; ++i) {
        if (value >= axis[i] && value <= axis[i + 1]) {
            double t = (value - axis[i]) / (axis[i + 1] - axis[i]);
            return {i, t};
        }
    }
    
    return {0, 0.0};
}

double VFPProdTable::interpolate5D(double thp, double rate, double wct,
                                    double gor, double alq) const {
    // Find brackets for each dimension
    auto [i_thp, t_thp] = findBracket(thp, thp_values_);
    auto [i_rate, t_rate] = findBracket(rate, rate_values_);
    auto [i_wct, t_wct] = findBracket(wct, wct_values_);
    auto [i_gor, t_gor] = findBracket(gor, gor_values_);
    auto [i_alq, t_alq] = findBracket(alq, alq_values_);
    
    // 5D linear interpolation (32 corner values)
    double result = 0.0;
    
    for (int d_thp = 0; d_thp <= 1; ++d_thp) {
        for (int d_rate = 0; d_rate <= 1; ++d_rate) {
            for (int d_wct = 0; d_wct <= 1; ++d_wct) {
                for (int d_gor = 0; d_gor <= 1; ++d_gor) {
                    for (int d_alq = 0; d_alq <= 1; ++d_alq) {
                        size_t idx = getIndex(
                            std::min(i_thp + d_thp, thp_values_.size() - 1),
                            std::min(i_rate + d_rate, rate_values_.size() - 1),
                            std::min(i_wct + d_wct, wct_values_.size() - 1),
                            std::min(i_gor + d_gor, gor_values_.size() - 1),
                            std::min(i_alq + d_alq, alq_values_.size() - 1)
                        );
                        
                        if (idx >= bhp_data_.size()) continue;
                        
                        double weight = 
                            (d_thp == 0 ? (1.0 - t_thp) : t_thp) *
                            (d_rate == 0 ? (1.0 - t_rate) : t_rate) *
                            (d_wct == 0 ? (1.0 - t_wct) : t_wct) *
                            (d_gor == 0 ? (1.0 - t_gor) : t_gor) *
                            (d_alq == 0 ? (1.0 - t_alq) : t_alq);
                        
                        result += weight * bhp_data_[idx];
                    }
                }
            }
        }
    }
    
    return result;
}

double VFPProdTable::calculateBHP(double thp, double flow_rate, double wct,
                                   double gor, double alq) const {
    if (bhp_data_.empty()) return thp;  // Fallback
    return interpolate5D(thp, flow_rate, wct, gor, alq);
}

double VFPProdTable::calculateFlowRate(double bhp, double thp, double wct,
                                        double gor, double alq) const {
    // Binary search to find flow rate that gives target BHP
    if (rate_values_.empty()) return 0.0;
    
    double rate_min = rate_values_.front();
    double rate_max = rate_values_.back();
    
    // Check if solution exists
    double bhp_min = calculateBHP(thp, rate_min, wct, gor, alq);
    double bhp_max = calculateBHP(thp, rate_max, wct, gor, alq);
    
    if (bhp < bhp_min || bhp > bhp_max) {
        // BHP out of table range
        return (bhp < bhp_min) ? 0.0 : rate_max;
    }
    
    // Bisection search
    double tol = 1e-6;
    int max_iter = 50;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        double rate_mid = 0.5 * (rate_min + rate_max);
        double bhp_mid = calculateBHP(thp, rate_mid, wct, gor, alq);
        
        if (std::abs(bhp_mid - bhp) < tol * std::abs(bhp)) {
            return rate_mid;
        }
        
        if (bhp_mid > bhp) {
            rate_max = rate_mid;
        } else {
            rate_min = rate_mid;
        }
    }
    
    return 0.5 * (rate_min + rate_max);
}

double VFPProdTable::optimizeGasLift(double bhp, double thp, double wct,
                                      double gor, double max_alq) const {
    if (alq_values_.empty() || alq_type_ == VFPArtificialLiftType::NONE) {
        return 0.0;
    }
    
    // Simple optimization: find ALQ that maximizes oil rate at given BHP
    double best_alq = 0.0;
    double best_rate = calculateFlowRate(bhp, thp, wct, gor, 0.0);
    
    for (const auto& alq : alq_values_) {
        if (alq > max_alq) break;
        
        double rate = calculateFlowRate(bhp, thp, wct, gor, alq);
        if (rate > best_rate) {
            best_rate = rate;
            best_alq = alq;
        }
    }
    
    return best_alq;
}

bool VFPProdTable::isValid() const {
    if (thp_values_.empty() || rate_values_.empty()) return false;
    
    size_t expected_size = thp_values_.size() * rate_values_.size() *
                          std::max(size_t(1), wct_values_.size()) *
                          std::max(size_t(1), gor_values_.size()) *
                          std::max(size_t(1), alq_values_.size());
    
    return bhp_data_.size() >= expected_size;
}

void VFPProdTable::printTable() const {
    std::cout << "\n=== VFPPROD Table " << table_number_ << " ===\n";
    std::cout << "Datum depth: " << datum_depth_ << " m\n";
    std::cout << "THP values: " << thp_values_.size() << "\n";
    std::cout << "Rate values: " << rate_values_.size() << "\n";
    std::cout << "WCT values: " << wct_values_.size() << "\n";
    std::cout << "GOR values: " << gor_values_.size() << "\n";
    std::cout << "ALQ values: " << alq_values_.size() << "\n";
    std::cout << "BHP data points: " << bhp_data_.size() << "\n";
}

// ============================================================================
// VFPInjTable Implementation
// ============================================================================

VFPInjTable::VFPInjTable(int table_num)
    : table_number_(table_num),
      datum_depth_(0.0),
      flow_type_(VFPFlowType::WATER),
      unit_system_(VFPUnitSystem::DEFAULT) {}

void VFPInjTable::setTHPValues(const std::vector<double>& thp) {
    thp_values_ = thp;
}

void VFPInjTable::setFlowRateValues(const std::vector<double>& rates) {
    rate_values_ = rates;
}

void VFPInjTable::setBHPData(const std::vector<double>& bhp) {
    bhp_data_ = bhp;
}

double VFPInjTable::interpolate2D(double thp, double rate) const {
    if (bhp_data_.empty() || thp_values_.empty() || rate_values_.empty()) {
        return thp;
    }
    
    // Find brackets
    size_t i_thp = 0, i_rate = 0;
    double t_thp = 0.0, t_rate = 0.0;
    
    // THP bracket
    if (thp <= thp_values_.front()) {
        i_thp = 0; t_thp = 0.0;
    } else if (thp >= thp_values_.back()) {
        i_thp = thp_values_.size() - 2; t_thp = 1.0;
    } else {
        for (size_t i = 0; i < thp_values_.size() - 1; ++i) {
            if (thp >= thp_values_[i] && thp <= thp_values_[i + 1]) {
                i_thp = i;
                t_thp = (thp - thp_values_[i]) / (thp_values_[i + 1] - thp_values_[i]);
                break;
            }
        }
    }
    
    // Rate bracket
    if (rate <= rate_values_.front()) {
        i_rate = 0; t_rate = 0.0;
    } else if (rate >= rate_values_.back()) {
        i_rate = rate_values_.size() - 2; t_rate = 1.0;
    } else {
        for (size_t i = 0; i < rate_values_.size() - 1; ++i) {
            if (rate >= rate_values_[i] && rate <= rate_values_[i + 1]) {
                i_rate = i;
                t_rate = (rate - rate_values_[i]) / (rate_values_[i + 1] - rate_values_[i]);
                break;
            }
        }
    }
    
    // Bilinear interpolation
    size_t n_rate = rate_values_.size();
    
    double bhp_00 = bhp_data_[i_thp * n_rate + i_rate];
    double bhp_01 = bhp_data_[i_thp * n_rate + std::min(i_rate + 1, n_rate - 1)];
    double bhp_10 = bhp_data_[std::min(i_thp + 1, thp_values_.size() - 1) * n_rate + i_rate];
    double bhp_11 = bhp_data_[std::min(i_thp + 1, thp_values_.size() - 1) * n_rate + 
                              std::min(i_rate + 1, n_rate - 1)];
    
    double bhp_0 = bhp_00 * (1.0 - t_rate) + bhp_01 * t_rate;
    double bhp_1 = bhp_10 * (1.0 - t_rate) + bhp_11 * t_rate;
    
    return bhp_0 * (1.0 - t_thp) + bhp_1 * t_thp;
}

double VFPInjTable::calculateBHP(double thp, double flow_rate) const {
    return interpolate2D(thp, flow_rate);
}

double VFPInjTable::calculateFlowRate(double bhp, double thp) const {
    if (rate_values_.empty()) return 0.0;
    
    double rate_min = rate_values_.front();
    double rate_max = rate_values_.back();
    
    // Bisection search
    double tol = 1e-6;
    int max_iter = 50;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        double rate_mid = 0.5 * (rate_min + rate_max);
        double bhp_mid = calculateBHP(thp, rate_mid);
        
        if (std::abs(bhp_mid - bhp) < tol * std::abs(bhp)) {
            return rate_mid;
        }
        
        if (bhp_mid > bhp) {
            rate_max = rate_mid;
        } else {
            rate_min = rate_mid;
        }
    }
    
    return 0.5 * (rate_min + rate_max);
}

bool VFPInjTable::isValid() const {
    return !thp_values_.empty() && !rate_values_.empty() &&
           bhp_data_.size() >= thp_values_.size() * rate_values_.size();
}

void VFPInjTable::printTable() const {
    std::cout << "\n=== VFPINJ Table " << table_number_ << " ===\n";
    std::cout << "Datum depth: " << datum_depth_ << " m\n";
    std::cout << "THP values: " << thp_values_.size() << "\n";
    std::cout << "Rate values: " << rate_values_.size() << "\n";
    std::cout << "BHP data points: " << bhp_data_.size() << "\n";
}

// ============================================================================
// VFPTableGenerator Implementation
// ============================================================================

VFPTableGenerator::VFPTableGenerator()
    : wellbore_radius_(0.1),
      roughness_(0.0001),
      total_depth_(2000.0),
      oil_density_(800.0), oil_viscosity_(0.005),
      water_density_(1000.0), water_viscosity_(0.001),
      gas_density_(1.0), gas_viscosity_(0.00001),
      surface_temp_(288.15), geothermal_grad_(0.03) {}

void VFPTableGenerator::setDeviation(const std::vector<double>& md,
                                      const std::vector<double>& inc) {
    md_ = md;
    inc_ = inc;
}

double VFPTableGenerator::calculateFrictionGradient(double rate, double density,
                                                     double viscosity) const {
    // Darcy-Weisbach friction
    double area = M_PI * wellbore_radius_ * wellbore_radius_;
    double velocity = rate / area;
    
    // Reynolds number
    double diameter = 2.0 * wellbore_radius_;
    double Re = density * std::abs(velocity) * diameter / viscosity;
    
    // Friction factor (Blasius for turbulent)
    double f;
    if (Re < 2300) {
        f = 64.0 / std::max(Re, 1.0);
    } else {
        f = 0.316 / std::pow(Re, 0.25);
    }
    
    // Friction pressure gradient (Pa/m)
    return f * density * velocity * std::abs(velocity) / (2.0 * diameter);
}

double VFPTableGenerator::calculateHydrostaticGradient(double density, double inc) const {
    return density * 9.81 * std::cos(inc * M_PI / 180.0);
}

double VFPTableGenerator::calculateTotalPressureDrop(double thp, double rate, double wct,
                                                      double gor, double alq) const {
    // Simplified calculation - in practice would use correlations like Beggs & Brill
    
    // Calculate mixture properties
    double fw = wct;
    double fo = 1.0 - wct;
    
    double liquid_rate = rate;  // Assuming oil rate
    double gas_rate = rate * gor + alq;  // Gas from GOR + gas lift
    
    double liquid_density = fo * oil_density_ + fw * water_density_;
    double liquid_viscosity = fo * oil_viscosity_ + fw * water_viscosity_;
    
    // Simplified mixture density
    double total_rate = liquid_rate + gas_rate;
    double liquid_holdup = 0.8;  // Simplified
    
    double mixture_density = liquid_holdup * liquid_density + 
                            (1.0 - liquid_holdup) * gas_density_;
    
    // Pressure drop
    double hydrostatic_dp = calculateHydrostaticGradient(mixture_density, 0.0) * total_depth_;
    double friction_dp = calculateFrictionGradient(total_rate, liquid_density, 
                                                    liquid_viscosity) * total_depth_;
    
    return hydrostatic_dp + friction_dp;
}

VFPProdTable VFPTableGenerator::generateProductionTable(
    const std::vector<double>& thp_range,
    const std::vector<double>& rate_range,
    const std::vector<double>& wct_range,
    const std::vector<double>& gor_range,
    const std::vector<double>& alq_range) const {
    
    VFPProdTable table(1);
    table.setDatum(total_depth_);
    table.setTHPValues(thp_range);
    table.setFlowRateValues(rate_range);
    table.setWaterCutValues(wct_range);
    table.setGORValues(gor_range);
    table.setALQValues(alq_range);
    
    // Calculate BHP for each combination
    for (double thp : thp_range) {
        for (double rate : rate_range) {
            for (double wct : wct_range) {
                for (double gor : gor_range) {
                    for (double alq : alq_range) {
                        double dp = calculateTotalPressureDrop(thp, rate, wct, gor, alq);
                        double bhp = thp + dp;
                        
                        VFPDataPoint point;
                        point.thp = thp;
                        point.flow_rate = rate;
                        point.wct = wct;
                        point.gor = gor;
                        point.alq = alq;
                        point.bhp = bhp;
                        
                        table.addDataPoint(point);
                    }
                }
            }
        }
    }
    
    table.buildTableFromPoints();
    return table;
}

VFPInjTable VFPTableGenerator::generateInjectionTable(
    const std::vector<double>& thp_range,
    const std::vector<double>& rate_range,
    const std::string& fluid) const {
    
    VFPInjTable table(1);
    table.setDatum(total_depth_);
    table.setTHPValues(thp_range);
    table.setFlowRateValues(rate_range);
    
    double density, viscosity;
    if (fluid == "WATER") {
        density = water_density_;
        viscosity = water_viscosity_;
        table.setFlowType(VFPFlowType::WATER);
    } else if (fluid == "GAS") {
        density = gas_density_;
        viscosity = gas_viscosity_;
        table.setFlowType(VFPFlowType::GAS);
    } else {
        density = oil_density_;
        viscosity = oil_viscosity_;
        table.setFlowType(VFPFlowType::OIL);
    }
    
    std::vector<double> bhp_data;
    
    for (double thp : thp_range) {
        for (double rate : rate_range) {
            double hydrostatic = calculateHydrostaticGradient(density, 0.0) * total_depth_;
            double friction = calculateFrictionGradient(rate, density, viscosity) * total_depth_;
            double bhp = thp + hydrostatic + friction;
            bhp_data.push_back(bhp);
        }
    }
    
    table.setBHPData(bhp_data);
    return table;
}

// ============================================================================
// GasLiftOptimizer Implementation
// ============================================================================

GasLiftOptimizer::GasLiftOptimizer()
    : vfp_table_(nullptr),
      oil_price_(50.0),
      gas_price_(3.0),
      gaslift_cost_(1.0),
      max_gaslift_(1e6),
      min_bhp_(0.0),
      min_thp_(0.0) {}

void GasLiftOptimizer::setVFPTable(const VFPProdTable& table) {
    vfp_table_ = &table;
}

double GasLiftOptimizer::calculateNetRevenue(double oil_rate, double gas_rate,
                                              double gaslift_rate) const {
    return oil_price_ * oil_rate + gas_price_ * gas_rate - gaslift_cost_ * gaslift_rate;
}

std::pair<double, double> GasLiftOptimizer::optimizeSingleWell(
    double thp, double wct, double gor, double reservoir_pressure) const {
    
    if (!vfp_table_) return {0.0, 0.0};
    
    double best_alq = 0.0;
    double best_revenue = -1e20;
    double best_oil_rate = 0.0;
    
    // Simple grid search
    for (double alq = 0.0; alq <= max_gaslift_; alq += max_gaslift_ / 20.0) {
        // Calculate achievable rate
        double oil_rate = vfp_table_->calculateFlowRate(reservoir_pressure, thp, 
                                                         wct, gor, alq);
        
        double gas_rate = oil_rate * gor;
        double revenue = calculateNetRevenue(oil_rate, gas_rate, alq);
        
        if (revenue > best_revenue) {
            best_revenue = revenue;
            best_alq = alq;
            best_oil_rate = oil_rate;
        }
    }
    
    return {best_alq, best_oil_rate};
}

std::vector<double> GasLiftOptimizer::optimizeField(
    const std::vector<std::tuple<const VFPProdTable*, double, double, double, double>>& wells,
    double total_gas_available) const {
    
    // Greedy allocation based on incremental return
    std::vector<double> alq_allocation(wells.size(), 0.0);
    double remaining_gas = total_gas_available;
    double step_size = total_gas_available / 100.0;
    
    while (remaining_gas > step_size) {
        double best_gain = 0.0;
        size_t best_well = 0;
        
        for (size_t w = 0; w < wells.size(); ++w) {
            const auto& [table, thp, wct, gor, pr] = wells[w];
            if (!table) continue;
            
            double current_rate = table->calculateFlowRate(pr, thp, wct, gor, 
                                                           alq_allocation[w]);
            double new_rate = table->calculateFlowRate(pr, thp, wct, gor,
                                                        alq_allocation[w] + step_size);
            
            double gain = (new_rate - current_rate) * oil_price_ - step_size * gaslift_cost_;
            
            if (gain > best_gain) {
                best_gain = gain;
                best_well = w;
            }
        }
        
        if (best_gain <= 0.0) break;  // No profitable allocation
        
        alq_allocation[best_well] += step_size;
        remaining_gas -= step_size;
    }
    
    return alq_allocation;
}

// ============================================================================
// VFPTableManager Implementation
// ============================================================================

VFPTableManager::VFPTableManager() {}

void VFPTableManager::addProductionTable(std::unique_ptr<VFPProdTable> table) {
    int num = table->getTableNumber();
    prod_tables_[num] = std::move(table);
}

void VFPTableManager::addInjectionTable(std::unique_ptr<VFPInjTable> table) {
    int num = table->getTableNumber();
    inj_tables_[num] = std::move(table);
}

const VFPProdTable* VFPTableManager::getProductionTable(int table_num) const {
    auto it = prod_tables_.find(table_num);
    return (it != prod_tables_.end()) ? it->second.get() : nullptr;
}

const VFPInjTable* VFPTableManager::getInjectionTable(int table_num) const {
    auto it = inj_tables_.find(table_num);
    return (it != inj_tables_.end()) ? it->second.get() : nullptr;
}

double VFPTableManager::calculateProductionBHP(int table_num, double thp, double rate,
                                                double wct, double gor, double alq) const {
    const auto* table = getProductionTable(table_num);
    if (!table) return thp;
    return table->calculateBHP(thp, rate, wct, gor, alq);
}

double VFPTableManager::calculateInjectionBHP(int table_num, double thp, double rate) const {
    const auto* table = getInjectionTable(table_num);
    if (!table) return thp;
    return table->calculateBHP(thp, rate);
}

double VFPTableManager::calculateProductionRate(int table_num, double bhp, double thp,
                                                 double wct, double gor, double alq) const {
    const auto* table = getProductionTable(table_num);
    if (!table) return 0.0;
    return table->calculateFlowRate(bhp, thp, wct, gor, alq);
}

double VFPTableManager::calculateInjectionRate(int table_num, double bhp, double thp) const {
    const auto* table = getInjectionTable(table_num);
    if (!table) return 0.0;
    return table->calculateFlowRate(bhp, thp);
}

void VFPTableManager::parseVFPPROD(const std::vector<std::string>& lines) {
    // Parse ECLIPSE VFPPROD format
    // This is a simplified parser - full implementation would handle all formats
    
    if (lines.empty()) return;
    
    auto table = std::make_unique<VFPProdTable>();
    
    // Parse header line for table number and flow type
    std::istringstream header(lines[0]);
    int table_num;
    header >> table_num;
    table->setTableNumber(table_num);
    
    // Parse remaining data...
    // (Full implementation would parse all axis values and BHP data)
    
    addProductionTable(std::move(table));
}

void VFPTableManager::parseVFPINJ(const std::vector<std::string>& lines) {
    if (lines.empty()) return;
    
    auto table = std::make_unique<VFPInjTable>();
    
    std::istringstream header(lines[0]);
    int table_num;
    header >> table_num;
    table->setTableNumber(table_num);
    
    addInjectionTable(std::move(table));
}

// ============================================================================
// Parse functions
// ============================================================================

VFPFlowType parseVFPFlowType(const std::string& str) {
    if (str == "OIL") return VFPFlowType::OIL;
    if (str == "LIQ") return VFPFlowType::LIQ;
    if (str == "GAS") return VFPFlowType::GAS;
    if (str == "WATER" || str == "WAT") return VFPFlowType::WATER;
    return VFPFlowType::OIL;
}

VFPWaterCutType parseVFPWaterCutType(const std::string& str) {
    if (str == "WCT") return VFPWaterCutType::WCT;
    if (str == "WGR") return VFPWaterCutType::WGR;
    return VFPWaterCutType::WCT;
}

VFPGasLiquidType parseVFPGasLiquidType(const std::string& str) {
    if (str == "GOR") return VFPGasLiquidType::GOR;
    if (str == "GLR") return VFPGasLiquidType::GLR;
    if (str == "OGR") return VFPGasLiquidType::OGR;
    return VFPGasLiquidType::GOR;
}

VFPArtificialLiftType parseVFPArtificialLiftType(const std::string& str) {
    if (str == "NONE" || str == "") return VFPArtificialLiftType::NONE;
    if (str == "GRAT" || str == "GASLIFT") return VFPArtificialLiftType::GASLIFT;
    if (str == "ESP") return VFPArtificialLiftType::ESP;
    if (str == "CHOKE") return VFPArtificialLiftType::CHOKE;
    return VFPArtificialLiftType::NONE;
}

} // namespace FSRM
