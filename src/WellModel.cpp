#include "WellModel.hpp"
#include <cmath>
#include <algorithm>

namespace ResSim {

// ============================================================================
// WellModel Base Class
// ============================================================================

WellModel::WellModel(const std::string& name, WellType type)
    : well_name(name), well_type(type), 
      control_mode(WellControlMode::RATE_CONTROL), target_value(0.0),
      reference_depth(0.0), use_wellbore_model(false),
      wellbore_inner_radius(0.1), wellbore_outer_radius(0.15),
      wellbore_roughness(1e-5), total_rate(0.0), bhp(1e5),
      cumulative_production(0.0) {
    
    constraints.max_rate = 1e10;
    constraints.min_rate = 0.0;
    constraints.max_bhp = 1e8;
    constraints.min_bhp = 1e5;
    constraints.max_thp = 1e8;
    constraints.min_thp = 0.0;
}

void WellModel::addCompletion(const WellCompletion& comp) {
    completions.push_back(comp);
}

void WellModel::setControl(WellControlMode mode, double target) {
    control_mode = mode;
    target_value = target;
}

void WellModel::setConstraints(const WellConstraints& constr) {
    constraints = constr;
}

void WellModel::setReferenceDepth(double depth) {
    reference_depth = depth;
}

double WellModel::computePeacemanRadius(double dx, double dy, 
                                        double kx, double ky) const {
    // Peaceman well model equivalent radius
    double r0 = 0.28 * std::sqrt(std::sqrt(ky/kx) * dx * dx + 
                                 std::sqrt(kx/ky) * dy * dy) /
                (std::pow(ky/kx, 0.25) + std::pow(kx/ky, 0.25));
    return r0;
}

double WellModel::computeWellIndex(int completion_idx, 
                                   double kx, double ky, double kz,
                                   double dx, double dy, double dz) const {
    if (completion_idx >= (int)completions.size()) return 0.0;
    
    const auto& comp = completions[completion_idx];
    
    // Horizontal well index for vertical well
    double k_h = std::sqrt(kx * ky);
    double r0 = computePeacemanRadius(dx, dy, kx, ky);
    double rw = comp.diameter / 2.0;
    
    // Well index: WI = 2*pi*k*h / (ln(r0/rw) + s)
    double height = dz;
    double WI = 2.0 * M_PI * k_h * height / 
                (std::log(r0 / rw) + comp.skin_factor);
    
    return WI;
}

double WellModel::computeBottomholePressure(double reservoir_pressure,
                                            double rate, int comp_idx) const {
    if (comp_idx >= (int)completions.size()) return reservoir_pressure;
    
    const auto& comp = completions[comp_idx];
    
    // Simplified: P_bhp = P_res - rate / WI
    double WI = comp.well_index;
    if (WI < 1e-10) WI = 1.0;
    
    return reservoir_pressure - rate / WI;
}

double WellModel::computeRate(double reservoir_pressure, 
                              double bhp, int comp_idx) const {
    if (comp_idx >= (int)completions.size()) return 0.0;
    
    const auto& comp = completions[comp_idx];
    
    // Rate = WI * (P_res - P_bhp)
    double rate = comp.well_index * (reservoir_pressure - bhp);
    
    // Apply constraints
    rate = std::max(constraints.min_rate, 
                   std::min(constraints.max_rate, rate));
    
    return rate;
}

void WellModel::computePhaseRates(const std::vector<double>& res_pressures,
                                   const std::vector<double>& saturations,
                                   double bhp,
                                   std::vector<double>& phase_rates) const {
    // Base implementation - override in derived classes
    phase_rates.resize(1);
    phase_rates[0] = 0.0;
    
    for (size_t i = 0; i < completions.size(); ++i) {
        if (i < res_pressures.size()) {
            phase_rates[0] += computeRate(res_pressures[i], bhp, i);
        }
    }
}

void WellModel::contributeToResidual(Vec F, Vec U, DM dm) const {
    // Add well contributions to residual
    // This would involve getting local values from U at well locations
    // and adding source/sink terms to F
}

void WellModel::contributeToJacobian(Mat J, Vec U, DM dm) const {
    // Add well Jacobian contributions
}

void WellModel::enableWellboreModel(bool enable) {
    use_wellbore_model = enable;
}

void WellModel::setWellboreGeometry(double inner_radius, double outer_radius) {
    wellbore_inner_radius = inner_radius;
    wellbore_outer_radius = outer_radius;
}

void WellModel::setWellboreRoughness(double roughness) {
    wellbore_roughness = roughness;
}

void WellModel::computeWellborePressureDrop(const std::vector<double>& phase_rates,
                                            double depth_interval,
                                            double& dp) const {
    // Drift flux model for multiphase flow
    // Simplified implementation
    
    double total_rate = 0.0;
    for (double rate : phase_rates) {
        total_rate += rate;
    }
    
    // Hydrostatic component
    double rho_avg = 1000.0;  // Average mixture density
    double g = 9.81;
    double dp_hydrostatic = rho_avg * g * depth_interval;
    
    // Frictional component (Darcy-Weisbach)
    double velocity = total_rate / (M_PI * wellbore_inner_radius * wellbore_inner_radius);
    double mu_avg = 0.001;  // Average viscosity
    double Re = rho_avg * velocity * 2.0 * wellbore_inner_radius / mu_avg;
    
    double friction_factor = 64.0 / Re;  // Laminar flow
    if (Re > 2300.0) {
        // Turbulent: Colebrook-White
        double relative_roughness = wellbore_roughness / (2.0 * wellbore_inner_radius);
        friction_factor = std::pow(1.0 / (-2.0 * std::log10(relative_roughness / 3.7 + 
                                    5.74 / std::pow(Re, 0.9))), 2.0);
    }
    
    double dp_friction = friction_factor * rho_avg * velocity * velocity * 
                        depth_interval / (2.0 * 2.0 * wellbore_inner_radius);
    
    dp = dp_hydrostatic + dp_friction;
}

void WellModel::updatePerformance(double dt, double rate) {
    total_rate = rate;
    cumulative_production += rate * dt;
}

// ============================================================================
// ProductionWell
// ============================================================================

ProductionWell::ProductionWell(const std::string& name)
    : WellModel(name, WellType::PRODUCER) {}

void ProductionWell::setTargetOilRate(double rate) {
    control_mode = WellControlMode::RATE_CONTROL;
    target_value = rate;
}

void ProductionWell::setTargetGasRate(double rate) {
    control_mode = WellControlMode::RATE_CONTROL;
    target_value = rate;
}

void ProductionWell::setTargetLiquidRate(double rate) {
    control_mode = WellControlMode::RATE_CONTROL;
    target_value = rate;
}

void ProductionWell::setTargetReservoirVoidage(double rate) {
    control_mode = WellControlMode::RESERVOIR_VOIDAGE;
    target_value = rate;
}

void ProductionWell::computePhaseRates(const std::vector<double>& res_pressures,
                                       const std::vector<double>& saturations,
                                       double bhp,
                                       std::vector<double>& phase_rates) const {
    // Compute oil, water, gas rates
    phase_rates.resize(3, 0.0);
    
    for (size_t i = 0; i < completions.size(); ++i) {
        if (i < res_pressures.size() && i < saturations.size() / 3) {
            double total_rate = computeRate(res_pressures[i], bhp, i);
            
            // Split by mobility
            double So = saturations[3*i];
            double Sw = saturations[3*i + 1];
            double Sg = saturations[3*i + 2];
            
            // Simplified - should use relative permeabilities
            double total_sat = So + Sw + Sg;
            if (total_sat > 0.0) {
                phase_rates[0] += total_rate * So / total_sat;  // Oil
                phase_rates[1] += total_rate * Sw / total_sat;  // Water
                phase_rates[2] += total_rate * Sg / total_sat;  // Gas
            }
        }
    }
}

// ============================================================================
// InjectionWell
// ============================================================================

InjectionWell::InjectionWell(const std::string& name)
    : WellModel(name, WellType::INJECTOR),
      injection_fluid("WATER"), injection_temperature(300.0) {}

void InjectionWell::setInjectionFluid(const std::string& fluid_type) {
    injection_fluid = fluid_type;
}

void InjectionWell::setTargetInjectionRate(double rate) {
    control_mode = WellControlMode::RATE_CONTROL;
    target_value = rate;
}

void InjectionWell::setInjectionTemperature(double temp) {
    injection_temperature = temp;
}

void InjectionWell::setInjectionComposition(const std::vector<double>& composition) {
    injection_composition = composition;
}

// ============================================================================
// HorizontalWell
// ============================================================================

HorizontalWell::HorizontalWell(const std::string& name, WellType type)
    : WellModel(name, type), use_icds(false) {}

void HorizontalWell::setWellPath(const std::vector<std::vector<double>>& path) {
    well_path = path;
}

void HorizontalWell::enableInflowControlDevices(bool enable) {
    use_icds = enable;
}

double HorizontalWell::computeWellIndex(int completion_idx,
                                        double kx, double ky, double kz,
                                        double dx, double dy, double dz) const {
    // Horizontal well model (Joshi, Babu-Odeh, etc.)
    
    if (completion_idx >= (int)completions.size()) return 0.0;
    const auto& comp = completions[completion_idx];
    
    // Simplified anisotropic horizontal well index
    double L = dz;  // Lateral length in this cell
    double h = dy;  // Reservoir thickness
    double rw = comp.diameter / 2.0;
    
    double k_v = kz;
    double k_h = std::sqrt(kx * ky);
    double anisotropy = std::sqrt(k_v / k_h);
    
    // Joshi equation
    double a = L / 2.0;
    double I = anisotropy * h / 2.0;
    
    double numerator = 2.0 * M_PI * std::sqrt(k_h * k_v) * L;
    double denominator = std::log((a + std::sqrt(a*a + I*I)) / rw) + 
                        (k_h / k_v) * std::log(h / (2.0 * rw));
    
    double WI = numerator / denominator;
    
    return WI;
}

// ============================================================================
// MultilateralWell
// ============================================================================

MultilateralWell::MultilateralWell(const std::string& name, WellType type)
    : WellModel(name, type) {}

void MultilateralWell::addLateral(int lateral_id, const std::vector<WellCompletion>& comps) {
    laterals[lateral_id] = comps;
    lateral_status[lateral_id] = true;
}

void MultilateralWell::setLateralControl(int lateral_id, bool is_open) {
    if (lateral_status.find(lateral_id) != lateral_status.end()) {
        lateral_status[lateral_id] = is_open;
    }
}

} // namespace ResSim
