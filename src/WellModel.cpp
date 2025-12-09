#include "WellModel.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace FSRM {

// Constants
constexpr double PI = M_PI;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;

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

// ============================================================================
// DirectionalWell
// ============================================================================

DirectionalWell::DirectionalWell(const std::string& name, WellType type)
    : WellModel(name, type), kickoff_depth(0.0), build_rate(0.0), 
      max_inclination(0.0) {
    // Initialize with surface point
    SurveyPoint surface;
    surface.measured_depth = 0.0;
    surface.inclination = 0.0;
    surface.azimuth = 0.0;
    surface.tvd = 0.0;
    surface.north_offset = 0.0;
    surface.east_offset = 0.0;
    surface.dogleg_severity = 0.0;
    survey_data.push_back(surface);
}

void DirectionalWell::addSurveyPoint(double md, double inc, double azi) {
    SurveyPoint point;
    point.measured_depth = md;
    point.inclination = inc;
    point.azimuth = azi;
    point.tvd = 0.0;  // Will be calculated
    point.north_offset = 0.0;
    point.east_offset = 0.0;
    point.dogleg_severity = 0.0;
    survey_data.push_back(point);
}

void DirectionalWell::addSurveyPoint(const SurveyPoint& point) {
    survey_data.push_back(point);
}

double DirectionalWell::computeDogleg(const SurveyPoint& p1, const SurveyPoint& p2) const {
    // Dogleg angle between two survey points using spherical trigonometry
    double inc1 = p1.inclination * DEG_TO_RAD;
    double inc2 = p2.inclination * DEG_TO_RAD;
    double azi1 = p1.azimuth * DEG_TO_RAD;
    double azi2 = p2.azimuth * DEG_TO_RAD;
    
    double cos_dogleg = std::cos(inc2 - inc1) - 
                       std::sin(inc1) * std::sin(inc2) * (1.0 - std::cos(azi2 - azi1));
    
    // Clamp to valid range to avoid numerical issues
    cos_dogleg = std::max(-1.0, std::min(1.0, cos_dogleg));
    
    return std::acos(cos_dogleg) * RAD_TO_DEG;
}

void DirectionalWell::computeMinimumCurvature(const SurveyPoint& p1, 
                                              const SurveyPoint& p2,
                                              SurveyPoint& p2_out) {
    // Minimum curvature method for computing TVD and offsets
    p2_out = p2;
    
    double md_diff = p2.measured_depth - p1.measured_depth;
    if (md_diff < 1e-10) {
        p2_out.tvd = p1.tvd;
        p2_out.north_offset = p1.north_offset;
        p2_out.east_offset = p1.east_offset;
        return;
    }
    
    double inc1 = p1.inclination * DEG_TO_RAD;
    double inc2 = p2.inclination * DEG_TO_RAD;
    double azi1 = p1.azimuth * DEG_TO_RAD;
    double azi2 = p2.azimuth * DEG_TO_RAD;
    
    // Compute dogleg angle
    double dogleg = computeDogleg(p1, p2) * DEG_TO_RAD;
    
    // Ratio factor
    double rf = 1.0;
    if (dogleg > 1e-6) {
        rf = 2.0 / dogleg * std::tan(dogleg / 2.0);
    }
    
    // Compute increments
    double delta_north = 0.5 * md_diff * (std::sin(inc1) * std::cos(azi1) + 
                                          std::sin(inc2) * std::cos(azi2)) * rf;
    double delta_east = 0.5 * md_diff * (std::sin(inc1) * std::sin(azi1) + 
                                         std::sin(inc2) * std::sin(azi2)) * rf;
    double delta_tvd = 0.5 * md_diff * (std::cos(inc1) + std::cos(inc2)) * rf;
    
    p2_out.tvd = p1.tvd + delta_tvd;
    p2_out.north_offset = p1.north_offset + delta_north;
    p2_out.east_offset = p1.east_offset + delta_east;
    
    // Dogleg severity in degrees per 30m
    p2_out.dogleg_severity = (dogleg * RAD_TO_DEG) * 30.0 / md_diff;
}

void DirectionalWell::calculatePositions() {
    if (survey_data.size() < 2) return;
    
    for (size_t i = 1; i < survey_data.size(); ++i) {
        computeMinimumCurvature(survey_data[i-1], survey_data[i], survey_data[i]);
    }
}

void DirectionalWell::calculateDoglegSeverity() {
    if (survey_data.size() < 2) return;
    
    for (size_t i = 1; i < survey_data.size(); ++i) {
        double dogleg = computeDogleg(survey_data[i-1], survey_data[i]);
        double md_diff = survey_data[i].measured_depth - survey_data[i-1].measured_depth;
        if (md_diff > 1e-10) {
            survey_data[i].dogleg_severity = dogleg * 30.0 / md_diff;
        }
    }
}

void DirectionalWell::computeTrajectory() {
    calculatePositions();
    calculateDoglegSeverity();
    
    trajectory.clear();
    for (size_t i = 1; i < survey_data.size(); ++i) {
        TrajectorySegment seg;
        seg.start = survey_data[i-1];
        seg.end = survey_data[i];
        seg.length = survey_data[i].measured_depth - survey_data[i-1].measured_depth;
        seg.avg_inclination = 0.5 * (survey_data[i].inclination + survey_data[i-1].inclination);
        seg.avg_azimuth = 0.5 * (survey_data[i].azimuth + survey_data[i-1].azimuth);
        trajectory.push_back(seg);
    }
}

std::vector<TrajectorySegment> DirectionalWell::getTrajectorySegments() const {
    return trajectory;
}

double DirectionalWell::getTotalMeasuredDepth() const {
    if (survey_data.empty()) return 0.0;
    return survey_data.back().measured_depth;
}

double DirectionalWell::getTrueVerticalDepth() const {
    if (survey_data.empty()) return 0.0;
    return survey_data.back().tvd;
}

double DirectionalWell::getMaxInclination() const {
    double max_inc = 0.0;
    for (const auto& point : survey_data) {
        max_inc = std::max(max_inc, point.inclination);
    }
    return max_inc;
}

SurveyPoint DirectionalWell::interpolateSurvey(double md) const {
    if (survey_data.empty()) {
        throw std::runtime_error("No survey data available");
    }
    
    if (md <= survey_data.front().measured_depth) {
        return survey_data.front();
    }
    if (md >= survey_data.back().measured_depth) {
        return survey_data.back();
    }
    
    // Find bracketing points
    for (size_t i = 1; i < survey_data.size(); ++i) {
        if (md <= survey_data[i].measured_depth) {
            // Linear interpolation
            const auto& p1 = survey_data[i-1];
            const auto& p2 = survey_data[i];
            double t = (md - p1.measured_depth) / (p2.measured_depth - p1.measured_depth);
            
            SurveyPoint interp;
            interp.measured_depth = md;
            interp.inclination = p1.inclination + t * (p2.inclination - p1.inclination);
            interp.azimuth = p1.azimuth + t * (p2.azimuth - p1.azimuth);
            interp.tvd = p1.tvd + t * (p2.tvd - p1.tvd);
            interp.north_offset = p1.north_offset + t * (p2.north_offset - p1.north_offset);
            interp.east_offset = p1.east_offset + t * (p2.east_offset - p1.east_offset);
            interp.dogleg_severity = p2.dogleg_severity;
            
            return interp;
        }
    }
    
    return survey_data.back();
}

void DirectionalWell::computeGridIntersections(int nx, int ny, int nz,
                                               double dx, double dy, double dz) {
    cell_intersection_lengths.clear();
    
    for (const auto& seg : trajectory) {
        // Simple ray-grid intersection
        // For each segment, compute which cells it passes through
        // This is a simplified version - full implementation would use 
        // a proper ray-tracing algorithm
        
        double x1 = seg.start.east_offset;
        double y1 = seg.start.north_offset;
        double z1 = seg.start.tvd;
        
        double x2 = seg.end.east_offset;
        double y2 = seg.end.north_offset;
        double z2 = seg.end.tvd;
        
        // Convert to grid indices (approximate)
        int i1 = static_cast<int>(x1 / dx);
        int j1 = static_cast<int>(y1 / dy);
        int k1 = static_cast<int>(z1 / dz);
        
        int i2 = static_cast<int>(x2 / dx);
        int j2 = static_cast<int>(y2 / dy);
        int k2 = static_cast<int>(z2 / dz);
        
        // Add all cells along the path (simplified - should use DDA or similar)
        for (int i = std::min(i1, i2); i <= std::max(i1, i2) && i < nx; ++i) {
            for (int j = std::min(j1, j2); j <= std::max(j1, j2) && j < ny; ++j) {
                for (int k = std::min(k1, k2); k <= std::max(k1, k2) && k < nz; ++k) {
                    if (i >= 0 && j >= 0 && k >= 0) {
                        int cell_id = i + j * nx + k * nx * ny;
                        cell_intersection_lengths[cell_id] += seg.length / 
                            ((std::abs(i2-i1)+1) * (std::abs(j2-j1)+1) * (std::abs(k2-k1)+1));
                    }
                }
            }
        }
    }
}

std::vector<int> DirectionalWell::getCellsIntersected() const {
    std::vector<int> cells;
    for (const auto& pair : cell_intersection_lengths) {
        cells.push_back(pair.first);
    }
    return cells;
}

double DirectionalWell::computeWellIndex(int completion_idx,
                                        double kx, double ky, double kz,
                                        double dx, double dy, double dz) const {
    // For directional wells, use anisotropic well index formulation
    if (completion_idx >= (int)completions.size()) return 0.0;
    
    const auto& comp = completions[completion_idx];
    
    // Get well angle at this completion
    double inclination = 0.0;
    if (!survey_data.empty()) {
        // Find survey data near this completion depth
        inclination = getMaxInclination();  // Simplified
    }
    
    double inc_rad = inclination * DEG_TO_RAD;
    double rw = comp.diameter / 2.0;
    
    // Anisotropic well index for deviated wells
    double kh = std::sqrt(kx * ky);
    double kv = kz;
    // double anisotropy = std::sqrt(kv / kh);  // Not used in this simplified implementation
    
    // Effective wellbore length in cell
    double L_eff = dz / std::cos(inc_rad);
    if (std::abs(std::cos(inc_rad)) < 1e-6) {
        L_eff = std::sqrt(dx*dx + dy*dy);
    }
    
    // Deviated well index (Peaceman for deviated wells)
    double re = 0.28 * std::sqrt(dx*dx + dy*dy);
    double WI = 2.0 * PI * std::sqrt(kh * kv) * L_eff / 
                (std::log(re / rw) + comp.skin_factor);
    
    return WI;
}

double DirectionalWell::computeDeviatedWellPI(double length, double angle,
                                              double kh, double kv,
                                              double h, double rw) const {
    // Productivity index for deviated well
    // double inc_rad = angle * DEG_TO_RAD;  // Not used in simplified Joshi equation
    double anisotropy = std::sqrt(kv / kh);
    (void)angle;  // Suppress unused parameter warning
    
    // Joshi equation for deviated wells
    double a = length / 2.0;
    double I = anisotropy * h / 2.0;
    
    if (a < 1e-10) return 0.0;
    
    double numerator = 2.0 * PI * std::sqrt(kh * kv) * length;
    double denominator = std::log((a + std::sqrt(a*a + I*I)) / rw);
    
    return numerator / denominator;
}

void DirectionalWell::setKickoffDepth(double depth) {
    kickoff_depth = depth;
}

void DirectionalWell::setBuildRate(double rate_deg_per_30m) {
    build_rate = rate_deg_per_30m;
}

void DirectionalWell::setMaxInclination(double max_inc) {
    max_inclination = max_inc;
}

void DirectionalWell::buildJCurve(double kickoff_md, double build_rate_val, 
                                  double target_inc) {
    // J-curve: vertical, then build to target inclination, then hold
    survey_data.clear();
    
    // Surface point
    addSurveyPoint(0.0, 0.0, 0.0);
    
    // Kickoff point
    if (kickoff_md > 0.0) {
        addSurveyPoint(kickoff_md, 0.0, 0.0);
    }
    
    // Build section
    double build_length = target_inc / build_rate_val * 30.0;
    int n_points = 10;
    for (int i = 1; i <= n_points; ++i) {
        double frac = static_cast<double>(i) / n_points;
        double md = kickoff_md + frac * build_length;
        double inc = frac * target_inc;
        addSurveyPoint(md, inc, 0.0);
    }
    
    kickoff_depth = kickoff_md;
    build_rate = build_rate_val;
    max_inclination = target_inc;
    
    computeTrajectory();
}

void DirectionalWell::buildSCurve(double kickoff_md, double build_rate_val,
                                  double target_inc, double hold_md, 
                                  double drop_rate) {
    // S-curve: vertical, build, hold, drop back to vertical
    survey_data.clear();
    
    addSurveyPoint(0.0, 0.0, 0.0);
    
    if (kickoff_md > 0.0) {
        addSurveyPoint(kickoff_md, 0.0, 0.0);
    }
    
    // Build section
    double build_length = target_inc / build_rate_val * 30.0;
    int n_build = 10;
    for (int i = 1; i <= n_build; ++i) {
        double frac = static_cast<double>(i) / n_build;
        double md = kickoff_md + frac * build_length;
        double inc = frac * target_inc;
        addSurveyPoint(md, inc, 0.0);
    }
    
    // Hold section
    double hold_start = kickoff_md + build_length;
    if (hold_md > 0.0) {
        addSurveyPoint(hold_start + hold_md, target_inc, 0.0);
    }
    
    // Drop section
    double drop_length = target_inc / drop_rate * 30.0;
    int n_drop = 10;
    for (int i = 1; i <= n_drop; ++i) {
        double frac = static_cast<double>(i) / n_drop;
        double md = hold_start + hold_md + frac * drop_length;
        double inc = target_inc * (1.0 - frac);
        addSurveyPoint(md, inc, 0.0);
    }
    
    kickoff_depth = kickoff_md;
    build_rate = build_rate_val;
    max_inclination = target_inc;
    
    computeTrajectory();
}

void DirectionalWell::buildSlantWell(double kickoff_md, double target_inc, 
                                     double target_azi) {
    // Slant well: build to target angle and azimuth, then hold
    survey_data.clear();
    
    addSurveyPoint(0.0, 0.0, 0.0);
    
    if (kickoff_md > 0.0) {
        addSurveyPoint(kickoff_md, 0.0, 0.0);
    }
    
    // Build to target with default build rate
    double default_build_rate = 2.0;  // degrees per 30m
    double build_length = target_inc / default_build_rate * 30.0;
    
    int n_points = 10;
    for (int i = 1; i <= n_points; ++i) {
        double frac = static_cast<double>(i) / n_points;
        double md = kickoff_md + frac * build_length;
        double inc = frac * target_inc;
        double azi = frac * target_azi;
        addSurveyPoint(md, inc, azi);
    }
    
    kickoff_depth = kickoff_md;
    build_rate = default_build_rate;
    max_inclination = target_inc;
    
    computeTrajectory();
}

} // namespace FSRM
