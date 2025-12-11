#include "WellboreHydraulics.hpp"
#include <stdexcept>
#include <limits>
#include <numeric>
#
namespace FSRM {
#
static constexpr double PI = 3.141592653589793238462643383279502884;
static constexpr double G = 9.80665; // m/s^2 (standard gravity)
#
// ============================================================================
// WellboreSegment
// ============================================================================
#
WellboreSegment::WellboreSegment()
    : start_md(0.0),
      end_md(0.0),
      start_tvd(0.0),
      end_tvd(0.0),
      inclination(0.0),
      azimuth(0.0),
      inner_diameter(0.1),
      outer_diameter(0.12),
      roughness(1e-5),
      thermal_conductivity(50.0),
      cement_conductivity(1.0),
      formation_conductivity(2.5),
      flow_area(0.0),
      hydraulic_diameter(0.0),
      wetted_perimeter(0.0) {
    calculateDerivedProperties();
}
#
void WellboreSegment::calculateDerivedProperties() {
    if (inner_diameter <= 0.0) {
        flow_area = 0.0;
        hydraulic_diameter = 0.0;
        wetted_perimeter = 0.0;
        return;
    }
    flow_area = PI * inner_diameter * inner_diameter / 4.0;
    hydraulic_diameter = inner_diameter;
    wetted_perimeter = PI * inner_diameter;
}
#
// ============================================================================
// WellboreCompletion
// ============================================================================
#
WellboreCompletion::WellboreCompletion()
    : type(CompletionType::CASED_CEMENTED),
      casing_id(0.2),
      casing_od(0.25),
      tubing_id(0.1),
      tubing_od(0.1143),
      annulus_area(0.0),
      has_tubing(true),
      tubing_setting_depth(0.0) {
    if (casing_id > 0.0 && tubing_od > 0.0) {
        annulus_area = PI * (casing_id * casing_id - tubing_od * tubing_od) / 4.0;
    }
}
#
// ============================================================================
// FrictionModel
// ============================================================================
#
FrictionModel::FrictionModel() : correlation_(Correlation::HAALAND) {}
#
static double clamp01(double x) {
    return std::max(0.0, std::min(1.0, x));
}
#
double FrictionModel::colebrookWhite(double reynolds, double rel_roughness) const {
    // Solve for Darcy friction factor using Colebrook-White:
    // 1/sqrt(f) = -2 log10( (eps/D)/3.7 + 2.51/(Re*sqrt(f)) )
    // Use a robust fixed-point iteration on 1/sqrt(f).
    if (reynolds <= 0.0) return std::numeric_limits<double>::infinity();
    rel_roughness = std::max(0.0, rel_roughness);
#
    // Initial guess: Haaland
    const double inv_sqrt_f0 = -1.8 * std::log10(std::pow(rel_roughness / 3.7, 1.11) + 6.9 / reynolds);
    double inv_sqrt_f = std::max(1e-6, inv_sqrt_f0);
#
    for (int it = 0; it < 25; ++it) {
        const double rhs = -2.0 * std::log10(rel_roughness / 3.7 + 2.51 / (reynolds * inv_sqrt_f));
        const double next = 0.5 * inv_sqrt_f + 0.5 * rhs; // damping
        if (std::abs(next - inv_sqrt_f) < 1e-10) {
            inv_sqrt_f = next;
            break;
        }
        inv_sqrt_f = next;
    }
#
    const double f_darcy = 1.0 / (inv_sqrt_f * inv_sqrt_f);
    return std::max(1e-6, f_darcy);
}
#
double FrictionModel::fanningFrictionFactor(double reynolds, double relative_roughness) const {
    // Return Fanning friction factor f_F (Darcy f_D = 4 f_F).
    if (!std::isfinite(reynolds) || reynolds <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
#
    // Laminar (pipe): f_F = 16 / Re
    if (reynolds < 2100.0) {
        return 16.0 / reynolds;
    }
#
    relative_roughness = std::max(0.0, relative_roughness);
#
    double f_darcy = 0.0;
    switch (correlation_) {
        case Correlation::COLEBROOK_WHITE: {
            f_darcy = colebrookWhite(reynolds, relative_roughness);
            break;
        }
        case Correlation::HAALAND: {
            const double inv_sqrt_f = -1.8 * std::log10(std::pow(relative_roughness / 3.7, 1.11) + 6.9 / reynolds);
            f_darcy = 1.0 / (inv_sqrt_f * inv_sqrt_f);
            break;
        }
        case Correlation::SWAMEE_JAIN: {
            // Swamee-Jain (Darcy)
            const double term = relative_roughness / 3.7 + 5.74 / std::pow(reynolds, 0.9);
            const double inv_sqrt_f = -2.0 * std::log10(term);
            f_darcy = 1.0 / (inv_sqrt_f * inv_sqrt_f);
            break;
        }
        case Correlation::CHEN: {
            // Chen (Darcy), explicit approximation
            const double A = -2.0 * std::log10(relative_roughness / 3.7065 - 5.0452 / reynolds *
                                               std::log10(std::pow(relative_roughness, 1.1098) / 2.8257 +
                                                          std::pow(5.8506 / reynolds, 0.8981)));
            f_darcy = 1.0 / (A * A);
            break;
        }
    }
#
    f_darcy = std::max(1e-6, f_darcy);
    return 0.25 * f_darcy;
}
#
double FrictionModel::newtonianFrictionGradient(double velocity, double density,
                                                double viscosity, double diameter,
                                                double roughness) const {
    if (diameter <= 0.0) return std::numeric_limits<double>::infinity();
    if (!std::isfinite(velocity) || !std::isfinite(density) || !std::isfinite(viscosity)) {
        return std::numeric_limits<double>::infinity();
    }
    if (density <= 0.0 || viscosity <= 0.0) return 0.0;
#
    const double Re = density * std::abs(velocity) * diameter / viscosity;
    const double rel_eps = std::max(0.0, roughness) / diameter;
    const double fF = fanningFrictionFactor(Re, rel_eps);
#
    // Using Fanning: dp/dL = 2 fF rho v^2 / D
    return 2.0 * fF * density * velocity * velocity / diameter;
}
#
double FrictionModel::powerLawFrictionGradient(double velocity, double density,
                                               double k_prime, double n_prime,
                                               double diameter, double roughness) const {
    // Generalized Newtonian power-law fluid.
    // Laminar: use Metzner-Reed generalized Reynolds number and f_F = 16 / Re_g.
    // Turbulent: use a simple Dodge-Metzner style approximation (explicit).
    if (diameter <= 0.0) return std::numeric_limits<double>::infinity();
    if (density <= 0.0 || k_prime <= 0.0 || n_prime <= 0.0) return 0.0;
#
    const double v = std::abs(velocity);
    if (v < 1e-16) return 0.0;
#
    const double n = n_prime;
    const double K = k_prime;
#
    // Generalized Reynolds number for power-law in pipe
    // Re_g = (rho * v^(2-n) * D^n) / (K * 8^(n-1) * ( (3n+1)/(4n) )^n )
    const double C = std::pow(8.0, n - 1.0) * std::pow((3.0 * n + 1.0) / (4.0 * n), n);
    const double Re_g = (density * std::pow(v, 2.0 - n) * std::pow(diameter, n)) / (K * C);
    const double rel_eps = std::max(0.0, roughness) / diameter;
#
    double fF = 0.0;
    if (Re_g < 2100.0) {
        fF = 16.0 / std::max(1e-12, Re_g);
    } else {
        // Explicit approximation for turbulent power-law (very common engineering form).
        // This is not a full-featured implementation, but is materially better than
        // treating power-law fluids as Newtonian with constant viscosity.
        //
        // Use a Blasius-like exponent that depends on n, with a roughness correction
        // blended in via Haaland on an "effective" Re.
        const double m = 0.25 / std::max(0.2, n); // heuristic
        const double Re_eff = std::max(1e-12, Re_g);
        // Start with smooth pipe (Darcy) and blend with roughness
        const double fD_smooth = 0.3164 / std::pow(Re_eff, m);
        const double inv_sqrt_f_rough = -1.8 * std::log10(std::pow(rel_eps / 3.7, 1.11) + 6.9 / Re_eff);
        const double fD_rough = 1.0 / (inv_sqrt_f_rough * inv_sqrt_f_rough);
        const double w = clamp01((rel_eps * 1e3)); // roughness weight
        const double fD = (1.0 - w) * fD_smooth + w * fD_rough;
        fF = 0.25 * std::max(1e-6, fD);
    }
#
    return 2.0 * fF * density * velocity * velocity / diameter;
}
#
double FrictionModel::applyDragReduction(double base_gradient, double drag_reduction_factor) const {
    // drag_reduction_factor in [0,1]: 0=no DR, 1=fully eliminated friction
    const double dr = clamp01(drag_reduction_factor);
    return base_gradient * (1.0 - dr);
}
#
// ============================================================================
// PerforationFriction
// ============================================================================
#
PerforationFriction::PerforationFriction() : erosion_coeff_(1e-11), min_cd_(0.6) {}
#
double PerforationFriction::calculatePressureDrop(double flow_rate, int num_perforations,
                                                  double perf_diameter,
                                                  double discharge_coefficient,
                                                  double fluid_density) const {
    if (num_perforations <= 0 || perf_diameter <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (fluid_density <= 0.0) return 0.0;
    const double Cd = std::max(min_cd_, discharge_coefficient);
    const double area = num_perforations * PI * perf_diameter * perf_diameter / 4.0;
    if (area <= 0.0) return std::numeric_limits<double>::infinity();
    // Orifice equation: ΔP = ρ Q^2 / (2 Cd^2 A^2)
    return fluid_density * flow_rate * flow_rate / (2.0 * Cd * Cd * area * area);
}
#
double PerforationFriction::calculateErosion(double flow_rate, int num_perforations,
                                            double current_diameter,
                                            double proppant_concentration,
                                            double proppant_hardness,
                                            double dt) const {
    if (flow_rate <= 0.0 || dt <= 0.0 || num_perforations <= 0 || current_diameter <= 0.0) {
        return current_diameter;
    }
    const double area = num_perforations * PI * current_diameter * current_diameter / 4.0;
    if (area <= 0.0) return current_diameter;
#
    const double v = flow_rate / area;
    const double hardness = std::max(0.1, proppant_hardness);
    const double erosion_rate = erosion_coeff_ * v * v * std::max(0.0, proppant_concentration) / hardness;
#
    return std::max(1e-6, current_diameter + erosion_rate * dt);
}
#
double PerforationFriction::calculateDischargeCoefficient(double reynolds, double l_over_d) const {
    // A simple smooth transition to typical Cd values.
    // For sharp-edged orifices Cd ~ 0.6-0.85 depending on Re and L/D.
    (void)l_over_d;
    if (!std::isfinite(reynolds) || reynolds <= 0.0) return min_cd_;
    const double cd = 0.62 + 0.2 * (1.0 - std::exp(-reynolds / 2e4));
    return std::max(min_cd_, std::min(0.95, cd));
}
#
// ============================================================================
// NearWellboreFriction
// ============================================================================
#
NearWellboreFriction::NearWellboreFriction() : tortuosity_coeff_(0.1) {}
#
double NearWellboreFriction::getTortuosityFactor(double perforation_phasing,
                                                 double stress_anisotropy) const {
    // Factor >= 1; increases with misalignment (phasing away from frac plane)
    // and with stress anisotropy (harder to initiate aligned flow paths).
    const double phase_rad = perforation_phasing * PI / 180.0;
    const double misalign = std::abs(std::sin(phase_rad));
    const double stress_scale = 1.0 + std::max(0.0, stress_anisotropy) / 25e6; // normalize to 25 MPa
    return 1.0 + tortuosity_coeff_ * misalign * stress_scale;
}
#
double NearWellboreFriction::calculateTortuosityPressureDrop(double flow_rate,
                                                             int num_clusters,
                                                             double fluid_viscosity,
                                                             double perforation_phasing,
                                                             double stress_anisotropy) const {
    if (num_clusters <= 0) return 0.0;
    if (flow_rate <= 0.0 || fluid_viscosity <= 0.0) return 0.0;
    const double tf = getTortuosityFactor(perforation_phasing, stress_anisotropy);
    // Empirical: dp scales with viscosity and per-cluster rate.
    const double q_cluster = flow_rate / static_cast<double>(num_clusters);
    return tf * fluid_viscosity * q_cluster * 1e8; // calibrated scale (Pa), user-tunable via tortuosity_coeff_
}
#
// ============================================================================
// WellboreTemperature
// ============================================================================
#
WellboreTemperature::WellboreTemperature()
    : surface_temp_(288.15),
      geothermal_gradient_(0.025),
      use_transient_(false),
      formation_diffusivity_(1e-6),
      injection_time_(0.0) {}
#
void WellboreTemperature::initialize(const std::vector<WellboreSegment>& segments,
                                     double surface_temp,
                                     double geothermal_gradient) {
    segments_ = segments;
    surface_temp_ = surface_temp;
    geothermal_gradient_ = geothermal_gradient;
#
    formation_temp_.resize(segments_.size(), surface_temp_);
    current_temp_.resize(segments_.size(), surface_temp_);
#
    for (size_t i = 0; i < segments_.size(); ++i) {
        const double tvd = 0.5 * (segments_[i].start_tvd + segments_[i].end_tvd);
        formation_temp_[i] = surface_temp_ + geothermal_gradient_ * tvd;
        current_temp_[i] = formation_temp_[i];
    }
#
    injection_time_ = 0.0;
}
#
double WellboreTemperature::rameyTimeFunction(double time, double radius) const {
    // Simple, bounded approximation. This is not a full Ramey implementation,
    // but provides a stable transient knob.
    if (time <= 0.0) return 0.0;
    const double alpha = std::max(1e-12, formation_diffusivity_);
    const double x = radius * radius / (4.0 * alpha * time);
    return std::exp(-x);
}
#
double WellboreTemperature::calculateHeatLoss(const WellboreSegment& segment,
                                              double fluid_temp,
                                              double flow_rate) const {
    // Very simple UA model (W): UA ~ 2π k_eff L / ln(r_out/r_in)
    (void)flow_rate;
    const double L = std::max(0.0, segment.end_md - segment.start_md);
    const double r_in = std::max(1e-6, segment.inner_diameter * 0.5);
    const double r_out = std::max(r_in * 1.1, segment.outer_diameter * 0.5);
    const double k_eff = std::max(0.1, segment.formation_conductivity);
    const double UA = 2.0 * PI * k_eff * L / std::log(r_out / r_in);
#
    // Driving ΔT to formation temperature at segment midpoint.
    const double tvd = 0.5 * (segment.start_tvd + segment.end_tvd);
    const double T_form = surface_temp_ + geothermal_gradient_ * tvd;
    return UA * (fluid_temp - T_form);
}
#
std::vector<double> WellboreTemperature::calculateInjectionProfile(double flow_rate,
                                                                   double injection_temp,
                                                                   double fluid_density,
                                                                   double fluid_cp,
                                                                   double time) {
    injection_time_ = time;
    if (segments_.empty()) return {};
    if (flow_rate <= 0.0 || fluid_density <= 0.0 || fluid_cp <= 0.0) {
        current_temp_ = formation_temp_;
        return current_temp_;
    }
#
    const double m_dot = fluid_density * flow_rate;
    double T = injection_temp;
    for (size_t i = 0; i < segments_.size(); ++i) {
        const auto& seg = segments_[i];
        const double L = std::max(0.0, seg.end_md - seg.start_md);
        if (L <= 0.0) {
            current_temp_[i] = T;
            continue;
        }
#
        const double tvd = 0.5 * (seg.start_tvd + seg.end_tvd);
        const double T_form = surface_temp_ + geothermal_gradient_ * tvd;
#
        // Exponential relaxation toward formation temperature.
        // Add an optional transient limiter.
        double k_relax = 1.0;
        if (use_transient_) {
            const double r = std::max(0.05, seg.outer_diameter * 0.5);
            k_relax = rameyTimeFunction(std::max(1e-6, time), r);
        }
        const double tau = std::max(1e-12, (m_dot * fluid_cp) / (2.0 * PI * std::max(0.1, seg.formation_conductivity)));
        const double a = std::exp(-k_relax * L / tau);
        T = T_form + (T - T_form) * a;
        current_temp_[i] = T;
    }
#
    return current_temp_;
}
#
double WellboreTemperature::getTemperatureAtDepth(double md) const {
    if (segments_.empty()) return surface_temp_;
    for (size_t i = 0; i < segments_.size(); ++i) {
        if (md >= segments_[i].start_md && md <= segments_[i].end_md) {
            return current_temp_.empty() ? formation_temp_[i] : current_temp_[i];
        }
    }
    return current_temp_.empty() ? formation_temp_.back() : current_temp_.back();
}
#
double WellboreTemperature::getBottomholeTemperature() const {
    if (current_temp_.empty()) return surface_temp_;
    return current_temp_.back();
}
#
// ============================================================================
// WellboreHydraulicsCalculator
// ============================================================================
#
WellboreHydraulicsCalculator::WellboreHydraulicsCalculator()
    : drag_reduction_(0.0),
      use_temperature_model_(false) {}
#
static double deg2rad(double deg) { return deg * PI / 180.0; }
#
void WellboreHydraulicsCalculator::buildFromSurvey(const std::vector<double>& md,
                                                   const std::vector<double>& inc,
                                                   const std::vector<double>& azi,
                                                   const std::vector<double>& diameters,
                                                   double roughness) {
    if (md.size() < 2 || inc.size() != md.size() || azi.size() != md.size() ||
        diameters.size() != md.size()) {
        throw std::invalid_argument("buildFromSurvey: md/inc/azi/diameters must have same size >= 2");
    }
#
    segments_.clear();
#
    double tvd = 0.0;
    for (size_t i = 1; i < md.size(); ++i) {
        WellboreSegment seg;
        seg.start_md = md[i - 1];
        seg.end_md = md[i];
        seg.inclination = 0.5 * (inc[i - 1] + inc[i]);
        seg.azimuth = 0.5 * (azi[i - 1] + azi[i]);
        seg.inner_diameter = 0.5 * (diameters[i - 1] + diameters[i]);
        seg.roughness = roughness;
#
        const double dmd = seg.end_md - seg.start_md;
        const double inc_rad = deg2rad(seg.inclination);
        const double dtvd = dmd * std::cos(inc_rad);
#
        seg.start_tvd = tvd;
        seg.end_tvd = tvd + dtvd;
        tvd = seg.end_tvd;
#
        seg.calculateDerivedProperties();
        segments_.push_back(seg);
    }
#
    if (use_temperature_model_ && temp_model_) {
        // Keep previously configured thermal settings.
        // (Call enableTemperatureModel again if you want different settings.)
    }
}
#
void WellboreHydraulicsCalculator::setCompletion(const WellboreCompletion& completion) {
    completion_ = completion;
}
#
void WellboreHydraulicsCalculator::addPerforations(const std::vector<double>& cluster_mds,
                                                   const std::vector<int>& num_perfs,
                                                   const std::vector<double>& perf_diameters) {
    if (cluster_mds.size() != num_perfs.size() || cluster_mds.size() != perf_diameters.size()) {
        throw std::invalid_argument("addPerforations: arrays must have same size");
    }
#
    perforations_.clear();
    for (size_t i = 0; i < cluster_mds.size(); ++i) {
        PerfCluster pc;
        pc.md = cluster_mds[i];
        pc.num_perfs = num_perfs[i];
        pc.diameter = perf_diameters[i];
        pc.discharge_coeff = 0.85;
        perforations_.push_back(pc);
    }
}
#
double WellboreHydraulicsCalculator::getTotalMeasuredDepth() const {
    if (segments_.empty()) return 0.0;
    return segments_.back().end_md;
}
#
double WellboreHydraulicsCalculator::getTotalVerticalDepth() const {
    if (segments_.empty()) return 0.0;
    return segments_.back().end_tvd;
}
#
int WellboreHydraulicsCalculator::getSegmentAtMD(double md) const {
    for (size_t i = 0; i < segments_.size(); ++i) {
        if (md >= segments_[i].start_md && md <= segments_[i].end_md) {
            return static_cast<int>(i);
        }
    }
    return segments_.empty() ? -1 : static_cast<int>(segments_.size() - 1);
}
#
double WellboreHydraulicsCalculator::getTVDAtMD(double md) const {
    if (segments_.empty()) return 0.0;
    const int idx = getSegmentAtMD(md);
    if (idx < 0) return 0.0;
#
    const auto& seg = segments_[static_cast<size_t>(idx)];
    const double L = std::max(1e-12, seg.end_md - seg.start_md);
    const double t = std::max(0.0, std::min(1.0, (md - seg.start_md) / L));
    return seg.start_tvd + t * (seg.end_tvd - seg.start_tvd);
}
#
double WellboreHydraulicsCalculator::calculateHydrostatic(double md1, double md2,
                                                          double fluid_density) const {
    if (fluid_density <= 0.0) return 0.0;
    const double tvd1 = getTVDAtMD(md1);
    const double tvd2 = getTVDAtMD(md2);
    const double dTVD = tvd2 - tvd1;
    return fluid_density * G * dTVD;
}
#
double WellboreHydraulicsCalculator::calculateFriction(double md1, double md2,
                                                       double flow_rate,
                                                       double fluid_density,
                                                       double fluid_viscosity) const {
    if (segments_.empty()) return 0.0;
    if (flow_rate <= 0.0 || fluid_density <= 0.0 || fluid_viscosity <= 0.0) return 0.0;
#
    const double md_lo = std::min(md1, md2);
    const double md_hi = std::max(md1, md2);
    double dp = 0.0;
#
    // Integrate segment-by-segment
    for (const auto& seg : segments_) {
        const double a = std::max(md_lo, seg.start_md);
        const double b = std::min(md_hi, seg.end_md);
        if (b <= a) continue;
#
        const double L = b - a;
        if (seg.flow_area <= 0.0 || seg.inner_diameter <= 0.0) continue;
#
        const double v = flow_rate / seg.flow_area;
        double grad = friction_model_.newtonianFrictionGradient(v, fluid_density, fluid_viscosity,
                                                                seg.inner_diameter, seg.roughness);
        grad = friction_model_.applyDragReduction(grad, drag_reduction_);
        dp += grad * L;
    }
#
    return dp;
}
#
std::vector<double> WellboreHydraulicsCalculator::calculatePressureProfile(double surface_pressure,
                                                                           double flow_rate,
                                                                           double fluid_density,
                                                                           double fluid_viscosity,
                                                                           double k_prime,
                                                                           double n_prime) {
    pressure_breakdown_.clear();
#
    if (segments_.empty()) {
        pressure_breakdown_["hydrostatic"] = 0.0;
        pressure_breakdown_["friction"] = 0.0;
        pressure_breakdown_["perforations"] = 0.0;
        pressure_breakdown_["near_wellbore"] = 0.0;
        return {};
    }
#
    std::vector<double> pressures;
    pressures.reserve(segments_.size());
#
    // Use either Newtonian or power-law friction.
    const bool use_powerlaw = (k_prime > 0.0 && n_prime > 0.0 && std::abs(n_prime - 1.0) > 1e-12);
#
    double p = surface_pressure;
    double dp_h = 0.0;
    double dp_f = 0.0;
#
    // Segment pressures at end of each segment
    for (const auto& seg : segments_) {
        const double L = std::max(0.0, seg.end_md - seg.start_md);
        const double dTVD = seg.end_tvd - seg.start_tvd;
#
        const double dph = fluid_density > 0.0 ? fluid_density * G * dTVD : 0.0;
        dp_h += dph;
#
        double dpf = 0.0;
        if (flow_rate > 0.0 && seg.flow_area > 0.0 && seg.inner_diameter > 0.0) {
            const double v = flow_rate / seg.flow_area;
            double grad = 0.0;
            if (use_powerlaw) {
                grad = friction_model_.powerLawFrictionGradient(v, fluid_density, k_prime, n_prime,
                                                               seg.inner_diameter, seg.roughness);
            } else {
                grad = friction_model_.newtonianFrictionGradient(v, fluid_density, fluid_viscosity,
                                                                 seg.inner_diameter, seg.roughness);
            }
            grad = friction_model_.applyDragReduction(grad, drag_reduction_);
            dpf = grad * L;
        }
        dp_f += dpf;
#
        p += dph + dpf;
        pressures.push_back(p);
    }
#
    // Perforation friction: assume total rate through all clusters
    double dp_perf = 0.0;
    for (auto& pc : perforations_) {
        // Update discharge coefficient (optional)
        // Assume perf length/diam ~ 5 if unknown
        const double Re_perf = (fluid_density > 0.0 && fluid_viscosity > 0.0 && pc.diameter > 0.0) ?
            (fluid_density * (flow_rate / std::max(1e-12, PI * pc.diameter * pc.diameter / 4.0)) * pc.diameter / fluid_viscosity) :
            0.0;
        pc.discharge_coeff = perf_friction_.calculateDischargeCoefficient(Re_perf, 5.0);
        dp_perf += perf_friction_.calculatePressureDrop(flow_rate, pc.num_perfs, pc.diameter,
                                                        pc.discharge_coeff, fluid_density);
    }
#
    // Near-wellbore tortuosity friction: use a conservative default phasing=60deg, stress anisotropy=0
    const int nclusters = std::max(1, static_cast<int>(perforations_.size()));
    const double dp_nwb = nwb_friction_.calculateTortuosityPressureDrop(flow_rate, nclusters,
                                                                        fluid_viscosity, 60.0, 0.0);
#
    // Apply end effects to the last pressure (bottomhole)
    if (!pressures.empty()) {
        pressures.back() += dp_perf + dp_nwb;
    }
#
    pressure_breakdown_["hydrostatic"] = dp_h;
    pressure_breakdown_["friction"] = dp_f;
    pressure_breakdown_["perforations"] = dp_perf;
    pressure_breakdown_["near_wellbore"] = dp_nwb;
#
    return pressures;
}
#
double WellboreHydraulicsCalculator::calculateSurfacePressure(double bhp,
                                                              double flow_rate,
                                                              double fluid_density,
                                                              double fluid_viscosity) {
    // Invert using a single forward evaluation because the model is monotone in surface pressure:
    // bhp = surface_pressure + dp_total  => surface_pressure = bhp - dp_total
    calculatePressureProfile(0.0, flow_rate, fluid_density, fluid_viscosity);
    const double dp_total = pressure_breakdown_["hydrostatic"] +
                            pressure_breakdown_["friction"] +
                            pressure_breakdown_["perforations"] +
                            pressure_breakdown_["near_wellbore"];
    return bhp - dp_total;
}
#
double WellboreHydraulicsCalculator::calculateBHP(double surface_pressure,
                                                  double flow_rate,
                                                  double fluid_density,
                                                  double fluid_viscosity) {
    const auto profile = calculatePressureProfile(surface_pressure, flow_rate,
                                                  fluid_density, fluid_viscosity);
    if (profile.empty()) return surface_pressure;
    return profile.back();
}
#
std::map<std::string, double> WellboreHydraulicsCalculator::getPressureBreakdown() const {
    return pressure_breakdown_;
}
#
void WellboreHydraulicsCalculator::updatePerforations(double flow_rate,
                                                      double proppant_concentration,
                                                      double dt) {
    // Erode each perf cluster independently (very simplified).
    for (auto& pc : perforations_) {
        pc.diameter = perf_friction_.calculateErosion(flow_rate, pc.num_perfs, pc.diameter,
                                                      proppant_concentration, /*hardness*/1.0, dt);
    }
}
#
void WellboreHydraulicsCalculator::enableTemperatureModel(double surface_temp, double geothermal_grad) {
    if (!temp_model_) temp_model_ = std::make_unique<WellboreTemperature>();
    use_temperature_model_ = true;
    temp_model_->initialize(segments_, surface_temp, geothermal_grad);
}
#
std::vector<double> WellboreHydraulicsCalculator::getTemperatureProfile() const {
    if (!use_temperature_model_ || !temp_model_) return {};
    std::vector<double> temps;
    temps.reserve(segments_.size());
    for (const auto& seg : segments_) {
        temps.push_back(temp_model_->getTemperatureAtDepth(0.5 * (seg.start_md + seg.end_md)));
    }
    return temps;
}
#
} // namespace FSRM

