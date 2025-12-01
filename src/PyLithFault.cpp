/**
 * @file PyLithFault.cpp
 * @brief Implementation of PyLith-style fault representation
 */

#include "PyLithFault.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace FSRM {

// =============================================================================
// Utility Functions
// =============================================================================

static double parseDouble(const std::map<std::string, std::string>& config,
                         const std::string& key, double default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        try {
            return std::stod(it->second);
        } catch (...) {
            return default_val;
        }
    }
    return default_val;
}

static std::string parseString(const std::map<std::string, std::string>& config,
                               const std::string& key, const std::string& default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        return it->second;
    }
    return default_val;
}

static bool parseBool(const std::map<std::string, std::string>& config,
                     const std::string& key, bool default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        std::string s = it->second;
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        return (s == "true" || s == "1" || s == "yes" || s == "on");
    }
    return default_val;
}

// Vector operations
static double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static std::array<double, 3> cross3(const std::array<double, 3>& a, 
                                    const std::array<double, 3>& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

static double norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static void normalize3(std::array<double, 3>& v) {
    double n = norm3(v);
    if (n > 1e-15) {
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
    }
}

// =============================================================================
// SlipTimeFn Implementations
// =============================================================================

void SlipTimeFnRamp::configure(const std::map<std::string, std::string>& config) {
    rise_time = parseDouble(config, "rise_time", 1.0);
}

void SlipTimeFnBrune::configure(const std::map<std::string, std::string>& config) {
    tau_c = parseDouble(config, "characteristic_time", 1.0);
    // Alternative: compute from rise time
    double tr = parseDouble(config, "rise_time", -1.0);
    if (tr > 0) {
        // For Brune, tau_c ≈ rise_time / 2.5
        tau_c = tr / 2.5;
    }
}

// Liu et al. slip time function
void SlipTimeFnLiu::computeNormalization() {
    // Normalization so that final slip is amplitude
    // For Liu function, cnorm ≈ 1.0 (already normalized)
    cnorm = 1.0;
}

double SlipTimeFnLiu::computeSlip(double t, double slip_time, double amplitude) const {
    double tau = t - slip_time;
    if (tau <= 0.0) return 0.0;
    if (tau >= rise_time) return amplitude;
    
    // Liu et al. (2006) regularized function
    double x = tau / rise_time;
    double slip_norm;
    
    if (x < 0.5) {
        // Cubic polynomial for first half
        slip_norm = 4.0 * x * x * x;
    } else {
        // Symmetric completion
        double y = 1.0 - x;
        slip_norm = 1.0 - 4.0 * y * y * y;
    }
    
    return amplitude * slip_norm * cnorm;
}

double SlipTimeFnLiu::computeSlipRate(double t, double slip_time, double amplitude) const {
    double tau = t - slip_time;
    if (tau <= 0.0 || tau >= rise_time) return 0.0;
    
    double x = tau / rise_time;
    double rate_norm;
    
    if (x < 0.5) {
        rate_norm = 12.0 * x * x;
    } else {
        double y = 1.0 - x;
        rate_norm = 12.0 * y * y;
    }
    
    return amplitude * rate_norm * cnorm / rise_time;
}

void SlipTimeFnLiu::configure(const std::map<std::string, std::string>& config) {
    rise_time = parseDouble(config, "rise_time", 1.0);
    computeNormalization();
}

void SlipTimeFnConstantRate::configure(const std::map<std::string, std::string>& config) {
    slip_rate = parseDouble(config, "slip_rate", 1e-9);
}

// Time history implementation
double SlipTimeFnTimeHistory::interpolate(double t, double slip_time) const {
    double tau = t - slip_time;
    if (tau <= 0.0 || time_points.empty()) return 0.0;
    
    // Find bracketing indices
    if (tau <= time_points.front()) return slip_values.front();
    if (tau >= time_points.back()) return slip_values.back();
    
    for (size_t i = 1; i < time_points.size(); ++i) {
        if (tau <= time_points[i]) {
            double t0 = time_points[i-1];
            double t1 = time_points[i];
            double v0 = slip_values[i-1];
            double v1 = slip_values[i];
            double frac = (tau - t0) / (t1 - t0);
            return v0 + frac * (v1 - v0);
        }
    }
    return slip_values.back();
}

double SlipTimeFnTimeHistory::computeSlip(double t, double slip_time, double amplitude) const {
    return amplitude * interpolate(t, slip_time);
}

double SlipTimeFnTimeHistory::computeSlipRate(double t, double slip_time, double amplitude) const {
    double dt = 1e-6;
    double s1 = computeSlip(t + dt, slip_time, amplitude);
    double s0 = computeSlip(t, slip_time, amplitude);
    return (s1 - s0) / dt;
}

void SlipTimeFnTimeHistory::configure(const std::map<std::string, std::string>& config) {
    // Load from file if specified
    std::string filename = parseString(config, "time_history_file", "");
    if (!filename.empty()) {
        std::ifstream file(filename);
        if (file.is_open()) {
            time_points.clear();
            slip_values.clear();
            double t, s;
            while (file >> t >> s) {
                time_points.push_back(t);
                slip_values.push_back(s);
            }
        }
    }
}

// =============================================================================
// Friction Model Implementations
// =============================================================================

void StaticFriction::configure(const std::map<std::string, std::string>& config) {
    coef = parseDouble(config, "friction_coefficient", 0.6);
}

double SlipWeakeningFriction::computeFriction(const DynamicFaultState& state) const {
    // state_variable holds cumulative slip
    double slip = state.state_variable;
    
    if (slip <= 0.0) {
        return static_coef;
    } else if (slip >= d0) {
        return dynamic_coef;
    } else {
        // Linear interpolation
        return static_coef - (static_coef - dynamic_coef) * slip / d0;
    }
}

double SlipWeakeningFriction::computeStateRate(const DynamicFaultState& state) const {
    // State is cumulative slip, rate is slip rate
    return state.slip_rate;
}

void SlipWeakeningFriction::updateState(DynamicFaultState& state, double dt) const {
    // Update cumulative slip
    state.state_variable += state.slip_rate * dt;
    state.friction = computeFriction(state);
}

void SlipWeakeningFriction::configure(const std::map<std::string, std::string>& config) {
    static_coef = parseDouble(config, "static_friction", 0.6);
    dynamic_coef = parseDouble(config, "dynamic_friction", 0.4);
    d0 = parseDouble(config, "critical_slip_distance", 0.4);
    cohesion = parseDouble(config, "cohesion", 0.0);
}

// Rate-and-state friction (aging law)
double RateStateFrictionAging::computeFriction(const DynamicFaultState& state) const {
    double V = std::max(state.slip_rate, linear_threshold);
    double theta = std::max(state.state_variable, 1e-20);
    
    // f = f0 + a*ln(V/V0) + b*ln(θ*V0/Dc)
    double f = f0 + a * std::log(V / v0) + b * std::log(theta * v0 / dc);
    
    return std::max(0.0, f);
}

double RateStateFrictionAging::computeStateRate(const DynamicFaultState& state) const {
    double V = std::max(std::abs(state.slip_rate), linear_threshold);
    double theta = std::max(state.state_variable, 1e-20);
    
    // Aging law: dθ/dt = 1 - V*θ/Dc
    return 1.0 - V * theta / dc;
}

void RateStateFrictionAging::updateState(DynamicFaultState& state, double dt) const {
    double dtheta = computeStateRate(state) * dt;
    state.state_variable = std::max(1e-20, state.state_variable + dtheta);
    state.friction = computeFriction(state);
}

void RateStateFrictionAging::configure(const std::map<std::string, std::string>& config) {
    a = parseDouble(config, "a", 0.01);
    b = parseDouble(config, "b", 0.015);
    dc = parseDouble(config, "dc", 0.01);
    f0 = parseDouble(config, "f0", 0.6);
    v0 = parseDouble(config, "v0", 1e-6);
    linear_threshold = parseDouble(config, "linear_threshold", 1e-12);
}

double RateStateFrictionAging::getCriticalStiffness(double sigma_n) const {
    // k_c = (b - a) * σ_n / D_c
    return (b - a) * sigma_n / dc;
}

double RateStateFrictionAging::getNucleationLength(double G, double sigma_n) const {
    // L_c = G * D_c / ((b - a) * σ_n)
    if (b <= a) return std::numeric_limits<double>::max();  // Stable
    return G * dc / ((b - a) * sigma_n);
}

// Rate-and-state friction (slip law)
double RateStateFrictionSlip::computeFriction(const DynamicFaultState& state) const {
    double V = std::max(state.slip_rate, linear_threshold);
    double theta = std::max(state.state_variable, 1e-20);
    
    double f = f0 + a * std::log(V / v0) + b * std::log(theta * v0 / dc);
    return std::max(0.0, f);
}

double RateStateFrictionSlip::computeStateRate(const DynamicFaultState& state) const {
    double V = std::max(std::abs(state.slip_rate), linear_threshold);
    double theta = std::max(state.state_variable, 1e-20);
    
    // Slip law: dθ/dt = -V*θ/Dc * ln(V*θ/Dc)
    double x = V * theta / dc;
    return -x * std::log(x) * dc / theta;
}

void RateStateFrictionSlip::updateState(DynamicFaultState& state, double dt) const {
    double dtheta = computeStateRate(state) * dt;
    state.state_variable = std::max(1e-20, state.state_variable + dtheta);
    state.friction = computeFriction(state);
}

void RateStateFrictionSlip::configure(const std::map<std::string, std::string>& config) {
    a = parseDouble(config, "a", 0.01);
    b = parseDouble(config, "b", 0.015);
    dc = parseDouble(config, "dc", 0.01);
    f0 = parseDouble(config, "f0", 0.6);
    v0 = parseDouble(config, "v0", 1e-6);
    linear_threshold = parseDouble(config, "linear_threshold", 1e-12);
}

// =============================================================================
// FaultCohesive Base Class Implementation
// =============================================================================

FaultCohesive::FaultCohesive(CohesiveFaultType type)
    : name("fault"), fault_type(type),
      shear_modulus(30e9), bulk_modulus(50e9) {}

void FaultCohesive::configure(const std::map<std::string, std::string>& config) {
    name = parseString(config, "name", "fault");
    shear_modulus = parseDouble(config, "shear_modulus", 30e9);
    bulk_modulus = parseDouble(config, "bulk_modulus", 50e9);
    
    // Configure output fields
    std::string output_str = parseString(config, "output_fields", "slip,traction");
    // Parse comma-separated list (simplified)
    output_fields.clear();
    output_fields.push_back(FaultOutputField::SLIP);
    output_fields.push_back(FaultOutputField::TRACTION);
}

void FaultCohesive::setFaultVertices(const std::vector<FaultVertex>& verts) {
    vertices = verts;
    slip_field.resize(vertices.size());
    traction_field.resize(vertices.size());
}

void FaultCohesive::setCohesiveCells(const std::vector<CohesiveCell>& cells_in) {
    cells = cells_in;
}

void FaultCohesive::createCohesiveCells(const std::vector<double>& node_coords,
                                        const std::vector<int>& fault_node_indices,
                                        double strike, double dip) {
    vertices.clear();
    cells.clear();
    
    // Compute fault orientation vectors
    double cos_s = std::cos(strike);
    double sin_s = std::sin(strike);
    double cos_d = std::cos(dip);
    double sin_d = std::sin(dip);
    
    std::array<double, 3> strike_dir = {sin_s, cos_s, 0.0};
    std::array<double, 3> dip_dir = {-cos_s * cos_d, sin_s * cos_d, -sin_d};
    std::array<double, 3> normal = {cos_s * sin_d, -sin_s * sin_d, cos_d};
    
    // Create vertices for each fault node
    for (size_t i = 0; i < fault_node_indices.size(); ++i) {
        int node_idx = fault_node_indices[i];
        
        FaultVertex v;
        v.vertex_id = static_cast<int>(i);
        v.vertex_negative = node_idx;
        v.vertex_positive = -1;  // Will be assigned when mesh is split
        v.vertex_lagrange = -1;
        
        v.coords[0] = node_coords[3 * node_idx];
        v.coords[1] = node_coords[3 * node_idx + 1];
        v.coords[2] = node_coords[3 * node_idx + 2];
        
        v.normal = normal;
        v.along_strike = strike_dir;
        v.up_dip = dip_dir;
        v.area = 1.0;  // Will be computed from cells
        
        vertices.push_back(v);
    }
    
    slip_field.resize(vertices.size());
    traction_field.resize(vertices.size());
}

void FaultCohesive::computeFaultOrientation(const FaultVertex& v,
                                           std::array<double, 3>& normal,
                                           std::array<double, 3>& strike_dir,
                                           std::array<double, 3>& dip_dir) const {
    normal = v.normal;
    strike_dir = v.along_strike;
    dip_dir = v.up_dip;
}

void FaultCohesive::transformToFaultCoords(const double* global_vector,
                                          const FaultVertex& v,
                                          double& ll_component,
                                          double& ud_component,
                                          double& normal_component) const {
    // Project global vector onto fault coordinate system
    std::array<double, 3> vec = {global_vector[0], global_vector[1], global_vector[2]};
    
    ll_component = dot3(vec, v.along_strike);
    ud_component = dot3(vec, v.up_dip);
    normal_component = dot3(vec, v.normal);
}

void FaultCohesive::transformToGlobalCoords(double ll, double ud, double normal,
                                           const FaultVertex& v,
                                           double* global_vector) const {
    // Transform from fault coords to global
    for (int i = 0; i < 3; ++i) {
        global_vector[i] = ll * v.along_strike[i] + 
                          ud * v.up_dip[i] + 
                          normal * v.normal[i];
    }
}

void FaultCohesive::writeOutput(const std::string& filename, double time) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open fault output file: " << filename << std::endl;
        return;
    }
    
    file << "# Fault output: " << name << "\n";
    file << "# Time: " << time << " s\n";
    file << "# Vertices: " << vertices.size() << "\n";
    file << "# Format: x, y, z, slip_ll, slip_ud, slip_n, traction_ll, traction_ud, traction_n\n";
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& v = vertices[i];
        file << std::scientific << std::setprecision(6)
             << v.coords[0] << " " << v.coords[1] << " " << v.coords[2] << " "
             << slip_field.slip_left_lateral[i] << " "
             << slip_field.slip_reverse[i] << " "
             << slip_field.slip_opening[i] << " "
             << traction_field.traction_shear_ll[i] << " "
             << traction_field.traction_shear_ud[i] << " "
             << traction_field.traction_normal[i] << "\n";
    }
    
    file.close();
}

// =============================================================================
// FaultCohesiveKin Implementation
// =============================================================================

FaultCohesiveKin::FaultCohesiveKin()
    : FaultCohesive(CohesiveFaultType::KINEMATIC), current_time(0.0) {
    slip_time_fn = std::make_unique<SlipTimeFnBrune>();
}

void FaultCohesiveKin::configure(const std::map<std::string, std::string>& config) {
    FaultCohesive::configure(config);
    
    // Set slip time function type
    std::string fn_type = parseString(config, "slip_time_fn", "brune");
    setSlipTimeFnType(parseSlipTimeFnType(fn_type));
    slip_time_fn->configure(config);
    
    // Check for uniform slip parameters
    double slip_ll = parseDouble(config, "final_slip_ll", 0.0);
    double slip_rev = parseDouble(config, "final_slip_reverse", 0.0);
    double slip_open = parseDouble(config, "final_slip_opening", 0.0);
    double slip_time = parseDouble(config, "slip_time", 0.0);
    double rise_time = parseDouble(config, "rise_time", 1.0);
    
    if (std::abs(slip_ll) > 0 || std::abs(slip_rev) > 0 || std::abs(slip_open) > 0) {
        // Will be applied in initialize() once vertices are set
    }
}

void FaultCohesiveKin::setSlipTimeFn(std::unique_ptr<SlipTimeFn> fn) {
    slip_time_fn = std::move(fn);
}

void FaultCohesiveKin::setSlipTimeFnType(SlipTimeFnType type) {
    slip_time_fn = createSlipTimeFn(type);
}

void FaultCohesiveKin::setSlipParams(const std::vector<KinematicSlipParams>& params) {
    kin_params = params;
}

void FaultCohesiveKin::setUniformSlip(double slip_ll, double slip_reverse, double slip_opening,
                                     double slip_time, double rise_time) {
    kin_params.resize(vertices.size());
    for (auto& p : kin_params) {
        p.final_slip_ll = slip_ll;
        p.final_slip_reverse = slip_reverse;
        p.final_slip_opening = slip_opening;
        p.slip_time = slip_time;
        p.rise_time = rise_time;
    }
}

void FaultCohesiveKin::setCircularRupture(double hypo_x, double hypo_y, double hypo_z,
                                         double rupture_velocity, double final_slip,
                                         double rise_time, double rake_angle) {
    kin_params.resize(vertices.size());
    
    double cos_rake = std::cos(rake_angle);
    double sin_rake = std::sin(rake_angle);
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& v = vertices[i];
        
        // Distance from hypocenter
        double dx = v.coords[0] - hypo_x;
        double dy = v.coords[1] - hypo_y;
        double dz = v.coords[2] - hypo_z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        kin_params[i].slip_time = dist / rupture_velocity;
        kin_params[i].rise_time = rise_time;
        kin_params[i].final_slip_ll = final_slip * cos_rake;
        kin_params[i].final_slip_reverse = final_slip * sin_rake;
        kin_params[i].final_slip_opening = 0.0;
    }
}

void FaultCohesiveKin::initialize() {
    if (kin_params.empty() && !vertices.empty()) {
        // Initialize with zero slip
        kin_params.resize(vertices.size());
    }
}

void FaultCohesiveKin::prestep(double t, double dt) {
    (void)dt;
    current_time = t;
}

void FaultCohesiveKin::computeResidual(const double* solution, double* residual) {
    // For kinematic faults, the residual enforces:
    // [[u]] = prescribed_slip
    // where [[u]] is the displacement jump across the fault
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& v = vertices[i];
        if (v.vertex_negative < 0 || v.vertex_positive < 0) continue;
        
        // Get prescribed slip at current time
        double slip_ll, slip_rev, slip_open;
        getSlipAtTime(current_time, i, slip_ll, slip_rev, slip_open);
        
        // Get actual displacement jump from solution
        double u_neg[3], u_pos[3];
        for (int j = 0; j < 3; ++j) {
            u_neg[j] = solution[3 * v.vertex_negative + j];
            u_pos[j] = solution[3 * v.vertex_positive + j];
        }
        
        // Displacement jump in global coords
        double du[3] = {u_pos[0] - u_neg[0], 
                       u_pos[1] - u_neg[1], 
                       u_pos[2] - u_neg[2]};
        
        // Transform to fault coords
        double actual_ll, actual_ud, actual_n;
        transformToFaultCoords(du, v, actual_ll, actual_ud, actual_n);
        
        // Residual is difference (constraint violation)
        double r_ll = actual_ll - slip_ll;
        double r_ud = actual_ud - slip_rev;
        double r_n = actual_n - slip_open;
        
        // Transform back to global for residual
        double r_global[3];
        transformToGlobalCoords(r_ll, r_ud, r_n, v, r_global);
        
        // Add to residual at Lagrange multiplier DOFs
        if (v.vertex_lagrange >= 0) {
            for (int j = 0; j < 3; ++j) {
                residual[3 * v.vertex_lagrange + j] += r_global[j] * v.area;
            }
        }
        
        // Store current slip in field
        slip_field.slip_left_lateral[i] = slip_ll;
        slip_field.slip_reverse[i] = slip_rev;
        slip_field.slip_opening[i] = slip_open;
    }
}

void FaultCohesiveKin::computeJacobian(const double* solution, double* jacobian) {
    (void)solution;
    
    // Jacobian entries for kinematic constraint:
    // ∂R_λ/∂u+ = T^T (identity in fault coords)
    // ∂R_λ/∂u- = -T^T
    // where T is the transformation matrix
    
    // Also need entries for momentum equation:
    // ∂R_u/∂λ = ±T (traction from Lagrange multiplier)
    
    // Implementation depends on matrix storage format
    (void)jacobian;
}

void FaultCohesiveKin::poststep(double t, double dt) {
    (void)dt;
    
    // Update slip rates
    for (size_t i = 0; i < vertices.size(); ++i) {
        double rate_ll, rate_rev, rate_open;
        getSlipRateAtTime(t, i, rate_ll, rate_rev, rate_open);
        
        slip_field.slip_rate_ll[i] = rate_ll;
        slip_field.slip_rate_reverse[i] = rate_rev;
        slip_field.slip_rate_opening[i] = rate_open;
    }
}

void FaultCohesiveKin::getSlipAtTime(double t, size_t vertex_idx,
                                    double& slip_ll, double& slip_rev, double& slip_open) const {
    if (vertex_idx >= kin_params.size()) {
        slip_ll = slip_rev = slip_open = 0.0;
        return;
    }
    
    const auto& p = kin_params[vertex_idx];
    
    // Configure slip time function rise time if needed
    if (auto* ramp = dynamic_cast<SlipTimeFnRamp*>(slip_time_fn.get())) {
        const_cast<SlipTimeFnRamp*>(ramp)->setRiseTime(p.rise_time);
    } else if (auto* liu = dynamic_cast<SlipTimeFnLiu*>(slip_time_fn.get())) {
        const_cast<SlipTimeFnLiu*>(liu)->setRiseTime(p.rise_time);
    } else if (auto* brune = dynamic_cast<SlipTimeFnBrune*>(slip_time_fn.get())) {
        const_cast<SlipTimeFnBrune*>(brune)->setCharacteristicTime(p.rise_time / 2.5);
    }
    
    slip_ll = slip_time_fn->computeSlip(t, p.slip_time, p.final_slip_ll);
    slip_rev = slip_time_fn->computeSlip(t, p.slip_time, p.final_slip_reverse);
    slip_open = slip_time_fn->computeSlip(t, p.slip_time, p.final_slip_opening);
}

void FaultCohesiveKin::getSlipRateAtTime(double t, size_t vertex_idx,
                                        double& rate_ll, double& rate_rev, double& rate_open) const {
    if (vertex_idx >= kin_params.size()) {
        rate_ll = rate_rev = rate_open = 0.0;
        return;
    }
    
    const auto& p = kin_params[vertex_idx];
    
    rate_ll = slip_time_fn->computeSlipRate(t, p.slip_time, p.final_slip_ll);
    rate_rev = slip_time_fn->computeSlipRate(t, p.slip_time, p.final_slip_reverse);
    rate_open = slip_time_fn->computeSlipRate(t, p.slip_time, p.final_slip_opening);
}

// =============================================================================
// FaultCohesiveDyn Implementation
// =============================================================================

FaultCohesiveDyn::FaultCohesiveDyn()
    : FaultCohesive(CohesiveFaultType::DYNAMIC),
      pert_type(TractionPerturbationType::UNIFORM),
      pert_amplitude(0.0), pert_duration(0.0), pert_width(1000.0),
      current_time(0.0), cumulative_frictional_work(0.0) {
    pert_center[0] = pert_center[1] = pert_center[2] = 0.0;
    friction_model = std::make_unique<SlipWeakeningFriction>();
}

void FaultCohesiveDyn::configure(const std::map<std::string, std::string>& config) {
    FaultCohesive::configure(config);
    
    // Set friction model
    std::string friction_type = parseString(config, "friction_model", "slip_weakening");
    setFrictionModelType(friction_type);
    friction_model->configure(config);
    
    // Initial traction
    double tau_ll = parseDouble(config, "initial_traction_ll", 0.0);
    double tau_ud = parseDouble(config, "initial_traction_ud", 0.0);
    double sigma_n = parseDouble(config, "initial_traction_normal", -50e6);
    setUniformInitialTraction(tau_ll, tau_ud, sigma_n);
    
    // Traction perturbation
    pert_amplitude = parseDouble(config, "perturbation_amplitude", 0.0);
    pert_duration = parseDouble(config, "perturbation_duration", 0.0);
    pert_center[0] = parseDouble(config, "perturbation_center_x", 0.0);
    pert_center[1] = parseDouble(config, "perturbation_center_y", 0.0);
    pert_center[2] = parseDouble(config, "perturbation_center_z", 0.0);
    pert_width = parseDouble(config, "perturbation_width", 1000.0);
    
    std::string pert_type_str = parseString(config, "perturbation_type", "gaussian");
    std::transform(pert_type_str.begin(), pert_type_str.end(), 
                  pert_type_str.begin(), ::tolower);
    if (pert_type_str == "uniform") {
        pert_type = TractionPerturbationType::UNIFORM;
    } else if (pert_type_str == "gaussian") {
        pert_type = TractionPerturbationType::GAUSSIAN;
    } else if (pert_type_str == "asperity") {
        pert_type = TractionPerturbationType::ASPERITY;
    }
}

void FaultCohesiveDyn::setFrictionModel(std::unique_ptr<FaultFrictionModel> model) {
    friction_model = std::move(model);
}

void FaultCohesiveDyn::setFrictionModelType(const std::string& type) {
    friction_model = createFaultFrictionModel(type);
}

void FaultCohesiveDyn::setInitialTraction(const FaultTractionField& traction) {
    initial_traction = traction;
}

void FaultCohesiveDyn::setUniformInitialTraction(double tau_ll, double tau_ud, double sigma_n) {
    initial_traction.resize(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        initial_traction.traction_shear_ll[i] = tau_ll;
        initial_traction.traction_shear_ud[i] = tau_ud;
        initial_traction.traction_normal[i] = sigma_n;
    }
}

void FaultCohesiveDyn::setInitialState(const std::vector<DynamicFaultState>& states) {
    dynamic_states = states;
}

void FaultCohesiveDyn::setTractionPerturbation(TractionPerturbationType type,
                                              double amplitude, double duration,
                                              double cx, double cy, double cz,
                                              double width) {
    pert_type = type;
    pert_amplitude = amplitude;
    pert_duration = duration;
    pert_center[0] = cx;
    pert_center[1] = cy;
    pert_center[2] = cz;
    pert_width = width;
}

void FaultCohesiveDyn::setTractionPerturbationFunction(
    std::function<void(double, double, double, double, double&)> fn) {
    pert_function = fn;
    pert_type = TractionPerturbationType::SPATIALLY_VARIABLE;
}

void FaultCohesiveDyn::initialize() {
    // Initialize dynamic states
    dynamic_states.resize(vertices.size());
    
    if (initial_traction.size() == 0) {
        initial_traction.resize(vertices.size());
    }
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        dynamic_states[i].state_variable = 1e6;  // Large θ initially
        dynamic_states[i].slip_rate = 0.0;
        dynamic_states[i].effective_normal = -initial_traction.traction_normal[i];
        dynamic_states[i].shear_stress = std::sqrt(
            initial_traction.traction_shear_ll[i] * initial_traction.traction_shear_ll[i] +
            initial_traction.traction_shear_ud[i] * initial_traction.traction_shear_ud[i]
        );
        dynamic_states[i].friction = friction_model->computeFriction(dynamic_states[i]);
        dynamic_states[i].is_locked = true;
        dynamic_states[i].has_ruptured = false;
        dynamic_states[i].rupture_time = -1.0;
    }
    
    cumulative_frictional_work = 0.0;
}

void FaultCohesiveDyn::prestep(double t, double dt) {
    current_time = t;
    
    // Apply traction perturbation
    if (pert_amplitude > 0 && t < pert_duration) {
        for (size_t i = 0; i < vertices.size(); ++i) {
            double dtau = computeTractionPerturbation(t, 
                vertices[i].coords[0], vertices[i].coords[1], vertices[i].coords[2]);
            
            // Add perturbation to shear traction (in strike direction)
            traction_field.traction_shear_ll[i] = 
                initial_traction.traction_shear_ll[i] + dtau;
        }
    }
}

double FaultCohesiveDyn::computeTractionPerturbation(double t, double x, double y, double z) const {
    if (t > pert_duration) return 0.0;
    
    // Time factor (smooth ramp up)
    double time_factor = (t < pert_duration) ? 
        0.5 * (1.0 - std::cos(M_PI * t / pert_duration)) : 1.0;
    
    // Spatial factor
    double dx = x - pert_center[0];
    double dy = y - pert_center[1];
    double dz = z - pert_center[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    
    double spatial_factor;
    switch (pert_type) {
        case TractionPerturbationType::UNIFORM:
            spatial_factor = 1.0;
            break;
        case TractionPerturbationType::GAUSSIAN:
            spatial_factor = std::exp(-r2 / (2.0 * pert_width * pert_width));
            break;
        case TractionPerturbationType::ASPERITY:
            spatial_factor = (r2 < pert_width * pert_width) ? 1.0 : 0.0;
            break;
        default:
            spatial_factor = 1.0;
    }
    
    if (pert_function) {
        double dtau;
        pert_function(t, x, y, z, dtau);
        return dtau;
    }
    
    return pert_amplitude * time_factor * spatial_factor;
}

void FaultCohesiveDyn::computeResidual(const double* solution, double* residual) {
    // For dynamic faults, the residual enforces:
    // τ = f(V, θ) × σ_n^eff  (friction constitutive equation)
    
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& v = vertices[i];
        if (v.vertex_negative < 0 || v.vertex_positive < 0) continue;
        
        auto& state = dynamic_states[i];
        
        // Get displacement jump from solution
        double u_neg[3], u_pos[3];
        for (int j = 0; j < 3; ++j) {
            u_neg[j] = solution[3 * v.vertex_negative + j];
            u_pos[j] = solution[3 * v.vertex_positive + j];
        }
        
        double du[3] = {u_pos[0] - u_neg[0], 
                       u_pos[1] - u_neg[1], 
                       u_pos[2] - u_neg[2]};
        
        // Transform to fault coords
        double slip_ll, slip_ud, slip_n;
        transformToFaultCoords(du, v, slip_ll, slip_ud, slip_n);
        
        // Get traction from Lagrange multiplier if available
        double traction_n = -state.effective_normal;  // Normal traction
        double traction_ll = traction_field.traction_shear_ll[i];
        double traction_ud = traction_field.traction_shear_ud[i];
        
        // Compute friction coefficient
        state.friction = friction_model->computeFriction(state);
        
        // Friction strength
        double tau_strength = state.friction * state.effective_normal;
        double tau_applied = std::sqrt(traction_ll * traction_ll + traction_ud * traction_ud);
        
        // Contact constraints
        if (slip_n > 0) {
            // Tensile opening - no normal traction
            // Residual: λ_n = 0
        } else {
            // Compression - enforce non-penetration
            // Residual: slip_n >= 0 (complementarity)
        }
        
        // Friction constraint
        if (tau_applied < tau_strength - 1e-6) {
            // Stick: slip rate should be zero
            state.is_locked = true;
        } else {
            // Slip: traction equals friction strength
            state.is_locked = false;
            if (!state.has_ruptured) {
                state.has_ruptured = true;
                state.rupture_time = current_time;
            }
        }
        
        // Update slip field
        slip_field.slip_left_lateral[i] = slip_ll;
        slip_field.slip_reverse[i] = slip_ud;
        slip_field.slip_opening[i] = slip_n;
    }
}

void FaultCohesiveDyn::computeJacobian(const double* solution, double* jacobian) {
    (void)solution;
    (void)jacobian;
    // Jacobian for friction equations - linearization of f(V, θ)
}

void FaultCohesiveDyn::poststep(double t, double dt) {
    // Update friction state
    updateFrictionState(dt);
    
    // Track energy
    for (size_t i = 0; i < vertices.size(); ++i) {
        const auto& v = vertices[i];
        const auto& state = dynamic_states[i];
        
        if (!state.is_locked) {
            double tau = std::sqrt(
                traction_field.traction_shear_ll[i] * traction_field.traction_shear_ll[i] +
                traction_field.traction_shear_ud[i] * traction_field.traction_shear_ud[i]
            );
            double slip_increment = state.slip_rate * dt;
            cumulative_frictional_work += tau * slip_increment * v.area;
        }
    }
}

void FaultCohesiveDyn::updateFrictionState(double dt) {
    for (size_t i = 0; i < dynamic_states.size(); ++i) {
        friction_model->updateState(dynamic_states[i], dt);
    }
}

bool FaultCohesiveDyn::hasRuptured() const {
    for (const auto& state : dynamic_states) {
        if (state.has_ruptured) return true;
    }
    return false;
}

double FaultCohesiveDyn::getRuptureTime() const {
    double min_time = std::numeric_limits<double>::max();
    for (const auto& state : dynamic_states) {
        if (state.has_ruptured && state.rupture_time < min_time) {
            min_time = state.rupture_time;
        }
    }
    return min_time;
}

double FaultCohesiveDyn::getMaxSlipRate() const {
    double max_rate = 0.0;
    for (const auto& state : dynamic_states) {
        max_rate = std::max(max_rate, state.slip_rate);
    }
    return max_rate;
}

double FaultCohesiveDyn::getMaxSlip() const {
    double max_slip = 0.0;
    for (size_t i = 0; i < slip_field.size(); ++i) {
        double slip = std::sqrt(
            slip_field.slip_left_lateral[i] * slip_field.slip_left_lateral[i] +
            slip_field.slip_reverse[i] * slip_field.slip_reverse[i]
        );
        max_slip = std::max(max_slip, slip);
    }
    return max_slip;
}

double FaultCohesiveDyn::getFrictionalWork() const {
    return cumulative_frictional_work;
}

double FaultCohesiveDyn::getFractureEnergy() const {
    // Fracture energy from slip-weakening: G_c = 0.5 * (f_s - f_d) * σ_n * D_c
    if (auto* sw = dynamic_cast<SlipWeakeningFriction*>(friction_model.get())) {
        double avg_sigma_n = 0.0;
        for (const auto& state : dynamic_states) {
            avg_sigma_n += state.effective_normal;
        }
        avg_sigma_n /= dynamic_states.size();
        
        return 0.5 * (sw->getStaticCoefficient() - sw->getDynamicCoefficient()) *
               avg_sigma_n * sw->getCriticalSlipDistance();
    }
    return 0.0;
}

double FaultCohesiveDyn::getRadiatedEnergy() const {
    // Simplified: E_r = E_total - E_friction - E_fracture
    // This is a rough estimate
    return 0.0;
}

// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<SlipTimeFn> createSlipTimeFn(SlipTimeFnType type) {
    switch (type) {
        case SlipTimeFnType::STEP:
            return std::make_unique<SlipTimeFnStep>();
        case SlipTimeFnType::RAMP:
            return std::make_unique<SlipTimeFnRamp>();
        case SlipTimeFnType::BRUNE:
            return std::make_unique<SlipTimeFnBrune>();
        case SlipTimeFnType::LIU:
            return std::make_unique<SlipTimeFnLiu>();
        case SlipTimeFnType::CONSTANT_RATE:
            return std::make_unique<SlipTimeFnConstantRate>();
        case SlipTimeFnType::TIME_HISTORY:
            return std::make_unique<SlipTimeFnTimeHistory>();
        default:
            return std::make_unique<SlipTimeFnBrune>();
    }
}

std::unique_ptr<SlipTimeFn> createSlipTimeFn(const std::string& type_str) {
    return createSlipTimeFn(parseSlipTimeFnType(type_str));
}

std::unique_ptr<FaultFrictionModel> createFaultFrictionModel(const std::string& type_str) {
    std::string s = type_str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    
    if (s == "static" || s == "constant") {
        return std::make_unique<StaticFriction>();
    } else if (s == "slip_weakening" || s == "slipweakening" || s == "sw") {
        return std::make_unique<SlipWeakeningFriction>();
    } else if (s == "rate_state_aging" || s == "rs_aging" || s == "aging") {
        return std::make_unique<RateStateFrictionAging>();
    } else if (s == "rate_state_slip" || s == "rs_slip" || s == "slip") {
        return std::make_unique<RateStateFrictionSlip>();
    }
    
    return std::make_unique<SlipWeakeningFriction>();
}

std::unique_ptr<FaultCohesive> createFaultCohesive(
    const std::map<std::string, std::string>& config) {
    
    std::string type_str = parseString(config, "fault_type", "kinematic");
    CohesiveFaultType type = parseCohesiveFaultType(type_str);
    
    std::unique_ptr<FaultCohesive> fault;
    
    switch (type) {
        case CohesiveFaultType::KINEMATIC:
            fault = std::make_unique<FaultCohesiveKin>();
            break;
        case CohesiveFaultType::DYNAMIC:
            fault = std::make_unique<FaultCohesiveDyn>();
            break;
        default:
            fault = std::make_unique<FaultCohesiveKin>();
    }
    
    fault->configure(config);
    return fault;
}

SlipTimeFnType parseSlipTimeFnType(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    
    if (s == "step" || s == "heaviside") {
        return SlipTimeFnType::STEP;
    } else if (s == "ramp" || s == "linear") {
        return SlipTimeFnType::RAMP;
    } else if (s == "brune" || s == "exponential") {
        return SlipTimeFnType::BRUNE;
    } else if (s == "liu" || s == "regularized") {
        return SlipTimeFnType::LIU;
    } else if (s == "constant_rate" || s == "creep") {
        return SlipTimeFnType::CONSTANT_RATE;
    } else if (s == "time_history" || s == "history") {
        return SlipTimeFnType::TIME_HISTORY;
    }
    
    return SlipTimeFnType::BRUNE;
}

CohesiveFaultType parseCohesiveFaultType(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    
    if (s == "kinematic" || s == "prescribed" || s == "kin") {
        return CohesiveFaultType::KINEMATIC;
    } else if (s == "dynamic" || s == "spontaneous" || s == "dyn") {
        return CohesiveFaultType::DYNAMIC;
    } else if (s == "impulse" || s == "greens") {
        return CohesiveFaultType::IMPULSE;
    } else if (s == "absorbing") {
        return CohesiveFaultType::ABSORBING;
    }
    
    return CohesiveFaultType::KINEMATIC;
}

// =============================================================================
// FaultSystem Implementation
// =============================================================================

FaultSystem::FaultSystem() : compute_interaction(false) {}

void FaultSystem::configure(const std::vector<std::map<std::string, std::string>>& configs) {
    faults.clear();
    for (const auto& config : configs) {
        addFault(createFaultCohesive(config));
    }
}

void FaultSystem::addFault(std::unique_ptr<FaultCohesive> fault) {
    faults.push_back(std::move(fault));
}

void FaultSystem::initialize() {
    for (auto& fault : faults) {
        fault->initialize();
    }
}

void FaultSystem::prestep(double t, double dt) {
    for (auto& fault : faults) {
        fault->prestep(t, dt);
    }
}

void FaultSystem::computeResidual(const double* solution, double* residual) {
    for (auto& fault : faults) {
        fault->computeResidual(solution, residual);
    }
}

void FaultSystem::computeJacobian(const double* solution, double* jacobian) {
    for (auto& fault : faults) {
        fault->computeJacobian(solution, jacobian);
    }
}

void FaultSystem::poststep(double t, double dt) {
    for (auto& fault : faults) {
        fault->poststep(t, dt);
    }
}

FaultCohesive* FaultSystem::getFault(const std::string& name) {
    for (auto& fault : faults) {
        if (fault->getName() == name) {
            return fault.get();
        }
    }
    return nullptr;
}

void FaultSystem::writeOutput(const std::string& base_filename, double time) {
    for (size_t i = 0; i < faults.size(); ++i) {
        std::stringstream ss;
        ss << base_filename << "_" << faults[i]->getName() << "_" 
           << std::setfill('0') << std::setw(6) << static_cast<int>(time) << ".dat";
        faults[i]->writeOutput(ss.str(), time);
    }
}

// =============================================================================
// FaultUtils Implementation
// =============================================================================

namespace FaultUtils {

double computeMomentMagnitude(double slip, double area, double shear_modulus) {
    double M0 = computeSeismicMoment(slip, area, shear_modulus);
    return (2.0/3.0) * (std::log10(M0) - 9.1);
}

double computeSeismicMoment(double slip, double area, double shear_modulus) {
    return shear_modulus * area * slip;
}

double computeStressDrop(double slip, double radius, double shear_modulus) {
    // For circular crack: Δσ = (7/16) * μ * D / r
    return (7.0/16.0) * shear_modulus * slip / radius;
}

double computeCornerFrequency(double moment, double stress_drop, 
                             double shear_velocity, double density) {
    // Brune corner frequency: f_c = 0.49 * β * (Δσ/M0)^(1/3)
    double beta = shear_velocity;
    return 0.49 * beta * std::pow(stress_drop / moment, 1.0/3.0);
}

double computeNucleationLength(double a, double b, double dc, 
                              double sigma_n, double shear_modulus) {
    // L_c = G * D_c / ((b - a) * σ_n)
    if (b <= a) return std::numeric_limits<double>::max();
    return shear_modulus * dc / ((b - a) * sigma_n);
}

double computeCriticalStiffness(double b_minus_a, double sigma_n, double dc) {
    // k_c = (b - a) * σ_n / D_c
    return b_minus_a * sigma_n / dc;
}

double estimateRiseTime(double slip, double slip_rate) {
    if (slip_rate < 1e-15) return std::numeric_limits<double>::max();
    return slip / slip_rate;
}

double computeRuptureVelocity(double shear_wave_velocity, 
                             double rayleigh_velocity,
                             bool supershear) {
    if (supershear) {
        // Supershear: V_r > V_s, typically ~1.4 * V_s
        return 1.4 * shear_wave_velocity;
    } else {
        // Sub-Rayleigh: V_r ≈ 0.9 * V_R
        return 0.9 * rayleigh_velocity;
    }
}

} // namespace FaultUtils

} // namespace FSRM
