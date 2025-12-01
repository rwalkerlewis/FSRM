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

// =============================================================================
// Multi-Physics Property Configuration
// =============================================================================

void FaultPermeabilityTensor::configure(const std::map<std::string, std::string>& config) {
    k_strike = parseDouble(config, "fault_permeability_strike", 1e-12);
    k_dip = parseDouble(config, "fault_permeability_dip", 1e-12);
    k_normal = parseDouble(config, "fault_permeability_normal", 1e-15);
    k_strike_dip = parseDouble(config, "fault_permeability_strike_dip", 0.0);
    k_strike_normal = parseDouble(config, "fault_permeability_strike_normal", 0.0);
    k_dip_normal = parseDouble(config, "fault_permeability_dip_normal", 0.0);
    
    // Store reference values
    k_strike_ref = k_strike;
    k_dip_ref = k_dip;
    k_normal_ref = k_normal;
    
    // Check for isotropic specification
    double k_iso = parseDouble(config, "fault_permeability", -1.0);
    if (k_iso > 0) {
        setIsotropic(k_iso);
    }
    
    // Check for principal specification
    double k_par = parseDouble(config, "fault_permeability_parallel", -1.0);
    double k_perp = parseDouble(config, "fault_permeability_perpendicular", -1.0);
    if (k_par > 0 && k_perp > 0) {
        setPrincipal(k_par, k_perp);
    }
}

void FaultHydraulicProperties::configure(const std::map<std::string, std::string>& config) {
    porosity = parseDouble(config, "fault_porosity", 0.1);
    porosity_ref = porosity;
    fault_thickness = parseDouble(config, "fault_thickness", 0.01);
    damage_zone_thickness = parseDouble(config, "damage_zone_thickness", 10.0);
    specific_storage = parseDouble(config, "fault_specific_storage", 1e-9);
    fluid_compressibility = parseDouble(config, "fluid_compressibility", 4.5e-10);
    solid_compressibility = parseDouble(config, "solid_compressibility", 1e-11);
    biot_coefficient = parseDouble(config, "fault_biot_coefficient", 0.9);
    fluid_viscosity = parseDouble(config, "fault_fluid_viscosity", 1e-3);
    fluid_density = parseDouble(config, "fault_fluid_density", 1000.0);
    
    permeability.configure(config);
}

void FaultThermalProperties::configure(const std::map<std::string, std::string>& config) {
    conductivity_parallel = parseDouble(config, "fault_thermal_conductivity_parallel", 3.0);
    conductivity_normal = parseDouble(config, "fault_thermal_conductivity_normal", 2.0);
    specific_heat_solid = parseDouble(config, "fault_specific_heat_solid", 900.0);
    specific_heat_fluid = parseDouble(config, "fault_specific_heat_fluid", 4186.0);
    thermal_expansion_solid = parseDouble(config, "fault_thermal_expansion_solid", 1e-5);
    thermal_expansion_fluid = parseDouble(config, "fault_thermal_expansion_fluid", 3e-4);
    reference_temperature = parseDouble(config, "fault_reference_temperature", 293.15);
    weakening_temperature = parseDouble(config, "fault_weakening_temperature", 1200.0);
    contact_diameter = parseDouble(config, "fault_contact_diameter", 10e-6);
    thermal_diffusivity = parseDouble(config, "fault_thermal_diffusivity", 1e-6);
}

void FaultChemicalProperties::configure(const std::map<std::string, std::string>& config) {
    quartz_fraction = parseDouble(config, "fault_quartz_fraction", 0.6);
    calcite_fraction = parseDouble(config, "fault_calcite_fraction", 0.2);
    clay_fraction = parseDouble(config, "fault_clay_fraction", 0.2);
    dissolution_rate = parseDouble(config, "fault_dissolution_rate", 1e-12);
    precipitation_rate = parseDouble(config, "fault_precipitation_rate", 1e-13);
    reactive_surface_area = parseDouble(config, "fault_reactive_surface_area", 1000.0);
    initial_ph = parseDouble(config, "fault_initial_ph", 7.0);
    co2_concentration = parseDouble(config, "fault_co2_concentration", 0.0);
}

void PermeabilityEvolutionParams::configure(const std::map<std::string, std::string>& config) {
    std::string model_str = parseString(config, "permeability_update_model", "constant");
    model = parsePermeabilityUpdateModel(model_str);
    
    slip_coefficient = parseDouble(config, "permeability_slip_coefficient", 10.0);
    dilation_coefficient = parseDouble(config, "permeability_dilation_coefficient", 100.0);
    stress_coefficient = parseDouble(config, "permeability_stress_coefficient", 0.1);
    stress_reference = parseDouble(config, "permeability_stress_reference", 50e6);
    shear_coefficient = parseDouble(config, "permeability_shear_coefficient", 50.0);
    damage_exponent = parseDouble(config, "permeability_damage_exponent", 3.0);
    k_min = parseDouble(config, "permeability_minimum", 1e-20);
    k_max = parseDouble(config, "permeability_maximum", 1e-10);
    enhancement_time = parseDouble(config, "permeability_enhancement_time", 1.0);
    recovery_time = parseDouble(config, "permeability_recovery_time", 3.15e7);
}

void PorosityEvolutionParams::configure(const std::map<std::string, std::string>& config) {
    std::string model_str = parseString(config, "porosity_update_model", "constant");
    model = parsePorosityUpdateModel(model_str);
    
    compaction_coefficient = parseDouble(config, "porosity_compaction_coefficient", 0.01);
    compaction_reference = parseDouble(config, "porosity_compaction_reference", 50e6);
    dilation_coefficient = parseDouble(config, "porosity_dilation_coefficient", 0.1);
    chemical_rate = parseDouble(config, "porosity_chemical_rate", 1e-12);
    phi_min = parseDouble(config, "porosity_minimum", 0.001);
    phi_max = parseDouble(config, "porosity_maximum", 0.5);
}

// =============================================================================
// DefaultFaultPropertyUpdater Implementation
// =============================================================================

DefaultFaultPropertyUpdater::DefaultFaultPropertyUpdater() {}

double DefaultFaultPropertyUpdater::slipDependentPermeability(double k0, double slip) const {
    // k = k0 * exp(alpha * slip)
    double factor = std::exp(perm_params.slip_coefficient * slip);
    return std::clamp(k0 * factor, perm_params.k_min, perm_params.k_max);
}

double DefaultFaultPropertyUpdater::dilationDependentPermeability(double k0, double dilation) const {
    // k = k0 * (1 + beta * dilation)
    // Only positive dilation (opening) enhances permeability
    double factor = 1.0 + perm_params.dilation_coefficient * std::max(0.0, dilation);
    return std::clamp(k0 * factor, perm_params.k_min, perm_params.k_max);
}

double DefaultFaultPropertyUpdater::stressDependentPermeability(double k0, double sigma_n_eff) const {
    // k = k0 * exp(-gamma * sigma_n_eff / sigma_ref)
    // Higher effective stress reduces permeability
    double factor = std::exp(-perm_params.stress_coefficient * sigma_n_eff / perm_params.stress_reference);
    return std::clamp(k0 * factor, perm_params.k_min, perm_params.k_max);
}

double DefaultFaultPropertyUpdater::shearEnhancedPermeability(double k0, double shear_strain) const {
    // k = k0 * (1 + delta * |gamma|)
    double factor = 1.0 + perm_params.shear_coefficient * std::abs(shear_strain);
    return std::clamp(k0 * factor, perm_params.k_min, perm_params.k_max);
}

double DefaultFaultPropertyUpdater::damageDependentPermeability(double k0, double damage) const {
    // k = k0 * (1 - D)^(-n) where D in [0, 1)
    double clamped_damage = std::clamp(damage, 0.0, 0.999);
    double factor = std::pow(1.0 - clamped_damage, -perm_params.damage_exponent);
    return std::clamp(k0 * factor, perm_params.k_min, perm_params.k_max);
}

FaultPermeabilityTensor DefaultFaultPropertyUpdater::updatePermeability(
    const FaultMultiPhysicsState& state,
    const FaultHydraulicProperties& props,
    double /*dt*/) {
    
    FaultPermeabilityTensor result = state.current_permeability;
    
    // Apply selected model(s)
    double scale_factor = 1.0;
    
    switch (perm_params.model) {
        case PermeabilityUpdateModel::CONSTANT:
            // No update
            return result;
            
        case PermeabilityUpdateModel::SLIP_DEPENDENT:
            scale_factor = std::exp(perm_params.slip_coefficient * state.slip_total);
            break;
            
        case PermeabilityUpdateModel::DILATION_DEPENDENT:
            scale_factor = 1.0 + perm_params.dilation_coefficient * 
                          std::max(0.0, state.volumetric_strain);
            break;
            
        case PermeabilityUpdateModel::STRESS_DEPENDENT:
            scale_factor = std::exp(-perm_params.stress_coefficient * 
                          state.effective_normal_stress / perm_params.stress_reference);
            break;
            
        case PermeabilityUpdateModel::SHEAR_ENHANCED:
            scale_factor = 1.0 + perm_params.shear_coefficient * 
                          std::abs(state.shear_strain);
            break;
            
        case PermeabilityUpdateModel::DAMAGE_DEPENDENT:
            scale_factor = std::pow(1.0 - std::clamp(state.damage, 0.0, 0.999),
                                   -perm_params.damage_exponent);
            break;
            
        case PermeabilityUpdateModel::COMBINED:
            // Apply all effects multiplicatively
            scale_factor *= std::exp(perm_params.slip_coefficient * state.slip_total);
            scale_factor *= (1.0 + perm_params.dilation_coefficient * 
                           std::max(0.0, state.volumetric_strain));
            scale_factor *= std::exp(-perm_params.stress_coefficient * 
                           state.effective_normal_stress / perm_params.stress_reference);
            break;
    }
    
    // Apply scale factor to all permeability components
    result.k_strike = std::clamp(result.k_strike_ref * scale_factor, 
                                perm_params.k_min, perm_params.k_max);
    result.k_dip = std::clamp(result.k_dip_ref * scale_factor,
                             perm_params.k_min, perm_params.k_max);
    result.k_normal = std::clamp(result.k_normal_ref * scale_factor,
                                perm_params.k_min, perm_params.k_max);
    
    return result;
}

double DefaultFaultPropertyUpdater::updatePorosity(
    const FaultMultiPhysicsState& state,
    const FaultHydraulicProperties& props,
    double /*dt*/) {
    
    double phi = props.porosity_ref;
    
    switch (poro_params.model) {
        case PorosityUpdateModel::CONSTANT:
            return props.porosity;
            
        case PorosityUpdateModel::COMPACTION:
            // φ = φ0 * exp(-α * Δσ'_mean / σ_ref)
            phi *= std::exp(-poro_params.compaction_coefficient * 
                           state.effective_normal_stress / poro_params.compaction_reference);
            break;
            
        case PorosityUpdateModel::DILATION:
            // φ = φ0 + β * Δε_v
            phi += poro_params.dilation_coefficient * state.volumetric_strain;
            break;
            
        case PorosityUpdateModel::CHEMICAL:
            // φ changes from dissolution/precipitation
            // This would integrate over time
            phi += poro_params.chemical_rate * state.reaction_rate;
            break;
            
        case PorosityUpdateModel::COMBINED:
            phi *= std::exp(-poro_params.compaction_coefficient * 
                           state.effective_normal_stress / poro_params.compaction_reference);
            phi += poro_params.dilation_coefficient * state.volumetric_strain;
            break;
    }
    
    return std::clamp(phi, poro_params.phi_min, poro_params.phi_max);
}

double DefaultFaultPropertyUpdater::updateFriction(
    const FaultMultiPhysicsState& state,
    double base_friction) {
    
    // Thermal weakening effect
    double T_ref = 293.15;
    double T_weak = 1200.0;  // Typical weakening temperature
    
    if (state.temperature > T_ref) {
        // Linear reduction in friction with temperature
        double reduction = (state.temperature - T_ref) / (T_weak - T_ref);
        reduction = std::clamp(reduction, 0.0, 0.9);  // Max 90% reduction
        return base_friction * (1.0 - reduction);
    }
    
    return base_friction;
}

// =============================================================================
// FaultCohesiveTHM Implementation
// =============================================================================

FaultCohesiveTHM::FaultCohesiveTHM() : 
    FaultCohesiveDyn(),
    coupling_mode(FaultCouplingMode::MECHANICAL_ONLY) {
    
    property_updater = std::make_shared<DefaultFaultPropertyUpdater>();
}

void FaultCohesiveTHM::configure(const std::map<std::string, std::string>& config) {
    // Configure base class
    FaultCohesiveDyn::configure(config);
    
    // Parse coupling mode
    std::string mode_str = parseString(config, "fault_coupling_mode", "mechanical");
    coupling_mode = parseFaultCouplingMode(mode_str);
    
    // Configure multi-physics properties
    if (coupling_mode != FaultCouplingMode::MECHANICAL_ONLY) {
        hydraulic_props.configure(config);
    }
    
    if (coupling_mode == FaultCouplingMode::THERMOELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        thermal_props.configure(config);
    }
    
    if (coupling_mode == FaultCouplingMode::THMC) {
        chemical_props.configure(config);
    }
    
    // Configure property evolution
    perm_evolution.configure(config);
    poro_evolution.configure(config);
    
    // Set up the default updater with parameters
    auto default_updater = std::dynamic_pointer_cast<DefaultFaultPropertyUpdater>(property_updater);
    if (default_updater) {
        default_updater->setPermeabilityParams(perm_evolution);
        default_updater->setPorosityParams(poro_evolution);
    }
}

void FaultCohesiveTHM::setHydraulicProperties(const FaultHydraulicProperties& props) {
    hydraulic_props = props;
}

void FaultCohesiveTHM::setThermalProperties(const FaultThermalProperties& props) {
    thermal_props = props;
}

void FaultCohesiveTHM::setChemicalProperties(const FaultChemicalProperties& props) {
    chemical_props = props;
}

void FaultCohesiveTHM::setHydraulicField(const FaultHydraulicField& field) {
    hydraulic_field = field;
}

void FaultCohesiveTHM::setThermalField(const FaultThermalField& field) {
    thermal_field = field;
}

void FaultCohesiveTHM::setPermeabilityEvolution(const PermeabilityEvolutionParams& params) {
    perm_evolution = params;
    auto default_updater = std::dynamic_pointer_cast<DefaultFaultPropertyUpdater>(property_updater);
    if (default_updater) {
        default_updater->setPermeabilityParams(params);
    }
}

void FaultCohesiveTHM::setPorosityEvolution(const PorosityEvolutionParams& params) {
    poro_evolution = params;
    auto default_updater = std::dynamic_pointer_cast<DefaultFaultPropertyUpdater>(property_updater);
    if (default_updater) {
        default_updater->setPorosityParams(params);
    }
}

void FaultCohesiveTHM::setPropertyUpdater(std::shared_ptr<FaultPropertyUpdateCallback> updater) {
    property_updater = updater;
}

void FaultCohesiveTHM::initializeMultiPhysicsState() {
    size_t num_vertices = vertices.size();
    mp_states.resize(num_vertices);
    
    // Initialize hydraulic field if needed
    if (coupling_mode == FaultCouplingMode::POROELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        hydraulic_field.resize(num_vertices);
        
        // Set initial permeability from properties
        for (size_t i = 0; i < num_vertices; ++i) {
            hydraulic_field.porosity[i] = hydraulic_props.porosity;
            hydraulic_field.permeability[i] = hydraulic_props.permeability;
            mp_states[i].current_porosity = hydraulic_props.porosity;
            mp_states[i].current_permeability = hydraulic_props.permeability;
        }
    }
    
    // Initialize thermal field if needed
    if (coupling_mode == FaultCouplingMode::THERMOELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        thermal_field.resize(num_vertices);
        
        for (size_t i = 0; i < num_vertices; ++i) {
            thermal_field.temperature[i] = thermal_props.reference_temperature;
            mp_states[i].temperature = thermal_props.reference_temperature;
        }
    }
}

void FaultCohesiveTHM::initialize() {
    // Initialize base class
    FaultCohesiveDyn::initialize();
    
    // Initialize multi-physics state
    initializeMultiPhysicsState();
}

void FaultCohesiveTHM::prestep(double t, double dt) {
    // Call base class prestep
    FaultCohesiveDyn::prestep(t, dt);
    
    // Update effective stress from pore pressure
    if (coupling_mode != FaultCouplingMode::MECHANICAL_ONLY) {
        updateEffectiveStress();
    }
    
    // Update thermal effects
    if (coupling_mode == FaultCouplingMode::THERMOELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        computeFrictionalHeating();
    }
}

void FaultCohesiveTHM::computeResidual(const double* solution, double* residual) {
    // Base mechanical residual
    FaultCohesiveDyn::computeResidual(solution, residual);
    
    // Add poroelastic coupling terms
    if (coupling_mode == FaultCouplingMode::POROELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        addPoroelasticResidual(solution, residual);
    }
    
    // Add thermal coupling terms
    if (coupling_mode == FaultCouplingMode::THERMOELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        addThermalResidual(solution, residual);
    }
}

void FaultCohesiveTHM::computeJacobian(const double* solution, double* jacobian) {
    // Base mechanical Jacobian
    FaultCohesiveDyn::computeJacobian(solution, jacobian);
    
    // Add coupling contributions to Jacobian
    // The coupling terms would add:
    // - ∂R_mech/∂p (pore pressure effect on effective stress)
    // - ∂R_flow/∂u (volume change effect on flow)
    // - ∂R_mech/∂T (thermal stress)
    // - ∂R_thermal/∂u (deformation effect on heat)
    
    // This is a placeholder - actual implementation would depend on
    // the specific FE discretization and DOF ordering
}

void FaultCohesiveTHM::poststep(double t, double dt) {
    // Update fault properties based on current state
    updateFaultProperties(dt);
    
    // Update thermal state
    if (coupling_mode == FaultCouplingMode::THERMOELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        updateHeatTransfer(dt);
        updateThermalPressurization(dt);
    }
    
    // Update fluid flow along fault
    if (coupling_mode == FaultCouplingMode::POROELASTIC ||
        coupling_mode == FaultCouplingMode::THERMO_HYDRO_MECHANICAL ||
        coupling_mode == FaultCouplingMode::THMC) {
        updateFluidFlow(dt);
    }
    
    // Call base class poststep
    FaultCohesiveDyn::poststep(t, dt);
}

void FaultCohesiveTHM::setPorePressureField(const std::vector<double>& pressure) {
    if (pressure.size() != hydraulic_field.pore_pressure.size()) {
        hydraulic_field.pore_pressure.resize(pressure.size());
    }
    hydraulic_field.pore_pressure = pressure;
    
    // Also update mp_states
    for (size_t i = 0; i < std::min(mp_states.size(), pressure.size()); ++i) {
        mp_states[i].pore_pressure = pressure[i];
    }
    
    updateEffectiveStress();
}

void FaultCohesiveTHM::setTemperatureField(const std::vector<double>& temperature) {
    if (temperature.size() != thermal_field.temperature.size()) {
        thermal_field.temperature.resize(temperature.size());
    }
    thermal_field.temperature = temperature;
    
    // Also update mp_states
    for (size_t i = 0; i < std::min(mp_states.size(), temperature.size()); ++i) {
        mp_states[i].temperature = temperature[i];
        mp_states[i].temperature_rise = temperature[i] - thermal_props.reference_temperature;
    }
}

void FaultCohesiveTHM::getFluidSourceTerms(std::vector<double>& source_minus,
                                          std::vector<double>& source_plus) const {
    size_t n = hydraulic_field.size();
    source_minus.resize(n, 0.0);
    source_plus.resize(n, 0.0);
    
    // Fluid flux across fault (normal direction) creates source/sink
    // for domains on either side
    for (size_t i = 0; i < n; ++i) {
        double flux_normal = hydraulic_field.fluid_flux_normal[i];
        // Flux from minus side to plus side
        source_minus[i] = -flux_normal;  // Sink for minus side
        source_plus[i] = flux_normal;    // Source for plus side
    }
}

void FaultCohesiveTHM::getHeatSourceTerms(std::vector<double>& heat_source) const {
    size_t n = thermal_field.size();
    heat_source.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        heat_source[i] = thermal_field.frictional_heating[i];
    }
}

const FaultMultiPhysicsState& FaultCohesiveTHM::getMultiPhysicsState(size_t vertex_idx) const {
    if (vertex_idx >= mp_states.size()) {
        static FaultMultiPhysicsState default_state;
        return default_state;
    }
    return mp_states[vertex_idx];
}

FaultPermeabilityTensor FaultCohesiveTHM::getPermeabilityAt(size_t vertex_idx) const {
    if (vertex_idx < hydraulic_field.permeability.size()) {
        return hydraulic_field.permeability[vertex_idx];
    }
    return hydraulic_props.permeability;
}

double FaultCohesiveTHM::getPorosityAt(size_t vertex_idx) const {
    if (vertex_idx < hydraulic_field.porosity.size()) {
        return hydraulic_field.porosity[vertex_idx];
    }
    return hydraulic_props.porosity;
}

void FaultCohesiveTHM::updateEffectiveStress() {
    for (size_t i = 0; i < mp_states.size(); ++i) {
        // Get total normal stress from dynamic state
        mp_states[i].normal_stress = dyn_states[i].normal_stress;
        
        // Get pore pressure
        double p = (i < hydraulic_field.pore_pressure.size()) ? 
                   hydraulic_field.pore_pressure[i] : 0.0;
        mp_states[i].pore_pressure = p;
        
        // Compute effective normal stress
        mp_states[i].effective_normal_stress = 
            FaultPoroelastic::effectiveNormalStress(
                mp_states[i].normal_stress, p, hydraulic_props.biot_coefficient);
        
        // Update the dynamic state effective stress for friction calculation
        // This is how pore pressure affects the fault behavior
        dyn_states[i].normal_stress = mp_states[i].effective_normal_stress;
    }
}

void FaultCohesiveTHM::updateThermalPressurization(double dt) {
    for (size_t i = 0; i < mp_states.size(); ++i) {
        if (mp_states[i].frictional_heating > 0 && dt > 0) {
            // Thermal pressurization coefficient
            double lambda_f = thermal_props.thermal_expansion_fluid;
            double beta_f = hydraulic_props.fluid_compressibility;
            double pressurization_coeff = lambda_f / beta_f;  // Pa/K
            
            // Compute pore pressure change from thermal pressurization
            double dp = FaultThermal::thermalPressurization(
                mp_states[i].slip_rate,
                mp_states[i].shear_stress,
                hydraulic_props.fault_thickness,
                thermal_props.thermal_diffusivity,
                hydraulic_props.permeability.k_normal / (hydraulic_props.fluid_viscosity * 
                    hydraulic_props.specific_storage),
                pressurization_coeff,
                dt);
            
            // Update pore pressure
            if (i < hydraulic_field.pore_pressure.size()) {
                hydraulic_field.pore_pressure[i] += dp;
                mp_states[i].pore_pressure = hydraulic_field.pore_pressure[i];
            }
        }
    }
}

void FaultCohesiveTHM::updateFluidFlow(double /*dt*/) {
    // Compute fluid flux along fault using Darcy's law
    // This is a simplified 1D flow along fault
    
    size_t n = mp_states.size();
    if (n < 2) return;
    
    // For each vertex, estimate flux from pressure gradient
    for (size_t i = 0; i < n; ++i) {
        double dp_ds = 0.0;  // Pressure gradient along strike
        double dp_dd = 0.0;  // Pressure gradient along dip
        
        // Estimate gradients (simplified finite difference)
        // In practice, this would use actual mesh connectivity
        if (i > 0 && i < n-1) {
            // Central difference (simplified - assumes uniform spacing)
            double dx = 100.0;  // Assume 100m spacing
            dp_ds = (hydraulic_field.pore_pressure[i+1] - 
                    hydraulic_field.pore_pressure[i-1]) / (2.0 * dx);
        }
        
        // Darcy velocity: q = -k/μ * (∇p - ρg)
        double k_s = hydraulic_field.permeability[i].k_strike;
        double k_d = hydraulic_field.permeability[i].k_dip;
        double mu = hydraulic_props.fluid_viscosity;
        
        hydraulic_field.fluid_flux_strike[i] = -k_s / mu * dp_ds;
        hydraulic_field.fluid_flux_dip[i] = -k_d / mu * 
            (dp_dd - hydraulic_props.fluid_density * 9.81);
        
        mp_states[i].fluid_flux_strike = hydraulic_field.fluid_flux_strike[i];
        mp_states[i].fluid_flux_dip = hydraulic_field.fluid_flux_dip[i];
    }
}

void FaultCohesiveTHM::updateHeatTransfer(double dt) {
    // Compute temperature evolution from frictional heating and diffusion
    
    size_t n = mp_states.size();
    if (n == 0) return;
    
    // Volumetric heat capacity
    double rho = 2700.0;  // Rock density kg/m³
    double cp = thermal_props.specific_heat_solid;
    double rho_cp = rho * cp;
    
    for (size_t i = 0; i < n; ++i) {
        // Adiabatic temperature rise from frictional heating
        double dT = FaultThermal::adiabaticTemperatureRise(
            mp_states[i].shear_stress,
            mp_states[i].slip_rate * dt,  // slip in this time step
            hydraulic_props.fault_thickness,
            rho_cp);
        
        // Update temperature
        if (i < thermal_field.temperature.size()) {
            thermal_field.temperature[i] += dT;
            mp_states[i].temperature = thermal_field.temperature[i];
            mp_states[i].temperature_rise = 
                thermal_field.temperature[i] - thermal_props.reference_temperature;
        }
    }
}

void FaultCohesiveTHM::updateFaultProperties(double dt) {
    if (!property_updater) return;
    
    for (size_t i = 0; i < mp_states.size(); ++i) {
        // Update mechanical state from dynamic state
        mp_states[i].slip_total = dyn_states[i].slip[0];  // Total slip magnitude
        mp_states[i].slip_rate = dyn_states[i].slip_rate[0];
        mp_states[i].shear_stress = dyn_states[i].shear_stress;
        
        // Update permeability using callback
        FaultPermeabilityTensor new_perm = property_updater->updatePermeability(
            mp_states[i], hydraulic_props, dt);
        
        if (i < hydraulic_field.permeability.size()) {
            hydraulic_field.permeability[i] = new_perm;
        }
        mp_states[i].current_permeability = new_perm;
        
        // Update porosity using callback
        double new_phi = property_updater->updatePorosity(
            mp_states[i], hydraulic_props, dt);
        
        if (i < hydraulic_field.porosity.size()) {
            hydraulic_field.porosity[i] = new_phi;
        }
        mp_states[i].current_porosity = new_phi;
        
        // Update friction using callback (e.g., thermal weakening)
        double base_friction = friction_model->getFriction();
        mp_states[i].current_friction = property_updater->updateFriction(
            mp_states[i], base_friction);
    }
}

void FaultCohesiveTHM::computeFrictionalHeating() {
    for (size_t i = 0; i < mp_states.size(); ++i) {
        // Q = τ * V (heat generation rate)
        double heating = FaultThermal::frictionalHeatingRate(
            mp_states[i].shear_stress,
            mp_states[i].slip_rate);
        
        if (i < thermal_field.frictional_heating.size()) {
            thermal_field.frictional_heating[i] = heating;
        }
        mp_states[i].frictional_heating = heating;
    }
}

void FaultCohesiveTHM::addPoroelasticResidual(const double* /*solution*/, double* /*residual*/) {
    // Add pore pressure contribution to fault traction balance
    // R_fault = T_cohesive - T_contact - α * p * n
    //
    // This modifies the normal traction by the effective stress principle
    // The actual implementation depends on the specific DOF ordering
    // and finite element discretization
}

void FaultCohesiveTHM::addThermalResidual(const double* /*solution*/, double* /*residual*/) {
    // Add thermal stress contribution
    // R_thermal = ∫ α_T * ΔT * K * ε dV
    //
    // For faults, this appears as:
    // - Modified fault strength (thermal weakening)
    // - Thermal expansion effects on fault opening
}

void FaultCohesiveTHM::addFlowResidual(const double* /*solution*/, double* /*residual*/) {
    // Add fluid flow residual for fault
    // R_flow = S * ∂p/∂t + ∇·q - Q
    //
    // For faults, we have 2D flow on fault surface plus
    // cross-fault flow acting as source/sink
}

// =============================================================================
// Poroelastic Helper Functions
// =============================================================================

namespace FaultPoroelastic {

double undrainedPressureChange(double shear_strain, double dilatancy_coeff,
                               double undrained_bulk) {
    // Δp = -B * K_u * ε_v
    // where ε_v = -dilatancy_coeff * γ (shear-induced dilatancy)
    double eps_v = -dilatancy_coeff * std::abs(shear_strain);
    return -undrained_bulk * eps_v;
}

double diffusivePressure(double p0, double p_boundary, 
                        double hydraulic_diffusivity,
                        double distance, double time) {
    if (time <= 0 || hydraulic_diffusivity <= 0) return p0;
    
    // Error function solution: p = p0 + (p_b - p0) * erfc(x / sqrt(4 * D * t))
    double arg = distance / std::sqrt(4.0 * hydraulic_diffusivity * time);
    return p0 + (p_boundary - p0) * std::erfc(arg);
}

double skemptonCoefficient(double biot, double bulk_drained,
                          double bulk_solid, double bulk_fluid,
                          double porosity) {
    // B = 1 / (1 + φ * K_d / K_f * (1 - K_d/K_s))
    double term1 = porosity * bulk_drained / bulk_fluid;
    double term2 = 1.0 - bulk_drained / bulk_solid;
    return 1.0 / (1.0 + term1 * term2);
}

} // namespace FaultPoroelastic

// =============================================================================
// Thermal Functions
// =============================================================================

namespace FaultThermal {

double thermalPressurization(double slip_rate, double shear_stress,
                            double fault_width, double thermal_diff,
                            double hydraulic_diff, double pressurization_coeff,
                            double dt) {
    if (slip_rate < 1e-15 || dt <= 0) return 0.0;
    
    // Heat generation rate
    double Q = frictionalHeatingRate(shear_stress, slip_rate) / fault_width;
    
    // Pore pressure change from thermal pressurization
    // dp/dt ≈ (λ_f / β) * (Q / ρ_c) * sqrt(π * α_th / α_hy)
    // Simplified model assuming slip rate >> diffusion rate
    double ratio = thermal_diff / hydraulic_diff;
    if (ratio > 0) {
        double dp = pressurization_coeff * Q / (2.7e6) * std::sqrt(M_PI * ratio) * dt;
        return dp;
    }
    return 0.0;
}

double adiabaticTemperatureRise(double shear_stress, double slip,
                               double fault_width, double heat_capacity) {
    if (fault_width <= 0 || heat_capacity <= 0) return 0.0;
    
    // ΔT = τ * δ / (ρ * c * w)
    // Energy per unit area = τ * δ
    // Distributed over fault width w
    return shear_stress * std::abs(slip) / (heat_capacity * fault_width);
}

double flashTemperature(double slip_rate, double shear_stress,
                       double contact_diameter, double thermal_diff,
                       double heat_capacity, double ambient_temp) {
    if (slip_rate < 1e-15) return ambient_temp;
    
    // Flash temperature at asperity contacts (Rice 2006)
    // T_flash = T_0 + (τ * V * D) / (2 * k * sqrt(π * α * D/V))
    double k = heat_capacity * thermal_diff;  // Thermal conductivity approx
    double contact_time = contact_diameter / slip_rate;
    double diffusion_length = std::sqrt(thermal_diff * contact_time);
    
    if (diffusion_length > 0) {
        double dT = shear_stress * slip_rate * contact_diameter / 
                   (2.0 * k * std::sqrt(M_PI) * diffusion_length);
        return ambient_temp + dT;
    }
    return ambient_temp;
}

double thermalWeakeningFriction(double base_friction, double temperature,
                               double weakening_temp, double residual_friction) {
    if (temperature < weakening_temp) {
        return base_friction;
    }
    // Exponential weakening above threshold
    double dT = temperature - weakening_temp;
    double decay_scale = 100.0;  // Temperature scale for decay
    return residual_friction + (base_friction - residual_friction) * 
           std::exp(-dT / decay_scale);
}

} // namespace FaultThermal

// =============================================================================
// Factory Functions for Multi-Physics
// =============================================================================

std::unique_ptr<FaultCohesiveTHM> createFaultCohesiveTHM(
    const std::map<std::string, std::string>& config) {
    
    auto fault = std::make_unique<FaultCohesiveTHM>();
    fault->configure(config);
    return fault;
}

FaultCouplingMode parseFaultCouplingMode(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    
    if (s == "mechanical" || s == "mechanical_only") {
        return FaultCouplingMode::MECHANICAL_ONLY;
    } else if (s == "poroelastic" || s == "hm") {
        return FaultCouplingMode::POROELASTIC;
    } else if (s == "thermoelastic" || s == "tm") {
        return FaultCouplingMode::THERMOELASTIC;
    } else if (s == "thm" || s == "thermo_hydro_mechanical") {
        return FaultCouplingMode::THERMO_HYDRO_MECHANICAL;
    } else if (s == "thmc") {
        return FaultCouplingMode::THMC;
    }
    return FaultCouplingMode::MECHANICAL_ONLY;
}

PermeabilityUpdateModel parsePermeabilityUpdateModel(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    
    if (s == "constant" || s == "none") {
        return PermeabilityUpdateModel::CONSTANT;
    } else if (s == "slip" || s == "slip_dependent") {
        return PermeabilityUpdateModel::SLIP_DEPENDENT;
    } else if (s == "dilation" || s == "dilation_dependent") {
        return PermeabilityUpdateModel::DILATION_DEPENDENT;
    } else if (s == "stress" || s == "stress_dependent") {
        return PermeabilityUpdateModel::STRESS_DEPENDENT;
    } else if (s == "shear" || s == "shear_enhanced") {
        return PermeabilityUpdateModel::SHEAR_ENHANCED;
    } else if (s == "damage" || s == "damage_dependent") {
        return PermeabilityUpdateModel::DAMAGE_DEPENDENT;
    } else if (s == "combined" || s == "all") {
        return PermeabilityUpdateModel::COMBINED;
    }
    return PermeabilityUpdateModel::CONSTANT;
}

PorosityUpdateModel parsePorosityUpdateModel(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    
    if (s == "constant" || s == "none") {
        return PorosityUpdateModel::CONSTANT;
    } else if (s == "compaction") {
        return PorosityUpdateModel::COMPACTION;
    } else if (s == "dilation") {
        return PorosityUpdateModel::DILATION;
    } else if (s == "chemical") {
        return PorosityUpdateModel::CHEMICAL;
    } else if (s == "combined" || s == "all") {
        return PorosityUpdateModel::COMBINED;
    }
    return PorosityUpdateModel::CONSTANT;
}

} // namespace FSRM
