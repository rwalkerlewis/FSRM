#include "FaultModel.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <iostream>

namespace FSRM {

// =============================================================================
// Utility functions
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

// =============================================================================
// FaultStressState Implementation
// =============================================================================

double FaultStressState::getMagnitude() const {
    if (moment < 1e10) return -999.0;  // Too small
    return (2.0/3.0) * (std::log10(moment) - 9.1);  // Kanamori (1977)
}

// =============================================================================
// CoulombParameters Implementation
// =============================================================================

void CoulombParameters::configure(const std::map<std::string, std::string>& config) {
    static_friction = parseDouble(config, "static_friction", 0.6);
    dynamic_friction = parseDouble(config, "dynamic_friction", 0.4);
    cohesion = parseDouble(config, "cohesion", 1e6);
    slip_weakening_Dc = parseDouble(config, "slip_weakening_dc", 0.01);
}

// =============================================================================
// RateStateParameters Implementation
// =============================================================================

void RateStateParameters::configure(const std::map<std::string, std::string>& config) {
    a = parseDouble(config, "rate_state_a", 0.010);
    b = parseDouble(config, "rate_state_b", 0.015);
    Dc = parseDouble(config, "rate_state_dc", 1e-4);
    V0 = parseDouble(config, "rate_state_v0", 1e-6);
    f0 = parseDouble(config, "rate_state_f0", 0.6);
    theta0 = parseDouble(config, "rate_state_theta0", Dc / V0);
}

// =============================================================================
// FaultGeometry Implementation
// =============================================================================

std::array<double, 3> FaultGeometry::getStrikeVector() const {
    // Strike vector in ENU coordinates
    return {std::sin(strike), std::cos(strike), 0.0};
}

std::array<double, 3> FaultGeometry::getDipVector() const {
    // Down-dip vector
    double cos_s = std::cos(strike);
    double sin_s = std::sin(strike);
    double cos_d = std::cos(dip);
    double sin_d = std::sin(dip);
    
    return {-cos_s * cos_d, sin_s * cos_d, -sin_d};
}

std::array<double, 3> FaultGeometry::getNormalVector() const {
    // Outward normal (pointing into hanging wall for dipping fault)
    double cos_s = std::cos(strike);
    double sin_s = std::sin(strike);
    double cos_d = std::cos(dip);
    double sin_d = std::sin(dip);
    
    return {cos_s * sin_d, -sin_s * sin_d, cos_d};
}

void FaultGeometry::resolveStress(double sxx, double syy, double szz,
                                  double sxy, double sxz, double syz,
                                  double& sigma_n, double& tau) const {
    // Get fault normal
    auto n = getNormalVector();
    
    // Traction vector t_i = σ_ij * n_j
    double tx = sxx * n[0] + sxy * n[1] + sxz * n[2];
    double ty = sxy * n[0] + syy * n[1] + syz * n[2];
    double tz = sxz * n[0] + syz * n[1] + szz * n[2];
    
    // Normal stress (positive = compression)
    sigma_n = tx * n[0] + ty * n[1] + tz * n[2];
    
    // Shear traction (tangent to fault)
    double tau_x = tx - sigma_n * n[0];
    double tau_y = ty - sigma_n * n[1];
    double tau_z = tz - sigma_n * n[2];
    
    // Shear stress magnitude
    tau = std::sqrt(tau_x * tau_x + tau_y * tau_y + tau_z * tau_z);
}

void FaultGeometry::configure(const std::map<std::string, std::string>& config) {
    x = parseDouble(config, "x", 0.0);
    y = parseDouble(config, "y", 0.0);
    z = parseDouble(config, "z", 0.0);
    
    // Convert from degrees if specified
    strike = parseDouble(config, "strike", 0.0);
    if (strike > 2.0 * M_PI) strike *= M_PI / 180.0;  // Convert from degrees
    
    dip = parseDouble(config, "dip", 60.0);
    if (dip > M_PI) dip *= M_PI / 180.0;
    
    rake = parseDouble(config, "rake", 0.0);
    if (std::abs(rake) > M_PI) rake *= M_PI / 180.0;
    
    length = parseDouble(config, "length", 1000.0);
    width = parseDouble(config, "width", 500.0);
}

// =============================================================================
// CoulombFriction Implementation
// =============================================================================

CoulombFriction::CoulombFriction()
    : FrictionModelBase(FrictionLaw::COULOMB),
      current_slip(0.0) {}

double CoulombFriction::getFriction(double slip_rate, double state_var,
                                    double sigma_n_eff) const {
    // Use state_var as cumulative slip for slip-weakening
    double slip = state_var;
    
    if (slip >= params.slip_weakening_Dc) {
        // Fully weakened - dynamic friction
        return params.dynamic_friction;
    } else if (slip <= 0.0) {
        // Not slipping - static friction
        return params.static_friction;
    } else {
        // Linear slip-weakening
        double f = params.static_friction - 
                  (params.static_friction - params.dynamic_friction) * 
                  slip / params.slip_weakening_Dc;
        return f;
    }
}

double CoulombFriction::getStateEvolutionRate(double slip_rate, double state_var) const {
    // For Coulomb with slip-weakening, state = cumulative slip
    return slip_rate;
}

void CoulombFriction::configure(const std::map<std::string, std::string>& config) {
    params.configure(config);
}

// =============================================================================
// RateStateFriction Implementation
// =============================================================================

RateStateFriction::RateStateFriction(EvolutionLaw law)
    : FrictionModelBase(law == EvolutionLaw::AGING ? 
                       FrictionLaw::RATE_STATE_AGING : FrictionLaw::RATE_STATE_SLIP),
      evolution_law(law) {}

double RateStateFriction::getFriction(double slip_rate, double state_var,
                                      double sigma_n_eff) const {
    // Rate-and-state friction: f = f0 + a*ln(V/V0) + b*ln(θ*V0/Dc)
    
    // Regularize slip rate for zero/negative values
    double V = std::max(slip_rate, 1e-20);
    double theta = std::max(state_var, 1e-10);
    
    double f = params.f0 + 
              params.a * std::log(V / params.V0) + 
              params.b * std::log(theta * params.V0 / params.Dc);
    
    // Ensure positive friction
    return std::max(0.0, f);
}

double RateStateFriction::getStateEvolutionRate(double slip_rate, double state_var) const {
    // Regularize
    double V = std::max(std::abs(slip_rate), 1e-20);
    double theta = std::max(state_var, 1e-10);
    
    switch (evolution_law) {
        case EvolutionLaw::AGING:
            // Aging law: dθ/dt = 1 - Vθ/Dc
            return 1.0 - V * theta / params.Dc;
            
        case EvolutionLaw::SLIP:
            // Slip law: dθ/dt = -Vθ/Dc * ln(Vθ/Dc)
            return -V * theta / params.Dc * std::log(V * theta / params.Dc);
            
        default:
            return 1.0 - V * theta / params.Dc;
    }
}

void RateStateFriction::configure(const std::map<std::string, std::string>& config) {
    params.configure(config);
    
    std::string law = parseString(config, "rate_state_law", "AGING");
    std::transform(law.begin(), law.end(), law.begin(), ::toupper);
    if (law == "SLIP") {
        evolution_law = EvolutionLaw::SLIP;
        friction_law = FrictionLaw::RATE_STATE_SLIP;
    } else {
        evolution_law = EvolutionLaw::AGING;
        friction_law = FrictionLaw::RATE_STATE_AGING;
    }
}

double RateStateFriction::getSteadyStateFriction(double slip_rate) const {
    // At steady state, θ = Dc/V
    double V = std::max(slip_rate, 1e-20);
    return params.f0 + (params.a - params.b) * std::log(V / params.V0);
}

double RateStateFriction::getCriticalStiffness(double sigma_n_eff) const {
    // For spring-slider stability: k_c = (b-a)*σ_n/Dc
    return (params.b - params.a) * sigma_n_eff / params.Dc;
}

// =============================================================================
// SeismicFaultModel Implementation
// =============================================================================

SeismicFaultModel::SeismicFaultModel()
    : name("fault"),
      current_slip(0.0),
      current_state_var(1e6),
      current_slip_rate(1e-12),
      b_value(1.0),
      nucleation_size(1.0),
      seismic_slip_rate(1e-3) {
    
    // Default to Coulomb friction
    friction = std::make_unique<CoulombFriction>();
}

void SeismicFaultModel::configure(const std::map<std::string, std::string>& config) {
    name = parseString(config, "name", "fault");
    
    geometry.configure(config);
    
    b_value = parseDouble(config, "b_value", 1.0);
    nucleation_size = parseDouble(config, "nucleation_size", 1.0);
    seismic_slip_rate = parseDouble(config, "seismic_slip_rate", 1e-3);
    
    // Set friction model
    std::string law_str = parseString(config, "friction_law", "COULOMB");
    FrictionLaw law = parseFrictionLaw(law_str);
    friction = createFrictionModel(law, config);
}

void SeismicFaultModel::setGeometry(const FaultGeometry& geom) {
    geometry = geom;
}

void SeismicFaultModel::setGeometry(double x, double y, double z,
                            double strike, double dip, double length, double width) {
    geometry.x = x;
    geometry.y = y;
    geometry.z = z;
    geometry.strike = strike;
    geometry.dip = dip;
    geometry.length = length;
    geometry.width = width;
}

void SeismicFaultModel::setFrictionModel(std::unique_ptr<FrictionModelBase> model) {
    friction = std::move(model);
}

void SeismicFaultModel::setFrictionLaw(FrictionLaw law) {
    friction = createFrictionModel(law);
}

FaultStressState SeismicFaultModel::computeStressState(double sxx, double syy, double szz,
                                                double sxy, double sxz, double syz,
                                                double pore_pressure) const {
    FaultStressState state;
    
    // Resolve stresses onto fault plane
    geometry.resolveStress(sxx, syy, szz, sxy, sxz, syz,
                          state.sigma_n, state.tau);
    
    // Effective normal stress
    state.sigma_n_eff = state.sigma_n - pore_pressure;
    
    // Ensure compression (positive)
    if (state.sigma_n_eff < 0.0) {
        state.sigma_n_eff = 0.0;  // Tensile opening
    }
    
    // Get friction coefficient
    double mu = friction->getFriction(current_slip_rate, current_state_var,
                                      state.sigma_n_eff);
    
    // Frictional strength
    if (auto coulomb = dynamic_cast<CoulombFriction*>(friction.get())) {
        state.tau_strength = coulomb->params.cohesion + mu * state.sigma_n_eff;
    } else {
        state.tau_strength = mu * state.sigma_n_eff;
    }
    
    // Coulomb Failure Function
    state.CFF = state.tau - state.tau_strength;
    
    // Current slip state
    state.slip = current_slip;
    state.slip_rate = current_slip_rate;
    state.state_variable = current_state_var;
    
    state.is_slipping = (state.CFF >= 0.0);
    state.is_seismic = (current_slip_rate >= seismic_slip_rate);
    
    // Seismic moment
    double G = 30e9;  // Approximate shear modulus (would need to be passed)
    state.moment = G * geometry.area() * current_slip;
    state.stress_drop = state.is_slipping ? state.tau - state.tau_strength : 0.0;
    
    return state;
}

FaultStressState SeismicFaultModel::updateSlip(const FaultStressState& current_state,
                                        double shear_modulus, double dt) {
    FaultStressState new_state = current_state;
    
    if (current_state.CFF <= 0.0) {
        // Elastic: no slip, state evolves
        current_slip_rate = 0.0;
        
        // State variable evolution (healing)
        double dtheta = friction->getStateEvolutionRate(0.0, current_state_var) * dt;
        current_state_var = std::max(1e-10, current_state_var + dtheta);
        
        new_state.slip_rate = 0.0;
        new_state.state_variable = current_state_var;
        new_state.is_slipping = false;
        return new_state;
    }
    
    // Slip is occurring
    new_state.is_slipping = true;
    
    // Solve for slip rate (simplified - full solution requires stiffness)
    // Approximation: excess stress drives slip
    double stress_excess = current_state.CFF;
    double radiation_damping = 0.5 * shear_modulus / 
                               std::sqrt(shear_modulus / 2500.0);  // Rough approximation
    
    current_slip_rate = std::max(1e-12, stress_excess / radiation_damping);
    
    // Update slip
    double dslip = current_slip_rate * dt;
    current_slip += dslip;
    
    // Update state variable
    double dtheta = friction->getStateEvolutionRate(current_slip_rate, current_state_var) * dt;
    current_state_var = std::max(1e-10, current_state_var + dtheta);
    
    // Update new state
    new_state.slip = current_slip;
    new_state.slip_rate = current_slip_rate;
    new_state.state_variable = current_state_var;
    new_state.is_seismic = (current_slip_rate >= seismic_slip_rate);
    new_state.moment = shear_modulus * geometry.area() * current_slip;
    
    return new_state;
}

bool SeismicFaultModel::checkNucleation(const FaultStressState& state) const {
    // Check if stress exceeds strength and patch is large enough
    if (state.CFF < 0.0) return false;
    
    // Simplified nucleation criterion
    // Full implementation would check patch size vs nucleation length
    return (state.tau / state.tau_strength) > 1.0;
}

SeismicEvent SeismicFaultModel::computeEvent(const FaultStressState& initial_state,
                                     const FaultStressState& final_state,
                                     double shear_modulus,
                                     double current_time) const {
    SeismicEvent event;
    
    event.time = current_time;
    event.slip = final_state.slip - initial_state.slip;
    event.area = geometry.area();
    event.moment = shear_modulus * event.area * event.slip;
    event.magnitude = (2.0/3.0) * (std::log10(event.moment) - 9.1);
    event.stress_drop = initial_state.tau - final_state.tau;
    
    // Rupture duration (simplified: circular crack)
    double rupture_velocity = 0.9 * std::sqrt(shear_modulus / 2500.0);
    double radius = std::sqrt(event.area / M_PI);
    event.duration = radius / rupture_velocity;
    
    // Location (fault center)
    event.x = geometry.x;
    event.y = geometry.y;
    event.z = geometry.z;
    
    // Hypocenter (approximate as center)
    event.hypo_x = geometry.x;
    event.hypo_y = geometry.y;
    event.hypo_z = geometry.z;
    
    event.type = final_state.is_seismic ? "earthquake" : "slow_slip";
    
    return event;
}

void SeismicFaultModel::computeCoulombStressChange(double slip, double receiver_x,
                                            double receiver_y, double receiver_z,
                                            double shear_modulus,
                                            double& dCFF) const {
    // Simplified Coulomb stress change calculation
    // Full implementation would use Okada (1992) equations
    
    // Distance from fault center to receiver
    double dx = receiver_x - geometry.x;
    double dy = receiver_y - geometry.y;
    double dz = receiver_z - geometry.z;
    double r = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    if (r < 1.0) r = 1.0;  // Avoid singularity
    
    // Simplified 1/r² decay
    double stress_scale = shear_modulus * slip / (geometry.length);
    dCFF = stress_scale * geometry.area() / (4.0 * M_PI * r * r);
    
    // Azimuthal dependence (simplified)
    double azimuth = std::atan2(dy, dx) - geometry.strike;
    dCFF *= (1.0 + 0.5 * std::cos(2.0 * azimuth));
}

void SeismicFaultModel::addEvent(const SeismicEvent& event) {
    events.push_back(event);
}

void SeismicFaultModel::clearCatalog() {
    events.clear();
}

double SeismicFaultModel::getSeismicityRate(double stressing_rate, double sigma_n_eff,
                                     double temperature) const {
    // Dieterich (1994) seismicity rate model
    // R/r = 1 / (1 - (Δτ/(A*σ))*γ)
    // where γ evolves with time and stressing
    
    auto rs = dynamic_cast<RateStateFriction*>(friction.get());
    if (!rs) {
        return 1.0;  // Background rate
    }
    
    double A = rs->params.a;
    double Asigma = A * sigma_n_eff;
    
    // Steady-state rate
    return stressing_rate / Asigma;
}

// =============================================================================
// FaultNetwork Implementation
// =============================================================================

FaultNetwork::FaultNetwork()
    : use_stress_transfer(false) {}

void FaultNetwork::addFault(std::unique_ptr<SeismicFaultModel> fault) {
    faults.push_back(std::move(fault));
}

void FaultNetwork::addFault(const std::map<std::string, std::string>& config) {
    auto fault = std::make_unique<SeismicFaultModel>();
    fault->configure(config);
    faults.push_back(std::move(fault));
}

void FaultNetwork::configure(const std::vector<std::map<std::string, std::string>>& fault_configs) {
    faults.clear();
    for (const auto& config : fault_configs) {
        addFault(config);
    }
}

void FaultNetwork::update(double sxx, double syy, double szz,
                         double sxy, double sxz, double syz,
                         double pore_pressure, double shear_modulus, double dt) {
    for (auto& fault : faults) {
        auto state = fault->computeStressState(sxx, syy, szz, sxy, sxz, syz, pore_pressure);
        auto new_state = fault->updateSlip(state, shear_modulus, dt);
        
        // Record event if significant slip occurred
        if (new_state.slip - state.slip > 1e-6 && new_state.is_seismic) {
            auto event = fault->computeEvent(state, new_state, shear_modulus, 0.0);
            fault->addEvent(event);
        }
    }
}

void FaultNetwork::updateWithCascade(double sxx, double syy, double szz,
                                     double sxy, double sxz, double syz,
                                     double pore_pressure, double shear_modulus, double dt) {
    // First pass: compute stress states
    std::vector<FaultStressState> states(faults.size());
    for (size_t i = 0; i < faults.size(); ++i) {
        states[i] = faults[i]->computeStressState(sxx, syy, szz, sxy, sxz, syz, pore_pressure);
    }
    
    // Update with stress transfer
    for (size_t i = 0; i < faults.size(); ++i) {
        auto new_state = faults[i]->updateSlip(states[i], shear_modulus, dt);
        
        // If slip occurred, compute stress transfer to other faults
        double dslip = new_state.slip - states[i].slip;
        if (dslip > 1e-8 && use_stress_transfer) {
            for (size_t j = 0; j < faults.size(); ++j) {
                if (i == j) continue;
                
                double dCFF;
                faults[i]->computeCoulombStressChange(dslip,
                    faults[j]->getGeometry().x,
                    faults[j]->getGeometry().y,
                    faults[j]->getGeometry().z,
                    shear_modulus, dCFF);
                
                // Modify receiver fault's stress state
                // (simplified - would need to track accumulated stress)
            }
        }
        
        // Record event
        if (dslip > 1e-6 && new_state.is_seismic) {
            auto event = faults[i]->computeEvent(states[i], new_state, shear_modulus, 0.0);
            faults[i]->addEvent(event);
        }
    }
}

SeismicFaultModel* FaultNetwork::getFault(size_t index) {
    if (index < faults.size()) {
        return faults[index].get();
    }
    return nullptr;
}

SeismicFaultModel* FaultNetwork::getFault(const std::string& name) {
    for (auto& fault : faults) {
        if (fault->name == name) {
            return fault.get();
        }
    }
    return nullptr;
}

std::vector<SeismicEvent> FaultNetwork::getAllEvents() const {
    std::vector<SeismicEvent> all_events;
    for (const auto& fault : faults) {
        const auto& fault_events = fault->getEventCatalog();
        all_events.insert(all_events.end(), fault_events.begin(), fault_events.end());
    }
    
    // Sort by time
    std::sort(all_events.begin(), all_events.end(),
              [](const SeismicEvent& a, const SeismicEvent& b) {
                  return a.time < b.time;
              });
    
    return all_events;
}

double FaultNetwork::getTotalMomentRate() const {
    double total_moment = 0.0;
    for (const auto& fault : faults) {
        for (const auto& event : fault->getEventCatalog()) {
            total_moment += event.moment;
        }
    }
    // Would need time span to compute rate
    return total_moment;
}

double FaultNetwork::getMagnitudeOfCompleteness() const {
    auto events = getAllEvents();
    return SeismicityAnalysis::estimateMc(events);
}

// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<FrictionModelBase> createFrictionModel(
    FrictionLaw law,
    const std::map<std::string, std::string>& config) {
    
    std::unique_ptr<FrictionModelBase> model;
    
    switch (law) {
        case FrictionLaw::COULOMB:
        case FrictionLaw::SLIP_WEAKENING:
            model = std::make_unique<CoulombFriction>();
            break;
        case FrictionLaw::RATE_STATE_AGING:
            model = std::make_unique<RateStateFriction>(RateStateFriction::EvolutionLaw::AGING);
            break;
        case FrictionLaw::RATE_STATE_SLIP:
            model = std::make_unique<RateStateFriction>(RateStateFriction::EvolutionLaw::SLIP);
            break;
        default:
            model = std::make_unique<CoulombFriction>();
    }
    
    if (!config.empty()) {
        model->configure(config);
    }
    
    return model;
}

FrictionLaw parseFrictionLaw(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    
    if (s == "COULOMB" || s == "MOHR_COULOMB") {
        return FrictionLaw::COULOMB;
    } else if (s == "RATE_STATE" || s == "RATE_STATE_AGING" || s == "RS_AGING") {
        return FrictionLaw::RATE_STATE_AGING;
    } else if (s == "RATE_STATE_SLIP" || s == "RS_SLIP") {
        return FrictionLaw::RATE_STATE_SLIP;
    } else if (s == "SLIP_WEAKENING" || s == "SW") {
        return FrictionLaw::SLIP_WEAKENING;
    } else if (s == "FLASH_HEATING") {
        return FrictionLaw::FLASH_HEATING;
    } else if (s == "THERMAL_PRESSURIZATION") {
        return FrictionLaw::THERMAL_PRESSURIZATION;
    }
    
    return FrictionLaw::COULOMB;
}

// =============================================================================
// SeismicityAnalysis Functions
// =============================================================================

namespace SeismicityAnalysis {

double estimateBValue(const std::vector<SeismicEvent>& events, double Mc) {
    if (events.empty()) return 1.0;
    
    // Filter events above Mc
    std::vector<double> mags;
    for (const auto& e : events) {
        if (e.magnitude >= Mc) {
            mags.push_back(e.magnitude);
        }
    }
    
    if (mags.empty()) return 1.0;
    
    // Maximum likelihood estimate: b = log10(e) / (M_mean - M_c)
    double sum = std::accumulate(mags.begin(), mags.end(), 0.0);
    double M_mean = sum / mags.size();
    
    if (M_mean <= Mc) return 1.0;
    
    return std::log10(M_E) / (M_mean - Mc);
}

double estimateMc(const std::vector<SeismicEvent>& events) {
    if (events.empty()) return 0.0;
    
    // Maximum curvature method (MAXC)
    // Find magnitude with maximum frequency in binned histogram
    
    // Create magnitude bins (0.1 magnitude wide)
    std::map<int, int> bins;
    for (const auto& e : events) {
        int bin = static_cast<int>(e.magnitude * 10);
        bins[bin]++;
    }
    
    if (bins.empty()) return 0.0;
    
    // Find maximum
    int max_count = 0;
    int max_bin = 0;
    for (const auto& [bin, count] : bins) {
        if (count > max_count) {
            max_count = count;
            max_bin = bin;
        }
    }
    
    // Mc = bin center + 0.2 (typical correction)
    return max_bin / 10.0 + 0.05 + 0.2;
}

double omoriRate(double t, double K, double c, double p) {
    // Modified Omori law: R(t) = K / (t + c)^p
    return K / std::pow(t + c, p);
}

ETASParams fitETAS(const std::vector<SeismicEvent>& events) {
    // Placeholder - full ETAS fitting requires MLE optimization
    ETASParams params;
    params.mu = 1.0;
    params.K = 0.1;
    params.alpha = 1.0;
    params.c = 0.01;
    params.p = 1.1;
    
    if (events.size() < 10) return params;
    
    // Simple estimates from data
    params.mu = static_cast<double>(events.size()) / 
                (events.back().time - events.front().time);
    
    return params;
}

double rateChange(double stress_change, double stressing_rate,
                  double Asigma, double temperature) {
    // Dieterich (1994) rate change model
    // R/r = 1 / (1 - γ*Δτ/(A*σ))
    // where γ = 1/(stress rate * t_a)
    
    if (std::abs(stressing_rate) < 1e-20) return 1.0;
    if (std::abs(Asigma) < 1e-10) return 1.0;
    
    double ta = 100.0;  // Aftershock duration (simplified)
    double gamma = 1.0 / (stressing_rate * ta);
    
    double denominator = 1.0 - gamma * stress_change / Asigma;
    
    if (denominator <= 0.0) {
        return 1000.0;  // Large rate increase
    }
    
    return 1.0 / denominator;
}

double coulombStressChange(double sxx, double syy, double sxy,
                           double receiver_strike, double receiver_dip,
                           double receiver_rake, double mu) {
    // Resolve 2D stress onto receiver fault orientation
    double phi = receiver_strike;
    
    // Normal stress on receiver
    double sigma_n = 0.5 * (sxx + syy) + 
                    0.5 * (sxx - syy) * std::cos(2.0 * phi) +
                    sxy * std::sin(2.0 * phi);
    
    // Shear stress on receiver (in rake direction)
    double tau = -0.5 * (sxx - syy) * std::sin(2.0 * phi) +
                 sxy * std::cos(2.0 * phi);
    
    // Coulomb stress change
    return tau - mu * sigma_n;
}

} // namespace SeismicityAnalysis

} // namespace FSRM
