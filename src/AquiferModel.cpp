#include "AquiferModel.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace FSRM {

// ============================================================================
// DimensionlessInflux Implementation
// ============================================================================

double DimensionlessInflux::pD_infinite(double tD) {
    // Van Everdingen-Hurst solution for infinite aquifer
    // Using polynomial approximation (accurate to <0.1% for tD > 0.01)
    
    if (tD <= 0.0) return 0.0;
    
    if (tD < 0.01) {
        // Early time: pD ≈ 2*sqrt(tD/π)
        return 2.0 * std::sqrt(tD / M_PI);
    } else if (tD < 200.0) {
        // Intermediate time: polynomial fit
        double sqrtTD = std::sqrt(tD);
        double logTD = std::log(tD);
        
        // Chatas approximation
        double pD = (2.0 * sqrtTD / std::sqrt(M_PI)) 
                  + (0.5 * logTD + 0.80907) * (1.0 - std::exp(-tD))
                  - 0.049 * std::exp(-tD);
        
        return pD;
    } else {
        // Late time: pD ≈ 0.5*ln(tD) + 0.80907
        return 0.5 * std::log(tD) + 0.80907;
    }
}

double DimensionlessInflux::pD_finite(double tD, double reD) {
    // Finite aquifer solution
    if (tD <= 0.0) return 0.0;
    
    double pD_inf = pD_infinite(tD);
    
    if (reD > 100.0) {
        // Nearly infinite aquifer
        return pD_inf;
    }
    
    // Boundary effect correction
    double tD_boundary = reD * reD / 4.0;
    
    if (tD < tD_boundary) {
        // Infinite-acting period
        return pD_inf;
    }
    
    // Pseudo-steady state correction
    double delta = (tD - tD_boundary) / (reD * reD);
    double pD_pss = pD_inf - delta * (1.0 - 1.0 / (reD * reD));
    
    return std::max(0.0, pD_pss);
}

double DimensionlessInflux::dpD_dtD(double tD) {
    // Derivative for Carter-Tracy superposition
    if (tD <= 0.0) return 0.0;
    
    if (tD < 0.01) {
        return 1.0 / std::sqrt(M_PI * tD);
    } else {
        // Numerical derivative
        double eps = tD * 1e-6;
        return (pD_infinite(tD + eps) - pD_infinite(tD - eps)) / (2.0 * eps);
    }
}

double DimensionlessInflux::pD_linear(double tD) {
    // Linear aquifer dimensionless influx
    if (tD <= 0.0) return 0.0;
    return 2.0 * std::sqrt(tD / M_PI);
}

// ============================================================================
// AquiferModelBase Implementation
// ============================================================================

AquiferModelBase::AquiferModelBase(int id, AquiferType type)
    : aquifer_id(id), aquifer_type(type),
      initial_pressure(1.0e7), current_pressure(1.0e7),
      cumulative_influx(0.0), time(0.0) {}

void AquiferModelBase::addConnection(const AquiferConnection& conn) {
    connections.push_back(conn);
}

// ============================================================================
// CarterTracyAquifer Implementation
// ============================================================================

CarterTracyAquifer::CarterTracyAquifer(int id)
    : AquiferModelBase(id, AquiferType::CARTER_TRACY),
      permeability(1e-13), porosity(0.2),
      total_compressibility(1e-9), water_viscosity(0.001),
      thickness(50.0), encroachment_angle(1.0),
      geometry(AquiferGeometry::RADIAL),
      inner_radius(1000.0), outer_radius(10000.0),
      reD(10.0), influx_constant(0.0), time_constant(0.0),
      num_steps(0) {
    computeConstants();
}

void CarterTracyAquifer::computeConstants() {
    // Influx constant B = 1.119 * φ * ct * re² * h * θ (field units)
    // In SI: B = φ * ct * re² * h * θ
    influx_constant = encroachment_angle * porosity * total_compressibility 
                    * inner_radius * inner_radius * thickness;
    
    // Time constant for dimensionless time
    // tD = k * t / (φ * μ * ct * ro²)
    time_constant = porosity * water_viscosity * total_compressibility 
                  * inner_radius * inner_radius / permeability;
    
    // Dimensionless radius
    reD = outer_radius / inner_radius;
}

void CarterTracyAquifer::setAquiferProperties(double perm_md, double phi,
                                              double ct, double mu,
                                              double h, double theta) {
    permeability = perm_md * 9.869233e-16;  // mD to m²
    porosity = phi;
    total_compressibility = ct;
    water_viscosity = mu;
    thickness = h;
    encroachment_angle = theta;  // Fraction of circle (e.g., 0.25 for quarter)
    
    computeConstants();
}

void CarterTracyAquifer::setAquiferRadius(double r_inner, double r_outer) {
    inner_radius = r_inner;
    outer_radius = r_outer;
    computeConstants();
}

void CarterTracyAquifer::setInitialPressure(double P0) {
    initial_pressure = P0;
    current_pressure = P0;
}

void CarterTracyAquifer::setGeometry(AquiferGeometry geom) {
    geometry = geom;
}

double CarterTracyAquifer::calculateInflux(double reservoir_pressure, double dt) {
    // Store pressure history for superposition
    pressure_history.push_back(reservoir_pressure);
    time_history.push_back(time + dt);
    num_steps++;
    
    return superpositionInflux(reservoir_pressure);
}

double CarterTracyAquifer::superpositionInflux(double current_pressure) {
    // Carter-Tracy method using superposition
    double We = 0.0;  // Cumulative influx
    
    if (num_steps == 0) return 0.0;
    
    // Superposition principle
    for (int j = 0; j < num_steps; ++j) {
        // Dimensionless time from step j to current
        double delta_t = time - (j > 0 ? time_history[j-1] : 0.0);
        double tD = delta_t / time_constant;
        
        // Pressure drop at step j
        double delta_p = (j == 0) 
            ? (initial_pressure - pressure_history[0])
            : (pressure_history[j-1] - pressure_history[j]);
        
        // Dimensionless influx function
        double pD;
        if (geometry == AquiferGeometry::LINEAR) {
            pD = DimensionlessInflux::pD_linear(tD);
        } else {
            pD = (reD > 100.0) 
               ? DimensionlessInflux::pD_infinite(tD)
               : DimensionlessInflux::pD_finite(tD, reD);
        }
        
        // Add contribution
        We += influx_constant * delta_p * pD;
    }
    
    // Calculate instantaneous rate
    double rate = (We - cumulative_influx) / 
                  (time_history.back() - (num_steps > 1 ? time_history[num_steps-2] : 0.0));
    
    return std::max(0.0, rate);
}

void CarterTracyAquifer::updateState(double influx_rate, double dt) {
    cumulative_influx += influx_rate * dt;
    time += dt;
    
    // Update aquifer pressure (for finite aquifer)
    if (reD < 100.0) {
        double pore_volume = M_PI * (outer_radius * outer_radius - inner_radius * inner_radius) 
                           * thickness * porosity * encroachment_angle;
        current_pressure = initial_pressure - cumulative_influx / (pore_volume * total_compressibility);
    }
}

void CarterTracyAquifer::reset() {
    cumulative_influx = 0.0;
    time = 0.0;
    current_pressure = initial_pressure;
    pressure_history.clear();
    time_history.clear();
    num_steps = 0;
}

void CarterTracyAquifer::configure(const std::map<std::string, std::string>& config) {
    // Parse configuration
    if (config.count("permeability")) 
        permeability = std::stod(config.at("permeability")) * 9.869233e-16;
    if (config.count("porosity")) 
        porosity = std::stod(config.at("porosity"));
    if (config.count("compressibility")) 
        total_compressibility = std::stod(config.at("compressibility"));
    if (config.count("viscosity")) 
        water_viscosity = std::stod(config.at("viscosity"));
    if (config.count("thickness")) 
        thickness = std::stod(config.at("thickness"));
    if (config.count("encroachment_angle")) 
        encroachment_angle = std::stod(config.at("encroachment_angle"));
    if (config.count("inner_radius")) 
        inner_radius = std::stod(config.at("inner_radius"));
    if (config.count("outer_radius")) 
        outer_radius = std::stod(config.at("outer_radius"));
    if (config.count("initial_pressure")) 
        setInitialPressure(std::stod(config.at("initial_pressure")));
    if (config.count("geometry"))
        geometry = parseAquiferGeometry(config.at("geometry"));
    
    computeConstants();
}

// ============================================================================
// FetkovichAquifer Implementation
// ============================================================================

FetkovichAquifer::FetkovichAquifer(int id)
    : AquiferModelBase(id, AquiferType::FETKOVICH),
      productivity_index(1e-6), initial_volume(1e9),
      remaining_volume(1e9), total_compressibility(1e-9) {}

double FetkovichAquifer::calculateInflux(double reservoir_pressure, double dt) {
    // Fetkovich model: Q = J * (Pa - Pr) * (Wei - We) / Wei
    // where Pa = current aquifer pressure, Pr = reservoir pressure
    
    // Remaining aquifer pressure based on cumulative production
    double pressure_ratio = remaining_volume / initial_volume;
    double aquifer_pressure = initial_pressure * pressure_ratio;
    
    if (aquifer_pressure <= reservoir_pressure) {
        return 0.0;  // No influx if aquifer pressure <= reservoir pressure
    }
    
    // Influx rate
    double rate = productivity_index * (aquifer_pressure - reservoir_pressure);
    
    // Limit by remaining volume
    double max_rate = remaining_volume / dt;
    return std::min(rate, max_rate);
}

void FetkovichAquifer::updateState(double influx_rate, double dt) {
    double influx_volume = influx_rate * dt;
    
    cumulative_influx += influx_volume;
    remaining_volume = initial_volume - cumulative_influx;
    remaining_volume = std::max(0.0, remaining_volume);
    
    time += dt;
    
    // Update pressure
    current_pressure = initial_pressure * (remaining_volume / initial_volume);
}

void FetkovichAquifer::reset() {
    cumulative_influx = 0.0;
    remaining_volume = initial_volume;
    time = 0.0;
    current_pressure = initial_pressure;
}

void FetkovichAquifer::setProductivityIndex(double J) {
    productivity_index = J;
}

void FetkovichAquifer::setInitialVolume(double V0) {
    initial_volume = V0;
    remaining_volume = V0;
}

void FetkovichAquifer::setCompressibility(double ct) {
    total_compressibility = ct;
}

void FetkovichAquifer::setInitialPressure(double P0) {
    initial_pressure = P0;
    current_pressure = P0;
}

void FetkovichAquifer::configure(const std::map<std::string, std::string>& config) {
    if (config.count("productivity_index"))
        productivity_index = std::stod(config.at("productivity_index"));
    if (config.count("initial_volume"))
        setInitialVolume(std::stod(config.at("initial_volume")));
    if (config.count("compressibility"))
        total_compressibility = std::stod(config.at("compressibility"));
    if (config.count("initial_pressure"))
        setInitialPressure(std::stod(config.at("initial_pressure")));
}

// ============================================================================
// NumericalAquifer Implementation
// ============================================================================

NumericalAquifer::NumericalAquifer(int id)
    : AquiferModelBase(id, AquiferType::NUMERICAL),
      water_density(1000.0), water_viscosity(0.001),
      water_compressibility(4.5e-10) {}

void NumericalAquifer::addCell(const AquiferCell& cell) {
    cells.push_back(cell);
}

void NumericalAquifer::setWaterProperties(double density, double viscosity, 
                                          double compressibility) {
    water_density = density;
    water_viscosity = viscosity;
    water_compressibility = compressibility;
}

double NumericalAquifer::getTotalPoreVolume() const {
    double total = 0.0;
    for (const auto& cell : cells) {
        total += cell.volume * cell.porosity;
    }
    return total;
}

double NumericalAquifer::calculateInflux(double reservoir_pressure, double dt) {
    if (cells.empty()) return 0.0;
    
    // Simple implementation: use average aquifer pressure
    double avg_pressure = 0.0;
    double total_pv = 0.0;
    
    for (const auto& cell : cells) {
        double pv = cell.volume * cell.porosity;
        avg_pressure += cell.pressure * pv;
        total_pv += pv;
    }
    
    if (total_pv > 0.0) {
        avg_pressure /= total_pv;
    }
    
    // Calculate influx based on total transmissibility
    double total_trans = 0.0;
    for (const auto& conn : connections) {
        total_trans += conn.trans_mult * conn.area;  // Simplified
    }
    
    // Simple linear relationship
    double avg_perm = 0.0;
    for (const auto& cell : cells) {
        avg_perm += cell.permeability;
    }
    avg_perm /= cells.size();
    
    double T = total_trans * avg_perm / water_viscosity;
    double rate = T * (avg_pressure - reservoir_pressure);
    
    return std::max(0.0, rate);
}

void NumericalAquifer::updateState(double influx_rate, double dt) {
    cumulative_influx += influx_rate * dt;
    time += dt;
    
    // Distribute pressure change among cells
    double total_pv = getTotalPoreVolume();
    if (total_pv > 0.0) {
        double dp = -influx_rate * dt / (total_pv * water_compressibility);
        for (auto& cell : cells) {
            cell.pressure += dp;
        }
    }
    
    // Update average pressure
    double avg_pressure = 0.0;
    for (const auto& cell : cells) {
        avg_pressure += cell.pressure * cell.volume * cell.porosity;
    }
    current_pressure = avg_pressure / total_pv;
}

void NumericalAquifer::reset() {
    cumulative_influx = 0.0;
    time = 0.0;
    current_pressure = initial_pressure;
    
    for (auto& cell : cells) {
        cell.pressure = initial_pressure;
    }
}

void NumericalAquifer::configure(const std::map<std::string, std::string>& config) {
    if (config.count("water_density"))
        water_density = std::stod(config.at("water_density"));
    if (config.count("water_viscosity"))
        water_viscosity = std::stod(config.at("water_viscosity"));
    if (config.count("water_compressibility"))
        water_compressibility = std::stod(config.at("water_compressibility"));
    if (config.count("initial_pressure"))
        initial_pressure = std::stod(config.at("initial_pressure"));
}

double NumericalAquifer::calculateTransmissibility(const AquiferCell& c1, 
                                                   const AquiferCell& c2,
                                                   double dx, double dy, 
                                                   double dz) const {
    // Harmonic average of permeability
    double k_avg = 2.0 * c1.permeability * c2.permeability 
                 / (c1.permeability + c2.permeability);
    
    // Simplified transmissibility
    double area = dy * dz;
    double dist = dx;
    
    return k_avg * area / (water_viscosity * dist);
}

// ============================================================================
// ConstantPressureAquifer Implementation
// ============================================================================

ConstantPressureAquifer::ConstantPressureAquifer(int id)
    : AquiferModelBase(id, AquiferType::CONSTANT_PRESSURE),
      transmissibility(1e-6), constant_pressure(1e7) {
    initial_pressure = constant_pressure;
    current_pressure = constant_pressure;
}

double ConstantPressureAquifer::calculateInflux(double reservoir_pressure, double dt) {
    double rate = transmissibility * (constant_pressure - reservoir_pressure);
    return std::max(0.0, rate);
}

void ConstantPressureAquifer::updateState(double influx_rate, double dt) {
    cumulative_influx += influx_rate * dt;
    time += dt;
    // Pressure stays constant
}

void ConstantPressureAquifer::reset() {
    cumulative_influx = 0.0;
    time = 0.0;
}

void ConstantPressureAquifer::setTransmissibility(double T) {
    transmissibility = T;
}

void ConstantPressureAquifer::setConstantPressure(double P) {
    constant_pressure = P;
    initial_pressure = P;
    current_pressure = P;
}

void ConstantPressureAquifer::configure(const std::map<std::string, std::string>& config) {
    if (config.count("transmissibility"))
        transmissibility = std::stod(config.at("transmissibility"));
    if (config.count("constant_pressure"))
        setConstantPressure(std::stod(config.at("constant_pressure")));
}

// ============================================================================
// ConstantFluxAquifer Implementation
// ============================================================================

ConstantFluxAquifer::ConstantFluxAquifer(int id)
    : AquiferModelBase(id, AquiferType::CONSTANT_FLUX),
      flux_rate(0.001) {}

double ConstantFluxAquifer::calculateInflux(double reservoir_pressure, double dt) {
    return flux_rate;
}

void ConstantFluxAquifer::updateState(double influx_rate, double dt) {
    cumulative_influx += influx_rate * dt;
    time += dt;
}

void ConstantFluxAquifer::reset() {
    cumulative_influx = 0.0;
    time = 0.0;
}

void ConstantFluxAquifer::setFluxRate(double Q) {
    flux_rate = Q;
}

void ConstantFluxAquifer::configure(const std::map<std::string, std::string>& config) {
    if (config.count("flux_rate"))
        flux_rate = std::stod(config.at("flux_rate"));
}

// ============================================================================
// AquiferManager Implementation
// ============================================================================

AquiferManager::AquiferManager() {}

void AquiferManager::addCarterTracyAquifer(int id, 
                                           const std::map<std::string, std::string>& config) {
    auto aquifer = std::make_unique<CarterTracyAquifer>(id);
    aquifer->configure(config);
    aquifer_index[id] = aquifers.size();
    aquifers.push_back(std::move(aquifer));
}

void AquiferManager::addFetkovichAquifer(int id, 
                                         const std::map<std::string, std::string>& config) {
    auto aquifer = std::make_unique<FetkovichAquifer>(id);
    aquifer->configure(config);
    aquifer_index[id] = aquifers.size();
    aquifers.push_back(std::move(aquifer));
}

void AquiferManager::addNumericalAquifer(int id, 
                                         const std::map<std::string, std::string>& config) {
    auto aquifer = std::make_unique<NumericalAquifer>(id);
    aquifer->configure(config);
    aquifer_index[id] = aquifers.size();
    aquifers.push_back(std::move(aquifer));
}

void AquiferManager::addConstantPressureAquifer(int id, double pressure, 
                                                double transmissibility) {
    auto aquifer = std::make_unique<ConstantPressureAquifer>(id);
    aquifer->setConstantPressure(pressure);
    aquifer->setTransmissibility(transmissibility);
    aquifer_index[id] = aquifers.size();
    aquifers.push_back(std::move(aquifer));
}

void AquiferManager::addConstantFluxAquifer(int id, double flux_rate) {
    auto aquifer = std::make_unique<ConstantFluxAquifer>(id);
    aquifer->setFluxRate(flux_rate);
    aquifer_index[id] = aquifers.size();
    aquifers.push_back(std::move(aquifer));
}

void AquiferManager::addConnection(const AquiferConnection& conn) {
    auto it = aquifer_index.find(conn.aquifer_id);
    if (it != aquifer_index.end()) {
        aquifers[it->second]->addConnection(conn);
        
        // Add to cell map
        for (int i = conn.i1; i <= conn.i2; ++i) {
            for (int j = conn.j1; j <= conn.j2; ++j) {
                for (int k = conn.k1; k <= conn.k2; ++k) {
                    cell_aquifer_map[{i, j, k}].push_back(conn.aquifer_id);
                }
            }
        }
    }
}

double AquiferManager::calculateTotalInflux(double reservoir_pressure, double dt) {
    double total = 0.0;
    for (auto& aquifer : aquifers) {
        total += aquifer->calculateInflux(reservoir_pressure, dt);
    }
    return total;
}

double AquiferManager::calculateCellInflux(int i, int j, int k,
                                           double cell_pressure, double dt) {
    auto it = cell_aquifer_map.find({i, j, k});
    if (it == cell_aquifer_map.end()) {
        return 0.0;
    }
    
    double total = 0.0;
    for (int aq_id : it->second) {
        auto aq_it = aquifer_index.find(aq_id);
        if (aq_it != aquifer_index.end()) {
            total += aquifers[aq_it->second]->calculateInflux(cell_pressure, dt);
        }
    }
    return total;
}

void AquiferManager::updateAll(double dt) {
    for (auto& aquifer : aquifers) {
        // Note: In real implementation, would track actual influx
        aquifer->updateState(0.0, dt);
    }
}

void AquiferManager::resetAll() {
    for (auto& aquifer : aquifers) {
        aquifer->reset();
    }
}

AquiferModelBase* AquiferManager::getAquifer(int id) {
    auto it = aquifer_index.find(id);
    if (it != aquifer_index.end()) {
        return aquifers[it->second].get();
    }
    return nullptr;
}

double AquiferManager::getTotalCumulativeInflux() const {
    double total = 0.0;
    for (const auto& aquifer : aquifers) {
        total += aquifer->getCumulativeInflux();
    }
    return total;
}

void AquiferManager::printSummary() const {
    std::cout << "\n=== Aquifer Summary ===\n";
    for (const auto& aquifer : aquifers) {
        std::cout << "Aquifer " << aquifer->getId() << ": ";
        switch (aquifer->getType()) {
            case AquiferType::CARTER_TRACY: std::cout << "Carter-Tracy"; break;
            case AquiferType::FETKOVICH: std::cout << "Fetkovich"; break;
            case AquiferType::NUMERICAL: std::cout << "Numerical"; break;
            case AquiferType::CONSTANT_PRESSURE: std::cout << "Constant Pressure"; break;
            case AquiferType::CONSTANT_FLUX: std::cout << "Constant Flux"; break;
        }
        std::cout << "\n  Initial P: " << aquifer->getInitialPressure() / 1e6 << " MPa"
                  << "\n  Current P: " << aquifer->getCurrentPressure() / 1e6 << " MPa"
                  << "\n  Cumulative Influx: " << aquifer->getCumulativeInflux() << " m³\n";
    }
}

// ============================================================================
// Factory Functions
// ============================================================================

std::unique_ptr<AquiferModelBase> createAquifer(AquiferType type, int id,
                                                const std::map<std::string, std::string>& config) {
    std::unique_ptr<AquiferModelBase> aquifer;
    
    switch (type) {
        case AquiferType::CARTER_TRACY:
            aquifer = std::make_unique<CarterTracyAquifer>(id);
            break;
        case AquiferType::FETKOVICH:
            aquifer = std::make_unique<FetkovichAquifer>(id);
            break;
        case AquiferType::NUMERICAL:
            aquifer = std::make_unique<NumericalAquifer>(id);
            break;
        case AquiferType::CONSTANT_PRESSURE:
            aquifer = std::make_unique<ConstantPressureAquifer>(id);
            break;
        case AquiferType::CONSTANT_FLUX:
            aquifer = std::make_unique<ConstantFluxAquifer>(id);
            break;
        default:
            throw std::invalid_argument("Unknown aquifer type");
    }
    
    aquifer->configure(config);
    return aquifer;
}

AquiferType parseAquiferType(const std::string& type_str) {
    if (type_str == "CARTER_TRACY" || type_str == "AQUCT") {
        return AquiferType::CARTER_TRACY;
    } else if (type_str == "FETKOVICH" || type_str == "AQUFETP") {
        return AquiferType::FETKOVICH;
    } else if (type_str == "NUMERICAL" || type_str == "AQUNUM") {
        return AquiferType::NUMERICAL;
    } else if (type_str == "CONSTANT_PRESSURE") {
        return AquiferType::CONSTANT_PRESSURE;
    } else if (type_str == "CONSTANT_FLUX") {
        return AquiferType::CONSTANT_FLUX;
    }
    
    throw std::invalid_argument("Unknown aquifer type: " + type_str);
}

AquiferGeometry parseAquiferGeometry(const std::string& geom_str) {
    if (geom_str == "RADIAL" || geom_str == "radial") {
        return AquiferGeometry::RADIAL;
    } else if (geom_str == "LINEAR" || geom_str == "linear") {
        return AquiferGeometry::LINEAR;
    } else if (geom_str == "BOTTOM" || geom_str == "bottom") {
        return AquiferGeometry::BOTTOM;
    } else if (geom_str == "EDGE" || geom_str == "edge") {
        return AquiferGeometry::EDGE;
    }
    
    return AquiferGeometry::RADIAL;  // Default
}

AquiferFace parseAquiferFace(const std::string& face_str) {
    if (face_str == "I-" || face_str == "I_MINUS") return AquiferFace::I_MINUS;
    if (face_str == "I+" || face_str == "I_PLUS") return AquiferFace::I_PLUS;
    if (face_str == "J-" || face_str == "J_MINUS") return AquiferFace::J_MINUS;
    if (face_str == "J+" || face_str == "J_PLUS") return AquiferFace::J_PLUS;
    if (face_str == "K-" || face_str == "K_MINUS") return AquiferFace::K_MINUS;
    if (face_str == "K+" || face_str == "K_PLUS") return AquiferFace::K_PLUS;
    
    return AquiferFace::K_MINUS;  // Default (bottom water)
}

} // namespace FSRM
