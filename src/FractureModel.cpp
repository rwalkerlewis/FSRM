#include "FractureModel.hpp"
#include "FractureNetwork.hpp"
#include <cmath>
#include <random>
#include <algorithm>

namespace FSRM {

// ============================================================================
// FractureModel Base Class
// ============================================================================

FractureModel::FractureModel(FractureType type)
    : frac_type(type), aperture(1e-4), permeability(1e-12), compressibility(1e-9) {}

void FractureModel::setAperture(double ap) {
    aperture = ap;
    // Cubic law: k_f = b^2/12
    permeability = aperture * aperture / 12.0;
}

void FractureModel::setPermeability(double perm) {
    permeability = perm;
}

void FractureModel::setCompressibility(double comp) {
    compressibility = comp;
}

// ============================================================================
// NaturalFractureNetwork
// ============================================================================

NaturalFractureNetwork::NaturalFractureNetwork()
    : FractureModel(FractureType::NATURAL),
      use_dual_porosity(false), shape_factor(0.0), matrix_block_size(1.0) {}

void NaturalFractureNetwork::setGeometry(const std::vector<double>& coords) {
    coordinates = coords;
}

void NaturalFractureNetwork::updateGeometry(double /* dt */) {
    // Natural fractures don't propagate in this model
}

void NaturalFractureNetwork::addFracture(const std::vector<double>& points, double ap) {
    Fracture frac;
    frac.points = points;
    frac.aperture = ap;
    frac.permeability = ap * ap / 12.0;  // Cubic law
    
    // Determine orientation
    if (points.size() >= 6) {
        double dx = points[3] - points[0];
        double dy = points[4] - points[1];
        double dz = points[5] - points[2];
        
        double abs_dx = std::abs(dx);
        double abs_dy = std::abs(dy);
        double abs_dz = std::abs(dz);
        
        if (abs_dx > abs_dy && abs_dx > abs_dz) frac.orientation = 0;
        else if (abs_dy > abs_dx && abs_dy > abs_dz) frac.orientation = 1;
        else frac.orientation = 2;
    }
    
    fractures.push_back(frac);
}

void NaturalFractureNetwork::generateStochasticNetwork(int num_fractures, 
                                                       double mean_length,
                                                       double std_length, 
                                                       int seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> length_dist(mean_length, std_length);
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> pos_dist(0.0, 1.0);
    
    for (int i = 0; i < num_fractures; ++i) {
        double length = std::max(0.1, length_dist(gen));
        double angle = angle_dist(gen);
        double x0 = pos_dist(gen) * 1000.0;  // Domain size
        double y0 = pos_dist(gen) * 1000.0;
        double z0 = pos_dist(gen) * 100.0;
        
        std::vector<double> points = {
            x0, y0, z0,
            x0 + length * std::cos(angle),
            y0 + length * std::sin(angle),
            z0
        };
        
        double ap = 1e-4 * (0.5 + 0.5 * pos_dist(gen));
        addFracture(points, ap);
    }
}

void NaturalFractureNetwork::importFromFractureNetwork(const FractureNetwork* network) {
    if (!network) return;
    
    fractures.clear();
    
    const auto& discrete_fractures = network->getFractures();
    
    for (const auto& dfrac : discrete_fractures) {
        Fracture frac;
        
        // Convert DiscreteFracture to internal representation
        // Center point and a point along the strike direction
        frac.points = {
            dfrac.center[0], dfrac.center[1], dfrac.center[2],
            dfrac.center[0] + dfrac.radius * dfrac.strike_dir[0],
            dfrac.center[1] + dfrac.radius * dfrac.strike_dir[1],
            dfrac.center[2] + dfrac.radius * dfrac.strike_dir[2]
        };
        
        frac.aperture = dfrac.aperture;
        frac.permeability = dfrac.permeability;
        
        // Determine dominant orientation from normal vector
        double abs_nx = std::abs(dfrac.normal[0]);
        double abs_ny = std::abs(dfrac.normal[1]);
        double abs_nz = std::abs(dfrac.normal[2]);
        
        if (abs_nx >= abs_ny && abs_nx >= abs_nz) {
            frac.orientation = 0;  // x-dominant (fracture in y-z plane)
        } else if (abs_ny >= abs_nx && abs_ny >= abs_nz) {
            frac.orientation = 1;  // y-dominant (fracture in x-z plane)
        } else {
            frac.orientation = 2;  // z-dominant (fracture in x-y plane)
        }
        
        fractures.push_back(frac);
    }
    
    // Update overall properties
    if (!fractures.empty()) {
        double total_aperture = 0.0;
        for (const auto& f : fractures) {
            total_aperture += f.aperture;
        }
        aperture = total_aperture / fractures.size();
        permeability = aperture * aperture / 12.0;  // Cubic law
    }
}

void NaturalFractureNetwork::enableDualPorosity(bool enable) {
    use_dual_porosity = enable;
}

void NaturalFractureNetwork::setShapeFactorModel(const std::string& model) {
    shape_factor = computeShapeFactor(model);
}

void NaturalFractureNetwork::setMatrixBlockSize(double size) {
    matrix_block_size = size;
}

double NaturalFractureNetwork::computeShapeFactor(const std::string& model) const {
    // Shape factor for fracture-matrix transfer
    if (model == "WARREN-ROOT") {
        // sigma = 4*n*(n+2) / L^2 where n is number of normal sets
        int n = 3;  // Three orthogonal sets
        return 4.0 * n * (n + 2) / (matrix_block_size * matrix_block_size);
    } else if (model == "KAZEMI") {
        // sigma = 4*(1/Lx^2 + 1/Ly^2 + 1/Lz^2)
        double L = matrix_block_size;
        return 12.0 / (L * L);
    }
    
    return 0.0;
}

void NaturalFractureNetwork::computeFluidExchange(Vec U, DM dm,
                                                  std::vector<double>& exchange_rates) const {
    (void)U; (void)dm;  // Suppress unused parameter warnings - stub implementation
    if (use_dual_porosity) {
        // Compute fracture-matrix transfer
        // Q = sigma * k_m/mu * (P_f - P_m)
        exchange_rates.clear();
        
        for (const auto& frac : fractures) {
            (void)frac;  // Suppress unused warning - placeholder
            // Simplified - would need to interpolate pressure from U
            double transfer_rate = shape_factor * 1e-15 / 0.001;  // k/mu
            exchange_rates.push_back(transfer_rate);
        }
    }
}

void NaturalFractureNetwork::contributeToResidual(Vec /* F */, Vec /* U */, DM /* dm */) const {
    // Add fracture flow contributions
}

void NaturalFractureNetwork::contributeToJacobian(Mat /* J */, Vec /* U */, DM /* dm */) const {
    // Add fracture Jacobian contributions
}

// ============================================================================
// HydraulicFractureModel
// ============================================================================

HydraulicFractureModel::HydraulicFractureModel()
    : FractureModel(FractureType::INDUCED_HYDRAULIC),
      fracture_model("P3D"), length(10.0), height(30.0), width(0.001),
      fracture_toughness(1e6), min_horizontal_stress(20e6), propagation_velocity(0.0),
      transport_proppant(false), proppant_diameter(0.0003), 
      proppant_density(2650.0), proppant_concentration(0.0),
      enable_leakoff_model(true), leakoff_coefficient(1e-7),
      youngs_modulus(10e9), poisson_ratio(0.25),
      fluid_density(1000.0), fluid_viscosity(0.001),
      proppant_pack_porosity(0.35) {}

void HydraulicFractureModel::setGeometry(const std::vector<double>& coords) {
    coordinates = coords;
    if (coords.size() >= 3) {
        fracture_center = {coords[0], coords[1], coords[2]};
    }
    if (coords.size() >= 6) {
        fracture_normal = {coords[3], coords[4], coords[5]};
    }
}

void HydraulicFractureModel::updateGeometry(double dt) {
    // Update fracture geometry based on propagation
    length += propagation_velocity * dt;
}

void HydraulicFractureModel::setFractureModel(const std::string& model) {
    fracture_model = model;
}

void HydraulicFractureModel::setPropagationCriteria(double Kc, double sigma_min) {
    fracture_toughness = Kc;
    min_horizontal_stress = sigma_min;
}

void HydraulicFractureModel::computePropagationDirection(const Vec /* stress_field */) {
    // Determine propagation based on stress field
    // Fracture propagates perpendicular to minimum principal stress
}

void HydraulicFractureModel::updateFractureGeometry(double pressure, double dt) {
    if (fracture_model == "PKN") {
        propagatePKN(pressure, dt);
    } else if (fracture_model == "KGD") {
        propagateKGD(pressure, dt);
    } else if (fracture_model == "P3D") {
        propagateP3D(pressure, dt);
    }
}

void HydraulicFractureModel::propagatePKN(double pressure, double dt) {
    // PKN (Perkins-Kern-Nordgren) model
    // Assumes height is constant, length >> height
    (void)dt;  // Used in full implementation for time-dependent propagation
    
    double net_pressure = pressure - min_horizontal_stress;
    if (net_pressure > 0) {
        // Simplified propagation criterion
        double K_I = net_pressure * std::sqrt(M_PI * length / 2.0);
        
        if (K_I > fracture_toughness) {
            propagation_velocity = 0.1 * (K_I - fracture_toughness) / fracture_toughness;
        } else {
            propagation_velocity = 0.0;
        }
        
        // Update width (PKN width profile)
        // Use user-configurable formation properties
        double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
        
        width = 4.0 * net_pressure * height / E_prime;
    }
}

void HydraulicFractureModel::propagateKGD(double pressure, double dt) {
    // KGD (Khristianovic-Geertsma-de Klerk) model
    // Assumes height >> length
    (void)dt;  // Used in full implementation for time-dependent propagation
    
    double net_pressure = pressure - min_horizontal_stress;
    if (net_pressure > 0) {
        // Use user-configurable formation properties
        double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
        
        width = 4.0 * net_pressure * length / E_prime;
        
        double K_I = net_pressure * std::sqrt(M_PI * length / 2.0);
        if (K_I > fracture_toughness) {
            propagation_velocity = 0.1 * (K_I - fracture_toughness) / fracture_toughness;
        }
    }
}

void HydraulicFractureModel::propagateP3D(double pressure, double dt) {
    // Pseudo-3D model
    // More sophisticated than PKN/KGD, allows variable height
    
    double net_pressure = pressure - min_horizontal_stress;
    if (net_pressure > 0) {
        // Simplified P3D
        double K_I = net_pressure * std::sqrt(M_PI * std::sqrt(length * height) / 2.0);
        
        if (K_I > fracture_toughness) {
            propagation_velocity = 0.05 * (K_I - fracture_toughness) / fracture_toughness;
            
            // Both length and height can grow
            length += propagation_velocity * dt;
            height += 0.5 * propagation_velocity * dt;
        }
    }
}

void HydraulicFractureModel::enableProppantTransport(bool enable) {
    transport_proppant = enable;
}

void HydraulicFractureModel::setProppantProperties(double diameter, double density, double concentration) {
    proppant_diameter = diameter;
    proppant_density = density;
    proppant_concentration = concentration;
}

void HydraulicFractureModel::computeProppantDistribution(const Vec velocity_field) {
    (void)velocity_field;  // Used in full implementation for advection
    
    if (!transport_proppant) return;
    
    // Solve transport equation for proppant
    // Simplified: uniform distribution for now
    proppant_distribution.resize(100, proppant_concentration);
    
    // Account for settling
    // Use user-configurable fluid properties
    double g = 9.81;
    double stokes_velocity = 2.0 * std::pow(proppant_diameter, 2.0) * 
                            (proppant_density - fluid_density) * g / (18.0 * fluid_viscosity);
    (void)stokes_velocity;  // Used in full implementation
    
    // Proppant accumulates at bottom
    for (size_t i = 0; i < proppant_distribution.size(); ++i) {
        double height_fraction = (double)i / proppant_distribution.size();
        proppant_distribution[i] *= (1.0 - 0.5 * height_fraction);
    }
}

void HydraulicFractureModel::enableLeakoff(bool enable) {
    enable_leakoff_model = enable;
}

void HydraulicFractureModel::setLeakoffCoefficient(double C_L) {
    leakoff_coefficient = C_L;
}

void HydraulicFractureModel::computeLeakoff(Vec /* U */, DM /* dm */, 
                                           std::vector<double>& leak_rates) const {
    if (!enable_leakoff_model) return;
    
    // Carter leak-off model: v_L = C_L / sqrt(t)
    // Total leak-off rate depends on fracture area and time
    
    double fracture_area = 2.0 * length * height;  // Both sides
    double avg_leak_rate = leakoff_coefficient * fracture_area;
    
    leak_rates.push_back(avg_leak_rate);
}

void HydraulicFractureModel::computeFluidExchange(Vec U, DM dm,
                                                  std::vector<double>& exchange_rates) const {
    // Combine leak-off and other exchange mechanisms
    std::vector<double> leak_rates;
    computeLeakoff(U, dm, leak_rates);
    
    exchange_rates = leak_rates;
}

void HydraulicFractureModel::contributeToResidual(Vec /* F */, Vec /* U */, DM /* dm */) const {
    // Add hydraulic fracture contributions to residual
}

void HydraulicFractureModel::contributeToJacobian(Mat /* J */, Vec /* U */, DM /* dm */) const {
    // Add hydraulic fracture Jacobian contributions
}

void HydraulicFractureModel::checkClosure() {
    // Check if fracture closes when pressure drops below closure stress
    // For now, assume fracture stays open due to proppant
}

void HydraulicFractureModel::setFormationProperties(double E, double nu) {
    youngs_modulus = E;
    poisson_ratio = nu;
}

void HydraulicFractureModel::setFluidProperties(double rho, double mu) {
    fluid_density = rho;
    fluid_viscosity = mu;
}

void HydraulicFractureModel::setProppantPackPorosity(double phi) {
    proppant_pack_porosity = phi;
}

void HydraulicFractureModel::computeEffectivePermeability(
    const std::vector<double>& proppant_conc, double& k_eff) const {
    
    (void)proppant_conc;  // Used in full implementation for concentration-dependent permeability
    
    // Carman-Kozeny equation for proppant pack
    // Use user-configurable proppant pack porosity
    double porosity_pack = proppant_pack_porosity;
    double particle_diameter = proppant_diameter;
    
    k_eff = porosity_pack * porosity_pack * porosity_pack * 
            particle_diameter * particle_diameter / 
            (180.0 * (1.0 - porosity_pack) * (1.0 - porosity_pack));
}

// ============================================================================
// FaultModel
// ============================================================================

FaultModel::FaultModel()
    : fault_length(1000.0), fault_width(500.0),
      static_friction(0.6), dynamic_friction(0.4),
      cohesion(1e6), dilation_angle(5.0 * M_PI / 180.0),
      use_rate_state(false), a_parameter(0.01), b_parameter(0.015),
      Dc_parameter(0.001), state_variable(1.0),
      permeability(1e-15), cumulative_slip(0.0),
      current_slip_mode(SlipMode::NO_SLIP) {}

void FaultModel::setFaultPlane(const std::vector<double>& strike,
                               const std::vector<double>& dip,
                               double length, double width) {
    fault_strike = strike;
    fault_dip = dip;
    fault_length = length;
    fault_width = width;
}

void FaultModel::setFrictionCoefficient(double static_mu, double dynamic_mu) {
    static_friction = static_mu;
    dynamic_friction = dynamic_mu;
}

void FaultModel::setCohesion(double c) {
    cohesion = c;
}

void FaultModel::setDilationAngle(double angle) {
    dilation_angle = angle;
}

void FaultModel::computeSlip(const Vec /* stress_field */, Vec /* displacement_field */) {
    // Compute fault slip based on stress state
    // This requires extracting stress at fault location from stress_field
}

void FaultModel::checkSlipCriteria(double shear_stress, double normal_stress,
                                   SlipMode& mode) const {
    // Mohr-Coulomb criterion: |tau| < c + mu * sigma_n
    
    double eff_normal_stress = normal_stress;  // Should include pore pressure effect
    
    double tau_max = cohesion + static_friction * eff_normal_stress;
    
    if (std::abs(shear_stress) < tau_max) {
        mode = SlipMode::NO_SLIP;
    } else {
        // Check if slip is stable or unstable
        if (use_rate_state) {
            // Rate-state friction can cause unstable slip
            if (b_parameter > a_parameter) {
                mode = SlipMode::SEISMIC;  // Velocity weakening
            } else {
                mode = SlipMode::STABLE_SLIP;  // Velocity strengthening
            }
        } else {
            mode = SlipMode::STICK_SLIP;
        }
    }
    
    // Don't assign to current_slip_mode in const function
    // current_slip_mode = mode;
}

void FaultModel::enableRateStateFriction(bool enable) {
    use_rate_state = enable;
}

void FaultModel::setRateStateParameters(double a, double b, double Dc) {
    a_parameter = a;
    b_parameter = b;
    Dc_parameter = Dc;
}

void FaultModel::updateStateVariable(double slip_rate, double dt) {
    if (!use_rate_state) return;
    
    // Aging law: dtheta/dt = 1 - V*theta/Dc
    double V0 = 1e-9;  // Reference velocity
    double dtheta_dt = 1.0 - slip_rate * state_variable / Dc_parameter;
    
    state_variable += dtheta_dt * dt;
    state_variable = std::max(state_variable, Dc_parameter / (100.0 * V0));
}

void FaultModel::computeEffectiveStress(double pore_pressure, double& eff_normal_stress) {
    // Effective stress: sigma' = sigma - alpha * P_p
    // Typically alpha = 1 for faults
    eff_normal_stress -= pore_pressure;
}

void FaultModel::computeMomentMagnitude(double slip_area, double slip_amount,
                                       double& magnitude) const {
    // Seismic moment: M0 = mu * A * D
    // where mu is shear modulus, A is rupture area, D is average slip
    
    double shear_modulus = 30e9;  // Pa
    double M0 = shear_modulus * slip_area * slip_amount;
    
    // Moment magnitude: Mw = (2/3) * log10(M0) - 6.07
    magnitude = (2.0 / 3.0) * std::log10(M0) - 6.07;
}

void FaultModel::updatePermeability(double slip_amount) {
    // Permeability can increase or decrease with slip
    // Typically increases initially due to dilatancy, then decreases due to gouge production
    
    cumulative_slip += slip_amount;
    
    // Simplified model
    double k0 = 1e-18;  // Initial fault permeability
    double k_max = 1e-14;  // Maximum due to dilatancy
    double slip_threshold = 0.01;  // meters
    
    if (cumulative_slip < slip_threshold) {
        // Dilatancy increases permeability
        permeability = k0 + (k_max - k0) * (cumulative_slip / slip_threshold);
    } else {
        // Gouge production reduces permeability
        permeability = k_max * std::exp(-(cumulative_slip - slip_threshold) / 0.1);
    }
}

} // namespace FSRM
