/*
 * Enhanced Geothermal System (EGS) Optimization with Stochastic Subsurface
 * 
 * This example demonstrates multi-objective optimization of an Enhanced Geothermal
 * System (EGS) considering stochastic subsurface representations. It implements:
 * 
 * 1. Stochastic subsurface modeling with geostatistical realizations
 *    - Spatially correlated permeability fields
 *    - Random fracture network generation
 *    - Thermal property uncertainty
 * 
 * 2. Multi-objective optimization
 *    - Maximize thermal energy extraction (MWth)
 *    - Maximize Net Present Value (NPV)
 *    - Minimize induced seismicity risk (Mmax)
 * 
 * 3. Monte Carlo uncertainty quantification
 *    - Propagate geological uncertainty through simulations
 *    - Compute confidence intervals on objectives
 *    - Risk-aware decision making
 * 
 * 4. Well placement and rate optimization
 *    - Injection/production well locations
 *    - Flow rates and pressures
 *    - Well spacing and configuration
 * 
 * 5. Pareto front analysis
 *    - Non-dominated solution identification
 *    - Trade-off visualization
 *    - Decision support metrics
 * 
 * Physical Model:
 * - Thermo-Hydro-Mechanical (THM) coupling
 * - Discrete fracture network (DFN) in stimulated reservoir
 * - Temperature-dependent fluid properties
 * - Induced seismicity on pre-existing faults
 * 
 * References:
 *   Tester et al. (2006) "The Future of Geothermal Energy" MIT Report
 *   McClure & Horne (2014) "Investigation of injection-induced seismicity" 
 *   Doe et al. (2014) "Discrete Fracture Network Simulations of EGS"
 * 
 * Usage:
 *   ./egs_optimization -c config/egs_optimization.config
 *   mpirun -np 4 ./egs_optimization -c config/egs_optimization.config
 * 
 *   # Quick mode with reduced Monte Carlo iterations
 *   ./egs_optimization -c config/egs_optimization.config -quick
 * 
 *   # Single realization evaluation
 *   ./egs_optimization -c config/egs_optimization.config -single
 * 
 *   # Specify number of Monte Carlo realizations
 *   ./egs_optimization -c config/egs_optimization.config -mc 50
 */

#include "Simulator.hpp"
#include "WellModel.hpp"
#include "FractureModel.hpp"
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
#include <numeric>
#include <functional>
#include <set>

//=============================================================================
// Local Structures for EGS Optimization
//=============================================================================
namespace EGSOpt {

// Forward declarations
struct StochasticRealization;
struct EGSScenario;

//-----------------------------------------------------------------------------
// 3D Point structure
//-----------------------------------------------------------------------------
struct Point3D {
    double x, y, z;
    
    Point3D(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
    
    double distance(const Point3D& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

//-----------------------------------------------------------------------------
// Discrete fracture for fracture network
//-----------------------------------------------------------------------------
struct DiscreteFracture {
    Point3D center;
    double strike;           // degrees from north
    double dip;              // degrees from horizontal
    double length;           // m - along-strike length
    double height;           // m - down-dip dimension
    double aperture;         // m - hydraulic aperture
    double transmissivity;   // m²/s
    
    // Computed properties
    double area() const { return length * height; }
    double permeability() const { return aperture * aperture / 12.0; }  // Cubic law
};

//-----------------------------------------------------------------------------
// Fault structure for seismicity modeling
//-----------------------------------------------------------------------------
struct Fault {
    Point3D center;
    double strike;
    double dip;
    double length;
    double width;
    double initial_friction;
    double a_minus_b;        // Rate-state friction parameter
    double critical_slip;    // Dc - characteristic slip distance
    double initial_shear;    // Initial shear stress
    double initial_normal;   // Initial normal stress
    
    // Seismicity state
    double cumulative_slip;
    double max_slip_rate;
    double seismic_moment;
    
    Fault() : cumulative_slip(0), max_slip_rate(0), seismic_moment(0) {}
    
    // Compute moment magnitude from slip
    double momentMagnitude(double slip, double shear_modulus = 30e9) const {
        double area = length * width;
        double M0 = shear_modulus * area * slip;
        return (std::log10(M0) - 9.1) / 1.5;
    }
};

//-----------------------------------------------------------------------------
// Stochastic subsurface realization
//-----------------------------------------------------------------------------
struct StochasticRealization {
    int realization_id;
    std::vector<double> permeability_field;   // Log-normal permeability
    std::vector<double> porosity_field;
    std::vector<double> thermal_conductivity;
    std::vector<DiscreteFracture> fractures;
    std::vector<Fault> faults;
    
    // Statistical parameters used to generate this realization
    double perm_mean;
    double perm_var;
    double perm_corr_length;
    double frac_density;
    double fault_density;
};

//-----------------------------------------------------------------------------
// Well configuration
//-----------------------------------------------------------------------------
struct WellConfig {
    Point3D location;
    double depth_top;        // Top of completion interval
    double depth_bottom;     // Bottom of completion interval
    bool is_injector;
    double flow_rate;        // m³/s (positive = injection, negative = production)
    double temperature;      // °C - injection temperature for injectors
    double bhp_limit;        // Pa - pressure limit
    
    std::string name;
};

//-----------------------------------------------------------------------------
// EGS development scenario
//-----------------------------------------------------------------------------
struct EGSScenario {
    std::string name;
    
    // Well configuration
    std::vector<WellConfig> wells;
    int num_injectors;
    int num_producers;
    double well_spacing;
    
    // Operating parameters
    double injection_rate_total;    // m³/s
    double injection_temperature;   // °C
    double production_bhp;          // Pa
    double operation_years;
    
    // Stimulation parameters
    double stimulated_volume;       // m³
    double enhanced_permeability;   // Average post-stimulation
    
    // Results (averaged over Monte Carlo realizations)
    double thermal_power_mean;      // MWth
    double thermal_power_std;
    double npv_mean;                // $M
    double npv_std;
    double max_magnitude_mean;      // Maximum induced Mw
    double max_magnitude_std;
    double thermal_breakthrough_years_mean;
    double lcoe_mean;               // $/MWh - Levelized Cost of Energy
    
    // Multi-objective results
    std::vector<double> pareto_thermal;
    std::vector<double> pareto_npv;
    std::vector<double> pareto_seismicity;
    
    EGSScenario() : 
        num_injectors(1), num_producers(1), well_spacing(500.0),
        injection_rate_total(0.05), injection_temperature(40.0),
        production_bhp(5.0e6), operation_years(30.0),
        stimulated_volume(0), enhanced_permeability(1e-14),
        thermal_power_mean(0), thermal_power_std(0),
        npv_mean(0), npv_std(0),
        max_magnitude_mean(0), max_magnitude_std(0),
        thermal_breakthrough_years_mean(0), lcoe_mean(0) {}
    
    void print(std::ostream& os) const {
        os << std::fixed << std::setprecision(2);
        os << "EGS Scenario: " << name << "\n";
        os << "  Configuration:\n";
        os << "    Injectors: " << num_injectors << ", Producers: " << num_producers << "\n";
        os << "    Well spacing: " << well_spacing << " m\n";
        os << "    Injection rate: " << injection_rate_total * 1000.0 << " L/s\n";
        os << "    Injection temp: " << injection_temperature << " °C\n";
        os << "    Operation period: " << operation_years << " years\n";
        os << "  Performance (mean ± std):\n";
        os << "    Thermal power: " << thermal_power_mean << " ± " << thermal_power_std << " MWth\n";
        os << "    NPV: $" << npv_mean << " ± " << npv_std << " M\n";
        os << "    Max induced Mw: " << max_magnitude_mean << " ± " << max_magnitude_std << "\n";
        os << "    Breakthrough: " << thermal_breakthrough_years_mean << " years\n";
        os << "    LCOE: $" << lcoe_mean << "/MWh\n";
    }
};

//-----------------------------------------------------------------------------
// Simulation result for a single realization
//-----------------------------------------------------------------------------
struct SimulationResult {
    int realization_id;
    
    // Time series data
    std::vector<double> times;           // years
    std::vector<double> thermal_power;   // MWth
    std::vector<double> production_temp; // °C
    std::vector<double> injection_pressure;  // MPa
    std::vector<double> cumulative_energy;   // MWh
    
    // Seismicity data
    std::vector<double> event_times;
    std::vector<double> event_magnitudes;
    double max_magnitude;
    int total_events;
    
    // Summary metrics
    double avg_thermal_power;
    double thermal_breakthrough_time;  // years to 10°C drop
    double cumulative_heat_extracted;  // MWth-years
    double npv;
    double irr;
    
    SimulationResult() : 
        max_magnitude(-5.0), total_events(0),
        avg_thermal_power(0), thermal_breakthrough_time(30.0),
        cumulative_heat_extracted(0), npv(0), irr(0) {}
};

} // namespace EGSOpt

static char help[] = 
    "Enhanced Geothermal System (EGS) Optimization with Stochastic Subsurface\n"
    "=========================================================================\n"
    "Multi-objective optimization of EGS considering geological uncertainty.\n\n"
    "Usage: mpirun -np N ./egs_optimization -c <config_file> [options]\n\n"
    "Options:\n"
    "  -c <file>           Configuration file (required)\n"
    "  -mc <N>             Number of Monte Carlo realizations (default: 20)\n"
    "  -quick              Quick mode (reduced MC iterations, coarser grid)\n"
    "  -single             Single realization evaluation\n"
    "  -seed <N>           Random seed for reproducibility\n"
    "  -output_dir <dir>   Output directory (default: output/egs_optimization)\n"
    "  -pareto             Generate Pareto front only\n"
    "  -sensitivity        Run sensitivity analysis\n"
    "\n";

namespace FSRM {

//=============================================================================
// Geostatistical Realization Generator
//=============================================================================
class GeostatisticalGenerator {
public:
    GeostatisticalGenerator(unsigned int seed = 42) : rng_(seed), seed_(seed) {}
    
    void setSeed(unsigned int seed) { 
        seed_ = seed;
        rng_.seed(seed); 
    }
    
    /**
     * @brief Generate spatially correlated log-normal permeability field
     * 
     * Uses Sequential Gaussian Simulation (simplified version)
     * 
     * @param nx, ny, nz Grid dimensions
     * @param dx, dy, dz Grid spacing (m)
     * @param mean_log_k Mean of log10(k)
     * @param var_log_k Variance of log10(k)
     * @param corr_length Correlation length (m)
     * @return Vector of permeability values
     */
    std::vector<double> generatePermeabilityField(
        int nx, int ny, int nz,
        double dx, double dy, double dz,
        double mean_log_k, double var_log_k, double corr_length) {
        
        int n_total = nx * ny * nz;
        std::vector<double> perm(n_total);
        std::vector<double> log_k(n_total);
        
        std::normal_distribution<double> normal(0.0, 1.0);
        
        // Generate correlated Gaussian field using FFT-based method (simplified)
        // For production code, use proper FFT-MA or SGS
        
        // First pass: uncorrelated
        for (int i = 0; i < n_total; ++i) {
            log_k[i] = normal(rng_);
        }
        
        // Apply spatial smoothing for correlation
        std::vector<double> smoothed(n_total, 0.0);
        double corr_cells = corr_length / std::min({dx, dy, dz});
        int filter_size = std::max(1, static_cast<int>(corr_cells));
        
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int idx = i + j * nx + k * nx * ny;
                    double sum = 0.0;
                    double weight_sum = 0.0;
                    
                    // Moving average filter
                    for (int dk = -filter_size; dk <= filter_size; ++dk) {
                        for (int dj = -filter_size; dj <= filter_size; ++dj) {
                            for (int di = -filter_size; di <= filter_size; ++di) {
                                int ii = std::clamp(i + di, 0, nx - 1);
                                int jj = std::clamp(j + dj, 0, ny - 1);
                                int kk = std::clamp(k + dk, 0, nz - 1);
                                int neighbor_idx = ii + jj * nx + kk * nx * ny;
                                
                                double dist = std::sqrt(di*di*dx*dx + dj*dj*dy*dy + dk*dk*dz*dz);
                                double weight = std::exp(-dist / corr_length);
                                
                                sum += log_k[neighbor_idx] * weight;
                                weight_sum += weight;
                            }
                        }
                    }
                    
                    smoothed[idx] = sum / weight_sum;
                }
            }
        }
        
        // Normalize and transform to target distribution
        double mean_s = 0.0, var_s = 0.0;
        for (double v : smoothed) mean_s += v;
        mean_s /= n_total;
        for (double v : smoothed) var_s += (v - mean_s) * (v - mean_s);
        var_s /= n_total;
        
        double std_s = std::sqrt(var_s);
        double std_target = std::sqrt(var_log_k);
        
        for (int i = 0; i < n_total; ++i) {
            double normalized = (smoothed[i] - mean_s) / std_s;
            double log_k_value = mean_log_k + normalized * std_target;
            perm[i] = std::pow(10.0, log_k_value);
        }
        
        return perm;
    }
    
    /**
     * @brief Generate discrete fracture network
     * 
     * @param domain_size Domain dimensions (m)
     * @param fracture_density P32 intensity (m²/m³)
     * @param mean_length Mean fracture length (m)
     * @param mean_height Mean fracture height (m)
     * @param mean_aperture Mean hydraulic aperture (m)
     * @param preferred_strike Preferred strike direction (degrees)
     * @param strike_std Standard deviation of strike (degrees)
     * @return Vector of fractures
     */
    std::vector<EGSOpt::DiscreteFracture> generateFractureNetwork(
        const EGSOpt::Point3D& domain_size,
        double fracture_density,
        double mean_length, double mean_height,
        double mean_aperture,
        double preferred_strike, double strike_std,
        double preferred_dip = 80.0, double dip_std = 10.0) {
        
        std::vector<EGSOpt::DiscreteFracture> fractures;
        
        double domain_volume = domain_size.x * domain_size.y * domain_size.z;
        double mean_area = mean_length * mean_height;
        int num_fractures = static_cast<int>(fracture_density * domain_volume / mean_area);
        
        std::uniform_real_distribution<double> uniform_x(0, domain_size.x);
        std::uniform_real_distribution<double> uniform_y(0, domain_size.y);
        std::uniform_real_distribution<double> uniform_z(-domain_size.z, 0);
        std::normal_distribution<double> length_dist(mean_length, mean_length * 0.3);
        std::normal_distribution<double> height_dist(mean_height, mean_height * 0.3);
        std::normal_distribution<double> aperture_dist(mean_aperture, mean_aperture * 0.5);
        std::normal_distribution<double> strike_dist(preferred_strike, strike_std);
        std::normal_distribution<double> dip_dist(preferred_dip, dip_std);
        
        for (int i = 0; i < num_fractures; ++i) {
            EGSOpt::DiscreteFracture frac;
            frac.center = EGSOpt::Point3D(uniform_x(rng_), uniform_y(rng_), uniform_z(rng_));
            frac.length = std::max(10.0, length_dist(rng_));
            frac.height = std::max(10.0, height_dist(rng_));
            frac.aperture = std::max(1e-5, aperture_dist(rng_));
            frac.strike = std::fmod(strike_dist(rng_) + 360.0, 360.0);
            frac.dip = std::clamp(dip_dist(rng_), 30.0, 90.0);
            frac.transmissivity = frac.aperture * frac.aperture * frac.aperture / 12.0;
            
            fractures.push_back(frac);
        }
        
        return fractures;
    }
    
    /**
     * @brief Generate fault network for seismicity modeling
     */
    std::vector<EGSOpt::Fault> generateFaultNetwork(
        const EGSOpt::Point3D& domain_size,
        double fault_density,
        double mean_length, double mean_width,
        double stress_gradient,  // Pa/m
        double pore_pressure_gradient) {
        
        std::vector<EGSOpt::Fault> faults;
        
        double domain_area = domain_size.x * domain_size.y;
        double mean_fault_area = mean_length * mean_width;
        int num_faults = std::max(1, static_cast<int>(fault_density * domain_area / mean_fault_area));
        
        std::uniform_real_distribution<double> uniform_x(0, domain_size.x);
        std::uniform_real_distribution<double> uniform_y(0, domain_size.y);
        std::uniform_real_distribution<double> uniform_z(-domain_size.z * 0.8, -domain_size.z * 0.2);
        std::normal_distribution<double> length_dist(mean_length, mean_length * 0.3);
        std::normal_distribution<double> width_dist(mean_width, mean_width * 0.3);
        std::uniform_real_distribution<double> strike_dist(0, 180);
        std::uniform_real_distribution<double> dip_dist(60, 90);
        std::uniform_real_distribution<double> friction_dist(0.5, 0.7);
        std::normal_distribution<double> a_b_dist(-0.004, 0.002);  // Velocity-weakening
        
        for (int i = 0; i < num_faults; ++i) {
            EGSOpt::Fault fault;
            fault.center = EGSOpt::Point3D(uniform_x(rng_), uniform_y(rng_), uniform_z(rng_));
            fault.length = std::max(100.0, length_dist(rng_));
            fault.width = std::max(50.0, width_dist(rng_));
            fault.strike = strike_dist(rng_);
            fault.dip = dip_dist(rng_);
            fault.initial_friction = friction_dist(rng_);
            fault.a_minus_b = a_b_dist(rng_);
            fault.critical_slip = 0.01;  // 1 cm
            
            // Stress state at fault depth
            double depth = -fault.center.z;
            fault.initial_normal = stress_gradient * depth;
            fault.initial_shear = 0.4 * fault.initial_normal;  // Critically stressed
            
            faults.push_back(fault);
        }
        
        return faults;
    }
    
    /**
     * @brief Generate complete stochastic realization
     */
    EGSOpt::StochasticRealization generateRealization(
        int realization_id,
        int nx, int ny, int nz,
        double dx, double dy, double dz,
        double mean_log_k, double var_log_k, double corr_length,
        double frac_density, double fault_density,
        const EGSOpt::Point3D& domain_size) {
        
        EGSOpt::StochasticRealization real;
        real.realization_id = realization_id;
        real.perm_mean = mean_log_k;
        real.perm_var = var_log_k;
        real.perm_corr_length = corr_length;
        real.frac_density = frac_density;
        real.fault_density = fault_density;
        
        // Generate permeability field
        real.permeability_field = generatePermeabilityField(
            nx, ny, nz, dx, dy, dz, mean_log_k, var_log_k, corr_length);
        
        // Generate porosity (correlated with permeability)
        real.porosity_field.resize(nx * ny * nz);
        for (size_t i = 0; i < real.permeability_field.size(); ++i) {
            // Kozeny-Carman-like relationship
            double log_k = std::log10(real.permeability_field[i]);
            double phi = 0.01 + 0.002 * (log_k + 15.0);  // Base porosity ~1%
            real.porosity_field[i] = std::clamp(phi, 0.005, 0.05);
        }
        
        // Generate thermal conductivity (slight variation)
        real.thermal_conductivity.resize(nx * ny * nz);
        std::normal_distribution<double> tc_dist(3.0, 0.3);  // W/(m·K) for granite
        for (size_t i = 0; i < real.thermal_conductivity.size(); ++i) {
            real.thermal_conductivity[i] = std::max(2.0, tc_dist(rng_));
        }
        
        // Generate fracture network
        real.fractures = generateFractureNetwork(
            domain_size, frac_density,
            100.0, 80.0,   // mean length, height
            5e-4,          // mean aperture
            45.0, 20.0);   // strike, strike_std
        
        // Generate fault network
        real.faults = generateFaultNetwork(
            domain_size, fault_density,
            500.0, 300.0,    // mean length, width
            26000.0,         // stress gradient (Pa/m)
            10000.0);        // pore pressure gradient
        
        return real;
    }
    
private:
    std::mt19937 rng_;
    unsigned int seed_;
};

//=============================================================================
// THM Simulator Wrapper
//=============================================================================
class EGSTHMSimulator {
public:
    EGSTHMSimulator() : 
        reservoir_temp_(200.0), injection_temp_(40.0),
        fluid_heat_capacity_(4180.0), fluid_density_(900.0) {}
    
    void setReservoirProperties(
        double reservoir_temp,      // °C
        double injection_temp,      // °C
        double rock_heat_capacity,  // J/(kg·K)
        double rock_density) {      // kg/m³
        
        reservoir_temp_ = reservoir_temp;
        injection_temp_ = injection_temp;
        rock_heat_capacity_ = rock_heat_capacity;
        rock_density_ = rock_density;
    }
    
    /**
     * @brief Run THM simulation for EGS scenario
     * 
     * Simplified thermal-hydraulic model with analytical/semi-analytical approach
     * For production, would call full numerical simulator
     */
    EGSOpt::SimulationResult runSimulation(
        const EGSOpt::EGSScenario& scenario,
        const EGSOpt::StochasticRealization& realization,
        double end_time_years,
        double dt_output_years = 0.5) {
        
        EGSOpt::SimulationResult result;
        result.realization_id = realization.realization_id;
        
        // Compute effective reservoir properties
        double avg_perm = computeAveragePermeability(realization);
        double avg_porosity = computeAverageProperty(realization.porosity_field);
        double avg_thermal_cond = computeAverageProperty(realization.thermal_conductivity);
        
        // Fracture-enhanced permeability
        double frac_perm = computeFractureEnhancement(realization, scenario);
        double effective_perm = avg_perm + frac_perm;
        
        // Compute thermal breakthrough time (simplified Gringarten model)
        double flow_rate = scenario.injection_rate_total;
        double well_spacing = scenario.well_spacing;
        
        // Thermal diffusivity
        double alpha = avg_thermal_cond / (rock_density_ * rock_heat_capacity_);
        
        // Characteristic time for thermal breakthrough
        double breakthrough_time = well_spacing * well_spacing / (4.0 * alpha);
        breakthrough_time /= (365.25 * 24.0 * 3600.0);  // Convert to years
        
        // Production temperature model (simplified)
        double t = 0.0;
        double initial_temp = reservoir_temp_;
        double min_temp = injection_temp_;
        
        while (t <= end_time_years) {
            result.times.push_back(t);
            
            // Temperature decline model
            double temp_ratio = 1.0 - std::pow(t / breakthrough_time, 0.5) * 
                               std::exp(-t / (2.0 * breakthrough_time));
            temp_ratio = std::max(0.1, temp_ratio);
            
            double prod_temp = min_temp + (initial_temp - min_temp) * temp_ratio;
            result.production_temp.push_back(prod_temp);
            
            // Thermal power
            double delta_T = prod_temp - scenario.injection_temperature;
            double thermal_power = fluid_density_ * fluid_heat_capacity_ * 
                                  flow_rate * delta_T / 1.0e6;  // MWth
            result.thermal_power.push_back(thermal_power);
            
            // Injection pressure (simplified)
            double viscosity = 0.0002;  // Pa·s at high temp
            double inj_pressure = computeInjectionPressure(
                flow_rate, effective_perm, viscosity, well_spacing);
            result.injection_pressure.push_back(inj_pressure / 1.0e6);  // MPa
            
            // Track thermal breakthrough
            if (prod_temp < (initial_temp - 10.0) && 
                result.thermal_breakthrough_time > t) {
                result.thermal_breakthrough_time = t;
            }
            
            t += dt_output_years;
        }
        
        // Compute seismic events
        computeSeismicity(scenario, realization, result);
        
        // Compute summary metrics
        result.avg_thermal_power = computeMean(result.thermal_power);
        result.cumulative_heat_extracted = result.avg_thermal_power * end_time_years;
        
        // Economic analysis
        result.npv = computeNPV(scenario, result);
        result.irr = computeIRR(scenario, result);
        
        return result;
    }
    
private:
    double reservoir_temp_;
    double injection_temp_;
    double fluid_heat_capacity_;
    double fluid_density_;
    double rock_heat_capacity_ = 790.0;
    double rock_density_ = 2650.0;
    
    double computeAveragePermeability(const EGSOpt::StochasticRealization& real) {
        // Geometric mean for permeability
        double log_sum = 0.0;
        for (double k : real.permeability_field) {
            log_sum += std::log(k);
        }
        return std::exp(log_sum / real.permeability_field.size());
    }
    
    double computeAverageProperty(const std::vector<double>& field) {
        return std::accumulate(field.begin(), field.end(), 0.0) / field.size();
    }
    
    double computeFractureEnhancement(
        const EGSOpt::StochasticRealization& real,
        const EGSOpt::EGSScenario& scenario) {
        
        // Sum fracture transmissivity contribution
        double total_trans = 0.0;
        for (const auto& frac : real.fractures) {
            // Weight by distance to wells
            double min_well_dist = 1e10;
            for (const auto& well : scenario.wells) {
                double dist = frac.center.distance(well.location);
                min_well_dist = std::min(min_well_dist, dist);
            }
            
            // Fractures near wells contribute more
            double weight = std::exp(-min_well_dist / 200.0);
            total_trans += frac.transmissivity * frac.area() * weight;
        }
        
        // Convert to equivalent permeability
        double reservoir_thickness = 500.0;  // m
        return total_trans / (scenario.stimulated_volume / reservoir_thickness);
    }
    
    double computeInjectionPressure(
        double flow_rate, double permeability, double viscosity, double spacing) {
        
        // Simplified radial flow
        double reservoir_thickness = 500.0;
        double wellbore_radius = 0.1;
        
        double dp = (viscosity * flow_rate / (2.0 * M_PI * permeability * reservoir_thickness)) *
                   std::log(spacing / 2.0 / wellbore_radius);
        
        return std::max(0.0, dp);
    }
    
    void computeSeismicity(
        const EGSOpt::EGSScenario& scenario,
        const EGSOpt::StochasticRealization& real,
        EGSOpt::SimulationResult& result) {
        
        result.max_magnitude = -5.0;
        result.total_events = 0;
        
        // Pressure perturbation from injection
        double flow_rate = scenario.injection_rate_total;
        double inj_pressure_increase = flow_rate * 1.0e9;  // Simplified
        
        for (const auto& fault : real.faults) {
            // Distance to nearest injector
            double min_dist = 1e10;
            for (const auto& well : scenario.wells) {
                if (well.is_injector) {
                    min_dist = std::min(min_dist, fault.center.distance(well.location));
                }
            }
            
            // Pressure perturbation at fault
            double dp = inj_pressure_increase * std::exp(-min_dist / 500.0);
            
            // Coulomb failure
            double effective_normal = fault.initial_normal - dp;
            double coulomb_stress = fault.initial_shear - 
                                   fault.initial_friction * effective_normal;
            
            if (coulomb_stress > 0) {
                // Fault slips
                double slip = coulomb_stress / 1.0e9;  // Simplified
                double Mw = fault.momentMagnitude(slip);
                
                if (Mw > 0) {
                    result.event_times.push_back(1.0);  // Simplified timing
                    result.event_magnitudes.push_back(Mw);
                    result.total_events++;
                    result.max_magnitude = std::max(result.max_magnitude, Mw);
                }
            }
        }
    }
    
    double computeMean(const std::vector<double>& data) {
        if (data.empty()) return 0.0;
        return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }
    
    double computeNPV(
        const EGSOpt::EGSScenario& scenario,
        const EGSOpt::SimulationResult& result) {
        
        double discount_rate = 0.08;
        double electricity_price = 80.0;  // $/MWh
        double capacity_factor = 0.90;
        double conversion_efficiency = 0.12;  // Thermal to electric
        
        // Capital costs
        double drilling_cost = 10.0e6 * (scenario.num_injectors + scenario.num_producers);
        double stimulation_cost = 5.0e6;
        double plant_cost = 2000.0 * result.avg_thermal_power * 1000.0 * conversion_efficiency;
        double capex = drilling_cost + stimulation_cost + plant_cost;
        
        // Annual revenue
        double annual_energy = result.avg_thermal_power * 8760.0 * 
                              capacity_factor * conversion_efficiency;  // MWh_e/year
        double annual_revenue = annual_energy * electricity_price;
        
        // Annual operating cost (5% of capex)
        double annual_opex = 0.05 * capex;
        
        // NPV calculation
        double npv = -capex;
        for (int year = 1; year <= static_cast<int>(scenario.operation_years); ++year) {
            double cash_flow = annual_revenue - annual_opex;
            // Account for temperature decline
            double temp_factor = result.production_temp[std::min(
                static_cast<size_t>(year * 2), result.production_temp.size() - 1)] / 
                result.production_temp[0];
            cash_flow *= temp_factor;
            
            npv += cash_flow / std::pow(1.0 + discount_rate, year);
        }
        
        return npv / 1.0e6;  // Convert to $M
    }
    
    double computeIRR(
        const EGSOpt::EGSScenario& /* scenario */,
        const EGSOpt::SimulationResult& result) {
        
        // Simplified IRR estimate
        if (result.npv > 0) {
            return 0.08 + 0.02 * result.npv / 10.0;
        }
        return 0.0;
    }
};

//=============================================================================
// Multi-Objective Optimizer
//=============================================================================
class MultiObjectiveOptimizer {
public:
    MultiObjectiveOptimizer() : num_mc_realizations_(20), rng_(42) {}
    
    void setMonteCarloIterations(int n) { num_mc_realizations_ = n; }
    void setSeed(unsigned int seed) { rng_.seed(seed); }
    
    /**
     * @brief Run optimization with Monte Carlo uncertainty quantification
     */
    std::pair<EGSOpt::EGSScenario, std::vector<EGSOpt::EGSScenario>>
    runOptimization(
        const ConfigReader& config,
        GeostatisticalGenerator& geo_gen,
        EGSTHMSimulator& simulator) {
        
        // Read optimization bounds from config
        double ws_min = config.getDouble("OPTIMIZATION", "well_spacing_min", 300.0);
        double ws_max = config.getDouble("OPTIMIZATION", "well_spacing_max", 800.0);
        double ws_step = config.getDouble("OPTIMIZATION", "well_spacing_step", 100.0);
        
        double rate_min = config.getDouble("OPTIMIZATION", "injection_rate_min", 0.02);
        double rate_max = config.getDouble("OPTIMIZATION", "injection_rate_max", 0.10);
        double rate_step = config.getDouble("OPTIMIZATION", "injection_rate_step", 0.02);
        
        // Grid dimensions for geostatistics
        int nx = config.getInt("GRID", "nx", 40);
        int ny = config.getInt("GRID", "ny", 40);
        int nz = config.getInt("GRID", "nz", 30);
        double Lx = config.getDouble("GRID", "Lx", 2000.0);
        double Ly = config.getDouble("GRID", "Ly", 2000.0);
        double Lz = config.getDouble("GRID", "Lz", 1500.0);
        double dx = Lx / nx, dy = Ly / ny, dz = Lz / nz;
        
        // Geostatistical parameters
        double mean_log_k = config.getDouble("STOCHASTIC", "mean_log_permeability", -15.0);
        double var_log_k = config.getDouble("STOCHASTIC", "variance_log_permeability", 1.0);
        double corr_length = config.getDouble("STOCHASTIC", "correlation_length", 100.0);
        double frac_density = config.getDouble("STOCHASTIC", "fracture_density_p32", 0.01);
        double fault_density = config.getDouble("STOCHASTIC", "fault_density", 0.001);
        
        EGSOpt::Point3D domain_size(Lx, Ly, Lz);
        
        // Reservoir properties
        double reservoir_temp = config.getDouble("THERMAL", "reservoir_temperature", 200.0);
        double injection_temp = config.getDouble("THERMAL", "injection_temperature", 40.0);
        simulator.setReservoirProperties(reservoir_temp, injection_temp, 790.0, 2650.0);
        
        std::vector<EGSOpt::EGSScenario> all_scenarios;
        EGSOpt::EGSScenario best_scenario;
        double best_utility = -1e30;
        
        std::cout << "Running multi-objective optimization with Monte Carlo UQ...\n";
        std::cout << "MC realizations: " << num_mc_realizations_ << "\n\n";
        
        // Grid search over design variables
        int scenario_count = 0;
        for (double ws = ws_min; ws <= ws_max; ws += ws_step) {
            for (double rate = rate_min; rate <= rate_max; rate += rate_step) {
                for (int n_inj = 1; n_inj <= 2; ++n_inj) {
                    for (int n_prod = 1; n_prod <= 2; ++n_prod) {
                        EGSOpt::EGSScenario scenario;
                        scenario.name = "ws" + std::to_string(int(ws)) + 
                                       "_rate" + std::to_string(int(rate*1000)) +
                                       "_nI" + std::to_string(n_inj) +
                                       "_nP" + std::to_string(n_prod);
                        scenario.well_spacing = ws;
                        scenario.injection_rate_total = rate;
                        scenario.num_injectors = n_inj;
                        scenario.num_producers = n_prod;
                        scenario.operation_years = config.getDouble("SIMULATION", "operation_years", 30.0);
                        scenario.injection_temperature = injection_temp;
                        
                        // Set up well locations
                        setupWellConfiguration(scenario, Lx, Ly, Lz);
                        
                        // Monte Carlo evaluation
                        std::vector<double> thermal_samples, npv_samples, mag_samples;
                        
                        for (int mc = 0; mc < num_mc_realizations_; ++mc) {
                            // Generate stochastic realization
                            auto realization = geo_gen.generateRealization(
                                mc, nx, ny, nz, dx, dy, dz,
                                mean_log_k, var_log_k, corr_length,
                                frac_density, fault_density, domain_size);
                            
                            // Run simulation
                            auto result = simulator.runSimulation(
                                scenario, realization, scenario.operation_years);
                            
                            thermal_samples.push_back(result.avg_thermal_power);
                            npv_samples.push_back(result.npv);
                            mag_samples.push_back(result.max_magnitude);
                        }
                        
                        // Compute statistics
                        scenario.thermal_power_mean = computeMean(thermal_samples);
                        scenario.thermal_power_std = computeStd(thermal_samples);
                        scenario.npv_mean = computeMean(npv_samples);
                        scenario.npv_std = computeStd(npv_samples);
                        scenario.max_magnitude_mean = computeMean(mag_samples);
                        scenario.max_magnitude_std = computeStd(mag_samples);
                        
                        // LCOE calculation
                        if (scenario.thermal_power_mean > 0) {
                            double annual_energy = scenario.thermal_power_mean * 8760 * 0.9 * 0.12;
                            double total_cost = 25.0e6 * (n_inj + n_prod) + 5.0e6;
                            scenario.lcoe_mean = total_cost / (annual_energy * scenario.operation_years);
                        }
                        
                        all_scenarios.push_back(scenario);
                        
                        // Multi-objective utility function
                        // Normalize objectives and combine
                        double utility = computeUtility(scenario);
                        if (utility > best_utility) {
                            best_utility = utility;
                            best_scenario = scenario;
                        }
                        
                        ++scenario_count;
                        if (scenario_count % 10 == 0) {
                            std::cout << "  Evaluated " << scenario_count << " scenarios...\n";
                        }
                    }
                }
            }
        }
        
        std::cout << "\nTotal scenarios evaluated: " << scenario_count << "\n";
        
        // Compute Pareto front
        computeParetoFront(all_scenarios);
        
        return {best_scenario, all_scenarios};
    }
    
    /**
     * @brief Identify Pareto-optimal solutions
     */
    void computeParetoFront(std::vector<EGSOpt::EGSScenario>& scenarios) {
        // For each scenario, check if it's dominated
        for (auto& s : scenarios) {
            bool dominated = false;
            for (const auto& other : scenarios) {
                if (&s == &other) continue;
                
                // Check if 'other' dominates 's'
                // Maximize thermal and NPV, minimize seismicity
                bool better_thermal = other.thermal_power_mean >= s.thermal_power_mean;
                bool better_npv = other.npv_mean >= s.npv_mean;
                bool better_seis = other.max_magnitude_mean <= s.max_magnitude_mean;
                
                bool strictly_better = 
                    (other.thermal_power_mean > s.thermal_power_mean) ||
                    (other.npv_mean > s.npv_mean) ||
                    (other.max_magnitude_mean < s.max_magnitude_mean);
                
                if (better_thermal && better_npv && better_seis && strictly_better) {
                    dominated = true;
                    break;
                }
            }
            
            if (!dominated) {
                s.pareto_thermal.push_back(s.thermal_power_mean);
                s.pareto_npv.push_back(s.npv_mean);
                s.pareto_seismicity.push_back(s.max_magnitude_mean);
            }
        }
    }
    
private:
    int num_mc_realizations_;
    std::mt19937 rng_;
    
    void setupWellConfiguration(EGSOpt::EGSScenario& scenario, 
                                double Lx, double Ly, double Lz) {
        scenario.wells.clear();
        
        double cx = Lx / 2.0;
        double cy = Ly / 2.0;
        double depth = Lz * 0.8;  // 80% of domain depth
        
        // Place injectors
        for (int i = 0; i < scenario.num_injectors; ++i) {
            EGSOpt::WellConfig well;
            well.is_injector = true;
            well.flow_rate = scenario.injection_rate_total / scenario.num_injectors;
            well.temperature = scenario.injection_temperature;
            well.depth_top = depth - 100;
            well.depth_bottom = depth + 100;
            
            // Position: offset from center
            double angle = 2.0 * M_PI * i / scenario.num_injectors;
            well.location = EGSOpt::Point3D(
                cx + scenario.well_spacing/2 * std::cos(angle),
                cy + scenario.well_spacing/2 * std::sin(angle),
                -depth);
            well.name = "INJ" + std::to_string(i + 1);
            
            scenario.wells.push_back(well);
        }
        
        // Place producers
        for (int i = 0; i < scenario.num_producers; ++i) {
            EGSOpt::WellConfig well;
            well.is_injector = false;
            well.flow_rate = -scenario.injection_rate_total / scenario.num_producers;
            well.depth_top = depth - 100;
            well.depth_bottom = depth + 100;
            well.bhp_limit = 5.0e6;  // 5 MPa
            
            // Position: opposite side from injectors
            double angle = M_PI + 2.0 * M_PI * i / scenario.num_producers;
            well.location = EGSOpt::Point3D(
                cx + scenario.well_spacing/2 * std::cos(angle),
                cy + scenario.well_spacing/2 * std::sin(angle),
                -depth);
            well.name = "PROD" + std::to_string(i + 1);
            
            scenario.wells.push_back(well);
        }
        
        // Stimulated volume estimate
        scenario.stimulated_volume = M_PI * std::pow(scenario.well_spacing, 2) * 200.0;
    }
    
    double computeMean(const std::vector<double>& data) {
        if (data.empty()) return 0.0;
        return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }
    
    double computeStd(const std::vector<double>& data) {
        if (data.size() < 2) return 0.0;
        double mean = computeMean(data);
        double sq_sum = 0.0;
        for (double v : data) {
            sq_sum += (v - mean) * (v - mean);
        }
        return std::sqrt(sq_sum / (data.size() - 1));
    }
    
    double computeUtility(const EGSOpt::EGSScenario& s) {
        // Multi-objective utility with risk adjustment
        // Weights: thermal=0.4, NPV=0.4, seismicity=0.2
        
        // Normalize to typical ranges
        double thermal_norm = s.thermal_power_mean / 10.0;  // 10 MWth target
        double npv_norm = (s.npv_mean + 50.0) / 100.0;      // -50 to +50 $M range
        double seis_norm = (3.0 - s.max_magnitude_mean) / 3.0;  // Mw 0-3 range
        
        // Risk penalty (penalize high std)
        double risk_penalty = 0.1 * (s.thermal_power_std / s.thermal_power_mean +
                                     s.npv_std / std::max(1.0, std::abs(s.npv_mean)));
        
        return 0.4 * thermal_norm + 0.4 * npv_norm + 0.2 * seis_norm - risk_penalty;
    }
};

//=============================================================================
// Output Writer
//=============================================================================
class EGSOutputWriter {
public:
    EGSOutputWriter(const std::string& output_dir) : output_dir_(output_dir) {}
    
    void writeResults(
        const EGSOpt::EGSScenario& best,
        const std::vector<EGSOpt::EGSScenario>& all_scenarios) {
        
        // Summary file
        std::ofstream summary(output_dir_ + "/egs_optimization_summary.txt");
        summary << "Enhanced Geothermal System (EGS) Optimization Results\n";
        summary << "======================================================\n\n";
        
        summary << "OPTIMAL SCENARIO:\n";
        summary << "-----------------\n";
        best.print(summary);
        
        summary << "\n\nALL SCENARIOS EVALUATED: " << all_scenarios.size() << "\n";
        summary << "========================\n\n";
        
        // Top scenarios by different objectives
        writeTopScenarios(summary, all_scenarios);
        
        summary.close();
        
        // CSV for analysis
        writeCSV(all_scenarios);
        
        // Pareto front data
        writeParetoFront(all_scenarios);
        
        // Gnuplot scripts
        writeGnuplotScripts();
    }
    
private:
    std::string output_dir_;
    
    void writeTopScenarios(std::ofstream& out, 
                          const std::vector<EGSOpt::EGSScenario>& scenarios) {
        auto sorted = scenarios;
        
        // Top by thermal power
        out << "TOP 5 BY THERMAL POWER:\n";
        std::sort(sorted.begin(), sorted.end(),
                 [](const auto& a, const auto& b) { 
                     return a.thermal_power_mean > b.thermal_power_mean; 
                 });
        writeScenarioTable(out, sorted, 5);
        
        // Top by NPV
        out << "\nTOP 5 BY NPV:\n";
        std::sort(sorted.begin(), sorted.end(),
                 [](const auto& a, const auto& b) { 
                     return a.npv_mean > b.npv_mean; 
                 });
        writeScenarioTable(out, sorted, 5);
        
        // Top by low seismicity
        out << "\nTOP 5 BY LOW SEISMICITY RISK:\n";
        std::sort(sorted.begin(), sorted.end(),
                 [](const auto& a, const auto& b) { 
                     return a.max_magnitude_mean < b.max_magnitude_mean; 
                 });
        writeScenarioTable(out, sorted, 5);
    }
    
    void writeScenarioTable(std::ofstream& out,
                           const std::vector<EGSOpt::EGSScenario>& scenarios,
                           int n) {
        out << std::setw(30) << "Scenario"
            << std::setw(12) << "WellSpc(m)"
            << std::setw(12) << "Rate(L/s)"
            << std::setw(8) << "nInj"
            << std::setw(8) << "nProd"
            << std::setw(12) << "Power(MW)"
            << std::setw(12) << "NPV($M)"
            << std::setw(10) << "Mmax"
            << "\n";
        
        for (int i = 0; i < std::min(n, (int)scenarios.size()); ++i) {
            const auto& s = scenarios[i];
            out << std::setw(30) << s.name
                << std::setw(12) << std::fixed << std::setprecision(0) << s.well_spacing
                << std::setw(12) << std::setprecision(1) << s.injection_rate_total * 1000.0
                << std::setw(8) << s.num_injectors
                << std::setw(8) << s.num_producers
                << std::setw(12) << std::setprecision(2) << s.thermal_power_mean
                << std::setw(12) << std::setprecision(1) << s.npv_mean
                << std::setw(10) << std::setprecision(2) << s.max_magnitude_mean
                << "\n";
        }
    }
    
    void writeCSV(const std::vector<EGSOpt::EGSScenario>& scenarios) {
        std::ofstream csv(output_dir_ + "/egs_optimization_results.csv");
        
        csv << "scenario,well_spacing_m,injection_rate_m3s,num_injectors,num_producers,"
            << "thermal_power_mean_MW,thermal_power_std_MW,"
            << "npv_mean_MM,npv_std_MM,"
            << "max_magnitude_mean,max_magnitude_std,"
            << "lcoe_per_MWh,is_pareto\n";
        
        for (const auto& s : scenarios) {
            bool is_pareto = !s.pareto_thermal.empty();
            csv << s.name << ","
                << s.well_spacing << ","
                << s.injection_rate_total << ","
                << s.num_injectors << ","
                << s.num_producers << ","
                << s.thermal_power_mean << ","
                << s.thermal_power_std << ","
                << s.npv_mean << ","
                << s.npv_std << ","
                << s.max_magnitude_mean << ","
                << s.max_magnitude_std << ","
                << s.lcoe_mean << ","
                << (is_pareto ? 1 : 0) << "\n";
        }
        csv.close();
    }
    
    void writeParetoFront(const std::vector<EGSOpt::EGSScenario>& scenarios) {
        std::ofstream pf(output_dir_ + "/pareto_front.dat");
        pf << "# Pareto front: thermal_power NPV max_magnitude\n";
        
        for (const auto& s : scenarios) {
            if (!s.pareto_thermal.empty()) {
                pf << s.thermal_power_mean << " " 
                   << s.npv_mean << " "
                   << s.max_magnitude_mean << " "
                   << "\"" << s.name << "\"\n";
            }
        }
        pf.close();
    }
    
    void writeGnuplotScripts() {
        // Pareto front visualization
        std::ofstream gp(output_dir_ + "/plot_pareto.gp");
        gp << R"(set terminal pngcairo enhanced size 1200,900 font 'Arial,12'
set output 'pareto_front.png'

set multiplot layout 2,2 title 'EGS Multi-Objective Optimization: Pareto Front' font ',16'

# Thermal vs NPV
set xlabel 'Thermal Power (MWth)'
set ylabel 'NPV ($M)'
set grid
plot 'pareto_front.dat' u 1:2 w p pt 7 ps 1.5 lc rgb 'red' title 'Pareto Optimal', \
     'egs_optimization_results.csv' u 6:8 skip 1 w p pt 6 ps 0.8 lc rgb 'gray' title 'All Scenarios'

# Thermal vs Seismicity
set xlabel 'Thermal Power (MWth)'
set ylabel 'Max Induced Magnitude'
plot 'pareto_front.dat' u 1:3 w p pt 7 ps 1.5 lc rgb 'red' title 'Pareto Optimal', \
     'egs_optimization_results.csv' u 6:10 skip 1 w p pt 6 ps 0.8 lc rgb 'gray' title 'All Scenarios'

# NPV vs Seismicity
set xlabel 'NPV ($M)'
set ylabel 'Max Induced Magnitude'
plot 'pareto_front.dat' u 2:3 w p pt 7 ps 1.5 lc rgb 'red' title 'Pareto Optimal', \
     'egs_optimization_results.csv' u 8:10 skip 1 w p pt 6 ps 0.8 lc rgb 'gray' title 'All Scenarios'

# 3D view
set xlabel 'Thermal Power (MWth)'
set ylabel 'NPV ($M)'
set zlabel 'Max Magnitude'
set view 60, 30
splot 'pareto_front.dat' u 1:2:3 w p pt 7 ps 1.5 lc rgb 'red' title 'Pareto Front'

unset multiplot
)";
        gp.close();
        
        // Sensitivity analysis plot
        std::ofstream gp2(output_dir_ + "/plot_sensitivity.gp");
        gp2 << R"(set terminal pngcairo enhanced size 1400,600 font 'Arial,12'
set output 'sensitivity_analysis.png'

set multiplot layout 1,3 title 'EGS Sensitivity Analysis' font ',16'

# Well spacing effect
set xlabel 'Well Spacing (m)'
set ylabel 'Thermal Power (MWth)'
set y2label 'NPV ($M)'
set y2tics
set grid
stats 'egs_optimization_results.csv' u 2:6 skip 1 nooutput
plot 'egs_optimization_results.csv' u 2:6 skip 1 w p pt 7 lc rgb 'blue' title 'Thermal Power', \
     'egs_optimization_results.csv' u 2:8 skip 1 w p pt 5 lc rgb 'green' axes x1y2 title 'NPV'

# Injection rate effect
set xlabel 'Injection Rate (m³/s)'
plot 'egs_optimization_results.csv' u 3:6 skip 1 w p pt 7 lc rgb 'blue' title 'Thermal Power', \
     'egs_optimization_results.csv' u 3:8 skip 1 w p pt 5 lc rgb 'green' axes x1y2 title 'NPV'

# Trade-off plot
unset y2label
unset y2tics
set xlabel 'Max Induced Magnitude'
set ylabel 'NPV ($M)'
set cblabel 'Thermal Power (MWth)'
set palette defined (0 'blue', 5 'yellow', 10 'red')
plot 'egs_optimization_results.csv' u 10:8:6 skip 1 w p pt 7 ps 1.2 palette title ''

unset multiplot
)";
        gp2.close();
    }
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
    char config_file[PETSC_MAX_PATH_LEN] = "config/egs_optimization.config";
    PetscOptionsGetString(nullptr, nullptr, "-c", config_file,
                         sizeof(config_file), nullptr);
    
    char output_dir[PETSC_MAX_PATH_LEN] = "output/egs_optimization";
    PetscOptionsGetString(nullptr, nullptr, "-output_dir", output_dir,
                         sizeof(output_dir), nullptr);
    
    PetscInt mc_iters = 20;
    PetscOptionsGetInt(nullptr, nullptr, "-mc", &mc_iters, nullptr);
    
    PetscBool quick_mode = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-quick", &quick_mode, nullptr);
    
    PetscBool single_mode = PETSC_FALSE;
    PetscOptionsGetBool(nullptr, nullptr, "-single", &single_mode, nullptr);
    
    PetscInt seed = 42;
    PetscOptionsGetInt(nullptr, nullptr, "-seed", &seed, nullptr);
    
    if (quick_mode) {
        mc_iters = 5;
    }
    if (single_mode) {
        mc_iters = 1;
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (rank == 0) {
        std::cout << "================================================================\n";
        std::cout << "  Enhanced Geothermal System (EGS) Optimization\n";
        std::cout << "  with Stochastic Subsurface and Multi-Objective Optimization\n";
        std::cout << "================================================================\n\n";
        std::cout << "Config file: " << config_file << "\n";
        std::cout << "Output directory: " << output_dir << "\n";
        std::cout << "Monte Carlo iterations: " << mc_iters << "\n";
        std::cout << "Random seed: " << seed << "\n\n";
    }
    
    // Read configuration
    FSRM::ConfigReader config;
    config.loadFile(config_file);
    
    // Initialize components
    FSRM::GeostatisticalGenerator geo_gen(seed);
    FSRM::EGSTHMSimulator simulator;
    FSRM::MultiObjectiveOptimizer optimizer;
    
    optimizer.setMonteCarloIterations(mc_iters);
    optimizer.setSeed(seed);
    
    // Run optimization
    auto [best_scenario, all_scenarios] = optimizer.runOptimization(config, geo_gen, simulator);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    // Output results
    if (rank == 0) {
        std::cout << "\n================================================================\n";
        std::cout << "  OPTIMIZATION RESULTS\n";
        std::cout << "================================================================\n\n";
        
        std::cout << "Scenarios evaluated: " << all_scenarios.size() << "\n";
        std::cout << "Computation time: " << duration.count() << " seconds\n\n";
        
        std::cout << "OPTIMAL EGS DEVELOPMENT:\n";
        std::cout << "------------------------\n";
        best_scenario.print(std::cout);
        
        // Count Pareto-optimal solutions
        int pareto_count = 0;
        for (const auto& s : all_scenarios) {
            if (!s.pareto_thermal.empty()) pareto_count++;
        }
        std::cout << "\nPareto-optimal solutions: " << pareto_count << "\n";
        
        // Create output directory and write results
        std::string mkdir_cmd = "mkdir -p " + std::string(output_dir);
        int ret = system(mkdir_cmd.c_str());
        (void)ret;
        
        FSRM::EGSOutputWriter writer(output_dir);
        writer.writeResults(best_scenario, all_scenarios);
        
        std::cout << "\n================================================================\n";
        std::cout << "  OUTPUT FILES\n";
        std::cout << "================================================================\n\n";
        std::cout << "Results written to: " << output_dir << "/\n";
        std::cout << "  - egs_optimization_summary.txt  (human-readable summary)\n";
        std::cout << "  - egs_optimization_results.csv  (all scenarios for analysis)\n";
        std::cout << "  - pareto_front.dat              (Pareto-optimal solutions)\n";
        std::cout << "  - plot_pareto.gp                (Gnuplot visualization)\n";
        std::cout << "  - plot_sensitivity.gp           (Sensitivity analysis)\n\n";
        
        std::cout << "KEY FINDINGS:\n";
        std::cout << "-------------\n";
        std::cout << "• Optimal well spacing: " << best_scenario.well_spacing << " m\n";
        std::cout << "• Optimal injection rate: " << best_scenario.injection_rate_total * 1000.0 << " L/s\n";
        std::cout << "• Well configuration: " << best_scenario.num_injectors << " injector(s), "
                 << best_scenario.num_producers << " producer(s)\n";
        std::cout << "• Expected thermal power: " << std::fixed << std::setprecision(1)
                 << best_scenario.thermal_power_mean << " ± " << best_scenario.thermal_power_std << " MWth\n";
        std::cout << "• Expected NPV: $" << best_scenario.npv_mean << " ± " << best_scenario.npv_std << " M\n";
        std::cout << "• Seismicity risk (Mmax): " << std::setprecision(2) 
                 << best_scenario.max_magnitude_mean << " ± " << best_scenario.max_magnitude_std << "\n";
        std::cout << "• LCOE: $" << std::setprecision(1) << best_scenario.lcoe_mean << "/MWh\n\n";
        
        // Recommendations
        std::cout << "RECOMMENDATIONS:\n";
        std::cout << "----------------\n";
        if (best_scenario.max_magnitude_mean > 2.0) {
            std::cout << "• Consider traffic light protocol for seismicity management\n";
            std::cout << "• Reduce injection rate if seismicity exceeds threshold\n";
        }
        if (best_scenario.well_spacing < 400.0) {
            std::cout << "• Tight well spacing - monitor for thermal short-circuiting\n";
        }
        if (best_scenario.npv_mean < 0) {
            std::cout << "• Negative NPV - consider longer operation or cost reduction\n";
        }
        if (best_scenario.thermal_power_mean > 5.0) {
            std::cout << "• Good thermal output - consider binary power plant\n";
        }
        
        std::cout << "\n================================================================\n";
        std::cout << "  EGS Optimization Complete\n";
        std::cout << "================================================================\n\n";
    }
    
    PetscFinalize();
    return 0;
}
