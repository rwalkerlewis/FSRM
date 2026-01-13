#ifndef THERMAL_PRESSURIZATION_HPP
#define THERMAL_PRESSURIZATION_HPP

#include "FaultModel.hpp"
#include <vector>
#include <array>

namespace FSRM {

/**
 * @brief Thermal pressurization model for dynamic rupture
 * 
 * As frictional heating occurs during slip, the temperature of pore fluids
 * increases, causing thermal expansion. In low-permeability fault zones,
 * this leads to pore pressure increase, which reduces effective normal stress
 * and causes dramatic fault weakening.
 * 
 * This is a critical process in earthquake mechanics that can explain:
 * - Low apparent friction in large earthquakes
 * - Efficiency of seismic radiation
 * - Rupture propagation at low stress drop
 * 
 * References:
 * - Sibson (1973) - original concept
 * - Lachenbruch (1980) - analytical solutions
 * - Rice (2006) - full theory
 * - Noda & Lapusta (2010) - numerical implementation
 */
class ThermalPressurization {
public:
    /**
     * @brief Material properties for thermal pressurization
     */
    struct MaterialProperties {
        // Thermal properties
        double thermal_diffusivity;        // α_th [m²/s]
        double specific_heat;              // ρc [J/(m³·K)]
        
        // Hydraulic properties
        double hydraulic_diffusivity;      // α_hy [m²/s]
        double permeability;               // k [m²]
        double fluid_viscosity;            // μ [Pa·s]
        double porosity;                   // φ [-]
        
        // Compressibilities
        double fluid_compressibility;      // β_f [1/Pa]
        double pore_compressibility;       // β_p [1/Pa]
        
        // Thermal pressurization parameter
        double Lambda;                     // Λ [Pa/K] = ∂p/∂T at constant volume
        
        // Shear zone properties
        double half_width;                 // w [m] - half-width of shearing zone
        
        // Initial conditions
        double initial_temperature;        // T₀ [K]
        double initial_pore_pressure;      // p₀ [Pa]
        
        MaterialProperties() :
            thermal_diffusivity(1.0e-6),       // Typical for rocks
            specific_heat(2.7e6),              // J/(m³·K) for granite
            hydraulic_diffusivity(1.0e-4),     // m²/s
            permeability(1e-20),               // Very low permeability
            fluid_viscosity(1e-3),             // Water at room temp
            porosity(0.01),                    // 1% porosity
            fluid_compressibility(4.5e-10),    // Water
            pore_compressibility(1e-10),       // Typical
            Lambda(0.1e6),                     // 0.1 MPa/K typical
            half_width(0.01),                  // 1 cm shear zone
            initial_temperature(293.15),       // 20°C
            initial_pore_pressure(0.0) {}      // Hydrostatic pressure subtracted
    };
    
    /**
     * @brief State variables for TP
     */
    struct State {
        double temperature;         // T [K]
        double pore_pressure;       // p [Pa]
        double cumulative_slip;     // δ [m]
        double slip_rate;           // V [m/s]
        double shear_stress;        // τ [Pa]
        double normal_stress;       // σ_n [Pa]
        
        State() : temperature(293.15), pore_pressure(0.0),
                 cumulative_slip(0.0), slip_rate(0.0),
                 shear_stress(0.0), normal_stress(0.0) {}
    };
    
    ThermalPressurization();
    
    // Configuration
    void setMaterialProperties(const MaterialProperties& props);
    void enable(bool enable) { enabled = enable; }
    bool isEnabled() const { return enabled; }
    
    // Initialize state
    void initialize(State& state) const;
    
    // Time integration of TP equations
    void update(State& state, double dt);
    
    // Compute change in effective normal stress due to TP
    double getEffectiveNormalStressChange(const State& state) const;
    
    // Get current pore pressure and temperature
    double getPorePressure(const State& state) const { return state.pore_pressure; }
    double getTemperature(const State& state) const { return state.temperature; }
    
    // Compute frictional heating rate
    double getFrictionalHeatingRate(const State& state) const;
    
    // Analytical solutions (for verification)
    struct AnalyticalSolution {
        double temperature;
        double pore_pressure;
    };
    AnalyticalSolution getAnalyticalSolution(double slip, double slip_rate,
                                             double shear_stress, double time) const;
    
    // Get material properties
    const MaterialProperties& getMaterialProperties() const { return props; }
    
private:
    bool enabled;
    MaterialProperties props;
    
    // Diffusion solver for 1D problem across fault zone
    // Temperature and pressure diffuse perpendicular to fault
    struct DiffusionSolver {
        int n_points;                       // Grid points across fault
        double dx;                          // Grid spacing
        std::vector<double> T;              // Temperature profile
        std::vector<double> p;              // Pressure profile
        
        DiffusionSolver(int n = 51) : n_points(n) {
            T.resize(n, 293.15);
            p.resize(n, 0.0);
        }
        
        void solve_heat_equation(double dt, double alpha, double source_term);
        void solve_diffusion_equation(double dt, double alpha, double source_term);
    };
    
    mutable DiffusionSolver solver;
    
    // Compute effective diffusivities
    double effectiveThermalDiffusivity() const;
    double effectiveHydraulicDiffusivity() const;
    
    // Source term for heat equation (frictional heating)
    double heatSourceTerm(const State& state) const;
    
    // Coupling between temperature and pressure
    double pressureFromTemperature(double dT) const;
    double temperatureFromPressure(double dp) const;
};

/**
 * @brief TP-proxy slip-weakening friction law
 * 
 * Herrera et al. (2024) developed a computationally efficient proxy
 * for thermal pressurization that captures the essential physics without
 * solving the full diffusion equations.
 * 
 * Friction law:
 * τ = -σ'_n * [μ_d + (μ_s - μ_d) / (1 + δ/d_c)^α]
 * 
 * where α is the TP proxy exponent (typically 1/3).
 */
class TPProxyFriction : public FrictionModelBase {
public:
    struct Parameters {
        double mu_s;           // Static friction
        double mu_d;           // Dynamic friction (with full TP effect)
        double d_c;            // Critical slip distance
        double alpha;          // TP proxy exponent (default 1/3)
        double cohesion;       // Cohesion
        
        Parameters() : mu_s(0.6), mu_d(0.1), d_c(0.4), alpha(1.0/3.0), cohesion(0.0) {}
    };
    
    TPProxyFriction();
    
    double getFriction(double slip_rate, double state_var,
                      double sigma_n_eff) const override;
    
    double getStateEvolutionRate(double slip_rate, double state_var) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    void setParameters(const Parameters& p) { params = p; }
    const Parameters& getParameters() const { return params; }
    
private:
    Parameters params;
};

/**
 * @brief Rate-and-state friction with full thermal pressurization
 * 
 * Combines rate-and-state friction with TP diffusion equations.
 * This provides the most physically realistic model but is computationally
 * expensive due to solving diffusion equations at every fault point.
 */
class RateStateWithTP : public RateStateFriction {
public:
    RateStateWithTP(EvolutionLaw law = EvolutionLaw::AGING);
    
    // Override to include TP effects
    double getFriction(double slip_rate, double state_var,
                      double sigma_n_eff) const override;
    
    // Enable/disable TP
    void enableThermalPressurization(bool enable);
    void setTPProperties(const ThermalPressurization::MaterialProperties& props);
    
    // Update TP state (call before getFriction)
    void updateTPState(double dt, double slip_rate, double shear_stress, double normal_stress);
    
    // Get current TP state
    const ThermalPressurization::State& getTPState() const { return tp_state; }
    
private:
    std::unique_ptr<ThermalPressurization> thermal_press;
    ThermalPressurization::State tp_state;
};

/**
 * @brief Multiple nucleation episodes for dynamic rupture
 * 
 * SeisSol supports multiple artificial nucleation episodes to trigger
 * rupture at different times and locations. This is useful for:
 * - Simulating earthquake sequences
 * - Testing rupture interaction
 * - Modeling triggered seismicity
 */
struct NucleationEpisode {
    // Timing
    double start_time;              // t_0 [s]
    double characteristic_time;     // τ [s] - duration of nucleation
    
    // Spatial extent (can be point, line, or area)
    double x_center, y_center, z_center;  // Center location
    double radius;                         // Nucleation patch radius
    
    // Stress perturbation (added to ambient stress)
    double delta_shear_stress;      // Δτ [Pa]
    double delta_normal_stress;     // Δσ_n [Pa]
    
    // Alternative: specify as traction in fault coordinates
    double T_strike;                // Strike-direction traction
    double T_dip;                   // Dip-direction traction
    double T_normal;                // Normal traction
    
    // Temporal function type
    enum class TimeFunction {
        LINEAR,                     // Linear ramp
        SMOOTHSTEP,                 // Smooth (C^1 continuous)
        GAUSSIAN,                   // Gaussian pulse
        REGULARIZED                 // Regularized (exp-based)
    };
    TimeFunction time_function;
    
    // Methods
    NucleationEpisode() :
        start_time(0.0), characteristic_time(0.1),
        x_center(0.0), y_center(0.0), z_center(0.0), radius(100.0),
        delta_shear_stress(1e6), delta_normal_stress(0.0),
        T_strike(0.0), T_dip(0.0), T_normal(0.0),
        time_function(TimeFunction::SMOOTHSTEP) {}
    
    // Evaluate time function at time t
    double evaluateTimeFunction(double t) const;
    
    // Check if point is within nucleation patch
    bool isInside(double x, double y, double z) const;
    
    // Get stress perturbation at given time and location
    void getStressPerturbation(double t, double x, double y, double z,
                              double& delta_tau, double& delta_sigma_n) const;
};

/**
 * @brief Manager for multiple nucleation episodes
 */
class MultipleNucleation {
public:
    MultipleNucleation();
    
    // Add nucleation episode
    void addEpisode(const NucleationEpisode& episode);
    void clearEpisodes() { episodes.clear(); }
    
    // Get total stress perturbation from all active episodes
    void getTotalPerturbation(double t, double x, double y, double z,
                             double& delta_tau, double& delta_sigma_n) const;
    
    // Check which episodes are active
    std::vector<int> getActiveEpisodes(double t) const;
    
    // Get episodes
    const std::vector<NucleationEpisode>& getEpisodes() const { return episodes; }
    size_t getNumEpisodes() const { return episodes.size(); }
    
private:
    std::vector<NucleationEpisode> episodes;
};

/**
 * @brief Configuration helper for TP in config files
 */
struct ThermalPressurizationConfig {
    bool enabled = false;
    
    // Material properties
    double thermal_diffusivity = 1.0e-6;
    double specific_heat = 2.7e6;
    double hydraulic_diffusivity = 1.0e-4;
    double Lambda = 0.1e6;
    double half_width_shear_zone = 0.01;
    
    // Initial conditions
    double initial_temperature = 293.15;
    double initial_pore_pressure = 0.0;
    
    // Parse from config file section
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Convert to MaterialProperties
    ThermalPressurization::MaterialProperties toMaterialProperties() const;
};

} // namespace FSRM

#endif // THERMAL_PRESSURIZATION_HPP
