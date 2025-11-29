#ifndef TRACER_MODEL_HPP
#define TRACER_MODEL_HPP

/**
 * @file TracerModel.hpp
 * @brief ECLIPSE-compatible tracer modeling
 * 
 * Implements passive and reactive tracer simulation used in ECLIPSE:
 * - Passive water/oil/gas tracers
 * - Partitioning tracers (SWCT analysis)
 * - Reactive tracers with first-order decay
 * - Adsorbing tracers
 * - Inter-well tracer testing
 * 
 * ECLIPSE Keywords: TRACER, TRACERS, TRNUMF, TRNUMS
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <tuple>
#include <cmath>
#include <functional>

namespace FSRM {

/**
 * @brief Type of tracer (what phase it travels with)
 */
enum class TracerPhase {
    WATER,      ///< Water-soluble tracer
    OIL,        ///< Oil-soluble tracer
    GAS,        ///< Gas-soluble tracer
    FREE        ///< Free tracer (doesn't partition)
};

/**
 * @brief Tracer behavior type
 */
enum class TracerBehavior {
    PASSIVE,        ///< No interaction with rock (inert)
    ADSORBING,      ///< Adsorbs to rock surface (Langmuir/linear)
    DECAYING,       ///< First-order radioactive decay
    PARTITIONING,   ///< Partitions between phases
    REACTIVE        ///< First-order reaction/transformation
};

/**
 * @brief Adsorption model type
 */
enum class AdsorptionModel {
    NONE,           ///< No adsorption
    LINEAR,         ///< Linear isotherm: q = Kd * C
    LANGMUIR,       ///< Langmuir: q = qmax * K * C / (1 + K * C)
    FREUNDLICH      ///< Freundlich: q = Kf * C^n
};

/**
 * @brief Tracer injection specification
 */
struct TracerInjection {
    std::string well_name;          ///< Injection well name
    double start_time;              ///< Start of injection (s)
    double end_time;                ///< End of injection (s)
    double concentration;           ///< Injection concentration (kg/m³ or mol/m³)
    double mass_rate;               ///< Mass injection rate (kg/s) - alternative to concentration
    bool use_mass_rate;             ///< Use mass rate instead of concentration
    
    TracerInjection() : start_time(0.0), end_time(1e10), 
                        concentration(1.0), mass_rate(0.0), use_mass_rate(false) {}
};

/**
 * @brief Properties for a single tracer
 */
struct TracerProperties {
    std::string name;               ///< Tracer identifier
    TracerPhase phase;              ///< Carrier phase
    TracerBehavior behavior;        ///< Interaction behavior
    
    // Physical properties
    double molecular_weight;        ///< MW (g/mol)
    double diffusion_coefficient;   ///< Molecular diffusion (m²/s)
    double dispersivity_L;          ///< Longitudinal dispersivity (m)
    double dispersivity_T;          ///< Transverse dispersivity (m)
    
    // Decay properties
    double half_life;               ///< Half-life for decay (s)
    double decay_constant;          ///< λ = ln(2) / half_life
    
    // Adsorption properties
    AdsorptionModel adsorption_model;
    double Kd;                      ///< Linear partition coefficient (m³/kg)
    double Kf;                      ///< Freundlich coefficient
    double freundlich_n;            ///< Freundlich exponent
    double langmuir_qmax;           ///< Langmuir max capacity (kg/kg)
    double langmuir_K;              ///< Langmuir equilibrium constant
    
    // Partitioning properties
    double partition_coeff_ow;      ///< Oil-water partition coefficient
    double partition_coeff_gw;      ///< Gas-water partition coefficient
    
    // Reactive properties
    double reaction_rate;           ///< First-order reaction rate (1/s)
    std::string reaction_product;   ///< Product tracer name (if any)
    
    // Default constructor
    TracerProperties();
    
    // Calculate retardation factor
    double getRetardationFactor(double porosity, double rock_density) const;
    
    // Calculate effective dispersion
    double getDispersion(double velocity, double porosity) const;
};

/**
 * @brief Tracer concentration history at a point
 */
struct TracerBreakthroughCurve {
    std::string tracer_name;
    std::string location_name;      ///< Well or observation point
    std::vector<double> times;
    std::vector<double> concentrations;
    std::vector<double> cumulative_mass;
    
    // Analysis results
    double breakthrough_time;       ///< First arrival time
    double peak_time;               ///< Time of max concentration
    double peak_concentration;
    double mean_residence_time;
    double variance;
    double recovery_fraction;       ///< Mass recovered / mass injected
    
    // Compute statistics from data
    void computeStatistics(double total_injected_mass);
};

/**
 * @brief Base class for tracer transport solver
 */
class TracerSolver {
public:
    TracerSolver();
    virtual ~TracerSolver() = default;
    
    // Add tracers
    void addTracer(const TracerProperties& props);
    void addInjection(const std::string& tracer_name, const TracerInjection& inj);
    
    // Grid setup
    void setGrid(int nx, int ny, int nz, double dx, double dy, double dz);
    void setPorosity(const std::vector<double>& phi);
    void setPermeability(const std::vector<double>& k);
    void setRockDensity(double rho);
    
    // Flow field (from main simulator)
    void setVelocityField(const std::vector<double>& vx,
                          const std::vector<double>& vy,
                          const std::vector<double>& vz);
    void setSaturationField(const std::vector<double>& Sw,
                            const std::vector<double>& So,
                            const std::vector<double>& Sg);
    
    // Time stepping
    virtual void step(double dt) = 0;
    
    // Access concentrations
    const std::vector<double>& getConcentration(const std::string& tracer_name) const;
    double getConcentration(const std::string& tracer_name, int i, int j, int k) const;
    
    // Observation points
    void addObservationWell(const std::string& name, int i, int j, int k);
    TracerBreakthroughCurve getBreakthroughCurve(const std::string& tracer_name,
                                                  const std::string& well_name) const;
    
    // Well registration for injection/production
    void registerWell(const std::string& name, int i, int j, int k);
    
    // Output
    void writeVTK(const std::string& filename, int step) const;
    void writeSummary(const std::string& filename) const;
    
protected:
    // Grid
    int nx_, ny_, nz_;
    double dx_, dy_, dz_;
    int ncells_;
    
    // Rock properties
    std::vector<double> porosity_;
    std::vector<double> permeability_;
    double rock_density_;
    
    // Flow field
    std::vector<double> vx_, vy_, vz_;
    std::vector<double> Sw_, So_, Sg_;
    
    // Tracers
    std::map<std::string, TracerProperties> tracer_props_;
    std::map<std::string, std::vector<double>> concentrations_;
    std::map<std::string, std::vector<TracerInjection>> injections_;
    
    // Observation points
    struct ObsPoint {
        std::string name;
        int i, j, k;
        std::map<std::string, TracerBreakthroughCurve> curves;
    };
    std::vector<ObsPoint> obs_points_;
    
    // Time tracking
    double time_;
    
    // Well locations for injection lookup
    std::map<std::string, std::tuple<int, int, int>> well_locations_;
    
    // Helper functions
    int idx(int i, int j, int k) const { return i + j * nx_ + k * nx_ * ny_; }
    void applyInjections(double t, double dt);
    void recordObservations();
    int findWellIndex(const std::string& well_name) const;
};

/**
 * @brief Explicit finite difference tracer solver
 * 
 * Upwind scheme for advection, central difference for dispersion.
 * Simple but requires small time steps for stability.
 */
class ExplicitTracerSolver : public TracerSolver {
public:
    ExplicitTracerSolver();
    
    void step(double dt) override;
    
    // CFL stability check
    double getMaxTimeStep() const;
    void setMaxCourant(double cfl) { max_courant_ = cfl; }
    
private:
    double max_courant_;
    
    void advectionStep(const std::string& tracer_name, double dt);
    void dispersionStep(const std::string& tracer_name, double dt);
    void reactionStep(const std::string& tracer_name, double dt);
};

/**
 * @brief Implicit tracer solver using operator splitting
 * 
 * More stable for larger time steps.
 */
class ImplicitTracerSolver : public TracerSolver {
public:
    ImplicitTracerSolver();
    
    void step(double dt) override;
    
    // Solver settings
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    
private:
    int max_iterations_;
    double tolerance_;
    
    // Linear solver for implicit step
    void solveAdvectionDiffusion(const std::string& tracer_name, double dt);
};

/**
 * @brief Streamline-based tracer solver
 * 
 * Most efficient for dominated advection problems.
 */
class StreamlineTracerSolver : public TracerSolver {
public:
    StreamlineTracerSolver();
    
    void step(double dt) override;
    
    // Streamline options
    void setNumStreamlines(int n) { num_streamlines_ = n; }
    void traceStreamlines();
    
private:
    int num_streamlines_;
    
    struct Streamline {
        std::vector<double> x, y, z;    // Coordinates
        std::vector<double> tof;        // Time of flight
        std::vector<int> cells;         // Cell indices
    };
    std::vector<Streamline> streamlines_;
    
    void solve1D(Streamline& sl, const std::string& tracer_name, double dt);
};

/**
 * @brief Tracer analysis utilities
 */
class TracerAnalysis {
public:
    /**
     * @brief Single-well chemical tracer (SWCT) interpretation
     * 
     * Estimates residual oil saturation from partitioning tracer test.
     * 
     * @param reactive_curve Partitioning tracer breakthrough
     * @param passive_curve Non-partitioning tracer breakthrough
     * @param partition_coeff Oil-water partition coefficient
     * @return Estimated residual oil saturation
     */
    static double estimateSorFromSWCT(const TracerBreakthroughCurve& reactive_curve,
                                       const TracerBreakthroughCurve& passive_curve,
                                       double partition_coeff);
    
    /**
     * @brief Inter-well connectivity from tracer test
     * 
     * @param curves Breakthrough curves at producers
     * @return Map of producer name to connectivity fraction
     */
    static std::map<std::string, double> computeConnectivity(
        const std::vector<TracerBreakthroughCurve>& curves);
    
    /**
     * @brief Swept volume estimation
     * 
     * @param curve Breakthrough curve
     * @param injection_rate Injection rate during test (m³/s)
     * @return Swept pore volume (m³)
     */
    static double estimateSweptVolume(const TracerBreakthroughCurve& curve,
                                       double injection_rate);
    
    /**
     * @brief Heterogeneity characterization
     * 
     * @param curve Breakthrough curve
     * @return Lorenz coefficient (0 = homogeneous, 1 = very heterogeneous)
     */
    static double computeLorenzCoefficient(const TracerBreakthroughCurve& curve);
    
    /**
     * @brief Fit analytical solution to breakthrough curve
     * 
     * Fits 1D advection-dispersion solution to estimate Peclet number.
     * 
     * @param curve Measured breakthrough curve
     * @param distance Distance from injector to producer (m)
     * @return Estimated Peclet number
     */
    static double fitPecletNumber(const TracerBreakthroughCurve& curve,
                                   double distance);
};

/**
 * @brief Tracer manager for multi-tracer simulations
 */
class TracerManager {
public:
    TracerManager();
    
    // Configuration
    void configure(const std::map<std::string, std::string>& config);
    
    // Setup
    void setGrid(int nx, int ny, int nz, double dx, double dy, double dz);
    void setPorosity(const std::vector<double>& phi);
    void setPermeability(const std::vector<double>& k);
    void setRockDensity(double rho);
    
    // Add tracers from ECLIPSE-style input
    void parseTracerKeywords(const std::vector<std::string>& keywords);
    
    // Add individual tracers
    void addWaterTracer(const std::string& name, 
                        TracerBehavior behavior = TracerBehavior::PASSIVE);
    void addOilTracer(const std::string& name,
                      TracerBehavior behavior = TracerBehavior::PASSIVE);
    void addGasTracer(const std::string& name,
                      TracerBehavior behavior = TracerBehavior::PASSIVE);
    void addPartitioningTracer(const std::string& name, double Kow);
    void addRadioactiveTracer(const std::string& name, double half_life_days);
    
    // Injections
    void addInjection(const std::string& tracer_name, 
                      const std::string& well_name,
                      double concentration,
                      double start_time = 0.0, 
                      double end_time = 1e10);
    void addSlugInjection(const std::string& tracer_name,
                          const std::string& well_name,
                          double mass,
                          double duration);
    
    // Simulation
    void updateFlowField(const std::vector<double>& vx,
                         const std::vector<double>& vy,
                         const std::vector<double>& vz,
                         const std::vector<double>& Sw,
                         const std::vector<double>& So,
                         const std::vector<double>& Sg);
    void step(double dt);
    
    // Results
    std::vector<double> getConcentration(const std::string& tracer_name) const;
    TracerBreakthroughCurve getBreakthrough(const std::string& tracer_name,
                                            const std::string& location) const;
    
    // Output
    void addObservationWell(const std::string& name, int i, int j, int k);
    void writeOutput(const std::string& prefix, int step) const;
    void writeSummary(const std::string& filename) const;
    
    // Analysis
    double estimateResidualOil(const std::string& partitioning_tracer,
                                const std::string& passive_tracer,
                                const std::string& well) const;
    
private:
    std::unique_ptr<TracerSolver> solver_;
    std::vector<TracerProperties> tracers_;
    double time_;
    bool initialized_;
};

// ============================================================================
// TracerProperties Implementation
// ============================================================================

inline TracerProperties::TracerProperties()
    : name("TRACER1"),
      phase(TracerPhase::WATER),
      behavior(TracerBehavior::PASSIVE),
      molecular_weight(100.0),
      diffusion_coefficient(1e-9),
      dispersivity_L(1.0),
      dispersivity_T(0.1),
      half_life(1e20),
      decay_constant(0.0),
      adsorption_model(AdsorptionModel::NONE),
      Kd(0.0),
      Kf(0.0),
      freundlich_n(1.0),
      langmuir_qmax(0.0),
      langmuir_K(0.0),
      partition_coeff_ow(0.0),
      partition_coeff_gw(0.0),
      reaction_rate(0.0) {}

inline double TracerProperties::getRetardationFactor(double porosity, 
                                                      double rock_density) const {
    // R = 1 + (ρb/φ) * Kd for linear adsorption
    if (adsorption_model == AdsorptionModel::NONE || Kd <= 0.0) {
        return 1.0;
    }
    
    double bulk_density = rock_density * (1.0 - porosity);
    return 1.0 + (bulk_density / porosity) * Kd;
}

inline double TracerProperties::getDispersion(double velocity, double porosity) const {
    // D = Dm/τ + αL * v
    // where τ is tortuosity factor (approximated as φ^(-1/3))
    double tortuosity = std::pow(porosity, -1.0/3.0);
    double D_eff = diffusion_coefficient / tortuosity;
    double D_mech = dispersivity_L * std::abs(velocity);
    
    return D_eff + D_mech;
}

} // namespace FSRM

#endif // TRACER_MODEL_HPP
