#ifndef MULTI_STAGE_FRACTURING_HPP
#define MULTI_STAGE_FRACTURING_HPP

/**
 * @file MultiStageFracturing.hpp
 * @brief Comprehensive multi-stage hydraulic fracturing system
 * 
 * ResFrac-equivalent capabilities:
 * - Multi-cluster, multi-stage fracturing
 * - Limited-entry perforation design
 * - Cluster efficiency and uniformity
 * - Ball/plug drop operations
 * - Stage-by-stage pumping schedules
 * - Stress shadowing between clusters
 * - Perforation erosion modeling
 * - Real-time treatment analysis
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <array>
#include <cmath>

namespace FSRM {

/**
 * @brief Perforation cluster parameters
 */
struct PerforationCluster {
    int cluster_id;                     ///< Cluster identifier
    double measured_depth;              ///< MD along wellbore (m)
    double true_vertical_depth;         ///< TVD (m)
    int num_perforations;               ///< Number of perfs in cluster
    double perforation_diameter;        ///< Initial perf diameter (m)
    double perforation_length;          ///< Penetration depth (m)
    double phasing;                     ///< Perforation phasing (degrees)
    double shot_density;                ///< Shots per meter
    
    // Perforation state
    double current_diameter;            ///< Current diameter after erosion (m)
    double discharge_coefficient;       ///< Cd (typically 0.7-0.9)
    double friction_multiplier;         ///< Additional friction factor
    bool is_open;                       ///< Whether perfs are open
    double total_flow_area;             ///< Total open flow area (m²)
    
    // Stress state at cluster
    double min_horizontal_stress;       ///< Shmin at cluster depth (Pa)
    double max_horizontal_stress;       ///< SHmax at cluster depth (Pa)
    double vertical_stress;             ///< Sv at cluster depth (Pa)
    double pore_pressure;               ///< Initial pore pressure (Pa)
    double stress_shadow_delta;         ///< Stress increase from nearby fracs (Pa)
    
    // Flow allocation
    double allocated_rate;              ///< Current rate allocation (m³/s)
    double cumulative_volume;           ///< Total volume pumped (m³)
    double proppant_placed;             ///< Total proppant placed (kg)
    
    PerforationCluster();
    
    // Calculate perforation pressure drop
    double calculatePerfFriction(double rate, double fluid_density) const;
    
    // Update erosion state
    void updateErosion(double rate, double proppant_conc, double dt);
    
    // Calculate breakdown pressure
    double calculateBreakdownPressure() const;
};

/**
 * @brief Treatment stage parameters
 */
struct TreatmentStage {
    int stage_id;                       ///< Stage number
    std::string stage_name;             ///< Descriptive name
    double start_md;                    ///< Start measured depth (m)
    double end_md;                      ///< End measured depth (m)
    
    // Clusters in this stage
    std::vector<int> cluster_ids;       ///< IDs of clusters in stage
    
    // Isolation
    enum class IsolationMethod {
        PLUG,                           ///< Composite/cast iron plug
        BALL_DROP,                      ///< Ball-activated sliding sleeve
        PACKER,                         ///< Mechanical packer
        DIVERTER,                       ///< Chemical/degradable diverter
        OPENHOLE                        ///< No isolation (openhole)
    };
    
    IsolationMethod isolation_method;
    double isolation_depth;             ///< Depth of isolation device (m)
    bool isolation_confirmed;           ///< Isolation verified
    
    // Stage state
    enum class StageState {
        PENDING,                        ///< Not yet treated
        ISOLATING,                      ///< Setting isolation
        PERFORATING,                    ///< Creating perforations
        PUMPING,                        ///< Active treatment
        FLUSH,                          ///< Flushing stage
        COMPLETED,                      ///< Treatment complete
        FLOWBACK                        ///< In flowback
    };
    
    StageState state;
    double start_time;                  ///< Treatment start time (s)
    double end_time;                    ///< Treatment end time (s)
    
    // Treatment parameters
    double total_volume;                ///< Total fluid volume (m³)
    double total_proppant;              ///< Total proppant mass (kg)
    double max_rate;                    ///< Maximum pump rate (m³/s)
    double avg_pressure;                ///< Average treating pressure (Pa)
    double max_pressure;                ///< Maximum treating pressure (Pa)
    double isip;                        ///< Instantaneous shut-in pressure (Pa)
    double closure_pressure;            ///< Measured closure pressure (Pa)
    
    // Stage results
    double fracture_half_length;        ///< Created half-length (m)
    double fracture_height;             ///< Created height (m)
    double propped_half_length;         ///< Propped half-length (m)
    double propped_height;              ///< Propped height (m)
    double average_conductivity;        ///< Average propped conductivity (mD·ft)
    double cluster_efficiency;          ///< Fraction of clusters taking fluid
    double srw_volume;                  ///< Stimulated reservoir volume (m³)
    
    TreatmentStage();
};

/**
 * @brief Pumping schedule step
 */
struct PumpScheduleStep {
    std::string fluid_name;             ///< Fluid system name
    double rate;                        ///< Pump rate (m³/s)
    double volume;                      ///< Volume to pump (m³)
    double duration;                    ///< Duration (s) - calculated from rate/volume
    double proppant_concentration;      ///< Proppant concentration (kg/m³)
    std::string proppant_type;          ///< Proppant type name
    bool is_pad;                        ///< Is this a pad stage?
    bool is_flush;                      ///< Is this a flush stage?
    double fluid_viscosity;             ///< Fluid viscosity at downhole (Pa·s)
    double fluid_density;               ///< Fluid density (kg/m³)
    
    // Additives
    std::vector<std::pair<std::string, double>> additives;  ///< Name, concentration
    
    PumpScheduleStep();
};

/**
 * @brief Fracturing fluid system
 */
struct FracFluid {
    std::string name;                   ///< Fluid system name
    
    enum class FluidType {
        SLICKWATER,                     ///< Low-viscosity friction reducer
        LINEAR_GEL,                     ///< Linear gel (guar/HPG)
        CROSSLINKED_GEL,               ///< Crosslinked gel (borate/zirconate)
        VISCOELASTIC_SURFACTANT,       ///< VES fluid
        FOAM,                           ///< N2 or CO2 foam
        ACID,                           ///< Acid fracturing fluid
        HYBRID                          ///< Combination system
    };
    
    FluidType type;
    double base_viscosity;              ///< Base viscosity at surface (Pa·s)
    double density;                     ///< Density (kg/m³)
    double n_prime;                     ///< Power-law index n'
    double k_prime;                     ///< Consistency index K' (Pa·s^n)
    
    // Temperature dependence
    double viscosity_temp_coeff;        ///< dμ/dT coefficient
    double breakdown_temp;              ///< Temperature for gel breakdown (K)
    double breakdown_time;              ///< Time to break at breakdown_temp (s)
    
    // Leakoff properties
    double leakoff_coefficient;         ///< Carter leakoff C_L (m/s^0.5)
    double spurt_loss;                  ///< Spurt loss coefficient (m)
    double wall_building_coeff;         ///< Filter cake coefficient
    
    // Friction properties
    double friction_reduction;          ///< Friction reduction factor (0-1)
    
    FracFluid();
    
    // Get viscosity at conditions
    double getViscosity(double temperature, double shear_rate, double time) const;
    
    // Get effective leakoff
    double getLeakoff(double time, double pressure_diff) const;
};

/**
 * @brief Proppant properties
 */
struct Proppant {
    std::string name;                   ///< Proppant name
    
    enum class ProppantType {
        SAND,                           ///< Natural sand
        RESIN_COATED,                  ///< Resin-coated sand
        CERAMIC_LWC,                   ///< Lightweight ceramic
        CERAMIC_IDC,                   ///< Intermediate density ceramic
        CERAMIC_HDC,                   ///< High-density ceramic
        BAUXITE                        ///< Sintered bauxite
    };
    
    ProppantType type;
    double diameter;                    ///< Mean diameter (m)
    double specific_gravity;            ///< Specific gravity
    double bulk_density;                ///< Bulk density (kg/m³)
    double absolute_density;            ///< Particle density (kg/m³)
    double porosity;                    ///< Pack porosity
    double sphericity;                  ///< Krumbein sphericity
    double roundness;                   ///< Krumbein roundness
    
    // Strength properties
    double crush_strength;              ///< Crush resistance stress (Pa)
    double fines_generated;             ///< Fines at crush stress (fraction)
    
    // Permeability data (stress-dependent)
    std::vector<double> stress_points;  ///< Closure stress points (Pa)
    std::vector<double> permeability;   ///< Permeability at stress (m²)
    std::vector<double> beta_factor;    ///< Non-Darcy beta factor
    
    Proppant();
    
    // Get settling velocity in fluid
    double getSettlingVelocity(double fluid_viscosity, double fluid_density) const;
    
    // Get permeability at closure stress
    double getPermeability(double closure_stress) const;
    
    // Get conductivity (k*w)
    double getConductivity(double closure_stress, double width) const;
    
    // Get pack porosity based on sphericity
    double getPackPorosity() const;
};

/**
 * @brief Limited-entry design calculator
 */
class LimitedEntryDesign {
public:
    LimitedEntryDesign();
    
    /**
     * @brief Design limited-entry completion
     * 
     * Calculates number of perforations per cluster to achieve
     * uniform fluid distribution.
     * 
     * @param clusters Cluster locations and stresses
     * @param target_rate Target total pump rate (m³/s)
     * @param fluid_density Fluid density (kg/m³)
     * @param target_perf_friction Target perforation friction (Pa)
     * @return Optimized perforation count per cluster
     */
    std::vector<int> designPerfCount(
        const std::vector<PerforationCluster>& clusters,
        double target_rate,
        double fluid_density,
        double target_perf_friction);
    
    /**
     * @brief Calculate expected cluster efficiency
     * 
     * @param clusters Cluster configuration
     * @param rate Total pump rate (m³/s)
     * @param fluid_props Fluid properties
     * @return Expected efficiency (0-1) and rate per cluster
     */
    std::pair<double, std::vector<double>> predictClusterEfficiency(
        const std::vector<PerforationCluster>& clusters,
        double rate,
        const FracFluid& fluid_props);
    
    // Settings
    void setMinPerfFriction(double dp) { min_perf_friction_ = dp; }
    void setMaxPerfFriction(double dp) { max_perf_friction_ = dp; }
    void setDischargeCoefficient(double cd) { discharge_coeff_ = cd; }
    
private:
    double min_perf_friction_;          ///< Minimum target perf friction (Pa)
    double max_perf_friction_;          ///< Maximum acceptable (Pa)
    double discharge_coeff_;            ///< Perforation Cd
    
    // Calculate flow distribution
    std::vector<double> calculateFlowDistribution(
        const std::vector<PerforationCluster>& clusters,
        double total_rate,
        double fluid_density);
};

/**
 * @brief Multi-stage fracturing manager
 */
class MultiStageFracManager {
public:
    MultiStageFracManager();
    
    // Well setup
    void setWellGeometry(const std::vector<double>& md,
                         const std::vector<double>& tvd,
                         const std::vector<double>& inclination,
                         const std::vector<double>& azimuth);
    
    void setWellboreProperties(double casing_id, double tubing_id,
                               double roughness);
    
    // Stress profile
    void setStressProfile(const std::vector<double>& depths,
                          const std::vector<double>& shmin,
                          const std::vector<double>& shmax,
                          const std::vector<double>& sv,
                          const std::vector<double>& pp);
    
    // Rock properties profile
    void setRockProfile(const std::vector<double>& depths,
                        const std::vector<double>& youngs_modulus,
                        const std::vector<double>& poisson_ratio,
                        const std::vector<double>& toughness);
    
    // Add clusters and stages
    void addCluster(const PerforationCluster& cluster);
    void addStage(const TreatmentStage& stage);
    
    // Fluid and proppant systems
    void addFluidSystem(const FracFluid& fluid);
    void addProppant(const Proppant& proppant);
    
    // Pump schedule
    void setPumpSchedule(int stage_id, const std::vector<PumpScheduleStep>& schedule);
    
    /**
     * @brief Simulate a treatment stage
     * 
     * @param stage_id Stage to simulate
     * @param dt Time step (s)
     * @return Simulation results
     */
    struct StageSimResult {
        double final_half_length;
        double final_height;
        double propped_length;
        double propped_height;
        std::vector<double> cluster_volumes;
        std::vector<double> cluster_proppant;
        std::vector<double> cluster_efficiency;
        double total_leakoff;
        double isip;
        double closure_pressure;
        double net_pressure;
        std::vector<double> pressure_history;
        std::vector<double> rate_history;
        std::vector<double> time_history;
    };
    
    StageSimResult simulateStage(int stage_id, double dt = 1.0);
    
    /**
     * @brief Simulate full treatment (all stages)
     * 
     * @param dt Time step (s)
     * @return Results for all stages
     */
    std::vector<StageSimResult> simulateFullTreatment(double dt = 1.0);
    
    /**
     * @brief Update stress shadowing between clusters
     * 
     * Calculates stress perturbation from existing fractures
     * and updates cluster stress states.
     */
    void updateStressShadowing();
    
    // Optimization
    struct OptimizationObjective {
        double weight_length;           ///< Weight for fracture length
        double weight_height;           ///< Weight for height containment
        double weight_uniformity;       ///< Weight for cluster uniformity
        double weight_cost;             ///< Weight for treatment cost
    };
    
    /**
     * @brief Optimize cluster spacing
     * 
     * @param num_clusters Number of clusters per stage
     * @param stage_length Total stage length (m)
     * @param objective Optimization objective weights
     * @return Optimal cluster positions
     */
    std::vector<double> optimizeClusterSpacing(
        int num_clusters,
        double stage_length,
        const OptimizationObjective& objective);
    
    /**
     * @brief Optimize pump schedule
     * 
     * @param stage_id Stage to optimize
     * @param total_proppant Target proppant mass (kg)
     * @param max_rate Maximum pump rate (m³/s)
     * @param objective Optimization objective
     * @return Optimized pump schedule
     */
    std::vector<PumpScheduleStep> optimizePumpSchedule(
        int stage_id,
        double total_proppant,
        double max_rate,
        const OptimizationObjective& objective);
    
    // Access results
    const std::vector<PerforationCluster>& getClusters() const { return clusters_; }
    const std::vector<TreatmentStage>& getStages() const { return stages_; }
    TreatmentStage* getStage(int stage_id);
    PerforationCluster* getCluster(int cluster_id);
    
    // Output
    void writeStageSummary(std::ostream& os, int stage_id) const;
    void writeTreatmentReport(const std::string& filename) const;
    void writeVTKOutput(const std::string& filename, int stage_id) const;
    
private:
    // Well geometry
    std::vector<double> well_md_;
    std::vector<double> well_tvd_;
    std::vector<double> well_inc_;
    std::vector<double> well_azi_;
    double casing_id_;
    double tubing_id_;
    double roughness_;
    
    // Stress profile
    std::vector<double> stress_depths_;
    std::vector<double> shmin_;
    std::vector<double> shmax_;
    std::vector<double> sv_;
    std::vector<double> pp_;
    
    // Rock properties
    std::vector<double> rock_depths_;
    std::vector<double> youngs_modulus_;
    std::vector<double> poisson_ratio_;
    std::vector<double> toughness_;
    
    // Clusters and stages
    std::vector<PerforationCluster> clusters_;
    std::vector<TreatmentStage> stages_;
    
    // Fluids and proppants
    std::map<std::string, FracFluid> fluids_;
    std::map<std::string, Proppant> proppants_;
    
    // Pump schedules
    std::map<int, std::vector<PumpScheduleStep>> pump_schedules_;
    
    // Helper functions
    double interpolateStress(double depth, const std::vector<double>& depths,
                             const std::vector<double>& values) const;
    double interpolateRock(double depth, const std::vector<double>& depths,
                           const std::vector<double>& values) const;
    
    // Stress shadowing calculation
    double calculateStressShadow(double x, double y, double z,
                                 const std::vector<std::array<double, 6>>& existing_fracs) const;
    
    // Limited entry design helper
    LimitedEntryDesign limited_entry_;
};

/**
 * @brief Perforation erosion model
 */
class PerforationErosion {
public:
    PerforationErosion();
    
    /**
     * @brief Calculate erosion rate
     * 
     * Based on API RP 19B methodology
     * 
     * @param velocity Fluid velocity through perf (m/s)
     * @param proppant_conc Proppant concentration (kg/m³)
     * @param proppant_hardness Proppant hardness factor
     * @param time Exposure time (s)
     * @return Diameter increase rate (m/s)
     */
    double calculateErosionRate(double velocity, double proppant_conc,
                                double proppant_hardness, double time) const;
    
    /**
     * @brief Update perforation geometry
     * 
     * @param cluster Cluster to update
     * @param rate Flow rate (m³/s)
     * @param proppant_conc Proppant concentration (kg/m³)
     * @param dt Time step (s)
     */
    void updatePerforation(PerforationCluster& cluster, double rate,
                          double proppant_conc, double dt);
    
    // Model parameters
    void setErosionCoefficient(double c) { erosion_coeff_ = c; }
    void setHardnessExponent(double n) { hardness_exp_ = n; }
    void setVelocityExponent(double n) { velocity_exp_ = n; }
    
private:
    double erosion_coeff_;
    double hardness_exp_;
    double velocity_exp_;
};

/**
 * @brief Ball/plug drop simulator
 */
class IsolationSimulator {
public:
    IsolationSimulator();
    
    /**
     * @brief Simulate ball drop
     * 
     * @param ball_diameter Ball diameter (m)
     * @param seat_diameter Seat diameter (m)
     * @param flow_rate Current flow rate (m³/s)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @param well_md Measured depths along well
     * @param well_inc Inclination angles
     * @return Time to seat and seat depth
     */
    std::pair<double, double> simulateBallDrop(
        double ball_diameter,
        double seat_diameter,
        double flow_rate,
        double fluid_viscosity,
        const std::vector<double>& well_md,
        const std::vector<double>& well_inc);
    
    /**
     * @brief Simulate plug setting
     * 
     * @param plug_depth Target plug depth (m)
     * @param setting_pressure Required setting pressure (Pa)
     * @param wellbore_pressure Current wellbore pressure (Pa)
     * @return Success and actual set depth
     */
    std::pair<bool, double> simulatePlugSet(
        double plug_depth,
        double setting_pressure,
        double wellbore_pressure);
    
    // Ball properties
    void setBallDensity(double rho) { ball_density_ = rho; }
    void setBallDragCoefficient(double cd) { ball_cd_ = cd; }
    
private:
    double ball_density_;
    double ball_cd_;
    
    double calculateBallVelocity(double diameter, double fluid_density,
                                 double fluid_viscosity, double inclination) const;
};

/**
 * @brief Parse multi-stage configuration from file
 */
void parseMultiStageConfig(const std::string& filename,
                           MultiStageFracManager& manager);

} // namespace FSRM

#endif // MULTI_STAGE_FRACTURING_HPP
