#ifndef DIVERSION_MODELING_HPP
#define DIVERSION_MODELING_HPP

/**
 * @file DiversionModeling.hpp
 * @brief Diversion modeling for hydraulic fracturing treatments
 * 
 * ResFrac-equivalent capabilities:
 * - Degradable diverter particles (fiber, particulate)
 * - Ball sealers
 * - Chemical diversion
 * - Limited-entry mechanical diversion
 * - Diverter degradation kinetics
 * - Fracture temporary plugging
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <cmath>

namespace FSRM {

/**
 * @brief Diverter particle properties
 */
struct DiverterParticle {
    std::string name;                   ///< Particle type name
    
    enum class DiverterType {
        FIBER,                          ///< Degradable fiber
        PARTICULATE,                    ///< Particulate (sand/calcium carbonate)
        BALL_SEALER,                    ///< Ball sealers
        BENZOIC_ACID,                   ///< Benzoic acid flakes
        PLA,                            ///< Polylactic acid particles
        WAX,                            ///< Wax particles
        ROCK_SALT                       ///< Rock salt particles
    };
    
    DiverterType type;
    double diameter;                    ///< Mean diameter (m)
    double density;                     ///< Particle density (kg/m³)
    double sphericity;                  ///< Sphericity factor
    
    // Degradation properties
    double degradation_rate;            ///< First-order degradation rate (1/s)
    double activation_energy;           ///< Arrhenius activation energy (J/mol)
    double half_life_ref;               ///< Half-life at reference temp (s)
    double temperature_ref;             ///< Reference temperature (K)
    
    // Sealing properties
    double permeability_reduction;      ///< Permeability reduction factor
    double bridging_width_ratio;        ///< Width/diameter for bridging
    double max_pressure_hold;           ///< Maximum sealing pressure (Pa)
    
    // Size distribution
    std::vector<double> size_bins;      ///< Size bin edges (m)
    std::vector<double> mass_fractions; ///< Mass fraction in each bin
    
    DiverterParticle();
    
    // Calculate degradation at temperature
    double getDegradationRate(double temperature) const;
    
    // Calculate half-life at temperature
    double getHalfLife(double temperature) const;
    
    // Calculate settling velocity
    double getSettlingVelocity(double fluid_viscosity, 
                               double fluid_density) const;
};

/**
 * @brief Diverter stage in pumping schedule
 */
struct DiverterStage {
    std::string diverter_name;          ///< Diverter type name
    double start_time;                  ///< Start of diverter stage (s)
    double end_time;                    ///< End of diverter stage (s)
    double concentration;               ///< Concentration (kg/m³)
    double volume;                      ///< Total volume with diverter (m³)
    double mass;                        ///< Total diverter mass (kg)
    
    // Carrier fluid
    std::string carrier_fluid;          ///< Carrier fluid name
    double carrier_viscosity;           ///< Carrier viscosity (Pa·s)
    
    DiverterStage();
};

/**
 * @brief Temporary sealing plug state
 */
struct SealPlug {
    double location;                    ///< Location in fracture (m from wellbore)
    double width;                       ///< Plug width (m)
    double permeability;                ///< Current plug permeability (m²)
    double initial_mass;                ///< Initial diverter mass (kg)
    double current_mass;                ///< Current mass after degradation (kg)
    double pressure_drop;               ///< Pressure drop across plug (Pa)
    double formation_time;              ///< Time when plug formed (s)
    bool is_bridged;                    ///< True if fully bridged
    
    SealPlug();
};

/**
 * @brief Diverter transport and sealing model
 */
class DiverterTransport {
public:
    DiverterTransport();
    
    /**
     * @brief Initialize fracture grid for diverter tracking
     * 
     * @param length Fracture half-length (m)
     * @param height Fracture height (m)
     * @param nx Grid points along length
     * @param nz Grid points along height
     */
    void initializeGrid(double length, double height, int nx, int nz);
    
    /**
     * @brief Add diverter particle type
     */
    void addDiverterType(const DiverterParticle& diverter);
    
    /**
     * @brief Set fracture width field
     */
    void setWidthField(const std::vector<std::vector<double>>& width);
    
    /**
     * @brief Set velocity field
     */
    void setVelocityField(const std::vector<std::vector<double>>& vx,
                          const std::vector<std::vector<double>>& vz);
    
    /**
     * @brief Set temperature field
     */
    void setTemperatureField(const std::vector<std::vector<double>>& temperature);
    
    /**
     * @brief Inject diverter at wellbore
     * 
     * @param diverter_name Diverter type name
     * @param concentration Concentration (kg/m³)
     * @param rate Injection rate (m³/s)
     * @param dt Time step (s)
     */
    void injectDiverter(const std::string& diverter_name,
                        double concentration,
                        double rate,
                        double dt);
    
    /**
     * @brief Advance diverter transport one time step
     * 
     * @param dt Time step (s)
     * @return True if successful
     */
    bool step(double dt);
    
    /**
     * @brief Check for bridging and seal formation
     */
    void checkBridging();
    
    /**
     * @brief Update diverter degradation
     * 
     * @param dt Time step (s)
     */
    void updateDegradation(double dt);
    
    /**
     * @brief Get diverter concentration field
     */
    const std::vector<std::vector<double>>& getConcentration(
        const std::string& diverter_name) const;
    
    /**
     * @brief Get total diverter concentration
     */
    std::vector<std::vector<double>> getTotalConcentration() const;
    
    /**
     * @brief Get seal plugs
     */
    const std::vector<SealPlug>& getSealPlugs() const { return plugs_; }
    
    /**
     * @brief Get effective permeability multiplier
     * 
     * Returns 1.0 where no diverter, < 1.0 where diverter present.
     */
    std::vector<std::vector<double>> getPermeabilityMultiplier() const;
    
    /**
     * @brief Calculate diversion pressure
     * 
     * Total pressure drop due to diverter seals.
     * 
     * @param flow_rate Current flow rate (m³/s)
     * @return Additional pressure drop (Pa)
     */
    double calculateDiversionPressure(double flow_rate) const;
    
    /**
     * @brief Get diversion efficiency
     * 
     * @return Fraction of perforation blocked
     */
    double getDiversionEfficiency() const;
    
    // Output
    void writeVTK(const std::string& filename) const;
    
private:
    // Grid
    double length_, height_;
    int nx_, nz_;
    double dx_, dz_;
    
    // Fields
    std::vector<std::vector<double>> width_;
    std::vector<std::vector<double>> vx_, vz_;
    std::vector<std::vector<double>> temperature_;
    
    // Diverter types and concentrations
    std::map<std::string, DiverterParticle> diverter_types_;
    std::map<std::string, std::vector<std::vector<double>>> concentrations_;
    
    // Seal plugs
    std::vector<SealPlug> plugs_;
    
    // Time tracking
    double current_time_;
    
    // Helpers
    bool checkLocalBridging(int i, int j, double total_conc) const;
    void formSealPlug(int i, int j, double total_mass);
};

/**
 * @brief Ball sealer model
 */
class BallSealerModel {
public:
    BallSealerModel();
    
    /**
     * @brief Ball sealer properties
     */
    struct BallSealer {
        double diameter;                ///< Ball diameter (m)
        double density;                 ///< Ball density (kg/m³)
        double seat_diameter;           ///< Seat/perforation diameter (m)
        bool is_buoyant;                ///< Buoyant or sinking ball
        double drag_coefficient;        ///< Cd for flow calculation
        
        BallSealer();
    };
    
    /**
     * @brief Add ball sealers
     * 
     * @param num_balls Number of balls to drop
     * @param ball Ball properties
     */
    void addBallSealers(int num_balls, const BallSealer& ball);
    
    /**
     * @brief Set perforation locations and flow rates
     * 
     * @param perf_depths Perforation depths (m TVD)
     * @param perf_rates Flow rates into each perf (m³/s)
     */
    void setPerforations(const std::vector<double>& perf_depths,
                         const std::vector<double>& perf_rates);
    
    /**
     * @brief Simulate ball drop
     * 
     * @param fluid_velocity Fluid velocity in wellbore (m/s)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @param fluid_density Fluid density (kg/m³)
     * @param dt Time step (s)
     * @return Indices of sealed perforations
     */
    std::vector<int> simulateDrop(double fluid_velocity,
                                   double fluid_viscosity,
                                   double fluid_density,
                                   double dt);
    
    /**
     * @brief Check which perforations are sealed
     */
    std::vector<bool> getSealedPerforations() const;
    
    /**
     * @brief Get ball positions
     */
    std::vector<double> getBallPositions() const { return ball_positions_; }
    
    /**
     * @brief Calculate diversion pressure for sealed perfs
     */
    double calculateSealPressure(int perf_index, double flow_rate) const;
    
private:
    std::vector<BallSealer> balls_;
    std::vector<double> ball_positions_;
    std::vector<int> ball_status_;      ///< -1=falling, perf_idx=seated
    
    std::vector<double> perf_depths_;
    std::vector<double> perf_rates_;
    std::vector<bool> perf_sealed_;
    
    double calculateBallVelocity(const BallSealer& ball,
                                  double fluid_velocity,
                                  double fluid_viscosity,
                                  double fluid_density) const;
    
    bool checkSeat(const BallSealer& ball, double position) const;
};

/**
 * @brief Chemical diversion model
 */
class ChemicalDiversion {
public:
    ChemicalDiversion();
    
    /**
     * @brief Chemical diverter properties
     */
    struct ChemicalDiverter {
        std::string name;
        
        enum class ChemicalType {
            VES,                        ///< Viscoelastic surfactant
            CROSSLINKED_GEL,           ///< Self-diverting gel
            IN_SITU_CROSSLINK,         ///< In-situ crosslinking
            FOAM                        ///< Foamed diverter
        };
        
        ChemicalType type;
        
        // Viscosity behavior
        double initial_viscosity;       ///< Initial viscosity (Pa·s)
        double peak_viscosity;          ///< Peak viscosity (Pa·s)
        double time_to_peak;            ///< Time to reach peak (s)
        double break_time;              ///< Time to break (s)
        
        // Temperature sensitivity
        double activation_temperature;  ///< Temperature for activation (K)
        double break_temperature;       ///< Temperature for breaking (K)
        
        // Shear sensitivity
        double shear_thinning_index;    ///< Power-law n
        
        ChemicalDiverter();
    };
    
    /**
     * @brief Add chemical diverter system
     */
    void addDiverterSystem(const ChemicalDiverter& diverter);
    
    /**
     * @brief Calculate effective viscosity
     * 
     * @param time Elapsed time since injection (s)
     * @param temperature Local temperature (K)
     * @param shear_rate Local shear rate (1/s)
     * @return Effective viscosity (Pa·s)
     */
    double calculateViscosity(double time, double temperature,
                               double shear_rate) const;
    
    /**
     * @brief Calculate resistance factor
     * 
     * RF = k_water / k_diverter = mobility ratio
     * 
     * @param time Elapsed time (s)
     * @param temperature Temperature (K)
     * @return Resistance factor
     */
    double getResistanceFactor(double time, double temperature) const;
    
    /**
     * @brief Check if diverter is broken
     * 
     * @param time Elapsed time (s)
     * @param temperature Temperature (K)
     * @return True if broken
     */
    bool isBroken(double time, double temperature) const;
    
private:
    std::vector<ChemicalDiverter> diverters_;
};

/**
 * @brief Integrated diversion manager
 */
class DiversionManager {
public:
    DiversionManager();
    
    // Configuration
    void enableParticulate(bool enable) { use_particulate_ = enable; }
    void enableBallSealers(bool enable) { use_ball_sealers_ = enable; }
    void enableChemical(bool enable) { use_chemical_ = enable; }
    
    // Setup
    void setFractureGeometry(double length, double height, int nx, int nz);
    void setWidthField(const std::vector<std::vector<double>>& width);
    void setVelocityField(const std::vector<std::vector<double>>& vx,
                          const std::vector<std::vector<double>>& vz);
    void setTemperatureField(const std::vector<std::vector<double>>& temp);
    
    // Add diverter systems
    void addParticulateDiverter(const DiverterParticle& diverter);
    void addBallSealers(int count, const BallSealerModel::BallSealer& ball);
    void addChemicalDiverter(const ChemicalDiversion::ChemicalDiverter& chemical);
    
    // Set perforation data
    void setPerforations(const std::vector<double>& depths,
                         const std::vector<double>& rates);
    
    /**
     * @brief Pump diverter stage
     * 
     * @param stage Diverter stage definition
     */
    void pumpDiverterStage(const DiverterStage& stage);
    
    /**
     * @brief Advance simulation
     * 
     * @param dt Time step (s)
     * @return Total diversion pressure (Pa)
     */
    double step(double dt);
    
    /**
     * @brief Get effective flow distribution multiplier
     * 
     * For each cluster, returns multiplier on flow rate due to diversion.
     * 
     * @return Vector of multipliers (0-1)
     */
    std::vector<double> getFlowDistributionMultiplier() const;
    
    /**
     * @brief Get current diversion state summary
     */
    struct DiversionSummary {
        double total_pressure_drop;
        int clusters_diverted;
        int clusters_open;
        double diverter_mass_remaining;
        double degradation_fraction;
    };
    
    DiversionSummary getSummary() const;
    
    // Output
    void writeReport(std::ostream& os) const;
    
private:
    bool use_particulate_;
    bool use_ball_sealers_;
    bool use_chemical_;
    
    std::unique_ptr<DiverterTransport> particulate_model_;
    std::unique_ptr<BallSealerModel> ball_model_;
    std::unique_ptr<ChemicalDiversion> chemical_model_;
    
    double current_time_;
};

} // namespace FSRM

#endif // DIVERSION_MODELING_HPP
