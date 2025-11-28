#ifndef CONFIG_READER_HPP
#define CONFIG_READER_HPP

#include "ReservoirSim.hpp"
#include "FluidModel.hpp"
#include "MaterialModel.hpp"
#include "FaultModel.hpp"
#include "UnitSystem.hpp"
#include <string>
#include <map>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <memory>

namespace FSRM {

/**
 * @brief Enhanced INI-style configuration reader
 * 
 * Supports all physics models and can configure the entire simulation
 * from a single text file without any code changes.
 */
class ConfigReader {
public:
    // =========================================================================
    // Nested Struct Definitions (must be before method declarations that use them)
    // =========================================================================
    
    struct DynamicsConfig {
        bool enable_dynamics;
        bool use_static_triggering;
        double trigger_threshold;
        double event_duration;
        double damping_alpha;
        double damping_beta;
        double quality_factor;
        bool enable_dynamic_permeability;
        double permeability_sensitivity;
        double permeability_recovery_time;
    };
    
    struct SeismicityConfig {
        bool enable_seismicity;
        std::string friction_law;
        double nucleation_size;
        double seismic_slip_rate;
        double b_value;
        bool enable_aftershocks;
        bool enable_stress_transfer;
    };
    
    struct OutputConfig {
        std::string format;           // VTK, HDF5, ECLIPSE
        std::string base_path;
        int frequency;
        bool write_pressure;
        bool write_displacement;
        bool write_stress;
        bool write_strain;
        bool write_velocity;
        bool write_permeability;
        bool write_temperature;
        bool write_saturation;
        bool write_fault_slip;
        bool write_seismic_catalog;
        
        // Output unit preferences
        std::map<std::string, std::string> output_units;
    };
    
    struct WellConfig {
        std::string name;
        std::string type;  // PRODUCER, INJECTOR
        int i, j, k;
        std::string control_mode;  // RATE, BHP, THP
        double target_value;
        double max_rate;
        double min_bhp;
        double diameter;
        double skin;
    };
    
    struct FractureConfig {
        std::string type;  // NATURAL, HYDRAULIC
        std::vector<double> location;
        double aperture;
        double permeability;
        double toughness;
        double energy;
        bool enable_propagation;
        bool enable_proppant;
    };
    
    struct FaultConfig {
        std::string name;
        std::vector<double> strike;
        std::vector<double> dip;
        double length;
        double width;
        double static_friction;
        double dynamic_friction;
        double cohesion;
        bool use_rate_state;
        double a_parameter;
        double b_parameter;
        double Dc_parameter;
        
        // Split node options
        bool use_split_nodes;
        std::string split_node_method;      // PENALTY, LAGRANGE, NITSCHE, etc.
        std::string traction_type;          // PRESCRIBED, FRICTION_DEPENDENT, COHESIVE_ZONE
        double penalty_normal;
        double penalty_tangent;
        double prescribed_traction_normal;
        double prescribed_traction_strike;
        double prescribed_traction_dip;
        double cohesive_strength;
        double critical_opening;
        double critical_slip;
        bool allow_separation;
    };
    
    struct ParticleConfig {
        double diameter;
        double density;
        double concentration;
        double diffusivity;
        bool enable_settling;
        bool enable_bridging;
    };
    
    struct BoundaryCondition {
        std::string type;  // DIRICHLET, NEUMANN, ROBIN
        std::string field; // PRESSURE, TEMPERATURE, DISPLACEMENT
        std::string location; // XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
        double value;
        double gradient;  // For Neumann
    };
    
    struct InitialCondition {
        std::string field;
        std::string distribution; // UNIFORM, GRADIENT, FILE
        double value;
        std::vector<double> gradient;
        std::string file;
    };
    
    struct ValidationResult {
        bool valid;
        std::vector<std::string> errors;
        std::vector<std::string> warnings;
    };
    
    // =========================================================================
    // Adaptive Mesh Refinement Configuration
    // =========================================================================
    
    struct AMRFeatureWell {
        std::string name;
        double x, y, z;
        double radius;
        int level;
    };
    
    struct AMRFeatureFault {
        std::string name;
        std::vector<std::array<double, 3>> trace;
        double width;
        int level;
    };
    
    struct AMRFeatureBox {
        std::string name;
        double xmin, xmax;
        double ymin, ymax;
        double zmin, zmax;
        int level;
    };
    
    struct AMRConfig {
        bool enabled;
        std::string method;           // plex_refine, plex_adapt, forest_p4est, pragmatic, mmg
        int adapt_every;              // Time steps between adaptations
        bool adapt_on_change;
        
        // Criterion
        std::string criterion_type;   // gradient, hessian, residual, jump, feature, combined
        double weight_pressure;
        double weight_saturation;
        double weight_velocity;
        double weight_temperature;
        
        // Strategy
        std::string strategy;         // fixed_fraction, threshold, equilibration
        double refine_fraction;
        double coarsen_fraction;
        double refine_threshold;
        double coarsen_threshold;
        
        // Limits
        int max_level;
        int min_level;
        int max_cells;
        int min_cells;
        double min_cell_size;
        double max_cell_size;
        
        // Quality
        double min_quality;
        double max_aspect_ratio;
        
        // Features
        bool preserve_boundaries;
        bool preserve_wells;
        bool preserve_faults;
        int buffer_layers;
        std::vector<AMRFeatureWell> wells;
        std::vector<AMRFeatureFault> faults;
        std::vector<AMRFeatureBox> boxes;
        
        // Transfer
        std::string transfer_method;  // interpolation, projection, injection, conservative
        bool conservative;
        
        // Load balancing
        bool balance_enabled;
        std::string balance_method;   // parmetis, ptscotch, simple
        double rebalance_threshold;
        
        // Output
        bool write_mesh;
        bool write_errors;
        std::string output_format;
        std::string output_prefix;
        
        // Default constructor with sensible defaults
        AMRConfig() : enabled(false), method("plex_adapt"), adapt_every(5),
                     adapt_on_change(true), criterion_type("gradient"),
                     weight_pressure(1.0), weight_saturation(1.0),
                     weight_velocity(0.5), weight_temperature(0.5),
                     strategy("fixed_fraction"), refine_fraction(0.2),
                     coarsen_fraction(0.05), refine_threshold(0.5),
                     coarsen_threshold(0.05), max_level(5), min_level(0),
                     max_cells(1000000), min_cells(100),
                     min_cell_size(1e-6), max_cell_size(1e6),
                     min_quality(0.1), max_aspect_ratio(10.0),
                     preserve_boundaries(true), preserve_wells(true),
                     preserve_faults(true), buffer_layers(1),
                     transfer_method("interpolation"), conservative(true),
                     balance_enabled(true), balance_method("parmetis"),
                     rebalance_threshold(0.2), write_mesh(true),
                     write_errors(true), output_format("vtk"),
                     output_prefix("amr_output") {}
    };
    
    // =========================================================================
    // Constructor/Destructor
    // =========================================================================
    
    ConfigReader();
    ~ConfigReader() = default;
    
    // Load configuration file
    bool loadFile(const std::string& filename);
    
    // =========================================================================
    // Basic Parsing Methods
    // =========================================================================
    
    bool parseSimulationConfig(SimulationConfig& config);
    bool parseGridConfig(GridConfig& config);
    bool parseMaterialProperties(std::vector<MaterialProperties>& props);
    bool parseFluidProperties(FluidProperties& props);
    
    // =========================================================================
    // Extended Parsing Methods for New Libraries
    // =========================================================================
    
    std::unique_ptr<FluidModelBase> parseFluidModel();
    std::vector<RockProperties> parseRockProperties();
    std::unique_ptr<MaterialModelBase> parseMaterialModel(const std::string& section);
    std::unique_ptr<FaultNetwork> parseFaultNetwork();
    std::vector<FaultConfig> parseFaultsExtended();
    std::vector<std::unique_ptr<SplitNodeFault>> parseSplitNodeFaults();
    
    bool parseDynamicsConfig(DynamicsConfig& config);
    bool parseSeismicityConfig(SeismicityConfig& config);
    bool parseOutputConfig(OutputConfig& config);
    
    std::vector<WellConfig> parseWells();
    std::vector<FractureConfig> parseFractures();
    std::vector<FaultConfig> parseFaults();
    bool parseParticleProperties(ParticleConfig& config);
    std::vector<BoundaryCondition> parseBoundaryConditions();
    std::vector<InitialCondition> parseInitialConditions();
    
    /**
     * @brief Parse Adaptive Mesh Refinement configuration
     * @param config Output AMR configuration
     * @return true if AMR section exists and was parsed successfully
     */
    bool parseAMRConfig(AMRConfig& config);
    
    // =========================================================================
    // Value Accessors
    // =========================================================================
    
    std::string getString(const std::string& section, const std::string& key, 
                         const std::string& default_val = "") const;
    int getInt(const std::string& section, const std::string& key, 
              int default_val = 0) const;
    double getDouble(const std::string& section, const std::string& key, 
                    double default_val = 0.0) const;
    bool getBool(const std::string& section, const std::string& key, 
                bool default_val = false) const;
    std::vector<double> getDoubleArray(const std::string& section, 
                                       const std::string& key) const;
    
    // =========================================================================
    // Unit-Aware Value Accessors (converts to SI base units)
    // =========================================================================
    
    /**
     * @brief Get double value with automatic unit conversion to SI
     * @param section Config section
     * @param key Config key
     * @param default_val Default value (in SI units)
     * @param default_unit Default unit if no unit specified in value
     * @return Value converted to SI base units (m, kg, s)
     */
    double getDoubleWithUnit(const std::string& section, const std::string& key,
                            double default_val = 0.0, 
                            const std::string& default_unit = "") const;
    
    /**
     * @brief Get array of doubles with automatic unit conversion
     */
    std::vector<double> getDoubleArrayWithUnit(const std::string& section,
                                               const std::string& key,
                                               const std::string& default_unit = "") const;
    
    // =========================================================================
    // Section/Key Query Methods
    // =========================================================================
    
    bool hasSection(const std::string& section) const;
    bool hasKey(const std::string& section, const std::string& key) const;
    std::vector<std::string> getSections() const;
    std::vector<std::string> getKeys(const std::string& section) const;
    
    // =========================================================================
    // Template Generation
    // =========================================================================
    
    static void generateTemplate(const std::string& filename);
    static void generateCompleteTemplate(const std::string& filename);
    
    // =========================================================================
    // Utility Methods
    // =========================================================================
    
    std::map<std::string, std::string> getSectionData(const std::string& section) const;
    std::vector<std::string> getSectionsMatching(const std::string& prefix) const;
    bool mergeFile(const std::string& filename);
    ValidationResult validate() const;
    
private:
    std::map<std::string, std::map<std::string, std::string>> data;
    UnitSystem unit_system_;
    
    std::string trim(const std::string& str) const;
    std::vector<std::string> split(const std::string& str, char delim) const;
    
    std::vector<double> parseDoubleArray(const std::string& value) const;
    std::vector<int> parseIntArray(const std::string& value) const;
    std::vector<std::string> parseStringArray(const std::string& value) const;
};

} // namespace FSRM

#endif // CONFIG_READER_HPP
