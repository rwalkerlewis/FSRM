#ifndef CONFIG_READER_HPP
#define CONFIG_READER_HPP

#include "ReservoirSim.hpp"
#include "FluidModel.hpp"
#include "MaterialModel.hpp"
#include "FaultModel.hpp"
#include <string>
#include <map>
#include <vector>
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
    
    bool parseDynamicsConfig(DynamicsConfig& config);
    bool parseSeismicityConfig(SeismicityConfig& config);
    bool parseOutputConfig(OutputConfig& config);
    
    std::vector<WellConfig> parseWells();
    std::vector<FractureConfig> parseFractures();
    std::vector<FaultConfig> parseFaults();
    bool parseParticleProperties(ParticleConfig& config);
    std::vector<BoundaryCondition> parseBoundaryConditions();
    std::vector<InitialCondition> parseInitialConditions();
    
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
    
    std::string trim(const std::string& str) const;
    std::vector<std::string> split(const std::string& str, char delim) const;
    
    std::vector<double> parseDoubleArray(const std::string& value) const;
    std::vector<int> parseIntArray(const std::string& value) const;
    std::vector<std::string> parseStringArray(const std::string& value) const;
};

} // namespace FSRM

#endif // CONFIG_READER_HPP
