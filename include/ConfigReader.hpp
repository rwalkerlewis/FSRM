#ifndef CONFIG_READER_HPP
#define CONFIG_READER_HPP

#include "ReservoirSim.hpp"
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

namespace ResSim {

// Simple INI-style configuration reader
class ConfigReader {
public:
    ConfigReader();
    ~ConfigReader() = default;
    
    // Load configuration file
    bool loadFile(const std::string& filename);
    
    // Parse different sections
    bool parseSimulationConfig(SimulationConfig& config);
    bool parseGridConfig(GridConfig& config);
    bool parseMaterialProperties(std::vector<MaterialProperties>& props);
    bool parseFluidProperties(FluidProperties& props);
    
    // Parse wells
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
    std::vector<WellConfig> parseWells();
    
    // Parse fractures
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
    std::vector<FractureConfig> parseFractures();
    
    // Parse faults
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
    std::vector<FaultConfig> parseFaults();
    
    // Parse particle/proppant properties
    struct ParticleConfig {
        double diameter;
        double density;
        double concentration;
        double diffusivity;
        bool enable_settling;
        bool enable_bridging;
    };
    bool parseParticleProperties(ParticleConfig& config);
    
    // Parse boundary conditions
    struct BoundaryCondition {
        std::string type;  // DIRICHLET, NEUMANN, ROBIN
        std::string field; // PRESSURE, TEMPERATURE, DISPLACEMENT
        std::string location; // XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
        double value;
        double gradient;  // For Neumann
    };
    std::vector<BoundaryCondition> parseBoundaryConditions();
    
    // Parse initial conditions
    struct InitialCondition {
        std::string field;
        std::string distribution; // UNIFORM, GRADIENT, FILE
        double value;
        std::vector<double> gradient;
        std::string file;
    };
    std::vector<InitialCondition> parseInitialConditions();
    
    // Get values by key
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
    
    // Check if section/key exists
    bool hasSection(const std::string& section) const;
    bool hasKey(const std::string& section, const std::string& key) const;
    
    // List all sections
    std::vector<std::string> getSections() const;
    
    // Get all keys in a section
    std::vector<std::string> getKeys(const std::string& section) const;
    
    // Generate template configuration file
    static void generateTemplate(const std::string& filename);
    
private:
    std::map<std::string, std::map<std::string, std::string>> data;
    
    std::string trim(const std::string& str) const;
    std::vector<std::string> split(const std::string& str, char delim) const;
};

} // namespace ResSim

#endif // CONFIG_READER_HPP
