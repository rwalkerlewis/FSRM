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
