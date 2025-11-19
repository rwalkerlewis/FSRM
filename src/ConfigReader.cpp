#include "ConfigReader.hpp"
#include <iostream>
#include <algorithm>
#include <cctype>

namespace ResSim {

ConfigReader::ConfigReader() {}

bool ConfigReader::loadFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open configuration file: " << filename << std::endl;
        return false;
    }
    
    std::string current_section;
    std::string line;
    int line_num = 0;
    
    while (std::getline(file, line)) {
        line_num++;
        line = trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == ';') {
            continue;
        }
        
        // Check for section header [SECTION]
        if (line[0] == '[' && line.back() == ']') {
            current_section = line.substr(1, line.length() - 2);
            current_section = trim(current_section);
            continue;
        }
        
        // Parse key = value
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) {
            std::cerr << "Warning: Invalid line " << line_num << ": " << line << std::endl;
            continue;
        }
        
        std::string key = trim(line.substr(0, eq_pos));
        std::string value = trim(line.substr(eq_pos + 1));
        
        // Remove inline comments
        size_t comment_pos = value.find('#');
        if (comment_pos != std::string::npos) {
            value = trim(value.substr(0, comment_pos));
        }
        
        if (current_section.empty()) {
            std::cerr << "Warning: Key without section at line " << line_num << std::endl;
            continue;
        }
        
        data[current_section][key] = value;
    }
    
    file.close();
    return true;
}

std::string ConfigReader::trim(const std::string& str) const {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> ConfigReader::split(const std::string& str, char delim) const {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string item;
    
    while (std::getline(ss, item, delim)) {
        item = trim(item);
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    
    return result;
}

std::string ConfigReader::getString(const std::string& section, const std::string& key,
                                   const std::string& default_val) const {
    auto sec_it = data.find(section);
    if (sec_it == data.end()) return default_val;
    
    auto key_it = sec_it->second.find(key);
    if (key_it == sec_it->second.end()) return default_val;
    
    return key_it->second;
}

int ConfigReader::getInt(const std::string& section, const std::string& key,
                        int default_val) const {
    std::string val = getString(section, key);
    if (val.empty()) return default_val;
    
    try {
        return std::stoi(val);
    } catch (...) {
        return default_val;
    }
}

double ConfigReader::getDouble(const std::string& section, const std::string& key,
                               double default_val) const {
    std::string val = getString(section, key);
    if (val.empty()) return default_val;
    
    try {
        return std::stod(val);
    } catch (...) {
        return default_val;
    }
}

bool ConfigReader::getBool(const std::string& section, const std::string& key,
                          bool default_val) const {
    std::string val = getString(section, key);
    if (val.empty()) return default_val;
    
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    
    if (val == "true" || val == "yes" || val == "1" || val == "on") return true;
    if (val == "false" || val == "no" || val == "0" || val == "off") return false;
    
    return default_val;
}

std::vector<double> ConfigReader::getDoubleArray(const std::string& section,
                                                 const std::string& key) const {
    std::vector<double> result;
    std::string val = getString(section, key);
    if (val.empty()) return result;
    
    auto tokens = split(val, ',');
    for (const auto& token : tokens) {
        try {
            result.push_back(std::stod(token));
        } catch (...) {
            std::cerr << "Warning: Cannot parse '" << token << "' as double" << std::endl;
        }
    }
    
    return result;
}

bool ConfigReader::hasSection(const std::string& section) const {
    return data.find(section) != data.end();
}

bool ConfigReader::hasKey(const std::string& section, const std::string& key) const {
    auto sec_it = data.find(section);
    if (sec_it == data.end()) return false;
    return sec_it->second.find(key) != sec_it->second.end();
}

std::vector<std::string> ConfigReader::getSections() const {
    std::vector<std::string> sections;
    for (const auto& pair : data) {
        sections.push_back(pair.first);
    }
    return sections;
}

std::vector<std::string> ConfigReader::getKeys(const std::string& section) const {
    std::vector<std::string> keys;
    auto sec_it = data.find(section);
    if (sec_it != data.end()) {
        for (const auto& pair : sec_it->second) {
            keys.push_back(pair.first);
        }
    }
    return keys;
}

bool ConfigReader::parseSimulationConfig(SimulationConfig& config) {
    if (!hasSection("SIMULATION")) return false;
    
    config.start_time = getDouble("SIMULATION", "start_time", 0.0);
    config.end_time = getDouble("SIMULATION", "end_time", 1.0);
    config.dt_initial = getDouble("SIMULATION", "dt_initial", 0.01);
    config.dt_min = getDouble("SIMULATION", "dt_min", 1e-10);
    config.dt_max = getDouble("SIMULATION", "dt_max", 1.0);
    
    config.max_timesteps = getInt("SIMULATION", "max_timesteps", 10000);
    config.output_frequency = getInt("SIMULATION", "output_frequency", 10);
    
    config.enable_adaptive_timestepping = getBool("SIMULATION", "adaptive_timestepping", true);
    config.enable_checkpointing = getBool("SIMULATION", "enable_checkpointing", true);
    config.checkpoint_frequency = getInt("SIMULATION", "checkpoint_frequency", 100);
    
    config.rtol = getDouble("SIMULATION", "rtol", 1e-6);
    config.atol = getDouble("SIMULATION", "atol", 1e-8);
    config.max_nonlinear_iterations = getInt("SIMULATION", "max_nonlinear_iterations", 50);
    config.max_linear_iterations = getInt("SIMULATION", "max_linear_iterations", 1000);
    
    // Physics flags
    config.enable_geomechanics = getBool("SIMULATION", "enable_geomechanics", false);
    config.enable_thermal = getBool("SIMULATION", "enable_thermal", false);
    config.enable_fractures = getBool("SIMULATION", "enable_fractures", false);
    config.enable_particle_transport = getBool("SIMULATION", "enable_particle_transport", false);
    config.enable_faults = getBool("SIMULATION", "enable_faults", false);
    config.enable_tidal_forces = getBool("SIMULATION", "enable_tidal_forces", false);
    
    // Fluid model
    std::string fluid_model = getString("SIMULATION", "fluid_model", "SINGLE_COMPONENT");
    if (fluid_model == "SINGLE_COMPONENT") config.fluid_model = FluidModelType::SINGLE_COMPONENT;
    else if (fluid_model == "BLACK_OIL") config.fluid_model = FluidModelType::BLACK_OIL;
    else if (fluid_model == "COMPOSITIONAL") config.fluid_model = FluidModelType::COMPOSITIONAL;
    
    // Solid model
    std::string solid_model = getString("SIMULATION", "solid_model", "ELASTIC");
    if (solid_model == "ELASTIC") config.solid_model = SolidModelType::ELASTIC;
    else if (solid_model == "VISCOELASTIC") config.solid_model = SolidModelType::VISCOELASTIC;
    else if (solid_model == "POROELASTIC") config.solid_model = SolidModelType::POROELASTIC;
    
    // Output format
    config.output_format = getString("SIMULATION", "output_format", "VTK");
    
    return true;
}

bool ConfigReader::parseGridConfig(GridConfig& config) {
    if (!hasSection("GRID")) return false;
    
    config.nx = getInt("GRID", "nx", 10);
    config.ny = getInt("GRID", "ny", 10);
    config.nz = getInt("GRID", "nz", 10);
    
    config.Lx = getDouble("GRID", "Lx", 1000.0);
    config.Ly = getDouble("GRID", "Ly", 1000.0);
    config.Lz = getDouble("GRID", "Lz", 100.0);
    
    config.use_unstructured = getBool("GRID", "use_unstructured", false);
    config.mesh_file = getString("GRID", "mesh_file", "");
    
    return true;
}

bool ConfigReader::parseMaterialProperties(std::vector<MaterialProperties>& props) {
    // Look for ROCK sections (ROCK1, ROCK2, etc.) or single ROCK section
    std::vector<std::string> rock_sections;
    
    for (const auto& section : getSections()) {
        if (section.find("ROCK") == 0) {
            rock_sections.push_back(section);
        }
    }
    
    if (rock_sections.empty()) return false;
    
    props.clear();
    
    for (const auto& section : rock_sections) {
        MaterialProperties mat;
        
        // Rock properties
        mat.porosity = getDouble(section, "porosity", 0.2);
        mat.permeability_x = getDouble(section, "permeability_x", 100.0) * 1e-15; // mD to m²
        mat.permeability_y = getDouble(section, "permeability_y", 100.0) * 1e-15;
        mat.permeability_z = getDouble(section, "permeability_z", 10.0) * 1e-15;
        mat.compressibility = getDouble(section, "compressibility", 1e-9); // 1/Pa
        
        // Mechanical properties
        mat.youngs_modulus = getDouble(section, "youngs_modulus", 10e9); // Pa
        mat.poisson_ratio = getDouble(section, "poisson_ratio", 0.25);
        mat.density = getDouble(section, "density", 2500.0); // kg/m³
        mat.biot_coefficient = getDouble(section, "biot_coefficient", 1.0);
        
        // Viscoelastic properties
        mat.relaxation_time = getDouble(section, "relaxation_time", 1e6); // seconds
        mat.viscosity = getDouble(section, "viscosity", 1e19); // Pa·s
        
        // Thermal properties
        mat.thermal_conductivity = getDouble(section, "thermal_conductivity", 2.5); // W/(m·K)
        mat.heat_capacity = getDouble(section, "heat_capacity", 900.0); // J/(kg·K)
        mat.thermal_expansion = getDouble(section, "thermal_expansion", 1e-5); // 1/K
        
        // Fracture properties
        mat.fracture_toughness = getDouble(section, "fracture_toughness", 1e6); // Pa·m^0.5
        mat.fracture_energy = getDouble(section, "fracture_energy", 100.0); // J/m²
        
        props.push_back(mat);
    }
    
    return true;
}

bool ConfigReader::parseFluidProperties(FluidProperties& props) {
    if (!hasSection("FLUID")) return false;
    
    // Single phase properties
    props.density = getDouble("FLUID", "density", 1000.0); // kg/m³
    props.viscosity = getDouble("FLUID", "viscosity", 0.001); // Pa·s
    props.compressibility = getDouble("FLUID", "compressibility", 1e-9); // 1/Pa
    
    // Black oil properties
    props.oil_density_std = getDouble("FLUID", "oil_density_std", 800.0);
    props.gas_density_std = getDouble("FLUID", "gas_density_std", 1.0);
    props.water_density_std = getDouble("FLUID", "water_density_std", 1000.0);
    
    props.oil_viscosity = getDouble("FLUID", "oil_viscosity", 0.005);
    props.gas_viscosity = getDouble("FLUID", "gas_viscosity", 0.00001);
    props.water_viscosity = getDouble("FLUID", "water_viscosity", 0.001);
    
    props.solution_GOR = getDouble("FLUID", "solution_GOR", 100.0); // scf/stb
    
    // Compositional properties
    props.component_mw = getDoubleArray("FLUID", "component_mw");
    props.component_Tc = getDoubleArray("FLUID", "component_Tc");
    props.component_Pc = getDoubleArray("FLUID", "component_Pc");
    props.component_omega = getDoubleArray("FLUID", "component_omega");
    
    return true;
}

void ConfigReader::generateTemplate(const std::string& filename) {
    std::ofstream file(filename);
    
    file << "# ReservoirSim Configuration File\n";
    file << "# All units in SI unless otherwise specified\n";
    file << "#\n";
    file << "# Lines starting with # or ; are comments\n";
    file << "# Format: key = value\n\n";
    
    file << "[SIMULATION]\n";
    file << "# Time parameters (seconds)\n";
    file << "start_time = 0.0\n";
    file << "end_time = 86400.0                    # 1 day\n";
    file << "dt_initial = 3600.0                   # 1 hour\n";
    file << "dt_min = 1.0\n";
    file << "dt_max = 86400.0\n";
    file << "max_timesteps = 10000\n\n";
    
    file << "# Output control\n";
    file << "output_frequency = 10\n";
    file << "output_format = VTK                   # VTK, HDF5, or ECLIPSE\n";
    file << "enable_checkpointing = true\n";
    file << "checkpoint_frequency = 100\n\n";
    
    file << "# Solver tolerances\n";
    file << "rtol = 1.0e-6                         # Relative tolerance\n";
    file << "atol = 1.0e-8                         # Absolute tolerance\n";
    file << "max_nonlinear_iterations = 50\n";
    file << "max_linear_iterations = 1000\n\n";
    
    file << "# Physics models\n";
    file << "fluid_model = SINGLE_COMPONENT        # SINGLE_COMPONENT, BLACK_OIL, COMPOSITIONAL\n";
    file << "solid_model = ELASTIC                 # ELASTIC, VISCOELASTIC, POROELASTIC\n\n";
    
    file << "# Physics flags\n";
    file << "enable_geomechanics = false\n";
    file << "enable_thermal = false\n";
    file << "enable_fractures = false\n";
    file << "enable_particle_transport = false\n";
    file << "enable_faults = false\n";
    file << "enable_tidal_forces = false\n";
    file << "adaptive_timestepping = true\n\n";
    
    file << "[GRID]\n";
    file << "# Grid dimensions\n";
    file << "nx = 20\n";
    file << "ny = 20\n";
    file << "nz = 5\n\n";
    
    file << "# Domain size (meters)\n";
    file << "Lx = 1000.0\n";
    file << "Ly = 1000.0\n";
    file << "Lz = 100.0\n\n";
    
    file << "# Mesh type\n";
    file << "use_unstructured = false\n";
    file << "mesh_file =                           # Path to mesh file if unstructured\n\n";
    
    file << "[ROCK]\n";
    file << "# Rock properties\n";
    file << "porosity = 0.20                       # Fraction (0-1)\n";
    file << "permeability_x = 100.0                # milliDarcy\n";
    file << "permeability_y = 100.0\n";
    file << "permeability_z = 10.0\n";
    file << "compressibility = 1.0e-9              # 1/Pa\n\n";
    
    file << "# Mechanical properties\n";
    file << "youngs_modulus = 10.0e9               # Pa (10 GPa)\n";
    file << "poisson_ratio = 0.25                  # Dimensionless\n";
    file << "density = 2500.0                      # kg/m³\n";
    file << "biot_coefficient = 1.0                # Dimensionless\n\n";
    
    file << "# Viscoelastic properties (if enabled)\n";
    file << "relaxation_time = 1.0e6               # seconds\n";
    file << "viscosity = 1.0e19                    # Pa·s\n\n";
    
    file << "# Thermal properties\n";
    file << "thermal_conductivity = 2.5            # W/(m·K)\n";
    file << "heat_capacity = 900.0                 # J/(kg·K)\n";
    file << "thermal_expansion = 1.0e-5            # 1/K\n\n";
    
    file << "# Fracture properties\n";
    file << "fracture_toughness = 1.0e6            # Pa·m^0.5\n";
    file << "fracture_energy = 100.0               # J/m²\n\n";
    
    file << "[FLUID]\n";
    file << "# Single-phase fluid properties\n";
    file << "density = 1000.0                      # kg/m³\n";
    file << "viscosity = 0.001                     # Pa·s (1 cP)\n";
    file << "compressibility = 1.0e-9              # 1/Pa\n\n";
    
    file << "# Black oil model properties\n";
    file << "oil_density_std = 800.0               # kg/m³ at standard conditions\n";
    file << "gas_density_std = 1.0\n";
    file << "water_density_std = 1000.0\n\n";
    
    file << "oil_viscosity = 0.005                 # Pa·s (5 cP)\n";
    file << "gas_viscosity = 0.00001               # Pa·s (0.01 cP)\n";
    file << "water_viscosity = 0.001               # Pa·s (1 cP)\n\n";
    
    file << "solution_GOR = 100.0                  # scf/stb\n\n";
    
    file << "# Compositional model (comma-separated arrays)\n";
    file << "component_mw = 16.04, 30.07, 44.10    # g/mol (CH4, C2H6, C3H8)\n";
    file << "component_Tc = 190.6, 305.4, 369.8    # K\n";
    file << "component_Pc = 4.60e6, 4.88e6, 4.25e6 # Pa\n";
    file << "component_omega = 0.011, 0.099, 0.152 # Acentric factor\n\n";
    
    file << "# Multiple rock types can be defined as [ROCK1], [ROCK2], etc.\n";
    file << "# Use Eclipse INCLUDE files for complex property distributions\n";
    
    file.close();
}

} // namespace ResSim
