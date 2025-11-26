#include "ConfigReader.hpp"
#include <iostream>
#include <algorithm>
#include <cctype>

namespace FSRM {

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

std::vector<ConfigReader::WellConfig> ConfigReader::parseWells() {
    std::vector<WellConfig> wells;
    
    for (const auto& section : getSections()) {
        if (section.find("WELL") == 0) {
            WellConfig well;
            well.name = getString(section, "name", section);
            well.type = getString(section, "type", "PRODUCER");
            well.i = getInt(section, "i", 0);
            well.j = getInt(section, "j", 0);
            well.k = getInt(section, "k", 0);
            well.control_mode = getString(section, "control_mode", "RATE");
            well.target_value = getDouble(section, "target_value", 100.0);
            well.max_rate = getDouble(section, "max_rate", 1000.0);
            well.min_bhp = getDouble(section, "min_bhp", 5e6);
            well.diameter = getDouble(section, "diameter", 0.2);
            well.skin = getDouble(section, "skin", 0.0);
            wells.push_back(well);
        }
    }
    
    return wells;
}

std::vector<ConfigReader::FractureConfig> ConfigReader::parseFractures() {
    std::vector<FractureConfig> fractures;
    
    for (const auto& section : getSections()) {
        if (section.find("FRACTURE") == 0) {
            FractureConfig frac;
            frac.type = getString(section, "type", "NATURAL");
            frac.location = getDoubleArray(section, "location");
            frac.aperture = getDouble(section, "aperture", 1e-4);
            frac.permeability = getDouble(section, "permeability", 1e-12);
            frac.toughness = getDouble(section, "toughness", 1e6);
            frac.energy = getDouble(section, "energy", 100.0);
            frac.enable_propagation = getBool(section, "enable_propagation", false);
            frac.enable_proppant = getBool(section, "enable_proppant", false);
            fractures.push_back(frac);
        }
    }
    
    return fractures;
}

std::vector<ConfigReader::FaultConfig> ConfigReader::parseFaults() {
    std::vector<FaultConfig> faults;
    
    for (const auto& section : getSections()) {
        if (section.find("FAULT") == 0) {
            FaultConfig fault;
            fault.name = getString(section, "name", section);
            fault.strike = getDoubleArray(section, "strike");
            fault.dip = getDoubleArray(section, "dip");
            fault.length = getDouble(section, "length", 1000.0);
            fault.width = getDouble(section, "width", 500.0);
            fault.static_friction = getDouble(section, "static_friction", 0.6);
            fault.dynamic_friction = getDouble(section, "dynamic_friction", 0.4);
            fault.cohesion = getDouble(section, "cohesion", 1e6);
            fault.use_rate_state = getBool(section, "use_rate_state", false);
            fault.a_parameter = getDouble(section, "a_parameter", 0.01);
            fault.b_parameter = getDouble(section, "b_parameter", 0.015);
            fault.Dc_parameter = getDouble(section, "Dc_parameter", 0.001);
            faults.push_back(fault);
        }
    }
    
    return faults;
}

bool ConfigReader::parseParticleProperties(ParticleConfig& config) {
    if (!hasSection("PARTICLE")) return false;
    
    config.diameter = getDouble("PARTICLE", "diameter", 0.0003);
    config.density = getDouble("PARTICLE", "density", 2650.0);
    config.concentration = getDouble("PARTICLE", "concentration", 0.0);
    config.diffusivity = getDouble("PARTICLE", "diffusivity", 1e-9);
    config.enable_settling = getBool("PARTICLE", "enable_settling", true);
    config.enable_bridging = getBool("PARTICLE", "enable_bridging", false);
    
    return true;
}

std::vector<ConfigReader::BoundaryCondition> ConfigReader::parseBoundaryConditions() {
    std::vector<BoundaryCondition> bcs;
    
    for (const auto& section : getSections()) {
        if (section.find("BC") == 0) {
            BoundaryCondition bc;
            bc.type = getString(section, "type", "DIRICHLET");
            bc.field = getString(section, "field", "PRESSURE");
            bc.location = getString(section, "location", "XMIN");
            bc.value = getDouble(section, "value", 0.0);
            bc.gradient = getDouble(section, "gradient", 0.0);
            bcs.push_back(bc);
        }
    }
    
    return bcs;
}

std::vector<ConfigReader::InitialCondition> ConfigReader::parseInitialConditions() {
    std::vector<InitialCondition> ics;
    
    for (const auto& section : getSections()) {
        if (section.find("IC") == 0) {
            InitialCondition ic;
            ic.field = getString(section, "field", "PRESSURE");
            ic.distribution = getString(section, "distribution", "UNIFORM");
            ic.value = getDouble(section, "value", 1e7);
            ic.gradient = getDoubleArray(section, "gradient");
            ic.file = getString(section, "file", "");
            ics.push_back(ic);
        }
    }
    
    return ics;
}

void ConfigReader::generateTemplate(const std::string& filename) {
    std::ofstream file(filename);
    
    file << "# FSRM Configuration File\n";
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
    
    file << "[WELL1]\n";
    file << "# Well configuration\n";
    file << "name = PROD1\n";
    file << "type = PRODUCER                       # PRODUCER, INJECTOR\n";
    file << "i = 10                                # Grid indices (0-based)\n";
    file << "j = 10\n";
    file << "k = 5\n";
    file << "control_mode = RATE                   # RATE, BHP, THP\n";
    file << "target_value = 0.01                   # m³/s (rate) or Pa (pressure)\n";
    file << "max_rate = 0.1                        # m³/s\n";
    file << "min_bhp = 5.0e6                       # Pa (5 MPa)\n";
    file << "diameter = 0.2                        # meters (8 inch)\n";
    file << "skin = 0.0                            # Skin factor\n\n";
    
    file << "[WELL2]\n";
    file << "name = INJ1\n";
    file << "type = INJECTOR\n";
    file << "i = 1\n";
    file << "j = 1\n";
    file << "k = 5\n";
    file << "control_mode = RATE\n";
    file << "target_value = 0.015                  # m³/s\n";
    file << "max_rate = 0.2\n";
    file << "min_bhp = 50.0e6                      # Pa (50 MPa max injection pressure)\n";
    file << "diameter = 0.15\n";
    file << "skin = 0.0\n\n";
    
    file << "[FRACTURE1]\n";
    file << "# Hydraulic or natural fracture\n";
    file << "type = HYDRAULIC                      # NATURAL, HYDRAULIC\n";
    file << "location = 500.0, 500.0, 50.0, 0.0, 1.0, 0.0  # x,y,z center + normal vector\n";
    file << "aperture = 0.001                      # meters (1 mm)\n";
    file << "permeability = 1.0e-10                # m² (from cubic law)\n";
    file << "toughness = 1.0e6                     # Pa·m^0.5\n";
    file << "energy = 100.0                        # J/m²\n";
    file << "enable_propagation = true\n";
    file << "enable_proppant = true\n\n";
    
    file << "[FAULT1]\n";
    file << "# Pre-existing fault\n";
    file << "name = MAIN_FAULT\n";
    file << "strike = 1.0, 0.0, 0.0                # Direction vector (E-W)\n";
    file << "dip = 0.0, 0.707, 0.707               # 45 degree dip\n";
    file << "length = 2000.0                       # meters\n";
    file << "width = 1500.0                        # meters\n";
    file << "static_friction = 0.6                 # Byerlee's law\n";
    file << "dynamic_friction = 0.4\n";
    file << "cohesion = 1.0e6                      # Pa (1 MPa)\n";
    file << "use_rate_state = true                 # Rate-and-state friction\n";
    file << "a_parameter = 0.010                   # Rate-state a\n";
    file << "b_parameter = 0.015                   # Rate-state b (b>a: velocity weakening)\n";
    file << "Dc_parameter = 0.001                  # Rate-state critical slip distance (m)\n\n";
    
    file << "[PARTICLE]\n";
    file << "# Proppant or tracer properties\n";
    file << "diameter = 0.0003                     # meters (300 microns, 40/70 mesh)\n";
    file << "density = 2650.0                      # kg/m³ (ceramic proppant)\n";
    file << "concentration = 0.0                   # kg/m³ (set to 0 initially)\n";
    file << "diffusivity = 1.0e-9                  # m²/s\n";
    file << "enable_settling = true                # Gravitational settling\n";
    file << "enable_bridging = false               # Proppant bridging in fractures\n\n";
    
    file << "[BC1]\n";
    file << "# Boundary condition 1\n";
    file << "type = DIRICHLET                      # DIRICHLET, NEUMANN, ROBIN\n";
    file << "field = PRESSURE                      # PRESSURE, TEMPERATURE, DISPLACEMENT\n";
    file << "location = XMIN                       # XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX\n";
    file << "value = 20.0e6                        # Pa (20 MPa)\n";
    file << "gradient = 0.0                        # For Neumann BC\n\n";
    
    file << "[BC2]\n";
    file << "type = NEUMANN\n";
    file << "field = DISPLACEMENT\n";
    file << "location = ZMAX\n";
    file << "value = 0.0\n";
    file << "gradient = 0.0                        # No stress at top\n\n";
    
    file << "[IC1]\n";
    file << "# Initial condition for pressure\n";
    file << "field = PRESSURE\n";
    file << "distribution = GRADIENT               # UNIFORM, GRADIENT, FILE\n";
    file << "value = 20.0e6                        # Pa (base value)\n";
    file << "gradient = 0.0, 0.0, 10000.0          # Pa/m (hydrostatic gradient in z)\n";
    file << "file =                                # Optional file for complex distributions\n\n";
    
    file << "[IC2]\n";
    file << "# Initial temperature\n";
    file << "field = TEMPERATURE\n";
    file << "distribution = GRADIENT\n";
    file << "value = 300.0                         # K (at top)\n";
    file << "gradient = 0.0, 0.0, 0.025            # K/m (geothermal gradient)\n";
    file << "file = \n\n";
    
    file << "# Multiple rock types can be defined as [ROCK1], [ROCK2], etc.\n";
    file << "# Multiple wells as [WELL1], [WELL2], [WELL3], etc.\n";
    file << "# Multiple fractures as [FRACTURE1], [FRACTURE2], etc.\n";
    file << "# Multiple faults as [FAULT1], [FAULT2], etc.\n";
    file << "# Multiple boundary conditions as [BC1], [BC2], etc.\n";
    file << "# Multiple initial conditions as [IC1], [IC2], etc.\n";
    
    file.close();
}

} // namespace FSRM
