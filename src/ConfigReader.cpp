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

// =============================================================================
// Unit-Aware Value Accessors
// =============================================================================

double ConfigReader::getDoubleWithUnit(const std::string& section, const std::string& key,
                                       double default_val, const std::string& default_unit) const {
    std::string val = getString(section, key);
    if (val.empty()) {
        return default_val;
    }
    
    // Try to parse as value with unit
    double parsed_value;
    std::string parsed_unit;
    
    if (unit_system_.parseValueWithUnit(val, parsed_value, parsed_unit)) {
        // Value has unit specified
        if (!parsed_unit.empty()) {
            try {
                return unit_system_.toBase(parsed_value, parsed_unit);
            } catch (const std::exception& e) {
                std::cerr << "Warning: Unit conversion error for [" << section 
                         << "]:" << key << " - " << e.what() << std::endl;
                return default_val;
            }
        } else if (!default_unit.empty()) {
            // No unit specified, use default unit
            try {
                return unit_system_.toBase(parsed_value, default_unit);
            } catch (const std::exception& e) {
                std::cerr << "Warning: Unit conversion error for [" << section 
                         << "]:" << key << " - " << e.what() << std::endl;
                return default_val;
            }
        } else {
            // No unit specified and no default, assume already in SI
            return parsed_value;
        }
    }
    
    // Fallback: try to parse as plain double
    try {
        double value = std::stod(val);
        if (!default_unit.empty()) {
            return unit_system_.toBase(value, default_unit);
        }
        return value;
    } catch (...) {
        return default_val;
    }
}

std::vector<double> ConfigReader::getDoubleArrayWithUnit(const std::string& section,
                                                         const std::string& key,
                                                         const std::string& default_unit) const {
    std::vector<double> result;
    std::string val = getString(section, key);
    if (val.empty()) return result;
    
    auto tokens = split(val, ',');
    for (const auto& token : tokens) {
        double parsed_value;
        std::string parsed_unit;
        
        if (unit_system_.parseValueWithUnit(token, parsed_value, parsed_unit)) {
            if (!parsed_unit.empty()) {
                try {
                    result.push_back(unit_system_.toBase(parsed_value, parsed_unit));
                } catch (const std::exception& e) {
                    std::cerr << "Warning: Cannot convert '" << token << "': " 
                             << e.what() << std::endl;
                }
            } else if (!default_unit.empty()) {
                try {
                    result.push_back(unit_system_.toBase(parsed_value, default_unit));
                } catch (const std::exception& e) {
                    std::cerr << "Warning: Cannot convert '" << token << "': " 
                             << e.what() << std::endl;
                }
            } else {
                result.push_back(parsed_value);
            }
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
    config.output_format = getString("SIMULATION", "output_format", "HDF5");
    
    return true;
}

bool ConfigReader::parseGridConfig(GridConfig& config) {
    if (!hasSection("GRID")) return false;
    
    // Grid dimensions (for structured grids)
    config.nx = getInt("GRID", "nx", 10);
    config.ny = getInt("GRID", "ny", 10);
    config.nz = getInt("GRID", "nz", 10);
    
    config.Lx = getDouble("GRID", "Lx", 1000.0);
    config.Ly = getDouble("GRID", "Ly", 1000.0);
    config.Lz = getDouble("GRID", "Lz", 100.0);
    
    // Domain origin
    config.origin_x = getDouble("GRID", "origin_x", 0.0);
    config.origin_y = getDouble("GRID", "origin_y", 0.0);
    config.origin_z = getDouble("GRID", "origin_z", 0.0);
    
    // Mesh type
    std::string mesh_type_str = getString("GRID", "mesh_type", "CARTESIAN");
    if (mesh_type_str == "CARTESIAN") config.mesh_type = MeshType::CARTESIAN;
    else if (mesh_type_str == "CORNER_POINT") config.mesh_type = MeshType::CORNER_POINT;
    else if (mesh_type_str == "GMSH") config.mesh_type = MeshType::GMSH;
    else if (mesh_type_str == "EXODUS") config.mesh_type = MeshType::EXODUS;
    else if (mesh_type_str == "CUSTOM") config.mesh_type = MeshType::CUSTOM;
    
    config.use_unstructured = getBool("GRID", "use_unstructured", false);
    config.mesh_file = getString("GRID", "mesh_file", "");
    
    // Auto-detect mesh type from file extension
    if (!config.mesh_file.empty() && config.mesh_type == MeshType::CARTESIAN) {
        std::string ext = config.mesh_file.substr(config.mesh_file.find_last_of('.') + 1);
        if (ext == "msh") {
            config.mesh_type = MeshType::GMSH;
            config.use_unstructured = true;
        } else if (ext == "exo" || ext == "e" || ext == "ex2") {
            config.mesh_type = MeshType::EXODUS;
            config.use_unstructured = true;
        }
    }
    
    // Gmsh-specific options
    config.gmsh_physical_volume = getString("GRID", "gmsh_physical_volume", "");
    
    std::string boundary_str = getString("GRID", "gmsh_boundaries", "");
    if (!boundary_str.empty()) {
        config.gmsh_boundaries = split(boundary_str, ',');
    }
    
    config.gmsh_refinement_level = getInt("GRID", "gmsh_refinement_level", 0);
    
    // Coordinate Reference System (CRS)
    config.input_crs = getString("GRID", "input_crs", "");
    config.model_crs = getString("GRID", "model_crs", "");
    config.use_local_coordinates = getBool("GRID", "use_local_coordinates", true);
    config.local_origin_x = getDouble("GRID", "local_origin_x", 0.0);
    config.local_origin_y = getDouble("GRID", "local_origin_y", 0.0);
    config.local_origin_z = getDouble("GRID", "local_origin_z", 0.0);
    config.auto_detect_utm = getBool("GRID", "auto_detect_utm", false);
    
    // Grid quality settings
    config.min_cell_volume = getDouble("GRID", "min_cell_volume", 1e-10);
    config.max_aspect_ratio = getDouble("GRID", "max_aspect_ratio", 100.0);
    
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
        
        // Rock properties (with automatic unit conversion)
        mat.porosity = getDouble(section, "porosity", 0.2);  // Dimensionless
        mat.permeability_x = getDoubleWithUnit(section, "permeability_x", 100.0e-15, "mD");
        mat.permeability_y = getDoubleWithUnit(section, "permeability_y", 100.0e-15, "mD");
        mat.permeability_z = getDoubleWithUnit(section, "permeability_z", 10.0e-15, "mD");
        mat.compressibility = getDoubleWithUnit(section, "compressibility", 1e-9, "1/Pa");
        
        // Mechanical properties
        mat.youngs_modulus = getDoubleWithUnit(section, "youngs_modulus", 10e9, "Pa");
        mat.poisson_ratio = getDouble(section, "poisson_ratio", 0.25);  // Dimensionless
        mat.density = getDoubleWithUnit(section, "density", 2500.0, "kg/m3");
        mat.biot_coefficient = getDouble(section, "biot_coefficient", 1.0);  // Dimensionless
        
        // Viscoelastic properties
        mat.relaxation_time = getDoubleWithUnit(section, "relaxation_time", 1e6, "s");
        mat.viscosity = getDoubleWithUnit(section, "viscosity", 1e19, "Pa-s");
        
        // Thermal properties
        mat.thermal_conductivity = getDoubleWithUnit(section, "thermal_conductivity", 2.5, "W/(m-K)");
        mat.heat_capacity = getDoubleWithUnit(section, "heat_capacity", 900.0, "J/(kg-K)");
        mat.thermal_expansion = getDoubleWithUnit(section, "thermal_expansion", 1e-5, "1/K");
        
        // Fracture properties
        mat.fracture_toughness = getDoubleWithUnit(section, "fracture_toughness", 1e6, "Pa-m0.5");
        mat.fracture_energy = getDoubleWithUnit(section, "fracture_energy", 100.0, "J/m2");
        
        props.push_back(mat);
    }
    
    return true;
}

bool ConfigReader::parseFluidProperties(FluidProperties& props) {
    if (!hasSection("FLUID")) return false;
    
    // Single phase properties (with automatic unit conversion)
    props.density = getDoubleWithUnit("FLUID", "density", 1000.0, "kg/m3");
    props.viscosity = getDoubleWithUnit("FLUID", "viscosity", 0.001, "Pa-s");
    props.compressibility = getDoubleWithUnit("FLUID", "compressibility", 1e-9, "1/Pa");
    
    // Black oil properties
    props.oil_density_std = getDoubleWithUnit("FLUID", "oil_density_std", 800.0, "kg/m3");
    props.gas_density_std = getDoubleWithUnit("FLUID", "gas_density_std", 1.0, "kg/m3");
    props.water_density_std = getDoubleWithUnit("FLUID", "water_density_std", 1000.0, "kg/m3");
    
    props.oil_viscosity = getDoubleWithUnit("FLUID", "oil_viscosity", 0.005, "Pa-s");
    props.gas_viscosity = getDoubleWithUnit("FLUID", "gas_viscosity", 0.00001, "Pa-s");
    props.water_viscosity = getDoubleWithUnit("FLUID", "water_viscosity", 0.001, "Pa-s");
    
    props.solution_GOR = getDouble("FLUID", "solution_GOR", 100.0); // Dimensionless ratio
    
    // Compositional properties
    props.component_mw = getDoubleArrayWithUnit("FLUID", "component_mw", "g");
    props.component_Tc = getDoubleArrayWithUnit("FLUID", "component_Tc", "K");
    props.component_Pc = getDoubleArrayWithUnit("FLUID", "component_Pc", "Pa");
    props.component_omega = getDoubleArray("FLUID", "component_omega");  // Dimensionless
    
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
            
            // Split node options
            fault.use_split_nodes = getBool(section, "use_split_nodes", false);
            fault.split_node_method = getString(section, "split_node_method", "PENALTY");
            fault.traction_type = getString(section, "traction_type", "FRICTION_DEPENDENT");
            fault.penalty_normal = getDouble(section, "penalty_normal", 1e12);
            fault.penalty_tangent = getDouble(section, "penalty_tangent", 1e10);
            fault.prescribed_traction_normal = getDouble(section, "prescribed_traction_normal", 0.0);
            fault.prescribed_traction_strike = getDouble(section, "prescribed_traction_strike", 0.0);
            fault.prescribed_traction_dip = getDouble(section, "prescribed_traction_dip", 0.0);
            fault.cohesive_strength = getDouble(section, "cohesive_strength", 5e6);
            fault.critical_opening = getDouble(section, "critical_opening", 1e-4);
            fault.critical_slip = getDouble(section, "critical_slip", 1e-3);
            fault.allow_separation = getBool(section, "allow_separation", true);
            
            faults.push_back(fault);
        }
    }
    
    return faults;
}

std::vector<std::unique_ptr<SplitNodeFault>> ConfigReader::parseSplitNodeFaults() {
    std::vector<std::unique_ptr<SplitNodeFault>> split_faults;
    
    for (const auto& section : getSections()) {
        if (section.find("FAULT") == 0 || section.find("SPLIT_NODE_FAULT") == 0) {
            // Check if split nodes are enabled
            if (!getBool(section, "use_split_nodes", false)) {
                continue;
            }
            
            // Get section data and create split node fault
            auto section_data = getSectionData(section);
            auto fault = createSplitNodeFault(section_data);
            
            if (fault) {
                split_faults.push_back(std::move(fault));
            }
        }
    }
    
    return split_faults;
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
    file << "output_format = HDF5                  # HDF5 (default), VTK (secondary), or ECLIPSE\n";
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

// =============================================================================
// Extended Parsing Methods
// =============================================================================

std::unique_ptr<FluidModelBase> ConfigReader::parseFluidModel() {
    std::string fluid_type = getString("FLUID", "type", "SINGLE_PHASE");
    auto section_data = getSectionData("FLUID");
    return createFluidModel(fluid_type, section_data);
}

std::vector<RockProperties> ConfigReader::parseRockProperties() {
    std::vector<RockProperties> rocks;
    
    auto rock_sections = getSectionsMatching("ROCK");
    for (const auto& section : rock_sections) {
        RockProperties rock;
        auto section_data = getSectionData(section);
        rock.configure(section_data);
        rock.name = getString(section, "name", section);
        rocks.push_back(std::move(rock));
    }
    
    return rocks;
}

std::unique_ptr<MaterialModelBase> ConfigReader::parseMaterialModel(const std::string& section) {
    std::string model_type = getString(section, "constitutive_model", "LINEAR_ELASTIC");
    auto section_data = getSectionData(section);
    return createMaterialModel(model_type, section_data);
}

std::unique_ptr<FaultNetwork> ConfigReader::parseFaultNetwork() {
    auto network = std::make_unique<FaultNetwork>();
    
    auto fault_sections = getSectionsMatching("FAULT");
    for (const auto& section : fault_sections) {
        auto section_data = getSectionData(section);
        network->addFault(section_data);
    }
    
    return network;
}

std::vector<ConfigReader::FaultConfig> ConfigReader::parseFaultsExtended() {
    // Keep backward compatibility with existing method
    return parseFaults();
}

bool ConfigReader::parseDynamicsConfig(DynamicsConfig& config) {
    if (!hasSection("DYNAMICS")) return false;
    
    config.enable_dynamics = getBool("DYNAMICS", "enable", false);
    config.use_static_triggering = getBool("DYNAMICS", "static_triggering", false);
    config.trigger_threshold = getDouble("DYNAMICS", "trigger_threshold", 1e6);
    config.event_duration = getDouble("DYNAMICS", "event_duration", 10.0);
    config.damping_alpha = getDouble("DYNAMICS", "damping_alpha", 0.01);
    config.damping_beta = getDouble("DYNAMICS", "damping_beta", 0.001);
    config.quality_factor = getDouble("DYNAMICS", "quality_factor", 100.0);
    config.enable_dynamic_permeability = getBool("DYNAMICS", "dynamic_permeability", false);
    config.permeability_sensitivity = getDouble("DYNAMICS", "permeability_sensitivity", 1.0);
    config.permeability_recovery_time = getDouble("DYNAMICS", "permeability_recovery_time", 100.0);
    
    return true;
}

bool ConfigReader::parseSeismicityConfig(SeismicityConfig& config) {
    if (!hasSection("SEISMICITY")) return false;
    
    config.enable_seismicity = getBool("SEISMICITY", "enable", false);
    config.friction_law = getString("SEISMICITY", "friction_law", "COULOMB");
    config.nucleation_size = getDouble("SEISMICITY", "nucleation_size", 1.0);
    config.seismic_slip_rate = getDouble("SEISMICITY", "seismic_slip_rate", 1e-3);
    config.b_value = getDouble("SEISMICITY", "b_value", 1.0);
    config.enable_aftershocks = getBool("SEISMICITY", "aftershocks", false);
    config.enable_stress_transfer = getBool("SEISMICITY", "stress_transfer", false);
    
    return true;
}

bool ConfigReader::parseOutputConfig(OutputConfig& config) {
    if (!hasSection("OUTPUT")) {
        // Use defaults
        config.format = "HDF5";
        config.base_path = "output";
        config.frequency = 10;
        config.write_pressure = true;
        config.write_displacement = true;
        config.write_stress = false;
        config.write_strain = false;
        config.write_velocity = false;
        config.write_permeability = false;
        config.write_temperature = false;
        config.write_saturation = false;
        config.write_fault_slip = false;
        config.write_seismic_catalog = false;
        return true;
    }
    
    config.format = getString("OUTPUT", "format", "VTK");
    config.base_path = getString("OUTPUT", "path", "output");
    config.frequency = getInt("OUTPUT", "frequency", 10);
    config.write_pressure = getBool("OUTPUT", "pressure", true);
    config.write_displacement = getBool("OUTPUT", "displacement", true);
    config.write_stress = getBool("OUTPUT", "stress", false);
    config.write_strain = getBool("OUTPUT", "strain", false);
    config.write_velocity = getBool("OUTPUT", "velocity", false);
    config.write_permeability = getBool("OUTPUT", "permeability", false);
    config.write_temperature = getBool("OUTPUT", "temperature", false);
    config.write_saturation = getBool("OUTPUT", "saturation", false);
    config.write_fault_slip = getBool("OUTPUT", "fault_slip", false);
    config.write_seismic_catalog = getBool("OUTPUT", "seismic_catalog", false);
    
    // Parse output unit preferences
    config.output_units["pressure"] = getString("OUTPUT", "pressure_unit", "MPa");
    config.output_units["displacement"] = getString("OUTPUT", "displacement_unit", "m");
    config.output_units["stress"] = getString("OUTPUT", "stress_unit", "MPa");
    config.output_units["permeability"] = getString("OUTPUT", "permeability_unit", "mD");
    config.output_units["temperature"] = getString("OUTPUT", "temperature_unit", "degC");
    config.output_units["density"] = getString("OUTPUT", "density_unit", "kg/m3");
    config.output_units["viscosity"] = getString("OUTPUT", "viscosity_unit", "cP");
    config.output_units["length"] = getString("OUTPUT", "length_unit", "m");
    config.output_units["time"] = getString("OUTPUT", "time_unit", "s");
    
    return true;
}

std::map<std::string, std::string> ConfigReader::getSectionData(const std::string& section) const {
    auto it = data.find(section);
    if (it != data.end()) {
        return it->second;
    }
    return {};
}

std::vector<std::string> ConfigReader::getSectionsMatching(const std::string& prefix) const {
    std::vector<std::string> result;
    for (const auto& pair : data) {
        if (pair.first.find(prefix) == 0) {
            result.push_back(pair.first);
        }
    }
    return result;
}

bool ConfigReader::mergeFile(const std::string& filename) {
    ConfigReader other;
    if (!other.loadFile(filename)) {
        return false;
    }
    
    // Merge data - other file values override existing
    for (const auto& section : other.data) {
        for (const auto& key_val : section.second) {
            data[section.first][key_val.first] = key_val.second;
        }
    }
    
    return true;
}

ConfigReader::ValidationResult ConfigReader::validate() const {
    ValidationResult result;
    result.valid = true;
    
    // Check for required sections
    if (!hasSection("SIMULATION")) {
        result.warnings.push_back("No [SIMULATION] section found - using defaults");
    }
    
    if (!hasSection("GRID")) {
        result.warnings.push_back("No [GRID] section found - using defaults");
    }
    
    // Check for at least one rock definition
    if (getSectionsMatching("ROCK").empty()) {
        result.warnings.push_back("No rock properties defined");
    }
    
    // Check for wells if fluid flow enabled
    auto wells = getSectionsMatching("WELL");
    if (wells.empty()) {
        result.warnings.push_back("No wells defined");
    }
    
    // Validate physics consistency
    bool enable_geo = hasKey("SIMULATION", "enable_geomechanics") && 
                     getBool("SIMULATION", "enable_geomechanics", false);
    bool enable_faults = hasKey("SIMULATION", "enable_faults") &&
                        getBool("SIMULATION", "enable_faults", false);
    
    if (enable_faults && !enable_geo) {
        result.warnings.push_back("Faults enabled but geomechanics disabled - enabling geomechanics");
    }
    
    // Check for valid parameter ranges
    if (hasSection("ROCK")) {
        double porosity = getDouble("ROCK", "porosity", 0.2);
        if (porosity < 0.0 || porosity > 1.0) {
            result.errors.push_back("Invalid porosity value (must be 0-1)");
            result.valid = false;
        }
        
        double nu = getDouble("ROCK", "poisson_ratio", 0.25);
        if (nu < -1.0 || nu >= 0.5) {
            result.errors.push_back("Invalid Poisson's ratio (must be -1 to 0.5)");
            result.valid = false;
        }
    }
    
    return result;
}

std::vector<double> ConfigReader::parseDoubleArray(const std::string& value) const {
    return getDoubleArray("", "");  // Already implemented, reuse via split
}

std::vector<int> ConfigReader::parseIntArray(const std::string& value) const {
    std::vector<int> result;
    auto tokens = split(value, ',');
    for (const auto& token : tokens) {
        try {
            result.push_back(std::stoi(token));
        } catch (...) {}
    }
    return result;
}

std::vector<std::string> ConfigReader::parseStringArray(const std::string& value) const {
    return split(value, ',');
}

bool ConfigReader::parseAMRConfig(AMRConfig& config) {
    // Main AMR section
    if (!hasSection("amr")) {
        return false;  // AMR not configured
    }
    
    config.enabled = getBool("amr", "enabled", false);
    config.method = getString("amr", "method", "plex_adapt");
    config.adapt_every = getInt("amr", "adapt_every", 5);
    config.adapt_on_change = getBool("amr", "adapt_on_change", true);
    
    // Criterion subsection
    config.criterion_type = getString("amr.criterion", "type", "gradient");
    config.weight_pressure = getDouble("amr.criterion", "weight_pressure", 1.0);
    config.weight_saturation = getDouble("amr.criterion", "weight_saturation", 1.0);
    config.weight_velocity = getDouble("amr.criterion", "weight_velocity", 0.5);
    config.weight_temperature = getDouble("amr.criterion", "weight_temperature", 0.5);
    
    // Strategy subsection
    config.strategy = getString("amr.strategy", "type", "fixed_fraction");
    config.refine_fraction = getDouble("amr.strategy", "refine_fraction", 0.2);
    config.coarsen_fraction = getDouble("amr.strategy", "coarsen_fraction", 0.05);
    config.refine_threshold = getDouble("amr.strategy", "refine_threshold", 0.5);
    config.coarsen_threshold = getDouble("amr.strategy", "coarsen_threshold", 0.05);
    
    // Limits subsection
    config.max_level = getInt("amr.limits", "max_refinement_level", 5);
    config.min_level = getInt("amr.limits", "min_refinement_level", 0);
    config.max_cells = getInt("amr.limits", "max_cells", 1000000);
    config.min_cells = getInt("amr.limits", "min_cells", 100);
    config.min_cell_size = getDouble("amr.limits", "min_cell_size", 1e-6);
    config.max_cell_size = getDouble("amr.limits", "max_cell_size", 1e6);
    
    // Quality subsection
    config.min_quality = getDouble("amr.quality", "min_quality", 0.1);
    config.max_aspect_ratio = getDouble("amr.quality", "max_aspect_ratio", 10.0);
    
    // Features subsection
    config.preserve_boundaries = getBool("amr.features", "preserve_boundaries", true);
    config.preserve_wells = getBool("amr.features", "preserve_wells", true);
    config.preserve_faults = getBool("amr.features", "preserve_faults", true);
    config.buffer_layers = getInt("amr.features", "buffer_layers", 1);
    
    // Parse well features (amr.features.wells.WELLNAME sections)
    auto well_sections = getSectionsMatching("amr.features.wells.");
    for (const auto& sec : well_sections) {
        AMRFeatureWell well;
        well.name = getString(sec, "name", "");
        well.x = getDouble(sec, "x", 0.0);
        well.y = getDouble(sec, "y", 0.0);
        well.z = getDouble(sec, "z", 0.0);
        well.radius = getDouble(sec, "radius", 10.0);
        well.level = getInt(sec, "level", 2);
        
        if (!well.name.empty() || well.radius > 0) {
            config.wells.push_back(well);
        }
    }
    
    // Parse fault features
    auto fault_sections = getSectionsMatching("amr.features.faults.");
    for (const auto& sec : fault_sections) {
        AMRFeatureFault fault;
        fault.name = getString(sec, "name", "");
        fault.width = getDouble(sec, "width", 10.0);
        fault.level = getInt(sec, "level", 2);
        
        // Parse trace points (format: x1,y1,z1;x2,y2,z2;...)
        std::string trace_str = getString(sec, "trace", "");
        if (!trace_str.empty()) {
            auto points = split(trace_str, ';');
            for (const auto& pt : points) {
                auto coords = parseDoubleArray(pt);
                if (coords.size() >= 3) {
                    fault.trace.push_back({coords[0], coords[1], coords[2]});
                }
            }
        }
        
        if (!fault.trace.empty()) {
            config.faults.push_back(fault);
        }
    }
    
    // Parse box features
    auto box_sections = getSectionsMatching("amr.features.boxes.");
    for (const auto& sec : box_sections) {
        AMRFeatureBox box;
        box.name = getString(sec, "name", "");
        box.xmin = getDouble(sec, "xmin", 0.0);
        box.xmax = getDouble(sec, "xmax", 0.0);
        box.ymin = getDouble(sec, "ymin", 0.0);
        box.ymax = getDouble(sec, "ymax", 0.0);
        box.zmin = getDouble(sec, "zmin", 0.0);
        box.zmax = getDouble(sec, "zmax", 0.0);
        box.level = getInt(sec, "level", 2);
        
        if (box.xmax > box.xmin || box.ymax > box.ymin || box.zmax > box.zmin) {
            config.boxes.push_back(box);
        }
    }
    
    // Transfer subsection
    config.transfer_method = getString("amr.transfer", "method", "interpolation");
    config.conservative = getBool("amr.transfer", "conservative", true);
    
    // Balance subsection
    config.balance_enabled = getBool("amr.balance", "enabled", true);
    config.balance_method = getString("amr.balance", "method", "parmetis");
    config.rebalance_threshold = getDouble("amr.balance", "rebalance_threshold", 0.2);
    
    // Output subsection
    config.write_mesh = getBool("amr.output", "write_mesh", true);
    config.write_errors = getBool("amr.output", "write_errors", true);
    config.output_format = getString("amr.output", "format", "vtk");
    config.output_prefix = getString("amr.output", "prefix", "amr_output");
    
    return true;
}

void ConfigReader::generateCompleteTemplate(const std::string& filename) {
    std::ofstream file(filename);
    
    file << "# ============================================================================\n";
    file << "# FSRM Complete Configuration Template\n";
    file << "# ============================================================================\n";
    file << "# This file contains ALL available configuration options.\n";
    file << "# Edit values as needed for your simulation.\n";
    file << "# All units in SI unless otherwise specified.\n";
    file << "#\n";
    file << "# Sections:\n";
    file << "#   [SIMULATION] - Time stepping, physics, solver settings\n";
    file << "#   [GRID] - Domain and mesh configuration\n";
    file << "#   [ROCK] - Material properties\n";
    file << "#   [FLUID] - Fluid properties\n";
    file << "#   [DYNAMICS] - Wave propagation and dynamic effects\n";
    file << "#   [SEISMICITY] - Fault mechanics and earthquakes\n";
    file << "#   [OUTPUT] - Output configuration\n";
    file << "#   [WELL*] - Well definitions\n";
    file << "#   [FRACTURE*] - Fracture definitions\n";
    file << "#   [FAULT*] - Fault definitions\n";
    file << "#   [BC*] - Boundary conditions\n";
    file << "#   [IC*] - Initial conditions\n";
    file << "# ============================================================================\n\n";
    
    file << "[SIMULATION]\n";
    file << "# Time parameters\n";
    file << "start_time = 0.0\n";
    file << "end_time = 86400.0                    # seconds (1 day)\n";
    file << "dt_initial = 3600.0                   # seconds (1 hour)\n";
    file << "dt_min = 0.001                        # seconds\n";
    file << "dt_max = 86400.0                      # seconds\n";
    file << "max_timesteps = 10000\n";
    file << "adaptive_timestepping = true\n\n";
    
    file << "# Solver settings\n";
    file << "rtol = 1.0e-6                         # Relative tolerance\n";
    file << "atol = 1.0e-8                         # Absolute tolerance\n";
    file << "max_nonlinear_iterations = 50\n";
    file << "max_linear_iterations = 1000\n\n";
    
    file << "# Physics models\n";
    file << "fluid_model = SINGLE_COMPONENT        # SINGLE_COMPONENT, BLACK_OIL, COMPOSITIONAL\n";
    file << "solid_model = ELASTIC                 # ELASTIC, VISCOELASTIC, POROELASTIC\n\n";
    
    file << "# Physics modules (enable/disable)\n";
    file << "enable_geomechanics = false\n";
    file << "enable_thermal = false\n";
    file << "enable_fractures = false\n";
    file << "enable_particle_transport = false\n";
    file << "enable_faults = false\n";
    file << "enable_tidal_forces = false\n";
    file << "enable_elastodynamics = false\n";
    file << "enable_poroelastodynamics = false\n\n";
    
    file << "# GPU settings\n";
    file << "use_gpu = false\n";
    file << "gpu_mode = CPU_FALLBACK               # CPU_ONLY, GPU_ONLY, HYBRID, CPU_FALLBACK\n";
    file << "gpu_device_id = 0\n";
    file << "gpu_memory_fraction = 0.8\n\n";
    
    file << "[GRID]\n";
    file << "nx = 20\n";
    file << "ny = 20\n";
    file << "nz = 5\n";
    file << "Lx = 1000.0                           # meters\n";
    file << "Ly = 1000.0\n";
    file << "Lz = 100.0\n";
    file << "use_unstructured = false\n";
    file << "mesh_file = \n\n";
    
    file << "[ROCK]\n";
    file << "# Constitutive model\n";
    file << "constitutive_model = LINEAR_ELASTIC   # LINEAR_ELASTIC, VISCOELASTIC, POROELASTIC, etc.\n\n";
    
    file << "# Flow properties\n";
    file << "porosity = 0.20\n";
    file << "permeability_x = 100.0                # milliDarcy\n";
    file << "permeability_y = 100.0\n";
    file << "permeability_z = 10.0\n";
    file << "compressibility = 1.0e-9              # 1/Pa\n\n";
    
    file << "# Elastic properties\n";
    file << "youngs_modulus = 10.0e9               # Pa\n";
    file << "poisson_ratio = 0.25\n";
    file << "density = 2500.0                      # kg/m³\n";
    file << "biot_coefficient = 1.0\n\n";
    
    file << "# Viscoelastic properties (if model = VISCOELASTIC)\n";
    file << "viscoelastic_type = MAXWELL           # MAXWELL, KELVIN_VOIGT, SLS\n";
    file << "relaxation_time = 1.0e6               # seconds\n";
    file << "viscosity = 1.0e19                    # Pa·s\n";
    file << "long_term_modulus = 5.0e9             # Pa (for SLS)\n\n";
    
    file << "# Poroelastic properties (if model = POROELASTIC)\n";
    file << "biot_modulus = 1.0e10                 # Pa\n";
    file << "grain_bulk_modulus = 35.0e9           # Pa\n";
    file << "fluid_bulk_modulus = 2.2e9            # Pa\n\n";
    
    file << "# Plasticity (if model = ELASTOPLASTIC)\n";
    file << "failure_criterion = MOHR_COULOMB      # MOHR_COULOMB, DRUCKER_PRAGER, VON_MISES\n";
    file << "cohesion = 1.0e6                      # Pa\n";
    file << "friction_angle = 30.0                 # degrees\n";
    file << "dilation_angle = 0.0                  # degrees\n";
    file << "tensile_strength = 0.5e6              # Pa\n\n";
    
    file << "# Wave properties\n";
    file << "p_wave_velocity = 5000.0              # m/s\n";
    file << "s_wave_velocity = 3000.0              # m/s\n";
    file << "quality_factor = 100.0                # Q (attenuation)\n";
    file << "damping_alpha = 0.01                  # Rayleigh mass damping\n";
    file << "damping_beta = 0.001                  # Rayleigh stiffness damping\n\n";
    
    file << "# Thermal properties\n";
    file << "thermal_conductivity = 2.5            # W/(m·K)\n";
    file << "heat_capacity = 900.0                 # J/(kg·K)\n";
    file << "thermal_expansion = 1.0e-5            # 1/K\n\n";
    
    file << "# Fracture properties\n";
    file << "fracture_toughness = 1.0e6            # Pa·m^0.5\n";
    file << "fracture_energy = 100.0               # J/m²\n\n";
    
    file << "# Permeability model\n";
    file << "permeability_model = CONSTANT         # CONSTANT, KOZENY_CARMAN, EXPONENTIAL\n";
    file << "permeability_stress_coeff = 1.0e-8    # For stress-dependent models\n\n";
    
    file << "[FLUID]\n";
    file << "type = SINGLE_PHASE                   # SINGLE_PHASE, BLACK_OIL, COMPOSITIONAL, BRINE, CO2\n\n";
    
    file << "# Single-phase properties\n";
    file << "density = 1000.0                      # kg/m³\n";
    file << "viscosity = 0.001                     # Pa·s\n";
    file << "compressibility = 1.0e-9             # 1/Pa\n";
    file << "viscosity_model = constant            # constant, exponential, arrhenius\n\n";
    
    file << "# Black oil properties\n";
    file << "oil_density_std = 850.0               # kg/m³\n";
    file << "oil_viscosity_dead = 0.005            # Pa·s\n";
    file << "oil_api = 35.0                        # API gravity\n";
    file << "gas_density_std = 1.0                 # kg/m³\n";
    file << "gas_gravity = 0.7                     # Specific gravity (air=1)\n";
    file << "water_density_std = 1020.0            # kg/m³\n";
    file << "water_viscosity = 0.0005              # Pa·s\n";
    file << "solution_gor = 150.0                  # sm³/sm³\n";
    file << "bubble_point_pressure = 15.0e6        # Pa\n";
    file << "pvt_correlation = STANDING            # STANDING, VASQUEZ_BEGGS\n\n";
    
    file << "# Relative permeability (Corey)\n";
    file << "swc = 0.2                             # Connate water saturation\n";
    file << "sor = 0.2                             # Residual oil saturation\n";
    file << "sgc = 0.05                            # Critical gas saturation\n";
    file << "corey_nw = 3.0                        # Water Corey exponent\n";
    file << "corey_no = 2.0                        # Oil Corey exponent\n";
    file << "corey_ng = 2.0                        # Gas Corey exponent\n\n";
    
    file << "# Compositional properties (comma-separated)\n";
    file << "eos = PENG_ROBINSON                   # PENG_ROBINSON, SRK\n";
    file << "component_mw = 16.04, 30.07, 44.10    # g/mol\n";
    file << "component_tc = 190.6, 305.4, 369.8    # K\n";
    file << "component_pc = 4.6e6, 4.88e6, 4.25e6  # Pa\n";
    file << "component_omega = 0.011, 0.099, 0.152\n";
    file << "composition = 0.6, 0.3, 0.1           # Initial mole fractions\n\n";
    
    file << "[DYNAMICS]\n";
    file << "enable = false                        # Enable dynamic wave effects\n";
    file << "static_triggering = false             # Static stress triggers dynamic event\n";
    file << "trigger_threshold = 1.0e6             # Pa (stress threshold)\n";
    file << "event_duration = 10.0                 # seconds\n";
    file << "damping_alpha = 0.01                  # Rayleigh mass damping\n";
    file << "damping_beta = 0.001                  # Rayleigh stiffness damping\n";
    file << "quality_factor = 100.0                # Attenuation Q\n";
    file << "dynamic_permeability = false          # Enable k changes from waves\n";
    file << "permeability_sensitivity = 1.0        # Sensitivity factor\n";
    file << "permeability_recovery_time = 100.0    # seconds\n\n";
    
    file << "[SEISMICITY]\n";
    file << "enable = false                        # Enable seismicity modeling\n";
    file << "friction_law = COULOMB                # COULOMB, RATE_STATE_AGING, RATE_STATE_SLIP\n";
    file << "nucleation_size = 1.0                 # meters\n";
    file << "seismic_slip_rate = 1.0e-3            # m/s (threshold for seismic slip)\n";
    file << "b_value = 1.0                         # Gutenberg-Richter b-value\n";
    file << "aftershocks = false                   # Model aftershock sequences\n";
    file << "stress_transfer = false               # Include Coulomb stress transfer\n\n";
    
    file << "[OUTPUT]\n";
    file << "format = VTK                          # VTK, HDF5, ECLIPSE\n";
    file << "path = output                         # Output directory\n";
    file << "frequency = 10                        # Write every N steps\n";
    file << "pressure = true                       # Write pressure field\n";
    file << "displacement = true                   # Write displacement field\n";
    file << "stress = false                        # Write stress tensor\n";
    file << "strain = false                        # Write strain tensor\n";
    file << "velocity = false                      # Write velocity field\n";
    file << "permeability = false                  # Write permeability field\n";
    file << "temperature = false                   # Write temperature field\n";
    file << "saturation = false                    # Write saturation fields\n";
    file << "fault_slip = false                    # Write fault slip history\n";
    file << "seismic_catalog = false               # Write earthquake catalog\n\n";
    
    file << "# ============================================================================\n";
    file << "# ADAPTIVE MESH REFINEMENT (AMR)\n";
    file << "# ============================================================================\n\n";
    
    file << "[amr]\n";
    file << "enabled = false                       # Enable AMR\n";
    file << "method = plex_adapt                   # plex_refine, plex_adapt, forest_p4est\n";
    file << "adapt_every = 5                       # Time steps between adaptations\n";
    file << "adapt_on_change = true                # Adapt when solution changes significantly\n\n";
    
    file << "[amr.criterion]\n";
    file << "type = gradient                       # gradient, hessian, residual, jump, feature, combined\n";
    file << "weight_pressure = 1.0                 # Weight for pressure error\n";
    file << "weight_saturation = 1.0               # Weight for saturation error\n";
    file << "weight_velocity = 0.5                 # Weight for velocity error\n";
    file << "weight_temperature = 0.5              # Weight for temperature error\n\n";
    
    file << "[amr.strategy]\n";
    file << "type = fixed_fraction                 # fixed_fraction, threshold, equilibration\n";
    file << "refine_fraction = 0.2                 # Fraction of cells to refine\n";
    file << "coarsen_fraction = 0.05               # Fraction of cells to coarsen\n";
    file << "refine_threshold = 0.5                # Threshold for refinement (threshold mode)\n";
    file << "coarsen_threshold = 0.05              # Threshold for coarsening (threshold mode)\n\n";
    
    file << "[amr.limits]\n";
    file << "max_refinement_level = 4              # Maximum refinement depth\n";
    file << "min_refinement_level = 0              # Minimum refinement level\n";
    file << "max_cells = 100000                    # Maximum number of cells\n";
    file << "min_cells = 100                       # Minimum number of cells\n";
    file << "min_cell_size = 1.0                   # Minimum cell dimension [m]\n";
    file << "max_cell_size = 100.0                 # Maximum cell dimension [m]\n\n";
    
    file << "[amr.quality]\n";
    file << "min_quality = 0.2                     # Minimum cell quality (0-1)\n";
    file << "max_aspect_ratio = 5.0                # Maximum aspect ratio\n\n";
    
    file << "[amr.features]\n";
    file << "preserve_boundaries = true            # Preserve boundary resolution\n";
    file << "preserve_wells = true                 # Preserve well resolution\n";
    file << "preserve_faults = true                # Preserve fault resolution\n";
    file << "buffer_layers = 2                     # Buffer layers around refined cells\n\n";
    
    file << "# Well refinement features (add more sections as needed)\n";
    file << "[amr.features.wells.PROD1]\n";
    file << "name = PROD1\n";
    file << "x = 50.0                              # Well X coordinate [m]\n";
    file << "y = 50.0                              # Well Y coordinate [m]\n";
    file << "z = 0.0                               # Well Z coordinate [m]\n";
    file << "radius = 20.0                         # Refinement radius [m]\n";
    file << "level = 3                             # Target refinement level\n\n";
    
    file << "# Box refinement regions\n";
    file << "[amr.features.boxes.region1]\n";
    file << "name = high_res_region\n";
    file << "xmin = 0.0\n";
    file << "xmax = 100.0\n";
    file << "ymin = 0.0\n";
    file << "ymax = 100.0\n";
    file << "zmin = 0.0\n";
    file << "zmax = 50.0\n";
    file << "level = 2\n\n";
    
    file << "[amr.transfer]\n";
    file << "method = interpolation                # interpolation, projection, injection, conservative\n";
    file << "conservative = true                   # Ensure mass conservation\n\n";
    
    file << "[amr.balance]\n";
    file << "enabled = true                        # Enable load balancing\n";
    file << "method = parmetis                     # parmetis, ptscotch, simple\n";
    file << "rebalance_threshold = 0.2             # Rebalance if imbalance > threshold\n\n";
    
    file << "[amr.output]\n";
    file << "write_mesh = true                     # Output adapted mesh\n";
    file << "write_errors = true                   # Output error field\n";
    file << "format = vtk                          # vtk, hdf5, exodusii\n";
    file << "prefix = amr_output                   # Output file prefix\n\n";
    
    file << "[WELL1]\n";
    file << "name = PROD1\n";
    file << "type = PRODUCER                       # PRODUCER, INJECTOR, OBSERVATION\n";
    file << "i = 10                                # Grid index X\n";
    file << "j = 10                                # Grid index Y\n";
    file << "k = 5                                 # Grid index Z\n";
    file << "control_mode = RATE                   # RATE, BHP, THP\n";
    file << "target_value = 0.01                   # m³/s or Pa\n";
    file << "max_rate = 0.1                        # m³/s\n";
    file << "min_bhp = 5.0e6                       # Pa\n";
    file << "diameter = 0.2                        # meters\n";
    file << "skin = 0.0\n\n";
    
    file << "[FAULT1]\n";
    file << "name = MAIN_FAULT\n";
    file << "x = 5000.0                            # Center X (meters)\n";
    file << "y = 0.0                               # Center Y\n";
    file << "z = 4000.0                            # Center Z (depth)\n";
    file << "strike = 0.0                          # Strike angle (degrees from N)\n";
    file << "dip = 60.0                            # Dip angle (degrees)\n";
    file << "length = 2000.0                       # Along-strike length (meters)\n";
    file << "width = 1500.0                        # Down-dip width (meters)\n";
    file << "friction_law = RATE_STATE_AGING      # COULOMB, RATE_STATE_AGING, etc.\n";
    file << "static_friction = 0.6\n";
    file << "dynamic_friction = 0.4\n";
    file << "cohesion = 1.0e6                      # Pa\n";
    file << "rate_state_a = 0.010                  # Direct effect\n";
    file << "rate_state_b = 0.015                  # Evolution effect (b>a = unstable)\n";
    file << "rate_state_dc = 1.0e-4                # Critical slip distance (meters)\n";
    file << "rate_state_v0 = 1.0e-6                # Reference velocity (m/s)\n";
    file << "rate_state_f0 = 0.6                   # Reference friction\n\n";
    
    file << "[FRACTURE1]\n";
    file << "type = HYDRAULIC                      # NATURAL, HYDRAULIC\n";
    file << "location = 500, 500, 50, 0, 1, 0      # x,y,z center + normal vector\n";
    file << "aperture = 0.001                      # meters\n";
    file << "permeability = 1.0e-10                # m²\n";
    file << "toughness = 1.0e6                     # Pa·m^0.5\n";
    file << "energy = 100.0                        # J/m²\n";
    file << "enable_propagation = true\n";
    file << "enable_proppant = true\n\n";
    
    file << "[BC1]\n";
    file << "type = DIRICHLET                      # DIRICHLET, NEUMANN, ROBIN\n";
    file << "field = PRESSURE                      # PRESSURE, TEMPERATURE, DISPLACEMENT\n";
    file << "location = XMIN                       # XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX\n";
    file << "value = 20.0e6                        # Pa\n\n";
    
    file << "[IC1]\n";
    file << "field = PRESSURE\n";
    file << "distribution = GRADIENT               # UNIFORM, GRADIENT, FILE\n";
    file << "value = 20.0e6                        # Base value (Pa)\n";
    file << "gradient = 0, 0, 10000                # Pa/m (hydrostatic)\n\n";
    
    file << "# ============================================================================\n";
    file << "# END OF TEMPLATE\n";
    file << "# ============================================================================\n";
    
    file.close();
}

} // namespace FSRM
