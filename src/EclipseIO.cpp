#include "EclipseIO.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>

namespace FSRM {

EclipseIO::EclipseIO() : nx(0), ny(0), nz(0), ncells(0) {}

EclipseIO::~EclipseIO() {}

bool EclipseIO::readDeckFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    EclipseKeyword current_keyword;
    bool reading_keyword = false;
    
    while (std::getline(file, line)) {
        // Remove comments
        size_t comment_pos = line.find("--");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) continue;
        
        // Check if this is a keyword
        if (std::isupper(line[0]) && line.find(' ') == std::string::npos && 
            line.back() != '/') {
            // Save previous keyword
            if (reading_keyword) {
                keyword_map[current_keyword.name] = keywords.size();
                keywords.push_back(current_keyword);
            }
            
            // Start new keyword
            current_keyword = EclipseKeyword();
            current_keyword.name = line;
            reading_keyword = true;
        } else if (reading_keyword) {
            // Add data to current keyword
            if (line == "/") {
                // End of keyword
                keyword_map[current_keyword.name] = keywords.size();
                keywords.push_back(current_keyword);
                reading_keyword = false;
            } else {
                current_keyword.data.push_back(line);
            }
        }
    }
    
    file.close();
    return true;
}

EclipseKeyword* EclipseIO::findKeyword(const std::string& name) {
    auto it = keyword_map.find(name);
    if (it != keyword_map.end()) {
        return &keywords[it->second];
    }
    return nullptr;
}

std::vector<std::string> EclipseIO::tokenize(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    
    while (iss >> token) {
        tokens.push_back(token);
    }
    
    return tokens;
}

double EclipseIO::parseValue(const std::string& token) {
    // Handle notation like 100*0.2 (repeat 100 times the value 0.2)
    size_t star_pos = token.find('*');
    if (star_pos != std::string::npos) {
        return std::stod(token.substr(star_pos + 1));
    }
    return std::stod(token);
}

void EclipseIO::expandArrayNotation(std::vector<std::string>& tokens) {
    std::vector<std::string> expanded;
    
    for (const auto& token : tokens) {
        size_t star_pos = token.find('*');
        if (star_pos != std::string::npos) {
            // Format: N*value
            int count = std::stoi(token.substr(0, star_pos));
            std::string value = token.substr(star_pos + 1);
            for (int i = 0; i < count; ++i) {
                expanded.push_back(value);
            }
        } else {
            expanded.push_back(token);
        }
    }
    
    tokens = expanded;
}

bool EclipseIO::parseRunspec() {
    // Parse DIMENS
    parseDimens();
    
    // Parse phase keywords
    parseOil();
    parseWater();
    parseGas();
    
    return true;
}

void EclipseIO::parseDimens() {
    auto kw = findKeyword("DIMENS");
    if (!kw) return;
    
    std::vector<std::string> all_tokens;
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        all_tokens.insert(all_tokens.end(), tokens.begin(), tokens.end());
    }
    
    if (all_tokens.size() >= 3) {
        nx = std::stoi(all_tokens[0]);
        ny = std::stoi(all_tokens[1]);
        nz = std::stoi(all_tokens[2]);
        ncells = nx * ny * nz;
    }
}

void EclipseIO::parseOil() {
    // Check if OIL keyword exists
}

void EclipseIO::parseWater() {
    // Check if WATER keyword exists
}

void EclipseIO::parseGas() {
    // Check if GAS keyword exists
}

void EclipseIO::parseDisgas() {}
void EclipseIO::parseVapoil() {}

bool EclipseIO::parseGrid() {
    parseDX();
    parseDY();
    parseDZ();
    parseTOPS();
    parsePoro();
    parsePermx();
    parsePermy();
    parsePermz();
    parseFaults();
    
    return true;
}

std::vector<double> EclipseIO::parseDataArray(const std::string& keyword_name, int expected_size) {
    std::vector<double> result;
    auto kw = findKeyword(keyword_name);
    if (!kw) {
        // Return default values
        return std::vector<double>(expected_size, 1.0);
    }
    
    std::vector<std::string> all_tokens;
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        all_tokens.insert(all_tokens.end(), tokens.begin(), tokens.end());
    }
    
    expandArrayNotation(all_tokens);
    
    for (const auto& token : all_tokens) {
        if (token != "/") {
            result.push_back(parseValue(token));
        }
    }
    
    return result;
}

std::vector<int> EclipseIO::parseIntArray(const std::string& keyword_name, int expected_size) {
    std::vector<int> result;
    auto kw = findKeyword(keyword_name);
    if (!kw) {
        return std::vector<int>(expected_size, 1);
    }
    
    std::vector<std::string> all_tokens;
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        all_tokens.insert(all_tokens.end(), tokens.begin(), tokens.end());
    }
    
    expandArrayNotation(all_tokens);
    
    for (const auto& token : all_tokens) {
        if (token != "/") {
            result.push_back(std::stoi(token));
        }
    }
    
    return result;
}

void EclipseIO::parseDX() {
    dx = parseDataArray("DX", ncells);
}

void EclipseIO::parseDY() {
    dy = parseDataArray("DY", ncells);
}

void EclipseIO::parseDZ() {
    dz = parseDataArray("DZ", ncells);
}

void EclipseIO::parseTOPS() {
    tops = parseDataArray("TOPS", nx * ny);
}

void EclipseIO::parsePoro() {
    poro = parseDataArray("PORO", ncells);
}

void EclipseIO::parsePermx() {
    permx = parseDataArray("PERMX", ncells);
}

void EclipseIO::parsePermy() {
    permy = parseDataArray("PERMY", ncells);
}

void EclipseIO::parsePermz() {
    permz = parseDataArray("PERMZ", ncells);
}

void EclipseIO::parseFaults() {
    // Parse FAULTS keyword for fault definitions
}

bool EclipseIO::parseEdit() {
    return true;
}

void EclipseIO::parseMultiply() {}
void EclipseIO::parseAdd() {}
void EclipseIO::parseBox() {}

bool EclipseIO::parseProps() {
    parsePVTO();
    parsePVTG();
    parsePVTW();
    parseRock();
    parseDensity();
    parseSWOF();
    parseSGOF();
    
    return true;
}

void EclipseIO::parsePVTO() {
    // Parse oil PVT table
}

void EclipseIO::parsePVTG() {
    // Parse gas PVT table
}

void EclipseIO::parsePVTW() {
    // Parse water PVT
    auto kw = findKeyword("PVTW");
    if (kw && !kw->data.empty()) {
        auto tokens = tokenize(kw->data[0]);
        if (tokens.size() >= 5) {
            pvtw_data.resize(5);
            for (int i = 0; i < 5; ++i) {
                pvtw_data[i] = parseValue(tokens[i]);
            }
        }
    }
}

void EclipseIO::parseRock() {
    // Parse rock compressibility
}

void EclipseIO::parseDensity() {
    // Parse fluid densities at standard conditions
}

void EclipseIO::parseSWOF() {
    // Parse water-oil relative permeability table
}

void EclipseIO::parseSGOF() {
    // Parse gas-oil relative permeability table
}

void EclipseIO::parseSWFN() {}
void EclipseIO::parseSOF3() {}

bool EclipseIO::parseRegions() {
    parseSATNUM();
    parsePVTNUM();
    parseEQLNUM();
    parseFIPNUM();
    
    return true;
}

void EclipseIO::parseSATNUM() {
    satnum = parseIntArray("SATNUM", ncells);
}

void EclipseIO::parsePVTNUM() {
    pvtnum = parseIntArray("PVTNUM", ncells);
}

void EclipseIO::parseEQLNUM() {
    eqlnum = parseIntArray("EQLNUM", ncells);
}

void EclipseIO::parseFIPNUM() {
    fipnum = parseIntArray("FIPNUM", ncells);
}

bool EclipseIO::parseSolution() {
    parseEquil();
    parsePressure();
    parseSWAT();
    parseSGAS();
    parseRS();
    parseRV();
    
    return true;
}

void EclipseIO::parseEquil() {
    // Parse equilibration data
}

void EclipseIO::parsePressure() {
    pressure = parseDataArray("PRESSURE", ncells);
}

void EclipseIO::parseSWAT() {
    swat = parseDataArray("SWAT", ncells);
}

void EclipseIO::parseSGAS() {
    sgas = parseDataArray("SGAS", ncells);
}

void EclipseIO::parseRS() {}
void EclipseIO::parseRV() {}

bool EclipseIO::parseSummary() {
    parseSummaryKeywords();
    return true;
}

void EclipseIO::parseSummaryKeywords() {
    // Parse summary output keywords (WOPR, WGPR, FOPT, etc.)
}

bool EclipseIO::parseSchedule() {
    parseWelspecs();
    parseCompdat();
    parseWconprod();
    parseWconinje();
    parseWeltarg();
    parseTstep();
    parseDates();
    
    return true;
}

void EclipseIO::parseWelspecs() {
    auto kw = findKeyword("WELSPECS");
    if (!kw) return;
    
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        if (tokens.size() >= 5 && tokens[0] != "/") {
            WellSpec well;
            well.name = tokens[0];
            // Remove quotes if present
            well.name.erase(std::remove(well.name.begin(), well.name.end(), '\''), well.name.end());
            well.i = std::stoi(tokens[2]);
            well.j = std::stoi(tokens[3]);
            well.ref_depth = parseValue(tokens[4]);
            
            // Determine well type from group name or default to producer
            well.type = WellType::PRODUCER;
            
            wells.push_back(well);
        }
    }
}

void EclipseIO::parseCompdat() {
    auto kw = findKeyword("COMPDAT");
    if (!kw) return;
    
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        if (tokens.size() >= 7 && tokens[0] != "/") {
            Completion comp;
            comp.well_name = tokens[0];
            comp.well_name.erase(std::remove(comp.well_name.begin(), comp.well_name.end(), '\''), 
                                comp.well_name.end());
            comp.i = std::stoi(tokens[1]);
            comp.j = std::stoi(tokens[2]);
            comp.k1 = std::stoi(tokens[3]);
            comp.k2 = std::stoi(tokens[4]);
            comp.open = (tokens[5] == "OPEN" || tokens[5] == "'OPEN'");
            // Parse more fields as needed
            
            completions.push_back(comp);
        }
    }
}

void EclipseIO::parseWconprod() {
    auto kw = findKeyword("WCONPROD");
    if (!kw) return;
    
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        if (tokens.size() >= 3 && tokens[0] != "/") {
            WellControl ctrl;
            ctrl.well_name = tokens[0];
            ctrl.well_name.erase(std::remove(ctrl.well_name.begin(), ctrl.well_name.end(), '\''),
                                ctrl.well_name.end());
            ctrl.control_mode = tokens[1];
            // Parse control values
            
            well_controls.push_back(ctrl);
        }
    }
}

void EclipseIO::parseWconinje() {
    // Parse injection well controls
}

void EclipseIO::parseWeltarg() {}

void EclipseIO::parseTstep() {
    auto kw = findKeyword("TSTEP");
    if (!kw) return;
    
    for (const auto& line : kw->data) {
        auto tokens = tokenize(line);
        for (const auto& token : tokens) {
            if (token != "/") {
                timesteps.push_back(parseValue(token));
            }
        }
    }
}

void EclipseIO::parseDates() {
    // Parse DATES keyword
}

bool EclipseIO::parseSection(const std::string& section_name) {
    if (section_name == "RUNSPEC") return parseRunspec();
    if (section_name == "GRID") return parseGrid();
    if (section_name == "EDIT") return parseEdit();
    if (section_name == "PROPS") return parseProps();
    if (section_name == "REGIONS") return parseRegions();
    if (section_name == "SOLUTION") return parseSolution();
    if (section_name == "SUMMARY") return parseSummary();
    if (section_name == "SCHEDULE") return parseSchedule();
    
    return false;
}

// Output functions
bool EclipseIO::writeRestartFile(const std::string& filename, int step) {
    // Write binary restart file in Eclipse format
    return true;
}

bool EclipseIO::writeSummaryFile(const std::string& filename) {
    // Write summary file
    return true;
}

bool EclipseIO::writeInitFile(const std::string& filename) {
    // Write INIT file
    return true;
}

// Data accessors
GridConfig EclipseIO::getGridConfig() const {
    GridConfig config;
    config.nx = nx;
    config.ny = ny;
    config.nz = nz;
    
    // Calculate domain size from DX, DY, DZ
    if (!dx.empty() && !dy.empty() && !dz.empty()) {
        config.Lx = 0.0;
        config.Ly = 0.0;
        config.Lz = 0.0;
        
        for (int i = 0; i < nx; ++i) config.Lx += dx[i];
        for (int j = 0; j < ny; ++j) config.Ly += dy[j * nx];
        for (int k = 0; k < nz; ++k) config.Lz += dz[k * nx * ny];
    }
    
    return config;
}

MaterialProperties EclipseIO::getMaterialProperties(int i, int j, int k) const {
    MaterialProperties props;
    
    int idx = i + j * nx + k * nx * ny;
    
    if (idx < ncells) {
        if (!poro.empty()) props.porosity = poro[idx];
        if (!permx.empty()) props.permeability_x = permx[idx];
        if (!permy.empty()) props.permeability_y = permy[idx];
        if (!permz.empty()) props.permeability_z = permz[idx];
    }
    
    return props;
}

FluidProperties EclipseIO::getFluidProperties() const {
    FluidProperties props;
    // Return fluid properties from PVT data
    return props;
}

SimulationConfig EclipseIO::getSimulationConfig() const {
    SimulationConfig config;
    
    // Set timesteps from TSTEP data
    if (!timesteps.empty()) {
        config.dt_initial = timesteps[0];
        config.end_time = 0.0;
        for (double dt : timesteps) {
            config.end_time += dt;
        }
    }
    
    return config;
}

} // namespace FSRM
