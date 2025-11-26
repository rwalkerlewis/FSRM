#ifndef ECLIPSE_IO_HPP
#define ECLIPSE_IO_HPP

#include "ReservoirSim.hpp"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

namespace FSRM {

struct EclipseKeyword {
    std::string name;
    std::vector<std::string> data;
};

class EclipseIO {
public:
    EclipseIO();
    ~EclipseIO();
    
    // Input functions
    bool readDeckFile(const std::string& filename);
    bool parseSection(const std::string& section_name);
    
    // RUNSPEC section
    bool parseRunspec();
    void parseDimens();
    void parseOil();
    void parseWater();
    void parseGas();
    void parseDisgas();
    void parseVapoil();
    
    // GRID section
    bool parseGrid();
    void parseDX();
    void parseDY();
    void parseDZ();
    void parseTOPS();
    void parsePoro();
    void parsePermx();
    void parsePermy();
    void parsePermz();
    void parseFaults();
    
    // EDIT section
    bool parseEdit();
    void parseMultiply();
    void parseAdd();
    void parseBox();
    
    // PROPS section
    bool parseProps();
    void parsePVTO();
    void parsePVTG();
    void parsePVTW();
    void parseRock();
    void parseDensity();
    void parseSWOF();
    void parseSGOF();
    void parseSWFN();
    void parseSOF3();
    
    // REGIONS section
    bool parseRegions();
    void parseSATNUM();
    void parsePVTNUM();
    void parseEQLNUM();
    void parseFIPNUM();
    
    // SOLUTION section
    bool parseSolution();
    void parseEquil();
    void parsePressure();
    void parseSWAT();
    void parseSGAS();
    void parseRS();
    void parseRV();
    
    // SUMMARY section
    bool parseSummary();
    void parseSummaryKeywords();
    
    // SCHEDULE section
    bool parseSchedule();
    void parseWelspecs();
    void parseCompdat();
    void parseWconprod();
    void parseWconinje();
    void parseWeltarg();
    void parseTstep();
    void parseDates();
    
    // Output functions
    bool writeRestartFile(const std::string& filename, int step);
    bool writeSummaryFile(const std::string& filename);
    bool writeInitFile(const std::string& filename);
    
    // Data accessors
    GridConfig getGridConfig() const;
    MaterialProperties getMaterialProperties(int i, int j, int k) const;
    FluidProperties getFluidProperties() const;
    SimulationConfig getSimulationConfig() const;
    
    // Grid data
    std::vector<double> getDX() const { return dx; }
    std::vector<double> getDY() const { return dy; }
    std::vector<double> getDZ() const { return dz; }
    std::vector<double> getTOPS() const { return tops; }
    std::vector<double> getPORO() const { return poro; }
    std::vector<double> getPERMX() const { return permx; }
    std::vector<double> getPERMY() const { return permy; }
    std::vector<double> getPERMZ() const { return permz; }
    
    // Well data
    struct WellSpec {
        std::string name;
        int i, j;
        double ref_depth;
        WellType type;
    };
    
    struct Completion {
        std::string well_name;
        int i, j, k1, k2;
        bool open;
        double well_index;
        double diameter;
        double skin;
    };
    
    struct WellControl {
        std::string well_name;
        std::string control_mode; // RATE, BHP, THP, etc.
        double target_rate;
        double target_pressure;
        double max_rate;
        double max_pressure;
    };
    
    std::vector<WellSpec> getWells() const { return wells; }
    std::vector<Completion> getCompletions() const { return completions; }
    std::vector<WellControl> getWellControls() const { return well_controls; }
    
private:
    std::vector<EclipseKeyword> keywords;
    std::map<std::string, size_t> keyword_map;
    
    // Grid dimensions
    int nx, ny, nz;
    int ncells;
    
    // Grid properties
    std::vector<double> dx, dy, dz, tops;
    std::vector<double> poro, permx, permy, permz;
    std::vector<int> actnum;
    
    // PVT data
    std::vector<std::vector<double>> pvto_data; // Oil PVT
    std::vector<std::vector<double>> pvtg_data; // Gas PVT
    std::vector<double> pvtw_data; // Water PVT
    std::vector<double> rock_data; // Rock compressibility
    
    // Relative permeability data
    std::vector<std::vector<double>> swof_data;
    std::vector<std::vector<double>> sgof_data;
    
    // Region data
    std::vector<int> satnum, pvtnum, eqlnum, fipnum;
    
    // Solution data
    std::vector<double> pressure, swat, sgas;
    
    // Wells
    std::vector<WellSpec> wells;
    std::vector<Completion> completions;
    std::vector<WellControl> well_controls;
    
    // Timesteps
    std::vector<double> timesteps;
    
    // Helper functions
    std::vector<double> parseDataArray(const std::string& keyword_name, int expected_size);
    std::vector<int> parseIntArray(const std::string& keyword_name, int expected_size);
    EclipseKeyword* findKeyword(const std::string& name);
    std::vector<std::string> tokenize(const std::string& line);
    double parseValue(const std::string& token);
    void expandArrayNotation(std::vector<std::string>& tokens);
};

} // namespace FSRM

#endif // ECLIPSE_IO_HPP
