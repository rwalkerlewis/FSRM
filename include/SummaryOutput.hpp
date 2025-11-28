#ifndef SUMMARY_OUTPUT_HPP
#define SUMMARY_OUTPUT_HPP

/**
 * @file SummaryOutput.hpp
 * @brief ECLIPSE-compatible summary output
 * 
 * Implements full ECLIPSE summary file format and keywords:
 * - Well summary vectors (WOPR, WWPR, WGPR, WBHP, etc.)
 * - Group summary vectors (GOPR, GWPR, GGPR, etc.)
 * - Field summary vectors (FOPR, FWPR, FGPR, etc.)
 * - Region summary vectors (ROPR, RWPR, etc.)
 * - Block summary vectors (BPRES, BSWAT, etc.)
 * - Connection summary vectors (COFR, CWFR, etc.)
 * - Aquifer summary vectors (AAQP, AAQT, etc.)
 * - Performance vectors (CPU, ELAPSED, NEWTON, etc.)
 * 
 * ECLIPSE Keywords: SUMMARY, RPTSUM
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <functional>

namespace FSRM {

/**
 * @brief Summary vector category
 */
enum class SummaryCategory {
    FIELD,          ///< Field total (F*)
    GROUP,          ///< Group total (G*)
    WELL,           ///< Well (W*)
    BLOCK,          ///< Grid block (B*)
    REGION,         ///< Region (R*)
    CONNECTION,     ///< Well connection (C*)
    SEGMENT,        ///< Multi-segment well (S*)
    AQUIFER,        ///< Aquifer (A*)
    NETWORK,        ///< Network node (N*)
    PERFORMANCE     ///< Performance/miscellaneous (no prefix)
};

/**
 * @brief Summary vector type (what it measures)
 */
enum class SummaryType {
    // Rates
    OIL_PRODUCTION_RATE,
    WATER_PRODUCTION_RATE,
    GAS_PRODUCTION_RATE,
    LIQUID_PRODUCTION_RATE,
    RESERVOIR_VOIDAGE_RATE,
    OIL_INJECTION_RATE,
    WATER_INJECTION_RATE,
    GAS_INJECTION_RATE,
    
    // Cumulative
    OIL_PRODUCTION_TOTAL,
    WATER_PRODUCTION_TOTAL,
    GAS_PRODUCTION_TOTAL,
    OIL_INJECTION_TOTAL,
    WATER_INJECTION_TOTAL,
    GAS_INJECTION_TOTAL,
    
    // Pressures
    BHP,
    THP,
    AVERAGE_PRESSURE,
    BLOCK_PRESSURE,
    
    // Saturations
    WATER_SATURATION,
    OIL_SATURATION,
    GAS_SATURATION,
    
    // Ratios
    WATER_CUT,
    GAS_OIL_RATIO,
    WATER_GAS_RATIO,
    
    // Well performance
    PRODUCTIVITY_INDEX,
    WELL_POTENTIAL,
    GAS_LIFT_RATE,
    
    // Aquifer
    AQUIFER_PRESSURE,
    AQUIFER_INFLUX_RATE,
    AQUIFER_INFLUX_TOTAL,
    
    // Miscellaneous
    TIME,
    TIMESTEP,
    NEWTON_ITERATIONS,
    LINEAR_ITERATIONS,
    CPU_TIME,
    ELAPSED_TIME,
    
    // Region totals
    REGION_PORE_VOLUME,
    REGION_HYDROCARBON_PV,
    REGION_PRESSURE_AVERAGE,
    
    // Network
    NODE_PRESSURE,
    BRANCH_FLOW_RATE
};

/**
 * @brief Single summary vector specification
 */
struct SummaryVector {
    std::string keyword;            ///< ECLIPSE keyword (e.g., "WOPR")
    SummaryCategory category;       ///< Category (WELL, FIELD, etc.)
    SummaryType type;               ///< What is measured
    std::string entity_name;        ///< Well/group/region name
    int entity_id = 0;              ///< Block i,j,k or region number
    int block_i = 0, block_j = 0, block_k = 0;  ///< For block vectors
    std::string unit;               ///< Output unit
    
    // Time series data
    std::vector<double> values;
    
    SummaryVector() = default;
    SummaryVector(const std::string& kw, SummaryCategory cat, SummaryType t,
                  const std::string& entity = "", const std::string& u = "")
        : keyword(kw), category(cat), type(t), entity_name(entity), unit(u) {}
};

/**
 * @brief Summary data for a single time step
 */
struct SummaryStep {
    double time;                    ///< Simulation time (days)
    double dt;                      ///< Time step size
    int report_step;                ///< Report step number
    int ministep;                   ///< Mini-step within report step
    std::map<std::string, double> values;  ///< Vector values
};

/**
 * @brief ECLIPSE summary specification manager
 * 
 * Parses SUMMARY section and manages what vectors to output.
 */
class SummarySpec {
public:
    SummarySpec();
    
    // Add vectors by keyword
    void addKeyword(const std::string& keyword);
    void addWellKeyword(const std::string& keyword, const std::string& well);
    void addGroupKeyword(const std::string& keyword, const std::string& group);
    void addBlockKeyword(const std::string& keyword, int i, int j, int k);
    void addRegionKeyword(const std::string& keyword, int region);
    void addAquiferKeyword(const std::string& keyword, int aquifer);
    void addConnectionKeyword(const std::string& keyword, const std::string& well,
                               int i, int j, int k);
    
    // Add common sets
    void addAllWellVectors(const std::string& well);
    void addAllFieldVectors();
    void addAllAquiferVectors(int aquifer);
    void addDefaultVectors();
    
    // Query
    bool hasVector(const std::string& key) const;
    std::vector<std::string> getAllKeywords() const;
    std::vector<SummaryVector> getVectors() const { return vectors_; }
    
    // Parse from ECLIPSE format
    void parse(const std::vector<std::string>& summary_lines);
    
    // Set well/group/region lists for expansion of wildcards
    void setWells(const std::vector<std::string>& wells) { wells_ = wells; }
    void setGroups(const std::vector<std::string>& groups) { groups_ = groups; }
    void setRegions(const std::vector<int>& regions) { regions_ = regions; }
    void setAquifers(const std::vector<int>& aquifers) { aquifers_ = aquifers; }
    
private:
    std::vector<SummaryVector> vectors_;
    std::set<std::string> keywords_;
    
    std::vector<std::string> wells_;
    std::vector<std::string> groups_;
    std::vector<int> regions_;
    std::vector<int> aquifers_;
    
    // Parse single keyword
    SummaryVector parseKeyword(const std::string& keyword) const;
    void expandWildcards(const std::string& pattern);
};

/**
 * @brief Summary data collector
 * 
 * Collects data from simulation and stores time series.
 */
class SummaryCollector {
public:
    SummaryCollector();
    
    void setSpec(const SummarySpec& spec);
    
    // Record data for current timestep
    void startStep(double time, double dt, int report_step);
    void endStep();
    
    // Well data
    void recordWellRate(const std::string& well, SummaryType type, double value);
    void recordWellPressure(const std::string& well, SummaryType type, double value);
    void recordWellCumulative(const std::string& well, SummaryType type, double value);
    
    // Group data
    void recordGroupRate(const std::string& group, SummaryType type, double value);
    void recordGroupCumulative(const std::string& group, SummaryType type, double value);
    
    // Field data
    void recordFieldRate(SummaryType type, double value);
    void recordFieldCumulative(SummaryType type, double value);
    void recordFieldPressure(double value);
    
    // Block data
    void recordBlockData(int i, int j, int k, SummaryType type, double value);
    
    // Region data  
    void recordRegionData(int region, SummaryType type, double value);
    
    // Aquifer data
    void recordAquiferData(int aquifer, SummaryType type, double value);
    
    // Performance data
    void recordNewtonIterations(int count);
    void recordLinearIterations(int count);
    void recordCPUTime(double seconds);
    void recordTimestepSize(double dt);
    
    // Access time series
    std::vector<double> getTimeSeries(const std::string& keyword) const;
    std::vector<double> getTimes() const;
    double getLatestValue(const std::string& keyword) const;
    
    // Get all steps
    const std::vector<SummaryStep>& getSteps() const { return steps_; }
    int getNumSteps() const { return static_cast<int>(steps_.size()); }
    
private:
    SummarySpec spec_;
    std::vector<SummaryStep> steps_;
    SummaryStep current_step_;
    bool collecting_ = false;
    
    // Cumulative counters
    std::map<std::string, double> cumulatives_;
    
    std::string makeKey(const std::string& keyword, const std::string& entity = "") const;
};

/**
 * @brief Summary file writer
 * 
 * Writes ECLIPSE-format summary files (.SMSPEC, .UNSMRY/.SMRY)
 */
class SummaryWriter {
public:
    SummaryWriter(const std::string& case_name);
    
    // Initialize with vectors to output
    void initialize(const SummarySpec& spec);
    
    // Write step data
    void writeStep(const SummaryStep& step);
    
    // Finalize files
    void finalize();
    
    // Output format options
    void setUnformatted(bool unformatted) { unformatted_ = unformatted; }
    void setBinaryOutput(bool binary) { binary_ = binary; }
    
    // Write ASCII summary (alternative)
    void writeASCII(const std::string& filename, 
                    const std::vector<SummaryStep>& steps) const;
    
    // Write CSV summary
    void writeCSV(const std::string& filename,
                  const std::vector<SummaryStep>& steps) const;
    
private:
    std::string case_name_;
    std::string smspec_file_;
    std::string smry_file_;
    
    bool unformatted_ = true;
    bool binary_ = false;
    
    std::vector<SummaryVector> vectors_;
    std::ofstream smry_stream_;
    
    void writeSmspec(const SummarySpec& spec);
    void writeBinaryRecord(const std::vector<float>& data);
};

/**
 * @brief Summary file reader (for restart)
 */
class SummaryReader {
public:
    SummaryReader(const std::string& case_name);
    
    bool read();
    
    // Access data
    std::vector<double> getTimeSeries(const std::string& keyword) const;
    std::vector<double> getTimes() const;
    
    // Get available keywords
    std::vector<std::string> getKeywords() const;
    
    // Get data at specific time
    double getValue(const std::string& keyword, double time) const;
    
    // Get all steps
    const std::vector<SummaryStep>& getSteps() const { return steps_; }
    
private:
    std::string case_name_;
    std::vector<SummaryVector> vectors_;
    std::vector<SummaryStep> steps_;
    
    bool readSmspec();
    bool readSmry();
};

/**
 * @brief Commonly used ECLIPSE summary keywords
 */
namespace SummaryKeywords {
    // Field vectors
    const std::string FOPR = "FOPR";     // Field Oil Production Rate
    const std::string FWPR = "FWPR";     // Field Water Production Rate
    const std::string FGPR = "FGPR";     // Field Gas Production Rate
    const std::string FLPR = "FLPR";     // Field Liquid Production Rate
    const std::string FOPT = "FOPT";     // Field Oil Production Total
    const std::string FWPT = "FWPT";     // Field Water Production Total
    const std::string FGPT = "FGPT";     // Field Gas Production Total
    const std::string FWIR = "FWIR";     // Field Water Injection Rate
    const std::string FGIR = "FGIR";     // Field Gas Injection Rate
    const std::string FWIT = "FWIT";     // Field Water Injection Total
    const std::string FGIT = "FGIT";     // Field Gas Injection Total
    const std::string FWCT = "FWCT";     // Field Water Cut
    const std::string FGOR = "FGOR";     // Field Gas-Oil Ratio
    const std::string FPR  = "FPR";      // Field Pressure (Average)
    
    // Well vectors  
    const std::string WOPR = "WOPR";     // Well Oil Production Rate
    const std::string WWPR = "WWPR";     // Well Water Production Rate
    const std::string WGPR = "WGPR";     // Well Gas Production Rate
    const std::string WLPR = "WLPR";     // Well Liquid Production Rate
    const std::string WOPT = "WOPT";     // Well Oil Production Total
    const std::string WWPT = "WWPT";     // Well Water Production Total
    const std::string WGPT = "WGPT";     // Well Gas Production Total
    const std::string WWIR = "WWIR";     // Well Water Injection Rate
    const std::string WGIR = "WGIR";     // Well Gas Injection Rate
    const std::string WWIT = "WWIT";     // Well Water Injection Total
    const std::string WGIT = "WGIT";     // Well Gas Injection Total
    const std::string WBHP = "WBHP";     // Well BHP
    const std::string WTHP = "WTHP";     // Well THP
    const std::string WWCT = "WWCT";     // Well Water Cut
    const std::string WGOR = "WGOR";     // Well Gas-Oil Ratio
    const std::string WPI  = "WPI";      // Well Productivity Index
    const std::string WSTAT = "WSTAT";   // Well Status
    const std::string WGLIR = "WGLIR";   // Well Gas Lift Rate
    
    // Group vectors
    const std::string GOPR = "GOPR";     // Group Oil Production Rate
    const std::string GWPR = "GWPR";     // Group Water Production Rate
    const std::string GGPR = "GGPR";     // Group Gas Production Rate
    const std::string GOPT = "GOPT";     // Group Oil Production Total
    const std::string GWPT = "GWPT";     // Group Water Production Total
    const std::string GGPT = "GGPT";     // Group Gas Production Total
    const std::string GWIR = "GWIR";     // Group Water Injection Rate
    const std::string GGIR = "GGIR";     // Group Gas Injection Rate
    
    // Block vectors
    const std::string BPR  = "BPR";      // Block Pressure
    const std::string BSWAT = "BSWAT";   // Block Water Saturation
    const std::string BSGAS = "BSGAS";   // Block Gas Saturation
    const std::string BOSAT = "BOSAT";   // Block Oil Saturation
    
    // Region vectors
    const std::string RPR  = "RPR";      // Region Pressure
    const std::string ROPR = "ROPR";     // Region Oil Production Rate
    const std::string RWPR = "RWPR";     // Region Water Production Rate
    const std::string RGPR = "RGPR";     // Region Gas Production Rate
    const std::string ROPT = "ROPT";     // Region Oil Production Total
    
    // Aquifer vectors
    const std::string AAQP = "AAQP";     // Aquifer Pressure
    const std::string AAQR = "AAQR";     // Aquifer Influx Rate
    const std::string AAQT = "AAQT";     // Aquifer Influx Total
    
    // Connection vectors
    const std::string COFR = "COFR";     // Connection Oil Flow Rate
    const std::string CWFR = "CWFR";     // Connection Water Flow Rate
    const std::string CGFR = "CGFR";     // Connection Gas Flow Rate
    
    // Performance vectors
    const std::string TIME = "TIME";     // Simulation Time
    const std::string TCPU = "TCPU";     // CPU Time
    const std::string TIMESTEP = "TIMESTEP"; // Time Step Size
    const std::string NEWTON = "NEWTON"; // Newton Iterations
    const std::string MLINEAR = "MLINEAR"; // Linear Iterations
}

/**
 * @brief Convenient wrapper for summary output
 */
class SummaryOutput {
public:
    SummaryOutput(const std::string& case_name);
    
    // Setup
    void setWells(const std::vector<std::string>& wells);
    void setGroups(const std::vector<std::string>& groups);
    void setRegions(const std::vector<int>& regions);
    void setAquifers(const std::vector<int>& aquifers);
    
    // Parse SUMMARY section
    void parse(const std::vector<std::string>& lines);
    
    // Add common outputs
    void addDefaultOutput();
    void addWellOutput(const std::string& well);
    void addAllWellsOutput();
    void addFieldOutput();
    void addAquiferOutput();
    
    // Initialize output files
    void initialize();
    
    // Record timestep data
    void beginStep(double time, double dt, int report_step);
    
    // Record various data types
    void recordWell(const std::string& well, 
                    double oil_rate, double water_rate, double gas_rate,
                    double bhp, double thp,
                    double oil_total, double water_total, double gas_total);
    
    void recordInjector(const std::string& well,
                        double water_rate, double gas_rate,
                        double bhp);
    
    void recordGroup(const std::string& group,
                     double oil_rate, double water_rate, double gas_rate);
    
    void recordField(double oil_rate, double water_rate, double gas_rate,
                     double oil_total, double water_total, double gas_total,
                     double avg_pressure);
    
    void recordBlock(int i, int j, int k, double pressure, 
                     double Sw, double So, double Sg);
    
    void recordRegion(int region, double oil_rate, double water_rate, double gas_rate);
    
    void recordAquifer(int aq, double pressure, double influx_rate, double influx_total);
    
    void recordPerformance(int newton_iter, int linear_iter, double cpu_time);
    
    void endStep();
    
    // Finalize and close files
    void finalize();
    
    // Write human-readable summary
    void writeTextSummary(const std::string& filename) const;
    void writeCSV(const std::string& filename) const;
    
    // Access collected data
    std::vector<double> getTimeSeries(const std::string& keyword) const;
    double getLatestValue(const std::string& keyword) const;
    
private:
    std::string case_name_;
    SummarySpec spec_;
    SummaryCollector collector_;
    std::unique_ptr<SummaryWriter> writer_;
    
    std::vector<std::string> wells_;
    std::vector<std::string> groups_;
    std::vector<int> regions_;
    std::vector<int> aquifers_;
    
    bool initialized_ = false;
};

// Parse helpers
SummaryCategory parseSummaryCategory(char prefix);
SummaryType parseSummaryType(const std::string& keyword);
std::string summaryTypeToString(SummaryType type);

} // namespace FSRM

#endif // SUMMARY_OUTPUT_HPP
