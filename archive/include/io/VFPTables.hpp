#ifndef VFP_TABLES_HPP
#define VFP_TABLES_HPP

/**
 * @file VFPTables.hpp
 * @brief ECLIPSE-compatible Vertical Flow Performance tables
 * 
 * Implements VFP tables for production and injection wells used in ECLIPSE:
 * - VFPPROD: Production well VFP tables
 * - VFPINJ: Injection well VFP tables
 * - Multi-dimensional interpolation
 * - Gas lift optimization support
 * 
 * VFP tables relate wellhead conditions (THP, GOR, WCT, etc.) to 
 * bottomhole pressure at various flow rates.
 * 
 * ECLIPSE Keywords: VFPPROD, VFPINJ
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>
#include <cmath>
#include <algorithm>

namespace FSRM {

/**
 * @brief Flow rate type for VFP tables
 */
enum class VFPFlowType {
    OIL,        ///< Oil rate as primary flow variable
    LIQ,        ///< Total liquid rate (oil + water)
    GAS,        ///< Gas rate as primary flow variable
    WATER       ///< Water rate as primary flow variable
};

/**
 * @brief Water cut definition type
 */
enum class VFPWaterCutType {
    WCT,        ///< Water cut = water rate / (oil + water rate)
    WGR         ///< Water-gas ratio = water rate / gas rate
};

/**
 * @brief Gas-liquid ratio type
 */
enum class VFPGasLiquidType {
    GOR,        ///< Gas-oil ratio
    GLR,        ///< Gas-liquid ratio (total)
    OGR         ///< Oil-gas ratio
};

/**
 * @brief Artificial lift type
 */
enum class VFPArtificialLiftType {
    NONE,       ///< No artificial lift
    GASLIFT,    ///< Gas lift (ALQ = gas lift rate)
    ESP,        ///< Electric submersible pump (ALQ = pump speed)
    CHOKE       ///< Surface choke (ALQ = choke setting)
};

/**
 * @brief Unit system for VFP tables
 */
enum class VFPUnitSystem {
    METRIC,     ///< Metric units (bar, sm³/day)
    FIELD,      ///< Field units (psi, STB/day)
    LAB,        ///< Lab units
    DEFAULT     ///< Inherit from simulation
};

/**
 * @brief Single VFP data point
 */
struct VFPDataPoint {
    double thp;         ///< Tubing head pressure
    double flow_rate;   ///< Flow rate
    double wct;         ///< Water cut or WGR
    double gor;         ///< GOR or GLR
    double alq;         ///< Artificial lift quantity
    double bhp;         ///< Calculated BHP
};

/**
 * @brief Production VFP table (VFPPROD)
 * 
 * 5-dimensional table: BHP = f(THP, flow_rate, WCT, GOR, ALQ)
 */
class VFPProdTable {
public:
    VFPProdTable(int table_num = 1);
    
    // Configuration
    void setTableNumber(int num) { table_number_ = num; }
    void setDatum(double depth) { datum_depth_ = depth; }
    void setFlowType(VFPFlowType type) { flow_type_ = type; }
    void setWaterCutType(VFPWaterCutType type) { wct_type_ = type; }
    void setGasLiquidType(VFPGasLiquidType type) { glr_type_ = type; }
    void setArtificialLiftType(VFPArtificialLiftType type) { alq_type_ = type; }
    void setUnitSystem(VFPUnitSystem units) { unit_system_ = units; }
    
    // Set axis values
    void setTHPValues(const std::vector<double>& thp);
    void setFlowRateValues(const std::vector<double>& rates);
    void setWaterCutValues(const std::vector<double>& wct);
    void setGORValues(const std::vector<double>& gor);
    void setALQValues(const std::vector<double>& alq);
    
    // Set BHP data (5D array flattened)
    void setBHPData(const std::vector<double>& bhp);
    
    // Alternative: set BHP from individual points
    void addDataPoint(const VFPDataPoint& point);
    void buildTableFromPoints();
    
    /**
     * @brief Calculate BHP from VFP table
     * 
     * Uses multi-dimensional linear interpolation.
     * 
     * @param thp Tubing head pressure
     * @param flow_rate Flow rate
     * @param wct Water cut or WGR
     * @param gor GOR or GLR
     * @param alq Artificial lift quantity
     * @return Calculated BHP
     */
    double calculateBHP(double thp, double flow_rate, double wct, 
                        double gor, double alq = 0.0) const;
    
    /**
     * @brief Calculate flow rate for given BHP and THP
     * 
     * Inverse lookup - finds flow rate that gives target BHP.
     * 
     * @param bhp Target bottomhole pressure
     * @param thp Tubing head pressure
     * @param wct Water cut
     * @param gor GOR
     * @param alq Artificial lift quantity
     * @return Flow rate (or 0 if not achievable)
     */
    double calculateFlowRate(double bhp, double thp, double wct,
                             double gor, double alq = 0.0) const;
    
    /**
     * @brief Find optimal ALQ for gas lift
     * 
     * @param bhp Target BHP
     * @param thp Tubing head pressure
     * @param wct Water cut
     * @param gor GOR
     * @param max_alq Maximum available gas lift
     * @return Optimal ALQ value
     */
    double optimizeGasLift(double bhp, double thp, double wct,
                           double gor, double max_alq) const;
    
    // Accessors
    int getTableNumber() const { return table_number_; }
    double getDatumDepth() const { return datum_depth_; }
    VFPFlowType getFlowType() const { return flow_type_; }
    
    const std::vector<double>& getTHPValues() const { return thp_values_; }
    const std::vector<double>& getFlowRateValues() const { return rate_values_; }
    const std::vector<double>& getWaterCutValues() const { return wct_values_; }
    const std::vector<double>& getGORValues() const { return gor_values_; }
    const std::vector<double>& getALQValues() const { return alq_values_; }
    
    // Validation
    bool isValid() const;
    void printTable() const;
    
private:
    int table_number_;
    double datum_depth_;
    
    // Types
    VFPFlowType flow_type_;
    VFPWaterCutType wct_type_;
    VFPGasLiquidType glr_type_;
    VFPArtificialLiftType alq_type_;
    VFPUnitSystem unit_system_;
    
    // Axis values
    std::vector<double> thp_values_;
    std::vector<double> rate_values_;
    std::vector<double> wct_values_;
    std::vector<double> gor_values_;
    std::vector<double> alq_values_;
    
    // BHP data: bhp_[thp_idx][rate_idx][wct_idx][gor_idx][alq_idx]
    // Stored as flattened 1D array
    std::vector<double> bhp_data_;
    
    // Data points for table building
    std::vector<VFPDataPoint> data_points_;
    
    // Helper functions
    size_t getIndex(size_t i_thp, size_t i_rate, size_t i_wct, 
                    size_t i_gor, size_t i_alq) const;
    
    std::pair<size_t, double> findBracket(double value, 
                                          const std::vector<double>& axis) const;
    
    double interpolate5D(double thp, double rate, double wct, 
                         double gor, double alq) const;
};

/**
 * @brief Injection VFP table (VFPINJ)
 * 
 * 2-dimensional table: BHP = f(THP, injection_rate)
 */
class VFPInjTable {
public:
    VFPInjTable(int table_num = 1);
    
    // Configuration
    void setTableNumber(int num) { table_number_ = num; }
    void setDatum(double depth) { datum_depth_ = depth; }
    void setFlowType(VFPFlowType type) { flow_type_ = type; }
    void setUnitSystem(VFPUnitSystem units) { unit_system_ = units; }
    
    // Set axis values
    void setTHPValues(const std::vector<double>& thp);
    void setFlowRateValues(const std::vector<double>& rates);
    
    // Set BHP data (2D array flattened: thp x rate)
    void setBHPData(const std::vector<double>& bhp);
    
    /**
     * @brief Calculate BHP from VFP table
     * 
     * @param thp Tubing head pressure
     * @param flow_rate Injection rate
     * @return Calculated BHP
     */
    double calculateBHP(double thp, double flow_rate) const;
    
    /**
     * @brief Calculate injection rate for given BHP and THP
     * 
     * @param bhp Target bottomhole pressure
     * @param thp Tubing head pressure
     * @return Injection rate (or 0 if not achievable)
     */
    double calculateFlowRate(double bhp, double thp) const;
    
    // Accessors
    int getTableNumber() const { return table_number_; }
    double getDatumDepth() const { return datum_depth_; }
    
    const std::vector<double>& getTHPValues() const { return thp_values_; }
    const std::vector<double>& getFlowRateValues() const { return rate_values_; }
    
    bool isValid() const;
    void printTable() const;
    
private:
    int table_number_;
    double datum_depth_;
    VFPFlowType flow_type_;
    VFPUnitSystem unit_system_;
    
    std::vector<double> thp_values_;
    std::vector<double> rate_values_;
    std::vector<double> bhp_data_;  // [thp_idx * n_rates + rate_idx]
    
    double interpolate2D(double thp, double rate) const;
};

/**
 * @brief Generate VFP table from wellbore model
 * 
 * Creates VFP table by running steady-state wellbore calculations.
 */
class VFPTableGenerator {
public:
    VFPTableGenerator();
    
    // Wellbore properties
    void setWellboreRadius(double r) { wellbore_radius_ = r; }
    void setRoughness(double eps) { roughness_ = eps; }
    void setTotalDepth(double tvd) { total_depth_ = tvd; }
    void setDeviation(const std::vector<double>& md, 
                      const std::vector<double>& inc);
    
    // Fluid properties
    void setOilProperties(double density, double viscosity);
    void setWaterProperties(double density, double viscosity);
    void setGasProperties(double density, double viscosity);
    
    // Pressure/temperature gradient
    void setSurfaceTemperature(double T) { surface_temp_ = T; }
    void setGeothermalGradient(double grad) { geothermal_grad_ = grad; }
    
    /**
     * @brief Generate production VFP table
     * 
     * @param thp_range THP values
     * @param rate_range Flow rates
     * @param wct_range Water cuts
     * @param gor_range GOR values
     * @param alq_range ALQ values (for gas lift)
     * @return Generated VFPPROD table
     */
    VFPProdTable generateProductionTable(
        const std::vector<double>& thp_range,
        const std::vector<double>& rate_range,
        const std::vector<double>& wct_range,
        const std::vector<double>& gor_range,
        const std::vector<double>& alq_range) const;
    
    /**
     * @brief Generate injection VFP table
     * 
     * @param thp_range THP values
     * @param rate_range Injection rates
     * @param fluid Injection fluid type
     * @return Generated VFPINJ table
     */
    VFPInjTable generateInjectionTable(
        const std::vector<double>& thp_range,
        const std::vector<double>& rate_range,
        const std::string& fluid = "WATER") const;
    
private:
    double wellbore_radius_;
    double roughness_;
    double total_depth_;
    
    std::vector<double> md_;    // Measured depth
    std::vector<double> inc_;   // Inclination
    
    double oil_density_, oil_viscosity_;
    double water_density_, water_viscosity_;
    double gas_density_, gas_viscosity_;
    
    double surface_temp_;
    double geothermal_grad_;
    
    // Pressure drop calculations
    double calculateFrictionGradient(double rate, double density, 
                                     double viscosity) const;
    double calculateHydrostaticGradient(double density, double inc) const;
    double calculateTotalPressureDrop(double thp, double rate, double wct,
                                      double gor, double alq) const;
};

/**
 * @brief Gas lift optimization using VFP tables
 */
class GasLiftOptimizer {
public:
    GasLiftOptimizer();
    
    // Set VFP table
    void setVFPTable(const VFPProdTable& table);
    
    // Economic parameters
    void setOilPrice(double price) { oil_price_ = price; }
    void setGasPrice(double price) { gas_price_ = price; }
    void setGasLiftCost(double cost) { gaslift_cost_ = cost; }
    
    // Constraints
    void setMaxGasLift(double max) { max_gaslift_ = max; }
    void setMinBHP(double min) { min_bhp_ = min; }
    void setMinTHP(double min) { min_thp_ = min; }
    
    /**
     * @brief Optimize gas lift for single well
     * 
     * @param thp Current THP
     * @param wct Current water cut
     * @param gor Current GOR
     * @param reservoir_pressure Reservoir pressure
     * @return Optimal gas lift rate and corresponding oil rate
     */
    std::pair<double, double> optimizeSingleWell(
        double thp, double wct, double gor,
        double reservoir_pressure) const;
    
    /**
     * @brief Optimize gas lift for multiple wells
     * 
     * Distributes available gas lift to maximize total revenue.
     * 
     * @param wells Vector of (VFP table, thp, wct, gor, Pr) tuples
     * @param total_gas_available Total available gas lift
     * @return Vector of optimal ALQ for each well
     */
    std::vector<double> optimizeField(
        const std::vector<std::tuple<const VFPProdTable*, double, double, double, double>>& wells,
        double total_gas_available) const;
    
private:
    const VFPProdTable* vfp_table_;
    
    double oil_price_;
    double gas_price_;
    double gaslift_cost_;
    double max_gaslift_;
    double min_bhp_;
    double min_thp_;
    
    // Calculate net revenue for given conditions
    double calculateNetRevenue(double oil_rate, double gas_rate, 
                               double gaslift_rate) const;
};

/**
 * @brief VFP table manager
 */
class VFPTableManager {
public:
    VFPTableManager();
    
    // Add tables
    void addProductionTable(std::unique_ptr<VFPProdTable> table);
    void addInjectionTable(std::unique_ptr<VFPInjTable> table);
    
    // Parse from ECLIPSE format
    void parseVFPPROD(const std::vector<std::string>& lines);
    void parseVFPINJ(const std::vector<std::string>& lines);
    
    // Access tables
    const VFPProdTable* getProductionTable(int table_num) const;
    const VFPInjTable* getInjectionTable(int table_num) const;
    
    // Calculate BHP for wells
    double calculateProductionBHP(int table_num, double thp, double rate,
                                  double wct, double gor, double alq = 0.0) const;
    double calculateInjectionBHP(int table_num, double thp, double rate) const;
    
    // Inverse calculations
    double calculateProductionRate(int table_num, double bhp, double thp,
                                   double wct, double gor, double alq = 0.0) const;
    double calculateInjectionRate(int table_num, double bhp, double thp) const;
    
    // Output
    void writeVFPPROD(const std::string& filename, int table_num) const;
    void writeVFPINJ(const std::string& filename, int table_num) const;
    
private:
    std::map<int, std::unique_ptr<VFPProdTable>> prod_tables_;
    std::map<int, std::unique_ptr<VFPInjTable>> inj_tables_;
};

/**
 * @brief Parse VFP flow type from string
 */
VFPFlowType parseVFPFlowType(const std::string& str);

/**
 * @brief Parse VFP water cut type from string
 */
VFPWaterCutType parseVFPWaterCutType(const std::string& str);

/**
 * @brief Parse VFP gas-liquid type from string
 */
VFPGasLiquidType parseVFPGasLiquidType(const std::string& str);

/**
 * @brief Parse VFP artificial lift type from string
 */
VFPArtificialLiftType parseVFPArtificialLiftType(const std::string& str);

} // namespace FSRM

#endif // VFP_TABLES_HPP
