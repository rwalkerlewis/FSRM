#ifndef GROUP_CONTROL_HPP
#define GROUP_CONTROL_HPP

/**
 * @file GroupControl.hpp
 * @brief ECLIPSE-compatible well group and field management
 * 
 * Implements hierarchical well group control used in ECLIPSE:
 * - GRUPTREE: Group hierarchy definition
 * - GCONPROD: Production group constraints
 * - GCONINJE: Injection group constraints  
 * - GCONSALE: Gas sales contracts
 * - GCONSUMP: Gas consumption
 * - GEFAC: Group efficiency factors
 * - Guide rate allocation
 * 
 * ECLIPSE Keywords: GRUPTREE, GCONPROD, GCONINJE, GCONSALE, GCONSUMP, GEFAC
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <set>
#include <functional>

namespace FSRM {

/**
 * @brief Group control mode for production
 */
enum class GroupProdControlMode {
    NONE,       ///< No group control
    ORAT,       ///< Target oil rate
    WRAT,       ///< Target water rate
    GRAT,       ///< Target gas rate
    LRAT,       ///< Target liquid rate
    RESV,       ///< Target reservoir voidage
    PRBL,       ///< Pressure balance
    FLD         ///< Field control (highest level)
};

/**
 * @brief Group control mode for injection
 */
enum class GroupInjControlMode {
    NONE,       ///< No group control
    RATE,       ///< Target injection rate
    RESV,       ///< Target reservoir voidage replacement
    REIN,       ///< Reinjection of produced gas/water
    VREP,       ///< Voidage replacement
    WGRA,       ///< Water-gas ratio
    FLD         ///< Field control
};

/**
 * @brief Injection phase
 */
enum class InjectionPhase {
    WATER,
    GAS,
    OIL
};

/**
 * @brief Guide rate type for allocation
 */
enum class GuideRateType {
    OIL,        ///< Based on oil potential
    WATER,      ///< Based on water potential
    GAS,        ///< Based on gas potential
    LIQ,        ///< Based on liquid potential
    COMB,       ///< Combined guide rate
    NONE,       ///< No guide rate (equal allocation)
    FORM        ///< Formula-based
};

/**
 * @brief Production constraints for a group
 */
struct GroupProdConstraints {
    GroupProdControlMode control_mode = GroupProdControlMode::NONE;
    double oil_rate_target = 1e20;      ///< Target oil rate (m³/s)
    double water_rate_target = 1e20;    ///< Target water rate (m³/s)
    double gas_rate_target = 1e20;      ///< Target gas rate (m³/s)
    double liquid_rate_target = 1e20;   ///< Target liquid rate (m³/s)
    double reservoir_voidage_target = 1e20;  ///< Target voidage (m³/s)
    double oil_rate_max = 1e20;         ///< Maximum oil rate
    double water_rate_max = 1e20;       ///< Maximum water rate
    double gas_rate_max = 1e20;         ///< Maximum gas rate
    double liquid_rate_max = 1e20;      ///< Maximum liquid rate
    double water_cut_max = 1.0;         ///< Maximum water cut
    double gas_oil_ratio_max = 1e20;    ///< Maximum GOR
    GuideRateType guide_rate_type = GuideRateType::OIL;
    bool guide_rate_defined = false;
    double guide_rate_phase_weighting = 1.0;
};

/**
 * @brief Injection constraints for a group
 */
struct GroupInjConstraints {
    GroupInjControlMode control_mode = GroupInjControlMode::NONE;
    InjectionPhase phase = InjectionPhase::WATER;
    double surface_rate_target = 1e20;  ///< Target surface injection rate
    double reservoir_rate_target = 1e20;///< Target reservoir injection rate
    double vrep_target = 1.0;           ///< Voidage replacement target fraction
    double reinj_fraction = 1.0;        ///< Reinjection fraction
    double surface_rate_max = 1e20;     ///< Maximum surface rate
    double reservoir_rate_max = 1e20;   ///< Maximum reservoir rate
    double bhp_max = 1e20;              ///< Maximum BHP for all wells
    GuideRateType guide_rate_type = GuideRateType::NONE;
};

/**
 * @brief Gas sales contract specification
 */
struct GasSalesContract {
    double sales_rate_max = 1e20;       ///< Maximum sales rate (m³/s)
    double sales_rate_min = 0.0;        ///< Minimum sales rate (m³/s)
    double sales_rate_target = 0.0;     ///< Target sales rate
    double price = 1.0;                 ///< Gas price
    bool has_contract = false;
};

/**
 * @brief Gas consumption specification
 */
struct GasConsumption {
    double fuel_rate = 0.0;             ///< Fuel gas consumption (m³/s)
    double import_rate = 0.0;           ///< Gas import rate
    bool has_consumption = false;
};

/**
 * @brief Well group node in hierarchy
 */
class WellGroup {
public:
    WellGroup(const std::string& name, const std::string& parent = "FIELD");
    
    // Hierarchy
    void setParent(const std::string& parent) { parent_name_ = parent; }
    void addChild(const std::string& child) { children_.insert(child); }
    void addWell(const std::string& well) { wells_.insert(well); }
    void removeWell(const std::string& well) { wells_.erase(well); }
    
    const std::string& getName() const { return name_; }
    const std::string& getParent() const { return parent_name_; }
    const std::set<std::string>& getChildren() const { return children_; }
    const std::set<std::string>& getWells() const { return wells_; }
    
    bool isField() const { return name_ == "FIELD"; }
    bool hasChildren() const { return !children_.empty(); }
    bool hasWells() const { return !wells_.empty(); }
    
    // Constraints
    void setProductionConstraints(const GroupProdConstraints& c) { prod_constraints_ = c; }
    void setInjectionConstraints(const GroupInjConstraints& c) { inj_constraints_ = c; }
    void setGasSalesContract(const GasSalesContract& c) { gas_sales_ = c; }
    void setGasConsumption(const GasConsumption& c) { gas_consumption_ = c; }
    void setEfficiencyFactor(double eff) { efficiency_factor_ = eff; }
    
    const GroupProdConstraints& getProductionConstraints() const { return prod_constraints_; }
    const GroupInjConstraints& getInjectionConstraints() const { return inj_constraints_; }
    const GasSalesContract& getGasSalesContract() const { return gas_sales_; }
    const GasConsumption& getGasConsumption() const { return gas_consumption_; }
    double getEfficiencyFactor() const { return efficiency_factor_; }
    
    // Current production
    void updateProduction(double oil, double water, double gas, double resv);
    double getOilRate() const { return current_oil_rate_; }
    double getWaterRate() const { return current_water_rate_; }
    double getGasRate() const { return current_gas_rate_; }
    double getLiquidRate() const { return current_oil_rate_ + current_water_rate_; }
    double getReservoirVoidage() const { return current_resv_rate_; }
    
    // Current injection
    void updateInjection(double rate, double resv);
    double getInjectionRate() const { return current_inj_rate_; }
    double getInjectionResvRate() const { return current_inj_resv_rate_; }
    
    // Cumulative production
    void updateCumulative(double dt);
    double getCumulativeOil() const { return cumulative_oil_; }
    double getCumulativeWater() const { return cumulative_water_; }
    double getCumulativeGas() const { return cumulative_gas_; }
    double getCumulativeInjection() const { return cumulative_injection_; }
    
private:
    std::string name_;
    std::string parent_name_;
    std::set<std::string> children_;
    std::set<std::string> wells_;
    
    GroupProdConstraints prod_constraints_;
    GroupInjConstraints inj_constraints_;
    GasSalesContract gas_sales_;
    GasConsumption gas_consumption_;
    double efficiency_factor_ = 1.0;
    
    // Current rates
    double current_oil_rate_ = 0.0;
    double current_water_rate_ = 0.0;
    double current_gas_rate_ = 0.0;
    double current_resv_rate_ = 0.0;
    double current_inj_rate_ = 0.0;
    double current_inj_resv_rate_ = 0.0;
    
    // Cumulative
    double cumulative_oil_ = 0.0;
    double cumulative_water_ = 0.0;
    double cumulative_gas_ = 0.0;
    double cumulative_injection_ = 0.0;
};

/**
 * @brief Guide rate calculation for allocation
 */
class GuideRateCalculator {
public:
    GuideRateCalculator();
    
    /**
     * @brief Calculate guide rates for wells in a group
     * 
     * @param group Group to calculate for
     * @param well_potentials Map of well name to potential rates (oil, water, gas)
     * @return Map of well name to guide rate fraction
     */
    std::map<std::string, double> calculateGuideRates(
        const WellGroup& group,
        const std::map<std::string, std::array<double, 3>>& well_potentials) const;
    
    /**
     * @brief Allocate group target to wells
     * 
     * @param target Total target rate to allocate
     * @param guide_rates Guide rate fractions for each well
     * @return Map of well name to allocated rate
     */
    std::map<std::string, double> allocateTarget(
        double target,
        const std::map<std::string, double>& guide_rates) const;
    
    // Settings
    void setDampingFactor(double damping) { damping_factor_ = damping; }
    void setMinimumRate(double min_rate) { minimum_rate_ = min_rate; }
    
private:
    double damping_factor_ = 1.0;
    double minimum_rate_ = 0.0;
};

/**
 * @brief Group control manager
 * 
 * Manages the complete well group hierarchy and enforces constraints.
 */
class GroupControlManager {
public:
    GroupControlManager();
    
    // Build hierarchy (GRUPTREE)
    void createGroup(const std::string& name, const std::string& parent = "FIELD");
    void assignWellToGroup(const std::string& well, const std::string& group);
    void removeWellFromGroup(const std::string& well, const std::string& group);
    
    // Set constraints (GCONPROD, GCONINJE)
    void setProductionConstraints(const std::string& group, 
                                   const GroupProdConstraints& constraints);
    void setInjectionConstraints(const std::string& group,
                                  const GroupInjConstraints& constraints);
    
    // Gas sales/consumption (GCONSALE, GCONSUMP)
    void setGasSalesContract(const std::string& group, const GasSalesContract& contract);
    void setGasConsumption(const std::string& group, const GasConsumption& consumption);
    
    // Efficiency factors (GEFAC)
    void setEfficiencyFactor(const std::string& group, double efficiency);
    
    // Access groups
    WellGroup* getGroup(const std::string& name);
    const WellGroup* getGroup(const std::string& name) const;
    WellGroup* getFieldGroup() { return getGroup("FIELD"); }
    
    // Get group containing a well
    std::string getWellGroup(const std::string& well) const;
    
    // Get all groups at a level
    std::vector<WellGroup*> getGroupsAtLevel(int level);
    
    /**
     * @brief Calculate effective constraints for a well
     * 
     * Applies hierarchy: Field -> Region -> Group -> Well
     * Returns the most restrictive constraint at each level.
     * 
     * @param well Well name
     * @return Effective production constraints
     */
    GroupProdConstraints getEffectiveConstraints(const std::string& well) const;
    
    /**
     * @brief Calculate required well rate adjustments
     * 
     * Compares current group rates to targets and calculates
     * how much each well needs to adjust.
     * 
     * @param well_rates Current well rates (oil, water, gas)
     * @return Map of well name to rate adjustment factor
     */
    std::map<std::string, double> calculateRateAdjustments(
        const std::map<std::string, std::array<double, 3>>& well_rates) const;
    
    /**
     * @brief Apply group controls to well targets
     * 
     * Given well potentials, returns allocated rates that satisfy
     * all group constraints.
     * 
     * @param well_potentials Well potential rates
     * @return Allocated rates for each well
     */
    std::map<std::string, std::array<double, 3>> applyGroupControls(
        const std::map<std::string, std::array<double, 3>>& well_potentials);
    
    /**
     * @brief Update group rates from well data
     * 
     * Aggregates well production/injection to group totals.
     * 
     * @param well_rates Current well rates
     */
    void updateGroupRates(const std::map<std::string, std::array<double, 4>>& well_rates);
    
    /**
     * @brief Check if any group constraint is violated
     * 
     * @return List of violated constraints (group name, constraint type)
     */
    std::vector<std::pair<std::string, std::string>> checkConstraintViolations() const;
    
    // Voidage replacement calculation
    double calculateVoidageReplacementRatio(const std::string& group) const;
    double calculateRequiredInjection(const std::string& group, double vrep_target) const;
    
    // Time stepping
    void updateCumulatives(double dt);
    
    // Output
    void printHierarchy() const;
    void printGroupSummary(const std::string& group) const;
    void writeGCON(std::ostream& os) const;
    
    // Parse ECLIPSE keywords
    void parseGRUPTREE(const std::vector<std::string>& lines);
    void parseGCONPROD(const std::vector<std::string>& lines);
    void parseGCONINJE(const std::vector<std::string>& lines);
    void parseGEFAC(const std::vector<std::string>& lines);
    
private:
    std::map<std::string, std::unique_ptr<WellGroup>> groups_;
    std::map<std::string, std::string> well_to_group_;  // Well -> immediate group
    
    GuideRateCalculator guide_calc_;
    
    // Helper functions
    void aggregateToParent(const std::string& group);
    void propagateConstraints();
    int getGroupLevel(const std::string& group) const;
    std::vector<std::string> getGroupPath(const std::string& group) const;
};

/**
 * @brief Network balancing for pressure-constrained systems
 * 
 * Handles scenarios where multiple groups share common pressure constraints
 * (e.g., common separator pressure, pipeline capacity).
 */
class NetworkBalancer {
public:
    NetworkBalancer();
    
    // Set up network
    void addNode(const std::string& name, double pressure);
    void addConnection(const std::string& from, const std::string& to, 
                       double conductivity);
    
    // Connect wells to nodes
    void connectWellToNode(const std::string& well, const std::string& node);
    
    // Pressure constraints
    void setNodePressure(const std::string& node, double pressure);
    void setMinNodePressure(const std::string& node, double min_pressure);
    void setMaxNodePressure(const std::string& node, double max_pressure);
    
    // Flow constraints
    void setConnectionCapacity(const std::string& from, const std::string& to,
                               double max_rate);
    
    /**
     * @brief Balance network flows
     * 
     * Given well potentials and node pressures, calculate achievable
     * well rates that satisfy network constraints.
     * 
     * @param well_potentials Well potential rates at reference pressure
     * @return Achievable well rates
     */
    std::map<std::string, double> balanceNetwork(
        const std::map<std::string, double>& well_potentials);
    
    // Get results
    double getNodePressure(const std::string& node) const;
    double getConnectionFlow(const std::string& from, const std::string& to) const;
    
private:
    struct NetworkNode {
        std::string name;
        double pressure;
        double min_pressure;
        double max_pressure;
        std::set<std::string> connected_wells;
    };
    
    struct NetworkConnection {
        std::string from;
        std::string to;
        double conductivity;
        double max_rate;
        double current_flow;
    };
    
    std::map<std::string, NetworkNode> nodes_;
    std::vector<NetworkConnection> connections_;
    std::map<std::string, std::string> well_nodes_;  // Well -> node
    
    // Newton iteration for network solution
    bool solveNetwork(const std::map<std::string, double>& well_potentials,
                      std::map<std::string, double>& node_pressures,
                      std::map<std::string, double>& well_rates);
};

/**
 * @brief Parse group production control mode from string
 */
GroupProdControlMode parseGroupProdControlMode(const std::string& str);

/**
 * @brief Parse group injection control mode from string
 */
GroupInjControlMode parseGroupInjControlMode(const std::string& str);

/**
 * @brief Parse guide rate type from string
 */
GuideRateType parseGuideRateType(const std::string& str);

/**
 * @brief Parse injection phase from string
 */
InjectionPhase parseInjectionPhase(const std::string& str);

} // namespace FSRM

#endif // GROUP_CONTROL_HPP
