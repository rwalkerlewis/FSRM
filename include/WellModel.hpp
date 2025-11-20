#ifndef WELL_MODEL_HPP
#define WELL_MODEL_HPP

#include "ReservoirSim.hpp"
#include <vector>
#include <string>

namespace ResSim {

struct WellCompletion {
    int i, j, k;  // Grid indices
    double well_index;
    double diameter;
    double skin_factor;
    bool is_open;
    double perforation_length;
};

struct WellConstraints {
    double max_rate;
    double min_rate;
    double max_bhp;
    double min_bhp;
    double max_thp;
    double min_thp;
};

enum class WellControlMode {
    RATE_CONTROL,
    BHP_CONTROL,
    THP_CONTROL,
    RESERVOIR_VOIDAGE
};

class WellModel {
public:
    WellModel(const std::string& name, WellType type);
    virtual ~WellModel() = default;
    
    // Setup
    void addCompletion(const WellCompletion& comp);
    void setControl(WellControlMode mode, double target);
    void setConstraints(const WellConstraints& constraints);
    void setReferenceDepth(double depth);
    
    // Well equations
    virtual double computeWellIndex(int completion_idx, 
                                   double kx, double ky, double kz,
                                   double dx, double dy, double dz) const;
    
    virtual double computeBottomholePressure(double reservoir_pressure,
                                            double rate, int comp_idx) const;
    
    virtual double computeRate(double reservoir_pressure, 
                              double bhp, int comp_idx) const;
    
    // Multi-phase wells
    virtual void computePhaseRates(const std::vector<double>& res_pressures,
                                  const std::vector<double>& saturations,
                                  double bhp,
                                  std::vector<double>& phase_rates) const;
    
    // Contribution to residual and Jacobian
    virtual void contributeToResidual(Vec F, Vec U, DM dm) const;
    virtual void contributeToJacobian(Mat J, Vec U, DM dm) const;
    
    // Wellbore model
    void enableWellboreModel(bool enable);
    void setWellboreGeometry(double inner_radius, double outer_radius);
    void setWellboreRoughness(double roughness);
    
    // Drift flux model for multiphase flow in wellbore
    void computeWellborePressureDrop(const std::vector<double>& phase_rates,
                                    double depth_interval,
                                    double& dp) const;
    
    // Getters
    std::string getName() const { return well_name; }
    WellType getType() const { return well_type; }
    WellControlMode getControlMode() const { return control_mode; }
    double getTargetValue() const { return target_value; }
    const std::vector<WellCompletion>& getCompletions() const { return completions; }
    
    // Performance metrics
    double getTotalRate() const { return total_rate; }
    double getBottomholePressure() const { return bhp; }
    double getCumulativeProduction() const { return cumulative_production; }
    
    void updatePerformance(double dt, double rate);
    
protected:
    std::string well_name;
    WellType well_type;
    WellControlMode control_mode;
    double target_value;
    
    std::vector<WellCompletion> completions;
    WellConstraints constraints;
    
    double reference_depth;
    
    // Wellbore properties
    bool use_wellbore_model;
    double wellbore_inner_radius;
    double wellbore_outer_radius;
    double wellbore_roughness;
    
    // Performance data
    double total_rate;
    double bhp;
    double cumulative_production;
    
    // Peaceman well index calculation
    double computePeacemanRadius(double dx, double dy, 
                                double kx, double ky) const;
};

class ProductionWell : public WellModel {
public:
    ProductionWell(const std::string& name);
    
    void setTargetOilRate(double rate);
    void setTargetGasRate(double rate);
    void setTargetLiquidRate(double rate);
    void setTargetReservoirVoidage(double rate);
    
    void computePhaseRates(const std::vector<double>& res_pressures,
                          const std::vector<double>& saturations,
                          double bhp,
                          std::vector<double>& phase_rates) const override;
};

class InjectionWell : public WellModel {
public:
    InjectionWell(const std::string& name);
    
    void setInjectionFluid(const std::string& fluid_type);  // WATER, GAS, STEAM, etc.
    void setTargetInjectionRate(double rate);
    void setInjectionTemperature(double temp);
    void setInjectionComposition(const std::vector<double>& composition);
    
    std::string getInjectionFluid() const { return injection_fluid; }
    double getInjectionTemperature() const { return injection_temperature; }
    
private:
    std::string injection_fluid;
    double injection_temperature;
    std::vector<double> injection_composition;
};

// Advanced well models
class HorizontalWell : public WellModel {
public:
    HorizontalWell(const std::string& name, WellType type);
    
    void setWellPath(const std::vector<std::vector<double>>& path);
    void enableInflowControlDevices(bool enable);
    
    double computeWellIndex(int completion_idx,
                           double kx, double ky, double kz,
                           double dx, double dy, double dz) const override;
    
private:
    std::vector<std::vector<double>> well_path;
    bool use_icds;
};

class MultilateralWell : public WellModel {
public:
    MultilateralWell(const std::string& name, WellType type);
    
    void addLateral(int lateral_id, const std::vector<WellCompletion>& comps);
    void setLateralControl(int lateral_id, bool is_open);
    
private:
    std::map<int, std::vector<WellCompletion>> laterals;
    std::map<int, bool> lateral_status;
};

} // namespace ResSim

#endif // WELL_MODEL_HPP
