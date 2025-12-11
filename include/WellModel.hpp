#ifndef WELL_MODEL_HPP
#define WELL_MODEL_HPP

#include "FSRM.hpp"
#include "WellboreHydraulics.hpp"
#include <vector>
#include <string>
#include <map>
#include <memory>

namespace FSRM {

struct WellCompletion {
    int i, j, k;  // Grid indices
    double well_index;
    double diameter;
    double skin_factor;
    bool is_open;
    double perforation_length;
};

// Directional well survey data point
struct SurveyPoint {
    double measured_depth;      // Measured depth along wellbore (m)
    double inclination;         // Inclination angle from vertical (degrees)
    double azimuth;             // Azimuth angle from North (degrees)
    double tvd;                 // True vertical depth (m)
    double north_offset;        // North coordinate offset (m)
    double east_offset;         // East coordinate offset (m)
    double dogleg_severity;     // Dogleg severity (degrees/30m)
};

// Directional well trajectory segment
struct TrajectorySegment {
    SurveyPoint start;
    SurveyPoint end;
    double length;              // Segment length (m)
    double avg_inclination;     // Average inclination (degrees)
    double avg_azimuth;         // Average azimuth (degrees)
    std::vector<int> grid_cells; // Grid cells intersected by segment
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
    
    // High-fidelity wellbore hydraulics (THP/BHP conversion and breakdown)
    struct WellboreHydraulicsOptions {
        double measured_depth = 0.0;                 // Total measured depth (m)
        double tubing_id = 0.1;                      // Tubing inner diameter (m)
        double roughness = 1.5e-5;                   // Roughness (m)
        double drag_reduction_factor = 0.0;          // 0..1
        
        // Fluid properties for the wellbore model
        double fluid_density = 1000.0;               // kg/m^3
        double fluid_viscosity = 0.001;              // Pa*s
        
        // Optional power-law rheology (if enabled)
        bool use_power_law = false;
        double k_prime = 0.0;                        // Pa*s^n
        double n_prime = 1.0;                        // -
        
        // Friction correlation (defaults to Haaland)
        FrictionModel::Correlation friction_correlation =
            FrictionModel::Correlation::HAALAND;
    };
    
    void configureWellboreHydraulics(const WellboreHydraulicsOptions& options);
    double computeBHPFromTHP(double thp, double rate) const;
    double computeTHPFromBHP(double bhp_in, double rate) const;
    std::map<std::string, double> computeWellborePressureBreakdown(double thp, double rate) const;
    
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
    
    // Optional high-fidelity wellbore hydraulics engine
    bool use_wellbore_hydraulics_ = false;
    WellboreHydraulicsOptions wellbore_hydraulics_options_;
    mutable std::unique_ptr<WellboreHydraulicsCalculator> wellbore_hydraulics_;
    
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

// Directional/Deviated well model
class DirectionalWell : public WellModel {
public:
    DirectionalWell(const std::string& name, WellType type);
    
    // Survey data management
    void addSurveyPoint(double md, double inc, double azi);
    void addSurveyPoint(const SurveyPoint& point);
    void computeTrajectory();  // Build trajectory from survey points
    
    // Trajectory calculations
    void calculateDoglegSeverity();
    void calculatePositions();  // Compute TVD and offsets from angles
    std::vector<TrajectorySegment> getTrajectorySegments() const;
    
    // Well geometry queries
    double getTotalMeasuredDepth() const;
    double getTrueVerticalDepth() const;
    double getMaxInclination() const;
    SurveyPoint interpolateSurvey(double md) const;
    
    // Grid intersection
    void computeGridIntersections(int nx, int ny, int nz,
                                 double dx, double dy, double dz);
    std::vector<int> getCellsIntersected() const;
    
    // Well index for deviated wells
    double computeWellIndex(int completion_idx,
                           double kx, double ky, double kz,
                           double dx, double dy, double dz) const override;
    
    // Deviated well productivity
    double computeDeviatedWellPI(double length, double angle,
                                 double kh, double kv,
                                 double h, double rw) const;
    
    // Kickoff and build section parameters
    void setKickoffDepth(double depth);
    void setBuildRate(double rate_deg_per_30m);  // Build rate in deg/30m
    void setMaxInclination(double max_inc);
    
    // Advanced trajectory types
    void buildSCurve(double kickoff_md, double build_rate,
                     double target_inc, double hold_md, double drop_rate);
    void buildJCurve(double kickoff_md, double build_rate, double target_inc);
    void buildSlantWell(double kickoff_md, double target_inc, double target_azi);
    
    // Getters
    const std::vector<SurveyPoint>& getSurveyData() const { return survey_data; }
    double getKickoffDepth() const { return kickoff_depth; }
    double getBuildRate() const { return build_rate; }
    
protected:
    std::vector<SurveyPoint> survey_data;
    std::vector<TrajectorySegment> trajectory;
    
    double kickoff_depth;       // Depth where well deviates from vertical
    double build_rate;          // Build rate (degrees per 30m)
    double max_inclination;     // Maximum inclination angle
    
    // Minimum curvature method for trajectory calculation
    void computeMinimumCurvature(const SurveyPoint& p1, const SurveyPoint& p2,
                                SurveyPoint& p2_out);
    
    // Dogleg calculation between two survey points
    double computeDogleg(const SurveyPoint& p1, const SurveyPoint& p2) const;
    
    // Intersected cells with weights
    std::map<int, double> cell_intersection_lengths;
};

} // namespace FSRM

#endif // WELL_MODEL_HPP
