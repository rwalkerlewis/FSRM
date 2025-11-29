#ifndef WELLBORE_HYDRAULICS_HPP
#define WELLBORE_HYDRAULICS_HPP

/**
 * @file WellboreHydraulics.hpp
 * @brief Comprehensive wellbore hydraulics for fracturing operations
 * 
 * ResFrac-equivalent capabilities:
 * - Friction pressure drop (Darcy, Fanning, turbulent)
 * - Perforation friction
 * - Near-wellbore tortuosity
 * - Temperature profile (wellbore cooling)
 * - Multiphase flow (during flowback)
 * - Wellbore storage effects
 * - Surface/downhole pressure conversion
 */

#include <vector>
#include <memory>
#include <string>
#include <map>
#include <array>
#include <functional>
#include <cmath>

namespace FSRM {

/**
 * @brief Wellbore segment properties
 */
struct WellboreSegment {
    double start_md;                    ///< Start measured depth (m)
    double end_md;                      ///< End measured depth (m)
    double start_tvd;                   ///< Start true vertical depth (m)
    double end_tvd;                     ///< End true vertical depth (m)
    double inclination;                 ///< Average inclination (degrees)
    double azimuth;                     ///< Average azimuth (degrees)
    
    // Pipe properties
    double inner_diameter;              ///< Inner diameter (m)
    double outer_diameter;              ///< Outer diameter (m)
    double roughness;                   ///< Surface roughness (m)
    
    // Thermal properties
    double thermal_conductivity;        ///< Pipe thermal conductivity (W/m/K)
    double cement_conductivity;         ///< Cement thermal conductivity (W/m/K)
    double formation_conductivity;      ///< Formation thermal conductivity (W/m/K)
    
    // Calculated properties
    double flow_area;                   ///< Flow cross-sectional area (m²)
    double hydraulic_diameter;          ///< Hydraulic diameter (m)
    double wetted_perimeter;            ///< Wetted perimeter (m)
    
    WellboreSegment();
    
    void calculateDerivedProperties();
};

/**
 * @brief Wellbore completion description
 */
struct WellboreCompletion {
    enum class CompletionType {
        CASED_CEMENTED,                 ///< Cased and cemented
        CASED_UNCEMENTED,              ///< Cased, uncemented (liner)
        OPENHOLE,                       ///< Openhole completion
        SLOTTED_LINER,                  ///< Slotted/preperforated liner
        GRAVEL_PACK                     ///< Gravel pack completion
    };
    
    CompletionType type;
    double casing_id;                   ///< Casing inner diameter (m)
    double casing_od;                   ///< Casing outer diameter (m)
    double tubing_id;                   ///< Tubing inner diameter (m)
    double tubing_od;                   ///< Tubing outer diameter (m)
    double annulus_area;                ///< Casing-tubing annulus area (m²)
    
    bool has_tubing;                    ///< Is tubing present
    double tubing_setting_depth;        ///< Tubing end MD (m)
    
    WellboreCompletion();
};

/**
 * @brief Friction pressure drop correlations
 */
class FrictionModel {
public:
    FrictionModel();
    
    /**
     * @brief Calculate Fanning friction factor
     * 
     * @param reynolds Reynolds number
     * @param relative_roughness ε/D
     * @return Fanning friction factor
     */
    double fanningFrictionFactor(double reynolds, 
                                  double relative_roughness) const;
    
    /**
     * @brief Calculate friction pressure gradient for Newtonian fluid
     * 
     * @param velocity Flow velocity (m/s)
     * @param density Fluid density (kg/m³)
     * @param viscosity Dynamic viscosity (Pa·s)
     * @param diameter Pipe diameter (m)
     * @param roughness Surface roughness (m)
     * @return Pressure gradient (Pa/m)
     */
    double newtonianFrictionGradient(double velocity, double density,
                                      double viscosity, double diameter,
                                      double roughness) const;
    
    /**
     * @brief Calculate friction for power-law fluid
     * 
     * @param velocity Flow velocity (m/s)
     * @param density Fluid density (kg/m³)
     * @param k_prime Consistency index K' (Pa·s^n)
     * @param n_prime Flow behavior index n'
     * @param diameter Pipe diameter (m)
     * @param roughness Surface roughness (m)
     * @return Pressure gradient (Pa/m)
     */
    double powerLawFrictionGradient(double velocity, double density,
                                     double k_prime, double n_prime,
                                     double diameter, double roughness) const;
    
    /**
     * @brief Calculate friction with drag reduction
     * 
     * For slickwater with friction reducer.
     * 
     * @param base_gradient Base friction gradient (Pa/m)
     * @param drag_reduction_factor Drag reduction (0-1)
     * @return Reduced friction gradient (Pa/m)
     */
    double applyDragReduction(double base_gradient,
                               double drag_reduction_factor) const;
    
    // Correlation selection
    enum class Correlation {
        COLEBROOK_WHITE,                ///< Implicit Colebrook-White
        HAALAND,                        ///< Explicit Haaland approximation
        SWAMEE_JAIN,                    ///< Swamee-Jain explicit
        CHEN                            ///< Chen explicit
    };
    
    void setCorrelation(Correlation corr) { correlation_ = corr; }
    
private:
    Correlation correlation_;
    
    // Colebrook-White iterative solution
    double colebrookWhite(double reynolds, double rel_roughness) const;
};

/**
 * @brief Perforation friction calculator
 */
class PerforationFriction {
public:
    PerforationFriction();
    
    /**
     * @brief Calculate perforation pressure drop
     * 
     * Using orifice equation: ΔP = ρQ²/(2Cd²A²)
     * 
     * @param flow_rate Total flow rate through perfs (m³/s)
     * @param num_perforations Number of perforations
     * @param perf_diameter Perforation diameter (m)
     * @param discharge_coefficient Cd (typically 0.7-0.9)
     * @param fluid_density Fluid density (kg/m³)
     * @return Perforation pressure drop (Pa)
     */
    double calculatePressureDrop(double flow_rate, int num_perforations,
                                  double perf_diameter,
                                  double discharge_coefficient,
                                  double fluid_density) const;
    
    /**
     * @brief Calculate perforation erosion
     * 
     * @param flow_rate Flow rate (m³/s)
     * @param num_perforations Number of perfs
     * @param current_diameter Current perf diameter (m)
     * @param proppant_concentration Proppant concentration (kg/m³)
     * @param proppant_hardness Relative hardness
     * @param dt Time step (s)
     * @return New perforation diameter (m)
     */
    double calculateErosion(double flow_rate, int num_perforations,
                            double current_diameter,
                            double proppant_concentration,
                            double proppant_hardness,
                            double dt) const;
    
    /**
     * @brief Calculate discharge coefficient
     * 
     * Cd varies with Reynolds number and entrance conditions.
     * 
     * @param reynolds Reynolds number through perforation
     * @param l_over_d Length/diameter ratio
     * @return Discharge coefficient
     */
    double calculateDischargeCoefficient(double reynolds, 
                                          double l_over_d) const;
    
    // Settings
    void setErosionCoefficient(double c) { erosion_coeff_ = c; }
    void setMinDischargeCoeff(double cd) { min_cd_ = cd; }
    
private:
    double erosion_coeff_;
    double min_cd_;
};

/**
 * @brief Near-wellbore friction (tortuosity)
 */
class NearWellboreFriction {
public:
    NearWellboreFriction();
    
    /**
     * @brief Calculate near-wellbore pressure drop
     * 
     * Accounts for tortuous path from perforations to fracture.
     * 
     * @param flow_rate Flow rate (m³/s)
     * @param num_clusters Number of clusters
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @param perforation_phasing Perforation phasing (degrees)
     * @param stress_anisotropy SHmax - Shmin (Pa)
     * @return Near-wellbore friction (Pa)
     */
    double calculateTortuosityPressureDrop(double flow_rate,
                                            int num_clusters,
                                            double fluid_viscosity,
                                            double perforation_phasing,
                                            double stress_anisotropy) const;
    
    /**
     * @brief Calculate tortuosity factor
     * 
     * @param perforation_phasing Degrees from preferred fracture plane
     * @param stress_anisotropy Stress difference (Pa)
     * @return Tortuosity factor (1.0 = no tortuosity)
     */
    double getTortuosityFactor(double perforation_phasing,
                                double stress_anisotropy) const;
    
    // Model parameters
    void setTortuosityCoefficient(double c) { tortuosity_coeff_ = c; }
    
private:
    double tortuosity_coeff_;
};

/**
 * @brief Wellbore temperature model
 */
class WellboreTemperature {
public:
    WellboreTemperature();
    
    /**
     * @brief Initialize temperature profile
     * 
     * @param segments Wellbore segments
     * @param surface_temp Surface temperature (K)
     * @param geothermal_gradient Geothermal gradient (K/m)
     */
    void initialize(const std::vector<WellboreSegment>& segments,
                    double surface_temp,
                    double geothermal_gradient);
    
    /**
     * @brief Calculate temperature profile during injection
     * 
     * @param flow_rate Injection rate (m³/s)
     * @param injection_temp Injection temperature (K)
     * @param fluid_density Fluid density (kg/m³)
     * @param fluid_cp Fluid specific heat (J/kg/K)
     * @param time Injection time (s)
     * @return Temperature at each segment (K)
     */
    std::vector<double> calculateInjectionProfile(
        double flow_rate,
        double injection_temp,
        double fluid_density,
        double fluid_cp,
        double time);
    
    /**
     * @brief Get temperature at depth
     * 
     * @param md Measured depth (m)
     * @return Temperature (K)
     */
    double getTemperatureAtDepth(double md) const;
    
    /**
     * @brief Get bottomhole temperature
     * 
     * @return BHT (K)
     */
    double getBottomholeTemperature() const;
    
    /**
     * @brief Calculate heat loss to formation
     * 
     * @param segment Wellbore segment
     * @param fluid_temp Fluid temperature (K)
     * @param flow_rate Flow rate (m³/s)
     * @return Heat loss rate (W)
     */
    double calculateHeatLoss(const WellboreSegment& segment,
                              double fluid_temp,
                              double flow_rate) const;
    
    // Transient model
    void setTransientModel(bool use_transient) { use_transient_ = use_transient; }
    void setFormationDiffusivity(double alpha) { formation_diffusivity_ = alpha; }
    
private:
    std::vector<WellboreSegment> segments_;
    std::vector<double> formation_temp_;
    std::vector<double> current_temp_;
    
    double surface_temp_;
    double geothermal_gradient_;
    
    bool use_transient_;
    double formation_diffusivity_;
    double injection_time_;
    
    // Ramey's method for wellbore heat transfer
    double rameyTimeFunction(double time, double radius) const;
};

/**
 * @brief Complete wellbore hydraulics calculator
 */
class WellboreHydraulicsCalculator {
public:
    WellboreHydraulicsCalculator();
    
    /**
     * @brief Build wellbore model from survey
     * 
     * @param md Measured depths (m)
     * @param inc Inclinations (degrees)
     * @param azi Azimuths (degrees)
     * @param diameters Pipe inner diameters at each depth (m)
     * @param roughness Surface roughness (m)
     */
    void buildFromSurvey(const std::vector<double>& md,
                         const std::vector<double>& inc,
                         const std::vector<double>& azi,
                         const std::vector<double>& diameters,
                         double roughness);
    
    /**
     * @brief Set completion details
     */
    void setCompletion(const WellboreCompletion& completion);
    
    /**
     * @brief Add perforation clusters
     * 
     * @param cluster_mds Measured depths of clusters
     * @param num_perfs Number of perforations per cluster
     * @param perf_diameters Perforation diameters (m)
     */
    void addPerforations(const std::vector<double>& cluster_mds,
                         const std::vector<int>& num_perfs,
                         const std::vector<double>& perf_diameters);
    
    /**
     * @brief Calculate pressure profile during pumping
     * 
     * @param surface_pressure Surface treating pressure (Pa)
     * @param flow_rate Pump rate (m³/s)
     * @param fluid_density Fluid density (kg/m³)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @param k_prime Power-law K' (Pa·s^n) - optional
     * @param n_prime Power-law n' - optional
     * @return Pressure at each segment depth (Pa)
     */
    std::vector<double> calculatePressureProfile(
        double surface_pressure,
        double flow_rate,
        double fluid_density,
        double fluid_viscosity,
        double k_prime = 0.0,
        double n_prime = 1.0);
    
    /**
     * @brief Calculate surface pressure from BHP
     * 
     * @param bhp Bottomhole pressure (Pa)
     * @param flow_rate Flow rate (m³/s)
     * @param fluid_density Fluid density (kg/m³)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @return Required surface pressure (Pa)
     */
    double calculateSurfacePressure(double bhp, double flow_rate,
                                     double fluid_density,
                                     double fluid_viscosity);
    
    /**
     * @brief Calculate BHP from surface pressure
     * 
     * @param surface_pressure Surface pressure (Pa)
     * @param flow_rate Flow rate (m³/s)
     * @param fluid_density Fluid density (kg/m³)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @return Bottomhole pressure (Pa)
     */
    double calculateBHP(double surface_pressure, double flow_rate,
                        double fluid_density, double fluid_viscosity);
    
    /**
     * @brief Get component pressure drops
     * 
     * @return Map of component name to pressure drop (Pa)
     */
    std::map<std::string, double> getPressureBreakdown() const;
    
    /**
     * @brief Calculate hydrostatic pressure
     * 
     * @param md1 Start depth (m)
     * @param md2 End depth (m)
     * @param fluid_density Fluid density (kg/m³)
     * @return Hydrostatic pressure difference (Pa)
     */
    double calculateHydrostatic(double md1, double md2,
                                 double fluid_density) const;
    
    /**
     * @brief Calculate friction pressure drop
     * 
     * @param md1 Start depth (m)
     * @param md2 End depth (m)
     * @param flow_rate Flow rate (m³/s)
     * @param fluid_density Fluid density (kg/m³)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @return Friction pressure drop (Pa)
     */
    double calculateFriction(double md1, double md2,
                              double flow_rate,
                              double fluid_density,
                              double fluid_viscosity) const;
    
    /**
     * @brief Update perforation state (erosion)
     * 
     * @param flow_rate Flow rate (m³/s)
     * @param proppant_concentration Proppant concentration (kg/m³)
     * @param dt Time step (s)
     */
    void updatePerforations(double flow_rate,
                            double proppant_concentration,
                            double dt);
    
    // Temperature model
    void enableTemperatureModel(double surface_temp, double geothermal_grad);
    std::vector<double> getTemperatureProfile() const;
    
    // Access
    const std::vector<WellboreSegment>& getSegments() const { return segments_; }
    double getTotalMeasuredDepth() const;
    double getTotalVerticalDepth() const;
    
    // Drag reduction
    void setDragReduction(double dr) { drag_reduction_ = dr; }
    
private:
    std::vector<WellboreSegment> segments_;
    WellboreCompletion completion_;
    
    // Perforations
    struct PerfCluster {
        double md;
        int num_perfs;
        double diameter;
        double discharge_coeff;
    };
    std::vector<PerfCluster> perforations_;
    
    // Models
    FrictionModel friction_model_;
    PerforationFriction perf_friction_;
    NearWellboreFriction nwb_friction_;
    std::unique_ptr<WellboreTemperature> temp_model_;
    
    // Settings
    double drag_reduction_;
    bool use_temperature_model_;
    
    // Results cache
    mutable std::map<std::string, double> pressure_breakdown_;
    
    // Helpers
    double getTVDAtMD(double md) const;
    int getSegmentAtMD(double md) const;
};

/**
 * @brief Multiphase wellbore flow (for production/flowback)
 */
class MultiphaseWellboreFlow {
public:
    MultiphaseWellboreFlow();
    
    /**
     * @brief Calculate pressure gradient for multiphase flow
     * 
     * Uses Beggs-Brill or similar correlation.
     * 
     * @param oil_rate Oil rate (m³/s)
     * @param water_rate Water rate (m³/s)
     * @param gas_rate Gas rate (m³/s at standard conditions)
     * @param pressure Local pressure (Pa)
     * @param temperature Local temperature (K)
     * @param diameter Pipe diameter (m)
     * @param inclination Pipe inclination (degrees from vertical)
     * @return Pressure gradient (Pa/m)
     */
    double calculatePressureGradient(double oil_rate, double water_rate,
                                      double gas_rate, double pressure,
                                      double temperature, double diameter,
                                      double inclination);
    
    /**
     * @brief Calculate liquid holdup
     * 
     * @param superficial_liquid_velocity (m/s)
     * @param superficial_gas_velocity (m/s)
     * @param liquid_density (kg/m³)
     * @param gas_density (kg/m³)
     * @param inclination (degrees)
     * @return Liquid holdup (0-1)
     */
    double calculateHoldup(double vsl, double vsg,
                           double rho_l, double rho_g,
                           double inclination);
    
    /**
     * @brief Identify flow regime
     */
    enum class FlowRegime {
        SEGREGATED,                     ///< Stratified/slug
        INTERMITTENT,                   ///< Slug/plug
        DISTRIBUTED,                    ///< Bubble/mist
        TRANSITION                      ///< Between regimes
    };
    
    FlowRegime identifyFlowRegime(double vsl, double vsg,
                                   double rho_l, double rho_g,
                                   double diameter, double inclination);
    
    // Correlation selection
    enum class Correlation {
        BEGGS_BRILL,                    ///< Beggs and Brill (1973)
        HAGEDORN_BROWN,                 ///< Hagedorn and Brown
        DUNS_ROS,                       ///< Duns and Ros
        MUKHERJEE_BRILL,                ///< Mukherjee and Brill
        GRAY                            ///< Gray correlation
    };
    
    void setCorrelation(Correlation corr) { correlation_ = corr; }
    
private:
    Correlation correlation_;
    
    // Beggs-Brill parameters
    double beggsbrillHoldup(double vsl, double vsg, double diameter,
                            double rho_l, double rho_g, double mu_l,
                            double sigma, double inclination);
};

} // namespace FSRM

#endif // WELLBORE_HYDRAULICS_HPP
