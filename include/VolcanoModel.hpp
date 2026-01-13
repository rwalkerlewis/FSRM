/**
 * @file VolcanoModel.hpp
 * @brief Fully coupled volcano modeling physics
 *
 * Implements comprehensive multi-physics models for volcanic systems:
 * - Magma chamber dynamics (pressurization, crystallization, volatile exsolution)
 * - Conduit flow (magma ascent, fragmentation, bubble nucleation/growth)
 * - Volcanic deformation (Mogi source, finite element ground deformation)
 * - Eruption dynamics (effusive lava flows, explosive plumes)
 * - Volcanic seismicity (VT earthquakes, LP events, volcanic tremor)
 * - Pyroclastic density currents (PDCs)
 * - Lahar/volcanic debris flows
 * - Volcanic gas emissions (SO2, CO2, H2O)
 * - Tephra/ash dispersal
 * - Coupling with infrasound for monitoring
 *
 * Supports realistic volcanic scenarios:
 * - Stratovolcanoes (Mt. St. Helens, Pinatubo, Merapi)
 * - Shield volcanoes (Kilauea, Mauna Loa)
 * - Calderas (Long Valley, Yellowstone, Campi Flegrei)
 * - Monogenetic fields (Auckland Volcanic Field)
 *
 * Scientific References:
 * - Sparks, R.S.J. (1978) - Magma chamber dynamics
 * - Woods (1988) - Eruption column modeling
 * - Mastin et al. (2009) - Ash dispersal (Ash3d)
 * - Macedonio et al. (2005) - Pyroclastic density currents
 * - Mogi (1958) - Volcanic deformation source
 * - Chouet (1996) - Volcanic seismicity
 *
 * @author FSRM Development Team
 */

#ifndef VOLCANO_MODEL_HPP
#define VOLCANO_MODEL_HPP

#include "FSRM.hpp"
#include "FluidModel.hpp"
#include "ThermalPressurization.hpp"
#include "AtmosphericInfrasound.hpp"
#include "SeismicSource.hpp"
#include <vector>
#include <array>
#include <string>
#include <map>
#include <memory>
#include <functional>
#include <cmath>
#include <complex>

namespace FSRM {

// Forward declarations
class PhysicsKernel;
class DiscontinuousGalerkin;

// =============================================================================
// Physical Constants for Volcanic Systems
// =============================================================================
namespace VolcanoConstants {
    constexpr double GRAVITY = 9.81;                    // m/s²
    constexpr double GAS_CONSTANT = 8.314;              // J/(mol·K)
    constexpr double STEFAN_BOLTZMANN = 5.67e-8;        // W/(m²·K⁴)
    constexpr double ATMOSPHERIC_PRESSURE = 101325.0;   // Pa
    
    // Magma properties
    constexpr double SILICATE_MELT_DENSITY = 2400.0;    // kg/m³ (basalt melt)
    constexpr double RHYOLITE_DENSITY = 2300.0;         // kg/m³
    constexpr double CRYSTAL_DENSITY = 2700.0;          // kg/m³
    constexpr double MAGMA_SPECIFIC_HEAT = 1200.0;      // J/(kg·K)
    constexpr double LATENT_HEAT_CRYSTALLIZATION = 4e5; // J/kg
    
    // Volatile properties
    constexpr double H2O_MOLECULAR_WEIGHT = 0.018;      // kg/mol
    constexpr double CO2_MOLECULAR_WEIGHT = 0.044;      // kg/mol
    constexpr double SO2_MOLECULAR_WEIGHT = 0.064;      // kg/mol
    
    // Earth properties
    constexpr double CRUSTAL_DENSITY = 2700.0;          // kg/m³
    constexpr double MANTLE_DENSITY = 3300.0;           // kg/m³
    constexpr double SHEAR_MODULUS_CRUST = 30e9;        // Pa
}

// =============================================================================
// Enumerations
// =============================================================================

/**
 * @brief Types of volcanic eruptions
 */
enum class EruptionType {
    EFFUSIVE,           // Lava flow dominated
    HAWAIIAN,           // Fire fountaining
    STROMBOLIAN,        // Regular mild explosions
    VULCANIAN,          // Short violent explosions
    SUBPLINIAN,         // Moderate sustained columns
    PLINIAN,            // Large sustained columns
    ULTRAPLINIAN,       // Extreme columns (>40 km)
    SURTSEYAN,          // Phreatomagmatic (water interaction)
    PHREATOPLINIAN,     // Large phreatomagmatic
    DOME_COLLAPSE,      // Dome failure -> PDC
    LATERAL_BLAST       // Directed blast (Mt. St. Helens 1980)
};

/**
 * @brief Types of magma composition
 */
enum class MagmaComposition {
    BASALT,             // Low silica (~50% SiO2), low viscosity
    BASALTIC_ANDESITE,  // ~53% SiO2
    ANDESITE,           // ~58% SiO2
    DACITE,             // ~65% SiO2
    RHYOLITE,           // High silica (~75% SiO2), high viscosity
    PHONOLITE,          // Alkaline, intermediate
    TRACHYTE,           // Alkaline, high silica
    CARBONATITE         // Carbonate magma (rare)
};

/**
 * @brief Types of volcanic seismicity
 */
enum class VolcanicSeismicEventType {
    VT,                 // Volcano-tectonic (brittle fracture)
    LP,                 // Long-period (fluid-filled cracks)
    VLP,                // Very long-period (conduit resonance)
    TREMOR,             // Continuous tremor
    HYBRID,             // Mixed VT + LP characteristics
    EXPLOSION_QUAKE,    // Source = explosion
    ROCKFALL,           // Dome collapse/rockfall
    LAHAR               // Debris flow signal
};

/**
 * @brief Volcanic deformation source types
 */
enum class DeformationSourceType {
    MOGI,               // Point spherical source (Mogi, 1958)
    MCTIGUE,            // Finite spherical cavity (McTigue, 1987)
    PENNY_CRACK,        // Horizontal penny-shaped crack
    PROLATE_SPHEROID,   // Prolate (cigar-shaped) source
    OBLATE_SPHEROID,    // Oblate (pancake-shaped) source
    RECTANGULAR_DIKE,   // Vertical dike (Okada)
    HORIZONTAL_SILL,    // Horizontal sill intrusion
    FINITE_ELEMENT      // Full FEM solution
};

/**
 * @brief Conduit flow regimes
 */
enum class ConduitFlowRegime {
    BUBBLY,             // Dispersed bubbles in melt
    SLUG,               // Large gas slugs
    CHURN,              // Transitional chaotic
    ANNULAR,            // Gas core, liquid annulus
    DISPERSED,          // Post-fragmentation (particles in gas)
    DENSE_SUSPENSION    // High particle concentration
};

// =============================================================================
// Magma Properties and Rheology
// =============================================================================

/**
 * @brief Magma composition and volatile content
 */
struct MagmaProperties {
    MagmaComposition composition = MagmaComposition::ANDESITE;
    
    // Major oxide composition (wt%)
    double SiO2 = 60.0;
    double Al2O3 = 17.0;
    double FeO = 6.0;
    double MgO = 3.0;
    double CaO = 6.0;
    double Na2O = 4.0;
    double K2O = 2.0;
    
    // Volatile content (wt%)
    double H2O_total = 4.0;         // Total water content
    double CO2_total = 0.1;         // Total CO2 content
    double SO2_total = 0.05;        // Total SO2 content
    double Cl_total = 0.1;          // Total chlorine
    
    // Physical properties
    double temperature = 1273.15;   // K (~1000°C)
    double crystal_fraction = 0.2; // Volume fraction
    double bubble_fraction = 0.0;   // Volume fraction (calculated)
    
    // Get melt viscosity using Giordano et al. (2008) model
    double getMeltViscosity() const;
    
    // Get bulk magma viscosity including crystals (Costa, 2005)
    double getBulkViscosity() const;
    
    // Get density accounting for crystals and bubbles
    double getDensity(double pressure) const;
    
    // Get dissolved volatile content at given pressure
    double getDissolvedH2O(double pressure) const;  // wt%
    double getDissolvedCO2(double pressure) const;  // wt%
    
    // Get saturation pressure for given volatile content
    double getSaturationPressure() const;
    
    // Get liquidus and solidus temperatures
    double getLiquidusTemperature() const;
    double getSolidusTemperature() const;
    
    // Crystallization rate (fraction per degree cooling)
    double getCrystallizationRate() const;
};

/**
 * @brief Magma rheology model
 * 
 * Implements temperature, crystal, and bubble-dependent viscosity
 * following the Herschel-Bulkley model for crystal-bearing magmas.
 */
class MagmaRheology {
public:
    MagmaRheology();
    
    void setMagmaProperties(const MagmaProperties& props);
    
    /**
     * @brief Get effective viscosity at given shear rate
     * 
     * For crystal-rich magmas, uses Herschel-Bulkley model:
     * τ = τ_y + K * γ̇^n
     */
    double getEffectiveViscosity(double shear_rate) const;
    
    /**
     * @brief Get yield stress (for crystal-bearing magmas)
     */
    double getYieldStress() const;
    
    /**
     * @brief Check if magma can flow (τ > τ_y)
     */
    bool canFlow(double shear_stress) const;
    
    /**
     * @brief Get relative viscosity factor from crystals
     */
    double getCrystalViscosityFactor() const;
    
    /**
     * @brief Get relative viscosity factor from bubbles
     */
    double getBubbleViscosityFactor() const;
    
private:
    MagmaProperties magma;
    double melt_viscosity;
    double yield_stress;
    double consistency_K;
    double power_law_n;
    
    void computeRheologicalParameters();
};

// =============================================================================
// Magma Chamber Model
// =============================================================================

/**
 * @brief Magma chamber geometry and state
 */
struct MagmaChamberGeometry {
    // Position (center)
    double x_center = 0.0;          // m
    double y_center = 0.0;          // m
    double depth = 5000.0;          // m below surface
    
    // Dimensions
    double volume = 10e9;           // m³ (10 km³)
    double semi_axis_a = 2000.0;    // m (horizontal)
    double semi_axis_b = 2000.0;    // m (horizontal)
    double semi_axis_c = 1000.0;    // m (vertical)
    
    // Shape type for deformation modeling
    DeformationSourceType shape = DeformationSourceType::OBLATE_SPHEROID;
    
    // Compute volume from axes
    double computeVolume() const {
        return (4.0/3.0) * M_PI * semi_axis_a * semi_axis_b * semi_axis_c;
    }
    
    // Effective radius for Mogi-type calculations
    double effectiveRadius() const {
        return std::pow(3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    }
};

/**
 * @brief Magma chamber thermodynamic state
 */
struct MagmaChamberState {
    double pressure = 200e6;        // Pa (200 MPa)
    double temperature = 1173.15;   // K (900°C)
    double volume = 10e9;           // m³
    double mass = 2.4e13;           // kg
    
    // Volatile content
    double dissolved_H2O = 3.5;     // wt%
    double dissolved_CO2 = 0.05;    // wt%
    double exsolved_gas_fraction = 0.05;  // Volume fraction
    
    // Crystal content
    double crystal_fraction = 0.3;  // Volume fraction
    
    // Melt fraction
    double melt_fraction() const { 
        return 1.0 - crystal_fraction - exsolved_gas_fraction; 
    }
    
    // Overpressure relative to lithostatic
    double overpressure = 10e6;     // Pa (10 MPa)
};

/**
 * @brief Magma chamber dynamics model
 * 
 * Models the evolution of a magma chamber including:
 * - Magma recharge from depth
 * - Cooling and crystallization
 * - Volatile exsolution (second boiling)
 * - Pressure evolution
 * - Chamber failure criteria
 */
class MagmaChamberModel {
public:
    MagmaChamberModel();
    
    void initialize(const MagmaChamberGeometry& geometry,
                   const MagmaProperties& magma,
                   const MagmaChamberState& initial_state);
    
    /**
     * @brief Update chamber state over time
     * 
     * Solves coupled ODEs for:
     * - dP/dt = f(recharge, crystallization, volatile exsolution, eruption)
     * - dT/dt = f(cooling, latent heat, recharge)
     * - dφ_x/dt = f(cooling rate, crystallization kinetics)
     */
    void update(double dt);
    
    // Recharge and eruption
    void setRechargeRate(double rate_kg_s);    // kg/s
    void setEruptionRate(double rate_kg_s);    // kg/s
    void setWallRockCooling(double heat_flux); // W/m²
    
    // Get current state
    const MagmaChamberState& getState() const { return state; }
    double getPressure() const { return state.pressure; }
    double getTemperature() const { return state.temperature; }
    double getOverpressure() const { return state.overpressure; }
    
    // Failure criteria
    bool isNearFailure(double tensile_strength = 10e6) const;
    double getFailureProbability() const;
    
    // Compute surface deformation from chamber pressurization
    void computeSurfaceDeformation(double x, double y,
                                  double& u_radial, double& u_vertical) const;
    
    // Get seismic moment rate from volume change
    double getSeismicMomentRate() const;
    
    // Magma compressibility
    double getMagmaCompressibility() const;
    
private:
    MagmaChamberGeometry geometry;
    MagmaProperties magma;
    MagmaChamberState state;
    
    double recharge_rate;
    double eruption_rate;
    double wall_heat_flux;
    
    // Equation of state
    double computePressureChange(double dV, double dm, double dT);
    double computeVolatileExsolution(double dP, double dT);
    double computeCrystallization(double dT, double dt);
    
    // Deformation Green's functions
    void mogiSource(double x, double y, double depth, double dV,
                   double& ur, double& uz) const;
    void mctigueSource(double x, double y, double depth, double radius,
                      double dP, double& ur, double& uz) const;
};

// =============================================================================
// Volcanic Conduit Flow Model
// =============================================================================

/**
 * @brief Conduit geometry
 */
struct ConduitGeometry {
    double depth = 5000.0;          // m (depth to chamber)
    double length = 5000.0;         // m
    double radius_base = 15.0;      // m (at chamber)
    double radius_vent = 10.0;      // m (at surface)
    
    // Optional: geometry variations
    bool has_constriction = false;
    double constriction_depth = 2000.0;
    double constriction_radius = 5.0;
    
    // Get radius at depth z (0 = vent, length = chamber)
    double radiusAt(double z) const {
        if (has_constriction) {
            // More complex geometry
            if (std::abs(z - constriction_depth) < 100.0) {
                return constriction_radius;
            }
        }
        // Linear taper
        return radius_vent + (radius_base - radius_vent) * z / length;
    }
    
    // Cross-sectional area
    double areaAt(double z) const {
        double r = radiusAt(z);
        return M_PI * r * r;
    }
};

/**
 * @brief Conduit flow state at a point
 */
struct ConduitFlowState {
    double z = 0.0;                 // Position along conduit (m from vent)
    double pressure = 1e5;          // Pa
    double temperature = 1273.0;    // K
    double velocity = 10.0;         // m/s
    double density = 2400.0;        // kg/m³
    
    // Phase fractions
    double melt_fraction = 0.7;
    double crystal_fraction = 0.2;
    double gas_fraction = 0.1;
    double bubble_number_density = 1e12;  // m⁻³
    double mean_bubble_radius = 1e-4;     // m
    
    // Dissolved volatiles
    double dissolved_H2O = 2.0;     // wt%
    double dissolved_CO2 = 0.01;    // wt%
    
    // Flow regime
    ConduitFlowRegime regime = ConduitFlowRegime::BUBBLY;
    
    // Is flow fragmented (gas-particle mixture)?
    bool isFragmented() const { return gas_fraction > 0.75; }
};

/**
 * @brief Steady-state conduit flow model
 * 
 * Solves coupled equations for magma ascent:
 * - Mass conservation: d(ρuA)/dz = 0
 * - Momentum: ρu du/dz = -dP/dz - ρg - f_wall
 * - Energy: dH/dz = ...
 * - Volatile exsolution and bubble growth
 * - Fragmentation criteria
 * 
 * Based on Mastin & Ghiorso (2000), Papale (2001)
 */
class ConduitFlowModel {
public:
    ConduitFlowModel();
    
    void initialize(const ConduitGeometry& geometry,
                   const MagmaProperties& magma,
                   double chamber_pressure);
    
    /**
     * @brief Solve for steady-state flow profile
     * 
     * Integrates from chamber to surface, determining:
     * - Fragmentation depth
     * - Exit velocity and mass flux
     * - Volatile exsolution profile
     */
    bool solvesteadyState();
    
    /**
     * @brief Solve time-dependent conduit flow
     */
    void update(double dt);
    
    // Boundary conditions
    void setChamberPressure(double P);
    void setVentPressure(double P);      // Usually atmospheric
    void setMassFluxTarget(double mdot); // kg/s
    
    // Get solution
    const std::vector<ConduitFlowState>& getFlowProfile() const { return profile; }
    double getMassFlux() const { return mass_flux; }
    double getExitVelocity() const;
    double getFragmentationDepth() const { return fragmentation_depth; }
    
    // Eruption intensity metrics
    double getDischargeRate() const { return mass_flux; }  // kg/s
    double getVolumeFlux() const;                          // m³/s (DRE)
    double getColumnHeight() const;                        // Estimated plume height
    
    // Check stability
    bool isChokedFlow() const;
    double getMaxMassFlux() const;
    
private:
    ConduitGeometry geometry;
    MagmaProperties magma;
    
    std::vector<ConduitFlowState> profile;
    double mass_flux;
    double fragmentation_depth;
    
    int n_nodes;
    double dz;
    
    // Physical models
    double computeViscosity(const ConduitFlowState& state) const;
    double computeWallFriction(const ConduitFlowState& state) const;
    double computeBubbleGrowthRate(const ConduitFlowState& state) const;
    bool checkFragmentation(const ConduitFlowState& state) const;
    
    // Volatile solubility
    double H2O_solubility(double P, double T) const;
    double CO2_solubility(double P, double T) const;
    
    // Equation of state for gas phase
    double gasPhase_density(double P, double T, double X_H2O, double X_CO2) const;
};

// =============================================================================
// Eruption Column and Plume Model
// =============================================================================

/**
 * @brief Eruption column parameters
 */
struct EruptionColumnParameters {
    // Vent conditions (from conduit model)
    double vent_radius = 50.0;          // m
    double exit_velocity = 200.0;       // m/s
    double exit_temperature = 1273.0;   // K
    double mass_flux = 1e7;             // kg/s
    double gas_fraction = 0.95;         // Volume fraction at vent
    double particle_size = 1e-3;        // m (mean diameter)
    
    // Atmosphere
    double atmospheric_temperature = 288.0;  // K (surface)
    double tropopause_height = 11000.0;      // m
    double stratosphere_temperature = 220.0; // K
    
    // Volcanic gas composition (mass fractions in gas phase)
    double X_H2O = 0.85;
    double X_CO2 = 0.10;
    double X_SO2 = 0.03;
    double X_other = 0.02;
};

/**
 * @brief Eruption column state
 */
struct ColumnState {
    double height = 0.0;            // m above vent
    double radius = 0.0;            // m
    double velocity = 0.0;          // m/s (vertical)
    double temperature = 0.0;       // K
    double density = 0.0;           // kg/m³ (bulk)
    double gas_fraction = 0.0;      // Volume fraction
    double particle_fraction = 0.0; // Volume fraction
    double entrained_air = 0.0;     // Mass fraction
    
    // Is the column buoyant?
    bool isBuoyant(double ambient_density) const {
        return density < ambient_density;
    }
};

/**
 * @brief Eruption column model (Woods, 1988; Bursik, 2001)
 * 
 * Models convective eruption columns including:
 * - Jet thrust region (momentum-driven)
 * - Buoyancy-driven rise with air entrainment
 * - Column collapse criteria
 * - Neutral buoyancy level (umbrella cloud)
 * 
 * Solves coupled ODEs:
 * - d(ρuA)/dz = ρ_a * u_e * 2πr  (mass entrainment)
 * - d(ρu²A)/dz = (ρ_a - ρ)gA    (momentum)
 * - d(ρuHaA)/dz = ρ_a * u_e * 2πr * H_a  (enthalpy)
 */
class EruptionColumnModel {
public:
    EruptionColumnModel();
    
    void initialize(const EruptionColumnParameters& params);
    
    /**
     * @brief Solve for column profile
     * 
     * @return True if column reaches neutral buoyancy, false if collapse
     */
    bool solve();
    
    // Get results
    const std::vector<ColumnState>& getProfile() const { return profile; }
    double getMaxHeight() const { return max_height; }
    double getNeutralBuoyancyHeight() const { return NBL_height; }
    double getUmbrellaCloudRadius() const;
    
    // Column collapse
    bool doesCollapse() const { return column_collapsed; }
    double getCollapseHeight() const { return collapse_height; }
    double getCollapseFlux() const;  // Mass flux into PDC
    
    // Tephra dispersal parameters
    double getSO2_emission_rate() const;  // kg/s
    double getAsh_emission_rate() const;  // kg/s
    
private:
    EruptionColumnParameters params;
    std::vector<ColumnState> profile;
    
    double max_height;
    double NBL_height;        // Neutral Buoyancy Level
    bool column_collapsed;
    double collapse_height;
    
    // Atmospheric model
    double atmosphericDensity(double z) const;
    double atmosphericTemperature(double z) const;
    double atmosphericPressure(double z) const;
    
    // Entrainment coefficient (Carazzo et al., 2008)
    double entrainmentCoefficient(double Ri) const;  // Richardson number dependent
    
    // Particle settling
    double particleSettlingVelocity(double diameter, double gas_density) const;
};

// =============================================================================
// Pyroclastic Density Current (PDC) Model
// =============================================================================

/**
 * @brief PDC types
 */
enum class PDCType {
    COLUMN_COLLAPSE,        // From eruption column collapse
    DOME_COLLAPSE,          // From lava dome failure
    DIRECTED_BLAST,         // Lateral blast (Mt. St. Helens)
    BOILING_OVER,           // Fountain collapse
    SURGE                   // Dilute turbulent current
};

/**
 * @brief PDC initial conditions
 */
struct PDCSourceConditions {
    PDCType type = PDCType::COLUMN_COLLAPSE;
    double mass_flux = 1e7;         // kg/s
    double initial_velocity = 50.0; // m/s
    double initial_temperature = 700.0;  // K
    double particle_concentration = 0.01; // Volume fraction
    double particle_size = 1e-3;    // m
    double source_height = 1000.0;  // m above vent
    double initial_radius = 500.0;  // m
    
    // Direction (for directed blasts)
    double azimuth = 0.0;           // degrees from N
    double inclination = 0.0;       // degrees from horizontal
};

/**
 * @brief PDC state at a point
 */
struct PDCState {
    double x, y, z;                 // Position (m)
    double u, v, w;                 // Velocity (m/s)
    double thickness = 0.0;         // m
    double temperature = 600.0;     // K
    double density = 10.0;          // kg/m³
    double particle_conc = 0.005;   // Volume fraction
    double velocity_magnitude() const { 
        return std::sqrt(u*u + v*v + w*w); 
    }
    double dynamic_pressure() const {
        return 0.5 * density * velocity_magnitude() * velocity_magnitude();
    }
};

/**
 * @brief Pyroclastic Density Current model
 * 
 * Models dilute to concentrated PDCs using:
 * - Shallow water + suspension equations (Denlinger & Iverson, 2001)
 * - Two-layer model (Valentine, 1998)
 * - Density stratification and air entrainment
 * - Sedimentation and erosion
 */
class PDCModel {
public:
    PDCModel();
    
    void initialize(const PDCSourceConditions& source);
    void setTopography(std::function<double(double, double)> elevation_func);
    
    /**
     * @brief Advance PDC simulation by one time step
     */
    void update(double dt);
    
    /**
     * @brief Run simulation until PDC halts
     */
    void runToCompletion(double max_time = 3600.0);
    
    // Get results
    double getRunoutDistance() const;
    double getMaxVelocity() const;
    double getCurrentExtent() const;
    
    // Hazard metrics at a point
    double getDynamicPressure(double x, double y) const;
    double getTemperature(double x, double y) const;
    double getDepositThickness(double x, double y) const;
    
    // Is point within PDC inundation zone?
    bool isInundated(double x, double y) const;
    
    // Get time series at a point
    std::vector<std::pair<double, PDCState>> getTimeSeriesAt(double x, double y) const;
    
private:
    PDCSourceConditions source;
    std::vector<PDCState> current_state;
    std::function<double(double, double)> topography;
    
    double current_time;
    double runout_distance;
    
    // Grid for solution
    int nx, ny;
    double dx, dy;
    double xmin, xmax, ymin, ymax;
    
    std::vector<double> deposit_thickness;
    std::vector<double> max_dynamic_pressure;
    std::vector<double> max_temperature;
    
    // Physics
    double computeDragForce(const PDCState& state) const;
    double computeEntrainment(const PDCState& state) const;
    double computeSedimentation(const PDCState& state) const;
    double computeErosion(const PDCState& state) const;
};

// =============================================================================
// Lava Flow Model
// =============================================================================

/**
 * @brief Lava flow source parameters
 */
struct LavaSourceParameters {
    double effusion_rate = 10.0;    // m³/s (DRE)
    double temperature = 1373.15;   // K (1100°C)
    double vent_x = 0.0;            // m
    double vent_y = 0.0;            // m
    double vent_z = 0.0;            // m (elevation)
    double initial_thickness = 2.0; // m
};

/**
 * @brief Lava flow state
 */
struct LavaFlowState {
    double thickness = 0.0;         // m
    double velocity_x = 0.0;        // m/s
    double velocity_y = 0.0;        // m/s
    double temperature = 1373.15;   // K
    double viscosity = 1e4;         // Pa·s
    double crust_thickness = 0.0;   // m
    double crystal_fraction = 0.3;
    
    // Is the flow mobile (not solidified)?
    bool isMobile() const { 
        return temperature > 1000.0 && crystal_fraction < 0.7; 
    }
};

/**
 * @brief Lava flow model
 * 
 * Models lava flow propagation using:
 * - Shallow water/Bingham fluid equations
 * - Temperature-dependent rheology
 * - Crust formation and insulation
 * - Topographic channeling
 * 
 * Based on Harris & Rowland (2001), Crisci et al. (2004)
 */
class LavaFlowModel {
public:
    LavaFlowModel();
    
    void initialize(const LavaSourceParameters& source,
                   const MagmaProperties& magma);
    void setTopography(std::function<double(double, double)> elevation_func);
    void setEffusionRate(double rate_m3_s);
    
    /**
     * @brief Update lava flow for one time step
     */
    void update(double dt);
    
    /**
     * @brief Run simulation for specified duration
     */
    void run(double duration);
    
    // Get results
    double getFlowArea() const;
    double getFlowLength() const;
    double getFlowVolume() const;
    
    // State at a point
    const LavaFlowState* getStateAt(double x, double y) const;
    double getThicknessAt(double x, double y) const;
    bool isInundated(double x, double y) const;
    
    // Thermal output
    double getRadiativeHeatFlux() const;  // W
    
private:
    LavaSourceParameters source;
    MagmaProperties magma;
    std::function<double(double, double)> topography;
    
    // Grid
    int nx, ny;
    double dx, dy;
    double xmin, xmax, ymin, ymax;
    std::vector<LavaFlowState> state;
    
    // Physics
    double computeViscosity(double T, double phi_crystal) const;
    double computeYieldStrength(double phi_crystal) const;
    double computeCoolingRate(const LavaFlowState& state) const;
    double computeCrustGrowth(const LavaFlowState& state, double dt) const;
    
    // Numerical methods
    void computeFluxes();
    void advectTemperature(double dt);
};

// =============================================================================
// Lahar (Volcanic Debris Flow) Model
// =============================================================================

/**
 * @brief Lahar types
 */
enum class LaharType {
    PRIMARY_HOT,        // From eruption (hot)
    PRIMARY_COLD,       // From eruption (rain-triggered)
    SECONDARY,          // Post-eruption remobilization
    JÖKULHLAUP,         // Glacier outburst flood
    DAM_BREAK           // Crater lake breach
};

/**
 * @brief Lahar source conditions
 */
struct LaharSourceConditions {
    LaharType type = LaharType::PRIMARY_COLD;
    double volume = 1e7;            // m³
    double peak_discharge = 1e4;    // m³/s
    double sediment_concentration = 0.5;  // Volume fraction
    double temperature = 300.0;     // K (can be hot for primary)
    double source_x = 0.0;          // m
    double source_y = 0.0;          // m
    double source_elevation = 2000.0; // m
    
    // Duration of source
    double duration = 600.0;        // s
    
    // Hydrograph shape
    enum class HydrographShape {
        TRIANGULAR,
        TRAPEZOIDAL,
        DAM_BREAK
    };
    HydrographShape shape = HydrographShape::TRIANGULAR;
};

/**
 * @brief Lahar state
 */
struct LaharState {
    double depth = 0.0;             // m
    double velocity_x = 0.0;        // m/s
    double velocity_y = 0.0;        // m/s
    double sediment_conc = 0.5;     // Volume fraction
    double temperature = 300.0;     // K
    
    double velocity_magnitude() const {
        return std::sqrt(velocity_x*velocity_x + velocity_y*velocity_y);
    }
    
    // Bulk density
    double density() const {
        return 1000.0 * (1.0 - sediment_conc) + 2650.0 * sediment_conc;
    }
    
    // Dynamic pressure
    double dynamic_pressure() const {
        return 0.5 * density() * velocity_magnitude() * velocity_magnitude();
    }
};

/**
 * @brief Lahar model
 * 
 * Models volcanic debris flows using:
 * - Shallow water equations with rheology
 * - Variable sediment concentration
 * - Erosion and deposition
 * - Bulking and debulking
 * 
 * Based on Iverson (1997), Manville et al. (1998)
 */
class LaharModel {
public:
    LaharModel();
    
    void initialize(const LaharSourceConditions& source);
    void setTopography(std::function<double(double, double)> elevation_func);
    void setChannel(std::function<double(double, double)> channel_width);
    
    /**
     * @brief Update lahar simulation
     */
    void update(double dt);
    
    /**
     * @brief Run simulation until lahar stops
     */
    void runToCompletion(double max_time = 7200.0);
    
    // Results
    double getRunoutDistance() const;
    double getCurrentVolume() const;
    double getPeakDischarge() const;
    
    // State at a point
    double getDepthAt(double x, double y) const;
    double getVelocityAt(double x, double y) const;
    bool isInundated(double x, double y) const;
    
    // Travel time to a point
    double getTravelTimeTo(double x, double y) const;
    
    // Maximum values at a point (for hazard mapping)
    double getMaxDepth(double x, double y) const;
    double getMaxVelocity(double x, double y) const;
    double getMaxDynamicPressure(double x, double y) const;
    
private:
    LaharSourceConditions source;
    std::function<double(double, double)> topography;
    std::function<double(double, double)> channel_width;
    
    // Grid
    int nx, ny;
    double dx, dy;
    std::vector<LaharState> state;
    
    double current_time;
    
    // Tracking
    std::vector<double> max_depth;
    std::vector<double> max_velocity;
    std::vector<double> arrival_time;
    
    // Physics
    double computeFriction(const LaharState& state) const;
    double computeErosion(const LaharState& state) const;
    double computeDeposition(const LaharState& state) const;
};

// =============================================================================
// Volcanic Deformation Model
// =============================================================================

/**
 * @brief Volcanic deformation source parameters
 */
struct VolcanicDeformationSource {
    DeformationSourceType type = DeformationSourceType::MOGI;
    
    // Location
    double x = 0.0;                 // m (horizontal)
    double y = 0.0;                 // m (horizontal)
    double depth = 5000.0;          // m (below surface)
    
    // Source strength
    double volume_change = 1e6;     // m³ (for Mogi)
    double pressure_change = 10e6;  // Pa (for McTigue)
    
    // Geometric parameters (for non-spherical sources)
    double radius = 1000.0;         // m (sphere/cavity radius)
    double semi_major = 2000.0;     // m (spheroid a)
    double semi_minor = 1000.0;     // m (spheroid b)
    double strike = 0.0;            // degrees (for dikes/sills)
    double dip = 90.0;              // degrees
    double length = 5000.0;         // m (dike/sill)
    double width = 2000.0;          // m (dike/sill)
    double opening = 1.0;           // m (dike/sill opening)
};

/**
 * @brief Volcanic deformation model
 * 
 * Computes surface displacements, tilts, and strains from
 * subsurface pressure/volume sources.
 * 
 * Supports multiple analytical source types and can be
 * coupled with FEM for complex geometries.
 */
class VolcanicDeformationModel {
public:
    VolcanicDeformationModel();
    
    void addSource(const VolcanicDeformationSource& source);
    void clearSources();
    
    // Material properties
    void setElasticProperties(double shear_modulus, double poisson_ratio);
    
    /**
     * @brief Compute displacement at a surface point
     * 
     * @param x, y Surface coordinates (m)
     * @param ux, uy, uz Output displacements (m)
     */
    void computeDisplacement(double x, double y,
                            double& ux, double& uy, double& uz) const;
    
    /**
     * @brief Compute tilt at a surface point
     * 
     * @param x, y Surface coordinates (m)
     * @param tilt_x, tilt_y Output tilts (radians)
     */
    void computeTilt(double x, double y,
                    double& tilt_x, double& tilt_y) const;
    
    /**
     * @brief Compute strain at a surface point
     */
    void computeStrain(double x, double y,
                      double& exx, double& eyy, double& exy) const;
    
    /**
     * @brief Compute gravity change at a surface point
     */
    double computeGravityChange(double x, double y) const;
    
    // Displacement field on a grid
    void computeDisplacementField(const std::vector<double>& x_coords,
                                 const std::vector<double>& y_coords,
                                 std::vector<double>& ux,
                                 std::vector<double>& uy,
                                 std::vector<double>& uz) const;
    
    // InSAR line-of-sight displacement
    double computeLOS_displacement(double x, double y,
                                  double look_azimuth, double incidence_angle) const;
    
private:
    std::vector<VolcanicDeformationSource> sources;
    double shear_modulus;
    double poisson_ratio;
    
    // Analytical solutions
    void mogiDisplacement(double x, double y, double depth, double dV,
                         double& ux, double& uy, double& uz) const;
    void mctigueDisplacement(double x, double y, double depth, double radius, double dP,
                            double& ux, double& uy, double& uz) const;
    void okataDikeDisplacement(double x, double y, const VolcanicDeformationSource& src,
                              double& ux, double& uy, double& uz) const;
    void okataSillDisplacement(double x, double y, const VolcanicDeformationSource& src,
                              double& ux, double& uy, double& uz) const;
    void spheroidDisplacement(double x, double y, const VolcanicDeformationSource& src,
                             double& ux, double& uy, double& uz) const;
};

// =============================================================================
// Volcanic Seismicity Model
// =============================================================================

/**
 * @brief Volcanic seismic event parameters
 */
struct VolcanicSeismicEvent {
    VolcanicSeismicEventType type = VolcanicSeismicEventType::LP;
    
    // Location
    double x = 0.0;                 // m
    double y = 0.0;                 // m
    double z = -2000.0;             // m (depth, negative = below surface)
    
    // Timing
    double origin_time = 0.0;       // s
    
    // Size
    double magnitude = 1.5;         // Local magnitude
    double seismic_moment = 1e12;   // N·m (for VT)
    
    // LP/VLP source
    double dominant_frequency = 1.5;     // Hz
    double quality_factor = 20.0;        // Q
    double crack_length = 100.0;         // m
    double crack_width = 50.0;           // m
    double pressure_transient = 1e5;     // Pa
    
    // Tremor
    double tremor_duration = 300.0;      // s
    double amplitude_modulation = 0.3;   // Fraction
};

/**
 * @brief Volcanic tremor source
 */
struct VolcanicTremorSource {
    double x = 0.0, y = 0.0, z = -1000.0;  // Location (m)
    
    // Spectral characteristics
    std::vector<double> peak_frequencies = {1.0, 2.0, 3.0};  // Hz
    std::vector<double> peak_amplitudes = {1.0, 0.5, 0.3};
    double bandwidth = 0.5;             // Hz
    
    // Temporal
    double duration = 600.0;            // s
    double amplitude_variation = 0.2;   // Coefficient of variation
    
    // Source mechanism
    enum class Mechanism {
        FLUID_OSCILLATION,      // Resonating fluid-filled crack
        BUBBLE_BURST,           // Repeated bubble bursts
        REPEATING_LP,           // Closely spaced LP events
        CONDUIT_OSCILLATION     // Organ-pipe mode
    };
    Mechanism mechanism = Mechanism::FLUID_OSCILLATION;
};

/**
 * @brief Volcanic seismicity model
 * 
 * Models volcanic earthquake sources including:
 * - VT earthquakes (double-couple)
 * - LP events (fluid-filled crack resonance)
 * - VLP events (conduit resonance)
 * - Volcanic tremor
 * - Explosion quakes
 * 
 * Based on Chouet (1996), Kumagai (2002)
 */
class VolcanicSeismicityModel {
public:
    VolcanicSeismicityModel();
    
    /**
     * @brief Generate VT event source time function
     */
    void generateVT_source(const VolcanicSeismicEvent& event,
                          std::vector<double>& times,
                          std::vector<double>& moment_rate);
    
    /**
     * @brief Generate LP event source time function
     * 
     * Uses crack resonance model (Chouet, 1986)
     */
    void generateLP_source(const VolcanicSeismicEvent& event,
                          std::vector<double>& times,
                          std::vector<std::array<double, 6>>& moment_tensor_rate);
    
    /**
     * @brief Generate VLP event source time function
     */
    void generateVLP_source(const VolcanicSeismicEvent& event,
                           std::vector<double>& times,
                           std::vector<std::array<double, 6>>& moment_tensor_rate);
    
    /**
     * @brief Generate tremor source
     */
    void generateTremor_source(const VolcanicTremorSource& source,
                              std::vector<double>& times,
                              std::vector<double>& amplitude);
    
    /**
     * @brief Generate explosion quake source
     */
    void generateExplosionQuake(const VolcanicSeismicEvent& event,
                               std::vector<double>& times,
                               std::vector<std::array<double, 6>>& moment_tensor_rate);
    
    // Crack resonance model (Chouet, 1986)
    struct CrackResonance {
        double length;           // Crack length (m)
        double width;            // Crack width (m)
        double stiffness;        // Crack stiffness (Pa/m)
        double fluid_density;    // Fluid density (kg/m³)
        double fluid_viscosity;  // Fluid viscosity (Pa·s)
        double fluid_velocity;   // Acoustic velocity (m/s)
        
        // Compute resonant frequency
        double fundamentalFrequency() const;
        double qualityFactor() const;
        
        // Mode shapes
        std::vector<double> getModeFrequencies(int n_modes) const;
    };
    
    void setCrackModel(const CrackResonance& crack);
    
private:
    CrackResonance crack_model;
    
    // Helper functions
    double decayingOscillation(double t, double f, double Q) const;
    void applySourceTimeFunction(double t, double rise_time, 
                                double& value, double& derivative) const;
};

// =============================================================================
// Volcanic Gas Emission Model
// =============================================================================

/**
 * @brief Volcanic gas species
 */
enum class VolcanicGas {
    H2O,            // Water vapor
    CO2,            // Carbon dioxide
    SO2,            // Sulfur dioxide
    H2S,            // Hydrogen sulfide
    HCl,            // Hydrogen chloride
    HF,             // Hydrogen fluoride
    CO,             // Carbon monoxide
    H2,             // Hydrogen
    He,             // Helium
    Rn              // Radon (radioactive tracer)
};

/**
 * @brief Gas emission source
 */
struct GasEmissionSource {
    double x = 0.0, y = 0.0, z = 0.0;  // Location (m)
    double emission_rate = 1000.0;      // kg/s (total)
    double temperature = 373.15;        // K (100°C)
    double velocity = 10.0;             // m/s (exit velocity)
    
    // Gas composition (mole fractions)
    std::map<VolcanicGas, double> composition = {
        {VolcanicGas::H2O, 0.90},
        {VolcanicGas::CO2, 0.05},
        {VolcanicGas::SO2, 0.03},
        {VolcanicGas::H2S, 0.01},
        {VolcanicGas::HCl, 0.01}
    };
};

/**
 * @brief Volcanic gas emission and dispersal model
 * 
 * Models volcanic degassing and plume dispersal:
 * - Source emission from vents, fumaroles, ground
 * - Atmospheric dispersion (Gaussian plume)
 * - Chemical transformations
 * - Deposition (wet and dry)
 */
class VolcanicGasModel {
public:
    VolcanicGasModel();
    
    void addSource(const GasEmissionSource& source);
    void clearSources();
    
    // Meteorological conditions
    void setWindSpeed(double speed_m_s);
    void setWindDirection(double direction_deg);  // From direction
    void setAtmosphericStability(char stability_class);  // A-F (Pasquill)
    void setMixingHeight(double height_m);
    
    /**
     * @brief Compute ground-level concentration at a point
     */
    double getConcentration(VolcanicGas gas, double x, double y) const;  // kg/m³
    
    /**
     * @brief Compute concentration in ppm
     */
    double getConcentrationPPM(VolcanicGas gas, double x, double y) const;
    
    /**
     * @brief Compute SO2 column density (for DOAS/satellite)
     */
    double getSO2_column_density(double x, double y) const;  // kg/m²
    
    /**
     * @brief Compute total SO2 emission rate
     */
    double getTotalSO2_emission() const;  // kg/s
    
    /**
     * @brief Compute total CO2 emission rate
     */
    double getTotalCO2_emission() const;  // kg/s
    
    // Hazard thresholds
    bool exceedsHealthLimit(VolcanicGas gas, double x, double y) const;
    
    // Get molecular weight
    static double getMolecularWeight(VolcanicGas gas);
    
private:
    std::vector<GasEmissionSource> sources;
    double wind_speed;
    double wind_direction;
    char stability_class;
    double mixing_height;
    
    // Gaussian dispersion parameters
    double sigmaY(double x) const;
    double sigmaZ(double x) const;
    
    // Plume rise
    double plumeRise(const GasEmissionSource& source, double x) const;
};

// =============================================================================
// Tephra/Ash Dispersal Model
// =============================================================================

/**
 * @brief Tephra particle characteristics
 */
struct TephraParticle {
    double diameter = 1e-3;         // m
    double density = 2500.0;        // kg/m³
    double shape_factor = 0.7;      // Sphericity
    
    // Settling velocity (Wilson & Huang, 1979)
    double settlingVelocity(double air_density, double air_viscosity) const;
    
    // Drag coefficient
    double dragCoefficient(double Re) const;
};

/**
 * @brief Tephra source parameters
 */
struct TephraSourceParameters {
    double x = 0.0, y = 0.0, z = 0.0;  // Vent location (m)
    double column_height = 15000.0;     // m above vent
    double mass_eruption_rate = 1e7;    // kg/s
    double vent_radius = 100.0;         // m
    
    // Particle size distribution (Gaussian in phi scale)
    double phi_mean = 2.0;              // phi = -log2(d_mm)
    double phi_stddev = 2.0;
    
    // Particle density
    double particle_density = 2500.0;   // kg/m³
};

/**
 * @brief Tephra dispersal model (Ash3d-style)
 * 
 * Models volcanic ash dispersal using:
 * - Advection by wind field
 * - Turbulent diffusion
 * - Gravitational settling
 * - Aggregation (optional)
 * - Ground deposition
 * 
 * Based on Folch et al. (2009), Mastin et al. (2009)
 */
class TephraDispersalModel {
public:
    TephraDispersalModel();
    
    void initialize(const TephraSourceParameters& source);
    
    // Wind field (can be 3D)
    void setWindField(std::function<void(double x, double y, double z,
                                        double& u, double& v, double& w)> wind);
    
    // Atmospheric properties
    void setAtmosphere(std::function<double(double z)> density,
                      std::function<double(double z)> viscosity);
    
    /**
     * @brief Run dispersal calculation
     */
    void run(double duration);
    
    /**
     * @brief Get deposit thickness at a point
     */
    double getDepositThickness(double x, double y) const;  // m
    
    /**
     * @brief Get deposit mass loading at a point
     */
    double getDepositLoading(double x, double y) const;  // kg/m²
    
    /**
     * @brief Get airborne ash concentration at a point
     */
    double getAirborneConcentration(double x, double y, double z) const;  // kg/m³
    
    // Hazard metrics
    bool exceedsAviationLimit(double x, double y, double z) const;  // >2 mg/m³
    double getVisibility(double x, double y, double z) const;       // m
    
    // Isopach contours
    void getIsopachContours(const std::vector<double>& thicknesses,
                           std::vector<std::vector<std::pair<double, double>>>& contours) const;
    
private:
    TephraSourceParameters source;
    
    // 3D grid
    int nx, ny, nz;
    double dx, dy, dz;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    
    // Particle size bins
    int n_size_bins;
    std::vector<TephraParticle> particles;
    std::vector<double> mass_fractions;
    
    // Concentration field (per size bin)
    std::vector<std::vector<double>> concentration;
    
    // Deposit (per size bin)
    std::vector<std::vector<double>> deposit;
    
    // Wind and atmosphere functions
    std::function<void(double, double, double, double&, double&, double&)> wind_field;
    std::function<double(double)> air_density;
    std::function<double(double)> air_viscosity;
    
    // Numerical methods
    void advect(double dt);
    void diffuse(double dt);
    void settle(double dt);
    void deposit_particles();
    
    // Initialize particle size distribution
    void initializeParticleSizes();
};

// =============================================================================
// Coupled Volcano System
// =============================================================================

/**
 * @brief Coupled volcanic system configuration
 */
struct VolcanoSystemConfig {
    // Physical components to enable
    bool enable_magma_chamber = true;
    bool enable_conduit_flow = true;
    bool enable_eruption_column = true;
    bool enable_pdc = false;
    bool enable_lava_flow = false;
    bool enable_lahar = false;
    bool enable_deformation = true;
    bool enable_seismicity = true;
    bool enable_gas_emission = true;
    bool enable_tephra = false;
    
    // Coupling
    bool couple_chamber_conduit = true;
    bool couple_conduit_column = true;
    bool couple_column_pdc = true;
    bool couple_chamber_deformation = true;
    bool couple_chamber_seismicity = true;
    
    // Time stepping
    double dt_min = 0.01;           // s
    double dt_max = 60.0;           // s
    double end_time = 3600.0;       // s
    
    // Output
    int output_interval = 100;
    std::string output_directory = "output/volcano";
};

/**
 * @brief Fully coupled volcanic system model
 * 
 * Integrates all volcanic physics models:
 * - Magma chamber → Conduit → Eruption column → PDC/Lava/Tephra
 * - Deformation and seismicity from chamber processes
 * - Gas emissions from all sources
 * 
 * Supports:
 * - Unrest scenarios (pre-eruption monitoring)
 * - Eruption scenarios (full eruption dynamics)
 * - Post-eruption scenarios (cooling, degassing)
 */
class CoupledVolcanoSystem {
public:
    CoupledVolcanoSystem();
    
    /**
     * @brief Initialize the volcanic system
     */
    void initialize(const VolcanoSystemConfig& config,
                   const MagmaChamberGeometry& chamber_geometry,
                   const MagmaProperties& magma,
                   const ConduitGeometry& conduit);
    
    /**
     * @brief Set topography for surface flows
     */
    void setTopography(std::function<double(double, double)> elevation);
    
    /**
     * @brief Set atmospheric conditions
     */
    void setAtmosphere(std::function<void(double, double&, double&, double&)> profile);
    
    /**
     * @brief Set magma recharge rate
     */
    void setRechargeRate(double rate_kg_s);
    
    /**
     * @brief Trigger an eruption
     */
    void triggerEruption(EruptionType type = EruptionType::SUBPLINIAN);
    
    /**
     * @brief Update the system by one time step
     */
    void update(double dt);
    
    /**
     * @brief Run simulation to end time
     */
    void run();
    
    // State queries
    bool isErupting() const { return eruption_active; }
    EruptionType getCurrentEruptionType() const { return current_eruption_type; }
    double getMassEruptionRate() const;
    double getColumnHeight() const;
    
    // Access component models
    const MagmaChamberModel& getChamber() const { return *chamber; }
    const ConduitFlowModel& getConduit() const { return *conduit; }
    const EruptionColumnModel& getColumn() const { return *column; }
    const VolcanicDeformationModel& getDeformation() const { return *deformation; }
    
    // Monitoring outputs
    struct MonitoringOutput {
        double time;
        double chamber_pressure;
        double chamber_overpressure;
        double radial_displacement;
        double vertical_displacement;
        double SO2_emission_rate;
        double seismicity_rate;
        double tremor_amplitude;
        double mass_eruption_rate;
        double column_height;
    };
    
    std::vector<MonitoringOutput> getMonitoringTimeSeries() const { return monitoring; }
    
    /**
     * @brief Write output files
     */
    void writeOutput(const std::string& output_dir) const;
    
private:
    VolcanoSystemConfig config;
    
    // Component models
    std::unique_ptr<MagmaChamberModel> chamber;
    std::unique_ptr<ConduitFlowModel> conduit;
    std::unique_ptr<EruptionColumnModel> column;
    std::unique_ptr<PDCModel> pdc;
    std::unique_ptr<LavaFlowModel> lava;
    std::unique_ptr<LaharModel> lahar;
    std::unique_ptr<VolcanicDeformationModel> deformation;
    std::unique_ptr<VolcanicSeismicityModel> seismicity;
    std::unique_ptr<VolcanicGasModel> gas;
    std::unique_ptr<TephraDispersalModel> tephra;
    
    // State
    double current_time;
    bool eruption_active;
    EruptionType current_eruption_type;
    
    // External conditions
    std::function<double(double, double)> topography;
    
    // Monitoring time series
    std::vector<MonitoringOutput> monitoring;
    
    // Coupling functions
    void updateChamberConduitCoupling();
    void updateConduitColumnCoupling();
    void updateColumnPDCCoupling();
    void updateDeformationFromChamber();
    void updateSeismicityFromChamber();
    void recordMonitoring();
};

// =============================================================================
// Pre-defined Volcanic Scenarios
// =============================================================================

namespace VolcanoScenarios {
    
    /**
     * @brief Mount St. Helens 1980-style lateral blast
     */
    VolcanoSystemConfig mountStHelens1980();
    MagmaChamberGeometry mountStHelens1980_chamber();
    MagmaProperties mountStHelens1980_magma();
    
    /**
     * @brief Pinatubo 1991-style Plinian eruption
     */
    VolcanoSystemConfig pinatubo1991();
    MagmaChamberGeometry pinatubo1991_chamber();
    MagmaProperties pinatubo1991_magma();
    
    /**
     * @brief Kilauea effusive eruption
     */
    VolcanoSystemConfig kilauea_effusive();
    MagmaChamberGeometry kilauea_chamber();
    MagmaProperties kilauea_magma();
    
    /**
     * @brief Stromboli persistent activity
     */
    VolcanoSystemConfig stromboli_persistent();
    MagmaChamberGeometry stromboli_chamber();
    MagmaProperties stromboli_magma();
    
    /**
     * @brief Yellowstone caldera unrest
     */
    VolcanoSystemConfig yellowstone_unrest();
    MagmaChamberGeometry yellowstone_chamber();
    MagmaProperties yellowstone_magma();
    
    /**
     * @brief Campi Flegrei unrest (bradyseism)
     */
    VolcanoSystemConfig campi_flegrei_unrest();
    MagmaChamberGeometry campi_flegrei_chamber();
    MagmaProperties campi_flegrei_magma();
    
    /**
     * @brief Merapi dome-forming eruption
     */
    VolcanoSystemConfig merapi_dome();
    MagmaChamberGeometry merapi_chamber();
    MagmaProperties merapi_magma();
    
    /**
     * @brief Generic VEI-4 eruption
     */
    VolcanoSystemConfig generic_VEI4();
    
    /**
     * @brief Generic VEI-6 eruption (Pinatubo-scale)
     */
    VolcanoSystemConfig generic_VEI6();
}

// =============================================================================
// Configuration Reader
// =============================================================================

/**
 * @brief Configuration reader for volcano models
 */
struct VolcanoConfigReader {
    /**
     * @brief Parse volcano system configuration from file
     */
    static bool parseVolcanoConfig(const std::map<std::string, std::string>& section,
                                  VolcanoSystemConfig& config);
    
    /**
     * @brief Parse magma chamber configuration
     */
    static bool parseChamberConfig(const std::map<std::string, std::string>& section,
                                  MagmaChamberGeometry& geometry,
                                  MagmaChamberState& state);
    
    /**
     * @brief Parse magma properties
     */
    static bool parseMagmaProperties(const std::map<std::string, std::string>& section,
                                    MagmaProperties& props);
    
    /**
     * @brief Parse conduit configuration
     */
    static bool parseConduitConfig(const std::map<std::string, std::string>& section,
                                  ConduitGeometry& geometry);
    
    /**
     * @brief Parse deformation source
     */
    static bool parseDeformationSource(const std::map<std::string, std::string>& section,
                                      VolcanicDeformationSource& source);
    
    /**
     * @brief Parse gas emission source
     */
    static bool parseGasSource(const std::map<std::string, std::string>& section,
                              GasEmissionSource& source);
};

// =============================================================================
// Utility Functions
// =============================================================================

namespace VolcanoUtils {
    
    /**
     * @brief Convert VEI to approximate parameters
     */
    void VEI_to_parameters(int VEI, double& volume_km3, double& column_height_km,
                          double& mass_eruption_rate);
    
    /**
     * @brief Estimate column height from mass eruption rate (Mastin et al., 2009)
     */
    double columnHeight_from_MER(double mass_eruption_rate);
    
    /**
     * @brief Estimate mass eruption rate from column height
     */
    double MER_from_columnHeight(double height_km);
    
    /**
     * @brief Estimate eruption duration from volume and MER
     */
    double estimateDuration(double volume_km3, double mass_eruption_rate);
    
    /**
     * @brief Calculate magma viscosity (Giordano et al., 2008)
     */
    double giordano2008_viscosity(double T_celsius, double SiO2, double TiO2,
                                 double Al2O3, double FeO, double MnO,
                                 double MgO, double CaO, double Na2O,
                                 double K2O, double P2O5, double H2O);
    
    /**
     * @brief Water solubility in silicate melt (Newman & Lowenstern, 2002)
     */
    double H2O_solubility(double P_MPa, double T_celsius, double SiO2_wt);
    
    /**
     * @brief CO2 solubility in silicate melt
     */
    double CO2_solubility(double P_MPa, double T_celsius, double SiO2_wt);
    
    /**
     * @brief Convert phi to diameter and back
     */
    double phi_to_mm(double phi);
    double mm_to_phi(double mm);
    
    /**
     * @brief Great circle distance (for regional effects)
     */
    double haversineDistance(double lat1, double lon1, double lat2, double lon2);
}

} // namespace FSRM

#endif // VOLCANO_MODEL_HPP
