/**
 * @file NuclearAirburstEffects.hpp
 * @brief Parameterised nuclear airburst effects calculator
 *
 * Implements all analytic / semi-empirical models needed to reproduce the
 * Berkeley 100 kT example (and any other airburst scenario) from a compact
 * set of input parameters.
 *
 * Physics models and references
 * ─────────────────────────────
 *   A  Blast overpressure         — piecewise fit calibrated to Glasstone-Dolan (1977)
 *   B  Dynamic (wind) pressure    — Rankine-Hugoniot relation
 *   C  Thermal radiation fluence  — inverse-square with atmospheric attenuation
 *   D  Prompt nuclear radiation   — gamma / neutron dose with exponential shielding
 *   E  Fallout (Gaussian plume)   — Way-Wigner decay, wind-rotated deposition
 *   F  EMP (E1 component)         — double-exponential waveform
 *   G  Fireball & shock           — Sedov-Taylor self-similar solution
 *   H  Ground-coupled seismics    — acoustic-impedance PGV, Murphy mb scaling
 *   I  Casualty / damage mapping  — annular-ring population model
 *
 * The class is entirely standalone (no PETSc dependency) so it can be used
 * from lightweight driver programs, unit tests, or the full FSRM simulator.
 *
 * @author FSRM Development Team
 */

#ifndef NUCLEAR_AIRBURST_EFFECTS_HPP
#define NUCLEAR_AIRBURST_EFFECTS_HPP

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

namespace FSRM
{

// =============================================================================
// Physical constants
// =============================================================================
namespace AirburstConstants
{
  constexpr double JOULES_PER_KT       = 4.184e12;    // J / kt TNT
  constexpr double STANDARD_PRESSURE   = 101.325;     // kPa
  constexpr double SEA_LEVEL_DENSITY   = 1.225;       // kg / m^3
  constexpr double STEFAN_BOLTZMANN    = 5.670374e-8;  // W / (m^2 K^4)
  constexpr double PI                  = 3.14159265358979323846;
  constexpr double EARTH_RADIUS_M      = 6371000.0;
  constexpr double DEG_TO_RAD          = PI / 180.0;
  constexpr double RAD_TO_DEG          = 180.0 / PI;
} // namespace AirburstConstants

// =============================================================================
// Input parameter structures
// =============================================================================

/**
 * @brief All parameters that define a single airburst scenario
 *
 * Every field has a sensible default so that only the values of interest
 * need to be overridden from a config file.
 */
struct AirburstParameters
{
  // --- Device ---
  double yield_kt         = 100.0;    ///< Yield (kilotons TNT)
  double burst_height_m   = 50.0;     ///< Height of burst AGL (m)
  double fission_fraction = 0.50;     ///< Fraction of yield from fission

  // --- Location ---
  double latitude         = 37.8716;  ///< Ground zero latitude (deg N)
  double longitude        = -122.2727;///< Ground zero longitude (deg E)
  std::string location_name = "UC Berkeley campus";

  // --- Atmosphere ---
  double thermal_attenuation_length_m = 18000.0;  ///< 1/e length for absorption (m)

  // --- Wind ---
  double wind_speed_surface_mps  = 5.0;   ///< Surface wind speed (m/s)
  double wind_speed_altitude_mps = 15.0;  ///< Wind at cloud height (m/s)
  double wind_direction_deg      = 315.0; ///< Meteorological (from), degrees

  // --- Ground ---
  double ground_density_kgm3   = 2700.0;  ///< Near-surface rock density (kg/m^3)
  double ground_vs_mps         = 3000.0;  ///< Near-surface shear-wave speed (m/s)
  double seismic_coupling_eff  = 0.001;   ///< Air-to-ground coupling efficiency

  // --- Population ---
  double population_density_urban = 4000.0; ///< ppl / km^2
  double population_density_sf    = 7000.0;
  double population_density_sub   = 2000.0;

  // --- Prompt radiation reference ---
  double prompt_dose_ref_gy      = 10.0;  ///< Gy at 1 km for 1 kt
  double prompt_atten_length_m   = 2500.0;///< 1/e length in air (m)

  // --- EMP ---
  double emp_e0_vm          = 25000.0;  ///< Reference E1 peak (V/m) at 100 kt
  double emp_decay_length_m = 30000.0;  ///< E1 spatial decay length (m)
  double emp_tau_rise_s     = 2.5e-9;   ///< E1 rise time (s)
  double emp_tau_decay_s    = 1.0e-6;   ///< E1 decay time (s)
};

// =============================================================================
// Damage / threshold tables
// =============================================================================

struct DamageThreshold
{
  std::string label;
  double value;        ///< threshold in the native unit (kPa, kJ/m^2, Gy …)
};

// =============================================================================
// Result structure returned by computeAll()
// =============================================================================

/**
 * @brief Results at a single evaluation point (ground range r from GZ)
 */
struct AirburstPointResult
{
  double ground_range_m              = 0.0;
  double slant_range_m               = 0.0;

  // Blast
  double peak_overpressure_kpa       = 0.0;
  double dynamic_pressure_kpa        = 0.0;
  double positive_phase_duration_s   = 0.0;

  // Thermal
  double thermal_fluence_kjm2        = 0.0;

  // Prompt radiation
  double prompt_dose_gy              = 0.0;

  // Fallout (requires x,y not just r)
  double fallout_activity_norm       = 0.0;
  double fallout_dose_rate_svh       = 0.0;

  // EMP
  double emp_e1_peak_vm              = 0.0;

  // Seismic
  double pgv_ms                      = 0.0;
  double pgv_cms                     = 0.0;
};

/**
 * @brief Scenario-global results (scalar quantities)
 */
struct AirburstScenarioResult
{
  double fireball_max_radius_m  = 0.0;
  double cloud_top_height_m     = 0.0;
  double total_energy_j         = 0.0;
  double seismic_mb             = 0.0;

  // Damage radii (m) for standard thresholds
  std::map<std::string, double> blast_radii;
  std::map<std::string, double> thermal_radii;
  std::map<std::string, double> radiation_radii;

  // Casualty estimates
  double estimated_fatalities   = 0.0;
  double estimated_injuries     = 0.0;
};

// =============================================================================
// Main calculator class
// =============================================================================

/**
 * @brief Generalised nuclear airburst effects calculator
 *
 * Usage:
 * @code
 *   AirburstParameters params;
 *   params.yield_kt = 100.0;
 *   params.burst_height_m = 50.0;
 *
 *   NuclearAirburstEffects calc(params);
 *
 *   // Point query
 *   double P = calc.peakOverpressure(5000.0);  // kPa at 5 km
 *
 *   // Full scenario
 *   AirburstScenarioResult summary = calc.computeScenarioSummary();
 *
 *   // CSV output
 *   calc.writeRadialProfile("profile.csv", 100.0, 50000.0, 500);
 * @endcode
 */
class NuclearAirburstEffects
{
public:
  // ──────────────────────────────────────────────────────────────────────────
  // Construction
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Construct with default parameters
   */
  NuclearAirburstEffects();

  /**
   * @brief Construct from explicit parameter struct
   */
  explicit NuclearAirburstEffects(const AirburstParameters& params);

  /**
   * @brief Construct from INI-style config file
   *
   * Reads [AIRBURST] section for all fields in AirburstParameters.
   */
  explicit NuclearAirburstEffects(const std::string& config_path);

  // ──────────────────────────────────────────────────────────────────────────
  // Parameter access
  // ──────────────────────────────────────────────────────────────────────────

  const AirburstParameters& parameters() const { return params_; }
  void setParameters(const AirburstParameters& p) { params_ = p; }

  // ──────────────────────────────────────────────────────────────────────────
  // A – Blast physics
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Peak static overpressure (kPa) at ground range r (m)
   *
   * Piecewise empirical fit calibrated to Glasstone-Dolan (1977) data
   * with Mach-stem enhancement for low burst angles.
   */
  double peakOverpressure(double ground_range_m) const;

  /**
   * @brief Dynamic (wind) pressure (kPa) from Rankine-Hugoniot
   */
  double dynamicPressure(double overpressure_kpa) const;

  /**
   * @brief Find ground range (m) where overpressure equals target_kpa
   */
  double blastRadiusForPressure(double target_kpa) const;

  // ──────────────────────────────────────────────────────────────────────────
  // B – Thermal radiation
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Thermal fluence (kJ / m^2) at ground range r (m)
   */
  double thermalFluence(double ground_range_m) const;

  /**
   * @brief Maximum fireball radius (m)
   */
  double fireballMaxRadius() const;

  // ──────────────────────────────────────────────────────────────────────────
  // C – Prompt nuclear radiation
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Prompt radiation dose (Gy) at ground range r (m)
   */
  double promptRadiationDose(double ground_range_m) const;

  // ──────────────────────────────────────────────────────────────────────────
  // D – Fallout
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Stabilised cloud top height (m)
   */
  double cloudTopHeight() const;

  /**
   * @brief Normalised fallout activity at (x, y) in metres relative to GZ
   *
   * Uses Gaussian plume model rotated by wind direction.
   */
  double falloutActivity(double x_m, double y_m) const;

  /**
   * @brief Dose rate (Sv/h) at H+1 hour from normalised activity
   */
  double falloutDoseRate1hr(double activity_norm) const;

  // ──────────────────────────────────────────────────────────────────────────
  // E – EMP
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Peak E1 electric field (V/m) at distance (m)
   */
  double empE1Peak(double distance_m) const;

  /**
   * @brief E1 waveform (V/m) at time t (s) and distance (m)
   */
  double empE1Waveform(double t_s, double distance_m) const;

  // ──────────────────────────────────────────────────────────────────────────
  // F – Fireball & shock evolution
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Fireball radius (m) at time t (s)
   */
  double fireballRadius(double t_s) const;

  /**
   * @brief Sedov-Taylor blast shock radius (m) at time t (s)
   */
  double shockRadius(double t_s) const;

  // ──────────────────────────────────────────────────────────────────────────
  // G – Ground-coupled seismics
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Peak ground velocity (m/s) from air-blast coupling
   */
  double peakGroundVelocity(double ground_range_m) const;

  /**
   * @brief Approximate body-wave magnitude from air-ground coupling
   */
  double seismicMagnitude() const;

  // ──────────────────────────────────────────────────────────────────────────
  // H – Coordinate utilities
  // ──────────────────────────────────────────────────────────────────────────

  double metersToDegreesLat(double m) const;
  double metersToDegreesLon(double m) const;

  // ──────────────────────────────────────────────────────────────────────────
  // I – Aggregate queries
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Evaluate all effects at a single ground range
   */
  AirburstPointResult evaluateAtRange(double ground_range_m) const;

  /**
   * @brief Evaluate all effects at (x, y) metres from GZ (includes fallout)
   */
  AirburstPointResult evaluateAtPoint(double x_m, double y_m) const;

  /**
   * @brief Compute scenario-wide summary (radii, magnitudes, casualties)
   */
  AirburstScenarioResult computeScenarioSummary() const;

  // ──────────────────────────────────────────────────────────────────────────
  // J – Output
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Write radial profile to CSV
   */
  void writeRadialProfile(const std::string& path,
                          double r_min_m, double r_max_m, int n_points) const;

  /**
   * @brief Write 2-D grid results to CSV (for external plotting)
   *
   * Columns: x_m, y_m, lon, lat, overpressure_kpa, thermal_kjm2,
   *          prompt_dose_gy, fallout_activity, pgv_cms
   */
  void writeGridCSV(const std::string& path,
                    double half_extent_m, int n_points) const;

  /**
   * @brief Print scenario summary to a stream
   */
  void printSummary(std::ostream& os = std::cout) const;

  // ──────────────────────────────────────────────────────────────────────────
  // K – Config I/O
  // ──────────────────────────────────────────────────────────────────────────

  /**
   * @brief Load parameters from INI-style config file
   */
  bool loadConfig(const std::string& path);

  /**
   * @brief Write current parameters as a template config file
   */
  void writeTemplateConfig(const std::string& path) const;

  // ──────────────────────────────────────────────────────────────────────────
  // Standard damage / threshold tables
  // ──────────────────────────────────────────────────────────────────────────

  static std::vector<DamageThreshold> blastDamageThresholds();
  static std::vector<DamageThreshold> thermalDamageThresholds();
  static std::vector<DamageThreshold> radiationDamageThresholds();

private:
  AirburstParameters params_;

  // Helper: slant range from ground range
  double slantRange(double ground_range_m) const;
};

} // namespace FSRM

#endif // NUCLEAR_AIRBURST_EFFECTS_HPP
