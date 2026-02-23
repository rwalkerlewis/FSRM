/**
 * @file NuclearAirburstEffects.cpp
 * @brief Parameterised nuclear airburst effects — full implementation
 *
 * Mirrors every physics formula from the Python model_berkeley_100kt_airburst.py
 * in a generic, config-driven C++ class.
 *
 * All calculations use SI base units internally:
 *   distance  → metres
 *   pressure  → kPa (blast) or Pa where noted
 *   energy    → Joules
 *   dose      → Gray
 *   field     → V/m
 *
 * References:
 *   Glasstone & Dolan (1977) — Effects of Nuclear Weapons
 *   Brode (1955)             — Numerical Solutions of Spherical Blast Waves
 *   Sedov (1959)             — Similarity and Dimensional Methods in Mechanics
 *   Way & Wigner (1948)      — Radioactive Decay of Fission Products
 */

#include "NuclearAirburstEffects.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace FSRM
{

namespace C = AirburstConstants;

// =============================================================================
// Helpers
// =============================================================================

static inline double clamp(double v, double lo, double hi)
{
  return std::max(lo, std::min(v, hi));
}

// Trim whitespace from both ends of a string
static std::string trim(const std::string& s)
{
  auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

// =============================================================================
// Construction
// =============================================================================

NuclearAirburstEffects::NuclearAirburstEffects()
  : params_()
{
}

NuclearAirburstEffects::NuclearAirburstEffects(const AirburstParameters& params)
  : params_(params)
{
}

NuclearAirburstEffects::NuclearAirburstEffects(const std::string& config_path)
  : params_()
{
  if (!loadConfig(config_path))
  {
    std::cerr << "NuclearAirburstEffects: failed to load config from "
              << config_path << ", using defaults\n";
  }
}

// =============================================================================
// Private helper
// =============================================================================

double NuclearAirburstEffects::slantRange(double ground_range_m) const
{
  return std::sqrt(ground_range_m * ground_range_m
                   + params_.burst_height_m * params_.burst_height_m);
}

// =============================================================================
// A — Blast physics
// =============================================================================

double NuclearAirburstEffects::peakOverpressure(double ground_range_m) const
{
  const double r   = std::max(ground_range_m, 1.0);
  const double R   = slantRange(r);
  const double W   = std::max(params_.yield_kt, 1e-12);
  const double W_kg = W * 1.0e6;   // kt → kg TNT equivalent
  const double Zbar = R / std::cbrt(W_kg);  // scaled distance (m / kg^{1/3})

  // Kinney-Graham (1985) free-air peak overpressure — continuous closed-form
  // fit to compiled experimental data:
  //   ΔP/P₀ = 808 [1 + (Z̄/4.5)²]
  //           / {√[1+(Z̄/0.048)²] · √[1+(Z̄/0.32)²] · √[1+(Z̄/1.35)²]}
  //
  // where Z̄ = R / m_charge^{1/3} (metres per kg^{1/3}).
  // Reference: Kinney & Graham, "Explosive Shocks in Air", 2nd ed., 1985.
  const double P0 = C::STANDARD_PRESSURE;  // 101.325 kPa
  const double num = 808.0 * (1.0 + std::pow(Zbar / 4.5, 2));
  const double d1  = std::sqrt(1.0 + std::pow(Zbar / 0.048, 2));
  const double d2  = std::sqrt(1.0 + std::pow(Zbar / 0.32,  2));
  const double d3  = std::sqrt(1.0 + std::pow(Zbar / 1.35,  2));
  double P_kpa = P0 * num / (d1 * d2 * d3);

  // Mach stem enhancement for low burst angles (< 40 deg from horizontal)
  const double angle = std::atan2(params_.burst_height_m, std::max(r, 1.0));
  const double mach_factor = (angle < 40.0 * C::DEG_TO_RAD) ? 1.8 : 1.0;

  return P_kpa * mach_factor;
}

double NuclearAirburstEffects::dynamicPressure(double overpressure_kpa) const
{
  // Rankine-Hugoniot: q = 2.5 ΔP² / (7 P₀ + ΔP)
  const double P  = std::max(overpressure_kpa, 0.0);
  const double P0 = C::STANDARD_PRESSURE;
  return 2.5 * P * P / (7.0 * P0 + P);
}

double NuclearAirburstEffects::blastRadiusForPressure(double target_kpa) const
{
  // Simple numerical search (bisection over 1–200 km)
  double lo = 1.0, hi = 200000.0;
  for (int i = 0; i < 60; ++i)
  {
    double mid = 0.5 * (lo + hi);
    double P   = peakOverpressure(mid);
    if (P > target_kpa)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

// =============================================================================
// B — Thermal radiation
// =============================================================================

double NuclearAirburstEffects::thermalFluence(double ground_range_m) const
{
  const double r  = std::max(ground_range_m, 1.0);
  const double R  = slantRange(r);
  const double Et = 0.35 * params_.yield_kt * C::JOULES_PER_KT;  // 35 % of yield
  const double Q  = Et / (4.0 * C::PI * R * R);                  // J / m^2
  const double tau = std::exp(-R / params_.thermal_attenuation_length_m);
  return Q * tau / 1000.0;  // kJ / m^2
}

double NuclearAirburstEffects::fireballMaxRadius() const
{
  // Glasstone & Dolan: R_max = 66 W^{0.4}  (metres, W in kt)
  return 66.0 * std::pow(std::max(params_.yield_kt, 1e-12), 0.4);
}

// =============================================================================
// C — Prompt nuclear radiation
// =============================================================================

double NuclearAirburstEffects::promptRadiationDose(double ground_range_m) const
{
  const double r    = std::max(ground_range_m, 1.0);
  const double R    = slantRange(r);
  const double Dref = params_.prompt_dose_ref_gy;
  const double W    = params_.yield_kt;
  const double lam  = params_.prompt_atten_length_m;
  return Dref * W * std::pow(1000.0 / R, 2.0) * std::exp(-R / lam);
}

// =============================================================================
// D — Fallout
// =============================================================================

double NuclearAirburstEffects::cloudTopHeight() const
{
  // GD 1977: H ≈ 2200 W^{0.4}, capped at 25 km
  return std::min(25000.0, 2200.0 * std::pow(std::max(params_.yield_kt, 1e-12), 0.4));
}

double NuclearAirburstEffects::falloutActivity(double x_m, double y_m) const
{
  const double theta = params_.wind_direction_deg * C::DEG_TO_RAD;
  const double dx_w  = std::cos(theta);
  const double dy_w  = std::sin(theta);

  const double downwind  = x_m * dx_w + y_m * dy_w;
  const double crosswind = std::abs(-x_m * dy_w + y_m * dx_w);

  if (downwind <= 0.0) return 0.0;

  const double H     = cloudTopHeight();
  const double t_fall = H / 2.0;
  const double drift  = params_.wind_speed_altitude_mps * t_fall;

  const double sig_a = 0.08 * std::max(downwind, 1.0) + 800.0;
  const double sig_c = 0.05 * std::max(downwind, 1.0) + 500.0;

  const double g1 = (downwind - drift) / sig_a;
  const double g2 = crosswind / sig_c;
  return std::exp(-0.5 * g1 * g1) * std::exp(-0.5 * g2 * g2);
}

double NuclearAirburstEffects::falloutDoseRate1hr(double activity_norm) const
{
  const double fission_yield_kt = params_.yield_kt * params_.fission_fraction;
  const double ref_rate = 50.0 * fission_yield_kt;  // Sv/h at hotspot
  return activity_norm * ref_rate;
}

// =============================================================================
// E — EMP
// =============================================================================

double NuclearAirburstEffects::empE1Peak(double distance_m) const
{
  const double E0 = params_.emp_e0_vm * std::sqrt(params_.yield_kt / 100.0);
  return E0 * std::exp(-distance_m / params_.emp_decay_length_m);
}

double NuclearAirburstEffects::empE1Waveform(double t_s, double distance_m) const
{
  if (t_s <= 0.0) return 0.0;
  const double Ep = empE1Peak(distance_m);
  return Ep * (std::exp(-t_s / params_.emp_tau_decay_s)
               - std::exp(-t_s / params_.emp_tau_rise_s));
}

// =============================================================================
// F — Fireball & shock evolution
// =============================================================================

double NuclearAirburstEffects::fireballRadius(double t_s) const
{
  if (t_s <= 0.0) return 0.0;
  const double Rmax   = fireballMaxRadius();
  const double t_form = 0.001 * std::pow(std::max(params_.yield_kt, 1e-12), 0.3);
  if (t_s < t_form)
    return Rmax * std::pow(t_s / t_form, 0.4);
  return Rmax;
}

double NuclearAirburstEffects::shockRadius(double t_s) const
{
  // Sedov-Taylor: R = 1.15 (E / ρ₀)^{0.2} t^{0.4}
  if (t_s <= 0.0) return 0.0;
  const double E   = params_.yield_kt * C::JOULES_PER_KT;
  const double rho = C::SEA_LEVEL_DENSITY;
  return 1.15 * std::pow(E / rho, 0.2) * std::pow(t_s, 0.4);
}

// =============================================================================
// G — Ground-coupled seismics
// =============================================================================

double NuclearAirburstEffects::peakGroundVelocity(double ground_range_m) const
{
  const double P = peakOverpressure(ground_range_m);  // kPa
  const double impedance = params_.ground_density_kgm3 * params_.ground_vs_mps;
  // P*1000 converts kPa → Pa;  factor 0.5 accounts for coupling efficiency
  return P * 1000.0 / impedance * 0.5;
}

double NuclearAirburstEffects::seismicMagnitude() const
{
  const double Weff = params_.yield_kt * params_.seismic_coupling_eff;
  return 4.0 + 0.75 * std::log10(std::max(Weff, 1e-6));
}

// =============================================================================
// H — Coordinate utilities
// =============================================================================

double NuclearAirburstEffects::metersToDegreesLat(double m) const
{
  return m / 111320.0;
}

double NuclearAirburstEffects::metersToDegreesLon(double m) const
{
  return m / (111320.0 * std::cos(params_.latitude * C::DEG_TO_RAD));
}

// =============================================================================
// I — Point / grid evaluation
// =============================================================================

AirburstPointResult NuclearAirburstEffects::evaluateAtRange(double ground_range_m) const
{
  AirburstPointResult res;
  res.ground_range_m          = ground_range_m;
  res.slant_range_m           = slantRange(ground_range_m);
  res.peak_overpressure_kpa   = peakOverpressure(ground_range_m);
  res.dynamic_pressure_kpa    = dynamicPressure(res.peak_overpressure_kpa);
  res.thermal_fluence_kjm2    = thermalFluence(ground_range_m);
  res.prompt_dose_gy          = promptRadiationDose(ground_range_m);
  res.emp_e1_peak_vm          = empE1Peak(ground_range_m);
  res.pgv_ms                  = peakGroundVelocity(ground_range_m);
  res.pgv_cms                 = res.pgv_ms * 100.0;
  return res;
}

AirburstPointResult NuclearAirburstEffects::evaluateAtPoint(double x_m, double y_m) const
{
  const double r = std::sqrt(x_m * x_m + y_m * y_m);
  AirburstPointResult res = evaluateAtRange(r);
  res.fallout_activity_norm = falloutActivity(x_m, y_m);
  res.fallout_dose_rate_svh = falloutDoseRate1hr(res.fallout_activity_norm);
  return res;
}

// =============================================================================
// Scenario summary
// =============================================================================

AirburstScenarioResult NuclearAirburstEffects::computeScenarioSummary() const
{
  AirburstScenarioResult res;
  res.fireball_max_radius_m = fireballMaxRadius();
  res.cloud_top_height_m    = cloudTopHeight();
  res.total_energy_j        = params_.yield_kt * C::JOULES_PER_KT;
  res.seismic_mb            = seismicMagnitude();

  // Blast radii
  for (const auto& th : blastDamageThresholds())
  {
    res.blast_radii[th.label] = blastRadiusForPressure(th.value);
  }

  // Thermal radii (numerical search)
  for (const auto& th : thermalDamageThresholds())
  {
    double rr = 10.0;
    while (rr < 60000.0 && thermalFluence(rr) > th.value) rr += 50.0;
    res.thermal_radii[th.label] = rr;
  }

  // Radiation radii
  for (const auto& th : radiationDamageThresholds())
  {
    double rr = 10.0;
    while (rr < 60000.0 && promptRadiationDose(rr) > th.value) rr += 50.0;
    res.radiation_radii[th.label] = rr;
  }

  // Casualty estimation (annular ring model)
  const double ring_edges[] = {0, 500, 1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000};
  const double fatality_rate[] = {1.0, 0.98, 0.90, 0.60, 0.35, 0.15, 0.05, 0.02, 0.005, 0.001};
  const double injury_rate[]   = {0.0, 0.02, 0.10, 0.35, 0.50, 0.45, 0.30, 0.15, 0.05, 0.01};
  constexpr int n_rings = 10;

  res.estimated_fatalities = 0.0;
  res.estimated_injuries   = 0.0;
  for (int i = 0; i < n_rings; ++i)
  {
    double area_km2 = C::PI * (std::pow(ring_edges[i + 1] / 1000.0, 2)
                                - std::pow(ring_edges[i] / 1000.0, 2));
    double pop = area_km2 * params_.population_density_urban;
    res.estimated_fatalities += pop * fatality_rate[i];
    res.estimated_injuries   += pop * injury_rate[i];
  }

  return res;
}

// =============================================================================
// Output
// =============================================================================

void NuclearAirburstEffects::writeRadialProfile(const std::string& path,
                                                 double r_min_m, double r_max_m,
                                                 int n_points) const
{
  std::ofstream f(path);
  if (!f.is_open())
  {
    std::cerr << "Cannot open " << path << " for writing\n";
    return;
  }

  f << "range_m,slant_m,overpressure_kpa,dynamic_kpa,thermal_kjm2,"
       "prompt_gy,emp_e1_vm,pgv_ms,pgv_cms\n";

  const double dr = (r_max_m - r_min_m) / std::max(n_points - 1, 1);
  for (int i = 0; i < n_points; ++i)
  {
    double r = r_min_m + i * dr;
    auto pt = evaluateAtRange(r);
    f << std::fixed << std::setprecision(2) << r << ","
      << pt.slant_range_m << ","
      << std::scientific << std::setprecision(6)
      << pt.peak_overpressure_kpa << ","
      << pt.dynamic_pressure_kpa << ","
      << pt.thermal_fluence_kjm2 << ","
      << pt.prompt_dose_gy << ","
      << pt.emp_e1_peak_vm << ","
      << pt.pgv_ms << ","
      << pt.pgv_cms << "\n";
  }
  f.close();
  std::cout << "  Wrote radial profile: " << path << "  (" << n_points << " points)\n";
}

void NuclearAirburstEffects::writeGridCSV(const std::string& path,
                                           double half_extent_m, int n_points) const
{
  std::ofstream f(path);
  if (!f.is_open())
  {
    std::cerr << "Cannot open " << path << " for writing\n";
    return;
  }

  f << "x_m,y_m,lon,lat,overpressure_kpa,dynamic_kpa,thermal_kjm2,"
       "prompt_gy,fallout_activity,fallout_dose_svh,pgv_cms\n";

  const double dx = 2.0 * half_extent_m / std::max(n_points - 1, 1);
  for (int j = 0; j < n_points; ++j)
  {
    double y = -half_extent_m + j * dx;
    for (int i = 0; i < n_points; ++i)
    {
      double x = -half_extent_m + i * dx;
      auto pt = evaluateAtPoint(x, y);
      double lon = params_.longitude + metersToDegreesLon(x);
      double lat = params_.latitude  + metersToDegreesLat(y);
      f << std::fixed << std::setprecision(1) << x << "," << y << ","
        << std::setprecision(6) << lon << "," << lat << ","
        << std::scientific << std::setprecision(4)
        << pt.peak_overpressure_kpa << ","
        << pt.dynamic_pressure_kpa << ","
        << pt.thermal_fluence_kjm2 << ","
        << pt.prompt_dose_gy << ","
        << pt.fallout_activity_norm << ","
        << pt.fallout_dose_rate_svh << ","
        << pt.pgv_cms << "\n";
    }
  }
  f.close();
  std::cout << "  Wrote 2D grid: " << path << "  ("
            << n_points << " x " << n_points << ")\n";
}

void NuclearAirburstEffects::printSummary(std::ostream& os) const
{
  auto summary = computeScenarioSummary();

  os << std::string(72, '=') << "\n"
     << "  Nuclear Airburst Effects Summary\n"
     << std::string(72, '=') << "\n\n"
     << "  Yield:            " << params_.yield_kt << " kt\n"
     << "  Height of burst:  " << params_.burst_height_m << " m AGL\n"
     << "  Location:         " << params_.latitude << " N, "
                               << params_.longitude << " E\n"
     << "                    " << params_.location_name << "\n\n"
     << "  Total energy:     " << std::scientific << std::setprecision(3)
                               << summary.total_energy_j << " J\n"
     << "  Fireball Rmax:    " << std::fixed << std::setprecision(0)
                               << summary.fireball_max_radius_m << " m\n"
     << "  Cloud top:        " << summary.cloud_top_height_m / 1000.0 << " km\n"
     << "  Seismic mb:       " << std::setprecision(1) << summary.seismic_mb << "\n\n";

  os << "  BLAST DAMAGE RADII:\n";
  for (const auto& [label, radius] : summary.blast_radii)
  {
    os << "    " << std::left << std::setw(46) << label
       << std::right << std::setprecision(2) << radius / 1000.0 << " km\n";
  }

  os << "\n  THERMAL:\n";
  for (const auto& [label, radius] : summary.thermal_radii)
  {
    os << "    " << std::left << std::setw(46) << label
       << std::right << std::setprecision(2) << radius / 1000.0 << " km\n";
  }

  os << "\n  PROMPT RADIATION:\n";
  for (const auto& [label, radius] : summary.radiation_radii)
  {
    os << "    " << std::left << std::setw(46) << label
       << std::right << std::setprecision(2) << radius / 1000.0 << " km\n";
  }

  os << "\n  PGV @ 5 km:       " << std::setprecision(2)
     << peakGroundVelocity(5000.0) * 100.0 << " cm/s\n"
     << "  PGV @ 10 km:      "
     << peakGroundVelocity(10000.0) * 100.0 << " cm/s\n\n"
     << "  ESTIMATED CASUALTIES (flat terrain, "
     << params_.population_density_urban << " ppl/km²):\n"
     << "    Fatalities:     ~" << std::fixed << std::setprecision(0)
     << summary.estimated_fatalities << "\n"
     << "    Injuries:       ~" << summary.estimated_injuries << "\n\n"
     << std::string(72, '=') << "\n";
}

// =============================================================================
// Config I/O
// =============================================================================

bool NuclearAirburstEffects::loadConfig(const std::string& path)
{
  std::ifstream file(path);
  if (!file.is_open()) return false;

  std::string section;
  std::string line;

  while (std::getline(file, line))
  {
    line = trim(line);
    if (line.empty() || line[0] == '#' || line[0] == ';') continue;

    // Section header
    if (line.front() == '[' && line.back() == ']')
    {
      section = trim(line.substr(1, line.size() - 2));
      continue;
    }

    auto eq = line.find('=');
    if (eq == std::string::npos) continue;

    std::string key = trim(line.substr(0, eq));
    std::string val = trim(line.substr(eq + 1));
    // Strip inline comments
    auto cpos = val.find('#');
    if (cpos != std::string::npos) val = trim(val.substr(0, cpos));

    // Map keys from [AIRBURST] section to params_
    if (section == "AIRBURST" || section == "NUCLEAR_AIRBURST"
        || section == "DEVICE" || section == "NUCLEAR_DEVICE")
    {
      try
      {
        if      (key == "yield_kt")         params_.yield_kt         = std::stod(val);
        else if (key == "burst_height_m")    params_.burst_height_m   = std::stod(val);
        else if (key == "fission_fraction")  params_.fission_fraction = std::stod(val);
        else if (key == "latitude")          params_.latitude         = std::stod(val);
        else if (key == "longitude")         params_.longitude        = std::stod(val);
        else if (key == "location_name")     params_.location_name    = val;
        else if (key == "thermal_attenuation_length_m")
          params_.thermal_attenuation_length_m = std::stod(val);
      }
      catch (...)
      {
        std::cerr << "Warning: cannot parse [" << section << "] " << key << " = " << val << "\n";
      }
    }
    else if (section == "WIND")
    {
      try
      {
        if      (key == "speed_surface_mps")  params_.wind_speed_surface_mps  = std::stod(val);
        else if (key == "speed_altitude_mps") params_.wind_speed_altitude_mps = std::stod(val);
        else if (key == "direction_deg")      params_.wind_direction_deg      = std::stod(val);
      }
      catch (...) {}
    }
    else if (section == "GROUND" || section == "SITE")
    {
      try
      {
        if      (key == "density_kgm3")       params_.ground_density_kgm3  = std::stod(val);
        else if (key == "vs_mps")             params_.ground_vs_mps        = std::stod(val);
        else if (key == "seismic_coupling")   params_.seismic_coupling_eff = std::stod(val);
      }
      catch (...) {}
    }
    else if (section == "POPULATION")
    {
      try
      {
        if      (key == "density_urban")      params_.population_density_urban = std::stod(val);
        else if (key == "density_suburban")    params_.population_density_sub   = std::stod(val);
      }
      catch (...) {}
    }
    else if (section == "RADIATION")
    {
      try
      {
        if      (key == "prompt_dose_ref_gy")   params_.prompt_dose_ref_gy    = std::stod(val);
        else if (key == "prompt_atten_length_m") params_.prompt_atten_length_m = std::stod(val);
      }
      catch (...) {}
    }
    else if (section == "EMP")
    {
      try
      {
        if      (key == "e0_vm")           params_.emp_e0_vm          = std::stod(val);
        else if (key == "decay_length_m")  params_.emp_decay_length_m = std::stod(val);
        else if (key == "tau_rise_s")      params_.emp_tau_rise_s     = std::stod(val);
        else if (key == "tau_decay_s")     params_.emp_tau_decay_s    = std::stod(val);
      }
      catch (...) {}
    }
  }
  return true;
}

void NuclearAirburstEffects::writeTemplateConfig(const std::string& path) const
{
  std::ofstream f(path);
  if (!f.is_open())
  {
    std::cerr << "Cannot write template to " << path << "\n";
    return;
  }

  f << "# ==============================================================================\n"
    << "# Nuclear Airburst Effects — Configuration File\n"
    << "# ==============================================================================\n"
    << "# Edit the parameters below to define your airburst scenario.\n"
    << "# All internal calculations use SI units.\n"
    << "# Lines beginning with # are comments.\n"
    << "#\n"
    << "# Usage:\n"
    << "#   nuclear_airburst_effects -c this_file.config\n"
    << "# ==============================================================================\n\n"

    << "[AIRBURST]\n"
    << "yield_kt          = " << params_.yield_kt         << "    # Weapon yield (kt TNT)\n"
    << "burst_height_m    = " << params_.burst_height_m   << "    # Height of burst AGL (m)\n"
    << "fission_fraction  = " << params_.fission_fraction  << "    # Fraction from fission\n"
    << "latitude          = " << params_.latitude          << "    # Ground zero (deg N)\n"
    << "longitude         = " << params_.longitude         << "    # Ground zero (deg E, negative = W)\n"
    << "location_name     = " << params_.location_name     << "\n"
    << "thermal_attenuation_length_m = " << params_.thermal_attenuation_length_m
                                         << "  # 1/e absorption (m)\n\n"

    << "[WIND]\n"
    << "speed_surface_mps  = " << params_.wind_speed_surface_mps  << "    # Surface wind (m/s)\n"
    << "speed_altitude_mps = " << params_.wind_speed_altitude_mps << "    # Wind at cloud height (m/s)\n"
    << "direction_deg      = " << params_.wind_direction_deg      << "   # From direction (meteorological, deg)\n\n"

    << "[GROUND]\n"
    << "density_kgm3      = " << params_.ground_density_kgm3  << "  # Near-surface density (kg/m^3)\n"
    << "vs_mps            = " << params_.ground_vs_mps         << "  # Shear-wave speed (m/s)\n"
    << "seismic_coupling  = " << params_.seismic_coupling_eff  << " # Air-to-ground coupling efficiency\n\n"

    << "[POPULATION]\n"
    << "density_urban     = " << params_.population_density_urban << "  # Urban (ppl/km^2)\n"
    << "density_suburban  = " << params_.population_density_sub   << "  # Suburban (ppl/km^2)\n\n"

    << "[RADIATION]\n"
    << "prompt_dose_ref_gy    = " << params_.prompt_dose_ref_gy    << "    # Ref dose Gy at 1 km for 1 kt\n"
    << "prompt_atten_length_m = " << params_.prompt_atten_length_m << "  # 1/e length in air (m)\n\n"

    << "[EMP]\n"
    << "e0_vm           = " << params_.emp_e0_vm          << "  # E1 ref peak at 100 kt (V/m)\n"
    << "decay_length_m  = " << params_.emp_decay_length_m << "  # Spatial decay (m)\n"
    << "tau_rise_s      = " << params_.emp_tau_rise_s     << "   # E1 rise time (s)\n"
    << "tau_decay_s     = " << params_.emp_tau_decay_s    << "   # E1 decay time (s)\n";

  f.close();
  std::cout << "  Wrote template config: " << path << "\n";
}

// =============================================================================
// Static threshold tables
// =============================================================================

std::vector<DamageThreshold> NuclearAirburstEffects::blastDamageThresholds()
{
  return {
    {"Total destruction (>140 kPa)",             140.0},
    {"Reinforced concrete destroyed (>35 kPa)",   35.0},
    {"Most buildings destroyed (>14 kPa)",        14.0},
    {"Moderate damage (>7 kPa)",                   7.0},
    {"Light damage, window breakage (>3.5 kPa)",   3.5},
    {"Glass breakage (>1 kPa)",                    1.0},
  };
}

std::vector<DamageThreshold> NuclearAirburstEffects::thermalDamageThresholds()
{
  return {
    {"Third-degree burns (>670 kJ/m²)",  670.0},
    {"Second-degree burns (>335 kJ/m²)", 335.0},
    {"First-degree burns (>125 kJ/m²)",  125.0},
    {"Pain threshold (~50 kJ/m²)",        50.0},
  };
}

std::vector<DamageThreshold> NuclearAirburstEffects::radiationDamageThresholds()
{
  return {
    {"Fatal (>10 Gy)",                     10.0},
    {"LD50/60 (~4.5 Gy)",                  4.5},
    {"Acute radiation sickness (>2 Gy)",    2.0},
    {"Mild symptoms (>0.5 Gy)",             0.5},
    {"Detectable (>0.01 Gy)",               0.01},
  };
}

} // namespace FSRM
