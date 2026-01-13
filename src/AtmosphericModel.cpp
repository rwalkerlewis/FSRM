/**
 * @file AtmosphericModel.cpp
 * @brief Implementation of comprehensive atmospheric models
 */

#include "AtmosphericModel.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace FSRM {

// =============================================================================
// AtmosphereUtils Implementation
// =============================================================================

namespace AtmosphereUtils {

double geometricToGeopotential(double geometric_alt) {
    const double r = AtmosphereConstants::EARTH_RADIUS;
    return (r * geometric_alt) / (r + geometric_alt);
}

double geopotentialToGeometric(double geopotential_alt) {
    const double r = AtmosphereConstants::EARTH_RADIUS;
    return (r * geopotential_alt) / (r - geopotential_alt);
}

double gravity(double altitude, double latitude) {
    // Gravity variation with altitude and latitude
    const double r = AtmosphereConstants::EARTH_RADIUS;
    
    // Latitude correction (WGS84 approximation)
    double sin_lat = std::sin(latitude * M_PI / 180.0);
    double g_lat = 9.7803253359 * (1.0 + 0.00193185265241 * sin_lat * sin_lat) /
                   std::sqrt(1.0 - 0.00669437999014 * sin_lat * sin_lat);
    
    // Altitude correction (free-air)
    double g = g_lat * std::pow(r / (r + altitude), 2);
    
    return g;
}

double barometricFormula(double altitude, double T0, double P0, double lapse_rate) {
    const double g = AtmosphereConstants::G0;
    const double R = AtmosphereConstants::R_AIR_DRY;
    
    if (std::abs(lapse_rate) < 1e-10) {
        // Isothermal layer
        return P0 * std::exp(-g * altitude / (R * T0));
    } else {
        // Non-isothermal layer
        double T = T0 + lapse_rate * altitude;
        return P0 * std::pow(T / T0, -g / (lapse_rate * R));
    }
}

double hypsometricEquation(double P1, double P2, double T_mean) {
    const double R = AtmosphereConstants::R_AIR_DRY;
    const double g = AtmosphereConstants::G0;
    return (R * T_mean / g) * std::log(P1 / P2);
}

double dynamicViscosity(double temperature) {
    // Sutherland's law
    const double mu0 = 1.716e-5;  // Pa·s at T0
    const double T0 = 273.15;     // K
    const double S = 110.4;       // Sutherland constant for air
    
    return mu0 * std::pow(temperature / T0, 1.5) * (T0 + S) / (temperature + S);
}

double thermalConductivity(double temperature) {
    // Sutherland's law for thermal conductivity
    const double k0 = 0.0241;     // W/(m·K) at T0
    const double T0 = 273.15;     // K
    const double S = 194.0;       // Sutherland constant
    
    return k0 * std::pow(temperature / T0, 1.5) * (T0 + S) / (temperature + S);
}

double diffusivity(double temperature, double pressure) {
    // Binary diffusion coefficient (air self-diffusion approximation)
    const double D0 = 1.8e-5;     // m²/s at STP
    const double P0 = 101325.0;   // Pa
    const double T0 = 273.15;     // K
    
    return D0 * std::pow(temperature / T0, 1.81) * (P0 / pressure);
}

double soundSpeed(double temperature, double gamma) {
    const double R = AtmosphereConstants::R_AIR_DRY;
    return std::sqrt(gamma * R * temperature);
}

double saturationVaporPressure(double temperature) {
    // Tetens formula (valid for water, T > 0°C)
    double T_C = temperature - 273.15;
    if (T_C >= 0.0) {
        return 610.78 * std::exp(17.27 * T_C / (T_C + 237.3));
    } else {
        // Use ice formula for T < 0°C
        return saturationVaporPressureIce(temperature);
    }
}

double saturationVaporPressureIce(double temperature) {
    // Tetens formula for ice
    double T_C = temperature - 273.15;
    return 610.78 * std::exp(21.875 * T_C / (T_C + 265.5));
}

double virtualTemperature(double temperature, double mixing_ratio) {
    // Virtual temperature accounts for moisture
    return temperature * (1.0 + 0.608 * mixing_ratio);
}

double potentialTemperature(double temperature, double pressure) {
    const double P0 = AtmosphereConstants::P0_STD;
    const double kappa = AtmosphereConstants::R_AIR_DRY / AtmosphereConstants::CP_AIR;
    return temperature * std::pow(P0 / pressure, kappa);
}

double equivalentPotentialTemperature(double temperature, double pressure, double mixing_ratio) {
    double theta = potentialTemperature(temperature, pressure);
    const double Lv = AtmosphereConstants::LATENT_HEAT_VAP;
    const double cp = AtmosphereConstants::CP_AIR;
    return theta * std::exp(Lv * mixing_ratio / (cp * temperature));
}

double bruntVaisalaFrequency(double dT_dz, double T, double g) {
    // N² = (g/T) * (dT/dz + g/cp)
    const double cp = AtmosphereConstants::CP_AIR;
    double N2 = (g / T) * (dT_dz + g / cp);
    if (N2 < 0) return -std::sqrt(-N2);  // Unstable
    return std::sqrt(N2);
}

double richardsonNumber(double N2, double dU_dz) {
    if (std::abs(dU_dz) < 1e-10) return 1e10;  // Very stable
    return N2 / (dU_dz * dU_dz);
}

char pasquillStabilityClass(double wind_speed, double solar_radiation, bool night) {
    // Simplified Pasquill-Gifford stability classification
    if (night) {
        if (wind_speed < 2.0) return 'F';
        if (wind_speed < 3.0) return 'E';
        return 'D';
    } else {
        if (solar_radiation > 700) {
            if (wind_speed < 2.0) return 'A';
            if (wind_speed < 3.0) return 'A';
            if (wind_speed < 5.0) return 'B';
            return 'C';
        } else if (solar_radiation > 350) {
            if (wind_speed < 2.0) return 'A';
            if (wind_speed < 3.0) return 'B';
            if (wind_speed < 5.0) return 'B';
            return 'C';
        } else if (solar_radiation > 0) {
            if (wind_speed < 2.0) return 'B';
            if (wind_speed < 5.0) return 'C';
            return 'D';
        }
        return 'D';
    }
}

} // namespace AtmosphereUtils

// =============================================================================
// AtmosphericSpecies Implementation
// =============================================================================

double AtmosphericSpecies::mixingRatio(double altitude) const {
    switch (profile_type) {
        case ProfileType::CONSTANT:
            return constant_mixing_ratio;
            
        case ProfileType::EXPONENTIAL:
            return constant_mixing_ratio * std::exp(-altitude / scale_height);
            
        case ProfileType::TABULATED: {
            if (altitudes.empty()) return 0.0;
            if (altitude <= altitudes.front()) return values.front();
            if (altitude >= altitudes.back()) return values.back();
            
            // Linear interpolation
            auto it = std::lower_bound(altitudes.begin(), altitudes.end(), altitude);
            size_t i = it - altitudes.begin();
            if (i == 0) i = 1;
            double t = (altitude - altitudes[i-1]) / (altitudes[i] - altitudes[i-1]);
            return values[i-1] + t * (values[i] - values[i-1]);
        }
            
        case ProfileType::FUNCTION:
            if (profile_function) return profile_function(altitude);
            return 0.0;
    }
    return 0.0;
}

double AtmosphericSpecies::numberDensity(double altitude, double total_number_density) const {
    return mixingRatio(altitude) * total_number_density;
}

double AtmosphericSpecies::massDensity(double altitude, double total_density) const {
    double mr = mixingRatio(altitude);
    return mr * molecular_weight / AtmosphereConstants::M_AIR_DRY * total_density;
}

// =============================================================================
// WindProfile Implementation
// =============================================================================

void WindProfile::getWind(double altitude, double& u, double& v, double& w) const {
    w = 0.0;  // Vertical wind usually zero in background profile
    
    switch (type) {
        case Type::NONE:
            u = v = 0.0;
            return;
            
        case Type::CONSTANT: {
            double dir_rad = reference_direction * M_PI / 180.0;
            u = -reference_speed * std::sin(dir_rad);
            v = -reference_speed * std::cos(dir_rad);
            return;
        }
            
        case Type::LOGARITHMIC: {
            if (altitude <= displacement_height + roughness_length) {
                u = v = 0.0;
                return;
            }
            double speed = (friction_velocity / 0.4) * 
                          std::log((altitude - displacement_height) / roughness_length);
            double dir_rad = reference_direction * M_PI / 180.0;
            u = -speed * std::sin(dir_rad);
            v = -speed * std::cos(dir_rad);
            return;
        }
            
        case Type::POWER_LAW: {
            double speed = reference_speed * 
                          std::pow(altitude / reference_height, power_exponent);
            double dir_rad = reference_direction * M_PI / 180.0;
            u = -speed * std::sin(dir_rad);
            v = -speed * std::cos(dir_rad);
            return;
        }
            
        case Type::EKMAN_SPIRAL: {
            double z_norm = altitude / boundary_layer_height;
            double gamma = M_PI / 4.0;  // Ekman parameter
            double exp_factor = std::exp(-gamma * z_norm);
            double speed = geostrophic_wind * std::sqrt(
                1.0 - 2.0 * exp_factor * std::cos(gamma * z_norm) + exp_factor * exp_factor);
            double angle_deviation = std::atan2(
                exp_factor * std::sin(gamma * z_norm),
                1.0 - exp_factor * std::cos(gamma * z_norm));
            double dir_rad = (reference_direction - angle_deviation * 180.0 / M_PI) * M_PI / 180.0;
            u = -speed * std::sin(dir_rad);
            v = -speed * std::cos(dir_rad);
            return;
        }
            
        case Type::TABULATED: {
            if (altitudes.empty()) {
                u = v = 0.0;
                return;
            }
            double speed, direction;
            if (altitude <= altitudes.front()) {
                speed = speeds.front();
                direction = directions.front();
            } else if (altitude >= altitudes.back()) {
                speed = speeds.back();
                direction = directions.back();
            } else {
                auto it = std::lower_bound(altitudes.begin(), altitudes.end(), altitude);
                size_t i = it - altitudes.begin();
                if (i == 0) i = 1;
                double t = (altitude - altitudes[i-1]) / (altitudes[i] - altitudes[i-1]);
                speed = speeds[i-1] + t * (speeds[i] - speeds[i-1]);
                direction = directions[i-1] + t * (directions[i] - directions[i-1]);
            }
            double dir_rad = direction * M_PI / 180.0;
            u = -speed * std::sin(dir_rad);
            v = -speed * std::cos(dir_rad);
            return;
        }
            
        case Type::FUNCTION:
            if (wind_function) {
                wind_function(altitude, u, v, w);
            } else {
                u = v = w = 0.0;
            }
            return;
    }
}

double WindProfile::getSpeed(double altitude) const {
    double u, v, w;
    getWind(altitude, u, v, w);
    return std::sqrt(u*u + v*v);
}

double WindProfile::getDirection(double altitude) const {
    double u, v, w;
    getWind(altitude, u, v, w);
    double dir = std::atan2(-u, -v) * 180.0 / M_PI;
    if (dir < 0) dir += 360.0;
    return dir;
}

// =============================================================================
// HumidityProfile Implementation
// =============================================================================

double HumidityProfile::getRelativeHumidity(double altitude) const {
    switch (type) {
        case Type::NONE:
            return 0.0;
            
        case Type::CONSTANT_RH:
            return constant_relative_humidity;
            
        case Type::EXPONENTIAL:
            return surface_relative_humidity * std::exp(-altitude / humidity_scale_height);
            
        case Type::TABULATED: {
            if (altitudes.empty()) return 0.0;
            if (altitude <= altitudes.front()) return relative_humidities.front();
            if (altitude >= altitudes.back()) return relative_humidities.back();
            
            auto it = std::lower_bound(altitudes.begin(), altitudes.end(), altitude);
            size_t i = it - altitudes.begin();
            if (i == 0) i = 1;
            double t = (altitude - altitudes[i-1]) / (altitudes[i] - altitudes[i-1]);
            return relative_humidities[i-1] + t * (relative_humidities[i] - relative_humidities[i-1]);
        }
            
        case Type::FUNCTION:
            if (humidity_function) return humidity_function(altitude);
            return 0.0;
    }
    return 0.0;
}

double HumidityProfile::getSpecificHumidity(double altitude, double temperature, double pressure) const {
    double rh = getRelativeHumidity(altitude);
    double e_s = saturationVaporPressure(temperature);
    double e = rh * e_s;
    double w = 0.622 * e / (pressure - e);  // Mixing ratio
    return w / (1.0 + w);  // Specific humidity
}

double HumidityProfile::getAbsoluteHumidity(double altitude, double temperature, double pressure) const {
    (void)pressure;  // Absolute humidity depends on vapor pressure and temperature, not total pressure
    double rh = getRelativeHumidity(altitude);
    double e_s = saturationVaporPressure(temperature);
    double e = rh * e_s;
    return e * AtmosphereConstants::M_H2O / (AtmosphereConstants::R_UNIVERSAL * temperature);
}

double HumidityProfile::getVaporPressure(double altitude, double temperature, double pressure) const {
    (void)pressure;  // Vapor pressure derived from saturation pressure and relative humidity
    double rh = getRelativeHumidity(altitude);
    return rh * saturationVaporPressure(temperature);
}

double HumidityProfile::getMixingRatio(double altitude, double temperature, double pressure) const {
    double e = getVaporPressure(altitude, temperature, pressure);
    return 0.622 * e / (pressure - e);
}

double HumidityProfile::saturationVaporPressure(double temperature) {
    return AtmosphereUtils::saturationVaporPressure(temperature);
}

double HumidityProfile::saturationMixingRatio(double temperature, double pressure) {
    double e_s = saturationVaporPressure(temperature);
    return 0.622 * e_s / (pressure - e_s);
}

double HumidityProfile::dewPoint(double temperature, double relative_humidity) {
    double T_C = temperature - 273.15;
    double a = 17.27, b = 237.3;
    double gamma = a * T_C / (b + T_C) + std::log(relative_humidity);
    return 273.15 + b * gamma / (a - gamma);
}

double HumidityProfile::wetBulbTemperature(double temperature, double relative_humidity, double pressure) {
    (void)pressure;  // Simplified psychrometric equation at standard pressure
    // Iterative solution using psychrometric equation
    double T_wb = temperature;
    
    for (int iter = 0; iter < 20; iter++) {
        double e_s_wb = saturationVaporPressure(T_wb);
        double e = relative_humidity * saturationVaporPressure(temperature);
        double T_wb_new = temperature - (e_s_wb - e) * AtmosphereConstants::LATENT_HEAT_VAP / 
                         (AtmosphereConstants::CP_AIR * 1000.0);
        if (std::abs(T_wb_new - T_wb) < 0.01) break;
        T_wb = T_wb_new;
    }
    return T_wb;
}

// =============================================================================
// TurbulenceProfile Implementation
// =============================================================================

double TurbulenceProfile::getTKE(double altitude) const {
    switch (type) {
        case Type::NONE:
            return 0.0;
            
        case Type::CONSTANT:
            return constant_tke;
            
        case Type::BOUNDARY_LAYER: {
            if (altitude > boundary_layer_height) return 0.01;  // Free atmosphere
            double z_norm = altitude / boundary_layer_height;
            // Typical convective boundary layer profile
            double tke = friction_velocity * friction_velocity * (1.0 - z_norm) * (1.0 - z_norm);
            if (convective_velocity > 0) {
                tke += 0.4 * convective_velocity * convective_velocity * 
                       std::pow(z_norm, 2.0/3.0) * std::pow(1.0 - 0.8 * z_norm, 2);
            }
            return tke;
        }
            
        case Type::TABULATED: {
            if (altitudes.empty()) return 0.0;
            if (altitude <= altitudes.front()) return tke_values.front();
            if (altitude >= altitudes.back()) return tke_values.back();
            
            auto it = std::lower_bound(altitudes.begin(), altitudes.end(), altitude);
            size_t i = it - altitudes.begin();
            if (i == 0) i = 1;
            double t = (altitude - altitudes[i-1]) / (altitudes[i] - altitudes[i-1]);
            return tke_values[i-1] + t * (tke_values[i] - tke_values[i-1]);
        }
            
        case Type::FUNCTION:
            return 0.0;  // Would need function pointer
    }
    return 0.0;
}

double TurbulenceProfile::getDissipation(double altitude) const {
    if (type == Type::CONSTANT) return constant_epsilon;
    if (type == Type::TABULATED && !epsilon_values.empty()) {
        if (altitude <= altitudes.front()) return epsilon_values.front();
        if (altitude >= altitudes.back()) return epsilon_values.back();
        
        auto it = std::lower_bound(altitudes.begin(), altitudes.end(), altitude);
        size_t i = it - altitudes.begin();
        if (i == 0) i = 1;
        double t = (altitude - altitudes[i-1]) / (altitudes[i] - altitudes[i-1]);
        return epsilon_values[i-1] + t * (epsilon_values[i] - epsilon_values[i-1]);
    }
    
    // Estimate from TKE using mixing length
    double tke = getTKE(altitude);
    double l = getMixingLength(altitude);
    if (l > 0 && tke > 0) {
        return std::pow(tke, 1.5) / l;
    }
    return 0.0;
}

double TurbulenceProfile::getEddyViscosity(double altitude) const {
    double tke = getTKE(altitude);
    double eps = getDissipation(altitude);
    if (eps > 0) {
        return 0.09 * tke * tke / eps;  // k-epsilon model: nu_t = C_mu * k² / epsilon
    }
    return 0.0;
}

double TurbulenceProfile::getEddyDiffusivity(double altitude) const {
    // Assume turbulent Prandtl number of ~0.85
    return getEddyViscosity(altitude) / 0.85;
}

double TurbulenceProfile::getMixingLength(double altitude) const {
    const double kappa = 0.4;  // von Karman constant
    
    if (altitude < boundary_layer_height) {
        // Blackadar formula
        double l_inf = 0.00027 * std::abs(AtmosphereConstants::G0 / 
                      (AtmosphereConstants::OMEGA_EARTH * std::sin(M_PI/4)));  // ~40m typical
        l_inf = std::min(l_inf, 200.0);  // Cap at 200m
        return kappa * altitude / (1.0 + kappa * altitude / l_inf);
    }
    return 50.0;  // Free atmosphere
}

// =============================================================================
// AtmosphericModelBase Implementation
// =============================================================================

double AtmosphericModelBase::meanMolecularWeight(double /*altitude*/) const {
    return AtmosphereConstants::M_AIR_DRY;
}

double AtmosphericModelBase::specificGasConstant(double altitude) const {
    double M = meanMolecularWeight(altitude);
    return AtmosphereConstants::R_UNIVERSAL * 1000.0 / M;  // J/(kg·K)
}

double AtmosphericModelBase::gamma(double /*altitude*/) const {
    return AtmosphereConstants::GAMMA_AIR;
}

double AtmosphericModelBase::soundSpeed(double altitude) const {
    return AtmosphereUtils::soundSpeed(temperature(altitude), gamma(altitude));
}

double AtmosphericModelBase::scaleHeight(double altitude) const {
    double T = temperature(altitude);
    double g = gravity(altitude);
    double R = specificGasConstant(altitude);
    return R * T / g;
}

double AtmosphericModelBase::numberDensity(double altitude) const {
    double p = pressure(altitude);
    double T = temperature(altitude);
    return p / (AtmosphereConstants::BOLTZMANN * T);
}

double AtmosphericModelBase::potentialTemperature(double altitude) const {
    return AtmosphereUtils::potentialTemperature(temperature(altitude), pressure(altitude));
}

double AtmosphericModelBase::dynamicViscosity(double altitude) const {
    return AtmosphereUtils::dynamicViscosity(temperature(altitude));
}

double AtmosphericModelBase::kinematicViscosity(double altitude) const {
    return dynamicViscosity(altitude) / density(altitude);
}

double AtmosphericModelBase::thermalConductivity(double altitude) const {
    return AtmosphereUtils::thermalConductivity(temperature(altitude));
}

double AtmosphericModelBase::meanFreePath(double altitude) const {
    double n = numberDensity(altitude);
    const double sigma = 3.65e-10;  // Effective molecular diameter for air
    return 1.0 / (std::sqrt(2.0) * M_PI * sigma * sigma * n);
}

double AtmosphericModelBase::geometricToGeopotential(double geometric_alt) const {
    return AtmosphereUtils::geometricToGeopotential(geometric_alt);
}

double AtmosphericModelBase::geopotentialToGeometric(double geopotential_alt) const {
    return AtmosphereUtils::geopotentialToGeometric(geopotential_alt);
}

double AtmosphericModelBase::gravity(double altitude) const {
    return AtmosphereUtils::gravity(altitude, 45.0);
}

AtmosphericState AtmosphericModelBase::getState(double altitude) const {
    AtmosphericState state;
    
    state.altitude = altitude;
    state.geopotential_altitude = geometricToGeopotential(altitude);
    state.latitude = 45.0;
    state.longitude = 0.0;
    
    state.temperature = temperature(altitude);
    state.pressure = pressure(altitude);
    state.density = density(altitude);
    state.potential_temperature = potentialTemperature(altitude);
    
    state.mean_molecular_weight = meanMolecularWeight(altitude);
    state.specific_gas_constant = specificGasConstant(altitude);
    state.gamma = gamma(altitude);
    
    // Moisture
    if (has_humidity_) {
        state.relative_humidity = humidity_profile_.getRelativeHumidity(altitude);
        state.specific_humidity = humidity_profile_.getSpecificHumidity(
            altitude, state.temperature, state.pressure);
        state.mixing_ratio = humidity_profile_.getMixingRatio(
            altitude, state.temperature, state.pressure);
        state.vapor_pressure = humidity_profile_.getVaporPressure(
            altitude, state.temperature, state.pressure);
        state.dew_point = HumidityProfile::dewPoint(state.temperature, state.relative_humidity);
        state.virtual_temperature = AtmosphereUtils::virtualTemperature(
            state.temperature, state.mixing_ratio);
    } else {
        state.relative_humidity = 0.0;
        state.specific_humidity = 0.0;
        state.mixing_ratio = 0.0;
        state.vapor_pressure = 0.0;
        state.dew_point = 0.0;
        state.virtual_temperature = state.temperature;
    }
    
    // Derived quantities
    state.sound_speed = soundSpeed(altitude);
    state.dynamic_viscosity = dynamicViscosity(altitude);
    state.kinematic_viscosity = kinematicViscosity(altitude);
    state.thermal_conductivity = thermalConductivity(altitude);
    state.mean_free_path = meanFreePath(altitude);
    
    // Wind
    if (has_wind_) {
        wind_profile_.getWind(altitude, state.wind_u, state.wind_v, state.wind_w);
        state.wind_speed = std::sqrt(state.wind_u * state.wind_u + state.wind_v * state.wind_v);
        state.wind_direction = wind_profile_.getDirection(altitude);
    } else {
        state.wind_u = state.wind_v = state.wind_w = 0.0;
        state.wind_speed = 0.0;
        state.wind_direction = 0.0;
    }
    
    // Turbulence
    if (has_turbulence_) {
        state.tke = turbulence_profile_.getTKE(altitude);
        state.epsilon = turbulence_profile_.getDissipation(altitude);
        state.eddy_viscosity = turbulence_profile_.getEddyViscosity(altitude);
    } else {
        state.tke = 0.0;
        state.epsilon = 0.0;
        state.eddy_viscosity = 0.0;
    }
    
    // Stability (simplified)
    state.brunt_vaisala_frequency = 0.01;  // Would need derivative
    state.richardson_number = 1.0;
    state.potential_vorticity = 0.0;
    
    // Radiation
    state.ozone_density = 0.0;
    state.aerosol_optical_depth = 0.0;
    
    // Layer
    state.layer_index = getLayerIndex(altitude);
    state.layer_name = getLayerName(altitude);
    
    return state;
}

int AtmosphericModelBase::getLayerIndex(double altitude) const {
    double h = geometricToGeopotential(altitude);
    for (size_t i = 0; i < layers_.size(); i++) {
        if (h >= layers_[i].base_altitude && h < layers_[i].top_altitude) {
            return static_cast<int>(i);
        }
    }
    return static_cast<int>(layers_.size()) - 1;
}

std::string AtmosphericModelBase::getLayerName(double altitude) const {
    int idx = getLayerIndex(altitude);
    if (idx >= 0 && idx < static_cast<int>(layers_.size())) {
        return layers_[idx].name;
    }
    return "Unknown";
}

const std::vector<AtmosphericLayer>& AtmosphericModelBase::getLayers() const {
    return layers_;
}

bool AtmosphericModelBase::isValidAltitude(double altitude) const {
    return altitude >= getMinAltitude() && altitude <= getMaxAltitude();
}

void AtmosphericModelBase::setWindProfile(const WindProfile& wind) {
    wind_profile_ = wind;
    has_wind_ = (wind.type != WindProfile::Type::NONE);
}

void AtmosphericModelBase::getWind(double altitude, double& u, double& v, double& w) const {
    if (has_wind_) {
        wind_profile_.getWind(altitude, u, v, w);
    } else {
        u = v = w = 0.0;
    }
}

void AtmosphericModelBase::setHumidityProfile(const HumidityProfile& humidity) {
    humidity_profile_ = humidity;
    has_humidity_ = (humidity.type != HumidityProfile::Type::NONE);
}

double AtmosphericModelBase::getRelativeHumidity(double altitude) const {
    if (has_humidity_) {
        return humidity_profile_.getRelativeHumidity(altitude);
    }
    return 0.0;
}

double AtmosphericModelBase::getSpecificHumidity(double altitude) const {
    if (has_humidity_) {
        return humidity_profile_.getSpecificHumidity(altitude, temperature(altitude), pressure(altitude));
    }
    return 0.0;
}

void AtmosphericModelBase::setTurbulenceProfile(const TurbulenceProfile& turbulence) {
    turbulence_profile_ = turbulence;
    has_turbulence_ = (turbulence.type != TurbulenceProfile::Type::NONE);
}

double AtmosphericModelBase::getTKE(double altitude) const {
    if (has_turbulence_) {
        return turbulence_profile_.getTKE(altitude);
    }
    return 0.0;
}

double AtmosphericModelBase::getEddyViscosity(double altitude) const {
    if (has_turbulence_) {
        return turbulence_profile_.getEddyViscosity(altitude);
    }
    return 0.0;
}

bool AtmosphericModelBase::loadFromConfig(const std::string& /*filename*/) {
    // To be implemented by derived classes
    return false;
}

bool AtmosphericModelBase::saveToConfig(const std::string& /*filename*/) const {
    // To be implemented by derived classes
    return false;
}

// =============================================================================
// USStandardAtmosphere1976 Implementation
// =============================================================================

USStandardAtmosphere1976::USStandardAtmosphere1976() {
    initializeLayers();
    computeBasePressures();
}

void USStandardAtmosphere1976::initializeLayers() {
    layers_.clear();
    
    // Troposphere (0-11 km)
    layers_.push_back({"Troposphere", 0.0, 11000.0, 288.15, -0.0065, 
                       AtmosphereConstants::P0_STD, AtmosphereConstants::RHO0_STD,
                       AtmosphereConstants::M_AIR_DRY, false});
    
    // Tropopause/Lower Stratosphere (11-20 km)
    layers_.push_back({"Tropopause", 11000.0, 20000.0, 216.65, 0.0, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, true});
    
    // Stratosphere (20-32 km)
    layers_.push_back({"Stratosphere", 20000.0, 32000.0, 216.65, 0.001, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
    
    // Upper Stratosphere (32-47 km)
    layers_.push_back({"Upper Stratosphere", 32000.0, 47000.0, 228.65, 0.0028, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
    
    // Stratopause (47-51 km)
    layers_.push_back({"Stratopause", 47000.0, 51000.0, 270.65, 0.0, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, true});
    
    // Mesosphere (51-71 km)
    layers_.push_back({"Mesosphere", 51000.0, 71000.0, 270.65, -0.0028, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
    
    // Upper Mesosphere (71-84.852 km)
    layers_.push_back({"Upper Mesosphere", 71000.0, 84852.0, 214.65, -0.002, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
}

void USStandardAtmosphere1976::computeBasePressures() {
    const double g = AtmosphereConstants::G0;
    const double R = AtmosphereConstants::R_AIR_DRY;
    
    p_bases_[0] = AtmosphereConstants::P0_STD;
    
    for (size_t i = 1; i < H_BOUNDS.size(); i++) {
        double dH = H_BOUNDS[i] - H_BOUNDS[i-1];
        double T_base = T_BASES[i-1];
        double lapse = LAPSE_RATES[i-1];
        
        if (std::abs(lapse) < 1e-10) {
            // Isothermal layer
            p_bases_[i] = p_bases_[i-1] * std::exp(-g * dH / (R * T_base));
        } else {
            // Non-isothermal layer
            double T_top = T_base + lapse * dH;
            p_bases_[i] = p_bases_[i-1] * std::pow(T_top / T_base, -g / (lapse * R));
        }
    }
}

double USStandardAtmosphere1976::temperature(double altitude) const {
    double H = geometricToGeopotential(altitude);
    
    // Find layer
    size_t layer = 0;
    for (size_t i = 0; i < H_BOUNDS.size() - 1; i++) {
        if (H >= H_BOUNDS[i] && H < H_BOUNDS[i+1]) {
            layer = i;
            break;
        }
        if (i == H_BOUNDS.size() - 2) layer = i;
    }
    
    if (H < 0) return T_BASES[0];
    if (H >= H_BOUNDS.back()) {
        if (extended_model_) {
            // Simplified extended model
            return 186.946 + (altitude - 86000.0) * 0.003;  // Warming in thermosphere
        }
        return T_BASES.back();
    }
    
    double dH = H - H_BOUNDS[layer];
    return T_BASES[layer] + LAPSE_RATES[layer] * dH;
}

double USStandardAtmosphere1976::pressure(double altitude) const {
    double H = geometricToGeopotential(altitude);
    
    if (H < 0) {
        // Below sea level extrapolation
        return AtmosphereConstants::P0_STD * std::exp(
            AtmosphereConstants::G0 * H / (AtmosphereConstants::R_AIR_DRY * T_BASES[0]));
    }
    
    if (H >= H_BOUNDS.back()) {
        if (extended_model_) {
            // Simplified exponential decay above 86 km
            double T = temperature(altitude);
            double scale_h = AtmosphereConstants::R_AIR_DRY * T / AtmosphereConstants::G0;
            return p_bases_.back() * std::exp(-(H - H_BOUNDS.back()) / scale_h);
        }
        return p_bases_.back();
    }
    
    // Find layer
    size_t layer = 0;
    for (size_t i = 0; i < H_BOUNDS.size() - 1; i++) {
        if (H >= H_BOUNDS[i] && H < H_BOUNDS[i+1]) {
            layer = i;
            break;
        }
    }
    
    const double g = AtmosphereConstants::G0;
    const double R = AtmosphereConstants::R_AIR_DRY;
    double dH = H - H_BOUNDS[layer];
    double T_base = T_BASES[layer];
    double lapse = LAPSE_RATES[layer];
    
    if (std::abs(lapse) < 1e-10) {
        return p_bases_[layer] * std::exp(-g * dH / (R * T_base));
    } else {
        double T = T_base + lapse * dH;
        return p_bases_[layer] * std::pow(T / T_base, -g / (lapse * R));
    }
}

double USStandardAtmosphere1976::density(double altitude) const {
    double T = temperature(altitude);
    double P = pressure(altitude);
    return P / (AtmosphereConstants::R_AIR_DRY * T);
}

double USStandardAtmosphere1976::meanMolecularWeight(double altitude) const {
    if (altitude < 86000.0) {
        return AtmosphereConstants::M_AIR_DRY;
    }
    
    // Above 86 km, composition changes
    if (extended_model_) {
        // Simplified model: molecular weight decreases due to atomic species
        double H = altitude / 1000.0;  // km
        if (H < 100) return 28.9644;
        if (H < 200) return 28.9644 - 0.05 * (H - 100);  // Linear decrease
        if (H < 500) return 23.9644 - 0.03 * (H - 200);
        return 16.0;  // Mostly atomic oxygen
    }
    return AtmosphereConstants::M_AIR_DRY;
}

// =============================================================================
// ICAOStandardAtmosphere Implementation
// =============================================================================

ICAOStandardAtmosphere::ICAOStandardAtmosphere() {
    initializeLayers();
}

void ICAOStandardAtmosphere::initializeLayers() {
    layers_.clear();
    
    // ICAO layers (similar to US Standard, but up to 80 km with different conventions)
    layers_.push_back({"Troposphere", 0.0, 11000.0, 288.15, -0.0065, 
                       AtmosphereConstants::P0_STD, AtmosphereConstants::RHO0_STD,
                       AtmosphereConstants::M_AIR_DRY, false});
    
    layers_.push_back({"Lower Stratosphere", 11000.0, 20000.0, 216.65, 0.0, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, true});
    
    layers_.push_back({"Middle Stratosphere", 20000.0, 32000.0, 216.65, 0.001, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
    
    layers_.push_back({"Upper Stratosphere", 32000.0, 47000.0, 228.65, 0.0028, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
    
    layers_.push_back({"Stratopause", 47000.0, 51000.0, 270.65, 0.0, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, true});
    
    layers_.push_back({"Lower Mesosphere", 51000.0, 71000.0, 270.65, -0.0028, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
    
    layers_.push_back({"Upper Mesosphere", 71000.0, 80000.0, 214.65, -0.002, 
                       0.0, 0.0, AtmosphereConstants::M_AIR_DRY, false});
}

double ICAOStandardAtmosphere::temperature(double altitude) const {
    // ICAO is essentially the same as US Standard up to 80 km
    double H = geometricToGeopotential(altitude);
    
    if (H < 0) return 288.15;
    if (H < 11000.0) return 288.15 - 0.0065 * H;
    if (H < 20000.0) return 216.65;
    if (H < 32000.0) return 216.65 + 0.001 * (H - 20000.0);
    if (H < 47000.0) return 228.65 + 0.0028 * (H - 32000.0);
    if (H < 51000.0) return 270.65;
    if (H < 71000.0) return 270.65 - 0.0028 * (H - 51000.0);
    if (H < 80000.0) return 214.65 - 0.002 * (H - 71000.0);
    return 196.65;
}

double ICAOStandardAtmosphere::pressure(double altitude) const {
    double H = geometricToGeopotential(altitude);
    const double g = AtmosphereConstants::G0;
    const double R = AtmosphereConstants::R_AIR_DRY;
    
    // Compute pressure layer by layer
    double P = AtmosphereConstants::P0_STD;
    double T_base, T_top, lapse, dH;
    
    // Layer boundaries and lapse rates
    const double H_bounds[] = {0, 11000, 20000, 32000, 47000, 51000, 71000, 80000};
    const double T_bases[] = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 196.65};
    const double lapse_rates[] = {-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0};
    
    if (H < 0) {
        return P * std::exp(g * (-H) / (R * 288.15));
    }
    
    for (int i = 0; i < 7; i++) {
        if (H <= H_bounds[i+1]) {
            dH = H - H_bounds[i];
            T_base = T_bases[i];
            lapse = lapse_rates[i];
            
            if (std::abs(lapse) < 1e-10) {
                return P * std::exp(-g * dH / (R * T_base));
            } else {
                T_top = T_base + lapse * dH;
                return P * std::pow(T_top / T_base, -g / (lapse * R));
            }
        }
        
        // Move to next layer
        dH = H_bounds[i+1] - H_bounds[i];
        T_base = T_bases[i];
        lapse = lapse_rates[i];
        
        if (std::abs(lapse) < 1e-10) {
            P *= std::exp(-g * dH / (R * T_base));
        } else {
            T_top = T_base + lapse * dH;
            P *= std::pow(T_top / T_base, -g / (lapse * R));
        }
    }
    
    // Above 80 km
    return P * std::exp(-g * (H - 80000.0) / (R * 196.65));
}

double ICAOStandardAtmosphere::density(double altitude) const {
    return pressure(altitude) / (AtmosphereConstants::R_AIR_DRY * temperature(altitude));
}

double ICAOStandardAtmosphere::pressureAltitude(double p) const {
    // Iterative solution (Newton-Raphson)
    double H = 0.0;
    for (int i = 0; i < 20; i++) {
        double P_calc = pressure(H);
        double dP_dH = -pressure(H) * AtmosphereConstants::G0 / 
                       (AtmosphereConstants::R_AIR_DRY * temperature(H));
        H = H - (P_calc - p) / dP_dH;
    }
    return H;
}

double ICAOStandardAtmosphere::densityAltitude(double rho) const {
    // Iterative solution
    double H = 0.0;
    for (int i = 0; i < 20; i++) {
        double rho_calc = density(H);
        if (std::abs(rho_calc - rho) < 1e-10) break;
        // Approximate gradient
        double drho_dH = (density(H + 100) - rho_calc) / 100.0;
        H = H - (rho_calc - rho) / drho_dH;
    }
    return H;
}

double ICAOStandardAtmosphere::indicatedAltitude(double p, double qnh) const {
    // QNH correction for altimeter setting
    // Returns indicated altitude given static pressure and QNH setting
    const double T0 = 288.15;
    const double L = 0.0065;
    const double g = AtmosphereConstants::G0;
    const double R = AtmosphereConstants::R_AIR_DRY;
    
    // ISA altitude from pressure
    double H_isa = (T0 / L) * (1.0 - std::pow(p / AtmosphereConstants::P0_STD, L * R / g));
    
    // QNH correction
    double H_qnh = (T0 / L) * (1.0 - std::pow(AtmosphereConstants::P0_STD / qnh, L * R / g));
    
    return H_isa - H_qnh;
}

// =============================================================================
// NRLMSISE00Atmosphere Implementation (Simplified)
// =============================================================================

NRLMSISE00Atmosphere::NRLMSISE00Atmosphere() {
    computeProfile();
}

void NRLMSISE00Atmosphere::computeProfile() {
    // Build altitude grid
    alt_grid_.clear();
    T_profile_.clear();
    rho_profile_.clear();
    
    // Generate grid from 0 to 1000 km
    for (double h = 0; h <= 1000000; h += 1000) {
        alt_grid_.push_back(h);
    }
    
    // Compute simplified profiles based on MSISE model behavior
    // This is a parameterized approximation, not the full model
    for (double h : alt_grid_) {
        double z = h / 1000.0;  // km
        
        // Temperature profile
        double T;
        if (z < 86) {
            // Use US Standard below 86 km
            USStandardAtmosphere1976 us_std;
            T = us_std.temperature(h);
        } else if (z < 120) {
            // Transition region
            T = 186.946 + (z - 86) * 2.5;
        } else {
            // Exospheric temperature (simplified Bates profile)
            double T_inf = 800.0 + 4.0 * (f107a_ - 70.0) + 
                          2.0 * (f107_ - f107a_) + 15.0 * ap_;
            T_inf = std::max(500.0, std::min(T_inf, 2500.0));
            double sigma = 0.02;  // Shape parameter
            T = T_inf - (T_inf - 350.0) * std::exp(-sigma * (z - 120));
        }
        T_profile_.push_back(T);
        
        // Density profile (simplified)
        double rho;
        if (z < 86) {
            USStandardAtmosphere1976 us_std;
            rho = us_std.density(h);
        } else {
            // Exponential with varying scale height
            double H_scale = AtmosphereConstants::R_AIR_DRY * T / AtmosphereConstants::G0;
            // Adjust for F10.7 and Ap
            double factor = 1.0 + 0.01 * (f107_ - 150.0) / 50.0;
            double rho_86 = 6.96e-6;  // kg/m³ at 86 km
            rho = factor * rho_86 * std::exp(-(h - 86000.0) / H_scale);
        }
        rho_profile_.push_back(rho);
    }
    
    profile_valid_ = true;
}

double NRLMSISE00Atmosphere::temperature(double altitude) const {
    if (!profile_valid_) const_cast<NRLMSISE00Atmosphere*>(this)->computeProfile();
    
    if (altitude < 0) return T_profile_.front();
    if (altitude >= alt_grid_.back()) return T_profile_.back();
    
    // Linear interpolation
    auto it = std::lower_bound(alt_grid_.begin(), alt_grid_.end(), altitude);
    size_t i = it - alt_grid_.begin();
    if (i == 0) return T_profile_.front();
    
    double t = (altitude - alt_grid_[i-1]) / (alt_grid_[i] - alt_grid_[i-1]);
    return T_profile_[i-1] + t * (T_profile_[i] - T_profile_[i-1]);
}

double NRLMSISE00Atmosphere::pressure(double altitude) const {
    // Compute from ideal gas law
    return density(altitude) * specificGasConstant(altitude) * temperature(altitude);
}

double NRLMSISE00Atmosphere::density(double altitude) const {
    if (!profile_valid_) const_cast<NRLMSISE00Atmosphere*>(this)->computeProfile();
    
    if (altitude < 0) return rho_profile_.front();
    if (altitude >= alt_grid_.back()) return rho_profile_.back();
    
    // Log-linear interpolation for density
    auto it = std::lower_bound(alt_grid_.begin(), alt_grid_.end(), altitude);
    size_t i = it - alt_grid_.begin();
    if (i == 0) return rho_profile_.front();
    
    double t = (altitude - alt_grid_[i-1]) / (alt_grid_[i] - alt_grid_[i-1]);
    double log_rho = std::log(rho_profile_[i-1]) + 
                     t * (std::log(rho_profile_[i]) - std::log(rho_profile_[i-1]));
    return std::exp(log_rho);
}

double NRLMSISE00Atmosphere::meanMolecularWeight(double altitude) const {
    double z = altitude / 1000.0;  // km
    
    if (z < 86) return AtmosphereConstants::M_AIR_DRY;
    if (z < 200) {
        // Transition from N2/O2 to atomic O dominated
        double f = (z - 86) / (200 - 86);
        return AtmosphereConstants::M_AIR_DRY * (1.0 - 0.45 * f);
    }
    if (z < 500) {
        // O, He, H mix
        double f = (z - 200) / (500 - 200);
        return 16.0 * (1.0 - 0.6 * f) + 4.0 * 0.4 * f;
    }
    // Dominated by H and He
    return 4.0;
}

double NRLMSISE00Atmosphere::N2_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 86) return 0.78084 * total;
    if (z < 200) {
        double f = std::exp(-(z - 86) / 20.0);
        return 0.78084 * total * f;
    }
    return 0.0;
}

double NRLMSISE00Atmosphere::O2_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 86) return 0.20946 * total;
    if (z < 150) {
        double f = std::exp(-(z - 86) / 15.0);
        return 0.20946 * total * f;
    }
    return 0.0;
}

double NRLMSISE00Atmosphere::O_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 86) return 0.0;
    if (z < 500) {
        // Atomic oxygen peaks around 200-300 km
        double peak = std::exp(-std::pow((z - 250) / 100, 2));
        return total * 0.8 * peak;
    }
    return total * 0.1;
}

double NRLMSISE00Atmosphere::He_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 200) return total * 5.24e-6;  // Constant mixing ratio
    // He becomes dominant at high altitude
    double f = std::min(1.0, (z - 200) / 500.0);
    return total * (5.24e-6 + f * 0.3);
}

double NRLMSISE00Atmosphere::H_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 500) return 0.0;
    // H dominant above 500 km
    double f = std::min(1.0, (z - 500) / 500.0);
    return total * f * 0.7;
}

double NRLMSISE00Atmosphere::Ar_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 100) return 0.00934 * total;
    return 0.00934 * total * std::exp(-(z - 100) / 10.0);
}

double NRLMSISE00Atmosphere::N_density(double altitude) const {
    double z = altitude / 1000.0;
    double total = density(altitude);
    
    if (z < 150) return 0.0;
    // Similar profile to atomic O but lower
    double peak = std::exp(-std::pow((z - 200) / 80, 2));
    return total * 0.05 * peak;
}

void NRLMSISE00Atmosphere::setF107(double f107) {
    f107_ = f107;
    profile_valid_ = false;
}

void NRLMSISE00Atmosphere::setF107Average(double f107a) {
    f107a_ = f107a;
    profile_valid_ = false;
}

void NRLMSISE00Atmosphere::setAp(double ap) {
    ap_ = ap;
    profile_valid_ = false;
}

void NRLMSISE00Atmosphere::setApArray(const std::array<double, 7>& ap_array) {
    ap_array_ = ap_array;
    ap_ = ap_array[0];
    profile_valid_ = false;
}

void NRLMSISE00Atmosphere::setDateTime(int year, int day_of_year, double seconds_of_day) {
    year_ = year;
    day_of_year_ = day_of_year;
    seconds_ = seconds_of_day;
    profile_valid_ = false;
}

void NRLMSISE00Atmosphere::setLocation(double latitude, double longitude) {
    latitude_ = latitude;
    longitude_ = longitude;
    profile_valid_ = false;
}

// =============================================================================
// ReferenceAtmosphere Implementation (AFGL/HITRAN profiles)
// =============================================================================

ReferenceAtmosphere::ReferenceAtmosphere(ReferenceAtmosphereType type) : type_(type) {
    initializeProfiles();
}

std::string ReferenceAtmosphere::getName() const {
    switch (type_) {
        case ReferenceAtmosphereType::TROPICAL: return "TROPICAL";
        case ReferenceAtmosphereType::MIDLATITUDE_SUMMER: return "MIDLATITUDE_SUMMER";
        case ReferenceAtmosphereType::MIDLATITUDE_WINTER: return "MIDLATITUDE_WINTER";
        case ReferenceAtmosphereType::SUBARCTIC_SUMMER: return "SUBARCTIC_SUMMER";
        case ReferenceAtmosphereType::SUBARCTIC_WINTER: return "SUBARCTIC_WINTER";
        case ReferenceAtmosphereType::US_STANDARD: return "US_STANDARD";
    }
    return "UNKNOWN";
}

std::string ReferenceAtmosphere::getDescription() const {
    switch (type_) {
        case ReferenceAtmosphereType::TROPICAL: 
            return "AFGL Tropical Atmosphere (15°N annual average)";
        case ReferenceAtmosphereType::MIDLATITUDE_SUMMER: 
            return "AFGL Mid-latitude Summer (45°N July)";
        case ReferenceAtmosphereType::MIDLATITUDE_WINTER: 
            return "AFGL Mid-latitude Winter (45°N January)";
        case ReferenceAtmosphereType::SUBARCTIC_SUMMER: 
            return "AFGL Subarctic Summer (60°N July)";
        case ReferenceAtmosphereType::SUBARCTIC_WINTER: 
            return "AFGL Subarctic Winter (60°N January)";
        case ReferenceAtmosphereType::US_STANDARD: 
            return "US Standard Atmosphere 1976";
    }
    return "Unknown Reference Atmosphere";
}

void ReferenceAtmosphere::initializeProfiles() {
    // Standard altitude grid (km)
    altitudes_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                  20, 21, 22, 23, 24, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
                  90, 95, 100, 105, 110, 115, 120};
    
    switch (type_) {
        case ReferenceAtmosphereType::TROPICAL:
            temperatures_ = {300.0, 294.0, 288.0, 284.0, 277.0, 270.0, 264.0, 257.0, 250.0, 
                            244.0, 237.0, 230.0, 224.0, 217.0, 210.0, 204.0, 197.0, 195.0, 
                            199.0, 203.0, 207.0, 211.0, 215.0, 217.0, 219.0, 221.0, 232.0, 
                            243.0, 254.0, 265.0, 270.0, 263.0, 247.0, 233.0, 219.0, 208.0,
                            198.0, 189.0, 186.0, 188.0, 195.0, 208.0, 240.0, 300.0, 400.0};
            pressures_ = {1013.0, 904.0, 805.0, 715.0, 633.0, 559.0, 492.0, 432.0, 378.0, 
                         329.0, 286.0, 247.0, 213.0, 182.0, 156.0, 132.0, 111.0, 93.7, 
                         78.9, 66.6, 56.5, 48.0, 40.9, 35.0, 30.0, 25.7, 12.2, 5.75, 
                         2.87, 1.49, 0.798, 0.425, 0.219, 0.109, 0.0522, 0.0240, 0.0105,
                         0.00446, 0.00184, 0.000760, 0.000320, 0.000145, 0.0000710, 
                         0.0000390, 0.0000230};
            h2o_ = {19000, 13000, 9300, 4700, 2200, 1500, 850, 470, 250, 120, 76, 64, 
                   56, 50, 45, 40, 34, 28, 20, 14, 9.2, 6.0, 4.0, 3.4, 3.0, 2.9, 2.7, 
                   2.4, 2.0, 1.6, 1.3, 1.0, 0.8, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 
                   0.15, 0.1, 0.08, 0.06, 0.04};
            break;
            
        case ReferenceAtmosphereType::MIDLATITUDE_SUMMER:
            temperatures_ = {294.0, 290.0, 285.0, 279.0, 273.0, 267.0, 261.0, 255.0, 248.0, 
                            242.0, 235.0, 229.0, 222.0, 216.0, 216.0, 216.0, 216.0, 216.0, 
                            216.0, 217.0, 218.0, 219.0, 220.0, 222.0, 223.0, 224.0, 234.0, 
                            245.0, 258.0, 270.0, 276.0, 269.0, 257.0, 240.0, 218.0, 202.0,
                            187.0, 177.0, 177.0, 184.0, 199.0, 222.0, 262.0, 330.0, 440.0};
            pressures_ = {1013.0, 902.0, 802.0, 710.0, 628.0, 554.0, 487.0, 426.0, 372.0, 
                         324.0, 281.0, 243.0, 209.0, 179.0, 153.0, 130.0, 111.0, 95.0, 
                         81.2, 69.5, 59.5, 51.0, 43.7, 37.6, 32.2, 27.7, 13.2, 6.52, 
                         3.33, 1.76, 0.951, 0.514, 0.272, 0.136, 0.0648, 0.0295, 0.0128,
                         0.00535, 0.00220, 0.000910, 0.000380, 0.000170, 0.0000830, 
                         0.0000440, 0.0000260};
            h2o_ = {14000, 9300, 5900, 3300, 1900, 1100, 640, 380, 210, 120, 76, 64, 
                   55, 45, 38, 32, 26, 21, 16, 12, 8.5, 5.4, 3.5, 2.4, 1.7, 1.2, 0.75, 
                   0.58, 0.45, 0.35, 0.27, 0.21, 0.16, 0.12, 0.095, 0.074, 0.058, 0.046,
                   0.036, 0.028, 0.022, 0.017, 0.013, 0.010, 0.008};
            break;
            
        case ReferenceAtmosphereType::MIDLATITUDE_WINTER:
            temperatures_ = {272.2, 268.7, 265.2, 261.7, 255.7, 249.7, 243.7, 237.7, 231.7, 
                            225.7, 219.7, 219.2, 218.7, 218.2, 217.7, 217.2, 216.7, 216.2, 
                            215.7, 215.2, 215.2, 215.2, 215.2, 215.2, 215.2, 215.2, 217.4, 
                            227.8, 243.2, 258.5, 265.7, 260.6, 250.0, 235.0, 215.0, 198.0,
                            183.0, 173.0, 173.0, 180.0, 193.0, 213.0, 250.0, 315.0, 420.0};
            pressures_ = {1018.0, 897.3, 789.7, 693.8, 608.1, 531.3, 462.7, 401.6, 347.3, 
                         299.2, 256.8, 219.9, 188.2, 161.0, 137.8, 117.8, 100.7, 86.1, 
                         73.5, 62.8, 53.7, 45.8, 39.1, 33.4, 28.6, 24.3, 11.1, 5.18, 
                         2.53, 1.29, 0.682, 0.362, 0.188, 0.0947, 0.0460, 0.0214, 0.00956,
                         0.00411, 0.00174, 0.000720, 0.000299, 0.000135, 0.0000660, 
                         0.0000354, 0.0000210};
            h2o_ = {3500, 2500, 1800, 1200, 660, 380, 220, 150, 94, 60, 38, 21, 12, 
                   7.6, 5.0, 3.4, 2.3, 1.6, 1.1, 0.76, 0.64, 0.54, 0.45, 0.38, 0.32, 
                   0.27, 0.21, 0.16, 0.12, 0.09, 0.072, 0.058, 0.046, 0.036, 0.028, 
                   0.022, 0.017, 0.013, 0.010, 0.0078, 0.0061, 0.0048, 0.0038, 0.0030, 0.0024};
            break;
            
        case ReferenceAtmosphereType::SUBARCTIC_SUMMER:
            temperatures_ = {287.0, 282.0, 276.0, 271.0, 266.0, 260.0, 253.0, 246.0, 239.0, 
                            232.0, 225.0, 225.0, 225.0, 225.0, 225.0, 225.0, 225.0, 225.0, 
                            225.0, 225.0, 225.0, 225.0, 225.0, 225.0, 226.0, 228.0, 235.0, 
                            247.0, 262.0, 274.0, 277.0, 268.0, 254.0, 238.0, 216.0, 200.0,
                            186.0, 177.0, 177.0, 183.0, 195.0, 218.0, 258.0, 325.0, 435.0};
            pressures_ = {1010.0, 896.0, 792.9, 700.0, 616.0, 541.0, 473.0, 413.0, 359.0, 
                         310.5, 267.7, 230.0, 197.7, 170.0, 146.0, 125.3, 107.5, 92.0, 
                         78.9, 67.6, 58.0, 49.8, 42.6, 36.5, 31.3, 26.8, 12.4, 6.00, 
                         3.05, 1.59, 0.854, 0.458, 0.239, 0.120, 0.0575, 0.0263, 0.0114,
                         0.00478, 0.00197, 0.000820, 0.000340, 0.000155, 0.0000760, 
                         0.0000405, 0.0000240};
            h2o_ = {8600, 5400, 3400, 2000, 1100, 640, 380, 210, 120, 76, 64, 55, 
                   48, 42, 36, 31, 27, 23, 19, 15, 11, 8.0, 5.6, 3.9, 2.7, 1.9, 1.0, 
                   0.66, 0.47, 0.35, 0.27, 0.21, 0.16, 0.12, 0.095, 0.074, 0.058, 0.046,
                   0.036, 0.028, 0.022, 0.017, 0.013, 0.010, 0.008};
            break;
            
        case ReferenceAtmosphereType::SUBARCTIC_WINTER:
            temperatures_ = {257.1, 259.1, 255.9, 252.7, 247.7, 240.9, 234.1, 227.3, 220.6, 
                            217.2, 217.2, 217.2, 217.2, 217.2, 217.2, 217.2, 216.6, 216.0, 
                            215.4, 214.8, 214.1, 213.6, 213.0, 212.4, 211.8, 211.2, 216.0, 
                            222.2, 234.7, 247.0, 259.3, 259.1, 250.0, 235.0, 213.0, 197.0,
                            183.0, 173.0, 173.0, 180.0, 193.0, 213.0, 250.0, 315.0, 420.0};
            pressures_ = {1013.0, 887.8, 777.5, 679.8, 593.2, 515.8, 446.7, 385.3, 330.8, 
                         282.9, 241.8, 206.7, 176.6, 151.0, 129.1, 110.3, 94.31, 80.58, 
                         68.82, 58.75, 50.14, 42.77, 36.47, 31.09, 26.49, 22.56, 10.20, 
                         4.701, 2.243, 1.113, 0.5719, 0.2990, 0.1520, 0.0756, 0.0364, 
                         0.0167, 0.00729, 0.00310, 0.00130, 0.000540, 0.000227, 0.000103,
                         0.0000510, 0.0000274, 0.0000163};
            h2o_ = {1200, 1200, 940, 680, 410, 200, 98, 54, 29, 10, 6.0, 5.0, 
                   4.8, 4.5, 4.2, 3.9, 3.4, 2.9, 2.4, 1.9, 1.5, 1.1, 0.84, 0.62, 
                   0.46, 0.34, 0.19, 0.13, 0.091, 0.068, 0.052, 0.041, 0.032, 0.024, 
                   0.019, 0.015, 0.011, 0.009, 0.007, 0.0056, 0.0044, 0.0035, 0.0028, 
                   0.0022, 0.0018};
            break;
            
        case ReferenceAtmosphereType::US_STANDARD:
        default:
            // US Standard 1976 profile
            temperatures_ = {288.2, 281.7, 275.2, 268.7, 262.2, 255.7, 249.2, 242.7, 236.2, 
                            229.7, 223.3, 216.8, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7, 
                            216.7, 216.7, 216.7, 217.6, 218.6, 219.6, 220.6, 221.6, 226.5, 
                            237.0, 251.0, 265.0, 270.7, 260.0, 245.0, 230.0, 210.0, 195.0,
                            180.0, 170.0, 170.0, 178.0, 192.0, 215.0, 255.0, 320.0, 430.0};
            pressures_ = {1013.25, 898.8, 795.0, 701.2, 616.6, 540.5, 472.2, 411.0, 356.5, 
                         308.0, 265.0, 227.0, 194.0, 166.0, 142.0, 121.0, 103.5, 88.5, 
                         75.7, 64.7, 55.3, 47.3, 40.5, 34.6, 29.6, 25.3, 12.0, 5.75, 
                         2.87, 1.49, 0.798, 0.425, 0.219, 0.109, 0.0522, 0.0240, 0.0105,
                         0.00446, 0.00184, 0.000760, 0.000320, 0.000145, 0.0000710, 
                         0.0000390, 0.0000230};
            h2o_ = {7700, 6200, 4700, 3200, 1900, 1100, 640, 380, 210, 120, 76, 64, 
                   55, 45, 38, 32, 26, 21, 16, 12, 8.5, 5.4, 3.5, 2.4, 1.7, 1.2, 0.75, 
                   0.58, 0.45, 0.35, 0.27, 0.21, 0.16, 0.12, 0.095, 0.074, 0.058, 0.046,
                   0.036, 0.028, 0.022, 0.017, 0.013, 0.010, 0.008};
            break;
    }
    
    // CO2 (constant ~420 ppm below 80 km, then decreasing)
    co2_.resize(altitudes_.size());
    for (size_t i = 0; i < altitudes_.size(); i++) {
        if (altitudes_[i] < 80) {
            co2_[i] = 420.0;
        } else {
            co2_[i] = 420.0 * std::exp(-(altitudes_[i] - 80) / 20.0);
        }
    }
    
    // O3 (simplified ozone profile - peaks around 20-25 km)
    o3_.resize(altitudes_.size());
    for (size_t i = 0; i < altitudes_.size(); i++) {
        double z = altitudes_[i];
        o3_[i] = 8.0 * std::exp(-std::pow((z - 22.0) / 8.0, 2));  // ppmv
        if (z > 60) o3_[i] *= std::exp(-(z - 60) / 10.0);
    }
}

double ReferenceAtmosphere::interpolateProfile(double altitude, 
                                               const std::vector<double>& alts,
                                               const std::vector<double>& vals) const {
    double z_km = altitude / 1000.0;  // Convert to km
    
    if (z_km <= alts.front()) return vals.front();
    if (z_km >= alts.back()) return vals.back();
    
    auto it = std::lower_bound(alts.begin(), alts.end(), z_km);
    size_t i = it - alts.begin();
    if (i == 0) i = 1;
    
    double t = (z_km - alts[i-1]) / (alts[i] - alts[i-1]);
    return vals[i-1] + t * (vals[i] - vals[i-1]);
}

double ReferenceAtmosphere::temperature(double altitude) const {
    return interpolateProfile(altitude, altitudes_, temperatures_);
}

double ReferenceAtmosphere::pressure(double altitude) const {
    // Pressures are in mb, convert to Pa
    return interpolateProfile(altitude, altitudes_, pressures_) * 100.0;
}

double ReferenceAtmosphere::density(double altitude) const {
    double T = temperature(altitude);
    double P = pressure(altitude);
    return P / (AtmosphereConstants::R_AIR_DRY * T);
}

double ReferenceAtmosphere::H2O_mixingRatio(double altitude) const {
    return interpolateProfile(altitude, altitudes_, h2o_) * 1e-6;  // ppmv to ratio
}

double ReferenceAtmosphere::CO2_mixingRatio(double altitude) const {
    return interpolateProfile(altitude, altitudes_, co2_) * 1e-6;
}

double ReferenceAtmosphere::O3_mixingRatio(double altitude) const {
    return interpolateProfile(altitude, altitudes_, o3_) * 1e-6;
}

double ReferenceAtmosphere::N2O_mixingRatio(double altitude) const {
    // Simplified N2O profile (~320 ppb at surface, decreasing above stratosphere)
    double z_km = altitude / 1000.0;
    if (z_km < 15) return 320e-9;
    if (z_km < 50) return 320e-9 * std::exp(-(z_km - 15) / 20.0);
    return 10e-9;
}

double ReferenceAtmosphere::CO_mixingRatio(double altitude) const {
    // Simplified CO profile (~100 ppb)
    double z_km = altitude / 1000.0;
    if (z_km < 10) return 100e-9;
    return 100e-9 * std::exp(-(z_km - 10) / 30.0);
}

double ReferenceAtmosphere::CH4_mixingRatio(double altitude) const {
    // Simplified CH4 profile (~1.9 ppm at surface)
    double z_km = altitude / 1000.0;
    if (z_km < 15) return 1.9e-6;
    if (z_km < 50) return 1.9e-6 * std::exp(-(z_km - 15) / 25.0);
    return 0.2e-6;
}

// =============================================================================
// TabulatedAtmosphere Implementation
// =============================================================================

TabulatedAtmosphere::TabulatedAtmosphere() {}

TabulatedAtmosphere::TabulatedAtmosphere(const std::string& config_file) {
    loadFromConfig(config_file);
}

double TabulatedAtmosphere::interpolate(double altitude, const std::vector<double>& values,
                                        bool use_log) const {
    if (altitudes_.empty() || values.empty()) {
        throw std::runtime_error("TabulatedAtmosphere: No data loaded");
    }
    
    // Handle extrapolation
    if (altitude < altitudes_.front()) {
        switch (extrap_type_) {
            case ExtrapolationType::CONSTANT:
                return values.front();
            case ExtrapolationType::LINEAR: {
                if (altitudes_.size() < 2) return values.front();
                double slope = (values[1] - values[0]) / (altitudes_[1] - altitudes_[0]);
                return values.front() + slope * (altitude - altitudes_.front());
            }
            case ExtrapolationType::EXPONENTIAL:
                if (use_log && values.front() > 0) {
                    double H = (altitudes_[1] - altitudes_[0]) / 
                              std::log(values[0] / values[1]);
                    return values.front() * std::exp((altitudes_.front() - altitude) / H);
                }
                return values.front();
            case ExtrapolationType::ERROR:
                throw std::out_of_range("Altitude below tabulated range");
        }
    }
    
    if (altitude > altitudes_.back()) {
        switch (extrap_type_) {
            case ExtrapolationType::CONSTANT:
                return values.back();
            case ExtrapolationType::LINEAR: {
                size_t n = altitudes_.size();
                if (n < 2) return values.back();
                double slope = (values[n-1] - values[n-2]) / (altitudes_[n-1] - altitudes_[n-2]);
                return values.back() + slope * (altitude - altitudes_.back());
            }
            case ExtrapolationType::EXPONENTIAL:
                if (use_log && values.back() > 0) {
                    size_t n = altitudes_.size();
                    double H = (altitudes_[n-1] - altitudes_[n-2]) / 
                              std::log(values[n-2] / values[n-1]);
                    return values.back() * std::exp(-(altitude - altitudes_.back()) / H);
                }
                return values.back();
            case ExtrapolationType::ERROR:
                throw std::out_of_range("Altitude above tabulated range");
        }
    }
    
    // Find bracketing indices
    auto it = std::lower_bound(altitudes_.begin(), altitudes_.end(), altitude);
    size_t i = it - altitudes_.begin();
    if (i == 0) i = 1;
    
    double t = (altitude - altitudes_[i-1]) / (altitudes_[i] - altitudes_[i-1]);
    
    switch (interp_type_) {
        case InterpolationType::LINEAR:
            return values[i-1] + t * (values[i] - values[i-1]);
            
        case InterpolationType::LOG_LINEAR:
            if (use_log && values[i-1] > 0 && values[i] > 0) {
                double log_v = std::log(values[i-1]) + t * (std::log(values[i]) - std::log(values[i-1]));
                return std::exp(log_v);
            }
            return values[i-1] + t * (values[i] - values[i-1]);
            
        case InterpolationType::CUBIC_SPLINE:
        case InterpolationType::PCHIP:
            // Simplified: fall back to linear for now
            return values[i-1] + t * (values[i] - values[i-1]);
    }
    
    return values[i-1] + t * (values[i] - values[i-1]);
}

double TabulatedAtmosphere::temperature(double altitude) const {
    return interpolate(altitude, temperatures_, false);
}

double TabulatedAtmosphere::pressure(double altitude) const {
    return interpolate(altitude, pressures_, true);
}

double TabulatedAtmosphere::density(double altitude) const {
    if (has_density_data_ && !densities_.empty()) {
        return interpolate(altitude, densities_, true);
    }
    // Compute from ideal gas law
    return pressure(altitude) / (AtmosphereConstants::R_AIR_DRY * temperature(altitude));
}

double TabulatedAtmosphere::getMinAltitude() const {
    return altitudes_.empty() ? 0.0 : altitudes_.front();
}

double TabulatedAtmosphere::getMaxAltitude() const {
    return altitudes_.empty() ? 86000.0 : altitudes_.back();
}

void TabulatedAtmosphere::setAltitudes(const std::vector<double>& altitudes) {
    altitudes_ = altitudes;
}

void TabulatedAtmosphere::setTemperatures(const std::vector<double>& temperatures) {
    temperatures_ = temperatures;
}

void TabulatedAtmosphere::setPressures(const std::vector<double>& pressures) {
    pressures_ = pressures;
}

void TabulatedAtmosphere::setDensities(const std::vector<double>& densities) {
    densities_ = densities;
    has_density_data_ = !densities.empty();
}

void TabulatedAtmosphere::addDataPoint(double altitude, double temperature, double pressure,
                                       double density) {
    // Insert in sorted order
    auto it = std::lower_bound(altitudes_.begin(), altitudes_.end(), altitude);
    size_t idx = it - altitudes_.begin();
    
    altitudes_.insert(altitudes_.begin() + idx, altitude);
    temperatures_.insert(temperatures_.begin() + idx, temperature);
    pressures_.insert(pressures_.begin() + idx, pressure);
    
    if (density > 0) {
        densities_.insert(densities_.begin() + idx, density);
        has_density_data_ = true;
    } else if (has_density_data_) {
        // Compute from ideal gas
        double rho = pressure / (AtmosphereConstants::R_AIR_DRY * temperature);
        densities_.insert(densities_.begin() + idx, rho);
    }
}

void TabulatedAtmosphere::clearData() {
    altitudes_.clear();
    temperatures_.clear();
    pressures_.clear();
    densities_.clear();
    has_density_data_ = false;
}

bool TabulatedAtmosphere::loadFromConfig(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    clearData();
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        double alt, T, P, rho = -1;
        
        if (iss >> alt >> T >> P) {
            iss >> rho;  // Optional
            addDataPoint(alt, T, P, rho);
        }
    }
    
    return !altitudes_.empty();
}

bool TabulatedAtmosphere::loadFromCSV(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    clearData();
    std::string line;
    bool header = true;
    
    while (std::getline(file, line)) {
        if (header) { header = false; continue; }  // Skip header
        
        std::istringstream iss(line);
        std::string token;
        std::vector<double> values;
        
        while (std::getline(iss, token, ',')) {
            values.push_back(std::stod(token));
        }
        
        if (values.size() >= 3) {
            double rho = (values.size() >= 4) ? values[3] : -1;
            addDataPoint(values[0], values[1], values[2], rho);
        }
    }
    
    return !altitudes_.empty();
}

bool TabulatedAtmosphere::loadFromJSON(const std::string& /*filename*/) {
    // JSON parsing would require external library or more code
    // For now, not implemented
    return false;
}

// =============================================================================
// CompositeAtmosphere Implementation
// =============================================================================

CompositeAtmosphere::CompositeAtmosphere() {}

std::string CompositeAtmosphere::getDescription() const {
    std::string desc = "Composite atmosphere with " + 
                       std::to_string(models_.size()) + " regions";
    return desc;
}

void CompositeAtmosphere::addModel(std::shared_ptr<AtmosphericModelBase> model,
                                   double min_alt, double max_alt, double blend_width) {
    ModelRegion region;
    region.model = model;
    region.min_alt = min_alt;
    region.max_alt = max_alt;
    region.blend_width = blend_width;
    models_.push_back(region);
}

void CompositeAtmosphere::clearModels() {
    models_.clear();
}

double CompositeAtmosphere::getMinAltitude() const {
    if (models_.empty()) return 0.0;
    double min_alt = models_[0].min_alt;
    for (const auto& m : models_) {
        min_alt = std::min(min_alt, m.min_alt);
    }
    return min_alt;
}

double CompositeAtmosphere::getMaxAltitude() const {
    if (models_.empty()) return 86000.0;
    double max_alt = models_[0].max_alt;
    for (const auto& m : models_) {
        max_alt = std::max(max_alt, m.max_alt);
    }
    return max_alt;
}

double CompositeAtmosphere::blendedValue(double altitude,
    std::function<double(const AtmosphericModelBase*, double)> getter) const {
    
    if (models_.empty()) {
        throw std::runtime_error("CompositeAtmosphere: No models added");
    }
    
    // Find applicable models
    std::vector<std::pair<double, double>> contributions;  // weight, value
    
    for (const auto& region : models_) {
        if (altitude < region.min_alt - region.blend_width ||
            altitude > region.max_alt + region.blend_width) {
            continue;
        }
        
        double weight = 1.0;
        
        // Blend at lower boundary
        if (altitude < region.min_alt + region.blend_width && region.blend_width > 0) {
            double t = (altitude - (region.min_alt - region.blend_width)) / 
                      (2 * region.blend_width);
            weight *= 0.5 * (1.0 + std::tanh(4.0 * (t - 0.5)));
        }
        
        // Blend at upper boundary
        if (altitude > region.max_alt - region.blend_width && region.blend_width > 0) {
            double t = ((region.max_alt + region.blend_width) - altitude) / 
                      (2 * region.blend_width);
            weight *= 0.5 * (1.0 + std::tanh(4.0 * (t - 0.5)));
        }
        
        if (weight > 1e-6) {
            double value = getter(region.model.get(), altitude);
            contributions.push_back({weight, value});
        }
    }
    
    if (contributions.empty()) {
        // Use first model as fallback
        return getter(models_[0].model.get(), altitude);
    }
    
    // Weighted average
    double total_weight = 0.0;
    double result = 0.0;
    for (const auto& [w, v] : contributions) {
        total_weight += w;
        result += w * v;
    }
    
    return result / total_weight;
}

double CompositeAtmosphere::temperature(double altitude) const {
    return blendedValue(altitude, [](const AtmosphericModelBase* m, double z) {
        return m->temperature(z);
    });
}

double CompositeAtmosphere::pressure(double altitude) const {
    return blendedValue(altitude, [](const AtmosphericModelBase* m, double z) {
        return m->pressure(z);
    });
}

double CompositeAtmosphere::density(double altitude) const {
    return blendedValue(altitude, [](const AtmosphericModelBase* m, double z) {
        return m->density(z);
    });
}

// =============================================================================
// PerturbedAtmosphere Implementation
// =============================================================================

PerturbedAtmosphere::PerturbedAtmosphere(std::shared_ptr<AtmosphericModelBase> base)
    : base_model_(base) {}

std::string PerturbedAtmosphere::getName() const {
    return "PERTURBED_" + base_model_->getName();
}

std::string PerturbedAtmosphere::getDescription() const {
    return "Perturbed " + base_model_->getDescription();
}

double PerturbedAtmosphere::temperature(double altitude) const {
    double T = base_model_->temperature(altitude);
    for (const auto& pert : T_perturbations_) {
        T += pert(altitude);
    }
    return T;
}

double PerturbedAtmosphere::pressure(double altitude) const {
    double P = base_model_->pressure(altitude);
    for (const auto& pert : P_perturbations_) {
        P += pert(altitude);
    }
    return std::max(P, 0.0);  // Ensure non-negative
}

double PerturbedAtmosphere::density(double altitude) const {
    double rho = base_model_->density(altitude);
    for (const auto& pert : rho_perturbations_) {
        rho += pert(altitude);
    }
    return std::max(rho, 0.0);  // Ensure non-negative
}

void PerturbedAtmosphere::addTemperaturePerturbation(std::function<double(double)> dT) {
    T_perturbations_.push_back(dT);
}

void PerturbedAtmosphere::addPressurePerturbation(std::function<double(double)> dP) {
    P_perturbations_.push_back(dP);
}

void PerturbedAtmosphere::addDensityPerturbation(std::function<double(double)> dRho) {
    rho_perturbations_.push_back(dRho);
}

void PerturbedAtmosphere::addGravityWave(double amplitude, double wavelength, double phase) {
    auto wave = [amplitude, wavelength, phase](double z) {
        return amplitude * std::sin(2.0 * M_PI * z / wavelength + phase);
    };
    T_perturbations_.push_back(wave);
}

void PerturbedAtmosphere::addLocalizedPerturbation(double amplitude, double center_alt, 
                                                   double width, const std::string& variable) {
    auto gaussian = [amplitude, center_alt, width](double z) {
        return amplitude * std::exp(-std::pow((z - center_alt) / width, 2));
    };
    
    if (variable == "temperature" || variable == "T") {
        T_perturbations_.push_back(gaussian);
    } else if (variable == "pressure" || variable == "P") {
        P_perturbations_.push_back(gaussian);
    } else if (variable == "density" || variable == "rho") {
        rho_perturbations_.push_back(gaussian);
    }
}

void PerturbedAtmosphere::clearPerturbations() {
    T_perturbations_.clear();
    P_perturbations_.clear();
    rho_perturbations_.clear();
}

// =============================================================================
// AtmosphereFactory Implementation
// =============================================================================

std::map<std::string, AtmosphereFactory::ModelCreator>& AtmosphereFactory::getRegistry() {
    static std::map<std::string, ModelCreator> registry;
    return registry;
}

std::shared_ptr<AtmosphericModelBase> AtmosphereFactory::create(const std::string& name) {
    // Check registry first
    auto& registry = getRegistry();
    auto it = registry.find(name);
    if (it != registry.end()) {
        return it->second();
    }
    
    // Built-in models
    if (name == "US_STANDARD_1976" || name == "US_STANDARD" || name == "US1976") {
        return std::make_shared<USStandardAtmosphere1976>();
    }
    if (name == "ICAO_STANDARD" || name == "ICAO" || name == "ISA") {
        return std::make_shared<ICAOStandardAtmosphere>();
    }
    if (name == "NRLMSISE00" || name == "MSISE00" || name == "MSIS") {
        return std::make_shared<NRLMSISE00Atmosphere>();
    }
    if (name == "TROPICAL") {
        return std::make_shared<ReferenceAtmosphere>(ReferenceAtmosphereType::TROPICAL);
    }
    if (name == "MIDLATITUDE_SUMMER" || name == "MLS") {
        return std::make_shared<ReferenceAtmosphere>(ReferenceAtmosphereType::MIDLATITUDE_SUMMER);
    }
    if (name == "MIDLATITUDE_WINTER" || name == "MLW") {
        return std::make_shared<ReferenceAtmosphere>(ReferenceAtmosphereType::MIDLATITUDE_WINTER);
    }
    if (name == "SUBARCTIC_SUMMER" || name == "SAS") {
        return std::make_shared<ReferenceAtmosphere>(ReferenceAtmosphereType::SUBARCTIC_SUMMER);
    }
    if (name == "SUBARCTIC_WINTER" || name == "SAW") {
        return std::make_shared<ReferenceAtmosphere>(ReferenceAtmosphereType::SUBARCTIC_WINTER);
    }
    if (name == "CUSTOM" || name == "TABULATED") {
        return std::make_shared<TabulatedAtmosphere>();
    }
    if (name == "COMPOSITE") {
        return std::make_shared<CompositeAtmosphere>();
    }
    
    // Default to US Standard
    return std::make_shared<USStandardAtmosphere1976>();
}

std::shared_ptr<AtmosphericModelBase> AtmosphereFactory::createFromConfig(const std::string& filename) {
    AtmosphereConfig config;
    if (config.parse(filename)) {
        return config.createModel();
    }
    return nullptr;
}

std::shared_ptr<AtmosphericModelBase> AtmosphereFactory::createReference(ReferenceAtmosphereType type) {
    return std::make_shared<ReferenceAtmosphere>(type);
}

void AtmosphereFactory::registerModel(const std::string& name, ModelCreator creator) {
    getRegistry()[name] = creator;
}

std::vector<std::string> AtmosphereFactory::getAvailableModels() {
    std::vector<std::string> models = {
        "US_STANDARD_1976",
        "ICAO_STANDARD",
        "NRLMSISE00",
        "TROPICAL",
        "MIDLATITUDE_SUMMER",
        "MIDLATITUDE_WINTER",
        "SUBARCTIC_SUMMER",
        "SUBARCTIC_WINTER",
        "CUSTOM",
        "COMPOSITE"
    };
    
    // Add registered custom models
    for (const auto& [name, creator] : getRegistry()) {
        models.push_back(name);
    }
    
    return models;
}

// =============================================================================
// AtmosphereConfig Implementation
// =============================================================================

bool AtmosphereConfig::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    std::map<std::string, std::string> config_map;
    std::string line;
    
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        size_t end = line.find_last_not_of(" \t");
        line = line.substr(start, end - start + 1);
        
        // Parse key = value
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;
        
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);
        
        // Trim key and value
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        
        config_map[key] = value;
    }
    
    return parse(config_map);
}

bool AtmosphereConfig::parse(const std::map<std::string, std::string>& config) {
    auto getStr = [&](const std::string& key, const std::string& def) -> std::string {
        auto it = config.find(key);
        return (it != config.end()) ? it->second : def;
    };
    
    auto getDouble = [&](const std::string& key, double def) -> double {
        auto it = config.find(key);
        return (it != config.end()) ? std::stod(it->second) : def;
    };
    
    auto getBool = [&](const std::string& key, bool def) -> bool {
        auto it = config.find(key);
        if (it == config.end()) return def;
        std::string v = it->second;
        return (v == "true" || v == "1" || v == "yes" || v == "on");
    };
    
    auto getInt = [&](const std::string& key, int def) -> int {
        auto it = config.find(key);
        return (it != config.end()) ? std::stoi(it->second) : def;
    };
    
    // Model selection
    model_type = getStr("atmosphere_model", "US_STANDARD_1976");
    custom_file = getStr("atmosphere_file", "");
    
    // Reference values
    sea_level_pressure = getDouble("sea_level_pressure", AtmosphereConstants::P0_STD);
    sea_level_temperature = getDouble("sea_level_temperature", AtmosphereConstants::T0_STD);
    sea_level_density = getDouble("sea_level_density", AtmosphereConstants::RHO0_STD);
    
    // Domain
    min_altitude = getDouble("min_altitude", 0.0);
    max_altitude = getDouble("max_altitude", 86000.0);
    
    // Wind
    enable_wind = getBool("enable_wind", false);
    wind_config_file = getStr("wind_config", "");
    if (enable_wind) {
        std::string wind_type = getStr("wind_type", "CONSTANT");
        if (wind_type == "NONE") wind_profile.type = WindProfile::Type::NONE;
        else if (wind_type == "CONSTANT") wind_profile.type = WindProfile::Type::CONSTANT;
        else if (wind_type == "LOGARITHMIC") wind_profile.type = WindProfile::Type::LOGARITHMIC;
        else if (wind_type == "POWER_LAW") wind_profile.type = WindProfile::Type::POWER_LAW;
        else if (wind_type == "EKMAN") wind_profile.type = WindProfile::Type::EKMAN_SPIRAL;
        
        wind_profile.reference_speed = getDouble("wind_speed", 0.0);
        wind_profile.reference_direction = getDouble("wind_direction", 0.0);
        wind_profile.reference_height = getDouble("wind_reference_height", 10.0);
        wind_profile.power_exponent = getDouble("wind_power_exponent", 0.14);
    }
    
    // Humidity
    enable_humidity = getBool("enable_humidity", false);
    humidity_config_file = getStr("humidity_config", "");
    if (enable_humidity) {
        std::string hum_type = getStr("humidity_type", "EXPONENTIAL");
        if (hum_type == "NONE") humidity_profile.type = HumidityProfile::Type::NONE;
        else if (hum_type == "CONSTANT") humidity_profile.type = HumidityProfile::Type::CONSTANT_RH;
        else if (hum_type == "EXPONENTIAL") humidity_profile.type = HumidityProfile::Type::EXPONENTIAL;
        
        humidity_profile.surface_relative_humidity = getDouble("surface_rh", 0.7);
        humidity_profile.humidity_scale_height = getDouble("humidity_scale_height", 2500.0);
        humidity_profile.constant_relative_humidity = getDouble("constant_rh", 0.5);
    }
    
    // Turbulence
    enable_turbulence = getBool("enable_turbulence", false);
    turbulence_config_file = getStr("turbulence_config", "");
    if (enable_turbulence) {
        std::string turb_type = getStr("turbulence_type", "BOUNDARY_LAYER");
        if (turb_type == "NONE") turbulence_profile.type = TurbulenceProfile::Type::NONE;
        else if (turb_type == "CONSTANT") turbulence_profile.type = TurbulenceProfile::Type::CONSTANT;
        else if (turb_type == "BOUNDARY_LAYER") turbulence_profile.type = TurbulenceProfile::Type::BOUNDARY_LAYER;
        
        turbulence_profile.boundary_layer_height = getDouble("pbl_height", 1000.0);
        turbulence_profile.friction_velocity = getDouble("friction_velocity", 0.3);
        turbulence_profile.convective_velocity = getDouble("convective_velocity", 0.0);
    }
    
    // NRLMSISE-00 specific
    f107 = getDouble("f107", 150.0);
    f107a = getDouble("f107a", 150.0);
    ap = getDouble("ap", 4.0);
    year = getInt("year", 2000);
    day_of_year = getInt("day_of_year", 172);
    local_time = getDouble("local_time", 12.0);
    latitude = getDouble("latitude", 45.0);
    longitude = getDouble("longitude", 0.0);
    
    return true;
}

std::shared_ptr<AtmosphericModelBase> AtmosphereConfig::createModel() const {
    std::shared_ptr<AtmosphericModelBase> model;
    
    if (model_type == "NRLMSISE00" || model_type == "MSISE00") {
        auto nrl = std::make_shared<NRLMSISE00Atmosphere>();
        nrl->setF107(f107);
        nrl->setF107Average(f107a);
        nrl->setAp(ap);
        nrl->setDateTime(year, day_of_year, local_time * 3600.0);
        nrl->setLocation(latitude, longitude);
        model = nrl;
    } else if (model_type == "CUSTOM" || model_type == "TABULATED") {
        auto tab = std::make_shared<TabulatedAtmosphere>();
        if (!custom_file.empty()) {
            tab->loadFromConfig(custom_file);
        }
        model = tab;
    } else {
        model = AtmosphereFactory::create(model_type);
    }
    
    if (!model) {
        model = std::make_shared<USStandardAtmosphere1976>();
    }
    
    // Apply profiles
    if (enable_wind) {
        model->setWindProfile(wind_profile);
    }
    
    if (enable_humidity) {
        model->setHumidityProfile(humidity_profile);
    }
    
    if (enable_turbulence) {
        model->setTurbulenceProfile(turbulence_profile);
    }
    
    return model;
}

} // namespace FSRM
