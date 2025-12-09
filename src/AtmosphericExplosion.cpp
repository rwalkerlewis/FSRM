/**
 * @file AtmosphericExplosion.cpp
 * @brief Implementation of atmospheric explosion physics with instabilities
 * 
 * See AtmosphericExplosion.hpp for detailed documentation.
 */

#include "AtmosphericExplosion.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <numeric>

namespace FSRM {

using namespace AtmosphericExplosionConstants;

// =============================================================================
// ExtendedAtmosphericModel Implementation
// =============================================================================

ExtendedAtmosphericModel::ExtendedAtmosphericModel() {
    setStandardAtmosphere();
}

void ExtendedAtmosphericModel::setStandardAtmosphere() {
    // US Standard Atmosphere 1976 breakpoints
    altitude_levels_ = {0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 86000.0};
    temperature_levels_ = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.95};
    lapse_rates_ = {-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002};
    
    // Default wind profile (logarithmic in boundary layer, then geostrophic)
    wind_func_ = [](double z, double& u, double& v) {
        const double u_ref = 10.0;  // 10 m/s reference
        const double z_ref = 10.0;  // 10 m reference height
        const double z0 = 0.1;      // Roughness length
        
        if (z < 1000.0) {
            // Logarithmic boundary layer
            u = u_ref * std::log(z / z0 + 1.0) / std::log(z_ref / z0 + 1.0);
        } else if (z < 10000.0) {
            // Ekman spiral transition
            double f = (z - 1000.0) / 9000.0;
            double u_bl = u_ref * std::log(1000.0 / z0 + 1.0) / std::log(z_ref / z0 + 1.0);
            u = u_bl * (1.0 + 0.5 * f);
        } else {
            // Geostrophic + jet stream structure
            double u_geo = 15.0;
            double jet_peak = std::exp(-std::pow((z - 10000.0) / 3000.0, 2));
            u = u_geo + 35.0 * jet_peak;
        }
        v = 0.0;
    };
    
    // Default humidity profile
    humidity_func_ = [](double z) {
        // Exponential decrease with height
        double RH_surface = 0.70;
        double scale_height = 2500.0;
        return RH_surface * std::exp(-z / scale_height);
    };
}

void ExtendedAtmosphericModel::loadFromProfile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open atmosphere profile: " << filename << std::endl;
        return;
    }
    
    altitude_levels_.clear();
    temperature_levels_.clear();
    lapse_rates_.clear();
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        double z, T;
        if (iss >> z >> T) {
            altitude_levels_.push_back(z);
            temperature_levels_.push_back(T);
        }
    }
    
    // Compute lapse rates between levels
    for (size_t i = 0; i < altitude_levels_.size() - 1; ++i) {
        double dT = temperature_levels_[i+1] - temperature_levels_[i];
        double dz = altitude_levels_[i+1] - altitude_levels_[i];
        lapse_rates_.push_back(dT / dz);
    }
}

void ExtendedAtmosphericModel::setCustomProfile(const std::vector<double>& altitudes,
                                                const std::vector<double>& temperatures,
                                                const std::vector<double>& /*pressures*/) {
    altitude_levels_ = altitudes;
    temperature_levels_ = temperatures;
    
    lapse_rates_.clear();
    for (size_t i = 0; i < altitude_levels_.size() - 1; ++i) {
        double dT = temperature_levels_[i+1] - temperature_levels_[i];
        double dz = altitude_levels_[i+1] - altitude_levels_[i];
        lapse_rates_.push_back(dT / dz);
    }
}

void ExtendedAtmosphericModel::setWindProfile(std::function<void(double z, double& u, double& v)> wind_func) {
    wind_func_ = wind_func;
}

void ExtendedAtmosphericModel::setLogWindProfile(double u_ref, double z_ref, double z0) {
    wind_func_ = [u_ref, z_ref, z0](double z, double& u, double& v) {
        if (z <= 0) z = 0.1;
        u = u_ref * std::log((z + z0) / z0) / std::log((z_ref + z0) / z0);
        v = 0.0;
    };
}

void ExtendedAtmosphericModel::setJetStream(double jet_altitude, double jet_speed, double jet_width) {
    has_jet_stream_ = true;
    jet_altitude_ = jet_altitude;
    jet_speed_ = jet_speed;
    jet_width_ = jet_width;
}

double ExtendedAtmosphericModel::temperature(double z) const {
    return interpolateTemperature(z);
}

double ExtendedAtmosphericModel::interpolateTemperature(double z) const {
    if (z <= altitude_levels_.front()) {
        return temperature_levels_.front();
    }
    if (z >= altitude_levels_.back()) {
        return temperature_levels_.back();
    }
    
    // Find layer
    size_t i = 0;
    while (i < altitude_levels_.size() - 1 && altitude_levels_[i+1] < z) {
        ++i;
    }
    
    double dz = z - altitude_levels_[i];
    return temperature_levels_[i] + lapse_rates_[i] * dz;
}

double ExtendedAtmosphericModel::pressure(double z) const {
    return computePressure(z);
}

double ExtendedAtmosphericModel::computePressure(double z) const {
    // Integrate hydrostatic equation through layers
    double P = SEA_LEVEL_PRESSURE;
    double T, z_base = 0.0;
    
    for (size_t i = 0; i < altitude_levels_.size() - 1; ++i) {
        double z_top = std::min(z, altitude_levels_[i+1]);
        double dz = z_top - z_base;
        
        if (dz <= 0) break;
        
        double L = lapse_rates_[i];
        double T_base = temperature_levels_[i];
        
        if (std::abs(L) < 1e-10) {
            // Isothermal layer
            P *= std::exp(-GRAVITY * dz / (R_AIR * T_base));
        } else {
            // Linear temperature layer
            double T_top = T_base + L * dz;
            P *= std::pow(T_top / T_base, -GRAVITY / (R_AIR * L));
        }
        
        z_base = z_top;
        if (z_top >= z) break;
    }
    
    return P;
}

double ExtendedAtmosphericModel::density(double z) const {
    double P = pressure(z);
    double T = temperature(z);
    return P / (R_AIR * T);
}

double ExtendedAtmosphericModel::soundSpeed(double z) const {
    double T = temperature(z);
    return std::sqrt(GAMMA_AIR * R_AIR * T);
}

double ExtendedAtmosphericModel::potentialTemperature(double z) const {
    double T = temperature(z);
    double P = pressure(z);
    return T * std::pow(SEA_LEVEL_PRESSURE / P, R_AIR / CP_AIR);
}

double ExtendedAtmosphericModel::buoyancyFrequency(double z) const {
    // Brunt-Väisälä frequency: N = sqrt(g/θ * dθ/dz)
    double dz = 100.0;  // Finite difference interval
    double theta1 = potentialTemperature(z - dz/2);
    double theta2 = potentialTemperature(z + dz/2);
    double dtheta_dz = (theta2 - theta1) / dz;
    
    double theta = potentialTemperature(z);
    double N2 = GRAVITY / theta * dtheta_dz;
    
    return N2 > 0 ? std::sqrt(N2) : 0.0;
}

double ExtendedAtmosphericModel::richardsonNumber(double z) const {
    double N = buoyancyFrequency(z);
    double shear = windShear(z);
    
    if (std::abs(shear) < 1e-10) return 1e10;
    return N * N / (shear * shear);
}

bool ExtendedAtmosphericModel::isStaticallyStable(double z) const {
    return buoyancyFrequency(z) > 0;
}

bool ExtendedAtmosphericModel::isDynamicallyStable(double z) const {
    return richardsonNumber(z) > 0.25;
}

void ExtendedAtmosphericModel::wind(double z, double& u, double& v) const {
    wind_func_(z, u, v);
    
    // Add jet stream if enabled
    if (has_jet_stream_) {
        double jet_factor = std::exp(-std::pow((z - jet_altitude_) / jet_width_, 2));
        u += jet_speed_ * jet_factor;
    }
}

double ExtendedAtmosphericModel::windShear(double z) const {
    double dz = 100.0;
    double u1, v1, u2, v2;
    wind(z - dz/2, u1, v1);
    wind(z + dz/2, u2, v2);
    
    double du_dz = (u2 - u1) / dz;
    double dv_dz = (v2 - v1) / dz;
    
    return std::sqrt(du_dz * du_dz + dv_dz * dv_dz);
}

double ExtendedAtmosphericModel::windDirection(double z) const {
    double u, v;
    wind(z, u, v);
    return std::atan2(u, v) * 180.0 / M_PI;
}

double ExtendedAtmosphericModel::tropopauseHeight() const {
    // Find height where temperature minimum occurs (end of troposphere)
    return TROPOPAUSE_HEIGHT;
}

double ExtendedAtmosphericModel::stratopauseHeight() const {
    return STRATOPAUSE_HEIGHT;
}

double ExtendedAtmosphericModel::boundaryLayerHeight() const {
    // Estimate from Richardson number profile
    for (double z = 100.0; z < 3000.0; z += 50.0) {
        if (richardsonNumber(z) > 0.25) {
            return z;
        }
    }
    return 1000.0;  // Default
}

void ExtendedAtmosphericModel::setHumidityProfile(std::function<double(double)> humidity_func) {
    humidity_func_ = humidity_func;
}

double ExtendedAtmosphericModel::relativeHumidity(double z) const {
    return humidity_func_(z);
}

double ExtendedAtmosphericModel::saturationVaporPressure(double T) const {
    // Clausius-Clapeyron approximation
    return 611.2 * std::exp(17.67 * (T - 273.15) / (T - 29.65));
}

double ExtendedAtmosphericModel::turbulentKineticEnergy(double z) const {
    // Simple parameterization based on shear
    double shear = windShear(z);
    double Ri = richardsonNumber(z);
    
    // TKE ~ (shear)² / stability_function
    double stability_func = std::max(0.1, 1.0 - Ri / 0.25);
    double TKE = 0.1 * shear * shear * stability_func * 100.0;  // m²/s²
    
    return TKE;
}

double ExtendedAtmosphericModel::eddyViscosity(double z) const {
    // Mixing length model
    double l = std::min(100.0, 0.4 * z);  // Mixing length
    double TKE = turbulentKineticEnergy(z);
    return 0.09 * std::sqrt(TKE) * l;
}

// =============================================================================
// MushroomCloudModel Implementation
// =============================================================================

MushroomCloudModel::MushroomCloudModel() {
    atmosphere_ = std::make_shared<ExtendedAtmosphericModel>();
}

void MushroomCloudModel::setYield(double yield_kt) {
    yield_kt_ = yield_kt;
}

void MushroomCloudModel::setBurstHeight(double height_m) {
    burst_height_ = height_m;
}

void MushroomCloudModel::setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm) {
    atmosphere_ = atm;
}

void MushroomCloudModel::setGroundProperties(double dust_availability, double soil_density) {
    dust_availability_ = dust_availability;
    soil_density_ = soil_density;
}

void MushroomCloudModel::enableRayleighTaylor(bool enable) {
    enable_rt_ = enable;
}

void MushroomCloudModel::enableKelvinHelmholtz(bool enable) {
    enable_kh_ = enable;
}

void MushroomCloudModel::enableEntrainment(bool enable) {
    enable_entrainment_ = enable;
}

void MushroomCloudModel::enableCondensation(bool enable) {
    enable_condensation_ = enable;
}

void MushroomCloudModel::initialize() {
    current_time_ = 0.0;
    current_phase_ = CloudPhase::FIREBALL;
    
    // Initial fireball radius (Glasstone-Dolan scaling)
    // R_max = 66 * W^0.4 meters for W in kt
    double R_max = maxFireballRadius();
    
    cloud_center_z_ = burst_height_;
    cloud_radius_ = R_max * 0.1;  // Start small
    cloud_temperature_ = initialCloudTemperature();
    cloud_mass_ = initialCloudMass();
    debris_mass_ = 0.0;
    rise_velocity_ = 0.0;
    
    // Ambient density at burst height
    double rho_ambient = atmosphere_->density(burst_height_);
    cloud_density_ = rho_ambient * 0.01;  // Hot, low density initially
    
    // Initialize stem
    stem_radius_ = 0.0;
    stem_top_z_ = burst_height_;
    
    // Initialize vortex
    vortex_strength_ = 0.0;
    vortex_core_radius_ = R_max * 0.2;
    
    // Initialize RT instability
    rt_state_.density_heavy = rho_ambient;
    rt_state_.density_light = cloud_density_;
    rt_state_.initial_amplitude = R_max * RT_AMPLITUDE_SEED;
    rt_state_.dominant_wavelength = R_max * 0.5;
    
    // Initialize KH instability
    kh_state_.velocity_upper = 0.0;
    kh_state_.velocity_lower = 0.0;
    kh_state_.density_upper = cloud_density_;
    kh_state_.density_lower = rho_ambient;
    kh_state_.shear_layer_thickness = cloud_radius_ * 0.1;
}

void MushroomCloudModel::step(double dt) {
    switch (current_phase_) {
        case CloudPhase::FIREBALL:
            updateFireballPhase(dt);
            break;
        case CloudPhase::RISE:
        case CloudPhase::TOROIDAL_FORMATION:
        case CloudPhase::STEM_DEVELOPMENT:
            updateRisePhase(dt);
            break;
        case CloudPhase::STABILIZATION:
            updateStabilizationPhase(dt);
            break;
        case CloudPhase::DISPERSION:
            updateDispersionPhase(dt);
            break;
    }
    
    current_time_ += dt;
}

void MushroomCloudModel::updateFireballPhase(double dt) {
    // Fireball expansion
    double R_max = maxFireballRadius();
    double t_form = fireballMaxTime();
    
    if (current_time_ < t_form) {
        // Rapid expansion
        cloud_radius_ = R_max * std::pow(current_time_ / t_form, 0.4);
        
        // Temperature decay
        cloud_temperature_ = initialCloudTemperature() * 
                            std::pow(initialFireballRadius() / cloud_radius_, 1.5);
        
        // Density from ideal gas
        double P = atmosphere_->pressure(cloud_center_z_);
        cloud_density_ = P / (R_AIR * cloud_temperature_);
        
    } else {
        // Transition to rise phase
        current_phase_ = CloudPhase::RISE;
        cloud_radius_ = R_max;
        
        // Begin rise with buoyancy
        rise_velocity_ = 10.0 * std::pow(yield_kt_, 0.25);
    }
}

void MushroomCloudModel::updateRisePhase(double dt) {
    // Buoyant rise of thermal
    
    // Buoyancy
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    double g = GRAVITY;
    double buoyancy = g * (rho_ambient - cloud_density_) / cloud_density_;
    
    // Drag
    double Cd = 0.5;  // Sphere drag coefficient
    double A = M_PI * cloud_radius_ * cloud_radius_;
    double drag = 0.5 * Cd * rho_ambient * rise_velocity_ * std::abs(rise_velocity_) * A / cloud_mass_;
    
    // Equation of motion
    double dv_dt = buoyancy - drag;
    rise_velocity_ += dv_dt * dt;
    cloud_center_z_ += rise_velocity_ * dt;
    
    // Entrainment
    if (enable_entrainment_) {
        computeEntrainment(dt);
    }
    
    // Expansion from entrainment and pressure decrease
    double P_ratio = atmosphere_->pressure(cloud_center_z_) / atmosphere_->pressure(cloud_center_z_ - rise_velocity_ * dt);
    cloud_radius_ *= std::pow(P_ratio, -1.0/3.0);
    
    // Temperature evolution (adiabatic + entrainment mixing)
    double T_ambient = atmosphere_->temperature(cloud_center_z_);
    double mixing_rate = ENTRAINMENT_COEFFICIENT * 2.0 * cloud_radius_ * rise_velocity_;
    double dT_entrain = mixing_rate * (T_ambient - cloud_temperature_) / (cloud_radius_ * cloud_radius_ * rise_velocity_ + 1e-10);
    double dT_adiab = -g / CP_AIR * rise_velocity_ * dt;
    cloud_temperature_ += dT_adiab + dT_entrain * dt;
    
    // Update density
    double P = atmosphere_->pressure(cloud_center_z_);
    cloud_density_ = P / (R_AIR * cloud_temperature_);
    
    // Vortex development
    vortex_strength_ += 2.0 * M_PI * cloud_radius_ * buoyancy * dt;
    
    // Update instabilities
    if (enable_rt_ || enable_kh_) {
        computeInstabilities(dt);
    }
    
    // Update stem
    updateStem(dt);
    
    // Phase transitions
    if (current_time_ > fireballMaxTime() * 5.0) {
        current_phase_ = CloudPhase::TOROIDAL_FORMATION;
    }
    
    // Check for stabilization
    double N = atmosphere_->buoyancyFrequency(cloud_center_z_);
    if (buoyancy < 0.01 * g || (N > 0 && rise_velocity_ < 1.0)) {
        current_phase_ = CloudPhase::STABILIZATION;
    }
}

void MushroomCloudModel::updateStabilizationPhase(double dt) {
    // Overshoot and oscillation
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    double buoyancy = GRAVITY * (rho_ambient - cloud_density_) / cloud_density_;
    
    // Damped oscillation
    double N = atmosphere_->buoyancyFrequency(cloud_center_z_);
    double damping = 0.1 * N;
    
    double dv_dt = buoyancy - damping * rise_velocity_;
    rise_velocity_ += dv_dt * dt;
    cloud_center_z_ += rise_velocity_ * dt;
    
    // Horizontal spreading
    cloud_radius_ += std::abs(rise_velocity_) * ENTRAINMENT_COEFFICIENT * dt;
    
    // Continue entrainment
    if (enable_entrainment_) {
        computeEntrainment(dt);
    }
    
    // Transition to dispersion
    if (std::abs(rise_velocity_) < 0.5 && std::abs(buoyancy) < 0.1) {
        current_phase_ = CloudPhase::DISPERSION;
    }
}

void MushroomCloudModel::updateDispersionPhase(double dt) {
    // Horizontal spreading and wind advection
    double u_wind, v_wind;
    atmosphere_->wind(cloud_center_z_, u_wind, v_wind);
    
    // Advection (not tracking x,y center here, but radius expansion)
    
    // Radial spreading from gravity current
    double g_reduced = GRAVITY * std::abs(atmosphere_->density(cloud_center_z_) - cloud_density_) / 
                      atmosphere_->density(cloud_center_z_);
    double H = cloud_radius_ * 0.5;  // Cloud thickness
    double spreading_velocity = std::sqrt(g_reduced * H);
    cloud_radius_ += spreading_velocity * dt;
    
    // Turbulent diffusion
    double Kh = atmosphere_->eddyViscosity(cloud_center_z_);
    cloud_radius_ += std::sqrt(2.0 * Kh * dt);
}

void MushroomCloudModel::computeEntrainment(double dt) {
    // Morton-Taylor-Turner entrainment
    double alpha = ENTRAINMENT_COEFFICIENT;
    
    // Entrainment velocity proportional to rise velocity
    double v_e = alpha * std::abs(rise_velocity_);
    
    // Mass entrainment rate
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    double A_surface = 4.0 * M_PI * cloud_radius_ * cloud_radius_;
    double dm_dt = rho_ambient * v_e * A_surface;
    
    // Update mass
    double dm = dm_dt * dt;
    double m_old = cloud_mass_;
    cloud_mass_ += dm;
    
    // Mix temperature
    double T_ambient = atmosphere_->temperature(cloud_center_z_);
    cloud_temperature_ = (m_old * cloud_temperature_ + dm * T_ambient) / cloud_mass_;
    
    // Update debris mass (only near ground)
    if (cloud_center_z_ < 5000.0 && burst_height_ < 100.0) {
        debris_mass_ += dm * dust_availability_ * 0.1;
    }
}

void MushroomCloudModel::computeInstabilities(double dt) {
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    
    // Rayleigh-Taylor at cloud top (dense ambient over light cloud)
    if (enable_rt_) {
        rt_state_.density_heavy = rho_ambient;
        rt_state_.density_light = cloud_density_;
        
        // Effective acceleration (includes buoyancy)
        double g_eff = GRAVITY + std::abs(rise_velocity_ * rise_velocity_ / cloud_radius_);
        
        // RT mixing zone growth
        if (rt_state_.atwood_number() > ATWOOD_NUMBER_THRESHOLD) {
            // Nonlinear growth: h = α A g t²
            double t_eff = current_time_ - fireballMaxTime();
            if (t_eff > 0) {
                double h_rt = rt_state_.mixing_zone_width(t_eff, g_eff);
                // Mixing increases effective cloud thickness
            }
        }
    }
    
    // Kelvin-Helmholtz at cloud boundary (shear with ambient)
    if (enable_kh_) {
        double u_ambient, v_ambient;
        atmosphere_->wind(cloud_center_z_, u_ambient, v_ambient);
        double wind_speed = std::sqrt(u_ambient * u_ambient + v_ambient * v_ambient);
        
        kh_state_.velocity_upper = wind_speed;
        kh_state_.velocity_lower = rise_velocity_;  // Cloud vertical velocity
        kh_state_.density_upper = rho_ambient;
        kh_state_.density_lower = cloud_density_;
        kh_state_.shear_layer_thickness = cloud_radius_ * 0.1;
        kh_state_.buoyancy_frequency = atmosphere_->buoyancyFrequency(cloud_center_z_);
        
        // KH rollup scale
        if (kh_state_.is_unstable()) {
            double lambda_kh = kh_state_.most_unstable_wavelength();
            // Rollup creates mixing at this scale
        }
    }
}

void MushroomCloudModel::updateStem(double dt) {
    // Stem develops from surface material drawn up by rising cloud
    
    if (burst_height_ < 100.0 && current_phase_ == CloudPhase::RISE) {
        // Updraft velocity at surface
        double v_up = rise_velocity_ * 0.3;  // Reduced velocity in stem
        
        // Stem radius expands with time
        double r_expansion_rate = v_up * 0.1;
        stem_radius_ = std::max(stem_radius_, cloud_radius_ * 0.3);
        stem_radius_ += r_expansion_rate * dt;
        
        // Stem height follows cloud base
        stem_top_z_ = std::max(0.0, cloud_center_z_ - cloud_radius_);
        
        // Entrain debris
        double rho_debris = soil_density_ * dust_availability_;
        double A_stem = M_PI * stem_radius_ * stem_radius_;
        double dm_debris = rho_debris * v_up * A_stem * 0.01 * dt;
        debris_mass_ += dm_debris;
    }
}

CloudPhase MushroomCloudModel::getCurrentPhase() const {
    return current_phase_;
}

double MushroomCloudModel::cloudTopHeight() const {
    return cloud_center_z_ + cloud_radius_;
}

double MushroomCloudModel::cloudBottomHeight() const {
    return std::max(0.0, cloud_center_z_ - cloud_radius_);
}

double MushroomCloudModel::cloudCenterHeight() const {
    return cloud_center_z_;
}

double MushroomCloudModel::cloudRadius() const {
    return cloud_radius_;
}

double MushroomCloudModel::stemRadius() const {
    return stem_radius_;
}

double MushroomCloudModel::stemHeight() const {
    return stem_top_z_;
}

double MushroomCloudModel::riseVelocity() const {
    return rise_velocity_;
}

double MushroomCloudModel::rotationalVelocity() const {
    // From vortex strength: v = Γ / (2πr)
    if (vortex_core_radius_ < 1.0) return 0.0;
    return vortex_strength_ / (2.0 * M_PI * vortex_core_radius_);
}

double MushroomCloudModel::entrainmentRate() const {
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    double v_e = ENTRAINMENT_COEFFICIENT * std::abs(rise_velocity_);
    double A_surface = 4.0 * M_PI * cloud_radius_ * cloud_radius_;
    return rho_ambient * v_e * A_surface;
}

double MushroomCloudModel::buoyancyFlux() const {
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    double g_prime = GRAVITY * (rho_ambient - cloud_density_) / rho_ambient;
    double V = 4.0/3.0 * M_PI * std::pow(cloud_radius_, 3);
    return g_prime * rise_velocity_ * V;
}

double MushroomCloudModel::cloudTemperature() const {
    return cloud_temperature_;
}

double MushroomCloudModel::cloudDensity() const {
    return cloud_density_;
}

double MushroomCloudModel::thermalExcess() const {
    return cloud_temperature_ - atmosphere_->temperature(cloud_center_z_);
}

double MushroomCloudModel::buoyancy() const {
    double rho_ambient = atmosphere_->density(cloud_center_z_);
    return GRAVITY * (rho_ambient - cloud_density_) / cloud_density_;
}

double MushroomCloudModel::rtMixingZoneWidth() const {
    if (!enable_rt_) return 0.0;
    
    double t_eff = current_time_ - fireballMaxTime();
    if (t_eff <= 0) return 0.0;
    
    double g_eff = GRAVITY;
    return rt_state_.mixing_zone_width(t_eff, g_eff);
}

double MushroomCloudModel::khRollupScale() const {
    if (!enable_kh_ || !kh_state_.is_unstable()) return 0.0;
    return kh_state_.most_unstable_wavelength();
}

double MushroomCloudModel::turbulentIntensity() const {
    // Estimate from instability amplitudes
    double rt_contrib = enable_rt_ ? rt_state_.alpha_bubble * std::abs(buoyancy()) : 0.0;
    double kh_contrib = enable_kh_ && kh_state_.is_unstable() ? 
                       0.1 * kh_state_.velocity_difference() : 0.0;
    return std::sqrt(rt_contrib * rt_contrib + kh_contrib * kh_contrib);
}

double MushroomCloudModel::totalMass() const {
    return cloud_mass_;
}

double MushroomCloudModel::debrisMass() const {
    return debris_mass_;
}

double MushroomCloudModel::fireballMaxTime() const {
    // Time to reach maximum fireball radius
    // τ ≈ 0.001 * W^0.4 seconds
    return 0.001 * std::pow(yield_kt_, 0.4);
}

double MushroomCloudModel::riseStartTime() const {
    return fireballMaxTime() * 1.5;
}

double MushroomCloudModel::stabilizationTime() const {
    // Empirical: t_stab ≈ 60 * W^0.25 seconds
    return 60.0 * std::pow(yield_kt_, 0.25);
}

double MushroomCloudModel::stabilizationHeight() const {
    // Depends on yield and atmospheric stability
    // Empirical: H ≈ 10 + 3.5 * W^0.4 km
    return (10.0 + 3.5 * std::pow(yield_kt_, 0.4)) * 1000.0;
}

double MushroomCloudModel::initialFireballRadius() const {
    // Initial radius ~ 15 * W^0.33 meters
    return 15.0 * std::pow(yield_kt_, 1.0/3.0);
}

double MushroomCloudModel::maxFireballRadius() const {
    // Maximum radius ~ 66 * W^0.4 meters (Glasstone-Dolan)
    return 66.0 * std::pow(yield_kt_, 0.4);
}

double MushroomCloudModel::initialCloudTemperature() const {
    // Initial surface temperature ~ 10^8 K (X-ray domain)
    // But drops rapidly; at max radius ~ 8000 K
    return 1.0e7;  // K (interior average)
}

double MushroomCloudModel::initialCloudMass() const {
    // Mass of air in initial fireball volume
    double R = initialFireballRadius();
    double V = 4.0/3.0 * M_PI * R * R * R;
    double rho = atmosphere_->density(burst_height_);
    return rho * V;
}

void MushroomCloudModel::getVerticalProfile(std::vector<double>& z,
                                           std::vector<double>& radius,
                                           std::vector<double>& temperature,
                                           std::vector<double>& density,
                                           std::vector<double>& velocity) const {
    int n_points = 50;
    double z_min = std::max(0.0, cloud_center_z_ - 2.0 * cloud_radius_);
    double z_max = cloud_center_z_ + 2.0 * cloud_radius_;
    
    z.resize(n_points);
    radius.resize(n_points);
    temperature.resize(n_points);
    density.resize(n_points);
    velocity.resize(n_points);
    
    for (int i = 0; i < n_points; ++i) {
        z[i] = z_min + (z_max - z_min) * i / (n_points - 1);
        
        double dz = z[i] - cloud_center_z_;
        double r_frac = std::abs(dz) / cloud_radius_;
        
        if (r_frac < 1.0) {
            // Inside cloud
            double r_local = cloud_radius_ * std::sqrt(1.0 - r_frac * r_frac);
            radius[i] = r_local;
            
            // Gaussian-ish profiles
            double core_factor = std::exp(-r_frac * r_frac * 2.0);
            temperature[i] = atmosphere_->temperature(z[i]) + 
                            (cloud_temperature_ - atmosphere_->temperature(cloud_center_z_)) * core_factor;
            density[i] = atmosphere_->density(z[i]) * (1.0 - 0.5 * core_factor);
            velocity[i] = rise_velocity_ * core_factor;
        } else {
            // Outside cloud
            radius[i] = 0.0;
            temperature[i] = atmosphere_->temperature(z[i]);
            density[i] = atmosphere_->density(z[i]);
            velocity[i] = 0.0;
        }
    }
}

// =============================================================================
// BuoyantPlumeModel Implementation
// =============================================================================

BuoyantPlumeModel::BuoyantPlumeModel() 
    : source_buoyancy_flux_(1e6),
      source_momentum_flux_(1e4),
      source_radius_(100.0),
      source_height_(0.0),
      source_temperature_(1000.0),
      entrainment_model_(EntrainmentModel::MORTON_TAYLOR_TURNER),
      neutral_height_(0.0),
      max_height_(0.0) {
    atmosphere_ = std::make_shared<ExtendedAtmosphericModel>();
}

void BuoyantPlumeModel::setSourceStrength(double buoyancy_flux) {
    source_buoyancy_flux_ = buoyancy_flux;
}

void BuoyantPlumeModel::setSourceMomentum(double momentum_flux) {
    source_momentum_flux_ = momentum_flux;
}

void BuoyantPlumeModel::setSourceRadius(double radius) {
    source_radius_ = radius;
}

void BuoyantPlumeModel::setSourceHeight(double height) {
    source_height_ = height;
}

void BuoyantPlumeModel::setSourceTemperature(double temperature) {
    source_temperature_ = temperature;
}

void BuoyantPlumeModel::setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm) {
    atmosphere_ = atm;
}

void BuoyantPlumeModel::setEntrainmentModel(EntrainmentModel model) {
    entrainment_model_ = model;
}

void BuoyantPlumeModel::setCrosswind(double speed, double direction) {
    crosswind_speed_ = speed;
    crosswind_direction_ = direction;
}

double BuoyantPlumeModel::entrainmentCoefficient(double Ri_local) const {
    // Entrainment coefficient depends on local Richardson number
    double alpha0 = ENTRAINMENT_COEFFICIENT;
    
    switch (entrainment_model_) {
        case EntrainmentModel::MORTON_TAYLOR_TURNER:
            return alpha0;
            
        case EntrainmentModel::PRIESTLEY_BALL:
            // Stronger entrainment for vigorous plumes
            return alpha0 * 1.5;
            
        case EntrainmentModel::BRIGGS:
            // Reduced entrainment with crosswind
            if (crosswind_speed_ > 0.1) {
                return alpha0 * 0.6;
            }
            return alpha0;
            
        default:
            return alpha0;
    }
}

void BuoyantPlumeModel::solve() {
    // Integrate plume equations from source to max height
    z_profile_.clear();
    radius_profile_.clear();
    velocity_profile_.clear();
    temperature_profile_.clear();
    x_trajectory_.clear();
    y_trajectory_.clear();
    
    // Initial state: [M, B, r, x, y]
    // M = w*r² (momentum flux / pi)
    // B = w*r²*g' (buoyancy flux / pi)
    // r = radius
    
    double z = source_height_;
    double r = source_radius_;
    double T = source_temperature_;
    double T_amb = atmosphere_->temperature(z);
    double rho_amb = atmosphere_->density(z);
    double g_prime = GRAVITY * (T - T_amb) / T_amb;
    
    // Initial velocity from buoyancy flux
    double w = std::pow(source_buoyancy_flux_ / (M_PI * r * r * g_prime + 1e-10), 1.0/3.0);
    if (w < 1.0) w = std::sqrt(source_momentum_flux_ / (M_PI * rho_amb * r * r + 1e-10));
    
    double x = 0.0, y = 0.0;
    double dz = 10.0;  // Integration step
    
    // Integrate upward
    while (z < 100000.0 && w > 0.1) {
        z_profile_.push_back(z);
        radius_profile_.push_back(r);
        velocity_profile_.push_back(w);
        temperature_profile_.push_back(T);
        x_trajectory_.push_back(x);
        y_trajectory_.push_back(y);
        
        // Atmospheric properties
        T_amb = atmosphere_->temperature(z);
        rho_amb = atmosphere_->density(z);
        double N = atmosphere_->buoyancyFrequency(z);
        
        // Buoyancy
        g_prime = GRAVITY * (T - T_amb) / T_amb;
        
        // Entrainment
        double alpha = entrainmentCoefficient(N * N * r / (w * w + 1e-10));
        double dr_dz = alpha;
        
        // Momentum equation
        double dw_dz = g_prime / w - 2.0 * alpha * w / r;
        
        // Heat equation (mixing with ambient)
        double dT_dz = -2.0 * alpha * (T - T_amb) / r - GRAVITY / CP_AIR;
        
        // Crosswind trajectory
        double dx_dz = crosswind_speed_ * std::sin(crosswind_direction_ * M_PI / 180.0) / w;
        double dy_dz = crosswind_speed_ * std::cos(crosswind_direction_ * M_PI / 180.0) / w;
        
        // Update
        z += dz;
        r += dr_dz * dz;
        w += dw_dz * dz;
        T += dT_dz * dz;
        x += dx_dz * dz;
        y += dy_dz * dz;
        
        // Check for neutral buoyancy
        if (g_prime < 0 && neutral_height_ == 0) {
            neutral_height_ = z;
        }
        
        // Prevent unphysical values
        if (r < 1.0) r = 1.0;
        if (T < T_amb) T = T_amb;
    }
    
    max_height_ = z;
    if (neutral_height_ == 0) neutral_height_ = max_height_;
}

double BuoyantPlumeModel::plumeRadius(double z) const {
    if (z_profile_.empty()) return 0.0;
    
    // Linear interpolation
    for (size_t i = 0; i < z_profile_.size() - 1; ++i) {
        if (z >= z_profile_[i] && z < z_profile_[i+1]) {
            double f = (z - z_profile_[i]) / (z_profile_[i+1] - z_profile_[i]);
            return radius_profile_[i] + f * (radius_profile_[i+1] - radius_profile_[i]);
        }
    }
    return radius_profile_.back();
}

double BuoyantPlumeModel::plumeVelocity(double z) const {
    if (z_profile_.empty()) return 0.0;
    
    for (size_t i = 0; i < z_profile_.size() - 1; ++i) {
        if (z >= z_profile_[i] && z < z_profile_[i+1]) {
            double f = (z - z_profile_[i]) / (z_profile_[i+1] - z_profile_[i]);
            return velocity_profile_[i] + f * (velocity_profile_[i+1] - velocity_profile_[i]);
        }
    }
    return 0.0;
}

double BuoyantPlumeModel::plumeTemperature(double z) const {
    if (z_profile_.empty()) return atmosphere_->temperature(z);
    
    for (size_t i = 0; i < z_profile_.size() - 1; ++i) {
        if (z >= z_profile_[i] && z < z_profile_[i+1]) {
            double f = (z - z_profile_[i]) / (z_profile_[i+1] - z_profile_[i]);
            return temperature_profile_[i] + f * (temperature_profile_[i+1] - temperature_profile_[i]);
        }
    }
    return atmosphere_->temperature(z);
}

double BuoyantPlumeModel::neutralBuoyancyHeight() const {
    return neutral_height_;
}

double BuoyantPlumeModel::maximumRiseHeight() const {
    return max_height_;
}

// =============================================================================
// TurbulentMixingModel Implementation
// =============================================================================

TurbulentMixingModel::TurbulentMixingModel() {}

double TurbulentMixingModel::mixingLength(double z, double boundary_layer_height) const {
    // von Karman mixing length with BL cap
    double l_vk = karman_constant_ * z;
    double l_max = 0.1 * boundary_layer_height;
    return std::min(l_vk, l_max);
}

double TurbulentMixingModel::smagorinskyViscosity(double strain_rate, double grid_size) const {
    return std::pow(smagorinsky_constant_ * grid_size, 2) * strain_rate;
}

double TurbulentMixingModel::turbulentPrandtl(double Ri) const {
    // Stability-dependent turbulent Prandtl number
    double Pr_t_neutral = TURBULENT_PRANDTL;
    
    if (Ri > 0) {
        // Stable stratification - increased Pr_t
        return Pr_t_neutral * (1.0 + 5.0 * Ri);
    } else {
        // Unstable - decreased Pr_t
        return Pr_t_neutral / (1.0 - 5.0 * Ri);
    }
}

double TurbulentMixingModel::entrainmentVelocity(double w, double Ri_bulk) const {
    // Entrainment velocity parameterization
    double E0 = ENTRAINMENT_COEFFICIENT;
    
    if (Ri_bulk > 0) {
        // Stable interface - reduced entrainment
        return E0 * w / (1.0 + Ri_bulk);
    } else {
        // Unstable - enhanced entrainment
        return E0 * w * (1.0 - Ri_bulk);
    }
}

double TurbulentMixingModel::mixingEfficiency(double Ri, double /*Pe*/) const {
    // Mixing efficiency Γ = B_flux / (B_flux + ε)
    // Typically 0.15-0.25 for stratified turbulence
    double Gamma_max = 0.2;
    
    if (Ri > 1.0) {
        return Gamma_max * std::exp(-(Ri - 1.0));
    } else if (Ri > 0) {
        return Gamma_max;
    } else {
        return Gamma_max * 1.5;  // Enhanced for unstable
    }
}

// =============================================================================
// AtmosphericExplosionSolver Implementation
// =============================================================================

AtmosphericExplosionSolver::AtmosphericExplosionSolver() {
    atmosphere_ = std::make_shared<ExtendedAtmosphericModel>();
}

void AtmosphericExplosionSolver::setYield(double yield_kt) {
    yield_kt_ = yield_kt;
}

void AtmosphericExplosionSolver::setBurstHeight(double height_m) {
    burst_height_ = height_m;
}

void AtmosphericExplosionSolver::setBurstLocation(double x, double y, double z) {
    burst_x_ = x;
    burst_y_ = y;
    burst_z_ = z;
}

void AtmosphericExplosionSolver::setFissionFraction(double fraction) {
    fission_fraction_ = fraction;
}

void AtmosphericExplosionSolver::setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm) {
    atmosphere_ = atm;
}

void AtmosphericExplosionSolver::setTopography(std::shared_ptr<TopographyModel> topo) {
    topography_ = topo;
}

void AtmosphericExplosionSolver::enableMushroomCloud(bool enable) {
    enable_cloud_ = enable;
}

void AtmosphericExplosionSolver::enableInstabilities(bool enable) {
    enable_instabilities_ = enable;
}

void AtmosphericExplosionSolver::enableDebrisTransport(bool enable) {
    enable_debris_ = enable;
}

void AtmosphericExplosionSolver::enableTurbulentMixing(bool enable) {
    enable_mixing_ = enable;
}

void AtmosphericExplosionSolver::enableRadiation(bool enable) {
    enable_radiation_ = enable;
}

void AtmosphericExplosionSolver::initialize() {
    if (enable_cloud_) {
        cloud_model_ = std::make_shared<MushroomCloudModel>();
        cloud_model_->setYield(yield_kt_);
        cloud_model_->setBurstHeight(burst_height_);
        cloud_model_->setAtmosphere(atmosphere_);
        cloud_model_->enableRayleighTaylor(enable_instabilities_);
        cloud_model_->enableKelvinHelmholtz(enable_instabilities_);
        cloud_model_->initialize();
    }
    
    if (enable_debris_) {
        debris_model_ = std::make_shared<DebrisTransportModel>();
        debris_model_->setAtmosphere(atmosphere_);
        if (cloud_model_) {
            debris_model_->setCloudModel(cloud_model_);
        }
        debris_model_->setDefaultNuclearDistribution(yield_kt_, burst_height_ < 100.0);
        debris_model_->initialize();
    }
    
    if (enable_mixing_) {
        mixing_model_ = std::make_shared<TurbulentMixingModel>();
    }
    
    current_time_ = 0.0;
    initialized_ = true;
    
    // Initialize history tracking
    time_history_.clear();
    rt_amplitude_history_.clear();
    kh_amplitude_history_.clear();
    cloud_height_history_.clear();
}

void AtmosphericExplosionSolver::step(double dt) {
    if (!initialized_) {
        initialize();
    }
    
    // Step cloud model
    if (enable_cloud_ && cloud_model_) {
        cloud_model_->step(dt);
    }
    
    // Step debris transport
    if (enable_debris_ && debris_model_) {
        debris_model_->step(dt);
    }
    
    // Coupled physics updates
    updateCoupledPhysics(dt);
    
    current_time_ += dt;
    
    // Record history
    recordHistory();
}

void AtmosphericExplosionSolver::run(double end_time, double dt) {
    if (!initialized_) {
        initialize();
    }
    
    while (current_time_ < end_time) {
        step(dt);
    }
}

void AtmosphericExplosionSolver::updateCoupledPhysics(double /*dt*/) {
    // Update coupling between models
    
    if (cloud_model_ && debris_model_) {
        // Update debris source from cloud
        double debris = cloud_model_->debrisMass();
        // debris_model updates from cloud state automatically via shared pointer
    }
    
    if (cloud_model_ && mixing_model_) {
        // Update turbulent parameters
        double intensity = cloud_model_->turbulentIntensity();
        // Could use intensity to modify mixing rates
    }
}

void AtmosphericExplosionSolver::recordHistory() {
    time_history_.push_back(current_time_);
    
    if (cloud_model_) {
        rt_amplitude_history_.push_back(cloud_model_->rtMixingZoneWidth());
        kh_amplitude_history_.push_back(cloud_model_->khRollupScale());
        cloud_height_history_.push_back(cloud_model_->cloudTopHeight());
    } else {
        rt_amplitude_history_.push_back(0.0);
        kh_amplitude_history_.push_back(0.0);
        cloud_height_history_.push_back(0.0);
    }
}

void AtmosphericExplosionSolver::writeOutput(const std::string& filename) const {
    std::ofstream out(filename);
    
    out << "# Atmospheric Explosion Simulation Output\n";
    out << "# Yield: " << yield_kt_ << " kt\n";
    out << "# Burst Height: " << burst_height_ << " m\n";
    out << "# Time: " << current_time_ << " s\n";
    out << "#\n";
    
    if (cloud_model_) {
        out << "# Cloud State:\n";
        out << "#   Center Height: " << cloud_model_->cloudCenterHeight() << " m\n";
        out << "#   Top Height: " << cloud_model_->cloudTopHeight() << " m\n";
        out << "#   Radius: " << cloud_model_->cloudRadius() << " m\n";
        out << "#   Rise Velocity: " << cloud_model_->riseVelocity() << " m/s\n";
        out << "#   Temperature: " << cloud_model_->cloudTemperature() << " K\n";
        out << "#   Total Mass: " << cloud_model_->totalMass() << " kg\n";
        out << "#   Debris Mass: " << cloud_model_->debrisMass() << " kg\n";
        out << "#   RT Mixing Zone: " << cloud_model_->rtMixingZoneWidth() << " m\n";
        out << "#   KH Rollup Scale: " << cloud_model_->khRollupScale() << " m\n";
    }
    
    out << "#\n# Time History:\n";
    out << "# time(s)  cloud_height(m)  rt_amplitude(m)  kh_scale(m)\n";
    
    for (size_t i = 0; i < time_history_.size(); ++i) {
        out << time_history_[i] << "  " 
            << cloud_height_history_[i] << "  "
            << rt_amplitude_history_[i] << "  "
            << kh_amplitude_history_[i] << "\n";
    }
}

void AtmosphericExplosionSolver::writeCloudProfile(const std::string& filename) const {
    if (!cloud_model_) return;
    
    std::ofstream out(filename);
    
    out << "# Cloud Vertical Profile at t = " << current_time_ << " s\n";
    out << "# z(m)  radius(m)  temperature(K)  density(kg/m³)  velocity(m/s)\n";
    
    std::vector<double> z, r, T, rho, w;
    cloud_model_->getVerticalProfile(z, r, T, rho, w);
    
    for (size_t i = 0; i < z.size(); ++i) {
        out << z[i] << "  " << r[i] << "  " << T[i] << "  " << rho[i] << "  " << w[i] << "\n";
    }
}

void AtmosphericExplosionSolver::getInstabilityGrowthHistory(std::vector<double>& times,
                                                             std::vector<double>& rt_amplitude,
                                                             std::vector<double>& kh_amplitude) const {
    times = time_history_;
    rt_amplitude = rt_amplitude_history_;
    kh_amplitude = kh_amplitude_history_;
}

// =============================================================================
// DebrisTransportModel Implementation
// =============================================================================

DebrisTransportModel::DebrisTransportModel() : total_debris_mass_(0.0) {
    atmosphere_ = std::make_shared<ExtendedAtmosphericModel>();
}

void DebrisTransportModel::setInitialDebrisMass(double mass) {
    total_debris_mass_ = mass;
}

void DebrisTransportModel::setAtmosphere(std::shared_ptr<ExtendedAtmosphericModel> atm) {
    atmosphere_ = atm;
}

void DebrisTransportModel::setCloudModel(std::shared_ptr<MushroomCloudModel> cloud) {
    cloud_ = cloud;
}

void DebrisTransportModel::setDefaultNuclearDistribution(double yield_kt, bool surface_burst) {
    size_classes_.clear();
    
    // Particle size classes based on nuclear fallout data
    // Median diameter ~100-200 μm for surface bursts, ~10-50 μm for airbursts
    
    double median = surface_burst ? 150e-6 : 30e-6;
    double gsd = 2.5;  // Geometric standard deviation
    
    // Create size classes spanning the distribution
    std::vector<double> diameters = {1e-6, 5e-6, 10e-6, 25e-6, 50e-6, 100e-6, 200e-6, 500e-6, 1000e-6};
    
    for (size_t i = 0; i < diameters.size() - 1; ++i) {
        ParticleSizeClass pc;
        pc.diameter_min = diameters[i];
        pc.diameter_max = diameters[i+1];
        pc.median_diameter = std::sqrt(diameters[i] * diameters[i+1]);
        pc.density = surface_burst ? 2500.0 : 1500.0;  // Soil vs condensate
        
        // Log-normal mass fraction
        double ln_d = std::log(pc.median_diameter);
        double ln_median = std::log(median);
        double ln_sigma = std::log(gsd);
        pc.mass_fraction = std::exp(-0.5 * std::pow((ln_d - ln_median) / ln_sigma, 2));
        
        // Settling velocity (Stokes for small, modified for large)
        double mu_air = 1.8e-5;  // Air viscosity
        double rho_air = 1.2;
        if (pc.median_diameter < 50e-6) {
            pc.settling_velocity = (pc.density - rho_air) * GRAVITY * pc.median_diameter * pc.median_diameter / (18.0 * mu_air);
        } else {
            // Larger particles - use drag coefficient
            double Re_est = 100.0;
            double Cd = 24.0 / Re_est + 6.0 / (1.0 + std::sqrt(Re_est)) + 0.4;
            pc.settling_velocity = std::sqrt(4.0 * pc.median_diameter * (pc.density - rho_air) * GRAVITY / (3.0 * Cd * rho_air));
        }
        
        // Radioactivity
        double fission_yield = yield_kt * 0.5;  // Assume 50% fission
        double total_activity = fission_yield * 1.45e23 * 3.7e10 / std::pow(86400.0, 1.2);  // Way-Wigner at 1 day
        pc.specific_activity = total_activity / (yield_kt * 1e6);  // Bq/kg (rough)
        pc.gamma_energy = 0.7;  // MeV average
        
        size_classes_.push_back(pc);
    }
    
    // Normalize mass fractions
    double total_frac = 0.0;
    for (auto& pc : size_classes_) total_frac += pc.mass_fraction;
    for (auto& pc : size_classes_) pc.mass_fraction /= total_frac;
    
    // Set total mass
    if (surface_burst) {
        total_debris_mass_ = yield_kt * 200.0 * 1000.0;  // ~200 tons per kt for surface
    } else {
        total_debris_mass_ = yield_kt * 10.0 * 1000.0;   // ~10 tons per kt for airburst
    }
}

void DebrisTransportModel::initialize() {
    // Set up height bins
    height_bins_.clear();
    for (double z = 0.0; z <= 50000.0; z += 500.0) {
        height_bins_.push_back(z);
    }
    
    // Initialize particle mass arrays
    particle_mass_.resize(size_classes_.size());
    for (size_t c = 0; c < size_classes_.size(); ++c) {
        particle_mass_[c].resize(height_bins_.size(), 0.0);
    }
    
    // Set up ground grid
    double grid_size = 200000.0;  // 200 km domain
    double dx = 1000.0;           // 1 km resolution
    int nx = static_cast<int>(grid_size / dx);
    
    x_grid_.resize(nx);
    y_grid_.resize(nx);
    for (int i = 0; i < nx; ++i) {
        x_grid_[i] = -grid_size/2 + i * dx;
        y_grid_[i] = -grid_size/2 + i * dx;
    }
    
    ground_deposition_.resize(nx, std::vector<double>(nx, 0.0));
    
    current_time_ = 0.0;
}

void DebrisTransportModel::step(double dt) {
    if (!cloud_) return;
    
    // Get cloud state
    double cloud_center = cloud_->cloudCenterHeight();
    double cloud_radius = cloud_->cloudRadius();
    double cloud_bottom = cloud_->cloudBottomHeight();
    double cloud_top = cloud_->cloudTopHeight();
    
    // Distribute particles from cloud
    if (current_time_ < 300.0) {  // Active release period
        double release_rate = total_debris_mass_ / 300.0;  // Release over 5 minutes
        
        for (size_t c = 0; c < size_classes_.size(); ++c) {
            double mass_release = release_rate * size_classes_[c].mass_fraction * dt;
            
            // Find height bins within cloud
            for (size_t z = 0; z < height_bins_.size(); ++z) {
                if (height_bins_[z] >= cloud_bottom && height_bins_[z] <= cloud_top) {
                    particle_mass_[c][z] += mass_release / ((cloud_top - cloud_bottom) / 500.0);
                }
            }
        }
    }
    
    // Settle particles
    for (size_t c = 0; c < size_classes_.size(); ++c) {
        double v_settle = computeSettlingVelocity(size_classes_[c], cloud_center);
        double dz_settle = v_settle * dt;
        
        // Move mass downward
        for (size_t z = 1; z < height_bins_.size(); ++z) {
            size_t z_src = z;
            while (z_src > 0 && (height_bins_[z] - height_bins_[z_src]) < dz_settle) {
                --z_src;
            }
            
            if (z_src < height_bins_.size() - 1 && z_src >= z) {
                particle_mass_[c][z] += particle_mass_[c][z_src] * (dz_settle / 500.0);
                particle_mass_[c][z_src] *= (1.0 - dz_settle / 500.0);
            }
        }
        
        // Ground deposition from lowest layer
        if (particle_mass_[c][0] > 0) {
            // Deposit in pattern based on wind
            double u, v;
            atmosphere_->wind(100.0, u, v);
            double wind_speed = std::sqrt(u*u + v*v);
            double wind_dir = std::atan2(u, v);
            
            // Simple downwind deposition
            double downwind_dist = wind_speed * current_time_ * 0.1;
            int ix = static_cast<int>((downwind_dist * std::sin(wind_dir) + x_grid_.back()) / 1000.0);
            int iy = static_cast<int>((downwind_dist * std::cos(wind_dir) + y_grid_.back()) / 1000.0);
            
            ix = std::max(0, std::min(ix, static_cast<int>(x_grid_.size()) - 1));
            iy = std::max(0, std::min(iy, static_cast<int>(y_grid_.size()) - 1));
            
            ground_deposition_[ix][iy] += particle_mass_[c][0] * v_settle * dt / (1000.0 * 1000.0);
            particle_mass_[c][0] *= std::exp(-v_settle * dt / 500.0);
        }
    }
    
    current_time_ += dt;
}

double DebrisTransportModel::computeSettlingVelocity(const ParticleSizeClass& pc, double z) const {
    double rho_air = atmosphere_->density(z);
    double T = atmosphere_->temperature(z);
    double mu = 1.458e-6 * std::pow(T, 1.5) / (T + 110.4);  // Sutherland's law
    
    double d = pc.median_diameter;
    double rho_p = pc.density;
    
    // Stokes settling for small particles
    if (d < 50e-6) {
        return (rho_p - rho_air) * GRAVITY * d * d / (18.0 * mu);
    }
    
    // Iterative solution for larger particles
    double v = pc.settling_velocity;  // Initial guess
    for (int iter = 0; iter < 10; ++iter) {
        double Re = rho_air * v * d / mu;
        double Cd;
        if (Re < 1.0) {
            Cd = 24.0 / Re;
        } else if (Re < 1000.0) {
            Cd = 24.0 / Re * (1.0 + 0.15 * std::pow(Re, 0.687));
        } else {
            Cd = 0.44;
        }
        v = std::sqrt(4.0 * d * (rho_p - rho_air) * GRAVITY / (3.0 * Cd * rho_air));
    }
    
    return v;
}

double DebrisTransportModel::totalSuspendedMass() const {
    double total = 0.0;
    for (const auto& class_mass : particle_mass_) {
        for (double m : class_mass) {
            total += m;
        }
    }
    return total;
}

// =============================================================================
// Configuration Parser Implementation
// =============================================================================

bool AtmosphericExplosionConfig::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open config file: " << filename << std::endl;
        return false;
    }
    
    std::map<std::string, std::string> config;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '[') continue;
        
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        
        std::string key = line.substr(0, eq);
        std::string value = line.substr(eq + 1);
        
        // Trim
        while (!key.empty() && (key.back() == ' ' || key.back() == '\t')) key.pop_back();
        while (!value.empty() && (value.front() == ' ' || value.front() == '\t')) value.erase(0, 1);
        
        config[key] = value;
    }
    
    return parse(config);
}

bool AtmosphericExplosionConfig::parse(const std::map<std::string, std::string>& config) {
    auto getDouble = [&](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end()) {
            try { return std::stod(it->second); } catch (...) {}
        }
        return def;
    };
    
    auto getBool = [&](const std::string& key, bool def) {
        auto it = config.find(key);
        if (it != config.end()) {
            return it->second == "true" || it->second == "1" || it->second == "yes";
        }
        return def;
    };
    
    yield_kt_ = getDouble("yield_kt", 100.0);
    burst_height_ = getDouble("burst_height", 0.0);
    burst_x_ = getDouble("location_x", 0.0);
    burst_y_ = getDouble("location_y", 0.0);
    burst_z_ = getDouble("location_z", burst_height_);
    fission_fraction_ = getDouble("fission_fraction", 0.5);
    
    enable_cloud_ = getBool("enable_mushroom_cloud", true);
    enable_instabilities_ = getBool("enable_instabilities", true);
    enable_debris_ = getBool("enable_debris_transport", true);
    enable_mixing_ = getBool("enable_turbulent_mixing", true);
    enable_radiation_ = getBool("enable_radiation", true);
    
    auto it = config.find("atmosphere_model");
    if (it != config.end()) atmosphere_model_ = it->second;
    
    configured_ = true;
    return true;
}

std::unique_ptr<AtmosphericExplosionSolver> AtmosphericExplosionConfig::createSolver() const {
    auto solver = std::make_unique<AtmosphericExplosionSolver>();
    
    solver->setYield(yield_kt_);
    solver->setBurstHeight(burst_height_);
    solver->setBurstLocation(burst_x_, burst_y_, burst_z_);
    solver->setFissionFraction(fission_fraction_);
    
    solver->enableMushroomCloud(enable_cloud_);
    solver->enableInstabilities(enable_instabilities_);
    solver->enableDebrisTransport(enable_debris_);
    solver->enableTurbulentMixing(enable_mixing_);
    solver->enableRadiation(enable_radiation_);
    
    auto atm = std::make_shared<ExtendedAtmosphericModel>();
    atm->setStandardAtmosphere();
    solver->setAtmosphere(atm);
    
    return solver;
}

} // namespace FSRM
