/**
 * @file VolcanoModel.cpp
 * @brief Implementation of fully coupled volcano modeling physics
 * 
 * Implements comprehensive multi-physics models for volcanic systems including
 * magma chamber dynamics, conduit flow, eruption columns, PDCs, lava flows,
 * lahars, volcanic deformation, seismicity, gas emissions, and tephra dispersal.
 * 
 * @author FSRM Development Team
 */

#include "VolcanoModel.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>

namespace FSRM {

// =============================================================================
// Physical Constants (local to this file)
// =============================================================================

static const double PI = M_PI;
static const double G = VolcanoConstants::GRAVITY;
static const double R_GAS = VolcanoConstants::GAS_CONSTANT;
static const double STEFAN_BOLTZMANN = VolcanoConstants::STEFAN_BOLTZMANN;
static const double P_ATM = VolcanoConstants::ATMOSPHERIC_PRESSURE;

// =============================================================================
// Utility Functions (VolcanoUtils namespace)
// =============================================================================

namespace VolcanoUtils {

void VEI_to_parameters(int VEI, double& volume_km3, double& column_height_km,
                       double& mass_eruption_rate) {
    // Based on Newhall & Self (1982) VEI scale
    switch (VEI) {
        case 0:
            volume_km3 = 1e-6;
            column_height_km = 0.1;
            mass_eruption_rate = 1e3;
            break;
        case 1:
            volume_km3 = 1e-5;
            column_height_km = 0.5;
            mass_eruption_rate = 1e4;
            break;
        case 2:
            volume_km3 = 1e-3;
            column_height_km = 1.0;
            mass_eruption_rate = 1e5;
            break;
        case 3:
            volume_km3 = 1e-2;
            column_height_km = 5.0;
            mass_eruption_rate = 1e6;
            break;
        case 4:
            volume_km3 = 0.1;
            column_height_km = 10.0;
            mass_eruption_rate = 1e7;
            break;
        case 5:
            volume_km3 = 1.0;
            column_height_km = 20.0;
            mass_eruption_rate = 1e8;
            break;
        case 6:
            volume_km3 = 10.0;
            column_height_km = 30.0;
            mass_eruption_rate = 1e9;
            break;
        case 7:
            volume_km3 = 100.0;
            column_height_km = 40.0;
            mass_eruption_rate = 1e10;
            break;
        case 8:
            volume_km3 = 1000.0;
            column_height_km = 50.0;
            mass_eruption_rate = 1e11;
            break;
        default:
            volume_km3 = 0.1;
            column_height_km = 10.0;
            mass_eruption_rate = 1e7;
    }
}

double columnHeight_from_MER(double mass_eruption_rate) {
    // Mastin et al. (2009): H = 2.0 * (MER)^0.241
    // H in km, MER in kg/s
    if (mass_eruption_rate <= 0) return 0.0;
    return 2.0 * std::pow(mass_eruption_rate, 0.241);
}

double MER_from_columnHeight(double height_km) {
    // Inverse of Mastin et al. (2009)
    if (height_km <= 0) return 0.0;
    return std::pow(height_km / 2.0, 1.0 / 0.241);
}

double estimateDuration(double volume_km3, double mass_eruption_rate) {
    // Convert volume to mass (assuming density ~2500 kg/m³)
    double mass_kg = volume_km3 * 1e9 * 2500.0;  // km³ -> m³ -> kg
    return mass_kg / mass_eruption_rate;  // seconds
}

double giordano2008_viscosity(double T_celsius, double SiO2, double TiO2,
                              double Al2O3, double FeO, double MnO,
                              double MgO, double CaO, double Na2O,
                              double K2O, double P2O5, double H2O) {
    // Giordano, Russell, Dingwell (2008) viscosity model
    // Valid for 700-1600°C and various melt compositions
    
    // Normalize to 100% anhydrous
    double total = SiO2 + TiO2 + Al2O3 + FeO + MnO + MgO + CaO + Na2O + K2O + P2O5;
    if (total <= 0) total = 100.0;
    
    double norm = 100.0 / total;
    double sio2 = SiO2 * norm;
    double tio2 = TiO2 * norm;
    double al2o3 = Al2O3 * norm;
    double feo = FeO * norm;
    double mgo = MgO * norm;
    double cao = CaO * norm;
    double na2o = Na2O * norm;
    double k2o = K2O * norm;
    
    // Model parameters (simplified VFT form)
    // log10(η) = A + B / (T - C)
    
    // A parameter
    double A = -4.55;
    
    // B parameter (composition-dependent with TiO2 and FeO effects)
    double b1 = 159.6 - 8.43 * (na2o + k2o);
    double b2 = 1.0 - 0.01 * sio2 + 0.005 * tio2;  // TiO2 increases viscosity slightly
    double b3 = 1.0 + 0.03 * (al2o3);
    double b4 = 1.0 - 0.02 * (cao + mgo) - 0.015 * feo;  // FeO reduces viscosity (network modifier)
    double b5 = 1.0 - 0.5 * H2O;  // Water effect
    
    double B = b1 * b2 * b3 * b4 * b5 * 100.0;
    
    // C parameter (glass transition)
    double C = 150.0 + 10.0 * sio2 / 100.0 - 50.0 * H2O;
    
    // Temperature in Kelvin
    double T_K = T_celsius + 273.15;
    
    // Avoid singularity near glass transition
    if (T_K <= C + 10) {
        T_K = C + 10;
    }
    
    double log_eta = A + B / (T_K - C);
    
    return std::pow(10.0, log_eta);
}

double H2O_solubility(double P_MPa, double T_celsius, double SiO2_wt) {
    // Newman & Lowenstern (2002) - simplified power law
    // H2O solubility in wt%
    
    if (P_MPa <= 0) return 0.0;
    
    // Coefficients depend on composition
    double a, b;
    if (SiO2_wt > 70) {
        // Rhyolite
        a = 0.4111;
        b = 0.5;
    } else if (SiO2_wt > 60) {
        // Dacite/Andesite
        a = 0.3539;
        b = 0.5;
    } else {
        // Basalt
        a = 0.2891;
        b = 0.5;
    }
    
    return a * std::pow(P_MPa, b);
}

double CO2_solubility(double P_MPa, double T_celsius, double SiO2_wt) {
    // CO2 solubility in wt% (much lower than H2O)
    
    if (P_MPa <= 0) return 0.0;
    
    // Henry's law approximation
    double K_H;
    if (SiO2_wt > 65) {
        K_H = 5e-6;  // Rhyolite
    } else if (SiO2_wt > 55) {
        K_H = 1e-5;  // Andesite
    } else {
        K_H = 2e-5;  // Basalt
    }
    
    return K_H * P_MPa * 100.0;  // Convert to wt%
}

double phi_to_mm(double phi) {
    // phi = -log2(d_mm)
    return std::pow(2.0, -phi);
}

double mm_to_phi(double mm) {
    if (mm <= 0) return 10.0;  // Very fine
    return -std::log2(mm);
}

double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
    // Great circle distance in meters
    const double R = 6371000.0;  // Earth radius in meters
    const double DEG_TO_RAD = PI / 180.0;
    
    double dLat = (lat2 - lat1) * DEG_TO_RAD;
    double dLon = (lon2 - lon1) * DEG_TO_RAD;
    
    double a = std::sin(dLat/2) * std::sin(dLat/2) +
               std::cos(lat1 * DEG_TO_RAD) * std::cos(lat2 * DEG_TO_RAD) *
               std::sin(dLon/2) * std::sin(dLon/2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1-a));
    
    return R * c;
}

} // namespace VolcanoUtils

// =============================================================================
// MagmaProperties Implementation
// =============================================================================

double MagmaProperties::getMeltViscosity() const {
    // Use Giordano et al. (2008) model
    double T_C = temperature - 273.15;  // Convert to Celsius
    
    return VolcanoUtils::giordano2008_viscosity(
        T_C, SiO2, 0.5,  // TiO2 assumed
        Al2O3, FeO, 0.1,  // MnO assumed
        MgO, CaO, Na2O, K2O, 0.1,  // P2O5 assumed
        H2O_total
    );
}

double MagmaProperties::getBulkViscosity() const {
    // Costa (2005) model for crystal-bearing magma
    double eta_melt = getMeltViscosity();
    
    if (crystal_fraction <= 0) return eta_melt;
    
    // Maximum packing fraction
    double phi_max = 0.6;
    
    // Relative viscosity
    double phi_star = crystal_fraction / phi_max;
    if (phi_star >= 1.0) {
        return eta_melt * 1e10;  // Effectively solid
    }
    
    // Costa (2005) parameterization
    double gamma = 2.5;
    double delta = 1.0;
    double B = 2.5;
    
    double erf_arg = std::sqrt(PI) / (2 * (1 - phi_star)) *
                     crystal_fraction * (1 + std::pow(crystal_fraction, gamma));
    double F = (1 - phi_star) * std::erf(erf_arg);
    
    double eta_r = (1 + std::pow(crystal_fraction, delta)) / 
                   std::pow(1 - F, B * phi_max);
    
    // Account for bubbles (reduce viscosity at low bubble fraction)
    if (bubble_fraction > 0 && bubble_fraction < 0.5) {
        double bubble_factor = 1.0 - 0.5 * bubble_fraction;
        eta_r *= bubble_factor;
    }
    
    return eta_melt * eta_r;
}

double MagmaProperties::getDensity(double pressure) const {
    // Mixture density with pressure-dependent compressibility
    double P_MPa = pressure / 1e6;
    
    // Melt density (temperature and composition dependent)
    double rho_melt;
    switch (composition) {
        case MagmaComposition::BASALT:
            rho_melt = 2700.0 - 0.2 * (temperature - 1473);
            break;
        case MagmaComposition::ANDESITE:
            rho_melt = 2600.0 - 0.2 * (temperature - 1373);
            break;
        case MagmaComposition::DACITE:
            rho_melt = 2500.0 - 0.2 * (temperature - 1273);
            break;
        case MagmaComposition::RHYOLITE:
            rho_melt = 2400.0 - 0.2 * (temperature - 1173);
            break;
        default:
            rho_melt = 2500.0;
    }
    
    // Apply pressure-dependent compressibility
    // β = (1/ρ)(dρ/dP) ≈ 1e-5 /MPa for silicate melts
    double beta = 1.0e-5;  // Compressibility coefficient (1/MPa)
    rho_melt *= (1.0 + beta * P_MPa);
    
    // Crystal density
    double rho_crystal = VolcanoConstants::CRYSTAL_DENSITY;
    
    // Gas density (ideal gas approximation)
    double MW_gas = VolcanoConstants::H2O_MOLECULAR_WEIGHT;  // Assume mostly H2O
    double rho_gas = pressure * MW_gas / (R_GAS * temperature);
    
    // Mixture
    double phi_melt = 1.0 - crystal_fraction - bubble_fraction;
    phi_melt = std::max(0.0, phi_melt);
    
    return phi_melt * rho_melt + 
           crystal_fraction * rho_crystal + 
           bubble_fraction * rho_gas;
}

double MagmaProperties::getDissolvedH2O(double pressure) const {
    double P_MPa = pressure / 1e6;
    double T_C = temperature - 273.15;
    return VolcanoUtils::H2O_solubility(P_MPa, T_C, SiO2);
}

double MagmaProperties::getDissolvedCO2(double pressure) const {
    double P_MPa = pressure / 1e6;
    double T_C = temperature - 273.15;
    return VolcanoUtils::CO2_solubility(P_MPa, T_C, SiO2);
}

double MagmaProperties::getSaturationPressure() const {
    // Iteratively find pressure where dissolved = total volatiles
    double P_low = 0.1e6;   // 0.1 MPa
    double P_high = 500e6;  // 500 MPa
    
    for (int iter = 0; iter < 50; ++iter) {
        double P_mid = 0.5 * (P_low + P_high);
        double dissolved = getDissolvedH2O(P_mid) + getDissolvedCO2(P_mid);
        double total = H2O_total + CO2_total;
        
        if (dissolved < total) {
            P_low = P_mid;
        } else {
            P_high = P_mid;
        }
        
        if (std::abs(P_high - P_low) < 1e5) break;
    }
    
    return 0.5 * (P_low + P_high);
}

double MagmaProperties::getLiquidusTemperature() const {
    // Simplified liquidus based on composition
    // Higher SiO2 = lower liquidus
    double T_base;
    switch (composition) {
        case MagmaComposition::BASALT:
            T_base = 1523.0;  // ~1250°C
            break;
        case MagmaComposition::BASALTIC_ANDESITE:
            T_base = 1473.0;
            break;
        case MagmaComposition::ANDESITE:
            T_base = 1423.0;
            break;
        case MagmaComposition::DACITE:
            T_base = 1323.0;
            break;
        case MagmaComposition::RHYOLITE:
            T_base = 1223.0;  // ~950°C
            break;
        default:
            T_base = 1373.0;
    }
    
    // Water lowers liquidus (~100°C per 4 wt% H2O)
    return T_base - 25.0 * H2O_total;
}

double MagmaProperties::getSolidusTemperature() const {
    // Solidus ~200-300°C below liquidus
    double T_liq = getLiquidusTemperature();
    return T_liq - 250.0 - 20.0 * H2O_total;  // Water also lowers solidus
}

double MagmaProperties::getCrystallizationRate() const {
    // dφ/dT - fraction crystallized per degree cooling
    double T_liq = getLiquidusTemperature();
    double T_sol = getSolidusTemperature();
    double dT = T_liq - T_sol;
    
    if (dT <= 0) return 0.0;
    
    // Linear crystallization as approximation
    return 1.0 / dT;
}

// =============================================================================
// MagmaRheology Implementation
// =============================================================================

MagmaRheology::MagmaRheology() 
    : melt_viscosity(1e4), yield_stress(0.0), consistency_K(1e4), power_law_n(1.0) {}

void MagmaRheology::setMagmaProperties(const MagmaProperties& props) {
    magma = props;
    computeRheologicalParameters();
}

void MagmaRheology::computeRheologicalParameters() {
    melt_viscosity = magma.getMeltViscosity();
    
    // Yield stress depends on crystal fraction (Caricchi et al., 2007)
    if (magma.crystal_fraction > 0.25) {
        // Yield stress emerges above ~25% crystals
        double phi_excess = magma.crystal_fraction - 0.25;
        yield_stress = 1e4 * std::pow(phi_excess / 0.35, 2.0);  // Pa
        yield_stress = std::min(yield_stress, 1e6);  // Cap at 1 MPa
    } else {
        yield_stress = 0.0;
    }
    
    // Herschel-Bulkley parameters
    consistency_K = magma.getBulkViscosity();
    
    // Power-law index (shear-thinning for crystal-rich magmas)
    if (magma.crystal_fraction > 0.3) {
        power_law_n = 0.7 - 0.3 * (magma.crystal_fraction - 0.3) / 0.3;
        power_law_n = std::max(power_law_n, 0.3);
    } else {
        power_law_n = 1.0;  // Newtonian
    }
}

double MagmaRheology::getEffectiveViscosity(double shear_rate) const {
    if (shear_rate <= 1e-10) shear_rate = 1e-10;
    
    if (yield_stress <= 0 && std::abs(power_law_n - 1.0) < 0.01) {
        // Newtonian
        return consistency_K;
    }
    
    // Herschel-Bulkley apparent viscosity
    // τ = τ_y + K * γ̇^n
    // η_app = τ / γ̇ = τ_y/γ̇ + K * γ̇^(n-1)
    
    double eta_app = yield_stress / shear_rate + 
                     consistency_K * std::pow(shear_rate, power_law_n - 1.0);
    
    // Regularize at low shear rates
    double eta_max = consistency_K * 1e6;
    return std::min(eta_app, eta_max);
}

double MagmaRheology::getYieldStress() const {
    return yield_stress;
}

bool MagmaRheology::canFlow(double shear_stress) const {
    return shear_stress > yield_stress;
}

double MagmaRheology::getCrystalViscosityFactor() const {
    if (magma.crystal_fraction <= 0) return 1.0;
    
    double phi_max = 0.6;
    double phi_star = magma.crystal_fraction / phi_max;
    
    if (phi_star >= 0.99) return 1e10;
    
    // Einstein-Roscoe relation
    return std::pow(1.0 - phi_star, -2.5);
}

double MagmaRheology::getBubbleViscosityFactor() const {
    if (magma.bubble_fraction <= 0) return 1.0;
    
    // Low bubble fraction: slight increase
    // High bubble fraction: decrease (capillary number dependent)
    if (magma.bubble_fraction < 0.3) {
        return 1.0 + magma.bubble_fraction;
    } else {
        return 1.3 * (1.0 - magma.bubble_fraction);
    }
}

// =============================================================================
// MagmaChamberModel Implementation
// =============================================================================

MagmaChamberModel::MagmaChamberModel()
    : recharge_rate(0.0), eruption_rate(0.0), wall_heat_flux(0.0) {}

void MagmaChamberModel::initialize(const MagmaChamberGeometry& geom,
                                    const MagmaProperties& mag,
                                    const MagmaChamberState& initial) {
    geometry = geom;
    magma = mag;
    state = initial;
    
    // Ensure consistency
    state.volume = geometry.computeVolume();
    state.mass = state.volume * magma.getDensity(state.pressure);
}

void MagmaChamberModel::setRechargeRate(double rate_kg_s) {
    recharge_rate = rate_kg_s;
}

void MagmaChamberModel::setEruptionRate(double rate_kg_s) {
    eruption_rate = rate_kg_s;
}

void MagmaChamberModel::setWallRockCooling(double heat_flux) {
    wall_heat_flux = heat_flux;
}

void MagmaChamberModel::update(double dt) {
    // Mass balance
    double dm = (recharge_rate - eruption_rate) * dt;
    state.mass += dm;
    
    // Cooling from wall rocks
    double surface_area = 4.0 * PI * std::pow(geometry.effectiveRadius(), 2);
    double dQ = -wall_heat_flux * surface_area * dt;
    
    // Temperature change
    double cp = VolcanoConstants::MAGMA_SPECIFIC_HEAT;
    double dT = dQ / (state.mass * cp);
    
    // Add latent heat from crystallization
    double dphi_x = computeCrystallization(dT, dt);
    double L = VolcanoConstants::LATENT_HEAT_CRYSTALLIZATION;
    dT += L * dphi_x / cp;
    
    state.temperature += dT;
    state.crystal_fraction += dphi_x;
    state.crystal_fraction = std::clamp(state.crystal_fraction, 0.0, 0.7);
    
    // Volatile exsolution
    double dX_gas = computeVolatileExsolution(0, dT);
    state.exsolved_gas_fraction += dX_gas;
    state.exsolved_gas_fraction = std::clamp(state.exsolved_gas_fraction, 0.0, 0.5);
    
    // Pressure change from mass, temperature, and volatile changes
    double dV = dm / magma.getDensity(state.pressure);  // Volume change from mass
    double dP = computePressureChange(dV, dm, dT);
    state.pressure += dP;
    
    // Update overpressure (relative to lithostatic)
    double P_lith = VolcanoConstants::CRUSTAL_DENSITY * G * geometry.depth;
    state.overpressure = state.pressure - P_lith;
}

double MagmaChamberModel::computePressureChange(double dV, double dm, double dT) {
    // Pressure from compressibility
    double beta = getMagmaCompressibility();
    
    // Thermal expansion
    double alpha_T = 1e-4;  // 1/K
    
    // Pressure change
    double dP = -(dV / state.volume) / beta;
    dP += alpha_T * dT / beta;
    
    // Volatile exsolution increases pressure
    dP += state.exsolved_gas_fraction * R_GAS * state.temperature / 
          (VolcanoConstants::H2O_MOLECULAR_WEIGHT * state.volume * beta);
    
    return dP;
}

double MagmaChamberModel::computeVolatileExsolution(double dP, double dT) {
    // Check if volatiles can exsolve
    double P_sat = magma.getSaturationPressure();
    
    if (state.pressure > P_sat) {
        return 0.0;  // Undersaturated
    }
    
    // Amount of exsolution depends on how far below saturation
    double undersaturation = (P_sat - state.pressure) / P_sat;
    double max_exsolution = magma.H2O_total - state.dissolved_H2O;
    
    return std::min(undersaturation * 0.1, max_exsolution * 0.01);
}

double MagmaChamberModel::computeCrystallization(double dT, double dt) {
    if (dT >= 0) return 0.0;  // Heating, no crystallization
    
    // Check if below liquidus
    double T_liq = magma.getLiquidusTemperature();
    if (state.temperature > T_liq) return 0.0;
    
    // Crystallization rate
    double rate = magma.getCrystallizationRate();
    double dphi = -rate * dT;  // dT is negative for cooling
    
    // Limit by remaining melt
    double max_crystals = 0.7 - state.crystal_fraction;
    return std::min(dphi, max_crystals);
}

bool MagmaChamberModel::isNearFailure(double tensile_strength) const {
    return state.overpressure > 0.8 * tensile_strength;
}

double MagmaChamberModel::getFailureProbability() const {
    // Simple probability based on overpressure
    double tensile_strength = 10e6;  // 10 MPa typical
    double ratio = state.overpressure / tensile_strength;
    
    if (ratio <= 0) return 0.0;
    if (ratio >= 1) return 1.0;
    
    // Sigmoid probability
    return 1.0 / (1.0 + std::exp(-10.0 * (ratio - 0.7)));
}

void MagmaChamberModel::computeSurfaceDeformation(double x, double y,
                                                   double& u_radial, double& u_vertical) const {
    // Compute volume change from overpressure
    double beta = getMagmaCompressibility();
    double dV = beta * state.overpressure * state.volume;
    
    mogiSource(x, y, geometry.depth, dV, u_radial, u_vertical);
}

void MagmaChamberModel::mogiSource(double x, double y, double depth, double dV,
                                    double& ur, double& uz) const {
    // Mogi (1958) point source deformation with elastic parameters
    double r = std::sqrt(x*x + y*y);
    double R = std::sqrt(r*r + depth*depth);
    
    double nu = 0.25;  // Poisson's ratio
    double mu = VolcanoConstants::SHEAR_MODULUS_CRUST;  // Shear modulus (Pa)
    
    // Use full elastic solution: C = (1-ν)/(πμ) * ΔP * V
    // For point source: C = (1-ν) * dV / π
    double C = (1 - nu) * dV / PI;
    
    // Scale by elastic moduli for stress-dependent effects
    // In reality, mu affects the pressure-volume relationship
    // For typical crustal conditions: mu ≈ 10-40 GPa
    double elastic_factor = 1.0;  // Could be: 1.0 / (1.0 + mu/1e10) for softer response
    
    ur = C * r / std::pow(R, 3) * elastic_factor;
    uz = C * depth / std::pow(R, 3) * elastic_factor;
}

void MagmaChamberModel::mctigueSource(double x, double y, double depth, double radius,
                                       double dP, double& ur, double& uz) const {
    // McTigue (1987) finite spherical source
    double r = std::sqrt(x*x + y*y);
    double R = std::sqrt(r*r + depth*depth);
    
    double nu = 0.25;
    double mu = VolcanoConstants::SHEAR_MODULUS_CRUST;
    
    double a = radius;
    double a_R = a / R;
    
    // First-order correction for finite source
    double C = (1 - nu) * a * a * a * dP / (mu * std::pow(R, 3));
    C *= (1 + a_R * a_R / 2.0);  // McTigue correction
    
    ur = C * r / R;
    uz = C * depth / R;
}

double MagmaChamberModel::getSeismicMomentRate() const {
    // Seismic moment rate from volume change rate
    double mu = VolcanoConstants::SHEAR_MODULUS_CRUST;
    double beta = getMagmaCompressibility();
    
    // dV/dt from pressure rate (approximate)
    double dP_dt = state.overpressure / 86400.0;  // Assume builds over 1 day
    double dV_dt = beta * dP_dt * state.volume;
    
    // M0_dot = mu * dV_dt (for volumetric source)
    return mu * std::abs(dV_dt);
}

double MagmaChamberModel::getMagmaCompressibility() const {
    // Effective compressibility of bubbly magma
    double beta_melt = 1e-10;  // 1/Pa for silicate melt
    double beta_gas = 1.0 / state.pressure;  // Ideal gas
    
    // Weighted average
    double phi_gas = state.exsolved_gas_fraction;
    return (1 - phi_gas) * beta_melt + phi_gas * beta_gas;
}

// =============================================================================
// ConduitFlowModel Implementation
// =============================================================================

ConduitFlowModel::ConduitFlowModel()
    : mass_flux(0.0), fragmentation_depth(0.0), n_nodes(100), dz(50.0) {}

void ConduitFlowModel::initialize(const ConduitGeometry& geom,
                                   const MagmaProperties& mag,
                                   double chamber_pressure) {
    geometry = geom;
    magma = mag;
    
    n_nodes = 100;
    dz = geometry.length / n_nodes;
    
    // Initialize profile
    profile.resize(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        profile[i].z = i * dz;
        profile[i].pressure = chamber_pressure * (1.0 - static_cast<double>(i) / n_nodes);
        profile[i].temperature = magma.temperature;
        profile[i].velocity = 1.0;
        profile[i].density = magma.getDensity(profile[i].pressure);
        profile[i].melt_fraction = 1.0 - magma.crystal_fraction;
        profile[i].crystal_fraction = magma.crystal_fraction;
        profile[i].gas_fraction = 0.0;
        profile[i].regime = ConduitFlowRegime::BUBBLY;
    }
    
    // Set boundary pressures
    profile[0].pressure = chamber_pressure;
    profile[n_nodes - 1].pressure = P_ATM;
}

void ConduitFlowModel::setChamberPressure(double P) {
    if (!profile.empty()) {
        profile[0].pressure = P;
    }
}

void ConduitFlowModel::setVentPressure(double P) {
    if (!profile.empty()) {
        profile.back().pressure = P;
    }
}

void ConduitFlowModel::setMassFluxTarget(double mdot) {
    mass_flux = mdot;
}

bool ConduitFlowModel::solvesteadyState() {
    // Iterative solution for steady-state conduit flow
    // Uses shooting method from chamber to surface
    
    double P_chamber = profile[0].pressure;
    double T = magma.temperature;
    
    // Initial guess for velocity at chamber
    double u_init = 1.0;  // m/s
    
    for (int iter = 0; iter < 50; ++iter) {
        double P = P_chamber;
        double u = u_init;
        double rho = magma.getDensity(P);
        
        fragmentation_depth = 0.0;
        bool fragmented = false;
        
        for (int i = 0; i < n_nodes; ++i) {
            double z = i * dz;
            double A = geometry.areaAt(z);
            
            // Volatile exsolution
            double H2O_sat = H2O_solubility(P, T);
            double H2O_ex = std::max(0.0, magma.H2O_total - H2O_sat);
            
            // Gas fraction
            double rho_gas = P * VolcanoConstants::H2O_MOLECULAR_WEIGHT / (R_GAS * T);
            double rho_melt = magma.getDensity(P);
            double phi_gas = H2O_ex * rho_melt / (rho_gas + H2O_ex * (rho_melt - rho_gas));
            phi_gas = std::clamp(phi_gas, 0.0, 0.99);
            
            // Bulk density
            rho = (1.0 - phi_gas) * rho_melt + phi_gas * rho_gas;
            
            // Store state
            profile[i].z = z;
            profile[i].pressure = P;
            profile[i].temperature = T;
            profile[i].velocity = u;
            profile[i].density = rho;
            profile[i].gas_fraction = phi_gas;
            profile[i].dissolved_H2O = std::min(H2O_sat, magma.H2O_total);
            
            // Check fragmentation
            if (!fragmented && checkFragmentation(profile[i])) {
                fragmented = true;
                fragmentation_depth = z;
                profile[i].regime = ConduitFlowRegime::DISPERSED;
            } else if (fragmented) {
                profile[i].regime = ConduitFlowRegime::DISPERSED;
            } else if (phi_gas > 0.3) {
                profile[i].regime = ConduitFlowRegime::CHURN;
            } else {
                profile[i].regime = ConduitFlowRegime::BUBBLY;
            }
            
            // Momentum equation: ρu du/dz = -dP/dz - ρg - f_wall - τ_viscous
            double eta = computeViscosity(profile[i]);
            double f_wall = computeWallFriction(profile[i]);
            
            // Viscous stress gradient (simplified for pipe flow)
            // τ_viscous ≈ 8ηu/R² for laminar pipe flow
            double R_conduit = std::sqrt(geometry.areaAt(i * dz) / PI);
            double tau_viscous = 8.0 * eta * u / (R_conduit * R_conduit);
            
            // Pressure gradient (including all terms)
            double dP_dz = -rho * G - f_wall - tau_viscous - rho * u * (u / dz) * 0.1;  // Full balance
            
            // Update for next node
            if (i < n_nodes - 1) {
                P += dP_dz * dz;
                P = std::max(P, P_ATM);
                
                // Mass conservation: ρuA = const
                double A_next = geometry.areaAt((i + 1) * dz);
                double rho_next = magma.getDensity(P);
                u = u * rho * A / (rho_next * A_next);
            }
        }
        
        // Check convergence (vent pressure should be ~atmospheric)
        double P_vent = profile.back().pressure;
        if (std::abs(P_vent - P_ATM) < 0.1 * P_ATM) {
            // Converged
            break;
        }
        
        // Adjust initial velocity
        if (P_vent > P_ATM) {
            u_init *= 1.1;  // Need faster flow
        } else {
            u_init *= 0.9;
        }
    }
    
    // Compute mass flux
    double A_vent = geometry.areaAt(geometry.length);
    mass_flux = profile.back().density * profile.back().velocity * A_vent;
    
    return true;
}

void ConduitFlowModel::update(double dt) {
    // Time-dependent update (simplified)
    solvesteadyState();
}

double ConduitFlowModel::getExitVelocity() const {
    if (profile.empty()) return 0.0;
    return profile.back().velocity;
}

double ConduitFlowModel::getVolumeFlux() const {
    if (profile.empty()) return 0.0;
    double rho_dre = 2500.0;  // Dense rock equivalent
    return mass_flux / rho_dre;
}

double ConduitFlowModel::getColumnHeight() const {
    return VolcanoUtils::columnHeight_from_MER(mass_flux);
}

bool ConduitFlowModel::isChokedFlow() const {
    // Check if flow is sonic at vent
    if (profile.empty()) return false;
    
    const auto& vent = profile.back();
    double c_sound = std::sqrt(1.4 * R_GAS * vent.temperature / 
                               VolcanoConstants::H2O_MOLECULAR_WEIGHT);
    
    return vent.velocity >= 0.9 * c_sound;
}

double ConduitFlowModel::getMaxMassFlux() const {
    // Estimate choked flow limit
    double A_vent = geometry.areaAt(geometry.length);
    double c_sound = std::sqrt(1.4 * R_GAS * magma.temperature / 
                               VolcanoConstants::H2O_MOLECULAR_WEIGHT);
    double rho_gas = P_ATM * VolcanoConstants::H2O_MOLECULAR_WEIGHT / 
                     (R_GAS * magma.temperature);
    
    return rho_gas * c_sound * A_vent;
}

double ConduitFlowModel::computeViscosity(const ConduitFlowState& state) const {
    MagmaProperties local_magma = magma;
    local_magma.temperature = state.temperature;
    local_magma.crystal_fraction = state.crystal_fraction;
    local_magma.bubble_fraction = state.gas_fraction;
    
    return local_magma.getBulkViscosity();
}

double ConduitFlowModel::computeWallFriction(const ConduitFlowState& state) const {
    double eta = computeViscosity(state);
    double r = geometry.radiusAt(state.z);
    
    // Poiseuille friction factor
    double Re = state.density * state.velocity * 2 * r / eta;
    double f;
    
    if (Re < 2300) {
        f = 64.0 / std::max(Re, 1.0);
    } else {
        f = 0.316 / std::pow(Re, 0.25);  // Blasius formula
    }
    
    return f * state.density * state.velocity * state.velocity / (4 * r);
}

double ConduitFlowModel::computeBubbleGrowthRate(const ConduitFlowState& state) const {
    // Simplified bubble growth rate
    double P_sat = magma.getSaturationPressure();
    double dP = P_sat - state.pressure;
    
    if (dP <= 0) return 0.0;
    
    double eta = computeViscosity(state);
    double r = state.mean_bubble_radius;
    
    // Rayleigh-Plesset approximation
    return r * dP / (4 * eta);
}

bool ConduitFlowModel::checkFragmentation(const ConduitFlowState& state) const {
    // Strain rate criterion (Papale, 1999)
    double eta = computeViscosity(state);
    double du_dz = state.velocity / dz;  // Approximate strain rate
    double tau = eta * du_dz;
    
    // Tensile strength of magma
    double sigma_t = 1e6;  // ~1 MPa
    
    // Also check gas fraction criterion
    bool strain_criterion = tau > sigma_t;
    bool gas_criterion = state.gas_fraction > 0.75;
    
    return strain_criterion || gas_criterion;
}

double ConduitFlowModel::H2O_solubility(double P, double T) const {
    double P_MPa = P / 1e6;
    double T_C = T - 273.15;
    return VolcanoUtils::H2O_solubility(P_MPa, T_C, magma.SiO2);
}

double ConduitFlowModel::CO2_solubility(double P, double T) const {
    double P_MPa = P / 1e6;
    double T_C = T - 273.15;
    return VolcanoUtils::CO2_solubility(P_MPa, T_C, magma.SiO2);
}

double ConduitFlowModel::gasPhase_density(double P, double T, double X_H2O, double X_CO2) const {
    double MW = X_H2O * VolcanoConstants::H2O_MOLECULAR_WEIGHT +
                X_CO2 * VolcanoConstants::CO2_MOLECULAR_WEIGHT;
    return P * MW / (R_GAS * T);
}

// =============================================================================
// EruptionColumnModel Implementation
// =============================================================================

EruptionColumnModel::EruptionColumnModel()
    : max_height(0.0), NBL_height(0.0), column_collapsed(false), collapse_height(0.0) {}

void EruptionColumnModel::initialize(const EruptionColumnParameters& p) {
    params = p;
    profile.clear();
    max_height = 0.0;
    NBL_height = 0.0;
    column_collapsed = false;
    collapse_height = 0.0;
}

bool EruptionColumnModel::solve() {
    // Solve eruption column dynamics using integral model (Woods, 1988)
    profile.clear();
    
    // Initial conditions at vent
    ColumnState vent;
    vent.height = 0.0;
    vent.radius = params.vent_radius;
    vent.velocity = params.exit_velocity;
    vent.temperature = params.exit_temperature;
    vent.gas_fraction = params.gas_fraction;
    vent.particle_fraction = 1.0 - params.gas_fraction;
    vent.entrained_air = 0.0;
    
    // Gas and particle densities
    double rho_gas = P_ATM * 0.029 / (R_GAS * params.exit_temperature);  // Air MW
    double rho_particle = 2500.0;  // kg/m³
    vent.density = vent.gas_fraction * rho_gas + vent.particle_fraction * rho_particle;
    
    profile.push_back(vent);
    
    // Integration parameters
    double dz = 100.0;  // m step size
    double z_max = 50000.0;  // Maximum 50 km
    
    ColumnState state = vent;
    bool reached_NBL = false;
    
    while (state.height < z_max) {
        double z = state.height;
        
        // Atmospheric properties
        double rho_atm = atmosphericDensity(z);
        double T_atm = atmosphericTemperature(z);
        double P_atm_z = atmosphericPressure(z);
        
        // Entrainment
        double Ri = G * (rho_atm - state.density) * state.radius / 
                    (rho_atm * state.velocity * state.velocity);
        double alpha = entrainmentCoefficient(Ri);
        
        // Mass entrainment rate
        double dm_dz = 2 * PI * state.radius * alpha * rho_atm * state.velocity;
        
        // Current mass flux
        double m_dot = state.density * state.velocity * PI * state.radius * state.radius;
        
        // Momentum equation
        double buoyancy = (rho_atm - state.density) * G * PI * state.radius * state.radius;
        double d_momentum_dz = buoyancy;
        
        // Energy equation (enthalpy)
        double cp_mix = 1000.0;  // J/(kg·K) for mixture
        double cp_air = 1004.0;
        
        // Update state
        ColumnState new_state;
        new_state.height = z + dz;
        
        // Mass increases due to entrainment
        double new_mass_flux = m_dot + dm_dz * dz;
        
        // Entrained air fraction
        double entrained_mass = dm_dz * dz;
        new_state.entrained_air = (state.entrained_air * m_dot + entrained_mass) / new_mass_flux;
        
        // Momentum and velocity
        double momentum = m_dot * state.velocity;
        momentum += d_momentum_dz * dz;
        new_state.velocity = momentum / new_mass_flux;
        
        // Check for column collapse (velocity becomes negative or too small)
        if (new_state.velocity < 1.0) {
            column_collapsed = true;
            collapse_height = z;
            break;
        }
        
        // Temperature from energy balance
        double T_mix = (state.temperature * m_dot * cp_mix + T_atm * entrained_mass * cp_air) /
                       (new_mass_flux * cp_mix);
        new_state.temperature = T_mix;
        
        // Density of mixture
        double rho_gas_new = P_atm_z * 0.029 / (R_GAS * new_state.temperature);
        new_state.gas_fraction = state.gas_fraction * (1.0 - new_state.entrained_air) + 
                                 new_state.entrained_air;
        new_state.particle_fraction = state.particle_fraction * (1.0 - new_state.entrained_air);
        
        new_state.density = new_state.gas_fraction * rho_gas_new + 
                           new_state.particle_fraction * rho_particle;
        
        // Radius from mass conservation
        new_state.radius = std::sqrt(new_mass_flux / 
                                     (PI * new_state.density * new_state.velocity));
        
        // Check for neutral buoyancy
        if (!reached_NBL && new_state.density < rho_atm && state.density >= rho_atm) {
            NBL_height = z;
            reached_NBL = true;
        }
        
        // Check if column has become neutrally buoyant and spreading
        if (reached_NBL && new_state.velocity < 10.0) {
            max_height = z;
            break;
        }
        
        profile.push_back(new_state);
        state = new_state;
        max_height = z;
    }
    
    // If we didn't collapse, we reached maximum height
    return !column_collapsed;
}

double EruptionColumnModel::getUmbrellaCloudRadius() const {
    // Sparks (1986) umbrella cloud spreading
    if (NBL_height <= 0) return 0.0;
    
    // Umbrella spreads as gravity current
    double N = 0.01;  // Brunt-Väisälä frequency (1/s)
    double Q = params.mass_flux;  // kg/s
    
    // Time since eruption start (assume 1 hour)
    double t = 3600.0;
    
    // Radius grows as t^(2/3)
    double R = 0.2 * std::pow(Q / (N * N * N), 0.25) * std::pow(N * t, 2.0/3.0);
    
    return std::min(R, 500000.0);  // Cap at 500 km
}

double EruptionColumnModel::getCollapseFlux() const {
    if (!column_collapsed) return 0.0;
    
    // Fraction of mass flux that collapses
    return params.mass_flux * 0.8;  // Most of the flux collapses
}

double EruptionColumnModel::getSO2_emission_rate() const {
    // SO2 in volcanic gas
    return params.mass_flux * params.gas_fraction * params.X_SO2;
}

double EruptionColumnModel::getAsh_emission_rate() const {
    // Ash (particles) in eruption column
    return params.mass_flux * (1.0 - params.gas_fraction);
}

double EruptionColumnModel::atmosphericDensity(double z) const {
    // Standard atmosphere approximation
    double T = atmosphericTemperature(z);
    double P = atmosphericPressure(z);
    double M_air = 0.029;  // kg/mol
    return P * M_air / (R_GAS * T);
}

double EruptionColumnModel::atmosphericTemperature(double z) const {
    if (z < params.tropopause_height) {
        // Troposphere: -6.5 K/km lapse rate
        return params.atmospheric_temperature - 0.0065 * z;
    } else {
        // Stratosphere: approximately isothermal
        return params.stratosphere_temperature;
    }
}

double EruptionColumnModel::atmosphericPressure(double z) const {
    // Barometric formula
    double T0 = params.atmospheric_temperature;
    double L = 0.0065;  // Lapse rate K/m
    double M = 0.029;   // Molar mass of air
    
    if (z < params.tropopause_height) {
        return P_ATM * std::pow(1 - L * z / T0, G * M / (R_GAS * L));
    } else {
        // Above tropopause
        double P_trop = P_ATM * std::pow(1 - L * params.tropopause_height / T0, 
                                         G * M / (R_GAS * L));
        double dz = z - params.tropopause_height;
        return P_trop * std::exp(-G * M * dz / (R_GAS * params.stratosphere_temperature));
    }
}

double EruptionColumnModel::entrainmentCoefficient(double Ri) const {
    // Carazzo et al. (2008) Richardson-dependent entrainment
    double alpha_0 = 0.09;  // Jet limit
    double alpha_inf = 0.12;  // Plume limit
    
    if (Ri <= 0) {
        return alpha_0;  // Momentum-driven jet
    } else {
        // Buoyancy-driven plume
        return alpha_0 + (alpha_inf - alpha_0) * Ri / (Ri + 0.1);
    }
}

double EruptionColumnModel::particleSettlingVelocity(double diameter, double gas_density) const {
    // Settling velocity for spherical particles
    double rho_p = 2500.0;  // Particle density
    double mu_g = 1.8e-5;   // Gas viscosity (Pa·s)
    
    // Reynolds number iteration
    double Cd = 0.44;  // Drag coefficient (turbulent)
    double w = std::sqrt(4 * G * diameter * (rho_p - gas_density) / (3 * Cd * gas_density));
    
    double Re = gas_density * w * diameter / mu_g;
    
    // Adjust drag coefficient based on Re
    if (Re < 1) {
        Cd = 24.0 / Re;  // Stokes
    } else if (Re < 1000) {
        Cd = 24.0 / Re * (1 + 0.15 * std::pow(Re, 0.687));
    }
    
    w = std::sqrt(4 * G * diameter * (rho_p - gas_density) / (3 * Cd * gas_density));
    
    return w;
}

// =============================================================================
// PDCModel Implementation
// =============================================================================

PDCModel::PDCModel()
    : current_time(0.0), runout_distance(0.0), 
      nx(200), ny(200), dx(100.0), dy(100.0),
      xmin(-10000), xmax(10000), ymin(-10000), ymax(10000) {}

void PDCModel::initialize(const PDCSourceConditions& src) {
    source = src;
    current_time = 0.0;
    runout_distance = 0.0;
    
    // Initialize grids
    int n = nx * ny;
    deposit_thickness.assign(n, 0.0);
    max_dynamic_pressure.assign(n, 0.0);
    max_temperature.assign(n, 0.0);
    
    // Initialize current state with source conditions
    current_state.clear();
    
    PDCState initial;
    initial.x = 0.0;
    initial.y = 0.0;
    initial.z = source.source_height;
    
    // Initial velocity based on source type
    double v0 = source.initial_velocity;
    double azimuth_rad = source.azimuth * PI / 180.0;
    double incl_rad = source.inclination * PI / 180.0;
    
    initial.u = v0 * std::cos(incl_rad) * std::sin(azimuth_rad);
    initial.v = v0 * std::cos(incl_rad) * std::cos(azimuth_rad);
    initial.w = -v0 * std::sin(incl_rad);  // Downward
    
    initial.thickness = source.initial_radius;
    initial.temperature = source.initial_temperature;
    initial.particle_conc = source.particle_concentration;
    
    // Density of mixture
    double rho_gas = P_ATM * 0.029 / (R_GAS * source.initial_temperature);
    double rho_particle = 2500.0;
    initial.density = (1 - source.particle_concentration) * rho_gas +
                     source.particle_concentration * rho_particle;
    
    current_state.push_back(initial);
}

void PDCModel::setTopography(std::function<double(double, double)> elevation_func) {
    topography = elevation_func;
}

void PDCModel::update(double dt) {
    if (current_state.empty()) return;
    
    std::vector<PDCState> new_state;
    
    for (auto& state : current_state) {
        // Get topography
        double z_ground = topography ? topography(state.x, state.y) : 0.0;
        
        // Check if PDC has reached ground
        if (state.z <= z_ground + state.thickness) {
            state.z = z_ground + state.thickness;
            state.w = 0.0;  // No vertical velocity at ground
        }
        
        // Atmospheric density
        double rho_atm = P_ATM * 0.029 / (R_GAS * 288.0);  // Assume 288 K ambient
        
        // Drag force
        double drag = computeDragForce(state);
        
        // Gravity
        double g_eff = G * (state.density - rho_atm) / state.density;
        
        // Entrainment
        double entrainment = computeEntrainment(state);
        
        // Sedimentation
        double sedimentation = computeSedimentation(state);
        
        // Update velocities
        double speed = state.velocity_magnitude();
        if (speed > 0.1) {
            state.u -= drag * state.u / speed * dt;
            state.v -= drag * state.v / speed * dt;
        }
        state.w -= g_eff * dt;
        
        // Update position
        state.x += state.u * dt;
        state.y += state.v * dt;
        state.z += state.w * dt;
        state.z = std::max(state.z, z_ground);
        
        // Update thickness (entrainment minus sedimentation)
        state.thickness += (entrainment - sedimentation) * dt;
        state.thickness = std::max(state.thickness, 0.1);
        
        // Update particle concentration
        state.particle_conc -= sedimentation * state.particle_conc / state.thickness * dt;
        state.particle_conc = std::max(state.particle_conc, 0.001);
        
        // Update density
        double rho_gas = P_ATM * 0.029 / (R_GAS * state.temperature);
        state.density = (1 - state.particle_conc) * rho_gas + 
                       state.particle_conc * 2500.0;
        
        // Cooling
        state.temperature -= 0.1 * dt;  // Simple cooling rate
        state.temperature = std::max(state.temperature, 300.0);
        
        // Track maximum values on grid
        int i = static_cast<int>((state.x - xmin) / dx);
        int j = static_cast<int>((state.y - ymin) / dy);
        
        if (i >= 0 && i < nx && j >= 0 && j < ny) {
            int idx = i + j * nx;
            double dp = state.dynamic_pressure();
            max_dynamic_pressure[idx] = std::max(max_dynamic_pressure[idx], dp);
            max_temperature[idx] = std::max(max_temperature[idx], state.temperature);
            deposit_thickness[idx] += sedimentation * dt;
        }
        
        // Update runout distance
        double dist = std::sqrt(state.x * state.x + state.y * state.y);
        runout_distance = std::max(runout_distance, dist);
        
        // Keep state if still active
        if (state.velocity_magnitude() > 0.5 && state.particle_conc > 0.001) {
            new_state.push_back(state);
        }
    }
    
    current_state = new_state;
    current_time += dt;
    
    // Add new material from source if still active
    if (current_time < 600.0) {  // 10 minutes of source activity
        double source_rate = source.mass_flux * dt / (source.initial_radius * source.initial_radius * 2500.0);
        
        // Apply source term to particles in source region
        // Note: PDC is particle-based, so we would add new particles here
        // For now, just update thickness in existing particles near source
        for (auto& p : current_state) {
            double dx = p.x;
            double dy = p.y;
            double r = std::sqrt(dx*dx + dy*dy);
            
            if (r < source.initial_radius) {
                // Add thickness proportional to source rate and distance from center
                double weight = std::exp(-r*r / (source.initial_radius * source.initial_radius));
                p.thickness += source_rate * weight;
                p.temperature = std::max(p.temperature, source.initial_temperature);
            }
        }
    }
}

void PDCModel::runToCompletion(double max_time) {
    double dt = 0.5;  // Time step
    
    while (current_time < max_time && !current_state.empty()) {
        update(dt);
    }
}

double PDCModel::getRunoutDistance() const {
    return runout_distance;
}

double PDCModel::getMaxVelocity() const {
    double max_v = 0.0;
    for (const auto& state : current_state) {
        max_v = std::max(max_v, state.velocity_magnitude());
    }
    return max_v;
}

double PDCModel::getCurrentExtent() const {
    double max_dist = 0.0;
    for (const auto& state : current_state) {
        double dist = std::sqrt(state.x * state.x + state.y * state.y);
        max_dist = std::max(max_dist, dist);
    }
    return max_dist;
}

double PDCModel::getDynamicPressure(double x, double y) const {
    int i = static_cast<int>((x - xmin) / dx);
    int j = static_cast<int>((y - ymin) / dy);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return max_dynamic_pressure[i + j * nx];
    }
    return 0.0;
}

double PDCModel::getTemperature(double x, double y) const {
    int i = static_cast<int>((x - xmin) / dx);
    int j = static_cast<int>((y - ymin) / dy);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return max_temperature[i + j * nx];
    }
    return 288.0;  // Ambient
}

double PDCModel::getDepositThickness(double x, double y) const {
    int i = static_cast<int>((x - xmin) / dx);
    int j = static_cast<int>((y - ymin) / dy);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return deposit_thickness[i + j * nx];
    }
    return 0.0;
}

bool PDCModel::isInundated(double x, double y) const {
    return getDynamicPressure(x, y) > 0.0;
}

std::vector<std::pair<double, PDCState>> PDCModel::getTimeSeriesAt(double x, double y) const {
    // This would require storing time history - returning empty for now
    return {};
}

double PDCModel::computeDragForce(const PDCState& state) const {
    // Drag from ground and air entrainment
    double speed = state.velocity_magnitude();
    if (speed < 0.1) return 0.0;
    
    // Bottom friction
    double Cd = 0.01;  // Drag coefficient
    return Cd * speed * speed / state.thickness;
}

double PDCModel::computeEntrainment(const PDCState& state) const {
    // Air entrainment rate
    double speed = state.velocity_magnitude();
    double alpha = 0.1;  // Entrainment coefficient
    return alpha * speed;
}

double PDCModel::computeSedimentation(const PDCState& state) const {
    // Particle settling
    double ws = 1.0;  // Settling velocity (m/s) for mean particle size
    return ws * state.particle_conc;
}

double PDCModel::computeErosion(const PDCState& state) const {
    // Substrate erosion (if velocity high enough)
    double tau = state.dynamic_pressure();
    double tau_crit = 100.0;  // Critical shear stress
    
    if (tau > tau_crit) {
        return 0.001 * (tau - tau_crit);  // Erosion rate
    }
    return 0.0;
}

// =============================================================================
// LavaFlowModel Implementation
// =============================================================================

LavaFlowModel::LavaFlowModel()
    : nx(200), ny(200), dx(10.0), dy(10.0),
      xmin(-1000), xmax(1000), ymin(-1000), ymax(1000) {}

void LavaFlowModel::initialize(const LavaSourceParameters& src,
                               const MagmaProperties& mag) {
    source = src;
    magma = mag;
    
    // Set up grid around vent
    xmin = src.vent_x - 5000;
    xmax = src.vent_x + 5000;
    ymin = src.vent_y - 5000;
    ymax = src.vent_y + 5000;
    
    dx = (xmax - xmin) / nx;
    dy = (ymax - ymin) / ny;
    
    // Initialize state grid
    state.resize(nx * ny);
    for (int i = 0; i < nx * ny; ++i) {
        state[i].thickness = 0.0;
        state[i].velocity_x = 0.0;
        state[i].velocity_y = 0.0;
        state[i].temperature = 288.0;  // Ambient
        state[i].viscosity = 1e10;
        state[i].crust_thickness = 0.0;
        state[i].crystal_fraction = 0.0;
    }
    
    // Initialize vent cell
    int i_vent = static_cast<int>((src.vent_x - xmin) / dx);
    int j_vent = static_cast<int>((src.vent_y - ymin) / dy);
    
    if (i_vent >= 0 && i_vent < nx && j_vent >= 0 && j_vent < ny) {
        int idx = i_vent + j_vent * nx;
        state[idx].thickness = src.initial_thickness;
        state[idx].temperature = src.temperature;
        state[idx].crystal_fraction = magma.crystal_fraction;
        state[idx].viscosity = computeViscosity(src.temperature, magma.crystal_fraction);
    }
}

void LavaFlowModel::setTopography(std::function<double(double, double)> elevation_func) {
    topography = elevation_func;
}

void LavaFlowModel::setEffusionRate(double rate_m3_s) {
    source.effusion_rate = rate_m3_s;
}

void LavaFlowModel::update(double dt) {
    // Add lava at vent
    int i_vent = static_cast<int>((source.vent_x - xmin) / dx);
    int j_vent = static_cast<int>((source.vent_y - ymin) / dy);
    
    if (i_vent >= 0 && i_vent < nx && j_vent >= 0 && j_vent < ny) {
        int idx = i_vent + j_vent * nx;
        double dh = source.effusion_rate * dt / (dx * dy);
        state[idx].thickness += dh;
        state[idx].temperature = source.temperature;
    }
    
    // Compute flow for each cell
    std::vector<LavaFlowState> new_state = state;
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = i + j * nx;
            
            if (state[idx].thickness < 0.01) continue;  // No lava
            if (!state[idx].isMobile()) continue;  // Solidified
            
            // Get topography gradients
            double x = xmin + i * dx;
            double y = ymin + j * dy;
            double z0 = topography ? topography(x, y) : 0.0;
            double z_xp = topography ? topography(x + dx, y) : 0.0;
            double z_xm = topography ? topography(x - dx, y) : 0.0;
            double z_yp = topography ? topography(x, y + dy) : 0.0;
            double z_ym = topography ? topography(x, y - dy) : 0.0;
            
            // Surface slope (topography + lava surface)
            double h0 = state[idx].thickness;
            double total_z0 = z0 + h0;
            
            double dz_dx = (z_xp - z_xm) / (2 * dx);
            double dz_dy = (z_yp - z_ym) / (2 * dy);
            
            // Also include lava surface slope
            int idx_xp = (i + 1) + j * nx;
            int idx_xm = (i - 1) + j * nx;
            int idx_yp = i + (j + 1) * nx;
            int idx_ym = i + (j - 1) * nx;
            
            // Get neighboring thicknesses
            double h_xp = state[idx_xp].thickness;
            double h_xm = state[idx_xm].thickness;
            double h_yp = state[idx_yp].thickness;
            double h_ym = state[idx_ym].thickness;
            
            double dh_dx = (h_xp - h_xm) / (2 * dx);
            double dh_dy = (h_yp - h_ym) / (2 * dy);
            
            // Total surface gradient
            double Sx = -(dz_dx + dh_dx);
            double Sy = -(dz_dy + dh_dy);
            
            // Viscosity (temperature and crystal dependent)
            double eta = computeViscosity(state[idx].temperature, state[idx].crystal_fraction);
            double tau_y = computeYieldStrength(state[idx].crystal_fraction);
            
            // Bingham flow velocity
            double S_mag = std::sqrt(Sx * Sx + Sy * Sy);
            double rho = 2500.0;  // Lava density
            
            // Critical height for yield stress
            double h_crit = tau_y / (rho * G * std::max(S_mag, 1e-6));
            
            double u_x = 0.0, u_y = 0.0;
            if (h0 > h_crit) {
                // Bingham flow velocity profile average
                double h_eff = h0 - h_crit;
                double K = rho * G * S_mag * h_eff * h_eff / (3 * eta);
                u_x = K * Sx / std::max(S_mag, 1e-6);
                u_y = K * Sy / std::max(S_mag, 1e-6);
            }
            
            new_state[idx].velocity_x = u_x;
            new_state[idx].velocity_y = u_y;
            
            // Mass flux out of cell
            double flux_x = h0 * u_x;
            double flux_y = h0 * u_y;
            
            // Update thicknesses
            double dh_dt = 0.0;
            
            // Outflow from this cell
            if (u_x > 0 && i < nx - 1) {
                dh_dt -= flux_x / dx;
                new_state[idx_xp].thickness += flux_x * dt / dx;
            } else if (u_x < 0 && i > 0) {
                dh_dt += flux_x / dx;
                new_state[idx_xm].thickness -= flux_x * dt / dx;
            }
            
            if (u_y > 0 && j < ny - 1) {
                dh_dt -= flux_y / dy;
                new_state[idx_yp].thickness += flux_y * dt / dy;
            } else if (u_y < 0 && j > 0) {
                dh_dt += flux_y / dy;
                new_state[idx_ym].thickness -= flux_y * dt / dy;
            }
            
            new_state[idx].thickness += dh_dt * dt;
            new_state[idx].thickness = std::max(new_state[idx].thickness, 0.0);
            
            // Cooling
            double cooling = computeCoolingRate(state[idx]);
            new_state[idx].temperature -= cooling * dt;
            new_state[idx].temperature = std::max(new_state[idx].temperature, 288.0);
            
            // Crystallization from cooling
            if (new_state[idx].temperature < magma.getLiquidusTemperature()) {
                double dT = cooling * dt;
                double dphi = magma.getCrystallizationRate() * dT;
                new_state[idx].crystal_fraction += dphi;
                new_state[idx].crystal_fraction = std::min(new_state[idx].crystal_fraction, 0.7);
            }
            
            // Crust growth
            double crust_growth = computeCrustGrowth(state[idx], dt);
            new_state[idx].crust_thickness += crust_growth;
            
            // Update viscosity
            new_state[idx].viscosity = computeViscosity(new_state[idx].temperature, 
                                                        new_state[idx].crystal_fraction);
        }
    }
    
    state = new_state;
}

void LavaFlowModel::run(double duration) {
    double dt = 1.0;  // 1 second time step
    double t = 0.0;
    
    while (t < duration) {
        update(dt);
        t += dt;
    }
}

double LavaFlowModel::getFlowArea() const {
    double area = 0.0;
    for (const auto& s : state) {
        if (s.thickness > 0.01) {
            area += dx * dy;
        }
    }
    return area;
}

double LavaFlowModel::getFlowLength() const {
    double max_dist = 0.0;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (state[i + j * nx].thickness > 0.01) {
                double x = xmin + i * dx - source.vent_x;
                double y = ymin + j * dy - source.vent_y;
                double dist = std::sqrt(x*x + y*y);
                max_dist = std::max(max_dist, dist);
            }
        }
    }
    return max_dist;
}

double LavaFlowModel::getFlowVolume() const {
    double volume = 0.0;
    for (const auto& s : state) {
        volume += s.thickness * dx * dy;
    }
    return volume;
}

const LavaFlowState* LavaFlowModel::getStateAt(double x, double y) const {
    int i = static_cast<int>((x - xmin) / dx);
    int j = static_cast<int>((y - ymin) / dy);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return &state[i + j * nx];
    }
    return nullptr;
}

double LavaFlowModel::getThicknessAt(double x, double y) const {
    const LavaFlowState* s = getStateAt(x, y);
    return s ? s->thickness : 0.0;
}

bool LavaFlowModel::isInundated(double x, double y) const {
    return getThicknessAt(x, y) > 0.01;
}

double LavaFlowModel::getRadiativeHeatFlux() const {
    double Q = 0.0;
    double epsilon = 0.9;  // Emissivity
    
    for (const auto& s : state) {
        if (s.thickness > 0.01) {
            // Stefan-Boltzmann radiation
            double T_surface = s.temperature - 200.0;  // Crust is cooler
            T_surface = std::max(T_surface, 288.0);
            double q = epsilon * STEFAN_BOLTZMANN * std::pow(T_surface, 4);
            Q += q * dx * dy;
        }
    }
    return Q;
}

double LavaFlowModel::computeViscosity(double T, double phi_crystal) const {
    MagmaProperties local = magma;
    local.temperature = T;
    local.crystal_fraction = phi_crystal;
    return local.getBulkViscosity();
}

double LavaFlowModel::computeYieldStrength(double phi_crystal) const {
    if (phi_crystal < 0.25) return 0.0;
    
    // Yield strength increases dramatically above 25% crystals
    double phi_excess = phi_crystal - 0.25;
    return 1e4 * std::pow(phi_excess / 0.35, 2.0);  // Pa
}

double LavaFlowModel::computeCoolingRate(const LavaFlowState& s) const {
    if (s.thickness < 0.01) return 0.0;
    
    double T_ambient = 288.0;
    double T_surface = s.temperature;
    
    // Radiative cooling
    double epsilon = 0.9;
    double q_rad = epsilon * STEFAN_BOLTZMANN * (std::pow(T_surface, 4) - std::pow(T_ambient, 4));
    
    // Convective cooling
    double h_conv = 10.0;  // W/(m²·K)
    double q_conv = h_conv * (T_surface - T_ambient);
    
    // Total heat loss
    double q_total = q_rad + q_conv;
    
    // Temperature decrease rate
    double rho = 2500.0;
    double cp = VolcanoConstants::MAGMA_SPECIFIC_HEAT;
    
    return q_total / (rho * cp * s.thickness);
}

double LavaFlowModel::computeCrustGrowth(const LavaFlowState& s, double dt) const {
    if (s.thickness < 0.01) return 0.0;
    
    // Crust grows by conduction
    double kappa = 1e-6;  // Thermal diffusivity (m²/s)
    double T_interior = s.temperature;
    double T_solidus = magma.getSolidusTemperature();
    
    if (T_interior < T_solidus) {
        // Already solid
        return s.thickness - s.crust_thickness;
    }
    
    // Crust grows as sqrt(kappa * t)
    double growth_rate = 0.5 * std::sqrt(kappa / std::max(dt, 1.0));
    return growth_rate * dt;
}

void LavaFlowModel::computeFluxes() {
    // Helper function for flux computation (used in update)
}

void LavaFlowModel::advectTemperature(double dt) {
    // Temperature advection with flow (used in update)
}

// =============================================================================
// LaharModel Implementation
// =============================================================================

LaharModel::LaharModel()
    : nx(200), ny(200), dx(50.0), dy(50.0), current_time(0.0) {}

void LaharModel::initialize(const LaharSourceConditions& src) {
    source = src;
    current_time = 0.0;
    
    // Set up grid centered on source
    double domain_size = 50000.0;  // 50 km domain
    nx = 200;
    ny = 200;
    dx = domain_size / nx;
    dy = domain_size / ny;
    
    // Initialize state grid
    state.resize(nx * ny);
    max_depth.resize(nx * ny, 0.0);
    max_velocity.resize(nx * ny, 0.0);
    arrival_time.resize(nx * ny, -1.0);
    
    for (auto& s : state) {
        s.depth = 0.0;
        s.velocity_x = 0.0;
        s.velocity_y = 0.0;
        s.sediment_conc = 0.0;
        s.temperature = 288.0;
    }
}

void LaharModel::setTopography(std::function<double(double, double)> elevation_func) {
    topography = elevation_func;
}

void LaharModel::setChannel(std::function<double(double, double)> width_func) {
    channel_width = width_func;
}

void LaharModel::update(double dt) {
    // Source term - add material at source location
    double source_rate = 0.0;
    if (current_time < source.duration) {
        // Hydrograph shape
        double t_peak = source.duration / 3.0;
        if (current_time < t_peak) {
            source_rate = source.peak_discharge * (current_time / t_peak);
        } else {
            source_rate = source.peak_discharge * (1.0 - (current_time - t_peak) / (source.duration - t_peak));
        }
    }
    
    // Find source cell
    double x0 = nx / 2 * dx;  // Source position x
    double y0 = ny / 2 * dy;  // Source position y
    int i_src = nx / 2;
    int j_src = ny / 2;
    
    if (source_rate > 0) {
        int idx = i_src + j_src * nx;
        double A_cell = dx * dy;
        
        // Distribute source over Gaussian region for smoother input
        double sigma = 2.0 * std::max(dx, dy);
        for (int j = std::max(0, j_src - 3); j < std::min(ny, j_src + 4); ++j) {
            for (int i = std::max(0, i_src - 3); i < std::min(nx, i_src + 4); ++i) {
                int idx_ij = i + j * nx;
                double x_ij = i * dx;
                double y_ij = j * dy;
                double r2 = (x_ij - x0) * (x_ij - x0) + (y_ij - y0) * (y_ij - y0);
                double weight = std::exp(-r2 / (2 * sigma * sigma));
                
                state[idx_ij].depth += source_rate * dt / A_cell * weight;
                state[idx_ij].sediment_conc = source.sediment_concentration;
                state[idx_ij].temperature = source.temperature;
            }
        }
    }
    
    // Shallow water solver
    std::vector<LaharState> new_state = state;
    
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            int idx = i + j * nx;
            
            if (state[idx].depth < 0.001) continue;
            
            // Get topography
            double x = i * dx;
            double y = j * dy;
            double z0 = topography ? topography(x, y) : source.source_elevation - 0.01 * std::sqrt(x*x + y*y);
            double z_xp = topography ? topography(x + dx, y) : z0 - 0.01 * dx;
            double z_xm = topography ? topography(x - dx, y) : z0 + 0.01 * dx;
            double z_yp = topography ? topography(x, y + dy) : z0 - 0.01 * dy;
            double z_ym = topography ? topography(x, y - dy) : z0 + 0.01 * dy;
            
            // Surface slopes
            double dz_dx = (z_xp - z_xm) / (2 * dx);
            double dz_dy = (z_yp - z_ym) / (2 * dy);
            
            // Gravity acceleration
            double g_x = -G * dz_dx;
            double g_y = -G * dz_dy;
            
            // Friction
            double friction = computeFriction(state[idx]);
            double speed = state[idx].velocity_magnitude();
            
            double f_x = 0.0, f_y = 0.0;
            if (speed > 0.01) {
                f_x = friction * state[idx].velocity_x / speed;
                f_y = friction * state[idx].velocity_y / speed;
            }
            
            // Update velocities
            new_state[idx].velocity_x = state[idx].velocity_x + (g_x - f_x) * dt;
            new_state[idx].velocity_y = state[idx].velocity_y + (g_y - f_y) * dt;
            
            // Limit velocity
            double max_v = 50.0;  // m/s
            double new_speed = new_state[idx].velocity_magnitude();
            if (new_speed > max_v) {
                new_state[idx].velocity_x *= max_v / new_speed;
                new_state[idx].velocity_y *= max_v / new_speed;
            }
            
            // Mass conservation - flux to neighbors
            double h = state[idx].depth;
            double u = new_state[idx].velocity_x;
            double v = new_state[idx].velocity_y;
            
            int idx_xp = (i + 1) + j * nx;
            int idx_xm = (i - 1) + j * nx;
            int idx_yp = i + (j + 1) * nx;
            int idx_ym = i + (j - 1) * nx;
            
            // Outflow
            double dh_dt = 0.0;
            if (u > 0) {
                double flux = h * u / dx;
                dh_dt -= flux;
                new_state[idx_xp].depth += flux * dt;
            } else if (u < 0) {
                double flux = -h * u / dx;
                dh_dt -= flux;
                new_state[idx_xm].depth += flux * dt;
            }
            
            if (v > 0) {
                double flux = h * v / dy;
                dh_dt -= flux;
                new_state[idx_yp].depth += flux * dt;
            } else if (v < 0) {
                double flux = -h * v / dy;
                dh_dt -= flux;
                new_state[idx_ym].depth += flux * dt;
            }
            
            new_state[idx].depth += dh_dt * dt;
            new_state[idx].depth = std::max(new_state[idx].depth, 0.0);
            
            // Erosion and deposition
            double erosion = computeErosion(state[idx]);
            double deposition = computeDeposition(state[idx]);
            
            new_state[idx].sediment_conc += (erosion - deposition) * dt;
            new_state[idx].sediment_conc = std::clamp(new_state[idx].sediment_conc, 0.0, 0.7);
            
            // Track maxima
            double v_mag = new_state[idx].velocity_magnitude();
            max_depth[idx] = std::max(max_depth[idx], new_state[idx].depth);
            max_velocity[idx] = std::max(max_velocity[idx], v_mag);
            
            if (arrival_time[idx] < 0 && new_state[idx].depth > 0.1) {
                arrival_time[idx] = current_time;
            }
        }
    }
    
    state = new_state;
    current_time += dt;
}

void LaharModel::runToCompletion(double max_time) {
    double dt = 0.5;  // Time step
    
    while (current_time < max_time) {
        update(dt);
        
        // Check if lahar has stopped (all velocities near zero)
        double max_v = 0.0;
        for (const auto& s : state) {
            max_v = std::max(max_v, s.velocity_magnitude());
        }
        if (max_v < 0.1) break;
    }
}

double LaharModel::getRunoutDistance() const {
    double max_dist = 0.0;
    double x0 = nx / 2 * dx;
    double y0 = ny / 2 * dy;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            if (max_depth[i + j * nx] > 0.1) {
                double x = i * dx - x0;
                double y = j * dy - y0;
                double dist = std::sqrt(x*x + y*y);
                max_dist = std::max(max_dist, dist);
            }
        }
    }
    return max_dist;
}

double LaharModel::getCurrentVolume() const {
    double volume = 0.0;
    for (const auto& s : state) {
        volume += s.depth * dx * dy;
    }
    return volume;
}

double LaharModel::getPeakDischarge() const {
    return source.peak_discharge;
}

double LaharModel::getDepthAt(double x, double y) const {
    int i = static_cast<int>(x / dx + nx / 2);
    int j = static_cast<int>(y / dy + ny / 2);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return state[i + j * nx].depth;
    }
    return 0.0;
}

double LaharModel::getVelocityAt(double x, double y) const {
    int i = static_cast<int>(x / dx + nx / 2);
    int j = static_cast<int>(y / dy + ny / 2);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return state[i + j * nx].velocity_magnitude();
    }
    return 0.0;
}

bool LaharModel::isInundated(double x, double y) const {
    return getMaxDepth(x, y) > 0.1;
}

double LaharModel::getTravelTimeTo(double x, double y) const {
    int i = static_cast<int>(x / dx + nx / 2);
    int j = static_cast<int>(y / dy + ny / 2);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return arrival_time[i + j * nx];
    }
    return -1.0;
}

double LaharModel::getMaxDepth(double x, double y) const {
    int i = static_cast<int>(x / dx + nx / 2);
    int j = static_cast<int>(y / dy + ny / 2);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return max_depth[i + j * nx];
    }
    return 0.0;
}

double LaharModel::getMaxVelocity(double x, double y) const {
    int i = static_cast<int>(x / dx + nx / 2);
    int j = static_cast<int>(y / dy + ny / 2);
    
    if (i >= 0 && i < nx && j >= 0 && j < ny) {
        return max_velocity[i + j * nx];
    }
    return 0.0;
}

double LaharModel::getMaxDynamicPressure(double x, double y) const {
    double h = getMaxDepth(x, y);
    double v = getMaxVelocity(x, y);
    double rho = 1800.0;  // Approximate lahar density
    double g = 9.81;
    
    // Total pressure = dynamic + hydrostatic components
    double P_dynamic = 0.5 * rho * v * v;
    double P_hydrostatic = rho * g * h;
    
    // Return total impact pressure on structures
    return P_dynamic + P_hydrostatic;
}

double LaharModel::computeFriction(const LaharState& s) const {
    if (s.depth < 0.01) return 0.0;
    
    // O'Brien friction model for debris flows
    double tau_y = 100.0;  // Yield stress (Pa)
    double K = 0.05;  // Viscous parameter
    double n_t = 0.1;  // Turbulent coefficient
    
    double rho = s.density();
    double v = s.velocity_magnitude();
    
    // Friction slope
    double Sf = tau_y / (rho * G * s.depth) + K * v / s.depth + n_t * n_t * v * v / std::pow(s.depth, 4.0/3.0);
    
    return G * Sf;
}

double LaharModel::computeErosion(const LaharState& s) const {
    if (s.depth < 0.01) return 0.0;
    
    double v = s.velocity_magnitude();
    double v_crit = 2.0;  // Critical velocity for erosion
    
    if (v > v_crit) {
        return 0.001 * (v - v_crit);  // Erosion rate (1/s)
    }
    return 0.0;
}

double LaharModel::computeDeposition(const LaharState& s) const {
    if (s.sediment_conc < 0.01) return 0.0;
    
    double v = s.velocity_magnitude();
    double v_crit = 1.0;  // Critical velocity for deposition
    
    if (v < v_crit) {
        return 0.01 * s.sediment_conc * (v_crit - v);  // Deposition rate
    }
    return 0.0;
}

// =============================================================================
// VolcanicDeformationModel Implementation
// =============================================================================

VolcanicDeformationModel::VolcanicDeformationModel()
    : shear_modulus(VolcanoConstants::SHEAR_MODULUS_CRUST), poisson_ratio(0.25) {}

void VolcanicDeformationModel::addSource(const VolcanicDeformationSource& source) {
    sources.push_back(source);
}

void VolcanicDeformationModel::clearSources() {
    sources.clear();
}

void VolcanicDeformationModel::setElasticProperties(double mu, double nu) {
    shear_modulus = mu;
    poisson_ratio = nu;
}

void VolcanicDeformationModel::computeDisplacement(double x, double y,
                                                    double& ux, double& uy, double& uz) const {
    ux = uy = uz = 0.0;
    
    for (const auto& src : sources) {
        double dux, duy, duz;
        
        switch (src.type) {
            case DeformationSourceType::MOGI:
                mogiDisplacement(x - src.x, y - src.y, src.depth, src.volume_change,
                                dux, duy, duz);
                break;
                
            case DeformationSourceType::MCTIGUE:
                mctigueDisplacement(x - src.x, y - src.y, src.depth, src.radius, 
                                   src.pressure_change, dux, duy, duz);
                break;
                
            case DeformationSourceType::RECTANGULAR_DIKE:
                okataDikeDisplacement(x, y, src, dux, duy, duz);
                break;
                
            case DeformationSourceType::HORIZONTAL_SILL:
                okataSillDisplacement(x, y, src, dux, duy, duz);
                break;
                
            case DeformationSourceType::PROLATE_SPHEROID:
            case DeformationSourceType::OBLATE_SPHEROID:
                spheroidDisplacement(x, y, src, dux, duy, duz);
                break;
                
            case DeformationSourceType::PENNY_CRACK:
                // Use sill approximation for penny crack
                okataSillDisplacement(x, y, src, dux, duy, duz);
                break;
                
            default:
                // Default to Mogi
                mogiDisplacement(x - src.x, y - src.y, src.depth, src.volume_change,
                                dux, duy, duz);
                break;
        }
        
        ux += dux;
        uy += duy;
        uz += duz;
    }
}

void VolcanicDeformationModel::computeTilt(double x, double y,
                                           double& tilt_x, double& tilt_y) const {
    // Compute tilt from displacement gradient
    double h = 10.0;  // Finite difference step
    
    double uz_xp, uz_xm, uz_yp, uz_ym;
    double ux_tmp, uy_tmp;
    
    computeDisplacement(x + h, y, ux_tmp, uy_tmp, uz_xp);
    computeDisplacement(x - h, y, ux_tmp, uy_tmp, uz_xm);
    computeDisplacement(x, y + h, ux_tmp, uy_tmp, uz_yp);
    computeDisplacement(x, y - h, ux_tmp, uy_tmp, uz_ym);
    
    tilt_x = (uz_xp - uz_xm) / (2 * h);
    tilt_y = (uz_yp - uz_ym) / (2 * h);
}

void VolcanicDeformationModel::computeStrain(double x, double y,
                                              double& exx, double& eyy, double& exy) const {
    // Compute strain from displacement gradient
    double h = 10.0;
    
    double ux_xp, uy_xp, uz_xp;
    double ux_xm, uy_xm, uz_xm;
    double ux_yp, uy_yp, uz_yp;
    double ux_ym, uy_ym, uz_ym;
    
    computeDisplacement(x + h, y, ux_xp, uy_xp, uz_xp);
    computeDisplacement(x - h, y, ux_xm, uy_xm, uz_xm);
    computeDisplacement(x, y + h, ux_yp, uy_yp, uz_yp);
    computeDisplacement(x, y - h, ux_ym, uy_ym, uz_ym);
    
    double dux_dx = (ux_xp - ux_xm) / (2 * h);
    double duy_dy = (uy_yp - uy_ym) / (2 * h);
    double dux_dy = (ux_yp - ux_ym) / (2 * h);
    double duy_dx = (uy_xp - uy_xm) / (2 * h);
    
    exx = dux_dx;
    eyy = duy_dy;
    exy = 0.5 * (dux_dy + duy_dx);
}

double VolcanicDeformationModel::computeGravityChange(double x, double y) const {
    // Gravity change from deformation (Bouguer approximation)
    // dg = -2*pi*G*rho*dz + free-air effect
    
    double ux, uy, uz;
    computeDisplacement(x, y, ux, uy, uz);
    
    double rho = VolcanoConstants::CRUSTAL_DENSITY;
    double G_grav = 6.674e-11;  // Gravitational constant
    
    // Free-air gradient: -0.3086 mGal/m
    double free_air = -0.3086e-5 * uz;  // m/s²
    
    // Bouguer correction: 2*pi*G*rho = 0.04193 mGal/m for 2670 kg/m³
    double bouguer = 2 * PI * G_grav * rho * uz;
    
    return free_air + bouguer;
}

void VolcanicDeformationModel::computeDisplacementField(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_coords,
    std::vector<double>& ux,
    std::vector<double>& uy,
    std::vector<double>& uz) const {
    
    size_t n = x_coords.size();
    ux.resize(n);
    uy.resize(n);
    uz.resize(n);
    
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        computeDisplacement(x_coords[i], y_coords[i], ux[i], uy[i], uz[i]);
    }
}

double VolcanicDeformationModel::computeLOS_displacement(
    double x, double y,
    double look_azimuth, double incidence_angle) const {
    
    double ux, uy, uz;
    computeDisplacement(x, y, ux, uy, uz);
    
    // Convert angles to radians
    double az_rad = look_azimuth * PI / 180.0;
    double inc_rad = incidence_angle * PI / 180.0;
    
    // Unit vector toward satellite
    double e_x = std::sin(inc_rad) * std::sin(az_rad);
    double e_y = std::sin(inc_rad) * std::cos(az_rad);
    double e_z = std::cos(inc_rad);
    
    // LOS displacement (positive = toward satellite)
    return -(ux * e_x + uy * e_y + uz * e_z);
}

void VolcanicDeformationModel::mogiDisplacement(double x, double y, double depth, double dV,
                                                 double& ux, double& uy, double& uz) const {
    // Mogi (1958) point source
    double r = std::sqrt(x*x + y*y);
    double R = std::sqrt(r*r + depth*depth);
    
    if (R < 1.0) R = 1.0;  // Avoid singularity
    
    double C = (1 - poisson_ratio) * dV / PI;
    
    double R3 = R * R * R;
    
    // Radial and vertical displacement
    double ur = C * r / R3;
    uz = C * depth / R3;
    
    // Convert radial to x, y components
    if (r > 1.0) {
        ux = ur * x / r;
        uy = ur * y / r;
    } else {
        ux = 0.0;
        uy = 0.0;
    }
}

void VolcanicDeformationModel::mctigueDisplacement(double x, double y, double depth, 
                                                    double radius, double dP,
                                                    double& ux, double& uy, double& uz) const {
    // McTigue (1987) finite spherical cavity
    double r = std::sqrt(x*x + y*y);
    double R = std::sqrt(r*r + depth*depth);
    
    if (R < radius) R = radius;  // Inside source
    
    double a = radius;
    double a_d = a / depth;  // Dimensionless source radius
    
    // Volume change equivalent
    double dV = PI * a * a * a * dP * (1 - poisson_ratio) / shear_modulus;
    
    // First-order term (Mogi point source)
    double C = (1 - poisson_ratio) * dV / PI;
    double R3 = R * R * R;
    
    double ur_mogi = C * r / R3;
    double uz_mogi = C * depth / R3;
    
    // McTigue (1987) finite sphere correction terms
    // Includes O(a/d)^2 terms for finite source effects
    double nu = poisson_ratio;
    double a_d2 = a_d * a_d;
    double a_d3 = a_d2 * a_d;
    
    // Radial correction factor
    double r_R = r / R;
    double d_R = depth / R;
    double r_R2 = r_R * r_R;
    
    // Second-order correction (McTigue eq. 5-6)
    double corr_r = 1.0 + a_d2 * (1.0 + nu) / 2.0 
                       + a_d3 * (3.0 * r_R2 - 1.0) * (1.0 + nu) / 4.0;
    double corr_z = 1.0 + a_d2 * (1.0 + nu) / 2.0 
                       - a_d3 * (3.0 * d_R * d_R - 1.0) * (1.0 + nu) / 4.0;
    
    double ur = ur_mogi * corr_r;
    uz = uz_mogi * corr_z;
    
    // Convert radial to x, y
    if (r > 1.0) {
        ux = ur * x / r;
        uy = ur * y / r;
    } else {
        ux = 0.0;
        uy = 0.0;
    }
}

void VolcanicDeformationModel::okataDikeDisplacement(double x, double y, 
                                                      const VolcanicDeformationSource& src,
                                                      double& ux, double& uy, double& uz) const {
    // Simplified Okada (1985) for vertical dike
    // Full implementation would require complete Okada equations
    
    double strike_rad = src.strike * PI / 180.0;
    double dip_rad = src.dip * PI / 180.0;
    
    // Rotate coordinates to fault-aligned system
    double cos_s = std::cos(strike_rad);
    double sin_s = std::sin(strike_rad);
    
    double x_rot = (x - src.x) * sin_s + (y - src.y) * cos_s;
    double y_rot = -(x - src.x) * cos_s + (y - src.y) * sin_s;
    
    // Distance from dike
    double d = src.depth + y_rot * std::sin(dip_rad);
    double R = std::sqrt(x_rot * x_rot + y_rot * y_rot + d * d);
    
    if (R < 1.0) R = 1.0;
    
    // Tensile opening
    double U = src.opening;
    double L = src.length;
    double W = src.width;
    
    // Approximate displacement (simplified)
    double factor = U * L * W / (4 * PI * R * R * R) * (1 - poisson_ratio);
    
    double uz_local = factor * d;
    double ux_local = factor * x_rot;
    double uy_local = factor * y_rot;
    
    // Rotate back
    ux = ux_local * sin_s - uy_local * cos_s;
    uy = ux_local * cos_s + uy_local * sin_s;
    uz = uz_local;
}

void VolcanicDeformationModel::okataSillDisplacement(double x, double y,
                                                      const VolcanicDeformationSource& src,
                                                      double& ux, double& uy, double& uz) const {
    // Horizontal sill (simplified)
    double dx = x - src.x;
    double dy = y - src.y;
    double r = std::sqrt(dx*dx + dy*dy);
    double d = src.depth;
    double R = std::sqrt(r*r + d*d);
    
    if (R < 1.0) R = 1.0;
    
    double U = src.opening;
    double L = src.length;
    double W = src.width;
    double A = L * W;
    
    // Displacement from opening sill
    double factor = U * A / (2 * PI * R * R * R) * (1 - poisson_ratio);
    
    uz = factor * d;
    double ur = factor * r;
    
    if (r > 1.0) {
        ux = ur * dx / r;
        uy = ur * dy / r;
    } else {
        ux = 0.0;
        uy = 0.0;
    }
}

void VolcanicDeformationModel::spheroidDisplacement(double x, double y,
                                                     const VolcanicDeformationSource& src,
                                                     double& ux, double& uy, double& uz) const {
    // Yang et al. (1988) spheroid source (simplified)
    // Uses Mogi as approximation with aspect ratio correction
    
    double a = src.semi_major;  // Horizontal
    double c = src.semi_minor;  // Vertical
    
    // Effective volume
    double V = (4.0/3.0) * PI * a * a * c;
    
    // Aspect ratio correction
    double aspect = c / a;
    double corr = 1.0;
    if (src.type == DeformationSourceType::PROLATE_SPHEROID) {
        corr = 1.0 + 0.3 * (aspect - 1.0);
    } else {
        corr = 1.0 - 0.3 * (1.0 - aspect);
    }
    
    // Use Mogi with corrected volume
    double dV = src.pressure_change * V * (1 - poisson_ratio) / shear_modulus * corr;
    
    mogiDisplacement(x - src.x, y - src.y, src.depth, dV, ux, uy, uz);
}

// =============================================================================
// VolcanicSeismicityModel Implementation
// =============================================================================

VolcanicSeismicityModel::VolcanicSeismicityModel() {
    // Default crack model parameters
    crack_model.length = 100.0;
    crack_model.width = 50.0;
    crack_model.stiffness = 1e9;
    crack_model.fluid_density = 2500.0;
    crack_model.fluid_viscosity = 100.0;
    crack_model.fluid_velocity = 1000.0;
}

void VolcanicSeismicityModel::generateVT_source(const VolcanicSeismicEvent& event,
                                                 std::vector<double>& times,
                                                 std::vector<double>& moment_rate) {
    // VT earthquake - standard double-couple source
    double M0 = event.seismic_moment;
    double rise_time = 0.1 + 0.1 * std::pow(M0 / 1e12, 1.0/3.0);  // Empirical scaling
    
    // Time vector
    int n_samples = 500;
    double dt = 0.01;
    times.resize(n_samples);
    moment_rate.resize(n_samples);
    
    for (int i = 0; i < n_samples; ++i) {
        times[i] = i * dt;
        double t = times[i];
        
        // Regularized source time function (triangle or smoothed boxcar)
        double value, deriv;
        applySourceTimeFunction(t, rise_time, value, deriv);
        moment_rate[i] = M0 * deriv;
    }
}

void VolcanicSeismicityModel::generateLP_source(const VolcanicSeismicEvent& event,
                                                 std::vector<double>& times,
                                                 std::vector<std::array<double, 6>>& moment_tensor_rate) {
    // LP event - resonating fluid-filled crack (Chouet, 1986)
    double f0 = event.dominant_frequency;
    double Q = event.quality_factor;
    double P0 = event.pressure_transient;
    
    // Crack dimensions
    double L = event.crack_length;
    double W = event.crack_width;
    double A = L * W;
    
    // Time vector
    int n_samples = 1000;
    double dt = 0.01;
    times.resize(n_samples);
    moment_tensor_rate.resize(n_samples);
    
    for (int i = 0; i < n_samples; ++i) {
        times[i] = i * dt;
        double t = times[i];
        
        // Decaying oscillation
        double amp = decayingOscillation(t, f0, Q);
        
        // Moment tensor for crack opening (CLVD + isotropic)
        // For horizontal crack: M_zz dominant
        double M_dot = P0 * A * amp * 2 * PI * f0;  // Rate of moment
        
        // [M_xx, M_yy, M_zz, M_xy, M_xz, M_yz]
        std::array<double, 6> mt;
        mt[0] = M_dot * 0.5;    // M_xx
        mt[1] = M_dot * 0.5;    // M_yy
        mt[2] = M_dot;          // M_zz (dominant for horizontal crack)
        mt[3] = 0.0;            // M_xy
        mt[4] = 0.0;            // M_xz
        mt[5] = 0.0;            // M_yz
        
        moment_tensor_rate[i] = mt;
    }
}

void VolcanicSeismicityModel::generateVLP_source(const VolcanicSeismicEvent& event,
                                                  std::vector<double>& times,
                                                  std::vector<std::array<double, 6>>& moment_tensor_rate) {
    // VLP event - very low frequency conduit resonance
    double f0 = std::max(0.01, event.dominant_frequency * 0.1);  // VLP is ~0.01-0.1 Hz
    double Q = 5.0;  // VLP has low Q
    
    int n_samples = 2000;
    double dt = 0.1;  // Longer time step for low frequency
    times.resize(n_samples);
    moment_tensor_rate.resize(n_samples);
    
    double M0 = event.seismic_moment * 100;  // VLP has larger moment
    
    for (int i = 0; i < n_samples; ++i) {
        times[i] = i * dt;
        double t = times[i];
        
        double amp = decayingOscillation(t, f0, Q);
        double M_dot = M0 * amp * 2 * PI * f0;
        
        // Volumetric source (conduit inflation/deflation)
        std::array<double, 6> mt;
        mt[0] = M_dot / 3.0;   // M_xx
        mt[1] = M_dot / 3.0;   // M_yy
        mt[2] = M_dot / 3.0;   // M_zz
        mt[3] = 0.0;
        mt[4] = 0.0;
        mt[5] = 0.0;
        
        moment_tensor_rate[i] = mt;
    }
}

void VolcanicSeismicityModel::generateTremor_source(const VolcanicTremorSource& source,
                                                     std::vector<double>& times,
                                                     std::vector<double>& amplitude) {
    // Volcanic tremor - continuous oscillation
    double dt = 0.01;
    int n_samples = static_cast<int>(source.duration / dt);
    times.resize(n_samples);
    amplitude.resize(n_samples);
    
    std::mt19937 rng(42);  // Random generator
    std::normal_distribution<double> noise(0.0, source.amplitude_variation);
    
    for (int i = 0; i < n_samples; ++i) {
        times[i] = i * dt;
        double t = times[i];
        
        double amp = 0.0;
        
        // Sum of multiple spectral peaks
        for (size_t k = 0; k < source.peak_frequencies.size(); ++k) {
            double f = source.peak_frequencies[k];
            double A = source.peak_amplitudes[k];
            
            // Narrow-band oscillation with random phase
            double phase = noise(rng) * 2 * PI;
            amp += A * std::sin(2 * PI * f * t + phase);
        }
        
        // Amplitude modulation
        double envelope = 1.0 + source.amplitude_variation * std::sin(0.1 * t);
        
        amplitude[i] = amp * envelope;
    }
}

void VolcanicSeismicityModel::generateExplosionQuake(const VolcanicSeismicEvent& event,
                                                      std::vector<double>& times,
                                                      std::vector<std::array<double, 6>>& moment_tensor_rate) {
    // Explosion quake - impulsive isotropic + single force
    int n_samples = 500;
    double dt = 0.01;
    times.resize(n_samples);
    moment_tensor_rate.resize(n_samples);
    
    double M0 = event.seismic_moment;
    double rise_time = 0.05;  // Very short
    
    for (int i = 0; i < n_samples; ++i) {
        times[i] = i * dt;
        double t = times[i];
        
        double value, deriv;
        applySourceTimeFunction(t, rise_time, value, deriv);
        double M_dot = M0 * deriv;
        
        // Isotropic explosion source
        std::array<double, 6> mt;
        mt[0] = M_dot / 3.0;   // M_xx
        mt[1] = M_dot / 3.0;   // M_yy
        mt[2] = M_dot / 3.0;   // M_zz
        mt[3] = 0.0;
        mt[4] = 0.0;
        mt[5] = 0.0;
        
        moment_tensor_rate[i] = mt;
    }
}

double VolcanicSeismicityModel::CrackResonance::fundamentalFrequency() const {
    // Chouet (1986) crack resonance frequency
    // f = c / (2 * L * sqrt(1 + K))
    // where K = stiffness * L / (rho * c^2 * d)
    
    double c = fluid_velocity;
    double K = stiffness * length / (fluid_density * c * c * width);
    
    return c / (2 * length * std::sqrt(1 + K));
}

double VolcanicSeismicityModel::CrackResonance::qualityFactor() const {
    // Q from viscous damping
    double c = fluid_velocity;
    double f0 = fundamentalFrequency();
    
    return fluid_density * c * width * f0 / (3 * fluid_viscosity);
}

std::vector<double> VolcanicSeismicityModel::CrackResonance::getModeFrequencies(int n_modes) const {
    std::vector<double> freqs(n_modes);
    double f0 = fundamentalFrequency();
    
    for (int n = 0; n < n_modes; ++n) {
        freqs[n] = f0 * (2 * n + 1);  // Odd harmonics for open crack
    }
    return freqs;
}

void VolcanicSeismicityModel::setCrackModel(const CrackResonance& crack) {
    crack_model = crack;
}

double VolcanicSeismicityModel::decayingOscillation(double t, double f, double Q) const {
    if (t < 0) return 0.0;
    
    double omega = 2 * PI * f;
    double decay = std::exp(-omega * t / (2 * Q));
    
    return decay * std::sin(omega * t);
}

void VolcanicSeismicityModel::applySourceTimeFunction(double t, double rise_time,
                                                       double& value, double& derivative) const {
    if (t < 0) {
        value = 0.0;
        derivative = 0.0;
    } else if (t < rise_time) {
        // Smooth ramp up (cosine taper)
        value = 0.5 * (1.0 - std::cos(PI * t / rise_time));
        derivative = 0.5 * PI / rise_time * std::sin(PI * t / rise_time);
    } else if (t < 2 * rise_time) {
        // Ramp down
        double t2 = t - rise_time;
        value = 0.5 * (1.0 + std::cos(PI * t2 / rise_time));
        derivative = -0.5 * PI / rise_time * std::sin(PI * t2 / rise_time);
    } else {
        value = 0.0;
        derivative = 0.0;
    }
}

// =============================================================================
// VolcanicGasModel Implementation
// =============================================================================

VolcanicGasModel::VolcanicGasModel()
    : wind_speed(5.0), wind_direction(0.0), stability_class('D'), mixing_height(1000.0) {}

void VolcanicGasModel::addSource(const GasEmissionSource& source) {
    sources.push_back(source);
}

void VolcanicGasModel::clearSources() {
    sources.clear();
}

void VolcanicGasModel::setWindSpeed(double speed_m_s) {
    wind_speed = std::max(0.5, speed_m_s);  // Minimum wind speed
}

void VolcanicGasModel::setWindDirection(double direction_deg) {
    wind_direction = direction_deg;
}

void VolcanicGasModel::setAtmosphericStability(char stability) {
    if (stability >= 'A' && stability <= 'F') {
        stability_class = stability;
    }
}

void VolcanicGasModel::setMixingHeight(double height_m) {
    mixing_height = std::max(100.0, height_m);
}

double VolcanicGasModel::getConcentration(VolcanicGas gas, double x, double y) const {
    double total_conc = 0.0;
    
    for (const auto& source : sources) {
        // Get emission rate for this gas
        double Q = source.emission_rate;  // Total emission kg/s
        auto it = source.composition.find(gas);
        if (it == source.composition.end()) continue;
        
        double mole_fraction = it->second;
        double MW = getMolecularWeight(gas);
        
        // Mass emission rate of this gas
        double Q_gas = Q * mole_fraction * MW / 
                       (mole_fraction * MW + (1.0 - mole_fraction) * 0.018);  // Normalize
        
        // Transform coordinates to plume-aligned system
        double wind_rad = wind_direction * PI / 180.0;
        double dx = x - source.x;
        double dy = y - source.y;
        
        // Downwind distance
        double x_plume = dx * std::cos(wind_rad) + dy * std::sin(wind_rad);
        // Crosswind distance
        double y_plume = -dx * std::sin(wind_rad) + dy * std::cos(wind_rad);
        
        if (x_plume <= 0) continue;  // Upwind
        
        // Plume rise
        double H = source.z + plumeRise(source, x_plume);
        
        // Dispersion coefficients
        double sy = sigmaY(x_plume);
        double sz = sigmaZ(x_plume);
        
        // Gaussian plume formula (ground reflection included)
        double exp_y = std::exp(-y_plume * y_plume / (2 * sy * sy));
        double exp_z1 = std::exp(-H * H / (2 * sz * sz));  // Direct
        double exp_z2 = std::exp(-H * H / (2 * sz * sz));  // Reflected
        
        double C = Q_gas / (2 * PI * wind_speed * sy * sz) * exp_y * (exp_z1 + exp_z2);
        
        total_conc += C;
    }
    
    return total_conc;  // kg/m³
}

double VolcanicGasModel::getConcentrationPPM(VolcanicGas gas, double x, double y) const {
    double C = getConcentration(gas, x, y);  // kg/m³
    double MW = getMolecularWeight(gas);
    
    // Convert to mol/m³
    double C_mol = C / MW;
    
    // Molar volume at STP
    double V_m = R_GAS * 298.0 / P_ATM;  // m³/mol
    
    // Volume fraction (ppm)
    return C_mol * V_m * 1e6;
}

double VolcanicGasModel::getSO2_column_density(double x, double y) const {
    // Integrate SO2 concentration vertically (simplified)
    // Assumes well-mixed in mixing layer
    
    double C = getConcentration(VolcanicGas::SO2, x, y);
    return C * mixing_height;  // kg/m²
}

double VolcanicGasModel::getTotalSO2_emission() const {
    double total = 0.0;
    for (const auto& source : sources) {
        auto it = source.composition.find(VolcanicGas::SO2);
        if (it != source.composition.end()) {
            double MW_SO2 = getMolecularWeight(VolcanicGas::SO2);
            total += source.emission_rate * it->second * MW_SO2 / 0.018;
        }
    }
    return total;
}

double VolcanicGasModel::getTotalCO2_emission() const {
    double total = 0.0;
    for (const auto& source : sources) {
        auto it = source.composition.find(VolcanicGas::CO2);
        if (it != source.composition.end()) {
            double MW_CO2 = getMolecularWeight(VolcanicGas::CO2);
            total += source.emission_rate * it->second * MW_CO2 / 0.018;
        }
    }
    return total;
}

bool VolcanicGasModel::exceedsHealthLimit(VolcanicGas gas, double x, double y) const {
    double ppm = getConcentrationPPM(gas, x, y);
    
    // WHO/OSHA exposure limits (simplified)
    switch (gas) {
        case VolcanicGas::SO2:
            return ppm > 2.0;  // 2 ppm short-term limit
        case VolcanicGas::CO2:
            return ppm > 5000.0;  // 0.5% = 5000 ppm
        case VolcanicGas::H2S:
            return ppm > 10.0;  // 10 ppm short-term
        case VolcanicGas::HCl:
            return ppm > 5.0;
        case VolcanicGas::HF:
            return ppm > 3.0;
        case VolcanicGas::CO:
            return ppm > 50.0;
        default:
            return false;
    }
}

double VolcanicGasModel::getMolecularWeight(VolcanicGas gas) {
    switch (gas) {
        case VolcanicGas::H2O: return 0.018;
        case VolcanicGas::CO2: return 0.044;
        case VolcanicGas::SO2: return 0.064;
        case VolcanicGas::H2S: return 0.034;
        case VolcanicGas::HCl: return 0.0365;
        case VolcanicGas::HF:  return 0.020;
        case VolcanicGas::CO:  return 0.028;
        case VolcanicGas::H2:  return 0.002;
        case VolcanicGas::He:  return 0.004;
        case VolcanicGas::Rn:  return 0.222;
        default: return 0.029;  // Air
    }
}

double VolcanicGasModel::sigmaY(double x) const {
    // Pasquill-Gifford dispersion coefficients (horizontal)
    // x in meters, x_km for far-field scaling
    double x_km = x / 1000.0;
    
    // Use power-law form for near-field (x < 10 km), adjust for far-field
    double sigma;
    if (x < 10000.0) {
        // Near-field: standard Pasquill-Gifford
        switch (stability_class) {
            case 'A':  // Very unstable
                sigma = 0.22 * x * std::pow(1 + 0.0001 * x, -0.5);
                break;
            case 'B':  // Unstable
                sigma = 0.16 * x * std::pow(1 + 0.0001 * x, -0.5);
                break;
            case 'C':  // Slightly unstable
                sigma = 0.11 * x * std::pow(1 + 0.0001 * x, -0.5);
                break;
            case 'D':  // Neutral
                sigma = 0.08 * x * std::pow(1 + 0.0001 * x, -0.5);
                break;
            case 'E':  // Slightly stable
                sigma = 0.06 * x * std::pow(1 + 0.0001 * x, -0.5);
                break;
            case 'F':  // Stable
                sigma = 0.04 * x * std::pow(1 + 0.0001 * x, -0.5);
                break;
            default:
                sigma = 0.08 * x * std::pow(1 + 0.0001 * x, -0.5);
        }
    } else {
        // Far-field: reduced growth rate (volcanic plumes disperse differently)
        double sigma_10km = sigmaY(10000.0);  // Value at 10 km
        double growth_exponent = (stability_class <= 'C') ? 0.7 : 0.5;  // Slower growth
        sigma = sigma_10km * std::pow(x_km / 10.0, growth_exponent);
    }
    
    return sigma;
}

double VolcanicGasModel::sigmaZ(double x) const {
    // Pasquill-Gifford dispersion coefficients (vertical)
    switch (stability_class) {
        case 'A':
            return 0.20 * x;
        case 'B':
            return 0.12 * x;
        case 'C':
            return 0.08 * x * std::pow(1 + 0.0002 * x, -0.5);
        case 'D':
            return 0.06 * x * std::pow(1 + 0.0015 * x, -0.5);
        case 'E':
            return 0.03 * x * std::pow(1 + 0.0003 * x, -1.0);
        case 'F':
            return 0.016 * x * std::pow(1 + 0.0003 * x, -1.0);
        default:
            return 0.06 * x * std::pow(1 + 0.0015 * x, -0.5);
    }
}

double VolcanicGasModel::plumeRise(const GasEmissionSource& source, double x) const {
    // Briggs plume rise equations
    double T_a = 288.0;  // Ambient temperature
    double dT = source.temperature - T_a;
    
    if (dT <= 0) {
        // No buoyancy
        return 0.0;
    }
    
    // Buoyancy flux
    double F = G * source.velocity * source.emission_rate * dT / 
               (PI * source.temperature * 1000.0);  // Approximate
    
    // Final rise (neutral conditions)
    double dH = 1.6 * std::pow(F, 1.0/3.0) * std::pow(x, 2.0/3.0) / wind_speed;
    
    // Limit by mixing height
    return std::min(dH, mixing_height - source.z);
}

// =============================================================================
// TephraParticle Implementation
// =============================================================================

double TephraParticle::settlingVelocity(double air_density, double air_viscosity) const {
    // Wilson & Huang (1979) settling velocity
    double d = diameter;
    double rho_p = density;
    double rho_a = air_density;
    double mu = air_viscosity;
    
    // Reynolds number iteration
    double Re = 1.0;
    double Cd, w;
    
    for (int iter = 0; iter < 10; ++iter) {
        Cd = dragCoefficient(Re);
        w = std::sqrt(4 * G * d * (rho_p - rho_a) / (3 * Cd * rho_a));
        Re = rho_a * w * d / mu;
    }
    
    // Apply shape factor correction
    return w * std::pow(shape_factor, 0.5);
}

double TephraParticle::dragCoefficient(double Re) const {
    // Drag coefficient for spheres (Morsi & Alexander, 1972)
    if (Re < 0.1) {
        return 24.0 / Re;  // Stokes regime
    } else if (Re < 1.0) {
        return 24.0 / Re * (1 + 0.15 * std::pow(Re, 0.687));
    } else if (Re < 500) {
        return 24.0 / Re * (1 + 0.15 * std::pow(Re, 0.687));
    } else if (Re < 2e5) {
        return 0.44;  // Newton regime
    } else {
        return 0.1;  // Turbulent
    }
}

// =============================================================================
// TephraDispersalModel Implementation
// =============================================================================

TephraDispersalModel::TephraDispersalModel()
    : nx(100), ny(100), nz(50), 
      dx(1000), dy(1000), dz(500),
      xmin(-50000), xmax(50000), ymin(-50000), ymax(50000), zmin(0), zmax(25000),
      n_size_bins(10) {}

void TephraDispersalModel::initialize(const TephraSourceParameters& src) {
    source = src;
    
    // Set up grid
    dx = (xmax - xmin) / nx;
    dy = (ymax - ymin) / ny;
    dz = (zmax - zmin) / nz;
    
    // Initialize particle size distribution
    initializeParticleSizes();
    
    // Initialize concentration and deposit arrays
    int n_2d = nx * ny;
    int n_3d = nx * ny * nz;
    
    concentration.resize(n_size_bins);
    deposit.resize(n_size_bins);
    
    for (int k = 0; k < n_size_bins; ++k) {
        concentration[k].assign(n_3d, 0.0);
        deposit[k].assign(n_2d, 0.0);
    }
    
    // Default wind field (uniform)
    wind_field = [](double x, double y, double z, double& u, double& v, double& w) {
        u = 10.0;  // 10 m/s eastward
        v = 0.0;
        w = 0.0;
    };
    
    // Default atmosphere
    air_density = [](double z) {
        return 1.225 * std::exp(-z / 8500.0);  // Standard atmosphere
    };
    
    air_viscosity = [](double z) {
        return 1.8e-5;  // Approximately constant
    };
}

void TephraDispersalModel::setWindField(
    std::function<void(double, double, double, double&, double&, double&)> wind) {
    wind_field = wind;
}

void TephraDispersalModel::setAtmosphere(std::function<double(double)> density,
                                          std::function<double(double)> viscosity) {
    air_density = density;
    air_viscosity = viscosity;
}

void TephraDispersalModel::run(double duration) {
    double dt = 10.0;  // Time step (seconds)
    double time = 0.0;
    
    while (time < duration) {
        // Add source
        int i_src = nx / 2;
        int j_src = ny / 2;
        int k_src = static_cast<int>((source.column_height - zmin) / dz);
        k_src = std::min(k_src, nz - 1);
        
        if (k_src >= 0) {
            int idx_3d = i_src + j_src * nx + k_src * nx * ny;
            
            for (int m = 0; m < n_size_bins; ++m) {
                double mass_rate = source.mass_eruption_rate * mass_fractions[m];
                double cell_volume = dx * dy * dz;
                concentration[m][idx_3d] += mass_rate * dt / cell_volume;
            }
        }
        
        // Advect
        advect(dt);
        
        // Diffuse
        diffuse(dt);
        
        // Settle
        settle(dt);
        
        // Deposit
        deposit_particles();
        
        time += dt;
    }
}

double TephraDispersalModel::getDepositThickness(double x, double y) const {
    // Sum deposits from all size bins
    double loading = getDepositLoading(x, y);  // kg/m²
    double bulk_density = 1000.0;  // kg/m³ for loose tephra
    return loading / bulk_density;
}

double TephraDispersalModel::getDepositLoading(double x, double y) const {
    int i = static_cast<int>((x - xmin) / dx);
    int j = static_cast<int>((y - ymin) / dy);
    
    if (i < 0 || i >= nx || j < 0 || j >= ny) return 0.0;
    
    int idx_2d = i + j * nx;
    
    double total = 0.0;
    for (int m = 0; m < n_size_bins; ++m) {
        total += deposit[m][idx_2d];
    }
    return total;
}

double TephraDispersalModel::getAirborneConcentration(double x, double y, double z) const {
    int i = static_cast<int>((x - xmin) / dx);
    int j = static_cast<int>((y - ymin) / dy);
    int k = static_cast<int>((z - zmin) / dz);
    
    if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return 0.0;
    
    int idx_3d = i + j * nx + k * nx * ny;
    
    double total = 0.0;
    for (int m = 0; m < n_size_bins; ++m) {
        total += concentration[m][idx_3d];
    }
    return total;
}

bool TephraDispersalModel::exceedsAviationLimit(double x, double y, double z) const {
    // Aviation limit: 2 mg/m³
    double C = getAirborneConcentration(x, y, z);
    return C > 2e-6;  // 2 mg/m³ = 2e-6 kg/m³
}

double TephraDispersalModel::getVisibility(double x, double y, double z) const {
    // Koschmieder equation: visibility = 3.912 / extinction
    double C = getAirborneConcentration(x, y, z);
    
    if (C < 1e-9) return 100000.0;  // 100 km clear
    
    // Mass extinction coefficient for volcanic ash ~ 1-5 m²/g
    double k_ext = 3.0 * 1000.0;  // m²/kg
    double extinction = k_ext * C;
    
    return 3.912 / extinction;
}

void TephraDispersalModel::getIsopachContours(
    const std::vector<double>& thicknesses,
    std::vector<std::vector<std::pair<double, double>>>& contours) const {
    
    // Simplified contour extraction
    contours.resize(thicknesses.size());
    
    for (size_t t = 0; t < thicknesses.size(); ++t) {
        double target = thicknesses[t];
        
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                double x = xmin + i * dx;
                double y = ymin + j * dy;
                double h = getDepositThickness(x, y);
                double h_xp = getDepositThickness(x + dx, y);
                double h_xm = getDepositThickness(x - dx, y);
                double h_yp = getDepositThickness(x, y + dy);
                double h_ym = getDepositThickness(x, y - dy);
                
                // Check for contour crossing
                if ((h >= target && (h_xp < target || h_xm < target || 
                                     h_yp < target || h_ym < target)) ||
                    (h < target && (h_xp >= target || h_xm >= target || 
                                    h_yp >= target || h_ym >= target))) {
                    contours[t].push_back({x, y});
                }
            }
        }
    }
}

void TephraDispersalModel::advect(double dt) {
    // First-order upwind advection
    for (int m = 0; m < n_size_bins; ++m) {
        std::vector<double> new_conc = concentration[m];
        
        for (int k = 1; k < nz - 1; ++k) {
            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {
                    double x = xmin + i * dx;
                    double y = ymin + j * dy;
                    double z = zmin + k * dz;
                    
                    double u, v, w;
                    wind_field(x, y, z, u, v, w);
                    
                    int idx = i + j * nx + k * nx * ny;
                    int idx_xm = (i - 1) + j * nx + k * nx * ny;
                    int idx_xp = (i + 1) + j * nx + k * nx * ny;
                    int idx_ym = i + (j - 1) * nx + k * nx * ny;
                    int idx_yp = i + (j + 1) * nx + k * nx * ny;
                    int idx_zm = i + j * nx + (k - 1) * nx * ny;
                    int idx_zp = i + j * nx + (k + 1) * nx * ny;
                    
                    // Upwind
                    double dC_dx = (u > 0) ? (concentration[m][idx] - concentration[m][idx_xm]) / dx
                                          : (concentration[m][idx_xp] - concentration[m][idx]) / dx;
                    double dC_dy = (v > 0) ? (concentration[m][idx] - concentration[m][idx_ym]) / dy
                                          : (concentration[m][idx_yp] - concentration[m][idx]) / dy;
                    double dC_dz = (w > 0) ? (concentration[m][idx] - concentration[m][idx_zm]) / dz
                                          : (concentration[m][idx_zp] - concentration[m][idx]) / dz;
                    
                    new_conc[idx] -= dt * (u * dC_dx + v * dC_dy + w * dC_dz);
                    new_conc[idx] = std::max(0.0, new_conc[idx]);
                }
            }
        }
        
        concentration[m] = new_conc;
    }
}

void TephraDispersalModel::diffuse(double dt) {
    // Simple diffusion
    double K_h = 1000.0;  // Horizontal diffusivity (m²/s)
    double K_v = 10.0;    // Vertical diffusivity (m²/s)
    
    for (int m = 0; m < n_size_bins; ++m) {
        std::vector<double> new_conc = concentration[m];
        
        for (int k = 1; k < nz - 1; ++k) {
            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {
                    int idx = i + j * nx + k * nx * ny;
                    int idx_xm = (i - 1) + j * nx + k * nx * ny;
                    int idx_xp = (i + 1) + j * nx + k * nx * ny;
                    int idx_ym = i + (j - 1) * nx + k * nx * ny;
                    int idx_yp = i + (j + 1) * nx + k * nx * ny;
                    int idx_zm = i + j * nx + (k - 1) * nx * ny;
                    int idx_zp = i + j * nx + (k + 1) * nx * ny;
                    
                    double laplacian_h = (concentration[m][idx_xp] + concentration[m][idx_xm] - 2 * concentration[m][idx]) / (dx * dx)
                                       + (concentration[m][idx_yp] + concentration[m][idx_ym] - 2 * concentration[m][idx]) / (dy * dy);
                    double laplacian_v = (concentration[m][idx_zp] + concentration[m][idx_zm] - 2 * concentration[m][idx]) / (dz * dz);
                    
                    new_conc[idx] += dt * (K_h * laplacian_h + K_v * laplacian_v);
                    new_conc[idx] = std::max(0.0, new_conc[idx]);
                }
            }
        }
        
        concentration[m] = new_conc;
    }
}

void TephraDispersalModel::settle(double dt) {
    // Gravitational settling
    for (int m = 0; m < n_size_bins; ++m) {
        std::vector<double> new_conc = concentration[m];
        
        for (int k = 0; k < nz; ++k) {
            double z = zmin + k * dz;
            double rho_a = air_density(z);
            double mu = air_viscosity(z);
            double ws = particles[m].settlingVelocity(rho_a, mu);
            
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int idx = i + j * nx + k * nx * ny;
                    
                    // Settling flux
                    double flux_out = ws * concentration[m][idx];
                    new_conc[idx] -= flux_out * dt / dz;
                    
                    // Add to cell below
                    if (k > 0) {
                        int idx_below = i + j * nx + (k - 1) * nx * ny;
                        new_conc[idx_below] += flux_out * dt / dz;
                    }
                    
                    new_conc[idx] = std::max(0.0, new_conc[idx]);
                }
            }
        }
        
        concentration[m] = new_conc;
    }
}

void TephraDispersalModel::deposit_particles() {
    // Ground deposition from lowest layer
    for (int m = 0; m < n_size_bins; ++m) {
        double z = zmin;
        double rho_a = air_density(z);
        double mu = air_viscosity(z);
        double ws = particles[m].settlingVelocity(rho_a, mu);
        
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx_3d = i + j * nx;  // Lowest layer (k=0)
                int idx_2d = i + j * nx;
                
                // Deposition flux
                double flux = ws * concentration[m][idx_3d];
                deposit[m][idx_2d] += flux * 1.0;  // Assume 1 second accumulation
            }
        }
    }
}

void TephraDispersalModel::initializeParticleSizes() {
    // Initialize particle size bins (phi scale)
    particles.resize(n_size_bins);
    mass_fractions.resize(n_size_bins);
    
    double phi_min = source.phi_mean - 3 * source.phi_stddev;
    double phi_max = source.phi_mean + 3 * source.phi_stddev;
    double dphi = (phi_max - phi_min) / n_size_bins;
    
    double total_mass = 0.0;
    
    for (int m = 0; m < n_size_bins; ++m) {
        double phi = phi_min + (m + 0.5) * dphi;
        double d_mm = VolcanoUtils::phi_to_mm(phi);
        
        particles[m].diameter = d_mm / 1000.0;  // Convert to meters
        particles[m].density = source.particle_density;
        particles[m].shape_factor = 0.7;  // Typical for volcanic ash
        
        // Gaussian mass distribution
        double z = (phi - source.phi_mean) / source.phi_stddev;
        mass_fractions[m] = std::exp(-0.5 * z * z);
        total_mass += mass_fractions[m];
    }
    
    // Normalize
    for (int m = 0; m < n_size_bins; ++m) {
        mass_fractions[m] /= total_mass;
    }
}

// =============================================================================
// CoupledVolcanoSystem Implementation
// =============================================================================

CoupledVolcanoSystem::CoupledVolcanoSystem()
    : current_time(0.0), eruption_active(false), 
      current_eruption_type(EruptionType::EFFUSIVE) {}

void CoupledVolcanoSystem::initialize(const VolcanoSystemConfig& cfg,
                                       const MagmaChamberGeometry& chamber_geometry,
                                       const MagmaProperties& mag,
                                       const ConduitGeometry& cond) {
    config = cfg;
    
    // Initialize component models based on configuration
    if (config.enable_magma_chamber) {
        chamber = std::make_unique<MagmaChamberModel>();
        MagmaChamberState initial_state;
        initial_state.pressure = 200e6;
        initial_state.temperature = mag.temperature;
        initial_state.crystal_fraction = mag.crystal_fraction;
        chamber->initialize(chamber_geometry, mag, initial_state);
    }
    
    if (config.enable_conduit_flow) {
        conduit = std::make_unique<ConduitFlowModel>();
        conduit->initialize(cond, mag, 200e6);
    }
    
    if (config.enable_eruption_column) {
        column = std::make_unique<EruptionColumnModel>();
    }
    
    if (config.enable_pdc) {
        pdc = std::make_unique<PDCModel>();
    }
    
    if (config.enable_lava_flow) {
        lava = std::make_unique<LavaFlowModel>();
        LavaSourceParameters lava_src;
        lava_src.vent_x = 0.0;
        lava_src.vent_y = 0.0;
        lava_src.temperature = mag.temperature;
        lava->initialize(lava_src, mag);
    }
    
    if (config.enable_lahar) {
        lahar = std::make_unique<LaharModel>();
    }
    
    if (config.enable_deformation) {
        deformation = std::make_unique<VolcanicDeformationModel>();
        
        // Add source from chamber
        VolcanicDeformationSource src;
        src.type = chamber_geometry.shape;
        src.x = chamber_geometry.x_center;
        src.y = chamber_geometry.y_center;
        src.depth = chamber_geometry.depth;
        src.radius = chamber_geometry.effectiveRadius();
        deformation->addSource(src);
    }
    
    if (config.enable_seismicity) {
        seismicity = std::make_unique<VolcanicSeismicityModel>();
    }
    
    if (config.enable_gas_emission) {
        gas = std::make_unique<VolcanicGasModel>();
    }
    
    if (config.enable_tephra) {
        tephra = std::make_unique<TephraDispersalModel>();
    }
    
    current_time = 0.0;
    monitoring.clear();
}

void CoupledVolcanoSystem::setTopography(std::function<double(double, double)> elevation) {
    topography = elevation;
    
    if (lava) lava->setTopography(elevation);
    if (pdc) pdc->setTopography(elevation);
    if (lahar) lahar->setTopography(elevation);
}

void CoupledVolcanoSystem::setAtmosphere(
    std::function<void(double, double&, double&, double&)> profile) {
    // Store atmosphere profile for column and tephra models
}

void CoupledVolcanoSystem::setRechargeRate(double rate_kg_s) {
    if (chamber) {
        chamber->setRechargeRate(rate_kg_s);
    }
}

void CoupledVolcanoSystem::triggerEruption(EruptionType type) {
    eruption_active = true;
    current_eruption_type = type;
    
    // Initialize eruption-specific models
    switch (type) {
        case EruptionType::PLINIAN:
        case EruptionType::SUBPLINIAN:
        case EruptionType::ULTRAPLINIAN:
            if (column && conduit) {
                EruptionColumnParameters params;
                params.mass_flux = conduit->getMassFlux();
                params.exit_velocity = conduit->getExitVelocity();
                column->initialize(params);
            }
            break;
            
        case EruptionType::EFFUSIVE:
        case EruptionType::HAWAIIAN:
            // Lava flow dominated
            if (lava) {
                lava->setEffusionRate(10.0);  // Default 10 m³/s
            }
            break;
            
        case EruptionType::DOME_COLLAPSE:
            // Initialize PDC from dome collapse
            if (pdc) {
                PDCSourceConditions src;
                src.type = PDCType::DOME_COLLAPSE;
                src.initial_velocity = 30.0;
                src.initial_temperature = 800.0;
                pdc->initialize(src);
            }
            break;
            
        default:
            break;
    }
}

void CoupledVolcanoSystem::update(double dt) {
    // Update chamber
    if (chamber && config.enable_magma_chamber) {
        if (eruption_active && conduit) {
            chamber->setEruptionRate(conduit->getMassFlux());
        }
        chamber->update(dt);
    }
    
    // Update conduit flow
    if (conduit && config.enable_conduit_flow) {
        if (config.couple_chamber_conduit && chamber) {
            conduit->setChamberPressure(chamber->getPressure());
        }
        conduit->update(dt);
    }
    
    // Update eruption column
    if (column && config.enable_eruption_column && eruption_active) {
        if (config.couple_conduit_column && conduit) {
            EruptionColumnParameters params;
            params.mass_flux = conduit->getMassFlux();
            params.exit_velocity = conduit->getExitVelocity();
            column->initialize(params);
        }
        column->solve();
        
        // Check for column collapse -> PDC
        if (config.couple_column_pdc && column->doesCollapse() && pdc) {
            PDCSourceConditions src;
            src.type = PDCType::COLUMN_COLLAPSE;
            src.mass_flux = column->getCollapseFlux();
            src.source_height = column->getCollapseHeight();
            pdc->initialize(src);
        }
    }
    
    // Update surface flows
    if (pdc && config.enable_pdc && eruption_active) {
        pdc->update(dt);
    }
    
    if (lava && config.enable_lava_flow && eruption_active) {
        lava->update(dt);
    }
    
    if (lahar && config.enable_lahar) {
        lahar->update(dt);
    }
    
    // Update deformation
    if (deformation && config.enable_deformation) {
        if (config.couple_chamber_deformation && chamber) {
            updateDeformationFromChamber();
        }
    }
    
    // Update gas emissions
    if (gas && config.enable_gas_emission) {
        // Update gas sources based on eruption state
    }
    
    // Update tephra dispersal
    if (tephra && config.enable_tephra && eruption_active) {
        if (column) {
            TephraSourceParameters src;
            src.column_height = column->getMaxHeight();
            src.mass_eruption_rate = column->getAsh_emission_rate();
            tephra->initialize(src);
        }
        tephra->run(dt);
    }
    
    // Record monitoring data
    recordMonitoring();
    
    current_time += dt;
}

void CoupledVolcanoSystem::run() {
    double dt = config.dt_min;
    
    while (current_time < config.end_time) {
        update(dt);
        
        // Adaptive time stepping (simplified)
        if (eruption_active) {
            dt = config.dt_min;
        } else {
            dt = config.dt_max;
        }
    }
}

double CoupledVolcanoSystem::getMassEruptionRate() const {
    if (conduit) {
        return conduit->getMassFlux();
    }
    return 0.0;
}

double CoupledVolcanoSystem::getColumnHeight() const {
    if (column) {
        return column->getMaxHeight();
    }
    return 0.0;
}

void CoupledVolcanoSystem::updateChamberConduitCoupling() {
    if (!chamber || !conduit) return;
    conduit->setChamberPressure(chamber->getPressure());
}

void CoupledVolcanoSystem::updateConduitColumnCoupling() {
    if (!conduit || !column) return;
    
    EruptionColumnParameters params;
    params.mass_flux = conduit->getMassFlux();
    params.exit_velocity = conduit->getExitVelocity();
    params.exit_temperature = 1273.0;  // Typical
    column->initialize(params);
}

void CoupledVolcanoSystem::updateColumnPDCCoupling() {
    if (!column || !pdc) return;
    
    if (column->doesCollapse()) {
        PDCSourceConditions src;
        src.type = PDCType::COLUMN_COLLAPSE;
        src.mass_flux = column->getCollapseFlux();
        src.source_height = column->getCollapseHeight();
        pdc->initialize(src);
    }
}

void CoupledVolcanoSystem::updateDeformationFromChamber() {
    if (!deformation || !chamber) return;
    
    // Update deformation source based on chamber overpressure
    deformation->clearSources();
    
    VolcanicDeformationSource src;
    src.type = DeformationSourceType::MOGI;
    src.depth = 5000.0;  // Default depth
    src.volume_change = chamber->getState().overpressure * chamber->getState().volume * 
                        chamber->getMagmaCompressibility();
    
    deformation->addSource(src);
}

void CoupledVolcanoSystem::updateSeismicityFromChamber() {
    // Generate seismic events based on chamber state
}

void CoupledVolcanoSystem::recordMonitoring() {
    MonitoringOutput output;
    output.time = current_time;
    
    if (chamber) {
        output.chamber_pressure = chamber->getPressure();
        output.chamber_overpressure = chamber->getOverpressure();
    }
    
    if (deformation) {
        deformation->computeDisplacement(1000.0, 0.0, 
            output.radial_displacement, 
            output.radial_displacement,  // Placeholder
            output.vertical_displacement);
    }
    
    if (gas) {
        output.SO2_emission_rate = gas->getTotalSO2_emission();
    }
    
    if (conduit) {
        output.mass_eruption_rate = conduit->getMassFlux();
    }
    
    if (column) {
        output.column_height = column->getMaxHeight();
    }
    
    monitoring.push_back(output);
}

void CoupledVolcanoSystem::writeOutput(const std::string& output_dir) const {
    // Write monitoring time series
    std::ofstream file(output_dir + "/monitoring_timeseries.dat");
    file << "# Volcano monitoring time series\n";
    file << "# time(s) P_chamber(Pa) overpressure(Pa) u_r(m) u_z(m) SO2(kg/s) MER(kg/s) H_col(m)\n";
    
    for (const auto& m : monitoring) {
        file << std::scientific << std::setprecision(6)
             << m.time << " "
             << m.chamber_pressure << " "
             << m.chamber_overpressure << " "
             << m.radial_displacement << " "
             << m.vertical_displacement << " "
             << m.SO2_emission_rate << " "
             << m.mass_eruption_rate << " "
             << m.column_height << "\n";
    }
}

// =============================================================================
// VolcanoScenarios Implementation
// =============================================================================

namespace VolcanoScenarios {

VolcanoSystemConfig mountStHelens1980() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_eruption_column = true;
    config.enable_pdc = true;
    config.enable_deformation = true;
    config.enable_seismicity = true;
    config.enable_gas_emission = true;
    config.enable_tephra = true;
    config.end_time = 36000.0;  // 10 hours
    return config;
}

MagmaChamberGeometry mountStHelens1980_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 7000.0;
    geom.volume = 0.5e9;  // 0.5 km³
    geom.semi_axis_a = 1000.0;
    geom.semi_axis_b = 1000.0;
    geom.semi_axis_c = 500.0;
    return geom;
}

MagmaProperties mountStHelens1980_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::DACITE;
    mag.SiO2 = 64.0;
    mag.temperature = 1123.0;
    mag.H2O_total = 4.6;
    mag.crystal_fraction = 0.3;
    return mag;
}

VolcanoSystemConfig pinatubo1991() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_eruption_column = true;
    config.enable_pdc = true;
    config.enable_deformation = true;
    config.enable_gas_emission = true;
    config.enable_tephra = true;
    config.end_time = 86400.0;  // 24 hours
    return config;
}

MagmaChamberGeometry pinatubo1991_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 8000.0;
    geom.volume = 5e9;  // 5 km³
    geom.semi_axis_a = 2000.0;
    geom.semi_axis_b = 2000.0;
    geom.semi_axis_c = 1000.0;
    return geom;
}

MagmaProperties pinatubo1991_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::DACITE;
    mag.SiO2 = 65.0;
    mag.temperature = 1053.0;
    mag.H2O_total = 6.0;
    mag.SO2_total = 0.15;
    mag.crystal_fraction = 0.25;
    return mag;
}

VolcanoSystemConfig kilauea_effusive() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_lava_flow = true;
    config.enable_deformation = true;
    config.enable_gas_emission = true;
    config.enable_eruption_column = false;
    config.enable_pdc = false;
    config.enable_tephra = false;
    config.end_time = 604800.0;  // 1 week
    return config;
}

MagmaChamberGeometry kilauea_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 3000.0;
    geom.volume = 1e9;
    geom.shape = DeformationSourceType::OBLATE_SPHEROID;
    return geom;
}

MagmaProperties kilauea_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::BASALT;
    mag.SiO2 = 50.0;
    mag.temperature = 1423.0;
    mag.H2O_total = 0.5;
    mag.crystal_fraction = 0.05;
    return mag;
}

VolcanoSystemConfig stromboli_persistent() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_gas_emission = true;
    config.enable_seismicity = true;
    config.end_time = 3600.0;
    return config;
}

MagmaChamberGeometry stromboli_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 2000.0;
    geom.volume = 0.1e9;
    return geom;
}

MagmaProperties stromboli_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::BASALT;
    mag.SiO2 = 49.0;
    mag.temperature = 1373.0;
    mag.H2O_total = 2.5;
    return mag;
}

VolcanoSystemConfig yellowstone_unrest() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_deformation = true;
    config.enable_seismicity = true;
    config.enable_gas_emission = true;
    config.enable_conduit_flow = false;
    config.enable_eruption_column = false;
    config.end_time = 31536000.0;  // 1 year
    return config;
}

MagmaChamberGeometry yellowstone_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 8000.0;
    geom.volume = 4000e9;  // 4000 km³
    geom.shape = DeformationSourceType::HORIZONTAL_SILL;
    return geom;
}

MagmaProperties yellowstone_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::RHYOLITE;
    mag.SiO2 = 77.0;
    mag.temperature = 1073.0;
    mag.H2O_total = 4.0;
    mag.crystal_fraction = 0.4;
    return mag;
}

VolcanoSystemConfig campi_flegrei_unrest() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_deformation = true;
    config.enable_seismicity = true;
    config.enable_gas_emission = true;
    config.end_time = 31536000.0;  // 1 year
    return config;
}

MagmaChamberGeometry campi_flegrei_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 4000.0;
    geom.volume = 10e9;
    geom.shape = DeformationSourceType::HORIZONTAL_SILL;
    return geom;
}

MagmaProperties campi_flegrei_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::PHONOLITE;
    mag.SiO2 = 60.0;
    mag.temperature = 1123.0;
    mag.H2O_total = 3.0;
    return mag;
}

VolcanoSystemConfig merapi_dome() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_pdc = true;
    config.enable_deformation = true;
    config.enable_seismicity = true;
    config.end_time = 86400.0;
    return config;
}

MagmaChamberGeometry merapi_chamber() {
    MagmaChamberGeometry geom;
    geom.depth = 5000.0;
    geom.volume = 0.5e9;
    return geom;
}

MagmaProperties merapi_magma() {
    MagmaProperties mag;
    mag.composition = MagmaComposition::ANDESITE;
    mag.SiO2 = 55.0;
    mag.temperature = 1173.0;
    mag.H2O_total = 3.0;
    mag.crystal_fraction = 0.4;
    return mag;
}

VolcanoSystemConfig generic_VEI4() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_eruption_column = true;
    config.enable_tephra = true;
    config.enable_deformation = true;
    config.end_time = 14400.0;  // 4 hours
    return config;
}

VolcanoSystemConfig generic_VEI6() {
    VolcanoSystemConfig config;
    config.enable_magma_chamber = true;
    config.enable_conduit_flow = true;
    config.enable_eruption_column = true;
    config.enable_pdc = true;
    config.enable_tephra = true;
    config.enable_deformation = true;
    config.enable_gas_emission = true;
    config.end_time = 86400.0;  // 24 hours
    return config;
}

} // namespace VolcanoScenarios

// =============================================================================
// VolcanoConfigReader Implementation
// =============================================================================

bool VolcanoConfigReader::parseVolcanoConfig(const std::map<std::string, std::string>& section,
                                              VolcanoSystemConfig& config) {
    auto getBool = [&section](const std::string& key, bool default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return it->second == "true" || it->second == "1" || it->second == "yes";
    };
    
    auto getDouble = [&section](const std::string& key, double default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return std::stod(it->second);
    };
    
    config.enable_magma_chamber = getBool("enable_magma_chamber", true);
    config.enable_conduit_flow = getBool("enable_conduit_flow", true);
    config.enable_eruption_column = getBool("enable_eruption_column", true);
    config.enable_pdc = getBool("enable_pdc", false);
    config.enable_lava_flow = getBool("enable_lava_flow", false);
    config.enable_lahar = getBool("enable_lahar", false);
    config.enable_deformation = getBool("enable_deformation", true);
    config.enable_seismicity = getBool("enable_seismicity", true);
    config.enable_gas_emission = getBool("enable_gas_emission", true);
    config.enable_tephra = getBool("enable_tephra", false);
    
    config.couple_chamber_conduit = getBool("couple_chamber_conduit", true);
    config.couple_conduit_column = getBool("couple_conduit_column", true);
    config.couple_column_pdc = getBool("couple_column_pdc", true);
    config.couple_chamber_deformation = getBool("couple_chamber_deformation", true);
    config.couple_chamber_seismicity = getBool("couple_chamber_seismicity", true);
    
    config.dt_min = getDouble("dt_min", 0.01);
    config.dt_max = getDouble("dt_max", 60.0);
    config.end_time = getDouble("end_time", 3600.0);
    config.output_interval = static_cast<int>(getDouble("output_interval", 100));
    
    auto it = section.find("output_directory");
    if (it != section.end()) {
        config.output_directory = it->second;
    }
    
    return true;
}

bool VolcanoConfigReader::parseChamberConfig(const std::map<std::string, std::string>& section,
                                              MagmaChamberGeometry& geometry,
                                              MagmaChamberState& state) {
    auto getDouble = [&section](const std::string& key, double default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return std::stod(it->second);
    };
    
    geometry.x_center = getDouble("x_center", 0.0);
    geometry.y_center = getDouble("y_center", 0.0);
    geometry.depth = getDouble("depth", 5000.0);
    geometry.volume = getDouble("volume", 10e9);
    geometry.semi_axis_a = getDouble("semi_axis_a", 2000.0);
    geometry.semi_axis_b = getDouble("semi_axis_b", 2000.0);
    geometry.semi_axis_c = getDouble("semi_axis_c", 1000.0);
    
    state.pressure = getDouble("pressure", 200e6);
    state.temperature = getDouble("temperature", 1173.15);
    state.crystal_fraction = getDouble("crystal_fraction", 0.3);
    state.dissolved_H2O = getDouble("dissolved_H2O", 3.5);
    state.overpressure = getDouble("overpressure", 10e6);
    
    // Parse shape type
    auto it = section.find("shape");
    if (it != section.end()) {
        if (it->second == "MOGI") geometry.shape = DeformationSourceType::MOGI;
        else if (it->second == "MCTIGUE") geometry.shape = DeformationSourceType::MCTIGUE;
        else if (it->second == "PROLATE") geometry.shape = DeformationSourceType::PROLATE_SPHEROID;
        else if (it->second == "OBLATE") geometry.shape = DeformationSourceType::OBLATE_SPHEROID;
        else if (it->second == "SILL") geometry.shape = DeformationSourceType::HORIZONTAL_SILL;
    }
    
    return true;
}

bool VolcanoConfigReader::parseMagmaProperties(const std::map<std::string, std::string>& section,
                                                MagmaProperties& props) {
    auto getDouble = [&section](const std::string& key, double default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return std::stod(it->second);
    };
    
    props.SiO2 = getDouble("SiO2", 60.0);
    props.Al2O3 = getDouble("Al2O3", 17.0);
    props.FeO = getDouble("FeO", 6.0);
    props.MgO = getDouble("MgO", 3.0);
    props.CaO = getDouble("CaO", 6.0);
    props.Na2O = getDouble("Na2O", 4.0);
    props.K2O = getDouble("K2O", 2.0);
    
    props.H2O_total = getDouble("H2O_total", 4.0);
    props.CO2_total = getDouble("CO2_total", 0.1);
    props.SO2_total = getDouble("SO2_total", 0.05);
    
    props.temperature = getDouble("temperature", 1273.15);
    props.crystal_fraction = getDouble("crystal_fraction", 0.2);
    
    // Parse composition enum
    auto it = section.find("composition");
    if (it != section.end()) {
        if (it->second == "BASALT") props.composition = MagmaComposition::BASALT;
        else if (it->second == "ANDESITE") props.composition = MagmaComposition::ANDESITE;
        else if (it->second == "DACITE") props.composition = MagmaComposition::DACITE;
        else if (it->second == "RHYOLITE") props.composition = MagmaComposition::RHYOLITE;
    }
    
    return true;
}

bool VolcanoConfigReader::parseConduitConfig(const std::map<std::string, std::string>& section,
                                              ConduitGeometry& geometry) {
    auto getDouble = [&section](const std::string& key, double default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return std::stod(it->second);
    };
    
    auto getBool = [&section](const std::string& key, bool default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return it->second == "true" || it->second == "1";
    };
    
    geometry.depth = getDouble("depth", 5000.0);
    geometry.length = getDouble("length", 5000.0);
    geometry.radius_base = getDouble("radius_base", 15.0);
    geometry.radius_vent = getDouble("radius_vent", 10.0);
    geometry.has_constriction = getBool("has_constriction", false);
    geometry.constriction_depth = getDouble("constriction_depth", 2000.0);
    geometry.constriction_radius = getDouble("constriction_radius", 5.0);
    
    return true;
}

bool VolcanoConfigReader::parseDeformationSource(const std::map<std::string, std::string>& section,
                                                  VolcanicDeformationSource& source) {
    auto getDouble = [&section](const std::string& key, double default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return std::stod(it->second);
    };
    
    source.x = getDouble("x", 0.0);
    source.y = getDouble("y", 0.0);
    source.depth = getDouble("depth", 5000.0);
    source.volume_change = getDouble("volume_change", 1e6);
    source.pressure_change = getDouble("pressure_change", 10e6);
    source.radius = getDouble("radius", 1000.0);
    
    // Parse source type
    auto it = section.find("type");
    if (it != section.end()) {
        if (it->second == "MOGI") source.type = DeformationSourceType::MOGI;
        else if (it->second == "MCTIGUE") source.type = DeformationSourceType::MCTIGUE;
        else if (it->second == "DIKE") source.type = DeformationSourceType::RECTANGULAR_DIKE;
        else if (it->second == "SILL") source.type = DeformationSourceType::HORIZONTAL_SILL;
    }
    
    return true;
}

bool VolcanoConfigReader::parseGasSource(const std::map<std::string, std::string>& section,
                                          GasEmissionSource& source) {
    auto getDouble = [&section](const std::string& key, double default_val) {
        auto it = section.find(key);
        if (it == section.end()) return default_val;
        return std::stod(it->second);
    };
    
    source.x = getDouble("x", 0.0);
    source.y = getDouble("y", 0.0);
    source.z = getDouble("z", 0.0);
    source.emission_rate = getDouble("emission_rate", 1000.0);
    source.temperature = getDouble("temperature", 373.15);
    source.velocity = getDouble("velocity", 10.0);
    
    // Parse gas composition
    source.composition[VolcanicGas::H2O] = getDouble("X_H2O", 0.90);
    source.composition[VolcanicGas::CO2] = getDouble("X_CO2", 0.05);
    source.composition[VolcanicGas::SO2] = getDouble("X_SO2", 0.03);
    source.composition[VolcanicGas::H2S] = getDouble("X_H2S", 0.01);
    source.composition[VolcanicGas::HCl] = getDouble("X_HCl", 0.01);
    
    return true;
}

} // namespace FSRM
