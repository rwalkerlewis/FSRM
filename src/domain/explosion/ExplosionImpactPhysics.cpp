/**
 * @file ExplosionImpactPhysics.cpp
 * @brief Implementations for explosion/impact source helpers
 *
 * This translation unit provides the missing link-time implementations for
 * the lightweight analytic/proxy models declared in ExplosionImpactPhysics.hpp.
 *
 * Note: These implementations are intentionally simplified/phenomenological.
 * They are primarily used for coupling (e.g., stress proxies) and example configs.
 */

#include "domain/explosion/ExplosionImpactPhysics.hpp"

#include <algorithm>
#include <cmath>

namespace FSRM {

// =============================================================================
// NuclearSourceParameters
// =============================================================================

double NuclearSourceParameters::cavityCoefficient(MediumType medium) {
    // Coefficients in m / kt^(1/3). See class docstring for citations.
    switch (medium) {
        case MediumType::GRANITE:  return 11.0;
        case MediumType::TUFF:     return 18.0;
        case MediumType::SALT:     return 16.0;
        case MediumType::ALLUVIUM: return 22.0;
        case MediumType::SHALE:    return 14.0;
        case MediumType::GENERIC:  return 12.0;
    }
    return 12.0;
}

double NuclearSourceParameters::cavity_radius(double rock_density) const {
    // Legacy single-coefficient overload: keep the historical C = 12
    // behaviour so existing call sites are unchanged.
    return cavity_radius(rock_density, MediumType::GENERIC);
}

double NuclearSourceParameters::cavity_radius(double rock_density,
                                              MediumType medium) const {
    // Empirical scaling for underground-contained shots:
    // Rc = C(medium) * W^(1/3) * (rho_ref/rho)^(1/3)
    //
    // The medium-aware coefficient distinguishes granite, tuff, salt,
    // alluvium, and shale -- single-medium coefficient hides ~3x cavity
    // size variation across realistic host rocks at fixed density.
    const double rho_ref = 2650.0;
    const double C = cavityCoefficient(medium);

    const double W = std::max(0.0, yield_kt);
    const double rho = (rock_density > 1.0) ? rock_density : rho_ref;
    return C * std::cbrt(W) * std::cbrt(rho_ref / rho);
}

double NuclearSourceParameters::crushed_zone_radius() const {
    // Educational scaling: crushed zone is a few cavity radii.
    return 3.0 * cavity_radius();
}

double NuclearSourceParameters::fractured_zone_radius() const {
    // Educational scaling: fractured zone extends ~5–15 cavity radii.
    return 10.0 * cavity_radius();
}

double NuclearSourceParameters::scalar_moment() const {
    // Legacy default: GENERIC medium (C = 12) so call sites that did
    // not pass a medium see the historical M0 value.
    return scalar_moment(MediumType::GENERIC);
}

double NuclearSourceParameters::scalar_moment(MediumType medium) const {
    // Cavity mechanics formula: M0 = 4*pi*rho*vp^2*Rc^3
    // where Rc is the medium-aware cavity radius from empirical scaling.
    const double W = std::max(1e-12, yield_kt);
    (void)W;
    double Rc = cavity_radius(2700.0, medium);
    double rho = 2700.0;
    double vp = 5500.0;
    return 4.0 * M_PI * rho * vp * vp * Rc * Rc * Rc;
}

double NuclearSourceParameters::body_wave_magnitude() const {
    // Murphy (1981) mb-yield relation: mb = 4.45 + 0.75 * log10(W)
    // For DPRK 2017 (250 kt): mb = 4.45 + 0.75*2.398 = 6.25 (observed: 6.3)
    const double W = std::max(1e-12, yield_kt);
    return 4.45 + 0.75 * std::log10(W);
}

double NuclearSourceParameters::surface_wave_magnitude() const {
    // Rough Ms scaling for larger yields; cap low yields sensibly.
    const double W = std::max(1e-12, yield_kt);
    return 3.8 + 1.00 * std::log10(W);
}

// =============================================================================
// SphericalCavitySource (analytic proxy)
// =============================================================================

SphericalCavitySource::SphericalCavitySource()
    : cavity_radius(1.0),
      density(2700.0),
      bulk_modulus(30e9),
      shear_modulus(20e9),
      initial_pressure(0.0),
      rise_time(0.05),
      p_velocity(5500.0),
      s_velocity(3200.0) {}

void SphericalCavitySource::setCavityRadius(double R) {
    cavity_radius = std::max(1e-6, R);
}

void SphericalCavitySource::setMediumProperties(double rho, double K, double G) {
    density = std::max(1.0, rho);
    bulk_modulus = std::max(1e6, K);
    shear_modulus = std::max(1e6, G);

    // Derived wave speeds
    p_velocity = std::sqrt((bulk_modulus + 4.0 / 3.0 * shear_modulus) / density);
    s_velocity = std::sqrt(shear_modulus / density);
}

void SphericalCavitySource::setOverpressure(double P0) {
    initial_pressure = P0;
}

void SphericalCavitySource::setRiseTime(double tau) {
    rise_time = std::max(1e-9, tau);
}

static inline double smoothRise(double t, double tau) {
    if (t <= 0.0) return 0.0;
    // Exponential rise to 1
    return 1.0 - std::exp(-t / std::max(1e-12, tau));
}

double SphericalCavitySource::displacement(double r, double t) const {
    // Quasi-static spherical cavity displacement with a simple retarded-time rise.
    const double rr = std::max(r, cavity_radius);
    const double tret = t - (rr - cavity_radius) / std::max(1e-6, p_velocity);
    const double P = initial_pressure * smoothRise(tret, rise_time);

    // Poisson's ratio from K,G
    const double nu = (3.0 * bulk_modulus - 2.0 * shear_modulus) /
                      (2.0 * (3.0 * bulk_modulus + shear_modulus));
    const double fac = (1.0 - 2.0 * nu) / std::max(1e-12, (1.0 - nu));

    // u(r) ~ (P a^3 / (4 G r^2)) * (1-2nu)/(1-nu)
    return (P * std::pow(cavity_radius, 3) / (4.0 * shear_modulus * rr * rr)) * fac;
}

double SphericalCavitySource::velocity(double r, double t) const {
    // Finite-difference derivative of displacement (robust for this proxy).
    const double dt = std::max(1e-6, 0.01 * rise_time);
    const double u1 = displacement(r, t);
    const double u0 = displacement(r, t - dt);
    return (u1 - u0) / dt;
}

void SphericalCavitySource::stress(double r, double t, double& sigma_rr, double& sigma_tt) const {
    // Thick-walled sphere / elastic cavity proxy:
    // σ_rr(r) = -P(t_r) * (a^3 / r^3)
    // σ_tt(r) = +0.5 P(t_r) * (a^3 / r^3)
    const double rr = std::max(r, cavity_radius);
    const double tret = t - (rr - cavity_radius) / std::max(1e-6, p_velocity);
    const double P = initial_pressure * smoothRise(tret, rise_time);
    const double s = std::pow(cavity_radius / rr, 3);

    sigma_rr = -P * s;
    sigma_tt = 0.5 * P * s;
}

double SphericalCavitySource::equivalentMoment() const {
    // Isotropic equivalent moment proxy (order-of-magnitude):
    // M ~ (4/3)π a^3 P
    return (4.0 / 3.0) * M_PI * std::pow(cavity_radius, 3) * initial_pressure;
}

// =============================================================================
// MuellerMurphySource (reduced displacement potential proxy)
// =============================================================================

MuellerMurphySource::MuellerMurphySource()
    : medium_type(NuclearSourceParameters::MediumType::GENERIC),
      density(2700.0),
      p_velocity(5000.0),
      s_velocity(3000.0),
      corner_frequency(1.0),
      scalar_moment(0.0),
      overshoot(1.0),
      rise_time(0.05) {}

void MuellerMurphySource::setParameters(const NuclearSourceParameters& params) {
    source_params = params;
    computeDerivedQuantities();
}

void MuellerMurphySource::setMediumProperties(double rho, double vp, double vs) {
    density = std::max(1.0, rho);
    p_velocity = std::max(1.0, vp);
    s_velocity = std::max(1.0, vs);
    computeDerivedQuantities();
}

void MuellerMurphySource::setMedium(
        NuclearSourceParameters::MediumType medium) {
    medium_type = medium;
    computeDerivedQuantities();
}

void MuellerMurphySource::computeDerivedQuantities() {
    // Use the medium-aware scalar_moment overload so a setMedium call
    // correctly rescales M0 via the cavity coefficient table. The
    // default medium_type is GENERIC, which preserves the historical
    // M0 value for callers that did not call setMedium.
    scalar_moment = source_params.scalar_moment(medium_type);
    // Patton (1988) corner frequency: fc = 3.0 * W^(-1/3) Hz for hard rock
    // with density correction for non-granite media (baseline rho = 2650 kg/m^3)
    const double W = std::max(1e-6, source_params.yield_kt);
    corner_frequency = 3.0 * std::pow(W, -1.0 / 3.0) *
                       std::pow(density / 2650.0, -1.0 / 3.0);
    rise_time = 0.55 / corner_frequency;
    overshoot = 1.0;
}

double MuellerMurphySource::momentRate(double t) const {
    if (t < 0.0 || scalar_moment <= 0.0) {
        return 0.0;
    }
    const double tau = std::max(1e-12, rise_time);
    return (scalar_moment / tau) * std::exp(-t / tau);
}

std::complex<double> MuellerMurphySource::rdp(double omega) const {
    const double w = std::max(1e-12, std::abs(omega));
    const double omega_c = 2.0 * M_PI * std::max(1e-12, corner_frequency);
    const double den = 1.0 + (w / omega_c) * (w / omega_c);
    return std::complex<double>(scalar_moment / den, 0.0);
}

double MuellerMurphySource::getCornerFrequency() const {
    return corner_frequency;
}

double MuellerMurphySource::getOvershoot() const {
    return overshoot;
}

void MuellerMurphySource::getMomentTensor(double t, double* M) const {
    const double mr = momentRate(t);
    const double scale = mr / std::max(1e-30, scalar_moment) * (scalar_moment / 3.0);
    M[0] = M[1] = M[2] = scale;
    M[3] = M[4] = M[5] = 0.0;
}

}  // namespace FSRM

