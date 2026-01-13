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

#include "ExplosionImpactPhysics.hpp"

#include <algorithm>
#include <cmath>

namespace FSRM {

// =============================================================================
// NuclearSourceParameters
// =============================================================================

double NuclearSourceParameters::cavity_radius(double rock_density) const {
    // Simple empirical scaling for underground-contained shots:
    // Rc ≈ C * W^(1/3) * (rho_ref/rho)^(1/3)
    //
    // Typical coefficients for hard rock are on the order of 8–15 m/kt^(1/3).
    // We use a conservative mid-range value.
    const double rho_ref = 2650.0;
    const double C = 12.0;  // meters per kt^(1/3)

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
    // Phenomenological scaling for underground explosions (order-of-magnitude):
    //
    // M0 ~ 10^(13.3) * W^(0.75)  [N·m] is a common ballpark for mb-yield relations.
    // Keep it simple and monotone.
    const double W = std::max(0.0, yield_kt);
    return 2.0e13 * std::pow(W, 0.75);
}

double NuclearSourceParameters::body_wave_magnitude() const {
    // Simple mb scaling (educational): mb ≈ A + B log10(W)
    // Typical B ~ 0.75–0.85; choose 0.80 and A to give mb~4.5 at 1 kt.
    const double W = std::max(1e-12, yield_kt);
    return 4.5 + 0.80 * std::log10(W);
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

}  // namespace FSRM

