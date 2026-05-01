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
      damping(1.0),
      overshoot_zero_factor(1.0),
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

void MuellerMurphySource::setOvershoot(double B) {
    overshoot = std::max(0.0, B);
    // Note: do not call computeDerivedQuantities here because that
    // method now resets nothing related to overshoot; it would
    // recompute scalar_moment and corner_frequency, which are
    // independent of B / zeta / k_B.
}

void MuellerMurphySource::setDamping(double zeta) {
    // Clamp into the physically meaningful damping range (0, 1].
    if (zeta <= 0.0) zeta = 1e-3;
    if (zeta > 1.0) zeta = 1.0;
    damping = zeta;
}

void MuellerMurphySource::setOvershootZeroFactor(double k_B) {
    overshoot_zero_factor = std::max(1e-3, k_B);
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
    // overshoot, damping, and overshoot_zero_factor preserved across
    // rebuilds so user-set RDP shape parameters survive a setMedium
    // or setMediumProperties call.
}

// Mueller-Murphy (1971) moment-rate function (Fourier pair of ::rdp).
//
// The frequency-domain RDP is the canonical Mueller-Murphy damped
// second-order resonator
//
//   Psi(omega) = M0 * (1 + i (B - 1) omega / omega_zero)
//                / (1 - (omega/omega_p)^2 + 2 i zeta omega / omega_p)
//
// The time-domain moment rate Mdot(t) is its inverse Fourier transform.
// For the default elastic-overshoot-free shape (B = 1) the (B - 1)
// numerator zero is gated off and the time-domain function is the
// impulse response of a damped harmonic oscillator scaled by omega_p^2:
//
//   zeta < 1 (underdamped):
//     Mdot(t) = (M0 omega_p^2 / omega_d) * exp(-zeta omega_p t) * sin(omega_d t)
//     omega_d = omega_p sqrt(1 - zeta^2)
//
//   zeta = 1 (critically damped, the default):
//     Mdot(t) = M0 omega_p^2 * t * exp(-omega_p t)
//     non-negative for all t >= 0; integrates to M0; peak at t = 1/omega_p
//     with value M0 omega_p / e
//
//   zeta > 1 (overdamped):
//     Mdot(t) = (M0 omega_p^2 / (2 omega_h)) *
//               (exp(-(zeta omega_p - omega_h) t) - exp(-(zeta omega_p + omega_h) t))
//     omega_h = omega_p sqrt(zeta^2 - 1)
//
// For the elastic-overshoot case (B != 1) the additional numerator
// factor (B - 1) i omega / omega_zero in the frequency domain is the
// Laplace transform of a time derivative scaled by (B - 1) / omega_zero.
// The overshoot moment rate is therefore
//
//   Mdot_overshoot(t) = Mdot_base(t) + (B - 1) / omega_zero * dMdot_base/dt
//
// where dMdot_base/dt is computed analytically per damping case.
//
// Default damping zeta = 1.0: critically damped is non-negative for all
// t and integrates to M0 analytically -- no FD-derivative tail, no
// Gibbs ripple from the ::rdp inversion. Underdamped values can be
// chosen via setDamping(zeta < 1) when a user wants to model the
// elastic overshoot peak (Stevens & Day 1985 figure 4) explicitly.
//
// References: Mueller & Murphy (1971) BSSA 61(6); Murphy (1981) DARPA
// report; Stevens & Day (1985) JGR 90, eq. 4 and 6.
double MuellerMurphySource::momentRate(double t) const {
    if (t < 0.0 || scalar_moment <= 0.0) {
        return 0.0;
    }
    const double omega_p = 2.0 * M_PI * std::max(1e-12, corner_frequency);
    const double zeta = damping;
    const double M0 = scalar_moment;
    const double op2 = omega_p * omega_p;

    auto base_and_deriv = [&](double t_eval, double& h, double& hp) {
        // h(t) = base impulse response of M0 * H(s) where
        //   H(s) = omega_p^2 / (s^2 + 2 zeta omega_p s + omega_p^2)
        // hp(t) = dh/dt evaluated analytically per damping case.
        if (zeta < 1.0 - 1e-9) {
            const double zw = zeta * omega_p;
            const double od = omega_p * std::sqrt(std::max(0.0, 1.0 - zeta * zeta));
            const double e_d = std::exp(-zw * t_eval);
            const double s_d = std::sin(od * t_eval);
            const double c_d = std::cos(od * t_eval);
            // h(t) = (M0 * omega_p^2 / omega_d) * exp(-zeta omega_p t) sin(omega_d t)
            h = (M0 * op2 / std::max(1e-30, od)) * e_d * s_d;
            // hp(t) = (M0 * omega_p^2 / omega_d) * exp(-zeta omega_p t)
            //         * (-zeta omega_p sin(omega_d t) + omega_d cos(omega_d t))
            hp = (M0 * op2 / std::max(1e-30, od)) * e_d *
                 (-zw * s_d + od * c_d);
        } else if (zeta > 1.0 + 1e-9) {
            const double zw = zeta * omega_p;
            const double oh = omega_p * std::sqrt(zeta * zeta - 1.0);
            const double a1 = zw - oh;
            const double a2 = zw + oh;
            const double e1 = std::exp(-a1 * t_eval);
            const double e2 = std::exp(-a2 * t_eval);
            // h(t) = (M0 * omega_p^2 / (2 omega_h)) * (exp(-a1 t) - exp(-a2 t))
            h = (M0 * op2 / (2.0 * std::max(1e-30, oh))) * (e1 - e2);
            // hp(t) = (M0 * omega_p^2 / (2 omega_h)) * (-a1 exp(-a1 t) + a2 exp(-a2 t))
            hp = (M0 * op2 / (2.0 * std::max(1e-30, oh))) *
                 (-a1 * e1 + a2 * e2);
        } else {
            // Critically damped: zeta == 1 (within +/- 1e-9 tolerance).
            const double e_c = std::exp(-omega_p * t_eval);
            // h(t) = M0 * omega_p^2 * t * exp(-omega_p t)
            h = M0 * op2 * t_eval * e_c;
            // hp(t) = M0 * omega_p^2 * (1 - omega_p t) * exp(-omega_p t)
            hp = M0 * op2 * (1.0 - omega_p * t_eval) * e_c;
        }
    };

    double h = 0.0;
    double hp = 0.0;
    base_and_deriv(t, h, hp);

    // Default (B = 1) gating: the (B - 1) numerator zero is suppressed
    // and the time-domain function reduces to the pure damped impulse
    // response. setOvershoot(B > 1) activates the elastic-overshoot
    // peak; the corresponding time-domain term is (B - 1) / omega_zero
    // times the analytic derivative of the base impulse response.
    if (std::abs(overshoot - 1.0) <= 1e-12) {
        return h;
    }
    const double omega_zero =
        omega_p * std::max(1e-9, overshoot_zero_factor);
    return h + (overshoot - 1.0) / omega_zero * hp;
}

// Mueller-Murphy (1971) RDP in the frequency domain.
//
//   Psi(omega) = M0 * (1 + i * (B - 1) * omega / omega_zero)
//                / (1 - (omega/omega_p)^2 + 2 i zeta omega / omega_p)
//
// where M0 is the steady-state seismic moment (low-frequency plateau),
// omega_p = 2 pi f_c is the corner frequency in radians, omega_zero =
// omega_p * k_B places the numerator zero, B is the elastic overshoot
// factor and zeta is the damping factor.
//
// Defaults (B = 1, zeta = 0.7, k_B = 1) suppress the numerator zero
// via the (B - 1) gating so the spectrum is the pure damped second-
// order low-pass with omega^-2 high-frequency decay -- this preserves
// the existing RDPSpectralShape and RDPLowFrequencyPlateauEqualsM0
// tests. For tamped granite or studies that resolve the elastic-
// overshoot peak, set B > 1 (and optionally tighten zeta) via
// setOvershoot / setDamping to activate the M&M numerator zero.
//
// Reference: Mueller & Murphy (1971) BSSA 61(6) figs. 3, 6;
// Stevens & Day (1985) JGR 90, eq. 4.
std::complex<double> MuellerMurphySource::rdp(double omega) const {
    const double w = std::max(1e-12, std::abs(omega));
    const double omega_p = 2.0 * M_PI * std::max(1e-12, corner_frequency);
    const double zeta = damping;
    const double k_B = std::max(1e-9, overshoot_zero_factor);
    const double omega_zero = omega_p * k_B;

    const double x = w / omega_p;
    using cd = std::complex<double>;
    // (B - 1) gating: B = 1 default leaves numerator = 1 so the
    // spectrum is the pure damped second-order low-pass. B > 1
    // activates the overshoot zero with strength (B - 1).
    const cd numerator(1.0, (overshoot - 1.0) * w / omega_zero);
    const cd denom(1.0 - x * x, 2.0 * zeta * x);
    return scalar_moment * (numerator / denom);
}

double MuellerMurphySource::getCornerFrequency() const {
    return corner_frequency;
}

double MuellerMurphySource::getOvershoot() const {
    return overshoot;
}

double MuellerMurphySource::getDamping() const {
    return damping;
}

void MuellerMurphySource::getMomentTensor(double t, double* M) const {
    const double mr = momentRate(t);
    const double scale = mr / std::max(1e-30, scalar_moment) * (scalar_moment / 3.0);
    M[0] = M[1] = M[2] = scale;
    M[3] = M[4] = M[5] = 0.0;
}

}  // namespace FSRM

