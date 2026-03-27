/**
 * @file CaprockIntegrity.cpp
 * @brief Caprock integrity assessment for CCS (Mohr–Coulomb, ΔCFS, permeability).
 */

#include "domain/geomechanics/CaprockIntegrity.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace FSRM {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kDefaultFaultAreaM2 = 1.0e6;

inline double degreesToRadians(double deg) { return deg * kDegToRad; }

// Symmetric 3×3 stress from Voigt {xx, yy, zz, xy, xz, yz}.
inline void stressTensorFromVoigt(const std::array<double, 6>& s, double S[3][3]) {
    S[0][0] = s[0];
    S[0][1] = s[3];
    S[0][2] = s[4];
    S[1][0] = s[3];
    S[1][1] = s[1];
    S[1][2] = s[5];
    S[2][0] = s[4];
    S[2][1] = s[5];
    S[2][2] = s[2];
}

// Eigenvalues of real symmetric 3×3: characteristic cubic solved by depressed cubic + Cardano /
// Viète trigonometric form (three real roots).
void symmetricEigenvaluesCardano(const double S[3][3], double lam[3]) {
    const double s11 = S[0][0];
    const double s12 = S[0][1];
    const double s13 = S[0][2];
    const double s22 = S[1][1];
    const double s23 = S[1][2];
    const double s33 = S[2][2];

    const double I1 = s11 + s22 + s33;
    const double I2 = (s11 * s22 - s12 * s12) + (s22 * s33 - s23 * s23) + (s11 * s33 - s13 * s13);
    const double I3 = s11 * (s22 * s33 - s23 * s23) - s12 * (s12 * s33 - s13 * s23)
                    + s13 * (s12 * s23 - s13 * s22);

    const double a = -I1;
    const double b = I2;
    const double c = -I3;

    const double p = b - a * a / 3.0;
    const double q = c + (2.0 * a * a * a - 9.0 * a * b) / 27.0;
    const double shift = -a / 3.0;

    if (std::fabs(p) < 1.0e-14) {
        const double y = std::cbrt(-q);
        lam[0] = lam[1] = lam[2] = y + shift;
        return;
    }

    if (p < 0.0) {
        const double m = 2.0 * std::sqrt(-p / 3.0);
        const double arg = std::clamp((3.0 * q) / (2.0 * p) * std::sqrt(-3.0 / p), -1.0, 1.0);
        const double theta = std::acos(arg);
        lam[0] = m * std::cos(theta / 3.0) + shift;
        lam[1] = m * std::cos((theta + 2.0 * kPi) / 3.0) + shift;
        lam[2] = m * std::cos((theta + 4.0 * kPi) / 3.0) + shift;
    } else {
        const double disc = 0.25 * q * q + p * p * p / 27.0;
        const double sqrt_disc = std::sqrt(std::max(0.0, disc));
        const double u = std::cbrt(-0.5 * q + sqrt_disc);
        const double v = std::cbrt(-0.5 * q - sqrt_disc);
        lam[0] = lam[1] = lam[2] = u + v + shift;
    }
}

void normalize3(double v[3]) {
    const double nrm = std::hypot(v[0], std::hypot(v[1], v[2]));
    if (nrm > 1.0e-300) {
        v[0] /= nrm;
        v[1] /= nrm;
        v[2] /= nrm;
    } else {
        v[0] = 0.0;
        v[1] = 0.0;
        v[2] = 1.0;
    }
}

// Unit eigenvector for eigenvalue lambda (symmetric S).
void eigenvectorSymmetric(const double S[3][3], double lambda, double e[3]) {
    const double A[3][3] = {{S[0][0] - lambda, S[0][1], S[0][2]},
                          {S[0][1], S[1][1] - lambda, S[1][2]},
                          {S[0][2], S[1][2], S[2][2] - lambda}};

    const double r0[3] = {A[0][0], A[0][1], A[0][2]};
    const double r1[3] = {A[1][0], A[1][1], A[1][2]};
    const double r2[3] = {A[2][0], A[2][1], A[2][2]};

    double c0[3] = {r0[1] * r1[2] - r0[2] * r1[1], r0[2] * r1[0] - r0[0] * r1[2],
                    r0[0] * r1[1] - r0[1] * r1[0]};
    double c1[3] = {r0[1] * r2[2] - r0[2] * r2[1], r0[2] * r2[0] - r0[0] * r2[2],
                    r0[0] * r2[1] - r0[1] * r2[0]};
    double c2[3] = {r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2],
                    r1[0] * r2[1] - r1[1] * r2[0]};

    double n0 = std::hypot(c0[0], std::hypot(c0[1], c0[2]));
    double n1 = std::hypot(c1[0], std::hypot(c1[1], c1[2]));
    double n2 = std::hypot(c2[0], std::hypot(c2[1], c2[2]));

    if (n0 >= n1 && n0 >= n2) {
        e[0] = c0[0] / n0;
        e[1] = c0[1] / n0;
        e[2] = c0[2] / n0;
    } else if (n1 >= n2) {
        e[0] = c1[0] / n1;
        e[1] = c1[1] / n1;
        e[2] = c1[2] / n1;
    } else if (n2 > 1.0e-300) {
        e[0] = c2[0] / n2;
        e[1] = c2[1] / n2;
        e[2] = c2[2] / n2;
    } else {
        e[0] = 1.0;
        e[1] = 0.0;
        e[2] = 0.0;
    }
}

inline double dot3(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Normal and shear magnitude on the plane with unit normal n using principal stresses and
// orthonormal eigenvectors (columns of R in row form: e_i[j] = component j of eigenvector i).
void normalShearFromPrincipal(const std::array<double, 3>& principal,
                              const double evec[3][3], const double n[3], double& sigma_n,
                              double& tau) {
    double n_p[3] = {dot3(n, evec[0]), dot3(n, evec[1]), dot3(n, evec[2])};

    sigma_n = principal[0] * n_p[0] * n_p[0] + principal[1] * n_p[1] * n_p[1]
            + principal[2] * n_p[2] * n_p[2];

    double t[3] = {principal[0] * n_p[0] * evec[0][0] + principal[1] * n_p[1] * evec[1][0]
                       + principal[2] * n_p[2] * evec[2][0],
                   principal[0] * n_p[0] * evec[0][1] + principal[1] * n_p[1] * evec[1][1]
                       + principal[2] * n_p[2] * evec[2][1],
                   principal[0] * n_p[0] * evec[0][2] + principal[1] * n_p[1] * evec[1][2]
                       + principal[2] * n_p[2] * evec[2][2]};

    const double tx = t[0] - sigma_n * n[0];
    const double ty = t[1] - sigma_n * n[1];
    const double tz = t[2] - sigma_n * n[2];
    tau = std::hypot(tx, std::hypot(ty, tz));
}

} // namespace

double CaprockIntegrity::mohrCoulombFOS(double sigma_n, double tau, double cohesion,
                                        double friction_angle) {
    if (tau <= 1.0e-300) {
        // No resolved shear: unlimited FOS; use a large finite value so callers/tests
        // using std::isfinite remain well-defined.
        return 1.0e200;
    }
    const double phi = degreesToRadians(friction_angle);
    const double strength = cohesion + sigma_n * std::tan(phi);
    return strength / tau;
}

double CaprockIntegrity::caprockPermeability(double k0, double sigma_eff, double sigma_eff_ref,
                                             double gamma, double enhancement_factor, bool failed) {
    if (failed) {
        return k0 * enhancement_factor;
    }
    return k0 * std::exp(gamma * (sigma_eff - sigma_eff_ref));
}

double CaprockIntegrity::deltaCFS(double delta_tau, double delta_sigma_n, double friction,
                                    double delta_pore_pressure) {
    return delta_tau - friction * (delta_sigma_n - delta_pore_pressure);
}

double CaprockIntegrity::estimateMagnitude(double fault_area) {
    if (fault_area <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (std::log10(fault_area) - 4.07) / 0.98;
}

CaprockIntegrity::IntegrityResult CaprockIntegrity::evaluate(
    const std::array<double, 6>& stress, double pore_pressure, double cohesion,
    double friction_angle, const std::array<double, 3>& fault_normal) const {
    IntegrityResult r{};
    double S[3][3];
    stressTensorFromVoigt(stress, S);

    double lam[3];
    symmetricEigenvaluesCardano(S, lam);

    double evec[3][3];
    for (int k = 0; k < 3; ++k) {
        eigenvectorSymmetric(S, lam[k], evec[k]);
        normalize3(evec[k]);
    }

    int ord[3] = {0, 1, 2};
    std::sort(ord, ord + 3, [&](int i, int j) { return lam[i] > lam[j]; });
    std::array<double, 3> evals{};
    double e_sorted[3][3];
    for (int i = 0; i < 3; ++i) {
        evals[i] = lam[ord[i]];
        std::memcpy(e_sorted[i], evec[ord[i]], 3 * sizeof(double));
    }

    double n[3] = {fault_normal[0], fault_normal[1], fault_normal[2]};
    normalize3(n);

    double sigma_n = 0.0;
    double tau = 0.0;
    normalShearFromPrincipal(evals, e_sorted, n, sigma_n, tau);

    const double sigma_n_eff = sigma_n - pore_pressure;
    r.mohr_coulomb_safety_factor = mohrCoulombFOS(sigma_n_eff, tau, cohesion, friction_angle);

    const double mu = std::tan(degreesToRadians(friction_angle));
    r.max_delta_CFS = deltaCFS(tau, sigma_n, mu, pore_pressure);

    r.failure_predicted = (r.mohr_coulomb_safety_factor < 1.0);

    const double phi = degreesToRadians(friction_angle);
    const double tan_phi = std::tan(phi);
    if (tan_phi > 1.0e-14) {
        r.critical_pressure = sigma_n - (tau - cohesion) / tan_phi;
    } else {
        r.critical_pressure = std::numeric_limits<double>::infinity();
    }

    r.estimated_magnitude = estimateMagnitude(kDefaultFaultAreaM2);

    return r;
}

void CaprockIntegrity::setCaprockProperties(double k0, double gamma, double enhancement) {
    k0_ = k0;
    gamma_ = gamma;
    enhancement_ = enhancement;
}

void CaprockIntegrity::setMohrCoulombParams(double cohesion, double friction_angle) {
    cohesion_ = cohesion;
    friction_angle_ = friction_angle;
}

} // namespace FSRM
