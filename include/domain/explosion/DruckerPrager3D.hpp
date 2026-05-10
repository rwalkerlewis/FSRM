/**
 * @file DruckerPrager3D.hpp
 * @brief Pass-13b (axis-1b physics) explicit 3D Drucker-Prager radial
 *        return for the Source3DBall solver.
 *
 * The 1D radial Lagrangian (pass-7..10) uses the spherical-symmetry
 * constraint s_tt = -s_rr / 2 to reduce the deviatoric tensor to a
 * single scalar. Pass-13b replaces that scalar return with a 6-component
 * Voigt tensor return on the unstructured 3D source-ball mesh.
 *
 * Algorithm (Simo & Hughes 1998, "Computational Inelasticity" sec 3.6,
 * boxed algorithm 3.5.1 specialised to the linear Drucker-Prager
 * surface):
 *
 *   1. Trial elastic predictor:
 *        sigma_trial = sigma + C : delta_eps
 *      with C the standard isotropic elastic stiffness defined by
 *      bulk modulus K and shear modulus G.
 *   2. Decompose: I1 = trace(sigma_trial),  s = sigma_trial - I1/3 * I.
 *   3. Yield function:
 *        f(sigma) = sqrt(J_2(s)) - alpha * I1 / 3 - k.
 *      With sqrt(J_2(s)) = sqrt(0.5 * s : s).
 *   4. If f(sigma_trial) <= 0: accept the elastic predictor unchanged.
 *   5. Else: explicit (non-iterative) radial return. The standard
 *      first-order DP return uses a single plastic multiplier
 *        delta_lambda = f_trial / (G + K * alpha^2)
 *      and projects
 *        sigma_new = sigma_trial - delta_lambda * (
 *                       G * s_trial / sqrt(J2_trial) +
 *                       K * alpha * I).
 *      The increment in plastic strain is
 *        delta_eps_p = delta_lambda * (
 *                       0.5 * s_trial / sqrt(J2_trial) +
 *                       (alpha / 3) * I).
 *      Equivalent plastic strain
 *        delta_eps_p_eq = delta_lambda * sqrt(2/3 + 2 * alpha^2 / 3).
 *
 * Pass-13b uses the explicit form (no inner Newton on delta_lambda) per
 * the design-doc rationale (CI tractability). Pass-15+ may switch to
 * an implicit return if accuracy demands.
 *
 * Spherical-symmetry reduction (used by the SphericalSymmetryReducesTo1D
 * unit gate): with sigma_trial isotropic plus a radial deviator
 * s_rr = 2*A, s_tt = s_phiphi = -A, the resulting sqrt(J_2) reduces to
 * sqrt(3) * |A|, and the pass-7..10 1D scalar return s_rr *= Y/sigma_eq
 * (sigma_eq = 1.5 * |s_rr| = 3 * |A|) gives the same projected stress
 * tensor under pure-Mises yield (alpha=0). The DP3D reproduces this
 * limit exactly.
 *
 * References:
 *   - Simo, J. C. and Hughes, T. J. R. (1998). Computational
 *     Inelasticity. Springer. Section 3.6 ("Pressure-dependent
 *     plasticity", Drucker-Prager radial return).
 *   - Drucker, D. C. and Prager, W. (1952). Soil mechanics and plastic
 *     analysis or limit design. Quarterly of Applied Math 10, pp 157-165.
 *   - Hoek, E. and Brown, E. T. (1980). Underground excavations in rock
 *     (medium-specific yield-surface parameters).
 */

#ifndef NEAR_FIELD_DRUCKER_PRAGER_3D_HPP
#define NEAR_FIELD_DRUCKER_PRAGER_3D_HPP

#include <array>
#include <cmath>
#include <string>

namespace FSRM {

/// Pass-13b Drucker-Prager parameter set.
/// Linear yield surface: sqrt(J2) = alpha_dp * (-I1/3) + k_dp.
/// Equivalently sqrt(J2) - alpha * I1 / 3 - k = 0 (sign convention:
/// compression negative for I1, so the alpha term is +alpha * P with
/// P = -I1/3 the mean compressive stress).
///
/// Per-medium presets are derived from the Hoek-Brown / Mohr-Coulomb
/// envelope for representative crustal rocks. The values are coarse
/// (factor-2 correct) and intended for the pass-13b validation gates;
/// pass-15+ axis-4 will refit against ANEOS.
struct DruckerPrager3DParameters
{
    /// Pressure-sensitivity coefficient. Dimensionless. Pure von Mises
    /// has alpha = 0 (used for unit gates).
    double alpha_dp = 0.3;

    /// Cohesive strength at zero confining pressure [Pa].
    double k_dp = 50.0e6;

    /// Diagnostic medium label.
    std::string medium_label = "default";
};

namespace DruckerPrager3DSets {

inline DruckerPrager3DParameters granite()
{
    DruckerPrager3DParameters p;
    p.alpha_dp = 0.30;
    p.k_dp = 70.0e6;
    p.medium_label = "GRANITE";
    return p;
}

inline DruckerPrager3DParameters tuff()
{
    DruckerPrager3DParameters p;
    p.alpha_dp = 0.22;
    p.k_dp = 30.0e6;
    p.medium_label = "TUFF";
    return p;
}

inline DruckerPrager3DParameters salt()
{
    DruckerPrager3DParameters p;
    p.alpha_dp = 0.10;
    p.k_dp = 20.0e6;
    p.medium_label = "SALT";
    return p;
}

inline DruckerPrager3DParameters alluvium()
{
    DruckerPrager3DParameters p;
    p.alpha_dp = 0.18;
    p.k_dp = 5.0e6;
    p.medium_label = "ALLUVIUM";
    return p;
}

/// Lookup by case-insensitive medium name. Falls back to granite when
/// the name is unrecognised. Used by the host's medium dispatch.
inline DruckerPrager3DParameters byName(const std::string& upper_name)
{
    if (upper_name == "GRANITE") return granite();
    if (upper_name == "TUFF")    return tuff();
    if (upper_name == "SALT")    return salt();
    if (upper_name == "ALLUVIUM") return alluvium();
    return granite();
}

}  // namespace DruckerPrager3DSets

/// Result of one radial-return projection on a single cell.
struct DruckerPrager3DReturnResult
{
    /// True if the trial state lay outside the yield surface and was
    /// projected. False for a pure elastic step (sigma_new = sigma_trial).
    bool yielded = false;
    /// Plastic multiplier from the linearised return [strain units].
    double delta_lambda = 0.0;
    /// Equivalent plastic strain increment from this step.
    double delta_eps_p_eq = 0.0;
};

/// Voigt convention used throughout the 3D source-ball solver:
///   index 0 = xx, 1 = yy, 2 = zz, 3 = xy, 4 = xz, 5 = yz.
/// Strain Voigt entries 3..5 use engineering shear strain
/// (gamma_ij = 2 * eps_ij). The radial-return below works in stress
/// Voigt and uses the standard isotropic elastic mapping
///   sigma_trial = K * trace(eps) * I + 2 G * dev(eps).
namespace Voigt6 {

inline double trace(const std::array<double, 6>& s) { return s[0] + s[1] + s[2]; }

inline std::array<double, 6> deviator(const std::array<double, 6>& s)
{
    const double p_iso = trace(s) / 3.0;
    return {s[0] - p_iso, s[1] - p_iso, s[2] - p_iso, s[3], s[4], s[5]};
}

inline double sqrtJ2(const std::array<double, 6>& dev)
{
    const double sq = 0.5 * (dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2])
                      + dev[3] * dev[3] + dev[4] * dev[4] + dev[5] * dev[5];
    return std::sqrt(sq > 0.0 ? sq : 0.0);
}

}  // namespace Voigt6

/// Apply the standard isotropic elastic stiffness to a small-strain
/// increment. eps Voigt indices 3..5 are engineering shear strains.
inline std::array<double, 6> applyIsotropicStiffness(double K, double G,
                                                     const std::array<double, 6>& eps)
{
    const double tr = eps[0] + eps[1] + eps[2];
    const double lam = K - (2.0 / 3.0) * G;
    std::array<double, 6> sigma;
    sigma[0] = lam * tr + 2.0 * G * eps[0];
    sigma[1] = lam * tr + 2.0 * G * eps[1];
    sigma[2] = lam * tr + 2.0 * G * eps[2];
    sigma[3] = G * eps[3];  // engineering shear: factor-2 absorbed already
    sigma[4] = G * eps[4];
    sigma[5] = G * eps[5];
    return sigma;
}

/// Pass-13b Drucker-Prager 3D explicit radial return.
///
/// Inputs:
///   sigma_in: 6-Voigt stress at start of the step [Pa].
///   delta_eps: 6-Voigt small-strain increment over dt (engineering
///              shear in 3..5).
///   K, G: bulk and shear modulus [Pa].
///   params: yield-surface parameters.
///
/// Outputs (modify in place):
///   sigma_out: stress after the radial return.
///   delta_eps_p_out: 6-Voigt plastic-strain increment (engineering
///                    shear in 3..5).
///
/// Returns the result struct describing whether the step yielded.
inline DruckerPrager3DReturnResult druckerPrager3DRadialReturn(
    const std::array<double, 6>& sigma_in,
    const std::array<double, 6>& delta_eps,
    double K, double G,
    const DruckerPrager3DParameters& params,
    std::array<double, 6>& sigma_out,
    std::array<double, 6>& delta_eps_p_out)
{
    DruckerPrager3DReturnResult res;

    const std::array<double, 6> dsigma_trial = applyIsotropicStiffness(K, G, delta_eps);
    std::array<double, 6> sigma_trial;
    for (int i = 0; i < 6; ++i) sigma_trial[i] = sigma_in[i] + dsigma_trial[i];

    const double I1_trial = Voigt6::trace(sigma_trial);
    const std::array<double, 6> s_trial = Voigt6::deviator(sigma_trial);
    const double sqrtJ2_trial = Voigt6::sqrtJ2(s_trial);

    const double f_trial = sqrtJ2_trial
                           - params.alpha_dp * I1_trial / 3.0
                           - params.k_dp;

    if (f_trial <= 0.0) {
        sigma_out = sigma_trial;
        for (int i = 0; i < 6; ++i) delta_eps_p_out[i] = 0.0;
        res.yielded = false;
        return res;
    }

    // Explicit linearised radial return (Simo & Hughes 1998 sec 3.6).
    // Plastic multiplier delta_lambda from the consistency condition
    // f(sigma_new) = 0 linearised about the trial state:
    //   delta_lambda = f_trial / (G + K * alpha^2).
    const double denom = G + K * params.alpha_dp * params.alpha_dp;
    const double dlambda = f_trial / (denom > 0.0 ? denom : 1.0);
    res.delta_lambda = dlambda;

    // Direction of the deviatoric return (unit tensor in deviatoric
    // space).
    const double inv_J2 = (sqrtJ2_trial > 0.0) ? 1.0 / sqrtJ2_trial : 0.0;
    std::array<double, 6> n_dev;
    n_dev[0] = 0.5 * s_trial[0] * inv_J2;
    n_dev[1] = 0.5 * s_trial[1] * inv_J2;
    n_dev[2] = 0.5 * s_trial[2] * inv_J2;
    n_dev[3] = s_trial[3] * inv_J2;
    n_dev[4] = s_trial[4] * inv_J2;
    n_dev[5] = s_trial[5] * inv_J2;

    // Stress projection.
    sigma_out[0] = sigma_trial[0] - dlambda * (2.0 * G * n_dev[0] + K * params.alpha_dp);
    sigma_out[1] = sigma_trial[1] - dlambda * (2.0 * G * n_dev[1] + K * params.alpha_dp);
    sigma_out[2] = sigma_trial[2] - dlambda * (2.0 * G * n_dev[2] + K * params.alpha_dp);
    sigma_out[3] = sigma_trial[3] - dlambda * (2.0 * G * n_dev[3]);
    sigma_out[4] = sigma_trial[4] - dlambda * (2.0 * G * n_dev[4]);
    sigma_out[5] = sigma_trial[5] - dlambda * (2.0 * G * n_dev[5]);

    // Plastic strain increment, Voigt with engineering shear in 3..5.
    delta_eps_p_out[0] = dlambda * (n_dev[0] + params.alpha_dp / 3.0);
    delta_eps_p_out[1] = dlambda * (n_dev[1] + params.alpha_dp / 3.0);
    delta_eps_p_out[2] = dlambda * (n_dev[2] + params.alpha_dp / 3.0);
    delta_eps_p_out[3] = dlambda * (2.0 * n_dev[3]);
    delta_eps_p_out[4] = dlambda * (2.0 * n_dev[4]);
    delta_eps_p_out[5] = dlambda * (2.0 * n_dev[5]);

    // Equivalent plastic strain. For pure deviatoric flow this is
    // delta_lambda * sqrt(2/3); the volumetric (alpha) contribution adds
    // (alpha/sqrt(3))^2 in quadrature.
    const double e2 = 2.0 / 3.0
                      + (2.0 / 3.0) * params.alpha_dp * params.alpha_dp;
    res.delta_eps_p_eq = dlambda * std::sqrt(e2 > 0.0 ? e2 : 0.0);
    res.yielded = true;
    return res;
}

}  // namespace FSRM

#endif  // NEAR_FIELD_DRUCKER_PRAGER_3D_HPP
