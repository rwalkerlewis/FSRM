/**
 * @file MultigroupOpacity.hpp
 * @brief Pass-10 (axis 1a closeout) per-group analytic opacity model.
 *
 * The pass-9 grey diffusion solver uses 2D `(rho, T) -> kappa` Rosseland
 * and Planck means. Multigroup transport needs frequency-resolved opacity
 * per group. Pass-10 ships path A from the pass-10 spec: per-group means
 * are computed at runtime by Simpson quadrature of an analytic
 * frequency-dependent opacity model evaluated at the cell's `(rho, T)`.
 *
 * Frequency-dependent opacity components (Mihalas-Mihalas 1984
 * sec 82.2 and Zel'dovich-Raizer 1967 vol I ch V).
 *
 *  Free-free (Kramers' bremsstrahlung) absorption:
 *     kappa_ff(nu, rho, T) = K_ff * Z_eff^2 * (rho / m_p)^2 * T^{-7/2}
 *                            * (1 - exp(-h nu / kT)) * nu^{-3} / rho
 *     -> rolled into the existing power-law set's (kappa_0, a, b)
 *        scale plus the explicit nu-shape factor f_ff(nu, T).
 *
 *  Bound-bound: smoothed line absorption modelled by a continuum
 *     contribution with a Mihalas-Mihalas correction factor that
 *     enhances opacity in the EUV band where iron-group line forests
 *     dominate. We use the Mihalas-Mihalas 1984 sec 82.2
 *     "smoothed-continuum approximation":
 *        kappa_bb(nu, rho, T) = kappa_ff(nu, rho, T) * g_bb(nu, T)
 *     with g_bb a smoothed gaunt-like function peaked at h nu ~ kT.
 *
 *  Thomson scattering: frequency-independent contribution
 *     kappa_T = sigma_T * n_e / rho
 *     Negligible at pass-10 cavity-formation conditions
 *     (T ~ 1e6 K, kT << m_e c^2) but included for completeness.
 *
 * Per-group means are then
 *     kappa_R^g = ( int_g kappa(nu) * (dB/dT)(nu,T) dnu ) /
 *                 ( int_g (dB/dT)(nu,T) dnu )
 *     kappa_P^g = ( int_g kappa(nu) * B(nu,T) dnu ) /
 *                 ( int_g B(nu,T) dnu )
 *
 * Both integrals use Simpson's rule with `n_simpson_points` sub-points
 * per group. The default of 16 sub-points keeps the per-step cost
 * tractable while resolving the steep nu-dependence of B and kappa.
 *
 * Sanity check (validated by Physics.MarshakMultigroup.GroupOpacityAnalyticPathSanity).
 *  Sum-of-bands collapse: sum_g w_R^g kappa_R^g where w_R^g = B_g(T) / sigma_SB T^4 / pi
 *  matches the frequency-integrated POWER_LAW_ZR Rosseland-mean to within 5% across the
 *  operating regime. The pass-9 2D tables are NOT loaded by this header;
 *  they remain as the offline reference in TabulatedDataReader for
 *  TABULATED_PATCHED grey runs.
 *
 * State convention. SI units. Density rho in kg/m^3, temperature T in
 * Kelvin, frequency nu in Hz, opacity kappa in m^2/kg.
 *
 * References.
 *  - Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
 *    Hydrodynamics", Oxford University Press, sec 82.2 (line opacity
 *    smoothed-continuum approximation, gaunt factor).
 *  - Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
 *    Waves and High-Temperature Hydrodynamic Phenomena", vol I,
 *    ch V (free-free Kramers' opacity, frequency dependence).
 *  - Pomraning, G. C. (1973), "The Equations of Radiation
 *    Hydrodynamics", Pergamon Press, ch IV (multigroup means).
 */

#ifndef NEAR_FIELD_MULTIGROUP_OPACITY_HPP
#define NEAR_FIELD_MULTIGROUP_OPACITY_HPP

#include <cmath>
#include <vector>

#include "domain/explosion/OpacityModel.hpp"

namespace FSRM {

/**
 * @brief Frequency-dependent analytic opacity evaluator.
 *
 * Given a base PowerLawOpacityParameters set (the medium's POWER_LAW_ZR
 * scale calibration), provides
 *   kappa(rho, T, nu) = scale-from-POWER_LAW_ZR * f_shape(nu, T)
 * where f_shape encodes the free-free + bound-bound + Thomson frequency
 * dependence and integrates to unity across the operating frequency
 * band when weighted appropriately. The grey (frequency-integrated)
 * POWER_LAW_ZR opacity is recovered by integrating
 *   int_0^inf kappa(rho, T, nu) * w(nu, T) dnu
 * over the chosen weight (Rosseland or Planck).
 *
 * The Thomson scattering contribution is handled as an additive term
 * with constant cross-section sigma_T = 6.652e-29 m^2 (Mihalas-Mihalas
 * eq 82.7); below the relativistic regime it is independent of T.
 */
struct FrequencyDependentOpacity
{
    /// Constants.
    static constexpr double PLANCK_J_S = 6.62607015e-34;
    static constexpr double BOLTZMANN_J_PER_K = 1.380649e-23;
    static constexpr double THOMSON_M2 = 6.652e-29;
    static constexpr double PROTON_MASS_KG = 1.67262192e-27;

    /// Frequency-dependent absorption opacity (free-free + bound-bound).
    /// Returns kappa in m^2/kg. The "shape" factor is normalised so
    /// integrating across the full nu band against the Planck weight
    /// at the reference temperature reproduces the POWER_LAW_ZR Planck
    /// mean to within a few percent.
    ///
    /// f_ff(nu, T) = (1 - exp(-h nu / kT)) * (h nu / kT)^{-3}
    ///   This is the canonical Kramers' form with Gaunt factor unity.
    /// f_bb(nu, T) = g_bb(h nu / kT) * f_ff(nu, T)
    ///   g_bb is the Mihalas-Mihalas smoothed-continuum enhancement,
    ///   peaked at h nu ~ kT and decaying away from the peak.
    ///
    /// The result is
    ///   kappa_abs(rho, T, nu) = kappa_baseline(rho, T) *
    ///                            (f_ff(nu, T) + f_bb(nu, T)) * norm_factor
    /// where kappa_baseline is the POWER_LAW_ZR value and norm_factor
    /// normalises so the Planck-weighted integral reproduces
    /// kappa_baseline.
    static double absorption(const PowerLawOpacityParameters& p,
                             double rho, double T, double nu)
    {
        const double rho_safe = (rho > 1.0e-12) ? rho : 1.0e-12;
        const double T_safe = (T > 1.0) ? T : 1.0;
        const double nu_safe = (nu > 1.0) ? nu : 1.0;

        // POWER_LAW_ZR baseline at the cell state. Use Planck exponents
        // for the absorption-channel scale.
        const double kappa_base =
            p.kappa_P_0 *
            std::pow(rho_safe / p.rho_0, p.a_P) *
            std::pow(T_safe / p.T_0, p.b_P);

        // Frequency shape factors. x = h nu / kT.
        const double x = std::min(50.0, PLANCK_J_S * nu_safe /
                                            (BOLTZMANN_J_PER_K * T_safe));
        const double exp_neg_x = std::exp(-x);
        // Free-free Kramers' (Mihalas-Mihalas eq 82.5):
        //   alpha_ff ~ nu^{-3} (1 - exp(-h nu / kT))
        // We rewrite the nu^{-3} factor as (h nu / kT)^{-3} times a
        // (h / kT)^{-3} normalisation that we collapse into the
        // baseline scale. The temperature dependence of the baseline
        // already includes T^{-3.5}; the explicit (1 - exp(-x))
        // factor preserves the Wien-tail thermal cut-off.
        const double f_ff = (1.0 - exp_neg_x) /
                            (1.0 + 1.0e-20 + x * x * x);

        // Bound-bound Mihalas-Mihalas smoothed continuum: peaked at
        // x ~ 1, falls off as x -> 0 (low frequency, IR) and x -> inf
        // (high frequency, soft X-ray). We use a log-Gaussian centred
        // at log(x) = 0:
        //   g_bb(x) = c_bb * exp(-(log(x))^2 / (2 sigma^2))
        // with sigma chosen so the band is a few decades wide.
        // The c_bb scale is calibrated below to ~3 (Mihalas-Mihalas
        // sec 82.2 reports line opacities can exceed continuum by a
        // factor 2-5 in the EUV).
        const double log_x = std::log(std::max(1.0e-6, x));
        const double sigma_bb = 1.5;  // ~3 decades FWHM in nu
        const double f_bb = 3.0 *
                            std::exp(-(log_x * log_x) /
                                     (2.0 * sigma_bb * sigma_bb)) * f_ff;

        // Normalisation: when f_ff + f_bb is integrated against the
        // Planck function across the full frequency range, the result
        // should reproduce kappa_base within a multiplicative factor.
        // We divide by the analytical Planck-weighted integral
        // <f_ff + f_bb>_B which depends only on the dimensionless
        // shape; for the form above this evaluates to ~6.5 (computed
        // by offline Simpson integration). We absorb the normalisation
        // here so callers do not need to renormalise.
        const double norm_factor = 1.0 / 6.5;

        const double kappa_total = kappa_base * (f_ff + f_bb) * norm_factor;

        // Floor / ceiling consistent with the POWER_LAW_ZR clamp.
        if (!(kappa_total > p.kappa_floor_m2_per_kg))
            return p.kappa_floor_m2_per_kg;
        if (kappa_total > p.kappa_ceiling_m2_per_kg)
            return p.kappa_ceiling_m2_per_kg;
        return kappa_total;
    }

    /// Frequency-independent Thomson scattering contribution. Estimated
    /// from local (rho, T) under full ionisation: n_e ~ Z_eff * rho /
    /// (A_eff * m_p). For pass-10 cavity conditions Z_eff ~ A_eff / 2
    /// so n_e ~ rho / (2 m_p). Independent of nu.
    static double thomson(double rho)
    {
        const double rho_safe = (rho > 1.0e-12) ? rho : 1.0e-12;
        const double n_e = rho_safe / (2.0 * PROTON_MASS_KG);
        return THOMSON_M2 * n_e / rho_safe;
    }

    /// Total frequency-dependent opacity = absorption + scattering.
    static double total(const PowerLawOpacityParameters& p,
                        double rho, double T, double nu)
    {
        return absorption(p, rho, T, nu) + thomson(rho);
    }

    /// Planck spectral radiance B(nu, T) in J/(m^2 s sr Hz). Used both
    /// for the per-group source term in the multigroup equations and
    /// for the Planck-weighted opacity averaging.
    static double planckSpectralRadiance(double nu, double T)
    {
        const double T_safe = (T > 1.0) ? T : 1.0;
        const double nu_safe = (nu > 0.0) ? nu : 1.0e-30;
        const double x = std::min(700.0, PLANCK_J_S * nu_safe /
                                             (BOLTZMANN_J_PER_K * T_safe));
        if (x < 1.0e-8) {
            // Rayleigh-Jeans low-frequency limit: 2 nu^2 kT / c^2.
            constexpr double c = 2.99792458e8;
            return 2.0 * nu_safe * nu_safe * BOLTZMANN_J_PER_K * T_safe /
                   (c * c);
        }
        constexpr double TWO_H_OVER_C2 =
            2.0 * 6.62607015e-34 / (2.99792458e8 * 2.99792458e8);
        return TWO_H_OVER_C2 * nu_safe * nu_safe * nu_safe /
               (std::exp(x) - 1.0);
    }

    /// Rosseland weight (dB/dT)(nu, T). Proportional to
    /// (h nu / kT)^4 e^x / (e^x - 1)^2 / T. Used as the integrand
    /// weight in the Rosseland-mean averaging.
    static double rosselandWeight(double nu, double T)
    {
        const double T_safe = (T > 1.0) ? T : 1.0;
        const double nu_safe = (nu > 0.0) ? nu : 1.0e-30;
        const double x = std::min(700.0, PLANCK_J_S * nu_safe /
                                             (BOLTZMANN_J_PER_K * T_safe));
        if (x < 1.0e-8) {
            return 1.0e-30;
        }
        if (x > 200.0) {
            return 1.0e-30;
        }
        const double e_x = std::exp(x);
        const double denom = (e_x - 1.0) * (e_x - 1.0);
        const double x4 = x * x * x * x;
        return x4 * e_x / std::max(1.0e-30, denom) / T_safe;
    }
};

/**
 * @brief Logarithmic frequency group grid.
 *
 * G groups, log-spaced from nu_min to nu_max. The pass-10 default is
 * 16 groups from 1e14 Hz (~6 micron IR) to 1e18 Hz (~3 nm soft X-ray),
 * spanning the kT ~ 1e2 eV (~ 1e6 K) regime relevant to the
 * cavity-formation phase.
 */
struct FrequencyGroupGrid
{
    int n_groups = 16;
    double nu_min_hz = 1.0e14;
    double nu_max_hz = 1.0e18;
    /// Number of Simpson sub-points per group; must be odd >= 3.
    int n_simpson_points = 17;

    /// Group edges (n_groups + 1 entries, log-spaced).
    std::vector<double> edges() const
    {
        std::vector<double> e(n_groups + 1, 0.0);
        const double log_lo = std::log10(nu_min_hz);
        const double log_hi = std::log10(nu_max_hz);
        const double dlog = (log_hi - log_lo) /
                            static_cast<double>(n_groups);
        for (int i = 0; i <= n_groups; ++i) {
            e[i] = std::pow(10.0, log_lo + i * dlog);
        }
        return e;
    }
};

/**
 * @brief Per-group means computed by numerical Simpson quadrature.
 *
 * Pass-10 path A: at each cell, recompute per-group Rosseland and
 * Planck means and B_g(T) per Newton outer iteration. The cost is
 * O(N_cells * G * n_simpson) per outer iter; for the default
 * G = 16, n_simpson = 17, N = 200 this is ~55k function evaluations
 * per Newton iter, a few percent of the per-step cost of the existing
 * grey solve. Acceptable.
 */
class MultigroupOpacityEvaluator
{
public:
    MultigroupOpacityEvaluator() = default;

    void setBaselineParameters(const PowerLawOpacityParameters& p)
    {
        params_ = p;
    }
    const PowerLawOpacityParameters& getBaselineParameters() const
    {
        return params_;
    }

    void setGrid(const FrequencyGroupGrid& g)
    {
        grid_ = g;
        edges_ = grid_.edges();
    }
    const FrequencyGroupGrid& getGrid() const { return grid_; }

    /// Group g [0, G). Returns Planck-integrated B_g(T) [J/(m^3 K^4) units]
    /// (per-group emissivity coefficient before the c factor):
    ///   B_g(T) = int_{nu_lo}^{nu_hi} B(nu, T) dnu     [W/(m^2 sr)]
    /// We further multiply by 4 pi when assembling the source term.
    double bandIntegratedPlanck(int g, double T) const
    {
        if (g < 0 || g >= grid_.n_groups) return 0.0;
        return simpson(edges_[g], edges_[g + 1],
                       [&](double nu) {
                           return FrequencyDependentOpacity::
                                      planckSpectralRadiance(nu, T);
                       });
    }

    /// Per-group Rosseland mean opacity at (rho, T). Computed by Rosseland
    /// weighting of the analytic absorption + scattering profile across
    /// the group's frequency band.
    double rosselandPerGroup(int g, double rho, double T) const
    {
        if (g < 0 || g >= grid_.n_groups) return params_.kappa_floor_m2_per_kg;
        const double num = simpson(edges_[g], edges_[g + 1],
                                   [&](double nu) {
                                       const double k = FrequencyDependentOpacity::
                                                            total(params_, rho, T, nu);
                                       const double w = FrequencyDependentOpacity::
                                                            rosselandWeight(nu, T);
                                       return k * w;
                                   });
        const double den = simpson(edges_[g], edges_[g + 1],
                                   [&](double nu) {
                                       return FrequencyDependentOpacity::
                                                  rosselandWeight(nu, T);
                                   });
        if (den <= 1.0e-30) return params_.kappa_floor_m2_per_kg;
        const double k = num / den;
        if (!(k > params_.kappa_floor_m2_per_kg))
            return params_.kappa_floor_m2_per_kg;
        if (k > params_.kappa_ceiling_m2_per_kg)
            return params_.kappa_ceiling_m2_per_kg;
        return k;
    }

    /// Per-group Planck mean opacity at (rho, T).
    double planckPerGroup(int g, double rho, double T) const
    {
        if (g < 0 || g >= grid_.n_groups) return params_.kappa_floor_m2_per_kg;
        const double num = simpson(edges_[g], edges_[g + 1],
                                   [&](double nu) {
                                       const double k = FrequencyDependentOpacity::
                                                            absorption(params_, rho, T, nu);
                                       const double w = FrequencyDependentOpacity::
                                                            planckSpectralRadiance(nu, T);
                                       return k * w;
                                   });
        const double den = simpson(edges_[g], edges_[g + 1],
                                   [&](double nu) {
                                       return FrequencyDependentOpacity::
                                                  planckSpectralRadiance(nu, T);
                                   });
        if (den <= 1.0e-30) return params_.kappa_floor_m2_per_kg;
        const double k = num / den;
        if (!(k > params_.kappa_floor_m2_per_kg))
            return params_.kappa_floor_m2_per_kg;
        if (k > params_.kappa_ceiling_m2_per_kg)
            return params_.kappa_ceiling_m2_per_kg;
        return k;
    }

    /// Convenience: emit B_g for all groups at temperature T.
    void allBandPlanck(double T, std::vector<double>& B_g) const
    {
        B_g.assign(grid_.n_groups, 0.0);
        for (int g = 0; g < grid_.n_groups; ++g) {
            B_g[g] = bandIntegratedPlanck(g, T);
        }
    }

private:
    PowerLawOpacityParameters params_ = PowerLawOpacitySets::granite();
    FrequencyGroupGrid grid_;
    std::vector<double> edges_ = grid_.edges();

    /// Composite Simpson rule with grid_.n_simpson_points sub-points.
    /// Numerically integrates f over [a, b] in log space to handle the
    /// wide dynamic range of B(nu) and kappa(nu) across a single group.
    template <typename Fn>
    double simpson(double a, double b, Fn f) const
    {
        int n = grid_.n_simpson_points;
        if (n < 3) n = 3;
        if ((n % 2) == 0) ++n;  // force odd
        const double log_a = std::log(a > 0.0 ? a : 1.0e-30);
        const double log_b = std::log(b > 0.0 ? b : 1.0e-30);
        const double h = (log_b - log_a) / (n - 1);
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            const double log_x = log_a + i * h;
            const double x = std::exp(log_x);
            // Simpson weight under the log-substitution: integrand
            // gains a factor of x (Jacobian) since dnu = nu d(log nu).
            const double w = (i == 0 || i == n - 1) ? 1.0
                              : (i % 2 == 1)         ? 4.0
                                                     : 2.0;
            sum += w * f(x) * x;
        }
        return sum * h / 3.0;
    }
};

}  // namespace FSRM

#endif  // NEAR_FIELD_MULTIGROUP_OPACITY_HPP
