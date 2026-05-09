/**
 * @file TillotsonEOS.hpp
 * @brief Tillotson equation of state for the host rock under
 *        nuclear-cavity hydrodynamic conditions.
 *
 * The Tillotson form (Tillotson 1962, "Metallic equations of state for
 * hypervelocity impact", General Atomic Report GA-3216) is an analytic
 * EOS that handles the four physical regimes a nuclear-cavity host
 * rock passes through during the early hydrodynamic phase: cold
 * compressed solid, cold expanded solid, hot expanded vapor, and the
 * mixed (partial vaporization) intermediate. It is calibrated against
 * shock-Hugoniot data and approaches the ideal-gas limit at high
 * specific internal energy.
 *
 * Pass-7 (axis 1 calibration follow-up) replaces the pass-6 ideal-gas
 * inner-cavity EOS with this Tillotson form for the host rock. The
 * inner cavity is no longer "detonation gas at gamma = 1.4"; it is
 * fully-vaporized rock plasma whose state at the radiation-to-
 * hydrodynamic transition time is set by a first-principles energy
 * partition (see RadialLagrangianSolver::solveCavityInitialState in
 * RadialLagrangian.cpp). This gets the (p, rho, e) curve into the
 * right physical regime: rock-plasma at temperatures four to six
 * orders of magnitude above chemical-detonation conditions, not
 * chemical-detonation products at energy-equivalent yield.
 *
 * The class is intentionally simple. All state-evaluation methods are
 * pure functions of (rho, e); no internal mutable state. Parameter
 * sets are static constexpr structs in this header so the class is
 * header-only-configurable.
 *
 * State convention: pressure p > 0 in compression (consistent with the
 * radial Lagrangian solver's gas-cell and Mie-Gruneisen conventions).
 * The Tillotson formulation in compression returns positive p when
 * either A * mu (cold pressure) or the thermal a + b/(1+e/(E0 eta^2))
 * term contributes positively.
 *
 * References.
 *  - Tillotson, J. H. (1962), "Metallic equations of state for
 *    hypervelocity impact", General Atomic Report GA-3216.
 *  - Melosh, H. J. (1989), "Impact Cratering: A Geologic Process",
 *    Oxford University Press, eqs 5.4.7-5.4.9 and Table A2.2 (the
 *    canonical reference for granite, basalt, and salt parameter
 *    sets used in cratering simulations).
 *  - Marsh, S. P. (ed., 1980), "LASL Shock Hugoniot Data", University
 *    of California Press (granite Hugoniot data used for the
 *    GraniteHugoniotCompression validation gate).
 *  - Trunin, R. F., Gudarenko, L. F., Zhernokletov, M. V., Simakov,
 *    G. V. (2001), "Experimental data on shock compressibility and
 *    adiabatic expansion of condensed substances", RFNC-VNIIEF
 *    (additional shock data for tuff and salt).
 *  - Carter, W. J. (1979), "Equation-of-state and shock-initiation
 *    investigations of NaCl", Los Alamos report (salt parameters).
 */

#ifndef NEAR_FIELD_TILLOTSON_EOS_HPP
#define NEAR_FIELD_TILLOTSON_EOS_HPP

#include <cmath>
#include <string>

namespace FSRM {

/**
 * @brief Tillotson EOS parameter set for a single host-rock medium.
 *
 * The 10-parameter Tillotson form. Units are SI throughout: density in
 * kg/m^3, specific internal energy in J/kg, pressures derived in Pa.
 *
 * Compressed and cold-expanded regime (rho >= rho_0 OR e < E_iv):
 *   p = (a + b / (1 + e / (E_0 eta^2))) rho e + A mu + B mu^2
 * with eta = rho / rho_0 and mu = eta - 1.
 *
 * Hot-expanded regime (rho < rho_0 AND e > E_cv):
 *   p = a rho e
 *     + ( b rho e / (1 + e / (E_0 eta^2))
 *         + A mu exp(-beta (rho_0/rho - 1)) ) exp(-alpha (rho_0/rho - 1)^2)
 *
 * Mixed regime (rho < rho_0 AND E_iv <= e <= E_cv):
 *   linear interpolation in e between the two expanded forms.
 */
struct TillotsonParameters
{
    double rho_0 = 2680.0;   ///< Reference (uncompressed solid) density [kg/m^3].
    double A = 1.8e10;       ///< Bulk modulus at p = 0 [Pa].
    double B = 1.8e10;       ///< Higher-order compression coefficient [Pa].
    double a = 0.5;          ///< Tillotson "a" coefficient (dimensionless).
    double b = 1.3;          ///< Tillotson "b" coefficient (dimensionless).
    double alpha = 5.0;      ///< Hot-expanded regime decay exponent.
    double beta = 5.0;       ///< Hot-expanded cold-pressure decay exponent.
    double E_0 = 1.6e7;      ///< Reference specific energy [J/kg].
    double E_iv = 3.5e6;     ///< Incipient-vaporization specific energy [J/kg].
    double E_cv = 1.8e7;     ///< Complete-vaporization specific energy [J/kg].
    double cv = 1.0e3;       ///< Specific heat at constant volume [J/(kg*K)].
    const char* name = "GENERIC";  ///< Diagnostic label for the parameter set.
};

/**
 * @brief Tabulated Tillotson parameter sets for the four media the
 *        historic-nuclear pipeline supports today.
 *
 * Each set bundles a citation in the comment block. Where a published
 * Tillotson set is not directly available (alluvium), the placeholder
 * is documented so it cannot be silently mistaken for production data.
 */
struct TillotsonParameterSets
{
    /// Granite. Melosh (1989) "Impact Cratering: A Geologic Process",
    /// Table A2.2 (the canonical cratering-literature parameter set).
    /// Validated against Marsh (1980) LASL Hugoniot data; the
    /// `Physics.TillotsonEOS.GraniteHugoniotCompression` gate locks
    /// this in to within 10 percent across 10-300 GPa.
    static constexpr TillotsonParameters granite()
    {
        TillotsonParameters p;
        p.rho_0 = 2680.0;
        p.A = 1.8e10;
        p.B = 1.8e10;
        p.a = 0.5;
        p.b = 1.3;
        p.alpha = 5.0;
        p.beta = 5.0;
        p.E_0 = 1.6e7;
        p.E_iv = 3.5e6;
        p.E_cv = 1.8e7;
        p.cv = 1.0e3;
        p.name = "GRANITE";
        return p;
    }

    /// Tuff. Derived from a Tillotson fit to Trunin et al. (2001)
    /// shock data for compacted volcanic tuff (rho_0 ~ 2000 kg/m^3,
    /// vp ~ 3500 m/s). Where the published Tillotson coefficients
    /// for tuff are not available the dimensionless "a", "b",
    /// "alpha", "beta" are taken from Melosh's basalt set (the
    /// closest-match volcanic glass) and the bulk-modulus parameters
    /// A, B and the energy scales are derived from the bulk sound
    /// speed and the latent-heat content per unit mass appropriate
    /// to a hydrated tuff matrix. The
    /// `Physics.TillotsonEOS.SaltAndAlluviumParameterSetSelfConsistency`
    /// gate exercises this set for thermodynamic self-consistency
    /// rather than against a primary Hugoniot fit.
    static constexpr TillotsonParameters tuff()
    {
        TillotsonParameters p;
        p.rho_0 = 2000.0;
        p.A = p.rho_0 * 3500.0 * 3500.0;
        p.B = p.A;
        p.a = 0.5;
        p.b = 1.3;
        p.alpha = 5.0;
        p.beta = 5.0;
        p.E_0 = 1.0e7;
        p.E_iv = 3.0e6;
        p.E_cv = 1.5e7;
        p.cv = 1.0e3;
        p.name = "TUFF";
        return p;
    }

    /// Salt (NaCl). Carter (1979) "Equation-of-state and
    /// shock-initiation investigations of NaCl", Los Alamos LA-7873;
    /// dimensionless coefficients from Melosh (1989) Table A2.2 NaCl
    /// row.
    static constexpr TillotsonParameters salt()
    {
        TillotsonParameters p;
        p.rho_0 = 2160.0;
        p.A = 2.5e10;
        p.B = 3.0e10;
        p.a = 0.5;
        p.b = 1.5;
        p.alpha = 5.0;
        p.beta = 5.0;
        p.E_0 = 1.5e7;
        p.E_iv = 1.5e6;
        p.E_cv = 8.0e6;
        p.cv = 0.85e3;
        p.name = "SALT";
        return p;
    }

    /// Alluvium. KNOWN GAP: no published Tillotson parameter set
    /// for alluvium-class soft sediments was located in the
    /// open literature. This set is granite's dimensionless
    /// constants scaled to alluvium's density (~ 1800 kg/m^3) and
    /// reduced bulk modulus (alluvium vp ~ 2400 m/s gives
    /// K ~ rho vp^2 ~ 1.0 GPa). The vaporization energy scales are
    /// scaled down with the bulk modulus so the qualitative
    /// vaporization-front behaviour is reasonable. Pass-8 should
    /// replace this with a fitted Hugoniot from the Carmichael
    /// (1989) Practical Handbook of Physical Properties or a
    /// purpose-built fit to Yucca Flat alluvium shock data.
    static constexpr TillotsonParameters alluvium()
    {
        TillotsonParameters p;
        p.rho_0 = 1800.0;
        p.A = p.rho_0 * 2400.0 * 2400.0;
        p.B = p.A;
        p.a = 0.5;
        p.b = 1.3;
        p.alpha = 5.0;
        p.beta = 5.0;
        p.E_0 = 5.0e6;
        p.E_iv = 1.0e6;
        p.E_cv = 5.0e6;
        p.cv = 1.2e3;
        p.name = "ALLUVIUM";
        return p;
    }

    /// Look up a parameter set by uppercase name. Returns granite()
    /// (the best-validated set) if the name is unrecognised; the
    /// caller is responsible for warning on the fall-through.
    static TillotsonParameters byName(const std::string& upper_name)
    {
        if (upper_name == "GRANITE") return granite();
        if (upper_name == "TUFF") return tuff();
        if (upper_name == "SALT") return salt();
        if (upper_name == "ALLUVIUM") return alluvium();
        return granite();
    }
};

/**
 * @brief Tillotson EOS evaluator. Pure-function methods on (rho, e).
 *
 * The class holds a single TillotsonParameters struct. All evaluations
 * branch internally on the (rho, e) regime per the Tillotson 1962
 * formulation.
 */
class TillotsonEOS
{
public:
    TillotsonEOS() = default;
    explicit TillotsonEOS(const TillotsonParameters& p) : params_(p) {}

    void setParameters(const TillotsonParameters& p) { params_ = p; }
    const TillotsonParameters& getParameters() const { return params_; }

    /// Pressure (compression positive) at (rho, e) [Pa]. Branches on
    /// regime per Tillotson 1962.
    double pressure(double rho, double e) const
    {
        const double rho_safe = (rho > 1.0e-12) ? rho : 1.0e-12;
        const double e_safe = (e > 0.0) ? e : 0.0;
        const double eta = rho_safe / params_.rho_0;
        const double mu = eta - 1.0;
        const bool compressed_or_cold =
            (rho_safe >= params_.rho_0) || (e_safe < params_.E_iv);
        const double p_compressed =
            pressureCompressed(rho_safe, e_safe, eta, mu);
        if (compressed_or_cold) {
            return p_compressed;
        }
        const double p_expanded_hot =
            pressureExpandedHot(rho_safe, e_safe, eta, mu);
        if (e_safe >= params_.E_cv) {
            return p_expanded_hot;
        }
        // Mixed regime: linear interpolation in e between cold-expanded
        // (which uses the compressed form) and hot-expanded.
        const double w = (e_safe - params_.E_iv) /
                         (params_.E_cv - params_.E_iv);
        return (1.0 - w) * p_compressed + w * p_expanded_hot;
    }

    /// Thermodynamic sound speed [m/s] at (rho, e). Computed from the
    /// standard relation
    ///   c^2 = (dp/drho)_e + (p / rho^2) (dp/de)_rho.
    /// The partial derivatives are taken numerically by central
    /// differences; the cost is dominated by four pressure() calls
    /// which are themselves O(1).
    double soundSpeed(double rho, double e) const
    {
        const double rho_safe = (rho > 1.0e-12) ? rho : 1.0e-12;
        const double e_safe = (e > 0.0) ? e : 0.0;
        const double drho = 1.0e-3 * rho_safe;
        const double de = 1.0e-3 * (e_safe + params_.E_0);
        const double p = pressure(rho_safe, e_safe);
        const double dp_drho = (pressure(rho_safe + drho, e_safe) -
                                pressure(rho_safe - drho, e_safe)) /
                               (2.0 * drho);
        const double dp_de = (pressure(rho_safe, e_safe + de) -
                              pressure(rho_safe, e_safe - de)) /
                             (2.0 * de);
        const double c2 = dp_drho + (p / (rho_safe * rho_safe)) * dp_de;
        return std::sqrt(c2 > 0.0 ? c2 : 0.0);
    }

    /// Approximate temperature [K] at (rho, e). Uses cv as the local
    /// reference: T ~ e / cv in the hot expanded regime (where the EOS
    /// is approximately ideal-gas-like with effective gamma = 1 + a),
    /// shifted by the latent-heat content in the cold and mixed
    /// regimes. This is a coarse approximation suitable for
    /// diagnostics; replace with a tabulated T(rho, e) when a SESAME-
    /// or QEOS-class table is plumbed.
    double temperature(double rho, double e) const
    {
        const double e_safe = (e > 0.0) ? e : 0.0;
        if (rho < params_.rho_0 && e_safe > params_.E_cv) {
            return e_safe / params_.cv;
        }
        if (rho < params_.rho_0 && e_safe > params_.E_iv) {
            return (e_safe - params_.E_iv) / params_.cv + 300.0;
        }
        return e_safe / params_.cv + 300.0;
    }

private:
    double pressureCompressed(double rho, double e,
                              double eta, double mu) const
    {
        const double eta2 = eta * eta;
        const double e_term = 1.0 + e / (params_.E_0 * eta2 + 1.0e-30);
        const double thermal_coef = params_.a + params_.b / e_term;
        const double thermal = thermal_coef * rho * e;
        const double cold = params_.A * mu + params_.B * mu * mu;
        return thermal + cold;
    }

    double pressureExpandedHot(double rho, double e,
                               double eta, double mu) const
    {
        const double eta2 = eta * eta;
        const double e_term = 1.0 + e / (params_.E_0 * eta2 + 1.0e-30);
        const double inv_eta_minus_1 = (params_.rho_0 / rho) - 1.0;
        const double exp_alpha =
            std::exp(-params_.alpha * inv_eta_minus_1 * inv_eta_minus_1);
        const double exp_beta = std::exp(-params_.beta * inv_eta_minus_1);
        const double thermal_iso = params_.a * rho * e;
        const double thermal_extra =
            (params_.b * rho * e / e_term + params_.A * mu * exp_beta) *
            exp_alpha;
        return thermal_iso + thermal_extra;
    }

    TillotsonParameters params_ = TillotsonParameterSets::granite();
};

} // namespace FSRM

#endif // NEAR_FIELD_TILLOTSON_EOS_HPP
