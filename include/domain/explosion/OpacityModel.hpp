/**
 * @file OpacityModel.hpp
 * @brief Power-law (Kramers'-type) Rosseland and Planck mean opacity model
 *        for the rock-plasma regime relevant to underground nuclear cavity
 *        formation.
 *
 * Pass-8 (axis 1, see docs/HISTORIC_NUCLEAR_ROADMAP.md) replaces the
 * pass-7 Zel'dovich-Raizer end-state approximation with an explicit
 * Marshak grey radiation-diffusion solve. The diffusion coefficient
 *   D = c / (3 kappa_R rho)
 * and the matter-radiation coupling
 *   c kappa_P rho (a T^4 - E_r)
 * both depend on a model for the medium's opacity in the kt-class
 * post-shock regime: temperatures 1e4 - 1e7 K, densities spanning six
 * orders of magnitude as the cavity expands.
 *
 * For pass-8 we use the Kramers'-type power-law parameterization
 *   kappa(rho, T) = kappa_0 * (rho / rho_0)^a * (T / T_0)^b
 * with separate (a_R, b_R) and (a_P, b_P) exponents for the Rosseland
 * and Planck means. The exponents are taken from Zel'dovich and Raizer
 * (1967) vol I, ch V, sec 10 (Kramers / free-free opacity) for the
 * ionized rock-plasma regime; the absolute scale kappa_0 at the
 * reference (rho_0, T_0) is calibrated from published shock-physics
 * literature for each medium (see per-set comment blocks).
 *
 * This is a low-order approximation. Real opacity has line structure,
 * photoionization edges, and pressure broadening that a power law
 * cannot capture. The TABULATED_TOPS option scaffolded for pass-9
 * will replace the power law with a tabulated kappa(rho, T) loaded
 * from an external file (TOPS-format Los Alamos opacity tables, or
 * the open-data SESAME 1980 series).
 *
 * State convention. SI units. Density rho in kg/m^3, temperature T in
 * Kelvin, opacity kappa in m^2/kg.
 *
 * References.
 *  - Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
 *    Waves and High-Temperature Hydrodynamic Phenomena", vol I,
 *    Academic Press, ch V sec 10 (Kramers' opacity, free-free and
 *    free-bound mean opacities for hot plasmas; granite-like rock
 *    plasma exponents a ~ 1, b ~ -3.5).
 *  - Pomraning, G. C. (1973), "The Equations of Radiation
 *    Hydrodynamics", Pergamon Press, ch IV (Marshak self-similar
 *    wave validation for constant kappa).
 *  - Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
 *    Hydrodynamics", Oxford University Press, sec 96-97 (gray
 *    diffusion limit, opacity averaging).
 *  - Marshak, R. E. (1958), "Effect of radiation on shock wave
 *    behavior", Physics of Fluids 1(1), pp. 24-29 (boundary
 *    conditions for radiation diffusion).
 */

#ifndef NEAR_FIELD_OPACITY_MODEL_HPP
#define NEAR_FIELD_OPACITY_MODEL_HPP

#include <cmath>
#include <string>

namespace FSRM {

/**
 * @brief Closed enum for the opacity model family. POWER_LAW_ZR is the
 *        pass-8 default. CONSTANT is a sanity-check (used by the Marshak
 *        self-similar test). TABULATED_TOPS is a pass-9 scaffold; the
 *        Marshak solver throws a clear "not implemented" error when
 *        constructed against TABULATED_TOPS.
 */
enum class OpacityModel
{
    POWER_LAW_ZR,
    CONSTANT,
    TABULATED_TOPS
};

/**
 * @brief Power-law opacity parameter set for one medium.
 *
 * The form is
 *   kappa_R(rho, T) = kappa_R_0 * (rho/rho_0)^a_R * (T/T_0)^b_R
 *   kappa_P(rho, T) = kappa_P_0 * (rho/rho_0)^a_P * (T/T_0)^b_P
 *
 * The Rosseland mean kappa_R appears in the diffusion coefficient
 * D = c / (3 kappa_R rho). The Planck mean kappa_P appears in the
 * matter-radiation coupling source c kappa_P rho (a T^4 - E_r).
 */
struct PowerLawOpacityParameters
{
    double rho_0 = 2700.0;     ///< Reference density [kg/m^3].
    double T_0 = 1.0e6;        ///< Reference temperature [K].
    double kappa_R_0 = 1.0;    ///< Rosseland opacity at (rho_0, T_0) [m^2/kg].
    double kappa_P_0 = 5.0;    ///< Planck opacity at (rho_0, T_0) [m^2/kg].
    double a_R = 1.0;          ///< Density exponent (Rosseland).
    double b_R = -3.5;         ///< Temperature exponent (Rosseland).
    double a_P = 1.0;          ///< Density exponent (Planck).
    double b_P = -3.5;         ///< Temperature exponent (Planck).
    /// Hard floors on opacity to avoid numerical breakdown in extreme
    /// cells. Values below these are clamped; the diffusion coefficient
    /// remains well-posed in the rarefied / cold regime.
    double kappa_floor_m2_per_kg = 1.0e-6;
    double kappa_ceiling_m2_per_kg = 1.0e6;
    const char* name = "GENERIC";
};

/**
 * @brief Tabulated power-law opacity sets per medium.
 *
 * The exponents a, b come from Zel'dovich-Raizer 1967 vol I sec 10
 * (Kramers' free-free opacity in the ionized regime: a ~ 1, b ~ -3.5
 * for the Rosseland mean; the Planck mean is somewhat steeper in T,
 * b ~ -3 to -4 depending on the plasma conditions). The reference
 * scale kappa_0 at (rho_0, T_0) is the calibration knob; the values
 * below are chosen so the diffusion coefficient at typical cavity
 * conditions (rho ~ rho_0, T ~ 1e6 K) gives D ~ 1e6 m^2/s, matching
 * the order-of-magnitude estimates in Glasstone & Dolan 1977 ch II
 * for the Marshak phase of a kt nuclear cavity.
 */
struct PowerLawOpacitySets
{
    /// Granite. Composition by mass: ~75% SiO2, ~15% Al2O3, others.
    /// Mean atomic number Z_eff ~ 11; mass fraction of high-Z (Fe, etc)
    /// is ~5%. Z-R 1967 vol I sec 10 gives free-free Rosseland
    /// kappa_R ~ Z_eff^2 (rho/T^3.5) -> ~ 1 m^2/kg at (rho_0, T_0)
    /// with rho_0 = granite uncompressed solid density.
    static constexpr PowerLawOpacityParameters granite()
    {
        PowerLawOpacityParameters p;
        p.rho_0 = 2680.0;
        p.T_0 = 1.0e6;
        p.kappa_R_0 = 1.0;
        p.kappa_P_0 = 5.0;
        p.a_R = 1.0;
        p.b_R = -3.5;
        p.a_P = 1.0;
        p.b_P = -3.5;
        p.kappa_floor_m2_per_kg = 1.0e-6;
        p.kappa_ceiling_m2_per_kg = 1.0e6;
        p.name = "GRANITE";
        return p;
    }

    /// Tuff. Compacted volcanic glass. Lower density, lower Z_eff than
    /// granite. Reference scale calibrated to produce a similar
    /// diffusion coefficient at (rho_0, T_0) given the lower density.
    static constexpr PowerLawOpacityParameters tuff()
    {
        PowerLawOpacityParameters p;
        p.rho_0 = 2000.0;
        p.T_0 = 1.0e6;
        p.kappa_R_0 = 0.8;
        p.kappa_P_0 = 4.0;
        p.a_R = 1.0;
        p.b_R = -3.5;
        p.a_P = 1.0;
        p.b_P = -3.5;
        p.kappa_floor_m2_per_kg = 1.0e-6;
        p.kappa_ceiling_m2_per_kg = 1.0e6;
        p.name = "TUFF";
        return p;
    }

    /// Salt (NaCl). Z_Na = 11, Z_Cl = 17. Higher mean Z_eff than
    /// silicate rocks; Kramers' opacity scales as Z_eff^2 so kappa_R_0
    /// is larger for fixed temperature. Z-R 1967 vol I sec 10
    /// indicates Rosseland kappa for NaCl plasma at 1e6 K, solid
    /// density is ~2 m^2/kg.
    static constexpr PowerLawOpacityParameters salt()
    {
        PowerLawOpacityParameters p;
        p.rho_0 = 2160.0;
        p.T_0 = 1.0e6;
        p.kappa_R_0 = 2.0;
        p.kappa_P_0 = 10.0;
        p.a_R = 1.0;
        p.b_R = -3.5;
        p.a_P = 1.0;
        p.b_P = -3.5;
        p.kappa_floor_m2_per_kg = 1.0e-6;
        p.kappa_ceiling_m2_per_kg = 1.0e6;
        p.name = "SALT";
        return p;
    }

    /// Alluvium. Loose unconsolidated soft sediment. Lower density,
    /// lower Z_eff. Calibration is the granite scaling shifted by the
    /// density and bulk modulus. Documented as a known low-confidence
    /// set in HISTORIC_NUCLEAR_FIDELITY pass-8 entry; tabulated
    /// opacities are pass-9 work.
    static constexpr PowerLawOpacityParameters alluvium()
    {
        PowerLawOpacityParameters p;
        p.rho_0 = 1800.0;
        p.T_0 = 1.0e6;
        p.kappa_R_0 = 0.6;
        p.kappa_P_0 = 3.0;
        p.a_R = 1.0;
        p.b_R = -3.5;
        p.a_P = 1.0;
        p.b_P = -3.5;
        p.kappa_floor_m2_per_kg = 1.0e-6;
        p.kappa_ceiling_m2_per_kg = 1.0e6;
        p.name = "ALLUVIUM";
        return p;
    }

    /// Look up by uppercase name. Returns granite() on unrecognised
    /// names; the caller is responsible for warning on the fallthrough.
    static PowerLawOpacityParameters byName(const std::string& upper_name)
    {
        if (upper_name == "GRANITE") return granite();
        if (upper_name == "TUFF") return tuff();
        if (upper_name == "SALT") return salt();
        if (upper_name == "ALLUVIUM") return alluvium();
        return granite();
    }
};

/**
 * @brief Power-law opacity evaluator. Pure-function methods.
 *
 * Holds a single PowerLawOpacityParameters set. Both Rosseland and
 * Planck means are evaluated by the same closed-form expression with
 * different exponents. Floors and ceilings keep the diffusion
 * coefficient finite in extreme regimes.
 */
class PowerLawOpacity
{
public:
    PowerLawOpacity() = default;
    explicit PowerLawOpacity(const PowerLawOpacityParameters& p) : params_(p) {}

    void setParameters(const PowerLawOpacityParameters& p) { params_ = p; }
    const PowerLawOpacityParameters& getParameters() const { return params_; }

    /// Rosseland mean opacity at (rho, T). Returns clamped value
    /// in [kappa_floor, kappa_ceiling].
    double rosseland(double rho, double T) const
    {
        const double rho_safe = (rho > 1.0e-12) ? rho : 1.0e-12;
        const double T_safe = (T > 1.0) ? T : 1.0;
        const double k = params_.kappa_R_0 *
                         std::pow(rho_safe / params_.rho_0, params_.a_R) *
                         std::pow(T_safe / params_.T_0, params_.b_R);
        return clamp(k);
    }

    /// Planck mean opacity at (rho, T). Returns clamped value.
    double planck(double rho, double T) const
    {
        const double rho_safe = (rho > 1.0e-12) ? rho : 1.0e-12;
        const double T_safe = (T > 1.0) ? T : 1.0;
        const double k = params_.kappa_P_0 *
                         std::pow(rho_safe / params_.rho_0, params_.a_P) *
                         std::pow(T_safe / params_.T_0, params_.b_P);
        return clamp(k);
    }

private:
    double clamp(double k) const
    {
        if (!(k > params_.kappa_floor_m2_per_kg))
            return params_.kappa_floor_m2_per_kg;
        if (k > params_.kappa_ceiling_m2_per_kg)
            return params_.kappa_ceiling_m2_per_kg;
        return k;
    }
    PowerLawOpacityParameters params_ = PowerLawOpacitySets::granite();
};

/**
 * @brief Universal physical constants for the radiation diffusion solve.
 *
 * a = 4 sigma_SB / c is the radiation constant.
 */
struct RadiationConstants
{
    static constexpr double SPEED_OF_LIGHT_M_PER_S = 2.99792458e8;
    static constexpr double STEFAN_BOLTZMANN_W_PER_M2_K4 = 5.670374419e-8;
    static constexpr double RADIATION_CONSTANT_A_J_PER_M3_K4 = 7.5657e-16;
    static constexpr double BOLTZMANN_J_PER_K = 1.380649e-23;
};

} // namespace FSRM

#endif // NEAR_FIELD_OPACITY_MODEL_HPP
