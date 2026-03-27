#ifndef CO2_BRINE_SOLUBILITY_HPP
#define CO2_BRINE_SOLUBILITY_HPP

/**
 * @file CO2BrineSolubility.hpp
 * @brief CO2-brine mutual solubility models for CCS simulation
 *
 * Implements:
 * - Duan & Sun (2003): CO2 solubility in brine as f(P, T, salinity)
 * - Spycher, Pruess & Ennis-King (2003/2005): Mutual solubility
 *   (CO2 dissolution in brine + H2O vaporization into CO2 phase)
 * - Activity coefficient model for dissolved CO2
 * - Henry's law with fugacity correction as simpler fallback
 *
 * References:
 *   Duan, Z. & Sun, R. (2003). An improved model calculating CO2 solubility
 *   in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to
 *   2000 bar. Chem. Geol. 193(3-4), 257-271.
 *
 *   Spycher, N., Pruess, K. & Ennis-King, J. (2003). CO2-H2O mixtures in the
 *   geological sequestration of CO2. I. Assessment and calculation of mutual
 *   solubilities from 12 to 100°C and up to 600 bar. Geochim. Cosmochim.
 *   Acta 67(16), 3015-3031.
 *
 * All inputs in SI: P [Pa], T [K], salinity [mass fraction NaCl]
 */

#include <utility>

namespace FSRM {

class CO2BrineSolubility {
public:
    // =========================================================================
    // Primary API
    // =========================================================================

    /**
     * @brief CO2 solubility in brine (Duan-Sun 2003)
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param salinity NaCl mass fraction (e.g. 0.035 for 3.5 wt%)
     * @return CO2 mole fraction in aqueous phase (xCO2)
     */
    static double co2SolubilityInBrine(double P, double T, double salinity);

    /**
     * @brief Mutual solubility via Spycher-Pruess (2003/2005)
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param salinity NaCl mass fraction
     * @return {xCO2_in_brine, yH2O_in_CO2_phase}
     */
    static std::pair<double, double> mutualSolubility(double P, double T,
                                                       double salinity);

    /**
     * @brief Activity coefficient of dissolved CO2 in brine (Duan-Sun)
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param salinity NaCl mass fraction
     * @return Activity coefficient (dimensionless, >1)
     */
    static double activityCoefficientCO2(double T, double P, double salinity);

    /**
     * @brief Henry's law constant for CO2 in water
     * @param T Temperature [K]
     * @param salinity NaCl mass fraction
     * @return Henry's constant [Pa] (fugacity / mole_fraction)
     */
    static double henryConstant(double T, double salinity);

    /**
     * @brief CO2 fugacity coefficient from Peng-Robinson EOS
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Fugacity coefficient (dimensionless)
     */
    static double co2FugacityCoefficient(double P, double T);

    /**
     * @brief Salinity correction factor for CO2 solubility
     * @param T Temperature [K]
     * @param P Pressure [Pa]
     * @param m_NaCl NaCl molality [mol/kg water]
     * @return Salinity reduction factor (multiply pure-water solubility by this)
     */
    static double salinityCorrection(double T, double P, double m_NaCl);

    // =========================================================================
    // Utility conversions
    // =========================================================================

    /** @brief Convert mass fraction to molality [mol NaCl / kg water] */
    static double massFractionToMolality(double w_NaCl);

    /** @brief Convert molality to mass fraction */
    static double molalityToMassFraction(double m_NaCl);

private:
    // =========================================================================
    // Duan-Sun parametric model
    // =========================================================================

    /**
     * @brief Duan-Sun Par function for CO2 chemical potential in pure water
     *
     * ln(mCO2) = Par(T,P) + 2λ*mNaCl + ζ*mNaCl²
     * Par is a polynomial in T, P, T² etc. (Table 2 of Duan-Sun)
     */
    static double computeParCO2(double T, double P);

    /**
     * @brief Duan-Sun interaction parameters λ and ζ
     * @return {lambda, zeta}
     */
    static std::pair<double, double> interactionParams(double T, double P);

    // =========================================================================
    // Spycher-Pruess equilibrium constants
    // =========================================================================

    /** @brief Equilibrium constant for CO2 dissolution (ln K0) */
    static double K0_CO2(double T);

    /** @brief Equilibrium constant for H2O vaporization (ln K0) */
    static double K0_H2O(double T);

    /** @brief Molar volume of CO2 in aqueous phase [m³/mol] */
    static double partialMolarVolumeCO2(double T);

    /** @brief Molar volume of H2O [m³/mol] */
    static double molarVolumeH2O(double T);
};

} // namespace FSRM

#endif // CO2_BRINE_SOLUBILITY_HPP
