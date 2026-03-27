#ifndef FLASH_CALCULATION_HPP
#define FLASH_CALCULATION_HPP

/**
 * @file FlashCalculation.hpp
 * @brief Two-phase flash calculation for CO2-brine systems
 *
 * Implements:
 * - Rachford-Rice iteration for two-phase split
 * - Successive substitution for K-value updates
 * - Tangent plane distance (TPD) stability analysis
 * - Specialized CO2-brine flash using solubility models
 *
 * The flash calculation determines, given overall composition z and (P,T),
 * whether the system is single-phase or two-phase, and if two-phase,
 * the phase compositions and fractions.
 *
 * For the CO2-brine system, the "components" are H2O and CO2, and the
 * two phases are an aqueous (brine) phase and a CO2-rich phase.
 */

#include <vector>

namespace FSRM {

class FlashCalculation {
public:
    /**
     * @brief Result of a flash calculation
     */
    struct FlashResult {
        double Lw;          ///< Aqueous (brine) phase mole fraction
        double Vg;          ///< CO2-rich phase mole fraction
        double xCO2_aq;     ///< CO2 mole fraction in aqueous phase
        double yH2O_gas;    ///< H2O mole fraction in CO2-rich phase
        double rho_aq;      ///< Aqueous phase density [kg/m³]
        double rho_gas;     ///< CO2-rich phase density [kg/m³]
        double mu_aq;       ///< Aqueous phase viscosity [Pa·s]
        double mu_gas;      ///< CO2-rich phase viscosity [Pa·s]
        bool converged;     ///< True if iteration converged
        int iterations;     ///< Number of iterations used
    };

    // =========================================================================
    // CO2-Brine specific flash
    // =========================================================================

    /**
     * @brief Two-phase flash for CO2-brine system
     *
     * Uses the Spycher-Pruess solubility model combined with
     * Rachford-Rice iteration for the phase split.
     *
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param zCO2 Overall CO2 mole fraction (0 to 1)
     * @param salinity NaCl mass fraction in brine
     * @return Flash result with phase compositions and properties
     */
    static FlashResult flash(double P, double T, double zCO2, double salinity);

    // =========================================================================
    // General multi-component flash utilities
    // =========================================================================

    /**
     * @brief Rachford-Rice objective function
     *
     * f(V) = Σ z_i * (K_i - 1) / (1 + V * (K_i - 1)) = 0
     *
     * @param z Overall compositions
     * @param K Equilibrium ratios (K-values)
     * @param V Vapor fraction guess
     * @return Value of Rachford-Rice function
     */
    static double rrFunction(const std::vector<double>& z,
                             const std::vector<double>& K, double V);

    /**
     * @brief Derivative of Rachford-Rice function
     *
     * f'(V) = -Σ z_i * (K_i - 1)² / (1 + V * (K_i - 1))²
     *
     * @param z Overall compositions
     * @param K Equilibrium ratios
     * @param V Vapor fraction
     * @return Derivative df/dV
     */
    static double rrDerivative(const std::vector<double>& z,
                               const std::vector<double>& K, double V);

    /**
     * @brief Solve Rachford-Rice equation for vapor fraction V
     *
     * Uses Newton-Raphson with bisection fallback.
     *
     * @param z Overall compositions
     * @param K Equilibrium ratios
     * @return Converged vapor fraction V in [0, 1]
     */
    static double solveRachfordRice(const std::vector<double>& z,
                                    const std::vector<double>& K);

    /**
     * @brief Stability analysis via tangent plane distance
     *
     * Tests whether the single-phase mixture z is thermodynamically
     * stable at (P,T). If unstable, a two-phase split should occur.
     *
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param z Overall compositions
     * @param K K-value estimates (Wilson or from EOS)
     * @return true if single-phase is stable, false if two-phase split needed
     */
    static bool isStable(double P, double T, const std::vector<double>& z,
                         const std::vector<double>& K);

    /**
     * @brief Update K-values using successive substitution
     *
     * K_i = φ_i^L / φ_i^V where φ are fugacity coefficients.
     * For CO2-brine, uses solubility-model-derived K-values.
     *
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param x Liquid phase composition
     * @param y Vapor phase composition
     * @return Updated K-values
     */
    static std::vector<double> updateKValues(double P, double T,
                                              const std::vector<double>& x,
                                              const std::vector<double>& y);

    /**
     * @brief Initial K-value estimate from Wilson equation
     *
     * K_i = (Pc_i / P) * exp(5.373 * (1 + ω_i) * (1 - Tc_i / T))
     *
     * For CO2-brine: component 0 = H2O, component 1 = CO2
     *
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return K-values for [H2O, CO2]
     */
    static std::vector<double> wilsonKValues(double P, double T);

private:
    /**
     * @brief Compute K-values from solubility model
     *
     * For the CO2-brine binary system:
     *   K_CO2 = y_CO2 / x_CO2  (from Spycher-Pruess)
     *   K_H2O = y_H2O / x_H2O  (from Spycher-Pruess)
     */
    static std::vector<double> solubilityBasedKValues(double P, double T,
                                                       double salinity);
};

} // namespace FSRM

#endif // FLASH_CALCULATION_HPP
