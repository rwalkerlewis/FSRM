#ifndef CO2_PROPERTIES_HPP
#define CO2_PROPERTIES_HPP

/**
 * @file CO2Properties.hpp
 * @brief CO2 equation of state and transport properties for CCS simulation
 *
 * Implements:
 * - Span & Wagner (1996) EOS for pure CO2 density and enthalpy
 * - Fenghour et al. (1998) viscosity correlation
 * - Vesovic et al. (1990) thermal conductivity (simplified)
 * - Peng-Robinson EOS as a faster alternative
 * - Phase boundary detection (liquid/gas/supercritical)
 *
 * Reference: Span, R. & Wagner, W. (1996). A New Equation of State for
 * Carbon Dioxide Covering the Fluid Region from the Triple-Point Temperature
 * to 1100 K at Pressures up to 800 MPa. J. Phys. Chem. Ref. Data, 25(6).
 *
 * All functions take (P [Pa], T [K]) and return properties in SI units.
 */

namespace FSRM {

class CO2Properties {
public:
    // =========================================================================
    // Critical constants for CO2
    // =========================================================================
    static constexpr double Tc   = 304.1282;    ///< Critical temperature [K]
    static constexpr double Pc   = 7.3773e6;    ///< Critical pressure [Pa]
    static constexpr double rhoc = 467.6;       ///< Critical density [kg/m³]
    static constexpr double Mw   = 0.0440098;   ///< Molar mass [kg/mol]
    static constexpr double R_specific = 188.9241; ///< Specific gas constant R/Mw [J/(kg·K)]
    static constexpr double R_universal = 8.31446; ///< Universal gas constant [J/(mol·K)]

    // Triple point
    static constexpr double T_triple = 216.592;  ///< Triple-point temperature [K]
    static constexpr double P_triple = 5.1795e5; ///< Triple-point pressure [Pa]

    // =========================================================================
    // Thermodynamic properties from Span-Wagner EOS
    // =========================================================================

    /**
     * @brief CO2 density from Span-Wagner EOS
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Density [kg/m³]
     */
    static double density(double P, double T);

    /**
     * @brief CO2 density from Peng-Robinson EOS (faster alternative)
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Density [kg/m³]
     */
    static double densityPR(double P, double T);

    /**
     * @brief CO2 specific enthalpy from Span-Wagner EOS
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Specific enthalpy [J/kg]
     */
    static double enthalpy(double P, double T);

    /**
     * @brief CO2 isothermal compressibility
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Compressibility [1/Pa]
     */
    static double compressibility(double P, double T);

    // =========================================================================
    // Transport properties
    // =========================================================================

    /**
     * @brief CO2 dynamic viscosity (Fenghour et al. 1998)
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Dynamic viscosity [Pa·s]
     */
    static double viscosity(double P, double T);

    /**
     * @brief CO2 thermal conductivity (Vesovic et al. 1990, simplified)
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @return Thermal conductivity [W/(m·K)]
     */
    static double thermalConductivity(double P, double T);

    // =========================================================================
    // Phase identification
    // =========================================================================

    /** @brief True if T > Tc and P > Pc */
    static bool isSupercritical(double P, double T);

    /** @brief True if below critical point and above saturation curve */
    static bool isLiquid(double P, double T);

    /** @brief True if below critical point and below saturation curve */
    static bool isGas(double P, double T);

    /**
     * @brief Saturation (vapor) pressure from Wagner equation
     * @param T Temperature [K], must be T_triple <= T <= Tc
     * @return Saturation pressure [Pa]
     */
    static double saturationPressure(double T);

private:
    // =========================================================================
    // Span-Wagner Helmholtz free energy (reduced form)
    // α(δ,τ) = α°(δ,τ) + αr(δ,τ)
    // where δ = ρ/ρc, τ = Tc/T
    // =========================================================================

    /// Ideal-gas part of dimensionless Helmholtz free energy
    static double alphaIdeal(double delta, double tau);

    /// Derivative dα°/dτ (for enthalpy)
    static double dalphaI_dtau(double delta, double tau);

    /// Residual part of dimensionless Helmholtz free energy
    static double alphaResidual(double delta, double tau);

    /// First derivative ∂αr/∂δ (for pressure)
    static double dalphaR_ddelta(double delta, double tau);

    /// Second derivative ∂²αr/∂δ² (for compressibility)
    static double d2alphaR_ddelta2(double delta, double tau);

    /// First derivative ∂αr/∂τ (for enthalpy)
    static double dalphaR_dtau(double delta, double tau);

    /// Cross derivative ∂²αr/(∂δ∂τ) (for speed of sound etc.)
    static double d2alphaR_ddeltadtau(double delta, double tau);

    // =========================================================================
    // Iterative density solver
    // =========================================================================

    /**
     * @brief Newton iteration to find density satisfying P = ρRT(1 + δ·∂αr/∂δ)
     * @param P Target pressure [Pa]
     * @param T Temperature [K]
     * @param rho_guess Initial density guess [kg/m³]
     * @return Converged density [kg/m³]
     */
    static double densityFromPressure(double P, double T, double rho_guess);

    /// Initial density guess based on phase region
    static double initialDensityGuess(double P, double T);
};

} // namespace FSRM

#endif // CO2_PROPERTIES_HPP
