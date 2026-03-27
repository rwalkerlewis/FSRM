/**
 * @file CO2BrineSolubility.cpp
 * @brief CO2-brine mutual solubility implementation
 *
 * Duan-Sun (2003) model for CO2 solubility in brine and
 * Spycher-Pruess (2003/2005) for mutual solubility.
 */

#include "domain/reservoir/CO2BrineSolubility.hpp"
#include "domain/reservoir/CO2Properties.hpp"
#include <cmath>
#include <algorithm>

namespace FSRM {

// ============================================================================
// Constants
// ============================================================================

namespace {
    static constexpr double R = 8.31446;           // Universal gas constant [J/(mol·K)]
    static constexpr double Mw_NaCl = 0.058443;    // NaCl molar mass [kg/mol]
    static constexpr double Mw_H2O  = 0.018015;    // H2O molar mass [kg/mol]
    static constexpr double Mw_CO2  = 0.044010;    // CO2 molar mass [kg/mol]
    static constexpr double P_ref   = 1.0e5;       // Reference pressure [Pa] (1 bar)
} // anonymous namespace

// ============================================================================
// Utility conversions
// ============================================================================

double CO2BrineSolubility::massFractionToMolality(double w_NaCl) {
    // m = w / (Mw_NaCl * (1 - w))
    if (w_NaCl <= 0.0) return 0.0;
    if (w_NaCl >= 1.0) return 1e10;
    return w_NaCl / (Mw_NaCl * (1.0 - w_NaCl));
}

double CO2BrineSolubility::molalityToMassFraction(double m_NaCl) {
    if (m_NaCl <= 0.0) return 0.0;
    return m_NaCl * Mw_NaCl / (1.0 + m_NaCl * Mw_NaCl);
}

// ============================================================================
// Duan-Sun (2003) model
// ============================================================================

double CO2BrineSolubility::computeParCO2(double T, double P) {
    // Par(T,P) = c1 + c2*T + c3/T + c4*T² + c5/(630-T) + c6*P
    //          + c7*P*ln(T) + c8*P/T + c9*P/(630-T) + c10*P²/(630-T)²
    //          + c11*T*ln(P)
    //
    // Coefficients from Duan & Sun (2003), Table 2 (Model A for T < 573 K)

    double P_bar = P / 1e5;  // Convert Pa to bar (Duan-Sun uses bar)

    // Coefficients for Model A (273–533 K, 0–2000 bar)
    static constexpr double c1  =  28.9447706;
    static constexpr double c2  = -0.0354581768;
    static constexpr double c3  = -4770.67077;
    static constexpr double c4  =  1.02782768e-5;
    static constexpr double c5  =  33.8126098;
    static constexpr double c6  =  9.04037140e-3;
    static constexpr double c7  = -1.14934031e-3;
    static constexpr double c8  = -0.307405726;
    static constexpr double c9  = -0.0907301486;
    static constexpr double c10 =  9.32713393e-4;
    static constexpr double c11 =  0.0;  // Not used in Model A

    double T_limit = std::min(T, 620.0);  // Avoid singularity at 630 K
    double inv_630_T = 1.0 / (630.0 - T_limit);

    double par = c1
               + c2 * T
               + c3 / T
               + c4 * T * T
               + c5 * inv_630_T
               + c6 * P_bar
               + c7 * P_bar * std::log(T)
               + c8 * P_bar / T
               + c9 * P_bar * inv_630_T
               + c10 * P_bar * P_bar * inv_630_T * inv_630_T
               + c11 * T * std::log(P_bar + 1e-30);

    return par;
}

std::pair<double, double> CO2BrineSolubility::interactionParams(double T, double P) {
    // λ(T,P) and ζ(T,P) for Na-CO2 and Na-Cl-CO2 interactions
    // From Duan & Sun (2003), Table 3

    double P_bar = P / 1e5;

    // λ(Na-CO2) coefficients
    static constexpr double d1  = -0.411370585;
    static constexpr double d2  =  6.07632013e-4;
    static constexpr double d3  =  97.5347708;
    static constexpr double d4  =  0.0;
    static constexpr double d5  =  0.0;
    static constexpr double d6  =  0.0;
    static constexpr double d7  = -0.0237622469;
    static constexpr double d8  =  0.0170656236;
    static constexpr double d9  =  0.0;
    static constexpr double d10 =  1.41335834e-5;
    static constexpr double d11 =  0.0;

    double T_limit = std::min(T, 620.0);
    double inv_630_T = 1.0 / (630.0 - T_limit);

    double lambda = d1 + d2 * T + d3 / T + d4 * T * T + d5 * inv_630_T
                  + d6 * P_bar + d7 * P_bar * std::log(T) + d8 * P_bar / T
                  + d9 * P_bar * inv_630_T + d10 * P_bar * P_bar * inv_630_T * inv_630_T
                  + d11 * T * std::log(P_bar + 1e-30);

    // ζ(Na-Cl-CO2) coefficients
    static constexpr double e1  =  3.36389723e-4;
    static constexpr double e2  = -1.98298980e-5;
    static constexpr double e3  =  0.0;
    static constexpr double e4  =  0.0;
    static constexpr double e5  =  0.0;
    static constexpr double e6  =  0.0;
    static constexpr double e7  =  2.12220830e-3;
    static constexpr double e8  = -5.24873303e-3;
    static constexpr double e9  =  0.0;
    static constexpr double e10 =  0.0;
    static constexpr double e11 =  0.0;

    double zeta = e1 + e2 * T + e3 / T + e4 * T * T + e5 * inv_630_T
                + e6 * P_bar + e7 * P_bar * std::log(T) + e8 * P_bar / T
                + e9 * P_bar * inv_630_T + e10 * P_bar * P_bar * inv_630_T * inv_630_T
                + e11 * T * std::log(P_bar + 1e-30);

    return {lambda, zeta};
}

double CO2BrineSolubility::co2SolubilityInBrine(double P, double T, double salinity) {
    // Duan-Sun model:
    // ln(m_CO2) = Par(T,P) + 2*λ*m_NaCl + ζ*m_NaCl²
    // where m_CO2 is molality of CO2 in solution [mol/kg_water]

    double m_NaCl = massFractionToMolality(salinity);

    double par = computeParCO2(T, P);
    auto [lambda, zeta] = interactionParams(T, P);

    double ln_m_CO2 = par + 2.0 * lambda * m_NaCl + zeta * m_NaCl * m_NaCl;
    double m_CO2 = std::exp(ln_m_CO2);

    // Convert molality to mole fraction
    // x_CO2 = m_CO2 / (m_CO2 + 1/Mw_H2O + m_NaCl)
    double n_water = 1.0 / Mw_H2O;      // Moles of water per kg water
    double x_CO2 = m_CO2 / (m_CO2 + n_water + m_NaCl);

    return std::max(0.0, std::min(x_CO2, 0.1));  // Physical bounds
}

double CO2BrineSolubility::activityCoefficientCO2(double T, double P, double salinity) {
    // γ_CO2 = exp(2*λ*m_NaCl + ζ*m_NaCl²) from Duan-Sun
    double m_NaCl = massFractionToMolality(salinity);
    auto [lambda, zeta] = interactionParams(T, P);

    return std::exp(2.0 * lambda * m_NaCl + zeta * m_NaCl * m_NaCl);
}

double CO2BrineSolubility::salinityCorrection(double T, double P, double m_NaCl) {
    // Ratio: solubility in brine / solubility in pure water
    auto [lambda, zeta] = interactionParams(T, P);
    return std::exp(2.0 * lambda * m_NaCl + zeta * m_NaCl * m_NaCl);
}

// ============================================================================
// Spycher-Pruess (2003/2005) mutual solubility
// ============================================================================

double CO2BrineSolubility::K0_CO2(double T) {
    // Equilibrium constant for CO2 dissolution
    // ln(K0_CO2) fit from Spycher & Pruess (2005), Eq. 5
    double T_C = T - 273.15;
    return 1.189 + 1.304e-2 * T_C - 5.446e-5 * T_C * T_C;
}

double CO2BrineSolubility::K0_H2O(double T) {
    // Equilibrium constant for H2O vaporization
    // From Spycher & Pruess (2003), fit to experimental data
    double T_C = T - 273.15;

    // ln(K0_H2O) expressed via Antoine-like fit
    // K0_H2O ≈ exp(a + b/T + c/T²) in appropriate units
    double log10_K = -2.209 + 3.097e-2 * T_C - 1.098e-4 * T_C * T_C
                   + 2.048e-7 * T_C * T_C * T_C;
    return std::pow(10.0, log10_K);
}

double CO2BrineSolubility::partialMolarVolumeCO2(double T) {
    // Partial molar volume of CO2 in aqueous solution [m³/mol]
    // From Spycher & Pruess (2003)
    double T_C = T - 273.15;
    return (32.6 + 3.413e-2 * T_C) * 1e-6;  // cm³/mol → m³/mol
}

double CO2BrineSolubility::molarVolumeH2O(double T) {
    // Molar volume of liquid H2O [m³/mol]
    double T_C = T - 273.15;
    // Simple fit: ~18.0 cm³/mol at 25°C, increasing with temperature
    return (18.0 + 0.0094 * T_C + 3.0e-5 * T_C * T_C) * 1e-6;
}

std::pair<double, double> CO2BrineSolubility::mutualSolubility(double P, double T,
                                                                double salinity) {
    // Spycher-Pruess (2003/2005) model for mutual solubility
    //
    // At equilibrium:
    //   y_H2O * P * φ_H2O = K0_H2O * a_H2O * exp(V_H2O * (P - P_ref) / (R*T))
    //   y_CO2 * P * φ_CO2 = K0_CO2 * a_CO2 * x_CO2 * exp(V_CO2 * (P - P_ref) / (R*T))
    //
    // where φ are fugacity coefficients, a are activities, K0 are equilibrium constants

    double P_bar = P / 1e5;
    (void)P_bar;

    // Fugacity coefficients from Peng-Robinson
    double phi_CO2 = co2FugacityCoefficient(P, T);

    // H2O fugacity coefficient (simplified — close to 1 for low P)
    // At high pressure, use Poynting correction
    double V_H2O = molarVolumeH2O(T);
    double phi_H2O_correction = std::exp(V_H2O * (P - P_ref) / (R * T));

    // K0 values
    double K_CO2 = K0_CO2(T);
    double K_H2O = K0_H2O(T);

    // Partial molar volume corrections
    double V_CO2 = partialMolarVolumeCO2(T);
    double poynting_CO2 = std::exp(V_CO2 * (P - P_ref) / (R * T));
    double poynting_H2O = phi_H2O_correction;

    // Activity of water in brine (Raoult's law approximation)
    double m_NaCl = massFractionToMolality(salinity);
    double x_NaCl_aq = m_NaCl / (m_NaCl + 1.0 / Mw_H2O);
    double a_H2O = 1.0 - x_NaCl_aq;  // Activity coefficient ≈ 1

    // Activity coefficient for CO2 in brine
    double gamma_CO2 = activityCoefficientCO2(T, P, salinity);

    // Solve for y_H2O (H2O mole fraction in CO2 phase)
    // y_H2O = K_H2O * a_H2O * poynting_H2O / (P * phi_H2O)
    // Approximate φ_H2O ≈ 1 for simplicity (valid at moderate P)
    double y_H2O_numer = K_H2O * a_H2O * poynting_H2O;
    double y_H2O_denom = P;
    double y_H2O = y_H2O_numer / y_H2O_denom;
    y_H2O = std::max(0.0, std::min(y_H2O, 0.1));  // Physical bound

    // y_CO2 = 1 - y_H2O (binary system)
    double y_CO2 = 1.0 - y_H2O;

    // Solve for x_CO2 (CO2 mole fraction in aqueous phase)
    // x_CO2 = y_CO2 * P * phi_CO2 / (K_CO2 * gamma_CO2 * poynting_CO2)
    double x_CO2_numer = y_CO2 * P * phi_CO2;
    double x_CO2_denom = K_CO2 * gamma_CO2 * poynting_CO2;

    double x_CO2 = (x_CO2_denom > 1e-30) ? x_CO2_numer / x_CO2_denom : 0.0;
    x_CO2 = std::max(0.0, std::min(x_CO2, 0.1));

    return {x_CO2, y_H2O};
}

// ============================================================================
// Henry's law fallback
// ============================================================================

double CO2BrineSolubility::henryConstant(double T, double salinity) {
    // Henry's constant for CO2 in water [Pa]
    // H = P_sat_ref * exp(A + B/T + C*ln(T))
    // Fit to experimental data (Carroll et al. 1991)

    double T_C = T - 273.15;
    (void)T_C;

    // ln(H/Pa) = A + B/T + C*ln(T) + D*T
    static constexpr double A_h = -6.8346;
    static constexpr double B_h =  1.2817e4;
    static constexpr double C_h = -3.7668e6;
    static constexpr double D_h =  2.997e8;

    double lnH = A_h + B_h / T + C_h / (T * T) + D_h / (T * T * T);
    double H_pure = std::exp(lnH) * 1e5;  // Convert from bar to Pa

    // Salinity correction (Sechenov coefficient)
    double m_NaCl = massFractionToMolality(salinity);
    double K_s = 0.11;  // Sechenov coefficient [kg/mol] (approximate for CO2-NaCl)
    double H_brine = H_pure * std::exp(K_s * m_NaCl);

    return H_brine;
}

// ============================================================================
// CO2 fugacity coefficient
// ============================================================================

double CO2BrineSolubility::co2FugacityCoefficient(double P, double T) {
    // Peng-Robinson EOS fugacity coefficient for pure CO2
    //
    // ln(φ) = (Z - 1) - ln(Z - B) - A/(2√2·B) * ln[(Z + (1+√2)B)/(Z + (1-√2)B)]

    static constexpr double omega_CO2 = 0.22394;
    static constexpr double Tc_CO2 = CO2Properties::Tc;
    static constexpr double Pc_CO2 = CO2Properties::Pc;

    double kappa = 0.37464 + 1.54226 * omega_CO2 - 0.26992 * omega_CO2 * omega_CO2;
    double alpha = std::pow(1.0 + kappa * (1.0 - std::sqrt(T / Tc_CO2)), 2);

    double a_PR = 0.45724 * R * R * Tc_CO2 * Tc_CO2 / Pc_CO2 * alpha;
    double b_PR = 0.07780 * R * Tc_CO2 / Pc_CO2;

    double A = a_PR * P / (R * R * T * T);
    double B = b_PR * P / (R * T);

    // Solve cubic for Z
    double c2 = -(1.0 - B);
    double c1 = A - 3.0 * B * B - 2.0 * B;
    double c0 = -(A * B - B * B - B * B * B);

    double p = c1 - c2 * c2 / 3.0;
    double q = c0 - c1 * c2 / 3.0 + 2.0 * c2 * c2 * c2 / 27.0;
    double D = q * q / 4.0 + p * p * p / 27.0;

    double Z;
    if (D > 0) {
        double sqrt_D = std::sqrt(D);
        double u = std::cbrt(-q / 2.0 + sqrt_D);
        double v = std::cbrt(-q / 2.0 - sqrt_D);
        Z = u + v - c2 / 3.0;
    } else {
        double theta = std::acos(-q / 2.0 * std::sqrt(-27.0 / (p * p * p)));
        double r_val = 2.0 * std::sqrt(-p / 3.0);
        double Z1 = r_val * std::cos(theta / 3.0) - c2 / 3.0;
        double Z2 = r_val * std::cos((theta + 2.0 * M_PI) / 3.0) - c2 / 3.0;
        double Z3 = r_val * std::cos((theta + 4.0 * M_PI) / 3.0) - c2 / 3.0;
        Z = std::max({Z1, Z2, Z3});  // Gas-phase root
    }
    Z = std::max(Z, B + 1e-10);

    // Fugacity coefficient
    double sqrt2 = std::sqrt(2.0);
    double arg_num = Z + (1.0 + sqrt2) * B;
    double arg_den = Z + (1.0 - sqrt2) * B;

    double ln_phi;
    if (arg_num > 0.0 && arg_den > 0.0) {
        ln_phi = (Z - 1.0) - std::log(Z - B)
               - A / (2.0 * sqrt2 * B) * std::log(arg_num / arg_den);
    } else {
        ln_phi = Z - 1.0 - std::log(Z);  // Simplified fallback
    }

    return std::exp(ln_phi);
}

} // namespace FSRM
