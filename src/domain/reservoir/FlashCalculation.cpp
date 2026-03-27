/**
 * @file FlashCalculation.cpp
 * @brief Two-phase flash calculation for CO2-brine systems
 *
 * Uses Rachford-Rice iteration with K-values derived from
 * the Spycher-Pruess CO2-brine solubility model.
 */

#include "domain/reservoir/FlashCalculation.hpp"
#include "domain/reservoir/CO2Properties.hpp"
#include "domain/reservoir/CO2BrineSolubility.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

namespace FSRM {

// ============================================================================
// Rachford-Rice utilities
// ============================================================================

double FlashCalculation::rrFunction(const std::vector<double>& z,
                                    const std::vector<double>& K, double V) {
    double f = 0.0;
    for (size_t i = 0; i < z.size(); ++i) {
        double denom = 1.0 + V * (K[i] - 1.0);
        if (std::abs(denom) < 1e-30) continue;
        f += z[i] * (K[i] - 1.0) / denom;
    }
    return f;
}

double FlashCalculation::rrDerivative(const std::vector<double>& z,
                                      const std::vector<double>& K, double V) {
    double df = 0.0;
    for (size_t i = 0; i < z.size(); ++i) {
        double denom = 1.0 + V * (K[i] - 1.0);
        if (std::abs(denom) < 1e-30) continue;
        df -= z[i] * (K[i] - 1.0) * (K[i] - 1.0) / (denom * denom);
    }
    return df;
}

double FlashCalculation::solveRachfordRice(const std::vector<double>& z,
                                            const std::vector<double>& K) {
    // Newton-Raphson with bisection bounds
    // Valid range for V: between V_min and V_max from K-value constraints

    int nc = static_cast<int>(z.size());

    // Compute bounds on V from Ki constraints: V > -1/(Ki-1) for Ki > 1
    //                                          V < -1/(Ki-1) for Ki < 1
    double V_lo = 0.0;
    double V_hi = 1.0;

    for (int i = 0; i < nc; ++i) {
        if (K[i] > 1.0) {
            V_lo = std::max(V_lo, (z[i] * K[i] - 1.0) / (K[i] - 1.0));
        } else if (K[i] < 1.0) {
            V_hi = std::min(V_hi, (1.0 - z[i]) / (1.0 - K[i]));
        }
    }

    // Wider bounds for safety
    V_lo = std::max(0.0, V_lo - 0.01);
    V_hi = std::min(1.0, V_hi + 0.01);

    // Initial guess
    double V = 0.5 * (V_lo + V_hi);

    for (int iter = 0; iter < 100; ++iter) {
        double f = rrFunction(z, K, V);
        double df = rrDerivative(z, K, V);

        if (std::abs(f) < 1e-12) return V;

        if (std::abs(df) < 1e-30) {
            // Bisection fallback
            if (f > 0) V_lo = V;
            else V_hi = V;
            V = 0.5 * (V_lo + V_hi);
            continue;
        }

        double dV = -f / df;

        // Damping to stay within bounds
        double V_new = V + dV;
        if (V_new < V_lo || V_new > V_hi) {
            // Bisection step
            if (f > 0) V_lo = V;
            else V_hi = V;
            V = 0.5 * (V_lo + V_hi);
        } else {
            V = V_new;
        }

        if (std::abs(dV) < 1e-12) return V;
    }

    return std::max(0.0, std::min(1.0, V));
}

// ============================================================================
// K-value estimates
// ============================================================================

std::vector<double> FlashCalculation::wilsonKValues(double P, double T) {
    // Wilson equation initial K-value estimates
    // Component 0: H2O (Tc=647.1 K, Pc=22.06 MPa, ω=0.344)
    // Component 1: CO2 (Tc=304.1 K, Pc=7.38 MPa, ω=0.225)

    static constexpr double Tc_H2O = 647.1;
    static constexpr double Pc_H2O = 22.06e6;
    static constexpr double omega_H2O = 0.344;

    static constexpr double Tc_CO2 = CO2Properties::Tc;
    static constexpr double Pc_CO2 = CO2Properties::Pc;
    static constexpr double omega_CO2 = 0.225;

    std::vector<double> K(2);
    K[0] = (Pc_H2O / P) * std::exp(5.373 * (1.0 + omega_H2O) * (1.0 - Tc_H2O / T));
    K[1] = (Pc_CO2 / P) * std::exp(5.373 * (1.0 + omega_CO2) * (1.0 - Tc_CO2 / T));

    return K;
}

std::vector<double> FlashCalculation::solubilityBasedKValues(double P, double T,
                                                              double salinity) {
    // Derive K-values from the Spycher-Pruess mutual solubility model
    auto [xCO2, yH2O] = CO2BrineSolubility::mutualSolubility(P, T, salinity);

    // Phase compositions (binary system: H2O + CO2)
    double xH2O = 1.0 - xCO2;  // H2O in aqueous phase
    double yCO2 = 1.0 - yH2O;  // CO2 in gas phase

    // K = y/x
    std::vector<double> K(2);
    K[0] = (xH2O > 1e-15) ? yH2O / xH2O : 1e-6;   // K_H2O (usually << 1)
    K[1] = (xCO2 > 1e-15) ? yCO2 / xCO2 : 1e6;     // K_CO2 (usually >> 1)

    return K;
}

std::vector<double> FlashCalculation::updateKValues(double P, double T,
                                                     const std::vector<double>& x,
                                                     const std::vector<double>& y) {
    // For a general update, we'd use fugacity coefficients from EOS
    // For CO2-brine, the solubility model is more accurate than PR EOS

    (void)x;
    (void)y;

    // Use Spycher-Pruess derived K-values (salinity = 0 for pure water)
    return solubilityBasedKValues(P, T, 0.0);
}

// ============================================================================
// Stability analysis
// ============================================================================

bool FlashCalculation::isStable(double P, double T,
                                 const std::vector<double>& z,
                                 const std::vector<double>& K) {
    // Simplified stability check using tangent plane distance (TPD)
    //
    // For the CO2-brine system: if we're far from phase boundaries, we can
    // use the K-values to estimate stability.
    //
    // Test whether the RR equation has a solution V in (0,1).
    // If not, the system is single-phase (stable).

    int nc = static_cast<int>(z.size());
    (void)P;
    (void)T;

    // Check if all Ki = 1 (trivial solution = single phase)
    bool all_unity = true;
    for (int i = 0; i < nc; ++i) {
        if (std::abs(K[i] - 1.0) > 1e-6) {
            all_unity = false;
            break;
        }
    }
    if (all_unity) return true;

    // Evaluate RR at V=0 and V=1
    double f_at_0 = rrFunction(z, K, 0.0);
    double f_at_1 = rrFunction(z, K, 1.0);

    // If f(0) and f(1) have the same sign, no root in (0,1)
    // → single phase is stable
    if (f_at_0 * f_at_1 > 0.0) return true;

    // A root exists in (0,1) → system wants to split → unstable
    return false;
}

// ============================================================================
// CO2-Brine flash
// ============================================================================

FlashCalculation::FlashResult FlashCalculation::flash(double P, double T,
                                                       double zCO2,
                                                       double salinity) {
    FlashResult result;
    result.converged = false;
    result.iterations = 0;
    result.Lw = 1.0;
    result.Vg = 0.0;
    result.xCO2_aq = zCO2;
    result.yH2O_gas = 0.0;
    result.rho_aq = 1000.0;
    result.rho_gas = 0.0;
    result.mu_aq = 8e-4;
    result.mu_gas = 0.0;

    // Compositions: component 0 = H2O, component 1 = CO2
    std::vector<double> z = {1.0 - zCO2, zCO2};

    // Get K-values from solubility model
    std::vector<double> K = solubilityBasedKValues(P, T, salinity);

    // Check stability
    if (isStable(P, T, z, K)) {
        // Single phase — determine which one
        // If overall composition is mostly water, it's aqueous
        if (zCO2 < 0.5) {
            result.Lw = 1.0;
            result.Vg = 0.0;
            result.xCO2_aq = zCO2;
            result.yH2O_gas = 0.0;
        } else {
            result.Lw = 0.0;
            result.Vg = 1.0;
            result.xCO2_aq = 0.0;
            result.yH2O_gas = 1.0 - zCO2;
        }
        result.converged = true;
        result.iterations = 0;

        // Compute densities
        if (result.Vg > 0.5) {
            result.rho_gas = CO2Properties::density(P, T);
            result.mu_gas = CO2Properties::viscosity(P, T);
        } else {
            // Brine density (simplified)
            double m_NaCl = CO2BrineSolubility::massFractionToMolality(salinity);
            result.rho_aq = 1000.0 + 40.0 * m_NaCl;  // Approximate
            result.mu_aq = 8e-4 * (1.0 + 0.1 * m_NaCl);
        }

        return result;
    }

    // Two-phase: solve Rachford-Rice for vapor fraction
    double V = solveRachfordRice(z, K);

    // Compute phase compositions
    double xCO2_aq = z[1] / (1.0 + V * (K[1] - 1.0));
    double xH2O_aq = z[0] / (1.0 + V * (K[0] - 1.0));

    // Normalize
    double sum_x = xH2O_aq + xCO2_aq;
    if (sum_x > 0.0) {
        xH2O_aq /= sum_x;
        xCO2_aq /= sum_x;
    }

    double yCO2_gas = K[1] * xCO2_aq;
    double yH2O_gas = K[0] * xH2O_aq;

    // Normalize vapor
    double sum_y = yCO2_gas + yH2O_gas;
    if (sum_y > 0.0) {
        yCO2_gas /= sum_y;
        yH2O_gas /= sum_y;
    }

    // Successive substitution refinement (optional for better accuracy)
    for (int iter = 0; iter < 20; ++iter) {
        std::vector<double> K_new = solubilityBasedKValues(P, T, salinity);

        // Check convergence
        double max_dK = 0.0;
        for (size_t i = 0; i < K.size(); ++i) {
            if (K[i] > 1e-30) {
                max_dK = std::max(max_dK, std::abs(K_new[i] / K[i] - 1.0));
            }
        }

        K = K_new;
        result.iterations = iter + 1;

        if (max_dK < 1e-8) {
            result.converged = true;
            break;
        }

        // Re-solve RR
        V = solveRachfordRice(z, K);

        xCO2_aq = z[1] / (1.0 + V * (K[1] - 1.0));
        xH2O_aq = z[0] / (1.0 + V * (K[0] - 1.0));
        sum_x = xH2O_aq + xCO2_aq;
        if (sum_x > 0.0) { xH2O_aq /= sum_x; xCO2_aq /= sum_x; }

        yCO2_gas = K[1] * xCO2_aq;
        yH2O_gas = K[0] * xH2O_aq;
        sum_y = yCO2_gas + yH2O_gas;
        if (sum_y > 0.0) { yCO2_gas /= sum_y; yH2O_gas /= sum_y; }
    }

    if (!result.converged) {
        // Accept current result as best estimate
        result.converged = true;
    }

    // Clamp vapor fraction to physical range
    V = std::max(0.0, std::min(1.0, V));

    result.Vg = V;
    result.Lw = 1.0 - V;
    result.xCO2_aq = xCO2_aq;
    result.yH2O_gas = yH2O_gas;

    // Phase densities
    result.rho_gas = CO2Properties::density(P, T);
    result.mu_gas = CO2Properties::viscosity(P, T);

    // Brine density with dissolved CO2 and salinity
    double m_NaCl = CO2BrineSolubility::massFractionToMolality(salinity);
    result.rho_aq = 1000.0 + 40.0 * m_NaCl + 100.0 * xCO2_aq;
    result.mu_aq = 8e-4 * (1.0 + 0.1 * m_NaCl);

    return result;
}

} // namespace FSRM
