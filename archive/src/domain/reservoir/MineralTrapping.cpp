/**
 * @file MineralTrapping.cpp
 * @brief Mineral trapping / water–rock interaction kinetics for CCS (simplified TST law)
 */

#include "domain/reservoir/MineralTrapping.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

namespace FSRM {

namespace {

constexpr double kMinTemperature_K = 1.0;
constexpr double kPHMin = 0.0;
constexpr double kPHMax = 14.0;
constexpr double kPhiMin = 0.001;
constexpr double kPhiMax = 0.5;
constexpr double kPorosityEps = 1.0e-9;
constexpr double kPermPhiUpper = 1.0 - 1.0e-6;

double clamp(double v, double lo, double hi) {
    return std::fmax(lo, std::fmin(hi, v));
}

bool isDissolutionReaction(const std::string& name) {
    std::string lower;
    lower.reserve(name.size());
    for (unsigned char c : name) {
        lower.push_back(static_cast<char>(std::tolower(c)));
    }
    return lower.find("dissolution") != std::string::npos;
}

} // namespace

void MineralTrapping::addReaction(const MineralReaction& rxn) {
    reactions_.push_back(rxn);
}

std::vector<double> MineralTrapping::computeRates(
    double T, double pH, const std::vector<double>& mineral_fractions) const {
    const std::size_t n = reactions_.size();
    std::vector<double> rates(n, 0.0);

    if (n == 0) {
        return rates;
    }
    if (mineral_fractions.size() != n) {
        return rates;
    }
    if (T < kMinTemperature_K || !std::isfinite(T)) {
        return rates;
    }

    pH = clamp(pH, kPHMin, kPHMax);
    const double Q = std::pow(10.0, -pH * 2.0);
    if (!std::isfinite(Q) || Q < 0.0) {
        return rates;
    }

    const double arrhenius_den = R_gas * T;
    if (arrhenius_den <= 0.0 || !std::isfinite(arrhenius_den)) {
        return rates;
    }

    for (std::size_t i = 0; i < n; ++i) {
        const MineralReaction& rx = reactions_[i];
        if (rx.reactive_surface_A < 0.0 || rx.rate_constant_k0 < 0.0 ||
            rx.equilibrium_constant_Keq <= 0.0 || !std::isfinite(rx.equilibrium_constant_Keq)) {
            rates[i] = 0.0;
            continue;
        }

        const double exp_term = std::exp(-rx.activation_energy_Ea / arrhenius_den);
        if (!std::isfinite(exp_term)) {
            rates[i] = 0.0;
            continue;
        }

        const double Q_over_K = Q / rx.equilibrium_constant_Keq;
        double driving = 1.0 - Q_over_K;
        if (!std::isfinite(driving)) {
            driving = 0.0;
        }
        // Limit extreme departure from equilibrium for numerical stability
        driving = clamp(driving, -1.0, 1.0);

        double r = rx.reactive_surface_A * rx.rate_constant_k0 * exp_term * driving;
        if (!std::isfinite(r)) {
            r = 0.0;
        }
        rates[i] = r;
    }

    return rates;
}

void MineralTrapping::updateMineralFractions(std::vector<double>& fractions, double dt, double T,
                                             double pH) const {
    const std::size_t n = reactions_.size();
    if (n == 0 || fractions.size() != n) {
        return;
    }
    if (dt < 0.0 || !std::isfinite(dt)) {
        return;
    }

    const std::vector<double> rates = computeRates(T, pH, fractions);

    for (std::size_t i = 0; i < n; ++i) {
        const MineralReaction& rx = reactions_[i];
        if (rx.molar_volume < 0.0 || !std::isfinite(rx.molar_volume)) {
            continue;
        }

        const int sign = isDissolutionReaction(rx.name) ? -1 : 1;
        double df = sign * rates[i] * rx.molar_volume * dt;
        if (!std::isfinite(df)) {
            continue;
        }

        fractions[i] += df;
        fractions[i] = clamp(fractions[i], 0.0, 1.0);
    }
}

double MineralTrapping::updatePorosity(double phi0, double delta_mineral_volume) {
    if (!std::isfinite(phi0) || !std::isfinite(delta_mineral_volume)) {
        return clamp(phi0, kPhiMin, kPhiMax);
    }
    double phi_new = phi0 - delta_mineral_volume;
    return clamp(phi_new, kPhiMin, kPhiMax);
}

double MineralTrapping::updatePermeability(double k0, double phi0, double phi_new) {
    if (!std::isfinite(k0) || !std::isfinite(phi0) || !std::isfinite(phi_new)) {
        return k0;
    }
    if (k0 < 0.0) {
        return k0;
    }

    double p0 = phi0;
    double p1 = phi_new;

    if (p0 <= kPorosityEps || p0 >= kPermPhiUpper) {
        return k0;
    }
    p1 = clamp(p1, kPorosityEps, kPermPhiUpper);

    const double one_m_p0 = 1.0 - p0;
    const double one_m_p1 = 1.0 - p1;
    if (one_m_p0 <= kPorosityEps || one_m_p1 <= kPorosityEps) {
        return k0;
    }

    const double ratio_phi = p1 / p0;
    const double ratio_one = one_m_p0 / one_m_p1;
    const double k_new = k0 * std::pow(ratio_phi, 3.0) * std::pow(ratio_one, 2.0);

    if (!std::isfinite(k_new) || k_new < 0.0) {
        return k0;
    }
    return k_new;
}

std::vector<MineralReaction> MineralTrapping::defaultCCSReactions() {
    std::vector<MineralReaction> out;
    out.reserve(5);

    out.push_back(MineralReaction{"Calcite precipitation",
                                  1.6e-9,
                                  41800.0,
                                  500.0,
                                  std::pow(10.0, -8.48),
                                  3.69e-5});
    out.push_back(MineralReaction{"Dawsonite precipitation",
                                  1e-10,
                                  62760.0,
                                  100.0,
                                  std::pow(10.0, -3.1),
                                  5.89e-5});
    out.push_back(MineralReaction{"Ankerite precipitation",
                                  1e-11,
                                  56000.0,
                                  50.0,
                                  std::pow(10.0, -7.2),
                                  6.63e-5});
    out.push_back(MineralReaction{"K-feldspar dissolution",
                                  1e-12,
                                  67830.0,
                                  200.0,
                                  1e20,
                                  1.09e-4});
    out.push_back(MineralReaction{"Chlorite dissolution",
                                  1e-13,
                                  88000.0,
                                  100.0,
                                  1e15,
                                  2.1e-4});
    return out;
}

} // namespace FSRM
