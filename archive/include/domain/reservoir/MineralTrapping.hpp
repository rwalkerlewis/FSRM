#ifndef MINERAL_TRAPPING_HPP
#define MINERAL_TRAPPING_HPP

#include <string>
#include <vector>
#include <cmath>

namespace FSRM {

struct MineralReaction {
    std::string name;
    double rate_constant_k0;       // mol/(m²·s)
    double activation_energy_Ea;   // J/mol
    double reactive_surface_A;     // m²/m³_rock
    double equilibrium_constant_Keq; // dimensionless
    double molar_volume;             // m³/mol
};

class MineralTrapping {
public:
    void addReaction(const MineralReaction& rxn);

    /**
     * @brief TST-style kinetic rate per reaction [mol/(m³_bulk·s)]
     * r = A * k0 * exp(-Ea/(R*T)) * (1 - Q/Keq), Q ≈ 10^(-2·pH)
     * @param mineral_fractions volume fractions (one per reaction); must match reaction count
     */
    std::vector<double> computeRates(double T, double pH,
                                     const std::vector<double>& mineral_fractions) const;

    /**
     * @brief Explicit Euler update of mineral volume fractions [dimensionless]
     * df_i = sign_i * r_i * Vm_i * dt (sign negative for reactions named with "dissolution")
     */
    void updateMineralFractions(std::vector<double>& fractions, double dt, double T,
                                double pH) const;

    /** Porosity after mineral volume change; clamped to [0.001, 0.5] */
    static double updatePorosity(double phi0, double delta_mineral_volume);

    /** Kozeny–Carman permeability ratio; k0, phi0, phi_new in consistent SI */
    static double updatePermeability(double k0, double phi0, double phi_new);

    static std::vector<MineralReaction> defaultCCSReactions();

private:
    std::vector<MineralReaction> reactions_;
    static constexpr double R_gas = 8.31446; // J/(mol·K)
};

} // namespace FSRM

#endif // MINERAL_TRAPPING_HPP
