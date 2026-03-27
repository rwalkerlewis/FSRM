#ifndef RESIDUAL_TRAPPING_HPP
#define RESIDUAL_TRAPPING_HPP

/**
 * @file ResidualTrapping.hpp
 * @brief Residual CO2 trapping (Land 1968) and Killough (1976) hysteretic gas rel-perm
 *
 * Saturation is dimensionless (SI consistent). Land coefficient C and Corey exponents
 * are dimensionless; krg_max is dimensionless (relative permeability).
 */

#include <vector>

namespace FSRM {

/**
 * @brief Residual trapping and hysteretic gas relative permeability for CCS
 */
class ResidualTrapping {
public:
    /// Land (1968): Sgr = Sg_max / (1 + C * Sg_max)
    static double landTrappedSaturation(double Sg_max, double C_land);

    void updateHistory(int cell, double Sg_current);
    double getMaxHistoricalSg(int cell) const;
    double getTrappedSg(int cell) const;

    /// Killough (1976) imbibition scanning curve (gas phase)
    double krg_imbibition(double Sg, int cell) const;
    /// Corey drainage gas relative permeability
    double krg_drainage(double Sg) const;

    void setLandCoefficient(double C);
    void setDrainageRelPerm(double Sgr_drain, double krg_max, double ng);
    void setNumCells(int n);

private:
    double C_land_{2.0};
    double Sgr_drain_{0.05};
    double krg_max_{1.0};
    double ng_{2.0};
    std::vector<double> Sg_max_history_;
};

} // namespace FSRM

#endif // RESIDUAL_TRAPPING_HPP
