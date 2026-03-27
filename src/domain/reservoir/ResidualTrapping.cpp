/**
 * @file ResidualTrapping.cpp
 * @brief Land residual trapping and Killough hysteretic krg for CO2
 */

#include "domain/reservoir/ResidualTrapping.hpp"

#include <algorithm>
#include <cmath>

namespace FSRM {

namespace {

constexpr double kSaturationTol = 1.0e-14;

} // namespace

double ResidualTrapping::landTrappedSaturation(double Sg_max, double C_land) {
    if (Sg_max <= 0.0) {
        return 0.0;
    }
    if (!std::isfinite(Sg_max) || !std::isfinite(C_land)) {
        return 0.0;
    }
    const double denom = 1.0 + C_land * Sg_max;
    if (denom <= 0.0) {
        return 0.0;
    }
    const double sgr = Sg_max / denom;
    return std::clamp(sgr, 0.0, Sg_max);
}

void ResidualTrapping::updateHistory(int cell, double Sg_current) {
    if (cell < 0 || cell >= static_cast<int>(Sg_max_history_.size())) {
        return;
    }
    const double s = std::max(0.0, Sg_current);
    const std::size_t i = static_cast<std::size_t>(cell);
    Sg_max_history_[i] = std::max(Sg_max_history_[i], s);
}

double ResidualTrapping::getMaxHistoricalSg(int cell) const {
    if (cell < 0 || cell >= static_cast<int>(Sg_max_history_.size())) {
        return 0.0;
    }
    return Sg_max_history_[static_cast<std::size_t>(cell)];
}

double ResidualTrapping::getTrappedSg(int cell) const {
    return landTrappedSaturation(getMaxHistoricalSg(cell), C_land_);
}

double ResidualTrapping::krg_drainage(double Sg) const {
    if (!std::isfinite(Sg)) {
        return 0.0;
    }
    if (Sg <= Sgr_drain_) {
        return 0.0;
    }
    const double denom = 1.0 - Sgr_drain_;
    if (denom <= kSaturationTol) {
        return 0.0;
    }
    double Se = (Sg - Sgr_drain_) / denom;
    Se = std::clamp(Se, 0.0, 1.0);
    if (Se <= 0.0) {
        return 0.0;
    }
    return krg_max_ * std::pow(Se, ng_);
}

double ResidualTrapping::krg_imbibition(double Sg, int cell) const {
    if (!std::isfinite(Sg)) {
        return 0.0;
    }
    const double kr_d = krg_drainage(Sg);
    const double Sg_max_hist = getMaxHistoricalSg(cell);

    if (Sg_max_hist <= kSaturationTol) {
        return std::max(0.0, kr_d);
    }

    const double Sgr_trapped = landTrappedSaturation(Sg_max_hist, C_land_);
    const double span = Sg_max_hist - Sgr_trapped;

    if (Sg <= Sgr_trapped + kSaturationTol) {
        return 0.0;
    }

    if (Sg >= Sg_max_hist - kSaturationTol) {
        return std::max(0.0, krg_drainage(Sg));
    }

    const double kr_at_max = krg_drainage(Sg_max_hist);

    if (span <= kSaturationTol) {
        return std::clamp(std::min(kr_at_max, kr_d), 0.0, kr_d);
    }
    const double ratio = (Sg - Sgr_trapped) / span;
    double kr_imb = kr_at_max * std::pow(ratio, ng_);
    kr_imb = std::clamp(kr_imb, 0.0, kr_d);
    return kr_imb;
}

void ResidualTrapping::setLandCoefficient(double C) {
    C_land_ = std::isfinite(C) ? C : C_land_;
}

void ResidualTrapping::setDrainageRelPerm(double Sgr_drain, double krg_max, double ng) {
    if (std::isfinite(Sgr_drain)) {
        Sgr_drain_ = std::clamp(Sgr_drain, 0.0, 1.0);
    }
    if (std::isfinite(krg_max) && krg_max >= 0.0) {
        krg_max_ = krg_max;
    }
    if (std::isfinite(ng) && ng > 0.0) {
        ng_ = ng;
    }
}

void ResidualTrapping::setNumCells(int n) {
    if (n < 0) {
        Sg_max_history_.clear();
        return;
    }
    Sg_max_history_.assign(static_cast<std::size_t>(n), 0.0);
}

} // namespace FSRM
