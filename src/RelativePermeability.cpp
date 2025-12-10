/**
 * @file RelativePermeability.cpp
 * @brief Implementation of relative permeability models
 */

#include "RelativePermeability.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace FSRM {

// =============================================================================
// TwoPhaseRelPerm Implementation
// =============================================================================

TwoPhaseRelPerm::TwoPhaseRelPerm() {
    ep_wetting_.S_connate = 0.2;
    ep_wetting_.S_max = 1.0;
    ep_wetting_.kr_max = 1.0;
    ep_nonwetting_.S_connate = 0.2;
    ep_nonwetting_.S_max = 1.0;
    ep_nonwetting_.kr_max = 1.0;
}

std::pair<double, double> TwoPhaseRelPerm::derivatives(double Sw) const {
    double eps = 1e-6;
    double kr_w_plus = kr_wetting(Sw + eps);
    double kr_w_minus = kr_wetting(Sw - eps);
    double kr_nw_plus = kr_nonwetting(Sw + eps);
    double kr_nw_minus = kr_nonwetting(Sw - eps);
    
    return {(kr_w_plus - kr_w_minus) / (2.0 * eps),
            (kr_nw_plus - kr_nw_minus) / (2.0 * eps)};
}

void TwoPhaseRelPerm::setEndPoints(const EndPoints& wet_ep, const EndPoints& nonwet_ep) {
    ep_wetting_ = wet_ep;
    ep_nonwetting_ = nonwet_ep;
}

double TwoPhaseRelPerm::normalizedSaturation(double Sw) const {
    double Swc = ep_wetting_.S_connate;
    double Sor = ep_nonwetting_.S_connate;
    double Se = (Sw - Swc) / (1.0 - Swc - Sor);
    return std::max(0.0, std::min(1.0, Se));
}

// =============================================================================
// CoreyRelPerm Implementation
// =============================================================================

CoreyRelPerm::CoreyRelPerm() : TwoPhaseRelPerm() {
}

double CoreyRelPerm::kr_wetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    return ep_wetting_.kr_max * std::pow(Se, n_wetting_);
}

double CoreyRelPerm::kr_nonwetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    return ep_nonwetting_.kr_max * std::pow(1.0 - Se, n_nonwetting_);
}

// =============================================================================
// BrooksCoreyRelPerm Implementation
// =============================================================================

BrooksCoreyRelPerm::BrooksCoreyRelPerm() : TwoPhaseRelPerm() {
}

double BrooksCoreyRelPerm::kr_wetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    return ep_wetting_.kr_max * std::pow(Se, (2.0 + 3.0*lambda_) / lambda_);
}

double BrooksCoreyRelPerm::kr_nonwetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    return ep_nonwetting_.kr_max * std::pow(1.0 - Se, 2.0) * 
           (1.0 - std::pow(Se, (2.0 + lambda_)/lambda_));
}

// =============================================================================
// VanGenuchtenRelPerm Implementation
// =============================================================================

VanGenuchtenRelPerm::VanGenuchtenRelPerm() : TwoPhaseRelPerm() {
}

double VanGenuchtenRelPerm::kr_wetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    double inner = 1.0 - std::pow(1.0 - std::pow(Se, 1.0/m_), m_);
    return ep_wetting_.kr_max * std::sqrt(Se) * inner * inner;
}

double VanGenuchtenRelPerm::kr_nonwetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    double Snw = 1.0 - Se;
    return ep_nonwetting_.kr_max * std::sqrt(Snw) * 
           std::pow(1.0 - std::pow(Se, 1.0/m_), 2.0 * m_);
}

void VanGenuchtenRelPerm::setParameters(double m, double n) {
    m_ = m;
    n_ = n;
}

// =============================================================================
// LETRelPerm Implementation
// =============================================================================

LETRelPerm::LETRelPerm() : TwoPhaseRelPerm() {
}

double LETRelPerm::kr_wetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    if (Se <= 0.0) return 0.0;
    if (Se >= 1.0) return ep_wetting_.kr_max;
    
    double num = std::pow(Se, L_wet_);
    double den = num + E_wet_ * std::pow(1.0 - Se, T_wet_);
    return ep_wetting_.kr_max * num / den;
}

double LETRelPerm::kr_nonwetting(double Sw) const {
    double Se = normalizedSaturation(Sw);
    double Snw = 1.0 - Se;
    if (Snw <= 0.0) return 0.0;
    if (Snw >= 1.0) return ep_nonwetting_.kr_max;
    
    double num = std::pow(Snw, L_nw_);
    double den = num + E_nw_ * std::pow(1.0 - Snw, T_nw_);
    return ep_nonwetting_.kr_max * num / den;
}

void LETRelPerm::setWettingLET(double L, double E, double T) {
    L_wet_ = L;
    E_wet_ = E;
    T_wet_ = T;
}

void LETRelPerm::setNonWettingLET(double L, double E, double T) {
    L_nw_ = L;
    E_nw_ = E;
    T_nw_ = T;
}

// =============================================================================
// TabularRelPerm Implementation
// =============================================================================

TabularRelPerm::TabularRelPerm() : TwoPhaseRelPerm() {
}

double TabularRelPerm::kr_wetting(double Sw) const {
    return interpolate(Sw, kr_w_);
}

double TabularRelPerm::kr_nonwetting(double Sw) const {
    return interpolate(Sw, kr_nw_);
}

void TabularRelPerm::setTable(const std::vector<double>& S,
                               const std::vector<double>& kr_w,
                               const std::vector<double>& kr_nw,
                               const std::vector<double>& Pc) {
    S_ = S;
    kr_w_ = kr_w;
    kr_nw_ = kr_nw;
    Pc_ = Pc;
}

double TabularRelPerm::getPc(double Sw) const {
    if (Pc_.empty()) return 0.0;
    return interpolate(Sw, Pc_);
}

void TabularRelPerm::parseFromSWOF(const std::vector<std::string>& data) {
    (void)data;
    // Placeholder for SWOF parsing
}

void TabularRelPerm::parseFromSGOF(const std::vector<std::string>& data) {
    (void)data;
    // Placeholder for SGOF parsing
}

double TabularRelPerm::interpolate(double S, const std::vector<double>& values) const {
    if (S_.empty() || values.empty()) return 0.0;
    if (S_.size() != values.size()) return 0.0;
    
    if (S <= S_.front()) return values.front();
    if (S >= S_.back()) return values.back();
    
    auto it = std::lower_bound(S_.begin(), S_.end(), S);
    size_t idx = it - S_.begin();
    if (idx == 0) idx = 1;
    
    double S0 = S_[idx - 1];
    double S1 = S_[idx];
    double v0 = values[idx - 1];
    double v1 = values[idx];
    
    double t = (S - S0) / (S1 - S0);
    return v0 + t * (v1 - v0);
}

// =============================================================================
// HysteresisRelPerm Implementation
// =============================================================================

HysteresisRelPerm::HysteresisRelPerm(std::shared_ptr<TwoPhaseRelPerm> drainage,
                                     std::shared_ptr<TwoPhaseRelPerm> imbibition,
                                     HysteresisModel model)
    : TwoPhaseRelPerm(), drainage_(drainage), imbibition_(imbibition), model_(model) {
}

double HysteresisRelPerm::kr_wetting(double Sw) const {
    if (state_.increasing && imbibition_) {
        return imbibition_->kr_wetting(Sw);
    }
    return drainage_ ? drainage_->kr_wetting(Sw) : 0.0;
}

double HysteresisRelPerm::kr_nonwetting(double Sw) const {
    if (state_.increasing && imbibition_) {
        return imbibition_->kr_nonwetting(Sw);
    }
    return drainage_ ? drainage_->kr_nonwetting(Sw) : 0.0;
}

void HysteresisRelPerm::updateState(double Sw) {
    state_.update(Sw);
}

double HysteresisRelPerm::getTrappedSaturation() const {
    // Land trapping: S_gt = S_max / (1 + C * S_max)
    double S_max = state_.S_max_historical;
    return S_max / (1.0 + land_parameter_ * S_max);
}

// =============================================================================
// ThreePhaseRelPerm Implementation
// =============================================================================

ThreePhaseRelPerm::ThreePhaseRelPerm() {
}

void ThreePhaseRelPerm::setWaterOilTable(std::shared_ptr<TwoPhaseRelPerm> table) {
    water_oil_ = table;
}

void ThreePhaseRelPerm::setGasOilTable(std::shared_ptr<TwoPhaseRelPerm> table) {
    gas_oil_ = table;
}

void ThreePhaseRelPerm::setGasWaterTable(std::shared_ptr<TwoPhaseRelPerm> table) {
    gas_water_ = table;
}

std::array<double, 3> ThreePhaseRelPerm::calculate(double Sw, double So, double Sg) const {
    std::array<double, 3> kr = {0.0, 0.0, 0.0};
    
    // Normalize saturations
    double S_total = Sw + So + Sg;
    if (S_total <= 0.0) return kr;
    
    Sw = Sw / S_total;
    So = So / S_total;
    Sg = Sg / S_total;
    
    // Water relative permeability from water-oil curve
    if (water_oil_) {
        kr[0] = water_oil_->kr_wetting(Sw);
    }
    
    // Gas relative permeability from gas-oil curve
    if (gas_oil_) {
        kr[2] = gas_oil_->kr_wetting(Sg);  // krg from gas-oil curve
    }
    
    // Oil relative permeability from three-phase model
    switch (model_) {
        case ThreePhaseModel::STONE_I:
            kr[1] = stoneI_oil(Sw, So, Sg);
            break;
        case ThreePhaseModel::STONE_II:
            kr[1] = stoneII_oil(Sw, So, Sg);
            break;
        case ThreePhaseModel::BAKER:
            kr[1] = baker_oil(Sw, So, Sg);
            break;
        default:
            kr[1] = stoneI_oil(Sw, So, Sg);
    }
    
    return kr;
}

std::pair<std::array<double, 3>, std::array<double, 6>> 
ThreePhaseRelPerm::calculateWithDerivatives(double Sw, double So, double Sg) const {
    std::array<double, 3> kr = calculate(Sw, So, Sg);
    std::array<double, 6> dkr = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    // Numerical derivatives (simplified)
    double eps = 1e-6;
    auto kr_dSw = calculate(Sw + eps, So - eps/2, Sg - eps/2);
    auto kr_dSo = calculate(Sw - eps/2, So + eps, Sg - eps/2);
    
    dkr[0] = (kr_dSw[0] - kr[0]) / eps;  // dkrw/dSw
    dkr[1] = (kr_dSo[0] - kr[0]) / eps;  // dkrw/dSo
    dkr[2] = (kr_dSw[1] - kr[1]) / eps;  // dkro/dSw
    dkr[3] = (kr_dSo[1] - kr[1]) / eps;  // dkro/dSo
    dkr[4] = (kr_dSw[2] - kr[2]) / eps;  // dkrg/dSw
    dkr[5] = (kr_dSo[2] - kr[2]) / eps;  // dkrg/dSo
    
    return {kr, dkr};
}

double ThreePhaseRelPerm::kr_water(double Sw, double So, double Sg) const {
    return calculate(Sw, So, Sg)[0];
}

double ThreePhaseRelPerm::kr_oil(double Sw, double So, double Sg) const {
    return calculate(Sw, So, Sg)[1];
}

double ThreePhaseRelPerm::kr_gas(double Sw, double So, double Sg) const {
    return calculate(Sw, So, Sg)[2];
}

double ThreePhaseRelPerm::stoneI_oil(double Sw, double So, double Sg) const {
    if (!water_oil_ || !gas_oil_) return 0.0;
    
    double Swc = 0.2;  // Simplified
    double Sor = 0.2;
    
    double So_e = So / (1.0 - Swc);
    if (So_e <= Sor / (1.0 - Swc)) return 0.0;
    
    // Oil relative permeability at connate water saturation
    double kro_cw = water_oil_->kr_nonwetting(Swc);
    
    // krow at current water saturation
    double krow = water_oil_->kr_nonwetting(Sw);
    
    // krog at current gas saturation (oil-gas curve)
    double krog = gas_oil_->kr_nonwetting(1.0 - Sg);
    
    // Stone I formula
    double alpha = krow / (kro_cw + 1e-10);
    double beta = krog / (kro_cw + 1e-10);
    
    double kro = So_e * alpha * beta / (1.0 - Swc + 1e-10);
    
    return std::max(0.0, std::min(1.0, kro));
}

double ThreePhaseRelPerm::stoneII_oil(double Sw, double So, double Sg) const {
    if (!water_oil_ || !gas_oil_) return 0.0;
    
    double Swc = 0.2;
    double kro_cw = water_oil_->kr_nonwetting(Swc);
    double krow = water_oil_->kr_nonwetting(Sw);
    double krog = gas_oil_->kr_nonwetting(1.0 - Sg);
    double krw = water_oil_->kr_wetting(Sw);
    double krg = gas_oil_->kr_wetting(Sg);
    
    // Stone II formula
    double kro = kro_cw * ((krow / (kro_cw + 1e-10) + krw) * 
                           (krog / (kro_cw + 1e-10) + krg) - (krw + krg));
    
    return std::max(0.0, std::min(1.0, kro));
}

double ThreePhaseRelPerm::baker_oil(double Sw, double So, double Sg) const {
    if (!water_oil_ || !gas_oil_) return 0.0;
    
    // Baker linear interpolation
    double krow = water_oil_->kr_nonwetting(Sw);
    double krog = gas_oil_->kr_nonwetting(1.0 - Sg);
    
    double Swn = Sw / (Sw + Sg + 1e-10);
    double Sgn = Sg / (Sw + Sg + 1e-10);
    
    double kro = Swn * krow + Sgn * krog;
    
    return std::max(0.0, std::min(1.0, kro));
}

// =============================================================================
// CapillaryPressure Implementation
// =============================================================================

CapillaryPressure::CapillaryPressure() {
}

void CapillaryPressure::setParameters(double lambda_or_m) {
    if (model_ == CapillaryModel::BROOKS_COREY) {
        lambda_ = lambda_or_m;
    } else if (model_ == CapillaryModel::VAN_GENUCHTEN) {
        m_vg_ = lambda_or_m;
    }
}

double CapillaryPressure::Pc_ow(double Sw) const {
    double Swc = 0.2;
    double Sor = 0.2;
    double Se = (Sw - Swc) / (1.0 - Swc - Sor);
    Se = std::max(1e-6, std::min(1.0 - 1e-6, Se));
    
    switch (model_) {
        case CapillaryModel::BROOKS_COREY:
            return entry_pressure_ * std::pow(Se, -1.0 / lambda_);
        case CapillaryModel::VAN_GENUCHTEN:
            return (1.0 / alpha_vg_) * std::pow(std::pow(Se, -1.0 / m_vg_) - 1.0, 1.0 - m_vg_);
        case CapillaryModel::TABULAR: {
            // Inline interpolation for tabular Pc
            if (S_table_.empty() || Pc_table_.empty()) return 0.0;
            if (Sw <= S_table_.front()) return Pc_table_.front();
            if (Sw >= S_table_.back()) return Pc_table_.back();
            
            auto it = std::lower_bound(S_table_.begin(), S_table_.end(), Sw);
            size_t idx = it - S_table_.begin();
            if (idx == 0) idx = 1;
            
            double S0 = S_table_[idx - 1];
            double S1 = S_table_[idx];
            double v0 = Pc_table_[idx - 1];
            double v1 = Pc_table_[idx];
            
            double t = (Sw - S0) / (S1 - S0);
            return v0 + t * (v1 - v0);
        }
        default:
            return entry_pressure_ * std::pow(Se, -1.0 / lambda_);
    }
}

double CapillaryPressure::Pc_go(double Sg) const {
    double Se = 1.0 - Sg;
    Se = std::max(1e-6, std::min(1.0 - 1e-6, Se));
    return entry_pressure_ * std::pow(Se, -1.0 / lambda_);
}

double CapillaryPressure::Pc_gw(double Sw) const {
    return Pc_ow(Sw) + Pc_go(1.0 - Sw);
}

void CapillaryPressure::setTabularPc(const std::vector<double>& S,
                                      const std::vector<double>& Pc) {
    S_table_ = S;
    Pc_table_ = Pc;
    model_ = CapillaryModel::TABULAR;
}

void CapillaryPressure::setJFunction(double k, double phi, double sigma, double theta) {
    (void)k;
    (void)phi;
    (void)sigma;
    (void)theta;
}

double CapillaryPressure::JFunction(double Sw) const {
    double Swc = 0.2;
    double Sor = 0.2;
    double Se = (Sw - Swc) / (1.0 - Swc - Sor);
    Se = std::max(1e-6, std::min(1.0 - 1e-6, Se));
    return std::pow(Se, -1.0 / lambda_);
}

double CapillaryPressure::PcFromJ(double Sw, double k, double phi) const {
    double sigma = 0.025;
    double theta = 0.0;
    double J = JFunction(Sw);
    return sigma * std::cos(theta) * std::sqrt(phi / k) * J;
}

} // namespace FSRM
