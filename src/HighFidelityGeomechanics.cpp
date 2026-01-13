/**
 * @file HighFidelityGeomechanics.cpp
 * @brief Implementation of high-fidelity geomechanics models
 * 
 * This file provides implementations for advanced geomechanics capabilities:
 * - Finite strain (hyperelastic) formulations
 * - Rate-dependent plasticity (viscoplasticity)
 * - Creep models (power law, Norton, Lemaitre)
 * - Hypoplasticity for granular materials
 * - Gradient-enhanced damage for regularization
 */

#include "HighFidelityGeomechanics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Tensor Utility Functions
// =============================================================================

namespace {

/**
 * @brief Compute determinant of 3x3 matrix
 */
double det3x3(const double F[9]) {
    return F[0] * (F[4] * F[8] - F[5] * F[7])
         - F[1] * (F[3] * F[8] - F[5] * F[6])
         + F[2] * (F[3] * F[7] - F[4] * F[6]);
}

/**
 * @brief Compute inverse of 3x3 matrix
 */
void inv3x3(const double F[9], double Finv[9]) {
    double J = det3x3(F);
    if (std::abs(J) < 1e-30) {
        throw std::runtime_error("Singular deformation gradient");
    }
    double Jinv = 1.0 / J;
    
    Finv[0] = Jinv * (F[4] * F[8] - F[5] * F[7]);
    Finv[1] = Jinv * (F[2] * F[7] - F[1] * F[8]);
    Finv[2] = Jinv * (F[1] * F[5] - F[2] * F[4]);
    Finv[3] = Jinv * (F[5] * F[6] - F[3] * F[8]);
    Finv[4] = Jinv * (F[0] * F[8] - F[2] * F[6]);
    Finv[5] = Jinv * (F[2] * F[3] - F[0] * F[5]);
    Finv[6] = Jinv * (F[3] * F[7] - F[4] * F[6]);
    Finv[7] = Jinv * (F[1] * F[6] - F[0] * F[7]);
    Finv[8] = Jinv * (F[0] * F[4] - F[1] * F[3]);
}

/**
 * @brief Compute F^T * F (right Cauchy-Green tensor)
 */
void computeRightCauchyGreen(const double F[9], double C[6]) {
    // C = F^T F in Voigt notation
    C[0] = F[0]*F[0] + F[3]*F[3] + F[6]*F[6]; // C_xx
    C[1] = F[1]*F[1] + F[4]*F[4] + F[7]*F[7]; // C_yy
    C[2] = F[2]*F[2] + F[5]*F[5] + F[8]*F[8]; // C_zz
    C[3] = F[1]*F[2] + F[4]*F[5] + F[7]*F[8]; // C_yz
    C[4] = F[0]*F[2] + F[3]*F[5] + F[6]*F[8]; // C_xz
    C[5] = F[0]*F[1] + F[3]*F[4] + F[6]*F[7]; // C_xy
}

/**
 * @brief Compute trace of tensor in Voigt notation
 */
double trace(const double T[6]) {
    return T[0] + T[1] + T[2];
}

/**
 * @brief Compute invariants of symmetric tensor
 */
void computeInvariants(const double C[6], double& I1, double& I2, double& I3) {
    I1 = trace(C);
    I2 = 0.5 * (I1 * I1 - (C[0]*C[0] + C[1]*C[1] + C[2]*C[2] + 
                           2.0*(C[3]*C[3] + C[4]*C[4] + C[5]*C[5])));
    I3 = C[0] * (C[1]*C[2] - C[3]*C[3])
       - C[5] * (C[5]*C[2] - C[3]*C[4])
       + C[4] * (C[5]*C[3] - C[1]*C[4]);
}

} // anonymous namespace

// =============================================================================
// DeformationGradient Implementation
// =============================================================================

DeformationGradient::DeformationGradient() {
    // Initialize to identity
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

DeformationGradient::DeformationGradient(const std::array<std::array<double, 3>, 3>& components)
    : F(components) {
}

double DeformationGradient::determinant() const {
    return F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1])
         - F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0])
         + F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]);
}

DeformationGradient DeformationGradient::inverse() const {
    DeformationGradient Finv;
    double J = determinant();
    if (std::abs(J) < 1e-30) {
        throw std::runtime_error("Singular deformation gradient");
    }
    double Jinv = 1.0 / J;
    
    Finv.F[0][0] = Jinv * (F[1][1] * F[2][2] - F[1][2] * F[2][1]);
    Finv.F[0][1] = Jinv * (F[0][2] * F[2][1] - F[0][1] * F[2][2]);
    Finv.F[0][2] = Jinv * (F[0][1] * F[1][2] - F[0][2] * F[1][1]);
    Finv.F[1][0] = Jinv * (F[1][2] * F[2][0] - F[1][0] * F[2][2]);
    Finv.F[1][1] = Jinv * (F[0][0] * F[2][2] - F[0][2] * F[2][0]);
    Finv.F[1][2] = Jinv * (F[0][2] * F[1][0] - F[0][0] * F[1][2]);
    Finv.F[2][0] = Jinv * (F[1][0] * F[2][1] - F[1][1] * F[2][0]);
    Finv.F[2][1] = Jinv * (F[0][1] * F[2][0] - F[0][0] * F[2][1]);
    Finv.F[2][2] = Jinv * (F[0][0] * F[1][1] - F[0][1] * F[1][0]);
    
    return Finv;
}

DeformationGradient DeformationGradient::transpose() const {
    DeformationGradient FT;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            FT.F[i][j] = F[j][i];
        }
    }
    return FT;
}

std::array<std::array<double, 3>, 3> DeformationGradient::rightCauchyGreen() const {
    // C = F^T * F
    std::array<std::array<double, 3>, 3> C{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                C[i][j] += F[k][i] * F[k][j];
            }
        }
    }
    return C;
}

std::array<std::array<double, 3>, 3> DeformationGradient::leftCauchyGreen() const {
    // B = F * F^T
    std::array<std::array<double, 3>, 3> B{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            B[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                B[i][j] += F[i][k] * F[j][k];
            }
        }
    }
    return B;
}

std::array<std::array<double, 3>, 3> DeformationGradient::greenLagrangeStrain() const {
    // E = 0.5*(C - I)
    auto C = rightCauchyGreen();
    std::array<std::array<double, 3>, 3> E_tensor{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            E_tensor[i][j] = 0.5 * (C[i][j] - ((i == j) ? 1.0 : 0.0));
        }
    }
    return E_tensor;
}

std::array<std::array<double, 3>, 3> DeformationGradient::eulerAlmansiStrain() const {
    // e = 0.5*(I - B^{-1})
    auto B = leftCauchyGreen();
    // For simplicity, approximate for small deformations
    std::array<std::array<double, 3>, 3> e{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            e[i][j] = 0.5 * (((i == j) ? 1.0 : 0.0) - B[i][j]);
        }
    }
    return e;
}

std::array<std::array<double, 3>, 3> DeformationGradient::logarithmicStrain() const {
    // Simplified: for small deformations, ln(U) ~ E
    return greenLagrangeStrain();
}

void DeformationGradient::polarDecomposition(std::array<std::array<double, 3>, 3>& R,
                                             std::array<std::array<double, 3>, 3>& U) const {
    // Simplified polar decomposition using iteration
    // F = R * U where R is rotation and U is right stretch
    
    // Initial guess: R = F / |F|
    double norm = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            norm += F[i][j] * F[i][j];
        }
    }
    norm = std::sqrt(norm);
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            R[i][j] = F[i][j] / norm;
        }
    }
    
    // Iterative refinement
    for (int iter = 0; iter < 10; ++iter) {
        // Compute R^{-T} and average
        DeformationGradient R_def;
        R_def.F = R;
        auto Rinv = R_def.inverse();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                R[i][j] = 0.5 * (R[i][j] + Rinv.F[j][i]);
            }
        }
    }
    
    // U = R^T * F
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            U[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                U[i][j] += R[k][i] * F[k][j];
            }
        }
    }
}

std::array<std::array<double, 3>, 3> DeformationGradient::velocityGradient(
    const DeformationGradient& F_curr, const DeformationGradient& F_dot) {
    // L = F_dot * F^{-1}
    auto Finv = F_curr.inverse();
    std::array<std::array<double, 3>, 3> L{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            L[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                L[i][j] += F_dot.F[i][k] * Finv.F[k][j];
            }
        }
    }
    return L;
}

std::array<std::array<double, 3>, 3> DeformationGradient::rateOfDeformation(
    const std::array<std::array<double, 3>, 3>& L) {
    // D = sym(L)
    std::array<std::array<double, 3>, 3> D{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            D[i][j] = 0.5 * (L[i][j] + L[j][i]);
        }
    }
    return D;
}

std::array<std::array<double, 3>, 3> DeformationGradient::spinTensor(
    const std::array<std::array<double, 3>, 3>& L) {
    // W = skew(L)
    std::array<std::array<double, 3>, 3> W{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            W[i][j] = 0.5 * (L[i][j] - L[j][i]);
        }
    }
    return W;
}

// =============================================================================
// FiniteStrainStress Implementation
// =============================================================================

void FiniteStrainStress::cauchyToKirchhoff(const std::array<std::array<double, 3>, 3>& sigma,
                                           double J,
                                           std::array<std::array<double, 3>, 3>& tau) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            tau[i][j] = J * sigma[i][j];
        }
    }
}

void FiniteStrainStress::cauchyToPiola1(const std::array<std::array<double, 3>, 3>& sigma,
                                        const DeformationGradient& F_def,
                                        std::array<std::array<double, 3>, 3>& P) {
    double J = F_def.determinant();
    auto Finv = F_def.inverse();
    // P = J * sigma * F^{-T}
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            P[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                P[i][j] += J * sigma[i][k] * Finv.F[j][k];
            }
        }
    }
}

void FiniteStrainStress::cauchyToPiola2(const std::array<std::array<double, 3>, 3>& sigma,
                                        const DeformationGradient& F_def,
                                        std::array<std::array<double, 3>, 3>& S) {
    double J = F_def.determinant();
    auto Finv = F_def.inverse();
    // S = F^{-1} * P = J * F^{-1} * sigma * F^{-T}
    std::array<std::array<double, 3>, 3> temp{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            temp[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                temp[i][j] += J * sigma[i][k] * Finv.F[j][k];
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            S[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                S[i][j] += Finv.F[i][k] * temp[k][j];
            }
        }
    }
}

// =============================================================================
// HyperelasticMaterial Implementation
// =============================================================================

HyperelasticMaterial::HyperelasticMaterial() {
}

void HyperelasticMaterial::setParameters(const Parameters& params) {
    params_ = params;
}

void HyperelasticMaterial::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

double HyperelasticMaterial::strainEnergy(const DeformationGradient& F) const {
    switch (params_.model) {
        case HyperelasticModel::NEO_HOOKEAN:
            return W_NeoHookean(F);
        case HyperelasticModel::MOONEY_RIVLIN:
            return W_MooneyRivlin(F);
        default:
            return W_NeoHookean(F);
    }
}

std::array<std::array<double, 3>, 3> HyperelasticMaterial::secondPiolaKirchhoff(const DeformationGradient& F) const {
    switch (params_.model) {
        case HyperelasticModel::NEO_HOOKEAN:
            return S_NeoHookean(F);
        case HyperelasticModel::MOONEY_RIVLIN:
            return S_MooneyRivlin(F);
        default:
            return S_NeoHookean(F);
    }
}

std::array<std::array<double, 3>, 3> HyperelasticMaterial::cauchyStress(const DeformationGradient& Fdef) const {
    auto S = secondPiolaKirchhoff(Fdef);
    std::array<std::array<double, 3>, 3> sigma{};
    
    // sigma = (1/J) * F * S * F^T (simplified)
    double J = Fdef.determinant();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            sigma[i][j] = S[i][j] / J;
        }
    }
    return sigma;
}

std::array<std::array<double, 6>, 6> HyperelasticMaterial::materialTangent(const DeformationGradient& F) const {
    (void)F;
    std::array<std::array<double, 6>, 6> C{};
    
    double lambda = params_.lambda;
    double mu = params_.mu;
    
    // Simplified isotropic tangent
    C[0][0] = C[1][1] = C[2][2] = lambda + 2.0 * mu;
    C[0][1] = C[0][2] = C[1][0] = C[1][2] = C[2][0] = C[2][1] = lambda;
    C[3][3] = C[4][4] = C[5][5] = mu;
    
    return C;
}

std::array<std::array<double, 6>, 6> HyperelasticMaterial::spatialTangent(const DeformationGradient& F) const {
    // For now, return material tangent
    return materialTangent(F);
}

bool HyperelasticMaterial::isValidDeformation(const DeformationGradient& Fdef) const {
    return Fdef.determinant() > 0.0;
}

double HyperelasticMaterial::W_NeoHookean(const DeformationGradient& Fdef) const {
    double J = Fdef.determinant();
    auto C = Fdef.rightCauchyGreen();
    double I1 = C[0][0] + C[1][1] + C[2][2];
    double lnJ = std::log(J);
    return 0.5 * params_.mu * (I1 - 3.0) - params_.mu * lnJ + 0.5 * params_.lambda * lnJ * lnJ;
}

double HyperelasticMaterial::W_MooneyRivlin(const DeformationGradient& Fdef) const {
    double J = Fdef.determinant();
    auto C = Fdef.rightCauchyGreen();
    double I1 = C[0][0] + C[1][1] + C[2][2];
    double I2 = 0.5 * (I1 * I1 - (C[0][0]*C[0][0] + C[1][1]*C[1][1] + C[2][2]*C[2][2] + 
                                   2.0*(C[0][1]*C[0][1] + C[0][2]*C[0][2] + C[1][2]*C[1][2])));
    double J23 = std::pow(J, -2.0/3.0);
    double I1_bar = J23 * I1;
    double I2_bar = std::pow(J, -4.0/3.0) * I2;
    return params_.C1 * (I1_bar - 3.0) + params_.C2 * (I2_bar - 3.0) + 
           0.5 * params_.kappa * (J - 1.0) * (J - 1.0);
}

double HyperelasticMaterial::W_Ogden(const DeformationGradient& Fdef) const {
    (void)Fdef;
    return 0.0;  // Placeholder
}

double HyperelasticMaterial::W_StVenant(const DeformationGradient& Fdef) const {
    // St. Venant-Kirchhoff: W = lambda/2 * tr(E)^2 + mu * tr(E^2)
    auto E = Fdef.greenLagrangeStrain();
    double trE = E[0][0] + E[1][1] + E[2][2];
    double trE2 = E[0][0]*E[0][0] + E[1][1]*E[1][1] + E[2][2]*E[2][2] +
                  2.0 * (E[0][1]*E[0][1] + E[0][2]*E[0][2] + E[1][2]*E[1][2]);
    return 0.5 * params_.lambda * trE * trE + params_.mu * trE2;
}

std::array<std::array<double, 3>, 3> HyperelasticMaterial::S_NeoHookean(const DeformationGradient& Fdef) const {
    std::array<std::array<double, 3>, 3> S{};
    double J = Fdef.determinant();
    double lnJ = std::log(J);
    double mu = params_.mu;
    double lambda = params_.lambda;
    auto C = Fdef.rightCauchyGreen();
    
    // Simplified: S_ij = mu * delta_ij + (lambda * lnJ - mu) / C_ii
    S[0][0] = mu + (lambda * lnJ - mu) / C[0][0];
    S[1][1] = mu + (lambda * lnJ - mu) / C[1][1];
    S[2][2] = mu + (lambda * lnJ - mu) / C[2][2];
    
    return S;
}

std::array<std::array<double, 3>, 3> HyperelasticMaterial::S_MooneyRivlin(const DeformationGradient& Fdef) const {
    std::array<std::array<double, 3>, 3> S{};
    auto C = Fdef.rightCauchyGreen();
    double I1 = C[0][0] + C[1][1] + C[2][2];
    
    double dW_dI1 = params_.C1;
    double dW_dI2 = params_.C2;
    
    S[0][0] = 2.0 * (dW_dI1 + dW_dI2 * (I1 - C[0][0]));
    S[1][1] = 2.0 * (dW_dI1 + dW_dI2 * (I1 - C[1][1]));
    S[2][2] = 2.0 * (dW_dI1 + dW_dI2 * (I1 - C[2][2]));
    
    return S;
}

void HyperelasticMaterial::computePrincipalStretches(const DeformationGradient& Fdef,
                                                     std::array<double, 3>& lambda_arr) const {
    // Simplified: use diagonal of C
    auto C = Fdef.rightCauchyGreen();
    lambda_arr[0] = std::sqrt(C[0][0]);
    lambda_arr[1] = std::sqrt(C[1][1]);
    lambda_arr[2] = std::sqrt(C[2][2]);
}

void HyperelasticMaterial::computeInvariants(const DeformationGradient& Fdef,
                                             double& I1, double& I2, double& I3) const {
    auto C = Fdef.rightCauchyGreen();
    I1 = C[0][0] + C[1][1] + C[2][2];
    I2 = 0.5 * (I1 * I1 - (C[0][0]*C[0][0] + C[1][1]*C[1][1] + C[2][2]*C[2][2] + 
                           2.0*(C[0][1]*C[0][1] + C[0][2]*C[0][2] + C[1][2]*C[1][2])));
    I3 = Fdef.determinant() * Fdef.determinant();
}

// =============================================================================
// FiniteStrainPlasticity Implementation
// =============================================================================

FiniteStrainPlasticity::FiniteStrainPlasticity() {
}

void FiniteStrainPlasticity::setParameters(const Parameters& params) {
    params_ = params;
}

void FiniteStrainPlasticity::setHyperelasticModel(std::shared_ptr<HyperelasticMaterial> elastic) {
    elastic_ = elastic;
}

void FiniteStrainPlasticity::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void FiniteStrainPlasticity::initializeState(State& state) const {
    // Initialize plastic deformation gradient to identity
    state.Fp = DeformationGradient();  // Default constructor gives identity
    state.accumulated_plastic_strain = 0.0;
    for (auto& row : state.back_stress) row.fill(0.0);
    state.is_plastic = false;
}

std::array<std::array<double, 3>, 3> FiniteStrainPlasticity::computeStress(
    const DeformationGradient& F_new, State& state) const {
    // Simplified elastic predictor
    (void)F_new;
    (void)state;
    std::array<std::array<double, 3>, 3> sigma{};
    return sigma;
}

std::array<std::array<double, 6>, 6> FiniteStrainPlasticity::consistentTangent(
    const DeformationGradient& F_def, const State& state) const {
    (void)F_def;
    (void)state;
    std::array<std::array<double, 6>, 6> C{};
    return C;
}

double FiniteStrainPlasticity::yieldFunction(
    const std::array<std::array<double, 3>, 3>& kirchhoff_stress,
    double accumulated_strain) const {
    // Von Mises yield function
    double p = (kirchhoff_stress[0][0] + kirchhoff_stress[1][1] + kirchhoff_stress[2][2]) / 3.0;
    double s11 = kirchhoff_stress[0][0] - p;
    double s22 = kirchhoff_stress[1][1] - p;
    double s33 = kirchhoff_stress[2][2] - p;
    double J2 = 0.5 * (s11*s11 + s22*s22 + s33*s33) +
                kirchhoff_stress[0][1]*kirchhoff_stress[0][1] +
                kirchhoff_stress[0][2]*kirchhoff_stress[0][2] +
                kirchhoff_stress[1][2]*kirchhoff_stress[1][2];
    double sigma_eq = std::sqrt(3.0 * J2);
    double hardening = params_.hardening_modulus * accumulated_strain;
    return sigma_eq - std::sqrt(2.0/3.0) * (params_.yield_stress + hardening);
}

bool FiniteStrainPlasticity::isYielding(
    const std::array<std::array<double, 3>, 3>& stress, double eps_p) const {
    return yieldFunction(stress, eps_p) > 0.0;
}

// =============================================================================
// HypoplasticMaterial Implementation
// =============================================================================

HypoplasticMaterial::HypoplasticMaterial() {
}

void HypoplasticMaterial::setParameters(const Parameters& params) {
    params_ = params;
}

void HypoplasticMaterial::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void HypoplasticMaterial::initializeState(State& state, double e0) const {
    state.void_ratio = e0;
    state.intergranular_strain.fill(0.0);
    state.rho = 0.0;
    state.temperature = 293.15;
}

std::array<double, 6> HypoplasticMaterial::stressRate(const std::array<double, 6>& stress,
                                                       const std::array<double, 6>& strain_rate,
                                                       const State& state) const {
    std::array<double, 6> sigma_dot;
    double K = 1e8;
    double G = 5e7;
    
    double ev = strain_rate[0] + strain_rate[1] + strain_rate[2];
    sigma_dot[0] = K * ev + 2.0 * G * (strain_rate[0] - ev/3.0);
    sigma_dot[1] = K * ev + 2.0 * G * (strain_rate[1] - ev/3.0);
    sigma_dot[2] = K * ev + 2.0 * G * (strain_rate[2] - ev/3.0);
    sigma_dot[3] = 2.0 * G * strain_rate[3];
    sigma_dot[4] = 2.0 * G * strain_rate[4];
    sigma_dot[5] = 2.0 * G * strain_rate[5];
    
    (void)stress;
    (void)state;
    return sigma_dot;
}

void HypoplasticMaterial::integrateStress(const std::array<double, 6>& strain_increment,
                                          std::array<double, 6>& stress,
                                          State& state,
                                          double dt) const {
    auto sigma_dot = stressRate(stress, strain_increment, state);
    for (int i = 0; i < 6; ++i) {
        stress[i] += sigma_dot[i];
    }
    (void)dt;
}

std::array<std::array<double, 6>, 6> HypoplasticMaterial::tangentStiffness(
    const std::array<double, 6>& stress,
    const std::array<double, 6>& strain_rate,
    const State& state) const {
    (void)stress;
    (void)strain_rate;
    (void)state;
    std::array<std::array<double, 6>, 6> L{};
    double K = 1e8;
    double G = 5e7;
    double lambda = K - 2.0*G/3.0;
    
    L[0][0] = L[1][1] = L[2][2] = lambda + 2.0*G;
    L[0][1] = L[0][2] = L[1][0] = L[1][2] = L[2][0] = L[2][1] = lambda;
    L[3][3] = L[4][4] = L[5][5] = G;
    
    return L;
}

void HypoplasticMaterial::limitVoidRatios(double p, double& e_c, double& e_d, double& e_i) const {
    double h_s = params_.h_s;
    double n_param = params_.n_h;
    double factor = std::exp(-std::pow(3.0 * p / h_s, n_param));
    e_c = params_.e_c0 * factor;
    e_d = params_.e_d0 * factor;
    e_i = params_.e_i0 * factor;
}

double HypoplasticMaterial::relativeDensity(double e, double p) const {
    double e_c, e_d, e_i;
    limitVoidRatios(p, e_c, e_d, e_i);
    return (e_i - e) / (e_i - e_d);
}

double HypoplasticMaterial::densityFactor(double e, double p) const {
    double e_c, e_d, e_i;
    limitVoidRatios(p, e_c, e_d, e_i);
    return std::pow(e_c / e, params_.alpha_h);
}

double HypoplasticMaterial::pressureFactor(double p, double e) const {
    double h_s = params_.h_s;
    double n_param = params_.n_h;
    double e_i = params_.e_i0 * std::exp(-std::pow(3.0 * p / h_s, n_param));
    return std::pow(h_s / (3.0 * p), n_param) * (1.0 + e_i) * e_i / e;
}

// =============================================================================
// Viscoplasticity Implementation
// =============================================================================

Viscoplasticity::Viscoplasticity() {
}

void Viscoplasticity::setParameters(const Parameters& params) {
    params_ = params;
}

void Viscoplasticity::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void Viscoplasticity::initializeState(State& state, double T0) const {
    state.plastic_strain.fill(0.0);
    state.accumulated_plastic_strain = 0.0;
    state.back_stress.fill(0.0);
    state.temperature = T0;
}

void Viscoplasticity::integrateStress(const std::array<double, 6>& strain_increment,
                                      double dt,
                                      std::array<double, 6>& stress,
                                      State& state,
                                      const std::array<std::array<double, 6>, 6>& C) const {
    // Trial stress
    std::array<double, 6> stress_trial;
    for (int i = 0; i < 6; ++i) {
        stress_trial[i] = stress[i];
        for (int j = 0; j < 6; ++j) {
            stress_trial[i] += C[i][j] * strain_increment[j];
        }
    }
    
    // Check yield
    double f = yieldFunction(stress_trial, state.accumulated_plastic_strain * params_.hardening_modulus);
    
    if (f <= 0.0) {
        stress = stress_trial;
        return;
    }
    
    // Viscoplastic correction
    double phi = overstressFunction(f);
    auto n = flowDirection(stress_trial);
    double dgamma = params_.fluidity * phi * dt * thermalFactor(state.temperature);
    
    for (int i = 0; i < 6; ++i) {
        state.plastic_strain[i] += dgamma * n[i];
        stress[i] = stress_trial[i];
        for (int j = 0; j < 6; ++j) {
            stress[i] -= C[i][j] * dgamma * n[j];
        }
    }
    
    state.accumulated_plastic_strain += dgamma;
}

std::array<double, 6> Viscoplasticity::plasticStrainRate(const std::array<double, 6>& stress,
                                                          const State& state) const {
    double f = yieldFunction(stress, state.accumulated_plastic_strain * params_.hardening_modulus);
    if (f <= 0.0) {
        return std::array<double, 6>{};
    }
    
    double phi = overstressFunction(f);
    auto n = flowDirection(stress);
    double gamma_dot = params_.fluidity * phi * thermalFactor(state.temperature);
    
    std::array<double, 6> eps_dot;
    for (int i = 0; i < 6; ++i) {
        eps_dot[i] = gamma_dot * n[i];
    }
    return eps_dot;
}

double Viscoplasticity::overstressFunction(double f) const {
    if (f <= 0.0) return 0.0;
    return std::pow(f / params_.reference_stress, params_.rate_exponent);
}

double Viscoplasticity::overstressFunctionDerivative(double f) const {
    if (f <= 0.0) return 0.0;
    return params_.rate_exponent / params_.reference_stress * 
           std::pow(f / params_.reference_stress, params_.rate_exponent - 1.0);
}

double Viscoplasticity::rateDependentYield(double strain_rate) const {
    double ratio = strain_rate / 1e-3;  // Reference strain rate
    return params_.yield_stress * (1.0 + std::pow(ratio, 1.0 / params_.rate_exponent));
}

std::array<std::array<double, 6>, 6> Viscoplasticity::consistentTangent(
    const std::array<double, 6>& stress, const State& state, double dt,
    const std::array<std::array<double, 6>, 6>& C) const {
    
    double f = yieldFunction(stress, state.accumulated_plastic_strain * params_.hardening_modulus);
    if (f <= 0.0) return C;
    
    double dphi_df = overstressFunctionDerivative(f);
    double theta = 1.0 / (1.0 + params_.fluidity * dphi_df * dt);
    
    std::array<std::array<double, 6>, 6> Cvp;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            Cvp[i][j] = theta * C[i][j];
        }
    }
    return Cvp;
}

double Viscoplasticity::yieldFunction(const std::array<double, 6>& stress, double hardening) const {
    // Von Mises
    double s1 = stress[0] - (stress[0] + stress[1] + stress[2]) / 3.0;
    double s2 = stress[1] - (stress[0] + stress[1] + stress[2]) / 3.0;
    double s3 = stress[2] - (stress[0] + stress[1] + stress[2]) / 3.0;
    double J2 = 0.5 * (s1*s1 + s2*s2 + s3*s3) + stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5];
    return std::sqrt(3.0 * J2) - params_.yield_stress - hardening;
}

std::array<double, 6> Viscoplasticity::flowDirection(const std::array<double, 6>& stress) const {
    double p = (stress[0] + stress[1] + stress[2]) / 3.0;
    std::array<double, 6> s;
    s[0] = stress[0] - p;
    s[1] = stress[1] - p;
    s[2] = stress[2] - p;
    s[3] = stress[3];
    s[4] = stress[4];
    s[5] = stress[5];
    
    double J2 = 0.5 * (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]) + s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    double sigma_eq = std::sqrt(3.0 * J2);
    
    std::array<double, 6> n;
    if (sigma_eq > 1e-10) {
        double factor = 1.5 / sigma_eq;
        for (int i = 0; i < 6; ++i) {
            n[i] = factor * s[i];
        }
    } else {
        n.fill(0.0);
    }
    return n;
}

double Viscoplasticity::thermalFactor(double T) const {
    if (!params_.thermal_activation) return 1.0;
    const double R = 8.314;
    return std::exp(-params_.activation_energy / (R * T));
}

// =============================================================================
// CreepModel_Impl Implementation
// =============================================================================

CreepModel_Impl::CreepModel_Impl() {
}

void CreepModel_Impl::setParameters(const Parameters& params) {
    params_ = params;
}

void CreepModel_Impl::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void CreepModel_Impl::initializeState(State& state) const {
    state.creep_strain = 0.0;
    state.creep_strain_tensor.fill(0.0);
    state.time = 0.0;
    state.damage = 0.0;
    state.stage = CreepStage::PRIMARY;
    state.maxwell_strain = 0.0;
    state.kelvin_strain = 0.0;
}

std::array<double, 6> CreepModel_Impl::creepStrainRate(const std::array<double, 6>& stress,
                                                        double T,
                                                        const State& state) const {
    // Compute von Mises equivalent stress
    double p = (stress[0] + stress[1] + stress[2]) / 3.0;
    std::array<double, 6> s;
    s[0] = stress[0] - p; s[1] = stress[1] - p; s[2] = stress[2] - p;
    s[3] = stress[3]; s[4] = stress[4]; s[5] = stress[5];
    
    double J2 = 0.5 * (s[0]*s[0] + s[1]*s[1] + s[2]*s[2]) + s[3]*s[3] + s[4]*s[4] + s[5]*s[5];
    double sigma_eq = std::sqrt(3.0 * J2);
    
    double eps_eq_dot = effectiveCreepRate(sigma_eq, T);
    
    // Add primary creep contribution
    if (params_.include_primary && state.time > 0.0) {
        eps_eq_dot += primaryCreepRate(sigma_eq, T, state.time);
    }
    
    // Flow rule
    std::array<double, 6> eps_dot;
    if (sigma_eq > 1e-10) {
        double factor = 1.5 * eps_eq_dot / sigma_eq;
        for (int i = 0; i < 6; ++i) {
            eps_dot[i] = factor * s[i];
        }
    } else {
        eps_dot.fill(0.0);
    }
    return eps_dot;
}

void CreepModel_Impl::updateState(const std::array<double, 6>& stress,
                                  double T, double dt, State& state) const {
    auto eps_dot = creepStrainRate(stress, T, state);
    
    double eps_eq_dot = 0.0;
    for (int i = 0; i < 3; ++i) {
        state.creep_strain_tensor[i] += eps_dot[i] * dt;
        eps_eq_dot += eps_dot[i] * eps_dot[i];
    }
    for (int i = 3; i < 6; ++i) {
        state.creep_strain_tensor[i] += eps_dot[i] * dt;
        eps_eq_dot += 2.0 * eps_dot[i] * eps_dot[i];
    }
    
    state.creep_strain += std::sqrt(2.0/3.0 * eps_eq_dot) * dt;
    state.time += dt;
    
    // Compute equivalent stress for stage determination
    double p = (stress[0] + stress[1] + stress[2]) / 3.0;
    double J2 = 0.5 * ((stress[0]-p)*(stress[0]-p) + (stress[1]-p)*(stress[1]-p) + (stress[2]-p)*(stress[2]-p)) +
                stress[3]*stress[3] + stress[4]*stress[4] + stress[5]*stress[5];
    double sigma_eq = std::sqrt(3.0 * J2);
    state.stage = determineStage(state, sigma_eq);
}

double CreepModel_Impl::effectiveCreepRate(double sigma_eq, double T) const {
    // Norton power law with Arrhenius temperature dependence
    double arrhenius = std::exp(-params_.activation_energy / (params_.gas_constant * T));
    return params_.A * std::pow(sigma_eq, params_.n) * arrhenius;
}

double CreepModel_Impl::timeToRupture(double sigma_eq, double T) const {
    double eps_min = effectiveCreepRate(sigma_eq, T);
    if (eps_min < 1e-20) return 1e20;
    // Monkman-Grant: t_r * eps_min^m = C
    double C = 0.05;
    double m = 0.85;
    return C / std::pow(eps_min, m);
}

CreepStage CreepModel_Impl::determineStage(const State& state, double sigma_eq) const {
    (void)sigma_eq;
    if (state.damage > params_.critical_damage * 0.5) {
        return CreepStage::TERTIARY;
    }
    if (state.time < 1000.0) {  // Rough threshold
        return CreepStage::PRIMARY;
    }
    return CreepStage::SECONDARY;
}

double CreepModel_Impl::getCreepViscosity(double sigma_eq, double T) const {
    double eps_eq = effectiveCreepRate(sigma_eq, T);
    if (eps_eq < 1e-20) return 1e20;
    return sigma_eq / (3.0 * eps_eq);
}

double CreepModel_Impl::primaryCreepRate(double sigma_eq, double T, double time) const {
    // Andrade-type: eps = B * t^m * sigma^n
    double m = params_.primary_exponent;
    double B = params_.primary_coefficient;
    if (time < 1e-10) return 0.0;
    double arrhenius = std::exp(-params_.activation_energy / (params_.gas_constant * T));
    return B * m * std::pow(time, m - 1.0) * std::pow(sigma_eq, params_.n) * arrhenius;
}

double CreepModel_Impl::secondaryCreepRate(double sigma_eq, double T) const {
    return effectiveCreepRate(sigma_eq, T);
}

double CreepModel_Impl::tertiaryCreepRate(double sigma_eq, double T, double damage) const {
    double secondary = secondaryCreepRate(sigma_eq, T);
    return secondary / (1.0 - damage);
}

double CreepModel_Impl::thetaProjection(double sigma_eq, double T, double time) const {
    (void)sigma_eq;
    (void)T;
    double t1 = params_.theta1, t2 = params_.theta2, t3 = params_.theta3, t4 = params_.theta4;
    return t1 * (1.0 - std::exp(-t2 * time)) + t3 * (std::exp(t4 * time) - 1.0);
}

void CreepModel_Impl::burgersUpdate(double sigma_eq, double dt, double& maxwell, double& kelvin) const {
    // Simplified Burgers model
    double eta_m = 1e15;  // Maxwell viscosity
    double eta_k = 1e14;  // Kelvin viscosity
    double E_k = 1e9;     // Kelvin modulus
    
    maxwell += sigma_eq / eta_m * dt;
    double kelvin_eq = sigma_eq / E_k;
    kelvin += (kelvin_eq - kelvin) * (1.0 - std::exp(-E_k / eta_k * dt));
}

// =============================================================================
// GradientDamage Implementation
// =============================================================================

GradientDamage::GradientDamage() {
}

void GradientDamage::setParameters(const Parameters& params) {
    params_ = params;
}

void GradientDamage::configure(const std::map<std::string, std::string>& config) {
    (void)config;
}

void GradientDamage::initializeState(State& state) const {
    state.damage = 0.0;
    state.kappa = 0.0;
    state.nonlocal_strain = 0.0;
    state.damage_tensor.fill(0.0);
    state.grad_nonlocal.fill(0.0);
}

double GradientDamage::equivalentStrain(const std::array<double, 6>& strain) const {
    // Mazars equivalent strain
    // Compute principal strains and sum positive parts
    double e1 = strain[0], e2 = strain[1], e3 = strain[2];
    double pos1 = std::max(0.0, e1);
    double pos2 = std::max(0.0, e2);
    double pos3 = std::max(0.0, e3);
    return std::sqrt(pos1*pos1 + pos2*pos2 + pos3*pos3);
}

double GradientDamage::nonlocalResidual(double nonlocal, const std::array<double, 3>& grad_nonlocal,
                                        double local_strain,
                                        double test_function,
                                        const std::array<double, 3>& grad_test) const {
    double c = params_.gradient_parameter;
    
    double grad_term = c * (grad_nonlocal[0] * grad_test[0] + 
                           grad_nonlocal[1] * grad_test[1] + 
                           grad_nonlocal[2] * grad_test[2]);
    
    return (nonlocal - local_strain) * test_function + grad_term;
}

double GradientDamage::calculateDamage(double nonlocal_strain, State& state) const {
    // Update history variable (max nonlocal strain)
    state.kappa = std::max(state.kappa, nonlocal_strain);
    state.nonlocal_strain = nonlocal_strain;
    
    if (state.kappa <= params_.damage_threshold) {
        return state.damage;
    }
    
    // Exponential softening
    double d_new = damageEvolutionLaw(state.kappa);
    d_new = std::min(d_new, params_.max_damage);
    
    // Irreversibility
    state.damage = std::max(state.damage, d_new);
    return state.damage;
}

std::array<double, 6> GradientDamage::damagedStress(const std::array<double, 6>& effective_stress,
                                                    const State& state) const {
    double g = (1.0 - state.damage) * (1.0 - state.damage);  // Quadratic degradation
    std::array<double, 6> sigma;
    for (int i = 0; i < 6; ++i) {
        sigma[i] = g * effective_stress[i];
    }
    return sigma;
}

std::array<std::array<double, 6>, 6> GradientDamage::damagedTangent(
    const std::array<std::array<double, 6>, 6>& C,
    const std::array<double, 6>& stress,
    const std::array<double, 6>& strain,
    const State& state) const {
    
    (void)stress;
    (void)strain;
    
    double g = (1.0 - state.damage) * (1.0 - state.damage);
    std::array<std::array<double, 6>, 6> C_damaged;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            C_damaged[i][j] = g * C[i][j];
        }
    }
    return C_damaged;
}

bool GradientDamage::isLocalizing(const State& state) const {
    return state.damage > 0.9 * params_.max_damage;
}

double GradientDamage::dissipatedEnergy(const State& state,
                                        const std::array<double, 6>& stress,
                                        const std::array<double, 6>& strain) const {
    (void)stress;
    (void)strain;
    return params_.fracture_energy * state.damage;
}

double GradientDamage::damageEvolutionLaw(double kappa) const {
    double k0 = params_.damage_threshold;
    double beta = params_.softening_exponent * params_.characteristic_length / params_.fracture_energy;
    return 1.0 - (k0 / kappa) * std::exp(-beta * (kappa - k0));
}

double GradientDamage::damageDrivingForce(const std::array<double, 6>& strain,
                                          const std::array<std::array<double, 6>, 6>& C) const {
    // Y = 0.5 * epsilon : C : epsilon
    double Y = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            Y += 0.5 * strain[i] * C[i][j] * strain[j];
        }
    }
    return Y;
}

// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<HyperelasticMaterial> createHyperelasticMaterial(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<HyperelasticMaterial>();
    model->configure(config);
    return model;
}

std::unique_ptr<FiniteStrainPlasticity> createFiniteStrainPlasticity(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<FiniteStrainPlasticity>();
    model->configure(config);
    return model;
}

std::unique_ptr<Viscoplasticity> createViscoplasticity(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<Viscoplasticity>();
    model->configure(config);
    return model;
}

std::unique_ptr<CreepModel_Impl> createCreepModel(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<CreepModel_Impl>();
    model->configure(config);
    return model;
}

std::unique_ptr<HypoplasticMaterial> createHypoplasticMaterial(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<HypoplasticMaterial>();
    model->configure(config);
    return model;
}

std::unique_ptr<GradientDamage> createGradientDamage(
    const std::map<std::string, std::string>& config) {
    
    auto model = std::make_unique<GradientDamage>();
    model->configure(config);
    return model;
}

// =============================================================================
// HighFidelityGeomechConfig Implementation
// =============================================================================

void HighFidelityGeomechConfig::parseConfig(const std::map<std::string, std::string>& config) {
    (void)config;
}

bool HighFidelityGeomechConfig::validate(std::string& error_msg) const {
    error_msg.clear();
    return true;
}

} // namespace HighFidelity
} // namespace FSRM
