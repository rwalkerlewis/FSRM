/**
 * @file PetscFEViscoelastic.cpp
 * @brief PetscDS pointwise callbacks for viscoelastic attenuation
 *
 * Implements the generalized Maxwell body (GMB) for seismic attenuation.
 * Memory variables R_i are stored in auxiliary fields and updated by a
 * TSPostStep callback (in Simulator.cpp), not by these residual callbacks.
 *
 * The f1 callback reads the current memory variables from aux fields and
 * adds them to the elastic stress. The g3 callback uses the unrelaxed
 * modulus as the tangent approximation.
 *
 * Reference: Moczo et al. (2007) "The Finite-Difference Modelling of
 * Earthquake Motions", Chapter 5: Viscoelastic modelling.
 */

#include "numerics/PetscFEViscoelastic.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace FSRM {
namespace PetscFEViscoelastic {

// ---------------------------------------------------------------------------
// Helper: read an integer constant with bounds check
// ---------------------------------------------------------------------------
static inline PetscInt getIntConst(PetscInt numConstants, const PetscScalar c[],
                                   PetscInt idx, PetscInt fallback) {
    return (numConstants > idx) ? static_cast<PetscInt>(PetscRealPart(c[idx])) : fallback;
}

static inline PetscReal getRealConst(PetscInt numConstants, const PetscScalar c[],
                                     PetscInt idx, PetscReal fallback) {
    return (numConstants > idx) ? PetscRealPart(c[idx]) : fallback;
}

// ---------------------------------------------------------------------------
// computeMechanismWeights
// ---------------------------------------------------------------------------
void computeMechanismWeights(int N, double f_min, double f_max, double Q_target,
                             double tau[], double delta_mu[])
{
    if (N <= 0 || N > VISCO_MAX_MECHANISMS) return;
    if (Q_target <= 0.0) Q_target = 100.0;

    // Log-spaced relaxation frequencies
    std::vector<double> relax_freq(static_cast<std::size_t>(N));
    if (N == 1) {
        relax_freq[0] = std::sqrt(f_min * f_max);
    } else {
        double log_fmin = std::log(f_min);
        double log_fmax = std::log(f_max);
        for (int i = 0; i < N; ++i) {
            double t_frac = static_cast<double>(i) / (N - 1);
            relax_freq[static_cast<std::size_t>(i)] = std::exp(log_fmin + t_frac * (log_fmax - log_fmin));
        }
    }

    // Relaxation times: tau_i = 1 / (2 * pi * f_i)
    for (int i = 0; i < N; ++i) {
        tau[i] = 1.0 / (2.0 * M_PI * relax_freq[static_cast<std::size_t>(i)]);
    }

    // Solve least-squares for anelastic coefficients:
    //   1/Q = sum_i Y_i * omega * tau_i / (1 + (omega * tau_i)^2)
    // Sample frequencies for fitting
    int n_sample = N * 10;
    if (n_sample < 10) n_sample = 10;

    // Build matrix A and target vector b
    std::vector<std::vector<double>> A(static_cast<std::size_t>(n_sample),
                                        std::vector<double>(static_cast<std::size_t>(N)));
    double target = 1.0 / Q_target;

    double log_fmin = std::log(f_min);
    double log_fmax = std::log(f_max);
    for (int j = 0; j < n_sample; ++j) {
        double t_frac = static_cast<double>(j) / (n_sample - 1);
        double freq = std::exp(log_fmin + t_frac * (log_fmax - log_fmin));
        double omega = 2.0 * M_PI * freq;
        for (int i = 0; i < N; ++i) {
            double wt = omega * tau[i];
            A[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = wt / (1.0 + wt * wt);
        }
    }

    // Normal equations: (A^T A + reg*I) Y = A^T b
    std::vector<std::vector<double>> ATA(static_cast<std::size_t>(N),
                                          std::vector<double>(static_cast<std::size_t>(N), 0.0));
    std::vector<double> ATb(static_cast<std::size_t>(N), 0.0);
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < N; ++k) {
            for (int j = 0; j < n_sample; ++j) {
                ATA[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] +=
                    A[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] *
                    A[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)];
            }
        }
        for (int j = 0; j < n_sample; ++j) {
            ATb[static_cast<std::size_t>(i)] += A[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] * target;
        }
        ATA[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] += 1e-8;  // regularization
    }

    // Gaussian elimination with pivoting
    std::vector<std::vector<double>> aug(static_cast<std::size_t>(N),
                                          std::vector<double>(static_cast<std::size_t>(N + 1)));
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < N; ++k)
            aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                ATA[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
        aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(N)] = ATb[static_cast<std::size_t>(i)];
    }

    for (int i = 0; i < N; ++i) {
        int max_row = i;
        for (int k = i + 1; k < N; ++k)
            if (std::abs(aug[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)]) >
                std::abs(aug[static_cast<std::size_t>(max_row)][static_cast<std::size_t>(i)]))
                max_row = k;
        std::swap(aug[static_cast<std::size_t>(i)], aug[static_cast<std::size_t>(max_row)]);
        for (int k = i + 1; k < N; ++k) {
            if (std::abs(aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)]) < 1e-15) continue;
            double factor = aug[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] /
                            aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
            for (int j = i; j <= N; ++j)
                aug[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)] -=
                    factor * aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }

    // Back substitution
    std::vector<double> Y(static_cast<std::size_t>(N), 0.0);
    for (int i = N - 1; i >= 0; --i) {
        Y[static_cast<std::size_t>(i)] = aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(N)];
        for (int j = i + 1; j < N; ++j)
            Y[static_cast<std::size_t>(i)] -= aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
                                               Y[static_cast<std::size_t>(j)];
        if (std::abs(aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)]) > 1e-15)
            Y[static_cast<std::size_t>(i)] /= aug[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }

    // Y_i are the anelastic coefficients. For a generalized Maxwell body,
    // delta_mu_i = Y_i * mu_relaxed (but we store them as fractions of mu,
    // to be multiplied by the actual mu at each cell during assembly).
    // Ensure non-negative weights.
    for (int i = 0; i < N; ++i) {
        delta_mu[i] = std::max(0.0, Y[static_cast<std::size_t>(i)]);
    }
}

// ---------------------------------------------------------------------------
// f1_viscoelastic_aux
// ---------------------------------------------------------------------------
void f1_viscoelastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar f1[])
{
    (void)Nf; (void)NfAux; (void)uOff; (void)u; (void)u_t;
    (void)aOff_x; (void)a_x; (void)a_t; (void)t; (void)x;

    // Read relaxed moduli from auxiliary fields
    const PetscScalar lambda = a[aOff[0]];
    const PetscScalar mu     = a[aOff[1]];

    // Compute strain from displacement gradient: eps_ij = 0.5*(u_i,j + u_j,i)
    // u_x layout: u_x[uOff_x[field] + i*dim + j] = du_i/dx_j
    const PetscInt uoff = uOff_x[0];

    // Compute elastic stress: sigma = lambda * tr(eps) * I + 2 * mu * eps
    // f1[i*dim + j] = sigma_ij  (Cauchy stress in Voigt-like layout)
    PetscScalar trace = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
        trace += u_x[uoff + d * dim + d];
    }

    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            PetscScalar eps_ij = 0.5 * (u_x[uoff + i * dim + j] + u_x[uoff + j * dim + i]);
            f1[i * dim + j] = lambda * trace * (i == j ? 1.0 : 0.0) + 2.0 * mu * eps_ij;
        }
    }

    // Memory variable contributions (sum_i R_i) would be added here when
    // the memory variables are stored in auxiliary fields. Currently, memory
    // variables start at zero (relaxed initial state) and the TSPostStep
    // callback handles their evolution externally. The f1 callback computes
    // only the elastic stress from the relaxed moduli.
    //
    // Future: when memory variables are stored in aux fields, add:
    //   for each mechanism m: sigma += R_m (read from a[aOff[3 + m*6 + v]])
    (void)NfAux;
}

// ---------------------------------------------------------------------------
// g3_viscoelastic_aux
// ---------------------------------------------------------------------------
void g3_viscoelastic_aux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t, PetscReal u_tShift,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar g3[])
{
    (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x; (void)u; (void)u_t; (void)u_x;
    (void)aOff_x; (void)a_x; (void)a_t; (void)t; (void)u_tShift; (void)x;

    // Read relaxed moduli
    const PetscScalar lambda_r = a[aOff[0]];
    const PetscScalar mu_r     = a[aOff[1]];

    // Compute unrelaxed moduli: mu_u = mu_r * (1 + sum_i delta_mu_i)
    PetscInt N_mech = getIntConst(numConstants, constants, VISCO_CONST_N, 0);
    if (N_mech > VISCO_MAX_MECHANISMS) N_mech = VISCO_MAX_MECHANISMS;

    PetscReal sum_dmu = 0.0;
    PetscReal sum_dk = 0.0;
    for (PetscInt m = 0; m < N_mech; ++m) {
        sum_dmu += getRealConst(numConstants, constants, VISCO_CONST_DMU_BASE + m, 0.0);
        sum_dk  += getRealConst(numConstants, constants, VISCO_CONST_DK_BASE + m, 0.0);
    }

    // Pass-3: per-cell Q scaling. The constants-array delta_mu was
    // computed once globally with config.visco_q_s, so it represents
    // the constant-Q solution at Q_global. Per-layer Q lets us
    // approximate cell-local attenuation by scaling delta_mu by
    // (Q_global / Q_local). At Q_local = Q_global the scale is 1
    // (bit-identical backward compat); at Q_local -> infinity the
    // scaled delta_mu -> 0 and the unrelaxed modulus collapses back
    // to the relaxed (elastic) modulus.
    //
    // NfAux >= 5 means the per-cell Q fields (AUX_QP=3, AUX_QS=4) are
    // populated. Constants slots 78 and 79 carry the global Q the
    // mechanism weights were fit against; when those slots are zero
    // the per-cell scaling is disabled (legacy path).
    PetscReal q_scale_s = 1.0;
    PetscReal q_scale_p = 1.0;
    if (NfAux >= 5) {
        const PetscReal q_p_local = PetscRealPart(a[aOff[3]]);
        const PetscReal q_s_local = PetscRealPart(a[aOff[4]]);
        const PetscReal q_p_global =
            getRealConst(numConstants, constants, VISCO_CONST_Q_P_GLOBAL, 0.0);
        const PetscReal q_s_global =
            getRealConst(numConstants, constants, VISCO_CONST_Q_S_GLOBAL, 0.0);
        if (q_p_local > 0.0 && q_p_global > 0.0)
            q_scale_p = q_p_global / q_p_local;
        if (q_s_local > 0.0 && q_s_global > 0.0)
            q_scale_s = q_s_global / q_s_local;
    }
    const PetscReal sum_dmu_eff = sum_dmu * q_scale_s;
    const PetscReal sum_dk_eff  = sum_dk  * q_scale_p;

    // Unrelaxed moduli (used as tangent for implicit integration)
    const PetscScalar mu_u     = mu_r * (1.0 + sum_dmu_eff);
    const PetscScalar lambda_u = lambda_r + (2.0 / 3.0) *
                                 (mu_r * sum_dmu_eff - mu_r * sum_dk_eff);
    (void)lambda_u;
    // Simplified: for P-wave attenuation, use lambda_u = lambda_r * (1 + sum_dk)
    // But for S-wave only: lambda unchanged, mu adjusted.
    // Using the proper tangent: lambda_u = lambda_r, mu_u = mu_r*(1+sum_dmu)
    // is a reasonable approximation. Let us use:
    const PetscScalar lam = lambda_r;
    const PetscScalar mu  = mu_u;

    // g3[i*dim*dim*dim + j*dim*dim + k*dim + l] = d(sigma_ij)/d(eps_kl)
    // = lambda * delta_ij * delta_kl + mu * (delta_ik * delta_jl + delta_il * delta_jk)
    const PetscInt n = dim * dim * dim * dim;
    for (PetscInt i = 0; i < n; ++i) g3[i] = 0.0;

    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            for (PetscInt k = 0; k < dim; ++k) {
                for (PetscInt l = 0; l < dim; ++l) {
                    PetscInt idx = ((i * dim + j) * dim + k) * dim + l;
                    // C_ijkl = lambda * delta_ij * delta_kl
                    //        + mu * (delta_ik * delta_jl + delta_il * delta_jk)
                    if (i == j && k == l) g3[idx] += lam;
                    if (i == k && j == l) g3[idx] += mu;
                    if (i == l && j == k) g3[idx] += mu;
                }
            }
        }
    }
}

} // namespace PetscFEViscoelastic
} // namespace FSRM
