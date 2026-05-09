/**
 * @file MultigroupRadiationDiffusion.cpp
 * @brief Pass-10 multigroup radiation-diffusion solver implementation.
 *        See MultigroupRadiationDiffusion.hpp for the design rationale
 *        and references.
 */

#include "domain/explosion/MultigroupRadiationDiffusion.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace FSRM {

namespace
{

constexpr double FOUR_PI = 12.5663706143591729539;
constexpr double SPEED_OF_LIGHT_M_PER_S =
    RadiationConstants::SPEED_OF_LIGHT_M_PER_S;
constexpr double RADIATION_CONSTANT_A =
    RadiationConstants::RADIATION_CONSTANT_A_J_PER_M3_K4;

double safeMax(double a, double b) { return a > b ? a : b; }

}  // namespace

double MultigroupRadiationDiffusionSolver::cellVolumeSpherical(double r_lo,
                                                               double r_hi)
{
    return (FOUR_PI / 3.0) * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
}

double MultigroupRadiationDiffusionSolver::faceAreaSpherical(double r)
{
    return FOUR_PI * r * r;
}

double MultigroupRadiationDiffusionSolver::harmonicMean(double a, double b)
{
    if (a <= 0.0 || b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

void MultigroupRadiationDiffusionSolver::setConfig(const Config& cfg)
{
    config_ = cfg;
    if (config_.group_grid.n_groups < 1) config_.group_grid.n_groups = 1;
    if (config_.group_grid.n_simpson_points < 3)
        config_.group_grid.n_simpson_points = 3;
    if ((config_.group_grid.n_simpson_points % 2) == 0)
        ++config_.group_grid.n_simpson_points;
    opacity_.setBaselineParameters(config_.opacity_params);
    opacity_.setGrid(config_.group_grid);
    G_ = config_.group_grid.n_groups;
    time_integrator_ =
        DiffusionTimeIntegrator::create(config_.time_integrator);
    prev_step_valid_ = false;
}

void MultigroupRadiationDiffusionSolver::initialize(int N)
{
    N_ = N > 0 ? N : 0;
    if (G_ == 0) G_ = config_.group_grid.n_groups;
    a_lower_.assign(N_, 0.0);
    b_diag_.assign(N_, 0.0);
    c_upper_.assign(N_, 0.0);
    rhs_.assign(N_, 0.0);
    x_.assign(N_, 0.0);
    kappa_R_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    kappa_P_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    D_face_.assign(static_cast<std::size_t>(N_ + 1), 0.0);
    T_iter_.assign(N_, 0.0);
    T_m_old_.assign(N_, 0.0);
    E_iter_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    E_r_old_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    e_int_old_.assign(N_, 0.0);
    B_g_iter_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    dBg_dT_iter_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    B_g_old_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    E_r_prev_step_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    E_r_step_n_.assign(static_cast<std::size_t>(N_) * G_, 0.0);
    prev_step_valid_ = false;
    if (!time_integrator_) {
        time_integrator_ =
            DiffusionTimeIntegrator::create(config_.time_integrator);
    }
}

void MultigroupRadiationDiffusionSolver::recomputeOpacitiesAndPlanck(
    const std::vector<double>& rho)
{
    // Per cell, per group: re-evaluate kappa_R^g, kappa_P^g, B_g, dBg/dT.
    // The B_g and dBg/dT use the matter-temperature iterate T_iter_[i].
    // CONSTANT override (kappa_constant_m2_per_kg > 0) replaces the
    // analytic per-group means with a uniform value; it does not affect
    // B_g(T) or dBg/dT (those drive the source term and must reflect
    // the true Planck integrals).
    const double dT_eps_factor = 0.001;  // 0.1% finite-difference for dBg/dT
    const bool use_const = config_.kappa_constant_m2_per_kg > 0.0;
    const double k_const = config_.kappa_constant_m2_per_kg;
    for (int i = 0; i < N_; ++i) {
        const double T = safeMax(1.0, T_iter_[i]);
        const double dT = safeMax(1.0, dT_eps_factor * T);
        for (int g = 0; g < G_; ++g) {
            const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
            if (use_const) {
                kappa_R_[k] = k_const;
                kappa_P_[k] = k_const;
            } else {
                kappa_R_[k] = opacity_.rosselandPerGroup(g, rho[i], T);
                kappa_P_[k] = opacity_.planckPerGroup(g, rho[i], T);
            }
            const double Bg_T = opacity_.bandIntegratedPlanck(g, T);
            const double Bg_Tp = opacity_.bandIntegratedPlanck(g, T + dT);
            B_g_iter_[k] = Bg_T;
            dBg_dT_iter_[k] = (Bg_Tp - Bg_T) / dT;
        }
    }
}

void MultigroupRadiationDiffusionSolver::updateFaceDiffusion(
    const std::vector<double>& rho, int g)
{
    if (N_ == 0) return;
    std::vector<double> D_cell(N_, 0.0);
    for (int i = 0; i < N_; ++i) {
        const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
        const double denom =
            3.0 * safeMax(1.0e-30, kappa_R_[k]) * safeMax(1.0e-30, rho[i]);
        D_cell[i] = SPEED_OF_LIGHT_M_PER_S / denom;
    }
    D_face_[0] = D_cell[0];
    for (int i = 1; i < N_; ++i) {
        D_face_[i] = harmonicMean(D_cell[i - 1], D_cell[i]);
    }
    D_face_[N_] = D_cell[N_ - 1];
}

void MultigroupRadiationDiffusionSolver::assembleTridiagonal(
    double dt, int g,
    const std::vector<double>& r_cell,
    const std::vector<double>& r_face,
    const std::vector<double>& rho)
{
    // Pass-11 generalisation: parameterise the per-cell-per-group
    // assembly by the DiffusionTimeIntegrator coefficients. BACKWARD_EULER
    // recovers the pass-10 form byte-for-byte. CRANK_NICOLSON splits
    // the operator and source 50/50 between n and n+1; BDF2 reads the
    // prior-prior step E^{n-1}.
    //
    // Note on the multigroup E_r_old_ cache: the pass-10 outer Newton
    // loop mutates E_r_old_ as a Picard relaxation cache between iters.
    // BE under pass-11 must reproduce that bit-exactly, so it continues
    // to read E_r_old_. CN and BDF2 require a stable E^n snapshot, which
    // this solver maintains in E_r_step_n_ (populated once per step()).
    const double four_pi = FOUR_PI;
    const auto coef = time_integrator_
                          ? time_integrator_->getCoefficients()
                          : BackwardEulerIntegrator().getCoefficients();
    const bool needs_explicit_diffusion =
        std::abs(coef.rhs_explicit_diffusion_factor) > 0.0;
    const bool needs_explicit_source =
        std::abs(coef.rhs_source_explicit_factor) > 0.0;
    const bool needs_prev_state =
        std::abs(coef.rhs_volume_nm1_factor) > 0.0;
    const bool bootstrap_to_be = needs_prev_state && !prev_step_valid_;
    const auto eff = bootstrap_to_be
                         ? BackwardEulerIntegrator().getCoefficients()
                         : coef;
    const bool is_be =
        time_integrator_ &&
        time_integrator_->kind() ==
            DiffusionTimeIntegratorKind::BACKWARD_EULER;
    const std::vector<double>& E_n_buf =
        (is_be || bootstrap_to_be) ? E_r_old_ : E_r_step_n_;

    for (int i = 0; i < N_; ++i) {
        const double r_lo = r_face[i];
        const double r_hi = r_face[i + 1];
        const double V_i = cellVolumeSpherical(r_lo, r_hi);

        double a_e = 0.0;
        if (i < N_ - 1) {
            const double A_e = faceAreaSpherical(r_hi);
            const double dr_e = safeMax(1.0e-12, r_cell[i + 1] - r_cell[i]);
            a_e = A_e * D_face_[i + 1] / dr_e;
        }

        double a_w = 0.0;
        if (i > 0) {
            const double A_w = faceAreaSpherical(r_lo);
            const double dr_w = safeMax(1.0e-12, r_cell[i] - r_cell[i - 1]);
            a_w = A_w * D_face_[i] / dr_w;
        }

        const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
        const double c_kp_rho_V =
            V_i * SPEED_OF_LIGHT_M_PER_S * kappa_P_[k] * rho[i];
        const double S_emit_implicit =
            (four_pi * B_g_iter_[k]) / SPEED_OF_LIGHT_M_PER_S;
        double S_emit_explicit = 0.0;
        if (needs_explicit_source) {
            S_emit_explicit = (four_pi * B_g_old_[k]) / SPEED_OF_LIGHT_M_PER_S;
        }

        a_lower_[i] = -dt * a_w * eff.lhs_implicit_factor;
        c_upper_[i] = -dt * a_e * eff.lhs_implicit_factor;
        b_diag_[i] = V_i * eff.lhs_volume_factor +
                     dt * (a_e + a_w + c_kp_rho_V) * eff.lhs_implicit_factor;

        const double E_old = E_n_buf[k];
        double rhs_cell = V_i * eff.rhs_volume_n_factor * E_old;
        if (needs_prev_state && !bootstrap_to_be) {
            rhs_cell += V_i * eff.rhs_volume_nm1_factor * E_r_prev_step_[k];
        }
        if (needs_explicit_diffusion) {
            const std::size_t kL = (i > 0)
                ? static_cast<std::size_t>(i - 1) * G_ + g
                : k;
            const std::size_t kR = (i < N_ - 1)
                ? static_cast<std::size_t>(i + 1) * G_ + g
                : k;
            const double E_left  = (i > 0)      ? E_n_buf[kL] : E_old;
            const double E_right = (i < N_ - 1) ? E_n_buf[kR] : E_old;
            const double explicit_diff =
                a_e * (E_right - E_old) -
                a_w * (E_old - E_left) -
                c_kp_rho_V * E_old;
            rhs_cell += dt * eff.rhs_explicit_diffusion_factor * explicit_diff;
        }
        rhs_cell += dt * c_kp_rho_V *
                    (eff.rhs_source_implicit_factor * S_emit_implicit +
                     eff.rhs_source_explicit_factor * S_emit_explicit);
        rhs_[i] = rhs_cell;
    }

    // Outer Marshak / extrapolation BC. Per-group ghost = 4 pi B_g(T_amb).
    if (N_ >= 2) {
        const double T_amb = config_.T_ambient_K;
        const double B_amb = opacity_.bandIntegratedPlanck(g, T_amb);
        const double E_outer = (FOUR_PI * B_amb) / SPEED_OF_LIGHT_M_PER_S;
        const double r_last = r_face[N_];
        const double A_last = faceAreaSpherical(r_last);
        const double dr_ghost = safeMax(1.0e-3, r_face[N_] - r_cell[N_ - 1]);
        const double D_ghost = D_face_[N_];
        const double a_ghost = A_last * D_ghost / dr_ghost;
        b_diag_[N_ - 1] += dt * a_ghost * eff.lhs_implicit_factor;
        rhs_[N_ - 1] += dt * a_ghost * E_outer *
                        (eff.lhs_implicit_factor +
                         eff.rhs_explicit_diffusion_factor);
        if (std::abs(eff.rhs_explicit_diffusion_factor) > 0.0) {
            const std::size_t kLast =
                static_cast<std::size_t>(N_ - 1) * G_ + g;
            rhs_[N_ - 1] -= dt * a_ghost *
                            eff.rhs_explicit_diffusion_factor *
                            E_n_buf[kLast];
        }
    }
}

void MultigroupRadiationDiffusionSolver::solveTridiagonal()
{
    if (N_ == 0) return;
    for (int i = 1; i < N_; ++i) {
        const double m = a_lower_[i] / safeMax(1.0e-30, b_diag_[i - 1]);
        b_diag_[i] -= m * c_upper_[i - 1];
        rhs_[i] -= m * rhs_[i - 1];
    }
    x_[N_ - 1] = rhs_[N_ - 1] / safeMax(1.0e-30, b_diag_[N_ - 1]);
    for (int i = N_ - 2; i >= 0; --i) {
        x_[i] = (rhs_[i] - c_upper_[i] * x_[i + 1]) /
                safeMax(1.0e-30, b_diag_[i]);
    }
}

void MultigroupRadiationDiffusionSolver::updateMatterTemperature(
    double dt, const std::vector<double>& rho, const TillotsonEOS& eos)
{
    (void)rho;
    // Per-cell linearised matter-energy update with multigroup source.
    //   rho cv (T^{n+1} - T^n)/dt
    //     = sum_g c kappa_P^g rho ( E_g^{n+1} - 4 pi B_g(T^{n+1}) / c )
    // Linearise B_g(T^{n+1}) ~ B_g(T_iter) + dB_g/dT (T^{n+1} - T_iter)
    // and solve for T^{n+1}:
    //   T^{n+1} (1/dt + sum_g (c kappa_P^g / cv) * 4 pi (dB_g/dT)/c)
    //     = T^n / dt + sum_g (c kappa_P^g / cv) * E_g^{n+1}
    //                  - sum_g (c kappa_P^g / cv) * 4 pi B_g(T_iter)/c
    //                  + sum_g (c kappa_P^g / cv) * 4 pi (dB_g/dT)/c * T_iter
    const double cv = safeMax(1.0e-3, eos.getParameters().cv);
    const double four_pi = FOUR_PI;
    const double c = SPEED_OF_LIGHT_M_PER_S;
    for (int i = 0; i < N_; ++i) {
        const double T_n = T_m_old_[i];
        const double T_iter = T_iter_[i];
        double sum_alpha = 0.0;     // T^{n+1} coefficient on LHS
        double sum_rhs = 0.0;       // RHS contribution from groups
        for (int g = 0; g < G_; ++g) {
            const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
            const double a_g = (c * kappa_P_[k] / cv);
            const double Bg = B_g_iter_[k];
            const double dBg = dBg_dT_iter_[k];
            const double E_g = E_iter_[k];
            const double Sg = (four_pi * Bg) / c;          // emission scale
            const double dSg = (four_pi * dBg) / c;        // d/dT of emission
            sum_alpha += a_g * dSg;
            sum_rhs += a_g * E_g - a_g * Sg + a_g * dSg * T_iter;
        }
        const double rhs_cell = T_n / dt + sum_rhs;
        const double denom = 1.0 / dt + sum_alpha;
        const double T_new = rhs_cell / safeMax(1.0e-30, denom);
        T_iter_[i] = T_new > 1.0 ? T_new : 1.0;
    }
}

MultigroupRadiationDiffusionSolver::StepResult
MultigroupRadiationDiffusionSolver::step(
    double dt, const std::vector<double>& r_cell,
    const std::vector<double>& r_face, const std::vector<double>& rho,
    std::vector<double>& e_int, std::vector<double>& E_r,
    std::vector<double>& T_m, const TillotsonEOS& eos,
    const std::vector<int>& is_gas)
{
    StepResult res;
    if (N_ == 0 || G_ == 0) {
        res.converged = false;
        return res;
    }
    if (static_cast<int>(E_r.size()) < N_ * G_) {
        throw std::runtime_error(
            "MultigroupRadiationDiffusionSolver::step: E_r vector size "
            "smaller than N*G; caller must allocate N*G entries.");
    }

    // Snapshot the old state. E_r_old_ is the pass-10 Picard relaxation
    // cache (mutated each Newton iter); E_r_step_n_ is the stable
    // pass-11 E^n snapshot used by the CN / BDF2 explicit terms.
    for (int i = 0; i < N_; ++i) {
        T_m_old_[i] = T_m[i] > 1.0 ? T_m[i] : config_.T_ambient_K;
        T_iter_[i] = T_m_old_[i];
        e_int_old_[i] = e_int[i];
        for (int g = 0; g < G_; ++g) {
            const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
            E_r_old_[k] = E_r[k];
            E_iter_[k] = E_r[k];
            E_r_step_n_[k] = E_r[k];
        }
    }

    // Pass-11: cache B_g(T_m_old) once for the Crank-Nicolson explicit
    // source weight. Backward-Euler and BDF2 ignore this buffer.
    const bool needs_explicit_source =
        time_integrator_ &&
        std::abs(time_integrator_->getCoefficients()
                 .rhs_source_explicit_factor) > 0.0;
    if (needs_explicit_source) {
        for (int i = 0; i < N_; ++i) {
            const double T_old = T_m_old_[i];
            for (int g = 0; g < G_; ++g) {
                const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
                B_g_old_[k] = opacity_.bandIntegratedPlanck(g, T_old);
            }
        }
    }

    int iter = 0;
    double residual = 1.0e30;
    for (; iter < config_.max_newton_iter; ++iter) {
        // Step 1: refresh per-cell, per-group opacities and Planck source
        // at the current T_iter.
        recomputeOpacitiesAndPlanck(rho);

        // Step 2: per group, assemble + solve the tridiagonal for E_g^{n+1}.
        for (int g = 0; g < G_; ++g) {
            updateFaceDiffusion(rho, g);
            assembleTridiagonal(dt, g, r_cell, r_face, rho);
            solveTridiagonal();
            for (int i = 0; i < N_; ++i) {
                const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
                E_iter_[k] = x_[i];
            }
        }

        // Step 3: update matter temperature from the new E_g.
        updateMatterTemperature(dt, rho, eos);

        // Step 4: residual on E (max relative change across all cell-group
        // entries from the previous iteration).
        residual = 0.0;
        for (int i = 0; i < N_; ++i) {
            for (int g = 0; g < G_; ++g) {
                const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
                const double Eold = E_r_old_[k];
                const double Enew = E_iter_[k];
                const double scale =
                    safeMax(1.0e-12, std::abs(Enew) + std::abs(Eold));
                const double rel = std::abs(Enew - Eold) / scale;
                if (rel > residual) residual = rel;
                E_r_old_[k] = Enew;  // for next-iter linearisation cache
            }
        }
        if (residual < config_.newton_tolerance) {
            ++iter;
            break;
        }
    }

    // Apply the converged state to the host vectors. Compute matter
    // energy change as cv (T_new - T_old).
    const double cv = safeMax(1.0e-3, eos.getParameters().cv);
    double matter_dE_total = 0.0;
    for (int i = 0; i < N_; ++i) {
        const double T_new = T_iter_[i];
        const double de = cv * (T_new - T_m_old_[i]);
        T_m[i] = T_new;
        e_int[i] += de;
        const double V_i = cellVolumeSpherical(r_face[i], r_face[i + 1]);
        matter_dE_total += rho[i] * de * V_i;
        for (int g = 0; g < G_; ++g) {
            const std::size_t k = static_cast<std::size_t>(i) * G_ + g;
            E_r[k] = E_iter_[k];
        }
        (void)is_gas;
    }

    // Diagnostic: radiation front index based on summed E_r_total per cell.
    // Front is the outermost cell where sum_g E_r^g exceeds front_factor *
    // a T_amb^4 (the same threshold the grey solver uses for the gray E_r).
    const double T_amb = config_.T_ambient_K;
    const double T_amb4 = T_amb * T_amb * T_amb * T_amb;
    const double E_floor = config_.front_factor * RADIATION_CONSTANT_A * T_amb4;
    int front_idx = 0;
    for (int i = N_ - 1; i >= 0; --i) {
        double E_sum = 0.0;
        for (int g = 0; g < G_; ++g) {
            E_sum += E_r[static_cast<std::size_t>(i) * G_ + g];
        }
        if (E_sum > E_floor) { front_idx = i; break; }
    }
    double t_diff = 0.0;
    if (front_idx >= 0 && front_idx < N_) {
        // Use the harmonic-average opacity across all groups (last group's
        // updateFaceDiffusion left D_face_ at the last group; we recompute
        // the cell's effective gray D from the per-group Rosseland mean).
        double sum_kappa_R = 0.0;
        for (int g = 0; g < G_; ++g) {
            sum_kappa_R +=
                kappa_R_[static_cast<std::size_t>(front_idx) * G_ + g];
        }
        const double kappa_R_eff = sum_kappa_R / G_;
        const double D_eff =
            SPEED_OF_LIGHT_M_PER_S /
            (3.0 * safeMax(1.0e-30, kappa_R_eff) * safeMax(1.0e-30, rho[front_idx]));
        const double dr =
            safeMax(1.0e-9, r_face[front_idx + 1] - r_face[front_idx]);
        t_diff = dr * dr / D_eff;
    }

    double E_r_total = 0.0;
    for (int i = 0; i < N_; ++i) {
        const double V_i = cellVolumeSpherical(r_face[i], r_face[i + 1]);
        for (int g = 0; g < G_; ++g) {
            E_r_total +=
                E_r[static_cast<std::size_t>(i) * G_ + g] * V_i;
        }
    }

    res.converged = (residual < config_.newton_tolerance);
    res.newton_iters = iter;
    res.residual_inf_norm = residual;
    res.radiation_front_index = front_idx;
    res.radiation_front_radius_m = (front_idx < N_) ? r_cell[front_idx] : 0.0;
    res.t_diff_at_front_s = t_diff;
    res.total_radiation_energy_J = E_r_total;
    res.total_matter_energy_change_J = matter_dE_total;

    // Pass-11: shift prior-step state for the next BDF2 advance. The
    // value that was E^n on this call becomes E^{n-1} on the next call.
    for (std::size_t k = 0; k < E_r_prev_step_.size(); ++k) {
        E_r_prev_step_[k] = E_r_step_n_[k];
    }
    if (time_integrator_ && time_integrator_->needsTwoPriorStates()) {
        prev_step_valid_ = true;
    }
    return res;
}

}  // namespace FSRM
