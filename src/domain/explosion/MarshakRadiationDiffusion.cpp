/**
 * @file MarshakRadiationDiffusion.cpp
 * @brief Implementation of the 1D radial grey radiation-diffusion solver.
 *        See MarshakRadiationDiffusion.hpp for the design rationale and
 *        references.
 */

#include "domain/explosion/MarshakRadiationDiffusion.hpp"

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

} // namespace

double MarshakRadiationDiffusionSolver::cellVolumeSpherical(double r_lo,
                                                            double r_hi)
{
    return (FOUR_PI / 3.0) * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
}

double MarshakRadiationDiffusionSolver::faceAreaSpherical(double r)
{
    return FOUR_PI * r * r;
}

double MarshakRadiationDiffusionSolver::harmonicMean(double a, double b)
{
    if (a <= 0.0 || b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

void MarshakRadiationDiffusionSolver::initialize(int N)
{
    if (config_.opacity_model == OpacityModel::TABULATED_FULL) {
        throw std::runtime_error(
            "MarshakRadiationDiffusionSolver: opacity_model=TABULATED_FULL "
            "is pass-10 work; pass-9 ships TABULATED_PATCHED only.");
    }
    if (config_.opacity_model == OpacityModel::TABULATED_TOPS) {
        // Pass-8 retained TABULATED_TOPS as a throwing stub; pass-9
        // promotes TABULATED_PATCHED to a working option but TOPS
        // remains a named-only legacy alias for any external code that
        // still references it.
        throw std::runtime_error(
            "MarshakRadiationDiffusionSolver: opacity_model=TABULATED_TOPS "
            "is the pass-8 legacy throwing stub; use TABULATED_PATCHED in "
            "pass-9.");
    }

    // Pass-9: lazy-load tabulated opacity tables when configured. Both
    // Rosseland and Planck are needed to apply the patch; if only one
    // path is configured, log a warning and degrade to power-law for
    // both means (the other has no analytic counterpart in our
    // power-law framework that the patch could blend with).
    if (!opacity_table_load_attempted_ &&
        config_.opacity_model == OpacityModel::TABULATED_PATCHED) {
        opacity_table_load_attempted_ = true;
        const bool both_paths =
            !config_.tabulated_opacity_rosseland_path.empty() &&
            !config_.tabulated_opacity_planck_path.empty();
        if (!both_paths) {
            std::fprintf(stderr,
                "MarshakRadiationDiffusionSolver: TABULATED_PATCHED "
                "requires both rosseland and planck table paths; one or "
                "both are empty. Falling back to POWER_LAW_ZR.\n");
        } else {
            std::string err_R, err_P;
            const bool rok = rosseland_table_.load(
                config_.tabulated_opacity_rosseland_path, err_R);
            const bool pok = planck_table_.load(
                config_.tabulated_opacity_planck_path, err_P);
            if (!rok) {
                std::fprintf(stderr,
                    "MarshakRadiationDiffusionSolver: Rosseland table "
                    "load FAILED: %s. Falling back to POWER_LAW_ZR for "
                    "Rosseland.\n", err_R.c_str());
            }
            if (!pok) {
                std::fprintf(stderr,
                    "MarshakRadiationDiffusionSolver: Planck table "
                    "load FAILED: %s. Falling back to POWER_LAW_ZR for "
                    "Planck.\n", err_P.c_str());
            }
        }
    }

    N_ = N > 0 ? N : 0;
    a_lower_.assign(N_, 0.0);
    b_diag_.assign(N_, 0.0);
    c_upper_.assign(N_, 0.0);
    rhs_.assign(N_, 0.0);
    x_.assign(N_, 0.0);
    T_iter_.assign(N_, 0.0);
    E_iter_.assign(N_, 0.0);
    kappa_R_.assign(N_, 0.0);
    kappa_P_.assign(N_, 0.0);
    D_face_.assign(N_ + 1, 0.0);
    e_int_old_.assign(N_, 0.0);
    E_r_old_.assign(N_, 0.0);
    T_m_old_.assign(N_, 0.0);
    power_law_.setParameters(config_.opacity_params);
}

void MarshakRadiationDiffusionSolver::evaluateOpacity(double rho, double T,
                                                     double& kappa_R,
                                                     double& kappa_P) const
{
    evaluateOpacityPatched(rho, T, kappa_R, kappa_P);
}

void MarshakRadiationDiffusionSolver::evaluateOpacityPatched(
    double rho, double T, double& kappa_R, double& kappa_P) const
{
    switch (config_.opacity_model) {
    case OpacityModel::CONSTANT: {
        const double k = config_.kappa_constant_m2_per_kg > 0.0
                             ? config_.kappa_constant_m2_per_kg
                             : 1.0;
        kappa_R = k;
        kappa_P = k;
        return;
    }
    case OpacityModel::POWER_LAW_ZR: {
        PowerLawOpacity ev(config_.opacity_params);
        kappa_R = ev.rosseland(rho, T);
        kappa_P = ev.planck(rho, T);
        return;
    }
    case OpacityModel::TABULATED_PATCHED: {
        // Pass-9 sin^2 patch in temperature. Z-R for T < lower; tabulated
        // for T > upper; smooth blend in between. Out-of-table queries
        // (rho or T outside the table coverage) fall back to Z-R with
        // the reader's one-time warning.
        PowerLawOpacity ev(config_.opacity_params);
        const double k_R_pl = ev.rosseland(rho, T);
        const double k_P_pl = ev.planck(rho, T);
        const double T_lo = config_.tabulated_blend_lower_k;
        const double T_hi = config_.tabulated_blend_upper_k;
        if (T <= T_lo || !rosseland_table_.isLoaded() ||
            !planck_table_.isLoaded()) {
            kappa_R = k_R_pl;
            kappa_P = k_P_pl;
            return;
        }
        const double k_R_tab_raw = rosseland_table_.evaluate(rho, T);
        const double k_P_tab_raw = planck_table_.evaluate(rho, T);
        if (std::isnan(k_R_tab_raw) || std::isnan(k_P_tab_raw)) {
            // Out-of-table coverage; reader has logged the warning.
            kappa_R = k_R_pl;
            kappa_P = k_P_pl;
            return;
        }
        if (T >= T_hi) {
            kappa_R = k_R_tab_raw;
            kappa_P = k_P_tab_raw;
            return;
        }
        // Sin^2 blend in T.
        const double t = (T - T_lo) / (T_hi - T_lo);
        const double s = std::sin(0.5 * 3.14159265358979323846 * t);
        const double w = s * s;
        kappa_R = (1.0 - w) * k_R_pl + w * k_R_tab_raw;
        kappa_P = (1.0 - w) * k_P_pl + w * k_P_tab_raw;
        return;
    }
    case OpacityModel::TABULATED_FULL:
    case OpacityModel::TABULATED_TOPS:
    default:
        // Should be caught by initialize(); zero opacity here is a
        // diagnostic for "fell through unexpectedly".
        kappa_R = 0.0;
        kappa_P = 0.0;
        return;
    }
}

void MarshakRadiationDiffusionSolver::updateOpacities(
    const std::vector<double>& rho, const std::vector<double>& T_iter)
{
    for (int i = 0; i < N_; ++i) {
        double kr = 0.0, kp = 0.0;
        evaluateOpacity(rho[i], T_iter[i], kr, kp);
        kappa_R_[i] = kr;
        kappa_P_[i] = kp;
    }
}

void MarshakRadiationDiffusionSolver::updateFaceDiffusionCoefficient(
    const std::vector<double>& rho)
{
    // Cell-wise D_i = c / (3 kappa_R_i rho_i). Face value D_{i+1/2}
    // uses harmonic mean for proper diffusive flux conservation.
    if (N_ == 0) return;
    std::vector<double> D_cell(N_, 0.0);
    for (int i = 0; i < N_; ++i) {
        const double denom = 3.0 * safeMax(1.0e-30, kappa_R_[i]) *
                             safeMax(1.0e-30, rho[i]);
        D_cell[i] = SPEED_OF_LIGHT_M_PER_S / denom;
    }
    D_face_[0] = D_cell[0];
    for (int i = 1; i < N_; ++i) {
        D_face_[i] = harmonicMean(D_cell[i - 1], D_cell[i]);
    }
    D_face_[N_] = D_cell[N_ - 1];
}

void MarshakRadiationDiffusionSolver::assembleTridiagonal(
    double dt, const std::vector<double>& r_cell,
    const std::vector<double>& r_face, const std::vector<double>& rho)
{
    // Backward-Euler in E_r with linearised source term.
    //
    // Per-cell volume V_i = (4/3) pi (r_face[i+1]^3 - r_face[i]^3).
    // Face flux from cell i to i+1:
    //   F_{i+1/2} = - A_{i+1/2} * D_{i+1/2} * (E_{i+1} - E_i) /
    //               (r_cell[i+1] - r_cell[i])
    // Conservative discretisation of the divergence: per-cell change
    //   V_i d E_i / dt = - (F_{i+1/2} - F_{i-1/2}) + V_i S_i
    // where the source S_i = c kappa_P rho (a T^4 - E_r) is linearised
    // using d(T^4)/dT = 4 T^3 around T_iter; here we use the simpler
    // Picard linearisation: hold T_iter fixed in the source, treat
    // E^{n+1} implicitly. The outer Newton loop iterates T_iter.
    //
    // The matter equation is per-cell so it does not contribute to
    // the tridiagonal stencil. Substituting the linearised T_m^{n+1}
    // from the matter equation into the source for E_r tightens the
    // coupling (Larsen 1988). Pass-8 ships the simpler Picard split:
    // E_r solve with frozen T_iter, then T_m update from the new E_r,
    // outer Newton on the residual. Adequate for the moderate-stiffness
    // regimes pass-8 targets.
    for (int i = 0; i < N_; ++i) {
        const double r_lo = r_face[i];
        const double r_hi = r_face[i + 1];
        const double V_i = cellVolumeSpherical(r_lo, r_hi);

        // East coefficient (link to cell i+1).
        double a_e = 0.0;
        if (i < N_ - 1) {
            const double A_e = faceAreaSpherical(r_hi);
            const double dr_e = safeMax(1.0e-12, r_cell[i + 1] - r_cell[i]);
            a_e = A_e * D_face_[i + 1] / dr_e;
        }

        // West coefficient (link to cell i-1).
        double a_w = 0.0;
        if (i > 0) {
            const double A_w = faceAreaSpherical(r_lo);
            const double dr_w = safeMax(1.0e-12, r_cell[i] - r_cell[i - 1]);
            a_w = A_w * D_face_[i] / dr_w;
        }

        // Source linearisation: V_i c kappa_P rho_i. Sink coefficient
        // contributes to the diagonal as +V_i c kappa_P rho_i dt; the
        // emission term contributes a T_iter^4 to the RHS multiplied
        // by V_i c kappa_P rho_i dt.
        const double c_kp_rho_V = V_i * SPEED_OF_LIGHT_M_PER_S *
                                  kappa_P_[i] * rho[i];
        const double T = T_iter_[i];
        const double T4 = T * T * T * T;

        // Backward-Euler: V_i (E^{n+1} - E^n)/dt = (-flux + source) at
        // n+1. Multiplying through by dt and packing:
        //   V_i E^{n+1} - dt (a_e (E_{i+1} - E_i) - a_w (E_i - E_{i-1}))
        //     - dt c_kp_rho_V (a T^4 - E^{n+1}) = V_i E^n.
        // -> diag = V_i + dt (a_e + a_w) + dt c_kp_rho_V
        //    sub  = -dt a_w
        //    sup  = -dt a_e
        //    rhs  = V_i E^n + dt c_kp_rho_V * a * T_iter^4.
        a_lower_[i] = -dt * a_w;
        c_upper_[i] = -dt * a_e;
        b_diag_[i] = V_i + dt * (a_e + a_w) + dt * c_kp_rho_V;
        rhs_[i] = V_i * E_r_old_[i] +
                  dt * c_kp_rho_V * RADIATION_CONSTANT_A * T4;
    }

    // Inner BC (i = 0): zero-flux symmetry. The west coefficient is
    // already zero because i > 0 is false. No additional change.

    // Outer BC (i = N_-1): Marshak / extrapolation -- E_r at the outer
    // ghost = a T_amb^4, which couples to the rightmost cell through
    // an additional flux term. To keep the matrix tridiagonal we add
    // a Dirichlet contribution to b_diag_ and rhs_ at the last cell.
    if (N_ >= 2) {
        const double r_last = r_face[N_];
        const double A_last = faceAreaSpherical(r_last);
        const double dr_ghost = safeMax(1.0e-3, r_face[N_] - r_cell[N_ - 1]);
        const double D_ghost = D_face_[N_];
        const double a_ghost = A_last * D_ghost / dr_ghost;
        const double T_amb4 = config_.T_ambient_K *
                              config_.T_ambient_K *
                              config_.T_ambient_K *
                              config_.T_ambient_K;
        const double E_r_outer = RADIATION_CONSTANT_A * T_amb4;
        b_diag_[N_ - 1] += dt * a_ghost;
        rhs_[N_ - 1] += dt * a_ghost * E_r_outer;
    }
}

void MarshakRadiationDiffusionSolver::solveTridiagonal()
{
    if (N_ == 0) return;
    // Standard Thomas algorithm: forward sweep, back substitution.
    // Rewrites b_diag_ and rhs_ in place.
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

void MarshakRadiationDiffusionSolver::updateMatterTemperature(
    double dt, const std::vector<double>& rho, const TillotsonEOS& eos)
{
    (void)rho;
    // Pass-8: cv comes from the Tillotson parameter set (per-medium
    // constant). A cv(rho, e) lookup is named as a pass-9 candidate.
    // Per-cell linearised matter-energy update. From
    //   rho cv (T^{n+1} - T^n) / dt = c kappa_P rho ( E_r^{n+1} - a T^{n+1}^4 )
    // linearise T^{n+1}^4 ~ T_iter^4 + 4 T_iter^3 (T^{n+1} - T_iter):
    //   rho cv (T^{n+1} - T^n) / dt = c kappa_P rho * E_r^{n+1}
    //                                - c kappa_P rho a (T_iter^4 +
    //                                  4 T_iter^3 (T^{n+1} - T_iter))
    // Solve for T^{n+1}:
    //   T^{n+1} (1/(dt) + (c kappa_P / cv) a 4 T_iter^3)
    //     = T^n / dt + (c kappa_P / cv) (E_r^{n+1} / rho_factor)
    //                - (c kappa_P / cv) a (T_iter^4 - 4 T_iter^4)
    //
    // For pass-8 we use cv from the configured Tillotson parameter
    // set; cv is (Tillotson) cv_J_per_kg_K. Densities differ per cell;
    // the (rho cv) factor cancels with the rho on the RHS in the
    // emission/absorption term.
    const double cv = safeMax(1.0e-3, eos.getParameters().cv);
    for (int i = 0; i < N_; ++i) {
        const double T_n = T_m_old_[i];
        const double T_iter = T_iter_[i];
        const double T3 = T_iter * T_iter * T_iter;
        const double T4 = T3 * T_iter;
        const double E_new = x_[i];
        const double a = RADIATION_CONSTANT_A;
        const double c = SPEED_OF_LIGHT_M_PER_S;
        const double kp = kappa_P_[i];

        // Per-cell scalar update. Note rho cancels since both sides
        // scale linearly in rho when dividing through.
        const double alpha = c * kp * a * 4.0 * T3 / cv;
        const double rhs_cell = T_n / dt + (c * kp / cv) * E_new
                              - (c * kp / cv) * a * (T4 - 4.0 * T4);
        const double denom = 1.0 / dt + alpha;
        const double T_new = rhs_cell / safeMax(1.0e-30, denom);
        T_iter_[i] = T_new > 1.0 ? T_new : 1.0;
    }
}

MarshakRadiationDiffusionSolver::StepResult
MarshakRadiationDiffusionSolver::step(
    double dt, const std::vector<double>& r_cell,
    const std::vector<double>& r_face, const std::vector<double>& rho,
    std::vector<double>& e_int, std::vector<double>& E_r,
    std::vector<double>& T_m, const TillotsonEOS& eos,
    const std::vector<int>& is_gas)
{
    StepResult res;
    if (N_ == 0) {
        res.converged = false;
        return res;
    }

    // Snapshot the old state for the backward-Euler RHS and matter
    // energy bookkeeping. We do not need the per-cell e_int_old beyond
    // the cell-wise difference; record it for the StepResult tally.
    for (int i = 0; i < N_; ++i) {
        E_r_old_[i] = E_r[i];
        T_m_old_[i] = T_m[i] > 1.0 ? T_m[i] : config_.T_ambient_K;
        e_int_old_[i] = e_int[i];
        T_iter_[i] = T_m_old_[i];
        E_iter_[i] = E_r_old_[i];
    }

    int iter = 0;
    double residual = 1.0e30;
    for (; iter < config_.max_newton_iter; ++iter) {
        // Step 1: update opacities at the current Newton iterate.
        updateOpacities(rho, T_iter_);

        // Step 2: update face diffusion coefficient.
        updateFaceDiffusionCoefficient(rho);

        // Step 3: assemble + solve the tridiagonal for E_r^{n+1}.
        assembleTridiagonal(dt, r_cell, r_face, rho);
        solveTridiagonal();

        // Step 4: update matter temperature from the new E_r.
        updateMatterTemperature(dt, rho, eos);

        // Step 5: compute residual on E_r.
        residual = 0.0;
        for (int i = 0; i < N_; ++i) {
            const double dE = std::abs(x_[i] - E_iter_[i]);
            const double scale = safeMax(1.0e-12, std::abs(x_[i]) +
                                                    std::abs(E_iter_[i]));
            const double rel = dE / scale;
            if (rel > residual) residual = rel;
            E_iter_[i] = x_[i];
        }
        if (residual < config_.newton_tolerance) {
            ++iter;
            break;
        }
    }

    // Apply the converged state to the host's vectors. The per-cell
    // matter-energy increment de = cv (T_new - T_old) [J/kg]. We do
    // not modify gas cells aggressively: the cavity gas already starts
    // far above the radiation-energy level so the radiation source
    // there is negligible; preserving e_int there keeps the host's
    // EOS update simple.
    const double cv = safeMax(1.0e-3, eos.getParameters().cv);
    double matter_dE_total = 0.0;
    for (int i = 0; i < N_; ++i) {
        E_r[i] = E_iter_[i];
        const double T_new = T_iter_[i];
        T_m[i] = T_new;
        const double de = cv * (T_new - T_m_old_[i]);
        e_int[i] += de;
        const double V_i = cellVolumeSpherical(r_face[i], r_face[i + 1]);
        matter_dE_total += rho[i] * de * V_i;
        (void)is_gas;
    }

    // Hand-off diagnostic: outermost cell with E_r above the ambient
    // floor. t_diff at that cell is (dr)^2 / D_face.
    const double T_amb4 = config_.T_ambient_K * config_.T_ambient_K *
                          config_.T_ambient_K * config_.T_ambient_K;
    const double E_floor = config_.front_factor *
                           RADIATION_CONSTANT_A * T_amb4;
    int front_idx = 0;
    for (int i = N_ - 1; i >= 0; --i) {
        if (E_r[i] > E_floor) { front_idx = i; break; }
    }
    double t_diff = 0.0;
    if (front_idx >= 0 && front_idx < N_) {
        const double dr = safeMax(1.0e-9, r_face[front_idx + 1] -
                                            r_face[front_idx]);
        const double D = safeMax(1.0e-30, D_face_[front_idx + 1]);
        t_diff = dr * dr / D;
    }

    // Total E_r * V over all cells.
    double E_r_total = 0.0;
    for (int i = 0; i < N_; ++i) {
        E_r_total += E_r[i] *
                     cellVolumeSpherical(r_face[i], r_face[i + 1]);
    }

    res.converged = (residual < config_.newton_tolerance);
    res.newton_iters = iter;
    res.residual_inf_norm = residual;
    res.radiation_front_index = front_idx;
    res.radiation_front_radius_m = (front_idx < N_) ? r_cell[front_idx] : 0.0;
    res.t_diff_at_front_s = t_diff;
    res.total_radiation_energy_J = E_r_total;
    res.total_matter_energy_change_J = matter_dE_total;
    return res;
}

} // namespace FSRM
