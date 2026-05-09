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
    if (config_.opacity_model == OpacityModel::TABULATED_TOPS) {
        // Pass-8 retained TABULATED_TOPS as a throwing stub; pass-9
        // promotes TABULATED_PATCHED to a working option but TOPS
        // remains a named-only legacy alias for any external code that
        // still references it. Pass-10 keeps the alias.
        throw std::runtime_error(
            "MarshakRadiationDiffusionSolver: opacity_model=TABULATED_TOPS "
            "is the pass-8 legacy throwing stub; use TABULATED_PATCHED "
            "(grey) or TABULATED_FULL (pass-10) instead.");
    }

    // Pass-9/10: lazy-load tabulated opacity tables when configured.
    // TABULATED_PATCHED uses a sin^2 blend with the Z-R baseline; pass-10
    // TABULATED_FULL skips the blend and uses the table everywhere with
    // a Z-R safety net only for out-of-table queries.
    if (!opacity_table_load_attempted_ &&
        (config_.opacity_model == OpacityModel::TABULATED_PATCHED ||
         config_.opacity_model == OpacityModel::TABULATED_FULL)) {
        opacity_table_load_attempted_ = true;
        const bool both_paths =
            !config_.tabulated_opacity_rosseland_path.empty() &&
            !config_.tabulated_opacity_planck_path.empty();
        const char* mode_name =
            config_.opacity_model == OpacityModel::TABULATED_FULL
                ? "TABULATED_FULL"
                : "TABULATED_PATCHED";
        if (!both_paths) {
            std::fprintf(stderr,
                "MarshakRadiationDiffusionSolver: %s requires both "
                "rosseland and planck table paths; one or both are empty. "
                "Falling back to POWER_LAW_ZR.\n", mode_name);
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
    E_r_prev_step_.assign(N_, 0.0);
    prev_step_valid_ = false;
    power_law_.setParameters(config_.opacity_params);
    if (!time_integrator_) {
        time_integrator_ =
            DiffusionTimeIntegrator::create(config_.time_integrator);
    }
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
    case OpacityModel::TABULATED_FULL: {
        // Pass-10 grey TABULATED_FULL: pure tabulated everywhere with
        // a POWER_LAW_ZR safety net for out-of-table queries. The
        // pass-9 sin^2 blend window is removed entirely so users get
        // exactly what the table says (one-time stderr warning logged
        // by the reader on first OOR hit).
        if (!rosseland_table_.isLoaded() || !planck_table_.isLoaded()) {
            PowerLawOpacity ev(config_.opacity_params);
            kappa_R = ev.rosseland(rho, T);
            kappa_P = ev.planck(rho, T);
            return;
        }
        const double k_R_tab = rosseland_table_.evaluate(rho, T);
        const double k_P_tab = planck_table_.evaluate(rho, T);
        if (std::isnan(k_R_tab) || std::isnan(k_P_tab)) {
            PowerLawOpacity ev(config_.opacity_params);
            kappa_R = ev.rosseland(rho, T);
            kappa_P = ev.planck(rho, T);
            return;
        }
        kappa_R = k_R_tab;
        kappa_P = k_P_tab;
        return;
    }
    case OpacityModel::TABULATED_TOPS:
    default:
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
    // Pass-11 generalisation of the per-cell tridiagonal assembly.
    // The DiffusionTimeIntegrator strategy supplies the seven scalars
    // that parameterise the BE / CN / BDF2 forms; this routine combines
    // them with the spatial-discretisation coefficients (V_i, a_e, a_w,
    // c_kp_rho_V) and the linearised matter source S_emit(T_iter) to
    // build the linear system.
    //
    // For BACKWARD_EULER (the pass-10 default) the coefficients reduce
    // to lhs_volume = 1, lhs_implicit = 1, rhs_volume_n = 1, all other
    // coefficients zero except rhs_source_implicit = 1 -- which is the
    // pass-10 assembly byte-for-byte.
    const auto coef = time_integrator_
                          ? time_integrator_->getCoefficients()
                          : BackwardEulerIntegrator().getCoefficients();
    const bool needs_explicit_diffusion =
        std::abs(coef.rhs_explicit_diffusion_factor) > 0.0;
    const bool needs_explicit_source =
        std::abs(coef.rhs_source_explicit_factor) > 0.0;
    const bool needs_prev_state =
        std::abs(coef.rhs_volume_nm1_factor) > 0.0;
    // BDF2 bootstrap: when E^{n-1} is requested but unavailable (first
    // step after initialisation), fall back to backward-Euler for this
    // step. The next step has prev_step_valid_ = true so BDF2 engages.
    const bool bootstrap_to_be = needs_prev_state && !prev_step_valid_;
    const auto eff = bootstrap_to_be
                         ? BackwardEulerIntegrator().getCoefficients()
                         : coef;

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

        const double c_kp_rho_V = V_i * SPEED_OF_LIGHT_M_PER_S *
                                  kappa_P_[i] * rho[i];
        const double T = T_iter_[i];
        const double T4 = T * T * T * T;
        const double S_emit_implicit = RADIATION_CONSTANT_A * T4;
        // Explicit (n) source uses T^n (T_m_old_), evaluated outside
        // the Newton loop because that snapshot is fixed for the step.
        double S_emit_explicit = 0.0;
        if (needs_explicit_source) {
            const double T_old = T_m_old_[i];
            const double T_old4 = T_old * T_old * T_old * T_old;
            S_emit_explicit = RADIATION_CONSTANT_A * T_old4;
        }

        a_lower_[i] = -dt * a_w * eff.lhs_implicit_factor;
        c_upper_[i] = -dt * a_e * eff.lhs_implicit_factor;
        b_diag_[i] = V_i * eff.lhs_volume_factor +
                     dt * (a_e + a_w + c_kp_rho_V) * eff.lhs_implicit_factor;

        double rhs_cell = V_i * eff.rhs_volume_n_factor * E_r_old_[i];
        if (needs_prev_state && !bootstrap_to_be) {
            rhs_cell += V_i * eff.rhs_volume_nm1_factor * E_r_prev_step_[i];
        }
        if (needs_explicit_diffusion) {
            const double E_left  = (i > 0)        ? E_r_old_[i - 1] : E_r_old_[i];
            const double E_right = (i < N_ - 1)   ? E_r_old_[i + 1] : E_r_old_[i];
            const double explicit_diff =
                a_e * (E_right - E_r_old_[i]) -
                a_w * (E_r_old_[i] - E_left) -
                c_kp_rho_V * E_r_old_[i];
            rhs_cell += dt * eff.rhs_explicit_diffusion_factor * explicit_diff;
        }
        rhs_cell += dt * c_kp_rho_V *
                    (eff.rhs_source_implicit_factor * S_emit_implicit +
                     eff.rhs_source_explicit_factor * S_emit_explicit);
        rhs_[i] = rhs_cell;
    }

    // Inner BC (i = 0): zero-flux symmetry. The west coefficient is
    // already zero because i > 0 is false. No additional change.

    // Outer BC (i = N_-1): Marshak / extrapolation -- E_r at the outer
    // ghost = a T_amb^4, which couples to the rightmost cell through
    // an additional flux term. To keep the matrix tridiagonal we add
    // a Dirichlet contribution to b_diag_ and rhs_ at the last cell.
    // Apply the integrator's lhs_implicit_factor on the diag; the RHS
    // ghost carries the matching explicit weight when the integrator
    // splits the diffusion operator (CN). BE/BDF2 leave it implicit.
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
        b_diag_[N_ - 1] += dt * a_ghost * eff.lhs_implicit_factor;
        rhs_[N_ - 1] += dt * a_ghost * E_r_outer *
                        (eff.lhs_implicit_factor +
                         eff.rhs_explicit_diffusion_factor);
        if (std::abs(eff.rhs_explicit_diffusion_factor) > 0.0) {
            rhs_[N_ - 1] -= dt * a_ghost *
                            eff.rhs_explicit_diffusion_factor *
                            E_r_old_[N_ - 1];
        }
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

    // Snapshot the old state for the (BE/CN/BDF2) RHS. The
    // explicit-source weighting under CN reads T_m_old_, fixed for the
    // step. The BDF2 RHS reads E^{n-1} from E_r_prev_step_; we shift
    // the prior state at the end of a successful step.
    // Stash E^n in a local buffer so we can promote it into E_r_prev_
    // after the Newton iteration converges.
    std::vector<double> E_r_step_input(N_, 0.0);
    for (int i = 0; i < N_; ++i) {
        E_r_old_[i] = E_r[i];
        E_r_step_input[i] = E_r[i];
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

    // Pass-11: shift prior-step state for the next BDF2 advance. After
    // a successful step at time n+1, the value that was E^n becomes
    // E^{n-1} for the next call. We always update the buffer (cheap)
    // but only mark it valid when the integrator that just ran needs
    // two prior states (BDF2). Other integrators ignore it.
    for (int i = 0; i < N_; ++i) {
        E_r_prev_step_[i] = E_r_step_input[i];
    }
    if (time_integrator_ && time_integrator_->needsTwoPriorStates()) {
        prev_step_valid_ = true;
    }
    return res;
}

} // namespace FSRM
