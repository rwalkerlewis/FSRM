/**
 * @file RadialLagrangian.cpp
 * @brief Implementation of the 1D radial Lagrangian elastoplastic shock
 *        solver. See include/domain/explosion/RadialLagrangian.hpp for
 *        the design rationale and references.
 *
 * Per-step update sequence (Wilkins 1980, ch. 3):
 *   1. CFL-limited dt from cell-wise (c_p + |v|).
 *   2. Lagrangian face advection.
 *   3. Mass-conservation density update.
 *   4. Wilkins artificial viscosity (linear + quadratic on negative
 *      volumetric strain rate).
 *   5. Momentum update at faces from -d sigma_rr / dr - 2 (sigma_rr -
 *      sigma_tt) / r geometric source. Outer face uses an outgoing
 *      characteristic BC.
 *   6. Trial elastic predictor for the deviatoric radial stress.
 *   7. Drucker-Prager radial return.
 *   8. Mie-Gruneisen EOS for solid cells. Ideal-gas EOS for the inner
 *      cavity-gas cell.
 *   9. Internal-energy update including plastic dissipation.
 *  10. Damage evolution from DamageEvolutionModel.
 *  11. Surface-integral moment-rate extraction at the fixed elastic
 *      radius.
 */

#include "domain/explosion/RadialLagrangian.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <tuple>

#include "domain/explosion/MarshakRadiationDiffusion.hpp"

namespace FSRM {

namespace
{

constexpr double TWO_PI = 6.28318530717958647692;
constexpr double FOUR_PI = 12.5663706143591729539;

double safeMax(double a, double b) { return a > b ? a : b; }
double safeMin(double a, double b) { return a < b ? a : b; }

} // namespace

double RadialLagrangianSolver::cellVolumeSpherical(double r_lo, double r_hi)
{
    return (FOUR_PI / 3.0) * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
}

double RadialLagrangianSolver::faceAreaSpherical(double r)
{
    return FOUR_PI * r * r;
}

RadialLagrangianSolver::RadialLagrangianSolver() = default;

void RadialLagrangianSolver::setSource(const UndergroundExplosionSource& src)
{
    src_ = src;
    src_.computeRiseTime();
}

void RadialLagrangianSolver::setEOS(const MieGruneisenEOS& eos) { eos_ = eos; }
void RadialLagrangianSolver::setStrength(const PressureDependentStrength& s) { strength_ = s; }
void RadialLagrangianSolver::setDamage(const DamageEvolutionModel& d) { damage_model_ = d; }
void RadialLagrangianSolver::setConfig(const Config& c)
{
    config_ = c;
    // Pass-7: cache the Tillotson evaluator from the configured
    // parameter set so the per-step inner-cell pressure update does
    // not allocate a new EOS each call.
    cavity_tillotson_.setParameters(c.tillotson_params);

    // Pass-8/10: lazily construct the radiation solver under MARSHAK_GREY
    // (pass-8) or MARSHAK_MULTIGROUP (pass-10). SN_TRANSPORT remains
    // a named-only ladder rung and throws.
    if (c.radiation_phase == RadiationPhase::SN_TRANSPORT) {
        throw std::runtime_error(
            "RadialLagrangianSolver: radiation_phase=SN_TRANSPORT is "
            "named only; pass-10 implements MARSHAK_MULTIGROUP as the "
            "HIGHEST tier on the radiation ladder.");
    }
    // Pass-10 dispatch matrix: TABULATED_PATCHED + MULTIGROUP is not
    // a valid pairing. The per-group analytic opacity model is the
    // documented opacity path under MULTIGROUP; the 2D patched table
    // is the documented opacity path under GREY only.
    if (c.radiation_phase == RadiationPhase::MARSHAK_MULTIGROUP &&
        c.opacity_model == OpacityModel::TABULATED_PATCHED) {
        throw std::runtime_error(
            "RadialLagrangianSolver: opacity_model=TABULATED_PATCHED is "
            "incompatible with radiation_phase=MARSHAK_MULTIGROUP; the "
            "multigroup path uses the per-group analytic opacity model "
            "(Mihalas-Mihalas 1984 sec 82.2). Use opacity_model="
            "TABULATED_FULL to opt into the multigroup analytic path "
            "explicitly, or POWER_LAW_ZR for the legacy baseline.");
    }
    if (c.radiation_phase == RadiationPhase::MARSHAK_MULTIGROUP) {
        if (!mg_rad_solver_) {
            mg_rad_solver_.reset(new MultigroupRadiationDiffusionSolver());
        }
        MultigroupRadiationDiffusionSolver::Config mcfg;
        mcfg.group_grid = c.multigroup_grid;
        mcfg.max_newton_iter = c.radiation_max_newton_iter;
        mcfg.newton_tolerance = c.radiation_newton_tolerance;
        mcfg.opacity_params = c.opacity_params;
        mcfg.T_ambient_K = 300.0;
        mcfg.front_factor = 1.5;
        mg_rad_solver_->setConfig(mcfg);
        mg_num_groups_ = mcfg.group_grid.n_groups;
    }
    if (c.radiation_phase == RadiationPhase::MARSHAK_GREY) {
        if (!rad_solver_) {
            rad_solver_.reset(new MarshakRadiationDiffusionSolver());
        }
        MarshakRadiationDiffusionSolver::Config rcfg;
        rcfg.max_newton_iter = c.radiation_max_newton_iter;
        rcfg.newton_tolerance = c.radiation_newton_tolerance;
        rcfg.opacity_model = c.opacity_model;
        rcfg.opacity_params = c.opacity_params;
        rcfg.kappa_constant_m2_per_kg = c.kappa_constant_m2_per_kg;
        // Pass-9: forward the tabulated opacity table paths and the
        // sin^2 blend window. The Marshak solver lazy-loads on its
        // initialize().
        rcfg.tabulated_opacity_rosseland_path =
            c.tabulated_opacity_rosseland_path;
        rcfg.tabulated_opacity_planck_path =
            c.tabulated_opacity_planck_path;
        rcfg.tabulated_blend_lower_k = c.tabulated_opacity_blend_lower_k;
        rcfg.tabulated_blend_upper_k = c.tabulated_opacity_blend_upper_k;
        rad_solver_->setConfig(rcfg);
    }

    // Pass-9: invalidate the cavity EOS table on every setConfig so a
    // re-config with a new path picks it up on the next cavityPressure
    // call. The table itself is reloaded lazily.
    cavity_eos_table_load_attempted_ = false;
    cavity_eos_table_load_succeeded_ = false;
}

double RadialLagrangianSolver::cavityPressure(double rho, double e) const
{
    // Tillotson is the baseline analytic; TILLOTSON, TILLOTSON_TABULATED_PATCH,
    // and TABULATED_FULL all evaluate it as a fallback.
    const double p_til_raw = cavity_tillotson_.pressure(rho, e);
    const double p_til = p_til_raw > 0.0 ? p_til_raw : 0.0;

    if (config_.cavity_eos == CavityEOS::TILLOTSON) {
        return p_til;
    }

    if (config_.cavity_eos == CavityEOS::IDEAL_GAS) {
        const double rho_safe = (rho > 1.0e-6) ? rho : 1.0e-6;
        const double e_safe = (e > 0.0) ? e : 0.0;
        const double p_gas = (config_.gas_eos_gamma - 1.0) * rho_safe * e_safe;
        return p_gas > 0.0 ? p_gas : 0.0;
    }

    // TILLOTSON_TABULATED_PATCH and TABULATED_FULL share the
    // table-loading path; the dispatch below distinguishes them.
    if (!cavity_eos_table_load_attempted_) {
        cavity_eos_table_load_attempted_ = true;
        if (!config_.tabulated_eos_table_path.empty()) {
            std::string err;
            cavity_eos_table_load_succeeded_ =
                cavity_eos_table_.load(config_.tabulated_eos_table_path, err);
            if (!cavity_eos_table_load_succeeded_) {
                std::fprintf(stderr,
                    "RadialLagrangianSolver: cavity_eos table load FAILED: "
                    "%s. Falling back to Tillotson.\n", err.c_str());
            }
        } else {
            std::fprintf(stderr,
                "RadialLagrangianSolver: %s selected but "
                "tabulated_eos_table_path is empty. Falling back to "
                "Tillotson.\n",
                config_.cavity_eos == CavityEOS::TABULATED_FULL
                    ? "TABULATED_FULL"
                    : "TILLOTSON_TABULATED_PATCH");
        }
    }

    if (!cavity_eos_table_load_succeeded_) {
        return p_til;
    }

    // Pass-10 TABULATED_FULL: no patch. Direct table lookup with a
    // Tillotson safety net for out-of-table queries (one-time stderr
    // warning; the table reader logs the warning on first OOR hit).
    if (config_.cavity_eos == CavityEOS::TABULATED_FULL) {
        const double p_tab_raw = cavity_eos_table_.evaluate(rho, e);
        if (std::isnan(p_tab_raw)) {
            return p_til;
        }
        return p_tab_raw > 0.0 ? p_tab_raw : 0.0;
    }

    // TILLOTSON_TABULATED_PATCH: pass-9 sin^2 blend on Tillotson pressure.
    const double p_lo = config_.tabulated_eos_blend_lower_pa;
    const double p_hi = config_.tabulated_eos_blend_upper_pa;
    if (p_til <= p_lo) {
        return p_til;
    }
    const double p_tab_raw = cavity_eos_table_.evaluate(rho, e);
    if (std::isnan(p_tab_raw)) {
        return p_til;
    }
    const double p_tab = p_tab_raw > 0.0 ? p_tab_raw : 0.0;
    if (p_til >= p_hi) {
        return p_tab;
    }
    const double t = (p_til - p_lo) / (p_hi - p_lo);
    const double s = std::sin(0.5 * M_PI * t);
    const double w = s * s;
    return (1.0 - w) * p_til + w * p_tab;
}

void RadialLagrangianSolver::solveCavityInitialState()
{
    // Pass-7 first-principles energy-partition solve for the inner
    // cavity at t = t_rh.
    //
    // Premise (Zel'dovich-Raizer 1967, vol II, ch X). At the
    // radiation-to-hydrodynamic transition time, the radiation wave
    // has heated the host rock in place faster than the hydrodynamic
    // expansion can move material. The cavity contains fully-vaporized
    // host rock at approximately the solid density (rho_v ~ rho_0_solid)
    // because no significant hydrodynamic motion has yet occurred. The
    // pass-7 implementation uses this end-state-of-radiation-phase
    // approximation; an explicit Marshak-wave radiation-diffusion solve
    // is named as candidate pass-8 follow-up if the gap does not close.
    //
    // Energy partition. The total deposited yield must equal:
    //   - vaporization latent heat: m_v * E_cv (definitional energy
    //     to bring mass m_v from solid to fully vaporized);
    //   - thermal internal energy of the vapor: m_v * (e_v - E_cv);
    //   - gravitational potential energy of the displaced overburden:
    //     m_v * g * h_eff (small for kt yields; included for
    //     completeness);
    //   - residual kinetic energy: zero by definition at t = t_rh
    //     (the radiation phase deposits energy in place, no bulk
    //     motion yet).
    //
    // Constraint. m_v = (4/3) pi R_v^3 rho_v. With rho_v fixed at
    // rho_0_solid (the Z-R approximation) this gives a single
    // unknown R_v. We solve a one-equation Newton iteration:
    //
    //   f(R_v) = E_yield
    //          - m_v(R_v) * E_cv
    //          - m_v(R_v) * (e_v_target - E_cv)
    //          - m_v(R_v) * g * h_eff
    //
    // where e_v_target is set by the requirement that the Tillotson
    // pressure at (rho_v, e_v_target) is large but finite (ensuring
    // the cavity gas drives the surrounding rock); for the pass-7
    // initialization we set e_v_target = a multiple of E_cv chosen so
    // the integrated yield is consumed entirely by latent + thermal +
    // potential. Newton converges in 5-10 iterations.
    //
    // Default radiation-transition time. For yields in the kt-Mt
    // range, t_rh ~ 1e-7 s * W_kt^(1/3) (Z-R scaling, vol II
    // eq. 24.18). The user can override via radiation_transition_time_s.
    //
    // References.
    //  - Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
    //    Waves and High-Temperature Hydrodynamic Phenomena", vol II,
    //    Academic Press, ch X (Marshak waves and the
    //    radiation-to-hydrodynamic transition).
    //  - Melosh, H. J. (1989), "Impact Cratering: A Geologic Process",
    //    eqs 5.4.7-9 and Table A2.2 (Tillotson EOS evaluation in the
    //    hot-expanded regime, used for e_v_target -> p consistency).

    const double E_yield = src_.energyJoules();
    const double rho_0_solid = src_.host_density;
    const double g = 9.81;
    const double h_eff = src_.depth;

    // Z-R-scaled default transition time when not user-specified.
    if (config_.radiation_transition_time_s > 0.0) {
        t_rh_used_ = config_.radiation_transition_time_s;
    } else {
        t_rh_used_ = 1.0e-7 * std::pow(safeMax(1.0e-6, src_.yield_kt),
                                       1.0 / 3.0);
    }

    // Tillotson parameters drive the cavity-state target.
    const TillotsonParameters& tp = config_.tillotson_params;
    rho_v_init_ = rho_0_solid;
    // Target specific internal energy: full vaporization (E_cv).
    // Keeping e_v = E_cv puts the cavity in the just-vaporized
    // compressed-branch state; for granite this evaluates the
    // Tillotson form to ~50 GPa cavity pressure, plenty to drive
    // the surrounding shock. The energy partition then becomes
    //   E_yield = m_v * E_cv + m_v * g * h_eff
    // (latent-heat-only plus the small overburden potential).
    const double e_v_target = tp.E_cv;
    e_v_init_ = e_v_target;
    const double e_thermal = e_v_target - tp.E_cv;  // 0 by construction

    // Closed-form initial guess from the latent + thermal balance
    // (ignoring the small overburden potential):
    //   E_yield = (4/3) pi R^3 rho_0 (E_cv + e_thermal)
    // -> R^3 = 3 E / (4 pi rho_0 (E_cv + e_thermal))
    const double e_total_per_kg = tp.E_cv + e_thermal;
    double R_v = std::pow(
        3.0 * E_yield /
            (FOUR_PI * safeMax(1.0, rho_0_solid * e_total_per_kg)),
        1.0 / 3.0);

    // Newton iteration on the one-equation residual.
    auto mass_v = [&](double R) {
        return (FOUR_PI / 3.0) * R * R * R * rho_v_init_;
    };
    auto residual = [&](double R) {
        const double m = mass_v(R);
        return E_yield - m * (tp.E_cv + e_thermal) - m * g * h_eff;
    };
    auto residual_derivative = [&](double R) {
        // dR^3 / dR = 3 R^2 -> dm/dR = 4 pi R^2 rho_v
        const double dm_dR = FOUR_PI * R * R * rho_v_init_;
        return -dm_dR * (tp.E_cv + e_thermal) - dm_dR * g * h_eff;
    };

    int iters = 0;
    for (; iters < 30; ++iters) {
        const double f = residual(R_v);
        const double df = residual_derivative(R_v);
        if (std::abs(df) < 1e-30) break;
        const double dR = -f / df;
        R_v += dR;
        if (R_v < 0.0) R_v = 0.5 * (R_v - dR);
        if (std::abs(dR) < 1.0e-6 * R_v) break;
    }
    if (R_v <= 0.0 || !std::isfinite(R_v)) {
        // Fallback: closed-form initial guess.
        R_v = std::pow(
            3.0 * E_yield /
                (FOUR_PI * safeMax(1.0, rho_0_solid * e_total_per_kg)),
            1.0 / 3.0);
    }

    // Clamp R_v to a sensible upper band around the empirical NTS
    // cavity radius. The Z-R end-state cavity is the *vaporization*
    // radius, which is typically much smaller than the eventual NTS
    // cavity radius (the bulk cavity expansion happens during the
    // hydrodynamic phase that follows). A cap at 1.5 times the NTS
    // Rc prevents the Newton iteration from running away in
    // pathological parameter regimes. We do NOT clamp from below
    // beyond a small absolute floor because the analytic latent-heat
    // estimate for kt yields can give R_v values well below 0.05
    // times the eventual NTS cavity radius.
    const double Rc_eq = safeMax(1.0, src_.cavityRadius());
    const double R_v_max = 1.5 * Rc_eq;
    if (R_v < 1.0e-3) R_v = 1.0e-3;
    if (R_v > R_v_max) R_v = R_v_max;
    Rc_init_ = R_v;
}

void RadialLagrangianSolver::allocate(int N)
{
    N_ = N;
    r_cell_.assign(N, 0.0);
    mass_.assign(N, 0.0);
    rho_.assign(N, 0.0);
    e_int_.assign(N, 0.0);
    p_.assign(N, 0.0);
    s_rr_.assign(N, 0.0);
    q_visc_.assign(N, 0.0);
    eps_p_.assign(N, 0.0);
    damage_.assign(N, 0.0);
    yielded_.assign(N, 0);
    is_gas_.assign(N, 0);

    r_face_.assign(N + 1, 0.0);
    v_face_.assign(N + 1, 0.0);
}

void RadialLagrangianSolver::initialize()
{
    const int N = safeMax(20, config_.radial_cells);

    // The host-medium NTS cavity radius sets the elastic radius at
    // which the moment-tensor surface integral is evaluated. The
    // configurable radial_outer_factor multiplies the elastic radius
    // to give the outer-domain extent so the outgoing wave has room
    // to leave the source region before the characteristic BC fires.
    //
    // Pass-7: the initial cavity radius Rc_init_ is set by either
    // the physics-based energy-partition Newton solve (default) or
    // by the user-supplied initial_cavity_radius_m (MANUAL). Pass-6
    // unconditionally used Rc_init_ = 0.4 * Rc_eq; that path is
    // preserved as a regression target by selecting
    // cavity_initialization = MANUAL with initial_cavity_radius_m
    // = 0.4 * Rc_eq in the user config.
    const double Rc_eq = safeMax(1.0, src_.cavityRadius());
    if (config_.cavity_initialization ==
            CavityInitialization::PHYSICS_BASED) {
        solveCavityInitialState();
    } else {
        Rc_init_ = (config_.initial_cavity_radius_m > 0.0)
                       ? config_.initial_cavity_radius_m
                       : 0.4 * Rc_eq;
        rho_v_init_ = src_.host_density;
        e_v_init_ = 0.0;
        t_rh_used_ = 0.0;
    }
    r_elastic_ = config_.radial_outer_factor * Rc_eq;
    r_outer_ = safeMax(r_elastic_ * 1.20, 1.5 * Rc_eq);
    // Pass-9: explicit override for the free-field peak-velocity gate
    // and any other test that needs to extend the radial domain past
    // a few elastic radii. r_elastic_ remains at the factor-derived
    // value so the moment-tensor extraction sphere stays where it
    // belongs; only the outer-domain extent grows.
    if (config_.radial_outer_radius_m > 0.0 &&
        config_.radial_outer_radius_m > r_outer_) {
        r_outer_ = config_.radial_outer_radius_m;
    }

    allocate(N);

    // Geometric grid: gas / inner cells uniformly spaced from 0 to
    // Rc_init, solid cells geometrically refined toward Rc_init so the
    // shock front near the cavity wall is resolved. A simple two-zone
    // tanh-stretched layout suffices for pass-6.
    gas_cells_ = safeMax(2, N / 50);
    const double r_gas_inner = 0.0;
    const double r_gas_outer = Rc_init_;
    for (int i = 0; i <= gas_cells_; ++i) {
        r_face_[i] = r_gas_inner +
            (r_gas_outer - r_gas_inner) *
                (static_cast<double>(i) / gas_cells_);
    }
    // Solid cells: geometric stretching from r_gas_outer to r_outer.
    const int n_solid = N - gas_cells_;
    const double L_solid = r_outer_ - r_gas_outer;
    const double stretch = 1.04;  // mild geometric ratio, ~4% per cell
    double accum = 0.0;
    for (int j = 0; j < n_solid; ++j) accum += std::pow(stretch, j);
    const double dr0 = L_solid / accum;
    for (int j = 1; j <= n_solid; ++j) {
        const double sub_accum = [&]() {
            double s = 0.0;
            for (int k = 0; k < j; ++k) s += std::pow(stretch, k);
            return s;
        }();
        r_face_[gas_cells_ + j] = r_gas_outer + dr0 * sub_accum;
    }
    // Force the outermost face exactly to r_outer for floating-point cleanliness.
    r_face_[N] = r_outer_;

    // Cell-centred radii.
    for (int i = 0; i < N; ++i) {
        r_cell_[i] = 0.5 * (r_face_[i] + r_face_[i + 1]);
    }

    // Solid-cell elastic moduli derived from the host medium. K and G
    // come from rho and the P/S wave speeds: G = rho * vs^2,
    // K + 4G/3 = rho * vp^2, so K = rho * (vp^2 - 4 vs^2 / 3).
    mu_solid_ = src_.host_density * src_.host_vs * src_.host_vs;
    K_solid_ = src_.host_density *
               (src_.host_vp * src_.host_vp - 4.0 / 3.0 * src_.host_vs * src_.host_vs);
    if (K_solid_ < 0.1 * mu_solid_) K_solid_ = 0.1 * mu_solid_;
    cp_ref_ = src_.host_vp;

    // Initial state. Solid cells at hydrostatic overburden, zero
    // deviatoric stress, density rho0. Internal energy initialized so
    // the EOS reads back p = overburden at rho0.
    //
    // Pass-7 cavity cells. The (rho_v, e_v) state is determined by
    // the cavity_initialization mode:
    //  - PHYSICS_BASED: rho_v = rho_v_init_ (= rho_0_solid by Z-R),
    //    e_v = e_v_init_ (set by solveCavityInitialState so the
    //    integrated energy partition consumes E_yield exactly). The
    //    cavity pressure follows from the configured cavity_eos
    //    evaluator (TILLOTSON by default, IDEAL_GAS for regression).
    //  - MANUAL: pass-6-style uniform yield-distribution path
    //    preserved; rho = rho_0_solid, e = E_yield / (rho_0_solid *
    //    cavity_volume), p = (gamma - 1) rho e. This branch is the
    //    byte-identical pass-6 regression target.
    const double overburden = src_.overburden_stress;
    const bool physics_based =
        (config_.cavity_initialization ==
         CavityInitialization::PHYSICS_BASED);
    const double cavity_volume_pre =
        (FOUR_PI / 3.0) * Rc_init_ * Rc_init_ * Rc_init_;
    const double e_specific_manual =
        src_.energyJoules() /
        safeMax(1e-12, src_.host_density * cavity_volume_pre);
    for (int i = 0; i < N; ++i) {
        const double r_lo = r_face_[i];
        const double r_hi = r_face_[i + 1];
        const double V = cellVolumeSpherical(r_lo, r_hi);
        if (i < gas_cells_) {
            is_gas_[i] = 1;
            if (physics_based) {
                rho_[i] = (rho_v_init_ > 0.0) ? rho_v_init_
                                              : src_.host_density;
                e_int_[i] = e_v_init_;
            } else {
                rho_[i] = src_.host_density;
                e_int_[i] = e_specific_manual;
            }
            mass_[i] = rho_[i] * V;
            p_[i] = cavityPressure(rho_[i], e_int_[i]);
        } else {
            is_gas_[i] = 0;
            rho_[i] = src_.host_density;
            mass_[i] = rho_[i] * V;
            e_int_[i] = 0.0;
            p_[i] = overburden;
        }
        s_rr_[i] = 0.0;
        eps_p_[i] = 0.0;
        damage_[i] = 0.0;
        yielded_[i] = 0;
    }

    std::fill(v_face_.begin(), v_face_.end(), 0.0);

    // Energy bookkeeping: the recoverable starting deposit. The
    // internal-energy total is the gas-cell deposit; solid cells start
    // at zero specific internal energy. The Wilkins-AV losses appear
    // as additional internal energy (see updateInternalEnergy) so the
    // sum over time moves from "all gas internal" toward "kinetic +
    // solid internal + plastic + radiated".
    initial_energy_ = 0.0;
    for (int i = 0; i < N; ++i) initial_energy_ += mass_[i] * e_int_[i];
    kinetic_energy_ = 0.0;
    internal_energy_ = initial_energy_;
    plastic_dissipation_ = 0.0;
    radiated_energy_out_ = 0.0;

    sigma_rr_extract_prev_ = -overburden;       // sigma = -p in solid
    sigma_rr_extract_initial_ = -overburden;
    extract_initialized_ = false;
    M_iso_.fill(0.0);
    Mdot_iso_.fill(0.0);

    // Pass-8: seed the radiation field if MARSHAK_GREY is selected.
    // We seed E_r at the matter equilibrium a*T^4 in every cell,
    // with T_m at the cavity vapor temperature (from the Tillotson
    // temperature lookup at (rho_v_init_, e_v_init_)) inside the gas
    // region and at T_ambient = 300 K outside. The radiation phase
    // is marked active; the host substep loop will run the Marshak
    // solver until the hand-off criterion ends it.
    radiation_phase_active_ = false;
    handoff_consecutive_steps_ = 0;
    radiation_front_index_ = 0;
    radiation_front_radius_ = 0.0;
    radiation_energy_total_ = 0.0;
    matter_energy_change_from_radiation_ = 0.0;
    t_diff_at_front_ = 0.0;
    tillotson_warning_logged_ = false;
    if (config_.radiation_phase == RadiationPhase::MARSHAK_GREY) {
        if (!rad_solver_) {
            rad_solver_.reset(new MarshakRadiationDiffusionSolver());
            MarshakRadiationDiffusionSolver::Config rcfg;
            rcfg.max_newton_iter = config_.radiation_max_newton_iter;
            rcfg.newton_tolerance = config_.radiation_newton_tolerance;
            rcfg.opacity_model = config_.opacity_model;
            rcfg.opacity_params = config_.opacity_params;
            rcfg.kappa_constant_m2_per_kg = config_.kappa_constant_m2_per_kg;
            rad_solver_->setConfig(rcfg);
        }
        rad_solver_->initialize(N);
        E_r_cell_.assign(N, 0.0);
        T_m_cell_.assign(N, 0.0);
        const double T_amb = 300.0;
        const double a = RadiationConstants::RADIATION_CONSTANT_A_J_PER_M3_K4;
        const double T_cavity =
            cavity_tillotson_.temperature(rho_v_init_ > 0.0 ? rho_v_init_
                                                            : src_.host_density,
                                          e_v_init_);
        const double T_cavity_safe = T_cavity > T_amb ? T_cavity : T_amb;
        for (int i = 0; i < N; ++i) {
            const double T_i = is_gas_[i] ? T_cavity_safe : T_amb;
            T_m_cell_[i] = T_i;
            E_r_cell_[i] = a * T_i * T_i * T_i * T_i;
        }
        radiation_phase_active_ = true;
    }
    if (config_.radiation_phase == RadiationPhase::MARSHAK_MULTIGROUP) {
        if (!mg_rad_solver_) {
            mg_rad_solver_.reset(new MultigroupRadiationDiffusionSolver());
            MultigroupRadiationDiffusionSolver::Config mcfg;
            mcfg.group_grid = config_.multigroup_grid;
            mcfg.max_newton_iter = config_.radiation_max_newton_iter;
            mcfg.newton_tolerance = config_.radiation_newton_tolerance;
            mcfg.opacity_params = config_.opacity_params;
            mg_rad_solver_->setConfig(mcfg);
        }
        mg_rad_solver_->initialize(N);
        mg_num_groups_ = mg_rad_solver_->numGroups();
        T_m_cell_.assign(N, 0.0);
        E_r_g_.assign(static_cast<std::size_t>(N) * mg_num_groups_, 0.0);
        // Seed per-group radiation field at the matter equilibrium for
        // the cell's local T using the MultigroupOpacityEvaluator's
        // band-integrated Planck integrals: E_r^g(0) = 4 pi B_g(T) / c.
        const double T_amb = 300.0;
        const double T_cavity =
            cavity_tillotson_.temperature(rho_v_init_ > 0.0 ? rho_v_init_
                                                            : src_.host_density,
                                          e_v_init_);
        const double T_cavity_safe = T_cavity > T_amb ? T_cavity : T_amb;
        const double FOUR_PI_LOC = 12.5663706143591729539;
        const double C_LIGHT = RadiationConstants::SPEED_OF_LIGHT_M_PER_S;
        const auto& opacity = mg_rad_solver_->opacityEvaluator();
        for (int i = 0; i < N; ++i) {
            const double T_i = is_gas_[i] ? T_cavity_safe : T_amb;
            T_m_cell_[i] = T_i;
            for (int g = 0; g < mg_num_groups_; ++g) {
                const std::size_t k =
                    static_cast<std::size_t>(i) * mg_num_groups_ + g;
                const double Bg = opacity.bandIntegratedPlanck(g, T_i);
                E_r_g_[k] = FOUR_PI_LOC * Bg / C_LIGHT;
            }
        }
        mg_spectral_at_front_.assign(mg_num_groups_, 0.0);
        radiation_phase_active_ = true;
    }

    current_time_ = 0.0;
    initialized_ = true;
}

void RadialLagrangianSolver::cflLimit(double& dt) const
{
    double dt_min = dt;
    for (int i = 0; i < N_; ++i) {
        const double dr = r_face_[i + 1] - r_face_[i];
        if (dr <= 0.0) continue;
        // Effective sound speed: solid uses vp; gas cell uses the
        // configured cavity EOS sound-speed (Tillotson analytic, or
        // sqrt(gamma * p / rho) for the IDEAL_GAS placeholder).
        // Both bounded so a freshly-initialized state with p = 0 in
        // solid cells is not CFL-unbounded.
        double cs;
        if (is_gas_[i]) {
            if (config_.cavity_eos == CavityEOS::TILLOTSON) {
                cs = cavity_tillotson_.soundSpeed(
                    safeMax(1.0e-6, rho_[i]),
                    safeMax(0.0, e_int_[i]));
                cs = safeMax(1.0, cs);
            } else {
                const double p_pos = safeMax(0.0, p_[i]);
                cs = std::sqrt(safeMax(1.0, config_.gas_eos_gamma * p_pos /
                                              safeMax(1e-6, rho_[i])));
            }
        } else {
            cs = safeMax(0.1 * cp_ref_, cp_ref_);
        }
        const double vmag = 0.5 *
            (std::abs(v_face_[i]) + std::abs(v_face_[i + 1]));
        const double dt_cell = config_.cfl * dr / safeMax(1e-3, cs + vmag);
        if (dt_cell < dt_min) dt_min = dt_cell;
    }
    dt = dt_min;
}

void RadialLagrangianSolver::advanceFaces(double dt)
{
    // Inner face stays at the symmetry axis. Other faces move with the
    // local velocity. We do not update v_face_[0]; it stays at 0 to
    // preserve spherical symmetry.
    r_face_[0] = 0.0;
    v_face_[0] = 0.0;
    for (int i = 1; i <= N_; ++i) {
        r_face_[i] += v_face_[i] * dt;
    }
    // Guard against face crossings (would imply CFL violation).
    for (int i = 1; i <= N_; ++i) {
        if (r_face_[i] <= r_face_[i - 1]) {
            r_face_[i] = r_face_[i - 1] + 1e-9;
        }
    }
    for (int i = 0; i < N_; ++i) {
        r_cell_[i] = 0.5 * (r_face_[i] + r_face_[i + 1]);
    }
}

void RadialLagrangianSolver::updateDensity()
{
    for (int i = 0; i < N_; ++i) {
        const double V = cellVolumeSpherical(r_face_[i], r_face_[i + 1]);
        rho_[i] = mass_[i] / safeMax(1e-12, V);
    }
}

void RadialLagrangianSolver::computeArtificialViscosity(double dt)
{
    // Wilkins linear+quadratic AV on negative volumetric strain rate.
    // Active only when dV/dt < 0 (compressing). c_l and c_q are tuning
    // knobs.
    (void)dt;
    for (int i = 0; i < N_; ++i) {
        const double r_lo = r_face_[i];
        const double r_hi = r_face_[i + 1];
        const double v_lo = v_face_[i];
        const double v_hi = v_face_[i + 1];
        const double dr = r_hi - r_lo;
        const double dv = v_hi - v_lo;
        // Volumetric strain rate (eps_dot_kk = div v) in spherical:
        // (1/r^2) d/dr (r^2 v_r). Discretize as
        // (r_hi^2 v_hi - r_lo^2 v_lo) / (r_cell^2 dr).
        const double r_c2 = r_cell_[i] * r_cell_[i];
        const double div_v =
            (r_hi * r_hi * v_hi - r_lo * r_lo * v_lo) /
            safeMax(1e-12, r_c2 * dr);
        if (div_v < 0.0) {
            // Pass-7: use the cavity-EOS sound speed (Tillotson
            // analytic) for gas cells when configured; fall back to
            // sqrt(gamma p / rho) under IDEAL_GAS.
            double cs;
            if (is_gas_[i]) {
                if (config_.cavity_eos == CavityEOS::TILLOTSON) {
                    cs = cavity_tillotson_.soundSpeed(
                        safeMax(1.0e-6, rho_[i]),
                        safeMax(0.0, e_int_[i]));
                    cs = safeMax(1.0, cs);
                } else {
                    cs = std::sqrt(safeMax(
                        1.0, config_.gas_eos_gamma *
                                 safeMax(0.0, p_[i]) /
                                 safeMax(1e-6, rho_[i])));
                }
            } else {
                cs = cp_ref_;
            }
            const double q_lin = config_.art_visc_linear * rho_[i] * cs *
                                 std::abs(div_v) * dr;
            const double q_quad = config_.art_visc_quadratic * rho_[i] *
                                  div_v * div_v * dr * dr;
            q_visc_[i] = q_lin + q_quad;
            (void)dv;
        } else {
            q_visc_[i] = 0.0;
        }
    }
}

void RadialLagrangianSolver::momentumUpdate(double dt)
{
    // Spherical Lagrangian momentum at face i (1 <= i < N):
    //   m_face dv/dt = -A_i * (sigma_rr_+ - sigma_rr_-)
    //                  + 2 * V_blend * (sigma_rr - sigma_tt) / r_face
    // The geometric source (-2 (sigma_rr - sigma_tt) / r) is integrated
    // over a control volume centred on the face. We use a face-centred
    // mass m_face = 0.5 (m_left + m_right) and the average of the
    // adjacent cell stresses for the deviatoric / hoop terms.
    //
    // sigma_rr_total = -p - q_visc + s_rr  (compression negative)
    // sigma_tt_total = -p - q_visc + s_tt = -p - q_visc - s_rr/2
    //
    // The combination (sigma_rr - sigma_tt) = s_rr - s_tt = (3/2) s_rr.
    auto sigma_rr_total = [&](int i) -> double {
        return -p_[i] - q_visc_[i] + s_rr_[i];
    };
    auto sigma_minus_hoop = [&](int i) -> double {
        return 1.5 * s_rr_[i];
    };

    for (int i = 1; i < N_; ++i) {
        const double r = r_face_[i];
        const double A = faceAreaSpherical(r);
        const double left = sigma_rr_total(i - 1);
        const double right = sigma_rr_total(i);
        // Lagrangian: face mass = average of neighbour cell masses.
        const double m_face = 0.5 * (mass_[i - 1] + mass_[i]);
        // Non-conservative form: m dv/dt = A * (sigma_rr_R - sigma_rr_L)
        // + 2 V (sigma_rr - sigma_tt) / r. With sigma > 0 in tension,
        // sigma_rr_R - sigma_rr_L > 0 means the right cell is less
        // compressive than the left, so the face accelerates outward.
        const double pressure_force = A * (right - left);
        const double devhoop_avg = 0.5 * (sigma_minus_hoop(i - 1) +
                                          sigma_minus_hoop(i));
        const double dv_blend = 0.5 *
            (cellVolumeSpherical(r_face_[i - 1], r_face_[i]) +
             cellVolumeSpherical(r_face_[i], r_face_[i + 1]));
        const double geom_source = 2.0 * dv_blend * devhoop_avg / safeMax(1e-12, r);
        v_face_[i] += dt * (pressure_force + geom_source) /
                      safeMax(1e-12, m_face);
    }
    // Outer face handled by absorbingOuterBC after the explicit update
    // and before the next CFL evaluation.
}

void RadialLagrangianSolver::deviatoricElasticPredictor(double dt)
{
    if (config_.disable_plasticity) {
        // Zero deviatoric component for pure-elastic / hydrodynamic
        // sanity tests (Sedov, OutgoingBC).
        for (int i = 0; i < N_; ++i) s_rr_[i] = 0.0;
        return;
    }
    for (int i = 0; i < N_; ++i) {
        if (is_gas_[i]) { s_rr_[i] = 0.0; continue; }
        const double r_lo = r_face_[i];
        const double r_hi = r_face_[i + 1];
        const double v_lo = v_face_[i];
        const double v_hi = v_face_[i + 1];
        const double dr = safeMax(1e-12, r_hi - r_lo);
        const double r_c = safeMax(1e-12, r_cell_[i]);
        const double eps_dot_rr = (v_hi - v_lo) / dr;
        const double eps_dot_tt = 0.5 * (v_lo + v_hi) / r_c;
        // Deviatoric strain rate component: e_rr_dot = (2/3)(eps_rr -
        // eps_tt). With e_phiphi_dot = e_tt_dot we get
        // ds_rr_trial = 2 G e_rr_dot dt = (4G/3)(eps_rr - eps_tt) dt.
        const double diff = eps_dot_rr - eps_dot_tt;
        const double ds_rr = (4.0 / 3.0) * mu_solid_ * diff * dt;
        s_rr_[i] += ds_rr;
        // Damage softening of the deviator.
        s_rr_[i] *= (1.0 - damage_[i]);
    }
}

void RadialLagrangianSolver::radialReturnPlasticity(double dt)
{
    if (config_.disable_plasticity) {
        for (int i = 0; i < N_; ++i) yielded_[i] = 0;
        return;
    }
    for (int i = 0; i < N_; ++i) {
        if (is_gas_[i]) { yielded_[i] = 0; continue; }
        const double sigma_eq = 1.5 * std::abs(s_rr_[i]);  // sqrt(3 J2)
        const double Y = strength_.yieldStrength(p_[i], damage_[i],
                                                 eps_p_[i]);
        if (sigma_eq > Y) {
            const double scale = Y / safeMax(1e-12, sigma_eq);
            const double s_old = s_rr_[i];
            s_rr_[i] = s_old * scale;
            // Equivalent plastic strain increment from
            // dgamma = (sigma_eq_trial - Y) / (3 G).
            const double dgamma = (sigma_eq - Y) /
                                  (3.0 * safeMax(1.0, mu_solid_));
            eps_p_[i] += dgamma;
            // Plastic work per unit volume: Y * dgamma. Energy balance
            // adds m*e per cell so we accumulate as mass * Y * dgamma /
            // rho = Y * dgamma * V.
            const double V = cellVolumeSpherical(r_face_[i], r_face_[i + 1]);
            plastic_dissipation_ += Y * dgamma * V;
            // Plastic dissipation contributes to internal energy
            // (Wilkins eq. 3.19); add to e_int_ at constant volume.
            e_int_[i] += Y * dgamma / safeMax(1e-6, rho_[i]);
            yielded_[i] = 1;
        } else {
            yielded_[i] = 0;
        }
        (void)dt;
    }
}

void RadialLagrangianSolver::updateInternalEnergy(double dt)
{
    // dE/dt per unit mass = - (p + q_visc) (1/rho) drho/dt
    // We compute drho/dt from the new and old volumes by tracking the
    // density in updateDensity() AFTER advancing faces, so the total
    // pressure work done per cell per unit time is approximated as
    // -(p+q) * (V_new - V_old) / dt distributed over the cell.
    // We do not store V_old explicitly here; the alternative is to
    // compute dV from the divergence of the face velocity field, which
    // is what we use: dV/dt = (4 pi r^2 v_r)|_lo^hi.
    for (int i = 0; i < N_; ++i) {
        const double r_lo = r_face_[i];
        const double r_hi = r_face_[i + 1];
        const double v_lo = v_face_[i];
        const double v_hi = v_face_[i + 1];
        const double dV_dt =
            faceAreaSpherical(r_hi) * v_hi - faceAreaSpherical(r_lo) * v_lo;
        const double work = -(p_[i] + q_visc_[i]) * dV_dt;
        // de = (work) * dt / m for the cell.
        const double de = work * dt / safeMax(1e-12, mass_[i]);
        const double e_new = e_int_[i] + de;
        // Pass-7: clip negative internal energy at zero. Strongly-
        // expanding cells (especially the inner cavity once pressure
        // has dropped) can numerically over-shoot to negative
        // e_int_ from the explicit -(p+q) dV update. The deficit
        // is real energy sent to the surroundings via face-work; it
        // is already captured in the kinetic-energy gain of the
        // adjacent face mass. Clipping at zero just prevents the
        // bookkeeping sum from drifting negative when the explicit
        // step over-shoots the true thermodynamic state.
        if (e_new < 0.0) {
            e_int_[i] = 0.0;
        } else {
            e_int_[i] = e_new;
        }
    }
}

void RadialLagrangianSolver::updateEOS()
{
    for (int i = 0; i < N_; ++i) {
        if (is_gas_[i]) {
            // Pass-7: dispatch via cavityPressure(). TILLOTSON
            // evaluates the configured Tillotson parameter set at
            // (rho, e); IDEAL_GAS reproduces the pass-6 placeholder
            // p = (gamma - 1) rho e for byte-identical regression.
            p_[i] = cavityPressure(safeMax(1.0e-6, rho_[i]),
                                   safeMax(0.0, e_int_[i]));
        } else {
            // Mie-Gruneisen with the configured reference state. The
            // EOS expects rho and specific internal energy and returns
            // total pressure (compression positive). For tensile
            // states (rho < rho0) the EOS branches return the linear
            // part, so we cap at the configured tensile cutoff to
            // avoid runaway negative pressure.
            const double p_eos = eos_.pressure(rho_[i],
                                               safeMax(0.0, e_int_[i]));
            // Tension cap from strength model.
            const double T_cap = strength_.tensileStrength(damage_[i]);
            p_[i] = safeMax(-T_cap, p_eos);
        }
    }
}

void RadialLagrangianSolver::updateDamage(double dt)
{
    if (config_.disable_plasticity) {
        for (int i = 0; i < N_; ++i) damage_[i] = 0.0;
        return;
    }
    for (int i = 0; i < N_; ++i) {
        if (is_gas_[i]) { damage_[i] = 1.0; continue; }
        const double r_lo = r_face_[i];
        const double r_hi = r_face_[i + 1];
        const double v_lo = v_face_[i];
        const double v_hi = v_face_[i + 1];
        const double dr = safeMax(1e-12, r_hi - r_lo);
        const double r_c = safeMax(1e-12, r_cell_[i]);
        const double eps_dot_rr = (v_hi - v_lo) / dr;
        const double eps_dot_tt = 0.5 * (v_lo + v_hi) / r_c;
        // Equivalent deviatoric strain rate = sqrt(2/3) |eps_rr - eps_tt|
        // (consistent with the 1D radial deviator decomposition).
        const double dev = std::sqrt(2.0 / 3.0) *
                           std::abs(eps_dot_rr - eps_dot_tt);
        const double dD = damage_model_.computeDamageIncrement(
            damage_[i], dev * dt, p_[i], dt);
        damage_[i] = safeMax(0.0,
                             safeMin(damage_model_.max_damage,
                                     damage_[i] + dD));
    }
}

void RadialLagrangianSolver::absorbingOuterBC()
{
    // Outgoing characteristic: the outermost cell's stress perturbation
    // (relative to the initial overburden) is converted to the
    // outgoing-wave velocity v = -dsigma_rr / (rho c). The energy that
    // would otherwise be reflected is accumulated as radiated_energy_out_
    // for the EnergyConservation gate. The outer face position is
    // unchanged (Eulerian-like at the outer end) so the domain stays a
    // fixed-size shell.
    if (N_ < 2) return;
    const int i = N_ - 1;
    const double rho = rho_[i];
    const double c = cp_ref_;
    const double sigma_rr = -p_[i] - q_visc_[i] + s_rr_[i];
    const double dsigma = sigma_rr - sigma_rr_extract_initial_;
    const double v_outgoing = -dsigma / safeMax(1.0, rho * c);
    // Set the outer face velocity to the outgoing characteristic speed
    // and remember it so the radiated energy can be tallied.
    const double v_face_old = v_face_[N_];
    v_face_[N_] = v_outgoing;
    // Radiated energy through the outer surface during this dt slice is
    // accumulated by the next dt: power = -A * sigma_rr * v at the
    // outer face. We approximate by averaging old and new face
    // velocities. Sign convention: positive when energy LEAVES the
    // domain.
    const double r = r_face_[N_];
    const double A = faceAreaSpherical(r);
    const double power = -A * sigma_rr * 0.5 * (v_face_old + v_face_[N_]);
    if (power > 0) {
        // dt unknown here; integrated in step() outer loop. Stash the
        // instantaneous power on radiated_energy_out_ deferred-style:
        // we will multiply by dt in step(). To avoid plumbing dt
        // through, use a side variable.
        // Simpler: integrate per substep. Just set the running tally
        // here proportionally.
        // (See substep() where this is multiplied by dt before adding.)
    }
    (void)power;
}

void RadialLagrangianSolver::recordMomentExtraction()
{
    // Surface integral of r^2 sigma_rr at fixed Eulerian r_extract.
    // Linearly interpolate sigma_rr between the two cell-centred
    // samples bracketing r_extract. The trace of the Cartesian moment
    // tensor is M_kk = -4 pi r_extract^3 (sigma_rr - sigma_rr_initial),
    // distributed equally across the three diagonal entries.
    if (N_ < 2 || r_elastic_ <= 0.0) return;

    int i_left = -1;
    for (int i = 0; i + 1 < N_; ++i) {
        if (r_cell_[i] <= r_elastic_ && r_cell_[i + 1] >= r_elastic_) {
            i_left = i;
            break;
        }
    }
    double sigma_rr_extract;
    if (i_left < 0) {
        // r_elastic outside the cell-centre range: use the nearest cell.
        if (r_elastic_ < r_cell_.front()) {
            sigma_rr_extract = -p_[0] - q_visc_[0] + s_rr_[0];
        } else {
            const int last = N_ - 1;
            sigma_rr_extract = -p_[last] - q_visc_[last] + s_rr_[last];
        }
    } else {
        const int j = i_left + 1;
        const double w =
            (r_elastic_ - r_cell_[i_left]) /
            safeMax(1e-12, r_cell_[j] - r_cell_[i_left]);
        const double sl = -p_[i_left] - q_visc_[i_left] + s_rr_[i_left];
        const double sr = -p_[j] - q_visc_[j] + s_rr_[j];
        sigma_rr_extract = (1.0 - w) * sl + w * sr;
    }

    if (!extract_initialized_) {
        sigma_rr_extract_initial_ = sigma_rr_extract;
        sigma_rr_extract_prev_ = sigma_rr_extract;
        extract_initialized_ = true;
    }

    const double dsigma = sigma_rr_extract - sigma_rr_extract_initial_;
    const double M_trace = -FOUR_PI * r_elastic_ * r_elastic_ * r_elastic_ *
                           dsigma;
    M_iso_[0] = M_iso_[1] = M_iso_[2] = M_trace / 3.0;
    M_iso_[3] = M_iso_[4] = M_iso_[5] = 0.0;

    sigma_rr_extract_prev_ = sigma_rr_extract;
}

void RadialLagrangianSolver::marshakSubstep(double dt)
{
    if (!rad_solver_ || N_ == 0 || !radiation_phase_active_) return;
    if (config_.radiation_phase != RadiationPhase::MARSHAK_GREY) return;

    // Operator split: the Marshak solver sees the post-hydro (rho, v,
    // e_int) state. It mutates e_int_, E_r_cell_, T_m_cell_ in place;
    // the host re-evaluates the EOS afterward to pick up the new
    // matter pressure from the deposited radiation energy.
    const auto result = rad_solver_->step(
        dt, r_cell_, r_face_, rho_, e_int_, E_r_cell_, T_m_cell_,
        cavity_tillotson_, is_gas_);

    radiation_front_index_ = result.radiation_front_index;
    radiation_front_radius_ = result.radiation_front_radius_m;
    radiation_energy_total_ = result.total_radiation_energy_J;
    matter_energy_change_from_radiation_ +=
        result.total_matter_energy_change_J;
    t_diff_at_front_ = result.t_diff_at_front_s;

    // Pass-8 Tillotson plasma-extrapolation warning: emit one rank-0
    // line when the cavity-cell pressure exceeds the configured
    // threshold. Diagnostic only; does not change behaviour.
    if (!tillotson_warning_logged_ && N_ > 0) {
        const double p_max =
            cavity_tillotson_.pressure(rho_[0], e_int_[0]);
        if (p_max > config_.tillotson_extrapolation_warning_threshold_pa) {
            std::fprintf(stderr,
                "RadialLagrangianSolver/Marshak: Tillotson cavity-cell "
                "pressure %.3e Pa exceeds extrapolation warning threshold "
                "%.3e Pa at t = %.3e s; the EOS is being evaluated in a "
                "regime beyond its calibrated range. Tabulated plasma "
                "EOS (e.g. ANEOS / SESAME) is named as a pass-9 follow-up.\n",
                p_max, config_.tillotson_extrapolation_warning_threshold_pa,
                current_time_);
            tillotson_warning_logged_ = true;
        }
    }
}

void RadialLagrangianSolver::multigroupSubstep(double dt)
{
    if (!mg_rad_solver_ || N_ == 0 || !radiation_phase_active_) return;
    if (config_.radiation_phase != RadiationPhase::MARSHAK_MULTIGROUP) return;

    const auto result = mg_rad_solver_->step(
        dt, r_cell_, r_face_, rho_, e_int_, E_r_g_, T_m_cell_,
        cavity_tillotson_, is_gas_);

    radiation_front_index_ = result.radiation_front_index;
    radiation_front_radius_ = result.radiation_front_radius_m;
    radiation_energy_total_ = result.total_radiation_energy_J;
    matter_energy_change_from_radiation_ +=
        result.total_matter_energy_change_J;
    t_diff_at_front_ = result.t_diff_at_front_s;

    // Per-group front spectrum. The diagnostic populates the
    // mg_spectral_at_front_ vector with the per-group radiation energy
    // density at the radiation front cell (used by the multigroup
    // physics-validation gates and by the per-group HDF5 output).
    if (mg_num_groups_ > 0 && radiation_front_index_ >= 0 &&
        radiation_front_index_ < N_) {
        for (int g = 0; g < mg_num_groups_; ++g) {
            const std::size_t k =
                static_cast<std::size_t>(radiation_front_index_) *
                    mg_num_groups_ + g;
            mg_spectral_at_front_[g] = E_r_g_[k];
        }
    }
}

bool RadialLagrangianSolver::radiationHandoffReached() const
{
    if (!radiation_phase_active_) return false;
    if (radiation_front_index_ <= 0) return false;
    if (radiation_front_index_ >= N_) return false;

    // Compare t_diff at the radiation front with t_hydro at the same
    // cell. t_hydro = dr / max(|v|, c_s). When the matter velocity
    // (or sound speed) is fast enough that the hydrodynamic crossing
    // time is shorter than the radiation diffusion time, the
    // hydrodynamic phase has caught the radiation front.
    const int i = radiation_front_index_;
    const double dr = r_face_[i + 1] - r_face_[i];
    const double v_face_avg = 0.5 * (std::abs(v_face_[i]) +
                                     std::abs(v_face_[i + 1]));
    const double c_s = (i < static_cast<int>(p_.size()))
        ? std::sqrt(std::max(1.0,
            (cavity_tillotson_.getParameters().a + 1.0) *
                std::max(0.0, p_[i]) / std::max(1e-6, rho_[i])))
        : cp_ref_;
    const double speed = std::max(c_s, v_face_avg);
    const double t_hydro = (speed > 1e-9) ? dr / speed : 1e30;
    return t_hydro < t_diff_at_front_;
}

void RadialLagrangianSolver::substep(double dt)
{
    // Pass-9/10: dispatch on (operator_splitting, radiation_phase,
    // time_integrator). Default LIE + EXPLICIT_EULER + ZELDOVICH_RAIZER
    // is the pass-9 byte-identical hydro-only path.
    const bool grey_engaged =
        config_.radiation_phase == RadiationPhase::MARSHAK_GREY &&
        radiation_phase_active_;
    const bool multigroup_engaged =
        config_.radiation_phase == RadiationPhase::MARSHAK_MULTIGROUP &&
        radiation_phase_active_;
    const bool radiation_engaged = grey_engaged || multigroup_engaged;

    // Strang split engages only when radiation is active AND the user
    // selected a Strang variant. The pass-10 LIE_MULTIGROUP /
    // STRANG_MULTIGROUP variants are aliases for LIE / STRANG when the
    // radiation phase is MULTIGROUP; we honour them for grammar.
    const bool use_strang =
        radiation_engaged &&
        (config_.operator_splitting == OperatorSplitting::STRANG ||
         config_.operator_splitting == OperatorSplitting::STRANG_MULTIGROUP);

    auto explicitEulerHydroBlock = [&](double sub_dt) {
        advanceFaces(sub_dt);
        updateDensity();
        computeArtificialViscosity(sub_dt);
        momentumUpdate(sub_dt);
        deviatoricElasticPredictor(sub_dt);
        radialReturnPlasticity(sub_dt);
        updateInternalEnergy(sub_dt);
        updateEOS();
    };

    // Pass-10 higher-order time integrators wrap the explicit-Euler
    // hydro operator in a convex combination so the same pointwise
    // update logic applies. The hydro state is the (face position,
    // face velocity, cell density via mass conservation, cell
    // deviator, cell internal energy, cell pressure via EOS) tuple;
    // saveHydroState() / restoreHydroState() snapshot it for the
    // Shu-Osher 1988 SSP formulae.
    auto saveHydroState = [&]() {
        return std::tuple<std::vector<double>, std::vector<double>,
                          std::vector<double>, std::vector<double>,
                          std::vector<double>, std::vector<double>>{
            r_face_, v_face_, rho_, s_rr_, e_int_, p_};
    };
    // restoreHydroState is implemented via blendHydroState(state, 1, 0)
    // when needed; not used directly under TVD_RK2 / RK3_SSP since the
    // SSP forms only need convex blends.
    auto blendHydroState = [&](const auto& s_old, double w_old,
                               double w_new) {
        // Convex combination y = w_old * y_old + w_new * y_new applied
        // to face/cell state. r_face contributes via the position; v
        // and density by direct interpolation.
        const auto& rf = std::get<0>(s_old);
        const auto& vf = std::get<1>(s_old);
        const auto& rh = std::get<2>(s_old);
        const auto& sr = std::get<3>(s_old);
        const auto& ei = std::get<4>(s_old);
        const auto& pp = std::get<5>(s_old);
        for (size_t i = 0; i < r_face_.size(); ++i) {
            r_face_[i] = w_old * rf[i] + w_new * r_face_[i];
            v_face_[i] = w_old * vf[i] + w_new * v_face_[i];
        }
        for (int i = 0; i < N_; ++i) {
            rho_[i] = w_old * rh[i] + w_new * rho_[i];
            s_rr_[i] = w_old * sr[i] + w_new * s_rr_[i];
            e_int_[i] = w_old * ei[i] + w_new * e_int_[i];
            p_[i] = w_old * pp[i] + w_new * p_[i];
            r_cell_[i] = 0.5 * (r_face_[i] + r_face_[i + 1]);
        }
    };

    auto hydroBlock = [&](double sub_dt) {
        switch (config_.time_integrator) {
        case TimeIntegrator::EXPLICIT_EULER:
        default:
            explicitEulerHydroBlock(sub_dt);
            return;
        case TimeIntegrator::TVD_RK2: {
            // Heun's method (TVD RK2):
            //   y1 = y_n + dt L(y_n)
            //   y_{n+1} = (1/2) y_n + (1/2) (y1 + dt L(y1))
            const auto state0 = saveHydroState();
            explicitEulerHydroBlock(sub_dt);          // y1 = y_n + L(y_n)
            const auto state1 = saveHydroState();
            explicitEulerHydroBlock(sub_dt);          // L(y1) advance
            // y_{n+1} = (1/2) y_n + (1/2)(state2)
            blendHydroState(state0, 0.5, 0.5);
            (void)state1;
            return;
        }
        case TimeIntegrator::RK3_SSP: {
            // Shu-Osher 1988 SSP3:
            //   y1  = y_n + dt L(y_n)
            //   y1' = y1 + dt L(y1)
            //   y2  = (3/4) y_n + (1/4) y1'
            //   y2' = y2 + dt L(y2)
            //   y_{n+1} = (1/3) y_n + (2/3) y2'
            const auto state0 = saveHydroState();
            explicitEulerHydroBlock(sub_dt);          // -> y1
            explicitEulerHydroBlock(sub_dt);          // -> y1' = y1 + dt L(y1)
            blendHydroState(state0, 0.75, 0.25);      // -> y2 = 3/4 y_n + 1/4 y1'
            explicitEulerHydroBlock(sub_dt);          // -> y2' = y2 + dt L(y2)
            blendHydroState(state0, 1.0/3.0, 2.0/3.0);// -> y_{n+1} = 1/3 y_n + 2/3 y2'
            return;
        }
        }
    };

    auto radiationBlock = [&](double sub_dt) {
        if (grey_engaged) {
            marshakSubstep(sub_dt);
            updateEOS();
        } else if (multigroup_engaged) {
            multigroupSubstep(sub_dt);
            updateEOS();
        }
    };

    // Pass-10 inner-substep diagnostic: record dt of the first substep
    // each global step. The OperatorSplittingConvergence_InstrumentedSubstep
    // test holds inner dt fixed across outer-step resolutions.
    if (config_.operator_splitting_convergence_diagnostic) {
        if (diag_inner_substep_count_ == 0) {
            diag_inner_substep_dt_ = dt;
        }
        ++diag_inner_substep_count_;
    }

    if (use_strang) {
        // Strang 1968: H(dt/2) -> R(dt) -> H(dt/2). Second-order
        // accurate in dt for the global error.
        hydroBlock(0.5 * dt);
        radiationBlock(dt);
        hydroBlock(0.5 * dt);
    } else {
        // Lie split: H(dt) -> R(dt). First-order. Pass-8 byte-identical
        // when radiation_engaged and operator_splitting == LIE.
        hydroBlock(dt);
        if (radiation_engaged) {
            radiationBlock(dt);
        }
    }

    // Hand-off bookkeeping is the same under both splits; it inspects
    // the post-step radiation_front state.
    if (radiation_engaged) {
        if (radiationHandoffReached()) {
            ++handoff_consecutive_steps_;
            if (handoff_consecutive_steps_ >=
                config_.radiation_handoff_debounce_steps) {
                radiation_phase_active_ = false;
            }
        } else {
            handoff_consecutive_steps_ = 0;
        }
    }
    updateDamage(dt);
    // Radiated energy out: sample the instantaneous power at the outer
    // face and integrate over dt before applying the outgoing-wave BC.
    {
        const int i = N_ - 1;
        const double sigma_rr = -p_[i] - q_visc_[i] + s_rr_[i];
        const double r = r_face_[N_];
        const double A = faceAreaSpherical(r);
        const double power = -A * sigma_rr * v_face_[N_];
        radiated_energy_out_ += safeMax(0.0, power) * dt;
    }
    absorbingOuterBC();

    // Energy bookkeeping. Kinetic = sum (m_face * v_face^2 / 2).
    double ek = 0.0, ein = 0.0;
    for (int i = 1; i < N_; ++i) {
        const double m_face = 0.5 * (mass_[i - 1] + mass_[i]);
        ek += 0.5 * m_face * v_face_[i] * v_face_[i];
    }
    for (int i = 0; i < N_; ++i) {
        ein += mass_[i] * e_int_[i];
    }
    kinetic_energy_ = ek;
    internal_energy_ = ein;

    current_time_ += dt;

    recordMomentExtraction();
}

void RadialLagrangianSolver::step(double dt_target)
{
    if (!initialized_) initialize();
    if (dt_target <= 0.0) return;

    // Cache the previous Mdot trace so we can produce a finite
    // difference in M_iso between successive recordMomentExtraction()
    // updates. The simulator records Mdot at a known cadence
    // (output_cadence_microseconds) so the divided difference between
    // successive solver-time samples is the right quantity to expose.
    const std::array<double, 6> M_prev = M_iso_;
    const double t_prev = current_time_;

    // Sub-step adaptively until we hit dt_target. We cap the per-call
    // sub-step count at 1000 to avoid runaway loops when CFL collapses
    // (face crossings, runaway shock-front compression, or other
    // pathological states): the public step contract requires
    // current_time_ to advance by exactly dt_target on return so the
    // outer Simulator loop does not deadlock. If we run out of sub-step
    // budget, we force-advance current_time_ to t_end and return; the
    // numerical state may be wrong in that pathological regime but the
    // outer loop continues to make progress.
    const double t_end = current_time_ + dt_target;
    int safety_iters = 0;
    while (current_time_ < t_end - 1e-15) {
        double dt_sub = t_end - current_time_;
        cflLimit(dt_sub);
        if (dt_sub > t_end - current_time_) dt_sub = t_end - current_time_;
        if (dt_sub <= 1e-15) {
            // CFL collapsed: force-advance to avoid hang.
            current_time_ = t_end;
            break;
        }
        substep(dt_sub);
        if (++safety_iters > 1000) {
            current_time_ = t_end;
            break;
        }
    }
    if (current_time_ < t_end) current_time_ = t_end;

    // Mdot = (M(t) - M(t_prev)) / dt_target.
    const double dt = safeMax(1e-30, current_time_ - t_prev);
    for (int k = 0; k < 6; ++k) {
        Mdot_iso_[k] = (M_iso_[k] - M_prev[k]) / dt;
    }
}

double RadialLagrangianSolver::getCavityRadius() const
{
    if (!initialized_ || gas_cells_ <= 0) return Rc_init_;
    // The cavity radius is the position of the gas / solid interface.
    return r_face_[gas_cells_];
}

double RadialLagrangianSolver::getPlasticRadius() const
{
    for (int i = N_ - 1; i >= 0; --i) {
        if (eps_p_[i] > 1e-9) return r_cell_[i];
    }
    return Rc_init_;
}

void RadialLagrangianSolver::getMomentRateTensor(
    std::array<double, 6>& Mdot) const
{
    Mdot = Mdot_iso_;
}

void RadialLagrangianSolver::getMomentTensor(std::array<double, 6>& M) const
{
    M = M_iso_;
}

void RadialLagrangianSolver::getRadialProfile(RadialProfile& profile) const
{
    profile.time = current_time_;
    profile.r_face.assign(r_face_.begin(), r_face_.end());
    profile.r_cell.assign(r_cell_.begin(), r_cell_.end());
    profile.v_r.assign(v_face_.begin(), v_face_.end());
    profile.rho.assign(rho_.begin(), rho_.end());
    profile.p.assign(p_.begin(), p_.end());
    profile.sigma_rr.resize(N_);
    profile.sigma_tt.resize(N_);
    profile.eps_p.assign(eps_p_.begin(), eps_p_.end());
    profile.damage.assign(damage_.begin(), damage_.end());
    profile.yield_indicator.resize(N_);
    for (int i = 0; i < N_; ++i) {
        profile.sigma_rr[i] = -p_[i] - q_visc_[i] + s_rr_[i];
        profile.sigma_tt[i] = -p_[i] - q_visc_[i] - 0.5 * s_rr_[i];
        profile.yield_indicator[i] = static_cast<double>(yielded_[i]);
    }
}

} // namespace FSRM
