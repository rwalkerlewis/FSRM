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
void RadialLagrangianSolver::setConfig(const Config& c) { config_ = c; }

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

    // The host-medium NTS cavity radius sets both the initial cavity
    // size and the elastic radius at which the moment-tensor surface
    // integral is evaluated. The configurable radial_outer_factor
    // multiplies the elastic radius to give the outer-domain extent so
    // the outgoing wave has room to leave the source region before the
    // characteristic BC fires.
    const double Rc_eq = safeMax(1.0, src_.cavityRadius());
    Rc_init_ = 0.4 * Rc_eq;
    r_elastic_ = config_.radial_outer_factor * Rc_eq;
    r_outer_ = safeMax(r_elastic_ * 1.20, 1.5 * Rc_eq);

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
    const double overburden = src_.overburden_stress;
    for (int i = 0; i < N; ++i) {
        const double r_lo = r_face_[i];
        const double r_hi = r_face_[i + 1];
        const double V = cellVolumeSpherical(r_lo, r_hi);
        if (i < gas_cells_) {
            is_gas_[i] = 1;
            // Detonation gas. Yield-derived total energy distributed
            // uniformly over the cavity volume. Pass-6 uses an ideal
            // gas EOS with the configurable adiabatic index; pass-7
            // should replace this with a JWL detonation-products EOS.
            // (See HISTORIC_NUCLEAR_ROADMAP.md axis 1, follow-up.)
            // Glasstone & Dolan (1977) tabulates the gas-cavity
            // partition for nuclear yields at ~25 percent of total
            // energy; for the radial-shock-driving role here the full
            // E_yield deposited as gas internal energy is a defensible
            // first approximation since plasticity / radiation losses
            // accumulate in the surrounding solid cells over the
            // simulated time.
            rho_[i] = src_.host_density;  // inertia continuity at t=0
            mass_[i] = rho_[i] * V;
            const double cavity_volume = (FOUR_PI / 3.0) * Rc_init_ * Rc_init_ * Rc_init_;
            const double e_specific = src_.energyJoules() /
                                      safeMax(1e-12, src_.host_density * cavity_volume);
            e_int_[i] = e_specific;
            p_[i] = (config_.gas_eos_gamma - 1.0) * rho_[i] * e_int_[i];
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

    current_time_ = 0.0;
    initialized_ = true;
}

void RadialLagrangianSolver::cflLimit(double& dt) const
{
    double dt_min = dt;
    for (int i = 0; i < N_; ++i) {
        const double dr = r_face_[i + 1] - r_face_[i];
        if (dr <= 0.0) continue;
        // Effective sound speed: solid uses vp; gas cell uses
        // sqrt(gamma * p / rho). Both bounded by 0.1 * cp_ref so a
        // freshly-initialized state with p = 0 in solid cells is not
        // CFL-unbounded.
        double cs;
        if (is_gas_[i]) {
            const double p_pos = safeMax(0.0, p_[i]);
            cs = std::sqrt(safeMax(1.0, config_.gas_eos_gamma * p_pos /
                                          safeMax(1e-6, rho_[i])));
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
            const double cs =
                is_gas_[i]
                    ? std::sqrt(safeMax(
                          1.0, config_.gas_eos_gamma * safeMax(0.0, p_[i]) /
                                   safeMax(1e-6, rho_[i])))
                    : cp_ref_;
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
        e_int_[i] += de;
    }
}

void RadialLagrangianSolver::updateEOS()
{
    for (int i = 0; i < N_; ++i) {
        if (is_gas_[i]) {
            // Ideal gas. p = (gamma - 1) rho e.
            const double p_gas = (config_.gas_eos_gamma - 1.0) *
                                 safeMax(1e-6, rho_[i]) *
                                 safeMax(0.0, e_int_[i]);
            p_[i] = safeMax(0.0, p_gas);
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

void RadialLagrangianSolver::substep(double dt)
{
    advanceFaces(dt);
    updateDensity();
    computeArtificialViscosity(dt);
    momentumUpdate(dt);
    deviatoricElasticPredictor(dt);
    radialReturnPlasticity(dt);
    updateInternalEnergy(dt);
    updateEOS();
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
