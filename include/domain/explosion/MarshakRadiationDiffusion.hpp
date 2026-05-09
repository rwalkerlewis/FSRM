/**
 * @file MarshakRadiationDiffusion.hpp
 * @brief 1D radial grey radiation-diffusion solver, coupled to matter
 *        via emission-absorption with a Newton iteration on T_m^4.
 *
 * Pass-8 (axis 1, see docs/HISTORIC_NUCLEAR_ROADMAP.md) replaces the
 * pass-7 Zel'dovich-Raizer end-state approximation with an explicit
 * numerical solve of the actual radiation-transport phase. The
 * radiation phase runs from t = 0 (yield deposition) to t = t_rh
 * (radiation-to-hydrodynamic transition); the host RadialLagrangianSolver
 * operator-splits per global timestep, advancing the hydro state with
 * the existing Wilkins-AV Lagrangian update and then calling this
 * solver to advance the (E_r, T_m) coupled system.
 *
 * Equations (grey radiation diffusion in 1D radial spherical symmetry):
 *
 *   dE_r/dt = (1/r^2) d/dr [ r^2 (c / (3 kappa_R rho)) dE_r/dr ]
 *           + c kappa_P rho ( a T_m^4 - E_r )
 *
 *   rho cv dT_m/dt = c kappa_P rho ( E_r - a T_m^4 )      (per-cell, no
 *                                                          spatial coupling
 *                                                          here; advection
 *                                                          is folded into
 *                                                          the host hydro
 *                                                          substep)
 *
 * with a = 4 sigma_SB / c the radiation constant.
 *
 * Discretisation:
 *
 *  - Spatial: cell-centred FVM on the host's Lagrangian radial mesh.
 *    The diffusion coefficient D = c / (3 kappa_R rho) is evaluated at
 *    cell faces by the harmonic mean of the two adjacent cell values
 *    (standard FVM choice for a diffusive flux).
 *
 *  - Temporal: backward-Euler in E_r. The nonlinear coupling term
 *    ( a T_m^4 - E_r ) is linearised around the current Newton iterate
 *    using d(T^4)/dT = 4 T^3 to write the matter-energy equation as
 *    a closed-form per-cell expression for T_m_new. We then assemble
 *    a tridiagonal system in E_r (the diffusion + linearised source)
 *    and solve it with the Thomas algorithm. Outer Newton iteration
 *    drives the residual on (E_r, T_m^4) below newton_tolerance.
 *
 *  - Boundary conditions:
 *     - inner face (r = 0): zero-flux symmetry boundary, dE_r/dr = 0,
 *       enforced by collapsing the first row of the tridiagonal.
 *     - outer face (r = r_outer): Marshak / extrapolation boundary,
 *       E_r = a T_m^4 at the outer cell with T_m at ambient (300 K).
 *       The radiation flux out is therefore ~ negligible until the
 *       front reaches the outer cell, at which point the host should
 *       have already triggered the radiation-to-hydro hand-off.
 *
 *  - Hand-off detection: returned in StepResult. The host computes
 *    t_diffusion = (dr)^2 rho / (D) and t_hydro = dr / max(|v|, c_s)
 *    at the radiation-front cell (the outermost cell with E_r >
 *    1.5 a T_amb^4) and triggers hand-off when t_hydro < t_diffusion
 *    for radiation_handoff_debounce_steps consecutive substeps.
 *
 * State convention. SI units throughout. E_r in J/m^3; T_m in Kelvin;
 * kappa in m^2/kg.
 *
 * References.
 *  - Pomraning, G. C. (1973), "The Equations of Radiation
 *    Hydrodynamics", Pergamon Press, ch IV (grey diffusion limit;
 *    self-similar Marshak wave validation target).
 *  - Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
 *    Hydrodynamics", Oxford University Press, sec 96-97 (operator
 *    splitting between matter and radiation; backward-Euler stability
 *    for stiff coupling).
 *  - Marshak, R. E. (1958), "Effect of radiation on shock wave
 *    behavior", Physics of Fluids 1(1), pp. 24-29 (boundary
 *    conditions).
 *  - Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
 *    Waves and High-Temperature Hydrodynamic Phenomena", vol I ch V
 *    (Kramers' opacity; vol II ch X (radiation-to-hydrodynamic
 *    transition criterion).
 *  - Larsen, E. W. (1988), "A grey transport acceleration method for
 *    time-dependent radiative transfer", J. Comp. Phys 78, pp. 459-480
 *    (linearised Newton on T^4 closure).
 */

#ifndef NEAR_FIELD_MARSHAK_RADIATION_DIFFUSION_HPP
#define NEAR_FIELD_MARSHAK_RADIATION_DIFFUSION_HPP

#include <vector>

#include "domain/explosion/OpacityModel.hpp"
#include "domain/explosion/TillotsonEOS.hpp"

namespace FSRM {

/**
 * @brief Grey radiation-diffusion solver for the 1D radial Lagrangian
 *        host. Decoupled from RadialLagrangianSolver: the host calls
 *        step() with the current cell-centred state and the solver
 *        mutates E_r, T_m, and e_int in place.
 */
class MarshakRadiationDiffusionSolver
{
public:
    struct Config
    {
        int max_newton_iter = 10;
        double newton_tolerance = 1.0e-6;
        OpacityModel opacity_model = OpacityModel::POWER_LAW_ZR;
        PowerLawOpacityParameters opacity_params =
            PowerLawOpacitySets::granite();
        /// CONSTANT-only: user-supplied opacity; used to exercise the
        /// Marshak self-similar physics gate independent of the Z-R
        /// power-law parameterization.
        double kappa_constant_m2_per_kg = 0.0;
        /// Ambient (background) matter temperature [K]. Outer-boundary
        /// reservoir temperature for the Marshak BC; also the cold-
        /// rock seed for cells outside the radiation front.
        double T_ambient_K = 300.0;
        /// Radiation-front detection threshold: a cell is considered
        /// "in the radiation front" when E_r > front_factor * a T_amb^4.
        double front_factor = 1.5;
    };

    /// Per-step diagnostic returned to the host RadialLagrangianSolver.
    struct StepResult
    {
        bool converged = true;
        int newton_iters = 0;
        double residual_inf_norm = 0.0;
        /// Outermost cell with E_r > front_factor * a T_amb^4.
        int radiation_front_index = 0;
        double radiation_front_radius_m = 0.0;
        /// t_diff = (dr)^2 / D evaluated at the radiation front. The
        /// host uses this to compare with t_hydro and trigger
        /// radiation-to-hydrodynamic hand-off.
        double t_diff_at_front_s = 0.0;
        /// Sum of E_r over all cells weighted by cell volume [J].
        double total_radiation_energy_J = 0.0;
        /// Net change in matter internal energy from the radiation
        /// substep alone [J]. Sign: positive when radiation deposits
        /// energy in the matter (cooling of the radiation field).
        double total_matter_energy_change_J = 0.0;
    };

    MarshakRadiationDiffusionSolver() = default;

    void setConfig(const Config& cfg) { config_ = cfg; }
    const Config& getConfig() const { return config_; }

    /// Build internal workspace for N cells. Idempotent.
    void initialize(int N);

    /// Advance the (E_r, T_m, e_int) coupled system by dt on the
    /// host's current radial mesh. Vectors are length N (cell-centred)
    /// or N+1 (face). On return:
    ///  - E_r is mutated to the post-step value.
    ///  - T_m is mutated to the post-step matter temperature.
    ///  - e_int is incremented by the per-cell matter-energy change
    ///    (in J/kg).
    /// The Tillotson EOS reference is used to look up the cell-wise
    /// effective heat capacity dT/de = 1 / (rho cv_eff). For the
    /// pass-8 implementation we use cv from the TillotsonParameters
    /// directly; an explicit cv(rho, e) lookup is named as a pass-9
    /// candidate.
    StepResult step(double dt,
                    const std::vector<double>& r_cell,
                    const std::vector<double>& r_face,
                    const std::vector<double>& rho,
                    std::vector<double>& e_int,
                    std::vector<double>& E_r,
                    std::vector<double>& T_m,
                    const TillotsonEOS& eos,
                    const std::vector<int>& is_gas);

    /// Diagnostic: opacity at (rho, T) under the configured model.
    /// Returns (kappa_R, kappa_P) in m^2/kg. Used by the
    /// OpacityRegimeCoverage gate.
    void evaluateOpacity(double rho, double T,
                         double& kappa_R, double& kappa_P) const;

private:
    Config config_;
    int N_ = 0;
    PowerLawOpacity power_law_;

    // Workspace for the tridiagonal system.
    std::vector<double> a_lower_;
    std::vector<double> b_diag_;
    std::vector<double> c_upper_;
    std::vector<double> rhs_;
    std::vector<double> x_;

    // Workspace for the Newton outer iteration.
    std::vector<double> T_iter_;
    std::vector<double> E_iter_;
    std::vector<double> kappa_R_;
    std::vector<double> kappa_P_;
    std::vector<double> D_face_;
    std::vector<double> e_int_old_;
    std::vector<double> E_r_old_;
    std::vector<double> T_m_old_;

    /// Recompute opacities per cell from rho and current T iterate.
    void updateOpacities(const std::vector<double>& rho,
                         const std::vector<double>& T_iter);

    /// Recompute the face diffusion coefficient D = c/(3 kappa_R rho)
    /// using the harmonic mean of the adjacent cell values.
    void updateFaceDiffusionCoefficient(const std::vector<double>& rho);

    /// Assemble the tridiagonal operator for E_r at the current Newton
    /// iterate. The linearised source term contributes -beta_i E_i to
    /// the diagonal and +beta_i a (4 T_iter_i^3 T_new_i - 3 T_iter_i^4)
    /// (from d(T^4)/dT = 4T^3) to the RHS, with beta = c kappa_P rho dt.
    /// The matter-energy equation is solved per-cell after the
    /// E_r solve and contributes to the RHS for the next Newton iter.
    void assembleTridiagonal(double dt,
                             const std::vector<double>& r_cell,
                             const std::vector<double>& r_face,
                             const std::vector<double>& rho);

    /// Thomas algorithm in-place. Reads a_lower_, b_diag_, c_upper_,
    /// rhs_; writes x_. All four work arrays are mutated.
    void solveTridiagonal();

    /// After the E_r tridiagonal solve, update T_m per cell from the
    /// linearised matter-energy balance. Uses T_iter as the
    /// linearisation point; iterating outer to convergence.
    void updateMatterTemperature(double dt,
                                 const std::vector<double>& rho,
                                 const TillotsonEOS& eos);

    /// Helpers
    static double cellVolumeSpherical(double r_lo, double r_hi);
    static double faceAreaSpherical(double r);
    static double harmonicMean(double a, double b);
};

} // namespace FSRM

#endif // NEAR_FIELD_MARSHAK_RADIATION_DIFFUSION_HPP
