/**
 * @file MultigroupRadiationDiffusion.hpp
 * @brief Pass-10 (axis 1a closeout) multigroup radiation-diffusion
 *        solver. Frequency-discretized counterpart to the pass-8 grey
 *        solver in MarshakRadiationDiffusion.hpp.
 *
 * Equations (multigroup radiation diffusion in 1D radial spherical
 * symmetry, for group g in [0, G)).
 *
 *   dE_r^g/dt = (1/r^2) d/dr [ r^2 (c / (3 kappa_R^g rho)) dE_r^g/dr ]
 *             + c kappa_P^g rho ( 4 pi B_g(T_m) / c - E_r^g )
 *
 *   rho cv dT_m/dt = sum_g c kappa_P^g rho ( E_r^g - 4 pi B_g(T_m) / c )
 *                                                       (per-cell, no
 *                                                        spatial
 *                                                        coupling)
 *
 * with B_g(T) = int_{nu_g}^{nu_{g+1}} B(nu, T) dnu the Planck-integrated
 * blackbody source per group. The matter energy equation sums emission
 * and absorption across all G groups so the matter temperature couples
 * every group to every other group through the Newton iteration, even
 * though the per-group spatial diffusion is uncoupled.
 *
 * Discretisation:
 *
 *  - Per-group spatial: identical to the grey solver, conservative
 *    cell-centred FVM with harmonic-mean face diffusion. G separate
 *    tridiagonal systems per Newton iteration.
 *
 *  - Temporal: backward-Euler in E_r^g, identical structure to the
 *    grey solve. The nonlinear B_g(T_m) source is linearised around
 *    the current Newton iterate using d B_g / dT (computed by Simpson
 *    quadrature). Outer Newton iterates until the relative residual on
 *    (E_r^g, T_m) all drop below newton_tolerance.
 *
 *  - Per-group opacity: computed lazily at runtime from
 *    MultigroupOpacityEvaluator (Simpson quadrature of the analytic
 *    Mihalas-Mihalas + Kramers + Thomson model). Per Newton iter the
 *    solver re-evaluates kappa_R^g(rho_i, T_iter_i) and
 *    kappa_P^g(rho_i, T_iter_i) for every cell and every group. The
 *    cost is dominated by this evaluation; default parameters
 *    (G = 16, n_simpson = 17, N = 200) are tractable in CI budgets.
 *
 *  - Boundary conditions: per-group Marshak BC at the outer face
 *    (E_r^g = 4 pi B_g(T_amb) / c at the outer ghost), zero-flux
 *    symmetry at the inner face.
 *
 * State convention. SI units. E_r^g in J/m^3, T_m in K, kappa^g in
 * m^2/kg, B_g(T) in W/(m^2 sr) (band-integrated; a 4 pi factor is
 * applied internally to convert to the source term coefficient).
 *
 * Residuals (used by the gates listed in the pass-10 spec).
 *  - Per-group E_r relative change between Newton iters.
 *  - Per-cell T_m relative change between Newton iters.
 *  - Total energy E_radiation_total + matter internal energy for the
 *    energy-conservation gate.
 *
 * References.
 *  - Pomraning, G. C. (1973), "The Equations of Radiation
 *    Hydrodynamics", Pergamon Press, ch IV (multigroup formulation,
 *    operator-split iteration).
 *  - Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
 *    Hydrodynamics", Oxford University Press, sec 80 (multigroup
 *    diffusion), sec 82.2 (line opacity smoothed continuum).
 *  - Larsen, E. W. (1988), "A grey transport acceleration method for
 *    time-dependent radiative transfer", J. Comp. Phys 78, pp 459-480
 *    (linearised Newton on T^4 closure; pass-8 generalisation).
 *  - Strang, G. (1968), "On the construction and comparison of
 *    difference schemes", SIAM J. Num. Anal. 5(3), pp 506-517 (operator
 *    splitting carried over from the grey solver).
 */

#ifndef NEAR_FIELD_MULTIGROUP_RADIATION_DIFFUSION_HPP
#define NEAR_FIELD_MULTIGROUP_RADIATION_DIFFUSION_HPP

#include <vector>

#include "domain/explosion/MultigroupOpacity.hpp"
#include "domain/explosion/OpacityModel.hpp"
#include "domain/explosion/TillotsonEOS.hpp"

namespace FSRM {

/**
 * @brief Multigroup radiation-diffusion solver in spherical 1D.
 *
 * Allocated and configured once per host RadialLagrangianSolver
 * lifecycle. Per-step the host calls step() with the current (rho,
 * e_int, E_r^g, T_m) state; the solver mutates E_r^g, T_m, and e_int
 * in place and returns a StepResult with diagnostics for the host's
 * radiation-to-hydrodynamic hand-off logic.
 */
class MultigroupRadiationDiffusionSolver
{
public:
    struct Config
    {
        FrequencyGroupGrid group_grid;          ///< G, nu range, simpson points.
        int max_newton_iter = 10;
        double newton_tolerance = 1.0e-6;
        /// Baseline POWER_LAW_ZR parameters (kappa_0, exponents) used
        /// to scale the analytic frequency-dependent opacity. Pass-10
        /// path A: per-group opacity is derived from this baseline by
        /// Simpson quadrature of the Mihalas-Mihalas shape.
        PowerLawOpacityParameters opacity_params =
            PowerLawOpacitySets::granite();
        /// Ambient (background) matter temperature [K]. Outer-boundary
        /// reservoir temperature for the Marshak BC; cold-rock seed
        /// for cells outside the radiation front.
        double T_ambient_K = 300.0;
        /// Radiation-front detection threshold: a cell is "in the
        /// radiation front" when sum_g E_r^g > front_factor * a T_amb^4.
        double front_factor = 1.5;
        /// Per-group CONSTANT opacity override [m^2/kg]. When > 0,
        /// every per-group Rosseland and Planck mean is set to this
        /// value, bypassing the analytic Mihalas-Mihalas + Kramers
        /// path. Used by the SelfSimilarPureRadiation gate to pin the
        /// per-group diffusion coefficient at a known constant so the
        /// analytic Marshak self-similar comparison applies cleanly.
        double kappa_constant_m2_per_kg = 0.0;
    };

    struct StepResult
    {
        bool converged = true;
        int newton_iters = 0;
        double residual_inf_norm = 0.0;
        int radiation_front_index = 0;
        double radiation_front_radius_m = 0.0;
        double t_diff_at_front_s = 0.0;
        double total_radiation_energy_J = 0.0;
        double total_matter_energy_change_J = 0.0;
    };

    MultigroupRadiationDiffusionSolver() = default;

    void setConfig(const Config& cfg);
    const Config& getConfig() const { return config_; }
    int numGroups() const { return config_.group_grid.n_groups; }

    /// Build internal workspace for N cells and G groups.
    /// Idempotent: re-init resets Newton workspace but does not mutate
    /// previously-stored E_r^g or T_m on the caller side.
    void initialize(int N);

    /// Advance the (E_r^g, T_m, e_int) coupled system by dt. The
    /// frequency-resolved radiation field E_r is laid out
    /// row-major [cell, group]: E_r[i*G + g]. The host vectors must
    /// already be sized: rho/e_int/T_m of size N, E_r of size N*G,
    /// is_gas of size N.
    StepResult step(double dt,
                    const std::vector<double>& r_cell,
                    const std::vector<double>& r_face,
                    const std::vector<double>& rho,
                    std::vector<double>& e_int,
                    std::vector<double>& E_r,
                    std::vector<double>& T_m,
                    const TillotsonEOS& eos,
                    const std::vector<int>& is_gas);

    /// Per-group accessors useful for tests and for the host
    /// per-group HDF5 output.
    const MultigroupOpacityEvaluator& opacityEvaluator() const
    {
        return opacity_;
    }
    /// Layout: E_r[i * G + g] in row-major form (cell major).
    static int idx(int cell, int group, int G) { return cell * G + group; }

private:
    Config config_;
    int N_ = 0;
    int G_ = 0;
    MultigroupOpacityEvaluator opacity_;

    // Workspace.
    std::vector<double> a_lower_, b_diag_, c_upper_, rhs_, x_;
    std::vector<double> kappa_R_, kappa_P_;   ///< per group, per cell, contiguous
    std::vector<double> D_face_;              ///< (N+1) per group, contiguous
    std::vector<double> T_iter_;
    std::vector<double> T_m_old_;
    std::vector<double> E_iter_;              ///< (N*G)
    std::vector<double> E_r_old_;             ///< (N*G)
    std::vector<double> e_int_old_;
    std::vector<double> B_g_iter_;            ///< (N*G) B_g at T_iter
    std::vector<double> dBg_dT_iter_;         ///< (N*G) dBg/dT at T_iter

    void recomputeOpacitiesAndPlanck(const std::vector<double>& rho);
    void updateFaceDiffusion(const std::vector<double>& rho, int g);
    void assembleTridiagonal(double dt, int g,
                             const std::vector<double>& r_cell,
                             const std::vector<double>& r_face,
                             const std::vector<double>& rho);
    void solveTridiagonal();
    void updateMatterTemperature(double dt,
                                 const std::vector<double>& rho,
                                 const TillotsonEOS& eos);

    static double cellVolumeSpherical(double r_lo, double r_hi);
    static double faceAreaSpherical(double r);
    static double harmonicMean(double a, double b);
};

}  // namespace FSRM

#endif  // NEAR_FIELD_MULTIGROUP_RADIATION_DIFFUSION_HPP
