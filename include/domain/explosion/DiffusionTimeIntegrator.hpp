/**
 * @file DiffusionTimeIntegrator.hpp
 * @brief Pass-11 (axis 1c) time-integrator strategy for the per-group
 *        radiation-diffusion solve. Replaces the pass-10 hardcoded
 *        backward-Euler assembly in MarshakRadiationDiffusion.cpp and
 *        MultigroupRadiationDiffusion.cpp with a strategy interface
 *        that the host solvers parameterise by a small set of
 *        time-stepping coefficients.
 *
 * The pass-10 multigroup and grey solvers assemble a per-cell-per-group
 * tridiagonal system of the form
 *
 *   diag = V_i + dt*(a_e + a_w + c_kp_rho_V)
 *   sub  = -dt*a_w
 *   sup  = -dt*a_e
 *   rhs  = V_i*E^n + dt*c_kp_rho_V*S_emit(T_iter)
 *
 * Pass-11 abstracts the per-step coefficients into a small Coefficients
 * struct. The host solvers compute V_i, a_e, a_w, c_kp_rho_V, S_emit
 * from the spatial discretisation and the current Newton iterate and
 * combine them with the integrator's coefficients to assemble the
 * actual per-cell linear system.
 *
 * Backward-Euler (pass-10 default, byte-identical):
 *
 *   diag = V_i + dt*(a_e + a_w + c_kp_rho_V)
 *   rhs  = V_i*E^n + dt*c_kp_rho_V*S_emit(T_iter)
 *
 * Crank-Nicolson (Larsen 1988 sec 3 time-centred linearisation, second-order):
 *
 *   diag = V_i + (dt/2)*(a_e + a_w + c_kp_rho_V)
 *   sub  = -(dt/2)*a_w
 *   sup  = -(dt/2)*a_e
 *   rhs  = V_i*E^n
 *        + (dt/2)*[a_e*(E^n_{i+1} - E^n_i) - a_w*(E^n_i - E^n_{i-1})
 *                  - c_kp_rho_V*E^n_i]
 *        + (dt/2)*c_kp_rho_V*[S_emit(T_iter) + S_emit(T^n)]
 *
 * BDF2 (second-order, A-stable, no spurious oscillations):
 *
 *   diag = (3/2)*V_i + dt*(a_e + a_w + c_kp_rho_V)
 *   sub  = -dt*a_w
 *   sup  = -dt*a_e
 *   rhs  = 2*V_i*E^n - (1/2)*V_i*E^{n-1} + dt*c_kp_rho_V*S_emit(T_iter)
 *
 * The per-step coefficients are exposed as a flat struct so the host
 * solvers can apply them inline without an additional virtual call per
 * cell. The "needsTwoPriorStates" accessor tells the host whether to
 * maintain an E^{n-1} ring buffer.
 *
 * References.
 *  - Strang, G. (1968), "On the construction and comparison of
 *    difference schemes", SIAM J. Num. Anal. 5(3), pp 506-517.
 *  - Larsen, E. W. (1988), "A grey transport acceleration method for
 *    time-dependent radiative transfer", J. Comp. Phys 78, pp 459-480
 *    (sec 3, time-centred linearisation of the matter coupling).
 *  - Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
 *    Hydrodynamics", Oxford University Press, sec 97 (operator splitting
 *    under implicit time stepping).
 *  - Hairer, E. and Wanner, G. (1996), "Solving Ordinary Differential
 *    Equations II: Stiff and Differential-Algebraic Problems", Springer.
 *    Sec V.1 (BDF formulae); sec IV.3 (Crank-Nicolson and theta-method).
 */

#ifndef NEAR_FIELD_DIFFUSION_TIME_INTEGRATOR_HPP
#define NEAR_FIELD_DIFFUSION_TIME_INTEGRATOR_HPP

#include <memory>

namespace FSRM {

/// Pass-11 axis-1c selector for the per-group diffusion solve.
/// BACKWARD_EULER is the pass-10 byte-identical default.
enum class DiffusionTimeIntegratorKind {
    BACKWARD_EULER,
    CRANK_NICOLSON,
    BDF2
};

/**
 * @brief Strategy interface for the per-group implicit time-integration
 *        of the radiation-diffusion solve. The host solvers query
 *        getCoefficients() once per time step and use the returned
 *        scalars to assemble the per-cell tridiagonal system.
 *
 * The Coefficients struct describes the canonical assembly:
 *
 *   diag = V_i*lhs_volume_factor
 *        + dt*(a_e + a_w + c_kp_rho_V)*lhs_implicit_factor
 *   sub  = -dt*a_w * lhs_implicit_factor
 *   sup  = -dt*a_e * lhs_implicit_factor
 *   rhs  = V_i*(rhs_volume_n_factor*E^n + rhs_volume_nm1_factor*E^{n-1})
 *        + dt*rhs_explicit_diffusion_factor*[a_e*(E^n_{i+1} - E^n_i)
 *                                            - a_w*(E^n_i - E^n_{i-1})
 *                                            - c_kp_rho_V*E^n_i]
 *        + dt*c_kp_rho_V*[rhs_source_implicit_factor*S(T_iter)
 *                         + rhs_source_explicit_factor*S(T^n)]
 *
 * Backward-Euler (1, 1, 1, 0, 0, 1, 0) recovers the pass-10 assembly
 * byte-for-byte. Crank-Nicolson (1, 1/2, 1, 0, 1/2, 1/2, 1/2) and
 * BDF2 (3/2, 1, 2, -1/2, 0, 1, 0) extend to second-order accuracy.
 */
class DiffusionTimeIntegrator
{
public:
    struct Coefficients
    {
        double lhs_volume_factor = 1.0;
        double lhs_implicit_factor = 1.0;
        double rhs_volume_n_factor = 1.0;
        double rhs_volume_nm1_factor = 0.0;
        double rhs_explicit_diffusion_factor = 0.0;
        double rhs_source_implicit_factor = 1.0;
        double rhs_source_explicit_factor = 0.0;
    };

    virtual ~DiffusionTimeIntegrator() = default;
    virtual Coefficients getCoefficients() const = 0;
    virtual bool needsTwoPriorStates() const = 0;
    virtual DiffusionTimeIntegratorKind kind() const = 0;
    virtual const char* name() const = 0;

    /// Factory: build the integrator selected by kind. Returns
    /// BackwardEulerIntegrator if kind is unrecognised so the host
    /// solvers degrade safely.
    static std::unique_ptr<DiffusionTimeIntegrator> create(
        DiffusionTimeIntegratorKind kind);
};

/// Pass-10 byte-identical (BACKWARD_EULER) integrator.
class BackwardEulerIntegrator final : public DiffusionTimeIntegrator
{
public:
    Coefficients getCoefficients() const override
    {
        Coefficients c;
        c.lhs_volume_factor = 1.0;
        c.lhs_implicit_factor = 1.0;
        c.rhs_volume_n_factor = 1.0;
        c.rhs_volume_nm1_factor = 0.0;
        c.rhs_explicit_diffusion_factor = 0.0;
        c.rhs_source_implicit_factor = 1.0;
        c.rhs_source_explicit_factor = 0.0;
        return c;
    }
    bool needsTwoPriorStates() const override { return false; }
    DiffusionTimeIntegratorKind kind() const override
    {
        return DiffusionTimeIntegratorKind::BACKWARD_EULER;
    }
    const char* name() const override { return "BACKWARD_EULER"; }
};

/// Crank-Nicolson (theta-method with theta = 1/2). Second-order
/// accurate, A-stable. Documented oscillation failure mode on stiff
/// initial conditions; gate the OscillationStability test against this.
class CrankNicolsonIntegrator final : public DiffusionTimeIntegrator
{
public:
    Coefficients getCoefficients() const override
    {
        Coefficients c;
        c.lhs_volume_factor = 1.0;
        c.lhs_implicit_factor = 0.5;
        c.rhs_volume_n_factor = 1.0;
        c.rhs_volume_nm1_factor = 0.0;
        c.rhs_explicit_diffusion_factor = 0.5;
        c.rhs_source_implicit_factor = 0.5;
        c.rhs_source_explicit_factor = 0.5;
        return c;
    }
    bool needsTwoPriorStates() const override { return false; }
    DiffusionTimeIntegratorKind kind() const override
    {
        return DiffusionTimeIntegratorKind::CRANK_NICOLSON;
    }
    const char* name() const override { return "CRANK_NICOLSON"; }
};

/// BDF2 second-order multistep. Requires E^{n-1}; the host solver
/// bootstraps the first step with backward-Euler and stores the prior
/// state for subsequent steps.
class BDF2Integrator final : public DiffusionTimeIntegrator
{
public:
    Coefficients getCoefficients() const override
    {
        Coefficients c;
        c.lhs_volume_factor = 1.5;
        c.lhs_implicit_factor = 1.0;
        c.rhs_volume_n_factor = 2.0;
        c.rhs_volume_nm1_factor = -0.5;
        c.rhs_explicit_diffusion_factor = 0.0;
        c.rhs_source_implicit_factor = 1.0;
        c.rhs_source_explicit_factor = 0.0;
        return c;
    }
    bool needsTwoPriorStates() const override { return true; }
    DiffusionTimeIntegratorKind kind() const override
    {
        return DiffusionTimeIntegratorKind::BDF2;
    }
    const char* name() const override { return "BDF2"; }
};

inline std::unique_ptr<DiffusionTimeIntegrator>
DiffusionTimeIntegrator::create(DiffusionTimeIntegratorKind kind)
{
    switch (kind) {
    case DiffusionTimeIntegratorKind::CRANK_NICOLSON:
        return std::unique_ptr<DiffusionTimeIntegrator>(
            new CrankNicolsonIntegrator);
    case DiffusionTimeIntegratorKind::BDF2:
        return std::unique_ptr<DiffusionTimeIntegrator>(new BDF2Integrator);
    case DiffusionTimeIntegratorKind::BACKWARD_EULER:
    default:
        return std::unique_ptr<DiffusionTimeIntegrator>(
            new BackwardEulerIntegrator);
    }
}

}  // namespace FSRM

#endif  // NEAR_FIELD_DIFFUSION_TIME_INTEGRATOR_HPP
