#ifndef HIGH_FIDELITY_NUMERICS_HPP
#define HIGH_FIDELITY_NUMERICS_HPP

/**
 * @file HighFidelityNumerics.hpp
 * @brief Advanced numerical methods for high-fidelity simulations
 * 
 * This file provides sophisticated numerical techniques that extend the basic
 * FE/FV methods used in the core solver. These are optional high-fidelity
 * extensions for complex simulations requiring enhanced accuracy or stability.
 * 
 * Features:
 * - Discontinuous Galerkin (DG) stabilization
 * - SUPG/GLS stabilization for advection-dominated flows
 * - Higher-order time integration schemes
 * - Adaptive mesh refinement (h-adaptivity)
 * - Polynomial adaptivity (p-adaptivity)
 * - hp-adaptivity
 * - Space-time formulations
 * - Iterative coupling strategies
 * - Newton-Krylov solvers with physics-based preconditioners
 * 
 * Integrates with PETSc infrastructure in FSRM.hpp.
 * 
 * @note These are advanced numerical options - basic simulations should use
 *       standard methods.
 */

#include "FSRM.hpp"
#include "PhysicsKernel.hpp"
#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <map>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Discontinuous Galerkin Methods
// =============================================================================

/**
 * @brief DG numerical flux type
 */
enum class DGFluxType {
    CENTRAL,             ///< Central (unstable for hyperbolic)
    UPWIND,              ///< Full upwind
    LAX_FRIEDRICHS,      ///< Local Lax-Friedrichs
    ROE,                 ///< Roe approximate Riemann solver
    HLLC,                ///< HLLC flux for Euler equations
    INTERIOR_PENALTY,    ///< IP/SIPG for diffusion
    NIPG,                ///< Non-symmetric interior penalty
    BASSI_REBAY          ///< BR1/BR2 lifting operators
};

/**
 * @brief DG penalty formulation
 */
enum class DGPenaltyType {
    SIPG,                ///< Symmetric Interior Penalty Galerkin
    NIPG,                ///< Non-symmetric IPG
    IIPG,                ///< Incomplete IPG
    BASSI_REBAY_1,       ///< BR1
    BASSI_REBAY_2,       ///< BR2 (compact stencil)
    LDG                  ///< Local DG
};

/**
 * @brief Discontinuous Galerkin stabilization
 * 
 * Provides DG formulations for:
 * - Advection (transport equations)
 * - Diffusion (with interior penalty)
 * - Advection-diffusion
 * - Wave propagation
 * 
 * Key features:
 * - Local conservation
 * - Handles discontinuities (shocks)
 * - High-order accuracy
 * - Flexible hp-refinement
 */
class DGStabilization {
public:
    struct Parameters {
        DGFluxType flux_type = DGFluxType::UPWIND;
        DGPenaltyType penalty_type = DGPenaltyType::SIPG;
        
        // Penalty parameter
        double sigma = 10.0;              ///< Penalty coefficient
        bool auto_penalty = true;          ///< Compute from element size
        double penalty_scaling = 1.0;      ///< Scaling factor
        
        // Polynomial degree
        int polynomial_order = 1;          ///< p = 1, 2, 3, ...
        
        // Limiting (for shocks)
        bool use_limiter = false;
        std::string limiter_type = "minmod"; ///< minmod, tvb, weno
        double tvb_constant = 0.0;         ///< TVB limiter constant
        
        // Stabilization for advection
        double upwind_fraction = 1.0;      ///< 0 = central, 1 = full upwind
    };
    
    DGStabilization();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate penalty parameter
     * 
     * σ = C · (p+1)² / h
     * 
     * where p is polynomial order, h is element size.
     */
    double calculatePenalty(double h, int p) const;
    
    /**
     * @brief Calculate numerical flux for advection
     * 
     * Upwind: F = a·u⁺ if a·n > 0, else a·u⁻
     * 
     * @param u_plus Value on + side of face
     * @param u_minus Value on - side of face
     * @param velocity Advection velocity
     * @param normal Face normal (pointing from - to +)
     * @return Numerical flux
     */
    double advectionFlux(double u_plus, double u_minus,
                        const std::array<double, 3>& velocity,
                        const std::array<double, 3>& normal) const;
    
    /**
     * @brief Calculate numerical flux for diffusion (SIPG)
     * 
     * {{∇u}}·n - σ/h·[[u]]
     * 
     * where {{·}} = average, [[·]] = jump
     */
    double diffusionFlux(double u_plus, double u_minus,
                        const std::array<double, 3>& grad_u_plus,
                        const std::array<double, 3>& grad_u_minus,
                        double diffusivity, double h,
                        const std::array<double, 3>& normal) const;
    
    /**
     * @brief Calculate Lax-Friedrichs flux
     * 
     * F = 0.5·(F⁺ + F⁻) - 0.5·λ·(u⁺ - u⁻)
     * 
     * where λ is max eigenvalue.
     */
    double laxFriedrichsFlux(double u_plus, double u_minus,
                            double flux_plus, double flux_minus,
                            double max_eigenvalue) const;
    
    /**
     * @brief Apply slope limiter
     * 
     * Limits polynomial coefficients to prevent oscillations.
     */
    void applyLimiter(std::vector<double>& solution,
                     const std::vector<std::array<double, 3>>& gradients,
                     const std::vector<double>& cell_sizes) const;
    
    /**
     * @brief Minmod function
     * 
     * minmod(a,b) = sign(a)·min(|a|,|b|) if sign(a)=sign(b), else 0
     */
    double minmod(double a, double b) const;
    double minmod(double a, double b, double c) const;
    
    /**
     * @brief TVB modified minmod
     */
    double tvbMinmod(double a, double b, double c, double h) const;
    
    /**
     * @brief Calculate face integrals for DG assembly
     * 
     * Returns contributions to elements sharing a face.
     */
    void faceIntegral(int face_id,
                     const std::vector<double>& u_plus,
                     const std::vector<double>& u_minus,
                     const std::array<double, 3>& normal,
                     double face_area, double h_plus, double h_minus,
                     std::vector<double>& residual_plus,
                     std::vector<double>& residual_minus) const;
    
private:
    Parameters params_;
};


// =============================================================================
// SUPG/GLS Stabilization
// =============================================================================

/**
 * @brief Stabilization method type
 */
enum class StabilizationMethod {
    NONE,
    SUPG,                ///< Streamline-Upwind Petrov-Galerkin
    GLS,                 ///< Galerkin Least-Squares
    PSPG,                ///< Pressure-Stabilized Petrov-Galerkin
    VMS,                 ///< Variational Multi-Scale
    SUBGRID_VISCOSITY,   ///< Artificial viscosity
    BUBBLE_ENRICHMENT    ///< Bubble function enrichment
};

/**
 * @brief Stabilization parameter formula
 */
enum class TauFormula {
    SHAKIB,              ///< Shakib et al. (common choice)
    TEZDUYAR,            ///< Tezduyar (for incompressible flow)
    CODINA,              ///< Codina (projection-based)
    FRANCA_VALENTIN,     ///< Franca & Valentin
    OPTIMAL              ///< Element-optimal τ
};

/**
 * @brief SUPG/GLS stabilization for advection-dominated problems
 * 
 * Standard Galerkin FEM produces oscillations when Pe > 1.
 * Stabilization adds artificial diffusion along streamlines.
 * 
 * SUPG: (v + τ·a·∇v, L(u)) = (v + τ·a·∇v, f)
 * 
 * where τ is the stabilization parameter, a is advection velocity,
 * and L is the differential operator.
 */
class SUPGStabilization {
public:
    struct Parameters {
        StabilizationMethod method = StabilizationMethod::SUPG;
        TauFormula tau_formula = TauFormula::SHAKIB;
        
        // Stabilization parameters
        double tau_multiplier = 1.0;       ///< Scale factor for τ
        double crosswind_diffusion = 0.0;  ///< Crosswind artificial diffusion
        
        // For transient problems
        bool include_transient = true;     ///< Include ∂u/∂t in τ
        double time_scale = 1.0;           ///< dt for transient τ
        
        // For reaction terms
        bool include_reaction = true;      ///< Include reaction in τ
        
        // Shock capturing
        bool shock_capturing = false;
        double shock_capture_coeff = 0.5;  ///< DC coefficient
    };
    
    SUPGStabilization();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate stabilization parameter τ
     * 
     * Shakib formula:
     * τ = (4/dt² + (2|a|/h)² + 9(4κ/h²)²)^(-1/2)
     * 
     * @param velocity Advection velocity
     * @param diffusivity Diffusion coefficient
     * @param reaction Reaction coefficient
     * @param h Element size
     * @param dt Time step (0 for steady-state)
     * @return Stabilization parameter τ
     */
    double calculateTau(const std::array<double, 3>& velocity,
                       double diffusivity, double reaction,
                       double h, double dt = 0.0) const;
    
    /**
     * @brief Calculate element Peclet number
     * 
     * Pe = |a|·h / (2·κ)
     */
    double pecletNumber(double velocity_magnitude, double h,
                       double diffusivity) const;
    
    /**
     * @brief Calculate optimal diffusivity for Pe >> 1
     * 
     * κ_art = 0.5·|a|·h·(coth(Pe) - 1/Pe)
     */
    double optimalDiffusivity(double velocity_magnitude, double h,
                             double diffusivity) const;
    
    /**
     * @brief Add SUPG stabilization term to residual
     * 
     * R_stab = τ·(a·∇v)·(L(u) - f)
     */
    void addStabilization(const std::array<double, 3>& velocity,
                         const std::array<double, 3>& test_gradient,
                         double residual_strong,
                         double tau,
                         double& stabilization_term) const;
    
    /**
     * @brief Add GLS stabilization term
     * 
     * R_GLS = τ·L*(v)·(L(u) - f)
     * 
     * where L* is adjoint operator.
     */
    void addGLSStabilization(double L_adjoint_v,
                            double residual_strong,
                            double tau,
                            double& stabilization_term) const;
    
    /**
     * @brief Calculate shock-capturing diffusivity
     * 
     * Adds isotropic diffusion in regions of high gradients.
     */
    double shockCapturingDiffusivity(const std::array<double, 3>& grad_u,
                                     double u, double h,
                                     double residual) const;
    
    /**
     * @brief Calculate PSPG term for pressure stabilization
     * 
     * For incompressible flow to avoid pressure oscillations.
     */
    double pspgTerm(const std::array<double, 3>& grad_p,
                   const std::array<double, 3>& grad_q,
                   double tau) const;
    
private:
    Parameters params_;
    
    double tau_Shakib(const std::array<double, 3>& velocity,
                     double diffusivity, double reaction,
                     double h, double dt) const;
    double tau_Tezduyar(const std::array<double, 3>& velocity,
                       double diffusivity, double h, double dt) const;
    double tau_Codina(const std::array<double, 3>& velocity,
                     double diffusivity, double reaction,
                     double h, double dt) const;
};


// =============================================================================
// Higher-Order Time Integration
// =============================================================================

/**
 * @brief Time integration scheme
 */
enum class TimeScheme {
    // Single-step methods
    FORWARD_EULER,       ///< Explicit, 1st order
    BACKWARD_EULER,      ///< Implicit, 1st order
    CRANK_NICOLSON,      ///< Implicit, 2nd order
    THETA_METHOD,        ///< Generalized θ-method
    
    // Runge-Kutta
    RK2,                 ///< Explicit RK2
    RK4,                 ///< Classic RK4
    SSPRK3,              ///< Strong stability preserving RK3
    DIRK2,               ///< Diagonally implicit RK2
    DIRK3,               ///< DIRK3 (A-stable)
    ESDIRK4,             ///< Explicit first stage DIRK4
    
    // Multi-step
    BDF2,                ///< 2nd order BDF
    BDF3,                ///< 3rd order BDF
    BDF4,                ///< 4th order BDF
    ADAMS_BASHFORTH_2,   ///< AB2 (explicit)
    ADAMS_MOULTON_2,     ///< AM2 (implicit)
    
    // Special purpose
    NEWMARK_BETA,        ///< Newmark-β for dynamics
    HHT_ALPHA,           ///< HHT-α (numerical damping)
    GENERALIZED_ALPHA,   ///< Generalized-α method
    BATHE                ///< Bathe composite scheme
};

/**
 * @brief Higher-order time integration schemes
 * 
 * Provides various time integration methods beyond the basic
 * Backward Euler used in the core solver.
 * 
 * Key considerations:
 * - Stability (A-stable, L-stable)
 * - Order of accuracy
 * - Numerical dissipation (high-frequency damping)
 * - Energy conservation
 */
class AdvancedTimeIntegration {
public:
    struct Parameters {
        TimeScheme scheme = TimeScheme::BACKWARD_EULER;
        
        // θ-method parameters
        double theta = 1.0;                ///< 0.5 = CN, 1.0 = BE
        
        // Newmark parameters
        double newmark_beta = 0.25;        ///< Typically 0.25 or 0.3025
        double newmark_gamma = 0.5;        ///< Typically 0.5 or 0.6
        
        // Generalized-α parameters
        double rho_infinity = 0.5;         ///< Spectral radius at infinity
        double alpha_m = 0.0;              ///< Computed from rho_inf
        double alpha_f = 0.0;              ///< Computed from rho_inf
        
        // Adaptive time stepping
        bool adaptive = false;
        double error_tolerance = 1e-3;     ///< Local truncation error
        double safety_factor = 0.9;        ///< dt scaling factor
        double min_dt_factor = 0.1;        ///< min(dt_new/dt_old)
        double max_dt_factor = 5.0;        ///< max(dt_new/dt_old)
    };
    
    AdvancedTimeIntegration();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Get order of accuracy
     */
    int getOrder() const;
    
    /**
     * @brief Check if scheme is A-stable
     */
    bool isAStable() const;
    
    /**
     * @brief Check if scheme is L-stable
     */
    bool isLStable() const;
    
    /**
     * @brief Get Butcher tableau for RK methods
     * 
     * Returns (A, b, c) where:
     * - A is the RK matrix (s × s)
     * - b is weights (s)
     * - c is nodes (s)
     */
    void getButcherTableau(std::vector<std::vector<double>>& A,
                          std::vector<double>& b,
                          std::vector<double>& c) const;
    
    /**
     * @brief Initialize for Newmark method
     * 
     * Computes initial accelerations from equilibrium.
     */
    void initializeNewmark(const std::vector<double>& u0,
                          const std::vector<double>& v0,
                          const std::function<std::vector<double>(
                              const std::vector<double>&,
                              const std::vector<double>&)>& compute_rhs,
                          std::vector<double>& a0) const;
    
    /**
     * @brief Newmark predictor step
     * 
     * ũ = u_n + dt·v_n + dt²·(0.5-β)·a_n
     * ṽ = v_n + dt·(1-γ)·a_n
     */
    void newmarkPredict(const std::vector<double>& u_n,
                       const std::vector<double>& v_n,
                       const std::vector<double>& a_n,
                       double dt,
                       std::vector<double>& u_pred,
                       std::vector<double>& v_pred) const;
    
    /**
     * @brief Newmark corrector step
     * 
     * u_{n+1} = ũ + β·dt²·a_{n+1}
     * v_{n+1} = ṽ + γ·dt·a_{n+1}
     */
    void newmarkCorrect(const std::vector<double>& u_pred,
                       const std::vector<double>& v_pred,
                       const std::vector<double>& a_new,
                       double dt,
                       std::vector<double>& u_new,
                       std::vector<double>& v_new) const;
    
    /**
     * @brief Get effective mass for Newmark
     * 
     * M_eff = M + γ·dt·C + β·dt²·K
     */
    void effectiveNewmarkMatrix(double dt, double mass, double damping,
                               double stiffness, double& M_eff) const;
    
    /**
     * @brief Generalized-α coefficients from ρ_∞
     * 
     * α_m = (2·ρ_∞ - 1) / (ρ_∞ + 1)
     * α_f = ρ_∞ / (ρ_∞ + 1)
     * γ = 0.5 - α_m + α_f
     * β = 0.25·(1 - α_m + α_f)²
     */
    void computeGeneralizedAlphaCoeffs();
    
    /**
     * @brief Estimate local truncation error
     * 
     * For adaptive time stepping.
     */
    double estimateError(const std::vector<double>& u_high,
                        const std::vector<double>& u_low) const;
    
    /**
     * @brief Compute new time step from error estimate
     */
    double computeNewTimestep(double dt_current, double error,
                             double tolerance) const;
    
    /**
     * @brief BDF coefficients
     * 
     * Returns (α_k, β) where:
     * Σ α_k·u_{n-k} = β·dt·f(u_{n+1})
     */
    void bdfCoefficients(int order, std::vector<double>& alpha,
                        double& beta) const;
    
private:
    Parameters params_;
};


// =============================================================================
// Adaptive Mesh Refinement (AMR)
// =============================================================================

/**
 * @brief Refinement indicator type
 */
enum class RefinementIndicator {
    GRADIENT,            ///< Solution gradient magnitude
    RESIDUAL,            ///< Element residual
    RECOVERY,            ///< ZZ error estimator
    HESSIAN,             ///< Second derivative
    JUMPS,               ///< Inter-element jumps
    GOAL_ORIENTED,       ///< Dual-weighted residual
    PHYSICS_BASED        ///< Problem-specific indicators
};

/**
 * @brief Refinement strategy
 */
enum class RefinementStrategy {
    FIXED_FRACTION,      ///< Refine/coarsen fixed % of elements
    FIXED_NUMBER,        ///< Refine/coarsen fixed number
    THRESHOLD,           ///< Refine if error > threshold
    EQUILIBRATION,       ///< Equilibrate error across elements
    KELLY                ///< Kelly error estimator
};

/**
 * @brief h-adaptivity manager
 * 
 * Handles adaptive mesh refinement based on error indicators.
 * 
 * Supports:
 * - Isotropic refinement (subdivide all directions)
 * - Anisotropic refinement (directional)
 * - Hanging nodes
 * - Conforming refinement
 */
class AdaptiveMeshRefinementManager {
public:
    struct Parameters {
        RefinementIndicator indicator = RefinementIndicator::GRADIENT;
        RefinementStrategy strategy = RefinementStrategy::FIXED_FRACTION;
        
        // Refinement fractions
        double refine_fraction = 0.3;      ///< Top 30% refined
        double coarsen_fraction = 0.1;     ///< Bottom 10% coarsened
        
        // Threshold-based
        double refine_threshold = 1.0;     ///< Error threshold for refinement
        double coarsen_threshold = 0.1;    ///< Error threshold for coarsening
        
        // Limits
        int max_refinement_level = 5;      ///< Maximum refinement depth
        int min_refinement_level = 0;      ///< Minimum refinement depth
        
        // Mesh quality
        double min_element_size = 1e-6;    ///< Minimum h
        double max_aspect_ratio = 10.0;    ///< Maximum aspect ratio
        
        // Conforming mesh
        bool ensure_conforming = true;     ///< Eliminate hanging nodes
        int max_level_difference = 1;      ///< Max level diff between neighbors
    };
    
    AdaptiveMeshRefinementManager();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate gradient-based indicator
     * 
     * η_K = h_K · ||∇u||_K
     */
    double gradientIndicator(const std::array<double, 3>& gradient,
                            double h) const;
    
    /**
     * @brief Calculate residual-based indicator
     * 
     * η_K² = h_K²·||R||²_K + h_K·Σ||J||²_∂K
     * 
     * where R is element residual and J is jump residual.
     */
    double residualIndicator(double element_residual,
                            const std::vector<double>& face_jumps,
                            double h) const;
    
    /**
     * @brief Zienkiewicz-Zhu recovery-based indicator
     * 
     * Computes error from difference between FE solution
     * and smoothed (recovered) solution.
     */
    double recoveryIndicator(const std::array<double, 3>& raw_gradient,
                            const std::array<double, 3>& recovered_gradient) const;
    
    /**
     * @brief Goal-oriented (dual-weighted) indicator
     * 
     * η_K = |R(u_h, z_h)| where z_h solves adjoint problem.
     */
    double goalOrientedIndicator(double primal_residual,
                                double dual_weight) const;
    
    /**
     * @brief Mark elements for refinement
     * 
     * @param indicators Error indicators per element
     * @return Refinement flags (-1 = coarsen, 0 = keep, +1 = refine)
     */
    std::vector<int> markElements(const std::vector<double>& indicators) const;
    
    /**
     * @brief Ensure mesh conformity
     * 
     * Propagates refinement to avoid hanging nodes.
     */
    void enforceConformity(std::vector<int>& refinement_flags,
                          const std::vector<std::vector<int>>& neighbors) const;
    
    /**
     * @brief Compute global error estimate
     * 
     * ||e||² ≈ Σ η_K²
     */
    double globalErrorEstimate(const std::vector<double>& indicators) const;
    
    /**
     * @brief Compute effectivity index
     * 
     * θ = estimated_error / true_error
     * 
     * Ideally close to 1.0.
     */
    double effectivityIndex(double estimated_error, double true_error) const;
    
private:
    Parameters params_;
    
    void fixedFractionStrategy(const std::vector<double>& indicators,
                              std::vector<int>& flags) const;
    void thresholdStrategy(const std::vector<double>& indicators,
                          std::vector<int>& flags) const;
};


// =============================================================================
// p-Adaptivity
// =============================================================================

/**
 * @brief p-adaptivity manager
 * 
 * Adapts polynomial order locally based on solution smoothness.
 * 
 * Higher order beneficial for smooth regions, lower order
 * for regions with discontinuities.
 */
class PAdaptivityManager {
public:
    struct Parameters {
        int min_order = 1;                 ///< Minimum polynomial order
        int max_order = 5;                 ///< Maximum polynomial order
        int default_order = 2;             ///< Default starting order
        
        // Smoothness indicators
        double smoothness_threshold_inc = 0.8;  ///< Increase p if smoother
        double smoothness_threshold_dec = 0.2;  ///< Decrease p if rougher
        
        // Error-based
        double p_refinement_threshold = 0.5;    ///< Error threshold to increase p
    };
    
    PAdaptivityManager();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate smoothness indicator
     * 
     * Based on decay of Legendre coefficients.
     * Fast decay = smooth = high p beneficial.
     */
    double smoothnessIndicator(const std::vector<double>& legendre_coeffs) const;
    
    /**
     * @brief Decide whether to increase/decrease polynomial order
     * 
     * @return -1 = decrease, 0 = keep, +1 = increase
     */
    int decideOrderChange(double smoothness, double error_indicator,
                         int current_order) const;
    
    /**
     * @brief Project solution from order p to order p+1
     */
    void increaseOrder(const std::vector<double>& coeffs_p,
                      std::vector<double>& coeffs_p1,
                      int p) const;
    
    /**
     * @brief Project solution from order p to order p-1
     */
    void decreaseOrder(const std::vector<double>& coeffs_p,
                      std::vector<double>& coeffs_p1,
                      int p) const;
    
private:
    Parameters params_;
};


// =============================================================================
// Physics-Based Preconditioning
// =============================================================================

/**
 * @brief Preconditioner type for coupled systems
 */
enum class CoupledPreconditioner {
    BLOCK_JACOBI,        ///< Block diagonal
    BLOCK_GAUSS_SEIDEL,  ///< Block triangular
    SCHUR_COMPLEMENT,    ///< Schur complement reduction
    FIELD_SPLIT,         ///< PETSc field split
    PHYSICS_BASED,       ///< Problem-specific
    MULTIGRID            ///< Algebraic multigrid
};

/**
 * @brief Physics-based preconditioning for coupled systems
 * 
 * Exploits problem structure for efficient preconditioning:
 * - Block structure for multi-physics
 * - Approximate Schur complements
 * - Sequential solves following physics
 */
class PhysicsBasedPreconditioner {
public:
    struct Parameters {
        CoupledPreconditioner type = CoupledPreconditioner::FIELD_SPLIT;
        
        // Block structure
        std::vector<int> block_sizes;          ///< DOFs per block
        std::vector<std::string> block_names;  ///< Names for debugging
        
        // Field split options
        std::string split_type = "additive";   ///< additive, multiplicative
        std::vector<std::string> block_solvers; ///< Solver per block
        
        // Schur complement
        std::string schur_precond = "selfp";   ///< Schur approximation
        double schur_scale = 1.0;              ///< Scaling factor
        
        // Multigrid
        int mg_levels = 3;
        std::string mg_smoother = "sor";
        int mg_smooth_sweeps = 2;
    };
    
    PhysicsBasedPreconditioner();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Apply block Jacobi preconditioner
     * 
     * Solves each block independently.
     */
    void applyBlockJacobi(const std::vector<std::vector<double>>& blocks,
                         const std::vector<double>& rhs,
                         std::vector<double>& solution) const;
    
    /**
     * @brief Apply block Gauss-Seidel
     * 
     * Sequential block solve with coupling.
     */
    void applyBlockGaussSeidel(const std::vector<std::vector<double>>& diag_blocks,
                              const std::vector<std::vector<double>>& off_diag,
                              const std::vector<double>& rhs,
                              std::vector<double>& solution) const;
    
    /**
     * @brief Compute Schur complement
     * 
     * For [A B; C D]x = [f; g]:
     * S = D - C·A⁻¹·B
     * y = S⁻¹·(g - C·A⁻¹·f)
     * x = A⁻¹·(f - B·y)
     */
    void applySchurComplement(const std::vector<std::vector<double>>& A,
                             const std::vector<std::vector<double>>& B,
                             const std::vector<std::vector<double>>& C,
                             const std::vector<std::vector<double>>& D,
                             const std::vector<double>& f,
                             const std::vector<double>& g,
                             std::vector<double>& x,
                             std::vector<double>& y) const;
    
    /**
     * @brief Create PETSc PCFieldSplit preconditioner
     * 
     * Configures PETSc for block preconditioning.
     */
    void setupFieldSplit(PC pc, const std::vector<IS>& index_sets) const;
    
private:
    Parameters params_;
};


// =============================================================================
// High-Fidelity Numerics Configuration
// =============================================================================

/**
 * @brief Configuration for high-fidelity numerics
 */
struct HighFidelityNumericsConfig {
    // Master enable
    bool enable_high_fidelity_numerics = false;
    
    // DG stabilization
    bool enable_dg = false;
    DGStabilization::Parameters dg_params;
    
    // SUPG/GLS stabilization
    bool enable_supg = false;
    SUPGStabilization::Parameters supg_params;
    
    // Time integration
    bool enable_advanced_time = false;
    AdvancedTimeIntegration::Parameters time_params;
    
    // AMR
    bool enable_amr = false;
    AdaptiveMeshRefinementManager::Parameters amr_params;
    
    // p-adaptivity
    bool enable_p_adapt = false;
    PAdaptivityManager::Parameters p_adapt_params;
    
    // Preconditioning
    bool enable_physics_precond = false;
    PhysicsBasedPreconditioner::Parameters precond_params;
    
    // Parse from config
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Validate
    bool validate(std::string& error_msg) const;
};

/**
 * @brief Factory functions
 */
std::unique_ptr<DGStabilization> createDGStabilization(
    const std::map<std::string, std::string>& config);

std::unique_ptr<SUPGStabilization> createSUPGStabilization(
    const std::map<std::string, std::string>& config);

std::unique_ptr<AdvancedTimeIntegration> createTimeIntegration(
    const std::map<std::string, std::string>& config);

std::unique_ptr<AdaptiveMeshRefinementManager> createAMRManager(
    const std::map<std::string, std::string>& config);

} // namespace HighFidelity
} // namespace FSRM

#endif // HIGH_FIDELITY_NUMERICS_HPP
