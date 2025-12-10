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
 * - SUPG/GLS stabilization for advection-dominated flows
 * - Higher-order time integration schemes
 * - Polynomial adaptivity (p-adaptivity)
 * - Physics-based preconditioning
 * 
 * @note For Discontinuous Galerkin methods, use DiscontinuousGalerkin.hpp
 *       which provides a comprehensive DG implementation with ADER time
 *       integration and local time stepping.
 * 
 * @note For Adaptive Mesh Refinement (h-adaptivity), use AdaptiveMeshRefinement.hpp
 *       which provides DMPlex/DMForest-based AMR with multiple error indicators.
 * 
 * Integrates with PETSc infrastructure in FSRM.hpp.
 */

#include "FSRM.hpp"
#include "PhysicsKernel.hpp"
#include "DiscontinuousGalerkin.hpp"      // For DG methods - use existing implementation
#include "AdaptiveMeshRefinement.hpp"     // For AMR - use existing implementation
#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <map>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// SUPG/GLS Stabilization (NEW - not in existing codebase)
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
     */
    void addGLSStabilization(double L_adjoint_v,
                            double residual_strong,
                            double tau,
                            double& stabilization_term) const;
    
    /**
     * @brief Calculate shock-capturing diffusivity
     */
    double shockCapturingDiffusivity(const std::array<double, 3>& grad_u,
                                     double u, double h,
                                     double residual) const;
    
    /**
     * @brief Calculate PSPG term for pressure stabilization
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
// Higher-Order Time Integration (NEW - extends basic time stepping)
// =============================================================================

/**
 * @brief Time integration scheme
 * 
 * @note For ADER time integration with DG, use ADERTimeIntegrator from
 *       DiscontinuousGalerkin.hpp which provides space-time DG methods.
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
    
    // Special purpose (for structural dynamics)
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
        double error_tolerance = 1e-3;
        double safety_factor = 0.9;
        double min_dt_factor = 0.1;
        double max_dt_factor = 5.0;
    };
    
    AdvancedTimeIntegration();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    int getOrder() const;
    bool isAStable() const;
    bool isLStable() const;
    
    /**
     * @brief Get Butcher tableau for RK methods
     */
    void getButcherTableau(std::vector<std::vector<double>>& A,
                          std::vector<double>& b,
                          std::vector<double>& c) const;
    
    /**
     * @brief Initialize for Newmark method
     */
    void initializeNewmark(const std::vector<double>& u0,
                          const std::vector<double>& v0,
                          const std::function<std::vector<double>(
                              const std::vector<double>&,
                              const std::vector<double>&)>& compute_rhs,
                          std::vector<double>& a0) const;
    
    /**
     * @brief Newmark predictor step
     */
    void newmarkPredict(const std::vector<double>& u_n,
                       const std::vector<double>& v_n,
                       const std::vector<double>& a_n,
                       double dt,
                       std::vector<double>& u_pred,
                       std::vector<double>& v_pred) const;
    
    /**
     * @brief Newmark corrector step
     */
    void newmarkCorrect(const std::vector<double>& u_pred,
                       const std::vector<double>& v_pred,
                       const std::vector<double>& a_new,
                       double dt,
                       std::vector<double>& u_new,
                       std::vector<double>& v_new) const;
    
    /**
     * @brief Get effective mass for Newmark
     */
    void effectiveNewmarkMatrix(double dt, double mass, double damping,
                               double stiffness, double& M_eff) const;
    
    /**
     * @brief Generalized-α coefficients from ρ_∞
     */
    void computeGeneralizedAlphaCoeffs();
    
    /**
     * @brief Estimate local truncation error
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
     */
    void bdfCoefficients(int order, std::vector<double>& alpha,
                        double& beta) const;
    
private:
    Parameters params_;
};


// =============================================================================
// p-Adaptivity (NEW - complements h-adaptivity in AdaptiveMeshRefinement.hpp)
// =============================================================================

/**
 * @brief p-adaptivity manager
 * 
 * Adapts polynomial order locally based on solution smoothness.
 * Higher order beneficial for smooth regions, lower order
 * for regions with discontinuities.
 * 
 * @note For h-adaptivity, use AdaptiveMeshRefinement class from
 *       AdaptiveMeshRefinement.hpp which provides comprehensive mesh refinement.
 */
class PAdaptivityManager {
public:
    struct Parameters {
        int min_order = 1;
        int max_order = 5;
        int default_order = 2;
        
        // Smoothness indicators
        double smoothness_threshold_inc = 0.8;
        double smoothness_threshold_dec = 0.2;
        
        // Error-based
        double p_refinement_threshold = 0.5;
    };
    
    PAdaptivityManager();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate smoothness indicator
     */
    double smoothnessIndicator(const std::vector<double>& legendre_coeffs) const;
    
    /**
     * @brief Decide whether to increase/decrease polynomial order
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
// Physics-Based Preconditioning (NEW - for coupled multi-physics)
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
        std::vector<int> block_sizes;
        std::vector<std::string> block_names;
        
        // Field split options
        std::string split_type = "additive";
        std::vector<std::string> block_solvers;
        
        // Schur complement
        std::string schur_precond = "selfp";
        double schur_scale = 1.0;
        
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
     */
    void applyBlockJacobi(const std::vector<std::vector<double>>& blocks,
                         const std::vector<double>& rhs,
                         std::vector<double>& solution) const;
    
    /**
     * @brief Apply block Gauss-Seidel
     */
    void applyBlockGaussSeidel(const std::vector<std::vector<double>>& diag_blocks,
                              const std::vector<std::vector<double>>& off_diag,
                              const std::vector<double>& rhs,
                              std::vector<double>& solution) const;
    
    /**
     * @brief Apply Schur complement method
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
 * 
 * @note For DG configuration, use DiscontinuousGalerkin class directly.
 * @note For AMR configuration, use AdaptiveMeshRefinement class directly.
 */
struct HighFidelityNumericsConfig {
    // Master enable
    bool enable_high_fidelity_numerics = false;
    
    // SUPG/GLS stabilization
    bool enable_supg = false;
    SUPGStabilization::Parameters supg_params;
    
    // Time integration
    bool enable_advanced_time = false;
    AdvancedTimeIntegration::Parameters time_params;
    
    // p-adaptivity
    bool enable_p_adapt = false;
    PAdaptivityManager::Parameters p_adapt_params;
    
    // Preconditioning
    bool enable_physics_precond = false;
    PhysicsBasedPreconditioner::Parameters precond_params;
    
    // References to existing implementations (for convenience)
    // DG: Use createDGSolver() from DiscontinuousGalerkin.hpp
    // AMR: Use AdaptiveMeshRefinement class from AdaptiveMeshRefinement.hpp
    
    void parseConfig(const std::map<std::string, std::string>& config);
    bool validate(std::string& error_msg) const;
};

/**
 * @brief Factory functions
 */
std::unique_ptr<SUPGStabilization> createSUPGStabilization(
    const std::map<std::string, std::string>& config);

std::unique_ptr<AdvancedTimeIntegration> createTimeIntegration(
    const std::map<std::string, std::string>& config);

// Note: For DG, use createDGSolver() from DiscontinuousGalerkin.hpp
// Note: For AMR, use AdaptiveMeshRefinement class from AdaptiveMeshRefinement.hpp

} // namespace HighFidelity
} // namespace FSRM

#endif // HIGH_FIDELITY_NUMERICS_HPP
