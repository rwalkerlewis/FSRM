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
 * **PETSc Integration:**
 * - Time integration uses PETSc TS (supports BDF, RK, Generalized-alpha, etc.)
 * - Preconditioning uses PETSc PC (FieldSplit, Schur, Multigrid)
 * - Vector/Matrix operations use PETSc Vec/Mat
 * 
 * @note For Discontinuous Galerkin methods, use DiscontinuousGalerkin.hpp
 * @note For Adaptive Mesh Refinement (h-adaptivity), use AdaptiveMeshRefinement.hpp
 */

#include "FSRM.hpp"
#include "PhysicsKernel.hpp"
#include "DiscontinuousGalerkin.hpp"
#include "AdaptiveMeshRefinement.hpp"

// PETSc headers for numerical methods
#include <petscts.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <map>

namespace FSRM {
namespace HighFidelity {

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
 */
class SUPGStabilization {
public:
    struct Parameters {
        StabilizationMethod method = StabilizationMethod::SUPG;
        TauFormula tau_formula = TauFormula::SHAKIB;
        
        double tau_multiplier = 1.0;
        double crosswind_diffusion = 0.0;
        bool include_transient = true;
        double time_scale = 1.0;
        bool include_reaction = true;
        bool shock_capturing = false;
        double shock_capture_coeff = 0.5;
    };
    
    SUPGStabilization();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate stabilization parameter τ using PETSc-compatible scalars
     */
    PetscReal calculateTau(const PetscReal velocity[3],
                          PetscReal diffusivity, PetscReal reaction,
                          PetscReal h, PetscReal dt = 0.0) const;
    
    PetscReal pecletNumber(PetscReal velocity_magnitude, PetscReal h,
                          PetscReal diffusivity) const;
    
    PetscReal optimalDiffusivity(PetscReal velocity_magnitude, PetscReal h,
                                 PetscReal diffusivity) const;
    
    /**
     * @brief Add SUPG stabilization term (for use in PETSc residual functions)
     */
    void addStabilization(const PetscReal velocity[3],
                         const PetscReal test_gradient[3],
                         PetscReal residual_strong,
                         PetscReal tau,
                         PetscReal* stabilization_term) const;
    
    void addGLSStabilization(PetscReal L_adjoint_v,
                            PetscReal residual_strong,
                            PetscReal tau,
                            PetscReal* stabilization_term) const;
    
    PetscReal shockCapturingDiffusivity(const PetscReal grad_u[3],
                                        PetscReal u, PetscReal h,
                                        PetscReal residual) const;
    
    PetscReal pspgTerm(const PetscReal grad_p[3],
                      const PetscReal grad_q[3],
                      PetscReal tau) const;
    
private:
    Parameters params_;
    
    PetscReal tau_Shakib(const PetscReal velocity[3],
                        PetscReal diffusivity, PetscReal reaction,
                        PetscReal h, PetscReal dt) const;
    PetscReal tau_Tezduyar(const PetscReal velocity[3],
                          PetscReal diffusivity, PetscReal h, PetscReal dt) const;
    PetscReal tau_Codina(const PetscReal velocity[3],
                        PetscReal diffusivity, PetscReal reaction,
                        PetscReal h, PetscReal dt) const;
};


// =============================================================================
// PETSc-Based Time Integration
// =============================================================================

/**
 * @brief Time integration scheme (maps to PETSc TS types)
 * 
 * These map directly to PETSc's TS implementations for efficiency and reliability.
 */
enum class TimeScheme {
    // PETSc TSTHETA
    BACKWARD_EULER,      ///< TSTHETA with theta=1.0
    CRANK_NICOLSON,      ///< TSTHETA with theta=0.5
    THETA_METHOD,        ///< TSTHETA with custom theta
    
    // PETSc explicit RK (TSRK)
    RK1FE,               ///< Forward Euler (TSRK1FE)
    RK2A,                ///< 2nd order RK (TSRK2A)
    RK3,                 ///< 3rd order RK (TSRK3)
    RK4,                 ///< Classic RK4 (TSRK4)
    RK5F,                ///< Fehlberg 5th order (TSRK5F)
    RK5DP,               ///< Dormand-Prince 5th order (TSRK5DP)
    RK5BS,               ///< Bogacki-Shampine 5th order
    SSPRK3,              ///< SSP RK3 (TSSSP)
    
    // PETSc implicit RK (TSARKIMEX)
    ARKIMEX1,            ///< ARKIMEX 1st order
    ARKIMEX2,            ///< ARKIMEX 2nd order (L-stable)
    ARKIMEX3,            ///< ARKIMEX 3rd order
    ARKIMEX4,            ///< ARKIMEX 4th order
    ARKIMEX5,            ///< ARKIMEX 5th order
    
    // PETSc BDF (TSBDF)
    BDF1,                ///< BDF1 (Backward Euler)
    BDF2,                ///< BDF2
    BDF3,                ///< BDF3
    BDF4,                ///< BDF4
    BDF5,                ///< BDF5
    BDF6,                ///< BDF6
    
    // PETSc Rosenbrock-W (TSROSW)
    ROSW2M,              ///< 2nd order Rosenbrock-W
    ROSW2P,              ///< 2nd order Rosenbrock-W (another)
    ROSSW3P,             ///< 3rd order Rosenbrock-W
    ROSW4L,              ///< 4th order Rosenbrock-W (L-stable)
    
    // PETSc Generalized-alpha (TSALPHA/TSALPHA2)
    GENERALIZED_ALPHA,   ///< Generalized-α for first-order systems
    GENERALIZED_ALPHA2,  ///< Generalized-α for second-order (structural dynamics)
    
    // PETSc Newmark (TSNEWMARK) - specifically for structural dynamics
    NEWMARK_BETA         ///< Newmark-β method
};

/**
 * @brief PETSc-based time integration wrapper
 * 
 * This class wraps PETSc's TS (Time Stepper) for high-fidelity time integration.
 * All heavy computation is done by PETSc, not manual loops.
 * 
 * PETSc TS provides:
 * - Automatic time step adaptation
 * - Error estimation
 * - Event handling
 * - Adjoint sensitivity
 */
class PETScTimeIntegration {
public:
    struct Parameters {
        TimeScheme scheme = TimeScheme::BDF2;
        
        // θ-method parameters (for TSTHETA)
        PetscReal theta = 1.0;
        
        // Newmark parameters (for TSNEWMARK)
        PetscReal newmark_beta = 0.25;
        PetscReal newmark_gamma = 0.5;
        
        // Generalized-α parameters
        PetscReal rho_infinity = 0.5;
        
        // Adaptive time stepping (built into PETSc TS)
        PetscBool adaptive = PETSC_TRUE;
        PetscReal atol = 1e-6;           ///< Absolute tolerance
        PetscReal rtol = 1e-6;           ///< Relative tolerance
        PetscReal dt_min = 1e-10;
        PetscReal dt_max = 1e10;
        PetscInt max_steps = 10000;
        
        // For ARKIMEX
        PetscBool fully_implicit = PETSC_FALSE;
    };
    
    PETScTimeIntegration(MPI_Comm comm);
    ~PETScTimeIntegration();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Create and configure PETSc TS object
     */
    PetscErrorCode createTS(DM dm);
    
    /**
     * @brief Set the RHS function for du/dt = F(t, u)
     */
    PetscErrorCode setRHSFunction(TSRHSFunction func, void* ctx);
    
    /**
     * @brief Set the RHS Jacobian for implicit methods
     */
    PetscErrorCode setRHSJacobian(TSRHSJacobian func, Mat J, Mat P, void* ctx);
    
    /**
     * @brief Set the IFunction for implicit form F(t, u, u_t) = 0
     */
    PetscErrorCode setIFunction(TSIFunction func, Vec F, void* ctx);
    
    /**
     * @brief Set the IJacobian dF/du + shift*dF/du_t
     */
    PetscErrorCode setIJacobian(TSIJacobian func, Mat J, Mat P, void* ctx);
    
    /**
     * @brief For second-order systems (structural dynamics)
     */
    PetscErrorCode setI2Function(TSI2Function func, Vec F, void* ctx);
    PetscErrorCode setI2Jacobian(TSI2Jacobian func, Mat J, Mat P, void* ctx);
    
    /**
     * @brief Solve from t0 to tf
     */
    PetscErrorCode solve(Vec u, PetscReal t0, PetscReal tf);
    
    /**
     * @brief Single time step
     */
    PetscErrorCode step(Vec u);
    
    /**
     * @brief Get current time and time step
     */
    PetscErrorCode getTime(PetscReal* t) const;
    PetscErrorCode getTimeStep(PetscReal* dt) const;
    
    /**
     * @brief Set time step (for manual control)
     */
    PetscErrorCode setTimeStep(PetscReal dt);
    
    /**
     * @brief Get underlying PETSc TS object for advanced configuration
     */
    TS getTS() { return ts_; }
    
    /**
     * @brief Get order of accuracy
     */
    PetscInt getOrder() const;
    
    /**
     * @brief Check stability properties
     */
    PetscBool isAStable() const;
    PetscBool isLStable() const;
    
private:
    MPI_Comm comm_;
    TS ts_;
    Parameters params_;
    PetscBool ts_created_;
    
    PetscErrorCode configureTSType();
    PetscErrorCode configureAdaptivity();
};


// =============================================================================
// p-Adaptivity (complements h-adaptivity in AdaptiveMeshRefinement.hpp)
// =============================================================================

/**
 * @brief p-adaptivity manager using PETSc infrastructure
 */
class PAdaptivityManager {
public:
    struct Parameters {
        PetscInt min_order = 1;
        PetscInt max_order = 5;
        PetscInt default_order = 2;
        
        PetscReal smoothness_threshold_inc = 0.8;
        PetscReal smoothness_threshold_dec = 0.2;
        PetscReal p_refinement_threshold = 0.5;
    };
    
    PAdaptivityManager();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate smoothness indicator using PETSc Vec operations
     */
    PetscErrorCode smoothnessIndicator(Vec legendre_coeffs, PetscReal* indicator) const;
    
    /**
     * @brief Decide order change based on smoothness and error
     * @return -1 = decrease, 0 = keep, +1 = increase
     */
    PetscInt decideOrderChange(PetscReal smoothness, PetscReal error_indicator,
                               PetscInt current_order) const;
    
    /**
     * @brief Project solution to different polynomial order using PETSc
     */
    PetscErrorCode changeOrder(DM dm_old, Vec u_old, 
                               PetscInt new_order, DM* dm_new, Vec* u_new) const;
    
private:
    Parameters params_;
};


// =============================================================================
// PETSc-Based Physics Preconditioning
// =============================================================================

/**
 * @brief Preconditioner type (maps to PETSc PC types)
 */
enum class PreconditionerType {
    // Block methods
    BJACOBI,             ///< PCBJACOBI - Block Jacobi
    ASM,                 ///< PCASM - Additive Schwarz
    GASM,                ///< PCGASM - Generalized Additive Schwarz
    
    // Field split (for multi-physics)
    FIELDSPLIT_ADDITIVE,      ///< Additive field split
    FIELDSPLIT_MULTIPLICATIVE, ///< Multiplicative field split
    FIELDSPLIT_SYMMETRIC,     ///< Symmetric multiplicative
    FIELDSPLIT_SCHUR,         ///< Schur complement
    
    // Multigrid
    MG,                  ///< PCMG - Geometric multigrid
    GAMG,                ///< PCGAMG - Algebraic multigrid
    ML,                  ///< PCML - ML algebraic multigrid
    HYPRE_BOOMERAMG,     ///< PCHYPRE with BoomerAMG
    
    // Direct (for sub-blocks)
    LU,                  ///< LU factorization
    CHOLESKY,            ///< Cholesky (SPD systems)
    
    // Specialized
    LSC,                 ///< PCLSC - Least squares commutator (Stokes)
    SIMPLE,              ///< SIMPLE-type (incompressible flow)
    USER                 ///< User-defined
};

/**
 * @brief Schur complement factorization type (for FIELDSPLIT_SCHUR)
 */
enum class SchurFactType {
    DIAG,                ///< PC_FIELDSPLIT_SCHUR_FACT_DIAG
    LOWER,               ///< PC_FIELDSPLIT_SCHUR_FACT_LOWER
    UPPER,               ///< PC_FIELDSPLIT_SCHUR_FACT_UPPER
    FULL                 ///< PC_FIELDSPLIT_SCHUR_FACT_FULL
};

/**
 * @brief Schur complement preconditioning type
 */
enum class SchurPreType {
    SELF,                ///< PC_FIELDSPLIT_SCHUR_PRE_SELF
    SELFP,               ///< PC_FIELDSPLIT_SCHUR_PRE_SELFP
    A11,                 ///< PC_FIELDSPLIT_SCHUR_PRE_A11
    USER,                ///< PC_FIELDSPLIT_SCHUR_PRE_USER
    FULL                 ///< PC_FIELDSPLIT_SCHUR_PRE_FULL
};

/**
 * @brief PETSc-based physics preconditioning for coupled systems
 * 
 * Uses PETSc's PC infrastructure for all preconditioning operations.
 * Supports block methods, field split, Schur complement, and multigrid.
 */
class PETScPhysicsPreconditioner {
public:
    struct FieldInfo {
        std::string name;
        IS index_set;              ///< PETSc Index Set for this field
        PetscInt num_components;
        std::string sub_pc_type;   ///< Preconditioner for this block
    };
    
    struct Parameters {
        PreconditionerType type = PreconditionerType::FIELDSPLIT_SCHUR;
        
        // Field split options
        std::vector<FieldInfo> fields;
        
        // Schur complement options
        SchurFactType schur_fact = SchurFactType::FULL;
        SchurPreType schur_pre = SchurPreType::SELFP;
        PetscReal schur_scale = 1.0;
        
        // Multigrid options
        PetscInt mg_levels = 3;
        std::string mg_smoother = "sor";
        PetscInt mg_smooth_down = 2;
        PetscInt mg_smooth_up = 2;
        
        // GAMG options
        PetscInt gamg_agg_nsmooths = 1;
        PetscReal gamg_threshold = 0.0;
    };
    
    PETScPhysicsPreconditioner(MPI_Comm comm);
    ~PETScPhysicsPreconditioner();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Setup preconditioner on KSP
     */
    PetscErrorCode setup(KSP ksp, Mat A);
    
    /**
     * @brief Add a field for field split
     */
    PetscErrorCode addField(const std::string& name, IS is, 
                            PetscInt num_components,
                            const std::string& sub_pc = "ilu");
    
    /**
     * @brief Create index sets from DM section
     */
    PetscErrorCode createFieldISFromDM(DM dm, const std::vector<std::string>& field_names);
    
    /**
     * @brief Configure field split preconditioner
     */
    PetscErrorCode setupFieldSplit(PC pc);
    
    /**
     * @brief Configure Schur complement
     */
    PetscErrorCode setupSchurComplement(PC pc);
    
    /**
     * @brief Configure multigrid
     */
    PetscErrorCode setupMultigrid(PC pc);
    
    /**
     * @brief Get underlying PC for advanced configuration
     */
    PetscErrorCode getPC(KSP ksp, PC* pc);
    
private:
    MPI_Comm comm_;
    Parameters params_;
    std::vector<IS> field_is_;
    
    PetscErrorCode setPCType(PC pc);
};


// =============================================================================
// High-Fidelity Numerics Configuration
// =============================================================================

/**
 * @brief Configuration for high-fidelity numerics
 */
struct HighFidelityNumericsConfig {
    bool enable_high_fidelity_numerics = false;
    
    // SUPG/GLS stabilization
    bool enable_supg = false;
    SUPGStabilization::Parameters supg_params;
    
    // PETSc time integration
    bool enable_advanced_time = false;
    PETScTimeIntegration::Parameters time_params;
    
    // p-adaptivity
    bool enable_p_adapt = false;
    PAdaptivityManager::Parameters p_adapt_params;
    
    // Physics-based preconditioning
    bool enable_physics_precond = false;
    PETScPhysicsPreconditioner::Parameters precond_params;
    
    void parseConfig(const std::map<std::string, std::string>& config);
    PetscErrorCode validate(std::string& error_msg) const;
};

/**
 * @brief Factory functions
 */
std::unique_ptr<SUPGStabilization> createSUPGStabilization(
    const std::map<std::string, std::string>& config);

std::unique_ptr<PETScTimeIntegration> createTimeIntegration(
    MPI_Comm comm, const std::map<std::string, std::string>& config);

std::unique_ptr<PETScPhysicsPreconditioner> createPhysicsPreconditioner(
    MPI_Comm comm, const std::map<std::string, std::string>& config);

// For DG: use createDGSolver() from DiscontinuousGalerkin.hpp
// For AMR: use AdaptiveMeshRefinement from AdaptiveMeshRefinement.hpp

} // namespace HighFidelity
} // namespace FSRM

#endif // HIGH_FIDELITY_NUMERICS_HPP
