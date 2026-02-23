#ifndef DISCONTINUOUS_GALERKIN_HPP
#define DISCONTINUOUS_GALERKIN_HPP

#include "FSRM.hpp"
#include <petscfe.h>
#include <vector>
#include <array>
#include <memory>
#include <functional>

namespace FSRM {

/**
 * @brief Polynomial order for DG method
 */
enum class DGOrder {
    O1 = 1,  // Linear (2nd order accurate)
    O2 = 2,  // Quadratic (3rd order accurate)
    O3 = 3,  // Cubic (4th order accurate)
    O4 = 4,  // Quartic (5th order accurate)
    O5 = 5,  // Quintic (6th order accurate)
    O6 = 6,  // Sextic (7th order accurate)
    O7 = 7,  // Septic (8th order accurate)
    O8 = 8,  // Octic (9th order accurate)
    O9 = 9,  // Nonic (10th order accurate)
    O10 = 10 // Decic (11th order accurate)
};

/**
 * @brief Flux computation method for DG interfaces
 */
enum class FluxMethod {
    GODUNOV,           // Godunov (exact Riemann solver)
    ROE,               // Roe's approximate Riemann solver
    HLL,               // Harten-Lax-van Leer
    HLLC,              // HLL-Contact
    RUSANOV,           // Rusanov (local Lax-Friedrichs)
    PENALTY,           // Penalty method (for diffusion)
    LDG,               // Local Discontinuous Galerkin
    BR2,               // Bassi-Rebay 2
    IP                 // Interior Penalty
};

/**
 * @brief Basis function type for DG
 */
enum class BasisType {
    LAGRANGE,          // Lagrange interpolation polynomials
    LEGENDRE,          // Legendre orthogonal polynomials
    DUBINER,           // Dubiner (hierarchical) basis for triangles/tets
    MODAL,             // Modal basis (integrated modes)
    NODAL              // Nodal basis (point values)
};

/**
 * @brief Reference element type
 */
enum class ElementType {
    INTERVAL,          // 1D interval [0,1]
    TRIANGLE,          // 2D triangle
    QUADRILATERAL,     // 2D quadrilateral
    TETRAHEDRON,       // 3D tetrahedron
    HEXAHEDRON,        // 3D hexahedron
    PRISM,             // 3D prism (wedge)
    PYRAMID            // 3D pyramid
};

/**
 * @brief Quadrature rule for numerical integration
 */
struct QuadratureRule {
    std::vector<std::array<double, 3>> points;  // Quadrature points in reference element
    std::vector<double> weights;                 // Quadrature weights
    int order;                                   // Exact integration order
    
    // Factory methods for common rules
    static QuadratureRule gaussLegendre(int n_points, int dim = 1);
    static QuadratureRule gaussLobatto(int n_points, int dim = 1);
    static QuadratureRule stroud(int order, ElementType elem_type);
    static QuadratureRule dunavantTriangle(int order);
    static QuadratureRule grundmannMoellerTet(int order);
};

/**
 * @brief Basis functions on reference element
 */
class BasisFunctions {
public:
    BasisFunctions(BasisType type, DGOrder order, ElementType elem_type);
    
    // Number of basis functions (degrees of freedom per element)
    int getNumBasis() const { return num_basis; }
    
    // Evaluate basis function i at point xi
    double evaluate(int i, const std::array<double, 3>& xi) const;
    
    // Evaluate gradient of basis function i at point xi
    std::array<double, 3> evaluateGradient(int i, const std::array<double, 3>& xi) const;
    
    // Get all basis functions at point xi
    void evaluateAll(const std::array<double, 3>& xi, std::vector<double>& phi) const;
    
    // Get all gradients at point xi
    void evaluateAllGradients(const std::array<double, 3>& xi, 
                             std::vector<std::array<double, 3>>& grad_phi) const;
    
    // Mass matrix (integral of phi_i * phi_j)
    const std::vector<double>& getMassMatrix() const { return mass_matrix; }
    
    // Stiffness matrix (integral of grad_phi_i · grad_phi_j)
    const std::vector<double>& getStiffnessMatrix() const { return stiffness_matrix; }
    
    // Differentiation matrix (for spectral methods)
    const std::vector<double>& getDifferentiationMatrix(int dir) const;
    
private:
    BasisType basis_type;
    DGOrder order;
    ElementType elem_type;
    int num_basis;
    
    // Precomputed matrices
    std::vector<double> mass_matrix;
    std::vector<double> stiffness_matrix;
    std::vector<std::vector<double>> diff_matrices;  // One per spatial dimension
    
    // Basis function implementations
    double lagrangeBasis(int i, const std::array<double, 3>& xi) const;
    double legendreBasis(int i, const std::array<double, 3>& xi) const;
    double dubinerBasis(int i, const std::array<double, 3>& xi) const;
    
    void precomputeMatrices();
};

/**
 * @brief Riemann solver for computing interface fluxes
 */
class RiemannSolver {
public:
    RiemannSolver(FluxMethod method);
    
    // Compute numerical flux at interface
    // uL, uR: left and right states
    // normal: outward normal vector
    // flux: output flux (modified in place)
    void computeFlux(const std::vector<double>& uL, 
                    const std::vector<double>& uR,
                    const std::array<double, 3>& normal,
                    std::vector<double>& flux) const;
    
    // Compute flux and wave speeds
    void computeFluxAndSpeeds(const std::vector<double>& uL,
                              const std::vector<double>& uR,
                              const std::array<double, 3>& normal,
                              std::vector<double>& flux,
                              double& sL, double& sR) const;
    
    // Set equation-specific flux function
    void setPhysicalFlux(std::function<void(const std::vector<double>&,
                                           const std::array<double, 3>&,
                                           std::vector<double>&)> func);
    
    // Set wave speed function
    void setWaveSpeed(std::function<double(const std::vector<double>&,
                                          const std::array<double, 3>&)> func);
    
private:
    FluxMethod method;
    
    // User-defined physics-specific functions
    std::function<void(const std::vector<double>&,
                      const std::array<double, 3>&,
                      std::vector<double>&)> physical_flux;
    std::function<double(const std::vector<double>&,
                        const std::array<double, 3>&)> wave_speed;
    
    // Specific Riemann solver implementations
    void godunovFlux(const std::vector<double>& uL, const std::vector<double>& uR,
                     const std::array<double, 3>& normal, std::vector<double>& flux) const;
    void roeFlux(const std::vector<double>& uL, const std::vector<double>& uR,
                 const std::array<double, 3>& normal, std::vector<double>& flux) const;
    void hllFlux(const std::vector<double>& uL, const std::vector<double>& uR,
                 const std::array<double, 3>& normal, std::vector<double>& flux) const;
    void hllcFlux(const std::vector<double>& uL, const std::vector<double>& uR,
                  const std::array<double, 3>& normal, std::vector<double>& flux) const;
    void rusanovFlux(const std::vector<double>& uL, const std::vector<double>& uR,
                     const std::array<double, 3>& normal, std::vector<double>& flux) const;
};

/**
 * @brief DG spatial discretization
 * 
 * This class implements the discontinuous Galerkin method for spatial
 * discretization, which is the foundation of SeisSol's high accuracy.
 * 
 * Key features:
 * - Arbitrary high-order accuracy (up to O10)
 * - Element-local operations (perfect for parallelization)
 * - Natural handling of discontinuities (shocks, faults)
 * - hp-adaptivity support
 */
class DiscontinuousGalerkin {
public:
    DiscontinuousGalerkin(MPI_Comm comm);
    ~DiscontinuousGalerkin();
    
    // Configuration
    void setOrder(DGOrder order);
    void setBasisType(BasisType type);
    void setFluxMethod(FluxMethod method);
    void setQuadratureOrder(int order);
    
    // Mesh setup
    PetscErrorCode setupMesh(DM dm);
    PetscErrorCode setupBasisFunctions();
    PetscErrorCode setupQuadratureRules();
    
    // Spatial operators
    PetscErrorCode assembleVolumeIntegrals(Vec U, Vec R);
    PetscErrorCode assembleSurfaceIntegrals(Vec U, Vec R);
    PetscErrorCode assembleBoundaryConditions(Vec U, Vec R, double t);
    
    // Complete spatial operator: dU/dt = L(U)
    PetscErrorCode spatialOperator(Vec U, Vec R, double t);
    
    // Jacobian assembly (for implicit methods)
    PetscErrorCode assembleJacobian(Vec U, Mat J, double shift);
    
    // Mass matrix operations
    PetscErrorCode applyMassMatrix(Vec X, Vec Y);
    PetscErrorCode applyInverseMassMatrix(Vec X, Vec Y);
    
    // Projection operations
    PetscErrorCode projectFunction(std::function<void(const double*, double*)> func, Vec U);
    PetscErrorCode projectDerivative(Vec U, int direction, Vec DU);
    
    // Limiters (for shock capturing)
    PetscErrorCode applyLimiter(Vec U);
    PetscErrorCode applyTVBLimiter(Vec U, double M);
    PetscErrorCode applyMinModLimiter(Vec U);
    
    // Adaptive mesh refinement indicators
    double computeErrorIndicator(Vec U, int element_id);
    PetscErrorCode computeErrorMap(Vec U, Vec error_map);
    
    // Accessors
    DGOrder getOrder() const { return order; }
    int getNumDofsPerElement() const;
    int getTotalDofs() const;
    
    // Get basis functions for element
    const BasisFunctions* getBasis(int element_id) const;
    
    // Get Riemann solver
    RiemannSolver* getRiemannSolver() { return riemann_solver.get(); }
    
private:
    MPI_Comm comm;
    int rank, size;
    
    // Configuration
    DGOrder order;
    BasisType basis_type;
    FluxMethod flux_method;
    int quadrature_order;
    
    // Mesh
    DM dm;
    int num_elements;
    int num_faces;
    int num_dofs_per_elem;
    int total_dofs;
    
    // Basis functions (may vary per element for hp-adaptivity)
    std::vector<std::unique_ptr<BasisFunctions>> basis_functions;
    
    // Quadrature rules
    std::unique_ptr<QuadratureRule> volume_quadrature;
    std::unique_ptr<QuadratureRule> surface_quadrature;
    
    // Riemann solver
    std::unique_ptr<RiemannSolver> riemann_solver;
    
    // Mass matrix (block diagonal, one block per element)
    Mat mass_matrix;
    Mat inverse_mass_matrix;
    
    // Boundary conditions
    std::map<int, std::function<void(const double*, double*, double)>> boundary_conditions;
    
    // Helper functions
    PetscErrorCode getElementDoFs(int elem_id, Vec U, double* u_elem);
    PetscErrorCode setElementDoFs(int elem_id, const double* u_elem, Vec U);
    PetscErrorCode getElementGeometry(int elem_id, std::vector<double>& coords);
    PetscErrorCode computeJacobian(int elem_id, const std::array<double, 3>& xi,
                                   std::array<std::array<double, 3>, 3>& J);
    
    // Flux computation at element interfaces
    PetscErrorCode computeInterfaceFlux(int face_id, Vec U, 
                                       std::vector<double>& flux);
    
    // Limiter implementations
    double detectTroubledCell(const double* u_elem, int elem_id);
    void applyLimiterToElement(double* u_elem, int elem_id, double indicator);
};

/**
 * @brief ADER time integration
 * 
 * ADER (Arbitrary high-order DERivatives) provides single-step time
 * integration matching the spatial accuracy of the DG method.
 * This is a key component of SeisSol's efficiency.
 * 
 * Key features:
 * - Single-step method (no multi-stage RK)
 * - Arbitrary high-order accuracy
 * - Optimal CFL number usage
 * - Predictor-corrector structure
 */
class ADERTimeIntegrator {
public:
    ADERTimeIntegrator(DiscontinuousGalerkin* dg, int order);
    
    // Single time step: U^{n+1} = ADER(U^n, dt)
    PetscErrorCode step(Vec U, double dt, double t);
    
    // Compute stable time step from CFL condition
    double computeTimeStep(Vec U, double cfl_number = 0.5);
    
    // Get order of accuracy
    int getOrder() const { return order; }
    
private:
    DiscontinuousGalerkin* dg;
    int order;
    
    // ADER predictor: compute space-time predictor q^*
    PetscErrorCode predictor(Vec U, Vec Q_star, double dt, double t);
    
    // ADER corrector: update U using predictor
    PetscErrorCode corrector(Vec U, Vec Q_star, double dt, double t);
    
    // Cauchy-Kowalevski procedure for time derivatives
    PetscErrorCode cauchyKowalevski(Vec U, std::vector<Vec>& time_derivatives, double t);
    
    // Space-time quadrature
    struct SpaceTimeQuadrature {
        std::vector<double> time_weights;
        std::vector<double> time_points;
        QuadratureRule space_rule;
    };
    SpaceTimeQuadrature space_time_quad;
    
    // Temporary vectors
    Vec Q_star;  // Space-time predictor
    std::vector<Vec> time_derivatives;  // U, dU/dt, d²U/dt², ...
};

/**
 * @brief Local Time Stepping (LTS) manager
 * 
 * LTS allows different elements to use different time steps based on
 * their local CFL condition, leading to significant speedups (5-10x)
 * for meshes with varying element sizes.
 * 
 * Key features:
 * - Clustered LTS (rate-2, rate-3, or arbitrary)
 * - Automatic cluster optimization
 * - Wiggle factor tuning
 * - Load balancing
 */
class LocalTimeStepping {
public:
    LocalTimeStepping(DiscontinuousGalerkin* dg);
    
    // Configuration
    void setRate(int rate);  // Rate-N LTS
    void setMaxClusters(int max_clusters);
    void enableWiggleFactor(bool enable, double min_factor = 0.51);
    void enableAutoMerge(bool enable, double performance_loss_threshold = 0.01);
    
    // Setup LTS clusters
    PetscErrorCode setupClusters(Vec U);
    PetscErrorCode optimizeClusters();
    
    // Time stepping with LTS
    PetscErrorCode step(Vec U, double dt_global, double t);
    
    // Analysis
    int getNumClusters() const { return num_clusters; }
    double getSpeedup() const;
    void printClusterStatistics() const;
    
    // Get cluster assignment for element
    int getElementCluster(int elem_id) const;
    
private:
    DiscontinuousGalerkin* dg;
    
    // Configuration
    int rate;
    int max_clusters;
    bool use_wiggle_factor;
    double wiggle_factor_min;
    bool auto_merge_clusters;
    double perf_loss_threshold;
    
    // Cluster data
    int num_clusters;
    std::vector<int> element_to_cluster;      // Cluster assignment for each element
    std::vector<std::vector<int>> cluster_elements;  // Elements in each cluster
    std::vector<double> cluster_dt;           // Time step for each cluster
    std::vector<int> cluster_update_frequency; // Update frequency relative to smallest cluster
    
    // Time tracking
    std::vector<double> cluster_time;         // Current time for each cluster
    std::vector<int> cluster_steps;           // Steps taken by each cluster
    
    // Helper functions
    PetscErrorCode computeElementTimeSteps(Vec U, std::vector<double>& dt_elements);
    PetscErrorCode clusterElements(const std::vector<double>& dt_elements);
    PetscErrorCode enforceMaxDifferenceProperty();
    double optimizeWiggleFactor(const std::vector<double>& dt_elements);
    double estimateCost(const std::vector<int>& element_to_cluster);
    void mergeClusters(int cluster1, int cluster2);
    
    // Update a single cluster
    PetscErrorCode updateCluster(int cluster_id, Vec U, double dt, double t);
    
    // Interface flux handling for LTS
    PetscErrorCode handleInterfacesBetweenClusters(Vec U, double t);
};

/**
 * @brief Factory function to create DG solver with ADER-LTS
 */
std::unique_ptr<DiscontinuousGalerkin> createDGSolver(
    MPI_Comm comm,
    DGOrder order = DGOrder::O3,
    BasisType basis = BasisType::NODAL,
    FluxMethod flux = FluxMethod::RUSANOV);

/**
 * @brief Factory function to create ADER time integrator
 */
std::unique_ptr<ADERTimeIntegrator> createADERIntegrator(
    DiscontinuousGalerkin* dg,
    int order = 4);

/**
 * @brief Factory function to create LTS manager
 */
std::unique_ptr<LocalTimeStepping> createLTS(
    DiscontinuousGalerkin* dg,
    int rate = 2);

// =============================================================================
// GPU-Accelerated DG Solver Interface
// =============================================================================

/**
 * @brief GPU execution context for the DG wave equation solver
 *
 * Manages device memory and kernel dispatch for the ADER-DG method
 * running on NVIDIA CUDA GPUs. The DG method is embarrassingly parallel
 * over elements, achieving 30-50× speedup on modern GPUs.
 *
 * ## Usage
 * ```cpp
 * DGGPUContext gpu_ctx;
 * gpu_ctx.initialize(dg_solver, n_elements, n_dof, N_VAR_ELASTIC);
 * gpu_ctx.uploadSolution(host_Q);
 * gpu_ctx.timeStep(dt);  // all DG operations on GPU
 * gpu_ctx.downloadSolution(host_Q);
 * ```
 *
 * ## Memory Layout
 * Element-major: data[elem * n_dof * n_var + dof * n_var + var]
 */
class DGGPUContext {
public:
    DGGPUContext();
    ~DGGPUContext();

    // Disable copy
    DGGPUContext(const DGGPUContext&) = delete;
    DGGPUContext& operator=(const DGGPUContext&) = delete;

    /// Check if GPU acceleration is available
    bool isAvailable() const { return gpu_available_; }

    /**
     * @brief Initialize GPU context and allocate device memory
     *
     * @param n_elements  Total number of DG elements
     * @param n_dof       Degrees of freedom per element
     * @param n_var       Conservative variables per DOF (9 for elastic)
     * @param n_interior_faces Number of interior faces
     * @param n_pml_elements   Number of PML boundary elements
     * @return true on success
     */
    bool initialize(int n_elements, int n_dof, int n_var,
                    int n_interior_faces = 0, int n_pml_elements = 0);

    /// Free all GPU memory
    void finalize();

    // =========================================================================
    // Data Transfer
    // =========================================================================

    /// Upload solution from host to GPU
    void uploadSolution(const double* host_Q);

    /// Download solution from GPU to host
    void downloadSolution(double* host_Q) const;

    /// Upload material properties
    void uploadMaterials(const void* host_materials, size_t bytes);

    /// Upload reference element matrices (stiffness, mass inverse)
    void uploadReferenceMatrices(const double* stiffness_x,
                                 const double* stiffness_y,
                                 const double* stiffness_z,
                                 const double* mass_inv,
                                 int n_dof);

    /// Upload mesh geometry (Jacobians, face connectivity, normals)
    void uploadGeometry(const double* inv_jacobians,
                        const double* det_jacobians,
                        const int* face_neighbors,
                        const double* face_normals,
                        const double* face_areas);

    // =========================================================================
    // DG Time Stepping on GPU
    // =========================================================================

    /**
     * @brief Execute one complete ADER-DG time step on GPU
     *
     * Performs: ADER predictor → volume integral → numerical flux →
     *          source injection → PML → solution update
     *
     * @param dt         Time step size
     * @param ader_order ADER predictor order (0 = forward Euler)
     */
    void timeStep(double dt, int ader_order = 0);

    /// Compute volume integrals on GPU
    void computeVolumeIntegral();

    /// Compute numerical fluxes on GPU
    void computeNumericalFlux();

    /// Apply source term
    void applySource(int source_elem, const double* moment_tensor,
                     double amplitude, const double* basis_at_source);

    /// Apply PML damping
    void applyPML(double dt);

    /// Update solution
    void updateSolution(double dt);

    /// Extract receiver data
    void extractReceivers(const int* receiver_elems,
                          const double* receiver_basis,
                          double* receiver_data,
                          int n_receivers);

    // =========================================================================
    // Diagnostics
    // =========================================================================

    /// Compute L2 norm of a solution component
    double computeNorm(int component) const;

    /// Get GPU memory usage in bytes
    size_t getMemoryUsage() const { return total_gpu_memory_; }

private:
    bool gpu_available_ = false;
    int n_elem_ = 0;
    int n_dof_ = 0;
    int n_var_ = 0;
    int n_faces_ = 0;
    int n_pml_ = 0;
    size_t total_gpu_memory_ = 0;

#ifdef USE_CUDA
    // Solution arrays (device)
    double* d_Q_ = nullptr;           ///< Solution [n_elem × n_dof × n_var]
    double* d_R_vol_ = nullptr;       ///< Volume residual
    double* d_R_surf_ = nullptr;      ///< Surface residual
    double* d_R_src_ = nullptr;       ///< Source residual
    double* d_Q_predicted_ = nullptr; ///< ADER predictor

    // Material properties (device)
    void* d_material_ = nullptr;

    // Reference element matrices (device)
    double* d_stiffness_x_ = nullptr;
    double* d_stiffness_y_ = nullptr;
    double* d_stiffness_z_ = nullptr;
    double* d_mass_inv_ = nullptr;

    // Mesh geometry (device)
    double* d_inv_jacobian_ = nullptr;
    double* d_det_jacobian_ = nullptr;
    int* d_face_neighbors_ = nullptr;
    double* d_face_normals_ = nullptr;
    double* d_face_area_ = nullptr;

    // PML (device)
    double* d_Q_pml_ = nullptr;
    int* d_pml_elem_ids_ = nullptr;
    double* d_pml_damping_ = nullptr;
#endif
};

} // namespace FSRM

#endif // DISCONTINUOUS_GALERKIN_HPP
