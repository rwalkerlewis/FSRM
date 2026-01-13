/**
 * @file BoundaryConditions.hpp
 * @brief Boundary conditions for seismic wave propagation
 * 
 * Implements:
 * - Perfectly Matched Layers (PML) - absorbing boundaries
 * - Clayton-Engquist absorbing boundary conditions
 * - Free surface boundary conditions
 * - Periodic boundary conditions
 * - Symmetry/anti-symmetry conditions
 */

#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <vector>
#include <array>
#include <memory>
#include <functional>
#include <map>
#include <string>
#include <complex>

namespace FSRM {

// Forward declarations
class DiscontinuousGalerkin;

/**
 * @brief Types of boundary conditions
 */
enum class BoundaryType {
    FREE_SURFACE,           // Traction-free surface (σ·n = 0)
    ABSORBING_PML,         // Perfectly Matched Layer
    ABSORBING_CE,          // Clayton-Engquist (paraxial approximation)
    ABSORBING_LYSMER,      // Lysmer-Kuhlemeyer dashpot
    DIRICHLET,             // Fixed displacement (u = 0)
    NEUMANN,               // Prescribed traction
    PERIODIC,              // Periodic boundary
    SYMMETRY,              // Symmetry plane
    ANTISYMMETRY,          // Anti-symmetry plane
    FAULT,                 // Internal fault boundary
    INTERFACE              // Material interface
};

/**
 * @brief Boundary face information
 */
struct BoundaryFace {
    int element_id;                    // Owner element
    int local_face_id;                 // Local face number (0-3 for tet)
    BoundaryType type;                 // Boundary condition type
    
    std::array<double, 3> centroid;    // Face centroid
    std::array<double, 3> normal;      // Outward normal
    double area;                       // Face area
    
    std::vector<int> nodes;            // Node indices on face
    std::vector<std::array<double, 3>> node_coords;  // Node coordinates
    
    // For matched faces (periodic, interface)
    int matched_element_id;
    int matched_face_id;
    
    // Material properties (for absorbing BC)
    double rho;                        // Density
    double vp;                         // P-wave velocity
    double vs;                         // S-wave velocity
    double impedance_p;                // P-wave impedance (ρ * vp)
    double impedance_s;                // S-wave impedance (ρ * vs)
};

/**
 * @brief Base class for boundary condition implementations
 */
class BoundaryConditionBase {
public:
    virtual ~BoundaryConditionBase() = default;
    
    /**
     * @brief Apply boundary condition to residual
     * @param face Boundary face information
     * @param u Solution at face nodes [ndof]
     * @param u_ext Ghost/external state [ndof]
     * @param t Current time
     */
    virtual void applyToResidual(const BoundaryFace& face,
                                const double* u, double* residual,
                                double t) const = 0;
    
    /**
     * @brief Compute boundary flux for DG
     * @param face Boundary face information
     * @param u_int Interior state
     * @param flux_n Normal flux [nvar]
     * @param t Current time
     */
    virtual void computeBoundaryFlux(const BoundaryFace& face,
                                     const double* u_int, double* flux_n,
                                     double t) const = 0;
    
    /**
     * @brief Get boundary state (ghost/external state)
     * @param face Boundary face information
     * @param u_int Interior state
     * @param u_ext External/ghost state (output)
     * @param t Current time
     */
    virtual void getBoundaryState(const BoundaryFace& face,
                                  const double* u_int, double* u_ext,
                                  double t) const = 0;
    
    /**
     * @brief Check if this BC requires auxiliary PDE
     */
    virtual bool requiresAuxiliaryPDE() const { return false; }
    
    /**
     * @brief Get number of state variables for this BC
     */
    virtual int getNumStateVariables() const { return 0; }
};

/**
 * @brief Free surface boundary condition
 * 
 * Implements traction-free surface: σ·n = 0
 * For velocity-stress formulation, reflects waves
 */
class FreeSurfaceBC : public BoundaryConditionBase {
public:
    FreeSurfaceBC();
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;
    
    /**
     * @brief Apply mirror condition for stress tensor
     * @param sigma_int Interior stress [6]: σxx, σyy, σzz, σxy, σxz, σyz
     * @param n Normal vector
     * @param sigma_ext Exterior stress (output)
     */
    void applyMirrorStress(const double* sigma_int, const double* n,
                          double* sigma_ext) const;
    
    /**
     * @brief Compute surface traction
     */
    void computeSurfaceTraction(const double* sigma, const double* n,
                               double* traction) const;

private:
    bool use_characteristic_decomposition;
};

/**
 * @brief PML (Perfectly Matched Layer) absorbing boundary
 * 
 * Based on Bérenger's split-field PML with:
 * - Complex frequency-shifted (CFS) formulation for long-time stability
 * - Multi-axial implementation (MPML)
 * - Optimized damping profile
 */
class PMLBoundary : public BoundaryConditionBase {
public:
    PMLBoundary();
    
    /**
     * @brief Configure PML parameters
     */
    void setThickness(double d);                    // PML thickness
    void setOrder(int m);                           // Polynomial order for damping profile
    void setReflectionCoefficient(double R);        // Target reflection coefficient
    void setCFSParameters(double alpha_max, double kappa_max);  // CFS parameters
    void setDomainBounds(double x_min, double x_max,
                        double y_min, double y_max,
                        double z_min, double z_max);
    void setMaterialProperties(double rho, double vp, double vs);
    
    /**
     * @brief Determine if point is in PML region
     */
    bool isInPML(double x, double y, double z) const;
    int getPMLDirection(double x, double y, double z) const;  // Returns bit flags
    
    /**
     * @brief Get damping profile at point
     */
    void getDampingCoefficients(double x, double y, double z,
                               double& d_x, double& d_y, double& d_z) const;
    
    /**
     * @brief Get CFS stretching functions
     */
    void getStretchingFunctions(double x, double y, double z,
                               std::complex<double>& s_x,
                               std::complex<double>& s_y,
                               std::complex<double>& s_z) const;
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;
    
    bool requiresAuxiliaryPDE() const override { return true; }
    int getNumStateVariables() const override;
    
    /**
     * @brief Initialize PML auxiliary variables
     */
    void initializeAuxiliaryVariables(int num_pml_elements);
    
    /**
     * @brief Get auxiliary variable RHS for PML evolution
     * @param elem_id Element in PML
     * @param u Solution state
     * @param aux Current auxiliary variables
     * @param aux_rhs RHS for auxiliary variables
     */
    void computeAuxiliaryRHS(int elem_id, const double* u,
                            const double* aux, double* aux_rhs);
    
    /**
     * @brief Update auxiliary variables
     */
    void updateAuxiliaryVariables(int elem_id, const double* aux_new);
    
    /**
     * @brief Get modified PML flux Jacobian
     */
    void getPMLJacobian(int elem_id, double* jacobian) const;

private:
    double thickness;                // PML thickness [m]
    int damping_order;              // Polynomial order (typically 2-4)
    double reflection_coef;          // Target reflection coefficient
    
    // CFS parameters
    double alpha_max;                // Max frequency shift
    double kappa_max;                // Max coordinate stretching
    
    // Domain bounds
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
    
    // Material properties
    double rho_ref;
    double vp_ref;
    double vs_ref;
    
    // Computed damping coefficient
    double d_max;                    // Maximum damping
    
    // Auxiliary variables for split-field PML
    std::vector<std::vector<double>> aux_vars;  // [element][var]
    std::vector<int> pml_element_map;  // Maps global element to PML element
    
    // Helper functions
    double dampingProfile(double dist) const;
    double alphaProfile(double dist) const;
    double kappaProfile(double dist) const;
};

/**
 * @brief Clayton-Engquist absorbing boundary condition
 * 
 * Paraxial approximation for absorbing waves at boundary.
 * First and second order implementations.
 */
class ClaytonEngquistBC : public BoundaryConditionBase {
public:
    ClaytonEngquistBC();
    
    /**
     * @brief Set approximation order (1 or 2)
     */
    void setOrder(int order);
    
    /**
     * @brief Set material properties at boundary
     */
    void setMaterialProperties(double rho, double vp, double vs);
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;
    
    /**
     * @brief Compute impedance matrix for traction BC
     * @param n Normal vector
     * @param Z Impedance matrix [3x3]
     */
    void computeImpedanceMatrix(const double* n, double* Z) const;
    
    /**
     * @brief Apply first-order paraxial approximation
     * Traction ≈ -Z · ∂u/∂t
     */
    void applyFirstOrder(const BoundaryFace& face,
                        const double* v, double* traction) const;
    
    /**
     * @brief Apply second-order paraxial approximation
     * Includes tangential derivatives
     */
    void applySecondOrder(const BoundaryFace& face,
                         const double* v, const double* dv_dt,
                         double* traction) const;

private:
    int order;                       // 1 or 2
    double rho;                      // Density
    double vp;                       // P-wave velocity
    double vs;                       // S-wave velocity
    double lambda;                   // Lamé parameter
    double mu;                       // Shear modulus
};

/**
 * @brief Lysmer-Kuhlemeyer dashpot boundary condition
 * 
 * Simple viscous absorbing boundary: t = -ρ·c·v
 */
class LysmerKuhlemeyer : public BoundaryConditionBase {
public:
    LysmerKuhlemeyer();
    
    void setMaterialProperties(double rho, double vp, double vs);
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;
    
    /**
     * @brief Compute absorbing traction
     * @param v Velocity at boundary
     * @param n Normal direction
     * @param t Traction (output)
     */
    void computeAbsorbingTraction(const double* v, const double* n,
                                  double* traction) const;

private:
    double rho;
    double vp;
    double vs;
    double impedance_p;
    double impedance_s;
};

/**
 * @brief Dirichlet (fixed) boundary condition
 */
class DirichletBC : public BoundaryConditionBase {
public:
    DirichletBC();
    
    /**
     * @brief Set prescribed displacement/velocity
     */
    void setPrescribedValue(std::function<void(double, double, double, double, double*)> func);
    
    /**
     * @brief Set to zero (homogeneous)
     */
    void setZero();
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;

private:
    std::function<void(double, double, double, double, double*)> prescribed_func;
};

/**
 * @brief Neumann (traction) boundary condition
 */
class NeumannBC : public BoundaryConditionBase {
public:
    NeumannBC();
    
    /**
     * @brief Set prescribed traction
     */
    void setPrescribedTraction(std::function<void(double, double, double, double, double*)> func);
    
    /**
     * @brief Set constant traction
     */
    void setConstantTraction(double tx, double ty, double tz);
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;

private:
    std::function<void(double, double, double, double, double*)> traction_func;
    bool is_constant;
    std::array<double, 3> constant_traction;
};

/**
 * @brief Periodic boundary condition
 */
class PeriodicBC : public BoundaryConditionBase {
public:
    PeriodicBC();
    
    /**
     * @brief Set periodic pair mapping
     */
    void setPeriodicMapping(const std::vector<std::pair<int, int>>& face_pairs);
    
    /**
     * @brief Set periodic direction and period
     */
    void setPeriodicDirection(int dir, double period);
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;

private:
    std::map<int, int> periodic_pairs;  // face -> matched face
    int periodic_dir;
    double period;
};

/**
 * @brief Symmetry boundary condition
 */
class SymmetryBC : public BoundaryConditionBase {
public:
    SymmetryBC(bool is_antisymmetric = false);
    
    void applyToResidual(const BoundaryFace& face,
                        const double* u, double* residual,
                        double t) const override;
    
    void computeBoundaryFlux(const BoundaryFace& face,
                            const double* u_int, double* flux_n,
                            double t) const override;
    
    void getBoundaryState(const BoundaryFace& face,
                         const double* u_int, double* u_ext,
                         double t) const override;
    
    /**
     * @brief Apply reflection to velocity
     */
    void reflectVelocity(const double* v_int, const double* n,
                        double* v_ext) const;
    
    /**
     * @brief Apply reflection to stress tensor
     */
    void reflectStress(const double* sigma_int, const double* n,
                      double* sigma_ext) const;

private:
    bool antisymmetric;  // True for anti-symmetry (opposite sign)
};

/**
 * @brief Boundary condition manager
 * 
 * Manages all boundary conditions for a domain
 */
class BoundaryConditionManager {
public:
    BoundaryConditionManager();
    ~BoundaryConditionManager();
    
    /**
     * @brief Initialize from DG discretization
     */
    void initialize(DiscontinuousGalerkin* dg);
    
    /**
     * @brief Add boundary condition to face group
     * @param tag Boundary tag from mesh
     * @param bc Boundary condition implementation
     */
    void addBoundaryCondition(int tag, std::unique_ptr<BoundaryConditionBase> bc);
    
    /**
     * @brief Set default boundary condition for untagged faces
     */
    void setDefaultBC(std::unique_ptr<BoundaryConditionBase> bc);
    
    /**
     * @brief Classify boundary faces
     * @param face_tags Tags from mesh for each boundary face
     */
    void classifyBoundaryFaces(const std::vector<int>& face_tags);
    
    /**
     * @brief Setup PML regions
     * @param pml_thickness Thickness of PML layer
     * @param directions Bit flags for PML directions (1=x, 2=y, 4=z)
     */
    void setupPML(double pml_thickness, int directions = 7);
    
    /**
     * @brief Get PML boundary object
     */
    PMLBoundary* getPML() { return pml.get(); }
    
    /**
     * @brief Apply boundary conditions to residual
     */
    void applyBoundaryConditions(double* residual, const double* solution,
                                double t);
    
    /**
     * @brief Compute boundary fluxes for DG
     */
    void computeBoundaryFluxes(const double* solution, double* flux,
                              double t);
    
    /**
     * @brief Get boundary faces of a given type
     */
    const std::vector<BoundaryFace>& getBoundaryFaces(BoundaryType type) const;
    
    /**
     * @brief Get all boundary faces
     */
    const std::vector<BoundaryFace>& getAllBoundaryFaces() const;
    
    /**
     * @brief Get boundary condition for face
     */
    BoundaryConditionBase* getBCForFace(int face_id) const;
    
    /**
     * @brief Update material properties on boundaries
     */
    void updateMaterialProperties();
    
    /**
     * @brief Check if any face has PML
     */
    bool hasPML() const { return pml != nullptr; }
    
    /**
     * @brief Get number of PML auxiliary variables
     */
    int getNumPMLAuxVariables() const;
    
    /**
     * @brief Evolve PML auxiliary variables
     */
    void evolvePMLAuxVariables(const double* solution, double dt);

private:
    DiscontinuousGalerkin* dg;
    
    // Boundary conditions by tag
    std::map<int, std::unique_ptr<BoundaryConditionBase>> bc_by_tag;
    std::unique_ptr<BoundaryConditionBase> default_bc;
    
    // Boundary faces organized by type
    std::map<BoundaryType, std::vector<BoundaryFace>> faces_by_type;
    std::vector<BoundaryFace> all_boundary_faces;
    
    // PML
    std::unique_ptr<PMLBoundary> pml;
    std::vector<int> pml_elements;
    
    // Face to BC mapping
    std::vector<BoundaryConditionBase*> face_bc;
};

/**
 * @brief Factory for creating boundary conditions from configuration
 */
class BoundaryConditionFactory {
public:
    /**
     * @brief Create boundary condition from type string
     */
    static std::unique_ptr<BoundaryConditionBase> create(const std::string& type);
    
    /**
     * @brief Create from configuration map
     */
    static std::unique_ptr<BoundaryConditionBase> createFromConfig(
        const std::map<std::string, std::string>& config,
        const std::string& prefix);
    
    /**
     * @brief Create free surface BC
     */
    static std::unique_ptr<FreeSurfaceBC> createFreeSurface();
    
    /**
     * @brief Create absorbing BC
     */
    static std::unique_ptr<BoundaryConditionBase> createAbsorbing(
        const std::string& method,
        double rho, double vp, double vs);
    
    /**
     * @brief Create PML BC
     */
    static std::unique_ptr<PMLBoundary> createPML(
        double thickness, double reflection_coef = 1e-6,
        double x_min = 0, double x_max = 0,
        double y_min = 0, double y_max = 0,
        double z_min = 0, double z_max = 0);
};

/**
 * @brief Configuration structure for boundary conditions
 */
struct BoundaryConditionConfig {
    // General settings
    std::map<int, std::string> bc_types;  // tag -> type name
    std::string default_bc_type;
    
    // Absorbing BC settings
    std::string absorbing_method;         // "pml", "clayton_engquist", "lysmer"
    int absorbing_order;                  // Order for CE
    
    // PML settings
    double pml_thickness;
    double pml_reflection_coef;
    int pml_damping_order;
    double pml_alpha_max;
    double pml_kappa_max;
    
    // Material properties (for absorbing BCs)
    double boundary_density;
    double boundary_vp;
    double boundary_vs;
    
    // Free surface settings
    bool use_characteristic_bc;
    
    /**
     * @brief Parse from configuration map
     */
    void parseConfig(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Setup boundary condition manager from config
     */
    void setupBoundaryConditions(BoundaryConditionManager& bcm) const;
};

} // namespace FSRM

#endif // BOUNDARY_CONDITIONS_HPP
