#ifndef POROELASTIC_SOLVER_HPP
#define POROELASTIC_SOLVER_HPP

#include "core/FSRM.hpp"
#include <petscts.h>
#include <petscdmplex.h>
#include <petscfe.h>
#include <vector>
#include <string>

namespace FSRM {

// Simple vector types
struct Vec3 { double x, y, z; };
struct Vec3i { int x, y, z; };

/**
 * @brief Poroelastic solver using PETSc TS with analytical Jacobians
 * 
 * Solves coupled fluid flow and geomechanics using:
 * - DMPlex for unstructured grid support (3D)
 * - PETSc Time Stepper (TS) for implicit time integration
 * - Hand-coded analytical Jacobian (no automatic differentiation)
 * - Fully implicit multiphase (water + CO2) pressure-saturation coupling
 * 
 * Fields (6 DOF per cell):
 *   0: Pressure (P)
 *   1: Water saturation (Sw)
 *   2: X-displacement (ux)
 *   3: Y-displacement (uy)
 *   4: Z-displacement (uz)
 *   5: Porosity (phi)
 */
class PoroelasticSolver {
public:
    struct WellData {
        std::string name;
        Vec3 position;
        bool is_injector;
        double target_rate;      // m³/day
        double target_bhp;       // Pa
        double current_bhp;      // Pa
        double oil_rate;         // m³/day
        double water_rate;       // m³/day
        double cumulative_oil;   // m³
        double cumulative_water; // m³
    };
    
    struct PhysicsParams {
        // Rock properties
        double porosity0 = 0.2;
        double permeability0 = 100e-15;   // m² (isotropic reference)
        double permeability_x = 100e-15;  // m² (x-direction)
        double permeability_y = 100e-15;  // m² (y-direction)
        double permeability_z = 10e-15;   // m² (z-direction, typically lower)
        double youngs_modulus = 10e9;     // Pa
        double poisson_ratio = 0.25;
        double biot_coefficient = 0.7;
        
        // Fluid properties
        double fluid_density = 1000.0;    // kg/m³
        double water_viscosity = 1e-3;    // Pa·s
        double oil_viscosity = 2e-3;      // Pa·s
        
        // Compressibilities
        double rock_compressibility = 1e-9;   // Pa⁻¹
        double fluid_compressibility = 1e-9;  // Pa⁻¹ (total compressibility)
        double water_compressibility = 4.5e-10; // Pa⁻¹
        double oil_compressibility = 1.5e-9;    // Pa⁻¹
        
        // Residual saturations for relative permeability
        double water_residual_saturation = 0.2;   // Swc/Swr
        double oil_residual_saturation = 0.2;     // Sor
        double gas_residual_saturation = 0.05;    // Sgc
        
        // Corey relative permeability model
        double corey_exponent = 4.0;              // Corey exponent for rel perm
        double corey_exponent_water = 4.0;        // Corey exponent for water
        double corey_exponent_oil = 2.0;          // Corey exponent for oil
        double kr_max_water = 1.0;                // Max water rel perm at Swr
        double kr_max_oil = 1.0;                  // Max oil rel perm at Sor
        
        // Initial conditions
        double initial_pressure = 30e6;   // Pa
        double initial_saturation = 0.2;  // Water saturation
    };
    
    PoroelasticSolver(MPI_Comm comm);
    ~PoroelasticSolver();
    
    // Setup
    void setDomain(Vec3 dims, Vec3i cells);
    void setPhysicsParams(const PhysicsParams& params);
    void addWell(const WellData& well);
    void initialize();
    
    // Solve
    void solve(double final_time, int num_steps);
    
    // Access results
    DM getDM() const { return dm_; }
    Vec getSolution() const { return solution_; }
    const std::vector<WellData>& getWells() const { return wells_; }
    
    // Extract specific fields (cell-centered values)
    void getPressure(std::vector<std::vector<double>>& P) const;
    void getSaturation(std::vector<std::vector<double>>& Sw) const;
    void getDisplacement(std::vector<std::vector<double>>& ux,
                        std::vector<std::vector<double>>& uz) const;
    void getDisplacement3D(std::vector<std::vector<double>>& ux,
                           std::vector<std::vector<double>>& uy,
                           std::vector<std::vector<double>>& uz) const;
    void getPorosity(std::vector<std::vector<double>>& phi) const;
    
private:
    MPI_Comm comm_;
    DM dm_;
    TS ts_;
    Vec solution_;
    PetscDS prob_;
    
    Vec3 domain_size_;
    Vec3i grid_cells_;
    PhysicsParams params_;
    std::vector<WellData> wells_;
    
    // DMPlex setup
    PetscErrorCode setupDMPlex();
    PetscErrorCode setupFields();
    PetscErrorCode setupTS();
    
    // Physics residual and Jacobian callbacks
    static PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void* ctx);
    static PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec X, Vec Xdot,
                                       PetscReal shift, Mat J, Mat P, void* ctx);
    
    // Pointwise residual (called for each cell)
    static void f0_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[],
                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[],
                           const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants,
                           const PetscScalar constants[], PetscScalar f0[]);
    
    static void f1_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                           const PetscInt uOff[], const PetscInt uOff_x[],
                           const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                           const PetscInt aOff[], const PetscInt aOff_x[],
                           const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                           PetscReal t, const PetscReal x[], PetscInt numConstants,
                           const PetscScalar constants[], PetscScalar f1[]);
    
    // Pointwise Jacobian (analytical derivatives)
    static void g0_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]);
    
    static void g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    // Saturation equation residual (field 1)
    static void f0_saturation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                             const PetscInt uOff[], const PetscInt uOff_x[],
                             const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                             const PetscInt aOff[], const PetscInt aOff_x[],
                             const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                             PetscReal t, const PetscReal x[], PetscInt numConstants,
                             const PetscScalar constants[], PetscScalar f0[]);
    
    static void f1_saturation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                             const PetscInt uOff[], const PetscInt uOff_x[],
                             const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                             const PetscInt aOff[], const PetscInt aOff_x[],
                             const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                             PetscReal t, const PetscReal x[], PetscInt numConstants,
                             const PetscScalar constants[], PetscScalar f1[]);
    
    // =========================================================================
    // Displacement equation residual and Jacobian (fields 2, 3, 4: ux, uy, uz)
    // Momentum equation: div(sigma) + alpha*grad(P) = 0  (quasi-static)
    // sigma_ij = lambda*eps_kk*delta_ij + 2*mu*eps_ij  (Hooke's law)
    // =========================================================================
    
    // f0 for displacement: body force + Biot coupling = 0  (no body force)
    // For quasi-static: f0_u = 0 (inertia term only if poroelastodynamics)
    static void f0_displacement(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                               const PetscInt uOff[], const PetscInt uOff_x[],
                               const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                               const PetscInt aOff[], const PetscInt aOff_x[],
                               const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                               PetscReal t, const PetscReal x[], PetscInt numConstants,
                               const PetscScalar constants[], PetscScalar f0[]);
    
    // f1 for displacement: stress tensor + Biot pressure coupling
    // f1[i*dim+j] = sigma_ij + alpha*P*delta_ij
    // where sigma_ij = lambda*eps_kk*delta_ij + 2*mu*eps_ij
    // eps_ij = 0.5*(du_i/dx_j + du_j/dx_i)
    // Note: For separate scalar fields (ux=field2, uy=field3, uz=field4),
    // the strain must be assembled from gradients of multiple fields.
    // PETSc calls f1 for each field separately, so field 2 (ux) gets
    // f1[j] = sigma_{0j} + alpha*P*delta_{0j}
    static void f1_displacement_x(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                  const PetscInt uOff[], const PetscInt uOff_x[],
                                  const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                  const PetscInt aOff[], const PetscInt aOff_x[],
                                  const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                  PetscReal t, const PetscReal x[], PetscInt numConstants,
                                  const PetscScalar constants[], PetscScalar f1[]);
    
    static void f1_displacement_y(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                  const PetscInt uOff[], const PetscInt uOff_x[],
                                  const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                  const PetscInt aOff[], const PetscInt aOff_x[],
                                  const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                  PetscReal t, const PetscReal x[], PetscInt numConstants,
                                  const PetscScalar constants[], PetscScalar f1[]);
    
    static void f1_displacement_z(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                  const PetscInt uOff[], const PetscInt uOff_x[],
                                  const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                  const PetscInt aOff[], const PetscInt aOff_x[],
                                  const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                  PetscReal t, const PetscReal x[], PetscInt numConstants,
                                  const PetscScalar constants[], PetscScalar f1[]);
    
    // Jacobian: g3 for (displacement_i, displacement_j) = elasticity stiffness
    // d(sigma_{ij})/d(du_k/dx_l) = C_{ijkl} = lambda*delta_{ij}*delta_{kl} + mu*(delta_{ik}*delta_{jl}+delta_{il}*delta_{jk})
    // Since fields are separate scalars, we need g3 for each (field_i, field_j) pair.
    // g3_uxux: d(f1_ux)/d(grad_ux) — row 0, col 0 of elasticity
    static void g3_uxux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    // Generic elasticity Jacobian for any (displacement_i, displacement_j) pair
    // The actual i,j are encoded via constants or the callback registration
    static void g3_uxuy(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uxuz(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uyux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uyuy(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uyuz(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uzux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uzuy(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g3_uzuz(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    // Biot coupling Jacobians
    // g1_ux_p: d(f1_ux)/d(P) = alpha * delta_{0j} (Biot coupling in momentum)
    // Actually, Biot coupling P→u goes via the stress: sigma' = sigma + alpha*P*I
    // So d(f1_ui)/d(P) = alpha * delta_{ij} — this is a g2 term (d(f1)/d(u))
    // but since P is field 0, we use g2 for (field_ui, field_P)
    static void g2_ux_p(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g2[]);
    
    static void g2_uy_p(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    static void g2_uz_p(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    // Biot coupling from displacement to pressure equation
    // d(f0_p)/d(u_dot) contributes alpha * shift through displacement divergence
    // This is already partially in g0_pp[2] and g0_pp[3] but needs proper g1 terms
    // g1_p_ux: d(f0_p)/d(grad_ux) = alpha * shift (x-component of div(u_dot))
    static void g1_p_ux(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[]);
    
    static void g1_p_uy(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[]);
    
    static void g1_p_uz(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[],
                       const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[],
                       const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, PetscReal u_tShift, const PetscReal x[],
                       PetscInt numConstants, const PetscScalar constants[], PetscScalar g1[]);
    
    // Cross-coupling Jacobians for pressure-saturation system
    static void g0_ps(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]);
    
    static void g0_ss(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]);
    
    static void g3_sp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]);
    
    // Boundary conditions
    void applyBoundaryConditions();
    
    // Well handling
    void updateWellPerformance(double time, double dt);
    
    // Helper to get cell index from coordinates
    PetscInt getCellIndex(PetscInt i, PetscInt k) const {
        return k * grid_cells_.x + i;
    }
};

} // namespace FSRM

#endif // POROELASTIC_SOLVER_HPP
