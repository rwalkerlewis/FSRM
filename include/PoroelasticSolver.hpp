#ifndef POROELASTIC_SOLVER_HPP
#define POROELASTIC_SOLVER_HPP

#include "ReservoirSim.hpp"
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
 * - DMPlex for unstructured grid support
 * - PETSc Time Stepper (TS) for implicit time integration
 * - Hand-coded analytical Jacobian (no automatic differentiation)
 * 
 * Fields (5 DOF per cell):
 *   0: Pressure (P)
 *   1: Water saturation (Sw)
 *   2: X-displacement (ux)
 *   3: Z-displacement (uz)
 *   4: Porosity (phi)
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
