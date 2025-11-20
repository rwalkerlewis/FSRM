#include "PoroelasticSolver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ResSim {

PoroelasticSolver::PoroelasticSolver(MPI_Comm comm) 
    : comm_(comm), dm_(nullptr), ts_(nullptr), solution_(nullptr), prob_(nullptr) {}

PoroelasticSolver::~PoroelasticSolver() {
    if (solution_) VecDestroy(&solution_);
    if (ts_) TSDestroy(&ts_);
    if (dm_) DMDestroy(&dm_);
}

void PoroelasticSolver::setDomain(Vec3 dims, Vec3i cells) {
    domain_size_ = dims;
    grid_cells_ = cells;
}

void PoroelasticSolver::setPhysicsParams(const PhysicsParams& params) {
    params_ = params;
}

void PoroelasticSolver::addWell(const WellData& well) {
    wells_.push_back(well);
}

PetscErrorCode PoroelasticSolver::setupDMPlex() {
    PetscErrorCode ierr;
    
    PetscFunctionBeginUser;
    
    // Create structured box mesh using DMPlex
    PetscInt faces[2] = {grid_cells_.x, grid_cells_.z};
    PetscReal lower[2] = {0.0, 0.0};
    PetscReal upper[2] = {domain_size_.x, domain_size_.z};
    
    ierr = DMPlexCreateBoxMesh(comm_, 2, PETSC_FALSE, faces, lower, upper,
                               nullptr, PETSC_TRUE, &dm_); CHKERRQ(ierr);
    
    // Set name for visualization
    ierr = PetscObjectSetName((PetscObject)dm_, "Poroelastic Mesh"); CHKERRQ(ierr);
    
    ierr = DMSetFromOptions(dm_); CHKERRQ(ierr);
    ierr = DMViewFromOptions(dm_, nullptr, "-dm_view"); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode PoroelasticSolver::setupFields() {
    PetscErrorCode ierr;
    PetscFE fe[5];
    
    PetscFunctionBeginUser;
    
    // Create 5 finite element fields (all P1 continuous)
    for (int i = 0; i < 5; ++i) {
        ierr = PetscFECreateDefault(comm_, 2, 1, PETSC_FALSE, nullptr, -1, &fe[i]); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe[i], 
                                 i == 0 ? "pressure" :
                                 i == 1 ? "saturation" :
                                 i == 2 ? "ux" :
                                 i == 3 ? "uz" : "porosity"); CHKERRQ(ierr);
    }
    
    // Add fields to DM
    for (int i = 0; i < 5; ++i) {
        ierr = DMAddField(dm_, nullptr, (PetscObject)fe[i]); CHKERRQ(ierr);
    }
    
    ierr = DMCreateDS(dm_); CHKERRQ(ierr);
    ierr = DMGetDS(dm_, &prob_); CHKERRQ(ierr);
    
    // Set up residual forms (f0 = time derivative, f1 = flux)
    ierr = PetscDSSetResidual(prob_, 0, f0_pressure, f1_pressure); CHKERRQ(ierr);
    
    // Set up analytical Jacobian forms
    // g0 = df0/du, g1 = df0/du_x, g2 = df1/du, g3 = df1/du_x
    ierr = PetscDSSetJacobian(prob_, 0, 0, g0_pp, nullptr, nullptr, g3_pp); CHKERRQ(ierr);
    
    // Cleanup FE objects
    for (int i = 0; i < 5; ++i) {
        ierr = PetscFEDestroy(&fe[i]); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

void PoroelasticSolver::initialize() {
    PetscErrorCode ierr;
    
    ierr = setupDMPlex(); CHKERRABORT(comm_, ierr);
    ierr = setupFields(); CHKERRABORT(comm_, ierr);
    
    ierr = DMCreateGlobalVector(dm_, &solution_); CHKERRABORT(comm_, ierr);
    ierr = PetscObjectSetName((PetscObject)solution_, "Solution"); CHKERRABORT(comm_, ierr);
    
    // Set initial conditions using cell-centered values
    Vec local;
    ierr = DMGetLocalVector(dm_, &local); CHKERRABORT(comm_, ierr);
    
    PetscScalar *array;
    ierr = VecGetArray(local, &array); CHKERRABORT(comm_, ierr);
    
    PetscInt cStart, cEnd;
    ierr = DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd); CHKERRABORT(comm_, ierr);
    
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscScalar centroid[2];
        PetscReal volume;
        ierr = DMPlexComputeCellGeometryFVM(dm_, c, &volume, centroid, nullptr); CHKERRABORT(comm_, ierr);
        
        double depth = centroid[1];
        
        // Get DOF offset for this cell
        PetscInt offset;
        ierr = DMPlexGetPointGlobal(dm_, c, &offset, nullptr); CHKERRABORT(comm_, ierr);
        
        if (offset >= 0) {
            // Hydrostatic pressure
            array[offset + 0] = params_.initial_pressure + params_.fluid_density * 9.81 * depth;
            array[offset + 1] = params_.initial_saturation;
            array[offset + 2] = 0.0;  // ux
            array[offset + 3] = 0.0;  // uz
            array[offset + 4] = params_.porosity0;
        }
    }
    
    ierr = VecRestoreArray(local, &array); CHKERRABORT(comm_, ierr);
    ierr = DMLocalToGlobal(dm_, local, INSERT_VALUES, solution_); CHKERRABORT(comm_, ierr);
    ierr = DMRestoreLocalVector(dm_, &local); CHKERRABORT(comm_, ierr);
}

PetscErrorCode PoroelasticSolver::setupTS() {
    PetscErrorCode ierr;
    
    PetscFunctionBeginUser;
    
    ierr = TSCreate(comm_, &ts_); CHKERRQ(ierr);
    ierr = TSSetDM(ts_, dm_); CHKERRQ(ierr);
    ierr = TSSetProblemType(ts_, TS_NONLINEAR); CHKERRQ(ierr);
    
    // Let DM provide the residual and Jacobian through the DS
    ierr = DMTSSetBoundaryLocal(dm_, nullptr, nullptr); CHKERRQ(ierr);
    ierr = DMTSSetIFunctionLocal(dm_, DMPlexTSComputeIFunctionFEM, nullptr); CHKERRQ(ierr);
    ierr = DMTSSetIJacobianLocal(dm_, DMPlexTSComputeIJacobianFEM, nullptr); CHKERRQ(ierr);
    
    ierr = TSSetType(ts_, TSBEULER); CHKERRQ(ierr);  // Backward Euler
    ierr = TSSetMaxSteps(ts_, 1000000); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts_, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    
    ierr = TSSetFromOptions(ts_); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

void PoroelasticSolver::solve(double final_time, int num_steps) {
    PetscErrorCode ierr;
    
    ierr = setupTS(); CHKERRABORT(comm_, ierr);
    
    double dt = final_time / num_steps;
    
    ierr = TSSetTimeStep(ts_, dt); CHKERRABORT(comm_, ierr);
    ierr = TSSetMaxTime(ts_, final_time); CHKERRABORT(comm_, ierr);
    
    ierr = TSSolve(ts_, solution_); CHKERRABORT(comm_, ierr);
    
    TSConvergedReason reason;
    ierr = TSGetConvergedReason(ts_, &reason); CHKERRABORT(comm_, ierr);
    
    PetscInt steps;
    ierr = TSGetStepNumber(ts_, &steps); CHKERRABORT(comm_, ierr);
    
    PetscReal ftime;
    ierr = TSGetTime(ts_, &ftime); CHKERRABORT(comm_, ierr);
    
    PetscPrintf(comm_, "Solver finished: %d steps, t=%.2e, reason=%s\n",
                steps, ftime, TSConvergedReasons[reason]);
}

// ============================================================================
// Pointwise residual functions for pressure equation
// ============================================================================

void PoroelasticSolver::f0_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[],
                                    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[],
                                    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], PetscInt numConstants,
                                    const PetscScalar constants[], PetscScalar f0[]) {
    // f0 = phi * ct * dP/dt
    // u[0] = P, u[4] = phi, u_t[0] = dP/dt
    
    const PetscScalar P     = u[uOff[0]];
    const PetscScalar Sw    = u[uOff[1]];
    const PetscScalar phi   = u[uOff[4]];
    const PetscScalar P_t   = u_t[uOff[0]];
    
    // Total compressibility (approximation)
    const PetscScalar ct = 1e-9;  // Should get from constants
    
    f0[0] = phi * ct * P_t;
}

void PoroelasticSolver::f1_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[],
                                    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[],
                                    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], PetscInt numConstants,
                                    const PetscScalar constants[], PetscScalar f1[]) {
    // f1 = -k/mu * grad(P)
    // u_x[0] = dP/dx, u_x[1] = dP/dz (for 2D)
    
    const PetscScalar P_x = u_x[uOff_x[0] + 0];
    const PetscScalar P_z = u_x[uOff_x[0] + 1];
    
    const PetscScalar k = 100e-15;   // Should get from constants
    const PetscScalar mu = 1e-3;
    const PetscScalar mob = k / mu;
    
    f1[0] = -mob * P_x;
    f1[1] = -mob * P_z;
}

// ============================================================================
// Analytical Jacobian functions
// ============================================================================

void PoroelasticSolver::g0_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]) {
    // g0 = df0/dP = phi * ct * shift (where shift = d()/dP_t)
    const PetscScalar phi = u[uOff[4]];
    const PetscScalar ct = 1e-9;
    
    g0[0] = phi * ct * u_tShift;
}

void PoroelasticSolver::g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]) {
    // g3 = df1/d(grad P) = -k/mu * I (identity tensor)
    const PetscScalar k = 100e-15;
    const PetscScalar mu = 1e-3;
    const PetscScalar mob = k / mu;
    
    // g3 is flattened 2x2 matrix: [dF_x/dP_x, dF_x/dP_z; dF_z/dP_x, dF_z/dP_z]
    g3[0] = -mob;  // dF_x/dP_x
    g3[1] = 0.0;   // dF_x/dP_z
    g3[2] = 0.0;   // dF_z/dP_x
    g3[3] = -mob;  // dF_z/dP_z
}

// ============================================================================
// Result extraction
// ============================================================================

void PoroelasticSolver::getPressure(std::vector<std::vector<double>>& P) const {
    P.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
    
    Vec local;
    DMGetLocalVector(dm_, &local);
    DMGlobalToLocal(dm_, solution_, INSERT_VALUES, local);
    
    PetscSection section;
    DMGetLocalSection(dm_, &section);
    
    PetscScalar *array;
    VecGetArray(local, &array);
    
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt dof, offset;
        PetscSectionGetDof(section, c, &dof);
        PetscSectionGetOffset(section, c, &offset);
        
        if (dof >= 5) {
            // Determine (i, k) from cell index
            PetscInt cell_idx = c - cStart;
            PetscInt i = cell_idx % grid_cells_.x;
            PetscInt k = cell_idx / grid_cells_.x;
            
            if (k < grid_cells_.z && i < grid_cells_.x) {
                P[k][i] = array[offset + 0];
            }
        }
    }
    
    VecRestoreArray(local, &array);
    DMRestoreLocalVector(dm_, &local);
}

void PoroelasticSolver::getSaturation(std::vector<std::vector<double>>& Sw) const {
    Sw.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
    
    Vec local;
    DMGetLocalVector(dm_, &local);
    DMGlobalToLocal(dm_, solution_, INSERT_VALUES, local);
    
    PetscSection section;
    DMGetLocalSection(dm_, &section);
    
    PetscScalar *array;
    VecGetArray(local, &array);
    
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt dof, offset;
        PetscSectionGetDof(section, c, &dof);
        PetscSectionGetOffset(section, c, &offset);
        
        if (dof >= 5) {
            PetscInt cell_idx = c - cStart;
            PetscInt i = cell_idx % grid_cells_.x;
            PetscInt k = cell_idx / grid_cells_.x;
            
            if (k < grid_cells_.z && i < grid_cells_.x) {
                Sw[k][i] = array[offset + 1];
            }
        }
    }
    
    VecRestoreArray(local, &array);
    DMRestoreLocalVector(dm_, &local);
}

void PoroelasticSolver::getDisplacement(std::vector<std::vector<double>>& ux,
                                       std::vector<std::vector<double>>& uz) const {
    ux.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
    uz.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
    
    Vec local;
    DMGetLocalVector(dm_, &local);
    DMGlobalToLocal(dm_, solution_, INSERT_VALUES, local);
    
    PetscSection section;
    DMGetLocalSection(dm_, &section);
    
    PetscScalar *array;
    VecGetArray(local, &array);
    
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt dof, offset;
        PetscSectionGetDof(section, c, &dof);
        PetscSectionGetOffset(section, c, &offset);
        
        if (dof >= 5) {
            PetscInt cell_idx = c - cStart;
            PetscInt i = cell_idx % grid_cells_.x;
            PetscInt k = cell_idx / grid_cells_.x;
            
            if (k < grid_cells_.z && i < grid_cells_.x) {
                ux[k][i] = array[offset + 2];
                uz[k][i] = array[offset + 3];
            }
        }
    }
    
    VecRestoreArray(local, &array);
    DMRestoreLocalVector(dm_, &local);
}void PoroelasticSolver::getPorosity(std::vector<std::vector<double>>& phi) const {
    phi.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
    
    Vec local;
    DMGetLocalVector(dm_, &local);
    DMGlobalToLocal(dm_, solution_, INSERT_VALUES, local);
    
    PetscSection section;
    DMGetLocalSection(dm_, &section);
    
    PetscScalar *array;
    VecGetArray(local, &array);
    
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    
    for (PetscInt c = cStart; c < cEnd; ++c) {
        PetscInt dof, offset;
        PetscSectionGetDof(section, c, &dof);
        PetscSectionGetOffset(section, c, &offset);
        
        if (dof >= 5) {
            PetscInt cell_idx = c - cStart;
            PetscInt i = cell_idx % grid_cells_.x;
            PetscInt k = cell_idx / grid_cells_.x;
            
            if (k < grid_cells_.z && i < grid_cells_.x) {
                phi[k][i] = array[offset + 4];
            }
        }
    }
    
    VecRestoreArray(local, &array);
    DMRestoreLocalVector(dm_, &local);
}

void PoroelasticSolver::applyBoundaryConditions() {
    // TODO: Implement boundary conditions using DMPlex labels
}

void PoroelasticSolver::updateWellPerformance(double time, double dt) {
    // TODO: Implement well source terms
}

PetscErrorCode PoroelasticSolver::FormIFunction(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void* ctx) {
    // DMPlex uses pointwise functions, this is just a wrapper
    return 0;
}

PetscErrorCode PoroelasticSolver::FormIJacobian(TS ts, PetscReal t, Vec X, Vec Xdot,
                                                PetscReal shift, Mat J, Mat P, void* ctx) {
    // DMPlex uses pointwise functions, this is just a wrapper
    return 0;
}

} // namespace ResSim
