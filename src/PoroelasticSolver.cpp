#include "PoroelasticSolver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace FSRM {

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
    
    // Set up constants array with user-configurable physics parameters
    // Constants layout:
    // [0] = total compressibility
    // [1] = Biot coefficient
    // [2] = permeability_x
    // [3] = permeability_z
    // [4] = fluid viscosity
    // [5] = water compressibility
    // [6] = oil compressibility
    // [7] = water residual saturation (Swr)
    // [8] = oil residual saturation (Sor)
    // [9] = Corey exponent
    // [10] = max relative permeability
    PetscScalar constants[11];
    constants[0] = params_.fluid_compressibility;
    constants[1] = params_.biot_coefficient;
    constants[2] = params_.permeability_x;
    constants[3] = params_.permeability_z;
    constants[4] = params_.water_viscosity;
    constants[5] = params_.water_compressibility;
    constants[6] = params_.oil_compressibility;
    constants[7] = params_.water_residual_saturation;
    constants[8] = params_.oil_residual_saturation;
    constants[9] = params_.corey_exponent;
    constants[10] = params_.kr_max_water;
    
    ierr = PetscDSSetConstants(prob_, 11, constants); CHKERRQ(ierr);
    
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
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    // f0 = phi * ct * dP/dt
    // u[0] = P, u[4] = phi, u_t[0] = dP/dt
    
    const PetscScalar P     = u[uOff[0]];
    const PetscScalar Sw    = u[uOff[1]];
    const PetscScalar phi   = u[uOff[4]];
    const PetscScalar P_t   = u_t[uOff[0]];
    (void)P; (void)Sw;  // Used in full implementation
    
    // Total compressibility from user-configurable constants
    // Constants layout: [0] = total compressibility
    const PetscScalar ct = (numConstants > 0) ? constants[0] : 1e-9;
    
    f0[0] = phi * ct * P_t;
}

void PoroelasticSolver::f1_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                    const PetscInt uOff[], const PetscInt uOff_x[],
                                    const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                    const PetscInt aOff[], const PetscInt aOff_x[],
                                    const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                    PetscReal t, const PetscReal x[], PetscInt numConstants,
                                    const PetscScalar constants[], PetscScalar f1[]) {
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)u; (void)u_t;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    // f1 = -k/mu * grad(P)
    // u_x[0] = dP/dx, u_x[1] = dP/dz (for 2D)
    
    const PetscScalar P_x = u_x[uOff_x[0] + 0];
    const PetscScalar P_z = u_x[uOff_x[0] + 1];
    
    // Get permeability and viscosity from user-configurable constants
    // Constants layout: [2] = permeability_x, [3] = permeability_z, [4] = viscosity
    const PetscScalar kx = (numConstants > 2) ? constants[2] : 100e-15;
    const PetscScalar kz = (numConstants > 3) ? constants[3] : 10e-15;
    const PetscScalar mu = (numConstants > 4) ? constants[4] : 1e-3;
    
    const PetscScalar mob_x = kx / mu;
    const PetscScalar mob_z = kz / mu;
    
    f1[0] = -mob_x * P_x;
    f1[1] = -mob_z * P_z;
}

// ============================================================================
// Analytical Jacobian functions (Full implementations)
// ============================================================================

void PoroelasticSolver::g0_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]) {
    // Full Jacobian g0 = ∂f0/∂u for pressure equation
    //
    // The pressure equation in coupled poroelasticity:
    // ∂/∂t(φ*ct*P + α*∇·u) + ∇·q = Q
    //
    // where:
    // - φ = porosity (field 4)
    // - ct = total compressibility
    // - P = pressure (field 0)
    // - α = Biot coefficient
    // - u = displacement (fields 2,3)
    // - q = -k/μ * kr * ∇P (Darcy flux)
    //
    // g0 Jacobian contributions:
    // ∂f0/∂P: φ * ct * shift (from storage term)
    // ∂f0/∂S: φ * ∂ct/∂S * shift + ∂kr/∂S contribution (from saturation-dependent properties)
    // ∂f0/∂u: α * div(δu) * shift (Biot coupling to displacement divergence)
    // ∂f0/∂φ: ct * P * shift (from porosity variation)
    
    (void)dim; (void)Nf; (void)NfAux;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    // Extract solution fields
    const PetscScalar P = u[uOff[0]];
    const PetscScalar S = u[uOff[1]];
    const PetscScalar phi = u[uOff[4]];
    
    // Material parameters from user-configurable constants
    // Constants layout:
    // [0] = total compressibility, [1] = Biot coefficient
    // [2] = permeability_x, [3] = permeability_z, [4] = viscosity
    // [5] = water compressibility, [6] = oil compressibility
    // [7] = water residual saturation (Swr), [8] = oil residual saturation (Sor)
    // [9] = Corey exponent, [10] = max relative permeability
    PetscScalar ct = (numConstants > 0) ? constants[0] : 1e-9;      // Total compressibility [1/Pa]
    PetscScalar alpha = (numConstants > 1) ? constants[1] : 1.0;    // Biot coefficient
    
    // Saturation-dependent compressibility (if two-phase)
    // ct_eff = S*cw + (1-S)*co where cw, co are phase compressibilities
    PetscScalar cw = (numConstants > 5) ? constants[5] : 4.5e-10;  // Water compressibility
    PetscScalar co = (numConstants > 6) ? constants[6] : 1.5e-9;   // Oil compressibility
    PetscScalar ct_eff = S * cw + (1.0 - S) * co;
    
    // ∂ct/∂S for saturation Jacobian contribution
    PetscScalar dct_dS = cw - co;
    
    // Initialize g0 to zero (Nf x Nf matrix, but we only fill pressure row)
    // For full coupled system, g0 would be 5x5
    
    // ∂f0/∂P: storage term Jacobian
    g0[0] = phi * ct_eff * u_tShift;
    
    // ∂f0/∂S: saturation coupling (through compressibility and relative permeability)
    if (Nf > 1) {
        // Get relative permeability parameters from constants
        PetscScalar Swr = (numConstants > 7) ? constants[7] : 0.2;   // Water residual saturation
        PetscScalar Sor = (numConstants > 8) ? constants[8] : 0.2;   // Oil residual saturation
        PetscScalar n_corey = (numConstants > 9) ? constants[9] : 4.0;  // Corey exponent
        PetscScalar kx = (numConstants > 2) ? constants[2] : 100e-15;   // Permeability
        PetscScalar mu = (numConstants > 4) ? constants[4] : 1e-3;      // Viscosity
        
        // Compressibility change with saturation
        g0[1] = phi * dct_dS * P * u_tShift;
        
        // Relative permeability derivative effect on flux
        // d(kr)/dS affects the transmissibility
        PetscScalar S_denom = 1.0 - Swr - Sor;
        PetscScalar Sw_norm = std::max(0.0, std::min(1.0, (S - Swr) / S_denom));
        PetscScalar dkr_dS = (Sw_norm > 0.0 && Sw_norm < 1.0) ? 
                             n_corey * std::pow(Sw_norm, n_corey - 1.0) / S_denom : 0.0;
        
        // Flux contribution (would be multiplied by pressure gradient magnitude)
        const PetscScalar dPdx = u_x[uOff_x[0]];
        const PetscScalar dPdz = u_x[uOff_x[0] + 1];
        PetscScalar grad_P_mag = std::sqrt(dPdx*dPdx + dPdz*dPdz);
        
        g0[1] += kx / mu * dkr_dS * grad_P_mag;
    }
    
    // ∂f0/∂ux, ∂f0/∂uz: Biot coupling (displacement divergence affects storage)
    // From storage term: ∂/∂t(α * ∇·u) → α * ∂(∇·u)/∂t
    // This contributes through the gradient terms, handled in g1
    if (Nf > 2) {
        g0[2] = alpha * u_tShift;  // Coupling to ∂ux/∂x (through divergence)
        g0[3] = alpha * u_tShift;  // Coupling to ∂uz/∂z (through divergence)
    }
    
    // ∂f0/∂φ: porosity variation effect
    if (Nf > 4) {
        g0[4] = ct_eff * P * u_tShift;
    }
    
    // Suppress remaining unused parameters
    (void)u_t;
}

void PoroelasticSolver::g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]) {
    // Full Jacobian g3 = ∂f1/∂(∇u) for pressure equation flux term
    //
    // The flux term in pressure equation:
    // f1 = q = -k/μ * kr(S) * ∇P (Darcy flux)
    //
    // g3 Jacobian contributions:
    // ∂q/∂(∇P) = -k/μ * kr * I (isotropic permeability)
    // For anisotropic k: ∂q_i/∂(∂P/∂x_j) = -k_ij/μ * kr
    //
    // Layout: g3[d1*dim + d2] = ∂(flux_d1)/∂(grad_d2 P)
    
    (void)Nf; (void)NfAux;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;
    (void)u_t; (void)u_x;
    
    // Extract saturation for relative permeability
    const PetscScalar S = u[uOff[1]];
    
    // Get all parameters from user-configurable constants
    // Constants layout:
    // [0] = total compressibility, [1] = Biot coefficient
    // [2] = permeability_x, [3] = permeability_z, [4] = viscosity
    // [5] = water compressibility, [6] = oil compressibility
    // [7] = water residual saturation (Swr), [8] = oil residual saturation (Sor)
    // [9] = Corey exponent, [10] = max relative permeability
    PetscScalar kx = (numConstants > 2) ? constants[2] : 100e-15;  // Permeability x [m²]
    PetscScalar kz = (numConstants > 3) ? constants[3] : 10e-15;   // Permeability z [m²]
    PetscScalar mu = (numConstants > 4) ? constants[4] : 1e-3;     // Viscosity [Pa·s]
    
    // Relative permeability (Corey model) with user-configurable parameters
    // kr = krmax * ((S - Swr)/(1 - Swr - Sor))^n
    PetscScalar Swr = (numConstants > 7) ? constants[7] : 0.2;     // Water residual saturation
    PetscScalar Sor = (numConstants > 8) ? constants[8] : 0.2;     // Oil residual saturation
    PetscScalar n_corey = (numConstants > 9) ? constants[9] : 4.0; // Corey exponent
    PetscScalar krmax = (numConstants > 10) ? constants[10] : 1.0; // Max relative permeability
    
    PetscScalar S_norm = std::max(0.0, std::min(1.0, (S - Swr) / (1.0 - Swr - Sor)));
    PetscScalar kr = krmax * std::pow(S_norm, n_corey);
    
    // Mobilities
    PetscScalar mob_x = kx * kr / mu;
    PetscScalar mob_z = kz * kr / mu;
    
    // Fill g3 tensor (dim × dim matrix)
    // For 2D: g3[0] = ∂q_x/∂(∂P/∂x), g3[1] = ∂q_x/∂(∂P/∂z)
    //         g3[2] = ∂q_z/∂(∂P/∂x), g3[3] = ∂q_z/∂(∂P/∂z)
    // For 3D: 3×3 matrix
    
    if (dim == 2) {
        g3[0] = mob_x;   // ∂q_x/∂(∂P/∂x) - main diagonal
        g3[1] = 0.0;     // ∂q_x/∂(∂P/∂z) - off-diagonal (zero for diagonal tensor)
        g3[2] = 0.0;     // ∂q_z/∂(∂P/∂x) - off-diagonal
        g3[3] = mob_z;   // ∂q_z/∂(∂P/∂z) - main diagonal
    } else if (dim == 3) {
        PetscScalar ky = kx;  // Assume isotropic in x-y
        PetscScalar mob_y = ky * kr / mu;
        
        g3[0] = mob_x;  g3[1] = 0.0;    g3[2] = 0.0;     // Row 0
        g3[3] = 0.0;    g3[4] = mob_y;  g3[5] = 0.0;     // Row 1
        g3[6] = 0.0;    g3[7] = 0.0;    g3[8] = mob_z;   // Row 2
    }
    
    // For coupled poroelasticity, would also include:
    // ∂q/∂(∇u) terms if permeability depends on strain (stress-dependent permeability)
    // k = k0 * exp(α_k * (σ_eff - σ_ref)) where σ_eff = σ_total - α*p
    // This would add off-diagonal blocks to the full coupled Jacobian
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

} // namespace FSRM
