#include "numerics/PoroelasticSolver.hpp"
#include "physics/MultiphaseFlowKernel.hpp"
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
    
    // Create 3D structured box mesh using DMPlex
    PetscInt faces[3] = {grid_cells_.x, grid_cells_.y, grid_cells_.z};
    PetscReal lower[3] = {0.0, 0.0, 0.0};
    PetscReal upper[3] = {domain_size_.x, domain_size_.y, domain_size_.z};
    
    ierr = DMPlexCreateBoxMesh(comm_, 3, PETSC_FALSE, faces, lower, upper,
                               nullptr, PETSC_TRUE, &dm_); CHKERRQ(ierr);
    
    // Set name for visualization
    ierr = PetscObjectSetName((PetscObject)dm_, "Poroelastic Mesh"); CHKERRQ(ierr);
    
    ierr = DMSetFromOptions(dm_); CHKERRQ(ierr);
    ierr = DMViewFromOptions(dm_, nullptr, "-dm_view"); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode PoroelasticSolver::setupFields() {
    PetscErrorCode ierr;
    PetscFE fe[6];
    
    PetscFunctionBeginUser;
    
    // Create 6 finite element fields (all P1 continuous, 3D)
    // Fields: 0=pressure, 1=saturation, 2=ux, 3=uy, 4=uz, 5=porosity
    const char* field_names[6] = {
        "pressure", "saturation", "ux", "uy", "uz", "porosity"
    };
    for (int i = 0; i < 6; ++i) {
        ierr = PetscFECreateDefault(comm_, 3, 1, PETSC_FALSE, nullptr, -1, &fe[i]); CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)fe[i], field_names[i]); CHKERRQ(ierr);
    }
    
    // Add fields to DM
    for (int i = 0; i < 6; ++i) {
        ierr = DMAddField(dm_, nullptr, (PetscObject)fe[i]); CHKERRQ(ierr);
    }
    
    ierr = DMCreateDS(dm_); CHKERRQ(ierr);
    ierr = DMGetDS(dm_, &prob_); CHKERRQ(ierr);
    
    // Set up constants array with physics parameters for multiphase flow
    // Layout matches MultiphaseFlowKernel constants:
    // [0]  = porosity
    // [1]  = permeability_x       [m²]
    // [2]  = permeability_y       [m²]
    // [3]  = permeability_z       [m²]
    // [4]  = water_viscosity      [Pa·s]
    // [5]  = co2/oil viscosity    [Pa·s]
    // [6]  = water_density        [kg/m³]
    // [7]  = co2/oil density      [kg/m³]
    // [8]  = water_compressibility [1/Pa]
    // [9]  = co2/oil compressibility [1/Pa]
    // [10] = Swr (water residual saturation)
    // [11] = Sor/Sgr (oil/CO2 residual saturation)
    // [12] = Corey exponent water
    // [13] = Corey exponent oil/CO2
    // [14] = max water relative permeability
    // [15] = max oil/CO2 relative permeability
    // [16] = gravity [m/s²]
    // [17] = Biot coefficient
    PetscScalar constants[18];
    constants[0]  = params_.porosity0;
    constants[1]  = params_.permeability_x;
    constants[2]  = params_.permeability_y;
    constants[3]  = params_.permeability_z;
    constants[4]  = params_.water_viscosity;
    constants[5]  = params_.oil_viscosity;
    constants[6]  = params_.fluid_density;
    constants[7]  = 700.0;  // CO2/oil density (default)
    constants[8]  = params_.water_compressibility;
    constants[9]  = params_.oil_compressibility;
    constants[10] = params_.water_residual_saturation;
    constants[11] = params_.oil_residual_saturation;
    constants[12] = params_.corey_exponent_water;
    constants[13] = params_.corey_exponent_oil;
    constants[14] = params_.kr_max_water;
    constants[15] = params_.kr_max_oil;
    constants[16] = 9.81;
    constants[17] = params_.biot_coefficient;
    
    ierr = PetscDSSetConstants(prob_, 18, constants); CHKERRQ(ierr);
    
    // =========================================================================
    // Pressure equation (field 0): residual and Jacobian
    // =========================================================================
    ierr = PetscDSSetResidual(prob_, 0, f0_pressure, f1_pressure); CHKERRQ(ierr);
    
    // (pressure, pressure) block
    ierr = PetscDSSetJacobian(prob_, 0, 0, g0_pp, nullptr, nullptr, g3_pp); CHKERRQ(ierr);
    
    // (pressure, saturation) cross-coupling block
    ierr = PetscDSSetJacobian(prob_, 0, 1, g0_ps, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    
    // =========================================================================
    // Saturation equation (field 1): residual and Jacobian
    // =========================================================================
    ierr = PetscDSSetResidual(prob_, 1, f0_saturation, f1_saturation); CHKERRQ(ierr);
    
    // (saturation, saturation) block
    ierr = PetscDSSetJacobian(prob_, 1, 1, g0_ss, nullptr, nullptr, nullptr); CHKERRQ(ierr);
    
    // (saturation, pressure) cross-coupling block — water mobility in gradient
    ierr = PetscDSSetJacobian(prob_, 1, 0, nullptr, nullptr, nullptr, g3_sp); CHKERRQ(ierr);
    
    // Cleanup FE objects
    for (int i = 0; i < 6; ++i) {
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
        PetscScalar centroid[3];
        PetscReal volume;
        ierr = DMPlexComputeCellGeometryFVM(dm_, c, &volume, centroid, nullptr); CHKERRABORT(comm_, ierr);
        
        double depth = centroid[2];  // z-coordinate is depth in 3D
        
        // Get DOF offset for this cell
        PetscInt offset;
        ierr = DMPlexGetPointGlobal(dm_, c, &offset, nullptr); CHKERRABORT(comm_, ierr);
        
        if (offset >= 0) {
            // Hydrostatic pressure
            array[offset + 0] = params_.initial_pressure + params_.fluid_density * 9.81 * depth;
            array[offset + 1] = params_.initial_saturation;
            array[offset + 2] = 0.0;  // ux
            array[offset + 3] = 0.0;  // uy
            array[offset + 4] = 0.0;  // uz
            array[offset + 5] = params_.porosity0;
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
    
    // f0 = phi * ct_eff * dP/dt
    // Fields: 0=P, 1=Sw, 2=ux, 3=uy, 4=uz, 5=porosity
    
    const PetscScalar Sw    = u[uOff[1]];
    const PetscScalar phi   = u[uOff[5]];
    const PetscScalar P_t   = u_t[uOff[0]];
    const PetscScalar S_t   = u_t[uOff[1]];
    
    // Constants: [0]=porosity, [8]=cw, [9]=cg, [6]=rho_w, [7]=rho_g
    const PetscScalar phi_c = (numConstants > 0) ? constants[0] : 0.15;
    const PetscScalar cw    = (numConstants > 8) ? constants[8] : 4e-10;
    const PetscScalar cg    = (numConstants > 9) ? constants[9] : 1e-8;
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar rho_g = (numConstants > 7) ? constants[7] : 700.0;
    
    // Use field porosity if available, else constant
    PetscScalar phi_eff = (PetscRealPart(phi) > 0.01) ? phi : phi_c;
    
    // Effective compressibility: ct = Sw*cw + Sg*cg
    PetscScalar Sg = 1.0 - Sw;
    PetscScalar ct_eff = Sw * cw + Sg * cg;
    
    // Accumulation with density-difference coupling
    f0[0] = phi_eff * ct_eff * P_t + phi_eff * (rho_w - rho_g) * S_t * 1e-9;
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
    
    // f1 = total Darcy flux = Σ_α k * kr_α/μ_α * (∇P - ρ_α*g*ê_z)
    // For 3D: u_x[uOff_x[0] + d] = dP/dx_d for d = 0,1,2
    
    // Constants layout (updated):
    // [1] = permeability_x, [2] = permeability_y, [3] = permeability_z
    // [4] = water_viscosity, [5] = co2_viscosity
    // [6] = water_density, [7] = co2_density
    // [10-15] = rel perm params, [16] = gravity
    const PetscScalar kx   = (numConstants > 1) ? constants[1]  : 200e-15;
    const PetscScalar ky   = (numConstants > 2) ? constants[2]  : 200e-15;
    const PetscScalar kz   = (numConstants > 3) ? constants[3]  : 100e-15;
    const PetscScalar mu_w = (numConstants > 4) ? constants[4]  : 8e-4;
    const PetscScalar mu_g = (numConstants > 5) ? constants[5]  : 5e-5;
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar rho_g = (numConstants > 7) ? constants[7] : 700.0;
    const PetscScalar grav  = (numConstants > 16) ? constants[16] : 9.81;
    
    const PetscScalar Sw = u[uOff[1]];
    
    // Corey relative permeability
    PetscScalar Swr = (numConstants > 10) ? constants[10] : 0.2;
    PetscScalar Sgr = (numConstants > 11) ? constants[11] : 0.2;
    PetscScalar nw  = (numConstants > 12) ? constants[12] : 4.0;
    PetscScalar ng  = (numConstants > 13) ? constants[13] : 2.0;
    PetscScalar krw_max = (numConstants > 14) ? constants[14] : 1.0;
    PetscScalar krg_max = (numConstants > 15) ? constants[15] : 1.0;
    
    PetscScalar denom = 1.0 - Swr - Sgr;
    PetscScalar Se_w = std::max(0.0, std::min(1.0, PetscRealPart((Sw - Swr) / denom)));
    PetscScalar Se_g = std::max(0.0, std::min(1.0, PetscRealPart((1.0 - Sw - Sgr) / denom)));
    PetscScalar krw = krw_max * std::pow(Se_w, PetscRealPart(nw));
    PetscScalar krg = krg_max * std::pow(Se_g, PetscRealPart(ng));
    
    PetscScalar k_vals[3] = {kx, ky, kz};
    
    for (PetscInt d = 0; d < dim; ++d) {
        PetscScalar dP_dd = u_x[uOff_x[0] + d];
        PetscScalar k_d = (d < 3) ? k_vals[d] : kx;
        
        PetscScalar mob_w = k_d * krw / mu_w;
        PetscScalar mob_g = k_d * krg / mu_g;
        
        // Gravity in z-direction (last dimension in 3D)
        PetscScalar grav_w = (d == dim - 1) ? rho_w * grav : 0.0;
        PetscScalar grav_g = (d == dim - 1) ? rho_g * grav : 0.0;
        
        f1[d] = mob_w * (dP_dd - grav_w) + mob_g * (dP_dd - grav_g);
    }
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
    (void)u_t;
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
    
    // Extract solution fields (6-field layout: P, Sw, ux, uy, uz, phi)
    const PetscScalar P = u[uOff[0]];
    const PetscScalar S = u[uOff[1]];
    const PetscScalar phi = u[uOff[5]];
    
    // Material parameters from constants (18-element layout)
    // [0]=porosity, [1-3]=permeability, [4-5]=viscosity, [6-7]=density
    // [8-9]=compressibility, [10-11]=residual sat, [12-13]=corey exp
    // [14-15]=kr_max, [16]=gravity, [17]=biot
    PetscScalar alpha = (numConstants > 17) ? constants[17] : 1.0;    // Biot coefficient
    
    // Saturation-dependent compressibility
    PetscScalar cw = (numConstants > 8) ? constants[8] : 4.5e-10;  // Water compressibility
    PetscScalar co = (numConstants > 9) ? constants[9] : 1.5e-9;   // CO2/oil compressibility
    PetscScalar ct_eff = S * cw + (1.0 - S) * co;
    
    // ∂ct/∂S for saturation Jacobian contribution
    PetscScalar dct_dS = cw - co;
    
    // Initialize g0 to zero (Nf x Nf matrix, but we only fill pressure row)
    // For full coupled system, g0 would be 5x5
    
    // ∂f0/∂P: storage term Jacobian
    g0[0] = phi * ct_eff * u_tShift;
    
    // ∂f0/∂S: saturation coupling (through compressibility and relative permeability)
    if (Nf > 1) {
        // Get relative permeability parameters from constants (18-element layout)
        PetscScalar Swr = (numConstants > 10) ? constants[10] : 0.2;   // Water residual saturation
        PetscScalar Sor = (numConstants > 11) ? constants[11] : 0.2;   // Oil/CO2 residual saturation
        PetscScalar n_corey = (numConstants > 12) ? constants[12] : 4.0;  // Corey exponent
        PetscScalar kx = (numConstants > 1) ? constants[1] : 200e-15;   // Permeability
        PetscScalar mu = (numConstants > 4) ? constants[4] : 8e-4;      // Viscosity
        
        // Compressibility change with saturation
        g0[1] = phi * dct_dS * P * u_tShift;
        
        // Relative permeability derivative effect on flux
        // d(kr)/dS affects the transmissibility
        PetscScalar S_denom = 1.0 - Swr - Sor;
        PetscScalar Sw_norm = std::max(0.0, std::min(1.0, (S - Swr) / S_denom));
        PetscScalar dkr_dS = (Sw_norm > 0.0 && Sw_norm < 1.0) ? 
                             n_corey * std::pow(Sw_norm, n_corey - 1.0) / S_denom : 0.0;
        
        // Flux contribution (pressure gradient magnitude, 3D)
        PetscScalar grad_P_mag_sq = 0.0;
        for (PetscInt dd = 0; dd < dim; ++dd) {
            PetscScalar dP_dd = u_x[uOff_x[0] + dd];
            grad_P_mag_sq += dP_dd * dP_dd;
        }
        PetscScalar grad_P_mag = std::sqrt(PetscRealPart(grad_P_mag_sq));
        
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
}

void PoroelasticSolver::g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]) {
    (void)uOff_x;
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
    
    // Constants layout (18-element):
    // [1]=perm_x, [2]=perm_y, [3]=perm_z, [4]=mu_w, [5]=mu_g
    // [10]=Swr, [11]=Sgr, [12]=nw, [13]=ng, [14]=krw_max, [15]=krg_max
    PetscScalar kx = (numConstants > 1) ? constants[1] : 200e-15;  // Permeability x [m²]
    PetscScalar kz = (numConstants > 3) ? constants[3] : 100e-15;  // Permeability z [m²]
    PetscScalar mu = (numConstants > 4) ? constants[4] : 8e-4;     // Water viscosity [Pa·s]
    PetscScalar mu_g = (numConstants > 5) ? constants[5] : 5e-5;   // CO2 viscosity [Pa·s]
    
    // Relative permeability (Corey model)
    PetscScalar Swr = (numConstants > 10) ? constants[10] : 0.2;
    PetscScalar Sor = (numConstants > 11) ? constants[11] : 0.2;
    PetscScalar n_corey = (numConstants > 12) ? constants[12] : 4.0;
    PetscScalar krmax = (numConstants > 14) ? constants[14] : 1.0;
    
    PetscScalar S_norm = std::max(0.0, std::min(1.0, (S - Swr) / (1.0 - Swr - Sor)));
    PetscScalar kr = krmax * std::pow(S_norm, n_corey);
    
    // Mobilities
    PetscScalar mob_x = kx * kr / mu;
    PetscScalar mob_z = kz * kr / mu;
    
    // Fill g3 tensor (dim × dim matrix)
    // For 2D: g3[0] = ∂q_x/∂(∂P/∂x), g3[1] = ∂q_x/∂(∂P/∂z)
    //         g3[2] = ∂q_z/∂(∂P/∂x), g3[3] = ∂q_z/∂(∂P/∂z)
    // For 3D: 3×3 matrix
    
    // CO2/gas relative permeability
    PetscScalar Sg = 1.0 - S;
    PetscScalar Se_g = std::max(0.0, std::min(1.0, PetscRealPart((Sg - Sor) / (1.0 - Swr - Sor))));
    PetscScalar n_corey_g = (numConstants > 13) ? constants[13] : 2.0;
    PetscScalar krg_max = (numConstants > 15) ? constants[15] : 1.0;
    PetscScalar krg = krg_max * std::pow(Se_g, PetscRealPart(n_corey_g));
    
    if (dim == 2) {
        PetscScalar mob_total_x = kx * kr / mu + kx * krg / mu_g;
        PetscScalar mob_total_z = kz * kr / mu + kz * krg / mu_g;
        g3[0] = mob_total_x;
        g3[1] = 0.0;
        g3[2] = 0.0;
        g3[3] = mob_total_z;
    } else if (dim == 3) {
        PetscScalar ky = (numConstants > 2) ? constants[2] : kx;
        
        // Total mobility = water + CO2 for each direction
        PetscScalar mob_total_x = kx * kr / mu + kx * krg / mu_g;
        PetscScalar mob_total_y = ky * kr / mu + ky * krg / mu_g;
        PetscScalar mob_total_z = kz * kr / mu + kz * krg / mu_g;
        
        g3[0] = mob_total_x;  g3[1] = 0.0;          g3[2] = 0.0;           // Row 0
        g3[3] = 0.0;          g3[4] = mob_total_y;  g3[5] = 0.0;           // Row 1
        g3[6] = 0.0;          g3[7] = 0.0;          g3[8] = mob_total_z;   // Row 2
    }
    
    // For coupled poroelasticity, would also include:
    // ∂q/∂(∇u) terms if permeability depends on strain (stress-dependent permeability)
    // k = k0 * exp(α_k * (σ_eff - σ_ref)) where σ_eff = σ_total - α*p
    // This would add off-diagonal blocks to the full coupled Jacobian
}

// ============================================================================
// Pointwise residual functions for saturation equation (field 1)
// ============================================================================

void PoroelasticSolver::f0_saturation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                      const PetscInt uOff[], const PetscInt uOff_x[],
                                      const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                      const PetscInt aOff[], const PetscInt aOff_x[],
                                      const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                      PetscReal t, const PetscReal x[], PetscInt numConstants,
                                      const PetscScalar constants[], PetscScalar f0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    const PetscScalar S_t = u_t[uOff[1]];
    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;
    
    // f0_s = φ * dSw/dt
    f0[0] = phi * S_t;
}

void PoroelasticSolver::f1_saturation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                      const PetscInt uOff[], const PetscInt uOff_x[],
                                      const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                                      const PetscInt aOff[], const PetscInt aOff_x[],
                                      const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                                      PetscReal t, const PetscReal x[], PetscInt numConstants,
                                      const PetscScalar constants[], PetscScalar f1[]) {
    (void)Nf; (void)NfAux;
    (void)u; (void)u_t;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    const PetscScalar Sw = u[uOff[1]];
    
    const PetscScalar kx   = (numConstants > 1) ? constants[1]  : 200e-15;
    const PetscScalar ky   = (numConstants > 2) ? constants[2]  : 200e-15;
    const PetscScalar kz   = (numConstants > 3) ? constants[3]  : 100e-15;
    const PetscScalar mu_w = (numConstants > 4) ? constants[4]  : 8e-4;
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar grav  = (numConstants > 16) ? constants[16] : 9.81;
    
    // Corey water relative permeability
    PetscScalar Swr = (numConstants > 10) ? constants[10] : 0.2;
    PetscScalar Sgr = (numConstants > 11) ? constants[11] : 0.2;
    PetscScalar nw  = (numConstants > 12) ? constants[12] : 4.0;
    PetscScalar krw_max = (numConstants > 14) ? constants[14] : 1.0;
    
    PetscScalar Se = std::max(0.0, std::min(1.0,
        PetscRealPart((Sw - Swr) / (1.0 - Swr - Sgr))));
    PetscScalar krw = krw_max * std::pow(Se, PetscRealPart(nw));
    
    PetscScalar k_vals[3] = {kx, ky, kz};
    
    for (PetscInt d = 0; d < dim; ++d) {
        PetscScalar dP_dd = u_x[uOff_x[0] + d];
        PetscScalar k_d = (d < 3) ? k_vals[d] : kx;
        PetscScalar mob_w = k_d * krw / mu_w;
        PetscScalar grav_w = (d == dim - 1) ? rho_w * grav : 0.0;
        
        f1[d] = mob_w * (dP_dd - grav_w);
    }
}

// ============================================================================
// Cross-coupling Jacobian callbacks
// ============================================================================

void PoroelasticSolver::g0_ps(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    const PetscScalar P_t  = u_t[uOff[0]];
    const PetscScalar phi  = (numConstants > 0) ? constants[0] : 0.15;
    const PetscScalar cw   = (numConstants > 8) ? constants[8] : 4e-10;
    const PetscScalar cg   = (numConstants > 9) ? constants[9] : 1e-8;
    const PetscScalar rho_w = (numConstants > 6) ? constants[6] : 1100.0;
    const PetscScalar rho_g = (numConstants > 7) ? constants[7] : 700.0;
    
    // ∂f0_p/∂Sw = φ*(cw - cg)*dP/dt + φ*(ρw - ρg)*1e-9*shift
    g0[0] = phi * (cw - cg) * P_t + phi * (rho_w - rho_g) * 1e-9 * u_tShift;
}

void PoroelasticSolver::g0_ss(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x; (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)x;
    
    const PetscScalar phi = (numConstants > 0) ? constants[0] : 0.15;
    
    // ∂f0_s/∂Sw = φ * shift
    g0[0] = phi * u_tShift;
}

void PoroelasticSolver::g3_sp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                              PetscReal t, PetscReal u_tShift, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[]) {
    (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
    (void)t; (void)u_tShift; (void)x;
    
    const PetscScalar Sw = u[uOff[1]];
    
    const PetscScalar kx   = (numConstants > 1) ? constants[1]  : 200e-15;
    const PetscScalar ky   = (numConstants > 2) ? constants[2]  : 200e-15;
    const PetscScalar kz   = (numConstants > 3) ? constants[3]  : 100e-15;
    const PetscScalar mu_w = (numConstants > 4) ? constants[4]  : 8e-4;
    
    PetscScalar Swr = (numConstants > 10) ? constants[10] : 0.2;
    PetscScalar Sgr = (numConstants > 11) ? constants[11] : 0.2;
    PetscScalar nw  = (numConstants > 12) ? constants[12] : 4.0;
    PetscScalar krw_max = (numConstants > 14) ? constants[14] : 1.0;
    
    PetscScalar Se = std::max(0.0, std::min(1.0,
        PetscRealPart((Sw - Swr) / (1.0 - Swr - Sgr))));
    PetscScalar krw = krw_max * std::pow(Se, PetscRealPart(nw));
    
    PetscScalar k_vals[3] = {kx, ky, kz};
    
    // ∂f1_s/∂(∇P) = water mobility tensor (diagonal)
    for (PetscInt i = 0; i < dim; ++i) {
        for (PetscInt j = 0; j < dim; ++j) {
            if (i == j) {
                PetscScalar k_d = (i < 3) ? k_vals[i] : kx;
                g3[i * dim + j] = k_d * krw / mu_w;
            } else {
                g3[i * dim + j] = 0.0;
            }
        }
    }
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
}

void PoroelasticSolver::getDisplacement3D(std::vector<std::vector<double>>& ux,
                                           std::vector<std::vector<double>>& uy,
                                           std::vector<std::vector<double>>& uz) const {
    ux.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
    uy.resize(grid_cells_.z, std::vector<double>(grid_cells_.x, 0.0));
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
        
        if (dof >= 6) {
            PetscInt cell_idx = c - cStart;
            PetscInt i = cell_idx % grid_cells_.x;
            PetscInt k = cell_idx / grid_cells_.x;
            
            if (k < grid_cells_.z && i < grid_cells_.x) {
                ux[k][i] = array[offset + 2];
                uy[k][i] = array[offset + 3];
                uz[k][i] = array[offset + 4];
            }
        }
    }
    
    VecRestoreArray(local, &array);
    DMRestoreLocalVector(dm_, &local);
}

void PoroelasticSolver::getPorosity(std::vector<std::vector<double>>& phi) const {
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
        
        if (dof >= 6) {
            PetscInt cell_idx = c - cStart;
            PetscInt i = cell_idx % grid_cells_.x;
            PetscInt k = cell_idx / grid_cells_.x;
            
            if (k < grid_cells_.z && i < grid_cells_.x) {
                phi[k][i] = array[offset + 5];
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
    (void)time;
    (void)dt;
    // TODO: Implement well source terms
}

PetscErrorCode PoroelasticSolver::FormIFunction(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void* ctx) {
    (void)ts;
    (void)t;
    (void)X;
    (void)Xdot;
    (void)F;
    (void)ctx;
    // DMPlex uses pointwise functions, this is just a wrapper
    return 0;
}

PetscErrorCode PoroelasticSolver::FormIJacobian(TS ts, PetscReal t, Vec X, Vec Xdot,
                                                PetscReal shift, Mat J, Mat P, void* ctx) {
    (void)ts;
    (void)t;
    (void)X;
    (void)Xdot;
    (void)shift;
    (void)J;
    (void)P;
    (void)ctx;
    // DMPlex uses pointwise functions, this is just a wrapper
    return 0;
}

} // namespace FSRM
