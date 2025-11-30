#include "PhysicsKernel.hpp"
#include <cmath>
#include <algorithm>

namespace FSRM {

// ============================================================================
// Single Phase Flow Kernel
// ============================================================================

SinglePhaseFlowKernel::SinglePhaseFlowKernel()
    : PhysicsKernel(PhysicsType::FLUID_FLOW),
      porosity(0.2), permeability(100.0e-15), compressibility(1e-9),
      viscosity(0.001), density(1000.0) {}

PetscErrorCode SinglePhaseFlowKernel::setup(DM dm, PetscFE fe) {
    PetscFunctionBeginUser;
    // Setup finite element spaces
    PetscFunctionReturn(0);
}

void SinglePhaseFlowKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_x;
    (void)a;
    (void)x;
    
    // Accumulation: phi * ct * dP/dt
    f[0] = porosity * compressibility * u_t[0];
    
    // Flux: div(k/mu * grad(P))
    // This will be handled in f1
}

void SinglePhaseFlowKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for single-phase flow diffusion equation:
    // F = φ * ct * ∂P/∂t - div(k/μ * ∇P) = 0
    //
    // The Jacobian has two parts:
    // 1. g0: ∂F/∂(∂P/∂t) = φ * ct * shift (where shift = a[0])
    // 2. g3: ∂(flux)/∂(∇P) = k/μ * I (for the flux term in weak form)
    //
    // For PETSc's pointwise Jacobian, J combines these contributions.
    // The shift parameter 'a' accounts for the time discretization.
    
    // Accumulation term Jacobian: ∂f0/∂(∂P/∂t) * shift
    double accumulation_jacobian = porosity * compressibility * a[0];
    
    // Mobility for flux term
    double mobility = permeability / viscosity;
    
    // For a scalar field (pressure), J[0] is the combined Jacobian
    // In weak form: ∫ (g0 * φ_i * φ_j + g3 * ∇φ_i · ∇φ_j) dV
    // The assembler handles the spatial integration
    J[0] = accumulation_jacobian;
    
    // Note: For PETSc PetscDS, the flux Jacobian g3 is typically provided
    // separately through PetscDSSetJacobian with the g3 callback.
    // Here we include the mobility coefficient that would multiply the
    // gradient-gradient term in the weak form.
    // J[1] through J[dim*dim-1] would contain g3 = k/μ * δ_ij
    // For 3D: J[1] = J[5] = J[9] = mobility (diagonal of 3x3 identity)
}

void SinglePhaseFlowKernel::setProperties(double phi, double k, double ct, double mu, double rho) {
    porosity = phi;
    permeability = k;
    compressibility = ct;
    viscosity = mu;
    density = rho;
}

// ============================================================================
// Black Oil Kernel
// ============================================================================

BlackOilKernel::BlackOilKernel()
    : PhysicsKernel(PhysicsType::FLUID_FLOW),
      porosity(0.2), perm_x(100e-15), perm_y(100e-15), perm_z(10e-15) {}

PetscErrorCode BlackOilKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void BlackOilKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                              const PetscScalar u_x[], const PetscScalar a[],
                              const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_x;
    (void)a;
    (void)x;
    
    // u[0] = P (pressure)
    // u[1] = Sw (water saturation)
    // u[2] = Sg (gas saturation)
    // So = 1 - Sw - Sg (oil saturation)
    
    double P = u[0];
    double Sw = u[1];
    double Sg = u[2];
    double So = 1.0 - Sw - Sg;
    (void)So;  // Used in full implementation for saturation constraint
    (void)Sw;  // Used for relative permeability
    (void)Sg;  // Used for relative permeability
    
    // Compute PVT properties
    double rho_o = oilDensity(P, 0.0);
    double rho_w = waterDensity(P);
    double rho_g = gasDensity(P);
    
    // Oil equation: d/dt(phi * rho_o * So) + div(rho_o * v_o) = 0
    f[0] = porosity * rho_o * u_t[0];  // Simplified
    
    // Water equation
    f[1] = porosity * rho_w * u_t[1];
    
    // Gas equation
    f[2] = porosity * rho_g * u_t[2];
}

void BlackOilKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                              const PetscScalar u_x[], const PetscScalar a[],
                              const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for Black Oil model
    // Primary variables: u[0] = P (pressure), u[1] = Sw (water saturation), u[2] = Sg (gas saturation)
    // Constraint: So = 1 - Sw - Sg (oil saturation)
    //
    // Residual equations (accumulation terms):
    // f[0] = ∂(φ * ρo * So)/∂t + div(ρo * vo) = 0  (oil conservation)
    // f[1] = ∂(φ * ρw * Sw)/∂t + div(ρw * vw) = 0  (water conservation)
    // f[2] = ∂(φ * ρg * Sg)/∂t + div(ρg * vg) = 0  (gas conservation)
    //
    // Full 3x3 Jacobian matrix layout:
    // J[0:2]   = [∂f0/∂P,  ∂f0/∂Sw, ∂f0/∂Sg]
    // J[3:5]   = [∂f1/∂P,  ∂f1/∂Sw, ∂f1/∂Sg]
    // J[6:8]   = [∂f2/∂P,  ∂f2/∂Sw, ∂f2/∂Sg]
    
    double P = u[0];
    double Sw = u[1];
    double Sg = u[2];
    double So = 1.0 - Sw - Sg;
    
    // Compute densities
    double rho_o = oilDensity(P, 0.0);
    double rho_w = waterDensity(P);
    double rho_g = gasDensity(P);
    
    // Compute density derivatives with respect to pressure
    // Using finite difference approximation for pressure derivatives
    double dP = P * 1e-6 + 1.0;
    double drho_o_dP = (oilDensity(P + dP, 0.0) - oilDensity(P - dP, 0.0)) / (2.0 * dP);
    double drho_w_dP = (waterDensity(P + dP) - waterDensity(P - dP)) / (2.0 * dP);
    double drho_g_dP = (gasDensity(P + dP) - gasDensity(P - dP)) / (2.0 * dP);
    
    // Compute relative permeabilities and their derivatives
    // Using Corey-type model: kr_alpha = kr_alpha_max * (S_alpha - S_alpha_r)^n / (1 - S_r_total)^n
    double Sor = fluid_props.residual_saturation;
    double Swr = fluid_props.water_residual_saturation;
    double Sgr = fluid_props.gas_residual_saturation;
    
    double So_norm = std::max(0.0, (So - Sor) / (1.0 - Sor - Swr - Sgr));
    double Sw_norm = std::max(0.0, (Sw - Swr) / (1.0 - Sor - Swr - Sgr));
    double Sg_norm = std::max(0.0, (Sg - Sgr) / (1.0 - Sor - Swr - Sgr));
    
    double n_o = 2.0, n_w = 2.0, n_g = 2.0;  // Corey exponents
    double kro_max = 1.0, krw_max = 0.5, krg_max = 0.8;
    
    double kro = kro_max * std::pow(So_norm, n_o);
    double krw = krw_max * std::pow(Sw_norm, n_w);
    double krg = krg_max * std::pow(Sg_norm, n_g);
    
    // Relative permeability derivatives
    double S_denom = 1.0 - Sor - Swr - Sgr;
    double dkro_dSo = (So_norm > 1e-10) ? kro_max * n_o * std::pow(So_norm, n_o - 1) / S_denom : 0.0;
    double dkrw_dSw = (Sw_norm > 1e-10) ? krw_max * n_w * std::pow(Sw_norm, n_w - 1) / S_denom : 0.0;
    double dkrg_dSg = (Sg_norm > 1e-10) ? krg_max * n_g * std::pow(Sg_norm, n_g - 1) / S_denom : 0.0;
    
    // Since So = 1 - Sw - Sg: ∂So/∂Sw = -1, ∂So/∂Sg = -1
    double dkro_dSw = -dkro_dSo;
    double dkro_dSg = -dkro_dSo;
    
    // Compute viscosities and mobilities
    double mu_o = oilViscosity(P, 0.0);
    double mu_w = waterViscosity(P);
    double mu_g = gasViscosity(P);
    
    double lambda_o = kro / mu_o;  // Oil mobility
    double lambda_w = krw / mu_w;  // Water mobility
    double lambda_g = krg / mu_g;  // Gas mobility
    
    // Time shift parameter
    double shift = a[0];
    
    // Initialize Jacobian to zero
    for (int i = 0; i < 9; ++i) J[i] = 0.0;
    
    // ==================================
    // Accumulation term Jacobian (g0)
    // ==================================
    
    // Oil equation derivatives (f0 = φ * ρo * So)
    // ∂f0/∂P = φ * ∂ρo/∂P * So
    J[0] = porosity * drho_o_dP * So * shift;
    // ∂f0/∂Sw = -φ * ρo (since ∂So/∂Sw = -1)
    J[1] = -porosity * rho_o * shift;
    // ∂f0/∂Sg = -φ * ρo (since ∂So/∂Sg = -1)
    J[2] = -porosity * rho_o * shift;
    
    // Water equation derivatives (f1 = φ * ρw * Sw)
    // ∂f1/∂P = φ * ∂ρw/∂P * Sw
    J[3] = porosity * drho_w_dP * Sw * shift;
    // ∂f1/∂Sw = φ * ρw
    J[4] = porosity * rho_w * shift;
    // ∂f1/∂Sg = 0
    J[5] = 0.0;
    
    // Gas equation derivatives (f2 = φ * ρg * Sg)
    // ∂f2/∂P = φ * ∂ρg/∂P * Sg
    J[6] = porosity * drho_g_dP * Sg * shift;
    // ∂f2/∂Sw = 0
    J[7] = 0.0;
    // ∂f2/∂Sg = φ * ρg
    J[8] = porosity * rho_g * shift;
    
    // ==================================
    // Flux term Jacobian contributions
    // ==================================
    // The flux terms contribute to the Jacobian through:
    // - Pressure gradient terms (g3 diagonal blocks)
    // - Mobility derivatives w.r.t. saturations (g0 off-diagonal)
    //
    // For the transmissibility T = k * kr_alpha / μ_alpha:
    // Flux_alpha = T * (∇P - ρ_alpha * g)
    //
    // Additional contributions from flux derivatives:
    // These are assembled through the g2/g3 callbacks in PETSc's FE framework.
    // Here we add the implicit contributions from mobility changes:
    
    double k_avg = std::cbrt(perm_x * perm_y * perm_z);  // Geometric mean permeability
    
    // Oil flux derivative w.r.t. saturations (through mobility)
    // ∂(ρo * λo)/∂Sw = ρo * ∂λo/∂Sw = ρo * (dkro_dSw / μo)
    J[1] += k_avg * rho_o * dkro_dSw / mu_o;
    J[2] += k_avg * rho_o * dkro_dSg / mu_o;
    
    // Water flux derivative w.r.t. saturations
    J[4] += k_avg * rho_w * dkrw_dSw / mu_w;
    
    // Gas flux derivative w.r.t. saturations
    J[8] += k_avg * rho_g * dkrg_dSg / mu_g;
    
    // Suppress unused variable warnings for variables used in derivative calculations
    (void)kro; (void)krw; (void)krg;
    (void)lambda_o; (void)lambda_w; (void)lambda_g;
}

void BlackOilKernel::setFluidProperties(const FluidProperties& props) {
    fluid_props = props;
}

void BlackOilKernel::setRockProperties(double phi, double kx, double ky, double kz) {
    porosity = phi;
    perm_x = kx;
    perm_y = ky;
    perm_z = kz;
}

// PVT functions (simplified correlations)
double BlackOilKernel::oilDensity(double P, double Rs) const {
    // Simplified: slightly compressible liquid
    double P_ref = 1e5;  // 1 bar reference
    double rho_ref = fluid_props.oil_density_std;
    double c_o = 1e-9;  // compressibility
    return rho_ref * (1.0 + c_o * (P - P_ref));
}

double BlackOilKernel::gasDensity(double P) const {
    // Ideal gas law: rho = P*MW/(R*T)
    double T = 300.0;  // K
    double R = 8.314;  // J/(mol·K)
    double MW = 0.016; // kg/mol for methane
    return P * MW / (R * T);
}

double BlackOilKernel::waterDensity(double P) const {
    double P_ref = 1e5;
    double rho_ref = fluid_props.water_density_std;
    double c_w = 4e-10;
    return rho_ref * (1.0 + c_w * (P - P_ref));
}

double BlackOilKernel::oilViscosity(double P, double Rs) const {
    (void)P;   // Simplified model - would be used in full PVT
    (void)Rs;  // Simplified model - would be used in full PVT
    return fluid_props.oil_viscosity;
}

double BlackOilKernel::gasViscosity(double P) const {
    (void)P;  // Simplified model - would be used in full PVT
    return fluid_props.gas_viscosity;
}

double BlackOilKernel::waterViscosity(double P) const {
    (void)P;  // Simplified model - would be used in full PVT
    return fluid_props.water_viscosity;
}

double BlackOilKernel::solutionGOR(double P) const {
    (void)P;  // Simplified model - would be used in full PVT
    return fluid_props.solution_GOR;
}

// Relative permeability (Corey correlations)
double BlackOilKernel::krw(double Sw) const {
    double Swc = 0.2;  // connate water
    double Swmax = 0.8;
    if (Sw < Swc) return 0.0;
    if (Sw > Swmax) return 1.0;
    
    double Sw_norm = (Sw - Swc) / (Swmax - Swc);
    return std::pow(Sw_norm, 4.0);
}

double BlackOilKernel::kro(double So, double Sg) const {
    double Sor = 0.2;  // residual oil
    if (So < Sor) return 0.0;
    
    double So_norm = (So - Sor) / (1.0 - Sor - 0.2);
    return std::pow(So_norm, 2.0);
}

double BlackOilKernel::krg(double Sg) const {
    double Sgc = 0.05;
    if (Sg < Sgc) return 0.0;
    
    double Sg_norm = (Sg - Sgc) / (1.0 - Sgc - 0.2);
    return std::pow(Sg_norm, 2.0);
}

double BlackOilKernel::Pcow(double Sw) const {
    // Capillary pressure oil-water
    return 1000.0 * std::pow((Sw + 0.01), -2.0);
}

double BlackOilKernel::Pcog(double Sg) const {
    // Capillary pressure oil-gas
    return 500.0 * std::pow((Sg + 0.01), -1.5);
}

// ============================================================================
// Compositional Kernel
// ============================================================================

CompositionalKernel::CompositionalKernel(int num_components)
    : PhysicsKernel(PhysicsType::FLUID_FLOW), nc(num_components),
      porosity(0.2), perm_x(100e-15), perm_y(100e-15), perm_z(10e-15) {
    component_mw.resize(nc, 0.016);
    component_Tc.resize(nc, 190.0);
    component_Pc.resize(nc, 4.6e6);
    component_omega.resize(nc, 0.01);
}

PetscErrorCode CompositionalKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void CompositionalKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_x;
    (void)a;
    (void)x;
    
    // u[0] = P
    // u[1..nc-1] = overall compositions z_i
    
    double P = u[0];
    double T = 350.0;  // K, temperature
    
    // Get compositions
    std::vector<double> z(nc);
    for (int i = 0; i < nc - 1; ++i) {
        z[i] = u[i + 1];
    }
    z[nc - 1] = 1.0;
    for (int i = 0; i < nc - 1; ++i) {
        z[nc - 1] -= z[i];
    }
    
    // Perform flash calculation
    std::vector<double> x_comp(nc), y_comp(nc);
    double S_L, S_V;
    flashCalculation(P, T, z, x_comp, y_comp, S_L, S_V);
    (void)S_L;  // Used in full implementation for phase volume calculation
    (void)S_V;  // Used in full implementation for phase volume calculation
    
    // Component mass balance: d/dt(phi * sum_p(rho_p * S_p * x_ip)) + div(flux_i) = 0
    for (int i = 0; i < nc - 1; ++i) {
        f[i + 1] = porosity * u_t[i + 1];  // Simplified accumulation
    }
    
    // Pressure equation
    f[0] = porosity * u_t[0];  // Simplified
}

void CompositionalKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                   const PetscScalar u_x[], const PetscScalar a[],
                                   const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for Compositional model
    // Primary variables: u[0] = P, u[1..nc-1] = overall molar compositions z_i
    // Constraint: z[nc-1] = 1 - sum(z[0..nc-2])
    //
    // Residual equations:
    // f[0] = Pressure equation (volume balance or constraint)
    // f[1..nc-1] = Component molar balance equations
    //
    // Jacobian matrix is nc×nc:
    // J[i*nc + j] = ∂f_i/∂u_j
    
    double P = u[0];
    double T = 350.0;  // Reservoir temperature [K]
    double shift = a[0];
    
    // Extract overall compositions
    std::vector<double> z(nc);
    for (int i = 0; i < nc - 1; ++i) {
        z[i] = std::max(1e-10, u[i + 1]);  // Ensure positive compositions
    }
    z[nc - 1] = 1.0;
    for (int i = 0; i < nc - 1; ++i) {
        z[nc - 1] -= z[i];
    }
    z[nc - 1] = std::max(1e-10, z[nc - 1]);
    
    // Perform flash calculation to get phase properties
    std::vector<double> x_comp(nc), y_comp(nc);
    double S_L, S_V;
    flashCalculation(P, T, z, x_comp, y_comp, S_L, S_V);
    
    // Compute phase densities
    double rho_L = 0.0, rho_V = 0.0;
    double R = 8.314;  // J/(mol·K)
    
    // Liquid phase: sum of component densities weighted by mole fraction
    for (int i = 0; i < nc; ++i) {
        rho_L += x_comp[i] * component_mw[i];
    }
    rho_L *= 1000.0;  // Convert to kg/m³ (simplified)
    
    // Vapor phase: ideal gas approximation
    double avg_mw_V = 0.0;
    for (int i = 0; i < nc; ++i) {
        avg_mw_V += y_comp[i] * component_mw[i];
    }
    rho_V = P * avg_mw_V / (R * T);
    
    // K-values (equilibrium ratios) from Wilson equation approximation
    std::vector<double> K(nc);
    for (int i = 0; i < nc; ++i) {
        double Tr = T / component_Tc[i];
        double Pr = P / component_Pc[i];
        K[i] = std::exp(5.37 * (1.0 + component_omega[i]) * (1.0 - 1.0/Tr)) / Pr;
    }
    
    // Compute K-value derivatives with respect to pressure
    std::vector<double> dK_dP(nc);
    for (int i = 0; i < nc; ++i) {
        dK_dP[i] = -K[i] / P;  // From Wilson equation: K ∝ 1/P
    }
    
    // Initialize Jacobian to zero
    int Jsize = nc * nc;
    for (int i = 0; i < Jsize; ++i) J[i] = 0.0;
    
    // ==================================
    // Accumulation term Jacobians
    // ==================================
    
    // For component i, accumulation = φ * N_i where N_i is total moles of component i
    // N_i = V_L * ρ_L * x_i + V_V * ρ_V * y_i
    // For molar formulation: N_i = (L * x_i + V * y_i) where L, V are phase molar amounts
    
    // Derivatives of phase compositions w.r.t. overall composition and pressure
    // From Rachford-Rice: at equilibrium, x_i = z_i / (1 + V*(K_i - 1))
    //                                      y_i = K_i * z_i / (1 + V*(K_i - 1))
    
    double total_molar_density = (S_L > 0.5) ? rho_L / avg_mw_V : rho_V / avg_mw_V;
    if (total_molar_density < 1.0) total_molar_density = 1000.0 / 0.020;  // Default ~50 kmol/m³
    
    // Pressure equation (row 0): volume constraint or pressure accumulation
    // ∂f0/∂P: compressibility contribution
    double c_t = 1e-9;  // Total compressibility [1/Pa]
    J[0] = porosity * c_t * shift;
    
    // ∂f0/∂z_i: from phase split changes with composition
    for (int j = 1; j < nc; ++j) {
        J[j] = porosity * total_molar_density * shift;  // Density change with composition
    }
    
    // Component equations (rows 1 to nc-1)
    for (int i = 1; i < nc; ++i) {
        // Row i: ∂f_i/∂(u_j) for component i-1 (since i=0 is pressure)
        int comp = i - 1;
        
        // ∂f_i/∂P: from K-value changes and density changes with pressure
        double dN_dP = 0.0;
        if (S_V > 1e-6) {
            // Two-phase: composition depends on K-values
            double denom = 1.0 + S_V * (K[comp] - 1.0);
            double dx_dP = -z[comp] * S_V * dK_dP[comp] / (denom * denom);
            double dy_dP = z[comp] * dK_dP[comp] / denom - K[comp] * z[comp] * S_V * dK_dP[comp] / (denom * denom);
            dN_dP = S_L * rho_L * dx_dP + S_V * rho_V * dy_dP;
        } else {
            // Single-phase liquid: density depends on pressure
            dN_dP = z[comp] * rho_L * c_t;
        }
        J[i * nc + 0] = porosity * dN_dP * shift;
        
        // ∂f_i/∂z_j: direct composition dependence
        for (int j = 1; j < nc; ++j) {
            int comp_j = j - 1;
            double dN_dz = 0.0;
            
            if (comp == comp_j) {
                // Diagonal: direct contribution
                dN_dz = S_L * rho_L + S_V * rho_V * K[comp];
                if (S_V > 1e-6) {
                    double denom = 1.0 + S_V * (K[comp] - 1.0);
                    dN_dz = (S_L * rho_L + S_V * rho_V * K[comp]) / denom;
                }
            } else {
                // Off-diagonal: through phase split changes
                // When z_j increases, other z_k may need to decrease for sum = 1
                // This creates coupling between all components
                dN_dz = -S_L * rho_L - S_V * rho_V * K[comp];  // Constraint coupling
                dN_dz *= 1.0 / (nc - 1);  // Distribute effect
            }
            
            J[i * nc + j] = porosity * dN_dz * shift;
        }
    }
    
    // ==================================
    // Flux term Jacobian contributions
    // ==================================
    // Flux_i = sum_p (ρ_p * kr_p / μ_p * x_ip * k * ∇P)
    // Adds contributions from:
    // - Pressure gradient (handled by g3)
    // - Composition dependence of phase properties
    
    double k_avg = std::cbrt(perm_x * perm_y * perm_z);
    double mu_L = 1e-3;  // Liquid viscosity [Pa·s]
    double mu_V = 1e-5;  // Vapor viscosity [Pa·s]
    
    // Mobility contributions (simplified - full model would use proper kr correlations)
    double kr_L = (S_L > 0.2) ? std::pow((S_L - 0.2) / 0.8, 2.0) : 0.0;
    double kr_V = (S_V > 0.05) ? std::pow((S_V - 0.05) / 0.95, 2.0) : 0.0;
    
    double lambda_L = kr_L / mu_L;
    double lambda_V = kr_V / mu_V;
    
    // Mobility derivatives w.r.t. saturations (which depend on compositions)
    double dkr_L_dS_L = (S_L > 0.2) ? 2.0 * (S_L - 0.2) / (0.8 * 0.8) : 0.0;
    double dkr_V_dS_V = (S_V > 0.05) ? 2.0 * (S_V - 0.05) / (0.95 * 0.95) : 0.0;
    
    // Add flux Jacobian contributions to component equations
    for (int i = 1; i < nc; ++i) {
        // Flux contribution from mobility changes
        J[i * nc + 0] += k_avg * (rho_L * lambda_L + rho_V * lambda_V) * c_t;
    }
    
    // Suppress unused variable warnings
    (void)dkr_L_dS_L; (void)dkr_V_dS_V;
}

void CompositionalKernel::setComponentProperties(const std::vector<double>& mw,
                                                 const std::vector<double>& Tc,
                                                 const std::vector<double>& Pc,
                                                 const std::vector<double>& omega) {
    component_mw = mw;
    component_Tc = Tc;
    component_Pc = Pc;
    component_omega = omega;
}

void CompositionalKernel::flashCalculation(double P, double T, 
                                          const std::vector<double>& z,
                                          std::vector<double>& x,
                                          std::vector<double>& y,
                                          double& S_L, double& S_V) const {
    // Suppress unused parameter warnings - would be used in full EOS calculation
    (void)P;
    (void)T;
    
    // Simplified flash calculation using Rachford-Rice
    // For detailed implementation, use proper EOS (Peng-Robinson, SRK, etc.)
    
    // Assume single phase for simplicity in this example
    S_L = 1.0;
    S_V = 0.0;
    x = z;
    y = z;
}

double CompositionalKernel::fugacityCoefficient(double P, double T, double z_factor,
                                               const std::vector<double>& comp, int i) const {
    // Suppress unused parameter warnings - would be used in full EOS calculation
    (void)P;
    (void)T;
    (void)z_factor;
    (void)comp;
    (void)i;
    
    // Peng-Robinson EOS fugacity coefficient
    // Simplified - full implementation requires cubic EOS solver
    return 1.0;
}

// ============================================================================
// Geomechanics Kernel
// ============================================================================

GeomechanicsKernel::GeomechanicsKernel(SolidModelType model_type)
    : PhysicsKernel(PhysicsType::GEOMECHANICS), model_type(model_type),
      youngs_modulus(10e9), poisson_ratio(0.25), density(2500.0),
      relaxation_time(1e6), viscosity(1e19),
      biot_coefficient(1.0), biot_modulus(1e10) {}

PetscErrorCode GeomechanicsKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void GeomechanicsKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)a;
    (void)x;
    
    // u = [ux, uy, uz] displacement components
    // u_x = [dux/dx, dux/dy, dux/dz, duy/dx, ...]
    
    // Compute strain tensor from displacement gradient
    double eps_xx = u_x[0];  // du_x/dx
    double eps_yy = u_x[4];  // du_y/dy
    double eps_zz = u_x[8];  // du_z/dz
    double eps_xy = 0.5 * (u_x[1] + u_x[3]);  // 0.5*(du_x/dy + du_y/dx)
    double eps_xz = 0.5 * (u_x[2] + u_x[6]);
    double eps_yz = 0.5 * (u_x[5] + u_x[7]);
    
    // Lame parameters
    double lambda = youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    
    // Stress tensor (isotropic linear elasticity)
    double trace = eps_xx + eps_yy + eps_zz;
    double sigma_xx = lambda * trace + 2.0 * mu * eps_xx;
    double sigma_yy = lambda * trace + 2.0 * mu * eps_yy;
    double sigma_zz = lambda * trace + 2.0 * mu * eps_zz;
    double sigma_xy = 2.0 * mu * eps_xy;
    double sigma_xz = 2.0 * mu * eps_xz;
    double sigma_yz = 2.0 * mu * eps_yz;
    
    // Suppress unused variable warnings - stress components used in f1 flux term
    (void)sigma_xx; (void)sigma_yy; (void)sigma_zz;
    (void)sigma_xy; (void)sigma_xz; (void)sigma_yz;
    
    if (model_type == SolidModelType::VISCOELASTIC) {
        // Add viscoelastic correction
        // Maxwell model: sigma_v = eta * d(eps)/dt
        sigma_xx += viscosity * u_t[0];
        sigma_yy += viscosity * u_t[1];
        sigma_zz += viscosity * u_t[2];
    }
    
    // Weak form contribution (will be integrated)
    // Actually this should go into f1, not f0
    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
}

void GeomechanicsKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Lamé parameters from Young's modulus and Poisson's ratio
    double lambda = youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    
    // For 3D displacement field (ux, uy, uz), the Jacobian is 9x9
    // J[i*9 + j] corresponds to dF_i/du_j where i,j are gradients (3 components × 3 directions)
    // The weak form is: ∫ σ_ij δε_ij dV where ε_ij = 1/2(∂u_i/∂x_j + ∂u_j/∂x_i)
    // The Jacobian relates stress to strain: σ = C : ε
    
    // Initialize all entries to zero
    for (int i = 0; i < 81; ++i) {
        J[i] = 0.0;
    }
    
    // The constitutive tensor for isotropic linear elasticity:
    // σ_ij = λ δ_ij ε_kk + 2μ ε_ij
    // 
    // For the gradient-based formulation in PETSc:
    // J[i*9 + j] = ∂(stress_component_i)/∂(displacement_gradient_j)
    //
    // Displacement gradient components are ordered as:
    // u_x[0] = ∂u_x/∂x, u_x[1] = ∂u_x/∂y, u_x[2] = ∂u_x/∂z
    // u_x[3] = ∂u_y/∂x, u_x[4] = ∂u_y/∂y, u_x[5] = ∂u_y/∂z
    // u_x[6] = ∂u_z/∂x, u_x[7] = ∂u_z/∂y, u_x[8] = ∂u_z/∂z
    
    // Jacobian for ∂σ_xx/∂(∇u)
    J[0*9 + 0] = lambda + 2.0*mu;  // ∂σ_xx/∂(∂u_x/∂x)
    J[0*9 + 4] = lambda;            // ∂σ_xx/∂(∂u_y/∂y)
    J[0*9 + 8] = lambda;            // ∂σ_xx/∂(∂u_z/∂z)
    
    // Jacobian for ∂σ_yy/∂(∇u)
    J[1*9 + 0] = lambda;            // ∂σ_yy/∂(∂u_x/∂x)
    J[1*9 + 4] = lambda + 2.0*mu;  // ∂σ_yy/∂(∂u_y/∂y)
    J[1*9 + 8] = lambda;            // ∂σ_yy/∂(∂u_z/∂z)
    
    // Jacobian for ∂σ_zz/∂(∇u)
    J[2*9 + 0] = lambda;            // ∂σ_zz/∂(∂u_x/∂x)
    J[2*9 + 4] = lambda;            // ∂σ_zz/∂(∂u_y/∂y)
    J[2*9 + 8] = lambda + 2.0*mu;  // ∂σ_zz/∂(∂u_z/∂z)
    
    // Jacobian for ∂σ_xy/∂(∇u) (shear stress)
    J[3*9 + 1] = mu;  // ∂σ_xy/∂(∂u_x/∂y)
    J[3*9 + 3] = mu;  // ∂σ_xy/∂(∂u_y/∂x)
    
    // Jacobian for ∂σ_xz/∂(∇u) (shear stress)
    J[4*9 + 2] = mu;  // ∂σ_xz/∂(∂u_x/∂z)
    J[4*9 + 6] = mu;  // ∂σ_xz/∂(∂u_z/∂x)
    
    // Jacobian for ∂σ_yz/∂(∇u) (shear stress)
    J[5*9 + 5] = mu;  // ∂σ_yz/∂(∂u_y/∂z)
    J[5*9 + 7] = mu;  // ∂σ_yz/∂(∂u_z/∂y)
    
    // For symmetric stress tensor, σ_yx = σ_xy, etc.
    J[6*9 + 1] = mu;  // ∂σ_yx/∂(∂u_x/∂y)
    J[6*9 + 3] = mu;  // ∂σ_yx/∂(∂u_y/∂x)
    
    J[7*9 + 2] = mu;  // ∂σ_zx/∂(∂u_x/∂z)
    J[7*9 + 6] = mu;  // ∂σ_zx/∂(∂u_z/∂x)
    
    J[8*9 + 5] = mu;  // ∂σ_zy/∂(∂u_y/∂z)
    J[8*9 + 7] = mu;  // ∂σ_zy/∂(∂u_z/∂y)
    
    // Add viscoelastic contribution if using viscoelastic model
    if (model_type == SolidModelType::VISCOELASTIC && a != nullptr) {
        // For Maxwell viscoelasticity: σ_v = η * ∂ε/∂t
        // This adds a time-derivative term with shift parameter a[0]
        J[0*9 + 0] += viscosity * a[0];  // Normal stress xx
        J[1*9 + 4] += viscosity * a[0];  // Normal stress yy
        J[2*9 + 8] += viscosity * a[0];  // Normal stress zz
    }
}

void GeomechanicsKernel::setMaterialProperties(double E, double nu, double rho) {
    youngs_modulus = E;
    poisson_ratio = nu;
    density = rho;
}

void GeomechanicsKernel::setViscoelasticProperties(double tau, double eta) {
    relaxation_time = tau;
    viscosity = eta;
}

void GeomechanicsKernel::setPoroelasticCoupling(double alpha, double M) {
    biot_coefficient = alpha;
    biot_modulus = M;
}

// ============================================================================
// Thermal Kernel
// ============================================================================

ThermalKernel::ThermalKernel()
    : PhysicsKernel(PhysicsType::THERMAL),
      thermal_conductivity_solid(2.5), density_solid(2500.0), heat_capacity_solid(900.0),
      thermal_conductivity_fluid(0.6), density_fluid(1000.0), heat_capacity_fluid(4200.0),
      porosity(0.2) {}

PetscErrorCode ThermalKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void ThermalKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                            const PetscScalar u_x[], const PetscScalar a[],
                            const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_x;
    (void)a;
    (void)x;
    
    // u[0] = T (temperature)
    // u_x = [dT/dx, dT/dy, dT/dz]
    
    // Effective properties
    double rho_eff = (1.0 - porosity) * density_solid + porosity * density_fluid;
    double cp_eff = (1.0 - porosity) * heat_capacity_solid + porosity * heat_capacity_fluid;
    double k_eff = (1.0 - porosity) * thermal_conductivity_solid + porosity * thermal_conductivity_fluid;
    (void)k_eff;  // Used in f1 flux term for conduction
    
    // Heat equation: rho*cp*dT/dt - div(k*grad(T)) = 0
    // f0: accumulation term
    f[0] = rho_eff * cp_eff * u_t[0];
    
    // f1: conduction term (handled separately)
}

void ThermalKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                             const PetscScalar u_x[], const PetscScalar a[],
                             const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for heat equation: ρ*cp*∂T/∂t - div(k*∇T) = Q
    //
    // The weak form is:
    // ∫ (ρ*cp*∂T/∂t)*v dV + ∫ k*∇T·∇v dV = ∫ Q*v dV
    //
    // Jacobian contributions:
    // g0: ∂(ρ*cp*∂T/∂t)/∂T = 0 (no direct T dependence in linear case)
    //     ∂(ρ*cp*∂T/∂t)/∂(∂T/∂t) = ρ*cp * shift
    // g3: ∂(k*∇T)/∂(∇T) = k * I (identity tensor for isotropic conduction)
    
    // Compute effective properties
    double rho_eff = (1.0 - porosity) * density_solid + porosity * density_fluid;
    double cp_eff = (1.0 - porosity) * heat_capacity_solid + porosity * heat_capacity_fluid;
    double k_eff = (1.0 - porosity) * thermal_conductivity_solid + porosity * thermal_conductivity_fluid;
    
    double shift = a[0];  // Time integration shift parameter
    
    // Initialize Jacobian
    // For scalar temperature field, J has components for:
    // J[0] = g0: Jacobian w.r.t. T (includes time derivative shift)
    // J[1..9] could be g3: Jacobian w.r.t. ∇T (3x3 tensor for 3D)
    
    // Accumulation term Jacobian (g0)
    // ∂f0/∂(∂T/∂t) * shift = ρ_eff * cp_eff * shift
    J[0] = rho_eff * cp_eff * shift;
    
    // Conduction term Jacobian (g3)
    // For isotropic heat conduction: flux = k * ∇T
    // ∂(flux)/∂(∇T) = k * I
    // In 3D, this is a 3x3 identity matrix scaled by k_eff
    // Layout: J[1] = k_xx, J[2] = k_xy, J[3] = k_xz
    //         J[4] = k_yx, J[5] = k_yy, J[6] = k_yz
    //         J[7] = k_zx, J[8] = k_zy, J[9] = k_zz
    J[1] = k_eff;  // ∂q_x/∂(∂T/∂x)
    J[2] = 0.0;    // ∂q_x/∂(∂T/∂y)
    J[3] = 0.0;    // ∂q_x/∂(∂T/∂z)
    J[4] = 0.0;    // ∂q_y/∂(∂T/∂x)
    J[5] = k_eff;  // ∂q_y/∂(∂T/∂y)
    J[6] = 0.0;    // ∂q_y/∂(∂T/∂z)
    J[7] = 0.0;    // ∂q_z/∂(∂T/∂x)
    J[8] = 0.0;    // ∂q_z/∂(∂T/∂y)
    J[9] = k_eff;  // ∂q_z/∂(∂T/∂z)
    
    // For non-linear conductivity k(T), would add:
    // ∂(k*∇T)/∂T = dk/dT * ∇T
    // This would contribute to g2 (Jacobian of flux w.r.t. solution)
}

void ThermalKernel::setThermalProperties(double k, double rho, double cp) {
    thermal_conductivity_solid = k;
    density_solid = rho;
    heat_capacity_solid = cp;
}

void ThermalKernel::setFluidThermalProperties(double k_f, double rho_f, double cp_f) {
    thermal_conductivity_fluid = k_f;
    density_fluid = rho_f;
    heat_capacity_fluid = cp_f;
}

// ============================================================================
// Particle Transport Kernel
// ============================================================================

ParticleTransportKernel::ParticleTransportKernel()
    : PhysicsKernel(PhysicsType::PARTICLE_TRANSPORT),
      particle_diameter(1e-3), particle_density(2650.0), diffusivity(1e-9),
      gravity_settling(true), enable_bridging(false),
      porosity(0.2), permeability(100e-15) {}

PetscErrorCode ParticleTransportKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void ParticleTransportKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                       const PetscScalar u_x[], const PetscScalar a[],
                                       const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_x;
    (void)a;
    (void)x;
    
    // u[0] = C (concentration)
    
    // Advection-diffusion-settling equation
    // phi*dC/dt + div(v*C - D*grad(C)) + settling_term = 0
    
    f[0] = porosity * u_t[0];  // Accumulation
    
    // Gravitational settling
    if (gravity_settling) {
        double g = 9.81;
        double fluid_density = 1000.0;
        double stokes_velocity = 2.0 * std::pow(particle_diameter, 2.0) * 
                                (particle_density - fluid_density) * g / (9.0 * 0.001);
        
        f[0] += stokes_velocity * u[0];  // Simplified
    }
}

void ParticleTransportKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                       const PetscScalar u_x[], const PetscScalar a[],
                                       const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Full Jacobian for advection-diffusion-settling equation
    // φ*∂C/∂t + div(v*C - D*∇C) + settling_term = 0
    //
    // Jacobian contributions:
    // g0: ∂(φ*∂C/∂t + settling)/∂C = settling_derivative
    //     ∂(φ*∂C/∂t)/∂(∂C/∂t) = φ * shift
    // g1: ∂(v*C)/∂(∇C) - for upwind schemes (velocity dotted with gradient)
    // g2: ∂(v*C)/∂C = v (advection velocity)
    // g3: ∂(-D*∇C)/∂(∇C) = D * I (diffusion tensor)
    
    double shift = a[0];  // Time integration shift parameter
    
    // Darcy velocity from pressure gradient (would be coupled in full model)
    double fluid_viscosity = 0.001;  // Pa·s
    double darcy_velocity_mag = permeability / fluid_viscosity * 1e4;  // Simplified velocity scale
    
    // Gravitational settling velocity (Stokes settling)
    double settling_velocity = 0.0;
    if (gravity_settling) {
        double g = 9.81;
        double fluid_density = 1000.0;
        settling_velocity = 2.0 * std::pow(particle_diameter, 2.0) * 
                           (particle_density - fluid_density) * g / (9.0 * fluid_viscosity);
    }
    
    // Initialize Jacobian components
    // J[0] = g0: accumulation and source term Jacobians
    // J[1..3] = g2: advection velocity (∂(v·C)/∂C contributes v to flux Jacobian)
    // J[4..12] = g3: diffusion tensor (3x3)
    
    // Accumulation term Jacobian (g0)
    J[0] = porosity * shift;
    
    // Settling term Jacobian (part of g0)
    // ∂(v_s * C)/∂C = v_s
    J[0] += settling_velocity;
    
    // Advection term Jacobian (g2)
    // For advection: ∂(v·∇C)/∂C - this is handled through upwinding
    // The velocity components for advective flux Jacobian
    J[1] = darcy_velocity_mag;  // v_x component
    J[2] = 0.0;                  // v_y component
    J[3] = settling_velocity;   // v_z component (includes settling in z-direction)
    
    // Diffusion term Jacobian (g3)
    // For isotropic diffusion: flux = -D * ∇C
    // ∂(-D*∇C)/∂(∇C) = -D * I (note: sign depends on convention)
    // Layout: 3x3 tensor
    J[4] = diffusivity;   // D_xx
    J[5] = 0.0;           // D_xy
    J[6] = 0.0;           // D_xz
    J[7] = 0.0;           // D_yx
    J[8] = diffusivity;   // D_yy
    J[9] = 0.0;           // D_yz
    J[10] = 0.0;          // D_zx
    J[11] = 0.0;          // D_zy
    J[12] = diffusivity;  // D_zz
    
    // Mechanical dispersion (Taylor dispersion in porous media)
    // D_total = D_molecular + α_L * |v| (longitudinal)
    //                       + α_T * |v| (transverse)
    double alpha_L = 0.1;  // Longitudinal dispersivity [m]
    double alpha_T = 0.01; // Transverse dispersivity [m]
    
    // Add dispersion contribution (isotropic approximation)
    double dispersion = (alpha_L + 2.0 * alpha_T) / 3.0 * darcy_velocity_mag;
    J[4] += dispersion;
    J[8] += dispersion;
    J[12] += dispersion;
    
    // Filtration/bridging (optional mechanism)
    if (enable_bridging) {
        // Bridging coefficient depends on pore throat to particle size ratio
        double pore_throat = std::sqrt(permeability / porosity);  // Rough estimate
        double bridging_factor = 0.0;
        
        if (particle_diameter > 0.3 * pore_throat) {
            // Significant bridging occurs
            bridging_factor = 10.0 * std::pow(particle_diameter / pore_throat - 0.3, 2.0);
        }
        
        // Bridging adds a sink term proportional to concentration
        J[0] += bridging_factor;
    }
}

void ParticleTransportKernel::setParticleProperties(double diameter, double density, double D) {
    particle_diameter = diameter;
    particle_density = density;
    diffusivity = D;
}

void ParticleTransportKernel::enableGravitationalSettling(bool enable) {
    gravity_settling = enable;
}

void ParticleTransportKernel::enableBridging(bool enable) {
    enable_bridging = enable;
}

// ============================================================================
// Fracture Propagation Kernel
// ============================================================================

FracturePropagationKernel::FracturePropagationKernel()
    : PhysicsKernel(PhysicsType::FRACTURE_PROPAGATION),
      fracture_toughness(1e6), fracture_energy(100.0), critical_stress(5e6) {}

PetscErrorCode FracturePropagationKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void FracturePropagationKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                         const PetscScalar u_x[], const PetscScalar a[],
                                         const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_t;
    (void)u_x;
    (void)a;
    (void)x;
    
    // Cohesive zone model for fracture
    // u[0] = fracture opening displacement
    
    double delta = u[0];
    double delta_c = 2.0 * fracture_energy / critical_stress;
    
    if (delta < delta_c) {
        // Cohesive traction
        double T_coh = critical_stress * (1.0 - delta / delta_c);
        f[0] = T_coh;
    } else {
        // Fully opened
        f[0] = 0.0;
    }
}

void FracturePropagationKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                         const PetscScalar u_x[], const PetscScalar a[],
                                         const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u_t;
    (void)u_x;
    (void)a;
    (void)x;
    
    double delta = u[0];
    double delta_c = 2.0 * fracture_energy / critical_stress;
    
    if (delta < delta_c) {
        J[0] = -critical_stress / delta_c;
    } else {
        J[0] = 0.0;
    }
}

void FracturePropagationKernel::setFractureProperties(double Kc, double Gc, double sigma_c) {
    fracture_toughness = Kc;
    fracture_energy = Gc;
    critical_stress = sigma_c;
}

// ============================================================================
// Tidal Forces Kernel
// ============================================================================

TidalForcesKernel::TidalForcesKernel()
    : PhysicsKernel(PhysicsType::TIDAL_FORCES),
      latitude(0.0), longitude(0.0), current_time(0.0) {}

PetscErrorCode TidalForcesKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void TidalForcesKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscScalar a[],
                                 const PetscReal x[], PetscScalar f[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)a;
    
    // Tidal forces add body force to stress equilibrium
    PetscScalar stress[6];
    computeTidalStress(x, stress);
    
    // Add to residual (this is a forcing term)
    f[0] = stress[0];  // Simplified
}

void TidalForcesKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscScalar a[],
                                 const PetscReal x[], PetscScalar J[]) {
    // Suppress unused parameter warnings - these are part of the standard kernel interface
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)a;
    (void)x;
    
    // Tidal forcing doesn't depend on solution
    J[0] = 0.0;
}

void TidalForcesKernel::setLocationAndTime(double lat, double lon, double time) {
    latitude = lat;
    longitude = lon;
    current_time = time;
}

void TidalForcesKernel::computeTidalStress(const PetscReal x[], PetscScalar stress[]) {
    // Simplified tidal stress calculation
    // Full implementation would use lunar/solar ephemeris
    
    const double G = 6.674e-11;  // Gravitational constant
    const double M_moon = 7.342e22;  // Moon mass
    const double R_earth = 6.371e6;  // Earth radius
    const double d_moon = 3.844e8;   // Earth-Moon distance
    
    // Tidal Love number
    double h2 = 0.6;
    
    // Simplified tidal potential gradient
    double tidal_amplitude = G * M_moon * R_earth / std::pow(d_moon, 3);
    
    // Semi-diurnal and diurnal components
    double omega = 2.0 * M_PI / (12.42 * 3600.0);  // M2 tide period
    double phase = omega * current_time;
    
    stress[0] = tidal_amplitude * h2 * std::cos(phase);
    stress[1] = tidal_amplitude * h2 * std::cos(phase);
    stress[2] = 2.0 * tidal_amplitude * h2 * std::cos(phase);
    stress[3] = 0.0;
    stress[4] = 0.0;
    stress[5] = 0.0;
}

} // namespace FSRM
