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
    
    // Jacobian of accumulation term
    J[0] = porosity * compressibility * a[0];  // a is shift parameter
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
    
    // Simplified Jacobian
    double P = u[0];
    double Sw = u[1];
    double Sg = u[2];
    (void)Sw;  // Used in full implementation
    (void)Sg;  // Used in full implementation
    
    double rho_o = oilDensity(P, 0.0);
    double rho_w = waterDensity(P);
    double rho_g = gasDensity(P);
    
    // Diagonal blocks
    J[0] = porosity * rho_o * a[0];
    J[4] = porosity * rho_w * a[0];
    J[8] = porosity * rho_g * a[0];
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
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)x;
    
    // Simplified diagonal Jacobian
    for (int i = 0; i < nc; ++i) {
        J[i * (nc + 1) + i] = porosity * a[0];
    }
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
    
    double rho_eff = (1.0 - porosity) * density_solid + porosity * density_fluid;
    double cp_eff = (1.0 - porosity) * heat_capacity_solid + porosity * heat_capacity_fluid;
    
    J[0] = rho_eff * cp_eff * a[0];
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
    
    J[0] = porosity * a[0];
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
