#include "MaterialModel.hpp"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace FSRM {

// =============================================================================
// Utility functions
// =============================================================================

static double parseDouble(const std::map<std::string, std::string>& config,
                         const std::string& key, double default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        try {
            return std::stod(it->second);
        } catch (...) {
            return default_val;
        }
    }
    return default_val;
}

static std::string parseString(const std::map<std::string, std::string>& config,
                               const std::string& key, const std::string& default_val) {
    auto it = config.find(key);
    if (it != config.end() && !it->second.empty()) {
        return it->second;
    }
    return default_val;
}

// =============================================================================
// StressTensor Implementation
// =============================================================================

double StressTensor::vonMises() const {
    double s_mean = mean();
    double dev_xx = xx - s_mean;
    double dev_yy = yy - s_mean;
    double dev_zz = zz - s_mean;
    
    return std::sqrt(1.5 * (dev_xx*dev_xx + dev_yy*dev_yy + dev_zz*dev_zz +
                           2.0 * (xy*xy + xz*xz + yz*yz)));
}

std::array<double, 3> StressTensor::principalStresses() const {
    // Solve characteristic equation: det(σ - λI) = 0
    // λ³ - I₁λ² + I₂λ - I₃ = 0
    
    double I1 = xx + yy + zz;
    double I2 = xx*yy + yy*zz + zz*xx - xy*xy - yz*yz - xz*xz;
    double I3 = xx*yy*zz + 2.0*xy*yz*xz - xx*yz*yz - yy*xz*xz - zz*xy*xy;
    
    // Cardano's formula for cubic
    double p = I2 - I1*I1/3.0;
    double q = 2.0*I1*I1*I1/27.0 - I1*I2/3.0 + I3;
    
    double discriminant = q*q/4.0 + p*p*p/27.0;
    
    std::array<double, 3> principal;
    
    if (discriminant < 0) {
        // Three distinct real roots
        double r = std::sqrt(-p*p*p/27.0);
        double theta = std::acos(-q/(2.0*r));
        double rc = std::cbrt(r);
        
        principal[0] = 2.0*rc*std::cos(theta/3.0) + I1/3.0;
        principal[1] = 2.0*rc*std::cos((theta + 2.0*M_PI)/3.0) + I1/3.0;
        principal[2] = 2.0*rc*std::cos((theta + 4.0*M_PI)/3.0) + I1/3.0;
    } else {
        // One real root (degenerate case)
        double u = std::cbrt(-q/2.0 + std::sqrt(discriminant));
        double v = std::cbrt(-q/2.0 - std::sqrt(discriminant));
        principal[0] = u + v + I1/3.0;
        principal[1] = principal[0];
        principal[2] = principal[0];
    }
    
    // Sort: σ₁ >= σ₂ >= σ₃
    std::sort(principal.begin(), principal.end(), std::greater<double>());
    
    return principal;
}

double StressTensor::maxPrincipal() const {
    auto p = principalStresses();
    return p[0];
}

double StressTensor::minPrincipal() const {
    auto p = principalStresses();
    return p[2];
}

double StressTensor::maxShear() const {
    auto p = principalStresses();
    return (p[0] - p[2]) / 2.0;
}

// =============================================================================
// StrainTensor Implementation
// =============================================================================

double StrainTensor::deviatoric() const {
    double e_mean = volumetric() / 3.0;
    double dev_xx = xx - e_mean;
    double dev_yy = yy - e_mean;
    double dev_zz = zz - e_mean;
    
    return std::sqrt(2.0/3.0 * (dev_xx*dev_xx + dev_yy*dev_yy + dev_zz*dev_zz +
                               0.5 * (xy*xy + xz*xz + yz*yz)));
}

// =============================================================================
// LinearElasticMaterial Implementation
// =============================================================================

LinearElasticMaterial::LinearElasticMaterial()
    : MaterialModelBase(ConstitutiveModel::LINEAR_ELASTIC),
      youngs_modulus(10e9),
      poisson_ratio(0.25),
      density(2500.0),
      porosity(0.2),
      E_P_coeff(0.0),
      E_T_coeff(0.0),
      nu_P_coeff(0.0),
      nu_T_coeff(0.0),
      P_reference(1e5),
      T_reference(293.15) {}

StressTensor LinearElasticMaterial::computeStress(const StrainTensor& strain,
                                                   double P_pore,
                                                   double T) const {
    // Get effective moduli
    double E = youngs_modulus + E_P_coeff * (P_pore - P_reference) + 
               E_T_coeff * (T - T_reference);
    double nu = poisson_ratio + nu_P_coeff * (P_pore - P_reference) + 
                nu_T_coeff * (T - T_reference);
    
    // Clamp Poisson's ratio
    nu = std::max(-1.0, std::min(0.499, nu));
    
    // Lamé parameters
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    // Compute stress: σ = λ*tr(ε)*I + 2μ*ε
    double trace = strain.volumetric();
    
    StressTensor stress;
    stress.xx = lambda * trace + 2.0 * mu * strain.xx;
    stress.yy = lambda * trace + 2.0 * mu * strain.yy;
    stress.zz = lambda * trace + 2.0 * mu * strain.zz;
    stress.xy = mu * strain.xy;  // Note: engineering shear strain
    stress.xz = mu * strain.xz;
    stress.yz = mu * strain.yz;
    
    return stress;
}

std::array<double, 36> LinearElasticMaterial::getStiffnessMatrix(double P, double T) const {
    double E = youngs_modulus + E_P_coeff * (P - P_reference) + 
               E_T_coeff * (T - T_reference);
    double nu = poisson_ratio;
    
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    std::array<double, 36> C{};
    
    // Voigt notation: [xx, yy, zz, yz, xz, xy]
    C[0] = lambda + 2*mu;  C[1] = lambda;        C[2] = lambda;
    C[6] = lambda;         C[7] = lambda + 2*mu; C[8] = lambda;
    C[12] = lambda;        C[13] = lambda;       C[14] = lambda + 2*mu;
    C[21] = mu;  // yz-yz
    C[28] = mu;  // xz-xz
    C[35] = mu;  // xy-xy
    
    return C;
}

void LinearElasticMaterial::configure(const std::map<std::string, std::string>& config) {
    youngs_modulus = parseDouble(config, "youngs_modulus", 10e9);
    poisson_ratio = parseDouble(config, "poisson_ratio", 0.25);
    density = parseDouble(config, "density", 2500.0);
    porosity = parseDouble(config, "porosity", 0.2);
    
    E_P_coeff = parseDouble(config, "youngs_modulus_pressure_coeff", 0.0);
    E_T_coeff = parseDouble(config, "youngs_modulus_temperature_coeff", 0.0);
    P_reference = parseDouble(config, "reference_pressure", 1e5);
    T_reference = parseDouble(config, "reference_temperature", 293.15);
}

void LinearElasticMaterial::setElasticModuli(double E, double nu) {
    youngs_modulus = E;
    poisson_ratio = nu;
}

void LinearElasticMaterial::setLameParameters(double lambda, double mu) {
    youngs_modulus = mu * (3.0*lambda + 2.0*mu) / (lambda + mu);
    poisson_ratio = lambda / (2.0*(lambda + mu));
}

void LinearElasticMaterial::setDensity(double rho) {
    density = rho;
}

void LinearElasticMaterial::setPorosity(double phi) {
    porosity = phi;
}

double LinearElasticMaterial::getBulkModulus() const {
    return youngs_modulus / (3.0 * (1.0 - 2.0*poisson_ratio));
}

double LinearElasticMaterial::getShearModulus() const {
    return youngs_modulus / (2.0 * (1.0 + poisson_ratio));
}

double LinearElasticMaterial::getLambda() const {
    return youngs_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0*poisson_ratio));
}

double LinearElasticMaterial::getPWaveVelocity() const {
    double lambda = getLambda();
    double mu = getShearModulus();
    return std::sqrt((lambda + 2.0*mu) / density);
}

double LinearElasticMaterial::getSWaveVelocity() const {
    return std::sqrt(getShearModulus() / density);
}

// =============================================================================
// ViscoelasticMaterial Implementation
// =============================================================================

ViscoelasticMaterial::ViscoelasticMaterial(ViscoType type)
    : MaterialModelBase(ConstitutiveModel::VISCOELASTIC_MAXWELL),
      visco_type(type),
      E_instantaneous(10e9),
      E_long_term(5e9),
      poisson_ratio(0.25),
      relaxation_time(1e6),
      viscosity(1e19),
      density(2500.0),
      porosity(0.2) {
    
    switch (type) {
        case ViscoType::MAXWELL:
            model_type = ConstitutiveModel::VISCOELASTIC_MAXWELL;
            break;
        case ViscoType::KELVIN_VOIGT:
            model_type = ConstitutiveModel::VISCOELASTIC_KELVIN;
            break;
        case ViscoType::STANDARD_LINEAR_SOLID:
            model_type = ConstitutiveModel::VISCOELASTIC_SLS;
            break;
    }
}

StressTensor ViscoelasticMaterial::computeStress(const StrainTensor& strain,
                                                  double P_pore,
                                                  double T) const {
    (void)P_pore;  // Reserved for poroelastic coupling
    (void)T;  // Reserved for thermoelastic effects
    // For instantaneous response, use elastic moduli
    double E = E_instantaneous;
    double nu = poisson_ratio;
    
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    double trace = strain.volumetric();
    
    StressTensor stress;
    stress.xx = lambda * trace + 2.0 * mu * strain.xx;
    stress.yy = lambda * trace + 2.0 * mu * strain.yy;
    stress.zz = lambda * trace + 2.0 * mu * strain.zz;
    stress.xy = mu * strain.xy;
    stress.xz = mu * strain.xz;
    stress.yz = mu * strain.yz;
    
    return stress;
}

StressTensor ViscoelasticMaterial::computeStressRate(const StrainTensor& strain,
                                                      const StrainTensor& strain_rate,
                                                      const StressTensor& current_stress,
                                                      double dt) const {
    (void)strain;  // Current state tracked in stress history
    // Compute stress increment
    StressTensor stress_rate;
    
    double G = E_instantaneous / (2.0 * (1.0 + poisson_ratio));
    
    switch (visco_type) {
        case ViscoType::MAXWELL: {
            // Maxwell: dσ/dt + σ/τ = G*dε/dt
            // Backward Euler: σ_new = (σ_old + G*dε*dt) / (1 + dt/τ)
            double factor = 1.0 / (1.0 + dt / relaxation_time);
            stress_rate.xx = (2.0*G*strain_rate.xx - current_stress.xx/relaxation_time) * factor;
            stress_rate.yy = (2.0*G*strain_rate.yy - current_stress.yy/relaxation_time) * factor;
            stress_rate.zz = (2.0*G*strain_rate.zz - current_stress.zz/relaxation_time) * factor;
            stress_rate.xy = (G*strain_rate.xy - current_stress.xy/relaxation_time) * factor;
            stress_rate.xz = (G*strain_rate.xz - current_stress.xz/relaxation_time) * factor;
            stress_rate.yz = (G*strain_rate.yz - current_stress.yz/relaxation_time) * factor;
            break;
        }
        case ViscoType::KELVIN_VOIGT: {
            // Kelvin-Voigt: σ = G*ε + η*dε/dt
            stress_rate.xx = 2.0*G*strain_rate.xx + viscosity*strain_rate.xx;
            stress_rate.yy = 2.0*G*strain_rate.yy + viscosity*strain_rate.yy;
            stress_rate.zz = 2.0*G*strain_rate.zz + viscosity*strain_rate.zz;
            stress_rate.xy = G*strain_rate.xy + viscosity*strain_rate.xy;
            stress_rate.xz = G*strain_rate.xz + viscosity*strain_rate.xz;
            stress_rate.yz = G*strain_rate.yz + viscosity*strain_rate.yz;
            break;
        }
        case ViscoType::STANDARD_LINEAR_SOLID: {
            // SLS: G_R + G_1*exp(-t/τ) where G_R is long-term modulus
            double G_inf = E_long_term / (2.0 * (1.0 + poisson_ratio));
            double G_1 = G - G_inf;
            double exp_factor = std::exp(-dt / relaxation_time);
            
            // More complex update needed for SLS
            stress_rate.xx = G_inf * strain_rate.xx + G_1 * exp_factor * strain_rate.xx;
            stress_rate.yy = G_inf * strain_rate.yy + G_1 * exp_factor * strain_rate.yy;
            stress_rate.zz = G_inf * strain_rate.zz + G_1 * exp_factor * strain_rate.zz;
            stress_rate.xy = G_inf * strain_rate.xy + G_1 * exp_factor * strain_rate.xy;
            stress_rate.xz = G_inf * strain_rate.xz + G_1 * exp_factor * strain_rate.xz;
            stress_rate.yz = G_inf * strain_rate.yz + G_1 * exp_factor * strain_rate.yz;
            break;
        }
    }
    
    return stress_rate;
}

std::array<double, 36> ViscoelasticMaterial::getStiffnessMatrix(double P, double T) const {
    (void)P;  // Reserved for pressure-dependent moduli
    (void)T;  // Reserved for temperature-dependent moduli
    // Return instantaneous stiffness
    double E = E_instantaneous;
    double nu = poisson_ratio;
    
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    double mu = E / (2.0 * (1.0 + nu));
    
    std::array<double, 36> C{};
    C[0] = lambda + 2*mu;  C[1] = lambda;        C[2] = lambda;
    C[6] = lambda;         C[7] = lambda + 2*mu; C[8] = lambda;
    C[12] = lambda;        C[13] = lambda;       C[14] = lambda + 2*mu;
    C[21] = mu;
    C[28] = mu;
    C[35] = mu;
    
    return C;
}

void ViscoelasticMaterial::configure(const std::map<std::string, std::string>& config) {
    E_instantaneous = parseDouble(config, "youngs_modulus", 10e9);
    E_long_term = parseDouble(config, "long_term_modulus", E_instantaneous * 0.5);
    poisson_ratio = parseDouble(config, "poisson_ratio", 0.25);
    relaxation_time = parseDouble(config, "relaxation_time", 1e6);
    viscosity = parseDouble(config, "viscosity", 1e19);
    density = parseDouble(config, "density", 2500.0);
    porosity = parseDouble(config, "porosity", 0.2);
    
    std::string type = parseString(config, "viscoelastic_type", "MAXWELL");
    if (type == "MAXWELL") visco_type = ViscoType::MAXWELL;
    else if (type == "KELVIN" || type == "KELVIN_VOIGT") visco_type = ViscoType::KELVIN_VOIGT;
    else if (type == "SLS" || type == "STANDARD_LINEAR_SOLID") visco_type = ViscoType::STANDARD_LINEAR_SOLID;
}

void ViscoelasticMaterial::setInstantaneousModuli(double E, double nu) {
    E_instantaneous = E;
    poisson_ratio = nu;
}

void ViscoelasticMaterial::setRelaxationTime(double tau) {
    relaxation_time = tau;
}

void ViscoelasticMaterial::setViscosity(double eta) {
    viscosity = eta;
}

void ViscoelasticMaterial::setLongTermModulus(double E_inf) {
    E_long_term = E_inf;
}

double ViscoelasticMaterial::getRelaxationModulus(double t) const {
    double G0 = E_instantaneous / (2.0 * (1.0 + poisson_ratio));
    
    switch (visco_type) {
        case ViscoType::MAXWELL:
            return G0 * std::exp(-t / relaxation_time);
        case ViscoType::KELVIN_VOIGT:
            return G0;  // Constant
        case ViscoType::STANDARD_LINEAR_SOLID: {
            double G_inf = E_long_term / (2.0 * (1.0 + poisson_ratio));
            return G_inf + (G0 - G_inf) * std::exp(-t / relaxation_time);
        }
    }
    return G0;
}

double ViscoelasticMaterial::getCreepCompliance(double t) const {
    double J0 = 2.0 * (1.0 + poisson_ratio) / E_instantaneous;
    
    switch (visco_type) {
        case ViscoType::MAXWELL:
            return J0 * (1.0 + t / relaxation_time);
        case ViscoType::KELVIN_VOIGT:
            return J0 * (1.0 - std::exp(-t / relaxation_time));
        case ViscoType::STANDARD_LINEAR_SOLID: {
            double J_inf = 2.0 * (1.0 + poisson_ratio) / E_long_term;
            return J_inf - (J_inf - J0) * std::exp(-t / relaxation_time);
        }
    }
    return J0;
}

// =============================================================================
// PoroelasticMaterial Implementation
// =============================================================================

PoroelasticMaterial::PoroelasticMaterial()
    : MaterialModelBase(ConstitutiveModel::POROELASTIC),
      youngs_modulus(10e9),
      poisson_ratio(0.25),
      bulk_modulus_drained(0.0),
      shear_modulus(0.0),
      biot_coefficient(1.0),
      biot_modulus(1e10),
      density_solid(2650.0),
      density_fluid(1000.0),
      bulk_modulus_solid(35e9),
      bulk_modulus_fluid(2.2e9),
      fluid_viscosity(0.001),
      porosity(0.2),
      permeability(100e-15) {
    
    // Compute derived quantities
    bulk_modulus_drained = youngs_modulus / (3.0 * (1.0 - 2.0*poisson_ratio));
    shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
}

StressTensor PoroelasticMaterial::computeStress(const StrainTensor& strain,
                                                 double P_pore,
                                                 double T) const {
    (void)T;  // Reserved for thermoporoelastic effects
    // Total stress = effective stress - α*P*I
    StressTensor effective = computeEffectiveStress(strain);
    return computeTotalStress(effective, P_pore);
}

StressTensor PoroelasticMaterial::computeEffectiveStress(const StrainTensor& strain) const {
    double lambda = bulk_modulus_drained - 2.0*shear_modulus/3.0;
    double mu = shear_modulus;
    
    double trace = strain.volumetric();
    
    StressTensor stress;
    stress.xx = lambda * trace + 2.0 * mu * strain.xx;
    stress.yy = lambda * trace + 2.0 * mu * strain.yy;
    stress.zz = lambda * trace + 2.0 * mu * strain.zz;
    stress.xy = mu * strain.xy;
    stress.xz = mu * strain.xz;
    stress.yz = mu * strain.yz;
    
    return stress;
}

StressTensor PoroelasticMaterial::computeTotalStress(const StressTensor& effective_stress,
                                                      double P_pore) const {
    StressTensor total = effective_stress;
    total.xx -= biot_coefficient * P_pore;
    total.yy -= biot_coefficient * P_pore;
    total.zz -= biot_coefficient * P_pore;
    return total;
}

std::array<double, 36> PoroelasticMaterial::getStiffnessMatrix(double P, double T) const {
    (void)P;  // Reserved for pressure-dependent moduli
    (void)T;  // Reserved for temperature-dependent moduli
    double lambda = bulk_modulus_drained - 2.0*shear_modulus/3.0;
    double mu = shear_modulus;
    
    std::array<double, 36> C{};
    C[0] = lambda + 2*mu;  C[1] = lambda;        C[2] = lambda;
    C[6] = lambda;         C[7] = lambda + 2*mu; C[8] = lambda;
    C[12] = lambda;        C[13] = lambda;       C[14] = lambda + 2*mu;
    C[21] = mu;
    C[28] = mu;
    C[35] = mu;
    
    return C;
}

void PoroelasticMaterial::configure(const std::map<std::string, std::string>& config) {
    youngs_modulus = parseDouble(config, "youngs_modulus", 10e9);
    poisson_ratio = parseDouble(config, "poisson_ratio", 0.25);
    biot_coefficient = parseDouble(config, "biot_coefficient", 1.0);
    biot_modulus = parseDouble(config, "biot_modulus", 1e10);
    density_solid = parseDouble(config, "density", 2650.0);
    density_fluid = parseDouble(config, "fluid_density", 1000.0);
    bulk_modulus_solid = parseDouble(config, "grain_bulk_modulus", 35e9);
    bulk_modulus_fluid = parseDouble(config, "fluid_bulk_modulus", 2.2e9);
    fluid_viscosity = parseDouble(config, "fluid_viscosity", 0.001);
    porosity = parseDouble(config, "porosity", 0.2);
    permeability = parseDouble(config, "permeability", 100.0) * 1e-15;  // mD to m²
    
    // Update derived quantities
    bulk_modulus_drained = youngs_modulus / (3.0 * (1.0 - 2.0*poisson_ratio));
    shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
}

double PoroelasticMaterial::getSkemptonCoefficient() const {
    // B = α*M / (K_d + α²*M)
    return biot_coefficient * biot_modulus / 
           (bulk_modulus_drained + biot_coefficient * biot_coefficient * biot_modulus);
}

double PoroelasticMaterial::getUndrainedBulkModulus() const {
    // K_u = K_d + α²*M
    return bulk_modulus_drained + biot_coefficient * biot_coefficient * biot_modulus;
}

double PoroelasticMaterial::getStorageCoefficient() const {
    // S = 1/M + α*(α-φ)/K_d
    return 1.0/biot_modulus + 
           biot_coefficient * (biot_coefficient - porosity) / bulk_modulus_drained;
}

double PoroelasticMaterial::getFastPWaveVelocity() const {
    double K_u = getUndrainedBulkModulus();
    double rho = (1.0 - porosity) * density_solid + porosity * density_fluid;
    return std::sqrt((K_u + 4.0*shear_modulus/3.0) / rho);
}

double PoroelasticMaterial::getSlowPWaveVelocity() const {
    // Biot slow wave (highly attenuated, diffusive at low frequency)
    double rho_f = density_fluid;
    double mobility = permeability / fluid_viscosity;
    
    // Characteristic frequency
    double omega_c = porosity / (mobility * rho_f);
    
    // At low frequency (quasi-static), slow wave is diffusive
    // Velocity estimate
    return std::sqrt(biot_modulus * porosity / rho_f) * 0.1;  // Order of magnitude
}

double PoroelasticMaterial::getSWaveVelocity() const {
    double rho = (1.0 - porosity) * density_solid + porosity * density_fluid;
    return std::sqrt(shear_modulus / rho);
}

void PoroelasticMaterial::setDrainedModuli(double E, double nu) {
    youngs_modulus = E;
    poisson_ratio = nu;
    bulk_modulus_drained = E / (3.0 * (1.0 - 2.0*nu));
    shear_modulus = E / (2.0 * (1.0 + nu));
}

void PoroelasticMaterial::setBiotParameters(double alpha, double M) {
    biot_coefficient = alpha;
    biot_modulus = M;
}

void PoroelasticMaterial::setFluidProperties(double Kf, double rho_f, double mu_f) {
    bulk_modulus_fluid = Kf;
    density_fluid = rho_f;
    fluid_viscosity = mu_f;
}

void PoroelasticMaterial::setPermeability(double k) {
    permeability = k;
}

void PoroelasticMaterial::setPorosity(double phi) {
    porosity = phi;
}

// =============================================================================
// ElastoplasticMaterial Implementation
// =============================================================================

ElastoplasticMaterial::ElastoplasticMaterial(FailureCriterion crit)
    : MaterialModelBase(ConstitutiveModel::ELASTOPLASTIC_MC),
      criterion(crit),
      youngs_modulus(10e9),
      poisson_ratio(0.25),
      density(2500.0),
      porosity(0.2),
      cohesion(1e6),
      friction_angle(30.0 * M_PI / 180.0),
      dilation_angle(0.0),
      tensile_strength(0.5e6),
      dp_d(0.0),
      dp_k(0.0),
      hardening_modulus(0.0),
      initial_yield_stress(1e7),
      accumulated_plastic_strain(0.0) {
    
    // Initialize Drucker-Prager parameters from Mohr-Coulomb
    double sin_phi = std::sin(friction_angle);
    double cos_phi = std::cos(friction_angle);
    
    // Outer circle matching (plane strain)
    dp_k = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    dp_d = 6.0 * cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
}

StressTensor ElastoplasticMaterial::computeStress(const StrainTensor& strain,
                                                   double P_pore,
                                                   double T) const {
    (void)P_pore;  // Reserved for poroelastic coupling
    (void)T;  // Reserved for thermoelastic effects
    // Elastic trial stress
    double lambda = youngs_modulus * poisson_ratio / 
                   ((1.0 + poisson_ratio) * (1.0 - 2.0*poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    
    double trace = strain.volumetric();
    
    StressTensor trial;
    trial.xx = lambda * trace + 2.0 * mu * strain.xx;
    trial.yy = lambda * trace + 2.0 * mu * strain.yy;
    trial.zz = lambda * trace + 2.0 * mu * strain.zz;
    trial.xy = mu * strain.xy;
    trial.xz = mu * strain.xz;
    trial.yz = mu * strain.yz;
    
    // Check yield
    double f = evaluateYieldFunction(trial);
    
    if (f <= 0.0) {
        return trial;  // Elastic
    }
    
    // Plastic correction (return mapping)
    double eps_p_eq = accumulated_plastic_strain;
    StressTensor correction = computePlasticCorrection(trial, eps_p_eq);
    
    // Updated stress
    StressTensor stress;
    stress.xx = trial.xx - correction.xx;
    stress.yy = trial.yy - correction.yy;
    stress.zz = trial.zz - correction.zz;
    stress.xy = trial.xy - correction.xy;
    stress.xz = trial.xz - correction.xz;
    stress.yz = trial.yz - correction.yz;
    
    return stress;
}

StressTensor ElastoplasticMaterial::computePlasticCorrection(const StressTensor& trial_stress,
                                                              double& plastic_strain_eq) const {
    StressTensor correction{0, 0, 0, 0, 0, 0};
    
    double f = evaluateYieldFunction(trial_stress);
    if (f <= 0.0) return correction;
    
    double G = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    double K = youngs_modulus / (3.0 * (1.0 - 2.0*poisson_ratio));
    
    if (criterion == FailureCriterion::DRUCKER_PRAGER || 
        criterion == FailureCriterion::MOHR_COULOMB) {
        
        // J2 and I1 invariants
        double I1 = trial_stress.xx + trial_stress.yy + trial_stress.zz;
        double s_mean = I1 / 3.0;
        
        double s_xx = trial_stress.xx - s_mean;
        double s_yy = trial_stress.yy - s_mean;
        double s_zz = trial_stress.zz - s_mean;
        
        double J2 = 0.5 * (s_xx*s_xx + s_yy*s_yy + s_zz*s_zz) +
                    trial_stress.xy*trial_stress.xy + 
                    trial_stress.xz*trial_stress.xz + 
                    trial_stress.yz*trial_stress.yz;
        
        double sqrtJ2 = std::sqrt(J2);
        
        // Plastic multiplier (radial return)
        double denom = G + K * dp_k * dp_k + hardening_modulus;
        double dlambda = f / denom;
        
        // Update plastic strain
        plastic_strain_eq += dlambda * std::sqrt(2.0/3.0);
        
        // Correction (deviatoric + volumetric)
        double factor = G * dlambda / sqrtJ2;
        correction.xx = factor * s_xx + K * dp_k * dlambda;
        correction.yy = factor * s_yy + K * dp_k * dlambda;
        correction.zz = factor * s_zz + K * dp_k * dlambda;
        correction.xy = factor * trial_stress.xy;
        correction.xz = factor * trial_stress.xz;
        correction.yz = factor * trial_stress.yz;
    }
    
    return correction;
}

std::array<double, 36> ElastoplasticMaterial::getStiffnessMatrix(double P, double T) const {
    (void)P;  // Reserved for pressure-dependent moduli
    (void)T;  // Reserved for temperature-dependent moduli
    // Return elastic stiffness (tangent would be different during plastic flow)
    double lambda = youngs_modulus * poisson_ratio / 
                   ((1.0 + poisson_ratio) * (1.0 - 2.0*poisson_ratio));
    double mu = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    
    std::array<double, 36> C{};
    C[0] = lambda + 2*mu;  C[1] = lambda;        C[2] = lambda;
    C[6] = lambda;         C[7] = lambda + 2*mu; C[8] = lambda;
    C[12] = lambda;        C[13] = lambda;       C[14] = lambda + 2*mu;
    C[21] = mu;
    C[28] = mu;
    C[35] = mu;
    
    return C;
}

void ElastoplasticMaterial::configure(const std::map<std::string, std::string>& config) {
    youngs_modulus = parseDouble(config, "youngs_modulus", 10e9);
    poisson_ratio = parseDouble(config, "poisson_ratio", 0.25);
    density = parseDouble(config, "density", 2500.0);
    porosity = parseDouble(config, "porosity", 0.2);
    
    cohesion = parseDouble(config, "cohesion", 1e6);
    friction_angle = parseDouble(config, "friction_angle", 30.0) * M_PI / 180.0;
    dilation_angle = parseDouble(config, "dilation_angle", 0.0) * M_PI / 180.0;
    tensile_strength = parseDouble(config, "tensile_strength", 0.5e6);
    
    hardening_modulus = parseDouble(config, "hardening_modulus", 0.0);
    
    // Update DP parameters
    double sin_phi = std::sin(friction_angle);
    double cos_phi = std::cos(friction_angle);
    dp_k = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    dp_d = 6.0 * cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    
    std::string crit = parseString(config, "failure_criterion", "MOHR_COULOMB");
    if (crit == "MOHR_COULOMB" || crit == "MC") criterion = FailureCriterion::MOHR_COULOMB;
    else if (crit == "DRUCKER_PRAGER" || crit == "DP") criterion = FailureCriterion::DRUCKER_PRAGER;
    else if (crit == "VON_MISES") criterion = FailureCriterion::VON_MISES;
}

double ElastoplasticMaterial::evaluateYieldFunction(const StressTensor& stress) const {
    double I1 = stress.xx + stress.yy + stress.zz;
    double s_mean = I1 / 3.0;
    
    double s_xx = stress.xx - s_mean;
    double s_yy = stress.yy - s_mean;
    double s_zz = stress.zz - s_mean;
    
    double J2 = 0.5 * (s_xx*s_xx + s_yy*s_yy + s_zz*s_zz) +
                stress.xy*stress.xy + stress.xz*stress.xz + stress.yz*stress.yz;
    
    switch (criterion) {
        case FailureCriterion::DRUCKER_PRAGER:
        case FailureCriterion::MOHR_COULOMB:
            // f = sqrt(J2) + k*I1 - d
            return std::sqrt(J2) + dp_k * I1 - dp_d;
            
        case FailureCriterion::VON_MISES:
            // f = sqrt(3*J2) - sigma_y
            return std::sqrt(3.0 * J2) - initial_yield_stress;
            
        default:
            return std::sqrt(J2) + dp_k * I1 - dp_d;
    }
}

bool ElastoplasticMaterial::isYielding(const StressTensor& stress) const {
    return evaluateYieldFunction(stress) > 0.0;
}

void ElastoplasticMaterial::setMohrCoulombParameters(double c, double phi, double psi) {
    cohesion = c;
    friction_angle = phi;
    dilation_angle = psi;
    
    double sin_phi = std::sin(friction_angle);
    double cos_phi = std::cos(friction_angle);
    dp_k = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
    dp_d = 6.0 * cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
}

void ElastoplasticMaterial::setDruckerPragerParameters(double d, double k) {
    dp_d = d;
    dp_k = k;
}

void ElastoplasticMaterial::setTensileStrength(double sigma_t) {
    tensile_strength = sigma_t;
}

void ElastoplasticMaterial::setHardeningParameters(double H, double sigma_y) {
    hardening_modulus = H;
    initial_yield_stress = sigma_y;
}

// =============================================================================
// AnisotropicMaterial Implementation
// =============================================================================

AnisotropicMaterial::AnisotropicMaterial(bool vertical_ti)
    : MaterialModelBase(vertical_ti ? ConstitutiveModel::ANISOTROPIC_VTI : 
                                     ConstitutiveModel::ANISOTROPIC_HTI),
      is_VTI(vertical_ti),
      C11(0), C13(0), C33(0), C44(0), C66(0),
      E_vertical(10e9),
      E_horizontal(15e9),
      G_vertical(4e9),
      nu_vh(0.25),
      nu_hh(0.20),
      epsilon(0.1),
      delta(0.05),
      gamma(0.08),
      density(2500.0),
      porosity(0.1) {
    
    // Initialize stiffness components from engineering constants
    double nu_hv = nu_vh * E_horizontal / E_vertical;
    double D = (1.0 - nu_hh * nu_hh) * E_vertical / E_horizontal - 2.0 * nu_vh * nu_hv - nu_vh * nu_vh;
    
    C33 = E_vertical * (1.0 - nu_hh * nu_hh * E_vertical / E_horizontal) / D;
    C11 = E_horizontal * (1.0 - nu_vh * nu_hv) / D;
    C13 = E_horizontal * nu_vh * (1.0 + nu_hh) / D;
    C44 = G_vertical;
    C66 = E_horizontal / (2.0 * (1.0 + nu_hh));
}

StressTensor AnisotropicMaterial::computeStress(const StrainTensor& strain,
                                                 double P_pore,
                                                 double T) const {
    (void)P_pore;  // Reserved for poroelastic coupling
    (void)T;  // Reserved for thermoelastic effects
    StressTensor stress;
    
    if (is_VTI) {
        // VTI: Symmetry axis is z (vertical)
        stress.xx = C11 * strain.xx + (C11 - 2.0*C66) * strain.yy + C13 * strain.zz;
        stress.yy = (C11 - 2.0*C66) * strain.xx + C11 * strain.yy + C13 * strain.zz;
        stress.zz = C13 * strain.xx + C13 * strain.yy + C33 * strain.zz;
        stress.yz = C44 * strain.yz;
        stress.xz = C44 * strain.xz;
        stress.xy = C66 * strain.xy;
    } else {
        // HTI: Symmetry axis is x (horizontal)
        stress.xx = C33 * strain.xx + C13 * strain.yy + C13 * strain.zz;
        stress.yy = C13 * strain.xx + C11 * strain.yy + (C11 - 2.0*C66) * strain.zz;
        stress.zz = C13 * strain.xx + (C11 - 2.0*C66) * strain.yy + C11 * strain.zz;
        stress.yz = C66 * strain.yz;
        stress.xz = C44 * strain.xz;
        stress.xy = C44 * strain.xy;
    }
    
    return stress;
}

std::array<double, 36> AnisotropicMaterial::getStiffnessMatrix(double P, double T) const {
    (void)P;  // Reserved for pressure-dependent moduli
    (void)T;  // Reserved for temperature-dependent moduli
    std::array<double, 36> C{};
    
    if (is_VTI) {
        C[0] = C11;      C[1] = C11-2*C66; C[2] = C13;
        C[6] = C11-2*C66; C[7] = C11;      C[8] = C13;
        C[12] = C13;      C[13] = C13;     C[14] = C33;
        C[21] = C44;
        C[28] = C44;
        C[35] = C66;
    } else {
        C[0] = C33;  C[1] = C13;       C[2] = C13;
        C[6] = C13;  C[7] = C11;       C[8] = C11-2*C66;
        C[12] = C13; C[13] = C11-2*C66; C[14] = C11;
        C[21] = C66;
        C[28] = C44;
        C[35] = C44;
    }
    
    return C;
}

void AnisotropicMaterial::configure(const std::map<std::string, std::string>& config) {
    E_vertical = parseDouble(config, "youngs_modulus_vertical", 10e9);
    E_horizontal = parseDouble(config, "youngs_modulus_horizontal", 15e9);
    G_vertical = parseDouble(config, "shear_modulus_vertical", 4e9);
    nu_vh = parseDouble(config, "poisson_ratio_vh", 0.25);
    nu_hh = parseDouble(config, "poisson_ratio_hh", 0.20);
    density = parseDouble(config, "density", 2500.0);
    porosity = parseDouble(config, "porosity", 0.1);
    
    // Thomsen parameters
    epsilon = parseDouble(config, "thomsen_epsilon", 0.1);
    delta = parseDouble(config, "thomsen_delta", 0.05);
    gamma = parseDouble(config, "thomsen_gamma", 0.08);
    
    // Recalculate stiffness
    double nu_hv = nu_vh * E_horizontal / E_vertical;
    double D = (1.0 - nu_hh * nu_hh) * E_vertical / E_horizontal - 2.0 * nu_vh * nu_hv - nu_vh * nu_vh;
    
    C33 = E_vertical * (1.0 - nu_hh * nu_hh * E_vertical / E_horizontal) / D;
    C11 = E_horizontal * (1.0 - nu_vh * nu_hv) / D;
    C13 = E_horizontal * nu_vh * (1.0 + nu_hh) / D;
    C44 = G_vertical;
    C66 = E_horizontal / (2.0 * (1.0 + nu_hh));
}

void AnisotropicMaterial::setThomsenParameters(double eps, double del, double gam) {
    epsilon = eps;
    delta = del;
    gamma = gam;
    
    // Update stiffness from Thomsen parameters
    double VP0 = std::sqrt(C33 / density);
    double VS0 = std::sqrt(C44 / density);
    
    C11 = C33 * (1.0 + 2.0 * epsilon);
    C66 = C44 * (1.0 + 2.0 * gamma);
    // C13 from delta (more complex relationship)
}

void AnisotropicMaterial::setVerticalModuli(double Ev, double Gv, double nuv) {
    E_vertical = Ev;
    G_vertical = Gv;
    nu_vh = nuv;
    C44 = Gv;
}

void AnisotropicMaterial::setHorizontalModuli(double Eh, double Gh, double nuh) {
    E_horizontal = Eh;
    nu_hh = nuh;
    C66 = Gh;
}

double AnisotropicMaterial::getVP0() const {
    return std::sqrt(C33 / density);
}

double AnisotropicMaterial::getVS0() const {
    return std::sqrt(C44 / density);
}

double AnisotropicMaterial::getVPh() const {
    return std::sqrt(C11 / density);
}

double AnisotropicMaterial::getVSh() const {
    return std::sqrt(C66 / density);
}

// =============================================================================
// ThermalProperties Implementation
// =============================================================================

void ThermalProperties::configure(const std::map<std::string, std::string>& config) {
    thermal_conductivity = parseDouble(config, "thermal_conductivity", 2.5);
    specific_heat = parseDouble(config, "heat_capacity", 900.0);
    thermal_expansion = parseDouble(config, "thermal_expansion", 1e-5);
    reference_temperature = parseDouble(config, "reference_temperature", 293.15);
    k_T_coeff = parseDouble(config, "conductivity_temperature_coeff", 0.0);
    cp_T_coeff = parseDouble(config, "specific_heat_temperature_coeff", 0.0);
}

double ThermalProperties::getConductivity(double T) const {
    return thermal_conductivity + k_T_coeff * (T - reference_temperature);
}

double ThermalProperties::getSpecificHeat(double T) const {
    return specific_heat + cp_T_coeff * (T - reference_temperature);
}

StressTensor ThermalProperties::getThermalStress(double E, double nu, double dT) const {
    // Thermal stress for constrained heating: σ = -E*α*ΔT / (1-2ν)
    double stress_val = -E * thermal_expansion * dT / (1.0 - 2.0*nu);
    
    StressTensor stress;
    stress.xx = stress_val;
    stress.yy = stress_val;
    stress.zz = stress_val;
    stress.xy = 0.0;
    stress.xz = 0.0;
    stress.yz = 0.0;
    
    return stress;
}

// =============================================================================
// PermeabilityModel Implementation
// =============================================================================

PermeabilityModel::PermeabilityModel(PermeabilityType t)
    : type(t),
      k_initial(100e-15),
      stress_coefficient(1e-8),
      porosity_exponent(3.0),
      grain_diameter(0.1e-3),
      strain_sensitivity(1.0),
      recovery_time(100.0),
      k_min(1e-18),
      k_max(1e-10) {}

void PermeabilityModel::configure(const std::map<std::string, std::string>& config) {
    k_initial = parseDouble(config, "permeability", 100.0) * 1e-15;  // mD to m²
    stress_coefficient = parseDouble(config, "permeability_stress_coeff", 1e-8);
    porosity_exponent = parseDouble(config, "permeability_porosity_exp", 3.0);
    grain_diameter = parseDouble(config, "grain_diameter", 0.1e-3);
    strain_sensitivity = parseDouble(config, "permeability_strain_sensitivity", 1.0);
    recovery_time = parseDouble(config, "permeability_recovery_time", 100.0);
    k_min = parseDouble(config, "permeability_min", 0.001) * 1e-15;
    k_max = parseDouble(config, "permeability_max", 10000.0) * 1e-15;
    
    std::string t = parseString(config, "permeability_model", "CONSTANT");
    if (t == "CONSTANT") type = PermeabilityType::CONSTANT;
    else if (t == "KOZENY_CARMAN") type = PermeabilityType::KOZENY_CARMAN;
    else if (t == "EXPONENTIAL") type = PermeabilityType::EXPONENTIAL;
    else if (t == "POWER_LAW") type = PermeabilityType::POWER_LAW;
    else if (t == "CUBIC_LAW") type = PermeabilityType::CUBIC_LAW;
    else if (t == "DYNAMIC") type = PermeabilityType::DYNAMIC;
}

double PermeabilityModel::getPermeability(double effective_stress, double porosity,
                                          double aperture) const {
    double k = k_initial;
    
    switch (type) {
        case PermeabilityType::CONSTANT:
            k = k_initial;
            break;
            
        case PermeabilityType::KOZENY_CARMAN:
            // k = d² * φ³ / (180 * (1-φ)²)
            k = grain_diameter * grain_diameter * std::pow(porosity, porosity_exponent) /
                (180.0 * std::pow(1.0 - porosity, 2.0));
            break;
            
        case PermeabilityType::EXPONENTIAL:
            // k = k0 * exp(-c * σ')
            k = k_initial * std::exp(-stress_coefficient * effective_stress);
            break;
            
        case PermeabilityType::POWER_LAW:
            // k = k0 * (1 + c*σ')^(-n)
            k = k_initial * std::pow(1.0 + stress_coefficient * effective_stress, -porosity_exponent);
            break;
            
        case PermeabilityType::CUBIC_LAW:
            // For fractures: k = w²/12
            if (aperture > 0.0) {
                k = aperture * aperture / 12.0;
            }
            break;
            
        case PermeabilityType::DYNAMIC:
            k = k_initial;  // Base value, modified by getDynamicPermeability
            break;
    }
    
    return std::max(k_min, std::min(k_max, k));
}

double PermeabilityModel::getDynamicPermeability(double k0, double volumetric_strain,
                                                  double shear_strain, double dt) const {
    // Dynamic permeability change from transient deformation
    double dk = strain_sensitivity * (std::abs(volumetric_strain) + 0.5 * std::abs(shear_strain));
    double k_enhanced = k0 * (1.0 + dk);
    
    // Recovery toward initial value
    if (recovery_time > 0.0) {
        double decay = std::exp(-dt / recovery_time);
        k_enhanced = k0 + (k_enhanced - k0) * decay;
    }
    
    return std::max(k_min, std::min(k_max, k_enhanced));
}

void PermeabilityModel::setInitialPermeability(double k0) {
    k_initial = k0;
}

void PermeabilityModel::setStressSensitivity(double coeff) {
    stress_coefficient = coeff;
}

void PermeabilityModel::setPorosityExponent(double n) {
    porosity_exponent = n;
}

void PermeabilityModel::setRecoveryTime(double tau) {
    recovery_time = tau;
}

void PermeabilityModel::setBounds(double kmin, double kmax) {
    k_min = kmin;
    k_max = kmax;
}

// =============================================================================
// FractureProperties Implementation
// =============================================================================

void FractureProperties::configure(const std::map<std::string, std::string>& config) {
    toughness_KIc = parseDouble(config, "fracture_toughness", 1e6);
    toughness_KIIc = parseDouble(config, "fracture_toughness_mode2", toughness_KIc * 1.5);
    toughness_KIIIc = parseDouble(config, "fracture_toughness_mode3", toughness_KIc);
    fracture_energy = parseDouble(config, "fracture_energy", 100.0);
    tensile_strength = parseDouble(config, "tensile_strength", 5e6);
    cohesive_strength = parseDouble(config, "cohesive_strength", tensile_strength);
    critical_opening = parseDouble(config, "critical_opening", 2e-5);
    process_zone_length = parseDouble(config, "process_zone_length", 0.01);
}

double FractureProperties::computeKIc(double E, double nu) const {
    // K_Ic = sqrt(G_c * E')  where E' = E/(1-ν²) for plane strain
    double E_prime = E / (1.0 - nu * nu);
    return std::sqrt(fracture_energy * E_prime);
}

// =============================================================================
// RockProperties Implementation
// =============================================================================

RockProperties::RockProperties()
    : name("default"),
      region_id(0) {
    mechanical = std::make_unique<LinearElasticMaterial>();
    permeability = std::make_unique<PermeabilityModel>();
}

void RockProperties::configure(const std::map<std::string, std::string>& config) {
    name = parseString(config, "name", "default");
    region_id = static_cast<int>(parseDouble(config, "region_id", 0));
    
    // Create appropriate mechanical model
    std::string mech_type = parseString(config, "constitutive_model", "LINEAR_ELASTIC");
    createMechanicalModel(mech_type, config);
    
    // Configure permeability
    permeability->configure(config);
    
    // Configure thermal and fracture properties
    thermal.configure(config);
    fracture.configure(config);
}

void RockProperties::createMechanicalModel(const std::string& type,
                                            const std::map<std::string, std::string>& config) {
    ConstitutiveModel model = parseConstitutiveModel(type);
    
    switch (model) {
        case ConstitutiveModel::LINEAR_ELASTIC:
            mechanical = std::make_unique<LinearElasticMaterial>();
            break;
        case ConstitutiveModel::VISCOELASTIC_MAXWELL:
        case ConstitutiveModel::VISCOELASTIC_KELVIN:
        case ConstitutiveModel::VISCOELASTIC_SLS:
            mechanical = std::make_unique<ViscoelasticMaterial>();
            break;
        case ConstitutiveModel::POROELASTIC:
            mechanical = std::make_unique<PoroelasticMaterial>();
            break;
        case ConstitutiveModel::ELASTOPLASTIC_MC:
        case ConstitutiveModel::ELASTOPLASTIC_DP:
            mechanical = std::make_unique<ElastoplasticMaterial>();
            break;
        case ConstitutiveModel::ANISOTROPIC_VTI:
            mechanical = std::make_unique<AnisotropicMaterial>(true);
            break;
        case ConstitutiveModel::ANISOTROPIC_HTI:
            mechanical = std::make_unique<AnisotropicMaterial>(false);
            break;
        default:
            mechanical = std::make_unique<LinearElasticMaterial>();
    }
    
    mechanical->configure(config);
}

double RockProperties::getDensity() const {
    return mechanical ? mechanical->getDensity() : 2500.0;
}

double RockProperties::getPorosity() const {
    return mechanical ? mechanical->getPorosity() : 0.2;
}

double RockProperties::getYoungsModulus() const {
    return mechanical ? mechanical->getYoungsModulus() : 10e9;
}

double RockProperties::getPoissonsRatio() const {
    return mechanical ? mechanical->getPoissonsRatio() : 0.25;
}

double RockProperties::getPermeability() const {
    return permeability ? permeability->getInitialPermeability() : 100e-15;
}

double RockProperties::getThermalConductivity() const {
    return thermal.thermal_conductivity;
}

// =============================================================================
// Factory Functions
// =============================================================================

std::unique_ptr<MaterialModelBase> createMaterialModel(
    const std::string& type_str,
    const std::map<std::string, std::string>& config) {
    
    ConstitutiveModel model = parseConstitutiveModel(type_str);
    std::unique_ptr<MaterialModelBase> mat;
    
    switch (model) {
        case ConstitutiveModel::LINEAR_ELASTIC:
            mat = std::make_unique<LinearElasticMaterial>();
            break;
        case ConstitutiveModel::VISCOELASTIC_MAXWELL:
            mat = std::make_unique<ViscoelasticMaterial>(ViscoelasticMaterial::ViscoType::MAXWELL);
            break;
        case ConstitutiveModel::VISCOELASTIC_KELVIN:
            mat = std::make_unique<ViscoelasticMaterial>(ViscoelasticMaterial::ViscoType::KELVIN_VOIGT);
            break;
        case ConstitutiveModel::VISCOELASTIC_SLS:
            mat = std::make_unique<ViscoelasticMaterial>(ViscoelasticMaterial::ViscoType::STANDARD_LINEAR_SOLID);
            break;
        case ConstitutiveModel::POROELASTIC:
            mat = std::make_unique<PoroelasticMaterial>();
            break;
        case ConstitutiveModel::ELASTOPLASTIC_MC:
            mat = std::make_unique<ElastoplasticMaterial>(FailureCriterion::MOHR_COULOMB);
            break;
        case ConstitutiveModel::ELASTOPLASTIC_DP:
            mat = std::make_unique<ElastoplasticMaterial>(FailureCriterion::DRUCKER_PRAGER);
            break;
        case ConstitutiveModel::ANISOTROPIC_VTI:
            mat = std::make_unique<AnisotropicMaterial>(true);
            break;
        case ConstitutiveModel::ANISOTROPIC_HTI:
            mat = std::make_unique<AnisotropicMaterial>(false);
            break;
        default:
            mat = std::make_unique<LinearElasticMaterial>();
    }
    
    mat->configure(config);
    return mat;
}

ConstitutiveModel parseConstitutiveModel(const std::string& type_str) {
    std::string s = type_str;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    
    if (s == "LINEAR_ELASTIC" || s == "ELASTIC" || s == "HOOKEAN") {
        return ConstitutiveModel::LINEAR_ELASTIC;
    } else if (s == "VISCOELASTIC" || s == "MAXWELL" || s == "VISCOELASTIC_MAXWELL") {
        return ConstitutiveModel::VISCOELASTIC_MAXWELL;
    } else if (s == "KELVIN" || s == "KELVIN_VOIGT" || s == "VISCOELASTIC_KELVIN") {
        return ConstitutiveModel::VISCOELASTIC_KELVIN;
    } else if (s == "SLS" || s == "STANDARD_LINEAR_SOLID" || s == "VISCOELASTIC_SLS") {
        return ConstitutiveModel::VISCOELASTIC_SLS;
    } else if (s == "POROELASTIC" || s == "BIOT") {
        return ConstitutiveModel::POROELASTIC;
    } else if (s == "ELASTOPLASTIC" || s == "MOHR_COULOMB" || s == "MC") {
        return ConstitutiveModel::ELASTOPLASTIC_MC;
    } else if (s == "DRUCKER_PRAGER" || s == "DP") {
        return ConstitutiveModel::ELASTOPLASTIC_DP;
    } else if (s == "VTI" || s == "ANISOTROPIC_VTI") {
        return ConstitutiveModel::ANISOTROPIC_VTI;
    } else if (s == "HTI" || s == "ANISOTROPIC_HTI") {
        return ConstitutiveModel::ANISOTROPIC_HTI;
    }
    
    return ConstitutiveModel::LINEAR_ELASTIC;
}

} // namespace FSRM
