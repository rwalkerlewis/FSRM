/**
 * @file ExplosionDamageKernels.cpp
 * @brief Implementations for explosion-related physics kernels (minimal FE coupling)
 */

#include "physics/ExplosionDamageKernels.hpp"

#include <algorithm>
#include <cmath>

namespace FSRM {

ExplosionSourceKernel::ExplosionSourceKernel()
    : PhysicsKernel(PhysicsType::EXPLOSION_SOURCE),
      yield_kt_(0.0),
      depth_(0.0),
      source_x_(0.0),
      source_y_(0.0),
      source_z_(0.0),
      density_(2700.0),
      p_velocity_(5000.0),
      s_velocity_(3000.0),
      scalar_moment_(0.0),
      corner_frequency_(1.0),
      rise_time_(0.05),
      overshoot_(1.0),
      tectonic_release_(false),
      tectonic_orientation_(0.0),
      tectonic_fraction_(0.0),
      source_model_(std::make_unique<MuellerMurphySource>()) {}

PetscErrorCode ExplosionSourceKernel::setup(DM dm, PetscFE fe) {
    (void)dm;
    (void)fe;
    PetscFunctionBeginUser;
    PetscFunctionReturn(0);
}

void ExplosionSourceKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar f[]) {
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)a;
    (void)x;
    for (int i = 0; i < 6; ++i) {
        f[i] = PetscScalar(0.0);
    }
}

void ExplosionSourceKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar J[]) {
    (void)u;
    (void)u_t;
    (void)u_x;
    (void)a;
    (void)x;
    for (int i = 0; i < 36; ++i) {
        J[i] = PetscScalar(0.0);
    }
}

void ExplosionSourceKernel::setExplosionParameters(double yield_kt, double depth_m) {
    yield_kt_ = yield_kt;
    depth_ = depth_m;
    NuclearSourceParameters p;
    p.yield_kt = yield_kt_;
    p.depth_of_burial = depth_;
    p.x = source_x_;
    p.y = source_y_;
    p.z = source_z_;
    source_model_->setParameters(p);
    source_model_->setMediumProperties(density_, p_velocity_, s_velocity_);
    computeDerivedQuantities();
}

void ExplosionSourceKernel::setMediumProperties(double rho, double vp, double vs) {
    density_ = rho;
    p_velocity_ = vp;
    s_velocity_ = vs;
    source_model_->setMediumProperties(density_, p_velocity_, s_velocity_);
    computeDerivedQuantities();
}

void ExplosionSourceKernel::setSourceLocation(double x, double y, double z) {
    source_x_ = x;
    source_y_ = y;
    source_z_ = z;
    NuclearSourceParameters p;
    p.yield_kt = yield_kt_;
    p.depth_of_burial = depth_;
    p.x = source_x_;
    p.y = source_y_;
    p.z = source_z_;
    source_model_->setParameters(p);
}

void ExplosionSourceKernel::setTectonicRelease(bool enable, double orientation, double fraction) {
    tectonic_release_ = enable;
    tectonic_orientation_ = orientation;
    tectonic_fraction_ = std::clamp(fraction, 0.0, 1.0);
}

double ExplosionSourceKernel::getScalarMoment() const {
    return scalar_moment_;
}

double ExplosionSourceKernel::getCornerFrequency() const {
    return corner_frequency_;
}

double ExplosionSourceKernel::getMagnitude() const {
    NuclearSourceParameters p;
    p.yield_kt = yield_kt_;
    return p.body_wave_magnitude();
}

double ExplosionSourceKernel::sourceTimeFunction(double t) const {
    if (t < 0.0) {
        return 0.0;
    }
    const double tau = std::max(1e-12, rise_time_);
    return std::exp(-t / tau);
}

void ExplosionSourceKernel::getMomentTensor(double t, double Mij[6]) const {
    const double s = sourceTimeFunction(t);
    const double m_iso = scalar_moment_ * s / 3.0;
    Mij[0] = Mij[1] = Mij[2] = m_iso;
    Mij[3] = Mij[4] = Mij[5] = 0.0;
    (void)t;
    (void)tectonic_release_;
    (void)tectonic_orientation_;
    (void)tectonic_fraction_;
}

void ExplosionSourceKernel::getPointForce(const PetscReal x[], double t, PetscScalar f[3]) const {
    (void)x;
    (void)t;
    f[0] = f[1] = f[2] = PetscScalar(0.0);
}

void ExplosionSourceKernel::computeDerivedQuantities() {
    NuclearSourceParameters p;
    p.yield_kt = yield_kt_;
    p.depth_of_burial = depth_;
    scalar_moment_ = p.scalar_moment();
    corner_frequency_ = source_model_->getCornerFrequency();
    rise_time_ = 0.55 / std::max(1e-12, corner_frequency_);
    overshoot_ = source_model_->getOvershoot();
}

// =============================================================================
// NearFieldDamageKernel stub implementation
// =============================================================================

NearFieldDamageKernel::NearFieldDamageKernel()
    : PhysicsKernel(PhysicsType::NEAR_FIELD_DAMAGE),
      yield_kt_(0.0), source_x_(0.0), source_y_(0.0), source_z_(0.0),
      tensile_strength_(10e6), compressive_strength_(100e6),
      damage_rate_coeff_(1.0), strain_threshold_(1e-3),
      perm_alpha_(1.0), perm_exponent_(3.0),
      cavity_radius_(0.0), crushed_radius_(0.0), fractured_radius_(0.0) {}

PetscErrorCode NearFieldDamageKernel::setup(DM dm, PetscFE fe) {
    (void)dm; (void)fe; PetscFunctionBeginUser; PetscFunctionReturn(0);
}

void NearFieldDamageKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar f[]) {
    (void)u; (void)u_t; (void)u_x; (void)a; (void)x;
    f[0] = PetscScalar(0.0);
}

void NearFieldDamageKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar J[]) {
    (void)u; (void)u_t; (void)u_x; (void)a; (void)x;
    J[0] = PetscScalar(0.0);
}

void NearFieldDamageKernel::setExplosionSource(double yield_kt, double x, double y, double z) {
    yield_kt_ = yield_kt; source_x_ = x; source_y_ = y; source_z_ = z;
    computeZoneRadii();
}

void NearFieldDamageKernel::setMaterialProperties(double ts, double cs) {
    tensile_strength_ = ts; compressive_strength_ = cs;
}

void NearFieldDamageKernel::setDamageParameters(double rate, double thresh) {
    damage_rate_coeff_ = rate; strain_threshold_ = thresh;
}

void NearFieldDamageKernel::setPermeabilityModel(double alpha, double n) {
    perm_alpha_ = alpha; perm_exponent_ = n;
}

double NearFieldDamageKernel::getCavityRadius() const { return cavity_radius_; }
double NearFieldDamageKernel::getCrushedZoneRadius() const { return crushed_radius_; }
double NearFieldDamageKernel::getFracturedZoneRadius() const { return fractured_radius_; }

double NearFieldDamageKernel::getEnhancedPermeability(double D, double k0) const {
    return k0 * std::pow(1.0 + perm_alpha_ * D, perm_exponent_);
}

NearFieldDamageKernel::DamageZone NearFieldDamageKernel::classifyZone(const PetscReal x[]) const {
    double r = std::sqrt((x[0]-source_x_)*(x[0]-source_x_) +
                         (x[1]-source_y_)*(x[1]-source_y_) +
                         (x[2]-source_z_)*(x[2]-source_z_));
    if (r < cavity_radius_) return DamageZone::CAVITY;
    if (r < crushed_radius_) return DamageZone::CRUSHED;
    if (r < fractured_radius_) return DamageZone::FRACTURED;
    return DamageZone::INTACT;
}

NearFieldDamageKernel::DamageZone NearFieldDamageKernel::classifyZoneFromDamage(double D) const {
    if (D > 0.9) return DamageZone::CAVITY;
    if (D > 0.5) return DamageZone::CRUSHED;
    if (D > 0.01) return DamageZone::FRACTURED;
    return DamageZone::INTACT;
}

void NearFieldDamageKernel::computeZoneRadii() {
    double W = yield_kt_ * 4.184e12; // joules
    cavity_radius_   = 0.05 * std::pow(W, 1.0/3.0);
    crushed_radius_  = 2.0 * cavity_radius_;
    fractured_radius_= 5.0 * cavity_radius_;
}

double NearFieldDamageKernel::damageRate(double D, double strain, double strain_rate) const {
    (void)strain_rate;
    if (strain < strain_threshold_ || D >= 1.0) return 0.0;
    return damage_rate_coeff_ * (1.0 - D) * (strain - strain_threshold_);
}

// =============================================================================
// HydrodynamicKernel stub implementation
// =============================================================================

HydrodynamicKernel::HydrodynamicKernel()
    : PhysicsKernel(PhysicsType::HYDRODYNAMIC),
      eos_type_(EOSType::IDEAL_GAS), gamma_(1.4),
      rho0_mg_(2700.0), c0_mg_(5000.0), s_mg_(1.5), Gamma0_mg_(2.0),
      rho0_t_(2700.0), A_t_(0), B_t_(0), E0_t_(0), a_t_(0), b_t_(0), alpha_t_(0), beta_t_(0),
      A_jwl_(0), B_jwl_(0), R1_jwl_(0), R2_jwl_(0), omega_jwl_(0),
      C_q_(2.0), C_l_(0.5),
      gravity_{0, 0, 0},
      has_source_(false), source_yield_kt_(0), source_x_(0), source_y_(0), source_z_(0), source_t0_(0) {}

PetscErrorCode HydrodynamicKernel::setup(DM dm, PetscFE fe) {
    (void)dm; (void)fe; PetscFunctionBeginUser; PetscFunctionReturn(0);
}

void HydrodynamicKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar f[]) {
    (void)u; (void)u_t; (void)u_x; (void)a; (void)x;
    for (int i = 0; i < 5; ++i) f[i] = PetscScalar(0.0);
}

void HydrodynamicKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                  const PetscScalar u_x[], const PetscScalar a[],
                                  const PetscReal x[], PetscScalar J[]) {
    (void)u; (void)u_t; (void)u_x; (void)a; (void)x;
    for (int i = 0; i < 25; ++i) J[i] = PetscScalar(0.0);
}

void HydrodynamicKernel::setEOS(EOSType type) { eos_type_ = type; }
void HydrodynamicKernel::setIdealGasParameters(double gamma) { gamma_ = gamma; }
void HydrodynamicKernel::setMieGruneisenParameters(double rho0, double c0, double s, double G0) {
    rho0_mg_ = rho0; c0_mg_ = c0; s_mg_ = s; Gamma0_mg_ = G0;
}
void HydrodynamicKernel::setTillotsonParameters(double rho0, double A, double B, double E0,
                                                  double a, double b, double alpha, double beta) {
    rho0_t_ = rho0; A_t_ = A; B_t_ = B; E0_t_ = E0;
    a_t_ = a; b_t_ = b; alpha_t_ = alpha; beta_t_ = beta;
}
void HydrodynamicKernel::setJWLParameters(double A, double B, double R1, double R2, double omega) {
    A_jwl_ = A; B_jwl_ = B; R1_jwl_ = R1; R2_jwl_ = R2; omega_jwl_ = omega;
}
void HydrodynamicKernel::setArtificialViscosity(double Cq, double Cl) { C_q_ = Cq; C_l_ = Cl; }
void HydrodynamicKernel::setExplosionSource(double yield, double x, double y, double z, double t0) {
    has_source_ = true; source_yield_kt_ = yield; source_x_ = x; source_y_ = y; source_z_ = z; source_t0_ = t0;
}
void HydrodynamicKernel::setGravity(double gx, double gy, double gz) {
    gravity_[0] = gx; gravity_[1] = gy; gravity_[2] = gz;
}

double HydrodynamicKernel::pressure(double rho, double e) const {
    switch (eos_type_) {
        case EOSType::IDEAL_GAS: return pressureIdealGas(rho, e);
        case EOSType::MIE_GRUNEISEN: return pressureMieGruneisen(rho, e);
        case EOSType::TILLOTSON: return pressureTillotson(rho, e);
        case EOSType::JWL: return pressureJWL(rho, e);
        default: return pressureIdealGas(rho, e);
    }
}
double HydrodynamicKernel::soundSpeed(double rho, double p) const {
    return std::sqrt(std::max(0.0, gamma_ * p / std::max(1e-30, rho)));
}
double HydrodynamicKernel::temperature(double /*rho*/, double e) const {
    return e / 1000.0; // placeholder
}

double HydrodynamicKernel::pressureIdealGas(double rho, double e) const {
    return (gamma_ - 1.0) * rho * e;
}
double HydrodynamicKernel::pressureMieGruneisen(double rho, double e) const {
    double eta = rho / rho0_mg_ - 1.0;
    double p_H = rho0_mg_ * c0_mg_ * c0_mg_ * eta * (1.0 + 0.5 * (s_mg_ - 1.0) * eta);
    return p_H + Gamma0_mg_ * rho * e;
}
double HydrodynamicKernel::pressureTillotson(double rho, double e) const {
    (void)rho; (void)e; return 0.0; // placeholder
}
double HydrodynamicKernel::pressureJWL(double rho, double e) const {
    (void)rho; (void)e; return 0.0; // placeholder
}

void HydrodynamicKernel::riemannFlux(const double UL[5], const double UR[5],
                                      const double n[3], double flux[5]) const {
    HLLCFlux(UL, UR, n, flux);
}
void HydrodynamicKernel::HLLCFlux(const double UL[5], const double UR[5],
                                    const double n[3], double flux[5]) const {
    (void)UL; (void)UR; (void)n;
    for (int i = 0; i < 5; ++i) flux[i] = 0.0;
}

// =============================================================================
// CraterFormationKernel stub implementation
// =============================================================================

CraterFormationKernel::CraterFormationKernel()
    : PhysicsKernel(PhysicsType::CRATER_FORMATION),
      impactor_diameter_(1.0), impactor_velocity_(1e4), impactor_density_(7800.0),
      impact_angle_(90.0),
      target_density_(2700.0), target_strength_(10e6), gravity_(9.81),
      Z_exponent_(2.7), flow_center_depth_(0.0),
      transient_diameter_(0.0), transient_depth_(0.0),
      final_diameter_(0.0), final_depth_(0.0), formation_time_(0.0) {}

PetscErrorCode CraterFormationKernel::setup(DM dm, PetscFE fe) {
    (void)dm; (void)fe; PetscFunctionBeginUser; PetscFunctionReturn(0);
}

void CraterFormationKernel::residual(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar f[]) {
    (void)u; (void)u_t; (void)u_x; (void)a; (void)x;
    for (int i = 0; i < 6; ++i) f[i] = PetscScalar(0.0);
}

void CraterFormationKernel::jacobian(const PetscScalar u[], const PetscScalar u_t[],
                                     const PetscScalar u_x[], const PetscScalar a[],
                                     const PetscReal x[], PetscScalar J[]) {
    (void)u; (void)u_t; (void)u_x; (void)a; (void)x;
    for (int i = 0; i < 36; ++i) J[i] = PetscScalar(0.0);
}

void CraterFormationKernel::setImpactor(double d, double v, double rho, double angle) {
    impactor_diameter_ = d; impactor_velocity_ = v; impactor_density_ = rho; impact_angle_ = angle;
    computeCraterDimensions();
}
void CraterFormationKernel::setTarget(double rho, double strength, double g) {
    target_density_ = rho; target_strength_ = strength; gravity_ = g;
    computeCraterDimensions();
}
void CraterFormationKernel::setZModelParameters(double Z, double depth) {
    Z_exponent_ = Z; flow_center_depth_ = depth;
}

double CraterFormationKernel::getTransientCraterDiameter() const { return transient_diameter_; }
double CraterFormationKernel::getTransientCraterDepth() const { return transient_depth_; }
double CraterFormationKernel::getFinalCraterDiameter() const { return final_diameter_; }
double CraterFormationKernel::getFinalCraterDepth() const { return final_depth_; }
double CraterFormationKernel::getCraterFormationTime() const { return formation_time_; }
bool CraterFormationKernel::isComplexCrater() const { return final_diameter_ > 3000.0; }

void CraterFormationKernel::getExcavationVelocity(const PetscReal x[], double t,
                                                    PetscScalar v[3]) const {
    (void)x; (void)t; v[0] = v[1] = v[2] = PetscScalar(0.0);
}

double CraterFormationKernel::getEjectaThickness(double r) const {
    if (r <= 0.0 || transient_diameter_ <= 0.0) return 0.0;
    double R = transient_diameter_ / 2.0;
    return 0.14 * R * std::pow(r / R, -3.0);
}

void CraterFormationKernel::getEjectaVelocity(double r, double& v, double& angle) const {
    (void)r; v = 0.0; angle = 45.0;
}

double CraterFormationKernel::getCraterDepth(double r, double t) const {
    (void)r; (void)t; return transient_depth_;
}

void CraterFormationKernel::computeCraterDimensions() {
    double a = impactor_diameter_ / 2.0;
    double KE = 0.5 * impactor_density_ * (4.0/3.0 * M_PI * a*a*a) *
                impactor_velocity_ * impactor_velocity_;
    transient_diameter_ = 1.161 * std::pow(KE / (target_density_ * gravity_), 1.0/3.26);
    transient_depth_ = transient_diameter_ / 3.0;
    final_diameter_ = 1.17 * transient_diameter_;
    final_depth_ = final_diameter_ / 5.0;
    formation_time_ = 0.8 * std::sqrt(transient_diameter_ / gravity_);
}

}  // namespace FSRM
