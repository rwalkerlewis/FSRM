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

}  // namespace FSRM
