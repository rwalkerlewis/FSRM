/**
 * @file SeismicSource.cpp
 * @brief Implementation of seismic source models
 */

#include "SeismicSource.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <numeric>

namespace FSRM {

// =============================================================================
// MomentTensor Implementation
// =============================================================================

MomentTensor MomentTensor::fromFaultGeometry(double strike, double dip, double rake, double M0) {
    return doubleCouple(strike, dip, rake, M0);
}

MomentTensor MomentTensor::doubleCouple(double strike, double dip, double rake, double M0) {
    // Convert angles to radians if in degrees
    double phi = strike;     // Strike
    double delta = dip;      // Dip
    double lambda = rake;    // Rake
    
    // Ensure radians
    if (std::abs(phi) > 2.0 * M_PI) phi *= M_PI / 180.0;
    if (std::abs(delta) > M_PI) delta *= M_PI / 180.0;
    if (std::abs(lambda) > M_PI) lambda *= M_PI / 180.0;
    
    // Trigonometric functions
    double sin_d = std::sin(delta);
    double cos_d = std::cos(delta);
    double sin_2d = std::sin(2.0 * delta);
    double cos_2d = std::cos(2.0 * delta);
    double sin_l = std::sin(lambda);
    double cos_l = std::cos(lambda);
    double sin_s = std::sin(phi);
    double cos_s = std::cos(phi);
    double sin_2s = std::sin(2.0 * phi);
    double cos_2s = std::cos(2.0 * phi);
    
    // Aki & Richards convention (for NED coordinate system)
    MomentTensor M;
    
    // Convert to ENU/XYZ if needed (using standard seismological convention)
    M.Mxx = -M0 * (sin_d * cos_l * sin_2s + sin_2d * sin_l * sin_s * sin_s);
    M.Myy =  M0 * (sin_d * cos_l * sin_2s - sin_2d * sin_l * cos_s * cos_s);
    M.Mzz = M0 * sin_2d * sin_l;
    M.Mxy = M0 * (sin_d * cos_l * cos_2s + 0.5 * sin_2d * sin_l * sin_2s);
    M.Mxz = -M0 * (cos_d * cos_l * cos_s + cos_2d * sin_l * sin_s);
    M.Myz = -M0 * (cos_d * cos_l * sin_s - cos_2d * sin_l * cos_s);
    
    return M;
}

MomentTensor MomentTensor::explosion(double M0) {
    MomentTensor M;
    M.Mxx = M0;
    M.Myy = M0;
    M.Mzz = M0;
    M.Mxy = 0;
    M.Mxz = 0;
    M.Myz = 0;
    return M;
}

MomentTensor MomentTensor::clvd(double M0, double strike, double dip) {
    // Compensated Linear Vector Dipole
    // Axis of symmetry along dip direction
    MomentTensor M;
    
    double sin_d = std::sin(dip);
    double cos_d = std::cos(dip);
    double sin_s = std::sin(strike);
    double cos_s = std::cos(strike);
    
    // CLVD has eigenvalues (1, 1, -2) or (-1, -1, 2)
    double n1 = cos_s * sin_d;  // Normal direction
    double n2 = sin_s * sin_d;
    double n3 = cos_d;
    
    M.Mxx = M0 * (1.0 - 3.0 * n1 * n1);
    M.Myy = M0 * (1.0 - 3.0 * n2 * n2);
    M.Mzz = M0 * (1.0 - 3.0 * n3 * n3);
    M.Mxy = -3.0 * M0 * n1 * n2;
    M.Mxz = -3.0 * M0 * n1 * n3;
    M.Myz = -3.0 * M0 * n2 * n3;
    
    return M;
}

double MomentTensor::scalarMoment() const {
    // M0 = sqrt(sum(M_ij^2) / 2)
    double sum = Mxx * Mxx + Myy * Myy + Mzz * Mzz +
                2.0 * (Mxy * Mxy + Mxz * Mxz + Myz * Myz);
    return std::sqrt(sum / 2.0);
}

double MomentTensor::magnitude() const {
    double M0 = scalarMoment();
    if (M0 < 1e6) return -10.0;  // Too small
    return (2.0/3.0) * (std::log10(M0) - 9.1);  // Kanamori (1977)
}

void MomentTensor::decompose(double& M0_iso, double& M0_clvd, double& M0_dc,
                             double& strike, double& dip, double& rake) const {
    // Decompose into isotropic + CLVD + double-couple
    
    // Isotropic part
    double trace = (Mxx + Myy + Mzz) / 3.0;
    M0_iso = std::abs(trace);
    
    // Deviatoric part
    double dev_xx = Mxx - trace;
    double dev_yy = Myy - trace;
    double dev_zz = Mzz - trace;
    
    // Eigenvalue decomposition of deviatoric part
    // For simplicity, assume small CLVD (pure DC)
    M0_dc = scalarMoment();
    M0_clvd = 0.0;
    
    // Extract fault plane (simplified)
    // Full implementation requires eigendecomposition
    strike = 0.0;
    dip = M_PI / 2.0;
    rake = 0.0;
}

void MomentTensor::toArray(double* M) const {
    M[0] = Mxx;
    M[1] = Myy;
    M[2] = Mzz;
    M[3] = Mxy;
    M[4] = Mxz;
    M[5] = Myz;
}

void MomentTensor::toMatrix(double* M) const {
    M[0] = Mxx;  M[1] = Mxy;  M[2] = Mxz;
    M[3] = Mxy;  M[4] = Myy;  M[5] = Myz;
    M[6] = Mxz;  M[7] = Myz;  M[8] = Mzz;
}

MomentTensor MomentTensor::operator+(const MomentTensor& other) const {
    return MomentTensor(Mxx + other.Mxx, Myy + other.Myy, Mzz + other.Mzz,
                       Mxy + other.Mxy, Mxz + other.Mxz, Myz + other.Myz);
}

MomentTensor MomentTensor::operator*(double scale) const {
    return MomentTensor(scale * Mxx, scale * Myy, scale * Mzz,
                       scale * Mxy, scale * Mxz, scale * Myz);
}

// =============================================================================
// SourceTimeFunctionEvaluator Implementation
// =============================================================================

SourceTimeFunctionEvaluator::SourceTimeFunctionEvaluator(SourceTimeFunction t)
    : type(t), duration(1.0), rise_time(0.5), onset_time(0.0),
      peak_time(0.5), frequency(1.0) {}

void SourceTimeFunctionEvaluator::setType(SourceTimeFunction t) {
    type = t;
}

void SourceTimeFunctionEvaluator::setDuration(double T) {
    duration = T;
}

void SourceTimeFunctionEvaluator::setRiseTime(double tau) {
    rise_time = tau;
}

void SourceTimeFunctionEvaluator::setOnsetTime(double t0) {
    onset_time = t0;
}

void SourceTimeFunctionEvaluator::setPeakTime(double tp) {
    peak_time = tp;
}

void SourceTimeFunctionEvaluator::setFrequency(double f0) {
    frequency = f0;
}

void SourceTimeFunctionEvaluator::setCustomFunction(std::function<double(double)> func) {
    custom_func = func;
    type = SourceTimeFunction::CUSTOM;
}

double SourceTimeFunctionEvaluator::momentRate(double t) const {
    double t_rel = t - onset_time;
    if (t_rel < 0) return 0.0;
    
    switch (type) {
        case SourceTimeFunction::GAUSSIAN:
            return gaussianRate(t_rel);
        case SourceTimeFunction::RICKER:
            return rickerRate(t_rel);
        case SourceTimeFunction::STEP:
            return stepRate(t_rel);
        case SourceTimeFunction::RAMP:
            return rampRate(t_rel);
        case SourceTimeFunction::TRIANGLE:
            return triangleRate(t_rel);
        case SourceTimeFunction::SMOOTHED_RAMP:
            return smoothedRampRate(t_rel);
        case SourceTimeFunction::YOFFE:
            return yoffeRate(t_rel);
        case SourceTimeFunction::BRUNE:
            return bruneRate(t_rel);
        case SourceTimeFunction::CUSTOM:
            return custom_func ? custom_func(t_rel) : 0.0;
        default:
            return gaussianRate(t_rel);
    }
}

double SourceTimeFunctionEvaluator::moment(double t) const {
    // Numerical integration of moment rate
    double t_rel = t - onset_time;
    if (t_rel <= 0) return 0.0;
    
    // Simpson's rule
    int n = 100;
    double dt = t_rel / n;
    double sum = momentRate(onset_time) + momentRate(t);
    
    for (int i = 1; i < n; ++i) {
        double ti = onset_time + i * dt;
        if (i % 2 == 0) {
            sum += 2.0 * momentRate(ti);
        } else {
            sum += 4.0 * momentRate(ti);
        }
    }
    
    return sum * dt / 3.0;
}

double SourceTimeFunctionEvaluator::normalizedRate(double t) const {
    double rate = momentRate(t);
    double max_rate = momentRate(onset_time + peak_time);
    if (std::abs(max_rate) < 1e-30) return 0.0;
    return rate / max_rate;
}

double SourceTimeFunctionEvaluator::getDominantFrequency() const {
    return frequency;
}

double SourceTimeFunctionEvaluator::getCornerFrequency() const {
    // Corner frequency fc ≈ 1 / (π * T)
    return 1.0 / (M_PI * duration);
}

double SourceTimeFunctionEvaluator::gaussianRate(double t) const {
    // Gaussian: G(t) = (1/σ√(2π)) * exp(-(t-t_p)²/(2σ²))
    double sigma = duration / 4.0;  // Full width at half max ≈ 2.35σ
    double t_center = peak_time;
    
    double x = (t - t_center) / sigma;
    return std::exp(-0.5 * x * x) / (sigma * std::sqrt(2.0 * M_PI));
}

double SourceTimeFunctionEvaluator::rickerRate(double t) const {
    // Ricker wavelet (Mexican hat)
    // R(t) = (1 - 2π²f₀²(t-t_p)²) * exp(-π²f₀²(t-t_p)²)
    double t_center = peak_time;
    double tau = t - t_center;
    double x = M_PI * frequency * tau;
    double x2 = x * x;
    
    return (1.0 - 2.0 * x2) * std::exp(-x2);
}

double SourceTimeFunctionEvaluator::stepRate(double t) const {
    // Heaviside step function (derivative is delta function)
    // Smoothed version
    double tau = duration / 10.0;
    return 0.5 * (1.0 + std::tanh(t / tau));
}

double SourceTimeFunctionEvaluator::rampRate(double t) const {
    // Linear ramp over rise_time
    if (t <= 0) return 0.0;
    if (t >= rise_time) return 1.0 / rise_time;  // Constant rate
    return t / (rise_time * rise_time);  // Increasing rate
}

double SourceTimeFunctionEvaluator::triangleRate(double t) const {
    // Triangular moment rate
    double T = duration;
    double T_half = T / 2.0;
    
    if (t < 0 || t > T) return 0.0;
    
    if (t <= T_half) {
        return 4.0 * t / (T * T);
    } else {
        return 4.0 * (T - t) / (T * T);
    }
}

double SourceTimeFunctionEvaluator::smoothedRampRate(double t) const {
    // Smoothed (regularized) ramp
    // f(t) = (t/T)² * (3 - 2t/T) for 0 ≤ t ≤ T
    double T = rise_time;
    
    if (t <= 0) return 0.0;
    if (t >= T) return 0.0;  // Rate drops to zero after ramp
    
    double x = t / T;
    double f = x * x * (3.0 - 2.0 * x);  // Moment (integral)
    
    // Derivative: moment rate
    return (6.0 * x * (1.0 - x)) / T;
}

double SourceTimeFunctionEvaluator::yoffeRate(double t) const {
    // Regularized Yoffe function (Tinti et al., 2005)
    // Used for slip rate in kinematic sources
    
    double Tr = rise_time;      // Total rise time
    double Ts = Tr * 0.3;       // Acceleration time (tunable)
    double Tacc = Ts;
    
    if (t < 0 || t > Tr) return 0.0;
    
    // Yoffe function with regularization
    double Cn = 2.0 / (M_PI * Tacc * (1.0 - Tacc / Tr));
    
    if (t < Tacc) {
        return Cn * std::sqrt(t / Tacc);
    } else {
        double x = (t - Tacc) / (Tr - Tacc);
        return Cn * std::sqrt(Tacc / t) * (1.0 - x);
    }
}

double SourceTimeFunctionEvaluator::bruneRate(double t) const {
    // Brune ω⁻² model
    // M'(t) = M0 * ω_c² * t * exp(-ω_c * t)
    
    double fc = getCornerFrequency();
    double omega_c = 2.0 * M_PI * fc;
    
    if (t < 0) return 0.0;
    
    return omega_c * omega_c * t * std::exp(-omega_c * t);
}

// =============================================================================
// PointSource Implementation
// =============================================================================

PointSource::PointSource()
    : name("point_source"), loc_x(0), loc_y(0), loc_z(0),
      scalar_moment(1e18), spatial_sigma(0),
      use_spatial_smoothing(false) {}

void PointSource::setLocation(double x, double y, double z) {
    loc_x = x;
    loc_y = y;
    loc_z = z;
}

void PointSource::setMomentTensor(const MomentTensor& M) {
    moment_tensor = M;
    scalar_moment = M.scalarMoment();
}

void PointSource::setSourceTimeFunction(const SourceTimeFunctionEvaluator& s) {
    stf = s;
}

void PointSource::setScalarMoment(double M0) {
    scalar_moment = M0;
}

void PointSource::setFromFault(double strike, double dip, double rake, double M0,
                               double x, double y, double z) {
    setLocation(x, y, z);
    setMomentTensor(MomentTensor::doubleCouple(strike, dip, rake, M0));
}

void PointSource::getLocation(double& x, double& y, double& z) const {
    x = loc_x;
    y = loc_y;
    z = loc_z;
}

std::array<double, 3> PointSource::getLocation() const {
    return {loc_x, loc_y, loc_z};
}

MomentTensor PointSource::getMomentTensor(double t) const {
    double rate = stf.momentRate(t);
    return moment_tensor * rate;
}

double PointSource::getScalarMoment(double t) const {
    return scalar_moment * stf.moment(t);
}

double PointSource::getMomentRate(double t) const {
    return scalar_moment * stf.momentRate(t);
}

void PointSource::getSourceTerm(double t, double* f) const {
    // Source term for elastic wave equation
    // f_i = -∂j[M_ij * s(t) * δ(x)]
    // In practice, distributed using spatial smoothing
    
    MomentTensor M = getMomentTensor(t);
    
    // For DG/FEM, this contributes to the right-hand side
    // f[0..8] corresponds to σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz, fx, fy, fz
    f[0] = M.Mxx;
    f[1] = M.Myy;
    f[2] = M.Mzz;
    f[3] = M.Mxy;
    f[4] = M.Mxz;
    f[5] = M.Myz;
    f[6] = 0.0;  // Single force x
    f[7] = 0.0;  // Single force y
    f[8] = 0.0;  // Single force z
}

void PointSource::setSpatialSmoothing(double sigma) {
    spatial_sigma = sigma;
    use_spatial_smoothing = (sigma > 0);
}

double PointSource::getSpatialWeight(double x, double y, double z) const {
    if (!use_spatial_smoothing) {
        // Delta function (1 if at source, 0 otherwise)
        double tol = 1e-6;
        if (std::abs(x - loc_x) < tol && 
            std::abs(y - loc_y) < tol &&
            std::abs(z - loc_z) < tol) {
            return 1.0;
        }
        return 0.0;
    }
    
    // Gaussian smoothing
    double dx = x - loc_x;
    double dy = y - loc_y;
    double dz = z - loc_z;
    double r2 = dx*dx + dy*dy + dz*dz;
    double sigma2 = spatial_sigma * spatial_sigma;
    
    return std::exp(-r2 / (2.0 * sigma2)) / 
           std::pow(2.0 * M_PI * sigma2, 1.5);
}

bool PointSource::isActive(double t) const {
    return t >= getOnsetTime() && t <= getEndTime();
}

double PointSource::getOnsetTime() const {
    return stf.momentRate(-1e30) > 0 ? -1e30 : 0.0;
}

double PointSource::getEndTime() const {
    // Estimate end time (when moment rate drops to ~0)
    return 10.0;  // Would depend on STF parameters
}

// =============================================================================
// KinematicSubfault Implementation
// =============================================================================

void KinematicSubfault::computeMomentTensor(double shear_modulus) {
    seismic_moment = shear_modulus * area * slip;
    moment_tensor = MomentTensor::doubleCouple(strike, dip, rake, seismic_moment);
}

double KinematicSubfault::getSlip(double t) const {
    double t_rel = t - rupture_time;
    if (t_rel <= 0) return 0.0;
    if (t_rel >= rise_time) return slip;
    
    // Linear or smoothed slip
    double x = t_rel / rise_time;
    return slip * x * x * (3.0 - 2.0 * x);  // Smooth ramp
}

double KinematicSubfault::getSlipRate(double t) const {
    double t_rel = t - rupture_time;
    if (t_rel <= 0 || t_rel >= rise_time) return 0.0;
    
    // Derivative of slip function
    double x = t_rel / rise_time;
    return slip * 6.0 * x * (1.0 - x) / rise_time;
}

// =============================================================================
// KinematicSource Implementation
// =============================================================================

KinematicSource::KinematicSource()
    : name("kinematic_source"), shear_modulus(30e9),
      hypo_x(0), hypo_y(0), hypo_z(0), spatial_sigma(0) {}

void KinematicSource::createRectangularFault(double center_x, double center_y, double center_z,
                                             double strike, double dip,
                                             double length, double width,
                                             int n_along_strike, int n_down_dip) {
    subfaults.clear();
    
    double dl = length / n_along_strike;
    double dw = width / n_down_dip;
    double area = dl * dw;
    
    // Unit vectors
    double sin_s = std::sin(strike);
    double cos_s = std::cos(strike);
    double sin_d = std::sin(dip);
    double cos_d = std::cos(dip);
    
    // Strike direction
    double sx = sin_s;
    double sy = cos_s;
    double sz = 0.0;
    
    // Down-dip direction
    double dx = -cos_s * cos_d;
    double dy = sin_s * cos_d;
    double dz = -sin_d;
    
    // Corner of fault (top-left)
    double x0 = center_x - 0.5 * length * sx - 0.5 * width * dx;
    double y0 = center_y - 0.5 * length * sy - 0.5 * width * dy;
    double z0 = center_z - 0.5 * length * sz - 0.5 * width * dz;
    
    for (int j = 0; j < n_down_dip; ++j) {
        for (int i = 0; i < n_along_strike; ++i) {
            KinematicSubfault sf;
            
            // Center of subfault
            sf.x = x0 + (i + 0.5) * dl * sx + (j + 0.5) * dw * dx;
            sf.y = y0 + (i + 0.5) * dl * sy + (j + 0.5) * dw * dy;
            sf.z = z0 + (i + 0.5) * dl * sz + (j + 0.5) * dw * dz;
            
            sf.length = dl;
            sf.width = dw;
            sf.area = area;
            sf.strike = strike;
            sf.dip = dip;
            
            sf.slip = 0.0;
            sf.rake = 0.0;
            sf.rupture_time = 0.0;
            sf.rise_time = 1.0;
            
            subfaults.push_back(sf);
        }
    }
    
    // Set hypocenter to center
    hypo_x = center_x;
    hypo_y = center_y;
    hypo_z = center_z;
}

void KinematicSource::loadFromSRF(const std::string& filename) {
    // Standard Rupture Format (SRF)
    // Used by SCEC and other codes
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open SRF file: " + filename);
    }
    
    subfaults.clear();
    
    std::string line;
    
    // Read header
    std::getline(file, line);  // Version
    
    // Read PLANE headers
    int n_planes;
    file >> n_planes;
    
    for (int p = 0; p < n_planes; ++p) {
        // PLANE header
        double elon, elat, nstk, ndip, len, wid;
        double stk, dip, dtop, shyp, dhyp;
        
        file >> elon >> elat >> nstk >> ndip >> len >> wid;
        file >> stk >> dip >> dtop >> shyp >> dhyp;
        
        int n_subfaults = static_cast<int>(nstk * ndip);
        
        // Read POINTS data
        std::string points_header;
        std::getline(file, line);  // Empty line
        std::getline(file, points_header);
        
        int npts;
        sscanf(points_header.c_str(), "POINTS %d", &npts);
        
        for (int i = 0; i < n_subfaults; ++i) {
            KinematicSubfault sf;
            
            double lon, lat, dep, stk_sf, dip_sf, area_sf, tinit, dt;
            double rake_sf, slip1, nt1, slip2, nt2, slip3, nt3;
            
            file >> lon >> lat >> dep >> stk_sf >> dip_sf >> area_sf;
            file >> tinit >> dt >> rake_sf >> slip1 >> nt1 >> slip2 >> nt2 >> slip3 >> nt3;
            
            // Convert to Cartesian (simplified - would use proper projection)
            sf.x = lon * 111000.0 * std::cos(lat * M_PI / 180.0);
            sf.y = lat * 111000.0;
            sf.z = -dep * 1000.0;
            
            sf.strike = stk_sf * M_PI / 180.0;
            sf.dip = dip_sf * M_PI / 180.0;
            sf.area = area_sf * 1e6;  // km² to m²
            
            sf.rupture_time = tinit;
            sf.rake = rake_sf * M_PI / 180.0;
            sf.slip = slip1 / 100.0;  // cm to m
            sf.rise_time = nt1 * dt;
            
            sf.computeMomentTensor(shear_modulus);
            
            subfaults.push_back(sf);
        }
    }
    
    file.close();
}

void KinematicSource::loadFromFSP(const std::string& filename) {
    // Finite-Source Model database format
    // Similar parsing to SRF
    loadFromSRF(filename);  // Placeholder
}

void KinematicSource::setUniformSlip(double slip, double rake) {
    for (auto& sf : subfaults) {
        sf.slip = slip;
        sf.rake = rake;
        sf.computeMomentTensor(shear_modulus);
    }
}

void KinematicSource::setCircularRupture(double hx, double hy, double hz,
                                         double rupture_velocity) {
    hypo_x = hx;
    hypo_y = hy;
    hypo_z = hz;
    
    for (auto& sf : subfaults) {
        double dx = sf.x - hypo_x;
        double dy = sf.y - hypo_y;
        double dz = sf.z - hypo_z;
        double r = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        sf.rupture_time = r / rupture_velocity;
    }
}

void KinematicSource::setSlipFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return;
    
    for (size_t i = 0; i < subfaults.size(); ++i) {
        double slip, rake;
        if (file >> slip >> rake) {
            subfaults[i].slip = slip;
            subfaults[i].rake = rake;
            subfaults[i].computeMomentTensor(shear_modulus);
        }
    }
    
    file.close();
}

void KinematicSource::setUniformRiseTime(double rise_time) {
    for (auto& sf : subfaults) {
        sf.rise_time = rise_time;
    }
}

void KinematicSource::setVariableRiseTime(
    std::function<double(double, double, double)> func) {
    for (auto& sf : subfaults) {
        sf.rise_time = func(sf.x, sf.y, sf.z);
    }
}

void KinematicSource::setShearModulus(double mu) {
    shear_modulus = mu;
    for (auto& sf : subfaults) {
        sf.computeMomentTensor(shear_modulus);
    }
}

double KinematicSource::getTotalMoment() const {
    double M0 = 0.0;
    for (const auto& sf : subfaults) {
        M0 += sf.seismic_moment;
    }
    return M0;
}

double KinematicSource::getMagnitude() const {
    double M0 = getTotalMoment();
    if (M0 < 1e6) return -10.0;
    return (2.0/3.0) * (std::log10(M0) - 9.1);
}

double KinematicSource::getMomentRate(double t) const {
    double rate = 0.0;
    for (const auto& sf : subfaults) {
        double sr = sf.getSlipRate(t);
        rate += shear_modulus * sf.area * sr;
    }
    return rate;
}

void KinematicSource::getSourceTerm(double x, double y, double z, double t,
                                    double* f) const {
    // Initialize
    for (int i = 0; i < 9; ++i) f[i] = 0.0;
    
    // Sum contributions from all subfaults
    for (const auto& sf : subfaults) {
        double slip_rate = sf.getSlipRate(t);
        if (slip_rate < 1e-20) continue;
        
        // Distance weight
        double dx = x - sf.x;
        double dy = y - sf.y;
        double dz = z - sf.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        
        double weight = 1.0;
        if (spatial_sigma > 0) {
            weight = std::exp(-r2 / (2.0 * spatial_sigma * spatial_sigma));
        } else {
            // Nearest subfault
            weight = (r2 < sf.area) ? 1.0 / sf.area : 0.0;
        }
        
        if (weight < 1e-10) continue;
        
        // Moment rate contribution
        double M_rate = shear_modulus * sf.area * slip_rate;
        MomentTensor M = sf.moment_tensor * (M_rate / sf.seismic_moment);
        
        f[0] += weight * M.Mxx;
        f[1] += weight * M.Myy;
        f[2] += weight * M.Mzz;
        f[3] += weight * M.Mxy;
        f[4] += weight * M.Mxz;
        f[5] += weight * M.Myz;
    }
}

double KinematicSource::getStartTime() const {
    double t_min = 1e30;
    for (const auto& sf : subfaults) {
        t_min = std::min(t_min, sf.rupture_time);
    }
    return t_min;
}

double KinematicSource::getEndTime() const {
    double t_max = 0.0;
    for (const auto& sf : subfaults) {
        t_max = std::max(t_max, sf.rupture_time + sf.rise_time);
    }
    return t_max;
}

// =============================================================================
// SingleForceSource Implementation
// =============================================================================

SingleForceSource::SingleForceSource()
    : loc_x(0), loc_y(0), loc_z(0),
      force_x(0), force_y(0), force_z(1.0),
      amplitude(1.0) {}

void SingleForceSource::setLocation(double x, double y, double z) {
    loc_x = x;
    loc_y = y;
    loc_z = z;
}

void SingleForceSource::setForce(double fx, double fy, double fz) {
    force_x = fx;
    force_y = fy;
    force_z = fz;
}

void SingleForceSource::setSourceTimeFunction(const SourceTimeFunctionEvaluator& s) {
    stf = s;
}

void SingleForceSource::setAmplitude(double F0) {
    amplitude = F0;
}

void SingleForceSource::getForce(double t, double& fx, double& fy, double& fz) const {
    double rate = stf.momentRate(t);
    fx = amplitude * force_x * rate;
    fy = amplitude * force_y * rate;
    fz = amplitude * force_z * rate;
}

// =============================================================================
// SeismicReceiver Implementation
// =============================================================================

SeismicReceiver::SeismicReceiver()
    : name("receiver"), loc_x(0), loc_y(0), loc_z(0),
      record_type(RecordType::VELOCITY), sampling_dt(0.01),
      last_record_time(-1e30) {}

void SeismicReceiver::setLocation(double x, double y, double z) {
    loc_x = x;
    loc_y = y;
    loc_z = z;
}

void SeismicReceiver::setName(const std::string& n) {
    name = n;
}

void SeismicReceiver::setRecordType(RecordType type) {
    record_type = type;
}

void SeismicReceiver::setSamplingRate(double dt) {
    sampling_dt = dt;
}

void SeismicReceiver::record(double t, const double* u) {
    if (t - last_record_time < sampling_dt - 1e-10) return;
    
    times.push_back(t);
    data.push_back({u[0], u[1], u[2]});
    last_record_time = t;
}

void SeismicReceiver::recordVelocity(double t, const double* v) {
    record(t, v);
}

void SeismicReceiver::recordStress(double t, const double* sigma) {
    if (t - last_record_time < sampling_dt - 1e-10) return;
    
    times.push_back(t);
    stress_data.push_back({sigma[0], sigma[1], sigma[2], 
                          sigma[3], sigma[4], sigma[5]});
    last_record_time = t;
}

void SeismicReceiver::writeASCII(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "# Receiver: " << name << std::endl;
    file << "# Location: " << loc_x << " " << loc_y << " " << loc_z << std::endl;
    file << "# Time  X  Y  Z" << std::endl;
    
    for (size_t i = 0; i < times.size(); ++i) {
        file << times[i] << " " 
             << data[i][0] << " " 
             << data[i][1] << " " 
             << data[i][2] << std::endl;
    }
    
    file.close();
}

void SeismicReceiver::writeSAC(const std::string& filename) const {
    // SAC format (binary) - simplified
    // Full implementation would write proper SAC headers
    
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) return;
    
    // Write as simple binary for now
    int npts = times.size();
    file.write(reinterpret_cast<const char*>(&npts), sizeof(int));
    file.write(reinterpret_cast<const char*>(&sampling_dt), sizeof(double));
    
    for (size_t i = 0; i < times.size(); ++i) {
        double val = data[i][0];  // X component
        file.write(reinterpret_cast<const char*>(&val), sizeof(double));
    }
    
    file.close();
}

void SeismicReceiver::writeHDF5(const std::string& filename) const {
#ifdef HDF5_FOUND
    // HDF5 output implementation
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        writeASCII(filename + ".txt");  // Fallback
        return;
    }
    
    // Write metadata as attributes
    hid_t attr_space = H5Screate(H5S_SCALAR);
    
    // Write receiver name
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, name.size() + 1);
    hid_t attr_id = H5Acreate(file_id, "name", str_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, str_type, name.c_str());
    H5Aclose(attr_id);
    H5Tclose(str_type);
    
    // Write location
    double location[3] = {loc_x, loc_y, loc_z};
    hsize_t loc_dims[1] = {3};
    hid_t loc_space = H5Screate_simple(1, loc_dims, nullptr);
    hid_t loc_dset = H5Dcreate(file_id, "location", H5T_NATIVE_DOUBLE, loc_space, 
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(loc_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, location);
    H5Dclose(loc_dset);
    H5Sclose(loc_space);
    
    // Write time series
    if (!times.empty()) {
        hsize_t time_dims[1] = {times.size()};
        hid_t time_space = H5Screate_simple(1, time_dims, nullptr);
        hid_t time_dset = H5Dcreate(file_id, "time", H5T_NATIVE_DOUBLE, time_space,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(time_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, times.data());
        H5Dclose(time_dset);
        H5Sclose(time_space);
    }
    
    // Write data (3-component velocity/displacement)
    if (!data.empty()) {
        hsize_t data_dims[2] = {data.size(), 3};
        hid_t data_space = H5Screate_simple(2, data_dims, nullptr);
        
        // Flatten data to contiguous array
        std::vector<double> flat_data;
        flat_data.reserve(data.size() * 3);
        for (const auto& d : data) {
            flat_data.push_back(d[0]);
            flat_data.push_back(d[1]);
            flat_data.push_back(d[2]);
        }
        
        hid_t data_dset = H5Dcreate(file_id, "data", H5T_NATIVE_DOUBLE, data_space,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(data_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat_data.data());
        H5Dclose(data_dset);
        H5Sclose(data_space);
    }
    
    H5Sclose(attr_space);
    H5Fclose(file_id);
#else
    // HDF5 not available - use ASCII format with .h5 extension replaced
    std::string ascii_filename = filename;
    size_t pos = ascii_filename.rfind(".h5");
    if (pos != std::string::npos) {
        ascii_filename = ascii_filename.substr(0, pos) + ".txt";
    } else {
        ascii_filename += ".txt";
    }
    writeASCII(ascii_filename);
#endif
}

double SeismicReceiver::getPeakValue(int component) const {
    double peak = 0.0;
    
    for (const auto& d : data) {
        if (component >= 0 && component < 3) {
            peak = std::max(peak, std::abs(d[component]));
        } else {
            // Vector magnitude
            double mag = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
            peak = std::max(peak, mag);
        }
    }
    
    return peak;
}

double SeismicReceiver::getPeakGroundVelocity() const {
    return getPeakValue(-1);
}

double SeismicReceiver::getPeakGroundAcceleration() const {
    // Compute acceleration from velocity by differentiation
    double peak = 0.0;
    
    for (size_t i = 1; i < data.size(); ++i) {
        double dt = times[i] - times[i-1];
        if (dt < 1e-15) continue;
        
        double ax = (data[i][0] - data[i-1][0]) / dt;
        double ay = (data[i][1] - data[i-1][1]) / dt;
        double az = (data[i][2] - data[i-1][2]) / dt;
        double a_mag = std::sqrt(ax*ax + ay*ay + az*az);
        peak = std::max(peak, a_mag);
    }
    
    return peak;
}

void SeismicReceiver::getLocation(double& x, double& y, double& z) const {
    x = loc_x;
    y = loc_y;
    z = loc_z;
}

// =============================================================================
// FaultReceiver Implementation
// =============================================================================

FaultReceiver::FaultReceiver()
    : name("fault_receiver"), loc_x(0), loc_y(0), loc_z(0),
      sampling_dt(0.001), last_record_time(-1e30) {}

void FaultReceiver::setLocation(double x, double y, double z) {
    loc_x = x;
    loc_y = y;
    loc_z = z;
}

void FaultReceiver::setName(const std::string& n) {
    name = n;
}

void FaultReceiver::setSamplingRate(double dt) {
    sampling_dt = dt;
}

void FaultReceiver::record(double t, double s, double sr,
                           double tn, double ts, double td,
                           double sv) {
    if (t - last_record_time < sampling_dt - 1e-10) return;
    
    times.push_back(t);
    slip.push_back(s);
    slip_rate.push_back(sr);
    traction_n.push_back(tn);
    traction_s.push_back(ts);
    traction_d.push_back(td);
    state_var.push_back(sv);
    last_record_time = t;
}

double FaultReceiver::getRuptureTime(double threshold) const {
    for (size_t i = 0; i < times.size(); ++i) {
        if (slip_rate[i] > threshold) {
            return times[i];
        }
    }
    return -1.0;  // Not ruptured
}

double FaultReceiver::getPeakSlipRate() const {
    if (slip_rate.empty()) return 0.0;
    return *std::max_element(slip_rate.begin(), slip_rate.end());
}

void FaultReceiver::writeASCII(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "# Fault Receiver: " << name << std::endl;
    file << "# Location: " << loc_x << " " << loc_y << " " << loc_z << std::endl;
    file << "# Time  Slip  SlipRate  Tn  Ts  Td  StateVar" << std::endl;
    
    for (size_t i = 0; i < times.size(); ++i) {
        file << times[i] << " " 
             << slip[i] << " " << slip_rate[i] << " "
             << traction_n[i] << " " << traction_s[i] << " " << traction_d[i] << " "
             << state_var[i] << std::endl;
    }
    
    file.close();
}

// =============================================================================
// SourceReceiverManager Implementation
// =============================================================================

SourceReceiverManager::SourceReceiverManager() {}

void SourceReceiverManager::addPointSource(std::unique_ptr<PointSource> source) {
    point_sources.push_back(std::move(source));
}

void SourceReceiverManager::addKinematicSource(std::unique_ptr<KinematicSource> source) {
    kinematic_sources.push_back(std::move(source));
}

void SourceReceiverManager::addSingleForceSource(std::unique_ptr<SingleForceSource> source) {
    force_sources.push_back(std::move(source));
}

void SourceReceiverManager::addReceiver(std::unique_ptr<SeismicReceiver> receiver) {
    receivers.push_back(std::move(receiver));
}

void SourceReceiverManager::addFaultReceiver(std::unique_ptr<FaultReceiver> receiver) {
    fault_receivers.push_back(std::move(receiver));
}

void SourceReceiverManager::addReceiverLine(double x1, double y1, double z1,
                                            double x2, double y2, double z2,
                                            int num_receivers, 
                                            const std::string& prefix) {
    for (int i = 0; i < num_receivers; ++i) {
        double t = static_cast<double>(i) / (num_receivers - 1);
        
        auto rec = std::make_unique<SeismicReceiver>();
        rec->setLocation(x1 + t * (x2 - x1),
                        y1 + t * (y2 - y1),
                        z1 + t * (z2 - z1));
        rec->setName(prefix + std::to_string(i));
        
        receivers.push_back(std::move(rec));
    }
}

void SourceReceiverManager::addReceiverGrid(double x_min, double x_max, 
                                            double y_min, double y_max,
                                            double z, int nx, int ny,
                                            const std::string& prefix) {
    double dx = (x_max - x_min) / (nx - 1);
    double dy = (y_max - y_min) / (ny - 1);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            auto rec = std::make_unique<SeismicReceiver>();
            rec->setLocation(x_min + i * dx, y_min + j * dy, z);
            rec->setName(prefix + std::to_string(j * nx + i));
            
            receivers.push_back(std::move(rec));
        }
    }
}

void SourceReceiverManager::getSourceTerm(double x, double y, double z, double t,
                                          double* f) const {
    // Initialize
    for (int i = 0; i < 9; ++i) f[i] = 0.0;
    
    double f_temp[9];
    
    // Point sources
    for (const auto& src : point_sources) {
        if (!src->isActive(t)) continue;
        
        double w = src->getSpatialWeight(x, y, z);
        if (w < 1e-15) continue;
        
        src->getSourceTerm(t, f_temp);
        for (int i = 0; i < 9; ++i) {
            f[i] += w * f_temp[i];
        }
    }
    
    // Kinematic sources
    for (const auto& src : kinematic_sources) {
        src->getSourceTerm(x, y, z, t, f_temp);
        for (int i = 0; i < 9; ++i) {
            f[i] += f_temp[i];
        }
    }
    
    // Single force sources
    for (const auto& src : force_sources) {
        double fx, fy, fz;
        src->getForce(t, fx, fy, fz);
        
        // Would need spatial weight
        f[6] += fx;
        f[7] += fy;
        f[8] += fz;
    }
}

void SourceReceiverManager::recordAll(double t,
    std::function<void(double, double, double, double*)> get_solution) {
    
    for (auto& rec : receivers) {
        double x, y, z;
        rec->getLocation(x, y, z);
        
        double u[3];
        get_solution(x, y, z, u);
        rec->record(t, u);
    }
}

void SourceReceiverManager::writeAllReceivers(const std::string& output_dir) const {
    for (const auto& rec : receivers) {
        rec->writeASCII(output_dir + "/" + rec->name + ".txt");
    }
    
    for (const auto& rec : fault_receivers) {
        rec->writeASCII(output_dir + "/" + rec->name + "_fault.txt");
    }
}

double SourceReceiverManager::getSourceStartTime() const {
    double t_min = 1e30;
    
    for (const auto& src : point_sources) {
        t_min = std::min(t_min, src->getOnsetTime());
    }
    
    for (const auto& src : kinematic_sources) {
        t_min = std::min(t_min, src->getStartTime());
    }
    
    return t_min;
}

double SourceReceiverManager::getSourceEndTime() const {
    double t_max = 0.0;
    
    for (const auto& src : point_sources) {
        t_max = std::max(t_max, src->getEndTime());
    }
    
    for (const auto& src : kinematic_sources) {
        t_max = std::max(t_max, src->getEndTime());
    }
    
    return t_max;
}

// =============================================================================
// Configuration Helpers
// =============================================================================

void SourceConfig::parseConfig(const std::map<std::string, std::string>& config) {
    auto get_double = [&](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end()) {
            try { return std::stod(it->second); }
            catch (...) { return def; }
        }
        return def;
    };
    
    auto get_string = [&](const std::string& key, const std::string& def) {
        auto it = config.find(key);
        return (it != config.end()) ? it->second : def;
    };
    
    auto get_bool = [&](const std::string& key, bool def) {
        auto it = config.find(key);
        if (it != config.end()) {
            return it->second == "true" || it->second == "1";
        }
        return def;
    };
    
    use_point_source = get_bool("use_point_source", false);
    source_time_function = get_string("source_time_function", "gaussian");
    source_duration = get_double("source_duration", 1.0);
    source_onset = get_double("source_onset_time", 0.0);
    source_x = get_double("source_x", 0.0);
    source_y = get_double("source_y", 0.0);
    source_z = get_double("source_z", 0.0);
    strike = get_double("source_strike", 0.0);
    dip = get_double("source_dip", 90.0);
    rake = get_double("source_rake", 0.0);
    scalar_moment = get_double("scalar_moment", 1e18);
    
    use_kinematic_source = get_bool("use_kinematic_source", false);
    kinematic_file = get_string("kinematic_source_file", "");
    rupture_velocity = get_double("rupture_velocity", 3000.0);
    rise_time = get_double("rise_time", 1.0);
}

std::unique_ptr<PointSource> SourceConfig::createPointSource() const {
    auto src = std::make_unique<PointSource>();
    
    src->setLocation(source_x, source_y, source_z);
    
    SourceTimeFunctionEvaluator stf;
    if (source_time_function == "gaussian") {
        stf.setType(SourceTimeFunction::GAUSSIAN);
    } else if (source_time_function == "ricker") {
        stf.setType(SourceTimeFunction::RICKER);
    } else if (source_time_function == "triangle") {
        stf.setType(SourceTimeFunction::TRIANGLE);
    } else if (source_time_function == "yoffe") {
        stf.setType(SourceTimeFunction::YOFFE);
    }
    stf.setDuration(source_duration);
    stf.setOnsetTime(source_onset);
    
    src->setSourceTimeFunction(stf);
    src->setFromFault(strike * M_PI / 180.0, dip * M_PI / 180.0, 
                     rake * M_PI / 180.0, scalar_moment,
                     source_x, source_y, source_z);
    
    return src;
}

std::unique_ptr<KinematicSource> SourceConfig::createKinematicSource() const {
    auto src = std::make_unique<KinematicSource>();
    
    if (!kinematic_file.empty()) {
        src->loadFromSRF(kinematic_file);
    }
    
    src->setUniformRiseTime(rise_time);
    
    return src;
}

void ReceiverConfig::parseConfig(const std::map<std::string, std::string>& config) {
    auto get_string = [&](const std::string& key, const std::string& def) {
        auto it = config.find(key);
        return (it != config.end()) ? it->second : def;
    };
    
    auto get_double = [&](const std::string& key, double def) {
        auto it = config.find(key);
        if (it != config.end()) {
            try { return std::stod(it->second); }
            catch (...) { return def; }
        }
        return def;
    };
    
    output_format = get_string("receiver_output_format", "ASCII");
    sampling_rate = 1.0 / get_double("receiver_sampling_dt", 0.01);
    receiver_file = get_string("receiver_file", "");
}

std::vector<std::unique_ptr<SeismicReceiver>> ReceiverConfig::createReceivers() const {
    std::vector<std::unique_ptr<SeismicReceiver>> recs;
    
    // Load from file if specified
    if (!receiver_file.empty()) {
        std::ifstream file(receiver_file);
        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue;
                
                std::istringstream iss(line);
                std::string name;
                double x, y, z;
                
                if (iss >> name >> x >> y >> z) {
                    auto rec = std::make_unique<SeismicReceiver>();
                    rec->setName(name);
                    rec->setLocation(x, y, z);
                    rec->setSamplingRate(1.0 / sampling_rate);
                    recs.push_back(std::move(rec));
                }
            }
            file.close();
        }
    }
    
    // Add inline receivers
    for (size_t i = 0; i < receiver_locations.size(); ++i) {
        auto rec = std::make_unique<SeismicReceiver>();
        rec->setLocation(receiver_locations[i][0],
                        receiver_locations[i][1],
                        receiver_locations[i][2]);
        if (i < receiver_names.size()) {
            rec->setName(receiver_names[i]);
        } else {
            rec->setName("receiver_" + std::to_string(i));
        }
        rec->setSamplingRate(1.0 / sampling_rate);
        recs.push_back(std::move(rec));
    }
    
    return recs;
}

} // namespace FSRM
