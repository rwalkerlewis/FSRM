/**
 * @file FractureNetwork.cpp
 * @brief Implementation of advanced discrete fracture network (DFN) algorithms
 * 
 * This file implements comprehensive fracture network modeling capabilities
 * including stochastic generation, connectivity analysis, and upscaling.
 */

#include "FractureNetwork.hpp"
#include <algorithm>
#include <numeric>
#include <queue>
#include <stack>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace FSRM {

// ============================================================================
// DiscreteFracture Implementation
// ============================================================================

DiscreteFracture::DiscreteFracture()
    : id(-1), set_id(-1),
      center({0.0, 0.0, 0.0}),
      normal({0.0, 0.0, 1.0}),
      strike_dir({1.0, 0.0, 0.0}),
      dip_dir({0.0, 1.0, 0.0}),
      radius(1.0), length(2.0), height(2.0),
      strike(0.0), dip(90.0),
      aperture(1e-4), permeability(0.0), transmissivity(0.0), storativity(1e-9),
      normal_stiffness(1e10), shear_stiffness(1e9),
      friction_angle(30.0), cohesion(0.0), dilation_angle(0.0),
      is_conductive(true), current_aperture(1e-4), cumulative_slip(0.0) {
    computePermeability();
}

void DiscreteFracture::computeVectors() {
    // Convert strike and dip angles (degrees) to direction vectors
    // Strike: clockwise from north (y-axis), Dip: angle from horizontal
    
    double strike_rad = strike * M_PI / 180.0;
    double dip_rad = dip * M_PI / 180.0;
    
    // Strike direction (horizontal, perpendicular to dip direction)
    strike_dir[0] = std::sin(strike_rad);
    strike_dir[1] = std::cos(strike_rad);
    strike_dir[2] = 0.0;
    
    // Dip direction (right-hand rule from strike)
    double dip_azimuth = strike_rad + M_PI / 2.0;
    dip_dir[0] = std::sin(dip_azimuth) * std::cos(dip_rad);
    dip_dir[1] = std::cos(dip_azimuth) * std::cos(dip_rad);
    dip_dir[2] = -std::sin(dip_rad);
    
    // Normal vector (cross product of strike and dip directions)
    normal[0] = strike_dir[1] * dip_dir[2] - strike_dir[2] * dip_dir[1];
    normal[1] = strike_dir[2] * dip_dir[0] - strike_dir[0] * dip_dir[2];
    normal[2] = strike_dir[0] * dip_dir[1] - strike_dir[1] * dip_dir[0];
    
    // Normalize
    double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    if (norm > 1e-10) {
        normal[0] /= norm;
        normal[1] /= norm;
        normal[2] /= norm;
    }
}

void DiscreteFracture::computePermeability() {
    // Cubic law: k = b²/12
    permeability = aperture * aperture / 12.0;
    transmissivity = permeability * aperture;  // T = k * b
}

double DiscreteFracture::computeArea() const {
    // For a disc: A = π * r²
    // For a rectangle: A = length * height
    if (radius > 0) {
        return M_PI * radius * radius;
    }
    return length * height;
}

double DiscreteFracture::computeAperture(double normal_stress, 
                                          double reference_aperture,
                                          double reference_stress) const {
    // Barton-Bandis hyperbolic closure model
    // b = b0 * (1 - σn / (σn + Kn * b0))
    // Simplified exponential model: b = b0 * exp(-σn / σref)
    
    if (normal_stress <= 0.0) {
        return reference_aperture;
    }
    
    double effective_stress = normal_stress;
    if (reference_stress > 0.0) {
        return reference_aperture * std::exp(-effective_stress / reference_stress);
    }
    
    // Using normal stiffness
    double closure = effective_stress / normal_stiffness;
    return std::max(reference_aperture - closure, reference_aperture * 0.01);
}

bool DiscreteFracture::containsPoint(const std::array<double, 3>& point, 
                                      double tolerance) const {
    // Check if point lies on fracture plane within tolerance
    
    // Vector from center to point
    double dx = point[0] - center[0];
    double dy = point[1] - center[1];
    double dz = point[2] - center[2];
    
    // Distance from plane
    double dist_to_plane = std::abs(dx * normal[0] + dy * normal[1] + dz * normal[2]);
    
    if (dist_to_plane > tolerance) {
        return false;
    }
    
    // Check if within fracture extent (using radius for disc)
    double dist_sq = dx*dx + dy*dy + dz*dz - dist_to_plane*dist_to_plane;
    
    if (radius > 0) {
        return dist_sq <= radius * radius;
    }
    
    // For rectangle, project onto strike and dip directions
    double strike_proj = dx * strike_dir[0] + dy * strike_dir[1] + dz * strike_dir[2];
    double dip_proj = dx * dip_dir[0] + dy * dip_dir[1] + dz * dip_dir[2];
    
    return std::abs(strike_proj) <= length / 2.0 && std::abs(dip_proj) <= height / 2.0;
}

std::vector<std::array<double, 3>> DiscreteFracture::getVertices(int n_sides) const {
    std::vector<std::array<double, 3>> vertices;
    
    if (radius > 0) {
        // Disc approximation with n_sides
        for (int i = 0; i < n_sides; ++i) {
            double angle = 2.0 * M_PI * i / n_sides;
            double local_x = radius * std::cos(angle);
            double local_y = radius * std::sin(angle);
            
            std::array<double, 3> vertex;
            vertex[0] = center[0] + local_x * strike_dir[0] + local_y * dip_dir[0];
            vertex[1] = center[1] + local_x * strike_dir[1] + local_y * dip_dir[1];
            vertex[2] = center[2] + local_x * strike_dir[2] + local_y * dip_dir[2];
            vertices.push_back(vertex);
        }
    } else {
        // Rectangle: 4 corners
        double half_l = length / 2.0;
        double half_h = height / 2.0;
        
        std::array<std::pair<double, double>, 4> corners = {{
            {-half_l, -half_h}, {half_l, -half_h},
            {half_l, half_h}, {-half_l, half_h}
        }};
        
        for (const auto& corner : corners) {
            std::array<double, 3> vertex;
            vertex[0] = center[0] + corner.first * strike_dir[0] + corner.second * dip_dir[0];
            vertex[1] = center[1] + corner.first * strike_dir[1] + corner.second * dip_dir[1];
            vertex[2] = center[2] + corner.first * strike_dir[2] + corner.second * dip_dir[2];
            vertices.push_back(vertex);
        }
    }
    
    return vertices;
}

// ============================================================================
// FractureIntensity Implementation
// ============================================================================

FractureIntensity::FractureIntensity()
    : P10(0.0), P20(0.0), P21(0.0), P30(0.0), P32(0.0), P33(0.0) {}

void FractureIntensity::computeFromP32(double p32, double mean_radius) {
    // P32 = fracture area per unit volume
    // Assuming circular fractures with mean radius r:
    // P32 = N * π * r² / V, so P30 = P32 / (π * r²)
    
    P32 = p32;
    double mean_area = M_PI * mean_radius * mean_radius;
    P30 = p32 / mean_area;
    
    // P33 = P32 * mean_aperture (volume fraction)
    // Approximate other measures
    P21 = p32 * 2.0 * mean_radius;  // Trace length ~ 2*r per fracture
    P10 = P30 * mean_radius;         // Linear density approximation
    P20 = P30 * mean_radius;         // Areal density approximation
}

void FractureIntensity::computeFromP21(double p21, double mean_trace_length) {
    // P21 = trace length per unit area (from window sampling)
    
    P21 = p21;
    
    // Approximate P32 (Mauldon's estimator)
    // P32 ≈ (π/2) * P21 for randomly oriented fractures
    P32 = (M_PI / 2.0) * p21;
    
    // Approximate other measures
    double mean_area = (M_PI / 4.0) * mean_trace_length * mean_trace_length;
    P30 = P32 / mean_area;
    P10 = p21 / mean_trace_length;
    P20 = P30 * mean_trace_length;
}

// ============================================================================
// BinghamParameters Implementation
// ============================================================================

BinghamParameters::BinghamParameters() : kappa1(0.0), kappa2(0.0) {
    // Default: principal axes aligned with coordinate system
    principal_axes[0] = {1.0, 0.0, 0.0};
    principal_axes[1] = {0.0, 1.0, 0.0};
    principal_axes[2] = {0.0, 0.0, 1.0};
}

// ============================================================================
// FractureIntersection Implementation
// ============================================================================

FractureIntersection::FractureIntersection()
    : fracture1_id(-1), fracture2_id(-1),
      point1({0.0, 0.0, 0.0}), point2({0.0, 0.0, 0.0}),
      length(0.0), transmissivity(0.0) {}

// ============================================================================
// FractureSet Implementation
// ============================================================================

FractureSet::FractureSet(int id, const std::string& name)
    : id_(id), name_(name),
      orientation_type_(OrientationDistribution::FISHER),
      fisher_params_(0.0, 90.0, 10.0),
      size_type_(SizeDistribution::LOGNORMAL),
      size_param1_(1.0), size_param2_(0.5), size_param3_(0.0),
      min_size_(0.1), max_size_(100.0),
      spatial_type_(SpatialDistribution::UNIFORM_POISSON),
      p32_(0.1), cluster_radius_(10.0), cluster_intensity_(5.0),
      aperture_type_(SizeDistribution::LOGNORMAL),
      aperture_param1_(1e-4), aperture_param2_(0.3),
      aperture_correlation_exp_(0.5),
      normal_stiffness_(1e10), shear_stiffness_(1e9),
      friction_angle_(30.0), cohesion_(0.0),
      termination_type_(TerminationType::NONE),
      termination_probability_(1.0) {}

void FractureSet::setOrientationDistribution(OrientationDistribution type) {
    orientation_type_ = type;
}

void FractureSet::setFisherParameters(const FisherParameters& params) {
    fisher_params_ = params;
}

void FractureSet::setBinghamParameters(const BinghamParameters& params) {
    bingham_params_ = params;
}

void FractureSet::setBootstrapData(const std::vector<std::pair<double, double>>& strike_dip_pairs) {
    bootstrap_data_ = strike_dip_pairs;
}

void FractureSet::setSizeDistribution(SizeDistribution type) {
    size_type_ = type;
}

void FractureSet::setSizeParameters(double param1, double param2, double param3) {
    size_param1_ = param1;
    size_param2_ = param2;
    size_param3_ = param3;
}

void FractureSet::setMinMaxSize(double min_size, double max_size) {
    min_size_ = min_size;
    max_size_ = max_size;
}

void FractureSet::setSpatialDistribution(SpatialDistribution type) {
    spatial_type_ = type;
}

void FractureSet::setIntensity(double p32) {
    p32_ = p32;
}

void FractureSet::setClusteringParameters(double cluster_radius, double cluster_intensity) {
    cluster_radius_ = cluster_radius;
    cluster_intensity_ = cluster_intensity;
}

void FractureSet::setApertureDistribution(SizeDistribution type) {
    aperture_type_ = type;
}

void FractureSet::setApertureParameters(double param1, double param2) {
    aperture_param1_ = param1;
    aperture_param2_ = param2;
}

void FractureSet::setApertureCorrelation(double exponent) {
    aperture_correlation_exp_ = exponent;
}

void FractureSet::setMechanicalProperties(double Kn, double Ks, double phi, double c) {
    normal_stiffness_ = Kn;
    shear_stiffness_ = Ks;
    friction_angle_ = phi;
    cohesion_ = c;
}

void FractureSet::setTerminationRule(TerminationType type, double probability) {
    termination_type_ = type;
    termination_probability_ = probability;
}

void FractureSet::addTerminatingSet(int set_id) {
    terminating_sets_.push_back(set_id);
}

std::pair<double, double> FractureSet::sampleFisherOrientation(std::mt19937& rng) const {
    // Fisher (von Mises-Fisher) distribution on sphere
    // Samples orientations clustered around mean pole with concentration kappa
    
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    double kappa = fisher_params_.kappa;
    
    // Sample theta (deviation from mean) using inverse CDF
    double u = uniform(rng);
    double cos_theta;
    
    if (kappa < 1e-6) {
        // Uniform distribution
        cos_theta = 2.0 * u - 1.0;
    } else {
        // Fisher distribution: cos(theta) = 1 + ln(u + (1-u)*exp(-2*kappa)) / kappa
        cos_theta = 1.0 + std::log(u + (1.0 - u) * std::exp(-2.0 * kappa)) / kappa;
    }
    
    double theta = std::acos(std::max(-1.0, std::min(1.0, cos_theta)));
    
    // Sample phi uniformly
    double phi = 2.0 * M_PI * uniform(rng);
    
    // Convert to strike/dip relative to mean
    // Apply rotation (simplified - rotate around mean pole)
    double strike = fisher_params_.mean_strike + (phi * 180.0 / M_PI);
    double dip = fisher_params_.mean_dip + (theta * 180.0 / M_PI) * std::cos(phi);
    
    // Normalize angles
    while (strike < 0.0) strike += 360.0;
    while (strike >= 360.0) strike -= 360.0;
    dip = std::max(0.0, std::min(90.0, dip));
    
    return {strike, dip};
}

std::pair<double, double> FractureSet::sampleBinghamOrientation(std::mt19937& rng) const {
    // Bingham distribution - more complex, approximated here
    // Supports both clustered and girdle distributions
    
    std::normal_distribution<double> normal(0.0, 1.0);
    
    // Generate random point on sphere weighted by Bingham distribution
    // Simplified: use principal axes with concentration parameters
    
    double x1 = normal(rng) / std::sqrt(1.0 + std::abs(bingham_params_.kappa1));
    double x2 = normal(rng) / std::sqrt(1.0 + std::abs(bingham_params_.kappa2));
    double x3 = normal(rng);
    
    // Normalize to unit sphere
    double norm = std::sqrt(x1*x1 + x2*x2 + x3*x3);
    x1 /= norm; x2 /= norm; x3 /= norm;
    
    // Transform to global coordinates using principal axes
    double gx = x1 * bingham_params_.principal_axes[0][0] + 
                x2 * bingham_params_.principal_axes[1][0] + 
                x3 * bingham_params_.principal_axes[2][0];
    double gy = x1 * bingham_params_.principal_axes[0][1] + 
                x2 * bingham_params_.principal_axes[1][1] + 
                x3 * bingham_params_.principal_axes[2][1];
    double gz = x1 * bingham_params_.principal_axes[0][2] + 
                x2 * bingham_params_.principal_axes[1][2] + 
                x3 * bingham_params_.principal_axes[2][2];
    
    // Convert to strike/dip
    double dip = std::acos(std::abs(gz)) * 180.0 / M_PI;
    double strike = std::atan2(gx, gy) * 180.0 / M_PI;
    
    while (strike < 0.0) strike += 360.0;
    
    return {strike, dip};
}

std::pair<double, double> FractureSet::sampleOrientation(std::mt19937& rng) const {
    switch (orientation_type_) {
        case OrientationDistribution::UNIFORM: {
            std::uniform_real_distribution<double> strike_dist(0.0, 360.0);
            std::uniform_real_distribution<double> cos_dip_dist(-1.0, 1.0);
            double strike = strike_dist(rng);
            double dip = std::acos(cos_dip_dist(rng)) * 180.0 / M_PI;
            return {strike, dip};
        }
        
        case OrientationDistribution::FISHER:
            return sampleFisherOrientation(rng);
        
        case OrientationDistribution::BINGHAM:
            return sampleBinghamOrientation(rng);
        
        case OrientationDistribution::KENT:
            // Kent distribution (Fisher-Bingham) - use Fisher as approximation
            return sampleFisherOrientation(rng);
        
        case OrientationDistribution::BOOTSTRAP: {
            if (bootstrap_data_.empty()) {
                return sampleFisherOrientation(rng);
            }
            std::uniform_int_distribution<size_t> idx_dist(0, bootstrap_data_.size() - 1);
            return bootstrap_data_[idx_dist(rng)];
        }
        
        default:
            return sampleFisherOrientation(rng);
    }
}

double FractureSet::samplePowerLaw(double x_min, double x_max, double alpha, 
                                    std::mt19937& rng) const {
    // Power-law distribution: P(X > x) ~ x^(-alpha)
    // Inverse CDF: x = x_min * (1 - u * (1 - (x_min/x_max)^(alpha-1)))^(-1/(alpha-1))
    
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(rng);
    
    if (std::abs(alpha - 1.0) < 1e-6) {
        // Log-uniform for alpha = 1
        return x_min * std::pow(x_max / x_min, u);
    }
    
    double term = 1.0 - std::pow(x_min / x_max, alpha - 1.0);
    return x_min * std::pow(1.0 - u * term, -1.0 / (alpha - 1.0));
}

double FractureSet::sampleSize(std::mt19937& rng) const {
    double size = 0.0;
    
    switch (size_type_) {
        case SizeDistribution::CONSTANT:
            size = size_param1_;
            break;
        
        case SizeDistribution::UNIFORM: {
            std::uniform_real_distribution<double> dist(size_param1_, size_param2_);
            size = dist(rng);
            break;
        }
        
        case SizeDistribution::NORMAL: {
            std::normal_distribution<double> dist(size_param1_, size_param2_);
            size = dist(rng);
            break;
        }
        
        case SizeDistribution::LOGNORMAL: {
            // param1 = mean of log, param2 = std of log
            std::lognormal_distribution<double> dist(std::log(size_param1_), size_param2_);
            size = dist(rng);
            break;
        }
        
        case SizeDistribution::POWER_LAW:
            // param1 = x_min, param2 = x_max, param3 = alpha (exponent)
            size = samplePowerLaw(size_param1_, size_param2_, size_param3_, rng);
            break;
        
        case SizeDistribution::EXPONENTIAL: {
            std::exponential_distribution<double> dist(1.0 / size_param1_);
            size = dist(rng);
            break;
        }
        
        case SizeDistribution::PARETO: {
            // Pareto: x = x_min * u^(-1/alpha)
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            double u = uniform(rng);
            size = size_param1_ * std::pow(u, -1.0 / size_param2_);
            break;
        }
        
        default:
            size = size_param1_;
    }
    
    // Enforce bounds
    return std::max(min_size_, std::min(max_size_, size));
}

double FractureSet::sampleAperture(double fracture_size, std::mt19937& rng) const {
    double base_aperture = 0.0;
    
    switch (aperture_type_) {
        case SizeDistribution::CONSTANT:
            base_aperture = aperture_param1_;
            break;
        
        case SizeDistribution::LOGNORMAL: {
            std::lognormal_distribution<double> dist(std::log(aperture_param1_), aperture_param2_);
            base_aperture = dist(rng);
            break;
        }
        
        case SizeDistribution::NORMAL: {
            std::normal_distribution<double> dist(aperture_param1_, aperture_param2_);
            base_aperture = std::max(1e-6, dist(rng));
            break;
        }
        
        default:
            base_aperture = aperture_param1_;
    }
    
    // Apply size-aperture correlation: b ~ L^exponent
    if (aperture_correlation_exp_ > 0.0 && fracture_size > 0.0) {
        double size_factor = std::pow(fracture_size / size_param1_, aperture_correlation_exp_);
        base_aperture *= size_factor;
    }
    
    return base_aperture;
}

int FractureSet::computeFractureCount(double volume) const {
    // Compute number of fractures from P32 intensity
    // P32 = total_area / volume
    // N = P32 * volume / mean_area
    
    double mean_radius = size_param1_;
    double mean_area = M_PI * mean_radius * mean_radius;
    
    int expected_count = static_cast<int>(p32_ * volume / mean_area);
    return std::max(1, expected_count);
}

std::vector<DiscreteFracture> FractureSet::generateFractures(
    const std::array<double, 3>& domain_min,
    const std::array<double, 3>& domain_max,
    std::mt19937& rng) const {
    
    std::vector<DiscreteFracture> fractures;
    
    // Compute domain volume
    double volume = (domain_max[0] - domain_min[0]) *
                    (domain_max[1] - domain_min[1]) *
                    (domain_max[2] - domain_min[2]);
    
    // Determine number of fractures based on intensity
    int n_fractures = computeFractureCount(volume);
    
    // Poisson variation
    std::poisson_distribution<int> poisson(n_fractures);
    n_fractures = poisson(rng);
    
    std::uniform_real_distribution<double> x_dist(domain_min[0], domain_max[0]);
    std::uniform_real_distribution<double> y_dist(domain_min[1], domain_max[1]);
    std::uniform_real_distribution<double> z_dist(domain_min[2], domain_max[2]);
    
    // Generate parent locations for clustered distributions
    std::vector<std::array<double, 3>> cluster_centers;
    if (spatial_type_ == SpatialDistribution::CLUSTERED_POISSON ||
        spatial_type_ == SpatialDistribution::PARENT_DAUGHTER) {
        
        int n_clusters = std::max(1, static_cast<int>(n_fractures / cluster_intensity_));
        for (int i = 0; i < n_clusters; ++i) {
            cluster_centers.push_back({x_dist(rng), y_dist(rng), z_dist(rng)});
        }
    }
    
    for (int i = 0; i < n_fractures; ++i) {
        DiscreteFracture frac;
        frac.id = i;
        frac.set_id = id_;
        
        // Sample orientation
        auto [strike, dip] = sampleOrientation(rng);
        frac.strike = strike;
        frac.dip = dip;
        frac.computeVectors();
        
        // Sample size
        double size = sampleSize(rng);
        frac.radius = size / 2.0;
        frac.length = size;
        frac.height = size;
        
        // Sample aperture
        frac.aperture = sampleAperture(size, rng);
        frac.current_aperture = frac.aperture;
        frac.computePermeability();
        
        // Sample location based on spatial distribution
        switch (spatial_type_) {
            case SpatialDistribution::UNIFORM_POISSON:
                frac.center = {x_dist(rng), y_dist(rng), z_dist(rng)};
                break;
            
            case SpatialDistribution::CLUSTERED_POISSON:
            case SpatialDistribution::PARENT_DAUGHTER: {
                // Pick a random cluster center
                std::uniform_int_distribution<size_t> cluster_dist(0, cluster_centers.size() - 1);
                auto& center = cluster_centers[cluster_dist(rng)];
                
                // Add offset within cluster radius
                std::normal_distribution<double> offset_dist(0.0, cluster_radius_ / 3.0);
                frac.center[0] = center[0] + offset_dist(rng);
                frac.center[1] = center[1] + offset_dist(rng);
                frac.center[2] = center[2] + offset_dist(rng);
                break;
            }
            
            case SpatialDistribution::FRACTAL: {
                // Box-fractal distribution (simplified)
                // Higher probability near existing fractures
                if (fractures.empty()) {
                    frac.center = {x_dist(rng), y_dist(rng), z_dist(rng)};
                } else {
                    // 50% chance to cluster near existing fracture
                    std::uniform_real_distribution<double> p(0.0, 1.0);
                    if (p(rng) < 0.5) {
                        std::uniform_int_distribution<size_t> idx(0, fractures.size() - 1);
                        auto& ref = fractures[idx(rng)];
                        std::normal_distribution<double> offset(0.0, ref.radius);
                        frac.center[0] = ref.center[0] + offset(rng);
                        frac.center[1] = ref.center[1] + offset(rng);
                        frac.center[2] = ref.center[2] + offset(rng);
                    } else {
                        frac.center = {x_dist(rng), y_dist(rng), z_dist(rng)};
                    }
                }
                break;
            }
            
            case SpatialDistribution::LEVY_FLIGHT: {
                // Lévy flight - heavy-tailed step sizes
                if (fractures.empty()) {
                    frac.center = {x_dist(rng), y_dist(rng), z_dist(rng)};
                } else {
                    auto& prev = fractures.back();
                    // Lévy distribution (approximated with power-law)
                    double step = samplePowerLaw(1.0, volume / 10.0, 1.5, rng);
                    std::uniform_real_distribution<double> angle(0.0, 2.0 * M_PI);
                    double phi = angle(rng);
                    double theta = angle(rng) / 2.0;
                    
                    frac.center[0] = prev.center[0] + step * std::sin(theta) * std::cos(phi);
                    frac.center[1] = prev.center[1] + step * std::sin(theta) * std::sin(phi);
                    frac.center[2] = prev.center[2] + step * std::cos(theta);
                }
                break;
            }
            
            case SpatialDistribution::STRATIGRAPHIC: {
                // Layer-controlled: uniform in x-y, constrained in z
                frac.center[0] = x_dist(rng);
                frac.center[1] = y_dist(rng);
                
                // Concentrate in layer (using cluster_radius as layer thickness)
                double layer_center = (domain_min[2] + domain_max[2]) / 2.0;
                std::normal_distribution<double> z_dist_layer(layer_center, cluster_radius_ / 3.0);
                frac.center[2] = z_dist_layer(rng);
                break;
            }
            
            default:
                frac.center = {x_dist(rng), y_dist(rng), z_dist(rng)};
        }
        
        // Set mechanical properties
        frac.normal_stiffness = normal_stiffness_;
        frac.shear_stiffness = shear_stiffness_;
        frac.friction_angle = friction_angle_;
        frac.cohesion = cohesion_;
        
        fractures.push_back(frac);
    }
    
    return fractures;
}

// ============================================================================
// FractureConnectivityGraph Implementation
// ============================================================================

FractureConnectivityGraph::FractureConnectivityGraph()
    : num_components_(0) {}

void FractureConnectivityGraph::addNode(int fracture_id) {
    if (adjacency_.find(fracture_id) == adjacency_.end()) {
        adjacency_[fracture_id] = std::set<std::pair<int, double>>();
    }
}

void FractureConnectivityGraph::addEdge(int frac1_id, int frac2_id, double weight) {
    addNode(frac1_id);
    addNode(frac2_id);
    
    adjacency_[frac1_id].insert({frac2_id, weight});
    adjacency_[frac2_id].insert({frac1_id, weight});
}

void FractureConnectivityGraph::buildFromFractures(const std::vector<DiscreteFracture>& fractures) {
    // Build graph from fracture connectivity information
    adjacency_.clear();
    component_id_.clear();
    
    for (const auto& frac : fractures) {
        addNode(frac.id);
        
        for (int connected_id : frac.connected_fractures) {
            // Weight based on transmissivity (if available)
            double weight = 1.0;
            addEdge(frac.id, connected_id, weight);
        }
    }
    
    computeComponents();
}

void FractureConnectivityGraph::dfs(int node, int component, std::set<int>& visited) {
    visited.insert(node);
    component_id_[node] = component;
    
    auto it = adjacency_.find(node);
    if (it != adjacency_.end()) {
        for (const auto& neighbor : it->second) {
            if (visited.find(neighbor.first) == visited.end()) {
                dfs(neighbor.first, component, visited);
            }
        }
    }
}

void FractureConnectivityGraph::computeComponents() {
    std::set<int> visited;
    num_components_ = 0;
    component_id_.clear();
    
    for (const auto& node_pair : adjacency_) {
        int node = node_pair.first;
        if (visited.find(node) == visited.end()) {
            dfs(node, num_components_, visited);
            num_components_++;
        }
    }
}

bool FractureConnectivityGraph::isConnected(int frac1_id, int frac2_id) const {
    auto it1 = component_id_.find(frac1_id);
    auto it2 = component_id_.find(frac2_id);
    
    if (it1 == component_id_.end() || it2 == component_id_.end()) {
        return false;
    }
    
    return it1->second == it2->second;
}

std::vector<int> FractureConnectivityGraph::getConnectedComponent(int fracture_id) const {
    std::vector<int> component;
    
    auto it = component_id_.find(fracture_id);
    if (it == component_id_.end()) {
        return component;
    }
    
    int target_component = it->second;
    
    for (const auto& pair : component_id_) {
        if (pair.second == target_component) {
            component.push_back(pair.first);
        }
    }
    
    return component;
}

std::vector<std::vector<int>> FractureConnectivityGraph::getAllComponents() const {
    std::vector<std::vector<int>> components(num_components_);
    
    for (const auto& pair : component_id_) {
        if (pair.second >= 0 && pair.second < num_components_) {
            components[pair.second].push_back(pair.first);
        }
    }
    
    return components;
}

int FractureConnectivityGraph::getLargestComponentSize() const {
    if (num_components_ == 0) return 0;
    
    std::vector<int> sizes(num_components_, 0);
    
    for (const auto& pair : component_id_) {
        if (pair.second >= 0 && pair.second < num_components_) {
            sizes[pair.second]++;
        }
    }
    
    return *std::max_element(sizes.begin(), sizes.end());
}

std::vector<int> FractureConnectivityGraph::findShortestPath(int from_id, int to_id) const {
    // Dijkstra's algorithm
    std::vector<int> path;
    
    if (adjacency_.find(from_id) == adjacency_.end() ||
        adjacency_.find(to_id) == adjacency_.end()) {
        return path;
    }
    
    if (!isConnected(from_id, to_id)) {
        return path;
    }
    
    std::map<int, double> dist;
    std::map<int, int> prev;
    
    // Priority queue: (distance, node)
    std::priority_queue<std::pair<double, int>, 
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> pq;
    
    // Initialize
    for (const auto& node : adjacency_) {
        dist[node.first] = std::numeric_limits<double>::infinity();
    }
    dist[from_id] = 0.0;
    pq.push({0.0, from_id});
    
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        
        if (u == to_id) break;
        
        if (d > dist[u]) continue;
        
        auto it = adjacency_.find(u);
        if (it != adjacency_.end()) {
            for (const auto& [v, weight] : it->second) {
                double alt = dist[u] + weight;
                if (alt < dist[v]) {
                    dist[v] = alt;
                    prev[v] = u;
                    pq.push({alt, v});
                }
            }
        }
    }
    
    // Reconstruct path
    if (prev.find(to_id) != prev.end() || from_id == to_id) {
        int current = to_id;
        while (current != from_id) {
            path.push_back(current);
            auto it = prev.find(current);
            if (it == prev.end()) break;
            current = it->second;
        }
        path.push_back(from_id);
        std::reverse(path.begin(), path.end());
    }
    
    return path;
}

double FractureConnectivityGraph::getPathLength(int from_id, int to_id) const {
    auto path = findShortestPath(from_id, to_id);
    
    if (path.size() < 2) {
        return (from_id == to_id) ? 0.0 : std::numeric_limits<double>::infinity();
    }
    
    double length = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto it = adjacency_.find(path[i]);
        if (it != adjacency_.end()) {
            for (const auto& [neighbor, weight] : it->second) {
                if (neighbor == path[i + 1]) {
                    length += weight;
                    break;
                }
            }
        }
    }
    
    return length;
}

std::vector<std::vector<int>> FractureConnectivityGraph::findAllPaths(
    int from_id, int to_id, int max_depth) const {
    
    std::vector<std::vector<int>> all_paths;
    
    if (adjacency_.find(from_id) == adjacency_.end() ||
        adjacency_.find(to_id) == adjacency_.end()) {
        return all_paths;
    }
    
    // DFS to find all paths
    std::vector<int> current_path;
    std::set<int> visited;
    
    std::function<void(int, int)> dfs_paths = [&](int current, int depth) {
        if (depth > max_depth) return;
        
        current_path.push_back(current);
        visited.insert(current);
        
        if (current == to_id) {
            all_paths.push_back(current_path);
        } else {
            auto it = adjacency_.find(current);
            if (it != adjacency_.end()) {
                for (const auto& [neighbor, weight] : it->second) {
                    if (visited.find(neighbor) == visited.end()) {
                        dfs_paths(neighbor, depth + 1);
                    }
                }
            }
        }
        
        current_path.pop_back();
        visited.erase(current);
    };
    
    dfs_paths(from_id, 0);
    
    return all_paths;
}

double FractureConnectivityGraph::getConnectivity() const {
    if (adjacency_.empty()) return 0.0;
    
    int largest = getLargestComponentSize();
    return static_cast<double>(largest) / adjacency_.size();
}

double FractureConnectivityGraph::getAverageClusteringCoefficient() const {
    if (adjacency_.empty()) return 0.0;
    
    double total_cc = 0.0;
    int count = 0;
    
    for (const auto& [node, neighbors] : adjacency_) {
        int k = neighbors.size();
        if (k < 2) continue;
        
        // Count triangles
        int triangles = 0;
        std::vector<int> neighbor_list;
        for (const auto& [n, w] : neighbors) {
            neighbor_list.push_back(n);
        }
        
        for (size_t i = 0; i < neighbor_list.size(); ++i) {
            for (size_t j = i + 1; j < neighbor_list.size(); ++j) {
                auto it = adjacency_.find(neighbor_list[i]);
                if (it != adjacency_.end()) {
                    for (const auto& [n, w] : it->second) {
                        if (n == neighbor_list[j]) {
                            triangles++;
                            break;
                        }
                    }
                }
            }
        }
        
        // Local clustering coefficient
        double cc = 2.0 * triangles / (k * (k - 1));
        total_cc += cc;
        count++;
    }
    
    return (count > 0) ? total_cc / count : 0.0;
}

double FractureConnectivityGraph::getAveragePathLength() const {
    if (adjacency_.size() < 2) return 0.0;
    
    double total_length = 0.0;
    int count = 0;
    
    std::vector<int> nodes;
    for (const auto& [node, neighbors] : adjacency_) {
        nodes.push_back(node);
    }
    
    // Sample pairs if too many nodes
    size_t max_samples = std::min(nodes.size() * nodes.size(), size_t(1000));
    
    for (size_t i = 0; i < nodes.size() && count < (int)max_samples; ++i) {
        for (size_t j = i + 1; j < nodes.size() && count < (int)max_samples; ++j) {
            double length = getPathLength(nodes[i], nodes[j]);
            if (length < std::numeric_limits<double>::infinity()) {
                total_length += length;
                count++;
            }
        }
    }
    
    return (count > 0) ? total_length / count : 0.0;
}

int FractureConnectivityGraph::getNumComponents() const {
    return num_components_;
}

FractureConnectivityGraph::PercolationResult FractureConnectivityGraph::analyzePercolation(
    const std::vector<DiscreteFracture>& fractures,
    const std::array<double, 3>& domain_min,
    const std::array<double, 3>& domain_max) const {
    
    PercolationResult result;
    result.percolates_x = false;
    result.percolates_y = false;
    result.percolates_z = false;
    result.percolation_probability = 0.0;
    
    if (fractures.empty()) return result;
    
    // Find fractures touching each boundary
    std::set<int> x_min_fracs, x_max_fracs;
    std::set<int> y_min_fracs, y_max_fracs;
    std::set<int> z_min_fracs, z_max_fracs;
    
    double tol = 0.01 * (domain_max[0] - domain_min[0]);
    
    for (const auto& frac : fractures) {
        auto vertices = frac.getVertices(8);
        
        for (const auto& v : vertices) {
            if (std::abs(v[0] - domain_min[0]) < tol) x_min_fracs.insert(frac.id);
            if (std::abs(v[0] - domain_max[0]) < tol) x_max_fracs.insert(frac.id);
            if (std::abs(v[1] - domain_min[1]) < tol) y_min_fracs.insert(frac.id);
            if (std::abs(v[1] - domain_max[1]) < tol) y_max_fracs.insert(frac.id);
            if (std::abs(v[2] - domain_min[2]) < tol) z_min_fracs.insert(frac.id);
            if (std::abs(v[2] - domain_max[2]) < tol) z_max_fracs.insert(frac.id);
        }
    }
    
    // Check percolation in each direction
    auto checkPercolation = [this](const std::set<int>& inlet, const std::set<int>& outlet) {
        for (int in_id : inlet) {
            for (int out_id : outlet) {
                if (isConnected(in_id, out_id)) {
                    return true;
                }
            }
        }
        return false;
    };
    
    result.percolates_x = checkPercolation(x_min_fracs, x_max_fracs);
    result.percolates_y = checkPercolation(y_min_fracs, y_max_fracs);
    result.percolates_z = checkPercolation(z_min_fracs, z_max_fracs);
    
    // Compute percolation probability as fraction of directions that percolate
    int perc_count = (result.percolates_x ? 1 : 0) + 
                     (result.percolates_y ? 1 : 0) + 
                     (result.percolates_z ? 1 : 0);
    result.percolation_probability = perc_count / 3.0;
    
    // Extract backbone fractures (those in percolating path)
    if (result.percolates_x || result.percolates_y || result.percolates_z) {
        std::set<int> backbone_set;
        
        auto addBackbone = [&](const std::set<int>& inlet, const std::set<int>& outlet) {
            for (int in_id : inlet) {
                for (int out_id : outlet) {
                    if (isConnected(in_id, out_id)) {
                        auto path = findShortestPath(in_id, out_id);
                        for (int id : path) {
                            backbone_set.insert(id);
                        }
                    }
                }
            }
        };
        
        if (result.percolates_x) addBackbone(x_min_fracs, x_max_fracs);
        if (result.percolates_y) addBackbone(y_min_fracs, y_max_fracs);
        if (result.percolates_z) addBackbone(z_min_fracs, z_max_fracs);
        
        result.backbone_fractures.assign(backbone_set.begin(), backbone_set.end());
    }
    
    return result;
}

std::vector<int> FractureConnectivityGraph::extractBackbone(int inlet_id, int outlet_id) const {
    // Find all fractures that lie on any shortest path between inlet and outlet
    std::vector<int> backbone;
    
    if (!isConnected(inlet_id, outlet_id)) {
        return backbone;
    }
    
    // Find shortest path length
    auto shortest = findShortestPath(inlet_id, outlet_id);
    double min_length = getPathLength(inlet_id, outlet_id);
    
    // Find all fractures that can be on a shortest path
    std::set<int> backbone_set;
    
    for (const auto& [node, neighbors] : adjacency_) {
        double d_from_inlet = getPathLength(inlet_id, node);
        double d_to_outlet = getPathLength(node, outlet_id);
        
        // Node is on backbone if d(inlet, node) + d(node, outlet) == d(inlet, outlet)
        if (std::abs(d_from_inlet + d_to_outlet - min_length) < 1e-10) {
            backbone_set.insert(node);
        }
    }
    
    backbone.assign(backbone_set.begin(), backbone_set.end());
    return backbone;
}

// ============================================================================
// FractureNetwork Implementation
// ============================================================================

FractureNetwork::FractureNetwork()
    : domain_min_({0.0, 0.0, 0.0}),
      domain_max_({100.0, 100.0, 100.0}),
      domain_volume_(1e6),
      rng_(42) {}

void FractureNetwork::setDomain(const std::array<double, 3>& min_corner,
                                 const std::array<double, 3>& max_corner) {
    domain_min_ = min_corner;
    domain_max_ = max_corner;
    domain_volume_ = (max_corner[0] - min_corner[0]) *
                     (max_corner[1] - min_corner[1]) *
                     (max_corner[2] - min_corner[2]);
}

void FractureNetwork::setDomainFromMesh(DM dm) {
    // Extract domain bounds from PETSc DM
    PetscReal gmin[3], gmax[3];
    DMGetBoundingBox(dm, gmin, gmax);
    
    domain_min_ = {gmin[0], gmin[1], gmin[2]};
    domain_max_ = {gmax[0], gmax[1], gmax[2]};
    domain_volume_ = (gmax[0] - gmin[0]) * (gmax[1] - gmin[1]) * (gmax[2] - gmin[2]);
}

void FractureNetwork::addFractureSet(const FractureSet& set) {
    sets_.push_back(set);
}

FractureSet* FractureNetwork::getFractureSet(int set_id) {
    for (auto& set : sets_) {
        if (set.getId() == set_id) {
            return &set;
        }
    }
    return nullptr;
}

void FractureNetwork::setRandomSeed(unsigned int seed) {
    rng_.seed(seed);
}

void FractureNetwork::generate() {
    fractures_.clear();
    int current_id = 0;
    
    for (auto& set : sets_) {
        auto set_fractures = set.generateFractures(domain_min_, domain_max_, rng_);
        
        // Renumber IDs
        for (auto& frac : set_fractures) {
            frac.id = current_id++;
        }
        
        fractures_.insert(fractures_.end(), set_fractures.begin(), set_fractures.end());
    }
    
    // Clip to domain
    clipToDomain();
    
    // Compute intersections
    computeIntersections();
    
    // Apply termination rules
    applyTerminationRules();
    
    // Build connectivity graph
    buildConnectivityGraph();
}

void FractureNetwork::generateWithIntensityControl(double target_p32, double tolerance) {
    // Iteratively adjust intensity to match target P32
    
    // Initial generation
    generate();
    
    FractureIntensity intensity = computeIntensity();
    double current_p32 = intensity.P32;
    
    int max_iterations = 10;
    int iter = 0;
    
    while (std::abs(current_p32 - target_p32) / target_p32 > tolerance && iter < max_iterations) {
        // Adjust set intensities
        double ratio = target_p32 / current_p32;
        
        for (auto& set : sets_) {
            // Approximate adjustment (this is simplified)
            double current_intensity = set.getId() < (int)sets_.size() ? 
                current_p32 / sets_.size() : 0.1;
            set.setIntensity(current_intensity * ratio);
        }
        
        // Regenerate
        generate();
        
        intensity = computeIntensity();
        current_p32 = intensity.P32;
        iter++;
    }
}

void FractureNetwork::clipToDomain() {
    // Remove fractures completely outside domain
    // and clip those partially outside
    
    auto it = fractures_.begin();
    while (it != fractures_.end()) {
        bool inside = true;
        
        // Check if center is inside
        for (int i = 0; i < 3; ++i) {
            if (it->center[i] < domain_min_[i] - it->radius ||
                it->center[i] > domain_max_[i] + it->radius) {
                inside = false;
                break;
            }
        }
        
        if (!inside) {
            it = fractures_.erase(it);
        } else {
            ++it;
        }
    }
}

void FractureNetwork::applyTerminationRules() {
    // Apply termination rules from fracture sets
    // This is a simplified implementation
    
    for ([[maybe_unused]] auto& set : sets_) {
        // Check termination type through properties
        // Full implementation would track which fractures terminate at intersections
    }
}

DiscreteFracture* FractureNetwork::getFracture(int id) {
    for (auto& frac : fractures_) {
        if (frac.id == id) {
            return &frac;
        }
    }
    return nullptr;
}

const DiscreteFracture* FractureNetwork::getFracture(int id) const {
    for (const auto& frac : fractures_) {
        if (frac.id == id) {
            return &frac;
        }
    }
    return nullptr;
}

bool FractureNetwork::detectIntersection(const DiscreteFracture& f1, 
                                          const DiscreteFracture& f2,
                                          FractureIntersection& intersection) const {
    // Detect if two fracture planes intersect within their extents
    
    // Compute intersection line direction (cross product of normals)
    std::array<double, 3> line_dir;
    line_dir[0] = f1.normal[1] * f2.normal[2] - f1.normal[2] * f2.normal[1];
    line_dir[1] = f1.normal[2] * f2.normal[0] - f1.normal[0] * f2.normal[2];
    line_dir[2] = f1.normal[0] * f2.normal[1] - f1.normal[1] * f2.normal[0];
    
    double line_len = std::sqrt(line_dir[0]*line_dir[0] + 
                                line_dir[1]*line_dir[1] + 
                                line_dir[2]*line_dir[2]);
    
    // Parallel planes don't intersect
    if (line_len < 1e-10) {
        return false;
    }
    
    // Normalize
    line_dir[0] /= line_len;
    line_dir[1] /= line_len;
    line_dir[2] /= line_len;
    
    // Find a point on the intersection line
    // Solve system: n1 · (p - c1) = 0 and n2 · (p - c2) = 0
    
    // Distance from plane 1 to origin
    double d1 = f1.normal[0] * f1.center[0] + 
                f1.normal[1] * f1.center[1] + 
                f1.normal[2] * f1.center[2];
    
    // Distance from plane 2 to origin
    double d2 = f2.normal[0] * f2.center[0] + 
                f2.normal[1] * f2.center[1] + 
                f2.normal[2] * f2.center[2];
    
    // Find point on intersection line closest to origin
    double n1n2 = f1.normal[0] * f2.normal[0] + 
                  f1.normal[1] * f2.normal[1] + 
                  f1.normal[2] * f2.normal[2];
    double denom = 1.0 - n1n2 * n1n2;
    
    if (std::abs(denom) < 1e-10) {
        return false;
    }
    
    double c1 = (d1 - d2 * n1n2) / denom;
    double c2 = (d2 - d1 * n1n2) / denom;
    
    std::array<double, 3> line_point;
    line_point[0] = c1 * f1.normal[0] + c2 * f2.normal[0];
    line_point[1] = c1 * f1.normal[1] + c2 * f2.normal[1];
    line_point[2] = c1 * f1.normal[2] + c2 * f2.normal[2];
    
    // Find intersection segment within both fracture extents
    // Using disc approximation
    
    auto findIntersectionWithDisc = [](const std::array<double, 3>& point,
                                        const std::array<double, 3>& dir,
                                        const DiscreteFracture& frac,
                                        double& t_min, double& t_max) {
        // Find t values where line intersects disc boundary
        // Project to fracture plane coordinates
        
        double dx = point[0] - frac.center[0];
        double dy = point[1] - frac.center[1];
        double dz = point[2] - frac.center[2];
        
        // Project onto plane (remove normal component)
        double n_dot = dx * frac.normal[0] + dy * frac.normal[1] + dz * frac.normal[2];
        dx -= n_dot * frac.normal[0];
        dy -= n_dot * frac.normal[1];
        dz -= n_dot * frac.normal[2];
        
        double dir_proj[3] = {dir[0], dir[1], dir[2]};
        double dir_n = dir[0] * frac.normal[0] + dir[1] * frac.normal[1] + dir[2] * frac.normal[2];
        dir_proj[0] -= dir_n * frac.normal[0];
        dir_proj[1] -= dir_n * frac.normal[1];
        dir_proj[2] -= dir_n * frac.normal[2];
        
        // Solve |p + t*d|² = r²
        double a = dir_proj[0]*dir_proj[0] + dir_proj[1]*dir_proj[1] + dir_proj[2]*dir_proj[2];
        double b = 2.0 * (dx*dir_proj[0] + dy*dir_proj[1] + dz*dir_proj[2]);
        double c = dx*dx + dy*dy + dz*dz - frac.radius * frac.radius;
        
        double disc = b*b - 4.0*a*c;
        
        if (disc < 0 || std::abs(a) < 1e-10) {
            return false;
        }
        
        t_min = (-b - std::sqrt(disc)) / (2.0 * a);
        t_max = (-b + std::sqrt(disc)) / (2.0 * a);
        
        return true;
    };
    
    double t1_min, t1_max, t2_min, t2_max;
    
    if (!findIntersectionWithDisc(line_point, line_dir, f1, t1_min, t1_max) ||
        !findIntersectionWithDisc(line_point, line_dir, f2, t2_min, t2_max)) {
        return false;
    }
    
    // Find overlap
    double t_start = std::max(t1_min, t2_min);
    double t_end = std::min(t1_max, t2_max);
    
    if (t_start >= t_end) {
        return false;
    }
    
    // Compute intersection endpoints
    intersection.fracture1_id = f1.id;
    intersection.fracture2_id = f2.id;
    
    intersection.point1[0] = line_point[0] + t_start * line_dir[0];
    intersection.point1[1] = line_point[1] + t_start * line_dir[1];
    intersection.point1[2] = line_point[2] + t_start * line_dir[2];
    
    intersection.point2[0] = line_point[0] + t_end * line_dir[0];
    intersection.point2[1] = line_point[1] + t_end * line_dir[1];
    intersection.point2[2] = line_point[2] + t_end * line_dir[2];
    
    double dx = intersection.point2[0] - intersection.point1[0];
    double dy = intersection.point2[1] - intersection.point1[1];
    double dz = intersection.point2[2] - intersection.point1[2];
    intersection.length = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Compute intersection transmissivity (harmonic mean)
    if (f1.transmissivity > 0 && f2.transmissivity > 0) {
        intersection.transmissivity = 2.0 * f1.transmissivity * f2.transmissivity / 
                                       (f1.transmissivity + f2.transmissivity);
    }
    
    return intersection.length > 1e-6;
}

void FractureNetwork::computeIntersections() {
    intersections_.clear();
    
    // O(n²) pairwise check - could be optimized with spatial indexing
    for (size_t i = 0; i < fractures_.size(); ++i) {
        for (size_t j = i + 1; j < fractures_.size(); ++j) {
            FractureIntersection intersection;
            if (detectIntersection(fractures_[i], fractures_[j], intersection)) {
                intersections_.push_back(intersection);
                
                // Update fracture connectivity
                fractures_[i].connected_fractures.push_back(fractures_[j].id);
                fractures_[j].connected_fractures.push_back(fractures_[i].id);
            }
        }
    }
}

void FractureNetwork::buildConnectivityGraph() {
    connectivity_ = FractureConnectivityGraph();
    
    // Add all fractures as nodes
    for (const auto& frac : fractures_) {
        connectivity_.addNode(frac.id);
    }
    
    // Add edges from intersections
    for (const auto& intersection : intersections_) {
        // Weight by intersection transmissivity
        double weight = 1.0 / std::max(intersection.transmissivity, 1e-20);
        connectivity_.addEdge(intersection.fracture1_id, intersection.fracture2_id, weight);
    }
    
    // Build from fractures will recompute components
    connectivity_.buildFromFractures(fractures_);
}

FractureIntensity FractureNetwork::computeIntensity() const {
    FractureIntensity intensity;
    
    if (fractures_.empty() || domain_volume_ <= 0) {
        return intensity;
    }
    
    double total_area = 0.0;
    double total_volume = 0.0;
    
    for (const auto& frac : fractures_) {
        double area = frac.computeArea();
        total_area += area;
        total_volume += area * frac.aperture;
    }
    
    intensity.P30 = fractures_.size() / domain_volume_;
    intensity.P32 = total_area / domain_volume_;
    intensity.P33 = total_volume / domain_volume_;
    
    // Approximate P10, P21 from P32
    double mean_radius = std::sqrt(total_area / (M_PI * fractures_.size()));
    intensity.P10 = intensity.P30 * 2.0 * mean_radius;
    intensity.P21 = intensity.P32 / mean_radius;
    intensity.P20 = intensity.P30 * mean_radius;
    
    return intensity;
}

double FractureNetwork::computePercolationThreshold(int n_realizations) {
    // Monte Carlo estimation of percolation threshold
    // Vary P32 and find critical value where percolation probability = 0.5
    
    if (sets_.empty()) {
        return 0.0;
    }
    
    double p32_min = 0.01;
    double p32_max = 1.0;
    double target_prob = 0.5;
    
    auto computePercolationProbability = [&](double p32) {
        int perc_count = 0;
        
        for (int i = 0; i < n_realizations; ++i) {
            // Temporarily modify intensity
            for (auto& set : sets_) {
                set.setIntensity(p32 / sets_.size());
            }
            
            // Generate and check percolation
            FractureNetwork temp_net;
            temp_net.setDomain(domain_min_, domain_max_);
            temp_net.setRandomSeed(i);
            for (const auto& set : sets_) {
                temp_net.addFractureSet(set);
            }
            temp_net.generate();
            
            auto result = temp_net.connectivity_.analyzePercolation(
                temp_net.fractures_, domain_min_, domain_max_);
            
            if (result.percolates_x || result.percolates_y || result.percolates_z) {
                perc_count++;
            }
        }
        
        return static_cast<double>(perc_count) / n_realizations;
    };
    
    // Binary search for threshold
    int max_iter = 20;
    for (int iter = 0; iter < max_iter; ++iter) {
        double p32_mid = (p32_min + p32_max) / 2.0;
        double prob = computePercolationProbability(p32_mid);
        
        if (prob < target_prob) {
            p32_min = p32_mid;
        } else {
            p32_max = p32_mid;
        }
        
        if (std::abs(prob - target_prob) < 0.05) {
            break;
        }
    }
    
    return (p32_min + p32_max) / 2.0;
}

FractureNetwork::MINCResult FractureNetwork::computeMINC(int n_continua) const {
    MINCResult result;
    
    if (n_continua < 2) n_continua = 2;
    
    result.volume_fractions.resize(n_continua);
    result.shape_factors.resize(n_continua);
    result.transmissibilities.resize(n_continua, std::vector<double>(n_continua, 0.0));
    
    // Estimate matrix block size from fracture spacing
    FractureIntensity intensity = computeIntensity();
    double fracture_spacing = (intensity.P32 > 0) ? 1.0 / intensity.P32 : 10.0;
    
    // Nested continua with decreasing volumes
    // Volume fractions: geometric series
    double total = 0.0;
    for (int i = 0; i < n_continua; ++i) {
        result.volume_fractions[i] = std::pow(0.5, i);
        total += result.volume_fractions[i];
    }
    for (int i = 0; i < n_continua; ++i) {
        result.volume_fractions[i] /= total;
    }
    
    // Shape factors for each continuum (Kazemi-type)
    // sigma = 4*(1/L1² + 1/L2² + 1/L3²) for cubic blocks
    for (int i = 0; i < n_continua; ++i) {
        double L_i = fracture_spacing * std::pow(0.5, i);  // Decreasing block size
        result.shape_factors[i] = 12.0 / (L_i * L_i);
    }
    
    // Transmissibilities between adjacent continua
    double mean_perm = 0.0;
    for (const auto& frac : fractures_) {
        mean_perm += frac.permeability;
    }
    if (!fractures_.empty()) {
        mean_perm /= fractures_.size();
    }
    
    for (int i = 0; i < n_continua - 1; ++i) {
        // T_i,i+1 = (k_i * A) / d
        double L_i = fracture_spacing * std::pow(0.5, i);
        double L_next = fracture_spacing * std::pow(0.5, i + 1);
        double d = (L_i + L_next) / 4.0;  // Distance between continuum centers
        
        result.transmissibilities[i][i + 1] = mean_perm / d;
        result.transmissibilities[i + 1][i] = result.transmissibilities[i][i + 1];
    }
    
    return result;
}

std::vector<FractureNetwork::EDFMConnection> FractureNetwork::computeEDFM(DM dm) const {
    std::vector<EDFMConnection> connections;
    
    // Get mesh information
    PetscInt dim;
    DMGetDimension(dm, &dim);
    
    // For each fracture, find intersecting matrix cells
    // This is a simplified implementation
    
    for (const auto& frac : fractures_) {
        // Create fracture-matrix connection
        EDFMConnection fm_conn;
        fm_conn.type = EDFMConnection::FRACTURE_MATRIX;
        fm_conn.id1 = frac.id;
        fm_conn.id2 = -1;  // Would be matrix cell ID
        
        // EDFM transmissibility: T_fm = k_m * A_f / d_f
        // where d_f is average normal distance
        double cell_size = std::pow(domain_volume_ / 1000.0, 1.0/3.0);  // Approximate
        fm_conn.distance = cell_size / 4.0;  // Average distance to fracture
        fm_conn.area = frac.computeArea();
        
        // Assume matrix permeability of 1e-15 m²
        double k_matrix = 1e-15;
        fm_conn.transmissibility = k_matrix * fm_conn.area / fm_conn.distance;
        
        connections.push_back(fm_conn);
    }
    
    // Add fracture-fracture connections from intersections
    for (const auto& intersection : intersections_) {
        EDFMConnection ff_conn;
        ff_conn.type = EDFMConnection::FRACTURE_FRACTURE;
        ff_conn.id1 = intersection.fracture1_id;
        ff_conn.id2 = intersection.fracture2_id;
        ff_conn.distance = 0.0;  // Fractures intersect
        ff_conn.area = intersection.length * std::min(
            getFracture(intersection.fracture1_id) ? 
                getFracture(intersection.fracture1_id)->aperture : 1e-4,
            getFracture(intersection.fracture2_id) ? 
                getFracture(intersection.fracture2_id)->aperture : 1e-4
        );
        ff_conn.transmissibility = intersection.transmissivity * intersection.length;
        
        connections.push_back(ff_conn);
    }
    
    return connections;
}

std::array<std::array<double, 3>, 3> FractureNetwork::computeOdaTensor() const {
    // Oda's crack tensor for equivalent permeability
    // k_ij = (1/V) * Σ (b³/12) * A * (δ_ij - n_i * n_j)
    
    std::array<std::array<double, 3>, 3> tensor = {{{0,0,0}, {0,0,0}, {0,0,0}}};
    
    if (fractures_.empty() || domain_volume_ <= 0) {
        return tensor;
    }
    
    for (const auto& frac : fractures_) {
        double area = frac.computeArea();
        double b3_12 = frac.aperture * frac.aperture * frac.aperture / 12.0;
        double factor = b3_12 * area / domain_volume_;
        
        // k_ij += factor * (δ_ij - n_i * n_j)
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double delta_ij = (i == j) ? 1.0 : 0.0;
                tensor[i][j] += factor * (delta_ij - frac.normal[i] * frac.normal[j]);
            }
        }
    }
    
    return tensor;
}

std::array<double, 3> FractureNetwork::computeSnowPermeability() const {
    // Snow's model: parallel plate assumption
    // k_i = (1/V) * Σ (b³/12) * A * cos²(θ_i)
    // where θ_i is angle between fracture normal and i-th axis
    
    std::array<double, 3> perm = {0.0, 0.0, 0.0};
    
    if (fractures_.empty() || domain_volume_ <= 0) {
        return perm;
    }
    
    for (const auto& frac : fractures_) {
        double area = frac.computeArea();
        double b3_12 = frac.aperture * frac.aperture * frac.aperture / 12.0;
        double factor = b3_12 * area / domain_volume_;
        
        // Contribution to each direction
        // Flow in x-direction: perpendicular to normal's x-component
        // k_x += factor * (1 - n_x²) = factor * (n_y² + n_z²)
        perm[0] += factor * (1.0 - frac.normal[0] * frac.normal[0]);
        perm[1] += factor * (1.0 - frac.normal[1] * frac.normal[1]);
        perm[2] += factor * (1.0 - frac.normal[2] * frac.normal[2]);
    }
    
    return perm;
}

void FractureNetwork::updateApertures(const std::vector<double>& stress_field) {
    // Update fracture apertures based on stress state
    
    if (stress_field.size() < fractures_.size() * 6) {
        return;  // Need full stress tensor for each fracture
    }
    
    for (size_t i = 0; i < fractures_.size(); ++i) {
        auto& frac = fractures_[i];
        
        // Extract stress tensor at fracture location (simplified: use index)
        size_t idx = i * 6;
        double sxx = stress_field[idx];
        double syy = stress_field[idx + 1];
        double szz = stress_field[idx + 2];
        double sxy = stress_field[idx + 3];
        double sxz = stress_field[idx + 4];
        double syz = stress_field[idx + 5];
        
        // Compute normal stress on fracture plane
        // σn = n^T · σ · n
        double sigma_n = sxx * frac.normal[0] * frac.normal[0] +
                         syy * frac.normal[1] * frac.normal[1] +
                         szz * frac.normal[2] * frac.normal[2] +
                         2.0 * sxy * frac.normal[0] * frac.normal[1] +
                         2.0 * sxz * frac.normal[0] * frac.normal[2] +
                         2.0 * syz * frac.normal[1] * frac.normal[2];
        
        // Update aperture using stress-dependent model
        frac.current_aperture = frac.computeAperture(sigma_n, frac.aperture, 10e6);
    }
}

void FractureNetwork::updatePermeabilities() {
    for (auto& frac : fractures_) {
        // Store original aperture for reference
        double original_aperture = frac.aperture;
        
        // Update aperture to current stress-dependent value
        frac.aperture = frac.current_aperture;
        frac.computePermeability();
        
        // Restore original (unstressed) aperture
        frac.aperture = original_aperture;
    }
}

void FractureNetwork::importFromFile(const std::string& filename, const std::string& format) {
    std::string fmt = format;
    
    // Auto-detect format from extension
    if (fmt == "auto") {
        size_t dot_pos = filename.rfind('.');
        if (dot_pos != std::string::npos) {
            fmt = filename.substr(dot_pos + 1);
        }
    }
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    fractures_.clear();
    
    if (fmt == "csv" || fmt == "txt") {
        // Simple CSV format: id,x,y,z,strike,dip,radius,aperture
        std::string line;
        std::getline(file, line);  // Skip header
        
        int id = 0;
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string token;
            
            DiscreteFracture frac;
            frac.id = id++;
            
            std::getline(ss, token, ','); frac.center[0] = std::stod(token);
            std::getline(ss, token, ','); frac.center[1] = std::stod(token);
            std::getline(ss, token, ','); frac.center[2] = std::stod(token);
            std::getline(ss, token, ','); frac.strike = std::stod(token);
            std::getline(ss, token, ','); frac.dip = std::stod(token);
            std::getline(ss, token, ','); frac.radius = std::stod(token);
            std::getline(ss, token, ','); frac.aperture = std::stod(token);
            
            frac.computeVectors();
            frac.computePermeability();
            frac.length = 2.0 * frac.radius;
            frac.height = 2.0 * frac.radius;
            
            fractures_.push_back(frac);
        }
    }
    
    file.close();
    
    // Compute intersections and connectivity
    computeIntersections();
    buildConnectivityGraph();
}

void FractureNetwork::exportToFile(const std::string& filename, const std::string& format) const {
    if (format == "vtk") {
        writeVTK(filename);
    } else if (format == "csv" || format == "txt") {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        file << "id,x,y,z,strike,dip,radius,aperture,permeability\n";
        
        for (const auto& frac : fractures_) {
            file << frac.id << ","
                 << frac.center[0] << "," << frac.center[1] << "," << frac.center[2] << ","
                 << frac.strike << "," << frac.dip << ","
                 << frac.radius << "," << frac.aperture << "," << frac.permeability << "\n";
        }
        
        file.close();
    } else if (format == "gmsh") {
        exportToGmsh(filename);
    }
}

void FractureNetwork::exportToGmsh(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << "// Gmsh geometry file for fracture network\n";
    file << "// Generated by FSRM FractureNetwork\n\n";
    
    file << "lc = 1.0;\n\n";
    
    int point_id = 1;
    int line_id = 1;
    int surface_id = 1;
    
    for (const auto& frac : fractures_) {
        auto vertices = frac.getVertices(8);
        
        file << "// Fracture " << frac.id << "\n";
        
        std::vector<int> point_ids;
        for (const auto& v : vertices) {
            file << "Point(" << point_id << ") = {" 
                 << v[0] << ", " << v[1] << ", " << v[2] << ", lc};\n";
            point_ids.push_back(point_id++);
        }
        
        std::vector<int> line_ids;
        for (size_t i = 0; i < point_ids.size(); ++i) {
            int p1 = point_ids[i];
            int p2 = point_ids[(i + 1) % point_ids.size()];
            file << "Line(" << line_id << ") = {" << p1 << ", " << p2 << "};\n";
            line_ids.push_back(line_id++);
        }
        
        file << "Line Loop(" << surface_id << ") = {";
        for (size_t i = 0; i < line_ids.size(); ++i) {
            if (i > 0) file << ", ";
            file << line_ids[i];
        }
        file << "};\n";
        
        file << "Plane Surface(" << surface_id << ") = {" << surface_id << "};\n\n";
        surface_id++;
    }
    
    file.close();
}

void FractureNetwork::writeVTK(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << "# vtk DataFile Version 3.0\n";
    file << "Fracture Network\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    
    // Count total vertices
    int total_vertices = 0;
    int n_sides = 8;
    for (const auto& frac : fractures_) {
        total_vertices += n_sides;
    }
    
    // Write points
    file << "POINTS " << total_vertices << " double\n";
    
    for (const auto& frac : fractures_) {
        (void)frac;  // Suppress unused warning - placeholder for full implementation
        auto vertices = frac.getVertices(n_sides);
        for (const auto& v : vertices) {
            file << v[0] << " " << v[1] << " " << v[2] << "\n";
        }
    }
    
    // Write polygons
    file << "POLYGONS " << fractures_.size() << " " 
         << fractures_.size() * (n_sides + 1) << "\n";
    
    int vertex_offset = 0;
    for (size_t i = 0; i < fractures_.size(); ++i) {
        file << n_sides;
        for (int j = 0; j < n_sides; ++j) {
            file << " " << (vertex_offset + j);
        }
        file << "\n";
        vertex_offset += n_sides;
    }
    
    // Write cell data
    file << "CELL_DATA " << fractures_.size() << "\n";
    
    file << "SCALARS aperture double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& frac : fractures_) {
        file << frac.aperture << "\n";
    }
    
    file << "SCALARS permeability double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& frac : fractures_) {
        file << frac.permeability << "\n";
    }
    
    file << "SCALARS set_id int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& frac : fractures_) {
        file << frac.set_id << "\n";
    }
    
    file << "VECTORS normal double\n";
    for (const auto& frac : fractures_) {
        file << frac.normal[0] << " " << frac.normal[1] << " " << frac.normal[2] << "\n";
    }
    
    file.close();
}

void FractureNetwork::writeFractureTraces(const std::string& filename, double z_level) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    file << "# Fracture traces at z = " << z_level << "\n";
    file << "# x1, y1, x2, y2, fracture_id, aperture\n";
    
    for (const auto& frac : fractures_) {
        // Check if fracture intersects z-plane
        double z_dist = z_level - frac.center[2];
        
        // Distance from center to plane along normal
        double d = std::abs(z_dist * frac.normal[2]);
        
        if (d > frac.radius) {
            continue;  // Fracture doesn't reach this z-level
        }
        
        // Compute trace (intersection with horizontal plane)
        // Simplified: approximate as line segment
        double trace_half_length = std::sqrt(frac.radius * frac.radius - d * d);
        
        // Trace direction is perpendicular to both fracture normal and z-axis
        double tx = frac.normal[1];
        double ty = -frac.normal[0];
        double t_norm = std::sqrt(tx*tx + ty*ty);
        
        if (t_norm < 1e-10) {
            // Horizontal fracture - draw as circle approximation
            tx = 1.0; ty = 0.0;
        } else {
            tx /= t_norm;
            ty /= t_norm;
        }
        
        // Trace center (projection of fracture center onto z-plane)
        double cx = frac.center[0] + z_dist * frac.normal[0] / frac.normal[2];
        double cy = frac.center[1] + z_dist * frac.normal[1] / frac.normal[2];
        
        double x1 = cx - trace_half_length * tx;
        double y1 = cy - trace_half_length * ty;
        double x2 = cx + trace_half_length * tx;
        double y2 = cy + trace_half_length * ty;
        
        file << x1 << ", " << y1 << ", " << x2 << ", " << y2 << ", "
             << frac.id << ", " << frac.aperture << "\n";
    }
    
    file.close();
}

// ============================================================================
// FracturePropagationModel Implementation
// ============================================================================

FracturePropagationModel::FracturePropagationModel()
    : youngs_modulus_(30e9),
      poisson_ratio_(0.25),
      fracture_toughness_(1e6),
      fluid_viscosity_(0.001),
      fluid_compressibility_(4.5e-10),
      far_field_stress_({{{0,0,0}, {0,0,0}, {0,0,0}}}) {}

void FracturePropagationModel::setRockProperties(double E, double nu, double KIC) {
    youngs_modulus_ = E;
    poisson_ratio_ = nu;
    fracture_toughness_ = KIC;
}

void FracturePropagationModel::setFluidProperties(double viscosity, double compressibility) {
    fluid_viscosity_ = viscosity;
    fluid_compressibility_ = compressibility;
}

void FracturePropagationModel::setFarFieldStress(double Sxx, double Syy, double Szz,
                                                  double Sxy, double Sxz, double Syz) {
    far_field_stress_[0][0] = Sxx;
    far_field_stress_[1][1] = Syy;
    far_field_stress_[2][2] = Szz;
    far_field_stress_[0][1] = far_field_stress_[1][0] = Sxy;
    far_field_stress_[0][2] = far_field_stress_[2][0] = Sxz;
    far_field_stress_[1][2] = far_field_stress_[2][1] = Syz;
}

double FracturePropagationModel::computeKI(const DiscreteFracture& frac, double net_pressure) const {
    // Mode I stress intensity factor for penny-shaped crack
    // K_I = 2 * p_net * sqrt(a / π)
    // where a is the crack radius
    
    double a = frac.radius;
    return 2.0 * net_pressure * std::sqrt(a / M_PI);
}

std::array<double, 3> FracturePropagationModel::computePropagationDirection(
    const TipState& state) const {
    
    // Maximum circumferential stress criterion
    // Propagation direction based on mixed-mode SIF
    
    std::array<double, 3> dir = state.propagation_direction;
    
    if (std::abs(state.K_II) > 1e-10) {
        // Angle of propagation: θ = 2 * arctan((K_I - sqrt(K_I² + 8*K_II²)) / (4*K_II))
        double theta = 2.0 * std::atan2(
            state.K_I - std::sqrt(state.K_I * state.K_I + 8.0 * state.K_II * state.K_II),
            4.0 * state.K_II);
        
        // Rotate propagation direction by theta (simplified 2D rotation)
        double cos_t = std::cos(theta);
        double sin_t = std::sin(theta);
        
        double new_x = dir[0] * cos_t - dir[1] * sin_t;
        double new_y = dir[0] * sin_t + dir[1] * cos_t;
        
        dir[0] = new_x;
        dir[1] = new_y;
        
        // Normalize
        double norm = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        if (norm > 1e-10) {
            dir[0] /= norm;
            dir[1] /= norm;
            dir[2] /= norm;
        }
    }
    
    return dir;
}

std::vector<FracturePropagationModel::TipState> FracturePropagationModel::computeTipStates(
    const FractureNetwork& network) const {
    
    std::vector<TipState> states;
    
    for (const auto& frac : network.getFractures()) {
        // Compute normal stress on fracture from far-field
        double sigma_n = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                sigma_n += far_field_stress_[i][j] * frac.normal[i] * frac.normal[j];
            }
        }
        
        // Net pressure (assume some internal pressure)
        double p_internal = sigma_n + 1e6;  // Would come from flow solution
        double net_pressure = p_internal - sigma_n;
        
        TipState state;
        state.fracture_id = frac.id;
        
        // Multiple tip positions around fracture perimeter
        // For simplicity, use one representative tip
        state.tip_position = frac.center;
        state.tip_position[0] += frac.radius * frac.strike_dir[0];
        state.tip_position[1] += frac.radius * frac.strike_dir[1];
        state.tip_position[2] += frac.radius * frac.strike_dir[2];
        
        // Propagation direction (radially outward)
        state.propagation_direction = frac.strike_dir;
        
        // Compute stress intensity factors
        state.K_I = computeKI(frac, net_pressure);
        
        // Mode II from shear stress (simplified)
        double tau = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                tau += far_field_stress_[i][j] * frac.normal[i] * frac.strike_dir[j];
            }
        }
        state.K_II = tau * std::sqrt(M_PI * frac.radius);
        
        // Mode III from anti-plane shear
        double tau_3 = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                tau_3 += far_field_stress_[i][j] * frac.normal[i] * frac.dip_dir[j];
            }
        }
        state.K_III = tau_3 * std::sqrt(M_PI * frac.radius);
        
        // Energy release rate: G = (K_I² + K_II²)/E' + K_III²/(2μ)
        double E_prime = youngs_modulus_ / (1.0 - poisson_ratio_ * poisson_ratio_);
        double mu = youngs_modulus_ / (2.0 * (1.0 + poisson_ratio_));
        
        state.energy_release_rate = (state.K_I * state.K_I + state.K_II * state.K_II) / E_prime +
                                    state.K_III * state.K_III / (2.0 * mu);
        
        // Propagation criterion: K_I > K_IC
        state.will_propagate = (state.K_I > fracture_toughness_);
        
        states.push_back(state);
    }
    
    return states;
}

void FracturePropagationModel::propagateNetwork(FractureNetwork& network, double dt) {
    auto states = computeTipStates(network);
    
    for (const auto& state : states) {
        if (!state.will_propagate) continue;
        
        DiscreteFracture* frac = network.getFracture(state.fracture_id);
        if (!frac) continue;
        
        // Compute propagation velocity
        // Paris law: da/dt = C * (K_I / K_IC - 1)^m
        double C = 1e-3;  // Paris law coefficient
        double m = 2.0;   // Paris law exponent
        
        double K_ratio = state.K_I / fracture_toughness_;
        if (K_ratio > 1.0) {
            double velocity = C * std::pow(K_ratio - 1.0, m);
            double da = velocity * dt;
            
            // Increase fracture radius
            frac->radius += da;
            frac->length = 2.0 * frac->radius;
            frac->height = 2.0 * frac->radius;
        }
    }
    
    // Recompute intersections after propagation
    network.computeIntersections();
    network.buildConnectivityGraph();
}

void FracturePropagationModel::computeStressShadowing(FractureNetwork& network) const {
    // Stress shadowing: nearby fractures reduce stress intensity
    // Uses superposition of stress fields
    
    auto& fractures = network.getFractures();
    
    for (size_t i = 0; i < fractures.size(); ++i) {
        auto& frac_i = fractures[i];
        
        double stress_reduction = 0.0;
        
        for (size_t j = 0; j < fractures.size(); ++j) {
            if (i == j) continue;
            
            auto& frac_j = fractures[j];
            
            // Distance between fracture centers
            double dx = frac_i.center[0] - frac_j.center[0];
            double dy = frac_i.center[1] - frac_j.center[1];
            double dz = frac_i.center[2] - frac_j.center[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            // Stress shadowing effect decreases with distance
            // and depends on relative orientation
            if (dist < 5.0 * frac_j.radius) {
                // Dot product of normals (parallel fractures shadow more)
                double n_dot = frac_i.normal[0] * frac_j.normal[0] +
                               frac_i.normal[1] * frac_j.normal[1] +
                               frac_i.normal[2] * frac_j.normal[2];
                
                double shadow_factor = std::abs(n_dot) * frac_j.radius / dist;
                stress_reduction += shadow_factor;
            }
        }
        
        // Apply stress reduction to aperture (proxy for reduced opening)
        double reduction_factor = std::exp(-stress_reduction * 0.1);
        frac_i.current_aperture = frac_i.aperture * reduction_factor;
    }
}

std::vector<std::pair<int, int>> FracturePropagationModel::detectCoalescence(
    const FractureNetwork& network, double threshold_distance) const {
    
    std::vector<std::pair<int, int>> coalescence_pairs;
    
    const auto& fractures = network.getFractures();
    
    for (size_t i = 0; i < fractures.size(); ++i) {
        for (size_t j = i + 1; j < fractures.size(); ++j) {
            const auto& f1 = fractures[i];
            const auto& f2 = fractures[j];
            
            // Find closest approach between fracture tips
            // Simplified: check distance between perimeters
            
            auto v1 = f1.getVertices(8);
            auto v2 = f2.getVertices(8);
            
            double min_dist = std::numeric_limits<double>::infinity();
            
            for (const auto& p1 : v1) {
                for (const auto& p2 : v2) {
                    double dx = p1[0] - p2[0];
                    double dy = p1[1] - p2[1];
                    double dz = p1[2] - p2[2];
                    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                    min_dist = std::min(min_dist, dist);
                }
            }
            
            if (min_dist < threshold_distance) {
                coalescence_pairs.push_back({f1.id, f2.id});
            }
        }
    }
    
    return coalescence_pairs;
}

// ============================================================================
// FractureFlowSimulator Implementation
// ============================================================================

FractureFlowSimulator::FractureFlowSimulator()
    : network_(nullptr),
      viscosity_(0.001),
      compressibility_(4.5e-10) {}

void FractureFlowSimulator::setNetwork(FractureNetwork* network) {
    network_ = network;
}

void FractureFlowSimulator::setFluidProperties(double viscosity, double compressibility) {
    viscosity_ = viscosity;
    compressibility_ = compressibility;
}

void FractureFlowSimulator::setInletPressure(int fracture_id, double pressure) {
    inlet_pressures_[fracture_id] = pressure;
}

void FractureFlowSimulator::setOutletPressure(int fracture_id, double pressure) {
    outlet_pressures_[fracture_id] = pressure;
}

void FractureFlowSimulator::setInletFlowRate(int fracture_id, double rate) {
    inlet_rates_[fracture_id] = rate;
}

FractureFlowSimulator::FlowSolution FractureFlowSimulator::solveSteadyState() {
    FlowSolution solution;
    
    if (!network_ || network_->getNumFractures() == 0) {
        return solution;
    }
    
    const auto& fractures = network_->getFractures();
    const auto& intersections = network_->getIntersections();
    size_t n = fractures.size();
    
    solution.pressures.resize(n, 0.0);
    solution.flow_rates.resize(n, 0.0);
    
    // Build system matrix: A * p = b
    // Using simple node-based formulation
    // Each fracture is a node, intersections provide connectivity
    
    // Simple Gauss-Seidel iteration
    std::vector<double> pressure(n, 0.0);
    std::vector<double> source(n, 0.0);
    std::vector<std::vector<std::pair<size_t, double>>> connections(n);
    
    // Initialize with boundary conditions
    for (const auto& [frac_id, p] : inlet_pressures_) {
        for (size_t i = 0; i < n; ++i) {
            if (fractures[i].id == frac_id) {
                pressure[i] = p;
                break;
            }
        }
    }
    
    for (const auto& [frac_id, p] : outlet_pressures_) {
        for (size_t i = 0; i < n; ++i) {
            if (fractures[i].id == frac_id) {
                pressure[i] = p;
                break;
            }
        }
    }
    
    // Build connectivity from intersections
    auto getId = [&fractures](int frac_id) -> size_t {
        for (size_t i = 0; i < fractures.size(); ++i) {
            if (fractures[i].id == frac_id) return i;
        }
        return 0;
    };
    
    for (const auto& inter : intersections) {
        size_t i = getId(inter.fracture1_id);
        size_t j = getId(inter.fracture2_id);
        
        // Transmissibility between fractures
        double T = inter.transmissivity / viscosity_;
        
        connections[i].push_back({j, T});
        connections[j].push_back({i, T});
    }
    
    // Iterative solution
    int max_iter = 1000;
    double tol = 1e-6;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_change = 0.0;
        
        for (size_t i = 0; i < n; ++i) {
            // Skip fixed pressure nodes
            bool is_fixed = false;
            for (const auto& [frac_id, p] : inlet_pressures_) {
                if (fractures[i].id == frac_id) { is_fixed = true; break; }
            }
            for (const auto& [frac_id, p] : outlet_pressures_) {
                if (fractures[i].id == frac_id) { is_fixed = true; break; }
            }
            
            if (is_fixed) continue;
            
            double sum_T = 0.0;
            double sum_Tp = 0.0;
            
            for (const auto& [j, T] : connections[i]) {
                sum_T += T;
                sum_Tp += T * pressure[j];
            }
            
            if (sum_T > 1e-20) {
                double new_p = (sum_Tp + source[i]) / sum_T;
                max_change = std::max(max_change, std::abs(new_p - pressure[i]));
                pressure[i] = new_p;
            }
        }
        
        if (max_change < tol) break;
    }
    
    // Compute flow rates
    solution.pressures = pressure;
    solution.total_inflow = 0.0;
    solution.total_outflow = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        double net_flow = 0.0;
        
        for (const auto& [j, T] : connections[i]) {
            net_flow += T * (pressure[j] - pressure[i]);
        }
        
        solution.flow_rates[i] = net_flow;
        
        // Check if inlet/outlet
        for (const auto& [frac_id, p] : inlet_pressures_) {
            if (fractures[i].id == frac_id) {
                solution.total_inflow += std::abs(net_flow);
            }
        }
        for (const auto& [frac_id, p] : outlet_pressures_) {
            if (fractures[i].id == frac_id) {
                solution.total_outflow += std::abs(net_flow);
            }
        }
    }
    
    // Effective transmissivity
    double dp = 0.0;
    for (const auto& [id1, p1] : inlet_pressures_) {
        for (const auto& [id2, p2] : outlet_pressures_) {
            dp = std::max(dp, std::abs(p1 - p2));
        }
    }
    
    if (dp > 1e-10) {
        solution.effective_transmissivity = solution.total_inflow / dp;
    }
    
    return solution;
}

std::array<double, 3> FractureFlowSimulator::computeEquivalentPermeability() {
    std::array<double, 3> perm = {0.0, 0.0, 0.0};
    
    if (!network_) return perm;
    
    const auto& fractures = network_->getFractures();
    if (fractures.empty()) return perm;
    
    // Get domain bounds
    std::array<double, 3> L = {0, 0, 0};
    std::array<double, 3> domain_min = {1e30, 1e30, 1e30};
    std::array<double, 3> domain_max = {-1e30, -1e30, -1e30};
    
    for (const auto& frac : fractures) {
        for (int i = 0; i < 3; ++i) {
            domain_min[i] = std::min(domain_min[i], frac.center[i] - frac.radius);
            domain_max[i] = std::max(domain_max[i], frac.center[i] + frac.radius);
        }
    }
    
    for (int i = 0; i < 3; ++i) {
        L[i] = domain_max[i] - domain_min[i];
    }
    
    // Compute permeability in each direction by flow simulation
    double dp = 1e6;  // 1 MPa pressure difference
    
    for (int dir = 0; dir < 3; ++dir) {
        inlet_pressures_.clear();
        outlet_pressures_.clear();
        
        // Find fractures at boundaries
        for (const auto& frac : fractures) {
            if (frac.center[dir] - frac.radius < domain_min[dir] + L[dir] * 0.1) {
                inlet_pressures_[frac.id] = dp;
            }
            if (frac.center[dir] + frac.radius > domain_max[dir] - L[dir] * 0.1) {
                outlet_pressures_[frac.id] = 0.0;
            }
        }
        
        if (inlet_pressures_.empty() || outlet_pressures_.empty()) {
            continue;
        }
        
        auto solution = solveSteadyState();
        
        // k = Q * μ * L / (A * ΔP)
        double A = 1.0;
        for (int i = 0; i < 3; ++i) {
            if (i != dir) A *= L[i];
        }
        
        perm[dir] = solution.total_inflow * viscosity_ * L[dir] / (A * dp);
    }
    
    return perm;
}

// ============================================================================
// DFNStatistics Implementation
// ============================================================================

std::vector<FisherParameters> DFNStatistics::clusterOrientations(
    const std::vector<DiscreteFracture>& fractures, int n_clusters) {
    
    std::vector<FisherParameters> clusters;
    
    if (fractures.empty() || n_clusters <= 0) {
        return clusters;
    }
    
    // Simple k-means clustering on orientation vectors
    
    // Convert orientations to unit vectors (pole to fracture plane)
    std::vector<std::array<double, 3>> poles;
    for (const auto& frac : fractures) {
        poles.push_back(frac.normal);
    }
    
    // Initialize cluster centers randomly
    std::vector<std::array<double, 3>> centers(n_clusters);
    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> dist(0, poles.size() - 1);
    
    for (int i = 0; i < n_clusters; ++i) {
        centers[i] = poles[dist(rng)];
    }
    
    // K-means iterations
    std::vector<int> assignments(poles.size());
    
    for (int iter = 0; iter < 100; ++iter) {
        // Assign to nearest cluster
        for (size_t i = 0; i < poles.size(); ++i) {
            double max_dot = -2.0;
            int best_cluster = 0;
            
            for (int c = 0; c < n_clusters; ++c) {
                double dot = std::abs(poles[i][0] * centers[c][0] +
                                      poles[i][1] * centers[c][1] +
                                      poles[i][2] * centers[c][2]);
                if (dot > max_dot) {
                    max_dot = dot;
                    best_cluster = c;
                }
            }
            
            assignments[i] = best_cluster;
        }
        
        // Update centers
        for (int c = 0; c < n_clusters; ++c) {
            std::array<double, 3> sum = {0, 0, 0};
            int count = 0;
            
            for (size_t i = 0; i < poles.size(); ++i) {
                if (assignments[i] == c) {
                    // Ensure consistent direction
                    double dot = poles[i][0] * centers[c][0] +
                                 poles[i][1] * centers[c][1] +
                                 poles[i][2] * centers[c][2];
                    double sign = (dot >= 0) ? 1.0 : -1.0;
                    
                    sum[0] += sign * poles[i][0];
                    sum[1] += sign * poles[i][1];
                    sum[2] += sign * poles[i][2];
                    count++;
                }
            }
            
            if (count > 0) {
                double norm = std::sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
                if (norm > 1e-10) {
                    centers[c][0] = sum[0] / norm;
                    centers[c][1] = sum[1] / norm;
                    centers[c][2] = sum[2] / norm;
                }
            }
        }
    }
    
    // Convert centers to Fisher parameters
    for (int c = 0; c < n_clusters; ++c) {
        FisherParameters params;
        
        // Convert normal to strike/dip
        double dip = std::acos(std::abs(centers[c][2])) * 180.0 / M_PI;
        double strike = std::atan2(centers[c][0], centers[c][1]) * 180.0 / M_PI;
        if (strike < 0) strike += 360.0;
        
        params.mean_strike = strike;
        params.mean_dip = dip;
        
        // Estimate kappa from cluster dispersion
        std::vector<std::pair<double, double>> cluster_orientations;
        for (size_t i = 0; i < poles.size(); ++i) {
            if (assignments[i] == c) {
                double frac_dip = std::acos(std::abs(poles[i][2])) * 180.0 / M_PI;
                double frac_strike = std::atan2(poles[i][0], poles[i][1]) * 180.0 / M_PI;
                if (frac_strike < 0) frac_strike += 360.0;
                cluster_orientations.push_back({frac_strike, frac_dip});
            }
        }
        
        params.kappa = computeFisherKappa(cluster_orientations);
        clusters.push_back(params);
    }
    
    return clusters;
}

double DFNStatistics::computeFisherKappa(
    const std::vector<std::pair<double, double>>& orientations) {
    
    if (orientations.size() < 2) {
        return 100.0;  // High concentration for single point
    }
    
    // Compute resultant length R
    double Rx = 0.0, Ry = 0.0, Rz = 0.0;
    
    for (const auto& [strike, dip] : orientations) {
        double strike_rad = strike * M_PI / 180.0;
        double dip_rad = dip * M_PI / 180.0;
        
        // Convert to unit vector (pole)
        double x = std::sin(dip_rad) * std::sin(strike_rad);
        double y = std::sin(dip_rad) * std::cos(strike_rad);
        double z = std::cos(dip_rad);
        
        Rx += x;
        Ry += y;
        Rz += z;
    }
    
    double R = std::sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
    double R_bar = R / orientations.size();
    
    // Estimate kappa (Fisher concentration parameter)
    // For large kappa: kappa ≈ (n-1) / (n - R)
    // size_t n = orientations.size();  // Not used in simplified estimate
    
    if (R_bar > 0.999) {
        return 1000.0;  // Very high concentration
    }
    
    // Best estimate: κ ≈ R̄(d - R̄²) / (1 - R̄²) where d is dimension (3 for sphere)
    double kappa = R_bar * (3.0 - R_bar * R_bar) / (1.0 - R_bar * R_bar);
    
    return std::max(0.1, kappa);
}

std::pair<SizeDistribution, std::vector<double>> DFNStatistics::fitSizeDistribution(
    const std::vector<double>& sizes) {
    
    if (sizes.empty()) {
        return {SizeDistribution::CONSTANT, {1.0}};
    }
    
    // Compute statistics
    double mean = 0.0, var = 0.0;
    double min_size = sizes[0], max_size = sizes[0];
    
    for (double s : sizes) {
        mean += s;
        min_size = std::min(min_size, s);
        max_size = std::max(max_size, s);
    }
    mean /= sizes.size();
    
    for (double s : sizes) {
        var += (s - mean) * (s - mean);
    }
    var /= sizes.size();
    double std_dev = std::sqrt(var);
    
    // Log statistics
    double log_mean = 0.0, log_var = 0.0;
    for (double s : sizes) {
        log_mean += std::log(s);
    }
    log_mean /= sizes.size();
    
    for (double s : sizes) {
        double diff = std::log(s) - log_mean;
        log_var += diff * diff;
    }
    log_var /= sizes.size();
    
    // Compare coefficient of variation to choose distribution
    double cv = std_dev / mean;
    double log_cv = std::sqrt(log_var);
    
    if (cv < 0.1) {
        // Low variation -> constant
        return {SizeDistribution::CONSTANT, {mean}};
    } else if (log_cv < cv * 0.8) {
        // Log-normal fits better
        return {SizeDistribution::LOGNORMAL, {std::exp(log_mean), std::sqrt(log_var)}};
    } else if (max_size / min_size > 100) {
        // Large range suggests power-law
        // Estimate exponent using MLE
        double alpha = 1.0 + sizes.size() / 
            std::accumulate(sizes.begin(), sizes.end(), 0.0,
                [min_size](double sum, double s) { return sum + std::log(s / min_size); });
        
        return {SizeDistribution::POWER_LAW, {min_size, max_size, alpha}};
    } else {
        // Default to log-normal
        return {SizeDistribution::LOGNORMAL, {std::exp(log_mean), std::sqrt(log_var)}};
    }
}

double DFNStatistics::computeFractalDimension(
    const std::vector<DiscreteFracture>& fractures,
    const std::array<double, 3>& domain_size) {
    
    if (fractures.empty()) return 0.0;
    
    // Box-counting method
    std::vector<double> box_sizes;
    std::vector<int> box_counts;
    
    double min_dim = std::min({domain_size[0], domain_size[1], domain_size[2]});
    
    for (double box_size = min_dim / 2.0; box_size >= min_dim / 64.0; box_size /= 2.0) {
        int nx = static_cast<int>(domain_size[0] / box_size) + 1;
        int ny = static_cast<int>(domain_size[1] / box_size) + 1;
        int nz = static_cast<int>(domain_size[2] / box_size) + 1;
        
        std::set<std::tuple<int, int, int>> occupied_boxes;
        
        for (const auto& frac : fractures) {
            int ix = static_cast<int>(frac.center[0] / box_size);
            int iy = static_cast<int>(frac.center[1] / box_size);
            int iz = static_cast<int>(frac.center[2] / box_size);
            
            // Also count neighboring boxes within fracture radius
            int r_boxes = static_cast<int>(frac.radius / box_size) + 1;
            
            for (int di = -r_boxes; di <= r_boxes; ++di) {
                for (int dj = -r_boxes; dj <= r_boxes; ++dj) {
                    for (int dk = -r_boxes; dk <= r_boxes; ++dk) {
                        int bx = ix + di;
                        int by = iy + dj;
                        int bz = iz + dk;
                        
                        if (bx >= 0 && bx < nx && by >= 0 && by < ny && bz >= 0 && bz < nz) {
                            occupied_boxes.insert({bx, by, bz});
                        }
                    }
                }
            }
        }
        
        box_sizes.push_back(box_size);
        box_counts.push_back(occupied_boxes.size());
    }
    
    // Linear regression on log-log plot: log(N) = D * log(1/ε) + C
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    int n = box_sizes.size();
    
    for (int i = 0; i < n; ++i) {
        double x = std::log(1.0 / box_sizes[i]);
        double y = std::log(box_counts[i]);
        
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_xx += x * x;
    }
    
    double D = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    
    return std::max(0.0, std::min(3.0, D));
}

std::vector<double> DFNStatistics::computeRipleyK(
    const std::vector<DiscreteFracture>& fractures,
    const std::vector<double>& distances) {
    
    std::vector<double> K_values;
    
    if (fractures.size() < 2 || distances.empty()) {
        return K_values;
    }
    
    // Estimate domain volume from fracture extent
    double x_min = 1e30, y_min = 1e30, z_min = 1e30;
    double x_max = -1e30, y_max = -1e30, z_max = -1e30;
    
    for (const auto& frac : fractures) {
        x_min = std::min(x_min, frac.center[0]);
        y_min = std::min(y_min, frac.center[1]);
        z_min = std::min(z_min, frac.center[2]);
        x_max = std::max(x_max, frac.center[0]);
        y_max = std::max(y_max, frac.center[1]);
        z_max = std::max(z_max, frac.center[2]);
    }
    
    double volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    // double lambda = fractures.size() / volume;  // Intensity - not used in current implementation
    
    for (double r : distances) {
        int count = 0;
        
        for (size_t i = 0; i < fractures.size(); ++i) {
            for (size_t j = i + 1; j < fractures.size(); ++j) {
                double dx = fractures[i].center[0] - fractures[j].center[0];
                double dy = fractures[i].center[1] - fractures[j].center[1];
                double dz = fractures[i].center[2] - fractures[j].center[2];
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist <= r) {
                    count += 2;  // Count both directions
                }
            }
        }
        
        // K(r) = V / n² * Σ I(d_ij ≤ r)
        double K = volume * count / (fractures.size() * fractures.size());
        K_values.push_back(K);
    }
    
    return K_values;
}

FractureIntensity DFNStatistics::computeIntensityFromScanline(
    const std::vector<DiscreteFracture>& fractures,
    const std::array<double, 3>& line_start,
    const std::array<double, 3>& line_end) {
    
    FractureIntensity intensity;
    
    // Line direction and length
    double dx = line_end[0] - line_start[0];
    double dy = line_end[1] - line_start[1];
    double dz = line_end[2] - line_start[2];
    double line_length = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    if (line_length < 1e-10) return intensity;
    
    dx /= line_length;
    dy /= line_length;
    dz /= line_length;
    
    int intersection_count = 0;
    
    for (const auto& frac : fractures) {
        // Check if scanline intersects fracture disc
        // Line: P(t) = start + t * dir
        // Plane: n · (P - center) = 0
        
        double denom = frac.normal[0] * dx + frac.normal[1] * dy + frac.normal[2] * dz;
        
        if (std::abs(denom) < 1e-10) {
            continue;  // Line parallel to fracture
        }
        
        double t = ((frac.center[0] - line_start[0]) * frac.normal[0] +
                    (frac.center[1] - line_start[1]) * frac.normal[1] +
                    (frac.center[2] - line_start[2]) * frac.normal[2]) / denom;
        
        if (t < 0 || t > line_length) {
            continue;  // Intersection outside line segment
        }
        
        // Intersection point
        double px = line_start[0] + t * dx;
        double py = line_start[1] + t * dy;
        double pz = line_start[2] + t * dz;
        
        // Check if within fracture radius
        double dist_x = px - frac.center[0];
        double dist_y = py - frac.center[1];
        double dist_z = pz - frac.center[2];
        double dist_sq = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
        
        if (dist_sq <= frac.radius * frac.radius) {
            intersection_count++;
        }
    }
    
    intensity.P10 = intersection_count / line_length;
    
    return intensity;
}

FractureIntensity DFNStatistics::computeIntensityFromWindow(
    const std::vector<DiscreteFracture>& fractures,
    const std::array<double, 3>& window_center,
    const std::array<double, 3>& window_size) {
    
    FractureIntensity intensity;
    
    double window_area = window_size[0] * window_size[1];
    double window_volume = window_size[0] * window_size[1] * window_size[2];
    
    if (window_area < 1e-10) return intensity;
    
    int count = 0;
    double total_trace_length = 0.0;
    double total_area = 0.0;
    
    for (const auto& frac : fractures) {
        // Check if fracture intersects window
        bool in_window = true;
        for (int i = 0; i < 3; ++i) {
            if (frac.center[i] + frac.radius < window_center[i] - window_size[i]/2 ||
                frac.center[i] - frac.radius > window_center[i] + window_size[i]/2) {
                in_window = false;
                break;
            }
        }
        
        if (!in_window) continue;
        
        count++;
        
        // Approximate trace length (intersection with horizontal plane at window center)
        double z_dist = std::abs(window_center[2] - frac.center[2]);
        if (z_dist < frac.radius * std::abs(frac.normal[2])) {
            double trace = 2.0 * std::sqrt(std::max(0.0, 
                frac.radius * frac.radius - z_dist * z_dist / 
                (frac.normal[2] * frac.normal[2] + 1e-10)));
            total_trace_length += trace;
        }
        
        total_area += frac.computeArea();
    }
    
    intensity.P20 = count / window_area;
    intensity.P21 = total_trace_length / window_area;
    intensity.P30 = count / window_volume;
    intensity.P32 = total_area / window_volume;
    
    return intensity;
}

// ============================================================================
// Configuration Parser
// ============================================================================

void parseDFNConfig(const std::string& filename, FractureNetwork& network) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open DFN config file: " + filename);
    }
    
    std::string line;
    FractureSet* current_set = nullptr;
    int set_id = 0;
    
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#' || line[0] == ';') {
            continue;
        }
        
        std::istringstream ss(line);
        std::string key;
        ss >> key;
        
        if (key == "[DOMAIN]") {
            double xmin, ymin, zmin, xmax, ymax, zmax;
            std::getline(file, line);
            std::istringstream dim_ss(line);
            dim_ss >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax;
            network.setDomain({xmin, ymin, zmin}, {xmax, ymax, zmax});
        }
        else if (key == "[FRACTURE_SET]") {
            std::string name;
            ss >> name;
            network.addFractureSet(FractureSet(set_id++, name));
            current_set = network.getFractureSet(set_id - 1);
        }
        else if (key == "orientation" && current_set) {
            std::string type;
            ss >> type;
            if (type == "FISHER") {
                double strike, dip, kappa;
                ss >> strike >> dip >> kappa;
                current_set->setOrientationDistribution(OrientationDistribution::FISHER);
                current_set->setFisherParameters(FisherParameters(strike, dip, kappa));
            } else if (type == "UNIFORM") {
                current_set->setOrientationDistribution(OrientationDistribution::UNIFORM);
            }
        }
        else if (key == "size" && current_set) {
            std::string type;
            ss >> type;
            if (type == "LOGNORMAL") {
                double mean, std;
                ss >> mean >> std;
                current_set->setSizeDistribution(SizeDistribution::LOGNORMAL);
                current_set->setSizeParameters(mean, std);
            } else if (type == "POWER_LAW") {
                double xmin, xmax, alpha;
                ss >> xmin >> xmax >> alpha;
                current_set->setSizeDistribution(SizeDistribution::POWER_LAW);
                current_set->setSizeParameters(xmin, xmax, alpha);
            }
        }
        else if (key == "intensity" && current_set) {
            double p32;
            ss >> p32;
            current_set->setIntensity(p32);
        }
        else if (key == "aperture" && current_set) {
            double ap;
            ss >> ap;
            current_set->setApertureParameters(ap, 0.3);
        }
        else if (key == "seed") {
            unsigned int seed;
            ss >> seed;
            network.setRandomSeed(seed);
        }
    }
    
    file.close();
}

} // namespace FSRM
