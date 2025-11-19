/*
 * Linear Elastic Fracture Mechanics (LEFM) Implementation
 * 
 * This file provides detailed LEFM calculations for fracture propagation
 */

#include "PhysicsKernel.hpp"
#include <cmath>
#include <algorithm>

namespace ResSim {
namespace LEFM {

// Compute Mode I stress intensity factor for pressurized crack
// Using Sneddon's solution for internal pressure
double computeStressIntensityFactor(double pressure, double crack_length, 
                                    double min_stress, double geometry_factor = 1.0) {
    // K_I = Y * (P - σ) * sqrt(π * a)
    // where:
    //   Y = geometry factor (1.0 for infinite body, varies for finite geometry)
    //   P = internal pressure in crack
    //   σ = confining stress perpendicular to crack
    //   a = crack half-length
    
    double net_pressure = pressure - min_stress;
    double a = crack_length / 2.0;  // half-length
    double K_I = geometry_factor * net_pressure * std::sqrt(M_PI * a);
    
    return K_I;
}

// Check if fracture should initiate
bool checkInitiation(double K_I, double K_Ic) {
    return K_I >= K_Ic;
}

// Compute propagation velocity based on K_I
// Using Paris-like law: v = C * ((K_I - K_Ic) / K_Ic)^m
double computePropagationVelocity(double K_I, double K_Ic, 
                                  double C_coeff = 1.0, double m_exp = 1.0) {
    if (K_I <= K_Ic) {
        return 0.0;  // No propagation if below toughness
    }
    
    double normalized_excess = (K_I - K_Ic) / K_Ic;
    double velocity = C_coeff * std::pow(normalized_excess, m_exp);
    
    // Limit maximum velocity (physical constraint)
    double v_max = 1000.0;  // 1000 m/s (~ shear wave speed)
    return std::min(velocity, v_max);
}

// Compute fracture width profile using Sneddon's solution
// For pressurized crack in infinite elastic medium
void computeWidthProfile(double crack_length, double pressure, double min_stress,
                        double youngs_modulus, double poisson_ratio,
                        int num_points, double* x_coords, double* widths) {
    // w(x) = 4*(1-ν²)/E * sqrt(a² - x²) * (P - σ)
    
    double a = crack_length / 2.0;
    double net_pressure = pressure - min_stress;
    double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    double coeff = 4.0 * net_pressure / E_prime;
    
    for (int i = 0; i < num_points; ++i) {
        double x = x_coords[i];
        
        if (std::abs(x) >= a) {
            widths[i] = 0.0;  // Outside crack
        } else {
            widths[i] = coeff * std::sqrt(a*a - x*x);
        }
    }
}

// Compute maximum width at crack center
double computeMaxWidth(double crack_length, double pressure, double min_stress,
                      double youngs_modulus, double poisson_ratio) {
    double a = crack_length / 2.0;
    double net_pressure = pressure - min_stress;
    double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    
    // w_max = 4*(1-ν²)/E * a * (P - σ)
    double w_max = 4.0 * net_pressure * a / E_prime;
    
    return w_max;
}

// Energy release rate (G) from stress intensity factor
double computeEnergyReleaseRate(double K_I, double youngs_modulus, 
                               double poisson_ratio) {
    // G = K_I² * (1 - ν²) / E
    double G = K_I * K_I * (1.0 - poisson_ratio * poisson_ratio) / youngs_modulus;
    return G;
}

// Check Irwin criterion: G >= G_c
bool checkIrwinCriterion(double G, double G_c) {
    return G >= G_c;
}

// Compute equivalent K_Ic from G_c
double computeKIcFromGc(double G_c, double youngs_modulus, double poisson_ratio) {
    // K_Ic = sqrt(E * G_c / (1 - ν²))
    double K_Ic = std::sqrt(youngs_modulus * G_c / (1.0 - poisson_ratio * poisson_ratio));
    return K_Ic;
}

// PKN model - plane strain width
double computePKNWidth(double net_pressure, double height, 
                      double youngs_modulus, double poisson_ratio) {
    // w = 4 * (P - σ) * H / E'
    double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    return 4.0 * net_pressure * height / E_prime;
}

// KGD model - plane strain width
double computeKGDWidth(double net_pressure, double length,
                      double youngs_modulus, double poisson_ratio) {
    // w = 4 * (P - σ) * L / E'
    double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    return 4.0 * net_pressure * length / E_prime;
}

// Radial/penny-shaped fracture
double computeRadialWidth(double net_pressure, double radius,
                         double youngs_modulus, double poisson_ratio) {
    // w_max = 2 * R * K_I / (E' * sqrt(π))
    double E_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    double K_I = net_pressure * std::sqrt(M_PI * radius);
    return 2.0 * radius * K_I / (E_prime * std::sqrt(M_PI));
}

// Fluid flow in fracture - cubic law
double computeFlowRate(double width, double viscosity, double pressure_gradient,
                      double fracture_permeability_factor = 1.0/12.0) {
    // q = -w³/(12*μ) * dP/dx  (Poiseuille flow between parallel plates)
    double permeability = fracture_permeability_factor * width * width;
    double q = -permeability * width / viscosity * pressure_gradient;
    return q;
}

// Carter leak-off
double computeLeakoff(double time, double leakoff_coeff) {
    // v_L = C_L / sqrt(t)
    if (time <= 0.0) return 0.0;
    return leakoff_coeff / std::sqrt(time);
}

// Tip velocity vs. K_I relationship (Adachi & Detournay, 2007)
// For toughness-dominated regime
double computeTipVelocity_Toughness(double K_I, double K_Ic, double Kprime) {
    // v = Kprime * (K_I/K_Ic - 1)
    // where Kprime is a material parameter
    if (K_I <= K_Ic) return 0.0;
    return Kprime * (K_I / K_Ic - 1.0);
}

// Tip velocity for viscosity-dominated regime
double computeTipVelocity_Viscosity(double net_pressure, double viscosity,
                                   double length, double Eprime) {
    // v ~ (P - σ)³ * L / (μ * E')
    return std::pow(net_pressure, 3.0) * length / (viscosity * Eprime);
}

// T-stress (secondary stress parameter) for constraint
double computeTStress(double applied_stress, double crack_length, double a_over_W) {
    // T-stress affects crack path stability
    // T = σ * f(a/W)
    // Simplified: T ≈ -σ for through crack
    return -applied_stress;
}

// J-integral (energy release rate in nonlinear elasticity)
double computeJIntegral(double K_I, double K_II, double youngs_modulus, 
                       double poisson_ratio) {
    // J = (K_I² + K_II²) * (1 - ν²) / E
    // For Mode I only: J = G
    double J = (K_I*K_I + K_II*K_II) * (1.0 - poisson_ratio*poisson_ratio) / youngs_modulus;
    return J;
}

// Determine propagation direction (maximum hoop stress criterion)
double computePropagationAngle(double K_I, double K_II) {
    // θ_0 = 2*atan((K_I - sqrt(K_I² + 8*K_II²))/(4*K_II))
    // For Mode I only (K_II = 0): θ_0 = 0 (straight ahead)
    
    if (std::abs(K_II) < 1e-6) {
        return 0.0;  // Pure Mode I
    }
    
    double term = (K_I - std::sqrt(K_I*K_I + 8.0*K_II*K_II)) / (4.0 * K_II);
    return 2.0 * std::atan(term);
}

// Geometry factor Y for various configurations
namespace GeometryFactors {
    // Through crack in infinite plate
    double infinite_plate() { return 1.0; }
    
    // Edge crack
    double edge_crack(double a, double W) {
        double ratio = a / W;
        // Approximate formula
        return 1.12 - 0.23*ratio + 10.6*ratio*ratio - 21.7*ratio*ratio*ratio;
    }
    
    // Center crack in finite width plate
    double center_crack_finite(double a, double W) {
        double ratio = a / W;
        return std::sqrt(1.0 / std::cos(M_PI * ratio / 2.0));
    }
    
    // Elliptical crack
    double elliptical(double a, double c) {
        // a = semi-major axis, c = semi-minor axis
        double phi = 0.0;  // Angle around ellipse
        return std::sqrt(std::cos(phi)*std::cos(phi) + 
                        (a/c)*(a/c)*std::sin(phi)*std::sin(phi));
    }
    
    // Penny-shaped crack
    double penny_shaped() { return 2.0 / std::sqrt(M_PI); }
}

} // namespace LEFM
} // namespace ResSim
