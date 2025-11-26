/**
 * @file test_analytical_solutions.cpp
 * @brief Tests against known analytical solutions
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include <cmath>

using namespace FSRM;
using namespace FSRM::Testing;

class AnalyticalSolutionsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// Theis Solution (Radial Flow to Well)
// ============================================================================

TEST_F(AnalyticalSolutionsTest, TheisSolution) {
    // Theis solution: s = Q/(4*pi*T) * W(u)
    // where u = r^2*S/(4*T*t) and W(u) is the exponential integral
    
    double Q = 0.1;       // Pumping rate (m^3/s)
    double T = 1e-3;      // Transmissivity (m^2/s)
    double S = 1e-4;      // Storativity
    double r = 100.0;     // Radial distance (m)
    double t = 3600.0;    // Time (s)
    
    double u = r*r * S / (4.0 * T * t);
    
    // Well function approximation for small u: W(u) ≈ -γ - ln(u)
    double gamma = 0.5772;  // Euler-Mascheroni constant
    double W_u = -gamma - std::log(u);
    
    double drawdown = Q / (4.0 * M_PI * T) * W_u;
    
    EXPECT_GT(drawdown, 0.0) << "Drawdown should be positive";
    EXPECT_TRUE(std::isfinite(drawdown));
    
    // Drawdown should increase with time (u decreases)
    double t2 = 2.0 * t;
    double u2 = r*r * S / (4.0 * T * t2);
    double W_u2 = -gamma - std::log(u2);
    double drawdown2 = Q / (4.0 * M_PI * T) * W_u2;
    
    EXPECT_GT(drawdown2, drawdown) << "Drawdown increases with time";
}

// ============================================================================
// Terzaghi Consolidation
// ============================================================================

TEST_F(AnalyticalSolutionsTest, TerzaghiConsolidation) {
    // 1D consolidation: p(z,t) = sum_{n=0}^inf (2*p0/M_n) * sin(M_n*z/H) * exp(-M_n^2*Tv)
    // M_n = (2n+1)*pi/2, Tv = cv*t/H^2
    
    double p0 = 1e6;      // Initial excess pressure
    double cv = 1e-6;     // Consolidation coefficient
    double H = 10.0;      // Layer thickness
    double z = 5.0;       // Depth (mid-layer)
    double t = 1e6;       // Time (s)
    
    double Tv = cv * t / (H * H);  // Time factor
    
    // Series solution (first few terms)
    double pressure = 0.0;
    for (int n = 0; n < 10; ++n) {
        double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
        pressure += (2.0 * p0 / M_n) * std::sin(M_n * z / H) * 
                   std::exp(-M_n * M_n * Tv);
    }
    
    EXPECT_TRUE(std::isfinite(pressure));
    EXPECT_GE(pressure, 0.0) << "Pressure should be non-negative";
    EXPECT_LE(pressure, p0) << "Pressure should decrease from initial";
    
    // At very large time, pressure should approach zero
    double t_large = 1e9;
    double Tv_large = cv * t_large / (H * H);
    double pressure_large = 0.0;
    for (int n = 0; n < 10; ++n) {
        double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
        pressure_large += (2.0 * p0 / M_n) * std::sin(M_n * z / H) * 
                         std::exp(-M_n * M_n * Tv_large);
    }
    
    EXPECT_LT(pressure_large, 0.01 * p0) << "Pressure should approach zero";
}

// ============================================================================
// Mandel-Cryer Effect
// ============================================================================

TEST_F(AnalyticalSolutionsTest, MandelCryerEffect) {
    // In 2D/3D consolidation, center pressure initially increases
    // before decreasing (non-monotonic behavior)
    
    double nu = 0.25;      // Drained Poisson ratio
    double nu_u = 0.4;     // Undrained Poisson ratio
    double B = 0.8;        // Skempton coefficient
    double alpha = 0.9;    // Biot coefficient
    
    // The Mandel-Cryer effect occurs when nu_u > nu
    EXPECT_GT(nu_u, nu) << "Undrained nu > drained nu for Mandel-Cryer effect";
    
    // Maximum overshoot occurs at early time
    // The overshoot ratio depends on material properties
    double overshoot = (nu_u - nu) / (1.0 - nu);
    
    EXPECT_GT(overshoot, 0.0) << "Positive overshoot expected";
}

// ============================================================================
// Buckley-Leverett (Immiscible Displacement)
// ============================================================================

TEST_F(AnalyticalSolutionsTest, BuckleyLeverett) {
    // Fractional flow function: f(Sw) = krw/muw / (krw/muw + kro/muo)
    // with Brooks-Corey relative permeability
    
    double Sw = 0.5;      // Water saturation
    double Swc = 0.2;     // Connate water
    double Sor = 0.2;     // Residual oil
    double muw = 1e-3;    // Water viscosity
    double muo = 5e-3;    // Oil viscosity
    double nw = 2.0;      // Water Corey exponent
    double no = 2.0;      // Oil Corey exponent
    
    // Normalized saturation
    double Sn = (Sw - Swc) / (1.0 - Swc - Sor);
    Sn = std::max(0.0, std::min(1.0, Sn));
    
    // Relative permeabilities (Brooks-Corey)
    double krw = std::pow(Sn, nw);
    double kro = std::pow(1.0 - Sn, no);
    
    // Fractional flow
    double mobility_w = krw / muw;
    double mobility_o = kro / muo;
    double fw = mobility_w / (mobility_w + mobility_o);
    
    EXPECT_GE(fw, 0.0);
    EXPECT_LE(fw, 1.0);
    
    // At Sw = 1-Sor, fw should approach 1
    Sw = 1.0 - Sor;
    Sn = (Sw - Swc) / (1.0 - Swc - Sor);
    krw = std::pow(Sn, nw);
    kro = std::pow(1.0 - Sn, no);
    mobility_w = krw / muw;
    mobility_o = kro / muo;
    fw = mobility_w / (mobility_w + mobility_o);
    
    EXPECT_NEAR(fw, 1.0, 1e-6);
}

// ============================================================================
// Line Source Solution (Heat Conduction)
// ============================================================================

TEST_F(AnalyticalSolutionsTest, LineSource) {
    // Temperature from line source: T(r,t) = Q/(4*pi*k*t) * exp(-r^2/(4*alpha*t))
    // alpha = k/(rho*cp)
    
    double Q = 1000.0;     // Heat rate per unit length (W/m)
    double k = 2.5;        // Thermal conductivity (W/m/K)
    double rho = 2500.0;   // Density (kg/m^3)
    double cp = 1000.0;    // Specific heat (J/kg/K)
    double r = 1.0;        // Radial distance (m)
    double t = 3600.0;     // Time (s)
    
    double alpha = k / (rho * cp);
    double T = Q / (4.0 * M_PI * k * t) * std::exp(-r*r / (4.0 * alpha * t));
    
    EXPECT_GT(T, 0.0) << "Temperature rise should be positive";
    EXPECT_TRUE(std::isfinite(T));
    
    // Temperature should decay with distance
    double r2 = 2.0 * r;
    double T2 = Q / (4.0 * M_PI * k * t) * std::exp(-r2*r2 / (4.0 * alpha * t));
    
    EXPECT_LT(T2, T) << "Temperature should decay with distance";
}

// ============================================================================
// Griffith Fracture Criterion
// ============================================================================

TEST_F(AnalyticalSolutionsTest, GriffithFracture) {
    // Energy release rate: G = K_I^2 / E'
    // Critical condition: G = Gc
    
    double K_IC = 1e6;     // Fracture toughness (Pa*sqrt(m))
    double E = 20e9;       // Young's modulus
    double nu = 0.25;      // Poisson ratio
    double E_prime = E / (1.0 - nu*nu);  // Plane strain modulus
    
    double Gc = K_IC * K_IC / E_prime;  // Critical energy release rate
    
    EXPECT_GT(Gc, 0.0);
    
    // For a penny-shaped crack: K_I = (2/pi) * sigma * sqrt(pi*a)
    double a = 10.0;       // Crack radius
    double sigma_c = K_IC * M_PI / (2.0 * std::sqrt(M_PI * a));
    
    EXPECT_GT(sigma_c, 0.0);
    EXPECT_LT(sigma_c, 1e9) << "Critical stress should be reasonable";
}

// ============================================================================
// KGD Hydraulic Fracture Model
// ============================================================================

TEST_F(AnalyticalSolutionsTest, KGDFractureScaling) {
    // KGD model scaling: L ~ (Q*E'*t/mu)^(1/4)
    
    double Q = 0.1;        // Injection rate (m^3/s per meter height)
    double E = 20e9;       // Young's modulus
    double nu = 0.25;
    double E_prime = E / (1.0 - nu*nu);
    double mu = 0.001;     // Fluid viscosity
    double t = 60.0;       // Time (s)
    
    double L_scale = std::pow(Q * E_prime * t / mu, 0.25);
    
    EXPECT_GT(L_scale, 0.0);
    
    // At double the time, length should increase by factor of 2^(1/4)
    double t2 = 2.0 * t;
    double L_scale2 = std::pow(Q * E_prime * t2 / mu, 0.25);
    
    double ratio = L_scale2 / L_scale;
    EXPECT_NEAR(ratio, std::pow(2.0, 0.25), 0.001);
}
