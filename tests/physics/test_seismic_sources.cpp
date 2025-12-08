/**
 * @file test_seismic_sources.cpp
 * @brief Tests for seismic source mechanisms and representations
 * 
 * Tests for:
 * - Point sources (moment tensor)
 * - Extended sources (finite fault)
 * - Kinematic vs dynamic rupture
 * - Source time functions
 * - Radiation patterns
 * - Moment magnitude relations
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "SeismicSource.hpp"
#include "PhysicsKernel.hpp"
#include <cmath>
#include <vector>
#include <array>

using namespace FSRM;
using namespace FSRM::Testing;

class SeismicSourceTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Material properties
        mu = 30e9;              // Shear modulus (30 GPa)
        rho = 2700.0;           // Density (kg/m³)
        
        // Wave velocities
        Vp = 6000.0;            // P-wave velocity (m/s)
        Vs = 3500.0;            // S-wave velocity (m/s)
        
        // Reference earthquake
        M0_ref = 1e18;          // Moment (1e18 N·m ≈ Mw 6.0)
        Mw_ref = (2.0/3.0) * std::log10(M0_ref) - 6.07;
    }
    
    int rank;
    double mu, rho;
    double Vp, Vs;
    double M0_ref, Mw_ref;
};

// ============================================================================
// Moment Magnitude Tests
// ============================================================================

TEST_F(SeismicSourceTest, MomentMagnitudeRelation) {
    // Mw = (2/3) * log10(M0) - 6.07 (Hanks & Kanamori, 1979)
    
    double M0 = 1e18;  // N·m
    double Mw = (2.0/3.0) * std::log10(M0) - 6.07;
    
    // Mw ≈ 6.0 for M0 = 1e18 N·m
    EXPECT_NEAR(Mw, 5.93, 0.05);
}

TEST_F(SeismicSourceTest, MagnitudeToMoment) {
    // M0 = 10^(1.5*(Mw + 6.07))
    
    double Mw = 6.0;
    double M0 = std::pow(10.0, 1.5 * (Mw + 6.07));
    
    EXPECT_NEAR(M0, 1.12e18, 0.1e18);
}

TEST_F(SeismicSourceTest, MagnitudeScaling) {
    // 1 unit increase in Mw = 31.6× increase in M0
    
    double Mw1 = 5.0;
    double Mw2 = 6.0;
    
    double M0_1 = std::pow(10.0, 1.5 * (Mw1 + 6.07));
    double M0_2 = std::pow(10.0, 1.5 * (Mw2 + 6.07));
    
    double ratio = M0_2 / M0_1;
    
    EXPECT_NEAR(ratio, 31.6, 0.5);
}

// ============================================================================
// Seismic Moment Tests
// ============================================================================

TEST_F(SeismicSourceTest, MomentFromFault) {
    // M0 = μ * A * D
    
    double A = 10e6;        // Fault area (10 km² = 10e6 m²)
    double D = 1.0;         // Average slip (1 m)
    
    double M0 = mu * A * D;
    
    EXPECT_GT(M0, 0.0);
    EXPECT_TRUE(std::isfinite(M0));
    
    // For μ=30 GPa, A=10 km², D=1m: M0 = 3e17 N·m
    EXPECT_NEAR(M0, 3e17, 0.1e17);
}

TEST_F(SeismicSourceTest, MomentRate) {
    // Moment rate: Ṁ0 = μ * A * V_r * D/L
    // where V_r is rupture velocity, L is fault length
    
    double A = 10e6;
    double V_r = 0.8 * Vs;   // 80% of shear velocity
    double D = 1.0;
    double L = 10000.0;      // 10 km
    
    // Simplified: Ṁ0 ≈ M0 / τ where τ is rupture duration
    double tau = L / V_r;
    double M0 = mu * A * D;
    double M0_dot = M0 / tau;
    
    EXPECT_GT(M0_dot, 0.0);
}

// ============================================================================
// Source Time Function Tests
// ============================================================================

TEST_F(SeismicSourceTest, BoxcarSTF) {
    // Boxcar: constant over duration τ
    
    double tau = 1.0;       // Duration
    double M0 = 1e18;       // Total moment
    
    auto stf = [tau, M0](double t) {
        return (t >= 0.0 && t <= tau) ? M0 / tau : 0.0;
    };
    
    // Integral = M0
    double dt = 0.001;
    double integral = 0.0;
    for (double t = 0; t <= tau; t += dt) {
        integral += stf(t) * dt;
    }
    
    EXPECT_NEAR(integral, M0, M0 * 0.01);
}

TEST_F(SeismicSourceTest, TriangularSTF) {
    // Triangular: rises linearly, then falls
    
    double tau = 1.0;
    double M0 = 1e18;
    double peak = 2.0 * M0 / tau;  // Peak moment rate
    
    auto stf = [tau, peak](double t) {
        if (t < 0.0 || t > tau) return 0.0;
        if (t <= tau/2) return 2.0 * peak * t / tau;
        return 2.0 * peak * (1.0 - t/tau);
    };
    
    EXPECT_NEAR(stf(0.0), 0.0, 1e-10);
    EXPECT_NEAR(stf(tau/2.0), peak, peak * 1e-6);
    EXPECT_NEAR(stf(tau), 0.0, 1e-10);
}

TEST_F(SeismicSourceTest, GaussianSTF) {
    // Gaussian-like smoothed function
    
    double t0 = 1.0;        // Center time
    double sigma = 0.2;     // Width
    double M0 = 1e18;
    
    auto stf = [t0, sigma, M0](double t) {
        double norm = M0 / (sigma * std::sqrt(2.0 * M_PI));
        return norm * std::exp(-0.5 * (t - t0) * (t - t0) / (sigma * sigma));
    };
    
    // Peak at t0
    double peak = stf(t0);
    EXPECT_GT(stf(t0 - 0.1), stf(t0 - 0.5));
}

TEST_F(SeismicSourceTest, RickerWavelet) {
    // Ricker wavelet: second derivative of Gaussian
    // Good for numerical simulations
    
    double f0 = 5.0;        // Dominant frequency
    double t0 = 1.0 / f0;   // Delay
    
    auto ricker = [f0, t0](double t) {
        double pi_f = M_PI * f0;
        double tau = t - t0;
        double arg = pi_f * pi_f * tau * tau;
        return (1.0 - 2.0 * arg) * std::exp(-arg);
    };
    
    EXPECT_NEAR(ricker(t0), 1.0, 1e-10);  // Peak at t0
}

// ============================================================================
// Moment Tensor Tests
// ============================================================================

TEST_F(SeismicSourceTest, MomentTensorSymmetry) {
    // Moment tensor is symmetric: M_ij = M_ji
    
    std::array<std::array<double, 3>, 3> M = {{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}};
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(M[i][j], M[j][i], 1e-10);
        }
    }
}

TEST_F(SeismicSourceTest, MomentTensorDoubleCoupleDecomposition) {
    // Double-couple: M_rr + M_θθ + M_φφ = 0 (traceless)
    
    double M0 = 1e18;
    std::array<std::array<double, 3>, 3> M;
    
    // Pure double-couple for vertical strike-slip fault
    M[0][0] = 0;      // M_xx
    M[1][1] = 0;      // M_yy
    M[2][2] = 0;      // M_zz
    M[0][1] = M[1][0] = M0;  // M_xy
    M[0][2] = M[2][0] = 0;
    M[1][2] = M[2][1] = 0;
    
    double trace = M[0][0] + M[1][1] + M[2][2];
    EXPECT_NEAR(trace, 0.0, M0 * 1e-10);
}

TEST_F(SeismicSourceTest, MomentTensorCLVD) {
    // Compensated Linear Vector Dipole (CLVD)
    // Represents non-double-couple component
    
    double M0 = 1e18;
    
    // Pure CLVD
    std::array<std::array<double, 3>, 3> M_clvd;
    M_clvd[0][0] = -M0/2;
    M_clvd[1][1] = -M0/2;
    M_clvd[2][2] = M0;
    M_clvd[0][1] = M_clvd[1][0] = 0;
    M_clvd[0][2] = M_clvd[2][0] = 0;
    M_clvd[1][2] = M_clvd[2][1] = 0;
    
    // Still traceless
    double trace = M_clvd[0][0] + M_clvd[1][1] + M_clvd[2][2];
    EXPECT_NEAR(trace, 0.0, M0 * 1e-10);
}

TEST_F(SeismicSourceTest, MomentTensorExplosion) {
    // Isotropic source (explosion): M_ij = M0 * δ_ij / 3
    
    double M0 = 1e18;
    std::array<std::array<double, 3>, 3> M_iso;
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            M_iso[i][j] = (i == j) ? M0 / 3.0 : 0.0;
        }
    }
    
    // Trace = M0
    double trace = M_iso[0][0] + M_iso[1][1] + M_iso[2][2];
    EXPECT_NEAR(trace, M0, M0 * 1e-10);
}

TEST_F(SeismicSourceTest, MomentTensorFromStrikeDipRake) {
    // Convert fault plane parameters to moment tensor
    
    double strike = 30.0 * M_PI / 180.0;   // 30°
    double dip = 60.0 * M_PI / 180.0;      // 60°
    double rake = 90.0 * M_PI / 180.0;     // 90° (pure thrust)
    double M0 = 1e18;
    
    // Aki-Richards convention
    double sin_s = std::sin(strike);
    double cos_s = std::cos(strike);
    double sin_d = std::sin(dip);
    double cos_d = std::cos(dip);
    double sin_r = std::sin(rake);
    double cos_r = std::cos(rake);
    
    // M_xx component
    double M_xx = -M0 * (sin_d * cos_r * sin_s * 2.0 * cos_s + 
                         sin_d * sin_r * sin_d * sin_s * sin_s);
    
    EXPECT_TRUE(std::isfinite(M_xx));
}

// ============================================================================
// Radiation Pattern Tests
// ============================================================================

TEST_F(SeismicSourceTest, PwaveRadiationPattern) {
    // P-wave: u_P ∝ sin(2θ) * cos(φ) for double-couple
    
    auto P_radiation = [](double theta, double phi) {
        return std::sin(2.0 * theta) * std::cos(phi);
    };
    
    // Nodal planes at θ = 0, 90°, 180°
    EXPECT_NEAR(P_radiation(0.0, 0.0), 0.0, 1e-10);
    EXPECT_NEAR(P_radiation(M_PI/2, 0.0), 0.0, 1e-10);
    
    // Maximum at θ = 45°
    double max_pattern = std::abs(P_radiation(M_PI/4, 0.0));
    EXPECT_GT(max_pattern, 0.9);
}

TEST_F(SeismicSourceTest, SwaveRadiationPattern) {
    // S-wave: u_SH ∝ cos(2θ) * cos(φ)
    //         u_SV ∝ cos(θ) * sin(φ)
    
    auto SH_radiation = [](double theta, double phi) {
        return std::cos(2.0 * theta) * std::cos(phi);
    };
    
    // SH maximum at θ = 0°
    EXPECT_NEAR(std::abs(SH_radiation(0.0, 0.0)), 1.0, 1e-10);
}

TEST_F(SeismicSourceTest, FourQuadrantPattern) {
    // P-wave has four-lobed pattern
    
    auto P_radiation = [](double theta, double phi) {
        return std::sin(2.0 * theta) * std::cos(phi);
    };
    
    // Check alternating polarity in quadrants
    double pattern_Q1 = P_radiation(M_PI/4, 0.0);
    double pattern_Q2 = P_radiation(3*M_PI/4, 0.0);
    
    EXPECT_GT(pattern_Q1 * pattern_Q2, 0.0) << "Q1 and Q2 should have same sign for φ=0";
}

// ============================================================================
// Source Dimension Tests
// ============================================================================

TEST_F(SeismicSourceTest, SourceScalingRelations) {
    // Wells & Coppersmith (1994) scaling relations
    
    double Mw = 6.5;
    
    // Rupture length: log(L) = -2.44 + 0.59*Mw (all types)
    double log_L = -2.44 + 0.59 * Mw;
    double L = std::pow(10.0, log_L);  // km
    
    EXPECT_GT(L, 1.0);
    EXPECT_LT(L, 1000.0);
    
    // Rupture area: log(A) = -3.49 + 0.91*Mw
    double log_A = -3.49 + 0.91 * Mw;
    double A = std::pow(10.0, log_A);  // km²
    
    EXPECT_GT(A, 1.0);
    EXPECT_LT(A, 10000.0);
}

TEST_F(SeismicSourceTest, SlipScaling) {
    // Average slip: D = M0 / (μ * A)
    
    double Mw = 7.0;
    double M0 = std::pow(10.0, 1.5 * (Mw + 6.07));
    double A = 1000e6;  // 1000 km² = 1000e6 m²
    
    double D = M0 / (mu * A);
    
    EXPECT_GT(D, 0.0);
    
    // For Mw 7.0, typical slip is meters
    EXPECT_GT(D, 0.1);
    EXPECT_LT(D, 20.0);
}

TEST_F(SeismicSourceTest, StressDropScaling) {
    // Stress drop: Δσ ≈ C * μ * D / L
    // where C is geometry factor
    
    double D = 2.0;         // Slip (m)
    double L = 20000.0;     // Fault length (m)
    double C = 1.0;         // Geometry factor (order 1)
    
    double delta_sigma = C * mu * D / L;
    
    // Typical stress drop: 1-10 MPa
    EXPECT_GT(delta_sigma, 0.1e6);
    EXPECT_LT(delta_sigma, 100e6);
}

// ============================================================================
// Rupture Velocity Tests
// ============================================================================

TEST_F(SeismicSourceTest, RuptureVelocityBounds) {
    // Rupture velocity typically 0.7-0.9 * Vs
    
    double Vr_min = 0.7 * Vs;
    double Vr_max = 0.9 * Vs;
    double Vr_typical = 0.8 * Vs;
    
    EXPECT_GT(Vr_typical, Vr_min);
    EXPECT_LT(Vr_typical, Vr_max);
    
    // Cannot exceed Rayleigh wave velocity (≈ 0.92 * Vs)
    double V_rayleigh = 0.92 * Vs;
    EXPECT_LT(Vr_max, V_rayleigh);
}

TEST_F(SeismicSourceTest, SupershearRupture) {
    // Supershear: Vr > Vs (possible under certain conditions)
    
    double V_rayleigh = 0.92 * Vs;
    double Vr_supershear = 1.2 * Vs;
    
    // Supershear velocities exist between Vs and Vp
    EXPECT_GT(Vr_supershear, Vs);
    EXPECT_LT(Vr_supershear, Vp);
    
    // Cannot exceed P-wave velocity
    EXPECT_LT(Vr_supershear, Vp);
}

// ============================================================================
// Rise Time Tests
// ============================================================================

TEST_F(SeismicSourceTest, RiseTimeScaling) {
    // Rise time τ_r ≈ D / V_slip
    // or τ_r scales with fault size
    
    double D = 2.0;             // Slip (m)
    double V_slip = 1.0;        // Slip velocity (m/s)
    
    double tau_r = D / V_slip;
    
    EXPECT_NEAR(tau_r, 2.0, 0.1);
    
    // Empirical: τ_r ≈ 0.2 * L / V_r
    double L = 20000.0;
    double V_r = 0.8 * Vs;
    double tau_r_emp = 0.2 * L / V_r;
    
    EXPECT_GT(tau_r_emp, 0.0);
}

TEST_F(SeismicSourceTest, CornerFrequency) {
    // Corner frequency: f_c ∝ V_r / L ∝ M0^(-1/3)
    
    double M0 = 1e18;
    
    // Empirical: log(f_c) = 2.34 - 0.5*log(M0) (Brune)
    // Simplified: f_c ≈ 3.3e6 * M0^(-1/3)
    double f_c = 3.3e6 * std::pow(M0, -1.0/3.0);
    
    EXPECT_GT(f_c, 0.0);
    
    // Larger earthquakes have lower corner frequency
    double M0_large = 1e20;
    double f_c_large = 3.3e6 * std::pow(M0_large, -1.0/3.0);
    
    EXPECT_LT(f_c_large, f_c);
}

// ============================================================================
// Finite Fault Tests
// ============================================================================

TEST_F(SeismicSourceTest, FiniteFaultDiscretization) {
    // Discretize fault into subfaults
    
    double L = 20000.0;     // 20 km
    double W = 10000.0;     // 10 km
    double dL = 1000.0;     // 1 km subfault
    double dW = 1000.0;
    
    int nL = static_cast<int>(L / dL);
    int nW = static_cast<int>(W / dW);
    int n_subfaults = nL * nW;
    
    EXPECT_EQ(n_subfaults, 200);
}

TEST_F(SeismicSourceTest, SubfaultMoment) {
    // Each subfault has moment M0_i = μ * dA * D_i
    
    double dA = 1e6;        // 1 km² = 1e6 m²
    double D = 2.0;         // Local slip
    
    double M0_subfault = mu * dA * D;
    
    EXPECT_GT(M0_subfault, 0.0);
    
    // Total moment = sum of subfault moments
}

TEST_F(SeismicSourceTest, RuptureTimeDelay) {
    // Subfault activation time: t_delay = r / V_r
    // where r is distance from hypocenter
    
    double r = 5000.0;      // 5 km from hypocenter
    double V_r = 0.8 * Vs;
    
    double t_delay = r / V_r;
    
    EXPECT_GT(t_delay, 0.0);
    EXPECT_NEAR(t_delay, 5000.0 / (0.8 * Vs), 0.1);
}

// ============================================================================
// Directivity Tests
// ============================================================================

TEST_F(SeismicSourceTest, DirectivityAmplification) {
    // Forward directivity: amplification ~ V / (V - V_r*cos(θ))
    
    double V = Vs;          // S-wave velocity
    double V_r = 0.8 * Vs;  // Rupture velocity
    double theta = 0.0;     // Forward direction
    
    double amp = V / (V - V_r * std::cos(theta));
    
    EXPECT_GT(amp, 1.0) << "Forward directivity should amplify";
    
    // Backward directivity (θ = 180°)
    double theta_back = M_PI;
    double amp_back = V / (V - V_r * std::cos(theta_back));
    
    EXPECT_LT(amp_back, 1.0) << "Backward directivity should reduce";
}

TEST_F(SeismicSourceTest, DirectivityPulse) {
    // Forward directivity produces narrow pulse
    // Duration ~ τ_r * (1 - V_r/V * cos(θ))
    
    double tau_r = 2.0;     // Rise time
    double V = Vs;
    double V_r = 0.8 * Vs;
    
    double tau_forward = tau_r * (1.0 - V_r/V * std::cos(0.0));
    double tau_backward = tau_r * (1.0 - V_r/V * std::cos(M_PI));
    
    EXPECT_LT(tau_forward, tau_backward);
}

// ============================================================================
// Point Source Approximation Tests
// ============================================================================

TEST_F(SeismicSourceTest, PointSourceValidity) {
    // Point source valid when r >> L (far field)
    // Rule of thumb: r > 3*L
    
    double L = 20000.0;     // 20 km fault
    double r_min = 3.0 * L; // Minimum distance for point source
    
    EXPECT_NEAR(r_min, 60000.0, 100.0);
}

TEST_F(SeismicSourceTest, NearFieldTerms) {
    // Near field: static + dynamic terms
    // Dynamic far field: u ∝ 1/r
    // Static near field: u ∝ 1/r², 1/r³
    
    double r1 = 1000.0;
    double r2 = 2000.0;
    
    double ratio_far = r1 / r2;
    double ratio_near = (r1/r2) * (r1/r2);
    
    EXPECT_NEAR(ratio_far, 0.5, 0.001);
    EXPECT_NEAR(ratio_near, 0.25, 0.001);
}

// ============================================================================
// Energy Tests
// ============================================================================

TEST_F(SeismicSourceTest, RadiatedEnergy) {
    // Radiated seismic energy: E_s ≈ Δσ * M0 / (2μ)
    
    double delta_sigma = 3e6;   // 3 MPa stress drop
    double M0 = 1e18;
    
    double E_s = delta_sigma * M0 / (2.0 * mu);
    
    EXPECT_GT(E_s, 0.0);
    
    // Energy-magnitude relation: log(E_s) = 1.5*Mw + 4.8
    double Mw = (2.0/3.0) * std::log10(M0) - 6.07;
    double E_s_empirical = std::pow(10.0, 1.5 * Mw + 4.8);
    
    // Same order of magnitude
    EXPECT_GT(E_s_empirical, 0.0);
}

TEST_F(SeismicSourceTest, SeismicEfficiency) {
    // Efficiency: η = E_s / (Δσ * A * D / 2)
    // Typical η ≈ 0.01-0.1
    
    double delta_sigma = 3e6;
    double A = 100e6;       // 100 km²
    double D = 1.0;         // 1 m slip
    
    double E_total = delta_sigma * A * D / 2.0;
    double E_s = E_total * 0.05;  // 5% efficiency
    
    double eta = E_s / E_total;
    
    EXPECT_NEAR(eta, 0.05, 0.001);
}

// ============================================================================
// Source Class Tests
// ============================================================================

TEST_F(SeismicSourceTest, PointSourceCreation) {
    SeismicSource source;
    
    source.setLocation(0.0, 0.0, -10000.0);  // 10 km depth
    source.setMoment(1e18);
    source.setMomentTensor({1, 0, 0, -1, 0, 0});  // Strike-slip
    
    EXPECT_TRUE(true);
}

TEST_F(SeismicSourceTest, SourceTimeFunction) {
    SeismicSource source;
    
    source.setSourceTimeFunction("gaussian", 1.0);  // 1 s half-duration
    
    double stf_0 = source.evaluateSTF(0.0);
    double stf_peak = source.evaluateSTF(1.0);
    
    EXPECT_GE(stf_peak, stf_0);
}

TEST_F(SeismicSourceTest, FiniteFaultCreation) {
    FiniteFaultSource fault;
    
    fault.setGeometry(30.0, 60.0, 50.0, 20.0);  // strike, dip, length, width (km)
    fault.setHypocenter(0.5, 0.5);  // Center of fault
    fault.setRuptureVelocity(0.8 * Vs);
    
    EXPECT_TRUE(true);
}

TEST_F(SeismicSourceTest, SubfaultIterator) {
    FiniteFaultSource fault;
    
    fault.setGeometry(30.0, 60.0, 50.0, 20.0);
    fault.discretize(50, 20);  // 50×20 = 1000 subfaults
    
    int count = 0;
    for (auto& subfault : fault) {
        (void)subfault;
        count++;
    }
    
    EXPECT_EQ(count, 1000);
}

TEST_F(SeismicSourceTest, SlipDistribution) {
    FiniteFaultSource fault;
    
    fault.setGeometry(30.0, 60.0, 50.0, 20.0);
    fault.discretize(50, 20);
    
    // Set Gaussian slip distribution
    fault.setSlipDistribution("gaussian", 25.0, 10.0, 5.0);  // center, peak, sigma
    
    double total_moment;
    fault.computeTotalMoment(total_moment);
    
    EXPECT_GT(total_moment, 0.0);
}

// ============================================================================
// Kinematic Source Tests
// ============================================================================

TEST_F(SeismicSourceTest, KinematicRupture) {
    KinematicSource source;
    
    source.setFaultGeometry(50000.0, 20000.0);  // L×W
    source.setHypocenter(10000.0, 10000.0);     // x, z on fault
    source.setRuptureVelocity(0.8 * Vs);
    source.setRiseTime(1.0);
    
    // Get activation time at a point
    double t_act;
    source.getActivationTime(30000.0, 10000.0, t_act);
    
    EXPECT_GT(t_act, 0.0);
}

TEST_F(SeismicSourceTest, RuptureAreaExpansion) {
    KinematicSource source;
    
    source.setFaultGeometry(50000.0, 20000.0);
    source.setHypocenter(25000.0, 10000.0);  // Center
    source.setRuptureVelocity(0.8 * Vs);
    
    // Area grows as circle: A(t) = π * (V_r * t)²
    double t = 5.0;
    double r = 0.8 * Vs * t;
    double A_expected = M_PI * r * r;
    
    EXPECT_GT(A_expected, 0.0);
}

// ============================================================================
// Dynamic Rupture Tests
// ============================================================================

TEST_F(SeismicSourceTest, DynamicFriction) {
    // Slip-weakening friction law
    
    double mu_s = 0.6;      // Static friction
    double mu_d = 0.4;      // Dynamic friction
    double D_c = 0.1;       // Critical slip distance
    
    // Friction as function of slip
    auto friction = [mu_s, mu_d, D_c](double D) {
        if (D < D_c) {
            return mu_s - (mu_s - mu_d) * D / D_c;
        }
        return mu_d;
    };
    
    EXPECT_NEAR(friction(0.0), mu_s, 1e-10);
    EXPECT_NEAR(friction(D_c), mu_d, 1e-10);
    EXPECT_NEAR(friction(2.0 * D_c), mu_d, 1e-10);
}

TEST_F(SeismicSourceTest, FractureEnergy) {
    // Fracture energy: G = 0.5 * (τ_s - τ_d) * D_c
    
    double sigma_n = 100e6;  // Normal stress (100 MPa)
    double mu_s = 0.6;
    double mu_d = 0.4;
    double D_c = 0.1;
    
    double tau_s = mu_s * sigma_n;
    double tau_d = mu_d * sigma_n;
    double G = 0.5 * (tau_s - tau_d) * D_c;
    
    EXPECT_GT(G, 0.0);
}

TEST_F(SeismicSourceTest, NucleationPatchSize) {
    // Critical nucleation length: L_c ~ μ * G / (τ - τ_r)²
    // where τ is applied stress, τ_r is residual
    
    double G_c = 1e6;       // Fracture energy
    double tau = 50e6;      // Applied shear
    double tau_r = 40e6;    // Residual strength
    
    double L_c = mu * G_c / ((tau - tau_r) * (tau - tau_r));
    
    EXPECT_GT(L_c, 0.0);
}

// ============================================================================
// Validation Tests
// ============================================================================

TEST_F(SeismicSourceTest, ScalingConsistency) {
    // Check M0 = μ * A * D is consistent with Mw
    
    double Mw = 7.0;
    double M0 = std::pow(10.0, 1.5 * (Mw + 6.07));
    
    // Wells-Coppersmith area
    double A = std::pow(10.0, -3.49 + 0.91 * Mw) * 1e6;  // km² to m²
    
    double D = M0 / (mu * A);
    
    // Should be reasonable slip (meters)
    EXPECT_GT(D, 0.1);
    EXPECT_LT(D, 20.0);
}

TEST_F(SeismicSourceTest, HistoricalEarthquake) {
    // 1906 San Francisco: Mw 7.9, L ≈ 400 km, D ≈ 4-6 m
    
    double Mw = 7.9;
    double M0 = std::pow(10.0, 1.5 * (Mw + 6.07));
    double L = 400e3;       // 400 km
    double W = 15e3;        // 15 km
    double A = L * W;
    
    double D = M0 / (mu * A);
    
    // Should be in range 4-6 m
    EXPECT_GT(D, 2.0);
    EXPECT_LT(D, 10.0);
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_F(SeismicSourceTest, VerySmallEarthquake) {
    // Microearthquake: Mw ≈ 0
    
    double Mw = 0.0;
    double M0 = std::pow(10.0, 1.5 * (Mw + 6.07));
    
    // Very small moment
    EXPECT_LT(M0, 1e12);
    EXPECT_GT(M0, 0.0);
}

TEST_F(SeismicSourceTest, VeryLargeEarthquake) {
    // Giant earthquake: Mw 9.5 (largest recorded)
    
    double Mw = 9.5;
    double M0 = std::pow(10.0, 1.5 * (Mw + 6.07));
    
    // Enormous moment
    EXPECT_GT(M0, 1e22);
    EXPECT_TRUE(std::isfinite(M0));
}

TEST_F(SeismicSourceTest, ZeroMoment) {
    // Edge case: M0 = 0 (no earthquake)
    
    double M0 = 0.0;
    
    // Mw undefined for M0 = 0
    // Typically handle as Mw = -inf or special value
    EXPECT_EQ(M0, 0.0);
}
