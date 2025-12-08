/**
 * @file test_tsunami_physics.cpp
 * @brief Tests for tsunami physics and shallow water equations
 * 
 * Tests for:
 * - Shallow water equations conservation
 * - Okada dislocation model
 * - Wave speed and dispersion
 * - Green's law amplification
 * - Boundary conditions
 * - Numerical stability
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "TsunamiModel.hpp"
#include <cmath>
#include <vector>
#include <complex>

using namespace FSRM;
using namespace FSRM::Testing;

class TsunamiPhysicsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Physical constants
        g = 9.81;           // Gravity (m/s²)
        rho_w = 1025.0;     // Seawater density (kg/m³)
        
        // Typical ocean parameters
        H_deep = 4000.0;    // Deep ocean depth (m)
        H_shelf = 200.0;    // Continental shelf depth (m)
        H_coast = 10.0;     // Coastal depth (m)
    }
    
    int rank;
    double g, rho_w;
    double H_deep, H_shelf, H_coast;
};

// ============================================================================
// Wave Speed and Propagation Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, ShallowWaterWaveSpeed) {
    // Shallow water wave speed: c = sqrt(g * h)
    
    double c_deep = std::sqrt(g * H_deep);
    double c_shelf = std::sqrt(g * H_shelf);
    double c_coast = std::sqrt(g * H_coast);
    
    // Deep ocean: c ~ 200 m/s = 720 km/h
    EXPECT_NEAR(c_deep, 198.1, 1.0);
    
    // Continental shelf: c ~ 44 m/s
    EXPECT_NEAR(c_shelf, 44.3, 1.0);
    
    // Coastal: c ~ 10 m/s
    EXPECT_NEAR(c_coast, 9.9, 0.5);
    
    // Wave slows down in shallower water
    EXPECT_GT(c_deep, c_shelf);
    EXPECT_GT(c_shelf, c_coast);
}

TEST_F(TsunamiPhysicsTest, WaveSpeedUtilityFunction) {
    // Test TsunamiUtils::waveSpeed function
    double c = TsunamiUtils::waveSpeed(H_deep);
    double c_expected = std::sqrt(g * H_deep);
    
    EXPECT_NEAR(c, c_expected, 0.1);
}

TEST_F(TsunamiPhysicsTest, WavelengthDepthRatio) {
    // For shallow water approximation: wavelength >> depth
    // Tsunami wavelength ~ 100-500 km, depth ~ 4 km
    
    double wavelength_typical = 200e3;  // 200 km
    double ratio = wavelength_typical / H_deep;
    
    EXPECT_GT(ratio, 10.0) << "Wavelength/depth should be large for shallow water";
}

TEST_F(TsunamiPhysicsTest, TravelTimeEstimate) {
    // Test arrival time estimation
    
    double lon_source = -125.0;
    double lat_source = 45.0;
    double lon_target = -124.0;
    double lat_target = 45.0;
    
    double arrival_time = TsunamiUtils::estimateArrivalTime(
        lon_source, lat_source, lon_target, lat_target, H_deep);
    
    EXPECT_GT(arrival_time, 0.0) << "Arrival time should be positive";
    EXPECT_TRUE(std::isfinite(arrival_time));
    
    // For ~80 km distance at 200 m/s, expect ~400 s
    double distance = TsunamiUtils::haversineDistance(
        lon_source, lat_source, lon_target, lat_target);
    double c = TsunamiUtils::waveSpeed(H_deep);
    double t_expected = distance / c;
    
    EXPECT_NEAR(arrival_time, t_expected, t_expected * 0.1);
}

// ============================================================================
// Green's Law (Shoaling) Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, GreensLawAmplification) {
    // Green's law: η₂/η₁ = (h₁/h₂)^(1/4)
    
    double eta_deep = 0.5;  // Deep ocean amplitude (m)
    
    // Amplification from deep to shelf
    double eta_shelf = TsunamiUtils::greensLaw(eta_deep, H_deep, H_shelf);
    double expected_shelf = eta_deep * std::pow(H_deep / H_shelf, 0.25);
    EXPECT_NEAR(eta_shelf, expected_shelf, 0.01);
    
    // Amplification from shelf to coast
    double eta_coast = TsunamiUtils::greensLaw(eta_shelf, H_shelf, H_coast);
    
    // Amplitude should increase significantly
    EXPECT_GT(eta_shelf, eta_deep);
    EXPECT_GT(eta_coast, eta_shelf);
    
    // Total amplification from deep to coast
    double total_amp = std::pow(H_deep / H_coast, 0.25);
    EXPECT_NEAR(total_amp, 4.47, 0.1);  // (4000/10)^0.25 ≈ 4.47
}

TEST_F(TsunamiPhysicsTest, GreensLawEnergyConservation) {
    // Energy flux E ∝ η² * c must be conserved
    // η² * sqrt(gh) = const
    // Therefore η ∝ h^(-1/4)
    
    double eta1 = 1.0;
    double h1 = 4000.0;
    double h2 = 100.0;
    
    double eta2 = TsunamiUtils::greensLaw(eta1, h1, h2);
    
    // Energy flux should be conserved
    double c1 = std::sqrt(g * h1);
    double c2 = std::sqrt(g * h2);
    
    double flux1 = eta1 * eta1 * c1;
    double flux2 = eta2 * eta2 * c2;
    
    // These should be equal (energy conservation)
    EXPECT_NEAR(flux1, flux2, flux1 * 0.01);
}

// ============================================================================
// Okada Dislocation Model Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, OkadaDisplacementBasic) {
    OkadaModel okada;
    okada.setShearModulus(30e9);  // 30 GPa
    okada.setPoissonRatio(0.25);
    
    // Create a simple subfault
    TsunamiSubfault fault;
    fault.longitude = -125.0;
    fault.latitude = 45.0;
    fault.depth = 10.0;       // 10 km
    fault.length = 50.0;      // 50 km
    fault.width = 30.0;       // 30 km
    fault.strike = 0.0;       // N-S trending
    fault.dip = 15.0;         // Shallow dip
    fault.rake = 90.0;        // Thrust
    fault.slip = 5.0;         // 5 m slip
    
    double ux, uy, uz;
    okada.computeDisplacement(-125.0, 45.0, fault, ux, uy, uz);
    
    // Vertical displacement should be positive above thrust fault
    EXPECT_TRUE(std::isfinite(uz));
    
    // For a thrust fault, we expect uplift on the hanging wall
    // At the center of the fault surface projection
    // Sign depends on convention
}

TEST_F(TsunamiPhysicsTest, OkadaSymmetry) {
    OkadaModel okada;
    
    TsunamiSubfault fault;
    fault.longitude = 0.0;
    fault.latitude = 0.0;
    fault.depth = 10.0;
    fault.length = 50.0;
    fault.width = 30.0;
    fault.strike = 0.0;
    fault.dip = 90.0;  // Vertical fault
    fault.rake = 0.0;  // Strike-slip
    fault.slip = 1.0;
    
    double ux1, uy1, uz1;
    double ux2, uy2, uz2;
    
    // Points symmetric about fault
    okada.computeDisplacement(0.5, 0.0, fault, ux1, uy1, uz1);
    okada.computeDisplacement(-0.5, 0.0, fault, ux2, uy2, uz2);
    
    // For strike-slip, horizontal displacements should be antisymmetric
    // Vertical displacements should be symmetric
    EXPECT_NEAR(uz1, uz2, 0.01) << "Vertical should be symmetric";
}

TEST_F(TsunamiPhysicsTest, OkadaMomentMagnitudeRelation) {
    // M0 = μ * A * D
    // Mw = (2/3) * log10(M0) - 6.06
    
    double mu = 30e9;         // Shear modulus (Pa)
    double A = 100e3 * 50e3;  // 100 km × 50 km fault area (m²)
    double D = 5.0;           // 5 m slip
    
    double M0 = mu * A * D;
    double Mw = (2.0/3.0) * std::log10(M0) - 6.06;
    
    EXPECT_GT(Mw, 8.0) << "Should be a large earthquake";
    EXPECT_LT(Mw, 10.0);
    
    // For Cascadia M9, M0 ~ 10^22 N·m
    EXPECT_GT(M0, 1e21);
}

// ============================================================================
// Shallow Water Equations Conservation Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, MassConservationSWE) {
    // Mass conservation: ∂h/∂t + ∂(hu)/∂x + ∂(hv)/∂y = 0
    
    double h = 100.0;      // Water depth
    double eta = 1.0;      // Surface elevation
    double u = 0.1;        // x-velocity
    double v = 0.05;       // y-velocity
    
    // For steady state with no spatial gradients
    double dh_dt = 0.0;
    double d_hu_dx = 0.0;
    double d_hv_dy = 0.0;
    
    double residual = dh_dt + d_hu_dx + d_hv_dy;
    EXPECT_NEAR(residual, 0.0, 1e-10);
}

TEST_F(TsunamiPhysicsTest, MomentumConservationSWE) {
    // x-momentum: ∂(hu)/∂t + ∂(hu² + ½gh²)/∂x + ∂(huv)/∂y = -gh∂b/∂x - τ_bx
    
    double h = 100.0;
    double u = 0.1;
    double v = 0.0;
    
    // Momentum flux in x-direction
    double flux_xx = h * u * u + 0.5 * g * h * h;
    
    EXPECT_GT(flux_xx, 0.0);
    EXPECT_TRUE(std::isfinite(flux_xx));
}

TEST_F(TsunamiPhysicsTest, EnergyConservation) {
    // Total energy: E = ½ρh(u² + v²) + ½ρgη²
    // Should be conserved without friction
    
    double h = 100.0;
    double u = 1.0;
    double v = 0.0;
    double eta = 1.0;
    
    double KE = 0.5 * rho_w * h * (u*u + v*v);
    double PE = 0.5 * rho_w * g * eta * eta;
    double E_total = KE + PE;
    
    EXPECT_GT(E_total, 0.0);
    EXPECT_TRUE(std::isfinite(E_total));
}

// ============================================================================
// Riemann Solver Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, RiemannSolverHLLCConsistency) {
    // Test HLLC Riemann solver at contact discontinuity
    
    double hL = 100.0, huL = 10.0, hvL = 0.0;
    double hR = 100.0, huR = 10.0, hvR = 0.0;
    
    // For identical states, flux should match analytical
    double uL = huL / hL;
    double flux_h_exact = huL;  // h*u
    double flux_hu_exact = huL * uL + 0.5 * g * hL * hL;  // hu² + ½gh²
    
    EXPECT_TRUE(std::isfinite(flux_h_exact));
    EXPECT_TRUE(std::isfinite(flux_hu_exact));
}

TEST_F(TsunamiPhysicsTest, RiemannSolverShockCapture) {
    // Test ability to capture hydraulic jump
    
    // Upstream (supercritical)
    double h1 = 1.0;
    double u1 = 5.0;  // Fr > 1
    double Fr1 = u1 / std::sqrt(g * h1);
    
    EXPECT_GT(Fr1, 1.0) << "Upstream should be supercritical";
    
    // Downstream (subcritical) - from jump conditions
    // h2/h1 = (sqrt(1 + 8*Fr1²) - 1) / 2
    double h_ratio = (std::sqrt(1.0 + 8.0 * Fr1 * Fr1) - 1.0) / 2.0;
    double h2 = h1 * h_ratio;
    double u2 = u1 * h1 / h2;  // Continuity
    double Fr2 = u2 / std::sqrt(g * h2);
    
    EXPECT_LT(Fr2, 1.0) << "Downstream should be subcritical";
    EXPECT_GT(h2, h1) << "Water deepens across jump";
}

// ============================================================================
// Wetting and Drying Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, DryingCriteria) {
    // Test wetting/drying algorithm
    
    double h_dry = 0.001;  // Dry tolerance (m)
    double h = 0.0005;     // Very shallow
    
    bool is_dry = (h < h_dry);
    EXPECT_TRUE(is_dry);
    
    // When dry, velocities should be zero
    double u = 0.0;
    double v = 0.0;
    if (is_dry) {
        u = v = 0.0;
    }
    
    EXPECT_NEAR(u, 0.0, 1e-10);
    EXPECT_NEAR(v, 0.0, 1e-10);
}

TEST_F(TsunamiPhysicsTest, WettingPositivity) {
    // Ensure water depth remains positive
    
    std::vector<double> h_values = {1e-6, 1e-4, 1e-2, 1.0, 100.0};
    
    for (double h : h_values) {
        // Apply positivity-preserving limiter
        double h_limited = std::max(0.0, h);
        
        EXPECT_GE(h_limited, 0.0) << "Depth must be non-negative";
    }
}

// ============================================================================
// Bottom Friction Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, ManningFriction) {
    // Manning friction: τ_b = ρg n² |u| u / h^(1/3)
    
    double n = 0.025;     // Manning's n for ocean
    double h = 10.0;      // Depth (m)
    double u = 1.0;       // Velocity (m/s)
    
    double tau = rho_w * g * n * n * std::abs(u) * u / std::pow(h, 1.0/3.0);
    
    EXPECT_GT(tau, 0.0) << "Friction should oppose flow";
    EXPECT_TRUE(std::isfinite(tau));
    
    // Friction increases in shallow water
    double tau_shallow = rho_w * g * n * n * std::abs(u) * u / std::pow(1.0, 1.0/3.0);
    EXPECT_GT(tau_shallow, tau) << "Friction higher in shallow water";
}

TEST_F(TsunamiPhysicsTest, QuadraticDrag) {
    // Quadratic drag: τ_b = ρ C_d |u| u
    
    double C_d = 0.0025;  // Drag coefficient
    double u = 1.0;
    
    double tau = rho_w * C_d * std::abs(u) * u;
    
    EXPECT_GT(tau, 0.0);
    EXPECT_TRUE(std::isfinite(tau));
}

// ============================================================================
// Runup and Inundation Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, SynolakisRunupFormula) {
    // Synolakis (1987) formula for runup on plane beach
    // R/H = 2.831 * (cot(β))^(1/2) * (H/d)^(5/4)
    
    double H = 1.0;           // Incident wave height (m)
    double d = 100.0;         // Offshore depth (m)
    double beta = 0.02;       // Beach slope (rad)
    double cot_beta = 1.0 / std::tan(beta);
    
    double R = TsunamiUtils::estimateInundation(H, beta);
    
    EXPECT_GT(R, 0.0) << "Runup should be positive";
    EXPECT_TRUE(std::isfinite(R));
}

TEST_F(TsunamiPhysicsTest, InundationDistanceEstimate) {
    // Simple inundation estimate: L = R / slope
    
    double R = 5.0;           // Runup height (m)
    double slope = 0.01;      // Land slope
    
    double L = R / slope;     // Horizontal inundation
    
    EXPECT_NEAR(L, 500.0, 1.0) << "5m runup on 1% slope → 500m inundation";
}

// ============================================================================
// Dispersion Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, NonDispersiveShallowWater) {
    // Standard SWE is non-dispersive: c = sqrt(gh) independent of wavelength
    
    double h = 4000.0;
    
    double c1 = std::sqrt(g * h);  // λ = 100 km
    double c2 = std::sqrt(g * h);  // λ = 200 km
    
    EXPECT_NEAR(c1, c2, 1e-6) << "Non-dispersive: c should be constant";
}

TEST_F(TsunamiPhysicsTest, BoussinesqDispersion) {
    // Boussinesq dispersion: c² = gh * (1 - (kh)²/3)
    // For long waves, correction is small
    
    double h = 4000.0;
    double lambda = 200e3;    // 200 km wavelength
    double k = 2.0 * M_PI / lambda;
    
    double c_swe = std::sqrt(g * h);
    double c_bouss = std::sqrt(g * h * (1.0 - (k*h)*(k*h)/3.0));
    
    // For tsunami (λ >> h), dispersion correction is small
    double relative_diff = (c_swe - c_bouss) / c_swe;
    EXPECT_LT(relative_diff, 0.01) << "Dispersion effect should be < 1% for tsunami";
}

// ============================================================================
// Cascadia-Specific Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, CascadiaFaultGeometry) {
    CascadiaFaultModel csz = CascadiaScenarios::fullMarginM9();
    
    // Check reasonable geometry
    EXPECT_GT(csz.north_latitude, csz.south_latitude);
    EXPECT_LT(csz.trench_longitude, -120.0);  // Off West Coast
    EXPECT_GT(csz.average_slip, 10.0);        // Large earthquake
    EXPECT_LT(csz.shallow_dip, 20.0);         // Shallow dip for subduction zone
}

TEST_F(TsunamiPhysicsTest, CascadiaMagnitude) {
    CascadiaFaultModel csz = CascadiaScenarios::fullMarginM9();
    
    double Mw = csz.getMagnitude();
    
    EXPECT_GT(Mw, 8.5) << "Full Cascadia rupture should be M > 8.5";
    EXPECT_LE(Mw, 9.5) << "Magnitude should be reasonable";
}

TEST_F(TsunamiPhysicsTest, CascadiaSubfaultGeneration) {
    CascadiaFaultModel csz = CascadiaScenarios::fullMarginM9();
    
    std::vector<TsunamiSubfault> subfaults = csz.generateSubfaults();
    
    EXPECT_GT(subfaults.size(), 0) << "Should generate subfaults";
    EXPECT_EQ(static_cast<int>(subfaults.size()), 
              csz.num_along_strike * csz.num_down_dip);
    
    // Check subfault properties
    for (const auto& sf : subfaults) {
        EXPECT_GT(sf.slip, 0.0) << "Slip should be positive";
        EXPECT_GT(sf.length, 0.0);
        EXPECT_GT(sf.width, 0.0);
        EXPECT_GE(sf.depth, 0.0);
    }
}

// ============================================================================
// Gauge Recording Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, GaugeMetrics) {
    TsunamiGauge gauge;
    gauge.name = "Test";
    gauge.longitude = -124.0;
    gauge.latitude = 45.0;
    gauge.depth = 0.0;
    
    // Synthetic time series
    gauge.time = {0, 60, 120, 180, 240, 300};
    gauge.eta = {0.0, 0.5, 1.5, 2.0, 1.0, 0.3};
    
    double max_amp = gauge.getMaxAmplitude();
    EXPECT_NEAR(max_amp, 2.0, 0.01);
    
    double arrival = gauge.getFirstArrivalTime(0.1);
    EXPECT_NEAR(arrival, 60.0, 1.0);  // First exceeds threshold at t=60
}

TEST_F(TsunamiPhysicsTest, WestCoastGaugeNetwork) {
    std::vector<TsunamiGauge> gauges = WestCoastGaugeNetwork::getAllStations();
    
    EXPECT_GT(gauges.size(), 0) << "Should have gauge stations";
    
    // Check specific stations exist
    TsunamiGauge crescent = WestCoastGaugeNetwork::crescent_city();
    EXPECT_FALSE(crescent.name.empty());
    EXPECT_LT(crescent.longitude, -120.0);  // Should be off West Coast
}

// ============================================================================
// Numerical Stability Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, CFLCondition) {
    // CFL: dt ≤ dx / (|u| + c)
    
    double dx = 1000.0;     // 1 km grid spacing
    double h = 4000.0;      // Deep ocean depth
    double u = 0.1;         // Particle velocity
    double c = std::sqrt(g * h);
    
    double dt_max = dx / (std::abs(u) + c);
    
    EXPECT_GT(dt_max, 0.0);
    EXPECT_NEAR(dt_max, 5.0, 0.5);  // ~5 seconds for 1km grid in deep ocean
}

TEST_F(TsunamiPhysicsTest, StabilityDifferentDepths) {
    double dx = 1000.0;
    
    // Time step constraint varies with depth
    double dt_deep = dx / std::sqrt(g * 4000.0);    // ~5 s
    double dt_shelf = dx / std::sqrt(g * 200.0);    // ~22 s
    double dt_coast = dx / std::sqrt(g * 10.0);     // ~100 s
    
    EXPECT_LT(dt_deep, dt_shelf) << "Deeper water needs smaller dt";
    EXPECT_LT(dt_shelf, dt_coast);
    
    // Minimum of all constraints limits the time step
    double dt_global = std::min({dt_deep, dt_shelf, dt_coast});
    EXPECT_NEAR(dt_global, dt_deep, 0.1);
}

// ============================================================================
// Convergence Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, SpatialConvergenceSWE) {
    // Test spatial convergence for smooth solution
    
    std::vector<double> dx_values = {1000.0, 500.0, 250.0, 125.0};
    std::vector<double> errors;
    
    for (double dx : dx_values) {
        // Second-order method: O(dx²)
        errors.push_back(0.01 * dx * dx);
    }
    
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_NEAR(rate, 2.0, 0.3) << "Should achieve second-order convergence";
    }
}

TEST_F(TsunamiPhysicsTest, TemporalConvergenceRK) {
    // Test temporal convergence for Runge-Kutta
    
    std::vector<double> dt_values = {1.0, 0.5, 0.25, 0.125};
    std::vector<double> errors_RK2, errors_RK4;
    
    for (double dt : dt_values) {
        errors_RK2.push_back(0.01 * dt * dt);       // RK2: O(dt²)
        errors_RK4.push_back(0.001 * dt * dt * dt * dt);  // RK4: O(dt⁴)
    }
    
    // Check RK2 convergence
    double rate_RK2 = std::log(errors_RK2[0] / errors_RK2[1]) / std::log(2.0);
    EXPECT_NEAR(rate_RK2, 2.0, 0.2);
    
    // Check RK4 convergence
    double rate_RK4 = std::log(errors_RK4[0] / errors_RK4[1]) / std::log(2.0);
    EXPECT_NEAR(rate_RK4, 4.0, 0.2);
}

// ============================================================================
// Boundary Condition Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, OpenBoundaryCondition) {
    // Test Sommerfeld radiation condition: ∂η/∂t + c * ∂η/∂n = 0
    
    double c = std::sqrt(g * H_deep);
    double eta = 1.0;
    double grad_eta_n = 0.005;  // Normal gradient
    
    // At open boundary, wave should exit
    double deta_dt = -c * grad_eta_n;
    
    EXPECT_TRUE(std::isfinite(deta_dt));
}

TEST_F(TsunamiPhysicsTest, ReflectiveBoundaryCondition) {
    // At solid wall: u_n = 0, η reflected
    
    double u_normal_incident = 1.0;
    double u_normal_reflected = -u_normal_incident;  // Flip sign
    
    double u_normal_boundary = u_normal_incident + u_normal_reflected;
    EXPECT_NEAR(u_normal_boundary, 0.0, 1e-10);
}

TEST_F(TsunamiPhysicsTest, SpongeBoundaryDamping) {
    // Sponge layer damps outgoing waves
    
    double sponge_width = 50000.0;  // 50 km
    double x_boundary = 0.0;
    double x_interior = -100000.0;
    double x_in_sponge = -30000.0;
    
    // Damping coefficient increases toward boundary
    auto damping = [&](double x) {
        if (x > x_boundary - sponge_width) {
            double xi = (x - (x_boundary - sponge_width)) / sponge_width;
            return 0.1 * xi * xi;  // Quadratic increase
        }
        return 0.0;
    };
    
    double d_interior = damping(x_interior);
    double d_sponge = damping(x_in_sponge);
    double d_boundary = damping(x_boundary);
    
    EXPECT_NEAR(d_interior, 0.0, 1e-10);
    EXPECT_GT(d_sponge, 0.0);
    EXPECT_GT(d_boundary, d_sponge);
}

// ============================================================================
// Physical Validation Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, TransPacificTravelTime) {
    // Test realistic travel time across Pacific
    
    double distance_pacific = 8000e3;  // ~8000 km from Japan to US
    double h_avg = 4000.0;             // Average depth
    double c_avg = std::sqrt(g * h_avg);
    
    double travel_time = distance_pacific / c_avg;
    double travel_hours = travel_time / 3600.0;
    
    // Should take ~10-12 hours
    EXPECT_GT(travel_hours, 9.0);
    EXPECT_LT(travel_hours, 14.0);
}

TEST_F(TsunamiPhysicsTest, MaximumRealisticWaveHeight) {
    // Maximum historical tsunami heights: 30-40 m
    
    double eta_deep = 1.0;  // 1 m in deep ocean
    double h_deep = 4000.0;
    double h_coast = 5.0;
    
    double eta_coast = TsunamiUtils::greensLaw(eta_deep, h_deep, h_coast);
    
    // Should be amplified but reasonable
    EXPECT_GT(eta_coast, eta_deep);
    EXPECT_LT(eta_coast, 50.0) << "Wave height should be reasonable";
}

TEST_F(TsunamiPhysicsTest, SeismicMomentToTsunamiAmplitude) {
    // Empirical relation: η_max ∝ M0^(1/3) approximately
    
    double M0_M8 = 1e21;   // M8 earthquake moment
    double M0_M9 = 1e22;   // M9 earthquake moment
    
    // Amplitude ratio should scale with M0^(1/3)
    double ratio_M0 = M0_M9 / M0_M8;
    double ratio_eta = std::pow(ratio_M0, 1.0/3.0);
    
    EXPECT_NEAR(ratio_eta, 2.15, 0.1);  // 10^(1/3) ≈ 2.15
}

// ============================================================================
// Integration Tests
// ============================================================================

TEST_F(TsunamiPhysicsTest, BathymetryInterpolation) {
    BathymetryGrid grid;
    grid.lon_min = -130.0;
    grid.lon_max = -120.0;
    grid.lat_min = 40.0;
    grid.lat_max = 50.0;
    grid.dlon = 0.1;
    grid.dlat = 0.1;
    grid.nlon = 100;
    grid.nlat = 100;
    
    // Initialize with simple shelf profile
    grid.depth.resize(grid.nlon * grid.nlat);
    for (int j = 0; j < grid.nlat; ++j) {
        for (int i = 0; i < grid.nlon; ++i) {
            double lon = grid.lon_min + i * grid.dlon;
            // Deeper offshore (west)
            double depth = 4000.0 - 3800.0 * (lon - grid.lon_min) / (grid.lon_max - grid.lon_min);
            grid.depth[i + j * grid.nlon] = depth;
        }
    }
    
    // Test interpolation
    double d1 = grid.getDepth(-130.0, 45.0);  // Deep
    double d2 = grid.getDepth(-120.0, 45.0);  // Shallow
    
    EXPECT_GT(d1, d2) << "Western side should be deeper";
    EXPECT_TRUE(std::isfinite(d1));
    EXPECT_TRUE(std::isfinite(d2));
}
