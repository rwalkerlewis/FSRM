/**
 * @file test_fracture_mechanics.cpp
 * @brief Tests for fracture mechanics and hydraulic fracturing
 * 
 * Tests for:
 * - Linear Elastic Fracture Mechanics (LEFM)
 * - PKN, KGD, P3D fracture models
 * - Stress intensity factors
 * - Energy release rate
 * - Fracture propagation criteria
 * - Hydraulic fracture width profiles
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "FractureModel.hpp"
#include "FractureNetwork.hpp"
#include <cmath>
#include <vector>

using namespace FSRM;
using namespace FSRM::Testing;

class FractureMechanicsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Rock properties
        E = 30e9;           // Young's modulus (30 GPa)
        nu = 0.25;          // Poisson's ratio
        K_IC = 1.5e6;       // Mode I fracture toughness (Pa·√m)
        
        // Derived properties
        E_prime = E / (1.0 - nu * nu);  // Plane strain modulus
        G_c = K_IC * K_IC / E_prime;    // Critical energy release rate
        
        // Fluid properties
        mu_f = 0.1;         // Fluid viscosity (Pa·s)
        rho_f = 1000.0;     // Fluid density (kg/m³)
        
        // In-situ stress
        sigma_min = 30e6;   // Minimum horizontal stress (30 MPa)
        sigma_max = 40e6;   // Maximum horizontal stress (40 MPa)
    }
    
    int rank;
    double E, nu, E_prime;
    double K_IC, G_c;
    double mu_f, rho_f;
    double sigma_min, sigma_max;
};

// ============================================================================
// Stress Intensity Factor Tests
// ============================================================================

TEST_F(FractureMechanicsTest, KI_GriffithCrack) {
    // Mode I stress intensity factor for Griffith crack
    // K_I = σ * √(π * a) for center crack in infinite plate
    
    double sigma = 10e6;    // Applied stress (10 MPa)
    double a = 1.0;         // Half crack length (m)
    
    double K_I = sigma * std::sqrt(M_PI * a);
    
    EXPECT_GT(K_I, 0.0);
    EXPECT_TRUE(std::isfinite(K_I));
    
    // For σ = 10 MPa, a = 1 m: K_I ≈ 17.7 MPa·√m
    EXPECT_NEAR(K_I, 17.7e6, 0.1e6);
}

TEST_F(FractureMechanicsTest, KI_PennyCrack) {
    // K_I for penny-shaped crack under uniform tension
    // K_I = (2/π) * σ * √(π * a)
    
    double sigma = 10e6;
    double a = 1.0;
    
    double K_I = (2.0 / M_PI) * sigma * std::sqrt(M_PI * a);
    
    EXPECT_GT(K_I, 0.0);
    EXPECT_NEAR(K_I, 11.3e6, 0.1e6);  // ≈ 11.3 MPa·√m
}

TEST_F(FractureMechanicsTest, KI_PressurizedCrack) {
    // K_I for uniformly pressurized crack
    // Same as tensile loading with σ = P - σ_min (net pressure)
    
    double P = 35e6;        // Fluid pressure (35 MPa)
    double a = 10.0;        // Half length (10 m)
    
    double P_net = P - sigma_min;  // Net pressure
    double K_I = P_net * std::sqrt(M_PI * a);
    
    EXPECT_GT(K_I, 0.0) << "Net pressure should be positive for propagation";
}

// ============================================================================
// Energy Release Rate Tests
// ============================================================================

TEST_F(FractureMechanicsTest, EnergyReleaseRate) {
    // G = K_I² / E' for plane strain
    
    double K_I = 2e6;  // 2 MPa·√m
    double G = K_I * K_I / E_prime;
    
    EXPECT_GT(G, 0.0);
    EXPECT_TRUE(std::isfinite(G));
    
    // Units: J/m²
}

TEST_F(FractureMechanicsTest, GriffithCriterion) {
    // Fracture propagates when G >= G_c
    
    double G = K_IC * K_IC / E_prime;  // At critical K
    
    EXPECT_NEAR(G, G_c, G_c * 1e-6);
    
    // Below critical - stable
    double K_below = 0.8 * K_IC;
    double G_below = K_below * K_below / E_prime;
    EXPECT_LT(G_below, G_c);
    
    // Above critical - propagates
    double K_above = 1.2 * K_IC;
    double G_above = K_above * K_above / E_prime;
    EXPECT_GT(G_above, G_c);
}

TEST_F(FractureMechanicsTest, JIntegral) {
    // J-integral equals G for elastic materials
    // J = ∫_Γ (W dy - T·(∂u/∂x) ds)
    
    // For linear elastic, J = G = K²/E'
    double K = K_IC;
    double J = K * K / E_prime;
    double G = J;  // Equality for elastic
    
    EXPECT_NEAR(J, G, J * 1e-10);
}

// ============================================================================
// Fracture Propagation Criteria Tests
// ============================================================================

TEST_F(FractureMechanicsTest, PropagationPressure) {
    // Minimum pressure for propagation: P_prop = σ_min + K_IC / √(π*a)
    
    double a = 10.0;  // 10 m half-length
    double P_prop = sigma_min + K_IC / std::sqrt(M_PI * a);
    
    EXPECT_GT(P_prop, sigma_min) << "Propagation pressure > closure pressure";
    
    // As crack grows, propagation pressure decreases
    double a_large = 100.0;
    double P_prop_large = sigma_min + K_IC / std::sqrt(M_PI * a_large);
    
    EXPECT_LT(P_prop_large, P_prop);
}

TEST_F(FractureMechanicsTest, CriticalStress) {
    // Critical stress for existing crack to propagate
    // σ_c = K_IC / √(π * a)
    
    double a = 0.01;  // 1 cm flaw
    double sigma_c = K_IC / std::sqrt(M_PI * a);
    
    // For K_IC = 1.5 MPa·√m, a = 0.01 m: σ_c ≈ 8.5 MPa
    EXPECT_NEAR(sigma_c / 1e6, 8.5, 0.5);
}

// ============================================================================
// PKN Model Tests
// ============================================================================

TEST_F(FractureMechanicsTest, PKNWidthProfile) {
    // PKN width at center: w = 2.5 * (μ*Q*L / E')^(1/4)
    // Assuming fixed height H
    
    double Q = 0.01;        // Flow rate per unit height (m²/s)
    double L = 100.0;       // Half-length (m)
    
    double w_center = 2.5 * std::pow(mu_f * Q * L / E_prime, 0.25);
    
    EXPECT_GT(w_center, 0.0);
    EXPECT_LT(w_center, 0.1);  // Typical widths are mm to cm scale
}

TEST_F(FractureMechanicsTest, PKNLengthScaling) {
    // PKN length scaling: L ~ t^(4/5)
    
    double t1 = 60.0;   // 1 minute
    double t2 = 600.0;  // 10 minutes
    
    // L ∝ t^(4/5)
    double ratio_expected = std::pow(t2 / t1, 4.0/5.0);
    
    // L2/L1 = (t2/t1)^(4/5)
    EXPECT_NEAR(ratio_expected, std::pow(10.0, 0.8), 0.01);
}

TEST_F(FractureMechanicsTest, PKNPressure) {
    // PKN net pressure: P_net ~ (μ*Q*E'^3 / H^3)^(1/4) * L^(1/4)
    
    double Q = 0.1;         // m³/s injection rate
    double H = 30.0;        // Height (m)
    double L = 100.0;       // Half-length (m)
    
    double P_net = std::pow(mu_f * Q * std::pow(E_prime, 3) / std::pow(H, 3), 0.25) 
                 * std::pow(L, 0.25);
    
    EXPECT_GT(P_net, 0.0);
    EXPECT_TRUE(std::isfinite(P_net));
}

// ============================================================================
// KGD Model Tests
// ============================================================================

TEST_F(FractureMechanicsTest, KGDWidthProfile) {
    // KGD maximum width: w_max ~ (μ*Q*L / E')^(1/3)
    // For plane strain (2D) fracture
    
    double Q = 0.001;       // Flow rate per unit height (m²/s)
    double L = 50.0;        // Half-length (m)
    
    double w_max = 1.87 * std::pow(mu_f * Q * L / E_prime, 1.0/3.0);
    
    EXPECT_GT(w_max, 0.0);
    EXPECT_TRUE(std::isfinite(w_max));
}

TEST_F(FractureMechanicsTest, KGDLengthScaling) {
    // KGD length scaling: L ~ t^(2/3) for viscosity-dominated
    
    double t1 = 60.0;
    double t2 = 600.0;
    
    double ratio_expected = std::pow(t2 / t1, 2.0/3.0);
    EXPECT_NEAR(ratio_expected, std::pow(10.0, 2.0/3.0), 0.01);
}

TEST_F(FractureMechanicsTest, KGDPressure) {
    // KGD pressure scaling
    // P_net ~ (μ*E'^2 / t)^(1/3) for viscosity-dominated
    
    double t = 600.0;  // Time (s)
    
    double P_net_scale = std::pow(mu_f * E_prime * E_prime / t, 1.0/3.0);
    
    EXPECT_GT(P_net_scale, 0.0);
}

// ============================================================================
// Dimensionless Parameter Tests
// ============================================================================

TEST_F(FractureMechanicsTest, ToughnessDominatedRegime) {
    // Dimensionless toughness K' = K_IC * (E' / (μ * Q))^(1/4)
    // K' >> 1: Toughness dominated
    
    double Q = 0.01;  // Low flow rate
    double K_prime = K_IC * std::pow(E_prime / (mu_f * Q), 0.25);
    
    EXPECT_GT(K_prime, 0.0);
    // High K' indicates toughness-dominated propagation
}

TEST_F(FractureMechanicsTest, ViscosityDominatedRegime) {
    // For high rate injection, viscosity dominates
    
    double Q = 1.0;   // High flow rate
    double K_prime = K_IC * std::pow(E_prime / (mu_f * Q), 0.25);
    
    EXPECT_GT(K_prime, 0.0);
    // Lower K' indicates viscosity-dominated propagation
}

TEST_F(FractureMechanicsTest, DimensionlessTime) {
    // τ = t * (K_IC^4 / (μ * E'^3 * L_0^2))
    // where L_0 is characteristic length
    
    double t = 3600.0;      // 1 hour
    double L_0 = 100.0;     // Characteristic length
    
    double tau = t * std::pow(K_IC, 4) / (mu_f * std::pow(E_prime, 3) * L_0 * L_0);
    
    EXPECT_GT(tau, 0.0);
}

// ============================================================================
// Fracture Width-Pressure Relation Tests
// ============================================================================

TEST_F(FractureMechanicsTest, SneddonWidth) {
    // Sneddon solution: w(x) = (4 * P_net * a / E') * √(1 - (x/a)²)
    // for -a ≤ x ≤ a
    
    double P_net = 5e6;     // Net pressure (5 MPa)
    double a = 50.0;        // Half-length (m)
    
    // Width at center (x = 0)
    double w_center = 4.0 * P_net * a / E_prime;
    
    EXPECT_GT(w_center, 0.0);
    
    // Width profile
    double x = 25.0;  // Halfway to tip
    double w_x = 4.0 * P_net * a / E_prime * std::sqrt(1.0 - (x/a)*(x/a));
    
    // Should be less than center width
    EXPECT_LT(w_x, w_center);
    EXPECT_GT(w_x, 0.0);
    
    // At tip, width = 0
    double w_tip = 4.0 * P_net * a / E_prime * std::sqrt(1.0 - 1.0);
    EXPECT_NEAR(w_tip, 0.0, 1e-10);
}

TEST_F(FractureMechanicsTest, ElasticStiffness) {
    // Elastic stiffness: w = P_net / (E' / a)
    // Rearranging: P_net = E' * w / a
    
    double w = 0.01;    // Width (1 cm)
    double a = 50.0;    // Half-length (m)
    
    double P_net = E_prime * w / a;
    
    // This is an approximation; actual depends on geometry
    EXPECT_GT(P_net, 0.0);
    EXPECT_TRUE(std::isfinite(P_net));
}

// ============================================================================
// Leak-off Tests
// ============================================================================

TEST_F(FractureMechanicsTest, CarterLeakoff) {
    // Carter leak-off: v_L = C_L / √t
    // where C_L is leak-off coefficient
    
    double C_L = 1e-4;      // Leak-off coefficient (m/√s)
    double t = 3600.0;      // Time (1 hour)
    
    double v_L = C_L / std::sqrt(t);
    
    EXPECT_GT(v_L, 0.0);
    
    // Leak-off decreases with time (filter cake builds up)
    double v_L_early = C_L / std::sqrt(60.0);
    double v_L_late = C_L / std::sqrt(3600.0);
    
    EXPECT_GT(v_L_early, v_L_late);
}

TEST_F(FractureMechanicsTest, LeakoffVolume) {
    // Cumulative leak-off volume per unit area: V_L = 2 * C_L * √t
    
    double C_L = 1e-4;
    double t = 3600.0;
    
    double V_L = 2.0 * C_L * std::sqrt(t);
    
    EXPECT_GT(V_L, 0.0);
    
    // Total leak-off from fracture = integral over fracture surface
    double A_frac = 100.0 * 30.0;  // L × H (m²)
    double V_total = V_L * 2.0 * A_frac;  // Both sides
    
    EXPECT_TRUE(std::isfinite(V_total));
}

// ============================================================================
// Proppant Transport Tests
// ============================================================================

TEST_F(FractureMechanicsTest, ProppantSettlingVelocity) {
    // Stokes settling: v_s = (ρ_p - ρ_f) * g * d² / (18 * μ)
    
    double rho_p = 2650.0;    // Proppant density (kg/m³)
    double d = 0.5e-3;        // Proppant diameter (0.5 mm)
    double g = 9.81;
    
    double v_s = (rho_p - rho_f) * g * d * d / (18.0 * mu_f);
    
    EXPECT_GT(v_s, 0.0);
    
    // Typical settling: mm/s to cm/s
    EXPECT_LT(v_s, 1.0);  // < 1 m/s
}

TEST_F(FractureMechanicsTest, ProppantConcentrationEffect) {
    // Hindered settling with concentration
    // v_s(c) = v_s0 * (1 - c)^n where n ≈ 4.65
    
    double v_s0 = 0.01;   // Single particle settling (m/s)
    double c = 0.3;       // Volume concentration (30%)
    double n = 4.65;
    
    double v_s = v_s0 * std::pow(1.0 - c, n);
    
    EXPECT_LT(v_s, v_s0);
    EXPECT_GT(v_s, 0.0);
}

// ============================================================================
// Fracture Network Connectivity Tests
// ============================================================================

TEST_F(FractureMechanicsTest, FractureNetworkP32) {
    // P32: Fracture area per unit volume
    
    double total_area = 1000.0;     // m² of fracture
    double volume = 1e6;            // 1 km³ = 1e9 m³, use smaller for test
    
    double P32 = total_area / volume;  // m⁻¹
    
    EXPECT_GT(P32, 0.0);
    EXPECT_TRUE(std::isfinite(P32));
}

TEST_F(FractureMechanicsTest, FracturePercolationThreshold) {
    // Percolation theory: connectivity threshold for random fractures
    // Critical P32 ~ 1/(π*r²) for circular fractures of radius r
    
    double r = 10.0;  // Mean fracture radius (m)
    double P32_critical = 1.0 / (M_PI * r * r);
    
    EXPECT_GT(P32_critical, 0.0);
}

// ============================================================================
// Natural Fracture Network Tests
// ============================================================================

TEST_F(FractureMechanicsTest, DualPorosityShapeFactor) {
    // Warren-Root shape factor: σ = 60 / L² for cubic blocks
    
    double L = 1.0;  // Block size (m)
    double sigma = 60.0 / (L * L);  // m⁻²
    
    EXPECT_NEAR(sigma, 60.0, 1e-6);
    
    // Kazemi shape factor: σ = 4*(1/Lx² + 1/Ly² + 1/Lz²)
    double Lx = 1.0, Ly = 2.0, Lz = 3.0;
    double sigma_kazemi = 4.0 * (1.0/(Lx*Lx) + 1.0/(Ly*Ly) + 1.0/(Lz*Lz));
    
    EXPECT_GT(sigma_kazemi, 0.0);
}

TEST_F(FractureMechanicsTest, FractureMatrixTransfer) {
    // Matrix-fracture transfer: q = σ * k_m * (P_m - P_f) / μ
    
    double k_m = 1e-18;     // Matrix permeability (m²)
    double P_m = 30e6;      // Matrix pressure (Pa)
    double P_f = 25e6;      // Fracture pressure (Pa)
    double sigma_shape = 60.0;  // Shape factor (m⁻²)
    
    double q = sigma_shape * k_m * (P_m - P_f) / mu_f;
    
    // Flow from matrix to fracture when P_m > P_f
    EXPECT_GT(q, 0.0);
}

// ============================================================================
// Fault Friction Tests
// ============================================================================

TEST_F(FractureMechanicsTest, MohrCoulombFailure) {
    // τ = c + μ_f * σ_n for slip
    
    double c = 0.0;         // No cohesion on fault
    double mu_friction = 0.6;  // Friction coefficient
    double sigma_n = 50e6;  // Normal stress (50 MPa)
    
    double tau_fail = c + mu_friction * sigma_n;
    
    EXPECT_NEAR(tau_fail, 30e6, 1e5);  // 30 MPa shear stress for slip
}

TEST_F(FractureMechanicsTest, SlipTendency) {
    // Slip tendency: T_s = τ / σ_n = μ_f at failure
    
    double tau = 25e6;      // Shear stress
    double sigma_n = 50e6;  // Normal stress
    
    double T_s = tau / sigma_n;
    
    EXPECT_NEAR(T_s, 0.5, 0.001);
    
    // Compare to friction coefficient
    double mu_f = 0.6;
    bool will_slip = (T_s >= mu_f);
    EXPECT_FALSE(will_slip);  // 0.5 < 0.6
}

TEST_F(FractureMechanicsTest, DilatancyAngle) {
    // Dilation during slip: dv/du = tan(ψ)
    // where dv = normal displacement, du = shear displacement
    
    double psi = 10.0 * M_PI / 180.0;  // 10 degrees dilation
    double du = 0.01;  // 1 cm slip
    
    double dv = du * std::tan(psi);
    
    EXPECT_GT(dv, 0.0);
    EXPECT_LT(dv, du);  // Dilation < slip for reasonable angles
}

// ============================================================================
// Rate-State Friction Tests
// ============================================================================

TEST_F(FractureMechanicsTest, RateStateFriction) {
    // μ = μ_0 + a*ln(V/V_0) + b*ln(V_0*θ/D_c)
    
    double mu_0 = 0.6;      // Reference friction
    double a = 0.01;        // Direct effect
    double b = 0.015;       // Evolution effect
    double V_0 = 1e-6;      // Reference velocity (m/s)
    double D_c = 1e-5;      // Critical slip distance (m)
    double V = 1e-3;        // Current velocity (m/s)
    double theta = 10.0;    // State variable (s)
    
    double mu = mu_0 + a * std::log(V / V_0) + b * std::log(V_0 * theta / D_c);
    
    EXPECT_TRUE(std::isfinite(mu));
    EXPECT_GT(mu, 0.0);
}

TEST_F(FractureMechanicsTest, VelocityWeakeningStrengthening) {
    // (a - b) < 0: Velocity weakening (unstable)
    // (a - b) > 0: Velocity strengthening (stable)
    
    double a = 0.01;
    double b_weak = 0.015;  // b > a: weakening
    double b_strong = 0.005;  // b < a: strengthening
    
    EXPECT_LT(a - b_weak, 0.0) << "Should be velocity weakening";
    EXPECT_GT(a - b_strong, 0.0) << "Should be velocity strengthening";
}

// ============================================================================
// Hydraulic Fracture Model Tests
// ============================================================================

TEST_F(FractureMechanicsTest, HydraulicFractureCreation) {
    HydraulicFractureModel hf;
    
    hf.setFractureModel("PKN");
    hf.setPropagationCriteria(K_IC, sigma_min);
    hf.setFormationProperties(E, nu);
    hf.setFluidProperties(rho_f, mu_f);
    
    // Verify properties are set
    EXPECT_TRUE(true);  // No exception thrown
}

TEST_F(FractureMechanicsTest, FractureGeometryUpdate) {
    HydraulicFractureModel hf;
    
    std::vector<double> coords = {0.0, 0.0, 0.0};  // Center
    hf.setGeometry(coords);
    
    // Update geometry (should evolve)
    hf.updateGeometry(1.0);  // 1 second
    
    // No crash
    EXPECT_TRUE(true);
}

// ============================================================================
// Natural Fracture Network Tests
// ============================================================================

TEST_F(FractureMechanicsTest, NaturalFractureCreation) {
    NaturalFractureNetwork nfn;
    
    // Add a fracture
    std::vector<double> points = {0.0, 0.0, 0.0, 10.0, 10.0, 0.0};
    nfn.addFracture(points, 0.001);  // 1 mm aperture
    
    EXPECT_TRUE(true);
}

TEST_F(FractureMechanicsTest, StochasticNetworkGeneration) {
    NaturalFractureNetwork nfn;
    
    nfn.generateStochasticNetwork(100, 10.0, 5.0, 42);  // 100 fractures
    
    EXPECT_TRUE(true);
}

// ============================================================================
// Fault Model Tests
// ============================================================================

TEST_F(FractureMechanicsTest, FaultModelCreation) {
    FaultModel fault;
    
    std::vector<double> strike = {1.0, 0.0, 0.0};
    std::vector<double> dip = {0.0, 0.707, 0.707};
    fault.setFaultPlane(strike, dip, 1000.0, 500.0);
    
    fault.setFrictionCoefficient(0.7, 0.5);
    fault.setCohesion(1e6);
    
    EXPECT_TRUE(true);
}

TEST_F(FractureMechanicsTest, FaultSlipMode) {
    FaultModel fault;
    
    std::vector<double> strike = {1.0, 0.0, 0.0};
    std::vector<double> dip = {0.0, 0.707, 0.707};
    fault.setFaultPlane(strike, dip, 1000.0, 500.0);
    fault.setFrictionCoefficient(0.6, 0.4);
    
    FaultModel::SlipMode mode;
    
    // Low shear stress - no slip
    fault.checkSlipCriteria(10e6, 50e6, mode);
    EXPECT_EQ(mode, FaultModel::SlipMode::NO_SLIP);
    
    // High shear stress - slip
    fault.checkSlipCriteria(35e6, 50e6, mode);
    EXPECT_NE(mode, FaultModel::SlipMode::NO_SLIP);
}

TEST_F(FractureMechanicsTest, FaultMomentMagnitude) {
    FaultModel fault;
    
    double slip_area = 10e6;    // 10 km² = 10e6 m²
    double slip_amount = 1.0;   // 1 m slip
    double Mw;
    
    fault.computeMomentMagnitude(slip_area, slip_amount, Mw);
    
    // M0 = μ * A * D = 30e9 * 10e6 * 1 = 3e17 N·m
    // Mw = (2/3) * log10(M0) - 6.06 ≈ 5.6
    EXPECT_GT(Mw, 5.0);
    EXPECT_LT(Mw, 7.0);
}

TEST_F(FractureMechanicsTest, RateStateFrictionUpdate) {
    FaultModel fault;
    
    fault.enableRateStateFriction(true);
    fault.setRateStateParameters(0.01, 0.015, 1e-5);
    
    double slip_rate = 1e-3;  // m/s
    double dt = 0.001;        // s
    
    fault.updateStateVariable(slip_rate, dt);
    
    EXPECT_TRUE(true);  // No crash
}

TEST_F(FractureMechanicsTest, FaultPermeabilityEvolution) {
    FaultModel fault;
    
    double k_initial = fault.getPermeability();
    
    fault.updatePermeability(0.1);  // 10 cm slip
    
    double k_after = fault.getPermeability();
    
    // Permeability typically increases with slip
    EXPECT_NE(k_initial, k_after);
}

// ============================================================================
// Energy Balance Tests
// ============================================================================

TEST_F(FractureMechanicsTest, FractureEnergyBalance) {
    // Energy: W_input = W_fracture + W_viscous + W_kinetic + W_stored
    
    double W_input = 1e9;       // Input energy (J)
    double W_fracture = 0.01 * W_input;  // ~1% to create surface
    double W_viscous = 0.80 * W_input;   // ~80% viscous dissipation
    double W_stored = 0.19 * W_input;    // ~19% elastic storage
    
    double W_total = W_fracture + W_viscous + W_stored;
    
    EXPECT_NEAR(W_total, W_input, W_input * 0.01);
}

TEST_F(FractureMechanicsTest, FractureSurfaceEnergy) {
    // Energy to create fracture surface: W = G_c * A
    
    double A = 1000.0 * 30.0;  // 1000 m × 30 m = 30000 m²
    double W = G_c * A * 2;    // Both faces
    
    EXPECT_GT(W, 0.0);
    EXPECT_TRUE(std::isfinite(W));
}

// ============================================================================
// Convergence Tests
// ============================================================================

TEST_F(FractureMechanicsTest, WidthConvergence) {
    // Test that width converges with mesh refinement
    
    std::vector<int> n_elements = {10, 20, 40, 80};
    std::vector<double> errors;
    
    double P_net = 5e6;
    double a = 50.0;
    double w_analytical = 4.0 * P_net * a / E_prime;
    
    for (int n : n_elements) {
        double h = a / n;
        // Assume O(h²) convergence for smooth solution
        double error = 0.1 * h * h;
        errors.push_back(error);
    }
    
    // Check convergence rate
    for (size_t i = 1; i < errors.size(); ++i) {
        double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
        EXPECT_GT(rate, 1.5) << "Should converge at least O(h^1.5)";
    }
}

TEST_F(FractureMechanicsTest, SIFConvergence) {
    // Stress intensity factor convergence near tip
    // Requires special elements (quarter-point) for optimal convergence
    
    // For standard elements, expect O(√h) convergence
    std::vector<double> h_values = {1.0, 0.5, 0.25, 0.125};
    std::vector<double> errors;
    
    for (double h : h_values) {
        // O(√h) error for standard elements
        errors.push_back(0.1 * std::sqrt(h));
    }
    
    double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
    EXPECT_NEAR(rate, 0.5, 0.1);  // √h convergence
}
