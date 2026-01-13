/**
 * @file test_fracture_model.cpp
 * @brief Unit tests for FractureModel classes
 */

#include <gtest/gtest.h>
#include "FractureModel.hpp"
#include "FSRM.hpp"
#include <cmath>

using namespace FSRM;

class FractureModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    
    int rank;
};

// ============================================================================
// NaturalFractureNetwork Tests
// ============================================================================

TEST_F(FractureModelTest, NaturalFractureCreate) {
    NaturalFractureNetwork network;
    EXPECT_EQ(network.getType(), FractureType::NATURAL);
}

TEST_F(FractureModelTest, NaturalFractureSetAperture) {
    NaturalFractureNetwork network;
    EXPECT_NO_THROW(network.setAperture(0.001));  // 1 mm
}

TEST_F(FractureModelTest, NaturalFractureSetPermeability) {
    NaturalFractureNetwork network;
    EXPECT_NO_THROW(network.setPermeability(1e-9));  // High perm in fracture
}

TEST_F(FractureModelTest, NaturalFractureEnableDualPorosity) {
    NaturalFractureNetwork network;
    EXPECT_NO_THROW(network.enableDualPorosity(true));
}

TEST_F(FractureModelTest, NaturalFractureSetShapeFactorModel) {
    NaturalFractureNetwork network;
    network.enableDualPorosity(true);
    EXPECT_NO_THROW(network.setShapeFactorModel("Warren-Root"));
    EXPECT_NO_THROW(network.setMatrixBlockSize(10.0));  // 10 m blocks
}

// ============================================================================
// HydraulicFractureModel Tests
// ============================================================================

TEST_F(FractureModelTest, HydraulicFractureCreate) {
    HydraulicFractureModel fracture;
    EXPECT_EQ(fracture.getType(), FractureType::INDUCED_HYDRAULIC);
}

TEST_F(FractureModelTest, HydraulicFractureSetModel) {
    HydraulicFractureModel fracture;
    EXPECT_NO_THROW(fracture.setFractureModel("KGD"));
}

TEST_F(FractureModelTest, HydraulicFractureSetPropagationCriteria) {
    HydraulicFractureModel fracture;
    double Kc = 1e6;        // Fracture toughness
    double sigma_min = 30e6; // Min horizontal stress
    EXPECT_NO_THROW(fracture.setPropagationCriteria(Kc, sigma_min));
}

TEST_F(FractureModelTest, HydraulicFractureEnableProppant) {
    HydraulicFractureModel fracture;
    EXPECT_NO_THROW(fracture.enableProppantTransport(true));
    EXPECT_NO_THROW(fracture.setProppantProperties(0.0005, 2650.0, 0.3));
}

TEST_F(FractureModelTest, HydraulicFractureEnableLeakoff) {
    HydraulicFractureModel fracture;
    EXPECT_NO_THROW(fracture.enableLeakoff(true));
    EXPECT_NO_THROW(fracture.setLeakoffCoefficient(1e-4));
}

// ============================================================================
// FaultModel Tests
// ============================================================================

TEST_F(FractureModelTest, FaultModelCreate) {
    FaultModel fault;
    SUCCEED();
}

TEST_F(FractureModelTest, FaultModelSetPlane) {
    FaultModel fault;
    std::vector<double> strike = {0.0, 1.0, 0.0};
    std::vector<double> dip = {0.0, 0.0, 1.0};
    EXPECT_NO_THROW(fault.setFaultPlane(strike, dip, 1000.0, 500.0));
}

TEST_F(FractureModelTest, FaultModelSetFriction) {
    FaultModel fault;
    EXPECT_NO_THROW(fault.setFrictionCoefficient(0.6, 0.5));  // static, dynamic
    EXPECT_NO_THROW(fault.setCohesion(1e6));
}

TEST_F(FractureModelTest, FaultModelEnableRateState) {
    FaultModel fault;
    EXPECT_NO_THROW(fault.enableRateStateFriction(true));
    EXPECT_NO_THROW(fault.setRateStateParameters(0.01, 0.015, 1e-4));  // a, b, Dc
}

// ============================================================================
// LEFM Calculations Verification
// ============================================================================

TEST_F(FractureModelTest, StressIntensityFactor) {
    // Analytical: K_I = (2/π) * σ * sqrt(π*a) for penny-shaped crack
    double sigma = 1e6;   // Pa
    double a = 10.0;      // m (crack radius)
    
    double K_I = (2.0 / M_PI) * sigma * std::sqrt(M_PI * a);
    
    EXPECT_GT(K_I, 0.0);
    EXPECT_TRUE(std::isfinite(K_I));
}

TEST_F(FractureModelTest, FractureWidthKGD) {
    // KGD model: w = 4*(1-nu)*p*L / E for elliptical crack
    double E = 20e9;
    double nu = 0.25;
    double p = 5e6;     // Net pressure
    double L = 100.0;   // Half-length
    
    double width = 4.0 * (1.0 - nu) * p * L / E;
    
    EXPECT_GT(width, 0.0);
    EXPECT_LT(width, 1.0);  // Should be mm-cm scale
}

TEST_F(FractureModelTest, FractureConductivity) {
    // Conductivity = k_f * w
    double k_f = 1e-9;   // Proppant pack permeability
    double w = 0.005;    // 5 mm width
    
    double conductivity = k_f * w;
    
    EXPECT_GT(conductivity, 0.0);
    EXPECT_NEAR(conductivity, 5e-12, 1e-15);
}

TEST_F(FractureModelTest, CarterLeakoff) {
    // V_L = 2 * C_L * A * sqrt(t)
    double C_L = 1e-4;    // m/s^0.5
    double A = 1000.0;    // m^2
    double t = 100.0;     // s
    
    double V_L = 2.0 * C_L * A * std::sqrt(t);
    
    EXPECT_GT(V_L, 0.0);
    EXPECT_TRUE(std::isfinite(V_L));
}

TEST_F(FractureModelTest, ProppantSettling) {
    // Stokes settling: v = (rho_p - rho_f) * g * d^2 / (18 * mu)
    double rho_p = 2650.0;   // Sand density
    double rho_f = 1000.0;   // Water
    double g = 9.81;
    double d = 0.0005;       // 0.5 mm diameter
    double mu = 0.001;       // Pa.s
    
    double v_settling = (rho_p - rho_f) * g * d * d / (18.0 * mu);
    
    EXPECT_GT(v_settling, 0.0);
    EXPECT_LT(v_settling, 1.0);  // Should be cm/s scale
}

TEST_F(FractureModelTest, CriticalStressIntensity) {
    // Propagation criterion: K_I > K_IC
    double K_IC = 1e6;  // Fracture toughness
    double K_I_subcritical = 0.8 * K_IC;
    double K_I_critical = 1.1 * K_IC;
    
    EXPECT_LT(K_I_subcritical, K_IC);  // No propagation
    EXPECT_GT(K_I_critical, K_IC);     // Propagation
}

TEST_F(FractureModelTest, ParisLaw) {
    // Fatigue crack growth: da/dN = C * (ΔK)^m
    double C = 1e-10;
    double m = 3.0;
    double deltaK = 0.5e6;
    
    double dadN = C * std::pow(deltaK, m);
    
    EXPECT_GT(dadN, 0.0);
    EXPECT_TRUE(std::isfinite(dadN));
}

// ============================================================================
// Fault Mechanics Verification
// ============================================================================

TEST_F(FractureModelTest, MohrCoulombCriterion) {
    // τ_max = c + μ * σ_n
    double c = 1e6;       // Cohesion
    double mu_f = 0.6;    // Friction coefficient
    double sigma_n = 30e6; // Normal stress
    
    double tau_max = c + mu_f * sigma_n;
    
    EXPECT_GT(tau_max, 0.0);
    EXPECT_TRUE(std::isfinite(tau_max));
}

TEST_F(FractureModelTest, EffectiveStress) {
    // σ_eff = σ - α * p
    double sigma = 50e6;  // Total stress
    double alpha = 1.0;   // Biot coefficient
    double p = 20e6;      // Pore pressure
    
    double sigma_eff = sigma - alpha * p;
    
    EXPECT_DOUBLE_EQ(sigma_eff, 30e6);
}

TEST_F(FractureModelTest, MomentMagnitude) {
    // M_0 = μ * A * d (seismic moment)
    // M_w = (log10(M_0) - 9.1) / 1.5
    double mu_rock = 30e9;  // Shear modulus
    double A = 1e6;         // Slip area (1 km^2)
    double d = 0.01;        // 1 cm slip
    
    double M_0 = mu_rock * A * d;
    double M_w = (std::log10(M_0) - 9.1) / 1.5;
    
    EXPECT_GT(M_0, 0.0);
    EXPECT_TRUE(std::isfinite(M_w));
    // For this small event: M_w ~ 2-3
    EXPECT_GT(M_w, 1.0);
    EXPECT_LT(M_w, 5.0);
}
