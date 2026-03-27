/**
 * @file test_caprock_integrity.cpp
 * @brief Unit tests for CaprockIntegrity
 */
#include <gtest/gtest.h>
#include "domain/geomechanics/CaprockIntegrity.hpp"
#include <cmath>
#include <mpi.h>

using namespace FSRM;

class CaprockIntegrityTest : public ::testing::Test {};

// Mohr-Coulomb FOS > 1 for stable conditions
TEST_F(CaprockIntegrityTest, StableConditionsFOS) {
    // Low shear stress, high normal stress
    double fos = CaprockIntegrity::mohrCoulombFOS(50e6, 5e6, 5e6, 30.0);
    EXPECT_GT(fos, 1.0);
}

// Mohr-Coulomb FOS < 1 for failure
TEST_F(CaprockIntegrityTest, FailureConditionsFOS) {
    // High shear stress, low normal stress
    double fos = CaprockIntegrity::mohrCoulombFOS(1e6, 20e6, 1e6, 30.0);
    EXPECT_LT(fos, 1.0);
}

// Permeability below failure
TEST_F(CaprockIntegrityTest, PermeabilityBelowFailure) {
    double k = CaprockIntegrity::caprockPermeability(
        1e-20, 50e6, 50e6, 1e-7, 1000.0, false);
    EXPECT_NEAR(k, 1e-20, 1e-22);  // sigma_eff = sigma_eff_ref → k = k0
}

// Permeability at failure (enhanced)
TEST_F(CaprockIntegrityTest, PermeabilityAtFailure) {
    double k = CaprockIntegrity::caprockPermeability(
        1e-20, 50e6, 50e6, 1e-7, 1000.0, true);
    EXPECT_NEAR(k, 1e-17, 1e-19);  // k0 * 1000
}

// Coulomb failure stress change
TEST_F(CaprockIntegrityTest, DeltaCFS) {
    // ΔCFS = Δτ - μ*(Δσn - ΔP)
    double dcfs = CaprockIntegrity::deltaCFS(5e6, 10e6, 0.6, 8e6);
    // = 5e6 - 0.6*(10e6 - 8e6) = 5e6 - 1.2e6 = 3.8e6
    EXPECT_NEAR(dcfs, 3.8e6, 1e3);
}

// Wells-Coppersmith magnitude
TEST_F(CaprockIntegrityTest, MagnitudeEstimate) {
    // M = (log10(A) - 4.07) / 0.98 for A in m²
    // A = 1e6 m²: M = (6 - 4.07)/0.98 ≈ 1.97
    double M = CaprockIntegrity::estimateMagnitude(1e6);
    EXPECT_NEAR(M, 1.97, 0.05);
}

// Full evaluation
TEST_F(CaprockIntegrityTest, FullEvaluation) {
    CaprockIntegrity cap;
    cap.setMohrCoulombParams(5e6, 30.0);
    cap.setCaprockProperties(1e-20, 1e-7, 1000.0);

    // Hydrostatic stress + some deviatoric
    std::array<double,6> stress = {-50e6, -50e6, -60e6, 5e6, 0, 0};
    double pore_pressure = 30e6;
    std::array<double,3> fault_normal = {0, 0, 1};

    auto result = cap.evaluate(stress, pore_pressure, 5e6, 30.0, fault_normal);
    EXPECT_TRUE(std::isfinite(result.mohr_coulomb_safety_factor));
    EXPECT_TRUE(std::isfinite(result.max_delta_CFS));
    EXPECT_TRUE(std::isfinite(result.estimated_magnitude));
}
