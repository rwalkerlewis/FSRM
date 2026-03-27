#include <gtest/gtest.h>
#include <mpi.h>
#include <petsc.h>

#include "domain/reservoir/FlashCalculation.hpp"
#include <cmath>
#include <vector>

using namespace FSRM;

class FlashCalculationTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

// Test 1: Pure water (zCO2=0) should be single-phase liquid
TEST_F(FlashCalculationTest, PureWaterSinglePhase) {
    auto result = FlashCalculation::flash(10e6, 350.0, 0.0, 0.035);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.Lw, 1.0, 0.01);
    EXPECT_NEAR(result.Vg, 0.0, 0.01);
}

// Test 2: Nearly pure CO2 should have large vapor fraction
TEST_F(FlashCalculationTest, HighCO2TwoPhase) {
    auto result = FlashCalculation::flash(10e6, 350.0, 0.5, 0.035);
    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.Vg, 0.0);
}

// Test 3: Flash converges
TEST_F(FlashCalculationTest, FlashConverges) {
    auto result = FlashCalculation::flash(15e6, 350.0, 0.3, 0.035);
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.iterations, 30);
}

// Test 4: Phase fractions sum to 1
TEST_F(FlashCalculationTest, PhaseFractionsSumToOne) {
    auto result = FlashCalculation::flash(10e6, 330.0, 0.4, 0.02);
    EXPECT_NEAR(result.Lw + result.Vg, 1.0, 1e-10);
}

// Test 5: Rachford-Rice function at solution is near zero
TEST_F(FlashCalculationTest, RachfordRiceAtSolution) {
    std::vector<double> z = {0.7, 0.3};
    std::vector<double> K = {0.001, 50.0};
    double V = FlashCalculation::solveRachfordRice(z, K);
    double f = FlashCalculation::rrFunction(z, K, V);
    EXPECT_NEAR(f, 0.0, 1e-10);
}

// Test 6: Wilson K-values are reasonable
TEST_F(FlashCalculationTest, WilsonKValues) {
    auto K = FlashCalculation::wilsonKValues(10e6, 350.0);
    EXPECT_EQ(K.size(), 2u);
    EXPECT_GT(K[1], 1.0);  // CO2 should prefer gas phase (K > 1)
}

// Test 7: Stability check for pure water
TEST_F(FlashCalculationTest, StabilityPureWater) {
    std::vector<double> z = {1.0, 0.0};
    auto K = FlashCalculation::wilsonKValues(10e6, 350.0);
    bool stable = FlashCalculation::isStable(10e6, 350.0, z, K);
    EXPECT_TRUE(stable);
}

// Test 8: Densities are physical
TEST_F(FlashCalculationTest, DensitiesPhysical) {
    auto result = FlashCalculation::flash(10e6, 350.0, 0.4, 0.035);
    if (result.Vg > 0.01) {
        EXPECT_GT(result.rho_gas, 1.0);
        EXPECT_LT(result.rho_gas, 2000.0);
    }
    if (result.Lw > 0.01) {
        EXPECT_GT(result.rho_aq, 500.0);
        EXPECT_LT(result.rho_aq, 2000.0);
    }
}
