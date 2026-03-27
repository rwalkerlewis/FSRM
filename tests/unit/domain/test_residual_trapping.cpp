#include <gtest/gtest.h>
#include <mpi.h>
#include <petsc.h>

#include "domain/reservoir/ResidualTrapping.hpp"
#include <cmath>

using namespace FSRM;

class ResidualTrappingTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

// Test 1: Land model basic formula
TEST_F(ResidualTrappingTest, LandModelBasic) {
    double C = 2.0;
    double Sg_max = 0.5;
    double Sgr = ResidualTrapping::landTrappedSaturation(Sg_max, C);
    // Sgr = 0.5 / (1 + 2*0.5) = 0.5/2.0 = 0.25
    EXPECT_NEAR(Sgr, 0.25, 1e-10);
}

// Test 2: Land model limit for large Sgi
TEST_F(ResidualTrappingTest, LandModelLargeLimit) {
    double C = 2.0;
    double Sgr = ResidualTrapping::landTrappedSaturation(100.0, C);
    // Sgr → 1/C = 0.5 for large Sgi
    EXPECT_NEAR(Sgr, 1.0/C, 0.01);
}

// Test 3: Land model zero input
TEST_F(ResidualTrappingTest, LandModelZero) {
    double Sgr = ResidualTrapping::landTrappedSaturation(0.0, 2.0);
    EXPECT_NEAR(Sgr, 0.0, 1e-15);
}

// Test 4: History tracking
TEST_F(ResidualTrappingTest, HistoryTracking) {
    ResidualTrapping trap;
    trap.setNumCells(10);
    trap.setLandCoefficient(2.0);

    trap.updateHistory(0, 0.3);
    EXPECT_NEAR(trap.getMaxHistoricalSg(0), 0.3, 1e-10);

    trap.updateHistory(0, 0.5);
    EXPECT_NEAR(trap.getMaxHistoricalSg(0), 0.5, 1e-10);

    trap.updateHistory(0, 0.2);  // Should NOT decrease
    EXPECT_NEAR(trap.getMaxHistoricalSg(0), 0.5, 1e-10);
}

// Test 5: Trapped saturation from history
TEST_F(ResidualTrappingTest, TrappedSgFromHistory) {
    ResidualTrapping trap;
    trap.setNumCells(5);
    trap.setLandCoefficient(2.0);

    trap.updateHistory(0, 0.5);
    double Sgr = trap.getTrappedSg(0);
    EXPECT_NEAR(Sgr, 0.25, 1e-10);  // 0.5/(1+2*0.5)
}

// Test 6: Drainage rel perm is zero at critical
TEST_F(ResidualTrappingTest, DrainageRelPermAtCritical) {
    ResidualTrapping trap;
    trap.setDrainageRelPerm(0.05, 1.0, 2.0);
    double kr = trap.krg_drainage(0.05);
    EXPECT_NEAR(kr, 0.0, 1e-10);
}

// Test 7: Drainage rel perm at high Sg
TEST_F(ResidualTrappingTest, DrainageRelPermHigh) {
    ResidualTrapping trap;
    trap.setDrainageRelPerm(0.05, 1.0, 2.0);
    double kr = trap.krg_drainage(0.8);
    EXPECT_GT(kr, 0.0);
    EXPECT_LE(kr, 1.0);
}

// Test 8: Imbibition kr bounded by drainage kr
TEST_F(ResidualTrappingTest, ImbibitionBounded) {
    ResidualTrapping trap;
    trap.setNumCells(5);
    trap.setLandCoefficient(2.0);
    trap.setDrainageRelPerm(0.05, 1.0, 2.0);

    trap.updateHistory(0, 0.6);
    double kr_imb = trap.krg_imbibition(0.4, 0);
    double kr_drain = trap.krg_drainage(0.4);
    EXPECT_LE(kr_imb, kr_drain + 1e-10);
}
