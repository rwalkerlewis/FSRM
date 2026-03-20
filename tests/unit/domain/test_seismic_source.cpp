/**
 * @file test_seismic_source.cpp
 * @brief Google Test smoke tests for seismic source utilities
 */

#include "domain/seismic/SeismicSource.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace FSRM;

class SeismicSourceTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

TEST_F(SeismicSourceTest, MomentTensorDoubleCoupleTrace) {
    const double M0 = 1e18;
    MomentTensor M = MomentTensor::doubleCouple(0.0, M_PI / 2.0, 0.0, M0);
    double trace = M.Mxx + M.Myy + M.Mzz;
    EXPECT_NEAR(trace, 0.0, 1e-6 * M0);
    EXPECT_GT(M.scalarMoment(), 0.0);
}

TEST_F(SeismicSourceTest, MomentTensorExplosionIsotropic) {
    const double M0 = 1e15;
    MomentTensor M = MomentTensor::explosion(M0);
    EXPECT_NEAR(M.Mxx, M0, 1e-9);
    EXPECT_NEAR(M.Myy, M0, 1e-9);
    EXPECT_NEAR(M.Mzz, M0, 1e-9);
    EXPECT_NEAR(M.Mxy, 0.0, 1e-9);
    EXPECT_NEAR(M.Mxz, 0.0, 1e-9);
    EXPECT_NEAR(M.Myz, 0.0, 1e-9);
}

TEST_F(SeismicSourceTest, SourceTimeFunctionGaussianReasonable) {
    SourceTimeFunctionEvaluator stf(SourceTimeFunction::GAUSSIAN);
    stf.setDuration(0.5);
    stf.setOnsetTime(0.0);
    double r0 = stf.momentRate(0.0);
    double r1 = stf.momentRate(10.0);
    EXPECT_TRUE(std::isfinite(r0));
    EXPECT_TRUE(std::isfinite(r1));
    EXPECT_GE(r0, 0.0);
    EXPECT_GE(r1, 0.0);
}

TEST_F(SeismicSourceTest, PointSourceInstantiation) {
    PointSource src;
    src.setLocation(1.0, 2.0, -500.0);
    double x = 0, y = 0, z = 0;
    src.getLocation(x, y, z);
    EXPECT_DOUBLE_EQ(x, 1.0);
    EXPECT_DOUBLE_EQ(y, 2.0);
    EXPECT_DOUBLE_EQ(z, -500.0);
}

TEST_F(SeismicSourceTest, SourceReceiverManagerSmoke) {
    SourceReceiverManager mgr;
    EXPECT_EQ(mgr.getNumPointSources(), 0u);
    EXPECT_EQ(mgr.getNumReceivers(), 0u);
}
