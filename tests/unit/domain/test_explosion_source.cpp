/**
 * @file test_explosion_source.cpp
 * @brief Smoke tests for ExplosionSourceKernel
 */

#include "physics/ExplosionDamageKernels.hpp"
#include "core/FSRM.hpp"
#include <array>
#include <cmath>
#include <gtest/gtest.h>

using namespace FSRM;

class ExplosionSourceKernelTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

TEST_F(ExplosionSourceKernelTest, Construction) {
    ExplosionSourceKernel kernel;
    EXPECT_EQ(kernel.getType(), PhysicsType::EXPLOSION_SOURCE);
}

TEST_F(ExplosionSourceKernelTest, FieldLayout) {
    ExplosionSourceKernel kernel;
    EXPECT_EQ(kernel.getNumFields(), 1);
    EXPECT_EQ(kernel.getNumComponents(0), 6);
}

TEST_F(ExplosionSourceKernelTest, SetExplosionParametersDerived) {
    ExplosionSourceKernel kernel;
    kernel.setMediumProperties(2700.0, 5000.0, 3000.0);
    kernel.setExplosionParameters(5.0, 400.0);
    EXPECT_GT(kernel.getScalarMoment(), 0.0);
    EXPECT_GT(kernel.getCornerFrequency(), 0.0);
    EXPECT_TRUE(std::isfinite(kernel.getMagnitude()));
}

TEST_F(ExplosionSourceKernelTest, SourceTimeFunctionBounded) {
    ExplosionSourceKernel kernel;
    kernel.setExplosionParameters(1.0, 200.0);
    EXPECT_NEAR(kernel.sourceTimeFunction(0.0), 1.0, 1e-9);
    EXPECT_GE(kernel.sourceTimeFunction(1.0), 0.0);
    EXPECT_LE(kernel.sourceTimeFunction(1.0), 1.0);
}

TEST_F(ExplosionSourceKernelTest, MomentTensorVoigt) {
    ExplosionSourceKernel kernel;
    kernel.setExplosionParameters(2.0, 300.0);
    std::array<double, 6> M{};
    kernel.getMomentTensor(0.0, M.data());
    EXPECT_TRUE(std::isfinite(M[0]));
    EXPECT_TRUE(std::isfinite(M[1]));
    EXPECT_TRUE(std::isfinite(M[2]));
}

// DPRK 2017 validation test using NuclearSourceParameters directly
#include "domain/explosion/ExplosionImpactPhysics.hpp"

TEST_F(ExplosionSourceKernelTest, DPRK2017MagnitudeCheck) {
    NuclearSourceParameters params;
    params.yield_kt = 250.0;
    params.depth_of_burial = 800.0;
    double mb = params.body_wave_magnitude();
    // Murphy (1981): mb = 4.45 + 0.75*log10(250) = 4.45 + 0.75*2.398 = 6.25
    // Observed DPRK 2017: mb ~ 6.3
    EXPECT_NEAR(mb, 6.25, 0.15);
}
