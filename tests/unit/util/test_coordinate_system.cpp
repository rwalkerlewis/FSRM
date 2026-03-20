/**
 * @file test_coordinate_system.cpp
 * @brief Smoke tests for CoordinateSystemManager
 */

#include "util/CoordinateSystem.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>

using namespace FSRM;

class CoordinateSystemTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

TEST_F(CoordinateSystemTest, ManagerConstruction) {
    CoordinateSystemManager mgr;
    EXPECT_FALSE(mgr.isConfigured());
}

TEST_F(CoordinateSystemTest, IdentityTransformSameEpsg) {
    CoordinateSystemManager mgr;
    ASSERT_TRUE(mgr.setInputCRS("EPSG:32610"));
    ASSERT_TRUE(mgr.setModelCRS("EPSG:32610"));
    ASSERT_TRUE(mgr.isConfigured());

    GeoPoint p{100.0, 200.0, 10.0};
    GeoPoint q = mgr.toModelCoords(p);
    EXPECT_DOUBLE_EQ(q.x, p.x);
    EXPECT_DOUBLE_EQ(q.y, p.y);
    EXPECT_DOUBLE_EQ(q.z, p.z);
}

TEST_F(CoordinateSystemTest, LocalOriginOffset) {
    CoordinateSystemManager mgr;
    mgr.setInputCRS("EPSG:32610");
    mgr.setModelCRS("EPSG:32610");
    mgr.setLocalOrigin(GeoPoint(10.0, 20.0, 0.0));
    GeoPoint p{110.0, 220.0, 5.0};
    GeoPoint q = mgr.toModelCoords(p);
    EXPECT_NEAR(q.x, 100.0, 1e-9);
    EXPECT_NEAR(q.y, 200.0, 1e-9);
    EXPECT_DOUBLE_EQ(q.z, 5.0);
}
