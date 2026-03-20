/**
 * @file test_well_model.cpp
 * @brief Unit tests for WellModel (smoke and basic API)
 */

#include "domain/reservoir/WellModel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace FSRM;

class WellModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

TEST_F(WellModelTest, ConstructWithWellTypeProducer) {
    WellModel well("W1", WellType::PRODUCER);
    EXPECT_EQ(well.getName(), "W1");
    EXPECT_EQ(well.getType(), WellType::PRODUCER);
}

TEST_F(WellModelTest, ConstructWithWellTypeInjector) {
    WellModel well("W2", WellType::INJECTOR);
    EXPECT_EQ(well.getType(), WellType::INJECTOR);
}

TEST_F(WellModelTest, ProductionWellTypeIsProducer) {
    ProductionWell well("PROD1");
    EXPECT_EQ(well.getType(), WellType::PRODUCER);
}

TEST_F(WellModelTest, InjectionWellTypeIsInjector) {
    InjectionWell well("INJ1");
    EXPECT_EQ(well.getType(), WellType::INJECTOR);
}

TEST_F(WellModelTest, ReferenceDepthAndWellboreRadii) {
    ProductionWell well("PROD1");
    well.setReferenceDepth(2500.0);
    well.setWellboreGeometry(0.108, 0.15);
    well.setWellboreRoughness(1e-5);
    well.enableWellboreModel(true);
    SUCCEED();
}

TEST_F(WellModelTest, AddCompletionSmoke) {
    ProductionWell well("PROD1");
    WellCompletion comp{};
    comp.i = 3;
    comp.j = 4;
    comp.k = 2;
    comp.well_index = 1e-12;
    comp.diameter = 0.2;
    comp.skin_factor = 0.0;
    comp.is_open = true;
    comp.perforation_length = 10.0;
    EXPECT_NO_THROW(well.addCompletion(comp));
    ASSERT_EQ(well.getCompletions().size(), 1u);
    EXPECT_DOUBLE_EQ(well.getCompletions()[0].diameter, 0.2);
}

TEST_F(WellModelTest, ComputeRateFinite) {
    ProductionWell well("PROD1");
    WellCompletion comp{};
    comp.i = 5;
    comp.j = 5;
    comp.k = 0;
    comp.well_index = 1e-12;
    comp.is_open = true;
    well.addCompletion(comp);
    double rate = well.computeRate(20e6, 10e6, 0);
    EXPECT_TRUE(std::isfinite(rate));
}

TEST_F(WellModelTest, WellIndexFinite) {
    ProductionWell well("PROD1");
    double WI = well.computeWellIndex(0, 100e-15, 100e-15, 10e-15, 100.0, 100.0, 10.0);
    EXPECT_TRUE(std::isfinite(WI));
}
