/**
 * @file test_ml_registry.cpp
 * @brief Smoke tests for MLModelRegistry
 */

#include "ml/MLModelRegistry.hpp"
#include "core/FSRM.hpp"
#include <algorithm>
#include <gtest/gtest.h>

using namespace FSRM;
using namespace FSRM::ML;

class MLModelRegistryTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
};

TEST_F(MLModelRegistryTest, SingletonConstructs) {
    MLModelRegistry& r = MLModelRegistry::getInstance();
    (void)r;
    SUCCEED();
}

TEST_F(MLModelRegistryTest, RegisterAndDiscover) {
    MLModelRegistry& r = MLModelRegistry::getInstance();
    const std::string name = "__gtest_ml_dummy__";
    r.registerModel<int>(name, ModelCategory::FNO, std::make_shared<int>(42));
    EXPECT_TRUE(r.hasModel(name));
    auto names = r.getAllModelNames();
    EXPECT_NE(std::find(names.begin(), names.end(), name), names.end());
    auto fno = r.getModelsByCategory(ModelCategory::FNO);
    EXPECT_NE(std::find(fno.begin(), fno.end(), name), fno.end());
    int* m = r.getModel<int>(name);
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(*m, 42);
    EXPECT_EQ(r.getModelStatus(name), ModelStatus::INITIALIZED);
}

TEST_F(MLModelRegistryTest, SelectBestModelEmptyCategory) {
    MLModelRegistry& r = MLModelRegistry::getInstance();
    std::string pick = r.selectBestModel(ModelCategory::MF_COKRIGING);
    EXPECT_TRUE(pick.empty());
}
