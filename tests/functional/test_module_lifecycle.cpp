/**
 * @file test_module_lifecycle.cpp
 * @brief Physics module lifecycle (configure through finalize)
 */

#include "core/PhysicsModuleRegistry.hpp"
#include "core/PhysicsModuleInterface.hpp"
#include "core/FSRM.hpp"
#include "core/ConfigReader.hpp"
#include <gtest/gtest.h>
#include <memory>
#include <vector>

using namespace FSRM;

namespace {

class TestLifecycleModule : public PhysicsModuleInterface {
public:
    mutable bool configured_ = false;
    mutable bool validated_ = false;
    mutable bool initialized_ = false;
    mutable bool finalized_ = false;

    std::string name() const override { return "test.lifecycle"; }
    std::string description() const override { return "Test lifecycle tracking module"; }

    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::FLUID_FLOW};
    }

    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader& /*config*/) override {
        configured_ = true;
        return 0;
    }

    PetscErrorCode validateConfig() const override {
        validated_ = true;
        return 0;
    }

    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override {
        return {};
    }

    PetscErrorCode initialize(DM dm, Vec solution) override {
        (void)dm;
        (void)solution;
        initialized_ = true;
        return 0;
    }

    PetscErrorCode finalize() override {
        finalized_ = true;
        return 0;
    }

    bool isEnabled(const SimulationConfig& /*config*/) const override { return true; }
};

} // namespace

class ModuleLifecycleTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(ModuleLifecycleTest, FullLifecycleOrder) {
    TestLifecycleModule mod;
    ConfigReader reader;
    DM dm = nullptr;
    Vec sol = nullptr;

    EXPECT_EQ(mod.configure(reader), 0);
    EXPECT_EQ(mod.validateConfig(), 0);
    auto kernels = mod.createKernels();
    EXPECT_TRUE(kernels.empty());

    EXPECT_EQ(mod.initialize(dm, sol), 0);
    EXPECT_EQ(mod.finalize(), 0);

    EXPECT_TRUE(mod.configured_);
    EXPECT_TRUE(mod.validated_);
    EXPECT_TRUE(mod.initialized_);
    EXPECT_TRUE(mod.finalized_);
}

TEST_F(ModuleLifecycleTest, CreateKernelsReturnsEmptyVector) {
    TestLifecycleModule mod;
    ConfigReader reader;
    ASSERT_EQ(mod.configure(reader), 0);
    auto kernels = mod.createKernels();
    EXPECT_EQ(kernels.size(), 0u);
    (void)rank;
}

TEST_F(ModuleLifecycleTest, PhysicsModuleRegistrySingleton) {
    PhysicsModuleRegistry& reg = PhysicsModuleRegistry::instance();
    (void)reg;
    SUCCEED();
}
