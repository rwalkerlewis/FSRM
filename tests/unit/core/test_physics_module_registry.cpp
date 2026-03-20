/**
 * @file test_physics_module_registry.cpp
 * @brief Unit tests for PhysicsModuleRegistry
 */

#include <gtest/gtest.h>
#include "core/PhysicsModuleRegistry.hpp"
#include "core/PhysicsModuleInterface.hpp"
#include "core/ConfigReader.hpp"
#include <mpi.h>
#include <petscerror.h>
#include <algorithm>
#include <vector>

using namespace FSRM;

namespace {

class MockModule : public PhysicsModuleInterface {
public:
    std::string name() const override { return "mock_module"; }
    std::string description() const override { return "Mock for testing"; }
    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::FLUID_FLOW};
    }
    std::vector<PhysicsType> requiredPhysics() const override { return {}; }

    PetscErrorCode configure(const ConfigReader&) override { return 0; }
    PetscErrorCode validateConfig() const override { return 0; }
    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override { return {}; }
    PetscErrorCode initialize(DM, Vec) override { return 0; }

    bool isEnabled(const SimulationConfig& config) const override {
        return std::find(config.enabled_modules.begin(), config.enabled_modules.end(),
                         "mock_module") != config.enabled_modules.end();
    }
};

class MockModuleWithDeps : public PhysicsModuleInterface {
public:
    std::string name() const override { return "mock_with_deps"; }
    std::string description() const override { return "Requires geomechanics"; }
    std::vector<PhysicsType> providedPhysics() const override {
        return {PhysicsType::FLUID_FLOW};
    }
    std::vector<PhysicsType> requiredPhysics() const override {
        return {PhysicsType::GEOMECHANICS};
    }

    PetscErrorCode configure(const ConfigReader&) override { return 0; }
    PetscErrorCode validateConfig() const override { return 0; }
    std::vector<std::shared_ptr<PhysicsKernel>> createKernels() override { return {}; }
    PetscErrorCode initialize(DM, Vec) override { return 0; }

    bool isEnabled(const SimulationConfig&) const override { return true; }
};

} // namespace

class PhysicsModuleRegistryTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        PhysicsModuleRegistry::instance().clear();
    }

    void TearDown() override { PhysicsModuleRegistry::instance().clear(); }

    int rank = 0;
};

TEST_F(PhysicsModuleRegistryTest, InitiallyEmptyAfterClear) {
    auto& reg = PhysicsModuleRegistry::instance();
    EXPECT_TRUE(reg.availableModules().empty());
}

TEST_F(PhysicsModuleRegistryTest, RegisterAndAvailableModules) {
    auto& reg = PhysicsModuleRegistry::instance();
    reg.registerModule("mock_module", []() {
        return std::make_unique<MockModule>();
    });
    auto names = reg.availableModules();
    ASSERT_EQ(names.size(), 1u);
    EXPECT_EQ(names[0], "mock_module");
}

TEST_F(PhysicsModuleRegistryTest, HasModuleAndCreateModule) {
    auto& reg = PhysicsModuleRegistry::instance();
    reg.registerModule("mock_module", []() {
        return std::make_unique<MockModule>();
    });
    EXPECT_TRUE(reg.hasModule("mock_module"));
    EXPECT_FALSE(reg.hasModule("unknown"));

    auto m = reg.createModule("mock_module");
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->name(), "mock_module");

    EXPECT_EQ(reg.createModule("not_registered"), nullptr);
}

TEST_F(PhysicsModuleRegistryTest, CreateEnabledModulesReturnsMock) {
    auto& reg = PhysicsModuleRegistry::instance();
    reg.registerModule("mock_module", []() {
        return std::make_unique<MockModule>();
    });
    SimulationConfig cfg;
    cfg.enabled_modules = {"mock_module"};
    auto mods = reg.createEnabledModules(cfg);
    ASSERT_EQ(mods.size(), 1u);
    ASSERT_NE(mods[0], nullptr);
    EXPECT_EQ(mods[0]->name(), "mock_module");
}

TEST_F(PhysicsModuleRegistryTest, DuplicateRegistrationDoesNotCrash) {
    auto& reg = PhysicsModuleRegistry::instance();
    reg.registerModule("mock_module", []() {
        return std::make_unique<MockModule>();
    });
    reg.registerModule("mock_module", []() {
        return std::make_unique<MockModule>();
    });
    EXPECT_EQ(reg.availableModules().size(), 1u);
    auto m = reg.createModule("mock_module");
    ASSERT_NE(m, nullptr);
}

TEST_F(PhysicsModuleRegistryTest, ValidateDependenciesPassesWithNoRequirements) {
    auto& reg = PhysicsModuleRegistry::instance();
    reg.registerModule("mock_module", []() {
        return std::make_unique<MockModule>();
    });
    auto m = reg.createModule("mock_module");
    ASSERT_NE(m, nullptr);
    std::vector<PhysicsModuleInterface*> ptrs = {m.get()};
    PetscErrorCode ierr = reg.validateDependencies(ptrs);
    EXPECT_EQ(ierr, 0);
}

TEST_F(PhysicsModuleRegistryTest, ValidateDependenciesFailsWhenRequiredPhysicsMissing) {
    auto& reg = PhysicsModuleRegistry::instance();
    reg.registerModule("mock_with_deps", []() {
        return std::make_unique<MockModuleWithDeps>();
    });
    auto m = reg.createModule("mock_with_deps");
    ASSERT_NE(m, nullptr);
    std::vector<PhysicsModuleInterface*> ptrs = {m.get()};
    PetscErrorCode ierr = reg.validateDependencies(ptrs);
    EXPECT_EQ(ierr, PETSC_ERR_ARG_INCOMP);
}
