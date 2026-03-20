/**
 * @file test_simulator_init.cpp
 * @brief Functional tests for Simulator construction and basic grid setup
 */

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "core/ConfigReader.hpp"
#include <gtest/gtest.h>

using namespace FSRM;

class SimulatorInitTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(SimulatorInitTest, ConstructFromPetscCommWorld) {
    Simulator sim(PETSC_COMM_WORLD);
    (void)sim;
    SUCCEED();
}

TEST_F(SimulatorInitTest, InitializeDefaultSimulationConfig) {
    Simulator sim(PETSC_COMM_WORLD);
    SimulationConfig cfg{};
    PetscErrorCode ierr = sim.initialize(cfg);
    EXPECT_EQ(ierr, 0);
}

TEST_F(SimulatorInitTest, SetupDMWithDefaultGridConfig) {
    Simulator sim(PETSC_COMM_WORLD);
    SimulationConfig cfg{};
    PetscErrorCode ierr = sim.initialize(cfg);
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupDM();
    EXPECT_EQ(ierr, 0);
}

TEST_F(SimulatorInitTest, ConfigReaderUsableWithSimulationDefaults) {
    ConfigReader reader;
    SimulationConfig parsed{};
    EXPECT_FALSE(reader.parseSimulationConfig(parsed));
}
