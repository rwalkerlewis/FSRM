/**
 * @file test_full_simulation.cpp
 * @brief Integration: minimal config on disk, full setup pipeline
 */

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include "core/ConfigReader.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>

using namespace FSRM;

class FullSimulationTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            std::ofstream cfg(config_path_);
            cfg << "[SIMULATION]\n";
            cfg << "fluid_model = SINGLE_COMPONENT\n";
            cfg << "max_timesteps = 2\n";
            cfg << "dt_initial = 0.01\n";
            cfg << "end_time = 0.02\n";
            cfg << "\n[GRID]\n";
            cfg << "nx = 4\n";
            cfg << "ny = 4\n";
            cfg << "nz = 2\n";
            cfg << "Lx = 1.0\n";
            cfg << "Ly = 1.0\n";
            cfg << "Lz = 1.0\n";
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }

    void TearDown() override {
        MPI_Barrier(PETSC_COMM_WORLD);
        if (rank == 0) {
            std::remove(config_path_);
        }
    }

    int rank = 0;
    static constexpr const char* config_path_ = "integration_min.ini";
};

TEST_F(FullSimulationTest, LoadConfigAndCompleteSetupStages) {
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(config_path_));
    GridConfig grid{};
    ASSERT_TRUE(reader.parseGridConfig(grid));
    EXPECT_EQ(grid.nx, 4);

    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_path_);
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupDM();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupFields();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupPhysics();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupTimeStepper();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupSolvers();
    EXPECT_EQ(ierr, 0);
}
