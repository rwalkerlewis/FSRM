/**
 * @file test_solver_convergence.cpp
 * @brief Smoke test: small single-phase grid initialization
 */

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <string>

using namespace FSRM;

class SolverConvergenceTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            std::ofstream cfg(config_path_);
            cfg << "[SIMULATION]\n";
            cfg << "fluid_model = SINGLE_COMPONENT\n";
            cfg << "max_timesteps = 2\n";
            cfg << "dt_initial = 0.01\n";
            cfg << "\n[GRID]\n";
            cfg << "nx = 4\n";
            cfg << "ny = 4\n";
            cfg << "nz = 1\n";
            cfg << "Lx = 1.0\n";
            cfg << "Ly = 1.0\n";
            cfg << "Lz = 1.0\n";
            cfg << "mesh_type = CARTESIAN\n";
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
    static constexpr const char* config_path_ = "test_solver_conv_small.ini";
};

TEST_F(SolverConvergenceTest, InitializeSmallSinglePhaseGrid) {
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_path_);
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupDM();
    EXPECT_EQ(ierr, 0);
}
