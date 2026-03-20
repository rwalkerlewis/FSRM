/**
 * @file test_coupled_physics.cpp
 * @brief Fluid flow + geomechanics setup smoke test
 */

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>

using namespace FSRM;

class CoupledPhysicsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if (rank == 0) {
            std::ofstream cfg(config_path_);
            cfg << "[SIMULATION]\n";
            cfg << "fluid_model = SINGLE_COMPONENT\n";
            cfg << "enable_geomechanics = true\n";
            cfg << "solid_model = ELASTIC\n";
            cfg << "max_timesteps = 1\n";
            cfg << "dt_initial = 0.01\n";
            cfg << "\n[GRID]\n";
            cfg << "nx = 3\n";
            cfg << "ny = 3\n";
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
    static constexpr const char* config_path_ = "coupled_physics_smoke.ini";
};

TEST_F(CoupledPhysicsTest, FluidAndGeomechanicsSetupDoesNotCrash) {
    Simulator sim(PETSC_COMM_WORLD);
    PetscErrorCode ierr = sim.initializeFromConfigFile(config_path_);
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupDM();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupFields();
    ASSERT_EQ(ierr, 0);
    ierr = sim.setupPhysics();
    EXPECT_EQ(ierr, 0);
}
