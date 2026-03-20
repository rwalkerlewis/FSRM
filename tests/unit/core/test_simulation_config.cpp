/**
 * @file test_simulation_config.cpp
 * @brief Unit tests for SimulationConfig defaults and module loading
 */

#include <gtest/gtest.h>
#include "core/FSRM.hpp"
#include "core/ConfigReader.hpp"
#include <fstream>
#include <cstdio>
#include <mpi.h>

using namespace FSRM;

class SimulationConfigTest : public ::testing::Test {
protected:
    void SetUp() override { MPI_Comm_rank(PETSC_COMM_WORLD, &rank); }
    int rank = 0;
};

TEST_F(SimulationConfigTest, DefaultValues) {
    SimulationConfig cfg;
    EXPECT_FALSE(cfg.enable_geomechanics);
    EXPECT_EQ(cfg.fluid_model, FluidModelType::SINGLE_COMPONENT);
    EXPECT_EQ(cfg.solid_model, SolidModelType::ELASTIC);
    EXPECT_TRUE(cfg.enabled_modules.empty());
    EXPECT_DOUBLE_EQ(cfg.start_time, 0.0);
    EXPECT_DOUBLE_EQ(cfg.end_time, 1.0);
}

TEST_F(SimulationConfigTest, LoadModulesSection) {
    const char* fname = "test_sim_cfg_modules.cfg";
    if (rank == 0) {
        std::ofstream c(fname);
        c << "[SIMULATION]\nend_time = 3.0\n";
        c << "[MODULES]\nenabled = elastodynamics\n";
        c.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(fname));
    SimulationConfig cfg;
    ASSERT_TRUE(reader.parseSimulationConfig(cfg));
    ASSERT_FALSE(cfg.enabled_modules.empty());
    EXPECT_EQ(cfg.enabled_modules[0], "elastodynamics");

    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0) {
        std::remove(fname);
    }
}

TEST_F(SimulationConfigTest, BackwardCompatEnableGeomechanics) {
    const char* fname = "test_sim_cfg_geomech.cfg";
    if (rank == 0) {
        std::ofstream c(fname);
        c << "[SIMULATION]\n";
        c << "fluid_model = NONE\n";
        c << "enable_geomechanics = true\n";
        c.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(fname));
    SimulationConfig cfg;
    ASSERT_TRUE(reader.parseSimulationConfig(cfg));
    bool found = false;
    for (const auto& m : cfg.enabled_modules) {
        if (m == "geomechanics") {
            found = true;
        }
    }
    EXPECT_TRUE(found);

    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0) {
        std::remove(fname);
    }
}
