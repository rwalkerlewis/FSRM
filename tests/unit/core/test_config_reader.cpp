/**
 * @file test_config_reader.cpp
 * @brief Unit tests for ConfigReader
 */

#include <gtest/gtest.h>
#include "core/ConfigReader.hpp"
#include "core/FSRM.hpp"
#include <fstream>
#include <cstdio>
#include <mpi.h>

using namespace FSRM;

class ConfigReaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        test_config_file = "test_config_reader_unit.cfg";
        if (rank == 0) {
            std::ofstream config(test_config_file);
            config << "[demo]\n";
            config << "message = hello_fsrm\n";
            config << "count = 7\n";
            config << "pi = 3.25\n";
            config << "\n[SIMULATION]\n";
            config << "start_time = 0.0\n";
            config << "end_time = 2.0\n";
            config << "fluid_model = SINGLE_COMPONENT\n";
            config << "\n[GRID]\n";
            config << "nx = 10\n";
            config << "ny = 20\n";
            config << "nz = 5\n";
            config << "Lx = 100.0\n";
            config << "Ly = 200.0\n";
            config << "Lz = 50.0\n";
            config << "mesh_type = CARTESIAN\n";
            config.close();
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }

    void TearDown() override {
        MPI_Barrier(PETSC_COMM_WORLD);
        if (rank == 0) {
            std::remove(test_config_file.c_str());
        }
    }

    std::string test_config_file;
    int rank = 0;
};

TEST_F(ConfigReaderTest, LoadConfigFile) {
    ConfigReader reader;
    EXPECT_TRUE(reader.loadFile(test_config_file));
}

TEST_F(ConfigReaderTest, ReadIntegerDoubleString) {
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(test_config_file));
    EXPECT_EQ(reader.getInt("demo", "count", 0), 7);
    EXPECT_DOUBLE_EQ(reader.getDouble("demo", "pi", 0.0), 3.25);
    EXPECT_EQ(reader.getString("demo", "message", ""), "hello_fsrm");
}

TEST_F(ConfigReaderTest, DefaultValuesForMissingKeys) {
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(test_config_file));
    EXPECT_EQ(reader.getInt("demo", "missing_int", 42), 42);
    EXPECT_DOUBLE_EQ(reader.getDouble("demo", "missing_double", 2.718), 2.718);
    EXPECT_EQ(reader.getString("demo", "missing_str", "default"), "default");
}

TEST_F(ConfigReaderTest, HasSectionAndKey) {
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(test_config_file));
    EXPECT_TRUE(reader.hasSection("demo"));
    EXPECT_TRUE(reader.hasSection("SIMULATION"));
    EXPECT_TRUE(reader.hasSection("GRID"));
    EXPECT_TRUE(reader.hasKey("demo", "message"));
    EXPECT_FALSE(reader.hasSection("no_such_section"));
    EXPECT_FALSE(reader.hasKey("demo", "no_such_key"));
}

TEST_F(ConfigReaderTest, ParseGridConfig) {
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(test_config_file));
    GridConfig grid;
    ASSERT_TRUE(reader.parseGridConfig(grid));
    EXPECT_EQ(grid.nx, 10);
    EXPECT_EQ(grid.ny, 20);
    EXPECT_EQ(grid.nz, 5);
    EXPECT_DOUBLE_EQ(grid.Lx, 100.0);
    EXPECT_DOUBLE_EQ(grid.Ly, 200.0);
    EXPECT_DOUBLE_EQ(grid.Lz, 50.0);
    EXPECT_EQ(grid.mesh_type, MeshType::CARTESIAN);
}

TEST_F(ConfigReaderTest, ParseSimulationConfig) {
    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(test_config_file));
    SimulationConfig sim;
    ASSERT_TRUE(reader.parseSimulationConfig(sim));
    EXPECT_DOUBLE_EQ(sim.start_time, 0.0);
    EXPECT_DOUBLE_EQ(sim.end_time, 2.0);
    EXPECT_EQ(sim.fluid_model, FluidModelType::SINGLE_COMPONENT);
}

TEST_F(ConfigReaderTest, ModulesSectionPopulatesEnabledModules) {
    const char* fname = "test_modules_section.cfg";
    if (rank == 0) {
        std::ofstream c(fname);
        c << "[SIMULATION]\nend_time = 1.0\n";
        c << "[MODULES]\nenabled = thermal, fluid_flow\n";
        c.close();
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    ConfigReader reader;
    ASSERT_TRUE(reader.loadFile(fname));
    SimulationConfig sim;
    ASSERT_TRUE(reader.parseSimulationConfig(sim));
    ASSERT_EQ(sim.enabled_modules.size(), 2u);
    EXPECT_EQ(sim.enabled_modules[0], "thermal");
    EXPECT_EQ(sim.enabled_modules[1], "fluid_flow");

    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0) {
        std::remove(fname);
    }
}

TEST_F(ConfigReaderTest, BackwardCompatGeomechanicsFlag) {
    const char* fname = "test_geomech_compat.cfg";
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
    SimulationConfig sim;
    ASSERT_TRUE(reader.parseSimulationConfig(sim));
    bool found = false;
    for (const auto& m : sim.enabled_modules) {
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
