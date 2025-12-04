/**
 * @file test_config_reader.cpp
 * @brief Unit tests for ConfigReader class
 */

#include <gtest/gtest.h>
#include "ConfigReader.hpp"
#include "FSRM.hpp"
#include <fstream>
#include <cstdio>

using namespace FSRM;

class ConfigReaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        test_config_file = "test_config_unit.config";
        
        if (rank == 0) {
            std::ofstream config(test_config_file);
            config << "[simulation]\n";
            config << "physics_type = FLUID_FLOW\n";
            config << "timesteps = 100\n";
            config << "dt = 86400.0\n";
            config << "\n[grid]\n";
            config << "nx = 10\n";
            config << "ny = 20\n";
            config << "nz = 5\n";
            config << "Lx = 100.0\n";
            config << "Ly = 200.0\n";
            config << "Lz = 50.0\n";
            config << "\n[material]\n";
            config << "porosity = 0.2\n";
            config << "permeability_x = 100.0\n";
            config << "output_frequency = 10\n";
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
    int rank;
};

TEST_F(ConfigReaderTest, LoadConfigFile) {
    ConfigReader reader;
    bool loaded = reader.loadFile(test_config_file);
    EXPECT_TRUE(loaded) << "Should load config file successfully";
}

TEST_F(ConfigReaderTest, ReadIntegerValues) {
    ConfigReader reader;
    reader.loadFile(test_config_file);
    
    EXPECT_EQ(reader.getInt("grid", "nx", 0), 10);
    EXPECT_EQ(reader.getInt("grid", "ny", 0), 20);
    EXPECT_EQ(reader.getInt("grid", "nz", 0), 5);
    EXPECT_EQ(reader.getInt("simulation", "timesteps", 0), 100);
}

TEST_F(ConfigReaderTest, ReadDoubleValues) {
    ConfigReader reader;
    reader.loadFile(test_config_file);
    
    EXPECT_DOUBLE_EQ(reader.getDouble("grid", "Lx", 0.0), 100.0);
    EXPECT_DOUBLE_EQ(reader.getDouble("grid", "Ly", 0.0), 200.0);
    EXPECT_DOUBLE_EQ(reader.getDouble("grid", "Lz", 0.0), 50.0);
    EXPECT_DOUBLE_EQ(reader.getDouble("material", "porosity", 0.0), 0.2);
}

TEST_F(ConfigReaderTest, DefaultValues) {
    ConfigReader reader;
    reader.loadFile(test_config_file);
    
    // Non-existent key with default
    EXPECT_EQ(reader.getInt("nonexistent", "key", 42), 42);
    EXPECT_DOUBLE_EQ(reader.getDouble("nonexistent", "key", 3.14), 3.14);
}

TEST_F(ConfigReaderTest, HasSection) {
    ConfigReader reader;
    reader.loadFile(test_config_file);
    
    EXPECT_TRUE(reader.hasSection("simulation"));
    EXPECT_TRUE(reader.hasSection("grid"));
    EXPECT_TRUE(reader.hasSection("material"));
    EXPECT_FALSE(reader.hasSection("nonexistent"));
}

TEST_F(ConfigReaderTest, HasKey) {
    ConfigReader reader;
    reader.loadFile(test_config_file);
    
    EXPECT_TRUE(reader.hasKey("grid", "nx"));
    EXPECT_TRUE(reader.hasKey("material", "porosity"));
    EXPECT_FALSE(reader.hasKey("grid", "nonexistent"));
}

TEST_F(ConfigReaderTest, ParseGridConfig) {
    ConfigReader reader;
    reader.loadFile(test_config_file);
    
    GridConfig grid;
    bool parsed = reader.parseGridConfig(grid);
    
    if (parsed) {
        EXPECT_EQ(grid.nx, 10);
        EXPECT_EQ(grid.ny, 20);
        EXPECT_EQ(grid.nz, 5);
        EXPECT_DOUBLE_EQ(grid.Lx, 100.0);
    }
}
