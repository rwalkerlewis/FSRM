/**
 * @file test_restart.cpp
 * @brief Checkpoint / persistence smoke tests
 */

#include "core/Simulator.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include <filesystem>

using namespace FSRM;

class RestartTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    void TearDown() override {
        MPI_Barrier(PETSC_COMM_WORLD);
        if (rank == 0) {
            std::filesystem::remove(summary_path_);
            std::remove("test_restart.config");
            std::remove("checkpoint.bin");
        }
    }

    int rank = 0;
    static constexpr const char* summary_path_ = "restart_smoke_out_SUMMARY.txt";
};

TEST_F(RestartTest, CheckpointConfigEnabled) {
    SimulationConfig config{};
    config.enable_checkpointing = true;
    config.checkpoint_frequency = 50;
    EXPECT_TRUE(config.enable_checkpointing);
    EXPECT_EQ(config.checkpoint_frequency, 50);
}

TEST_F(RestartTest, WriteCheckpointReturnsSuccess) {
    Simulator sim(PETSC_COMM_WORLD);
    SimulationConfig cfg{};
    cfg.output_file = "restart_smoke_out";
    PetscErrorCode ierr = sim.initialize(cfg);
    ASSERT_EQ(ierr, 0);
    ierr = sim.writeCheckpoint(0);
    EXPECT_EQ(ierr, 0);
}

TEST_F(RestartTest, SummaryFileExistsAfterWriteSummary) {
    Simulator sim(PETSC_COMM_WORLD);
    SimulationConfig cfg{};
    cfg.output_file = "restart_smoke_out";
    PetscErrorCode ierr = sim.initialize(cfg);
    ASSERT_EQ(ierr, 0);
    ierr = sim.writeCheckpoint(0);
    ASSERT_EQ(ierr, 0);
    ierr = sim.writeSummary();
    ASSERT_EQ(ierr, 0);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0) {
        EXPECT_TRUE(std::filesystem::exists(summary_path_));
    }
}

TEST_F(RestartTest, OutputFormatOptions) {
    SimulationConfig config{};
    config.output_format = "ECLIPSE";
    EXPECT_EQ(config.output_format, "ECLIPSE");
    config.output_format = "VTK";
    EXPECT_EQ(config.output_format, "VTK");
    config.output_format = "HDF5";
    EXPECT_EQ(config.output_format, "HDF5");
}
