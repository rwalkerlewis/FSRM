/**
 * @file test_scaling.cpp
 * @brief Parallel scaling tests
 */

#include <gtest/gtest.h>
#include "ReservoirSim.hpp"
#include <chrono>
#include <vector>
#include <cmath>

using namespace FSRM;

class ScalingTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

TEST_F(ScalingTest, MPIRankAndSize) {
    EXPECT_GE(rank, 0);
    EXPECT_GE(size, 1);
    EXPECT_LT(rank, size);
    
    if (rank == 0) {
        std::cout << "\nRunning with " << size << " MPI processes\n";
    }
}

TEST_F(ScalingTest, MPIBroadcast) {
    double value = 0.0;
    
    if (rank == 0) {
        value = 42.0;
    }
    
    MPI_Bcast(&value, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    
    EXPECT_DOUBLE_EQ(value, 42.0);
}

TEST_F(ScalingTest, MPIReduce) {
    double local_value = rank + 1.0;  // Each rank contributes its rank+1
    double global_sum = 0.0;
    
    MPI_Reduce(&local_value, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
    
    if (rank == 0) {
        double expected_sum = size * (size + 1) / 2.0;
        EXPECT_DOUBLE_EQ(global_sum, expected_sum);
    }
}

TEST_F(ScalingTest, MPIAllReduce) {
    double local_value = 1.0;
    double global_sum = 0.0;
    
    MPI_Allreduce(&local_value, &global_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    
    EXPECT_DOUBLE_EQ(global_sum, (double)size);
}

TEST_F(ScalingTest, DomainDecomposition) {
    // Simulate domain decomposition for a 3D grid
    int total_nx = 100;
    int total_ny = 100;
    int total_nz = 50;
    
    // Simple 1D decomposition in x-direction
    int local_nx = total_nx / size;
    int remainder = total_nx % size;
    
    if (rank < remainder) {
        local_nx++;
    }
    
    int local_ny = total_ny;
    int local_nz = total_nz;
    
    int local_cells = local_nx * local_ny * local_nz;
    int total_cells = total_nx * total_ny * total_nz;
    
    // Verify all cells are accounted for
    int sum_cells = 0;
    MPI_Reduce(&local_cells, &sum_cells, 1, MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD);
    
    if (rank == 0) {
        EXPECT_EQ(sum_cells, total_cells);
        std::cout << "Total cells: " << total_cells << "\n";
        std::cout << "Cells per process (avg): " << total_cells / size << "\n";
    }
}

TEST_F(ScalingTest, LoadBalance) {
    // Each process does work proportional to local data
    int work_per_cell = 100;  // Operations per cell
    int local_cells = (1000 * 1000) / size;  // 1M cells total
    
    double local_result = 0.0;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < local_cells; ++i) {
        for (int j = 0; j < work_per_cell; ++j) {
            local_result += std::sin(i * 0.001) * std::cos(j * 0.001);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double local_time = duration.count();
    
    double min_time, max_time, avg_time;
    MPI_Reduce(&local_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, PETSC_COMM_WORLD);
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD);
    MPI_Reduce(&local_time, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
    
    if (rank == 0) {
        avg_time /= size;
        double imbalance = (max_time - min_time) / avg_time;
        
        std::cout << "Load balance test:\n";
        std::cout << "  Min time: " << min_time << " ms\n";
        std::cout << "  Max time: " << max_time << " ms\n";
        std::cout << "  Avg time: " << avg_time << " ms\n";
        std::cout << "  Imbalance: " << imbalance * 100 << "%\n";
        
        // Imbalance should be small for uniform work distribution
        EXPECT_LT(imbalance, 0.3);
    }
    
    EXPECT_TRUE(std::isfinite(local_result));
}

TEST_F(ScalingTest, CommunicationOverhead) {
    // Measure point-to-point communication
    const int msg_size = 10000;
    std::vector<double> send_buf(msg_size, 1.0);
    std::vector<double> recv_buf(msg_size, 0.0);
    
    MPI_Barrier(PETSC_COMM_WORLD);
    auto start = std::chrono::high_resolution_clock::now();
    
    // Ring communication pattern
    int dest = (rank + 1) % size;
    int source = (rank - 1 + size) % size;
    
    MPI_Request req[2];
    MPI_Isend(send_buf.data(), msg_size, MPI_DOUBLE, dest, 0, PETSC_COMM_WORLD, &req[0]);
    MPI_Irecv(recv_buf.data(), msg_size, MPI_DOUBLE, source, 0, PETSC_COMM_WORLD, &req[1]);
    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    
    MPI_Barrier(PETSC_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "Ring communication (" << msg_size << " doubles): " 
                  << duration.count() << " μs\n";
    }
    
    // Verify data received
    EXPECT_DOUBLE_EQ(recv_buf[0], 1.0);
}

TEST_F(ScalingTest, GlobalReductionPerformance) {
    std::vector<double> local_data(10000);
    for (size_t i = 0; i < local_data.size(); ++i) {
        local_data[i] = rank * 1000 + i;
    }
    
    double local_sum = 0.0;
    for (double v : local_data) {
        local_sum += v;
    }
    
    double global_sum = 0.0;
    
    MPI_Barrier(PETSC_COMM_WORLD);
    auto start = std::chrono::high_resolution_clock::now();
    
    const int num_reductions = 100;
    for (int i = 0; i < num_reductions; ++i) {
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    }
    
    MPI_Barrier(PETSC_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_reduce = duration.count() / (double)num_reductions;
    
    if (rank == 0) {
        std::cout << "Allreduce (1 double, " << size << " procs): " 
                  << time_per_reduce << " μs\n";
    }
    
    EXPECT_GT(global_sum, 0.0);
}
