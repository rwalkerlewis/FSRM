/**
 * @file test_benchmarks.cpp
 * @brief Performance benchmark tests
 */

#include <gtest/gtest.h>
#include "Simulator.hpp"
#include "PhysicsKernel.hpp"
#include "FSRM.hpp"
#include <chrono>
#include <vector>
#include <cmath>

using namespace FSRM;

class BenchmarkTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

// ============================================================================
// Kernel Performance Tests
// ============================================================================

TEST_F(BenchmarkTest, SinglePhaseKernelPerformance) {
    SinglePhaseFlowKernel kernel;
    kernel.setProperties(0.2, 100e-15, 1e-9, 1e-3, 1000.0);
    
    PetscScalar u[1] = {10e6};
    PetscScalar u_t[1] = {1000.0};
    PetscScalar u_x[3] = {1e5, 0.0, 0.0};
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1];
    
    const int num_iterations = 100000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u_x[0] = 1e5 + i * 0.001;
        kernel.residual(u, u_t, u_x, a, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "\nSingle phase kernel: " << time_per_eval << " μs/eval\n";
    }
    
    EXPECT_LT(time_per_eval, 100.0) << "Kernel should be fast";
}

TEST_F(BenchmarkTest, GeomechanicsKernelPerformance) {
    GeomechanicsKernel kernel(SolidModelType::ELASTIC);
    kernel.setMaterialProperties(10e9, 0.25, 2500.0);
    
    PetscScalar u[3] = {0.001, 0.0, 0.0};
    PetscScalar u_t[3] = {0.0, 0.0, 0.0};
    PetscScalar u_x[9] = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    PetscScalar a[3] = {0.0, 0.0, 0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[3];
    
    const int num_iterations = 100000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u_x[0] = 0.001 + i * 1e-9;
        kernel.residual(u, u_t, u_x, a, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Geomechanics kernel: " << time_per_eval << " μs/eval\n";
    }
    
    EXPECT_LT(time_per_eval, 100.0);
}

// ============================================================================
// Mathematical Operations Performance
// ============================================================================

TEST_F(BenchmarkTest, MatrixVectorMultiply) {
    // 3x3 matrix-vector multiply performance
    const int N = 3;
    double A[N][N] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    double x_vec[N] = {1, 2, 3};
    double y[N];
    
    const int num_iterations = 1000000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < num_iterations; ++iter) {
        for (int i = 0; i < N; ++i) {
            y[i] = 0.0;
            for (int j = 0; j < N; ++j) {
                y[i] += A[i][j] * x_vec[j];
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_mult = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "3x3 matrix-vector: " << time_per_mult << " ns/mult\n";
    }
    
    EXPECT_LT(time_per_mult, 1000.0);  // Should be < 1 μs
}

TEST_F(BenchmarkTest, LameParameterCalculation) {
    const int num_iterations = 1000000;
    double E = 10e9;
    double nu = 0.25;
    double lambda, mu;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        double E_var = E + i * 1.0;
        lambda = E_var * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu = E_var / (2.0 * (1.0 + nu));
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Lame parameters: " << time_per_calc << " ns/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 500.0);
    EXPECT_GT(lambda, 0.0);  // Use values to prevent optimization
    EXPECT_GT(mu, 0.0);
}

TEST_F(BenchmarkTest, WaveSpeedCalculation) {
    const int num_iterations = 1000000;
    double rho = 2500.0;
    double E = 50e9;
    double nu = 0.25;
    double V_p, V_s;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        double E_var = E + i * 1.0;
        double lambda = E_var * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double mu = E_var / (2.0 * (1.0 + nu));
        V_p = std::sqrt((lambda + 2.0*mu) / rho);
        V_s = std::sqrt(mu / rho);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Wave speeds: " << time_per_calc << " ns/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 1000.0);
    EXPECT_GT(V_p, 0.0);
    EXPECT_GT(V_s, 0.0);
}

// ============================================================================
// Memory Access Performance
// ============================================================================

TEST_F(BenchmarkTest, ArrayAccessPattern) {
    const int N = 1000;
    std::vector<double> data(N * N, 1.0);
    double sum = 0.0;
    
    // Row-major access (good cache behavior)
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            sum += data[i * N + j];
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration_row = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Column-major access (poor cache behavior)
    sum = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            sum += data[i * N + j];
        }
    }
    end = std::chrono::high_resolution_clock::now();
    
    auto duration_col = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "Row-major access: " << duration_row.count() << " μs\n";
        std::cout << "Column-major access: " << duration_col.count() << " μs\n";
    }
    
    // Row-major should typically be faster or equal
    EXPECT_LE(duration_row.count(), duration_col.count() * 2);
    EXPECT_GT(sum, 0.0);  // Use sum to prevent optimization
}

// ============================================================================
// DOF Counting Tests
// ============================================================================

TEST_F(BenchmarkTest, DOFCalculation) {
    int nx = 100, ny = 100, nz = 50;
    
    // Single phase flow: 1 DOF per cell
    int dofs_single_phase = nx * ny * nz;
    EXPECT_EQ(dofs_single_phase, 500000);
    
    // Poroelasticity: 4 DOFs per cell (p, ux, uy, uz)
    int dofs_poroelasticity = nx * ny * nz * 4;
    EXPECT_EQ(dofs_poroelasticity, 2000000);
    
    // Black oil: 3 DOFs per cell (po, sw, sg)
    int dofs_black_oil = nx * ny * nz * 3;
    EXPECT_EQ(dofs_black_oil, 1500000);
    
    if (rank == 0) {
        std::cout << "DOFs for 100x100x50 grid:\n";
        std::cout << "  Single phase: " << dofs_single_phase << "\n";
        std::cout << "  Poroelasticity: " << dofs_poroelasticity << "\n";
        std::cout << "  Black oil: " << dofs_black_oil << "\n";
    }
}
