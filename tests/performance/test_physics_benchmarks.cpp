/**
 * @file test_physics_benchmarks.cpp
 * @brief Physics-specific performance benchmarks
 * 
 * Benchmarks for different physics models:
 * - Poroelasticity
 * - Fracture mechanics
 * - Wave propagation
 * - Two-phase flow
 * - Thermal diffusion
 */

#include <gtest/gtest.h>
#include "PhysicsKernel.hpp"
#include "PoroelasticSolver.hpp"
#include "FractureModel.hpp"
#include "TwoPhaseFlow.hpp"
#include <chrono>
#include <vector>
#include <cmath>

using namespace FSRM;

class PhysicsBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

// ============================================================================
// Poroelasticity Benchmarks
// ============================================================================

TEST_F(PhysicsBenchmark, PoroelasticKernelPerformance) {
    PoroelasticKernel kernel;
    
    // Set material properties
    double phi = 0.2;              // Porosity
    double k = 1e-13;              // Permeability (m²)
    double K_s = 36e9;             // Solid bulk modulus (Pa)
    double K_f = 2.2e9;            // Fluid bulk modulus (Pa)
    double mu_s = 14e9;            // Shear modulus (Pa)
    double mu_f = 0.001;           // Fluid viscosity (Pa·s)
    double rho_s = 2650.0;         // Solid density (kg/m³)
    double rho_f = 1000.0;         // Fluid density (kg/m³)
    
    kernel.setProperties(phi, k, K_s, K_f, mu_s, mu_f, rho_s, rho_f);
    
    // State variables: [p, ux, uy, uz]
    PetscScalar u[4] = {1e7, 0.001, 0.0, 0.0};
    PetscScalar u_t[4] = {1000.0, 0.0, 0.0, 0.0};
    PetscScalar u_x[12] = {1e5, 0.0, 0.0,
                           0.0, 0.001, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0};
    PetscScalar a[4] = {0.0, 0.0, 0.0, 0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[4];
    
    const int num_iterations = 50000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u[0] = 1e7 + i * 100.0;
        kernel.residual(u, u_t, u_x, a, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "\n=== Poroelasticity Benchmark ===\n";
        std::cout << "Kernel evaluation: " << time_per_eval << " μs/eval\n";
        std::cout << "Throughput: " << 1.0e6 / time_per_eval << " eval/s\n";
    }
    
    EXPECT_LT(time_per_eval, 200.0);
}

TEST_F(PhysicsBenchmark, BiotCoefficientCalculation) {
    const int num_iterations = 1000000;
    double K_s = 36e9;
    double K_f = 2.2e9;
    double phi = 0.2;
    double alpha;  // Biot coefficient
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        double phi_var = phi + i * 1e-9;
        double K_d = K_s * (1.0 - phi_var);  // Drained bulk modulus (simplified)
        alpha = 1.0 - K_d / K_s;
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Biot coefficient: " << time_per_calc << " ns/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 100.0);
    EXPECT_GT(alpha, 0.0);
}

// ============================================================================
// Fracture Mechanics Benchmarks
// ============================================================================

TEST_F(PhysicsBenchmark, FractureGrowthCalculation) {
    FractureModel fracture;
    fracture.setProperties(1.0e6, 100.0, 0.001);  // K_IC, E, nu
    
    const int num_iterations = 10000;
    double K_I = 0.8e6;  // Stress intensity factor
    double crack_tip[3] = {0.0, 0.0, 0.0};
    double normal[3] = {1.0, 0.0, 0.0};
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        K_I += i * 100.0;
        double delta_a = fracture.computePropagationIncrement(K_I, crack_tip, normal);
        crack_tip[0] += delta_a;
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "\n=== Fracture Mechanics Benchmark ===\n";
        std::cout << "Propagation calc: " << time_per_calc << " μs/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 50.0);
}

TEST_F(PhysicsBenchmark, StressIntensityFactor) {
    const int num_iterations = 100000;
    double sigma = 10e6;  // Applied stress
    double a = 0.1;       // Crack length
    double K_I;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        double a_var = a + i * 1e-7;
        K_I = sigma * std::sqrt(M_PI * a_var);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Stress intensity: " << time_per_calc << " ns/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 200.0);
    EXPECT_GT(K_I, 0.0);
}

// ============================================================================
// Wave Propagation Benchmarks
// ============================================================================

TEST_F(PhysicsBenchmark, ElasticWaveKernelPerformance) {
    ElastodynamicsKernel kernel;
    kernel.setMaterialProperties(50e9, 0.25, 2700.0);  // E, nu, rho
    
    // State variables: [ux, uy, uz]
    PetscScalar u[3] = {0.001, 0.0, 0.0};
    PetscScalar u_t[3] = {1.0, 0.0, 0.0};
    PetscScalar u_tt[3] = {0.0, 0.0, 0.0};
    PetscScalar u_x[9] = {0.001, 0.0, 0.0,
                          0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[3];
    
    const int num_iterations = 50000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u_t[0] = 1.0 + i * 0.0001;
        kernel.residualElastodynamics(u, u_t, u_tt, u_x, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "\n=== Wave Propagation Benchmark ===\n";
        std::cout << "Elastodynamics: " << time_per_eval << " μs/eval\n";
    }
    
    EXPECT_LT(time_per_eval, 150.0);
}

TEST_F(PhysicsBenchmark, PoroelasticWaveKernelPerformance) {
    PoroelastodynamicsKernel kernel;
    kernel.setProperties(0.2, 1e-13, 36e9, 2.2e9, 14e9, 0.001, 2650.0, 1000.0);
    
    // State variables: [p, ux, uy, uz]
    PetscScalar u[4] = {1e7, 0.001, 0.0, 0.0};
    PetscScalar u_t[4] = {1000.0, 1.0, 0.0, 0.0};
    PetscScalar u_tt[4] = {0.0, 0.0, 0.0, 0.0};
    PetscScalar u_x[12];
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[4];
    
    const int num_iterations = 30000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u_t[0] = 1000.0 + i * 1.0;
        kernel.residualPoroelastodynamics(u, u_t, u_tt, u_x, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Poroelastodynamics: " << time_per_eval << " μs/eval\n";
    }
    
    EXPECT_LT(time_per_eval, 250.0);
}

// ============================================================================
// Two-Phase Flow Benchmarks
// ============================================================================

TEST_F(PhysicsBenchmark, TwoPhaseFlowKernelPerformance) {
    TwoPhaseFlowKernel kernel;
    kernel.setProperties(0.2, 100e-15, 1e-9, 0.003, 0.0005, 850.0, 1000.0);
    
    // State variables: [p, Sw]
    PetscScalar u[2] = {1e7, 0.5};
    PetscScalar u_t[2] = {1000.0, 0.01};
    PetscScalar u_x[6] = {1e5, 0.0, 0.0, 0.0, 0.0, 0.0};
    PetscScalar a[2] = {0.0, 0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[2];
    
    const int num_iterations = 50000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u[1] = 0.3 + (i % 5000) * 0.0001;  // Vary saturation
        kernel.residual(u, u_t, u_x, a, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "\n=== Two-Phase Flow Benchmark ===\n";
        std::cout << "Kernel evaluation: " << time_per_eval << " μs/eval\n";
    }
    
    EXPECT_LT(time_per_eval, 100.0);
}

TEST_F(PhysicsBenchmark, RelativePermeabilityCalculation) {
    const int num_iterations = 1000000;
    double Sw, kr_w, kr_o;
    
    // Corey model parameters
    double Sw_c = 0.2;   // Connate water
    double Sor = 0.2;    // Residual oil
    double n_w = 2.0;    // Water exponent
    double n_o = 2.0;    // Oil exponent
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        Sw = 0.2 + (i % 10000) * 0.00006;  // Range [0.2, 0.8]
        double S_star = (Sw - Sw_c) / (1.0 - Sw_c - Sor);
        kr_w = std::pow(S_star, n_w);
        kr_o = std::pow(1.0 - S_star, n_o);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Rel. perm (Corey): " << time_per_calc << " ns/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 100.0);
    EXPECT_GE(kr_w, 0.0);
    EXPECT_LE(kr_o, 1.0);
}

TEST_F(PhysicsBenchmark, CapillaryPressureCalculation) {
    const int num_iterations = 1000000;
    double Sw, Pc;
    
    // Brooks-Corey model
    double Pd = 1000.0;   // Entry pressure (Pa)
    double lambda = 2.0;  // Pore size distribution
    double Sw_c = 0.2;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        Sw = 0.2 + (i % 10000) * 0.00006;
        double S_e = (Sw - Sw_c) / (1.0 - Sw_c);
        Pc = Pd * std::pow(S_e, -1.0/lambda);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "Capillary pressure: " << time_per_calc << " ns/calc\n";
    }
    
    EXPECT_LT(time_per_calc, 150.0);
    EXPECT_GT(Pc, 0.0);
}

// ============================================================================
// Thermal Benchmarks
// ============================================================================

TEST_F(PhysicsBenchmark, ThermalDiffusionKernelPerformance) {
    ThermalKernel kernel;
    kernel.setProperties(2.5, 800.0, 2650.0);  // k_thermal, cp, rho
    
    PetscScalar u[1] = {350.0};  // Temperature in K
    PetscScalar u_t[1] = {0.1};
    PetscScalar u_x[3] = {10.0, 0.0, 0.0};  // Temp gradient
    PetscScalar a[1] = {0.0};
    PetscReal x[3] = {0.0, 0.0, 0.0};
    PetscScalar f[1];
    
    const int num_iterations = 100000;
    
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        u[0] = 300.0 + i * 0.001;
        kernel.residual(u, u_t, u_x, a, x, f);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_iterations;
    
    if (rank == 0) {
        std::cout << "\n=== Thermal Benchmark ===\n";
        std::cout << "Heat diffusion: " << time_per_eval << " μs/eval\n";
    }
    
    EXPECT_LT(time_per_eval, 50.0);
}

// ============================================================================
// Grid Size Scalability Tests
// ============================================================================

TEST_F(PhysicsBenchmark, GridSizeScalability) {
    std::vector<int> grid_sizes = {10, 20, 40, 80, 100};
    
    if (rank == 0) {
        std::cout << "\n=== Grid Size Scalability ===\n";
        std::cout << "Grid Size | DOFs | Time/iter (ms) | DOFs/sec\n";
        std::cout << "----------|------|----------------|----------\n";
    }
    
    for (int n : grid_sizes) {
        int num_cells = n * n * n;
        int dofs_per_cell = 4;  // Poroelasticity
        int total_dofs = num_cells * dofs_per_cell;
        
        // Simulate work proportional to DOFs
        double result = 0.0;
        auto start = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < total_dofs; ++i) {
            result += std::sin(i * 0.001) * std::cos(i * 0.001);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        double throughput = total_dofs / (duration.count() / 1.0e6);
        
        if (rank == 0) {
            printf("%3dx%3dx%3d | %7d | %14.3f | %10.0f\n", 
                   n, n, n, total_dofs, time_ms, throughput);
        }
        
        EXPECT_GT(result, 0.0);  // Use result
    }
}
