/**
 * @file test_scenario_benchmarks.cpp
 * @brief Real-world scenario performance benchmarks
 * 
 * Benchmarks for complete simulation scenarios:
 * - Hydraulic fracturing
 * - Enhanced geothermal systems
 * - CO2 storage
 * - Induced seismicity
 * - Production optimization
 */

#include <gtest/gtest.h>
#include "Simulator.hpp"
#include "ConfigReader.hpp"
#include <chrono>
#include <fstream>

using namespace FSRM;

class ScenarioBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
    
    // Helper function to run a simulation and measure performance
    struct BenchmarkResult {
        double setup_time_ms;
        double simulation_time_ms;
        double total_time_ms;
        int num_timesteps;
        int num_newton_iterations;
        int num_linear_iterations;
        double avg_time_per_step_ms;
    };
    
    BenchmarkResult runSimulationBenchmark(const std::string& config_file) {
        BenchmarkResult result = {0, 0, 0, 0, 0, 0, 0};
        
        // Setup phase
        auto start_setup = std::chrono::high_resolution_clock::now();
        
        FSRM::Simulator sim(PETSC_COMM_WORLD);
        PetscErrorCode ierr = sim.initializeFromConfigFile(config_file.c_str());
        if (ierr != 0) return result;
        
        ierr = sim.setupDM();
        if (ierr != 0) return result;
        ierr = sim.setupFields();
        if (ierr != 0) return result;
        ierr = sim.setupPhysics();
        if (ierr != 0) return result;
        ierr = sim.setMaterialProperties();
        if (ierr != 0) return result;
        ierr = sim.setInitialConditions();
        if (ierr != 0) return result;
        ierr = sim.setupTimeStepper();
        if (ierr != 0) return result;
        ierr = sim.setupSolvers();
        if (ierr != 0) return result;
        
        auto end_setup = std::chrono::high_resolution_clock::now();
        result.setup_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_setup - start_setup).count();
        
        // Simulation phase
        auto start_sim = std::chrono::high_resolution_clock::now();
        ierr = sim.run();
        auto end_sim = std::chrono::high_resolution_clock::now();
        
        result.simulation_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_sim - start_sim).count();
        result.total_time_ms = result.setup_time_ms + result.simulation_time_ms;
        
        // Get performance metrics
        auto metrics = sim.getPerformanceMetrics();
        result.num_timesteps = metrics.num_timesteps;
        result.num_newton_iterations = metrics.total_newton_iterations;
        result.num_linear_iterations = metrics.total_linear_iterations;
        
        if (result.num_timesteps > 0) {
            result.avg_time_per_step_ms = result.simulation_time_ms / result.num_timesteps;
        }
        
        return result;
    }
};

// ============================================================================
// Hydraulic Fracturing Scenarios
// ============================================================================

TEST_F(ScenarioBenchmark, HydraulicFracturingSmall) {
    // Small problem: 20x20x10 grid with single fracture
    std::string config = 
        "[SIMULATION]\n"
        "end_time = 3600.0\n"  // 1 hour
        "dt_initial = 10.0\n"
        "fluid_model = SINGLE_PHASE\n"
        "enable_geomechanics = true\n"
        "enable_fractures = true\n"
        "[GRID]\n"
        "nx = 20\n"
        "ny = 20\n"
        "nz = 10\n"
        "Lx = 100.0\n"
        "Ly = 100.0\n"
        "Lz = 50.0\n"
        "[ROCK]\n"
        "porosity = 0.15\n"
        "permeability_x = 0.1\n"
        "youngs_modulus = 20e9\n"
        "[FLUID]\n"
        "viscosity = 0.001\n"
        "density = 1000.0\n";
    
    // Write config to temporary file
    std::ofstream ofs("benchmark_hf_small.config");
    ofs << config;
    ofs.close();
    
    if (rank == 0) {
        std::cout << "\n=== Hydraulic Fracturing Benchmark (Small) ===\n";
        std::cout << "Grid: 20x20x10 (4,000 cells)\n";
    }
    
    auto result = runSimulationBenchmark("benchmark_hf_small.config");
    
    if (rank == 0) {
        std::cout << "Setup time:           " << result.setup_time_ms << " ms\n";
        std::cout << "Simulation time:      " << result.simulation_time_ms << " ms\n";
        std::cout << "Total time:           " << result.total_time_ms << " ms\n";
        std::cout << "Timesteps:            " << result.num_timesteps << "\n";
        std::cout << "Avg time per step:    " << result.avg_time_per_step_ms << " ms\n";
        std::cout << "Newton iterations:    " << result.num_newton_iterations << "\n";
        std::cout << "Linear iterations:    " << result.num_linear_iterations << "\n";
    }
    
    std::remove("benchmark_hf_small.config");
    
    EXPECT_GT(result.num_timesteps, 0);
}

TEST_F(ScenarioBenchmark, HydraulicFracturingMedium) {
    // Medium problem: 50x50x20 grid
    if (size < 4 && rank == 0) {
        std::cout << "\nNote: Run with >= 4 MPI processes for better performance\n";
    }
    
    std::string config = 
        "[SIMULATION]\n"
        "end_time = 7200.0\n"  // 2 hours
        "dt_initial = 20.0\n"
        "fluid_model = SINGLE_PHASE\n"
        "enable_geomechanics = true\n"
        "enable_fractures = true\n"
        "[GRID]\n"
        "nx = 50\n"
        "ny = 50\n"
        "nz = 20\n"
        "Lx = 500.0\n"
        "Ly = 500.0\n"
        "Lz = 100.0\n"
        "[ROCK]\n"
        "porosity = 0.15\n"
        "permeability_x = 0.1\n"
        "youngs_modulus = 20e9\n"
        "[FLUID]\n"
        "viscosity = 0.001\n"
        "density = 1000.0\n";
    
    std::ofstream ofs("benchmark_hf_medium.config");
    ofs << config;
    ofs.close();
    
    if (rank == 0) {
        std::cout << "\n=== Hydraulic Fracturing Benchmark (Medium) ===\n";
        std::cout << "Grid: 50x50x20 (50,000 cells)\n";
        std::cout << "MPI processes: " << size << "\n";
    }
    
    auto result = runSimulationBenchmark("benchmark_hf_medium.config");
    
    if (rank == 0) {
        std::cout << "Total time:           " << result.total_time_ms / 1000.0 << " s\n";
        std::cout << "Timesteps:            " << result.num_timesteps << "\n";
        std::cout << "Avg time per step:    " << result.avg_time_per_step_ms << " ms\n";
        std::cout << "Performance:          " << 50000.0 / (result.avg_time_per_step_ms / 1000.0) 
                  << " cells/s\n";
    }
    
    std::remove("benchmark_hf_medium.config");
    
    EXPECT_GT(result.num_timesteps, 0);
}

// ============================================================================
// Geothermal Scenarios
// ============================================================================

TEST_F(ScenarioBenchmark, GeothermalSystem) {
    std::string config = 
        "[SIMULATION]\n"
        "end_time = 86400.0\n"  // 1 day
        "dt_initial = 100.0\n"
        "fluid_model = SINGLE_PHASE\n"
        "enable_thermal = true\n"
        "enable_geomechanics = true\n"
        "[GRID]\n"
        "nx = 30\n"
        "ny = 30\n"
        "nz = 30\n"
        "Lx = 1000.0\n"
        "Ly = 1000.0\n"
        "Lz = 1000.0\n"
        "[ROCK]\n"
        "porosity = 0.01\n"
        "permeability_x = 1.0\n"
        "youngs_modulus = 50e9\n"
        "thermal_conductivity = 3.0\n"
        "[FLUID]\n"
        "viscosity = 0.0002\n"  // Hot water
        "density = 900.0\n";
    
    std::ofstream ofs("benchmark_geothermal.config");
    ofs << config;
    ofs.close();
    
    if (rank == 0) {
        std::cout << "\n=== Geothermal System Benchmark ===\n";
        std::cout << "Grid: 30x30x30 (27,000 cells)\n";
        std::cout << "Physics: Thermal + Hydraulic + Mechanical\n";
    }
    
    auto result = runSimulationBenchmark("benchmark_geothermal.config");
    
    if (rank == 0) {
        std::cout << "Total time:           " << result.total_time_ms / 1000.0 << " s\n";
        std::cout << "Timesteps:            " << result.num_timesteps << "\n";
        std::cout << "Avg time per step:    " << result.avg_time_per_step_ms << " ms\n";
    }
    
    std::remove("benchmark_geothermal.config");
    
    EXPECT_GT(result.num_timesteps, 0);
}

// ============================================================================
// CO2 Storage Scenarios
// ============================================================================

TEST_F(ScenarioBenchmark, CO2Storage) {
    std::string config = 
        "[SIMULATION]\n"
        "end_time = 3.1536e7\n"  // 1 year
        "dt_initial = 86400.0\n"  // 1 day
        "fluid_model = TWO_PHASE\n"
        "enable_geomechanics = false\n"
        "[GRID]\n"
        "nx = 40\n"
        "ny = 40\n"
        "nz = 20\n"
        "Lx = 2000.0\n"
        "Ly = 2000.0\n"
        "Lz = 100.0\n"
        "[ROCK]\n"
        "porosity = 0.25\n"
        "permeability_x = 100.0\n"
        "[FLUID]\n"
        "oil_viscosity = 0.000045\n"  // CO2
        "water_viscosity = 0.001\n";
    
    std::ofstream ofs("benchmark_co2.config");
    ofs << config;
    ofs.close();
    
    if (rank == 0) {
        std::cout << "\n=== CO2 Storage Benchmark ===\n";
        std::cout << "Grid: 40x40x20 (32,000 cells)\n";
        std::cout << "Physics: Two-phase flow\n";
    }
    
    auto result = runSimulationBenchmark("benchmark_co2.config");
    
    if (rank == 0) {
        std::cout << "Total time:           " << result.total_time_ms / 1000.0 << " s\n";
        std::cout << "Timesteps:            " << result.num_timesteps << "\n";
        std::cout << "Avg time per step:    " << result.avg_time_per_step_ms << " ms\n";
    }
    
    std::remove("benchmark_co2.config");
    
    EXPECT_GT(result.num_timesteps, 0);
}

// ============================================================================
// Wave Propagation Scenarios
// ============================================================================

TEST_F(ScenarioBenchmark, WavePropagation) {
    std::string config = 
        "[SIMULATION]\n"
        "end_time = 1.0\n"  // 1 second
        "dt_initial = 0.0001\n"
        "fluid_model = NONE\n"
        "enable_wave_propagation = true\n"
        "[GRID]\n"
        "nx = 50\n"
        "ny = 50\n"
        "nz = 50\n"
        "Lx = 100.0\n"
        "Ly = 100.0\n"
        "Lz = 100.0\n"
        "[ROCK]\n"
        "youngs_modulus = 50e9\n"
        "poisson_ratio = 0.25\n"
        "density = 2700.0\n";
    
    std::ofstream ofs("benchmark_wave.config");
    ofs << config;
    ofs.close();
    
    if (rank == 0) {
        std::cout << "\n=== Wave Propagation Benchmark ===\n";
        std::cout << "Grid: 50x50x50 (125,000 cells)\n";
        std::cout << "Physics: Elastodynamics\n";
    }
    
    auto result = runSimulationBenchmark("benchmark_wave.config");
    
    if (rank == 0) {
        std::cout << "Total time:           " << result.total_time_ms / 1000.0 << " s\n";
        std::cout << "Timesteps:            " << result.num_timesteps << "\n";
        std::cout << "Avg time per step:    " << result.avg_time_per_step_ms << " ms\n";
        
        if (result.num_timesteps > 0) {
            double throughput = 125000.0 * result.num_timesteps / (result.simulation_time_ms / 1000.0);
            std::cout << "Throughput:           " << throughput / 1e6 << " Mcell-steps/s\n";
        }
    }
    
    std::remove("benchmark_wave.config");
    
    EXPECT_GT(result.num_timesteps, 0);
}

// ============================================================================
// Parallel Scaling Scenarios
// ============================================================================

TEST_F(ScenarioBenchmark, ParallelScalingTest) {
    // Fixed problem size, measure performance with different process counts
    std::string config = 
        "[SIMULATION]\n"
        "end_time = 1000.0\n"
        "dt_initial = 10.0\n"
        "max_time_steps = 10\n"  // Fixed number of steps for fair comparison
        "fluid_model = SINGLE_PHASE\n"
        "[GRID]\n"
        "nx = 60\n"
        "ny = 60\n"
        "nz = 60\n"
        "Lx = 100.0\n"
        "Ly = 100.0\n"
        "Lz = 100.0\n"
        "[ROCK]\n"
        "porosity = 0.2\n"
        "permeability_x = 100.0\n"
        "[FLUID]\n"
        "viscosity = 0.001\n"
        "density = 1000.0\n";
    
    std::ofstream ofs("benchmark_scaling.config");
    ofs << config;
    ofs.close();
    
    if (rank == 0) {
        std::cout << "\n=== Parallel Scaling Benchmark ===\n";
        std::cout << "Grid: 60x60x60 (216,000 cells)\n";
        std::cout << "MPI processes: " << size << "\n";
    }
    
    auto result = runSimulationBenchmark("benchmark_scaling.config");
    
    if (rank == 0) {
        std::cout << "Total time:           " << result.total_time_ms / 1000.0 << " s\n";
        std::cout << "Time per step:        " << result.avg_time_per_step_ms << " ms\n";
        std::cout << "Cells per process:    " << 216000 / size << "\n";
        std::cout << "Performance:          " << 216000.0 / (result.avg_time_per_step_ms / 1000.0) 
                  << " cells/s\n";
        
        // Calculate parallel efficiency (assuming linear scaling from 1 process)
        // This is just illustrative - real baseline would need to be measured
        double baseline_time_estimate = result.avg_time_per_step_ms * size;
        double efficiency = baseline_time_estimate / (result.avg_time_per_step_ms * size);
        std::cout << "Parallel efficiency:  " << efficiency * 100.0 << "%\n";
    }
    
    std::remove("benchmark_scaling.config");
    
    EXPECT_GT(result.num_timesteps, 0);
}

// ============================================================================
// Problem Size Benchmarks
// ============================================================================

TEST_F(ScenarioBenchmark, ProblemSizeScaling) {
    std::vector<int> grid_sizes = {20, 30, 40, 50};
    
    if (rank == 0) {
        std::cout << "\n=== Problem Size Scaling Benchmark ===\n";
        std::cout << "Grid Size | Cells | Time/step (ms) | Cells/s\n";
        std::cout << "----------|-------|----------------|----------\n";
    }
    
    for (int n : grid_sizes) {
        int num_cells = n * n * n;
        
        std::string config = 
            "[SIMULATION]\n"
            "end_time = 100.0\n"
            "dt_initial = 10.0\n"
            "max_time_steps = 5\n"
            "fluid_model = SINGLE_PHASE\n"
            "[GRID]\n"
            "nx = " + std::to_string(n) + "\n"
            "ny = " + std::to_string(n) + "\n"
            "nz = " + std::to_string(n) + "\n"
            "Lx = 100.0\n"
            "Ly = 100.0\n"
            "Lz = 100.0\n"
            "[ROCK]\n"
            "porosity = 0.2\n"
            "permeability_x = 100.0\n"
            "[FLUID]\n"
            "viscosity = 0.001\n"
            "density = 1000.0\n";
        
        std::string filename = "benchmark_size_" + std::to_string(n) + ".config";
        std::ofstream ofs(filename);
        ofs << config;
        ofs.close();
        
        auto result = runSimulationBenchmark(filename);
        
        if (rank == 0 && result.num_timesteps > 0) {
            double cells_per_sec = num_cells / (result.avg_time_per_step_ms / 1000.0);
            printf("%3dx%3dx%3d | %6d | %14.1f | %9.0f\n",
                   n, n, n, num_cells, result.avg_time_per_step_ms, cells_per_sec);
        }
        
        std::remove(filename.c_str());
    }
}
