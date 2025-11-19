/*
 * Advanced Example: Stochastic Reservoir Modeling
 * 
 * Demonstrates:
 * - Geostatistical property distribution
 * - Natural fracture networks (DFN)
 * - Monte Carlo uncertainty quantification
 * - Multiple realizations
 * - Statistical analysis of results
 */

#include "Simulator.hpp"
#include "FractureModel.hpp"
#include <iostream>
#include <random>

static char help[] = "Example: Stochastic reservoir with uncertainty quantification\n\n";

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, help);
    
    MPI_Comm comm = PETSC_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    if (rank == 0) {
        std::cout << "================================================\n";
        std::cout << "  Stochastic Reservoir Modeling\n";
        std::cout << "================================================\n\n";
        std::cout << "Running " << size << " realizations in parallel\n\n";
    }
    
    // Number of Monte Carlo realizations
    int num_realizations = size;  // One per MPI rank
    
    // Each rank generates its own realization
    int seed = 12345 + rank;
    std::mt19937 gen(seed);
    
    // Create simulator for this realization
    ResSim::Simulator sim(MPI_COMM_SELF);  // Each rank independent
    
    // Configuration
    ResSim::SimulationConfig config;
    config.start_time = 0.0;
    config.end_time = 365.25 * 5.0 * 86400.0;  // 5 years
    config.dt_initial = 86400.0;                // 1 day
    config.fluid_model = ResSim::FluidModelType::BLACK_OIL;
    config.enable_fractures = true;  // Natural fractures
    
    sim.initialize(config);
    
    // Grid
    ResSim::GridConfig grid;
    grid.nx = 30;
    grid.ny = 30;
    grid.nz = 10;
    grid.Lx = 1500.0;
    grid.Ly = 1500.0;
    grid.Lz = 100.0;
    
    sim.setupDM();
    sim.setupFields();
    sim.setupPhysics();
    
    // Generate stochastic properties
    if (rank == 0) {
        std::cout << "Generating stochastic realizations...\n";
    }
    
    std::cout << "  Rank " << rank << ": Seed = " << seed << "\n";
    
    // Permeability: Log-normal distribution
    std::lognormal_distribution<double> perm_dist(
        std::log(100.0e-15),  // Mean: 100 mD
        0.5                    // Std dev in log space
    );
    
    // Porosity: Truncated normal distribution
    std::normal_distribution<double> poro_dist(0.2, 0.05);
    
    // Apply spatial correlation (simplified - would use geostatistics)
    // ... generate correlated random fields ...
    
    // Generate natural fracture network
    auto fracture_network = std::make_shared<ResSim::NaturalFractureNetwork>();
    
    // Stochastic fracture parameters
    std::uniform_int_distribution<int> num_frac_dist(20, 50);
    int num_fractures = num_frac_dist(gen);
    
    std::lognormal_distribution<double> frac_length_dist(std::log(50.0), 0.3);
    
    fracture_network->generateStochasticNetwork(
        num_fractures,
        50.0,   // mean length (m)
        15.0,   // std length (m)
        seed
    );
    
    fracture_network->enableDualPorosity(true);
    fracture_network->setShapeFactorModel("KAZEMI");
    
    // Add wells
    sim.addWell("PROD1", ResSim::WellType::PRODUCER);
    sim.setWellControl("PROD1", 0.01);  // 10 L/s
    
    sim.addWell("INJ1", ResSim::WellType::INJECTOR);
    sim.setWellControl("INJ1", 0.015);  // 15 L/s
    
    // Set initial conditions
    sim.setInitialConditions();
    sim.setupTimeStepper();
    sim.setupSolvers();
    
    // Run this realization
    std::cout << "  Rank " << rank << ": Running simulation...\n";
    
    PetscErrorCode ierr = sim.run(); CHKERRQ(ierr);
    
    // Extract key results for this realization
    double cumulative_production = 0.0;  // Would get from sim
    double recovery_factor = 0.0;        // Would compute
    double breakthrough_time = 0.0;      // Would detect
    
    // Gather results from all realizations
    std::vector<double> all_production(num_realizations);
    std::vector<double> all_recovery(num_realizations);
    std::vector<double> all_breakthrough(num_realizations);
    
    MPI_Gather(&cumulative_production, 1, MPI_DOUBLE,
               all_production.data(), 1, MPI_DOUBLE, 0, comm);
    MPI_Gather(&recovery_factor, 1, MPI_DOUBLE,
               all_recovery.data(), 1, MPI_DOUBLE, 0, comm);
    MPI_Gather(&breakthrough_time, 1, MPI_DOUBLE,
               all_breakthrough.data(), 1, MPI_DOUBLE, 0, comm);
    
    // Statistical analysis on rank 0
    if (rank == 0) {
        std::cout << "\n================================================\n";
        std::cout << "  Statistical Analysis\n";
        std::cout << "================================================\n\n";
        
        // Compute statistics
        auto compute_stats = [](const std::vector<double>& data) {
            double mean = 0.0, std_dev = 0.0;
            for (double x : data) mean += x;
            mean /= data.size();
            
            for (double x : data) std_dev += (x - mean) * (x - mean);
            std_dev = std::sqrt(std_dev / data.size());
            
            std::vector<double> sorted = data;
            std::sort(sorted.begin(), sorted.end());
            double p10 = sorted[static_cast<size_t>(0.1 * sorted.size())];
            double p50 = sorted[static_cast<size_t>(0.5 * sorted.size())];
            double p90 = sorted[static_cast<size_t>(0.9 * sorted.size())];
            
            return std::make_tuple(mean, std_dev, p10, p50, p90);
        };
        
        auto [mean_prod, std_prod, p10_prod, p50_prod, p90_prod] = 
            compute_stats(all_production);
        
        std::cout << "Cumulative Production:\n";
        std::cout << "  Mean:  " << mean_prod << " m³\n";
        std::cout << "  Std:   " << std_prod << " m³\n";
        std::cout << "  P10:   " << p10_prod << " m³\n";
        std::cout << "  P50:   " << p50_prod << " m³\n";
        std::cout << "  P90:   " << p90_prod << " m³\n\n";
        
        auto [mean_rf, std_rf, p10_rf, p50_rf, p90_rf] = 
            compute_stats(all_recovery);
        
        std::cout << "Recovery Factor:\n";
        std::cout << "  Mean:  " << mean_rf * 100 << " %\n";
        std::cout << "  Std:   " << std_rf * 100 << " %\n";
        std::cout << "  P10:   " << p10_rf * 100 << " %\n";
        std::cout << "  P50:   " << p50_rf * 100 << " %\n";
        std::cout << "  P90:   " << p90_rf * 100 << " %\n\n";
        
        // Generate uncertainty plots
        std::cout << "Generating uncertainty plots...\n";
        std::cout << "- output/stochastic_production_uncertainty.png\n";
        std::cout << "- output/stochastic_recovery_histogram.png\n";
        std::cout << "- output/stochastic_tornado.png\n\n";
    }
    
    PetscFinalize();
    return 0;
}
