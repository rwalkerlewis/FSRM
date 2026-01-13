/**
 * @file test_thermal_eor_benchmarks.cpp
 * @brief Thermal and Enhanced Oil Recovery (EOR) benchmarks
 * 
 * Benchmarks for:
 * - Steam injection (SAGD, cyclic steam)
 * - Hot water flooding
 * - In-situ combustion
 * - Polymer flooding
 * - Surfactant flooding
 * - CO2 flooding
 */

#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <vector>

class ThermalEORBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

// ============================================================================
// Steam Flooding - Marx-Langenheim Model
// ============================================================================

TEST_F(ThermalEORBenchmark, SteamFloodingMarxLangenheim) {
    // Marx-Langenheim analytical solution for steam zone growth
    
    double T_s = 250.0 + 273.15;  // Steam temperature (K)
    double T_r = 50.0 + 273.15;   // Reservoir temperature (K)
    double h = 10.0;              // Reservoir thickness (m)
    double M_s = 100.0;           // Steam injection rate (kg/s)
    double L_s = 2.26e6;          // Latent heat of steam (J/kg)
    double rho_r = 2200.0;        // Rock density (kg/m³)
    double c_r = 1000.0;          // Rock heat capacity (J/kg·K)
    double phi = 0.3;             // Porosity
    double rho_f = 800.0;         // Fluid density (kg/m³)
    double c_f = 2000.0;          // Fluid heat capacity (J/kg·K)
    
    // Volumetric heat capacity
    double M_v = rho_r * c_r * (1.0 - phi) + rho_f * c_f * phi;
    
    // Thermal diffusivity of overburden
    double alpha = 1.0e-6;        // m²/s
    
    if (rank == 0) {
        std::cout << "\n=== Steam Flooding (Marx-Langenheim) ===\n";
        std::cout << "T_steam = " << T_s - 273.15 << " °C\n";
        std::cout << "Injection rate = " << M_s << " kg/s\n\n";
        std::cout << "Time (days) | Steam zone area (m²) | Heat efficiency\n";
        std::cout << "------------|----------------------|----------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> times = {10, 30, 100, 365, 1000};  // days
    
    for (double t_days : times) {
        double t = t_days * 86400.0;  // Convert to seconds
        
        // Dimensionless time
        double t_D = 4.0 * alpha * t / (h * h);
        
        // Heat loss function (complementary error function)
        double G_t = std::sqrt(t_D / M_PI) * std::exp(-1.0 / (4.0 * t_D)) - 
                     0.5 * std::erfc(1.0 / (2.0 * std::sqrt(t_D)));
        
        // Heat efficiency
        double E_h = 1.0 - G_t;
        
        // Steam zone volume
        double Q_inj = M_s * L_s * t * E_h;  // Net heat injected
        double V_s = Q_inj / (M_v * (T_s - T_r));
        
        // Steam zone area (assuming uniform height h)
        double A_s = V_s / h;
        
        if (rank == 0) {
            printf("%11.0f | %20.1f | %14.3f\n", t_days, A_s, E_h);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Heat efficiency decreases with time due to losses\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// SAGD (Steam-Assisted Gravity Drainage)
// ============================================================================

TEST_F(ThermalEORBenchmark, SAGDPerformance) {
    // SAGD production rate estimation
    
    double k = 5000.0e-15;        // Permeability (m²)
    double phi = 0.35;            // Porosity
    double h = 30.0;              // Pay zone thickness (m)
    double mu_o_hot = 0.001;      // Hot oil viscosity (Pa·s)
    double mu_o_cold = 1.0;       // Cold oil viscosity (Pa·s)
    double Delta_rho = 150.0;     // Density difference (kg/m³)
    double g = 9.81;              // Gravity (m/s²)
    double L = 1000.0;            // Well pair length (m)
    
    if (rank == 0) {
        std::cout << "\n=== SAGD Performance Benchmark ===\n";
        std::cout << "Steam-Assisted Gravity Drainage\n";
        std::cout << "Permeability = " << k*1e15 << " mD\n\n";
        std::cout << "Chamber height (m) | Production rate (m³/day) | Oil rate (bbl/day)\n";
        std::cout << "-------------------|--------------------------|--------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Chamber growth over time
    for (double h_c = 5.0; h_c <= h; h_c += 5.0) {
        // Butler's theory: drainage rate
        double q = k * Delta_rho * g * h_c * h_c * L / (2.0 * mu_o_hot);
        
        // Convert to m³/day
        double q_day = q * 86400.0;
        
        // Convert to bbl/day
        double q_bbl = q_day * 6.289;  // 1 m³ = 6.289 bbl
        
        if (rank == 0) {
            printf("%18.1f | %24.2f | %18.1f\n", h_c, q_day, q_bbl);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Production increases with chamber height\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}

// ============================================================================
// Polymer Flooding
// ============================================================================

TEST_F(ThermalEORBenchmark, PolymerFloodingViscosity) {
    // Polymer solution viscosity models
    
    if (rank == 0) {
        std::cout << "\n=== Polymer Flooding Viscosity ===\n\n";
        std::cout << "Conc (ppm) | Flory-Huggins | Carreau | Power-law | Shear rate effect\n";
        std::cout << "-----------|---------------|---------|-----------|------------------\n";
    }
    
    double mu_w = 0.001;          // Water viscosity (Pa·s)
    double gamma_dot = 10.0;      // Shear rate (1/s)
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> concentrations = {100, 500, 1000, 2000, 5000};  // ppm
    
    for (double C : concentrations) {
        // Flory-Huggins model (dilute solution)
        double intrinsic_visc = 2.5;  // dL/g
        double mu_FH = mu_w * (1.0 + intrinsic_visc * C / 1000.0);
        
        // Carreau model (shear-thinning)
        double lambda = 1.0;  // Time constant
        double n = 0.5;       // Power-law index
        double mu_inf = mu_w;
        double mu_0 = mu_FH;
        double mu_Carreau = mu_inf + (mu_0 - mu_inf) * 
                           std::pow(1.0 + std::pow(lambda * gamma_dot, 2.0), (n - 1.0) / 2.0);
        
        // Power-law model
        double K = 0.01;      // Consistency index
        double mu_PL = K * std::pow(gamma_dot, n - 1.0);
        
        // Shear rate effect factor
        double shear_factor = mu_Carreau / mu_0;
        
        if (rank == 0) {
            printf("%10.0f | %13.4f | %7.4f | %9.4f | %16.3f\n",
                   C, mu_FH*1000, mu_Carreau*1000, mu_PL*1000, shear_factor);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nViscosity in cP\n";
        std::cout << "Computation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// Surfactant Flooding - IFT Reduction
// ============================================================================

TEST_F(ThermalEORBenchmark, SurfactantIFTReduction) {
    // Interfacial tension reduction with surfactant
    
    double sigma_0 = 0.030;       // Initial IFT (N/m)
    double CMC = 100.0;           // Critical micelle concentration (ppm)
    
    if (rank == 0) {
        std::cout << "\n=== Surfactant IFT Reduction ===\n";
        std::cout << "Initial IFT = " << sigma_0*1000 << " mN/m\n\n";
        std::cout << "Conc (ppm) | IFT (mN/m) | Capillary number | Residual oil mobilization\n";
        std::cout << "-----------|------------|------------------|-------------------------\n";
    }
    
    double u = 1.0e-5;            // Darcy velocity (m/s)
    double mu = 0.001;            // Viscosity (Pa·s)
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> concentrations = {0, 50, 100, 500, 1000, 5000};
    
    for (double C : concentrations) {
        // IFT reduction (log-linear below CMC, plateau above)
        double sigma;
        if (C < CMC) {
            sigma = sigma_0 * std::exp(-5.0 * C / CMC);
        } else {
            sigma = sigma_0 * 1.0e-3;  // Ultra-low IFT
        }
        
        // Capillary number
        double N_c = mu * u / sigma;
        
        // Residual oil mobilization (correlation)
        double S_or_initial = 0.30;
        double S_or_reduced = S_or_initial / (1.0 + 10.0 * N_c);
        double mobilization = (S_or_initial - S_or_reduced) / S_or_initial;
        
        if (rank == 0) {
            printf("%10.0f | %10.3f | %16.2e | %23.2f%%\n",
                   C, sigma*1000, N_c, mobilization*100);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Ultra-low IFT needed to mobilize residual oil\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// CO2 Flooding - Minimum Miscibility Pressure
// ============================================================================

TEST_F(ThermalEORBenchmark, CO2MiscibilityPressure) {
    // Estimate minimum miscibility pressure (MMP)
    
    if (rank == 0) {
        std::cout << "\n=== CO2 Miscibility Pressure ===\n\n";
        std::cout << "Temperature (°C) | Reservoir °API | MMP (psi) | MMP (MPa)\n";
        std::cout << "-----------------|----------------|-----------|----------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> temperatures = {30, 50, 70, 90, 110};
    std::vector<double> api_gravities = {25, 30, 35, 40};
    
    for (double T_C : temperatures) {
        for (double API : api_gravities) {
            // Correlation for MMP (Alston et al.)
            // MMP (psi) = a * T + b * API + c
            double a = 18.0;
            double b = -5.0;
            double c = 1200.0;
            
            double MMP_psi = a * T_C + b * API + c;
            double MMP_MPa = MMP_psi * 0.00689476;  // psi to MPa
            
            if (rank == 0 && static_cast<int>(T_C) % 40 == 30) {
                printf("%16.0f | %14.0f | %9.0f | %8.1f\n",
                       T_C, API, MMP_psi, MMP_MPa);
            }
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: MMP increases with temperature, decreases with API\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Thermal Conduction in Porous Media
// ============================================================================

TEST_F(ThermalEORBenchmark, ThermalConductivityEffective) {
    // Effective thermal conductivity models
    
    if (rank == 0) {
        std::cout << "\n=== Effective Thermal Conductivity ===\n\n";
        std::cout << "Porosity | Arithmetic | Geometric | Harmonic | Maxwell\n";
        std::cout << "---------|------------|-----------|----------|---------\n";
    }
    
    double k_s = 2.5;             // Solid conductivity (W/m·K)
    double k_f = 0.6;             // Fluid conductivity (W/m·K)
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (double phi = 0.1; phi <= 0.4; phi += 0.05) {
        // Arithmetic mean (parallel)
        double k_arith = phi * k_f + (1.0 - phi) * k_s;
        
        // Geometric mean
        double k_geom = std::pow(k_f, phi) * std::pow(k_s, 1.0 - phi);
        
        // Harmonic mean (series)
        double k_harm = 1.0 / (phi / k_f + (1.0 - phi) / k_s);
        
        // Maxwell model (spherical inclusions)
        double k_maxwell = k_s * (2.0 * k_s + k_f - 2.0 * phi * (k_s - k_f)) /
                                 (2.0 * k_s + k_f + phi * (k_s - k_f));
        
        if (rank == 0) {
            printf("%8.2f | %10.3f | %9.3f | %8.3f | %7.3f\n",
                   phi, k_arith, k_geom, k_harm, k_maxwell);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nConductivity in W/m·K\n";
        std::cout << "Computation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}

// ============================================================================
// Cyclic Steam Stimulation (CSS)
// ============================================================================

TEST_F(ThermalEORBenchmark, CyclicSteamStimulation) {
    // CSS cycle performance
    
    double V_inj = 1000.0;        // Steam injection volume (m³)
    double t_soak = 5.0;          // Soak time (days)
    double h = 20.0;              // Pay zone thickness (m)
    double phi = 0.35;            // Porosity
    double S_o = 0.80;            // Oil saturation
    double B_o = 1.2;             // Oil formation volume factor
    
    if (rank == 0) {
        std::cout << "\n=== Cyclic Steam Stimulation ===\n";
        std::cout << "Huff-and-puff cycles\n\n";
        std::cout << "Cycle | Steam (m³) | Soak (days) | Oil produced (m³) | SOR\n";
        std::cout << "------|------------|-------------|-------------------|-----\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    double total_steam = 0.0;
    double total_oil = 0.0;
    
    for (int cycle = 1; cycle <= 10; ++cycle) {
        // Heated zone radius (simplified)
        double r_h = 10.0 * std::sqrt(cycle);
        
        // Recoverable oil
        double V_heated = M_PI * r_h * r_h * h;
        double N_oil = V_heated * phi * S_o / B_o;
        
        // Production declines with cycles
        double recovery_factor = 0.3 * std::exp(-0.15 * cycle);
        double N_produced = N_oil * recovery_factor;
        
        total_steam += V_inj;
        total_oil += N_produced;
        
        // Steam-oil ratio (SOR)
        double SOR = total_steam / total_oil;
        
        if (rank == 0) {
            printf("%5d | %10.0f | %11.0f | %17.1f | %4.1f\n",
                   cycle, V_inj, t_soak, N_produced, SOR);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nTotal steam: " << total_steam << " m³\n";
        std::cout << "Total oil: " << total_oil << " m³\n";
        std::cout << "Overall SOR: " << total_steam / total_oil << "\n";
        std::cout << "Computation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// In-Situ Combustion
// ============================================================================

TEST_F(ThermalEORBenchmark, InSituCombustion) {
    // In-situ combustion front propagation
    
    double m_fuel = 50.0;         // Fuel concentration (kg/m³)
    double H_c = 40.0e6;          // Heat of combustion (J/kg)
    double rho_r = 2200.0;        // Rock density (kg/m³)
    double c_r = 1000.0;          // Rock heat capacity (J/kg·K)
    double phi = 0.25;            // Porosity
    double air_flux = 50.0;       // Air flux (m³/m²·day)
    
    if (rank == 0) {
        std::cout << "\n=== In-Situ Combustion ===\n";
        std::cout << "Combustion front propagation\n\n";
        std::cout << "Time (days) | Front position (m) | Temperature (°C) | Air consumed\n";
        std::cout << "------------|-------------------|------------------|-------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Volumetric heat capacity
    double M_v = rho_r * c_r * (1.0 - phi);
    
    // Combustion temperature rise
    double Delta_T = m_fuel * H_c / M_v;
    double T_max = 500.0;  // Maximum combustion temperature (°C)
    
    std::vector<double> times = {10, 30, 100, 365};
    
    for (double t_days : times) {
        // Air requirement for combustion (stoichiometric)
        double air_per_fuel = 12.0;  // m³ air / kg fuel
        
        // Front velocity
        double v_f = air_flux / air_per_fuel / m_fuel;  // m/day
        
        // Front position
        double x_f = v_f * t_days;
        
        // Air consumed
        double air_consumed = air_flux * t_days;
        
        if (rank == 0) {
            printf("%11.0f | %17.2f | %16.0f | %11.0f\n",
                   t_days, x_f, T_max, air_consumed);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}
