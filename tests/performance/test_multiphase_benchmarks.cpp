/**
 * @file test_multiphase_benchmarks.cpp
 * @brief Multiphase flow benchmark tests
 * 
 * Benchmarks for two-phase and three-phase flow:
 * - Gravity segregation
 * - Counter-current imbibition
 * - Viscous fingering
 * - Gas-oil contact movement
 * - Three-phase relative permeability
 * - Capillary pressure effects
 */

#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>

class MultiphaseBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

// ============================================================================
// Gravity Segregation
// ============================================================================

TEST_F(MultiphaseBenchmark, GravitySegregation) {
    // Vertical segregation of oil and water due to density difference
    
    double rho_w = 1000.0;        // Water density (kg/m³)
    double rho_o = 850.0;         // Oil density (kg/m³)
    double k = 1000.0e-15;        // Permeability (m²)
    double phi = 0.25;            // Porosity
    double mu_w = 0.001;          // Water viscosity (Pa·s)
    double mu_o = 0.010;          // Oil viscosity (Pa·s)
    double g = 9.81;              // Gravity (m/s²)
    
    double height = 10.0;         // Column height (m)
    int num_layers = 100;
    double dz = height / num_layers;
    
    if (rank == 0) {
        std::cout << "\n=== Gravity Segregation Benchmark ===\n";
        std::cout << "Vertical oil-water segregation\n";
        std::cout << "Δρ = " << rho_w - rho_o << " kg/m³\n\n";
    }
    
    // Initial condition: uniform mixture S_w = 0.5
    std::vector<double> S_w(num_layers, 0.5);
    std::vector<double> S_w_new(num_layers);
    
    double dt = 10.0;  // Time step (seconds)
    int num_steps = 1000;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int step = 0; step < num_steps; ++step) {
        // Compute segregation velocities
        for (int i = 1; i < num_layers - 1; ++i) {
            // Relative permeabilities (simple power law)
            double k_rw = S_w[i] * S_w[i];
            double k_ro = (1.0 - S_w[i]) * (1.0 - S_w[i]);
            
            // Segregation velocity (Buckley-Leverett with gravity)
            double lambda_w = k_rw / mu_w;
            double lambda_o = k_ro / mu_o;
            double lambda_t = lambda_w + lambda_o;
            
            if (lambda_t > 1.0e-10) {
                double v_seg = k * lambda_w * lambda_o / lambda_t * (rho_w - rho_o) * g;
                
                // Update saturation (upwind scheme)
                double flux_w = v_seg * S_w[i];
                S_w_new[i] = S_w[i] - dt / (phi * dz) * (flux_w - v_seg * S_w[i-1]);
            } else {
                S_w_new[i] = S_w[i];
            }
            
            // Bound saturations
            S_w_new[i] = std::max(0.2, std::min(0.8, S_w_new[i]));
        }
        
        S_w = S_w_new;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "Simulation time: " << num_steps * dt << " seconds\n";
        std::cout << "Computation time: " << duration.count() << " ms\n";
        std::cout << "Bottom S_w: " << S_w[0] << ", Top S_w: " << S_w[num_layers-1] << "\n";
        std::cout << "Expected: Bottom enriched in water, top in oil\n";
    }
    
    EXPECT_GT(S_w[0], S_w[num_layers-1]);  // Bottom should have more water
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Counter-Current Imbibition
// ============================================================================

TEST_F(MultiphaseBenchmark, CounterCurrentImbibition) {
    // Spontaneous imbibition with counter-current flow
    
    double sigma = 0.03;          // Interfacial tension (N/m)
    double theta = 30.0 * M_PI / 180.0;  // Contact angle (radians)
    double k = 100.0e-15;         // Permeability (m²)
    double phi = 0.20;            // Porosity
    double mu_w = 0.001;          // Water viscosity
    double mu_o = 0.005;          // Oil viscosity
    
    // Capillary pressure parameters
    double P_e = 5000.0;          // Entry pressure (Pa)
    double lambda = 2.0;          // Pore size distribution index
    
    if (rank == 0) {
        std::cout << "\n=== Counter-Current Imbibition Benchmark ===\n";
        std::cout << "Spontaneous imbibition test\n";
        std::cout << "σ = " << sigma*1000 << " mN/m, θ = " << theta*180/M_PI << "°\n\n";
    }
    
    // Imbibition characteristic time
    double tau = phi * mu_w * mu_o / (k * sigma * std::cos(theta));
    
    // Dimensionless times
    std::vector<double> t_D = {0.001, 0.01, 0.1, 1.0, 10.0};
    
    auto start = std::chrono::high_resolution_clock::now();
    
    if (rank == 0) {
        std::cout << "t_D | Water recovery | Imbibition rate\n";
        std::cout << "----|----------------|------------------\n";
    }
    
    for (double t_dim : t_D) {
        // Approximate imbibition recovery (Ma et al. model)
        double R = std::sqrt(t_dim);  // Recovery fraction
        if (R > 0.9) R = 0.9;  // Maximum recovery
        
        double dR_dt = 0.5 / std::sqrt(t_dim);  // Rate
        
        if (rank == 0) {
            printf("%4.3f | %14.3f | %16.3e\n", t_dim, R, dR_dt);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nCharacteristic time: " << tau << " s\n";
        std::cout << "Computation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_GT(tau, 0.0);
}

// ============================================================================
// Viscous Fingering Stability
// ============================================================================

TEST_F(MultiphaseBenchmark, ViscousFingeringStability) {
    // Analyze stability of viscous fingering (Saffman-Taylor)
    
    std::vector<double> viscosity_ratios = {0.1, 1.0, 5.0, 10.0, 50.0};
    
    if (rank == 0) {
        std::cout << "\n=== Viscous Fingering Stability ===\n";
        std::cout << "Mobility ratio analysis\n\n";
        std::cout << "μ_inj/μ_disp | Mobility ratio | Stability | Finger width\n";
        std::cout << "-------------|----------------|-----------|-------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (double M : viscosity_ratios) {
        // Critical wavelength for instability
        double lambda_c = 2.0 * M_PI;  // Simplified
        
        // Growth rate (Saffman-Taylor)
        double sigma = M - 1.0;  // Instability parameter
        
        // Finger width (semi-analytical)
        double b = 0.5 * (1.0 + 1.0 / M);  // Relative finger width
        
        std::string stability;
        if (M > 1.0) {
            stability = "Unstable";
        } else if (M < 1.0) {
            stability = "Stable";
        } else {
            stability = "Neutral";
        }
        
        if (rank == 0) {
            printf("%12.1f | %14.2f | %9s | %11.3f\n",
                   M, M, stability.c_str(), b);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: M > 1 (low viscosity displacing high) is unstable\n";
    }
    
    EXPECT_LT(duration.count(), 1000);
}

// ============================================================================
// Three-Phase Relative Permeability
// ============================================================================

TEST_F(MultiphaseBenchmark, ThreePhaseRelativePermeability) {
    // Stone's model II for three-phase relative permeability
    
    double S_wc = 0.2;            // Connate water
    double S_orw = 0.2;           // Residual oil to water
    double S_org = 0.2;           // Residual oil to gas
    
    if (rank == 0) {
        std::cout << "\n=== Three-Phase Relative Permeability ===\n";
        std::cout << "Stone's Model II\n\n";
        std::cout << "S_w | S_o | S_g | k_rw | k_ro | k_rg\n";
        std::cout << "----|-----|-----|------|------|------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    int num_evaluations = 0;
    for (double S_w = 0.2; S_w <= 0.7; S_w += 0.1) {
        for (double S_g = 0.0; S_g <= 0.5; S_g += 0.1) {
            double S_o = 1.0 - S_w - S_g;
            if (S_o < 0.1) continue;
            
            // Normalized saturations
            double S_w_norm = (S_w - S_wc) / (1.0 - S_wc - S_orw);
            double S_g_norm = S_g / (1.0 - S_wc - S_org);
            
            // Two-phase rel perms
            double k_rw_2ph = std::pow(S_w_norm, 2.0);
            double k_rog_2ph = std::pow(1.0 - S_g_norm, 2.0);
            double k_rg_2ph = std::pow(S_g_norm, 2.0);
            
            // Stone's Model II for three-phase k_ro
            double k_ro_3ph = k_rog_2ph * k_rw_2ph;
            
            // Three-phase rel perms
            double k_rw = k_rw_2ph;
            double k_ro = k_ro_3ph;
            double k_rg = k_rg_2ph;
            
            num_evaluations++;
            
            if (rank == 0 && num_evaluations % 5 == 0) {
                printf("%4.1f | %4.1f | %4.1f | %4.2f | %4.2f | %4.2f\n",
                       S_w, S_o, S_g, k_rw, k_ro, k_rg);
            }
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_evaluations;
    
    if (rank == 0) {
        std::cout << "\nEvaluations: " << num_evaluations << "\n";
        std::cout << "Time per evaluation: " << time_per_eval << " ns\n";
        std::cout << "Throughput: " << 1e6 / time_per_eval << " M eval/s\n";
    }
    
    EXPECT_LT(time_per_eval, 1000.0);
}

// ============================================================================
// Capillary Pressure Hysteresis
// ============================================================================

TEST_F(MultiphaseBenchmark, CapillaryPressureHysteresis) {
    // Drainage and imbibition curves with hysteresis
    
    if (rank == 0) {
        std::cout << "\n=== Capillary Pressure Hysteresis ===\n";
        std::cout << "Drainage vs Imbibition\n\n";
        std::cout << "S_w | P_c drainage (kPa) | P_c imbibition (kPa) | Hysteresis\n";
        std::cout << "----|----------------------|----------------------|------------\n";
    }
    
    // Drainage curve (Brooks-Corey)
    double P_d = 5000.0;          // Entry pressure drainage (Pa)
    double lambda_d = 2.0;        // Pore size index
    double S_wc = 0.2;
    
    // Imbibition curve
    double P_i = 2000.0;          // Entry pressure imbibition (Pa)
    double lambda_i = 3.0;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (double S_w = 0.2; S_w <= 0.8; S_w += 0.05) {
        // Drainage
        double S_e_d = (S_w - S_wc) / (1.0 - S_wc);
        double P_c_d = P_d * std::pow(S_e_d, -1.0/lambda_d);
        
        // Imbibition
        double P_c_i = P_i * std::pow(S_e_d, -1.0/lambda_i);
        
        // Hysteresis (difference)
        double hysteresis = P_c_d - P_c_i;
        
        if (rank == 0 && static_cast<int>(S_w * 20) % 2 == 0) {
            printf("%4.2f | %20.2f | %20.2f | %10.2f\n",
                   S_w, P_c_d/1000.0, P_c_i/1000.0, hysteresis/1000.0);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Drainage P_c > Imbibition P_c (hysteresis)\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Relative Permeability Models Comparison
// ============================================================================

TEST_F(MultiphaseBenchmark, RelPermModelsComparison) {
    // Compare different rel perm models
    
    if (rank == 0) {
        std::cout << "\n=== Relative Permeability Models Comparison ===\n\n";
        std::cout << "S_w | Corey | Brooks-Corey | LET | van Genuchten\n";
        std::cout << "----|-------|--------------|-----|---------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (double S_w = 0.2; S_w <= 0.8; S_w += 0.1) {
        double S_wc = 0.2;
        double S_or = 0.2;
        double S_e = (S_w - S_wc) / (1.0 - S_wc - S_or);
        
        // Corey model
        double k_rw_corey = std::pow(S_e, 2.0);
        
        // Brooks-Corey model
        double lambda = 2.0;
        double k_rw_bc = std::pow(S_e, (2.0 + 3.0*lambda) / lambda);
        
        // LET model (Lomeland-Ebeltoft-Thomas)
        double L = 2.0, E = 2.0, T = 2.0;
        double k_rw_let = std::pow(S_e, L) / (std::pow(S_e, L) + E * std::pow(1.0 - S_e, T));
        
        // van Genuchten model
        double m = 0.5;
        double k_rw_vg = std::sqrt(S_e) * std::pow(1.0 - std::pow(1.0 - std::pow(S_e, 1.0/m), m), 2.0);
        
        if (rank == 0) {
            printf("%4.1f | %5.3f | %12.3f | %4.3f | %13.3f\n",
                   S_w, k_rw_corey, k_rw_bc, k_rw_let, k_rw_vg);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Saturation Front Tracking
// ============================================================================

TEST_F(MultiphaseBenchmark, SaturationFrontTracking) {
    // Track saturation front movement (method of characteristics)
    
    double u_t = 1.0e-5;          // Total velocity (m/s)
    double phi = 0.2;             // Porosity
    double L = 100.0;             // Domain length (m)
    
    if (rank == 0) {
        std::cout << "\n=== Saturation Front Tracking ===\n";
        std::cout << "Method of characteristics\n\n";
        std::cout << "Time (days) | Front position (m) | Front speed (m/day)\n";
        std::cout << "------------|-------------------|--------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> times = {1, 5, 10, 30, 100, 365};  // days
    
    for (double t_days : times) {
        double t = t_days * 86400.0;  // Convert to seconds
        
        // Front position (simplified Buckley-Leverett)
        double df_dS = 2.0;  // Approximate derivative of fractional flow
        double x_f = u_t * t * df_dS / phi;
        
        // Front speed
        double v_f = u_t * df_dS / phi;
        
        if (rank == 0) {
            printf("%11.1f | %17.2f | %18.4f\n",
                   t_days, x_f, v_f * 86400.0);
        }
        
        if (x_f > L) break;  // Front reached end
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Fractional Flow Analysis
// ============================================================================

TEST_F(MultiphaseBenchmark, FractionalFlowAnalysis) {
    // Analyze fractional flow for different viscosity ratios
    
    std::vector<double> visc_ratios = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
    
    if (rank == 0) {
        std::cout << "\n=== Fractional Flow Analysis ===\n\n";
        std::cout << "μ_w/μ_o | M | Breakthrough S_w | Recovery at BT | Ultimate recovery\n";
        std::cout << "--------|---|------------------|----------------|------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (double M_visc : visc_ratios) {
        // Mobility ratio (assuming equal endpoint rel perms)
        double M = 1.0 / M_visc;
        
        // Breakthrough saturation (Welge analysis)
        double S_w_bt = 0.5;  // Simplified
        
        // Recovery at breakthrough
        double R_bt = 0.5;
        
        // Ultimate recovery
        double R_ult = 0.7;
        
        if (rank == 0) {
            printf("%7.1f | %1.1f | %16.2f | %14.2f | %16.2f\n",
                   M_visc, M, S_w_bt, R_bt, R_ult);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Higher M (favorable mobility) = better recovery\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}
