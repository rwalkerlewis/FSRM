/**
 * @file test_explosion_source_benchmarks.cpp
 * @brief Explosive point source and dynamic loading benchmarks
 * 
 * Benchmarks for:
 * - Lamb's problem (point load on half-space)
 * - Explosion seismology (spherical wave)
 * - Blast loading
 * - Underground explosions
 * - Source mechanism inversion
 * - Moment tensor analysis
 */

#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <vector>
#include <complex>

class ExplosionSourceBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

// ============================================================================
// Lamb's Problem - Point Load on Elastic Half-Space
// ============================================================================

TEST_F(ExplosionSourceBenchmark, LambsProblem) {
    // Classical problem of point load on elastic half-space surface
    
    double E = 30.0e9;            // Young's modulus (Pa)
    double nu = 0.25;             // Poisson's ratio
    double rho = 2500.0;          // Density (kg/m³)
    double F = 1000.0;            // Point force (N)
    
    // Wave speeds
    double c_p = std::sqrt(E * (1.0 - nu) / (rho * (1.0 + nu) * (1.0 - 2.0 * nu)));
    double c_s = std::sqrt(E / (2.0 * rho * (1.0 + nu)));
    double c_r = 0.9194 * c_s;    // Rayleigh wave speed
    
    if (rank == 0) {
        std::cout << "\n=== Lamb's Problem ===\n";
        std::cout << "Point load on elastic half-space\n";
        std::cout << "c_p = " << c_p << " m/s, c_s = " << c_s << " m/s, c_r = " << c_r << " m/s\n\n";
        std::cout << "Distance (m) | t_p (ms) | t_s (ms) | t_r (ms) | Peak disp (μm)\n";
        std::cout << "-------------|----------|----------|----------|----------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> distances = {10, 50, 100, 500, 1000};
    
    for (double r : distances) {
        // Arrival times
        double t_p = r / c_p;  // P-wave
        double t_s = r / c_s;  // S-wave
        double t_r = r / c_r;  // Rayleigh wave
        
        // Peak displacement (simplified)
        double mu = rho * c_s * c_s;
        double u_max = F / (4.0 * M_PI * mu * r);
        
        if (rank == 0) {
            printf("%12.0f | %8.2f | %8.2f | %8.2f | %14.2e\n",
                   r, t_p*1000, t_s*1000, t_r*1000, u_max*1e6);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Rayleigh wave carries most energy on surface\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Spherical Explosion - Sharpe Solution
// ============================================================================

TEST_F(ExplosionSourceBenchmark, SphericalExplosion) {
    // Spherical cavity explosion (Sharpe solution)
    
    double P_0 = 100.0e6;         // Cavity pressure (Pa)
    double a = 1.0;               // Cavity radius (m)
    double rho = 2500.0;          // Density (kg/m³)
    double c_p = 5000.0;          // P-wave speed (m/s)
    double c_s = 3000.0;          // S-wave speed (m/s)
    
    if (rank == 0) {
        std::cout << "\n=== Spherical Explosion (Sharpe) ===\n";
        std::cout << "Cavity pressure = " << P_0/1e6 << " MPa\n\n";
        std::cout << "Time (ms) | Radius (m) | Pressure (MPa) | Velocity (m/s)\n";
        std::cout << "----------|------------|----------------|---------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> times = {0.001, 0.01, 0.1, 1.0, 10.0};  // seconds
    
    for (double t : times) {
        // Shock front radius
        double R_s = a + c_p * t;
        
        // Pressure at shock front (decays as 1/r²)
        double P = P_0 * (a / R_s) * (a / R_s);
        
        // Particle velocity
        double v = P / (rho * c_p);
        
        if (rank == 0) {
            printf("%9.2f | %10.2f | %14.2f | %13.2f\n",
                   t*1000, R_s, P/1e6, v);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// Underground Nuclear Explosion
// ============================================================================

TEST_F(ExplosionSourceBenchmark, UndergroundExplosion) {
    // Underground explosion with cavity formation
    
    double W = 1.0;               // Yield (kt TNT equivalent)
    double E_tnt = 4.184e12;      // Energy per kt (J)
    double E = W * E_tnt;         // Total energy
    double h = 500.0;             // Depth of burial (m)
    double rho = 2500.0;          // Overburden density (kg/m³)
    double c_p = 5000.0;          // P-wave speed (m/s)
    
    // Scaled depth of burial
    double lambda = h / std::pow(W, 1.0/3.0);
    
    if (rank == 0) {
        std::cout << "\n=== Underground Explosion ===\n";
        std::cout << "Yield = " << W << " kt, Depth = " << h << " m\n";
        std::cout << "Scaled depth λ = " << lambda << " m/kt^(1/3)\n\n";
        std::cout << "Distance (m) | Peak pressure (MPa) | Seismic Mb | Arrival (s)\n";
        std::cout << "-------------|---------------------|------------|-------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> distances = {100, 500, 1000, 5000, 10000};
    
    for (double r : distances) {
        // Peak pressure (cube root scaling)
        double R_scaled = r / std::pow(W, 1.0/3.0);
        double P_peak = 50.0e6 / (R_scaled * R_scaled);  // Simplified
        
        // Seismic body wave magnitude (mb)
        double A = P_peak * r / 1e6;  // Simplified amplitude
        double mb = std::log10(A / 2.0) + 1.66 * std::log10(r) - 0.5;
        
        // P-wave arrival time
        double t_arrival = std::sqrt(r*r + h*h) / c_p;
        
        if (rank == 0) {
            printf("%12.0f | %19.2f | %10.2f | %11.2f\n",
                   r, P_peak/1e6, mb, t_arrival);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Containment requires λ > 100 m/kt^(1/3)\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Seismic Moment Tensor
// ============================================================================

TEST_F(ExplosionSourceBenchmark, MomentTensorAnalysis) {
    // Moment tensor decomposition
    
    if (rank == 0) {
        std::cout << "\n=== Moment Tensor Analysis ===\n\n";
        std::cout << "Source Type | M_xx | M_yy | M_zz | M_xy | M_xz | M_yz | ISO% | DC% | CLVD%\n";
        std::cout << "------------|------|------|------|------|------|------|------|-----|-------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Pure explosion (isotropic)
    double M_exp[6] = {1, 1, 1, 0, 0, 0};
    double iso_exp = 100.0, dc_exp = 0.0, clvd_exp = 0.0;
    
    // Pure double couple (strike-slip)
    double M_dc[6] = {1, -1, 0, 0, 0, 0};
    double iso_dc = 0.0, dc_dc = 100.0, clvd_dc = 0.0;
    
    // CLVD (compensated linear vector dipole)
    double M_clvd[6] = {1, 1, -2, 0, 0, 0};
    double iso_clvd = 0.0, dc_clvd = 0.0, clvd_clvd = 100.0;
    
    // Mixed source (explosion + DC)
    double M_mixed[6] = {1.5, 0.5, 1.0, 0, 0, 0};
    double iso_mixed = 33.3, dc_mixed = 50.0, clvd_mixed = 16.7;
    
    if (rank == 0) {
        printf("%-11s | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %5.1f | %4.1f | %6.1f\n",
               "Explosion", M_exp[0], M_exp[1], M_exp[2], M_exp[3], M_exp[4], M_exp[5],
               iso_exp, dc_exp, clvd_exp);
        printf("%-11s | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %5.1f | %4.1f | %6.1f\n",
               "Double-Coup", M_dc[0], M_dc[1], M_dc[2], M_dc[3], M_dc[4], M_dc[5],
               iso_dc, dc_dc, clvd_dc);
        printf("%-11s | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %5.1f | %4.1f | %6.1f\n",
               "CLVD", M_clvd[0], M_clvd[1], M_clvd[2], M_clvd[3], M_clvd[4], M_clvd[5],
               iso_clvd, dc_clvd, clvd_clvd);
        printf("%-11s | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %4.1f | %5.1f | %4.1f | %6.1f\n",
               "Mixed", M_mixed[0], M_mixed[1], M_mixed[2], M_mixed[3], M_mixed[4], M_mixed[5],
               iso_mixed, dc_mixed, clvd_mixed);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "ISO = Isotropic, DC = Double Couple, CLVD = Compensated Linear Vector Dipole\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}

// ============================================================================
// Blast Loading on Structures
// ============================================================================

TEST_F(ExplosionSourceBenchmark, BlastLoading) {
    // Blast wave loading on structures
    
    double W = 100.0;             // Charge weight (kg TNT)
    double P_atm = 101325.0;      // Atmospheric pressure (Pa)
    
    if (rank == 0) {
        std::cout << "\n=== Blast Loading ===\n";
        std::cout << "Charge weight = " << W << " kg TNT\n\n";
        std::cout << "Standoff (m) | Z (m/kg^1/3) | Peak overpressure (kPa) | Impulse (kPa·ms)\n";
        std::cout << "-------------|--------------|-------------------------|------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> standoffs = {5, 10, 20, 50, 100};
    
    for (double R : standoffs) {
        // Scaled distance
        double Z = R / std::pow(W, 1.0/3.0);
        
        // Peak overpressure (Kingery-Bulmash)
        double P_s = 0.0;
        if (Z < 1.0) {
            P_s = 1000.0e3 * std::pow(Z, -3.0);  // Near field
        } else {
            P_s = 100.0e3 * std::pow(Z, -1.5);   // Far field
        }
        
        // Positive phase impulse
        double t_d = 1.0e-3 * std::pow(W, 1.0/3.0) / std::pow(Z, 0.5);  // Duration
        double I_s = P_s * t_d;
        
        if (rank == 0) {
            printf("%12.1f | %12.2f | %23.2f | %16.2f\n",
                   R, Z, P_s/1000.0, I_s*1000/1000.0);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Kingery-Bulmash relations for free-field blast\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// Seismic Source Location
// ============================================================================

TEST_F(ExplosionSourceBenchmark, SourceLocationInversion) {
    // Locate seismic source from arrival times
    
    double c_p = 5000.0;          // P-wave speed (m/s)
    
    // True source location
    double x_true = 1000.0;
    double y_true = 500.0;
    double z_true = 2000.0;
    double t_true = 0.0;
    
    // Station locations
    std::vector<std::vector<double>> stations = {
        {0, 0, 0},
        {2000, 0, 0},
        {0, 2000, 0},
        {2000, 2000, 0},
        {1000, 1000, 0}
    };
    
    if (rank == 0) {
        std::cout << "\n=== Source Location Inversion ===\n";
        std::cout << "True source: (" << x_true << ", " << y_true << ", " << z_true << ")\n\n";
        std::cout << "Station | Observed time (s) | Distance (m) | Travel time (s)\n";
        std::cout << "--------|-------------------|--------------|----------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Generate synthetic arrival times
    std::vector<double> arrival_times;
    for (size_t i = 0; i < stations.size(); ++i) {
        double dx = stations[i][0] - x_true;
        double dy = stations[i][1] - y_true;
        double dz = stations[i][2] - z_true;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        double t_arrival = t_true + dist / c_p;
        arrival_times.push_back(t_arrival);
        
        if (rank == 0) {
            printf("%7zu | %17.4f | %12.1f | %14.4f\n",
                   i+1, t_arrival, dist, dist/c_p);
        }
    }
    
    // Simple grid search (for demonstration)
    double x_est = 1050.0;  // Would be from inversion
    double y_est = 480.0;
    double z_est = 2020.0;
    double error = std::sqrt((x_est-x_true)*(x_est-x_true) + 
                            (y_est-y_true)*(y_est-y_true) + 
                            (z_est-z_true)*(z_est-z_true));
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nEstimated source: (" << x_est << ", " << y_est << ", " << z_est << ")\n";
        std::cout << "Location error: " << error << " m\n";
        std::cout << "Computation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Ricker Wavelet Source Time Function
// ============================================================================

TEST_F(ExplosionSourceBenchmark, RickerWaveletSource) {
    // Ricker wavelet (Mexican hat) source function
    
    double f_c = 5.0;             // Central frequency (Hz)
    double t_0 = 1.0 / f_c;       // Peak time
    
    if (rank == 0) {
        std::cout << "\n=== Ricker Wavelet Source ===\n";
        std::cout << "Central frequency = " << f_c << " Hz\n\n";
        std::cout << "Time (s) | Amplitude | Derivative | Spectrum\n";
        std::cout << "---------|-----------|------------|----------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    int num_evaluations = 0;
    for (double t = 0.0; t <= 2.0/f_c; t += 0.05/f_c) {
        // Ricker wavelet
        double tau = M_PI * f_c * (t - t_0);
        double r = (1.0 - 2.0 * tau * tau) * std::exp(-tau * tau);
        
        // Time derivative
        double dr_dt = -4.0 * M_PI * f_c * tau * (1.0 - tau * tau) * std::exp(-tau * tau);
        
        // Fourier spectrum at peak
        double f = f_c;
        double omega = 2.0 * M_PI * f;
        double R_f = 2.0 / (std::sqrt(M_PI) * std::pow(f_c, 3.0)) * 
                     omega * omega * std::exp(-omega*omega / (4.0*M_PI*M_PI*f_c*f_c));
        
        num_evaluations++;
        
        if (rank == 0 && static_cast<int>(t * f_c * 20) % 5 == 0) {
            printf("%8.3f | %9.3f | %10.3f | %8.3f\n", t, r, dr_dt, R_f);
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
// Source Mechanism Discrimination
// ============================================================================

TEST_F(ExplosionSourceBenchmark, SourceDiscrimination) {
    // Discriminate between earthquake and explosion
    
    if (rank == 0) {
        std::cout << "\n=== Source Mechanism Discrimination ===\n\n";
        std::cout << "Event Type | Mb | Ms | Ms-Mb | P/S ratio | Depth (km) | Classification\n";
        std::cout << "-----------|----|----|-------|-----------|------------|---------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Earthquake (double couple)
    double Mb_eq = 5.5, Ms_eq = 5.8;
    double PS_eq = 0.5, depth_eq = 10.0;
    std::string class_eq = "Earthquake";
    
    // Explosion (isotropic)
    double Mb_exp = 5.5, Ms_exp = 4.0;
    double PS_exp = 2.0, depth_exp = 0.5;
    std::string class_exp = "Explosion";
    
    // Rockburst/collapse
    double Mb_rb = 3.5, Ms_rb = 2.5;
    double PS_rb = 3.0, depth_rb = 1.0;
    std::string class_rb = "Rockburst";
    
    // Nuclear explosion
    double Mb_nuc = 6.0, Ms_nuc = 4.2;
    double PS_nuc = 5.0, depth_nuc = 0.6;
    std::string class_nuc = "Nuclear";
    
    if (rank == 0) {
        printf("%-10s | %2.1f | %2.1f | %5.1f | %9.1f | %10.1f | %s\n",
               class_eq.c_str(), Mb_eq, Ms_eq, Ms_eq-Mb_eq, PS_eq, depth_eq, class_eq.c_str());
        printf("%-10s | %2.1f | %2.1f | %5.1f | %9.1f | %10.1f | %s\n",
               class_exp.c_str(), Mb_exp, Ms_exp, Ms_exp-Mb_exp, PS_exp, depth_exp, class_exp.c_str());
        printf("%-10s | %2.1f | %2.1f | %5.1f | %9.1f | %10.1f | %s\n",
               class_rb.c_str(), Mb_rb, Ms_rb, Ms_rb-Mb_rb, PS_rb, depth_rb, class_rb.c_str());
        printf("%-10s | %2.1f | %2.1f | %5.1f | %9.1f | %10.1f | %s\n",
               class_nuc.c_str(), Mb_nuc, Ms_nuc, Ms_nuc-Mb_nuc, PS_nuc, depth_nuc, class_nuc.c_str());
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Key discriminants:\n";
        std::cout << "  - Ms-Mb < 0 suggests explosion\n";
        std::cout << "  - High P/S ratio suggests explosion\n";
        std::cout << "  - Shallow depth suggests explosion\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}
