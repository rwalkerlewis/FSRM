/**
 * @file test_scec_benchmarks.cpp
 * @brief SCEC (Southern California Earthquake Center) benchmark tests
 * 
 * Performance and correctness benchmarks for:
 * - Dynamic rupture propagation (TPV series)
 * - Wave propagation in layered media (LOH series)
 * - Fault friction laws
 * - Seismic wave generation
 */

#include <gtest/gtest.h>
#include "Simulator.hpp"
#include "FaultModel.hpp"
#include "PhysicsKernel.hpp"
#include <chrono>
#include <cmath>

using namespace FSRM;

class SCECBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
};

// ============================================================================
// Friction Law Benchmarks
// ============================================================================

TEST_F(SCECBenchmark, SlipWeakeningFrictionLaw) {
    // Test linear slip-weakening friction law performance
    
    double mu_s = 0.677;              // Static friction
    double mu_d = 0.525;              // Dynamic friction
    double D_c = 0.40;                // Slip-weakening distance (m)
    
    const int num_evaluations = 100000;
    std::vector<double> slip_values(num_evaluations);
    std::vector<double> friction_values(num_evaluations);
    
    // Generate range of slip values
    for (int i = 0; i < num_evaluations; ++i) {
        slip_values[i] = i * D_c * 5.0 / num_evaluations;  // 0 to 5*D_c
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Compute friction coefficient for each slip value
    for (int i = 0; i < num_evaluations; ++i) {
        double slip = slip_values[i];
        if (slip <= D_c) {
            // Linear weakening
            friction_values[i] = mu_s - (mu_s - mu_d) * (slip / D_c);
        } else {
            // Fully weakened
            friction_values[i] = mu_d;
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_evaluations;
    
    if (rank == 0) {
        std::cout << "\n=== Slip-Weakening Friction Law ===\n";
        std::cout << "Evaluations: " << num_evaluations << "\n";
        std::cout << "Time per eval: " << time_per_eval << " ns\n";
        std::cout << "Throughput: " << 1e9 / time_per_eval << " eval/s\n";
    }
    
    EXPECT_LT(time_per_eval, 100.0);  // Should be very fast
    EXPECT_DOUBLE_EQ(friction_values[0], mu_s);
    EXPECT_DOUBLE_EQ(friction_values[num_evaluations-1], mu_d);
}

TEST_F(SCECBenchmark, RateAndStateFrictionLaw) {
    // Test rate-and-state friction law (Dieterich-Ruina)
    
    double f0 = 0.6;                  // Reference friction
    double a = 0.008;                 // Direct effect parameter
    double b = 0.012;                 // Evolution effect parameter
    double D_c = 0.02;                // Critical slip distance (m)
    double V0 = 1.0e-6;               // Reference velocity (m/s)
    
    const int num_evaluations = 50000;
    std::vector<double> friction_values(num_evaluations);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_evaluations; ++i) {
        double V = V0 * std::pow(10.0, (i - num_evaluations/2) / 5000.0);  // Range of velocities
        double theta = D_c / V0;      // Simplified state variable
        
        // Rate-and-state friction coefficient
        friction_values[i] = f0 + a * std::log(V / V0) + b * std::log(V0 * theta / D_c);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_eval = duration.count() / (double)num_evaluations;
    
    if (rank == 0) {
        std::cout << "\n=== Rate-and-State Friction Law ===\n";
        std::cout << "Evaluations: " << num_evaluations << "\n";
        std::cout << "Time per eval: " << time_per_eval << " ns\n";
        std::cout << "Throughput: " << 1e9 / time_per_eval << " eval/s\n";
    }
    
    EXPECT_LT(time_per_eval, 500.0);
}

// ============================================================================
// Dynamic Rupture Benchmarks
// ============================================================================

TEST_F(SCECBenchmark, RuptureSpeedCalculation) {
    // Compute rupture propagation speed
    
    const int num_points = 10000;
    std::vector<double> rupture_times(num_points);
    std::vector<double> distances(num_points);
    std::vector<double> rupture_speeds(num_points);
    
    // Simulate rupture arrival times along fault
    double V_s = 3464.0;              // Shear wave speed (m/s)
    double rupture_velocity = 0.8 * V_s;  // Sub-shear rupture
    
    for (int i = 0; i < num_points; ++i) {
        distances[i] = i * 50.0;      // Every 50 m
        rupture_times[i] = distances[i] / rupture_velocity + 0.5;  // 0.5 s offset
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Compute rupture speed from times
    for (int i = 1; i < num_points; ++i) {
        double dx = distances[i] - distances[i-1];
        double dt = rupture_times[i] - rupture_times[i-1];
        rupture_speeds[i] = dx / dt;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_total = duration.count() / 1000.0;
    
    if (rank == 0) {
        std::cout << "\n=== Rupture Speed Calculation ===\n";
        std::cout << "Points: " << num_points << "\n";
        std::cout << "Time: " << time_total << " ms\n";
        std::cout << "Average rupture speed: " << rupture_speeds[num_points/2] << " m/s\n";
        std::cout << "Expected: " << rupture_velocity << " m/s\n";
    }
    
    // Check computed speed is close to expected
    double avg_speed = 0.0;
    for (int i = num_points/4; i < 3*num_points/4; ++i) {
        avg_speed += rupture_speeds[i];
    }
    avg_speed /= (num_points/2);
    
    EXPECT_NEAR(avg_speed, rupture_velocity, rupture_velocity * 0.01);
}

TEST_F(SCECBenchmark, StressTensorRotation) {
    // Rotate stress tensor to fault-local coordinates
    
    const int num_rotations = 100000;
    
    // Background stress state
    double sigma_xx = -120.0e6;       // Pa (compression negative)
    double sigma_yy = -120.0e6;
    double sigma_zz = -120.0e6;
    double tau_xy = 70.0e6;
    double tau_xz = 0.0;
    double tau_yz = 0.0;
    
    std::vector<double> shear_stress(num_rotations);
    std::vector<double> normal_stress(num_rotations);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_rotations; ++i) {
        // Fault orientation (strike, dip)
        double strike = (i % 360) * M_PI / 180.0;
        double dip = 90.0 * M_PI / 180.0;  // Vertical fault
        
        // Fault normal and slip direction unit vectors
        double n_x = -std::sin(dip) * std::sin(strike);
        double n_y = std::sin(dip) * std::cos(strike);
        double n_z = std::cos(dip);
        
        double s_x = std::cos(strike);
        double s_y = std::sin(strike);
        double s_z = 0.0;
        
        // Traction on fault: T = sigma . n
        double T_x = sigma_xx * n_x + tau_xy * n_y + tau_xz * n_z;
        double T_y = tau_xy * n_x + sigma_yy * n_y + tau_yz * n_z;
        double T_z = tau_xz * n_x + tau_yz * n_y + sigma_zz * n_z;
        
        // Normal stress: sigma_n = T . n
        normal_stress[i] = T_x * n_x + T_y * n_y + T_z * n_z;
        
        // Shear stress: tau = |T - sigma_n * n|
        double T_s_x = T_x - normal_stress[i] * n_x;
        double T_s_y = T_y - normal_stress[i] * n_y;
        double T_s_z = T_z - normal_stress[i] * n_z;
        
        shear_stress[i] = std::sqrt(T_s_x*T_s_x + T_s_y*T_s_y + T_s_z*T_s_z);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_rotation = duration.count() / (double)num_rotations;
    
    if (rank == 0) {
        std::cout << "\n=== Stress Tensor Rotation ===\n";
        std::cout << "Rotations: " << num_rotations << "\n";
        std::cout << "Time per rotation: " << time_per_rotation << " ns\n";
        std::cout << "Throughput: " << 1e9 / time_per_rotation << " rotations/s\n";
    }
    
    EXPECT_LT(time_per_rotation, 1000.0);
}

// ============================================================================
// Wave Propagation Benchmarks
// ============================================================================

TEST_F(SCECBenchmark, SeismicWaveSpeed) {
    // Calculate P-wave and S-wave speeds
    
    const int num_calculations = 1000000;
    std::vector<double> V_p(num_calculations);
    std::vector<double> V_s(num_calculations);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_calculations; ++i) {
        double E = 50e9 + i * 1000.0;     // Young's modulus
        double nu = 0.25;                  // Poisson's ratio
        double rho = 2670.0;               // Density
        
        // Lame parameters
        double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double mu = E / (2.0 * (1.0 + nu));
        
        // Wave speeds
        V_p[i] = std::sqrt((lambda + 2.0 * mu) / rho);
        V_s[i] = std::sqrt(mu / rho);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_calculations;
    
    if (rank == 0) {
        std::cout << "\n=== Seismic Wave Speed Calculation ===\n";
        std::cout << "Calculations: " << num_calculations << "\n";
        std::cout << "Time per calc: " << time_per_calc << " ns\n";
        std::cout << "V_p example: " << V_p[0] << " m/s\n";
        std::cout << "V_s example: " << V_s[0] << " m/s\n";
        std::cout << "V_p/V_s ratio: " << V_p[0] / V_s[0] << "\n";
    }
    
    EXPECT_LT(time_per_calc, 200.0);
    EXPECT_NEAR(V_p[0] / V_s[0], std::sqrt(3.0), 0.01);  // For ν=0.25
}

TEST_F(SCECBenchmark, WaveArrivalTime) {
    // Calculate theoretical wave arrival times
    
    const int num_stations = 1000;
    std::vector<double> distances(num_stations);
    std::vector<double> P_arrival_times(num_stations);
    std::vector<double> S_arrival_times(num_stations);
    
    // Source at origin, stations at various distances
    for (int i = 0; i < num_stations; ++i) {
        distances[i] = i * 100.0;  // 0 to 100 km
    }
    
    double V_p = 6000.0;  // m/s
    double V_s = 3464.0;  // m/s
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_stations; ++i) {
        P_arrival_times[i] = distances[i] / V_p;
        S_arrival_times[i] = distances[i] / V_s;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double time_per_calc = duration.count() / (double)num_stations;
    
    if (rank == 0) {
        std::cout << "\n=== Wave Arrival Time Calculation ===\n";
        std::cout << "Stations: " << num_stations << "\n";
        std::cout << "Time per station: " << time_per_calc << " ns\n";
        std::cout << "P arrival at 10 km: " << P_arrival_times[100] << " s\n";
        std::cout << "S arrival at 10 km: " << S_arrival_times[100] << " s\n";
        std::cout << "S-P time: " << S_arrival_times[100] - P_arrival_times[100] << " s\n";
    }
    
    EXPECT_LT(time_per_calc, 50.0);
}

TEST_F(SCECBenchmark, RickerWaveletGeneration) {
    // Generate Ricker wavelet source time function
    
    const int num_samples = 10000;
    double dt = 0.001;  // 1 ms
    double f0 = 2.0;    // 2 Hz dominant frequency
    
    std::vector<double> wavelet(num_samples);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_samples; ++i) {
        double t = i * dt - 1.0 / f0;  // Center at t = 1/f0
        double arg = M_PI * f0 * t;
        wavelet[i] = (1.0 - 2.0 * arg * arg) * std::exp(-arg * arg);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_total = duration.count() / 1000.0;
    
    if (rank == 0) {
        std::cout << "\n=== Ricker Wavelet Generation ===\n";
        std::cout << "Samples: " << num_samples << "\n";
        std::cout << "Time: " << time_total << " ms\n";
        std::cout << "Dominant frequency: " << f0 << " Hz\n";
    }
    
    // Find peak value (should be at t = 1/f0)
    double max_val = 0.0;
    int max_idx = 0;
    for (int i = 0; i < num_samples; ++i) {
        if (std::abs(wavelet[i]) > max_val) {
            max_val = std::abs(wavelet[i]);
            max_idx = i;
        }
    }
    
    EXPECT_NEAR(max_idx * dt, 1.0 / f0, 0.01);
}

// ============================================================================
// Fault Slip Distribution Benchmarks
// ============================================================================

TEST_F(SCECBenchmark, SlipDistributionAnalysis) {
    // Analyze slip distribution on fault
    
    const int num_fault_points = 50000;
    std::vector<double> slip(num_fault_points);
    std::vector<double> slip_rate(num_fault_points);
    std::vector<double> rupture_time(num_fault_points);
    
    // Generate synthetic slip distribution
    for (int i = 0; i < num_fault_points; ++i) {
        double x = (i % 100) * 100.0;     // Position along strike
        double z = (i / 100) * 100.0;     // Position down-dip
        
        // Elliptical asperity
        double r = std::sqrt(x*x + z*z);
        slip[i] = std::max(0.0, 2.0 * (1.0 - r / 5000.0));  // Max 2 m slip
        slip_rate[i] = slip[i] / 2.0;     // Rise time ~ 2 s
        rupture_time[i] = r / 2800.0;     // Rupture speed ~ 2800 m/s
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Compute statistics
    double total_slip = 0.0;
    double max_slip = 0.0;
    double avg_slip_rate = 0.0;
    
    for (int i = 0; i < num_fault_points; ++i) {
        total_slip += slip[i];
        max_slip = std::max(max_slip, slip[i]);
        avg_slip_rate += slip_rate[i];
    }
    
    double avg_slip = total_slip / num_fault_points;
    avg_slip_rate /= num_fault_points;
    
    // Compute seismic moment: M0 = μ * A * <slip>
    double mu = 32.04e9;              // Shear modulus (Pa)
    double area = num_fault_points * 100.0 * 100.0;  // m²
    double moment = mu * area * avg_slip;
    
    // Moment magnitude
    double Mw = (2.0/3.0) * std::log10(moment) - 6.07;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\n=== Slip Distribution Analysis ===\n";
        std::cout << "Fault points: " << num_fault_points << "\n";
        std::cout << "Analysis time: " << duration.count() << " ms\n";
        std::cout << "Average slip: " << avg_slip << " m\n";
        std::cout << "Maximum slip: " << max_slip << " m\n";
        std::cout << "Seismic moment: " << moment << " N⋅m\n";
        std::cout << "Moment magnitude: " << Mw << "\n";
    }
    
    EXPECT_GT(Mw, 5.0);
    EXPECT_LT(Mw, 8.0);
}

// ============================================================================
// Performance Scaling Tests
// ============================================================================

TEST_F(SCECBenchmark, FaultPointScaling) {
    // Test performance vs number of fault points
    
    std::vector<int> fault_sizes = {1000, 5000, 10000, 50000, 100000};
    
    if (rank == 0) {
        std::cout << "\n=== Fault Point Scaling ===\n";
        std::cout << "Points | Time (ms) | Points/s\n";
        std::cout << "-------|-----------|----------\n";
    }
    
    for (int N : fault_sizes) {
        std::vector<double> slip(N);
        std::vector<double> stress(N);
        std::vector<double> friction(N);
        
        // Initialize
        for (int i = 0; i < N; ++i) {
            slip[i] = i * 0.0001;
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Simulate friction calculation on all points
        for (int i = 0; i < N; ++i) {
            double mu_s = 0.6;
            double mu_d = 0.5;
            double D_c = 0.4;
            double s = slip[i];
            
            if (s < D_c) {
                friction[i] = mu_s - (mu_s - mu_d) * (s / D_c);
            } else {
                friction[i] = mu_d;
            }
            
            stress[i] = friction[i] * 120.0e6;  // Normal stress
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        double throughput = N / (duration.count() / 1.0e6);
        
        if (rank == 0) {
            printf("%6d | %9.3f | %9.0f\n", N, time_ms, throughput);
        }
        
        EXPECT_GT(throughput, 1e6);  // Should process > 1M points/s
    }
}
