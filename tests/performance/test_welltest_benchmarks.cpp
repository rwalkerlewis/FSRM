/**
 * @file test_welltest_benchmarks.cpp
 * @brief Well test analysis benchmarks
 * 
 * Benchmarks for:
 * - Pressure buildup/drawdown analysis
 * - Derivative analysis
 * - Wellbore storage effects
 * - Skin factor determination
 * - Reservoir boundaries detection
 * - Interference testing
 */

#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <vector>

class WellTestBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    int rank, size;
    
    // Exponential integral (used in well testing)
    double Ei(double x) {
        if (x < 1.0e-10) return 1.0e10;
        if (x > 10.0) return 0.0;
        
        const double EULER = 0.5772156649;
        if (x < 1.0) {
            double sum = -EULER - std::log(x);
            double term = x;
            for (int n = 1; n < 100; ++n) {
                sum += term / (n * std::tgamma(n + 1));
                term *= -x;
                if (std::abs(term) < 1.0e-10) break;
            }
            return sum;
        }
        
        double a = 1.0, b = x + 1.0, c = 1.0e30, d = 1.0 / b;
        double h = d;
        for (int n = 1; n < 100; ++n) {
            double an = -n * n;
            b += 2.0;
            d = 1.0 / (an * d + b);
            c = b + an / c;
            double delta = c * d;
            h *= delta;
            if (std::abs(delta - 1.0) < 1.0e-10) break;
        }
        return h * std::exp(-x);
    }
};

// ============================================================================
// Pressure Drawdown Analysis
// ============================================================================

TEST_F(WellTestBenchmark, PressureDrawdownAnalysis) {
    // Analyze pressure drawdown test
    
    double Q = 0.01;              // Flow rate (m³/s)
    double k = 100.0e-15;         // Permeability (m²)
    double h = 10.0;              // Thickness (m)
    double phi = 0.2;             // Porosity
    double mu = 0.001;            // Viscosity (Pa·s)
    double c_t = 1.0e-9;          // Compressibility (1/Pa)
    double r_w = 0.1;             // Wellbore radius (m)
    double S = 0.0;               // Skin factor
    
    double T = k * h / mu;        // Transmissivity
    double S_storage = phi * c_t * h;  // Storativity
    
    if (rank == 0) {
        std::cout << "\n=== Pressure Drawdown Analysis ===\n";
        std::cout << "Q = " << Q*86400 << " m³/day, k = " << k*1e15 << " mD\n\n";
        std::cout << "Time (hr) | Pressure drop (kPa) | Derivative | Flow regime\n";
        std::cout << "----------|---------------------|------------|-------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> times = {0.01, 0.1, 1, 10, 100};  // hours
    
    for (double t_hr : times) {
        double t = t_hr * 3600.0;  // Convert to seconds
        
        // Dimensionless time
        double t_D = T * t / (r_w * r_w * S_storage);
        
        // Line source solution (Theis)
        double p_D = 0.5 * Ei(1.0 / (4.0 * t_D));
        
        // Pressure drop (including skin)
        double Delta_p = Q * mu / (4.0 * M_PI * k * h) * (p_D + S);
        
        // Derivative (for flow regime identification)
        double derivative = Q * mu / (4.0 * M_PI * k * h * t);
        
        // Identify flow regime
        std::string regime;
        if (t_D < 100) {
            regime = "Radial";
        } else if (t_D < 1000) {
            regime = "Transition";
        } else {
            regime = "Boundary";
        }
        
        if (rank == 0) {
            printf("%9.2f | %19.2f | %10.2e | %11s\n",
                   t_hr, Delta_p/1000.0, derivative, regime.c_str());
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
// Pressure Buildup Analysis
// ============================================================================

TEST_F(WellTestBenchmark, PressureBuildupAnalysis) {
    // Homer plot analysis for buildup test
    
    double Q = 0.01;              // Production rate (m³/s)
    double k = 100.0e-15;         // Permeability
    double h = 10.0;              // Thickness
    double mu = 0.001;            // Viscosity
    double t_p = 86400.0;         // Production time (s) = 1 day
    
    if (rank == 0) {
        std::cout << "\n=== Pressure Buildup Analysis ===\n";
        std::cout << "Horner plot method\n\n";
        std::cout << "Δt (hr) | Horner time | Pressure (kPa) | Extrapolated P*\n";
        std::cout << "--------|-------------|----------------|------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> delta_t = {0.1, 1, 10, 100, 1000};  // hours
    
    double m = Q * mu / (4.0 * M_PI * k * h);  // Slope in Horner plot
    
    for (double dt_hr : delta_t) {
        double dt = dt_hr * 3600.0;  // Convert to seconds
        
        // Horner time ratio
        double t_ratio = (t_p + dt) / dt;
        
        // Pressure (buildup approximation)
        double p_ws = m * std::log(t_ratio);
        
        // Extrapolated initial pressure (at infinite shut-in)
        double p_star = p_ws + m * std::log(t_ratio);
        
        if (rank == 0) {
            printf("%7.1f | %11.2f | %14.2f | %16.2f\n",
                   dt_hr, t_ratio, p_ws/1000.0, p_star/1000.0);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Plot pressure vs log(Horner time) for straight line\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// Wellbore Storage and Skin Effects
// ============================================================================

TEST_F(WellTestBenchmark, WellboreStorageSkin) {
    // Calculate wellbore storage coefficient and skin effects
    
    double r_w = 0.1;             // Wellbore radius (m)
    double h = 10.0;              // Thickness (m)
    double c_f = 1.0e-9;          // Fluid compressibility (1/Pa)
    
    // Wellbore storage coefficient
    double C_w = M_PI * r_w * r_w * h * c_f;
    
    if (rank == 0) {
        std::cout << "\n=== Wellbore Storage and Skin ===\n";
        std::cout << "Wellbore storage = " << C_w << " m³/Pa\n\n";
        std::cout << "Skin | Pressure drop (kPa) | Productivity | Flow efficiency\n";
        std::cout << "-----|---------------------|--------------|----------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    double k = 100.0e-15;
    double mu = 0.001;
    double Q = 0.01;
    
    std::vector<double> skin_values = {-5, -2, 0, 2, 5, 10, 20};
    
    for (double S : skin_values) {
        // Additional pressure drop due to skin
        double Delta_p_skin = Q * mu * S / (2.0 * M_PI * k * h);
        
        // Productivity index ratio
        double J_ratio = 1.0 / (1.0 + S);
        
        // Flow efficiency
        double FE = 100.0 * J_ratio;
        
        if (rank == 0) {
            printf("%4.0f | %19.2f | %12.2f | %14.1f%%\n",
                   S, Delta_p_skin/1000.0, J_ratio, FE);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Negative skin = stimulated, positive = damaged\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}

// ============================================================================
// Reservoir Boundary Detection
// ============================================================================

TEST_F(WellTestBenchmark, ReservoirBoundaryDetection) {
    // Detect reservoir boundaries from pressure response
    
    double k = 100.0e-15;
    double h = 10.0;
    double phi = 0.2;
    double mu = 0.001;
    double c_t = 1.0e-9;
    double Q = 0.01;
    
    double T = k * h / mu;
    double S = phi * c_t * h;
    
    if (rank == 0) {
        std::cout << "\n=== Reservoir Boundary Detection ===\n\n";
        std::cout << "Boundary type | Distance (m) | Time to detect (hr) | Signature\n";
        std::cout << "--------------|--------------|---------------------|------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> distances = {100, 500, 1000, 2000};
    
    for (double L : distances) {
        // Time for boundary to affect well
        double t_boundary = 0.25 * L * L * S / T;
        
        // Convert to hours
        double t_hr = t_boundary / 3600.0;
        
        // Boundary effects on derivative
        std::string signature;
        if (L < 500) {
            signature = "Early slope change";
        } else {
            signature = "Late doubling";
        }
        
        if (rank == 0) {
            printf("%-13s | %12.0f | %19.2f | %s\n",
                   "No-flow", L, t_hr, signature.c_str());
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Derivative doubles at no-flow boundary\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}

// ============================================================================
// Interference Testing
// ============================================================================

TEST_F(WellTestBenchmark, InterferenceTesting) {
    // Multi-well interference test analysis
    
    double Q = 0.01;              // Production rate at active well
    double k = 100.0e-15;
    double h = 10.0;
    double phi = 0.2;
    double mu = 0.001;
    double c_t = 1.0e-9;
    
    double T = k * h / mu;
    double S = phi * c_t * h;
    
    if (rank == 0) {
        std::cout << "\n=== Interference Testing ===\n";
        std::cout << "Observation well response\n\n";
        std::cout << "Distance (m) | Time (hr) | Pressure response (kPa) | Arrival time\n";
        std::cout << "-------------|-----------|-------------------------|-------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> distances = {50, 100, 200, 500};
    std::vector<double> times = {1, 10, 100};  // hours
    
    for (double r : distances) {
        // Pressure response arrival time
        double t_arrival = r * r * S / (4.0 * T);
        
        for (double t_hr : times) {
            double t = t_hr * 3600.0;
            
            if (t > t_arrival) {
                // Theis solution at distance r
                double u = r * r * S / (4.0 * T * t);
                double W_u = Ei(u);
                double Delta_p = Q * mu / (4.0 * M_PI * k * h) * W_u;
                
                if (rank == 0) {
                    printf("%12.0f | %9.0f | %23.3f | %11.2f hr\n",
                           r, t_hr, Delta_p/1000.0, t_arrival/3600.0);
                }
            }
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
// Type Curve Matching
// ============================================================================

TEST_F(WellTestBenchmark, TypeCurveMatching) {
    // Generate type curves for different reservoir models
    
    if (rank == 0) {
        std::cout << "\n=== Type Curve Generation ===\n";
        std::cout << "Dimensionless pressure vs dimensionless time\n\n";
        std::cout << "t_D | p_D (Infinite) | p_D (Closed) | p_D (Constant P)\n";
        std::cout << "----|----------------|--------------|------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> t_D_values;
    for (double t_D = 0.01; t_D <= 1000.0; t_D *= 2.0) {
        t_D_values.push_back(t_D);
    }
    
    for (double t_D : t_D_values) {
        // Infinite reservoir (line source)
        double p_D_inf = 0.5 * Ei(1.0 / (4.0 * t_D));
        
        // Closed reservoir (pseudo-steady state)
        double p_D_closed = 2.0 * t_D + std::log(4.0 * t_D);
        
        // Constant pressure boundary
        double p_D_const = 0.5 * (1.0 - std::exp(-4.0 * t_D));
        
        if (rank == 0 && static_cast<int>(std::log10(t_D) * 2) % 2 == 0) {
            printf("%1.2e | %14.3f | %12.3f | %16.3f\n",
                   t_D, p_D_inf, p_D_closed, p_D_const);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Match field data to type curves for reservoir ID\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Fractured Well Analysis
// ============================================================================

TEST_F(WellTestBenchmark, FracturedWellAnalysis) {
    // Infinite conductivity vertical fracture
    
    double x_f = 50.0;            // Fracture half-length (m)
    double k = 10.0e-15;          // Formation permeability (m²)
    double k_f = 1000.0e-15;      // Fracture permeability (m²)
    double w_f = 0.005;           // Fracture width (m)
    double h = 10.0;              // Thickness (m)
    
    // Fracture conductivity
    double F_cD = k_f * w_f / (k * x_f);
    
    if (rank == 0) {
        std::cout << "\n=== Fractured Well Analysis ===\n";
        std::cout << "Vertical hydraulic fracture\n";
        std::cout << "F_cD = " << F_cD << "\n\n";
        std::cout << "Time (hr) | Flow regime | Productivity enhancement\n";
        std::cout << "----------|-------------|-------------------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<double> times = {0.01, 0.1, 1, 10, 100};
    
    for (double t_hr : times) {
        std::string regime;
        double PI_ratio;
        
        if (t_hr < 0.1) {
            regime = "Fracture linear";
            PI_ratio = 5.0;
        } else if (t_hr < 10.0) {
            regime = "Bilinear";
            PI_ratio = 3.0;
        } else {
            regime = "Pseudo-radial";
            PI_ratio = 1.5;
        }
        
        if (rank == 0) {
            printf("%9.2f | %-15s | %23.1f\n", t_hr, regime.c_str(), PI_ratio);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Note: Fractures significantly enhance productivity\n";
    }
    
    EXPECT_LT(duration.count(), 2000);
}
