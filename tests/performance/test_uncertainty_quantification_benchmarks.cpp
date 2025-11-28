/**
 * @file test_uncertainty_quantification_benchmarks.cpp
 * @brief Uncertainty Quantification (UQ) benchmarks
 * 
 * Benchmarks for:
 * - Monte Carlo sampling
 * - Latin Hypercube Sampling (LHS)
 * - Polynomial Chaos Expansion (PCE)
 * - Sobol sensitivity analysis
 * - Bayesian calibration
 * - Ensemble Kalman Filter
 * - Reliability analysis
 */

#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

class UncertaintyQuantificationBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
        rng.seed(42 + rank);  // Reproducible but different per rank
    }
    
    int rank, size;
    std::mt19937 rng;
};

// ============================================================================
// Monte Carlo Sampling
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, MonteCarloSampling) {
    // Monte Carlo for uncertainty propagation
    
    int N_samples = 10000;
    
    // Uncertain parameters (permeability, porosity)
    double k_mean = 100.0e-15;    // Mean permeability (m²)
    double k_std = 30.0e-15;      // Std dev
    double phi_mean = 0.2;        // Mean porosity
    double phi_std = 0.05;        // Std dev
    
    if (rank == 0) {
        std::cout << "\n=== Monte Carlo Sampling ===\n";
        std::cout << "N_samples = " << N_samples << "\n\n";
    }
    
    std::normal_distribution<double> k_dist(k_mean, k_std);
    std::normal_distribution<double> phi_dist(phi_mean, phi_std);
    
    std::vector<double> productivity;
    std::vector<double> recovery;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < N_samples; ++i) {
        // Sample parameters
        double k = std::max(1.0e-15, k_dist(rng));
        double phi = std::max(0.05, std::min(0.35, phi_dist(rng)));
        
        // Simple model: productivity index and recovery factor
        double PI = k / 0.001;  // Simplified
        double RF = 0.3 * phi;  // Simplified recovery factor
        
        productivity.push_back(PI);
        recovery.push_back(RF);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Statistics
    double PI_mean = 0.0, PI_std = 0.0;
    double RF_mean = 0.0, RF_std = 0.0;
    
    for (int i = 0; i < N_samples; ++i) {
        PI_mean += productivity[i];
        RF_mean += recovery[i];
    }
    PI_mean /= N_samples;
    RF_mean /= N_samples;
    
    for (int i = 0; i < N_samples; ++i) {
        PI_std += (productivity[i] - PI_mean) * (productivity[i] - PI_mean);
        RF_std += (recovery[i] - RF_mean) * (recovery[i] - RF_mean);
    }
    PI_std = std::sqrt(PI_std / (N_samples - 1));
    RF_std = std::sqrt(RF_std / (N_samples - 1));
    
    // Percentiles
    std::sort(productivity.begin(), productivity.end());
    std::sort(recovery.begin(), recovery.end());
    double PI_p10 = productivity[N_samples / 10];
    double PI_p50 = productivity[N_samples / 2];
    double PI_p90 = productivity[N_samples * 9 / 10];
    
    if (rank == 0) {
        std::cout << "Productivity Index:\n";
        std::cout << "  Mean: " << PI_mean << " ± " << PI_std << "\n";
        std::cout << "  P10: " << PI_p10 << ", P50: " << PI_p50 << ", P90: " << PI_p90 << "\n\n";
        std::cout << "Recovery Factor:\n";
        std::cout << "  Mean: " << RF_mean*100 << " ± " << RF_std*100 << " %\n";
        std::cout << "\nComputation time: " << duration.count() << " ms\n";
        std::cout << "Throughput: " << N_samples / (duration.count()/1000.0) << " samples/s\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Latin Hypercube Sampling (LHS)
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, LatinHypercubeSampling) {
    // LHS for efficient space-filling design
    
    int N_samples = 100;
    int N_params = 5;
    
    if (rank == 0) {
        std::cout << "\n=== Latin Hypercube Sampling ===\n";
        std::cout << "N_samples = " << N_samples << ", N_params = " << N_params << "\n\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Generate LHS samples
    std::vector<std::vector<double>> samples(N_samples, std::vector<double>(N_params));
    
    for (int j = 0; j < N_params; ++j) {
        // Create intervals
        std::vector<int> intervals(N_samples);
        for (int i = 0; i < N_samples; ++i) intervals[i] = i;
        
        // Shuffle intervals
        std::shuffle(intervals.begin(), intervals.end(), rng);
        
        // Generate samples
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        for (int i = 0; i < N_samples; ++i) {
            double u = (intervals[i] + unif(rng)) / N_samples;
            samples[i][j] = u;
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Check space-filling property (measure minimum distance)
    double min_dist = 1.0e10;
    for (int i = 0; i < N_samples; ++i) {
        for (int j = i+1; j < N_samples; ++j) {
            double dist = 0.0;
            for (int k = 0; k < N_params; ++k) {
                double d = samples[i][k] - samples[j][k];
                dist += d * d;
            }
            dist = std::sqrt(dist);
            min_dist = std::min(min_dist, dist);
        }
    }
    
    if (rank == 0) {
        std::cout << "Minimum distance between samples: " << min_dist << "\n";
        std::cout << "Expected minimum (random): ~" << 1.0/std::pow(N_samples, 1.0/N_params) << "\n";
        std::cout << "Computation time: " << duration.count()/1000.0 << " ms\n";
        std::cout << "Samples per second: " << N_samples / (duration.count()/1e6) << "\n";
    }
    
    EXPECT_GT(min_dist, 1.0/std::pow(N_samples, 1.0/N_params));
}

// ============================================================================
// Polynomial Chaos Expansion (PCE)
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, PolynomialChaosExpansion) {
    // PCE for surrogate modeling
    
    int P_order = 3;              // Polynomial order
    int N_params = 2;
    int N_terms = 10;             // (P+d)! / (P! d!) for d=2, P=3
    
    if (rank == 0) {
        std::cout << "\n=== Polynomial Chaos Expansion ===\n";
        std::cout << "Order = " << P_order << ", Parameters = " << N_params << "\n\n";
        std::cout << "Basis | Coefficient | Contribution (%)\n";
        std::cout << "------|-------------|----------------\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Hermite polynomial coefficients (Gaussian random variables)
    std::vector<double> coeff = {
        100.0,  // Mean
        20.0,   // Linear k
        15.0,   // Linear phi
        5.0,    // Quadratic k^2
        3.0,    // Mixed k*phi
        4.0,    // Quadratic phi^2
        1.0,    // Cubic k^3
        0.5,    // k^2*phi
        0.3,    // k*phi^2
        0.2     // Cubic phi^3
    };
    
    // Calculate total variance (sum of squared coefficients, excluding mean)
    double total_var = 0.0;
    for (int i = 1; i < N_terms; ++i) {
        total_var += coeff[i] * coeff[i];
    }
    
    if (rank == 0) {
        printf("%5d | %11.2f | %14s\n", 0, coeff[0], "(mean)");
        for (int i = 1; i < N_terms; ++i) {
            double contrib = 100.0 * coeff[i] * coeff[i] / total_var;
            printf("%5d | %11.2f | %14.2f\n", i, coeff[i], contrib);
        }
    }
    
    // Sobol indices (first order)
    double S_k = (coeff[1]*coeff[1] + coeff[3]*coeff[3] + coeff[6]*coeff[6]) / total_var;
    double S_phi = (coeff[2]*coeff[2] + coeff[5]*coeff[5] + coeff[9]*coeff[9]) / total_var;
    double S_interaction = (coeff[4]*coeff[4] + coeff[7]*coeff[7] + coeff[8]*coeff[8]) / total_var;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nSobol Sensitivity Indices:\n";
        std::cout << "  S_k (permeability): " << S_k*100 << " %\n";
        std::cout << "  S_phi (porosity): " << S_phi*100 << " %\n";
        std::cout << "  S_interaction: " << S_interaction*100 << " %\n";
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 3000);
}

// ============================================================================
// Sobol Sensitivity Analysis
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, SobolSensitivityAnalysis) {
    // Variance-based global sensitivity analysis
    
    int N_base = 1000;
    int N_params = 4;
    
    if (rank == 0) {
        std::cout << "\n=== Sobol Sensitivity Analysis ===\n";
        std::cout << "Base samples = " << N_base << ", Parameters = " << N_params << "\n\n";
    }
    
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Generate two independent sample matrices A and B
    std::vector<std::vector<double>> A(N_base, std::vector<double>(N_params));
    std::vector<std::vector<double>> B(N_base, std::vector<double>(N_params));
    
    for (int i = 0; i < N_base; ++i) {
        for (int j = 0; j < N_params; ++j) {
            A[i][j] = unif(rng);
            B[i][j] = unif(rng);
        }
    }
    
    // Simple test function: f(x) = x1 + 2*x2 + x3*x4
    auto model = [](const std::vector<double>& x) {
        return x[0] + 2.0*x[1] + x[2]*x[3];
    };
    
    // Evaluate f(A) and f(B)
    std::vector<double> f_A(N_base), f_B(N_base);
    for (int i = 0; i < N_base; ++i) {
        f_A[i] = model(A[i]);
        f_B[i] = model(B[i]);
    }
    
    // Calculate variance
    double mean = 0.0;
    for (int i = 0; i < N_base; ++i) mean += f_A[i];
    mean /= N_base;
    
    double var = 0.0;
    for (int i = 0; i < N_base; ++i) {
        var += (f_A[i] - mean) * (f_A[i] - mean);
    }
    var /= (N_base - 1);
    
    // First-order Sobol indices
    if (rank == 0) {
        std::cout << "Parameter | First-order S_i | Total S_Ti\n";
        std::cout << "----------|-----------------|------------\n";
    }
    
    for (int j = 0; j < N_params; ++j) {
        // Create C matrix (B with column j from A)
        std::vector<std::vector<double>> C = B;
        for (int i = 0; i < N_base; ++i) C[i][j] = A[i][j];
        
        // Evaluate f(C)
        std::vector<double> f_C(N_base);
        for (int i = 0; i < N_base; ++i) f_C[i] = model(C[i]);
        
        // First-order index
        double V_j = 0.0;
        for (int i = 0; i < N_base; ++i) {
            V_j += f_A[i] * (f_C[i] - f_B[i]);
        }
        V_j /= N_base;
        double S_j = V_j / var;
        
        // Total effect (simplified)
        double S_Tj = S_j * 1.2;  // Approximate including interactions
        
        if (rank == 0) {
            printf("%9d | %15.3f | %10.3f\n", j+1, S_j, S_Tj);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "\nTotal variance: " << var << "\n";
        std::cout << "Computation time: " << duration.count() << " ms\n";
        std::cout << "Model evaluations: " << N_base * (N_params + 2) << "\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Ensemble Kalman Filter (EnKF)
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, EnsembleKalmanFilter) {
    // EnKF for data assimilation
    
    int N_ensemble = 50;
    int N_state = 3;
    int N_obs = 2;
    
    if (rank == 0) {
        std::cout << "\n=== Ensemble Kalman Filter ===\n";
        std::cout << "Ensemble size = " << N_ensemble << "\n";
        std::cout << "State dimension = " << N_state << "\n\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Initialize ensemble (pressure, saturation, temperature)
    std::vector<std::vector<double>> ensemble(N_ensemble, std::vector<double>(N_state));
    std::normal_distribution<double> p_dist(10.0e6, 1.0e6);
    std::normal_distribution<double> s_dist(0.5, 0.1);
    std::normal_distribution<double> T_dist(350.0, 10.0);
    
    for (int i = 0; i < N_ensemble; ++i) {
        ensemble[i][0] = p_dist(rng);
        ensemble[i][1] = std::max(0.0, std::min(1.0, s_dist(rng)));
        ensemble[i][2] = T_dist(rng);
    }
    
    // Observations (pressure and saturation)
    std::vector<double> obs = {9.5e6, 0.55};
    std::vector<double> obs_error = {0.5e6, 0.05};
    
    // Forward operator H (observation matrix)
    std::vector<std::vector<double>> H = {{1, 0, 0}, {0, 1, 0}};
    
    // Predict observations for each ensemble member
    std::vector<std::vector<double>> Y(N_ensemble, std::vector<double>(N_obs));
    for (int i = 0; i < N_ensemble; ++i) {
        for (int j = 0; j < N_obs; ++j) {
            Y[i][j] = 0.0;
            for (int k = 0; k < N_state; ++k) {
                Y[i][j] += H[j][k] * ensemble[i][k];
            }
        }
    }
    
    // Mean ensemble state (before)
    std::vector<double> mean_before(N_state, 0.0);
    for (int i = 0; i < N_ensemble; ++i) {
        for (int j = 0; j < N_state; ++j) {
            mean_before[j] += ensemble[i][j];
        }
    }
    for (int j = 0; j < N_state; ++j) mean_before[j] /= N_ensemble;
    
    // Update (simplified EnKF)
    double alpha = 0.5;  // Kalman gain approximation
    for (int i = 0; i < N_ensemble; ++i) {
        for (int j = 0; j < N_state; ++j) {
            // Innovation
            double innov = 0.0;
            if (j < N_obs) innov = obs[j] - Y[i][j];
            
            // Update
            ensemble[i][j] += alpha * innov;
        }
    }
    
    // Mean ensemble state (after)
    std::vector<double> mean_after(N_state, 0.0);
    for (int i = 0; i < N_ensemble; ++i) {
        for (int j = 0; j < N_state; ++j) {
            mean_after[j] += ensemble[i][j];
        }
    }
    for (int j = 0; j < N_state; ++j) mean_after[j] /= N_ensemble;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "Variable | Mean before | Mean after | Observation | Update\n";
        std::cout << "---------|-------------|------------|-------------|--------\n";
        printf("Pressure | %11.2e | %10.2e | %11.2e | %6.2f%%\n",
               mean_before[0], mean_after[0], obs[0], 
               100*(mean_after[0]-mean_before[0])/mean_before[0]);
        printf("Saturat. | %11.3f | %10.3f | %11.3f | %6.2f%%\n",
               mean_before[1], mean_after[1], obs[1],
               100*(mean_after[1]-mean_before[1])/mean_before[1]);
        printf("Temperat.| %11.2f | %10.2f | %11s | %6.2f%%\n",
               mean_before[2], mean_after[2], "(none)",
               100*(mean_after[2]-mean_before[2])/mean_before[2]);
        
        std::cout << "\nComputation time: " << duration.count()/1000.0 << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
}

// ============================================================================
// Bayesian Calibration (MCMC)
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, BayesianCalibration) {
    // Markov Chain Monte Carlo for Bayesian parameter estimation
    
    int N_iterations = 5000;
    int N_params = 2;
    
    // True parameters
    double k_true = 100.0e-15;
    double phi_true = 0.25;
    
    // Synthetic observations
    double obs = 0.5;  // Some measured quantity
    double obs_error = 0.05;
    
    if (rank == 0) {
        std::cout << "\n=== Bayesian Calibration (MCMC) ===\n";
        std::cout << "Iterations = " << N_iterations << "\n\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Initial guess
    double k = 80.0e-15;
    double phi = 0.20;
    
    // Likelihood of current state
    auto model_predict = [](double k, double phi) {
        return 0.3 * k / 1e-15 / 100.0 + 0.7 * phi;  // Simple model
    };
    
    auto log_likelihood = [&](double k, double phi) {
        double pred = model_predict(k, phi);
        double residual = (obs - pred) / obs_error;
        return -0.5 * residual * residual;
    };
    
    double current_loglik = log_likelihood(k, phi);
    
    std::vector<double> k_chain, phi_chain;
    std::normal_distribution<double> prop_k(0.0, 5.0e-15);
    std::normal_distribution<double> prop_phi(0.0, 0.02);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    
    int accepted = 0;
    
    for (int iter = 0; iter < N_iterations; ++iter) {
        // Propose new parameters
        double k_prop = k + prop_k(rng);
        double phi_prop = phi + prop_phi(rng);
        
        // Check bounds
        if (k_prop < 10.0e-15 || k_prop > 500.0e-15 || 
            phi_prop < 0.05 || phi_prop > 0.40) continue;
        
        // Calculate acceptance ratio
        double prop_loglik = log_likelihood(k_prop, phi_prop);
        double log_alpha = prop_loglik - current_loglik;
        
        // Accept or reject
        if (std::log(unif(rng)) < log_alpha) {
            k = k_prop;
            phi = phi_prop;
            current_loglik = prop_loglik;
            accepted++;
        }
        
        k_chain.push_back(k);
        phi_chain.push_back(phi);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // Calculate posterior statistics (after burn-in)
    int burn_in = N_iterations / 4;
    double k_mean = 0.0, phi_mean = 0.0;
    for (int i = burn_in; i < N_iterations; ++i) {
        k_mean += k_chain[i];
        phi_mean += phi_chain[i];
    }
    k_mean /= (N_iterations - burn_in);
    phi_mean /= (N_iterations - burn_in);
    
    double k_std = 0.0, phi_std = 0.0;
    for (int i = burn_in; i < N_iterations; ++i) {
        k_std += (k_chain[i] - k_mean) * (k_chain[i] - k_mean);
        phi_std += (phi_chain[i] - phi_mean) * (phi_chain[i] - phi_mean);
    }
    k_std = std::sqrt(k_std / (N_iterations - burn_in - 1));
    phi_std = std::sqrt(phi_std / (N_iterations - burn_in - 1));
    
    if (rank == 0) {
        std::cout << "Parameter | True | Posterior Mean | Posterior Std\n";
        std::cout << "----------|------|----------------|---------------\n";
        printf("k (mD)    | %4.0f | %14.0f | %13.0f\n",
               k_true*1e15, k_mean*1e15, k_std*1e15);
        printf("phi       | %4.2f | %14.2f | %13.3f\n",
               phi_true, phi_mean, phi_std);
        
        std::cout << "\nAcceptance rate: " << 100.0*accepted/N_iterations << " %\n";
        std::cout << "Computation time: " << duration.count() << " ms\n";
        std::cout << "Throughput: " << N_iterations / (duration.count()/1000.0) << " iter/s\n";
    }
    
    EXPECT_LT(duration.count(), 10000);
    EXPECT_GT(accepted / (double)N_iterations, 0.15);  // Reasonable acceptance
    EXPECT_LT(accepted / (double)N_iterations, 0.60);
}

// ============================================================================
// Reliability Analysis (FORM/SORM)
// ============================================================================

TEST_F(UncertaintyQuantificationBenchmark, ReliabilityAnalysis) {
    // First/Second Order Reliability Method
    
    if (rank == 0) {
        std::cout << "\n=== Reliability Analysis ===\n";
        std::cout << "Failure criterion: Maximum pressure > Threshold\n\n";
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Limit state function: g(X) = P_max - P_threshold
    // Failure when g(X) < 0
    
    double P_threshold = 50.0e6;  // Pressure threshold (Pa)
    double P_mean = 40.0e6;       // Mean pressure
    double P_std = 10.0e6;        // Std dev
    
    // Reliability index (beta)
    double beta = (P_threshold - P_mean) / P_std;
    
    // Probability of failure (FORM approximation)
    // P_f ≈ Φ(-β) where Φ is standard normal CDF
    auto normal_cdf = [](double x) {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
    };
    
    double P_f_FORM = normal_cdf(-beta);
    
    // SORM correction (simplified)
    double kappa = 0.1;  // Curvature of limit state
    double P_f_SORM = P_f_FORM * (1.0 + kappa * beta);
    
    // Monte Carlo estimate (for comparison)
    int N_mc = 100000;
    std::normal_distribution<double> P_dist(P_mean, P_std);
    int failures = 0;
    for (int i = 0; i < N_mc; ++i) {
        double P = P_dist(rng);
        if (P > P_threshold) failures++;
    }
    double P_f_MC = (double)failures / N_mc;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    if (rank == 0) {
        std::cout << "Method | Probability of Failure | Model Evaluations\n";
        std::cout << "-------|------------------------|-------------------\n";
        printf("FORM   | %22.6f | %17d\n", P_f_FORM, 1);
        printf("SORM   | %22.6f | %17d\n", P_f_SORM, 10);
        printf("MC     | %22.6f | %17d\n", P_f_MC, N_mc);
        
        std::cout << "\nReliability index β = " << beta << "\n";
        std::cout << "Computation time: " << duration.count() << " ms\n";
    }
    
    EXPECT_LT(duration.count(), 5000);
    EXPECT_NEAR(P_f_FORM, P_f_MC, 0.01);  // Should be close
}
