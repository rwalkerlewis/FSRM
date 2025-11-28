/**
 * @file test_gpu_benchmarks.cpp
 * @brief GPU performance benchmarks
 * 
 * Benchmarks comparing CPU vs GPU performance for:
 * - Memory bandwidth
 * - Kernel execution
 * - Data transfers
 * - Different problem sizes
 */

#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <petsc.h>

#ifdef ENABLE_CUDA
#include "GPUManager.hpp"
#include "GPUKernels.cuh"
#include "PhysicsKernel_GPU.hpp"
using namespace FSRM;
#endif

class GPUBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
        
#ifdef ENABLE_CUDA
        // Initialize GPU if available
        gpu_available = GPUManager::getInstance().isAvailable();
        if (gpu_available && rank == 0) {
            std::cout << "\nGPU detected: " << GPUManager::getInstance().getDeviceName() << "\n";
            std::cout << "GPU Memory: " << GPUManager::getInstance().getTotalMemoryGB() << " GB\n\n";
        }
#else
        gpu_available = false;
        if (rank == 0) {
            std::cout << "\nGPU support not enabled in this build\n";
            std::cout << "Rebuild with -DENABLE_CUDA=ON to run GPU benchmarks\n\n";
        }
#endif
    }
    
    int rank, size;
    bool gpu_available = false;
};

#ifdef ENABLE_CUDA

// ============================================================================
// Memory Bandwidth Benchmarks
// ============================================================================

TEST_F(GPUBenchmark, MemoryBandwidthCopyToDevice) {
    if (!gpu_available) {
        GTEST_SKIP() << "GPU not available";
    }
    
    std::vector<size_t> sizes = {
        1024 * 1024,           // 1 MB
        10 * 1024 * 1024,      // 10 MB
        100 * 1024 * 1024,     // 100 MB
        1024 * 1024 * 1024     // 1 GB
    };
    
    if (rank == 0) {
        std::cout << "=== GPU Memory Bandwidth (Host to Device) ===\n";
        std::cout << "Size (MB) | Time (ms) | Bandwidth (GB/s)\n";
        std::cout << "----------|-----------|------------------\n";
    }
    
    for (size_t size : sizes) {
        std::vector<double> host_data(size / sizeof(double));
        for (size_t i = 0; i < host_data.size(); ++i) {
            host_data[i] = i * 0.001;
        }
        
        double* device_data = nullptr;
        cudaMalloc(&device_data, size);
        
        // Warm-up
        cudaMemcpy(device_data, host_data.data(), size, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        
        // Benchmark
        auto start = std::chrono::high_resolution_clock::now();
        cudaMemcpy(device_data, host_data.data(), size, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        double bandwidth_gbps = (size / 1e9) / (duration.count() / 1e6);
        
        if (rank == 0) {
            printf("%9.1f | %9.3f | %16.2f\n", 
                   size / (1024.0 * 1024.0), time_ms, bandwidth_gbps);
        }
        
        cudaFree(device_data);
        
        EXPECT_GT(bandwidth_gbps, 1.0);  // Should be at least 1 GB/s
    }
}

TEST_F(GPUBenchmark, MemoryBandwidthCopyToHost) {
    if (!gpu_available) {
        GTEST_SKIP() << "GPU not available";
    }
    
    std::vector<size_t> sizes = {
        1024 * 1024,           // 1 MB
        10 * 1024 * 1024,      // 10 MB
        100 * 1024 * 1024,     // 100 MB
        1024 * 1024 * 1024     // 1 GB
    };
    
    if (rank == 0) {
        std::cout << "\n=== GPU Memory Bandwidth (Device to Host) ===\n";
        std::cout << "Size (MB) | Time (ms) | Bandwidth (GB/s)\n";
        std::cout << "----------|-----------|------------------\n";
    }
    
    for (size_t size : sizes) {
        std::vector<double> host_data(size / sizeof(double));
        double* device_data = nullptr;
        cudaMalloc(&device_data, size);
        
        // Warm-up
        cudaMemcpy(host_data.data(), device_data, size, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        
        // Benchmark
        auto start = std::chrono::high_resolution_clock::now();
        cudaMemcpy(host_data.data(), device_data, size, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        double bandwidth_gbps = (size / 1e9) / (duration.count() / 1e6);
        
        if (rank == 0) {
            printf("%9.1f | %9.3f | %16.2f\n", 
                   size / (1024.0 * 1024.0), time_ms, bandwidth_gbps);
        }
        
        cudaFree(device_data);
        
        EXPECT_GT(bandwidth_gbps, 1.0);
    }
}

// ============================================================================
// Kernel Performance Benchmarks
// ============================================================================

TEST_F(GPUBenchmark, VectorAdditionPerformance) {
    if (!gpu_available) {
        GTEST_SKIP() << "GPU not available";
    }
    
    std::vector<int> sizes = {1000, 10000, 100000, 1000000, 10000000};
    
    if (rank == 0) {
        std::cout << "\n=== Vector Addition Performance ===\n";
        std::cout << "Size | CPU (ms) | GPU (ms) | Speedup\n";
        std::cout << "-----|----------|----------|--------\n";
    }
    
    for (int N : sizes) {
        std::vector<double> a(N), b(N), c_cpu(N), c_gpu(N);
        
        // Initialize
        for (int i = 0; i < N; ++i) {
            a[i] = i * 0.001;
            b[i] = i * 0.002;
        }
        
        // CPU version
        auto start_cpu = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            c_cpu[i] = a[i] + b[i];
        }
        auto end_cpu = std::chrono::high_resolution_clock::now();
        double time_cpu_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_cpu - start_cpu).count() / 1000.0;
        
        // GPU version
        double *d_a, *d_b, *d_c;
        cudaMalloc(&d_a, N * sizeof(double));
        cudaMalloc(&d_b, N * sizeof(double));
        cudaMalloc(&d_c, N * sizeof(double));
        
        cudaMemcpy(d_a, a.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_b, b.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        
        auto start_gpu = std::chrono::high_resolution_clock::now();
        
        int blockSize = 256;
        int numBlocks = (N + blockSize - 1) / blockSize;
        vectorAdd<<<numBlocks, blockSize>>>(d_a, d_b, d_c, N);
        cudaDeviceSynchronize();
        
        auto end_gpu = std::chrono::high_resolution_clock::now();
        double time_gpu_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_gpu - start_gpu).count() / 1000.0;
        
        cudaMemcpy(c_gpu.data(), d_c, N * sizeof(double), cudaMemcpyDeviceToHost);
        
        double speedup = time_cpu_ms / time_gpu_ms;
        
        if (rank == 0) {
            printf("%5d | %8.3f | %8.3f | %6.2fx\n", N, time_cpu_ms, time_gpu_ms, speedup);
        }
        
        cudaFree(d_a);
        cudaFree(d_b);
        cudaFree(d_c);
        
        // Verify correctness
        double max_error = 0.0;
        for (int i = 0; i < N; ++i) {
            max_error = std::max(max_error, std::abs(c_cpu[i] - c_gpu[i]));
        }
        EXPECT_LT(max_error, 1e-10);
    }
}

TEST_F(GPUBenchmark, SinglePhaseFlowKernelGPU) {
    if (!gpu_available) {
        GTEST_SKIP() << "GPU not available";
    }
    
    std::vector<int> grid_sizes = {10, 20, 40, 80};
    
    if (rank == 0) {
        std::cout << "\n=== Single Phase Flow Kernel (GPU) ===\n";
        std::cout << "Grid | Cells | CPU (ms) | GPU (ms) | Speedup\n";
        std::cout << "-----|-------|----------|----------|--------\n";
    }
    
    for (int n : grid_sizes) {
        int num_cells = n * n * n;
        
        // Create test data
        std::vector<double> pressure(num_cells);
        std::vector<double> pressure_old(num_cells);
        std::vector<double> residual_cpu(num_cells);
        std::vector<double> residual_gpu(num_cells);
        
        for (int i = 0; i < num_cells; ++i) {
            pressure[i] = 1e7 + i * 1000.0;
            pressure_old[i] = 1e7 + (i - 1) * 1000.0;
        }
        
        // Material properties
        double phi = 0.2;
        double k = 100e-15;
        double mu = 0.001;
        double rho = 1000.0;
        double c_f = 1e-9;
        double dt = 1.0;
        
        // CPU version
        auto start_cpu = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_cells; ++i) {
            double accumulation = phi * c_f * (pressure[i] - pressure_old[i]) / dt;
            double flux = k / mu * 1e5;  // Simplified flux
            residual_cpu[i] = accumulation - flux;
        }
        auto end_cpu = std::chrono::high_resolution_clock::now();
        double time_cpu_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_cpu - start_cpu).count() / 1000.0;
        
        // GPU version
        PhysicsKernelGPU kernel_gpu;
        kernel_gpu.setSinglePhaseProperties(phi, k, c_f, mu, rho);
        
        auto start_gpu = std::chrono::high_resolution_clock::now();
        kernel_gpu.computeSinglePhaseResidual(
            pressure.data(), pressure_old.data(), residual_gpu.data(), 
            num_cells, dt);
        auto end_gpu = std::chrono::high_resolution_clock::now();
        double time_gpu_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_gpu - start_gpu).count() / 1000.0;
        
        double speedup = time_cpu_ms / time_gpu_ms;
        
        if (rank == 0) {
            printf("%3dx%3dx%3d | %6d | %8.3f | %8.3f | %6.2fx\n", 
                   n, n, n, num_cells, time_cpu_ms, time_gpu_ms, speedup);
        }
        
        EXPECT_GT(speedup, 0.5);  // GPU should be competitive
    }
}

TEST_F(GPUBenchmark, PoroelasticKernelGPU) {
    if (!gpu_available) {
        GTEST_SKIP() << "GPU not available";
    }
    
    std::vector<int> grid_sizes = {10, 20, 40};
    
    if (rank == 0) {
        std::cout << "\n=== Poroelastic Kernel (GPU) ===\n";
        std::cout << "Grid | Cells | DOFs | CPU (ms) | GPU (ms) | Speedup\n";
        std::cout << "-----|-------|------|----------|----------|--------\n";
    }
    
    for (int n : grid_sizes) {
        int num_cells = n * n * n;
        int dofs = num_cells * 4;  // p, ux, uy, uz
        
        std::vector<double> state(dofs);
        std::vector<double> residual_cpu(dofs);
        std::vector<double> residual_gpu(dofs);
        
        // Initialize
        for (int i = 0; i < dofs; ++i) {
            state[i] = (i % 4 == 0) ? 1e7 : 0.001;
        }
        
        // Material properties
        double phi = 0.2;
        double k = 1e-13;
        double K_s = 36e9;
        double K_f = 2.2e9;
        double mu_s = 14e9;
        double mu_f = 0.001;
        double rho_s = 2650.0;
        double rho_f = 1000.0;
        
        // CPU version (simplified)
        auto start_cpu = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_cells; ++i) {
            int idx = i * 4;
            // Simplified poroelastic residual
            residual_cpu[idx] = phi * (state[idx] - 1e7) * 1e-9;  // Pressure
            residual_cpu[idx+1] = mu_s * state[idx+1];            // ux
            residual_cpu[idx+2] = mu_s * state[idx+2];            // uy
            residual_cpu[idx+3] = mu_s * state[idx+3];            // uz
        }
        auto end_cpu = std::chrono::high_resolution_clock::now();
        double time_cpu_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_cpu - start_cpu).count() / 1000.0;
        
        // GPU version
        PhysicsKernelGPU kernel_gpu;
        kernel_gpu.setPoroelasticProperties(phi, k, K_s, K_f, mu_s, mu_f, rho_s, rho_f);
        
        auto start_gpu = std::chrono::high_resolution_clock::now();
        kernel_gpu.computePoroelasticResidual(
            state.data(), residual_gpu.data(), num_cells);
        auto end_gpu = std::chrono::high_resolution_clock::now();
        double time_gpu_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_gpu - start_gpu).count() / 1000.0;
        
        double speedup = time_cpu_ms / time_gpu_ms;
        
        if (rank == 0) {
            printf("%3dx%3dx%3d | %6d | %5d | %8.3f | %8.3f | %6.2fx\n", 
                   n, n, n, num_cells, dofs, time_cpu_ms, time_gpu_ms, speedup);
        }
        
        EXPECT_GT(speedup, 0.5);
    }
}

// ============================================================================
// GPU Scaling Benchmarks
// ============================================================================

TEST_F(GPUBenchmark, GPUStrongScaling) {
    if (!gpu_available) {
        GTEST_SKIP() << "GPU not available";
    }
    
    int num_cells = 1000000;  // 1M cells
    int num_iterations = 100;
    
    std::vector<int> block_sizes = {32, 64, 128, 256, 512, 1024};
    
    if (rank == 0) {
        std::cout << "\n=== GPU Strong Scaling (1M cells) ===\n";
        std::cout << "Block Size | Time (ms) | Throughput (Gcells/s)\n";
        std::cout << "-----------|-----------|----------------------\n";
    }
    
    std::vector<double> data(num_cells);
    for (int i = 0; i < num_cells; ++i) {
        data[i] = i * 0.001;
    }
    
    double* d_data;
    cudaMalloc(&d_data, num_cells * sizeof(double));
    cudaMemcpy(d_data, data.data(), num_cells * sizeof(double), cudaMemcpyHostToDevice);
    
    for (int blockSize : block_sizes) {
        int numBlocks = (num_cells + blockSize - 1) / blockSize;
        
        // Warm-up
        simpleKernel<<<numBlocks, blockSize>>>(d_data, num_cells);
        cudaDeviceSynchronize();
        
        // Benchmark
        auto start = std::chrono::high_resolution_clock::now();
        for (int iter = 0; iter < num_iterations; ++iter) {
            simpleKernel<<<numBlocks, blockSize>>>(d_data, num_cells);
        }
        cudaDeviceSynchronize();
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / (1000.0 * num_iterations);
        double throughput = num_cells / (time_ms * 1e6);  // Gcells/s
        
        if (rank == 0) {
            printf("%10d | %9.3f | %20.2f\n", blockSize, time_ms, throughput);
        }
    }
    
    cudaFree(d_data);
}

#else

// Placeholder tests when GPU is not available
TEST_F(GPUBenchmark, GPUNotAvailable) {
    GTEST_SKIP() << "GPU support not enabled. Rebuild with -DENABLE_CUDA=ON";
}

#endif  // ENABLE_CUDA
