/**
 * @file test_memory_io_benchmarks.cpp
 * @brief Memory and I/O performance benchmarks
 * 
 * Benchmarks for:
 * - Memory allocation/deallocation
 * - Memory access patterns
 * - Cache efficiency
 * - File I/O operations
 * - HDF5 I/O performance
 */

#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <mpi.h>
#include <petsc.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

class MemoryIOBenchmark : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        MPI_Comm_size(PETSC_COMM_WORLD, &size);
    }
    
    void TearDown() override {
        // Clean up test files
        if (rank == 0) {
            std::remove("benchmark_test.dat");
            std::remove("benchmark_test.h5");
        }
    }
    
    int rank, size;
};

// ============================================================================
// Memory Allocation Benchmarks
// ============================================================================

TEST_F(MemoryIOBenchmark, MemoryAllocationPerformance) {
    std::vector<size_t> sizes = {
        1024,              // 1 KB
        1024 * 1024,       // 1 MB
        10 * 1024 * 1024,  // 10 MB
        100 * 1024 * 1024  // 100 MB
    };
    
    if (rank == 0) {
        std::cout << "\n=== Memory Allocation Performance ===\n";
        std::cout << "Size (MB) | Alloc (ms) | Dealloc (ms) | Total (ms)\n";
        std::cout << "----------|------------|--------------|------------\n";
    }
    
    for (size_t size : sizes) {
        size_t num_doubles = size / sizeof(double);
        
        // Allocation
        auto start_alloc = std::chrono::high_resolution_clock::now();
        double* data = new double[num_doubles];
        auto end_alloc = std::chrono::high_resolution_clock::now();
        
        // Initialize to prevent optimization
        for (size_t i = 0; i < std::min(num_doubles, size_t(1000)); ++i) {
            data[i] = i * 0.001;
        }
        
        // Deallocation
        auto start_dealloc = std::chrono::high_resolution_clock::now();
        delete[] data;
        auto end_dealloc = std::chrono::high_resolution_clock::now();
        
        double alloc_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_alloc - start_alloc).count() / 1000.0;
        double dealloc_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_dealloc - start_dealloc).count() / 1000.0;
        double total_ms = alloc_ms + dealloc_ms;
        
        if (rank == 0) {
            printf("%9.1f | %10.3f | %12.3f | %10.3f\n",
                   size / (1024.0 * 1024.0), alloc_ms, dealloc_ms, total_ms);
        }
    }
}

TEST_F(MemoryIOBenchmark, VectorReallocationPerformance) {
    if (rank == 0) {
        std::cout << "\n=== Vector Reallocation Performance ===\n";
        std::cout << "Final Size | Reserve | No Reserve | Speedup\n";
        std::cout << "-----------|---------|------------|--------\n";
    }
    
    std::vector<int> sizes = {1000, 10000, 100000, 1000000};
    
    for (int final_size : sizes) {
        // With reserve
        auto start_reserve = std::chrono::high_resolution_clock::now();
        std::vector<double> vec_reserve;
        vec_reserve.reserve(final_size);
        for (int i = 0; i < final_size; ++i) {
            vec_reserve.push_back(i * 0.001);
        }
        auto end_reserve = std::chrono::high_resolution_clock::now();
        
        // Without reserve
        auto start_no_reserve = std::chrono::high_resolution_clock::now();
        std::vector<double> vec_no_reserve;
        for (int i = 0; i < final_size; ++i) {
            vec_no_reserve.push_back(i * 0.001);
        }
        auto end_no_reserve = std::chrono::high_resolution_clock::now();
        
        double time_reserve = std::chrono::duration_cast<std::chrono::microseconds>(
            end_reserve - start_reserve).count() / 1000.0;
        double time_no_reserve = std::chrono::duration_cast<std::chrono::microseconds>(
            end_no_reserve - start_no_reserve).count() / 1000.0;
        double speedup = time_no_reserve / time_reserve;
        
        if (rank == 0) {
            printf("%10d | %7.3f | %10.3f | %6.2fx\n",
                   final_size, time_reserve, time_no_reserve, speedup);
        }
        
        EXPECT_GT(speedup, 1.0);  // Reserve should be faster
    }
}

// ============================================================================
// Cache Performance Benchmarks
// ============================================================================

TEST_F(MemoryIOBenchmark, CacheLinePerformance) {
    const int N = 10000000;  // 10M elements
    std::vector<double> data(N, 1.0);
    
    if (rank == 0) {
        std::cout << "\n=== Cache Performance ===\n";
        std::cout << "Stride | Time (ms) | Bandwidth (GB/s)\n";
        std::cout << "-------|-----------|------------------\n";
    }
    
    std::vector<int> strides = {1, 2, 4, 8, 16, 32, 64};
    
    for (int stride : strides) {
        double sum = 0.0;
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; i += stride) {
            sum += data[i];
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        double bytes_accessed = (N / stride) * sizeof(double);
        double bandwidth_gbps = (bytes_accessed / 1e9) / (duration.count() / 1e6);
        
        if (rank == 0) {
            printf("%6d | %9.3f | %16.2f\n", stride, time_ms, bandwidth_gbps);
        }
        
        EXPECT_GT(sum, 0.0);  // Use sum to prevent optimization
    }
}

TEST_F(MemoryIOBenchmark, MatrixAccessPattern) {
    const int N = 2000;
    std::vector<double> matrix(N * N, 1.0);
    
    if (rank == 0) {
        std::cout << "\n=== Matrix Access Patterns (2000x2000) ===\n";
        std::cout << "Pattern | Time (ms) | Bandwidth (GB/s)\n";
        std::cout << "--------|-----------|------------------\n";
    }
    
    // Row-major access
    double sum = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            sum += matrix[i * N + j];
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_row = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_row = duration_row.count() / 1000.0;
    double bytes_row = N * N * sizeof(double);
    double bw_row = (bytes_row / 1e9) / (duration_row.count() / 1e6);
    
    // Column-major access
    sum = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            sum += matrix[i * N + j];
        }
    }
    end = std::chrono::high_resolution_clock::now();
    auto duration_col = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_col = duration_col.count() / 1000.0;
    double bw_col = (bytes_row / 1e9) / (duration_col.count() / 1e6);
    
    // Random access
    std::random_device rd;
    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::uniform_int_distribution<> dis(0, N*N - 1);
    
    sum = 0.0;
    start = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < N * N; ++k) {
        int idx = dis(gen);
        sum += matrix[idx];
    }
    end = std::chrono::high_resolution_clock::now();
    auto duration_rand = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    double time_rand = duration_rand.count() / 1000.0;
    double bw_rand = (bytes_row / 1e9) / (duration_rand.count() / 1e6);
    
    if (rank == 0) {
        printf("Row-maj | %9.3f | %16.2f\n", time_row, bw_row);
        printf("Col-maj | %9.3f | %16.2f\n", time_col, bw_col);
        printf("Random  | %9.3f | %16.2f\n", time_rand, bw_rand);
    }
    
    EXPECT_GT(sum, 0.0);
}

// ============================================================================
// File I/O Benchmarks
// ============================================================================

TEST_F(MemoryIOBenchmark, BinaryFileIOPerformance) {
    std::vector<size_t> sizes = {
        1024 * 1024,       // 1 MB
        10 * 1024 * 1024,  // 10 MB
        100 * 1024 * 1024  // 100 MB
    };
    
    if (rank == 0) {
        std::cout << "\n=== Binary File I/O Performance ===\n";
        std::cout << "Size (MB) | Write (ms) | Read (ms) | Write BW (MB/s) | Read BW (MB/s)\n";
        std::cout << "----------|------------|-----------|-----------------|----------------\n";
    }
    
    for (size_t size : sizes) {
        size_t num_doubles = size / sizeof(double);
        std::vector<double> write_data(num_doubles);
        std::vector<double> read_data(num_doubles);
        
        // Initialize write data
        for (size_t i = 0; i < num_doubles; ++i) {
            write_data[i] = i * 0.001;
        }
        
        // Write performance
        auto start_write = std::chrono::high_resolution_clock::now();
        {
            std::ofstream ofs("benchmark_test.dat", std::ios::binary);
            ofs.write(reinterpret_cast<const char*>(write_data.data()), size);
        }
        auto end_write = std::chrono::high_resolution_clock::now();
        
        // Read performance
        auto start_read = std::chrono::high_resolution_clock::now();
        {
            std::ifstream ifs("benchmark_test.dat", std::ios::binary);
            ifs.read(reinterpret_cast<char*>(read_data.data()), size);
        }
        auto end_read = std::chrono::high_resolution_clock::now();
        
        double write_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_write - start_write).count() / 1000.0;
        double read_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_read - start_read).count() / 1000.0;
        
        double write_bw = (size / (1024.0 * 1024.0)) / (write_ms / 1000.0);
        double read_bw = (size / (1024.0 * 1024.0)) / (read_ms / 1000.0);
        
        if (rank == 0) {
            printf("%9.1f | %10.3f | %9.3f | %15.1f | %14.1f\n",
                   size / (1024.0 * 1024.0), write_ms, read_ms, write_bw, read_bw);
        }
        
        // Verify correctness
        double max_error = 0.0;
        for (size_t i = 0; i < std::min(num_doubles, size_t(1000)); ++i) {
            max_error = std::max(max_error, std::abs(write_data[i] - read_data[i]));
        }
        EXPECT_LT(max_error, 1e-10);
    }
}

#ifdef HAVE_HDF5

TEST_F(MemoryIOBenchmark, HDF5PerformanceSingleProcess) {
    if (size > 1) {
        GTEST_SKIP() << "HDF5 single process test - run with 1 MPI process";
    }
    
    std::vector<int> sizes = {100, 200, 400};  // Grid sizes
    
    if (rank == 0) {
        std::cout << "\n=== HDF5 I/O Performance ===\n";
        std::cout << "Grid | Cells | Write (ms) | Read (ms) | File Size (MB)\n";
        std::cout << "-----|-------|------------|-----------|----------------\n";
    }
    
    for (int n : sizes) {
        int num_cells = n * n * n;
        std::vector<double> data(num_cells);
        std::vector<double> read_data(num_cells);
        
        // Initialize
        for (int i = 0; i < num_cells; ++i) {
            data[i] = i * 0.001;
        }
        
        // Write
        auto start_write = std::chrono::high_resolution_clock::now();
        {
            hid_t file_id = H5Fcreate("benchmark_test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            
            hsize_t dims[1] = {static_cast<hsize_t>(num_cells)};
            hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
            hid_t dataset_id = H5Dcreate2(file_id, "/pressure", H5T_NATIVE_DOUBLE,
                                         dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
            
            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
            H5Fclose(file_id);
        }
        auto end_write = std::chrono::high_resolution_clock::now();
        
        // Read
        auto start_read = std::chrono::high_resolution_clock::now();
        {
            hid_t file_id = H5Fopen("benchmark_test.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
            hid_t dataset_id = H5Dopen2(file_id, "/pressure", H5P_DEFAULT);
            
            H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_data.data());
            
            H5Dclose(dataset_id);
            H5Fclose(file_id);
        }
        auto end_read = std::chrono::high_resolution_clock::now();
        
        double write_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_write - start_write).count() / 1000.0;
        double read_ms = std::chrono::duration_cast<std::chrono::microseconds>(
            end_read - start_read).count() / 1000.0;
        
        // Get file size
        std::ifstream file("benchmark_test.h5", std::ios::binary | std::ios::ate);
        double file_size_mb = file.tellg() / (1024.0 * 1024.0);
        file.close();
        
        if (rank == 0) {
            printf("%3dx%3dx%3d | %6d | %10.3f | %9.3f | %14.2f\n",
                   n, n, n, num_cells, write_ms, read_ms, file_size_mb);
        }
        
        // Verify
        double max_error = 0.0;
        for (int i = 0; i < std::min(num_cells, 1000); ++i) {
            max_error = std::max(max_error, std::abs(data[i] - read_data[i]));
        }
        EXPECT_LT(max_error, 1e-10);
    }
}

#endif  // HAVE_HDF5

// ============================================================================
// Memory Bandwidth Benchmarks
// ============================================================================

TEST_F(MemoryIOBenchmark, MemoryBandwidthCopy) {
    std::vector<size_t> sizes = {
        1024 * 1024,       // 1 MB
        10 * 1024 * 1024,  // 10 MB
        100 * 1024 * 1024  // 100 MB
    };
    
    if (rank == 0) {
        std::cout << "\n=== Memory Copy Bandwidth ===\n";
        std::cout << "Size (MB) | Time (ms) | Bandwidth (GB/s)\n";
        std::cout << "----------|-----------|------------------\n";
    }
    
    for (size_t size : sizes) {
        size_t num_doubles = size / sizeof(double);
        std::vector<double> src(num_doubles, 1.0);
        std::vector<double> dst(num_doubles);
        
        auto start = std::chrono::high_resolution_clock::now();
        std::copy(src.begin(), src.end(), dst.begin());
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        double bandwidth_gbps = (size / 1e9) / (duration.count() / 1e6);
        
        if (rank == 0) {
            printf("%9.1f | %9.3f | %16.2f\n",
                   size / (1024.0 * 1024.0), time_ms, bandwidth_gbps);
        }
        
        EXPECT_DOUBLE_EQ(dst[0], 1.0);
    }
}

TEST_F(MemoryIOBenchmark, MemoryBandwidthStreamTriad) {
    // Classic STREAM Triad: a[i] = b[i] + c[i] * scalar
    std::vector<size_t> sizes = {
        1024 * 1024,       // 1 MB
        10 * 1024 * 1024,  // 10 MB
        100 * 1024 * 1024  // 100 MB
    };
    
    if (rank == 0) {
        std::cout << "\n=== STREAM Triad Benchmark ===\n";
        std::cout << "Size (MB) | Time (ms) | Bandwidth (GB/s)\n";
        std::cout << "----------|-----------|------------------\n";
    }
    
    for (size_t size : sizes) {
        size_t num_doubles = size / sizeof(double);
        std::vector<double> a(num_doubles);
        std::vector<double> b(num_doubles, 1.0);
        std::vector<double> c(num_doubles, 2.0);
        double scalar = 3.0;
        
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < num_doubles; ++i) {
            a[i] = b[i] + c[i] * scalar;
        }
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        // 3 arrays * size bytes (2 reads + 1 write)
        double bandwidth_gbps = (3.0 * size / 1e9) / (duration.count() / 1e6);
        
        if (rank == 0) {
            printf("%9.1f | %9.3f | %16.2f\n",
                   size / (1024.0 * 1024.0), time_ms, bandwidth_gbps);
        }
        
        EXPECT_DOUBLE_EQ(a[0], 7.0);
    }
}
