#include "GPUManager.hpp"
#include <iostream>
#include <stdexcept>
#include <cstring>

#ifdef USE_CUDA
#include <cuda_runtime.h>

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error: " << cudaGetErrorString(err) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            throw std::runtime_error(cudaGetErrorString(err)); \
        } \
    } while(0)

#define CUDA_CHECK_LAST() \
    do { \
        cudaError_t err = cudaGetLastError(); \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA Error: " << cudaGetErrorString(err) \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        } \
    } while(0)

#endif

namespace ResSim {

GPUManager::GPUManager() 
    : gpu_available(false), 
      backend(GPUBackend::NONE),
      current_device(-1),
      device_count(0) {
#ifdef USE_CUDA
    CUDA_CHECK(cudaEventCreate(&start_event));
    CUDA_CHECK(cudaEventCreate(&stop_event));
#endif
}

GPUManager::~GPUManager() {
#ifdef USE_CUDA
    if (gpu_available) {
        cudaEventDestroy(start_event);
        cudaEventDestroy(stop_event);
    }
#endif
}

GPUManager& GPUManager::getInstance() {
    static GPUManager instance;
    return instance;
}

bool GPUManager::initialize() {
#ifdef USE_CUDA
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0) {
        std::cerr << "No CUDA devices found or CUDA error: " 
                  << cudaGetErrorString(err) << std::endl;
        gpu_available = false;
        backend = GPUBackend::NONE;
        return false;
    }
    
    backend = GPUBackend::CUDA;
    gpu_available = true;
    
    // Query device properties
    devices.clear();
    for (int i = 0; i < device_count; ++i) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
        
        GPUDeviceInfo info;
        info.device_id = i;
        info.name = prop.name;
        info.total_memory = prop.totalGlobalMem;
        info.compute_capability_major = prop.major;
        info.compute_capability_minor = prop.minor;
        info.multiprocessor_count = prop.multiProcessorCount;
        info.max_threads_per_block = prop.maxThreadsPerBlock;
        info.unified_memory_support = (prop.unifiedAddressing != 0);
        
        // Get free memory
        size_t free_mem, total_mem;
        CUDA_CHECK(cudaSetDevice(i));
        CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
        info.free_memory = free_mem;
        
        devices.push_back(info);
    }
    
    // Set default device (device 0)
    if (device_count > 0) {
        setDevice(0);
    }
    
    return true;
    
#elif defined(USE_HIP)
    // HIP implementation would go here
    backend = GPUBackend::HIP;
    gpu_available = true;
    return true;
    
#else
    gpu_available = false;
    backend = GPUBackend::NONE;
    return false;
#endif
}

int GPUManager::getDeviceCount() const {
    return device_count;
}

bool GPUManager::setDevice(int device_id) {
#ifdef USE_CUDA
    if (device_id < 0 || device_id >= device_count) {
        std::cerr << "Invalid device ID: " << device_id << std::endl;
        return false;
    }
    
    CUDA_CHECK(cudaSetDevice(device_id));
    current_device = device_id;
    return true;
#else
    return false;
#endif
}

int GPUManager::getCurrentDevice() const {
    return current_device;
}

GPUDeviceInfo GPUManager::getDeviceInfo(int device_id) const {
    if (device_id < 0 || device_id >= static_cast<int>(devices.size())) {
        throw std::out_of_range("Invalid device ID");
    }
    return devices[device_id];
}

void GPUManager::printDeviceInfo() const {
    if (!gpu_available) {
        std::cout << "No GPU devices available" << std::endl;
        return;
    }
    
    std::cout << "\n=== GPU Device Information ===" << std::endl;
    std::cout << "Backend: " << (backend == GPUBackend::CUDA ? "CUDA" : 
                                 backend == GPUBackend::HIP ? "HIP" : "NONE") << std::endl;
    std::cout << "Number of devices: " << device_count << std::endl;
    
    for (const auto& dev : devices) {
        std::cout << "\nDevice " << dev.device_id << ": " << dev.name << std::endl;
        std::cout << "  Compute Capability: " << dev.compute_capability_major 
                  << "." << dev.compute_capability_minor << std::endl;
        std::cout << "  Total Memory: " << dev.total_memory / (1024*1024) << " MB" << std::endl;
        std::cout << "  Free Memory: " << dev.free_memory / (1024*1024) << " MB" << std::endl;
        std::cout << "  Multiprocessors: " << dev.multiprocessor_count << std::endl;
        std::cout << "  Max Threads/Block: " << dev.max_threads_per_block << std::endl;
        std::cout << "  Unified Memory: " << (dev.unified_memory_support ? "Yes" : "No") << std::endl;
    }
    std::cout << "==============================\n" << std::endl;
}

void* GPUManager::allocateDevice(size_t size) {
#ifdef USE_CUDA
    void* ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&ptr, size));
    return ptr;
#else
    return nullptr;
#endif
}

void GPUManager::freeDevice(void* ptr) {
#ifdef USE_CUDA
    if (ptr) {
        CUDA_CHECK(cudaFree(ptr));
    }
#endif
}

void GPUManager::copyHostToDevice(void* dst, const void* src, size_t size) {
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
#endif
}

void GPUManager::copyDeviceToHost(void* dst, const void* src, size_t size) {
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost));
#endif
}

void GPUManager::copyDeviceToDevice(void* dst, const void* src, size_t size) {
#ifdef USE_CUDA
    CUDA_CHECK(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
#endif
}

void GPUManager::synchronize() {
#ifdef USE_CUDA
    CUDA_CHECK(cudaDeviceSynchronize());
#endif
}

void GPUManager::synchronizeStream(void* stream) {
#ifdef USE_CUDA
    CUDA_CHECK(cudaStreamSynchronize(static_cast<cudaStream_t>(stream)));
#endif
}

std::string GPUManager::getLastError() const {
#ifdef USE_CUDA
    cudaError_t err = cudaGetLastError();
    return std::string(cudaGetErrorString(err));
#else
    return "No GPU support compiled";
#endif
}

void GPUManager::checkError(const char* msg) const {
#ifdef USE_CUDA
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error after " << msg << ": " 
                  << cudaGetErrorString(err) << std::endl;
    }
#endif
}

void GPUManager::startTimer() {
#ifdef USE_CUDA
    CUDA_CHECK(cudaEventRecord(start_event));
#endif
}

float GPUManager::stopTimer() {
#ifdef USE_CUDA
    CUDA_CHECK(cudaEventRecord(stop_event));
    CUDA_CHECK(cudaEventSynchronize(stop_event));
    
    float elapsed_ms = 0.0f;
    CUDA_CHECK(cudaEventElapsedTime(&elapsed_ms, start_event, stop_event));
    return elapsed_ms;
#else
    return 0.0f;
#endif
}

} // namespace ResSim
