#ifndef GPU_MANAGER_HPP
#define GPU_MANAGER_HPP

#include <string>
#include <vector>
#include <memory>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

namespace FSRM {

enum class GPUBackend {
    NONE,
    CUDA,
    HIP
};

struct GPUDeviceInfo {
    int device_id;
    std::string name;
    size_t total_memory;
    size_t free_memory;
    int compute_capability_major;
    int compute_capability_minor;
    int multiprocessor_count;
    int max_threads_per_block;
    bool unified_memory_support;
};

class GPUManager {
public:
    static GPUManager& getInstance();
    
    // Initialization
    bool initialize();
    bool isAvailable() const { return gpu_available; }
    GPUBackend getBackend() const { return backend; }
    
    // Device management
    int getDeviceCount() const;
    bool setDevice(int device_id);
    int getCurrentDevice() const;
    GPUDeviceInfo getDeviceInfo(int device_id) const;
    void printDeviceInfo() const;
    
    // Memory management
    void* allocateDevice(size_t size);
    void freeDevice(void* ptr);
    void copyHostToDevice(void* dst, const void* src, size_t size);
    void copyDeviceToHost(void* dst, const void* src, size_t size);
    void copyDeviceToDevice(void* dst, const void* src, size_t size);
    
    // Synchronization
    void synchronize();
    void synchronizeStream(void* stream);
    
    // Error handling
    std::string getLastError() const;
    void checkError(const char* msg) const;
    
    // Performance
    void startTimer();
    float stopTimer(); // returns elapsed time in ms
    
    // Destructor
    ~GPUManager();
    
private:
    GPUManager();
    GPUManager(const GPUManager&) = delete;
    GPUManager& operator=(const GPUManager&) = delete;
    
    bool gpu_available;
    GPUBackend backend;
    int current_device;
    int device_count;
    std::vector<GPUDeviceInfo> devices;
    
#ifdef USE_CUDA
    cudaEvent_t start_event, stop_event;
#endif
};

// RAII wrapper for GPU memory
template<typename T>
class GPUArray {
public:
    GPUArray(size_t n) : size_(n), data_(nullptr) {
        data_ = static_cast<T*>(GPUManager::getInstance().allocateDevice(n * sizeof(T)));
    }
    
    ~GPUArray() {
        if (data_) {
            GPUManager::getInstance().freeDevice(data_);
        }
    }
    
    // Delete copy
    GPUArray(const GPUArray&) = delete;
    GPUArray& operator=(const GPUArray&) = delete;
    
    // Move semantics
    GPUArray(GPUArray&& other) noexcept : size_(other.size_), data_(other.data_) {
        other.data_ = nullptr;
        other.size_ = 0;
    }
    
    T* data() { return data_; }
    const T* data() const { return data_; }
    size_t size() const { return size_; }
    
    void copyFromHost(const T* host_data) {
        GPUManager::getInstance().copyHostToDevice(data_, host_data, size_ * sizeof(T));
    }
    
    void copyToHost(T* host_data) const {
        GPUManager::getInstance().copyDeviceToHost(host_data, data_, size_ * sizeof(T));
    }
    
private:
    size_t size_;
    T* data_;
};

} // namespace FSRM

#endif // GPU_MANAGER_HPP
