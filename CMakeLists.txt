cmake_minimum_required(VERSION 3.15)
project(FSRM VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Use g++ if available
if(EXISTS "/usr/bin/g++")
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
endif()

# Options
option(ENABLE_TESTING "Enable testing" ON)
option(ENABLE_MPI "Enable MPI support" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(ENABLE_CUDA "Enable CUDA GPU support" OFF)
option(ENABLE_HIP "Enable AMD ROCm/HIP GPU support" OFF)

# Find required packages
# Set MPI paths if they exist
if(EXISTS "/usr/lib/x86_64-linux-gnu/openmpi")
    set(MPI_HOME "/usr/lib/x86_64-linux-gnu/openmpi")
    set(MPI_C_INCLUDE_PATH "/usr/lib/x86_64-linux-gnu/openmpi/include")
    set(MPI_CXX_INCLUDE_PATH "/usr/lib/x86_64-linux-gnu/openmpi/include")
endif()
set(MPI_C_COMPILER mpicc)
set(MPI_CXX_COMPILER mpicxx)
find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

# Find PETSc
pkg_check_modules(PETSC REQUIRED PETSc>=3.15)

# Find HDF5 for output (use pkg-config)
pkg_check_modules(HDF5 REQUIRED hdf5-openmpi)

# GPU support
if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit)
    if(CUDAToolkit_FOUND)
        message(STATUS "CUDA enabled: ${CUDAToolkit_VERSION}")
        set(CMAKE_CUDA_STANDARD 14)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        add_definitions(-DUSE_CUDA)
        
        # Set CUDA architectures (adjust based on your GPU)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set(CMAKE_CUDA_ARCHITECTURES 60 70 75 80 86)
        endif()
    else()
        message(WARNING "CUDA requested but not found. Building CPU-only version.")
        set(ENABLE_CUDA OFF)
    endif()
endif()

if(ENABLE_HIP)
    find_package(hip)
    if(hip_FOUND)
        message(STATUS "HIP/ROCm enabled")
        add_definitions(-DUSE_HIP)
    else()
        message(WARNING "HIP requested but not found.")
        set(ENABLE_HIP OFF)
    endif()
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
    ${PETSC_INCLUDE_DIRS}
    ${MPI_CXX_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# Link directories
link_directories(
    ${PETSC_LIBRARY_DIRS}
)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.hpp")

# GPU source files
if(ENABLE_CUDA)
    file(GLOB_RECURSE CUDA_SOURCES "src/*.cu")
    set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CUDA)
    list(APPEND SOURCES ${CUDA_SOURCES})
endif()

# Create library
add_library(fsrmlib SHARED ${SOURCES})
target_link_libraries(fsrmlib 
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${HDF5_LDFLAGS}
    m
)

# Link CUDA libraries if enabled
if(ENABLE_CUDA AND CUDAToolkit_FOUND)
    target_link_libraries(fsrmlib CUDA::cudart CUDA::cublas CUDA::cusparse)
    set_target_properties(fsrmlib PROPERTIES 
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )
endif()

if(ENABLE_HIP AND hip_FOUND)
    target_link_libraries(fsrmlib hip::host)
endif()

# Main executable
add_executable(fsrm src/main.cpp)
target_link_libraries(fsrm fsrmlib)

# Testing
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Install
install(TARGETS fsrm fsrmlib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
