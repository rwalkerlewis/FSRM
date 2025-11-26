cmake_minimum_required(VERSION 3.15)

# Set compilers before project() call if not already set
if(NOT DEFINED CMAKE_C_COMPILER AND EXISTS "/usr/bin/gcc")
    set(CMAKE_C_COMPILER "/usr/bin/gcc")
endif()
if(NOT DEFINED CMAKE_CXX_COMPILER AND EXISTS "/usr/bin/g++")
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
endif()

project(FSRM VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(ENABLE_TESTING "Enable testing" ON)
option(ENABLE_MPI "Enable MPI support" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(ENABLE_CUDA "Enable CUDA GPU support" OFF)
option(ENABLE_HIP "Enable AMD ROCm/HIP GPU support" OFF)

# Find required packages
# Set MPI compiler wrappers
if(NOT DEFINED MPI_C_COMPILER)
    find_program(MPI_C_COMPILER NAMES mpicc)
endif()
if(NOT DEFINED MPI_CXX_COMPILER)
    find_program(MPI_CXX_COMPILER NAMES mpicxx mpic++)
endif()

find_package(MPI REQUIRED)
find_package(PkgConfig)

# Find PETSc with multiple fallback methods
set(PETSC_FOUND FALSE)

# Method 1: Try pkg-config
if(PkgConfig_FOUND)
    pkg_check_modules(PETSC PETSc)
    if(PETSC_FOUND)
        message(STATUS "Found PETSc via pkg-config")
    endif()
endif()

# Method 2: Try environment variables
if(NOT PETSC_FOUND AND DEFINED ENV{PETSC_DIR})
    set(PETSC_DIR $ENV{PETSC_DIR})
    message(STATUS "Using PETSC_DIR from environment: ${PETSC_DIR}")
    
    if(DEFINED ENV{PETSC_ARCH})
        set(PETSC_ARCH $ENV{PETSC_ARCH})
        set(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include" "${PETSC_DIR}/${PETSC_ARCH}/include")
        set(PETSC_LIBRARY_DIRS "${PETSC_DIR}/${PETSC_ARCH}/lib")
        message(STATUS "Using PETSC_ARCH: ${PETSC_ARCH}")
    else()
        set(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include")
        set(PETSC_LIBRARY_DIRS "${PETSC_DIR}/lib")
    endif()
    
    set(PETSC_LIBRARIES "petsc")
    set(PETSC_FOUND TRUE)
    message(STATUS "PETSc configured from environment variables")
endif()

# Method 3: Search common locations
if(NOT PETSC_FOUND)
    message(STATUS "Searching for PETSc in common locations...")
    
    # Try to find petsc.h
    find_path(PETSC_INCLUDE_DIR 
        NAMES petsc.h petsc/petsc.h
        PATHS 
            /usr/include/petsc
            /usr/include
            /usr/local/include/petsc
            /usr/local/include
            /usr/lib/petsc/include
            /usr/lib/petscdir/petsc-3.18/x86_64-linux-gnu-real/include
            /usr/lib/petscdir/petsc-3.19/x86_64-linux-gnu-real/include
            /usr/lib/petscdir/petsc-3.20/x86_64-linux-gnu-real/include
        PATH_SUFFIXES petsc
        DOC "PETSc include directory"
    )
    
    # Try to find libpetsc
    find_library(PETSC_LIBRARY
        NAMES petsc
        PATHS
            /usr/lib
            /usr/lib/x86_64-linux-gnu
            /usr/local/lib
            /usr/lib/petscdir/petsc-3.18/x86_64-linux-gnu-real/lib
            /usr/lib/petscdir/petsc-3.19/x86_64-linux-gnu-real/lib
            /usr/lib/petscdir/petsc-3.20/x86_64-linux-gnu-real/lib
        DOC "PETSc library"
    )
    
    if(PETSC_INCLUDE_DIR AND PETSC_LIBRARY)
        set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR})
        set(PETSC_LIBRARIES ${PETSC_LIBRARY})
        set(PETSC_FOUND TRUE)
        message(STATUS "Found PETSc:")
        message(STATUS "  Include: ${PETSC_INCLUDE_DIR}")
        message(STATUS "  Library: ${PETSC_LIBRARY}")
    endif()
endif()

if(NOT PETSC_FOUND)
    message(FATAL_ERROR "PETSc not found. Please install petsc-dev or set PETSC_DIR environment variable.")
endif()

# Find HDF5 (make it optional)
if(PkgConfig_FOUND)
    pkg_check_modules(HDF5 hdf5-openmpi)
endif()

if(NOT HDF5_FOUND)
    find_package(HDF5 COMPONENTS C)
    if(HDF5_FOUND)
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
        set(HDF5_LIBRARIES ${HDF5_LIBRARIES})
        set(HDF5_LDFLAGS ${HDF5_LIBRARIES})
        message(STATUS "Found HDF5 via find_package")
    else()
        message(WARNING "HDF5 not found - continuing without it")
        set(HDF5_INCLUDE_DIRS "")
        set(HDF5_LIBRARIES "")
        set(HDF5_LDFLAGS "")
    endif()
endif()

# GPU support
if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit)
    if(CUDAToolkit_FOUND)
        message(STATUS "CUDA enabled: ${CUDAToolkit_VERSION}")
        set(CMAKE_CUDA_STANDARD 14)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        add_definitions(-DUSE_CUDA)
        
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
if(PETSC_LIBRARY_DIRS)
    link_directories(${PETSC_LIBRARY_DIRS})
endif()

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

# Print configuration summary
message(STATUS "")
message(STATUS "FSRM Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "  CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  MPI enabled: YES")
message(STATUS "  CUDA enabled: ${ENABLE_CUDA}")
message(STATUS "  HIP enabled: ${ENABLE_HIP}")
message(STATUS "  Testing enabled: ${ENABLE_TESTING}")
message(STATUS "  Examples enabled: ${BUILD_EXAMPLES}")
message(STATUS "  PETSc found: ${PETSC_FOUND}")
if(HDF5_FOUND)
    message(STATUS "  HDF5 found: YES")
else()
    message(STATUS "  HDF5 found: NO (optional)")
endif()
message(STATUS "")
