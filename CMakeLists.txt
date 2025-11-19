cmake_minimum_required(VERSION 3.15)
project(ReservoirSim VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(ENABLE_TESTING "Enable testing" ON)
option(ENABLE_MPI "Enable MPI support" ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Find required packages
find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

# Find PETSc
pkg_check_modules(PETSC REQUIRED PETSc>=3.15)

# Find HDF5 for output
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
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

# Create library
add_library(reservoirlib SHARED ${SOURCES})
target_link_libraries(reservoirlib 
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    m
)

# Main executable
add_executable(reservoirsim src/main.cpp)
target_link_libraries(reservoirsim reservoirlib)

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
install(TARGETS reservoirsim reservoirlib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
