# Examples CMakeLists
# All examples now run through config-driven executable with different config files

# Helper function to create examples
function(add_reservoir_example NAME SOURCE)
    add_executable(${NAME} ${SOURCE})
    target_link_libraries(${NAME} fsrmlib)
    
    # Install example data files
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${NAME}.DATA)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${NAME}.DATA 
                      ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.DATA COPYONLY)
    endif()
endfunction()

# Main config-driven executable - ALL examples run through this
add_reservoir_example(simulator ex_config_driven.cpp)

# SPE benchmarks (kept for validation/comparison)
add_reservoir_example(spe1 spe1.cpp)
add_reservoir_example(spe3 spe3.cpp)
add_reservoir_example(spe9 spe9.cpp)
add_reservoir_example(spe10 spe10.cpp)

# Copy all config files to build directory
file(GLOB CONFIG_FILES "${CMAKE_SOURCE_DIR}/config/*.config")
foreach(CONFIG_FILE ${CONFIG_FILES})
    get_filename_component(FILENAME ${CONFIG_FILE} NAME)
    configure_file(${CONFIG_FILE} ${CMAKE_CURRENT_BINARY_DIR}/config/${FILENAME} COPYONLY)
endforeach()

# Copy data files
file(GLOB DATA_FILES "*.DATA" "*.INC")
foreach(DATA_FILE ${DATA_FILES})
    get_filename_component(FILENAME ${DATA_FILE} NAME)
    configure_file(${DATA_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} COPYONLY)
endforeach()

# Display usage message
message(STATUS "===========================================================")
message(STATUS "Examples Build Configuration")
message(STATUS "===========================================================")
message(STATUS "Building config-driven simulator executable")
message(STATUS "")
message(STATUS "Usage:")
message(STATUS "  ./simulator -c config/<example>.config")
message(STATUS "")
message(STATUS "Available example configs:")
message(STATUS "  - config/single_phase.config")
message(STATUS "  - config/hydraulic_fracturing.config")
message(STATUS "  - config/lefm_fracture_growth.config")
message(STATUS "  - config/induced_seismicity.config")
message(STATUS "  - config/wave_propagation.config")
message(STATUS "  - config/geothermal.config")
message(STATUS "  - config/co2_storage.config")
message(STATUS "  - config/shale_reservoir.config")
message(STATUS "  - config/default.config")
message(STATUS "")
message(STATUS "SPE Benchmark configs:")
message(STATUS "  - config/spe1_benchmark.config")
message(STATUS "  - config/spe3_benchmark.config")
message(STATUS "  - config/spe9_benchmark.config")
message(STATUS "  - config/spe10_benchmark.config")
message(STATUS "")
message(STATUS "All examples now use HDF5 output by default")
message(STATUS "Use VTK for legacy visualization if needed")
message(STATUS "===========================================================")
