# Examples CMakeLists

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

# Basic examples
add_reservoir_example(ex01_single_phase ex01_single_phase.cpp)
# add_reservoir_example(ex02_black_oil ex02_black_oil.cpp)  # TODO: Create
# add_reservoir_example(ex03_compositional ex03_compositional.cpp)  # TODO: Create
# add_reservoir_example(ex04_geomechanics ex04_geomechanics.cpp)  # TODO: Create

# SPE benchmarks
add_reservoir_example(spe1 spe1.cpp)
# add_reservoir_example(spe3 spe3.cpp)  # TODO: Create
# add_reservoir_example(spe9 spe9.cpp)  # TODO: Create
# add_reservoir_example(spe10 spe10.cpp)  # TODO: Create

# Advanced examples
add_reservoir_example(ex_hydraulic_fracturing ex_hydraulic_fracturing.cpp)
add_reservoir_example(ex_induced_seismicity ex_induced_seismicity.cpp)
# add_reservoir_example(ex_thermal_recovery ex_thermal_recovery.cpp)  # TODO: Create
add_reservoir_example(ex_stochastic_reservoir ex_stochastic_reservoir.cpp)
add_reservoir_example(ex_config_driven ex_config_driven.cpp)
add_reservoir_example(ex_lefm_fracture_growth ex_lefm_fracture_growth.cpp)
add_reservoir_example(ex_wave_propagation ex_wave_propagation.cpp)
add_reservoir_example(ex_buckley_leverett_2d ex_buckley_leverett_2d.cpp)
add_reservoir_example(ex_coupled_reservoir_2d ex_coupled_reservoir_2d.cpp)
add_reservoir_example(ex_reservoir_2d_vertical ex_reservoir_2d_vertical.cpp)
add_reservoir_example(ex_reservoir_2d_vertical_enhanced ex_reservoir_2d_vertical_enhanced.cpp)
add_reservoir_example(ex_gnuplot_solver ex_gnuplot_solver.cpp)

# Copy data files
file(GLOB DATA_FILES "*.DATA" "*.INC")
foreach(DATA_FILE ${DATA_FILES})
    get_filename_component(FILENAME ${DATA_FILE} NAME)
    configure_file(${DATA_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} COPYONLY)
endforeach()
