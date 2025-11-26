# Tests CMakeLists
# Organized test suite with unit, functional, physics, integration, and performance tests

find_package(GTest)

# ============================================================================
# Define Test Source Files by Category
# ============================================================================

# Unit Tests - Test individual components in isolation
set(UNIT_TEST_SOURCES
    unit/test_config_reader.cpp
    unit/test_physics_kernels.cpp
    unit/test_well_model.cpp
    unit/test_fracture_model.cpp
    unit/test_eclipse_io.cpp
    unit/test_fluid_model.cpp
    unit/test_material_model.cpp
)

# Functional Tests - Test features working correctly
set(FUNCTIONAL_TEST_SOURCES
    functional/test_simulator_init.cpp
    functional/test_solver_convergence.cpp
    functional/test_boundary_conditions.cpp
    functional/test_well_operations.cpp
)

# Physics Tests (MMS) - Method of Manufactured Solutions tests
set(PHYSICS_TEST_SOURCES
    physics/test_mms_diffusion.cpp
    physics/test_mms_elasticity.cpp
    physics/test_mms_wave_propagation.cpp
    physics/test_analytical_solutions.cpp
)

# Integration Tests - Full-scale end-to-end tests
set(INTEGRATION_TEST_SOURCES
    integration/test_full_simulation.cpp
    integration/test_restart.cpp
)

# Performance Tests - Benchmarks and scaling tests
set(PERFORMANCE_TEST_SOURCES
    performance/test_benchmarks.cpp
    performance/test_scaling.cpp
)

# Legacy Tests (kept for backward compatibility)
set(LEGACY_TEST_SOURCES
    test_mms_simple.cpp
    test_physics_detailed.cpp
)

# ============================================================================
# Build Test Executables
# ============================================================================

if(GTest_FOUND)
    # Main comprehensive test executable
    add_executable(run_all_tests
        test_main_gtest.cpp
        ${UNIT_TEST_SOURCES}
        ${FUNCTIONAL_TEST_SOURCES}
        ${PHYSICS_TEST_SOURCES}
        ${INTEGRATION_TEST_SOURCES}
        ${PERFORMANCE_TEST_SOURCES}
    )
    
    target_link_libraries(run_all_tests
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        GTest::GTest
    )
    
    target_include_directories(run_all_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    # Category-specific executables for faster focused testing
    
    # Unit tests executable
    add_executable(run_unit_tests
        test_main_gtest.cpp
        ${UNIT_TEST_SOURCES}
    )
    target_link_libraries(run_unit_tests
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        GTest::GTest
    )
    target_include_directories(run_unit_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    # Functional tests executable
    add_executable(run_functional_tests
        test_main_gtest.cpp
        ${FUNCTIONAL_TEST_SOURCES}
    )
    target_link_libraries(run_functional_tests
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        GTest::GTest
    )
    target_include_directories(run_functional_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    # Physics/MMS tests executable
    add_executable(run_physics_tests
        test_main_gtest.cpp
        ${PHYSICS_TEST_SOURCES}
    )
    target_link_libraries(run_physics_tests
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        GTest::GTest
    )
    target_include_directories(run_physics_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    # Integration tests executable
    add_executable(run_integration_tests
        test_main_gtest.cpp
        ${INTEGRATION_TEST_SOURCES}
    )
    target_link_libraries(run_integration_tests
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        GTest::GTest
    )
    target_include_directories(run_integration_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    # Performance tests executable
    add_executable(run_performance_tests
        test_main_gtest.cpp
        ${PERFORMANCE_TEST_SOURCES}
    )
    target_link_libraries(run_performance_tests
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
        GTest::GTest
    )
    target_include_directories(run_performance_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    # ============================================================================
    # Register Tests with CTest
    # ============================================================================
    
    # Quick unit tests (run first, fast feedback)
    add_test(NAME Unit.ConfigReader 
        COMMAND run_unit_tests --gtest_filter=ConfigReaderTest.*)
    add_test(NAME Unit.PhysicsKernels 
        COMMAND run_unit_tests --gtest_filter=*KernelTest.*)
    add_test(NAME Unit.WellModel 
        COMMAND run_unit_tests --gtest_filter=WellModelTest.*)
    add_test(NAME Unit.FractureModel 
        COMMAND run_unit_tests --gtest_filter=*FractureModelTest.*:*FractureNetworkTest.*)
    add_test(NAME Unit.EclipseIO 
        COMMAND run_unit_tests --gtest_filter=EclipseIOTest.*)
    add_test(NAME Unit.FluidModel
        COMMAND run_unit_tests --gtest_filter=FluidModelTest.*)
    add_test(NAME Unit.MaterialModel
        COMMAND run_unit_tests --gtest_filter=MaterialModelTest.*)
    
    # Functional tests
    add_test(NAME Functional.SimulatorInit 
        COMMAND run_functional_tests --gtest_filter=SimulatorInitTest.*)
    add_test(NAME Functional.SolverConvergence 
        COMMAND run_functional_tests --gtest_filter=SolverConvergenceTest.*)
    add_test(NAME Functional.BoundaryConditions 
        COMMAND run_functional_tests --gtest_filter=BoundaryConditionTest.*)
    add_test(NAME Functional.WellOperations 
        COMMAND run_functional_tests --gtest_filter=WellOperationsTest.*)
    
    # Physics/MMS tests
    add_test(NAME Physics.MMS.Diffusion 
        COMMAND run_physics_tests --gtest_filter=MMSDiffusionTest.*)
    add_test(NAME Physics.MMS.Elasticity 
        COMMAND run_physics_tests --gtest_filter=MMSElasticityTest.*)
    add_test(NAME Physics.MMS.WavePropagation 
        COMMAND run_physics_tests --gtest_filter=MMSWavePropagationTest.*)
    add_test(NAME Physics.AnalyticalSolutions 
        COMMAND run_physics_tests --gtest_filter=AnalyticalSolutionsTest.*)
    
    # Integration tests (longer running)
    add_test(NAME Integration.FullSimulation 
        COMMAND run_integration_tests --gtest_filter=FullSimulationTest.*)
    add_test(NAME Integration.Restart 
        COMMAND run_integration_tests --gtest_filter=RestartTest.*)
    
    # Performance tests
    add_test(NAME Performance.Benchmarks 
        COMMAND run_performance_tests --gtest_filter=BenchmarkTest.*)
    add_test(NAME Performance.Scaling 
        COMMAND run_performance_tests --gtest_filter=ScalingTest.*)
    
    # ============================================================================
    # Set Test Properties (Timeouts, Labels)
    # ============================================================================
    
    # Unit tests - should be fast
    set_tests_properties(
        Unit.ConfigReader
        Unit.PhysicsKernels
        Unit.WellModel
        Unit.FractureModel
        Unit.EclipseIO
        Unit.FluidModel
        Unit.MaterialModel
        PROPERTIES 
            TIMEOUT 60
            LABELS "unit;quick"
    )
    
    # Functional tests - moderate time
    set_tests_properties(
        Functional.SimulatorInit
        Functional.SolverConvergence
        Functional.BoundaryConditions
        Functional.WellOperations
        PROPERTIES 
            TIMEOUT 120
            LABELS "functional"
    )
    
    # Physics tests - moderate time
    set_tests_properties(
        Physics.MMS.Diffusion
        Physics.MMS.Elasticity
        Physics.MMS.WavePropagation
        Physics.AnalyticalSolutions
        PROPERTIES 
            TIMEOUT 180
            LABELS "physics;mms"
    )
    
    # Integration tests - longer time
    set_tests_properties(
        Integration.FullSimulation
        Integration.Restart
        PROPERTIES 
            TIMEOUT 600
            LABELS "integration;slow"
    )
    
    # Performance tests - longest time
    set_tests_properties(
        Performance.Benchmarks
        Performance.Scaling
        PROPERTIES 
            TIMEOUT 900
            LABELS "performance;slow"
    )
    
else()
    message(WARNING "GTest not found. Building legacy test executable only.")
    
    # Legacy test executable (no GTest)
    add_executable(run_tests
        test_main.cpp
        ${LEGACY_TEST_SOURCES}
    )
    
    target_link_libraries(run_tests 
        fsrmlib
        ${PETSC_LIBRARIES}
        ${MPI_CXX_LIBRARIES}
    )
    
    target_include_directories(run_tests PRIVATE 
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
    )
    
    add_test(NAME LegacyTests COMMAND run_tests)
    set_tests_properties(LegacyTests PROPERTIES TIMEOUT 300)
endif()

# ============================================================================
# Custom Test Targets
# ============================================================================

add_custom_target(test_quick
    COMMAND ${CMAKE_CTEST_COMMAND} -L "quick" --output-on-failure
    DEPENDS run_unit_tests
    COMMENT "Running quick unit tests"
)

add_custom_target(test_unit
    COMMAND ${CMAKE_CTEST_COMMAND} -L "unit" --output-on-failure
    DEPENDS run_unit_tests
    COMMENT "Running all unit tests"
)

add_custom_target(test_functional
    COMMAND ${CMAKE_CTEST_COMMAND} -L "functional" --output-on-failure
    DEPENDS run_functional_tests
    COMMENT "Running functional tests"
)

add_custom_target(test_physics
    COMMAND ${CMAKE_CTEST_COMMAND} -L "physics" --output-on-failure
    DEPENDS run_physics_tests
    COMMENT "Running physics/MMS tests"
)

add_custom_target(test_integration
    COMMAND ${CMAKE_CTEST_COMMAND} -L "integration" --output-on-failure
    DEPENDS run_integration_tests
    COMMENT "Running integration tests"
)

add_custom_target(test_performance
    COMMAND ${CMAKE_CTEST_COMMAND} -L "performance" --output-on-failure
    DEPENDS run_performance_tests
    COMMENT "Running performance tests"
)

add_custom_target(test_all
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS run_all_tests
    COMMENT "Running all tests"
)

# ============================================================================
# Install Test Data
# ============================================================================

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
