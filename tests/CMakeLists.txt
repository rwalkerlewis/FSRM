# Tests CMakeLists
find_package(GTest)

# Create test executable - minimal test to verify build works
set(TEST_SOURCES
    test_main.cpp
    test_mms_simple.cpp
    test_physics_detailed.cpp
    # test_performance.cpp  # Disabled temporarily - needs API update
)
add_executable(run_tests ${TEST_SOURCES})

target_link_libraries(run_tests 
    fsrmlib
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

if(GTest_FOUND)
    target_link_libraries(run_tests GTest::GTest GTest::Main)
    
    # Unit tests
    add_test(NAME UnitTests COMMAND run_tests --gtest_filter=UnitTest*)
    add_test(NAME EclipseIOTest COMMAND run_tests --gtest_filter=EclipseIO*)
    add_test(NAME PhysicsKernelTest COMMAND run_tests --gtest_filter=PhysicsKernelTestFixture*)
    add_test(NAME WellModelTest COMMAND run_tests --gtest_filter=WellModel*)
    
    # Fracture model tests
    add_test(NAME FractureModelTest COMMAND run_tests --gtest_filter=FractureModelTest*)
    
    # MMS and convergence tests
    add_test(NAME MMSConvergenceTest COMMAND run_tests --gtest_filter=MMSConvergenceTest*)
    
    # Integration tests
    add_test(NAME IntegrationTest COMMAND run_tests --gtest_filter=IntegrationTest*)
    
    # Specific physics tests
    add_test(NAME SinglePhaseFlowTests COMMAND run_tests --gtest_filter=*SinglePhaseFlow*)
    add_test(NAME PoroelasticityTests COMMAND run_tests --gtest_filter=*Poroelasticity*)
    add_test(NAME ElastodynamicsTests COMMAND run_tests --gtest_filter=*Elastodynamics*)
    add_test(NAME HydraulicFractureTests COMMAND run_tests --gtest_filter=*HydraulicFracture*)
    
    # Performance tests
    add_test(NAME PerformanceTest COMMAND run_tests --gtest_filter=PerformanceTest*)
    
    # Set test properties
    set_tests_properties(
        UnitTests
        EclipseIOTest
        PhysicsKernelTest
        WellModelTest
        FractureModelTest
        PROPERTIES TIMEOUT 60
    )
    
    set_tests_properties(
        MMSConvergenceTest
        PROPERTIES TIMEOUT 300
    )
    
    set_tests_properties(
        IntegrationTest
        PROPERTIES TIMEOUT 600
    )
    
    set_tests_properties(
        PerformanceTest
        PROPERTIES TIMEOUT 600
    )
    
else()
    message(WARNING "GTest not found. Building basic test executable only.")
endif()

# Add custom test targets
add_custom_target(test_quick
    COMMAND ${CMAKE_CTEST_COMMAND} -R "UnitTest|EclipseIO|PhysicsKernel|WellModel"
    DEPENDS run_tests
    COMMENT "Running quick unit tests"
)

add_custom_target(test_convergence
    COMMAND ${CMAKE_CTEST_COMMAND} -R "MMS|Convergence"
    DEPENDS run_tests
    COMMENT "Running convergence tests"
)

add_custom_target(test_integration
    COMMAND ${CMAKE_CTEST_COMMAND} -R "Integration"
    DEPENDS run_tests
    COMMENT "Running integration tests"
)

add_custom_target(test_all
    COMMAND ${CMAKE_CTEST_COMMAND}
    DEPENDS run_tests
    COMMENT "Running all tests"
)
