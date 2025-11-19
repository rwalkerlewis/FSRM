# Tests CMakeLists
find_package(GTest)

# Create test executable
file(GLOB TEST_SOURCES "*.cpp")
add_executable(run_tests ${TEST_SOURCES})

target_link_libraries(run_tests 
    reservoirlib
    ${PETSC_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

if(GTest_FOUND)
    target_link_libraries(run_tests GTest::GTest GTest::Main)
endif()

# Add tests
add_test(NAME UnitTests COMMAND run_tests)
add_test(NAME EclipseIOTest COMMAND run_tests --gtest_filter=EclipseIO*)
add_test(NAME PhysicsKernelTest COMMAND run_tests --gtest_filter=PhysicsKernel*)
add_test(NAME WellModelTest COMMAND run_tests --gtest_filter=WellModel*)

# MMS convergence tests
add_test(NAME MMSConvergenceTest COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_mms_tests.sh)
