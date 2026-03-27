# =============================================================================
# FSRM Test Suite — CMakeLists.txt
# =============================================================================
# Organized by: unit, functional, physics_validation, integration,
#               performance, experimental
# =============================================================================

find_package(GTest)

# =============================================================================
# Source File Lists by Category
# =============================================================================

set(UNIT_TEST_SOURCES
    # Core
    unit/core/test_config_reader.cpp
    unit/core/test_physics_module_registry.cpp
    unit/core/test_simulation_config.cpp
    # Physics kernels
    unit/physics/test_physics_kernel_base.cpp
    unit/physics/test_single_phase_flow.cpp
    unit/physics/test_geomechanics.cpp
    unit/physics/test_thermal.cpp
    unit/physics/test_elastodynamics.cpp
    unit/physics/test_poroelastodynamics.cpp
    unit/physics/test_boundary_conditions.cpp
    # Domain
    unit/domain/test_well_model.cpp
    unit/domain/test_fracture_model.cpp
    unit/domain/test_material_model.cpp
    unit/domain/test_seismic_source.cpp
    unit/domain/test_explosion_source.cpp
    unit/domain/test_co2_properties.cpp
    unit/domain/test_flash_calculation.cpp
    unit/domain/test_residual_trapping.cpp
    # I/O
    unit/io/test_eclipse_io.cpp
    # ML
    unit/ml/test_ml_registry.cpp
    # Utilities
    unit/util/test_unit_system.cpp
    unit/util/test_coordinate_system.cpp
)

set(FUNCTIONAL_TEST_SOURCES
    functional/test_simulator_init.cpp
    functional/test_module_lifecycle.cpp
    functional/test_solver_convergence.cpp
    functional/test_well_operations.cpp
)

set(PHYSICS_VALIDATION_TEST_SOURCES
    physics_validation/test_mms_diffusion.cpp
    physics_validation/test_mms_elasticity.cpp
    physics_validation/test_mms_wave_propagation.cpp
    physics_validation/test_analytical_solutions.cpp
    physics_validation/test_thermoporoelastic.cpp
)

set(INTEGRATION_TEST_SOURCES
    integration/test_full_simulation.cpp
    integration/test_restart.cpp
    integration/test_coupled_physics.cpp
)

set(PERFORMANCE_TEST_SOURCES
    performance/test_benchmarks.cpp
    performance/test_scaling.cpp
    performance/test_memory.cpp
)

set(EXPERIMENTAL_TEST_SOURCES
    experimental/test_neural_stubs.cpp
)

# =============================================================================
# Build Test Executables
# =============================================================================

if(GTest_FOUND)

    # ---------- Helper macro for test executables ----------
    macro(fsrm_add_test_executable TARGET_NAME)
        add_executable(${TARGET_NAME} test_main_gtest.cpp ${ARGN})
        target_link_libraries(${TARGET_NAME}
            fsrmlib
            ${PETSC_LIBRARIES}
            ${MPI_CXX_LIBRARIES}
            GTest::GTest
        )
        target_include_directories(${TARGET_NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/include
            ${CMAKE_SOURCE_DIR}/tests
        )
        target_compile_definitions(${TARGET_NAME} PRIVATE FSRM_TEST_NO_MAIN)
    endmacro()

    # ---------- Per-category executables ----------
    fsrm_add_test_executable(run_unit_tests          ${UNIT_TEST_SOURCES})
    fsrm_add_test_executable(run_functional_tests     ${FUNCTIONAL_TEST_SOURCES})
    fsrm_add_test_executable(run_physics_tests        ${PHYSICS_VALIDATION_TEST_SOURCES})
    fsrm_add_test_executable(run_integration_tests    ${INTEGRATION_TEST_SOURCES})
    fsrm_add_test_executable(run_performance_tests    ${PERFORMANCE_TEST_SOURCES})
    fsrm_add_test_executable(run_experimental_tests   ${EXPERIMENTAL_TEST_SOURCES})

    # ---------- All-in-one executable ----------
    fsrm_add_test_executable(run_all_tests
        ${UNIT_TEST_SOURCES}
        ${FUNCTIONAL_TEST_SOURCES}
        ${PHYSICS_VALIDATION_TEST_SOURCES}
        ${INTEGRATION_TEST_SOURCES}
        ${PERFORMANCE_TEST_SOURCES}
        ${EXPERIMENTAL_TEST_SOURCES}
    )

    # =========================================================================
    # CTest Registrations
    # =========================================================================

    # --- Unit tests ---
    add_test(NAME Unit.ConfigReader           COMMAND run_unit_tests --gtest_filter=ConfigReaderTest.*)
    add_test(NAME Unit.PhysicsModuleRegistry  COMMAND run_unit_tests --gtest_filter=PhysicsModuleRegistryTest.*)
    add_test(NAME Unit.SimulationConfig       COMMAND run_unit_tests --gtest_filter=SimulationConfigTest.*)
    add_test(NAME Unit.PhysicsKernelBase      COMMAND run_unit_tests --gtest_filter=PhysicsKernelBaseTest.*)
    add_test(NAME Unit.SinglePhaseFlow        COMMAND run_unit_tests --gtest_filter=SinglePhaseFlowKernelTest.*)
    add_test(NAME Unit.Geomechanics           COMMAND run_unit_tests --gtest_filter=GeomechanicsKernelTest.*)
    add_test(NAME Unit.Thermal                COMMAND run_unit_tests --gtest_filter=ThermalKernelTest.*)
    add_test(NAME Unit.Elastodynamics         COMMAND run_unit_tests --gtest_filter=ElastodynamicsKernelTest.*)
    add_test(NAME Unit.Poroelastodynamics     COMMAND run_unit_tests --gtest_filter=PoroelastodynamicsKernelTest.*)
    add_test(NAME Unit.BoundaryConditions     COMMAND run_unit_tests --gtest_filter=*BoundaryCondition*.*)
    add_test(NAME Unit.WellModel              COMMAND run_unit_tests --gtest_filter=WellModelTest.*)
    add_test(NAME Unit.FractureModel          COMMAND run_unit_tests --gtest_filter=*FractureModel*.*)
    add_test(NAME Unit.MaterialModel          COMMAND run_unit_tests --gtest_filter=MaterialModelTest.*)
    add_test(NAME Unit.SeismicSource          COMMAND run_unit_tests --gtest_filter=*SeismicSource*.*)
    add_test(NAME Unit.ExplosionSource        COMMAND run_unit_tests --gtest_filter=ExplosionSourceKernelTest.*)
    add_test(NAME Unit.EclipseIO              COMMAND run_unit_tests --gtest_filter=EclipseIOTest.*)
    add_test(NAME Unit.MLRegistry             COMMAND run_unit_tests --gtest_filter=MLModelRegistryTest.*)
    add_test(NAME Unit.UnitSystem             COMMAND run_unit_tests --gtest_filter=UnitSystemTest.*)
    add_test(NAME Unit.CoordinateSystem       COMMAND run_unit_tests --gtest_filter=CoordinateSystemTest.*)
    add_test(NAME Unit.CO2Properties          COMMAND run_unit_tests --gtest_filter=CO2PropertiesTest.*)
    add_test(NAME Unit.FlashCalculation       COMMAND run_unit_tests --gtest_filter=FlashCalculationTest.*)
    add_test(NAME Unit.ResidualTrapping       COMMAND run_unit_tests --gtest_filter=ResidualTrappingTest.*)

    # --- Functional tests ---
    add_test(NAME Functional.SimulatorInit    COMMAND run_functional_tests --gtest_filter=SimulatorInitTest.*)
    add_test(NAME Functional.ModuleLifecycle  COMMAND run_functional_tests --gtest_filter=ModuleLifecycleTest.*)
    add_test(NAME Functional.SolverConvergence COMMAND run_functional_tests --gtest_filter=SolverConvergenceTest.*)
    add_test(NAME Functional.WellOperations   COMMAND run_functional_tests --gtest_filter=WellOperationsTest.*)

    # --- Physics validation (MMS) tests ---
    add_test(NAME Physics.MMS.Diffusion       COMMAND run_physics_tests --gtest_filter=MMSDiffusionTest.*)
    add_test(NAME Physics.MMS.Elasticity      COMMAND run_physics_tests --gtest_filter=MMSElasticityTest.*)
    add_test(NAME Physics.MMS.WavePropagation COMMAND run_physics_tests --gtest_filter=MMSWavePropagationTest.*)
    add_test(NAME Physics.AnalyticalSolutions COMMAND run_physics_tests --gtest_filter=AnalyticalSolutionsTest.*)
    add_test(NAME Physics.Thermoporoelastic   COMMAND run_physics_tests --gtest_filter=ThermoporoelasticConsolidationTest.*)

    # --- Integration tests ---
    add_test(NAME Integration.FullSimulation  COMMAND run_integration_tests --gtest_filter=FullSimulationTest.*)
    add_test(NAME Integration.Restart         COMMAND run_integration_tests --gtest_filter=RestartTest.*)
    add_test(NAME Integration.CoupledPhysics  COMMAND run_integration_tests --gtest_filter=CoupledPhysicsTest.*)

    # --- Performance tests ---
    add_test(NAME Performance.Benchmarks      COMMAND run_performance_tests --gtest_filter=BenchmarkTest.*)
    add_test(NAME Performance.Scaling         COMMAND run_performance_tests --gtest_filter=ScalingTest.*)
    add_test(NAME Performance.Memory          COMMAND run_performance_tests --gtest_filter=MemoryStabilityTest.*)

    # --- Experimental tests ---
    add_test(NAME Experimental.NeuralStubs    COMMAND run_experimental_tests --gtest_filter=ExperimentalNeuralTest.*)

    # =========================================================================
    # Test Properties (Timeouts and Labels)
    # =========================================================================

    # Unit tests — fast
    set_tests_properties(
        Unit.ConfigReader Unit.PhysicsModuleRegistry Unit.SimulationConfig
        Unit.PhysicsKernelBase Unit.SinglePhaseFlow Unit.Geomechanics
        Unit.Thermal Unit.Elastodynamics Unit.Poroelastodynamics
        Unit.BoundaryConditions Unit.WellModel Unit.FractureModel
        Unit.MaterialModel Unit.SeismicSource Unit.ExplosionSource
        Unit.EclipseIO Unit.MLRegistry Unit.UnitSystem Unit.CoordinateSystem
        Unit.CO2Properties Unit.FlashCalculation Unit.ResidualTrapping
        PROPERTIES TIMEOUT 120 LABELS "unit"
    )

    # Functional tests
    set_tests_properties(
        Functional.SimulatorInit Functional.ModuleLifecycle
        Functional.SolverConvergence Functional.WellOperations
        PROPERTIES TIMEOUT 180 LABELS "functional"
    )

    # Physics validation tests
    set_tests_properties(
        Physics.MMS.Diffusion Physics.MMS.Elasticity
        Physics.MMS.WavePropagation Physics.AnalyticalSolutions
        Physics.Thermoporoelastic
        PROPERTIES TIMEOUT 300 LABELS "physics_validation"
    )

    # Integration tests
    set_tests_properties(
        Integration.FullSimulation Integration.Restart
        Integration.CoupledPhysics
        PROPERTIES TIMEOUT 600 LABELS "integration"
    )

    # Performance tests
    set_tests_properties(
        Performance.Benchmarks Performance.Scaling Performance.Memory
        PROPERTIES TIMEOUT 900 LABELS "performance"
    )

    # Experimental tests
    set_tests_properties(
        Experimental.NeuralStubs
        PROPERTIES TIMEOUT 120 LABELS "experimental"
    )

else()
    message(WARNING "GTest not found. Skipping test suite.")
endif()

# =========================================================================
# Custom Make Targets
# =========================================================================

if(GTest_FOUND)
    add_custom_target(test_unit
        COMMAND ${CMAKE_CTEST_COMMAND} -L "unit" --output-on-failure
        DEPENDS run_unit_tests
        COMMENT "Running unit tests")

    add_custom_target(test_functional
        COMMAND ${CMAKE_CTEST_COMMAND} -L "functional" --output-on-failure
        DEPENDS run_functional_tests
        COMMENT "Running functional tests")

    add_custom_target(test_physics
        COMMAND ${CMAKE_CTEST_COMMAND} -L "physics_validation" --output-on-failure
        DEPENDS run_physics_tests
        COMMENT "Running physics validation tests")

    add_custom_target(test_integration
        COMMAND ${CMAKE_CTEST_COMMAND} -L "integration" --output-on-failure
        DEPENDS run_integration_tests
        COMMENT "Running integration tests")

    add_custom_target(test_performance
        COMMAND ${CMAKE_CTEST_COMMAND} -L "performance" --output-on-failure
        DEPENDS run_performance_tests
        COMMENT "Running performance tests")

    add_custom_target(test_all
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS run_all_tests
        COMMENT "Running all tests")
endif()

# =========================================================================
# Copy Test Data
# =========================================================================
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
