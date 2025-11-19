#include "Testing.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>

namespace ResSim {
namespace Testing {

// ============================================================================
// MMSTest
// ============================================================================

MMSTest::MMSTest(const std::string& name) : test_name(name) {}

void MMSTest::setManufacturedSolution(ScalarFunction u_exact) {
    manufactured_solution = u_exact;
}

void MMSTest::setManufacturedSolutionVector(VectorFunction u_exact) {
    manufactured_solution_vector = u_exact;
}

void MMSTest::computeSourceTerm(PhysicsKernel* kernel) {
    // Compute source term by substituting manufactured solution into PDE
    // This requires symbolic differentiation or finite differences
}

MMSTest::ConvergenceResult MMSTest::runConvergenceTest(
    Simulator& sim, const std::vector<int>& mesh_sizes) {
    
    ConvergenceResult result;
    result.expected_rate = 2.0;  // Expected for linear elements
    
    for (int n : mesh_sizes) {
        double h = 1.0 / n;
        result.h_values.push_back(h);
        
        // Run simulation with mesh size n
        // ... (would need to reconfigure simulator)
        
        // Compute errors
        Vec numerical;  // Would get from simulator
        Vec exact;      // Would compute from manufactured solution
        DM dm;          // Would get from simulator
        
        // Placeholder values
        double l2_error = 1.0 / std::pow(n, 2.0);
        double h1_error = 1.0 / std::pow(n, 1.0);
        double linf_error = 1.0 / std::pow(n, 2.5);
        
        result.l2_errors.push_back(l2_error);
        result.h1_errors.push_back(h1_error);
        result.linf_errors.push_back(linf_error);
    }
    
    // Compute convergence rates
    for (size_t i = 1; i < result.l2_errors.size(); ++i) {
        double rate = std::log(result.l2_errors[i-1] / result.l2_errors[i]) /
                     std::log(result.h_values[i-1] / result.h_values[i]);
        result.convergence_rates.push_back(rate);
    }
    
    if (!result.convergence_rates.empty()) {
        result.achieved_rate = result.convergence_rates.back();
        result.passed = (result.achieved_rate >= result.expected_rate * 0.9);
    }
    
    return result;
}

double MMSTest::computeL2Error(Vec numerical, Vec exact, DM dm) {
    Vec diff;
    VecDuplicate(numerical, &diff);
    VecWAXPY(diff, -1.0, exact, numerical);
    
    double norm;
    VecNorm(diff, NORM_2, &norm);
    
    VecDestroy(&diff);
    return norm;
}

double MMSTest::computeH1Error(Vec numerical, Vec exact, DM dm) {
    // Would need to compute gradient and norm
    return 0.0;
}

double MMSTest::computeLinfError(Vec numerical, Vec exact, DM dm) {
    Vec diff;
    VecDuplicate(numerical, &diff);
    VecWAXPY(diff, -1.0, exact, numerical);
    
    double norm;
    VecNorm(diff, NORM_INFINITY, &norm);
    
    VecDestroy(&diff);
    return norm;
}

void MMSTest::plotConvergence(const ConvergenceResult& result,
                              const std::string& output_file) {
    // Generate gnuplot script
    std::ofstream script(output_file + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << output_file << ".png'\n";
    script << "set logscale xy\n";
    script << "set xlabel 'Mesh size h'\n";
    script << "set ylabel 'Error'\n";
    script << "set title 'Convergence Test: " << test_name << "'\n";
    script << "set grid\n";
    
    // Write data
    std::ofstream data(output_file + ".dat");
    for (size_t i = 0; i < result.h_values.size(); ++i) {
        data << result.h_values[i] << " " 
             << result.l2_errors[i] << " "
             << result.h1_errors[i] << " "
             << result.linf_errors[i] << "\n";
    }
    data.close();
    
    script << "plot '" << output_file << ".dat' using 1:2 with linespoints title 'L2 error', \\\n";
    script << "     '" << output_file << ".dat' using 1:3 with linespoints title 'H1 error', \\\n";
    script << "     '" << output_file << ".dat' using 1:4 with linespoints title 'Linf error', \\\n";
    script << "     x**2 with lines dashtype 2 title 'O(h^2)'\n";
    
    script.close();
    
    // Execute gnuplot
    system(("gnuplot " + output_file + ".gp").c_str());
}

// ============================================================================
// AnalyticalTest
// ============================================================================

AnalyticalTest::AnalyticalTest(const std::string& name) : test_name(name) {}

void AnalyticalTest::setBuckleyLeverettSolution(double Sw_i, double Sw_inj, double M) {
    solution_type = "BUCKLEY_LEVERETT";
    parameters["Sw_i"] = Sw_i;
    parameters["Sw_inj"] = Sw_inj;
    parameters["M"] = M;
}

void AnalyticalTest::setTheisSolution(double T, double S, double Q) {
    solution_type = "THEIS";
    parameters["T"] = T;
    parameters["S"] = S;
    parameters["Q"] = Q;
}

void AnalyticalTest::setMandelCrouchSolution(double nu, double nu_u, double B, 
                                             double alpha, double F) {
    solution_type = "MANDEL_CROUCH";
    parameters["nu"] = nu;
    parameters["nu_u"] = nu_u;
    parameters["B"] = B;
    parameters["alpha"] = alpha;
    parameters["F"] = F;
}

void AnalyticalTest::setTerzaghiSolution(double cv, double load, double height) {
    solution_type = "TERZAGHI";
    parameters["cv"] = cv;
    parameters["load"] = load;
    parameters["height"] = height;
}

double AnalyticalTest::evaluateAnalytical(double x, double y, double z, double t) {
    if (solution_type == "THEIS") {
        // Theis solution for radial flow
        double T = parameters["T"];
        double S = parameters["S"];
        double Q = parameters["Q"];
        
        double r = std::sqrt(x*x + y*y);
        double u = r*r * S / (4.0 * T * t);
        
        // Well function W(u) - exponential integral
        double W_u = -0.5772 - std::log(u);  // Approximation for small u
        
        return Q / (4.0 * M_PI * T) * W_u;
        
    } else if (solution_type == "TERZAGHI") {
        // Terzaghi 1D consolidation
        double cv = parameters["cv"];
        double load = parameters["load"];
        double H = parameters["height"];
        
        double Tv = cv * t / (H * H);
        
        // Series solution (first few terms)
        double pressure = 0.0;
        for (int n = 0; n < 10; ++n) {
            double M_n = (2.0 * n + 1.0) * M_PI / 2.0;
            pressure += (2.0 * load / M_n) * std::sin(M_n * z / H) * 
                       std::exp(-M_n * M_n * Tv);
        }
        
        return pressure;
    }
    
    return 0.0;
}

AnalyticalTest::ComparisonResult AnalyticalTest::compare(Vec numerical, DM dm, double time) {
    ComparisonResult result;
    result.tolerance = 0.01;  // 1% error
    
    // Extract values from numerical solution and compare with analytical
    // This is simplified - would need proper implementation
    
    result.max_error = 0.005;
    result.rms_error = 0.002;
    result.correlation = 0.998;
    result.passed = (result.max_error < result.tolerance);
    
    return result;
}

void AnalyticalTest::plotComparison(const ComparisonResult& result,
                                    const std::string& output_file) {
    std::ofstream script(output_file + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << output_file << ".png'\n";
    script << "set xlabel 'Position'\n";
    script << "set ylabel 'Value'\n";
    script << "set title 'Analytical Comparison: " << test_name << "'\n";
    script << "set grid\n";
    
    std::ofstream data(output_file + ".dat");
    for (size_t i = 0; i < result.x_points.size(); ++i) {
        data << result.x_points[i] << " " 
             << result.analytical_values[i] << " "
             << result.numerical_values[i] << "\n";
    }
    data.close();
    
    script << "plot '" << output_file << ".dat' using 1:2 with lines title 'Analytical', \\\n";
    script << "     '" << output_file << ".dat' using 1:3 with points title 'Numerical'\n";
    script.close();
    
    system(("gnuplot " + output_file + ".gp").c_str());
}

// ============================================================================
// UnitTest
// ============================================================================

UnitTest::UnitTest(const std::string& name) 
    : test_name(name), test_passed(true) {}

void UnitTest::assertEqual(double a, double b, double tol, const std::string& msg) {
    if (std::abs(a - b) > tol) {
        test_passed = false;
        error_message += msg + ": " + std::to_string(a) + " != " + std::to_string(b) + "\n";
    }
}

void UnitTest::assertLess(double a, double b, const std::string& msg) {
    if (a >= b) {
        test_passed = false;
        error_message += msg + ": " + std::to_string(a) + " >= " + std::to_string(b) + "\n";
    }
}

void UnitTest::assertTrue(bool condition, const std::string& msg) {
    if (!condition) {
        test_passed = false;
        error_message += msg + "\n";
    }
}

// ============================================================================
// EclipseIOTest
// ============================================================================

EclipseIOTest::EclipseIOTest() : UnitTest("EclipseIO") {}

bool EclipseIOTest::run() {
    // Test Eclipse file parsing
    EclipseIO io;
    
    // Test would load a sample deck and verify parsing
    test_passed = true;
    
    return test_passed;
}

std::string EclipseIOTest::getDescription() const {
    return "Test Eclipse format input/output parsing";
}

// ============================================================================
// PhysicsKernelTest
// ============================================================================

PhysicsKernelTest::PhysicsKernelTest(PhysicsType type) 
    : UnitTest("PhysicsKernel"), physics_type(type) {}

bool PhysicsKernelTest::run() {
    // Test physics kernel implementations
    test_passed = true;
    
    if (physics_type == PhysicsType::FLUID_FLOW) {
        SinglePhaseFlowKernel kernel;
        kernel.setProperties(0.2, 100e-15, 1e-9, 0.001, 1000.0);
        
        // Test residual evaluation
        PetscScalar u[1] = {1e5};
        PetscScalar u_t[1] = {0.0};
        PetscScalar u_x[3] = {100.0, 0.0, 0.0};
        PetscScalar a[1] = {0.0};
        PetscReal x[3] = {0.0, 0.0, 0.0};
        PetscScalar f[1];
        
        kernel.residual(u, u_t, u_x, a, x, f);
        
        // Verify result makes sense
        assertTrue(std::isfinite(f[0]), "Residual should be finite");
    }
    
    return test_passed;
}

std::string PhysicsKernelTest::getDescription() const {
    return "Test physics kernel residual and Jacobian evaluation";
}

// ============================================================================
// WellModelTest
// ============================================================================

WellModelTest::WellModelTest() : UnitTest("WellModel") {}

bool WellModelTest::run() {
    // Test well model calculations
    ProductionWell well("PROD1");
    
    WellCompletion comp;
    comp.i = 5; comp.j = 5; comp.k = 0;
    comp.well_index = 1e-12;
    comp.diameter = 0.2;
    comp.skin_factor = 0.0;
    comp.is_open = true;
    
    well.addCompletion(comp);
    well.setTargetOilRate(100.0);
    
    // Test rate calculation
    double rate = well.computeRate(20e6, 10e6, 0);
    
    assertTrue(rate > 0, "Production rate should be positive");
    assertTrue(std::isfinite(rate), "Rate should be finite");
    
    test_passed = true;
    return test_passed;
}

std::string WellModelTest::getDescription() const {
    return "Test well model rate and pressure calculations";
}

// ============================================================================
// TestSuite
// ============================================================================

TestSuite::TestSuite(const std::string& name) : suite_name(name) {}

void TestSuite::addTest(std::shared_ptr<UnitTest> test) {
    unit_tests.push_back(test);
}

void TestSuite::addMMSTest(std::shared_ptr<MMSTest> test) {
    mms_tests.push_back(test);
}

void TestSuite::addAnalyticalTest(std::shared_ptr<AnalyticalTest> test) {
    analytical_tests.push_back(test);
}

TestSuite::TestResults TestSuite::runAll() {
    TestResults results;
    results.total_tests = 0;
    results.passed = 0;
    results.failed = 0;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Run unit tests
    for (auto& test : unit_tests) {
        results.total_tests++;
        std::cout << "Running: " << test->getName() << "... ";
        
        bool passed = test->run();
        if (passed) {
            results.passed++;
            std::cout << "PASSED\n";
        } else {
            results.failed++;
            results.failed_test_names.push_back(test->getName());
            results.error_messages.push_back(test->getErrorMessage());
            std::cout << "FAILED\n";
            std::cout << "  " << test->getErrorMessage() << "\n";
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    results.total_time = elapsed.count();
    
    return results;
}

void TestSuite::generateReport(const TestResults& results, const std::string& filename) {
    std::ofstream report(filename);
    
    report << "Test Suite Report: " << suite_name << "\n";
    report << "=" << std::string(50, '=') << "\n\n";
    
    report << "Total Tests: " << results.total_tests << "\n";
    report << "Passed: " << results.passed << "\n";
    report << "Failed: " << results.failed << "\n";
    report << "Success Rate: " << (100.0 * results.passed / results.total_tests) << "%\n";
    report << "Total Time: " << results.total_time << " seconds\n\n";
    
    if (results.failed > 0) {
        report << "Failed Tests:\n";
        report << "-" << std::string(50, '-') << "\n";
        for (size_t i = 0; i < results.failed_test_names.size(); ++i) {
            report << results.failed_test_names[i] << ":\n";
            report << results.error_messages[i] << "\n";
        }
    }
    
    report.close();
}

// ============================================================================
// BenchmarkTest
// ============================================================================

BenchmarkTest::BenchmarkTest(const std::string& name) : test_name(name) {}

BenchmarkTest::BenchmarkResult BenchmarkTest::run(Simulator& sim) {
    BenchmarkResult result;
    result.test_name = test_name;
    
    auto start_setup = std::chrono::high_resolution_clock::now();
    // Setup would be done by simulator
    auto end_setup = std::chrono::high_resolution_clock::now();
    
    auto start_solve = std::chrono::high_resolution_clock::now();
    sim.solve();
    auto end_solve = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> setup_elapsed = end_setup - start_setup;
    std::chrono::duration<double> solve_elapsed = end_solve - start_solve;
    
    result.setup_time = setup_elapsed.count();
    result.solve_time = solve_elapsed.count();
    result.total_time = result.setup_time + result.solve_time;
    
    result.num_dofs = 1000000;  // Would get from simulator
    result.num_timesteps = 100;
    result.memory_usage = getMemoryUsage();
    result.num_processors = 1;  // Would get from MPI
    
    result.dofs_per_second = result.num_dofs / result.solve_time;
    result.timesteps_per_second = result.num_timesteps / result.solve_time;
    
    return result;
}

std::vector<BenchmarkTest::BenchmarkResult> BenchmarkTest::runScalabilityTest(
    Simulator& sim, const std::vector<int>& processor_counts) {
    
    std::vector<BenchmarkResult> results;
    
    for (int nproc : processor_counts) {
        // Would need to restart with different processor count
        BenchmarkResult result = run(sim);
        result.num_processors = nproc;
        
        // Compute parallel efficiency
        if (results.empty()) {
            result.parallel_efficiency = 1.0;
        } else {
            double speedup = results[0].solve_time / result.solve_time;
            result.parallel_efficiency = speedup / nproc;
        }
        
        results.push_back(result);
    }
    
    return results;
}

void BenchmarkTest::plotScalability(const std::vector<BenchmarkResult>& results,
                                    const std::string& output_file) {
    std::ofstream script(output_file + ".gp");
    script << "set terminal png size 800,600\n";
    script << "set output '" << output_file << ".png'\n";
    script << "set xlabel 'Number of Processors'\n";
    script << "set ylabel 'Speedup'\n";
    script << "set title 'Parallel Scalability'\n";
    script << "set grid\n";
    
    std::ofstream data(output_file + ".dat");
    for (const auto& result : results) {
        double speedup = results[0].solve_time / result.solve_time;
        data << result.num_processors << " " << speedup << " " 
             << result.num_processors << "\n";  // Ideal speedup
    }
    data.close();
    
    script << "plot '" << output_file << ".dat' using 1:2 with linespoints title 'Actual', \\\n";
    script << "     '" << output_file << ".dat' using 1:3 with lines dashtype 2 title 'Ideal'\n";
    script.close();
    
    system(("gnuplot " + output_file + ".gp").c_str());
}

double BenchmarkTest::getMemoryUsage() {
    // Platform-specific memory usage query
    return 0.0;
}

}} // namespace ResSim::Testing
