#ifndef TESTING_HPP
#define TESTING_HPP

#include "ReservoirSim.hpp"
#include "Simulator.hpp"
#include <functional>
#include <vector>
#include <string>

namespace FSRM {
namespace Testing {

// Method of Manufactured Solutions (MMS)
class MMSTest {
public:
    MMSTest(const std::string& name);
    
    // Define manufactured solution
    using ScalarFunction = std::function<double(double, double, double, double)>;
    using VectorFunction = std::function<void(double, double, double, double, double*)>;
    
    void setManufacturedSolution(ScalarFunction u_exact);
    void setManufacturedSolutionVector(VectorFunction u_exact);
    
    // Compute source term from manufactured solution
    void computeSourceTerm(PhysicsKernel* kernel);
    
    // Run convergence test
    struct ConvergenceResult {
        std::vector<double> h_values;      // mesh sizes
        std::vector<double> l2_errors;     // L2 errors
        std::vector<double> h1_errors;     // H1 errors
        std::vector<double> linf_errors;   // L-infinity errors
        std::vector<double> convergence_rates;
        bool passed;
        double expected_rate;
        double achieved_rate;
    };
    
    ConvergenceResult runConvergenceTest(Simulator& sim,
                                        const std::vector<int>& mesh_sizes);
    
    // Compute errors
    double computeL2Error(Vec numerical, Vec exact, DM dm);
    double computeH1Error(Vec numerical, Vec exact, DM dm);
    double computeLinfError(Vec numerical, Vec exact, DM dm);
    
    // Visualization
    void plotConvergence(const ConvergenceResult& result,
                        const std::string& output_file);
    
private:
    std::string test_name;
    ScalarFunction manufactured_solution;
    VectorFunction manufactured_solution_vector;
    ScalarFunction source_term;
};

// Analytical solution comparison
class AnalyticalTest {
public:
    AnalyticalTest(const std::string& name);
    
    // Common analytical solutions
    void setBuckleyLeverettSolution(double Sw_i, double Sw_inj,
                                   double mobility_ratio);
    void setTheisSolution(double transmissivity, double storativity,
                         double pumping_rate);
    void setMandelCrouchSolution(double nu, double nu_u, double B,
                                double alpha, double F);
    void setTerzaghiSolution(double cv, double load, double height);
    
    // Compare with numerical solution
    struct ComparisonResult {
        double max_error;
        double rms_error;
        double correlation;
        std::vector<double> x_points;
        std::vector<double> analytical_values;
        std::vector<double> numerical_values;
        bool passed;
        double tolerance;
    };
    
    ComparisonResult compare(Vec numerical, DM dm, double time);
    
    // Visualization
    void plotComparison(const ComparisonResult& result,
                       const std::string& output_file);
    
private:
    std::string test_name;
    std::string solution_type;
    std::map<std::string, double> parameters;
    
    double evaluateAnalytical(double x, double y, double z, double t);
};

// Unit test framework
class UnitTest {
public:
    UnitTest(const std::string& name);
    virtual ~UnitTest() = default;
    
    virtual bool run() = 0;
    virtual std::string getDescription() const = 0;
    
    std::string getName() const { return test_name; }
    bool passed() const { return test_passed; }
    std::string getErrorMessage() const { return error_message; }
    
protected:
    std::string test_name;
    bool test_passed;
    std::string error_message;
    
    // Assertion helpers
    void assertEqual(double a, double b, double tol, const std::string& msg);
    void assertLess(double a, double b, const std::string& msg);
    void assertTrue(bool condition, const std::string& msg);
};

// Specific unit tests
class EclipseIOTest : public UnitTest {
public:
    EclipseIOTest();
    bool run() override;
    std::string getDescription() const override;
};

class PhysicsKernelTest : public UnitTest {
public:
    PhysicsKernelTest(PhysicsType type);
    bool run() override;
    std::string getDescription() const override;
    
private:
    PhysicsType physics_type;
};

class WellModelTest : public UnitTest {
public:
    WellModelTest();
    bool run() override;
    std::string getDescription() const override;
};

// Test suite manager
class TestSuite {
public:
    TestSuite(const std::string& name);
    
    void addTest(std::shared_ptr<UnitTest> test);
    void addMMSTest(std::shared_ptr<MMSTest> test);
    void addAnalyticalTest(std::shared_ptr<AnalyticalTest> test);
    
    struct TestResults {
        int total_tests;
        int passed;
        int failed;
        std::vector<std::string> failed_test_names;
        std::vector<std::string> error_messages;
        double total_time;
        
        // Detailed per-test information
        std::vector<std::string> test_names;
        std::vector<std::string> test_descriptions;
        std::vector<double> test_times;
        std::vector<bool> test_passed;
    };
    
    TestResults runAll();
    void generateReport(const TestResults& results, const std::string& filename);
    
private:
    std::string suite_name;
    std::vector<std::shared_ptr<UnitTest>> unit_tests;
    std::vector<std::shared_ptr<MMSTest>> mms_tests;
    std::vector<std::shared_ptr<AnalyticalTest>> analytical_tests;
};

// Benchmark tests
class BenchmarkTest {
public:
    BenchmarkTest(const std::string& name);
    
    struct BenchmarkResult {
        std::string test_name;
        double setup_time;
        double solve_time;
        double total_time;
        int num_dofs;
        int num_timesteps;
        double memory_usage;
        int num_processors;
        
        // Performance metrics
        double dofs_per_second;
        double timesteps_per_second;
        double parallel_efficiency;
    };
    
    BenchmarkResult run(Simulator& sim);
    
    // Scalability tests
    std::vector<BenchmarkResult> runScalabilityTest(Simulator& sim,
                                                    const std::vector<int>& processor_counts);
    
    void plotScalability(const std::vector<BenchmarkResult>& results,
                        const std::string& output_file);
    
private:
    std::string test_name;
    double getMemoryUsage();
};

}} // namespace FSRM::Testing

#endif // TESTING_HPP
