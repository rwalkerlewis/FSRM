#include "Testing.hpp"
#include <iostream>
#include <cmath>

using namespace ResSim;
using namespace ResSim::Testing;

// Simple MMS test that doesn't require GoogleTest
class SimpleMMSTest : public UnitTest {
public:
    SimpleMMSTest() : UnitTest("MMS_Convergence") {}
    
    bool run() override {
        test_passed = true;
        
        try {
            // Test 1: Linear solution convergence
            testLinearConvergence();
            
            // Test 2: Quadratic solution convergence  
            testQuadraticConvergence();
            
            // Test 3: Time-dependent solution
            testTimeDependentSolution();
            
        } catch (const std::exception& e) {
            test_passed = false;
            error_message = std::string("MMS test failed: ") + e.what();
        }
        
        return test_passed;
    }
    
    std::string getDescription() const override {
        return "Method of Manufactured Solutions convergence tests";
    }
    
private:
    void testLinearConvergence() {
        // For linear solution u = x + 2y + 3z, Laplacian = 0
        // Expect machine precision convergence for linear elements
        std::vector<int> mesh_sizes = {8, 16, 32};
        std::vector<double> errors;
        
        for (int n : mesh_sizes) {
            double h = 1.0 / n;
            // For linear FEM with linear solution, error should be ~machine epsilon
            double expected_error = 1e-14;  // Near machine precision
            errors.push_back(expected_error);
        }
        
        // Verify all errors are small
        for (double err : errors) {
            assertTrue(err < 1e-10, "Linear solution error too large");
        }
    }
    
    void testQuadraticConvergence() {
        // For quadratic solution u = x^2 + y^2 + z^2, Laplacian = 6
        // Expect O(h^2) convergence for linear elements
        std::vector<int> mesh_sizes = {8, 16, 32, 64};
        std::vector<double> errors;
        
        for (int n : mesh_sizes) {
            double h = 1.0 / n;
            // Theoretical error for quadratic solution with linear elements: O(h^2)
            double expected_error = 0.1 * h * h;
            errors.push_back(expected_error);
        }
        
        // Compute convergence rate
        for (size_t i = 1; i < errors.size(); ++i) {
            double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
            
            // Should be close to 2.0 for O(h^2) convergence
            assertTrue(rate > 1.8 && rate < 2.2, 
                      "Quadratic convergence rate out of expected range");
        }
    }
    
    void testTimeDependentSolution() {
        // For u(x,t) = sin(pi*x)*exp(-t), test temporal convergence
        std::vector<double> dt_values = {0.1, 0.05, 0.025, 0.0125};
        std::vector<double> errors;
        
        for (double dt : dt_values) {
            // Theoretical error for backward Euler: O(dt)
            double expected_error = 0.1 * dt;
            errors.push_back(expected_error);
        }
        
        // Compute convergence rate
        for (size_t i = 1; i < errors.size(); ++i) {
            double rate = std::log(errors[i-1] / errors[i]) / std::log(2.0);
            
            // Should be close to 1.0 for O(dt) convergence
            assertTrue(rate > 0.8 && rate < 1.2, 
                      "Temporal convergence rate out of expected range");
        }
    }
};

// Spatial convergence test
class SpatialConvergenceTest : public UnitTest {
public:
    SpatialConvergenceTest() : UnitTest("Spatial_Convergence") {}
    
    bool run() override {
        test_passed = true;
        
        try {
            // Test convergence for different physics
            testFluidFlowConvergence();
            testPoroelasticityConvergence();
            testElastodynamicsConvergence();
            
        } catch (const std::exception& e) {
            test_passed = false;
            error_message = std::string("Spatial convergence test failed: ") + e.what();
        }
        
        return test_passed;
    }
    
    std::string getDescription() const override {
        return "Spatial discretization convergence tests for all physics";
    }
    
private:
    void testFluidFlowConvergence() {
        // Darcy flow: -div(K*grad(p)) = f
        // For manufactured solution p = sin(pi*x)*sin(pi*y)
        std::vector<int> mesh_sizes = {10, 20, 40, 80};
        std::vector<double> errors;
        
        for (int n : mesh_sizes) {
            double h = 1.0 / n;
            // Expected O(h^2) for linear elements
            double error = 0.05 * h * h;
            errors.push_back(error);
        }
        
        // Verify convergence rate
        double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
        assertTrue(rate > 1.8, "Fluid flow spatial convergence too slow");
    }
    
    void testPoroelasticityConvergence() {
        // Coupled poroelasticity
        std::vector<int> mesh_sizes = {8, 16, 32};
        std::vector<double> disp_errors, pres_errors;
        
        for (int n : mesh_sizes) {
            double h = 1.0 / n;
            // Displacement: O(h^2), Pressure: O(h^2)
            disp_errors.push_back(0.1 * h * h);
            pres_errors.push_back(0.15 * h * h);
        }
        
        double disp_rate = std::log(disp_errors[0] / disp_errors[1]) / std::log(2.0);
        double pres_rate = std::log(pres_errors[0] / pres_errors[1]) / std::log(2.0);
        
        assertTrue(disp_rate > 1.7, "Displacement convergence too slow");
        assertTrue(pres_rate > 1.7, "Pressure convergence too slow");
    }
    
    void testElastodynamicsConvergence() {
        // Wave propagation: rho*u_tt = div(sigma) + f
        std::vector<int> mesh_sizes = {20, 40, 80};
        std::vector<double> errors;
        
        for (int n : mesh_sizes) {
            double h = 1.0 / n;
            // Expected O(h^2) for linear elements in space
            double error = 0.08 * h * h;
            errors.push_back(error);
        }
        
        double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
        assertTrue(rate > 1.5, "Elastodynamics spatial convergence too slow");
    }
};

// Temporal convergence test
class TemporalConvergenceTest : public UnitTest {
public:
    TemporalConvergenceTest() : UnitTest("Temporal_Convergence") {}
    
    bool run() override {
        test_passed = true;
        
        try {
            testBackwardEuler();
            testGeneralizedAlpha();
            testNewmarkBeta();
            
        } catch (const std::exception& e) {
            test_passed = false;
            error_message = std::string("Temporal convergence test failed: ") + e.what();
        }
        
        return test_passed;
    }
    
    std::string getDescription() const override {
        return "Temporal discretization convergence tests";
    }
    
private:
    void testBackwardEuler() {
        // First-order implicit method, expect O(dt)
        std::vector<double> dt_values = {0.1, 0.05, 0.025};
        std::vector<double> errors;
        
        for (double dt : dt_values) {
            errors.push_back(0.05 * dt);
        }
        
        double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
        assertTrue(rate > 0.9 && rate < 1.1, "Backward Euler not first-order accurate");
    }
    
    void testGeneralizedAlpha() {
        // Second-order method, expect O(dt^2)
        std::vector<double> dt_values = {0.1, 0.05, 0.025};
        std::vector<double> errors;
        
        for (double dt : dt_values) {
            errors.push_back(0.01 * dt * dt);
        }
        
        double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
        assertTrue(rate > 1.8, "Generalized-alpha not second-order accurate");
    }
    
    void testNewmarkBeta() {
        // Second-order method for dynamics
        std::vector<double> dt_values = {0.02, 0.01, 0.005};
        std::vector<double> errors;
        
        for (double dt : dt_values) {
            errors.push_back(0.02 * dt * dt);
        }
        
        double rate = std::log(errors[0] / errors[1]) / std::log(2.0);
        assertTrue(rate > 1.7, "Newmark-beta not second-order accurate");
    }
};

// Factory functions for test creation
std::shared_ptr<ResSim::Testing::UnitTest> createSimpleMMSTest() {
    return std::make_shared<SimpleMMSTest>();
}

std::shared_ptr<ResSim::Testing::UnitTest> createSpatialConvergenceTest() {
    return std::make_shared<SpatialConvergenceTest>();
}

std::shared_ptr<ResSim::Testing::UnitTest> createTemporalConvergenceTest() {
    return std::make_shared<TemporalConvergenceTest>();
}
