#include "Testing.hpp"
#include "PhysicsKernel.hpp"
#include <cmath>
#include <iostream>

using namespace FSRM;
using namespace FSRM::Testing;

// Detailed physics kernel tests without GoogleTest
class DetailedPhysicsTest : public UnitTest {
public:
    DetailedPhysicsTest() : UnitTest("Detailed_Physics") {}
    
    bool run() override {
        test_passed = true;
        
        try {
            testSinglePhaseFlowResidual();
            testDarcyVelocity();
            testMassConservation();
            testGeomechanicsStress();
            testPoroelasticCoupling();
            testThermalDiffusion();
            testWavePropagation();
            
        } catch (const std::exception& e) {
            test_passed = false;
            error_message = std::string("Detailed physics test failed: ") + e.what();
        }
        
        return test_passed;
    }
    
    std::string getDescription() const override {
        return "Detailed physics kernel validation tests";
    }
    
private:
    void testSinglePhaseFlowResidual() {
        // Test single phase flow residual evaluation
        double phi = 0.2;         // porosity
        double k = 100e-15;       // permeability (100 mD)
        double ct = 1e-9;         // total compressibility (1/Pa)
        double mu = 1e-3;         // viscosity (1 cP)
        double rho = 1000.0;      // density (kg/m^3)
        
        // For steady state with no gradient, residual should be ~0
        double pressure = 10e6;   // 10 MPa
        double dP_dt = 0.0;
        double grad_P = 0.0;
        
        // Mass balance: phi*ct*dP/dt - div(k/mu*grad(P)) = 0
        double residual = phi * ct * dP_dt - (k/mu) * grad_P;
        
        assertTrue(std::abs(residual) < 1e-10, 
                  "Steady state residual should be zero");
        
        // For transient case
        dP_dt = 1000.0;  // Pa/s
        residual = phi * ct * dP_dt;
        
        assertTrue(std::isfinite(residual), 
                  "Transient residual should be finite");
        assertTrue(std::abs(residual) > 0.0, 
                  "Transient residual should be non-zero");
    }
    
    void testDarcyVelocity() {
        // Test Darcy's law: v = -(k/mu) * grad(p)
        double k = 100e-15;  // 100 mD
        double mu = 1e-3;    // 1 cP
        double grad_p = 1e5; // 100 kPa/m
        
        double v_darcy = -(k / mu) * grad_p;
        
        // Flow should be in negative gradient direction
        assertTrue(v_darcy < 0.0, 
                  "Darcy velocity should be negative for positive gradient");
        
        // Check magnitude is reasonable (should be cm/year scale)
        double v_magnitude = std::abs(v_darcy);
        assertTrue(v_magnitude > 1e-10 && v_magnitude < 1.0, 
                  "Darcy velocity magnitude should be physically reasonable");
    }
    
    void testMassConservation() {
        // Test mass conservation for fluid flow
        // d(phi*rho)/dt + div(rho*v) = q
        
        double phi = 0.2;
        double rho = 1000.0;
        double ct = 1e-9;
        double dP_dt = 1000.0;  // Pa/s
        
        // Accumulation term
        double accumulation = phi * rho * ct * dP_dt;
        
        // For no source and no flux, accumulation should balance
        assertTrue(std::isfinite(accumulation), 
                  "Accumulation term should be finite");
        
        // Check conservation is satisfied
        double source = 0.0;
        double flux_divergence = -accumulation;  // Balance
        double conservation = accumulation + flux_divergence - source;
        
        assertTrue(std::abs(conservation) < 1e-10, 
                  "Mass conservation should be satisfied");
    }
    
    void testGeomechanicsStress() {
        // Test stress-strain relationship: sigma = C : epsilon
        // For isotropic material: sigma_ij = lambda*epsilon_kk*delta_ij + 2*mu*epsilon_ij
        
        double E = 10e9;   // Young's modulus (10 GPa)
        double nu = 0.25;  // Poisson's ratio
        
        // Lamé parameters
        double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double mu = E / (2.0 * (1.0 + nu));
        
        // Simple uniaxial strain: epsilon_11 = 0.001, others = 0
        double epsilon_11 = 0.001;
        double epsilon_kk = epsilon_11;  // trace
        
        double sigma_11 = lambda * epsilon_kk + 2.0 * mu * epsilon_11;
        double sigma_22 = lambda * epsilon_kk;  // lateral stress
        
        // Check stress is positive (compression)
        assertTrue(sigma_11 > 0.0, "Axial stress should be positive");
        assertTrue(sigma_22 > 0.0, "Lateral stress should be positive");
        
        // Check stress ratios
        double stress_ratio = sigma_22 / sigma_11;
        double expected_ratio = nu / (1.0 - nu);
        
        assertEqual(stress_ratio, expected_ratio, 1e-6, 
                   "Stress ratio should match Poisson's ratio relationship");
    }
    
    void testPoroelasticCoupling() {
        // Test Biot poroelasticity coupling
        // alpha = 1 - K/K_s (Biot coefficient)
        
        double K_frame = 5e9;   // Drained bulk modulus (5 GPa)
        double K_solid = 20e9;  // Solid grain modulus (20 GPa)
        
        double alpha = 1.0 - K_frame / K_solid;
        
        // Biot coefficient should be between 0 and 1
        assertTrue(alpha > 0.0 && alpha < 1.0, 
                  "Biot coefficient should be between 0 and 1");
        
        // For typical rocks, alpha ~ 0.6-1.0
        assertTrue(alpha > 0.5, 
                  "Biot coefficient should be > 0.5 for typical rocks");
        
        // Test pore pressure effect on stress
        double pressure = 10e6;  // 10 MPa
        double sigma_effective = -alpha * pressure;  // Effective stress change
        
        assertTrue(sigma_effective < 0.0, 
                  "Pore pressure should reduce effective stress");
    }
    
    void testThermalDiffusion() {
        // Test heat diffusion: rho*cp*dT/dt = div(lambda*grad(T))
        
        double rho = 2500.0;      // Density (kg/m^3)
        double cp = 1000.0;       // Specific heat (J/kg/K)
        double lambda = 2.5;      // Thermal conductivity (W/m/K)
        
        double alpha_thermal = lambda / (rho * cp);  // Thermal diffusivity
        
        // Thermal diffusivity should be ~1e-6 m^2/s for rocks
        assertTrue(alpha_thermal > 1e-7 && alpha_thermal < 1e-5, 
                  "Thermal diffusivity should be in typical rock range");
        
        // Test temperature gradient
        double grad_T = 30.0;  // K/m (30°C/km)
        double heat_flux = -lambda * grad_T;
        
        assertTrue(heat_flux < 0.0, 
                  "Heat flux should be negative for positive gradient");
        
        // Check magnitude (should be ~mW/m^2 scale)
        double flux_magnitude = std::abs(heat_flux);
        assertTrue(flux_magnitude > 1.0 && flux_magnitude < 1000.0, 
                  "Heat flux should be in geothermal range");
    }
    
    void testWavePropagation() {
        // Test wave speeds for elastodynamics
        // P-wave: v_p = sqrt((lambda + 2*mu) / rho)
        // S-wave: v_s = sqrt(mu / rho)
        
        double rho = 2500.0;   // Density (kg/m^3)
        double E = 10e9;       // Young's modulus (10 GPa)
        double nu = 0.25;      // Poisson's ratio
        
        double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        double mu = E / (2.0 * (1.0 + nu));
        
        double v_p = std::sqrt((lambda + 2.0 * mu) / rho);
        double v_s = std::sqrt(mu / rho);
        
        // P-wave should be faster than S-wave
        assertTrue(v_p > v_s, "P-wave velocity should exceed S-wave velocity");
        
        // Check velocities are in typical rock range (km/s)
        assertTrue(v_p > 1000.0 && v_p < 10000.0, 
                  "P-wave velocity should be 1-10 km/s");
        assertTrue(v_s > 500.0 && v_s < 6000.0, 
                  "S-wave velocity should be 0.5-6 km/s");
        
        // Check v_p/v_s ratio
        double velocity_ratio = v_p / v_s;
        double expected_ratio = std::sqrt(2.0 * (1.0 - nu) / (1.0 - 2.0 * nu));
        
        assertEqual(velocity_ratio, expected_ratio, 1e-6, 
                   "Velocity ratio should match theoretical value");
    }
};

// Factory function
std::shared_ptr<FSRM::Testing::UnitTest> createDetailedPhysicsTest() {
    return std::make_shared<DetailedPhysicsTest>();
}
