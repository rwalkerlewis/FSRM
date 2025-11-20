#include "Testing.hpp"
#include "FractureModel.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace ResSim;
using namespace ResSim::Testing;

// Test fixture for fracture model tests
class FractureModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize MPI/PETSc if needed
    }
    
    void TearDown() override {
        // Cleanup
    }
};

// Test LEFM stress intensity factor calculation
TEST_F(FractureModelTest, StressIntensityFactor) {
    LEFMFractureModel fracture;
    
    // Mode I loading on a penny-shaped crack
    double radius = 10.0;  // meters
    double pressure = 1e6;  // Pa
    
    double K_I = fracture.computeStressIntensityFactor(
        FractureMode::MODE_I, radius, pressure);
    
    // Analytical solution: K_I = (2/π) * σ * sqrt(π*a)
    double K_I_analytical = (2.0 / M_PI) * pressure * std::sqrt(M_PI * radius);
    
    EXPECT_NEAR(K_I, K_I_analytical, 1e3) 
        << "Mode I stress intensity factor should match analytical solution";
}

// Test fracture propagation criterion
TEST_F(FractureModelTest, PropagationCriterion) {
    LEFMFractureModel fracture;
    
    // Set material properties
    double K_IC = 1e6;  // Fracture toughness (Pa·m^0.5)
    fracture.setFractureToughness(K_IC);
    
    // Test sub-critical loading
    double K_I_subcritical = 0.8 * K_IC;
    EXPECT_FALSE(fracture.willPropagate(K_I_subcritical, 0.0, 0.0))
        << "Fracture should not propagate below critical stress intensity";
    
    // Test critical loading
    double K_I_critical = 1.1 * K_IC;
    EXPECT_TRUE(fracture.willPropagate(K_I_critical, 0.0, 0.0))
        << "Fracture should propagate above critical stress intensity";
}

// Test fracture growth rate
TEST_F(FractureModelTest, GrowthRate) {
    LEFMFractureModel fracture;
    
    double K_IC = 1e6;
    fracture.setFractureToughness(K_IC);
    
    // Paris law parameters
    double C = 1e-10;  // Material constant
    double m = 3.0;    // Paris exponent
    fracture.setParisLawParameters(C, m);
    
    double K_I = 1.5 * K_IC;
    double deltaK = K_I - 0.5 * K_IC;  // Stress intensity range
    
    double growth_rate = fracture.computeGrowthRate(deltaK);
    
    // Paris law: da/dN = C * (ΔK)^m
    double expected_rate = C * std::pow(deltaK, m);
    
    EXPECT_NEAR(growth_rate, expected_rate, 1e-13)
        << "Growth rate should follow Paris law";
}

// Test hydraulic fracture width calculation
TEST_F(FractureModelTest, HydraulicFractureWidth) {
    LEFMFractureModel fracture;
    
    // Set mechanical properties
    double E = 20e9;     // Young's modulus (Pa)
    double nu = 0.25;    // Poisson's ratio
    fracture.setMechanicalProperties(E, nu);
    
    double pressure = 5e6;  // Pa
    double length = 100.0;   // m
    
    double width = fracture.computeFractureWidth(pressure, length);
    
    // Plane strain modulus
    double E_prime = E / (1.0 - nu * nu);
    
    // Analytical solution for elliptical crack
    double width_analytical = 4.0 * (1.0 - nu) * pressure * length / E;
    
    EXPECT_NEAR(width, width_analytical, 1.0)
        << "Fracture width should match analytical solution";
}

// Test fracture conductivity
TEST_F(FractureModelTest, FractureConductivity) {
    LEFMFractureModel fracture;
    
    double width = 0.001;      // m (1 mm)
    double permeability = 1e-9;  // m^2 (proppant pack)
    
    fracture.setProppantPermeability(permeability);
    
    double conductivity = fracture.computeConductivity(width);
    
    // Conductivity = k * w
    double expected_conductivity = permeability * width;
    
    EXPECT_NEAR(conductivity, expected_conductivity, 1e-15)
        << "Fracture conductivity calculation";
}

// Test fracture network connectivity
TEST_F(FractureModelTest, NetworkConnectivity) {
    FractureNetwork network;
    
    // Create a simple network
    Fracture frac1;
    frac1.center = {0.0, 0.0, 0.0};
    frac1.strike = 0.0;
    frac1.dip = 90.0;
    frac1.length = 50.0;
    frac1.height = 30.0;
    network.addFracture(frac1);
    
    Fracture frac2;
    frac2.center = {25.0, 0.0, 0.0};
    frac2.strike = 90.0;
    frac2.dip = 90.0;
    frac2.length = 50.0;
    frac2.height = 30.0;
    network.addFracture(frac2);
    
    network.computeConnectivity();
    
    // These fractures should intersect
    EXPECT_TRUE(network.areConnected(0, 1))
        << "Perpendicular fractures should intersect";
    
    double percolation = network.computePercolation();
    EXPECT_GT(percolation, 0.0)
        << "Connected network should have non-zero percolation";
}

// Test fracture aperture evolution with stress
TEST_F(FractureModelTest, ApertureStressEvolution) {
    LEFMFractureModel fracture;
    
    double initial_aperture = 0.0001;  // m (0.1 mm)
    double normal_stress = 10e6;        // Pa
    double stiffness = 100e9;           // Pa/m
    
    fracture.setInitialAperture(initial_aperture);
    fracture.setNormalStiffness(stiffness);
    
    double aperture = fracture.computeApertureUnderStress(normal_stress);
    
    // Linear elastic closure
    double closure = normal_stress / stiffness;
    double expected_aperture = std::max(0.0, initial_aperture - closure);
    
    EXPECT_NEAR(aperture, expected_aperture, 1e-10)
        << "Aperture should decrease with normal stress";
}

// Test fluid-driven fracture propagation
TEST_F(FractureModelTest, FluidDrivenPropagation) {
    LEFMFractureModel fracture;
    
    // Set up hydraulic fracturing scenario
    double E = 20e9;
    double nu = 0.25;
    double K_IC = 1e6;
    double sigma_h = 30e6;  // Minimum horizontal stress
    
    fracture.setMechanicalProperties(E, nu);
    fracture.setFractureToughness(K_IC);
    fracture.setInSituStress(sigma_h, sigma_h, 50e6);
    
    // Injection parameters
    double injection_rate = 0.1;  // m^3/s
    double fluid_viscosity = 0.001;  // Pa·s
    double dt = 10.0;  // s
    
    double initial_length = 10.0;
    double new_length = fracture.propagateHydraulicFracture(
        initial_length, injection_rate, fluid_viscosity, dt);
    
    EXPECT_GT(new_length, initial_length)
        << "Fracture should grow during injection";
    
    EXPECT_LT(new_length, initial_length + injection_rate * dt / 0.001)
        << "Growth should be physically reasonable";
}

// Test KGD model for hydraulic fracture
TEST_F(FractureModelTest, KGDModel) {
    LEFMFractureModel fracture;
    
    double E = 20e9;
    double nu = 0.25;
    double K_IC = 1e6;
    double mu = 0.001;
    double Q = 0.1;
    double t = 60.0;  // 1 minute
    
    fracture.setMechanicalProperties(E, nu);
    fracture.setFractureToughness(K_IC);
    
    auto result = fracture.solveKGDModel(Q, mu, t);
    
    // KGD scaling: L ~ (Q * E' * t / mu)^(1/4)
    double E_prime = E / (1.0 - nu * nu);
    double L_scale = std::pow(Q * E_prime * t / mu, 0.25);
    
    EXPECT_GT(result.length, 0.0) << "Fracture length should be positive";
    EXPECT_GT(result.width, 0.0) << "Fracture width should be positive";
    EXPECT_NEAR(result.length / L_scale, 1.0, 0.5)
        << "Length should scale correctly with KGD model";
}

// Test fracture intersection geometry
TEST_F(FractureModelTest, IntersectionGeometry) {
    // Fracture 1: vertical, E-W oriented
    Fracture frac1;
    frac1.center = {0.0, 0.0, 0.0};
    frac1.strike = 90.0;  // E-W
    frac1.dip = 90.0;
    frac1.length = 100.0;
    frac1.height = 50.0;
    
    // Fracture 2: vertical, N-S oriented
    Fracture frac2;
    frac2.center = {0.0, 0.0, 0.0};
    frac2.strike = 0.0;  // N-S
    frac2.dip = 90.0;
    frac2.length = 100.0;
    frac2.height = 50.0;
    
    FractureNetwork network;
    auto intersection = network.computeIntersection(frac1, frac2);
    
    EXPECT_TRUE(intersection.exists) << "Perpendicular fractures should intersect";
    EXPECT_NEAR(intersection.length, 50.0, 1.0)
        << "Intersection length should equal minimum height";
}

// Test proppant transport in fracture
TEST_F(FractureModelTest, ProppantTransport) {
    LEFMFractureModel fracture;
    
    double width = 0.005;  // m (5 mm)
    double length = 100.0;  // m
    double fluid_velocity = 1.0;  // m/s
    double proppant_diameter = 0.0005;  // m (0.5 mm)
    double proppant_density = 2650.0;  // kg/m^3 (sand)
    double fluid_density = 1000.0;  // kg/m^3
    double fluid_viscosity = 0.001;  // Pa·s
    
    fracture.setFractureGeometry(width, length);
    fracture.setProppantProperties(proppant_diameter, proppant_density);
    
    double settling_velocity = fracture.computeSettlingVelocity(
        proppant_diameter, proppant_density, fluid_density, fluid_viscosity);
    
    // Stokes settling velocity
    double g = 9.81;
    double v_stokes = (proppant_density - fluid_density) * g * 
                      proppant_diameter * proppant_diameter / (18.0 * fluid_viscosity);
    
    EXPECT_NEAR(settling_velocity, v_stokes, v_stokes * 0.1)
        << "Settling velocity should match Stokes law";
    
    // Test proppant placement
    double placement_efficiency = fracture.computePlacementEfficiency(
        fluid_velocity, settling_velocity, length, width);
    
    EXPECT_GT(placement_efficiency, 0.0) << "Should have some proppant placement";
    EXPECT_LE(placement_efficiency, 1.0) << "Efficiency should be <= 100%";
}

// Test leak-off during hydraulic fracturing
TEST_F(FractureModelTest, CarterLeakOff) {
    LEFMFractureModel fracture;
    
    double C_L = 1e-4;  // Leak-off coefficient (m/s^0.5)
    double t = 100.0;    // s
    double area = 1000.0;  // m^2
    
    fracture.setLeakOffCoefficient(C_L);
    
    double leak_off_volume = fracture.computeLeakOff(area, t);
    
    // Carter leak-off: V_L = 2 * C_L * A * sqrt(t)
    double expected_volume = 2.0 * C_L * area * std::sqrt(t);
    
    EXPECT_NEAR(leak_off_volume, expected_volume, 1e-6)
        << "Leak-off should follow Carter model";
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    
    int result = RUN_ALL_TESTS();
    
    PetscFinalize();
    MPI_Finalize();
    return result;
}
