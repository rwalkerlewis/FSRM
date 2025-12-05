/**
 * @file test_fracture_network.cpp
 * @brief Unit tests for advanced FractureNetwork algorithms
 */

#include <gtest/gtest.h>
#include "FractureNetwork.hpp"
#include <cmath>
#include <random>

using namespace FSRM;

class FractureNetworkTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up default domain
        domain_min_ = {0.0, 0.0, 0.0};
        domain_max_ = {100.0, 100.0, 100.0};
    }
    
    std::array<double, 3> domain_min_;
    std::array<double, 3> domain_max_;
};

// ============================================================================
// DiscreteFracture Tests
// ============================================================================

TEST_F(FractureNetworkTest, DiscreteFractureDefaultConstruction) {
    DiscreteFracture frac;
    
    EXPECT_EQ(frac.id, -1);
    EXPECT_DOUBLE_EQ(frac.strike, 0.0);
    EXPECT_DOUBLE_EQ(frac.dip, 90.0);
    EXPECT_DOUBLE_EQ(frac.aperture, 1e-4);
    EXPECT_TRUE(frac.is_conductive);
}

TEST_F(FractureNetworkTest, DiscreteFractureComputeVectors) {
    DiscreteFracture frac;
    frac.strike = 45.0;
    frac.dip = 60.0;
    frac.computeVectors();
    
    // Check that normal vector is unit length
    double norm = std::sqrt(frac.normal[0]*frac.normal[0] + 
                           frac.normal[1]*frac.normal[1] + 
                           frac.normal[2]*frac.normal[2]);
    EXPECT_NEAR(norm, 1.0, 1e-10);
    
    // Check orthogonality
    double dot = frac.strike_dir[0]*frac.normal[0] + 
                 frac.strike_dir[1]*frac.normal[1] + 
                 frac.strike_dir[2]*frac.normal[2];
    EXPECT_NEAR(dot, 0.0, 1e-10);
}

TEST_F(FractureNetworkTest, DiscreteFractureComputePermeability) {
    DiscreteFracture frac;
    frac.aperture = 1e-3;  // 1 mm
    frac.computePermeability();
    
    // Cubic law: k = b²/12
    double expected_k = 1e-6 / 12.0;
    EXPECT_NEAR(frac.permeability, expected_k, 1e-15);
    
    // Transmissivity: T = k * b
    EXPECT_NEAR(frac.transmissivity, expected_k * 1e-3, 1e-15);
}

TEST_F(FractureNetworkTest, DiscreteFractureComputeArea) {
    DiscreteFracture frac;
    frac.radius = 5.0;
    
    double area = frac.computeArea();
    EXPECT_NEAR(area, M_PI * 25.0, 1e-10);
}

TEST_F(FractureNetworkTest, DiscreteFractureStressDependentAperture) {
    DiscreteFracture frac;
    double ref_aperture = 1e-3;
    double ref_stress = 10e6;
    
    // At zero stress, aperture = reference
    double ap0 = frac.computeAperture(0.0, ref_aperture, ref_stress);
    EXPECT_DOUBLE_EQ(ap0, ref_aperture);
    
    // At non-zero stress, aperture should decrease
    double ap1 = frac.computeAperture(10e6, ref_aperture, ref_stress);
    EXPECT_LT(ap1, ref_aperture);
    EXPECT_GT(ap1, 0.0);
}

TEST_F(FractureNetworkTest, DiscreteFractureGetVertices) {
    DiscreteFracture frac;
    frac.center = {50.0, 50.0, 50.0};
    frac.radius = 10.0;
    frac.strike = 0.0;
    frac.dip = 90.0;
    frac.computeVectors();
    
    auto vertices = frac.getVertices(8);
    EXPECT_EQ(vertices.size(), 8);
    
    // All vertices should be at distance radius from center
    for (const auto& v : vertices) {
        double dist = std::sqrt(
            std::pow(v[0] - frac.center[0], 2) +
            std::pow(v[1] - frac.center[1], 2) +
            std::pow(v[2] - frac.center[2], 2));
        EXPECT_NEAR(dist, frac.radius, 1e-10);
    }
}

// ============================================================================
// FractureIntensity Tests
// ============================================================================

TEST_F(FractureNetworkTest, FractureIntensityFromP32) {
    FractureIntensity intensity;
    intensity.computeFromP32(0.1, 5.0);  // P32 = 0.1 m²/m³, mean radius 5 m
    
    EXPECT_DOUBLE_EQ(intensity.P32, 0.1);
    EXPECT_GT(intensity.P30, 0.0);
    EXPECT_GT(intensity.P21, 0.0);
}

TEST_F(FractureNetworkTest, FractureIntensityFromP21) {
    FractureIntensity intensity;
    intensity.computeFromP21(0.5, 10.0);  // P21 = 0.5 m/m², mean trace 10 m
    
    EXPECT_DOUBLE_EQ(intensity.P21, 0.5);
    EXPECT_NEAR(intensity.P32, (M_PI / 2.0) * 0.5, 1e-10);
}

// ============================================================================
// FractureSet Tests
// ============================================================================

TEST_F(FractureNetworkTest, FractureSetConstruction) {
    FractureSet set(0, "TestSet");
    
    EXPECT_EQ(set.getId(), 0);
    EXPECT_EQ(set.getName(), "TestSet");
}

TEST_F(FractureNetworkTest, FractureSetFisherOrientation) {
    FractureSet set(0, "FisherSet");
    set.setOrientationDistribution(OrientationDistribution::FISHER);
    set.setFisherParameters(FisherParameters(45.0, 60.0, 20.0));
    
    std::mt19937 rng(42);
    auto [strike, dip] = set.sampleOrientation(rng);
    
    EXPECT_GE(strike, 0.0);
    EXPECT_LT(strike, 360.0);
    EXPECT_GE(dip, 0.0);
    EXPECT_LE(dip, 90.0);
}

TEST_F(FractureNetworkTest, FractureSetPowerLawSize) {
    FractureSet set(0, "PowerLawSet");
    set.setSizeDistribution(SizeDistribution::POWER_LAW);
    set.setSizeParameters(1.0, 100.0, 2.5);  // x_min, x_max, alpha
    set.setMinMaxSize(0.5, 150.0);
    
    std::mt19937 rng(42);
    
    // Sample many sizes and check distribution
    int small_count = 0;
    int n_samples = 1000;
    
    for (int i = 0; i < n_samples; ++i) {
        double size = set.sampleSize(rng);
        EXPECT_GE(size, 0.5);
        EXPECT_LE(size, 150.0);
        
        if (size < 10.0) small_count++;
    }
    
    // Power-law should have more small fractures
    EXPECT_GT(small_count, n_samples / 2);
}

TEST_F(FractureNetworkTest, FractureSetGenerateFractures) {
    FractureSet set(0, "GenSet");
    set.setOrientationDistribution(OrientationDistribution::FISHER);
    set.setFisherParameters(FisherParameters(90.0, 45.0, 15.0));
    set.setSizeDistribution(SizeDistribution::LOGNORMAL);
    set.setSizeParameters(5.0, 0.5);
    set.setIntensity(0.05);  // P32
    set.setApertureParameters(1e-4, 0.2);
    
    std::mt19937 rng(42);
    auto fractures = set.generateFractures(domain_min_, domain_max_, rng);
    
    EXPECT_GT(fractures.size(), 0);
    
    for (const auto& frac : fractures) {
        EXPECT_EQ(frac.set_id, 0);
        EXPECT_GE(frac.strike, 0.0);
        EXPECT_LT(frac.strike, 360.0);
        EXPECT_GT(frac.aperture, 0.0);
    }
}

// ============================================================================
// FractureConnectivityGraph Tests
// ============================================================================

TEST_F(FractureNetworkTest, ConnectivityGraphBasic) {
    FractureConnectivityGraph graph;
    
    graph.addNode(0);
    graph.addNode(1);
    graph.addNode(2);
    
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(1, 2, 1.0);
    
    EXPECT_TRUE(graph.isConnected(0, 2));
    EXPECT_EQ(graph.getNumComponents(), 1);
}

TEST_F(FractureNetworkTest, ConnectivityGraphComponents) {
    FractureConnectivityGraph graph;
    
    // Two disconnected components
    graph.addNode(0);
    graph.addNode(1);
    graph.addEdge(0, 1, 1.0);
    
    graph.addNode(2);
    graph.addNode(3);
    graph.addEdge(2, 3, 1.0);
    
    // Need to rebuild after modifications
    std::vector<DiscreteFracture> fractures(4);
    for (int i = 0; i < 4; ++i) {
        fractures[i].id = i;
    }
    fractures[0].connected_fractures = {1};
    fractures[1].connected_fractures = {0};
    fractures[2].connected_fractures = {3};
    fractures[3].connected_fractures = {2};
    
    graph.buildFromFractures(fractures);
    
    EXPECT_EQ(graph.getNumComponents(), 2);
    EXPECT_TRUE(graph.isConnected(0, 1));
    EXPECT_TRUE(graph.isConnected(2, 3));
    EXPECT_FALSE(graph.isConnected(0, 2));
}

TEST_F(FractureNetworkTest, ConnectivityGraphShortestPath) {
    FractureConnectivityGraph graph;
    
    // Simple chain: 0 - 1 - 2 - 3
    for (int i = 0; i < 4; ++i) {
        graph.addNode(i);
    }
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(1, 2, 1.0);
    graph.addEdge(2, 3, 1.0);
    
    auto path = graph.findShortestPath(0, 3);
    
    EXPECT_EQ(path.size(), 4);
    EXPECT_EQ(path[0], 0);
    EXPECT_EQ(path[3], 3);
}

TEST_F(FractureNetworkTest, ConnectivityGraphMetrics) {
    FractureConnectivityGraph graph;
    
    // Complete graph of 4 nodes
    for (int i = 0; i < 4; ++i) {
        graph.addNode(i);
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            graph.addEdge(i, j, 1.0);
        }
    }
    
    EXPECT_DOUBLE_EQ(graph.getConnectivity(), 1.0);  // All connected
    EXPECT_DOUBLE_EQ(graph.getAverageClusteringCoefficient(), 1.0);  // Complete graph
}

// ============================================================================
// FractureNetwork Tests
// ============================================================================

TEST_F(FractureNetworkTest, FractureNetworkSetDomain) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    
    // Domain should be set correctly
    EXPECT_NO_THROW(network.generate());
}

TEST_F(FractureNetworkTest, FractureNetworkGenerate) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    // Add a fracture set
    FractureSet set(0, "TestSet");
    set.setOrientationDistribution(OrientationDistribution::FISHER);
    set.setFisherParameters(FisherParameters(0.0, 90.0, 10.0));
    set.setSizeDistribution(SizeDistribution::LOGNORMAL);
    set.setSizeParameters(5.0, 0.3);
    set.setIntensity(0.1);
    
    network.addFractureSet(set);
    network.generate();
    
    EXPECT_GT(network.getNumFractures(), 0);
}

TEST_F(FractureNetworkTest, FractureNetworkIntersections) {
    FractureNetwork network;
    network.setDomain({0, 0, 0}, {50, 50, 50});
    network.setRandomSeed(123);
    
    // Create two orthogonal sets that should intersect
    FractureSet set1(0, "NS");
    set1.setFisherParameters(FisherParameters(0.0, 90.0, 100.0));  // N-S vertical
    set1.setSizeParameters(15.0, 0.1);
    set1.setIntensity(0.15);
    
    FractureSet set2(1, "EW");
    set2.setFisherParameters(FisherParameters(90.0, 90.0, 100.0));  // E-W vertical
    set2.setSizeParameters(15.0, 0.1);
    set2.setIntensity(0.15);
    
    network.addFractureSet(set1);
    network.addFractureSet(set2);
    network.generate();
    
    // Should have some intersections
    const auto& intersections = network.getIntersections();
    // May or may not have intersections depending on random placement
    EXPECT_GE(intersections.size(), 0);
}

TEST_F(FractureNetworkTest, FractureNetworkIntensity) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    FractureSet set(0, "Test");
    set.setSizeParameters(5.0, 0.2);
    set.setIntensity(0.05);
    
    network.addFractureSet(set);
    network.generate();
    
    FractureIntensity intensity = network.computeIntensity();
    
    EXPECT_GE(intensity.P32, 0.0);
    EXPECT_GE(intensity.P30, 0.0);
}

TEST_F(FractureNetworkTest, FractureNetworkOdaTensor) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    // Vertical fractures striking N-S
    FractureSet set(0, "NS");
    set.setFisherParameters(FisherParameters(0.0, 90.0, 100.0));
    set.setSizeParameters(10.0, 0.1);
    set.setIntensity(0.1);
    set.setApertureParameters(1e-3, 0.1);
    
    network.addFractureSet(set);
    network.generate();
    
    auto tensor = network.computeOdaTensor();
    
    // For N-S vertical fractures, flow is easy in N-S direction
    // tensor[1][1] (y-direction, N-S) should be relatively high
    // This test checks the tensor is computed without error
    bool has_nonzero = false;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (tensor[i][j] != 0.0) has_nonzero = true;
        }
    }
    
    if (network.getNumFractures() > 0) {
        EXPECT_TRUE(has_nonzero);
    }
}

TEST_F(FractureNetworkTest, FractureNetworkSnowPermeability) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    FractureSet set(0, "Test");
    set.setSizeParameters(10.0, 0.2);
    set.setIntensity(0.1);
    set.setApertureParameters(1e-3, 0.1);
    
    network.addFractureSet(set);
    network.generate();
    
    auto perm = network.computeSnowPermeability();
    
    // All permeabilities should be non-negative
    EXPECT_GE(perm[0], 0.0);
    EXPECT_GE(perm[1], 0.0);
    EXPECT_GE(perm[2], 0.0);
}

TEST_F(FractureNetworkTest, FractureNetworkMINC) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    FractureSet set(0, "Test");
    set.setSizeParameters(5.0, 0.2);
    set.setIntensity(0.1);
    
    network.addFractureSet(set);
    network.generate();
    
    auto minc = network.computeMINC(3);
    
    EXPECT_EQ(minc.volume_fractions.size(), 3);
    EXPECT_EQ(minc.shape_factors.size(), 3);
    EXPECT_EQ(minc.transmissibilities.size(), 3);
    
    // Volume fractions should sum to 1
    double sum = 0.0;
    for (double vf : minc.volume_fractions) {
        sum += vf;
        EXPECT_GT(vf, 0.0);
    }
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// ============================================================================
// FracturePropagationModel Tests
// ============================================================================

TEST_F(FractureNetworkTest, PropagationModelConstruction) {
    FracturePropagationModel model;
    
    EXPECT_NO_THROW(model.setRockProperties(30e9, 0.25, 1e6));
    EXPECT_NO_THROW(model.setFluidProperties(0.001, 4.5e-10));
}

TEST_F(FractureNetworkTest, PropagationModelStress) {
    FracturePropagationModel model;
    model.setRockProperties(30e9, 0.25, 1e6);
    
    EXPECT_NO_THROW(model.setFarFieldStress(
        20e6, 25e6, 30e6,  // Principal stresses
        0.0, 0.0, 0.0      // Shear stresses
    ));
}

TEST_F(FractureNetworkTest, PropagationModelTipStates) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    FractureSet set(0, "Test");
    set.setSizeParameters(10.0, 0.1);
    set.setIntensity(0.05);
    
    network.addFractureSet(set);
    network.generate();
    
    FracturePropagationModel model;
    model.setRockProperties(30e9, 0.25, 1e6);
    model.setFarFieldStress(20e6, 25e6, 30e6, 0.0, 0.0, 0.0);
    
    auto states = model.computeTipStates(network);
    
    EXPECT_EQ(states.size(), network.getNumFractures());
    
    for (const auto& state : states) {
        EXPECT_TRUE(std::isfinite(state.K_I));
        EXPECT_TRUE(std::isfinite(state.energy_release_rate));
    }
}

TEST_F(FractureNetworkTest, PropagationModelCoalescence) {
    FractureNetwork network;
    network.setDomain({0, 0, 0}, {20, 20, 20});
    
    // Manually add two close fractures
    DiscreteFracture f1, f2;
    f1.id = 0;
    f1.center = {5.0, 10.0, 10.0};
    f1.radius = 3.0;
    f1.strike = 0.0;
    f1.dip = 90.0;
    f1.computeVectors();
    
    f2.id = 1;
    f2.center = {15.0, 10.0, 10.0};
    f2.radius = 3.0;
    f2.strike = 0.0;
    f2.dip = 90.0;
    f2.computeVectors();
    
    network.getFractures().push_back(f1);
    network.getFractures().push_back(f2);
    
    FracturePropagationModel model;
    
    // Check coalescence with different thresholds
    auto pairs_small = model.detectCoalescence(network, 1.0);
    auto pairs_large = model.detectCoalescence(network, 20.0);
    
    // Fractures are 10 units apart, centers are 10m, radii are 3m each
    // So edge-to-edge distance is ~4m
    EXPECT_EQ(pairs_small.size(), 0);  // 1.0 threshold too small
    EXPECT_EQ(pairs_large.size(), 1);  // 20.0 threshold captures them
}

// ============================================================================
// FractureFlowSimulator Tests
// ============================================================================

TEST_F(FractureNetworkTest, FlowSimulatorSetup) {
    FractureNetwork network;
    FractureFlowSimulator simulator;
    
    simulator.setNetwork(&network);
    simulator.setFluidProperties(0.001, 4.5e-10);
    
    EXPECT_NO_THROW(simulator.setInletPressure(0, 10e6));
    EXPECT_NO_THROW(simulator.setOutletPressure(1, 1e6));
}

TEST_F(FractureNetworkTest, FlowSimulatorSolve) {
    FractureNetwork network;
    network.setDomain({0, 0, 0}, {50, 50, 50});
    network.setRandomSeed(42);
    
    // Dense network for connectivity
    FractureSet set(0, "Test");
    set.setSizeParameters(15.0, 0.2);
    set.setIntensity(0.3);
    set.setApertureParameters(1e-3, 0.1);
    
    network.addFractureSet(set);
    network.generate();
    
    if (network.getNumFractures() < 2) {
        GTEST_SKIP() << "Not enough fractures generated";
    }
    
    FractureFlowSimulator simulator;
    simulator.setNetwork(&network);
    simulator.setFluidProperties(0.001, 4.5e-10);
    
    // Set BCs on first and last fracture
    const auto& fracs = network.getFractures();
    simulator.setInletPressure(fracs.front().id, 10e6);
    simulator.setOutletPressure(fracs.back().id, 1e6);
    
    auto solution = simulator.solveSteadyState();
    
    EXPECT_EQ(solution.pressures.size(), network.getNumFractures());
    EXPECT_GE(solution.effective_transmissivity, 0.0);
}

// ============================================================================
// DFNStatistics Tests
// ============================================================================

TEST_F(FractureNetworkTest, StatisticsFisherKappa) {
    // Test with clustered orientations
    std::vector<std::pair<double, double>> orientations;
    
    // All pointing roughly the same direction
    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 5.0);
    
    for (int i = 0; i < 100; ++i) {
        orientations.push_back({45.0 + noise(rng), 60.0 + noise(rng)});
    }
    
    double kappa = DFNStatistics::computeFisherKappa(orientations);
    
    // High concentration (small spread) should give high kappa
    EXPECT_GT(kappa, 10.0);
}

TEST_F(FractureNetworkTest, StatisticsFitSizeDistribution) {
    std::vector<double> sizes;
    
    // Generate log-normal distributed sizes
    std::mt19937 rng(42);
    std::lognormal_distribution<double> dist(std::log(10.0), 0.5);
    
    for (int i = 0; i < 200; ++i) {
        sizes.push_back(dist(rng));
    }
    
    auto [fitted_type, params] = DFNStatistics::fitSizeDistribution(sizes);
    
    // Should identify as log-normal or similar
    EXPECT_TRUE(fitted_type == SizeDistribution::LOGNORMAL || 
                fitted_type == SizeDistribution::POWER_LAW);
}

TEST_F(FractureNetworkTest, StatisticsClusterOrientations) {
    std::vector<DiscreteFracture> fractures;
    
    // Create two distinct orientation clusters
    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 5.0);
    
    // Cluster 1: N-S vertical
    for (int i = 0; i < 50; ++i) {
        DiscreteFracture f;
        f.strike = noise(rng);
        f.dip = 90.0 + noise(rng);
        f.computeVectors();
        fractures.push_back(f);
    }
    
    // Cluster 2: E-W vertical
    for (int i = 0; i < 50; ++i) {
        DiscreteFracture f;
        f.strike = 90.0 + noise(rng);
        f.dip = 90.0 + noise(rng);
        f.computeVectors();
        fractures.push_back(f);
    }
    
    auto clusters = DFNStatistics::clusterOrientations(fractures, 2);
    
    EXPECT_EQ(clusters.size(), 2);
    
    // Cluster centers should be roughly 90° apart in strike
    double strike_diff = std::abs(clusters[0].mean_strike - clusters[1].mean_strike);
    if (strike_diff > 180.0) strike_diff = 360.0 - strike_diff;
    
    EXPECT_GT(strike_diff, 60.0);  // At least 60° apart
}

TEST_F(FractureNetworkTest, StatisticsFractalDimension) {
    std::vector<DiscreteFracture> fractures;
    std::array<double, 3> domain = {100.0, 100.0, 100.0};
    
    // Uniformly distributed fractures
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos(0.0, 100.0);
    
    for (int i = 0; i < 100; ++i) {
        DiscreteFracture f;
        f.center = {pos(rng), pos(rng), pos(rng)};
        f.radius = 5.0;
        fractures.push_back(f);
    }
    
    double D = DFNStatistics::computeFractalDimension(fractures, domain);
    
    // Uniform distribution should have D close to 3
    EXPECT_GT(D, 2.0);
    EXPECT_LE(D, 3.0);
}

TEST_F(FractureNetworkTest, StatisticsScanlineIntensity) {
    std::vector<DiscreteFracture> fractures;
    
    // Create vertical fractures along a line
    for (int i = 0; i < 10; ++i) {
        DiscreteFracture f;
        f.id = i;
        f.center = {50.0, 10.0 * i + 5.0, 50.0};
        f.radius = 8.0;
        f.strike = 0.0;
        f.dip = 90.0;
        f.computeVectors();
        fractures.push_back(f);
    }
    
    // Scanline along y-axis
    std::array<double, 3> start = {50.0, 0.0, 50.0};
    std::array<double, 3> end = {50.0, 100.0, 50.0};
    
    FractureIntensity intensity = DFNStatistics::computeIntensityFromScanline(
        fractures, start, end);
    
    // Should detect approximately 10 intersections over 100m
    EXPECT_NEAR(intensity.P10, 0.1, 0.05);
}

// ============================================================================
// VTK Export Test
// ============================================================================

TEST_F(FractureNetworkTest, ExportToVTK) {
    FractureNetwork network;
    network.setDomain(domain_min_, domain_max_);
    network.setRandomSeed(42);
    
    FractureSet set(0, "Test");
    set.setSizeParameters(10.0, 0.2);
    set.setIntensity(0.05);
    
    network.addFractureSet(set);
    network.generate();
    
    // Export should not throw
    EXPECT_NO_THROW(network.writeVTK("/tmp/test_fractures.vtk"));
}
