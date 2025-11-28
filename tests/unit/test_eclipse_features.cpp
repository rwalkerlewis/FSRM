/**
 * @file test_eclipse_features.cpp
 * @brief Unit tests for ECLIPSE-compatible features
 * 
 * Tests for:
 * - Aquifer models (Carter-Tracy, Fetkovich)
 * - Tracer simulation
 * - VFP tables
 * - Group control
 * - Relative permeability (hysteresis, 3-phase)
 * - Summary output
 */

#include <gtest/gtest.h>
#include <cmath>
#include <memory>

// Include ECLIPSE feature headers
#include "AquiferModel.hpp"
#include "TracerModel.hpp"
#include "VFPTables.hpp"
#include "GroupControl.hpp"
#include "RelativePermeability.hpp"
#include "SummaryOutput.hpp"

using namespace FSRM;

// ============================================================================
// Aquifer Model Tests
// ============================================================================

class AquiferModelTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(AquiferModelTest, CarterTracyInitialization) {
    CarterTracyAquifer aquifer(1);
    
    EXPECT_EQ(aquifer.getId(), 1);
    EXPECT_EQ(aquifer.getType(), AquiferType::CARTER_TRACY);
    EXPECT_GT(aquifer.getInfluxConstant(), 0.0);
}

TEST_F(AquiferModelTest, CarterTracyInfluxCalculation) {
    CarterTracyAquifer aquifer(1);
    
    // Set typical aquifer properties
    aquifer.setAquiferProperties(
        100.0,      // permeability (mD)
        0.2,        // porosity
        1e-9,       // compressibility (1/Pa)
        0.001,      // viscosity (Pa·s)
        50.0,       // thickness (m)
        0.5         // encroachment angle (half circle)
    );
    aquifer.setAquiferRadius(1000.0, 10000.0);
    aquifer.setInitialPressure(30e6);  // 30 MPa
    
    // Calculate influx when reservoir pressure drops
    double reservoir_pressure = 25e6;  // 25 MPa (5 MPa drop)
    double dt = 86400.0;  // 1 day
    
    double influx = aquifer.calculateInflux(reservoir_pressure, dt);
    
    // Influx should be positive (water entering reservoir)
    EXPECT_GT(influx, 0.0);
}

TEST_F(AquiferModelTest, FetkovichInfluxCalculation) {
    FetkovichAquifer aquifer(2);
    
    aquifer.setProductivityIndex(1e-6);  // m³/s/Pa
    aquifer.setInitialVolume(1e9);       // 1 billion m³
    aquifer.setInitialPressure(30e6);
    
    // Test influx
    double reservoir_pressure = 25e6;
    double dt = 86400.0;
    
    double influx = aquifer.calculateInflux(reservoir_pressure, dt);
    
    EXPECT_GT(influx, 0.0);
    
    // Update and check cumulative
    aquifer.updateState(influx, dt);
    EXPECT_GT(aquifer.getCumulativeInflux(), 0.0);
}

TEST_F(AquiferModelTest, ConstantPressureAquifer) {
    ConstantPressureAquifer aquifer(3);
    
    aquifer.setConstantPressure(30e6);
    aquifer.setTransmissibility(1e-6);
    
    // Test influx at different reservoir pressures
    double influx1 = aquifer.calculateInflux(25e6, 86400.0);
    double influx2 = aquifer.calculateInflux(20e6, 86400.0);
    
    // Higher pressure drop should give higher influx
    EXPECT_GT(influx2, influx1);
}

TEST_F(AquiferModelTest, AquiferManager) {
    AquiferManager manager;
    
    // Add different aquifer types
    std::map<std::string, std::string> ct_config = {
        {"permeability", "100.0"},
        {"porosity", "0.2"},
        {"initial_pressure", "30e6"}
    };
    manager.addCarterTracyAquifer(1, ct_config);
    
    std::map<std::string, std::string> fet_config = {
        {"productivity_index", "1e-6"},
        {"initial_volume", "1e9"},
        {"initial_pressure", "30e6"}
    };
    manager.addFetkovichAquifer(2, fet_config);
    
    // Test retrieval
    auto* aq1 = manager.getAquifer(1);
    auto* aq2 = manager.getAquifer(2);
    
    EXPECT_NE(aq1, nullptr);
    EXPECT_NE(aq2, nullptr);
    EXPECT_EQ(aq1->getType(), AquiferType::CARTER_TRACY);
    EXPECT_EQ(aq2->getType(), AquiferType::FETKOVICH);
}

// ============================================================================
// Tracer Model Tests
// ============================================================================

class TracerModelTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(TracerModelTest, TracerPropertiesDefaults) {
    TracerProperties props;
    
    EXPECT_EQ(props.phase, TracerPhase::WATER);
    EXPECT_EQ(props.behavior, TracerBehavior::PASSIVE);
    EXPECT_GT(props.molecular_weight, 0.0);
    EXPECT_GT(props.diffusion_coefficient, 0.0);
}

TEST_F(TracerModelTest, TracerRetardationFactor) {
    TracerProperties props;
    props.adsorption_model = AdsorptionModel::LINEAR;
    props.Kd = 0.001;  // m³/kg
    
    double porosity = 0.2;
    double rock_density = 2650.0;
    
    double R = props.getRetardationFactor(porosity, rock_density);
    
    // R = 1 + (ρb/φ) * Kd
    // With φ=0.2, ρ=2650, Kd=0.001:
    // R = 1 + (2650*0.8/0.2) * 0.001 ≈ 11.6
    EXPECT_GT(R, 1.0);
    EXPECT_NEAR(R, 11.6, 1.0);
}

TEST_F(TracerModelTest, TracerDispersion) {
    TracerProperties props;
    props.diffusion_coefficient = 1e-9;
    props.dispersivity_L = 1.0;
    
    double velocity = 1e-5;  // m/s
    double porosity = 0.2;
    
    double D = props.getDispersion(velocity, porosity);
    
    EXPECT_GT(D, props.diffusion_coefficient);  // Mechanical dispersion adds to molecular
}

TEST_F(TracerModelTest, BreakthroughCurveStatistics) {
    TracerBreakthroughCurve curve;
    curve.tracer_name = "TEST";
    curve.location_name = "PROD1";
    
    // Create synthetic breakthrough curve
    for (int i = 0; i < 100; ++i) {
        double t = i * 86400.0;  // Time in seconds
        double c = std::exp(-std::pow((t - 50*86400.0) / (10*86400.0), 2));
        curve.times.push_back(t);
        curve.concentrations.push_back(c);
    }
    
    curve.computeStatistics(100.0);  // 100 kg injected
    
    EXPECT_GT(curve.mean_residence_time, 0.0);
    EXPECT_GT(curve.peak_concentration, 0.0);
    EXPECT_GE(curve.breakthrough_time, 0.0);
}

TEST_F(TracerModelTest, TracerManager) {
    TracerManager manager;
    manager.setGrid(10, 10, 1, 10.0, 10.0, 1.0);
    
    manager.addWaterTracer("TRACER1", TracerBehavior::PASSIVE);
    manager.addPartitioningTracer("TRACER2", 2.0);
    manager.addRadioactiveTracer("TRACER3", 100.0);  // 100 day half-life
    
    // This should work without errors
    EXPECT_NO_THROW(manager.addObservationWell("OBS1", 5, 5, 0));
}

// ============================================================================
// VFP Tables Tests
// ============================================================================

class VFPTablesTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(VFPTablesTest, VFPProdTableCreation) {
    VFPProdTable table(1);
    
    table.setTHPValues({1e6, 2e6, 3e6});
    table.setFlowRateValues({0.001, 0.01, 0.1});
    table.setWaterCutValues({0.0, 0.5, 0.9});
    table.setGORValues({50.0, 100.0, 200.0});
    table.setALQValues({0.0});  // No gas lift
    
    EXPECT_EQ(table.getTableNumber(), 1);
    EXPECT_EQ(table.getTHPValues().size(), 3);
    EXPECT_EQ(table.getFlowRateValues().size(), 3);
}

TEST_F(VFPTablesTest, VFPTableGenerator) {
    VFPTableGenerator generator;
    
    generator.setTotalDepth(2000.0);
    generator.setWellboreRadius(0.1);
    generator.setOilProperties(800.0, 0.005);
    generator.setWaterProperties(1000.0, 0.001);
    generator.setGasProperties(1.0, 0.00001);
    
    std::vector<double> thp_range = {1e6, 2e6, 3e6};
    std::vector<double> rate_range = {0.001, 0.01, 0.1};
    std::vector<double> wct_range = {0.0, 0.5, 0.9};
    std::vector<double> gor_range = {50.0, 100.0};
    std::vector<double> alq_range = {0.0};
    
    auto table = generator.generateProductionTable(
        thp_range, rate_range, wct_range, gor_range, alq_range);
    
    EXPECT_TRUE(table.isValid());
}

TEST_F(VFPTablesTest, VFPBHPCalculation) {
    VFPTableGenerator generator;
    generator.setTotalDepth(2000.0);
    generator.setWellboreRadius(0.1);
    generator.setOilProperties(800.0, 0.005);
    generator.setWaterProperties(1000.0, 0.001);
    
    auto table = generator.generateProductionTable(
        {1e6, 2e6, 3e6},
        {0.001, 0.01, 0.1},
        {0.0, 0.5},
        {50.0, 100.0},
        {0.0}
    );
    
    // BHP should increase with THP
    double bhp1 = table.calculateBHP(1e6, 0.01, 0.0, 50.0, 0.0);
    double bhp2 = table.calculateBHP(2e6, 0.01, 0.0, 50.0, 0.0);
    
    EXPECT_GT(bhp2, bhp1);
    
    // BHP should be higher than THP (due to hydrostatic)
    EXPECT_GT(bhp1, 1e6);
}

TEST_F(VFPTablesTest, VFPInjTable) {
    VFPInjTable table(1);
    
    table.setTHPValues({1e6, 2e6, 3e6});
    table.setFlowRateValues({0.001, 0.01, 0.1});
    
    // Simple linear BHP data
    std::vector<double> bhp_data;
    for (double thp : {1e6, 2e6, 3e6}) {
        for (double rate : {0.001, 0.01, 0.1}) {
            bhp_data.push_back(thp + 15e6 + rate * 1e8);  // Simple model
        }
    }
    table.setBHPData(bhp_data);
    
    EXPECT_TRUE(table.isValid());
    
    double bhp = table.calculateBHP(2e6, 0.01);
    EXPECT_GT(bhp, 2e6);  // BHP > THP for injection
}

// ============================================================================
// Group Control Tests
// ============================================================================

class GroupControlTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(GroupControlTest, GroupCreation) {
    WellGroup group("GROUP1", "FIELD");
    
    EXPECT_EQ(group.getName(), "GROUP1");
    EXPECT_EQ(group.getParent(), "FIELD");
    EXPECT_FALSE(group.isField());
}

TEST_F(GroupControlTest, GroupHierarchy) {
    GroupControlManager manager;
    
    manager.createGroup("PLATFORM_A", "FIELD");
    manager.createGroup("REGION_NORTH", "PLATFORM_A");
    
    manager.assignWellToGroup("WELL1", "REGION_NORTH");
    manager.assignWellToGroup("WELL2", "REGION_NORTH");
    
    EXPECT_EQ(manager.getWellGroup("WELL1"), "REGION_NORTH");
    
    auto* group = manager.getGroup("REGION_NORTH");
    ASSERT_NE(group, nullptr);
    EXPECT_EQ(group->getWells().size(), 2);
}

TEST_F(GroupControlTest, ProductionConstraints) {
    GroupControlManager manager;
    manager.createGroup("PROD_GROUP", "FIELD");
    
    GroupProdConstraints constraints;
    constraints.control_mode = GroupProdControlMode::ORAT;
    constraints.oil_rate_target = 0.1;  // m³/s
    constraints.oil_rate_max = 0.15;
    constraints.water_cut_max = 0.9;
    
    manager.setProductionConstraints("PROD_GROUP", constraints);
    
    auto* group = manager.getGroup("PROD_GROUP");
    ASSERT_NE(group, nullptr);
    
    const auto& c = group->getProductionConstraints();
    EXPECT_EQ(c.control_mode, GroupProdControlMode::ORAT);
    EXPECT_NEAR(c.oil_rate_target, 0.1, 1e-10);
}

TEST_F(GroupControlTest, GuideRateCalculation) {
    GuideRateCalculator calc;
    
    WellGroup group("TEST", "FIELD");
    group.addWell("W1");
    group.addWell("W2");
    group.addWell("W3");
    
    GroupProdConstraints constraints;
    constraints.guide_rate_type = GuideRateType::OIL;
    constraints.guide_rate_defined = true;
    group.setProductionConstraints(constraints);
    
    // Well potentials: {oil, water, gas}
    std::map<std::string, std::array<double, 3>> potentials = {
        {"W1", {0.1, 0.01, 10.0}},
        {"W2", {0.2, 0.02, 20.0}},
        {"W3", {0.05, 0.005, 5.0}}
    };
    
    auto guide_rates = calc.calculateGuideRates(group, potentials);
    
    // Total should sum to 1.0
    double total = 0.0;
    for (const auto& [well, rate] : guide_rates) {
        total += rate;
        EXPECT_GT(rate, 0.0);
    }
    EXPECT_NEAR(total, 1.0, 1e-10);
    
    // W2 should have highest guide rate (highest oil potential)
    EXPECT_GT(guide_rates["W2"], guide_rates["W1"]);
    EXPECT_GT(guide_rates["W2"], guide_rates["W3"]);
}

// ============================================================================
// Relative Permeability Tests
// ============================================================================

class RelativePermeabilityTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(RelativePermeabilityTest, CoreyModel) {
    CoreyRelPerm corey;
    
    EndPoints wet_ep, nw_ep;
    wet_ep.S_connate = 0.2;
    wet_ep.S_max = 1.0;
    wet_ep.kr_max = 1.0;
    
    nw_ep.S_connate = 0.0;
    nw_ep.S_max = 0.8;
    nw_ep.kr_max = 1.0;
    
    corey.setEndPoints(wet_ep, nw_ep);
    corey.setWettingExponent(2.0);
    corey.setNonWettingExponent(2.0);
    
    // At connate water, kr_w should be 0
    EXPECT_NEAR(corey.kr_wetting(0.2), 0.0, 1e-10);
    
    // At max water saturation, kr_w should be max
    EXPECT_NEAR(corey.kr_wetting(1.0), 1.0, 1e-10);
    
    // At connate water, kr_nw should be max
    EXPECT_NEAR(corey.kr_nonwetting(0.2), 1.0, 1e-10);
}

TEST_F(RelativePermeabilityTest, ThreePhaseStoneI) {
    auto water_oil = std::make_shared<CoreyRelPerm>();
    auto gas_oil = std::make_shared<CoreyRelPerm>();
    
    // Set up two-phase curves
    EndPoints wet_ep, nw_ep;
    wet_ep.S_connate = 0.2;
    wet_ep.kr_max = 1.0;
    nw_ep.S_connate = 0.2;
    nw_ep.kr_max = 1.0;
    
    water_oil->setEndPoints(wet_ep, nw_ep);
    gas_oil->setEndPoints(wet_ep, nw_ep);
    
    ThreePhaseRelPerm three_phase;
    three_phase.setWaterOilTable(water_oil);
    three_phase.setGasOilTable(gas_oil);
    three_phase.setThreePhaseModel(ThreePhaseModel::STONE_I);
    
    // Test at single phase (all oil)
    auto kr = three_phase.calculate(0.2, 0.8, 0.0);
    
    EXPECT_NEAR(kr[0], 0.0, 1e-10);  // kr_w at connate
    EXPECT_GT(kr[1], 0.0);           // kr_o should be positive
    EXPECT_NEAR(kr[2], 0.0, 1e-10);  // kr_g at zero gas
}

TEST_F(RelativePermeabilityTest, HysteresisKillough) {
    auto drainage = std::make_shared<CoreyRelPerm>();
    
    EndPoints wet_ep, nw_ep;
    wet_ep.S_connate = 0.2;
    wet_ep.kr_max = 1.0;
    nw_ep.S_connate = 0.2;
    nw_ep.kr_max = 1.0;
    
    drainage->setEndPoints(wet_ep, nw_ep);
    
    HysteresisRelPerm hysteresis(drainage, nullptr, HysteresisModel::KILLOUGH);
    hysteresis.setKilloughCurvature(0.1);
    hysteresis.setLandParameter(2.0);
    
    // Drainage scan
    hysteresis.updateState(0.3);
    double kr1 = hysteresis.kr_wetting(0.3);
    
    hysteresis.updateState(0.5);
    double kr2 = hysteresis.kr_wetting(0.5);
    
    // Imbibition scan (decreasing saturation)
    hysteresis.updateState(0.4);
    double kr3 = hysteresis.kr_wetting(0.4);
    
    // Basic sanity checks
    EXPECT_GT(kr2, kr1);  // Higher saturation = higher kr
    EXPECT_GE(kr3, 0.0);  // Non-negative
}

TEST_F(RelativePermeabilityTest, TabularRelPerm) {
    TabularRelPerm table;
    
    std::vector<double> S = {0.2, 0.4, 0.6, 0.8, 1.0};
    std::vector<double> krw = {0.0, 0.05, 0.2, 0.5, 1.0};
    std::vector<double> kro = {1.0, 0.6, 0.2, 0.05, 0.0};
    
    table.setTable(S, krw, kro);
    
    // Test interpolation
    EXPECT_NEAR(table.kr_wetting(0.2), 0.0, 1e-10);
    EXPECT_NEAR(table.kr_wetting(1.0), 1.0, 1e-10);
    EXPECT_NEAR(table.kr_nonwetting(0.2), 1.0, 1e-10);
    EXPECT_NEAR(table.kr_nonwetting(1.0), 0.0, 1e-10);
    
    // Interpolated value
    double krw_mid = table.kr_wetting(0.5);
    EXPECT_GT(krw_mid, 0.05);
    EXPECT_LT(krw_mid, 0.2);
}

// ============================================================================
// Summary Output Tests
// ============================================================================

class SummaryOutputTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(SummaryOutputTest, SummaryVectorCreation) {
    SummaryVector vec("WOPR", SummaryCategory::WELL, 
                      SummaryType::OIL_PRODUCTION_RATE, "PROD1", "m3/s");
    
    EXPECT_EQ(vec.keyword, "WOPR");
    EXPECT_EQ(vec.category, SummaryCategory::WELL);
    EXPECT_EQ(vec.entity_name, "PROD1");
}

TEST_F(SummaryOutputTest, SummarySpec) {
    SummarySpec spec;
    spec.setWells({"PROD1", "PROD2", "INJ1"});
    
    spec.addKeyword("FOPR");
    spec.addWellKeyword("WOPR", "PROD1");
    spec.addWellKeyword("WBHP", "PROD1");
    spec.addFieldVectors();
    
    EXPECT_TRUE(spec.hasVector("FOPR"));
    
    auto keywords = spec.getAllKeywords();
    EXPECT_GT(keywords.size(), 0);
}

TEST_F(SummaryOutputTest, SummaryCollector) {
    SummarySpec spec;
    spec.addKeyword("FOPR");
    spec.addKeyword("TIME");
    
    SummaryCollector collector;
    collector.setSpec(spec);
    
    // Record some steps
    for (int i = 0; i < 5; ++i) {
        collector.startStep(i * 86400.0, 86400.0, i);
        collector.recordFieldRate(SummaryType::OIL_PRODUCTION_RATE, 0.1 + i * 0.01);
        collector.endStep();
    }
    
    EXPECT_EQ(collector.getNumSteps(), 5);
    
    auto times = collector.getTimes();
    EXPECT_EQ(times.size(), 5);
    EXPECT_NEAR(times[0], 0.0, 1e-10);
    EXPECT_NEAR(times[4], 4 * 86400.0, 1e-10);
}

TEST_F(SummaryOutputTest, SummaryOutput) {
    SummaryOutput output("TEST_CASE");
    
    output.setWells({"PROD1", "PROD2", "INJ1"});
    output.addDefaultOutput();
    output.addWellOutput("PROD1");
    
    // Don't initialize (would create files) - just test setup
    EXPECT_NO_THROW(output.addFieldOutput());
}

// ============================================================================
// Integration Tests
// ============================================================================

class EclipseFeaturesIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(EclipseFeaturesIntegrationTest, AquiferWithSummary) {
    // Create aquifer
    FetkovichAquifer aquifer(1);
    aquifer.setProductivityIndex(1e-6);
    aquifer.setInitialVolume(1e9);
    aquifer.setInitialPressure(30e6);
    
    // Create summary collector
    SummarySpec spec;
    spec.addAquiferKeyword("AAQP", 1);
    spec.addAquiferKeyword("AAQR", 1);
    spec.addAquiferKeyword("AAQT", 1);
    
    SummaryCollector collector;
    collector.setSpec(spec);
    
    // Simulate
    double reservoir_pressure = 28e6;
    double dt = 86400.0;
    
    for (int step = 0; step < 10; ++step) {
        collector.startStep(step * dt, dt, step);
        
        double influx = aquifer.calculateInflux(reservoir_pressure, dt);
        aquifer.updateState(influx, dt);
        
        collector.recordAquiferData(1, SummaryType::AQUIFER_PRESSURE, 
                                    aquifer.getCurrentPressure());
        collector.recordAquiferData(1, SummaryType::AQUIFER_INFLUX_RATE, influx);
        collector.recordAquiferData(1, SummaryType::AQUIFER_INFLUX_TOTAL, 
                                    aquifer.getCumulativeInflux());
        
        collector.endStep();
    }
    
    EXPECT_EQ(collector.getNumSteps(), 10);
    EXPECT_GT(aquifer.getCumulativeInflux(), 0.0);
}

TEST_F(EclipseFeaturesIntegrationTest, GroupControlWithVFP) {
    // Create group hierarchy
    GroupControlManager groups;
    groups.createGroup("PLATFORM", "FIELD");
    groups.assignWellToGroup("PROD1", "PLATFORM");
    groups.assignWellToGroup("PROD2", "PLATFORM");
    
    // Set group constraints
    GroupProdConstraints constraints;
    constraints.control_mode = GroupProdControlMode::ORAT;
    constraints.oil_rate_target = 0.1;
    groups.setProductionConstraints("PLATFORM", constraints);
    
    // Create VFP table
    VFPTableGenerator gen;
    gen.setTotalDepth(2000.0);
    auto vfp = gen.generateProductionTable(
        {1e6, 2e6}, {0.01, 0.1}, {0.0, 0.5}, {50.0}, {0.0});
    
    VFPTableManager vfp_manager;
    vfp_manager.addProductionTable(std::make_unique<VFPProdTable>(vfp));
    
    // Use VFP with group control
    double bhp = vfp_manager.calculateProductionBHP(1, 1.5e6, 0.05, 0.2, 75.0, 0.0);
    EXPECT_GT(bhp, 1.5e6);
}

// Main function
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
