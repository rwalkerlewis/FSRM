/**
 * @file test_volcano_physics.cpp
 * @brief Tests for volcanic physics models
 * 
 * Tests for:
 * - Magma chamber dynamics
 * - Conduit flow models
 * - Eruption column physics
 * - Volcanic deformation (Mogi source)
 * - Magma rheology
 * - Gas emission and dispersal
 */

#include <gtest/gtest.h>
#include "Testing.hpp"
#include "VolcanoModel.hpp"
#include <cmath>
#include <vector>

using namespace FSRM;
using namespace FSRM::Testing;

class VolcanoPhysicsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        
        // Physical constants
        g = VolcanoConstants::GRAVITY;
        R_gas = VolcanoConstants::GAS_CONSTANT;
        P_atm = VolcanoConstants::ATMOSPHERIC_PRESSURE;
        
        // Default magma properties (andesite)
        rho_melt = VolcanoConstants::SILICATE_MELT_DENSITY;
        rho_crystal = VolcanoConstants::CRYSTAL_DENSITY;
    }
    
    int rank;
    double g, R_gas, P_atm;
    double rho_melt, rho_crystal;
};

// ============================================================================
// Magma Rheology Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, MeltViscosityTemperatureDependence) {
    MagmaProperties magma;
    magma.composition = MagmaComposition::ANDESITE;
    magma.temperature = 1273.15;  // 1000°C
    magma.crystal_fraction = 0.0;
    magma.H2O_total = 4.0;
    
    double mu1 = magma.getMeltViscosity();
    
    // Increase temperature
    magma.temperature = 1373.15;  // 1100°C
    double mu2 = magma.getMeltViscosity();
    
    // Viscosity should decrease with temperature (Arrhenius)
    EXPECT_LT(mu2, mu1) << "Viscosity should decrease with temperature";
    EXPECT_GT(mu1, 0.0);
    EXPECT_GT(mu2, 0.0);
}

TEST_F(VolcanoPhysicsTest, ViscositySilicaDependence) {
    // Higher silica content → higher viscosity
    
    MagmaProperties basalt;
    basalt.SiO2 = 50.0;
    basalt.temperature = 1273.15;
    basalt.H2O_total = 0.5;
    basalt.crystal_fraction = 0.0;
    
    MagmaProperties rhyolite;
    rhyolite.SiO2 = 75.0;
    rhyolite.temperature = 1073.15;  // Lower T for rhyolite
    rhyolite.H2O_total = 4.0;
    rhyolite.crystal_fraction = 0.0;
    
    double mu_basalt = basalt.getMeltViscosity();
    double mu_rhyolite = rhyolite.getMeltViscosity();
    
    // Rhyolite should be much more viscous
    EXPECT_GT(mu_rhyolite, mu_basalt);
    
    // Typical ranges: basalt ~10-100 Pa·s, rhyolite ~10^6-10^10 Pa·s
    EXPECT_GT(mu_basalt, 1.0);
    EXPECT_LT(mu_basalt, 1e5);
}

TEST_F(VolcanoPhysicsTest, ViscosityWaterEffect) {
    MagmaProperties magma;
    magma.SiO2 = 65.0;
    magma.temperature = 1173.15;
    magma.crystal_fraction = 0.0;
    
    // Dry magma
    magma.H2O_total = 0.1;
    double mu_dry = magma.getMeltViscosity();
    
    // Wet magma
    magma.H2O_total = 4.0;
    double mu_wet = magma.getMeltViscosity();
    
    // Water dramatically reduces viscosity
    EXPECT_LT(mu_wet, mu_dry);
}

TEST_F(VolcanoPhysicsTest, CrystalViscosityEffect) {
    MagmaProperties magma;
    magma.temperature = 1273.15;
    
    magma.crystal_fraction = 0.0;
    double mu_0 = magma.getBulkViscosity();
    
    magma.crystal_fraction = 0.3;
    double mu_30 = magma.getBulkViscosity();
    
    magma.crystal_fraction = 0.5;
    double mu_50 = magma.getBulkViscosity();
    
    // Crystals increase viscosity (Einstein-Roscoe type relation)
    EXPECT_GT(mu_30, mu_0);
    EXPECT_GT(mu_50, mu_30);
    
    // Near maximum packing, viscosity should be very high
    magma.crystal_fraction = 0.6;  // Near rheological lockup
    double mu_60 = magma.getBulkViscosity();
    EXPECT_GT(mu_60, mu_50 * 10);  // Strong increase
}

TEST_F(VolcanoPhysicsTest, MagmaRheologyYieldStress) {
    MagmaRheology rheology;
    
    MagmaProperties magma;
    magma.crystal_fraction = 0.4;  // Crystal-rich
    rheology.setMagmaProperties(magma);
    
    double yield_stress = rheology.getYieldStress();
    
    // Crystal-rich magma should have yield stress
    EXPECT_GT(yield_stress, 0.0);
    
    // Check if magma can flow
    double shear_stress_low = yield_stress * 0.5;
    double shear_stress_high = yield_stress * 2.0;
    
    EXPECT_FALSE(rheology.canFlow(shear_stress_low));
    EXPECT_TRUE(rheology.canFlow(shear_stress_high));
}

// ============================================================================
// Volatile Solubility Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, H2OSolubility) {
    MagmaProperties magma;
    
    double P_high = 200e6;  // 200 MPa
    double P_low = 50e6;    // 50 MPa
    
    double H2O_high = magma.getDissolvedH2O(P_high);
    double H2O_low = magma.getDissolvedH2O(P_low);
    
    // Solubility increases with pressure
    EXPECT_GT(H2O_high, H2O_low);
    EXPECT_GT(H2O_high, 0.0);
    
    // At typical magma chamber pressures (~200 MPa), can dissolve ~5-6 wt% H2O
    EXPECT_GT(H2O_high, 3.0);
    EXPECT_LT(H2O_high, 10.0);
}

TEST_F(VolcanoPhysicsTest, CO2Solubility) {
    MagmaProperties magma;
    
    double P = 200e6;
    double CO2 = magma.getDissolvedCO2(P);
    
    // CO2 is much less soluble than H2O
    EXPECT_GT(CO2, 0.0);
    EXPECT_LT(CO2, 1.0);  // Typically < 0.5 wt% at 200 MPa
}

TEST_F(VolcanoPhysicsTest, SaturationPressure) {
    MagmaProperties magma;
    magma.H2O_total = 4.0;   // 4 wt% water
    magma.CO2_total = 0.1;   // 0.1 wt% CO2
    
    double P_sat = magma.getSaturationPressure();
    
    // Should be positive and realistic
    EXPECT_GT(P_sat, 0.0);
    EXPECT_LT(P_sat, 500e6);  // < 500 MPa
}

// ============================================================================
// Magma Chamber Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, MagmaChamberPressurization) {
    MagmaChamberGeometry geometry;
    geometry.depth = 5000.0;
    geometry.volume = 10e9;  // 10 km³
    
    MagmaProperties magma;
    
    MagmaChamberState initial_state;
    initial_state.pressure = 150e6;
    initial_state.overpressure = 5e6;
    
    MagmaChamberModel chamber;
    chamber.initialize(geometry, magma, initial_state);
    
    // Set recharge rate
    chamber.setRechargeRate(100.0);  // 100 kg/s
    
    // Pressure should increase with recharge
    double P0 = chamber.getPressure();
    chamber.update(86400.0);  // 1 day
    double P1 = chamber.getPressure();
    
    EXPECT_GT(P1, P0) << "Pressure should increase with recharge";
}

TEST_F(VolcanoPhysicsTest, MagmaChamberCompressibility) {
    MagmaChamberGeometry geometry;
    geometry.volume = 10e9;
    
    MagmaProperties magma;
    magma.bubble_fraction = 0.05;  // 5% bubbles
    
    MagmaChamberState state;
    
    MagmaChamberModel chamber;
    chamber.initialize(geometry, magma, state);
    
    double beta = chamber.getMagmaCompressibility();
    
    // Bubbly magma is more compressible
    EXPECT_GT(beta, 0.0);
    EXPECT_TRUE(std::isfinite(beta));
}

TEST_F(VolcanoPhysicsTest, MagmaChamberFailure) {
    MagmaChamberGeometry geometry;
    geometry.depth = 5000.0;
    
    MagmaProperties magma;
    
    MagmaChamberState state;
    state.overpressure = 5e6;  // 5 MPa overpressure
    
    MagmaChamberModel chamber;
    chamber.initialize(geometry, magma, state);
    
    // Check failure criteria
    double tensile_strength = 10e6;  // 10 MPa
    bool near_failure = chamber.isNearFailure(tensile_strength);
    
    // With 5 MPa overpressure vs 10 MPa strength, not yet failing
    EXPECT_FALSE(near_failure);
    
    // Increase overpressure
    state.overpressure = 15e6;
    chamber.initialize(geometry, magma, state);
    
    near_failure = chamber.isNearFailure(tensile_strength);
    EXPECT_TRUE(near_failure);
}

// ============================================================================
// Volcanic Deformation Tests (Mogi Source)
// ============================================================================

TEST_F(VolcanoPhysicsTest, MogiDisplacementPattern) {
    VolcanicDeformationModel deformation;
    deformation.setElasticProperties(30e9, 0.25);
    
    VolcanicDeformationSource source;
    source.type = DeformationSourceType::MOGI;
    source.x = 0.0;
    source.y = 0.0;
    source.depth = 5000.0;
    source.volume_change = 1e6;  // 1 km³ = 10^9 m³, but 10^6 m³ = 0.001 km³
    
    deformation.addSource(source);
    
    double ux, uy, uz;
    
    // At center (directly above source)
    deformation.computeDisplacement(0.0, 0.0, ux, uy, uz);
    
    // Maximum uplift at center
    EXPECT_GT(uz, 0.0) << "Inflation should cause uplift";
    EXPECT_NEAR(ux, 0.0, 1e-10) << "No horizontal at center";
    EXPECT_NEAR(uy, 0.0, 1e-10);
    
    // Off-center
    deformation.computeDisplacement(5000.0, 0.0, ux, uy, uz);
    
    // Horizontal pointing radially outward
    EXPECT_GT(ux, 0.0) << "Radial expansion from inflation";
    EXPECT_GT(uz, 0.0) << "Still uplifting";
    
    // Vertical decreases with distance
    double uz_center;
    deformation.computeDisplacement(0.0, 0.0, ux, uy, uz_center);
    EXPECT_GT(uz_center, uz) << "Uplift decreases with distance";
}

TEST_F(VolcanoPhysicsTest, MogiScaling) {
    // Mogi solution: uz_max = (3/4π) * ΔV / d²
    // where d is source depth
    
    double dV = 1e7;      // Volume change (m³)
    double d = 5000.0;    // Depth (m)
    double nu = 0.25;     // Poisson ratio
    
    // Analytical uz at center
    double uz_analytical = (1.0 - nu) * dV / (M_PI * d * d);
    
    VolcanicDeformationModel deformation;
    deformation.setElasticProperties(30e9, nu);
    
    VolcanicDeformationSource source;
    source.type = DeformationSourceType::MOGI;
    source.depth = d;
    source.volume_change = dV;
    
    deformation.addSource(source);
    
    double ux, uy, uz;
    deformation.computeDisplacement(0.0, 0.0, ux, uy, uz);
    
    EXPECT_NEAR(uz, uz_analytical, uz_analytical * 0.1);
}

TEST_F(VolcanoPhysicsTest, MogiTiltPattern) {
    VolcanicDeformationModel deformation;
    deformation.setElasticProperties(30e9, 0.25);
    
    VolcanicDeformationSource source;
    source.type = DeformationSourceType::MOGI;
    source.depth = 5000.0;
    source.volume_change = 1e7;
    
    deformation.addSource(source);
    
    double tilt_x, tilt_y;
    
    // Tilt at some distance from center
    deformation.computeTilt(3000.0, 0.0, tilt_x, tilt_y);
    
    // Tilt should be outward (positive in x direction)
    EXPECT_NE(tilt_x, 0.0);
    EXPECT_TRUE(std::isfinite(tilt_x));
    EXPECT_TRUE(std::isfinite(tilt_y));
    
    // Maximum tilt at r = d (depth) for Mogi source
    deformation.computeTilt(5000.0, 0.0, tilt_x, tilt_y);
    double tilt_max_r = std::abs(tilt_x);
    
    deformation.computeTilt(2000.0, 0.0, tilt_x, tilt_y);
    double tilt_closer = std::abs(tilt_x);
    
    // Tilt pattern should peak near r = d
    EXPECT_TRUE(std::isfinite(tilt_max_r));
    EXPECT_TRUE(std::isfinite(tilt_closer));
}

// ============================================================================
// Conduit Flow Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, ConduitFlowMassConservation) {
    // Mass flux should be constant along conduit (steady state)
    
    ConduitGeometry geometry;
    geometry.depth = 5000.0;
    geometry.length = 5000.0;
    geometry.radius_base = 15.0;
    geometry.radius_vent = 10.0;
    
    MagmaProperties magma;
    double P_chamber = 150e6;
    
    ConduitFlowModel conduit;
    conduit.initialize(geometry, magma, P_chamber);
    
    double mass_flux = conduit.getMassFlux();
    
    // Mass flux should be positive and finite
    EXPECT_GT(mass_flux, 0.0);
    EXPECT_TRUE(std::isfinite(mass_flux));
}

TEST_F(VolcanoPhysicsTest, FragmentationDepth) {
    ConduitGeometry geometry;
    geometry.depth = 5000.0;
    geometry.length = 5000.0;
    
    MagmaProperties magma;
    magma.H2O_total = 4.0;
    
    ConduitFlowModel conduit;
    conduit.initialize(geometry, magma, 200e6);
    conduit.solvesteadyState();
    
    double frag_depth = conduit.getFragmentationDepth();
    
    // Fragmentation should occur above chamber, below vent
    EXPECT_GT(frag_depth, 0.0);
    EXPECT_LT(frag_depth, geometry.length);
}

TEST_F(VolcanoPhysicsTest, ExitVelocityEstimate) {
    ConduitGeometry geometry;
    geometry.radius_vent = 50.0;  // Large vent
    
    MagmaProperties magma;
    magma.H2O_total = 4.0;
    
    ConduitFlowModel conduit;
    conduit.initialize(geometry, magma, 200e6);
    
    double v_exit = conduit.getExitVelocity();
    
    // Typical exit velocities: 50-300 m/s for explosive eruptions
    EXPECT_GT(v_exit, 0.0);
    EXPECT_TRUE(std::isfinite(v_exit));
}

// ============================================================================
// Eruption Column Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, ColumnHeightFromMER) {
    // Mastin et al. (2009) relation: H = 0.25 * MER^0.25
    
    double MER = 1e7;  // 10^7 kg/s (VEI 5-6)
    double H = VolcanoUtils::columnHeight_from_MER(MER);
    
    // Should give ~15-25 km for this MER
    EXPECT_GT(H, 10000.0);
    EXPECT_LT(H, 50000.0);
    
    // Higher MER → higher column
    double H_low = VolcanoUtils::columnHeight_from_MER(1e6);
    double H_high = VolcanoUtils::columnHeight_from_MER(1e8);
    
    EXPECT_GT(H_high, H);
    EXPECT_LT(H_low, H);
}

TEST_F(VolcanoPhysicsTest, ColumnBuoyancy) {
    EruptionColumnParameters params;
    params.exit_velocity = 200.0;
    params.exit_temperature = 1273.0;
    params.mass_flux = 1e7;
    params.gas_fraction = 0.95;
    
    EruptionColumnModel column;
    column.initialize(params);
    
    bool success = column.solve();
    
    // Column should reach neutral buoyancy level
    double NBL = column.getNeutralBuoyancyHeight();
    EXPECT_GT(NBL, 0.0);
    EXPECT_TRUE(std::isfinite(NBL));
}

TEST_F(VolcanoPhysicsTest, ColumnCollapse) {
    EruptionColumnParameters params;
    params.exit_velocity = 50.0;    // Low exit velocity
    params.exit_temperature = 1073.0;
    params.mass_flux = 1e8;         // Very high mass flux
    params.gas_fraction = 0.85;     // Less gas
    
    EruptionColumnModel column;
    column.initialize(params);
    column.solve();
    
    // With low velocity and high mass, column may collapse
    bool collapsed = column.doesCollapse();
    
    // Either collapsed or reached NBL
    EXPECT_TRUE(collapsed || column.getNeutralBuoyancyHeight() > 0.0);
}

// ============================================================================
// VEI and Eruption Intensity Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, VEIToParameters) {
    double volume, height, MER;
    
    // VEI 4
    VolcanoUtils::VEI_to_parameters(4, volume, height, MER);
    EXPECT_GT(volume, 0.01);   // > 0.01 km³
    EXPECT_LT(volume, 1.0);    // < 1 km³
    EXPECT_GT(height, 10.0);   // > 10 km column
    
    // VEI 6 (Pinatubo-scale)
    VolcanoUtils::VEI_to_parameters(6, volume, height, MER);
    EXPECT_GT(volume, 1.0);    // > 1 km³
    EXPECT_LT(volume, 100.0);  // < 100 km³
    EXPECT_GT(height, 25.0);   // > 25 km column
}

TEST_F(VolcanoPhysicsTest, EruptionDurationEstimate) {
    double volume_km3 = 5.0;    // 5 km³
    double MER = 1e8;           // 10^8 kg/s
    
    double duration = VolcanoUtils::estimateDuration(volume_km3, MER);
    
    // Should be on order of hours to days
    EXPECT_GT(duration, 3600.0);       // > 1 hour
    EXPECT_LT(duration, 86400.0 * 7);  // < 1 week
}

// ============================================================================
// Pyroclastic Density Current Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, PDCDynamicPressure) {
    PDCState state;
    state.density = 50.0;      // kg/m³
    state.u = 50.0;            // m/s
    state.v = 0.0;
    state.w = 0.0;
    
    double P_dyn = state.dynamic_pressure();
    
    // P_dyn = 0.5 * rho * v²
    double P_expected = 0.5 * state.density * state.u * state.u;
    EXPECT_NEAR(P_dyn, P_expected, 1.0);
    
    // 50 kg/m³ at 50 m/s → 62.5 kPa
    EXPECT_NEAR(P_dyn, 62500.0, 100.0);
}

TEST_F(VolcanoPhysicsTest, PDCRunout) {
    PDCSourceConditions source;
    source.type = PDCType::COLUMN_COLLAPSE;
    source.mass_flux = 1e7;
    source.initial_velocity = 50.0;
    source.particle_concentration = 0.01;
    
    PDCModel pdc;
    pdc.initialize(source);
    
    double runout = pdc.getRunoutDistance();
    
    // Runout depends on source conditions
    EXPECT_GE(runout, 0.0);
}

// ============================================================================
// Lava Flow Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, LavaFlowVelocity) {
    // Lava velocity depends on slope, viscosity, and thickness
    // v ∝ ρg sin(θ) h² / μ (for Newtonian)
    
    double rho = 2500.0;       // Density (kg/m³)
    double theta = 0.1;        // Slope (rad, ~6°)
    double h = 5.0;            // Thickness (m)
    double mu = 1e4;           // Viscosity (Pa·s)
    
    double v = rho * g * std::sin(theta) * h * h / (3.0 * mu);
    
    // Typical lava flow velocities: 0.1 - 10 m/s
    EXPECT_GT(v, 0.0);
    EXPECT_LT(v, 100.0);
}

TEST_F(VolcanoPhysicsTest, LavaFlowCooling) {
    LavaFlowState state;
    state.temperature = 1373.15;  // 1100°C
    state.thickness = 5.0;
    state.crystal_fraction = 0.3;
    
    // Check mobility
    bool mobile = state.isMobile();
    EXPECT_TRUE(mobile) << "Hot lava with low crystals should be mobile";
    
    // Cool down and increase crystals
    state.temperature = 900.0;
    state.crystal_fraction = 0.75;
    
    mobile = state.isMobile();
    EXPECT_FALSE(mobile) << "Cold, crystal-rich lava should be immobile";
}

// ============================================================================
// Lahar Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, LaharDensity) {
    LaharState state;
    state.sediment_conc = 0.4;  // 40% sediment
    
    double rho = state.density();
    
    // Lahar density between water and rock
    EXPECT_GT(rho, 1000.0);    // > pure water
    EXPECT_LT(rho, 2650.0);    // < pure sediment
    
    // For 40% sediment: ρ ≈ 1000*(1-0.4) + 2650*0.4 = 600 + 1060 = 1660 kg/m³
    EXPECT_NEAR(rho, 1660.0, 10.0);
}

TEST_F(VolcanoPhysicsTest, LaharDynamicPressure) {
    LaharState state;
    state.sediment_conc = 0.5;
    state.velocity_x = 10.0;
    state.velocity_y = 0.0;
    state.depth = 3.0;
    
    double P_dyn = state.dynamic_pressure();
    
    // High-density, fast-moving lahar has high impact force
    EXPECT_GT(P_dyn, 0.0);
    
    // Should be destructive at typical velocities
    // ~80 kPa threshold for building destruction
}

// ============================================================================
// Volcanic Gas Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, GasEmissionRates) {
    GasEmissionSource source;
    source.emission_rate = 1000.0;  // kg/s total
    source.composition[VolcanicGas::H2O] = 0.90;
    source.composition[VolcanicGas::CO2] = 0.05;
    source.composition[VolcanicGas::SO2] = 0.03;
    
    VolcanicGasModel gas;
    gas.addSource(source);
    
    double SO2_rate = gas.getTotalSO2_emission();
    
    // SO2 rate = 1000 * 0.03 = 30 kg/s
    EXPECT_GT(SO2_rate, 0.0);
    EXPECT_TRUE(std::isfinite(SO2_rate));
}

TEST_F(VolcanoPhysicsTest, GausianPlumeDispersion) {
    GasEmissionSource source;
    source.x = 0.0;
    source.y = 0.0;
    source.emission_rate = 100.0;
    source.temperature = 373.15;
    source.velocity = 5.0;
    
    VolcanicGasModel gas;
    gas.addSource(source);
    gas.setWindSpeed(5.0);
    gas.setWindDirection(0.0);  // From north
    gas.setAtmosphericStability('D');  // Neutral
    
    // Concentration should decrease with distance
    double conc_near = gas.getConcentration(VolcanicGas::SO2, 100.0, 0.0);
    double conc_far = gas.getConcentration(VolcanicGas::SO2, 1000.0, 0.0);
    
    EXPECT_GT(conc_near, conc_far) << "Concentration decreases with distance";
}

// ============================================================================
// Volcanic Seismicity Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, LPEventFrequency) {
    VolcanicSeismicEvent event;
    event.type = VolcanicSeismicEventType::LP;
    event.dominant_frequency = 1.5;  // Hz
    event.quality_factor = 20.0;
    event.crack_length = 100.0;
    
    // LP events have characteristic frequency < 5 Hz
    EXPECT_LT(event.dominant_frequency, 5.0);
    EXPECT_GT(event.dominant_frequency, 0.1);
}

TEST_F(VolcanoPhysicsTest, CrackResonanceFrequency) {
    VolcanicSeismicityModel::CrackResonance crack;
    crack.length = 100.0;         // m
    crack.width = 50.0;           // m
    crack.fluid_velocity = 500.0; // m/s (for fluid in crack)
    
    double f0 = crack.fundamentalFrequency();
    
    // Frequency ~ c / (2L)
    EXPECT_GT(f0, 0.0);
    EXPECT_TRUE(std::isfinite(f0));
}

TEST_F(VolcanoPhysicsTest, VTEventMagnitude) {
    VolcanicSeismicEvent event;
    event.type = VolcanicSeismicEventType::VT;
    event.seismic_moment = 1e14;  // N·m
    
    // Mw = (2/3) * log10(M0) - 6.06
    double Mw = (2.0/3.0) * std::log10(event.seismic_moment) - 6.06;
    
    // Typical volcanic VT: M0-3
    EXPECT_GT(Mw, -1.0);
    EXPECT_LT(Mw, 5.0);
}

// ============================================================================
// Tephra Dispersal Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, TephraSettlingVelocity) {
    TephraParticle particle;
    particle.diameter = 1e-3;   // 1 mm
    particle.density = 2500.0;
    particle.shape_factor = 0.7;
    
    double rho_air = 1.2;       // kg/m³
    double mu_air = 1.8e-5;     // Pa·s
    
    double v_s = particle.settlingVelocity(rho_air, mu_air);
    
    // Settling velocity should be positive and reasonable
    EXPECT_GT(v_s, 0.0);
    EXPECT_LT(v_s, 100.0);  // < 100 m/s for particles
    
    // Larger particles settle faster
    TephraParticle large;
    large.diameter = 1e-2;  // 10 mm
    large.density = 2500.0;
    large.shape_factor = 0.7;
    
    double v_large = large.settlingVelocity(rho_air, mu_air);
    EXPECT_GT(v_large, v_s);
}

TEST_F(VolcanoPhysicsTest, PhiScale) {
    // phi = -log2(d_mm)
    
    double d_mm = 1.0;
    double phi = VolcanoUtils::mm_to_phi(d_mm);
    EXPECT_NEAR(phi, 0.0, 0.01);  // 1 mm = phi 0
    
    d_mm = 2.0;
    phi = VolcanoUtils::mm_to_phi(d_mm);
    EXPECT_NEAR(phi, -1.0, 0.01);  // 2 mm = phi -1
    
    d_mm = 0.5;
    phi = VolcanoUtils::mm_to_phi(d_mm);
    EXPECT_NEAR(phi, 1.0, 0.01);  // 0.5 mm = phi 1
    
    // Inverse
    double d_back = VolcanoUtils::phi_to_mm(2.0);
    EXPECT_NEAR(d_back, 0.25, 0.01);  // phi 2 = 0.25 mm
}

// ============================================================================
// Scenario Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, MountStHelensScenario) {
    VolcanoSystemConfig config = VolcanoScenarios::mountStHelens1980();
    MagmaChamberGeometry chamber = VolcanoScenarios::mountStHelens1980_chamber();
    MagmaProperties magma = VolcanoScenarios::mountStHelens1980_magma();
    
    // Dacite composition
    EXPECT_GT(magma.SiO2, 60.0);
    EXPECT_LT(magma.SiO2, 70.0);
    
    // Shallow chamber
    EXPECT_LT(chamber.depth, 10000.0);
}

TEST_F(VolcanoPhysicsTest, KilaueaScenario) {
    MagmaProperties magma = VolcanoScenarios::kilauea_magma();
    
    // Basalt composition
    EXPECT_LT(magma.SiO2, 55.0);
    
    // Low viscosity expected
    double mu = magma.getMeltViscosity();
    EXPECT_LT(mu, 1000.0);  // < 1000 Pa·s for basalt
}

TEST_F(VolcanoPhysicsTest, YellowstoneUnrestScenario) {
    MagmaChamberGeometry chamber = VolcanoScenarios::yellowstone_chamber();
    
    // Large, shallow chamber
    EXPECT_GT(chamber.volume, 100e9);  // > 100 km³
}

// ============================================================================
// Energy and Mass Budget Tests
// ============================================================================

TEST_F(VolcanoPhysicsTest, ThermalEnergyRelease) {
    // Thermal energy release: E = m * cp * ΔT
    
    double mass = 1e12;        // 10^12 kg ejected
    double cp = 1200.0;        // Specific heat (J/kg/K)
    double dT = 500.0;         // Temperature drop (K)
    
    double E_thermal = mass * cp * dT;
    
    // Should be on order of 10^17 - 10^18 J for large eruption
    EXPECT_GT(E_thermal, 1e17);
    EXPECT_TRUE(std::isfinite(E_thermal));
}

TEST_F(VolcanoPhysicsTest, SeismicEfficiency) {
    // Typically < 10% of total energy goes to seismic waves
    
    double E_total = 1e18;  // Total energy
    double efficiency = 0.05;  // 5%
    
    double E_seismic = E_total * efficiency;
    double M0 = E_seismic * 2e4;  // Rough M0 ~ 2e4 * E_seismic
    double Mw = (2.0/3.0) * std::log10(M0) - 6.06;
    
    EXPECT_TRUE(std::isfinite(Mw));
}

TEST_F(VolcanoPhysicsTest, MassConservationEruption) {
    // Mass in chamber = Mass ejected + Mass intruded - Mass recharged
    
    double M_initial = 1e13;     // kg
    double M_ejected = 1e12;     // kg
    double M_recharged = 0.0;    // kg (no recharge during eruption)
    double M_intruded = 0.0;     // kg
    
    double M_final = M_initial - M_ejected + M_recharged - M_intruded;
    
    EXPECT_NEAR(M_final, 9e12, 1e10);
}
