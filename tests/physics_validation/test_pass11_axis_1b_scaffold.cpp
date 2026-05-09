/**
 * @file test_pass11_axis_1b_scaffold.cpp
 * @brief Pass-11 (axis 1b) scaffold regression gates. Verifies:
 *  1. cavity_geometry = THREE_DIMENSIONAL throws on
 *     RadialLagrangianSolver::setConfig with a clear "pass-12 work"
 *     diagnostic message.
 *  2. cavity_geometry = SPHERICAL is the byte-identical default that
 *     does not exercise the scaffold code.
 *  3. The Source3DBall factory throws on construction.
 */

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

#include "domain/explosion/NearFieldExplosion.hpp"
#include "domain/explosion/RadialLagrangian.hpp"
#include "domain/explosion/Source3DBall.hpp"
#include "domain/explosion/TillotsonEOS.hpp"

using FSRM::DamageEvolutionModel;
using FSRM::MieGruneisenEOS;
using FSRM::PressureDependentStrength;
using FSRM::RadialLagrangianSolver;
using FSRM::Source3DBall;
using FSRM::Source3DBallConfig;
using FSRM::TillotsonParameterSets;
using FSRM::UndergroundExplosionSource;
using FSRM::makeSource3DBall;

class Pass11Axis1bScaffoldTest : public ::testing::Test
{
};

// =========================================================================
// 1. cavity_geometry = THREE_DIMENSIONAL throws.
// =========================================================================
TEST_F(Pass11Axis1bScaffoldTest, ThreeDimensionalCavityGeometryThrows)
{
    UndergroundExplosionSource src;
    src.yield_kt = 1.0;
    src.depth = 500.0;
    src.location = {0.0, 0.0, -500.0};
    src.host_density = 2700.0;
    src.host_vp = 5500.0;
    src.host_vs = 3200.0;
    src.host_porosity = 0.005;
    src.overburden_stress = src.host_density * 9.81 * src.depth;

    RadialLagrangianSolver solver;
    solver.setSource(src);
    MieGruneisenEOS eos;
    eos.rho0 = src.host_density;
    eos.c0 = src.host_vp;
    solver.setEOS(eos);
    PressureDependentStrength strength;
    solver.setStrength(strength);
    DamageEvolutionModel damage;
    solver.setDamage(damage);

    RadialLagrangianSolver::Config cfg;
    cfg.cavity_geometry =
        RadialLagrangianSolver::CavityGeometry::THREE_DIMENSIONAL;
    cfg.tillotson_params = TillotsonParameterSets::granite();

    // Throw must contain a clear "pass-12 work" or "axis-1b" diagnostic.
    bool caught = false;
    std::string what;
    try {
        solver.setConfig(cfg);
    } catch (const std::runtime_error& e) {
        caught = true;
        what = e.what();
    }
    EXPECT_TRUE(caught)
        << "RadialLagrangianSolver must throw when cavity_geometry = "
           "THREE_DIMENSIONAL is selected (pass-11 scaffold).";
    EXPECT_NE(what.find("pass-12"), std::string::npos)
        << "Throw message should reference pass-12: " << what;
    EXPECT_NE(what.find("axis-1b"), std::string::npos)
        << "Throw message should reference axis-1b: " << what;
}

// =========================================================================
// 2. SPHERICAL (default) does not throw.
// =========================================================================
TEST_F(Pass11Axis1bScaffoldTest, SphericalCavityGeometryDoesNotThrow)
{
    UndergroundExplosionSource src;
    src.yield_kt = 1.0;
    src.depth = 500.0;
    src.location = {0.0, 0.0, -500.0};
    src.host_density = 2700.0;
    src.host_vp = 5500.0;
    src.host_vs = 3200.0;
    src.host_porosity = 0.005;
    src.overburden_stress = src.host_density * 9.81 * src.depth;

    RadialLagrangianSolver solver;
    solver.setSource(src);
    MieGruneisenEOS eos;
    eos.rho0 = src.host_density;
    eos.c0 = src.host_vp;
    solver.setEOS(eos);
    PressureDependentStrength strength;
    solver.setStrength(strength);
    DamageEvolutionModel damage;
    solver.setDamage(damage);

    RadialLagrangianSolver::Config cfg;
    cfg.cavity_geometry =
        RadialLagrangianSolver::CavityGeometry::SPHERICAL;
    cfg.tillotson_params = TillotsonParameterSets::granite();

    EXPECT_NO_THROW(solver.setConfig(cfg))
        << "SPHERICAL cavity_geometry is the pass-10 default and must not throw.";
}

// =========================================================================
// 3. Source3DBall factory throws.
// =========================================================================
TEST_F(Pass11Axis1bScaffoldTest, Source3DBallFactoryThrows)
{
    Source3DBallConfig cfg;
    cfg.outer_radius_m = 500.0;
    cfg.mesh_cell_size_m = 5.0;
    cfg.asymmetric_overburden = true;
    cfg.constitutive = Source3DBallConfig::Constitutive::DRUCKER_PRAGER_3D;

    bool caught = false;
    std::string what;
    try {
        auto ball = makeSource3DBall(cfg);
        (void)ball;
    } catch (const std::runtime_error& e) {
        caught = true;
        what = e.what();
    }
    EXPECT_TRUE(caught) << "makeSource3DBall must throw in pass-11.";
    EXPECT_NE(what.find("pass-12"), std::string::npos)
        << "Throw message should reference pass-12: " << what;
}
