/**
 * @file test_pass11_axis_1b_scaffold.cpp
 * @brief Pass-11 / pass-13a regression gates for the axis-1b scaffold
 *        and foundation slice. Verifies:
 *
 *  1. cavity_geometry = THREE_DIMENSIONAL throws on
 *     RadialLagrangianSolver::setConfig. Pass-13a foundation kept
 *     this throw because the host-side delegation from
 *     RadialLagrangianSolver into Source3DBallImpl is pass-13b/c work;
 *     the throw message now references pass-13 explicitly.
 *  2. cavity_geometry = SPHERICAL is the byte-identical default that
 *     does not exercise the scaffold code.
 *  3. Pass-13a foundation: makeSource3DBall() no longer throws and
 *     returns a Source3DBallImpl instance whose name() reports the
 *     foundation-skeleton tag. Replaces the pass-11
 *     Source3DBallFactoryThrows assertion.
 */

#include <gtest/gtest.h>

#include <memory>
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
           "THREE_DIMENSIONAL is selected (pass-13b/c host delegation "
           "still pending; foundation slice ships only the leaf "
           "Source3DBallImpl).";
    EXPECT_NE(what.find("pass-13"), std::string::npos)
        << "Throw message should reference pass-13: " << what;
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
// 3. Pass-13a foundation: makeSource3DBall() returns a Source3DBallImpl.
//    Replaces the pass-11 Source3DBallFactoryThrows assertion.
// =========================================================================
TEST_F(Pass11Axis1bScaffoldTest, Source3DBallFactoryReturnsImpl)
{
    Source3DBallConfig cfg;
    cfg.outer_radius_m = 500.0;
    cfg.mesh_cell_size_m = 5.0;
    cfg.asymmetric_overburden = true;
    cfg.constitutive = Source3DBallConfig::Constitutive::DRUCKER_PRAGER_3D;

    std::unique_ptr<Source3DBall> ball;
    EXPECT_NO_THROW({ ball = makeSource3DBall(cfg); })
        << "Pass-13a foundation: makeSource3DBall must construct without "
           "throwing. The pass-11 throw has been replaced with a "
           "Source3DBallImpl instance.";
    ASSERT_NE(ball, nullptr);
    const std::string nm = ball->name();
    EXPECT_NE(nm.find("pass-13a"), std::string::npos)
        << "Source3DBallImpl::name() should tag the pass-13a foundation "
           "skeleton: " << nm;
}
