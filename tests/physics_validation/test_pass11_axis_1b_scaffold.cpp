/**
 * @file test_pass11_axis_1b_scaffold.cpp
 * @brief Pass-11 / pass-13a / pass-13b axis-1b regression gates.
 *
 * Pass-13b updates the 13a throw gate to a positive delegation gate:
 *  1. cavity_geometry = THREE_DIMENSIONAL no longer throws on
 *     setConfig (pass-13a kept the throw with a "pass-13b/c follow-on"
 *     message; pass-13b lands the host-side delegation into
 *     Source3DBallImpl).
 *  2. After successful setConfig under THREE_DIMENSIONAL the host's
 *     name() reports the pass-13b 3D impl tag.
 *  3. SPHERICAL (default) does not throw and reports the legacy 1D
 *     name. Backward-compat for the 32 historic-event integration
 *     tests.
 *  4. makeSource3DBall() returns a Source3DBallImpl whose name()
 *     reports the pass-13b physics tag (pass-13a foundation tag
 *     bumped on this branch).
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

namespace
{

UndergroundExplosionSource makeTestSource()
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
    return src;
}

}  // namespace

// =========================================================================
// 1. Pass-13b: cavity_geometry = THREE_DIMENSIONAL delegates to
//    Source3DBallImpl. The pass-11/13a throw is gone.
// =========================================================================
TEST_F(Pass11Axis1bScaffoldTest, ThreeDimensionalCavityGeometryDelegates)
{
    auto src = makeTestSource();
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
    // Foundation no-op mode: empty mesh path so the 3D impl
    // initialises without touching the disk. This is the host-
    // delegation gate; the full mesh-loaded path is exercised by the
    // SphericalCellLevelEquivalenceVs1D physics-validation gate.
    cfg.source_3d_ball.mesh_path = "";

    EXPECT_NO_THROW(solver.setConfig(cfg))
        << "Pass-13c: cavity_geometry = THREE_DIMENSIONAL must no "
           "longer throw. The host constructs a Source3DBallImpl and "
           "delegates step / moment-tensor accessors to it.";

    const std::string nm = solver.name();
    EXPECT_NE(nm.find("Source3DBallImpl"), std::string::npos)
        << "name() should report the 3D impl tag when delegation is "
           "active: " << nm;
    EXPECT_NE(nm.find("pass13c"), std::string::npos)
        << "name() should reference pass-13c: " << nm;

    // Foundation no-op mode (empty mesh_path) leaves M / Mdot at zero
    // because no surface-integral set exists. The full integral is
    // exercised by Source3DBallSurfaceIntegral.* unit gates and the
    // Salmon3DWithOverburdenAtMPI4 integration gate.
    std::array<double, 6> M{};
    std::array<double, 6> Mdot{};
    solver.getMomentTensor(M);
    solver.getMomentRateTensor(Mdot);
    for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(M[i], 0.0)
            << "Foundation no-op mode (empty mesh_path) must leave M=0; "
               "M[" << i << "] = " << M[i];
        EXPECT_EQ(Mdot[i], 0.0);
    }
}

// =========================================================================
// 2. SPHERICAL (default) does not throw and reports the legacy 1D name.
// =========================================================================
TEST_F(Pass11Axis1bScaffoldTest, SphericalCavityGeometryDoesNotThrow)
{
    auto src = makeTestSource();
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
        << "SPHERICAL cavity_geometry is the pass-10 default and must "
           "not throw.";
    const std::string nm = solver.name();
    EXPECT_EQ(nm.find("Source3DBallImpl"), std::string::npos)
        << "Under SPHERICAL the host should report the legacy 1D "
           "name without the 3D impl tag: " << nm;
}

// =========================================================================
// 3. Pass-13b: makeSource3DBall() returns a Source3DBallImpl with the
//    pass-13b physics tag.
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
        << "Pass-13c: makeSource3DBall must construct without "
           "throwing. The pass-11 throw was replaced in pass-13a; "
           "pass-13c bumps the implementation tag to the validation "
           "milestone.";
    ASSERT_NE(ball, nullptr);
    const std::string nm = ball->name();
    EXPECT_NE(nm.find("pass13c"), std::string::npos)
        << "Source3DBallImpl::name() should tag the pass-13c "
           "validation implementation: " << nm;
}
