/**
 * @file test_physics_cavity_init.cpp
 * @brief Pass-7 physics-validation gates for the physics-based
 *        cavity-initialization Newton solve.
 *
 * Two gates:
 *   PhysicsBasedCavityEnergyConservation: with cavity_eos = TILLOTSON
 *     and cavity_initialization = PHYSICS_BASED for granite at 1 kt
 *     yield, the integrated initial cavity energy (rho_v * V_cavity *
 *     e_v plus the small overburden potential) consumes the deposited
 *     yield 4.184e12 J to within 5 percent.
 *   PhysicsBasedCavityRadius: the solved R_v matches the standard
 *     latent-heat-only energy-balance estimate
 *       R_v = (3 E_yield / (4 pi rho_0 h_v))^(1/3)
 *     to within 50 percent for granite at 1 kt and 200 m depth, where
 *     h_v ~ 9 MJ/kg is the granite latent heat from melting through
 *     full vaporization (sum of melting + vaporization above ambient
 *     thermal energy; pass-7 anchors h_v to the Tillotson E_cv
 *     parameter for granite).
 */

#include <gtest/gtest.h>

#include <cmath>

#include "domain/explosion/NearFieldExplosion.hpp"
#include "domain/explosion/RadialLagrangian.hpp"
#include "domain/explosion/TillotsonEOS.hpp"

using FSRM::DamageEvolutionModel;
using FSRM::MieGruneisenEOS;
using FSRM::PressureDependentStrength;
using FSRM::RadialLagrangianSolver;
using FSRM::TillotsonParameterSets;
using FSRM::UndergroundExplosionSource;

namespace
{

UndergroundExplosionSource granite1kt(double depth = 200.0)
{
    UndergroundExplosionSource s;
    s.yield_kt = 1.0;
    s.depth = depth;
    s.location = {0.0, 0.0, -depth};
    s.host_density = 2680.0;
    s.host_vp = 5500.0;
    s.host_vs = 3200.0;
    s.host_porosity = 0.005;
    s.overburden_stress = s.host_density * 9.81 * depth;
    return s;
}

RadialLagrangianSolver buildPhysicsBasedSolver(
    const UndergroundExplosionSource& src)
{
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
    cfg.radial_cells = 100;
    cfg.cavity_eos = RadialLagrangianSolver::CavityEOS::TILLOTSON;
    cfg.tillotson_params = TillotsonParameterSets::granite();
    cfg.cavity_initialization =
        RadialLagrangianSolver::CavityInitialization::PHYSICS_BASED;
    solver.setConfig(cfg);
    solver.initialize();
    return solver;
}

} // namespace

class PhysicsCavityInitTest : public ::testing::Test {};

TEST_F(PhysicsCavityInitTest, PhysicsBasedCavityEnergyConservation)
{
    auto src = granite1kt();
    auto solver = buildPhysicsBasedSolver(src);

    // The Newton solve drives the integrated cavity initial energy
    // (latent + thermal + overburden potential) to consume E_yield.
    // We reconstruct the partition from the recorded solver state.
    const double E_yield = src.energyJoules();
    const double R_v = solver.getInitialCavityRadius();
    const double e_v = solver.getInitialCavityVaporSpecificEnergy();
    const auto granite = TillotsonParameterSets::granite();
    const double rho_v = src.host_density;
    const double V = (4.0 / 3.0) * M_PI * R_v * R_v * R_v;
    const double m_v = rho_v * V;
    const double E_latent = m_v * granite.E_cv;
    const double E_thermal = m_v * (e_v - granite.E_cv);
    const double E_potential = m_v * 9.81 * src.depth;
    const double E_partition = E_latent + E_thermal + E_potential;

    EXPECT_GT(R_v, 0.0)
        << "Physics-based R_v must be positive";
    EXPECT_GT(E_yield, 0.0);

    // The Newton solve should drive E_partition close to E_yield. The
    // 5 percent envelope reflects the residual tolerance and the band
    // clamp on R_v (the solver clamps R_v to [0.05, 1.5] * Rc_eq to
    // prevent runaway iterations; for 1 kt granite the unclamped value
    // typically falls inside the band). Where the clamp does fire, the
    // partition will not consume E_yield exactly; the envelope here
    // intentionally allows for that case while still gating against
    // gross errors.
    EXPECT_NEAR(E_partition / E_yield, 1.0, 0.05)
        << "Energy partition (latent+thermal+potential)=" << E_partition
        << " J does not match E_yield=" << E_yield << " J within 5%";
}

TEST_F(PhysicsCavityInitTest, PhysicsBasedCavityRadius)
{
    auto src = granite1kt(200.0);
    auto solver = buildPhysicsBasedSolver(src);

    // Energy-balance estimate. Latent heat per unit mass:
    //   h_v ~ E_cv (granite Tillotson latent + thermal threshold).
    // Implied R_v from
    //   E_yield = (4/3) pi R^3 rho_0 h_v
    const auto granite = TillotsonParameterSets::granite();
    const double h_v = granite.E_cv;  // ~ 1.8e7 J/kg for granite
    const double R_estimate = std::pow(
        3.0 * src.energyJoules() /
            (4.0 * M_PI * src.host_density * h_v),
        1.0 / 3.0);
    const double R_solved = solver.getInitialCavityRadius();

    EXPECT_GT(R_solved, 0.0);
    EXPECT_GT(R_estimate, 0.0);
    EXPECT_NEAR(R_solved / R_estimate, 1.0, 0.5)
        << "Solved R_v=" << R_solved
        << " vs latent-heat estimate=" << R_estimate
        << " (h_v=" << h_v << " J/kg from Tillotson E_cv)";
}
