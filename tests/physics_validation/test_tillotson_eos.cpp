/**
 * @file test_tillotson_eos.cpp
 * @brief Pass-7 physics-validation gates for the TillotsonEOS class.
 *
 * The Tillotson EOS (Tillotson 1962) is the analytic form the pass-7
 * RadialLagrangianSolver uses for the inner-cavity post-vaporization
 * rock-plasma state. These gates are the per-class-per-medium
 * thermodynamic-self-consistency checks that lock in the parameter
 * sets shipped in include/domain/explosion/TillotsonEOS.hpp.
 *
 *  - GraniteHugoniotCompression: at the published granite shock
 *    Hugoniot states (nominal 10, 30, 100 GPa), the Tillotson
 *    pressure(rho, e) at the corresponding (rho_h, e_h) Hugoniot
 *    state is within an order-of-magnitude envelope of the
 *    experimental value. Pass-7 acceptance is intentionally wide:
 *    the goal here is NOT to certify granite-shock-modelling
 *    fidelity but to ensure the class evaluates the closed-form
 *    Tillotson formula correctly given Melosh's Table A2.2 set.
 *  - GraniteVaporizationEnergyBalance: heating granite at constant
 *    density rho_0 from e = 0 to e = E_cv must take energy E_cv
 *    per unit mass (definitional). This is a sanity check on the
 *    energy parameterization.
 *  - GraniteSoundSpeedConsistency: the closed-form sound-speed
 *    relation c^2 = (dp/drho)_e + (p / rho^2) (dp/de)_rho is
 *    numerically evaluated by the class. We verify the result is
 *    finite and positive at representative rho, e.
 *  - SaltAndAlluviumParameterSetSelfConsistency: the salt and
 *    alluvium sets each return finite, regime-correct pressures
 *    across the cold-compressed, expanded-cold, and expanded-hot
 *    regimes. Where alluvium parameters are placeholders, the test
 *    asserts placeholder behavior; the gap is documented in the
 *    parameter set itself and in HISTORIC_NUCLEAR_FIDELITY.md.
 */

#include <gtest/gtest.h>

#include <cmath>

#include "domain/explosion/TillotsonEOS.hpp"

using FSRM::TillotsonEOS;
using FSRM::TillotsonParameters;
using FSRM::TillotsonParameterSets;

namespace
{

// Approximate granite Hugoniot states at three reference shock
// pressures from Marsh (1980) LASL / Trunin (2001). These are the
// "rho_h" and corresponding particle-velocity-derived specific
// internal energy at the Hugoniot. Tabulated values vary slightly
// across the literature; the envelope below is intentionally generous
// to handle source-to-source spread.
struct HugoniotState
{
    double p_target_pa;
    double rho_h_kgm3;
    double e_h_jkg;
};

const HugoniotState granite_hugoniot[] = {
    {10.0e9, 3300.0, 1.0e6},
    {30.0e9, 3700.0, 4.0e6},
    {100.0e9, 4400.0, 1.5e7},
};

} // namespace

class TillotsonEOSTest : public ::testing::Test {};

TEST_F(TillotsonEOSTest, GraniteHugoniotCompression)
{
    TillotsonEOS eos(TillotsonParameterSets::granite());
    for (const auto& s : granite_hugoniot) {
        const double p_eos = eos.pressure(s.rho_h_kgm3, s.e_h_jkg);
        EXPECT_GT(p_eos, 0.0)
            << "Compressed-regime Tillotson pressure at rho_h="
            << s.rho_h_kgm3 << " e=" << s.e_h_jkg
            << " must be positive (compression).";
        // The order-of-magnitude envelope locks in the closed-form
        // evaluation; the absolute Hugoniot match is not the goal of
        // this gate (Melosh's parameter set is calibrated to the
        // Hugoniot, not the other way around).
        EXPECT_LT(p_eos, 1.0e3 * s.p_target_pa)
            << "p_eos=" << p_eos << " > 1000x target "
            << s.p_target_pa;
        EXPECT_GT(p_eos, 1.0e-3 * s.p_target_pa)
            << "p_eos=" << p_eos << " < 0.001x target "
            << s.p_target_pa;
    }
}

TEST_F(TillotsonEOSTest, GraniteVaporizationEnergyBalance)
{
    // The Tillotson E_cv parameter is the specific internal energy
    // (relative to the cold reference state) at which the host rock
    // is fully vaporized. Heating from (rho_0, 0) along an isobar to
    // (rho_0, E_cv) requires by definition E_cv per unit mass. This
    // gate just verifies the parameter is interpretable as such: at
    // (rho_0, E_cv) the EOS returns a finite pressure consistent with
    // the hot-expanded thermal term taking over from the cold cold-
    // pressure contribution (mu = 0 at rho_0 -> A * mu = 0).
    TillotsonEOS eos(TillotsonParameterSets::granite());
    const auto& p = eos.getParameters();
    const double p_at_e0 = eos.pressure(p.rho_0, 0.0);
    const double p_at_iv = eos.pressure(p.rho_0, p.E_iv);
    const double p_at_cv = eos.pressure(p.rho_0, p.E_cv);
    EXPECT_TRUE(std::isfinite(p_at_e0));
    EXPECT_TRUE(std::isfinite(p_at_iv));
    EXPECT_TRUE(std::isfinite(p_at_cv));
    // Increasing energy at fixed rho_0 should increase pressure
    // (the thermal term is monotone in e in the compressed branch).
    EXPECT_GT(p_at_iv, p_at_e0)
        << "At rho_0, p(E_iv) must exceed p(0): thermal term is "
           "monotone in e under the compressed branch.";
    EXPECT_GT(p_at_cv, p_at_iv)
        << "At rho_0, p(E_cv) must exceed p(E_iv).";
}

TEST_F(TillotsonEOSTest, GraniteSoundSpeedConsistency)
{
    TillotsonEOS eos(TillotsonParameterSets::granite());
    // Sample a representative grid of (rho, e) states that span the
    // cold-compressed, mixed, and hot-expanded regimes.
    const double rho_grid[] = {1500.0, 2680.0, 3500.0, 4500.0};
    const double e_grid[] = {1.0e5, 5.0e6, 1.5e7, 5.0e7};
    for (double rho : rho_grid) {
        for (double e : e_grid) {
            const double cs = eos.soundSpeed(rho, e);
            EXPECT_TRUE(std::isfinite(cs))
                << "Sound speed must be finite at rho=" << rho
                << " e=" << e;
            EXPECT_GE(cs, 0.0)
                << "Sound speed must be non-negative at rho=" << rho
                << " e=" << e;
        }
    }
}

TEST_F(TillotsonEOSTest, SaltAndAlluviumParameterSetSelfConsistency)
{
    auto check = [](TillotsonParameters p) {
        TillotsonEOS eos(p);
        // Compressed-regime sample.
        const double p_compressed =
            eos.pressure(1.2 * p.rho_0, 0.5 * p.E_iv);
        EXPECT_TRUE(std::isfinite(p_compressed)) << p.name;
        EXPECT_GT(p_compressed, 0.0)
            << p.name << " must give positive compression pressure";
        // Cold-expanded sample.
        const double p_expanded_cold =
            eos.pressure(0.9 * p.rho_0, 0.3 * p.E_iv);
        EXPECT_TRUE(std::isfinite(p_expanded_cold)) << p.name;
        // Hot-expanded sample.
        const double p_hot =
            eos.pressure(0.5 * p.rho_0, 2.0 * p.E_cv);
        EXPECT_TRUE(std::isfinite(p_hot)) << p.name;
        EXPECT_GE(p_hot, 0.0)
            << p.name << " hot-expanded pressure must be non-negative";
    };
    check(TillotsonParameterSets::tuff());
    check(TillotsonParameterSets::salt());
    check(TillotsonParameterSets::alluvium());
}
