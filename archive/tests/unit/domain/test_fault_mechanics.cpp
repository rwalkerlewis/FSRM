/**
 * @file test_fault_mechanics.cpp
 * @brief Unit tests for fault friction and Coulomb failure (FaultModel.hpp / FaultMechanics.cpp)
 */

#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <petscsys.h>

#include "domain/geomechanics/FaultModel.hpp"

using namespace FSRM;

class FaultMechanicsTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

// Coulomb slip-weakening: linear mu from mu_s at slip=0 to mu_d at slip >= Dc
TEST_F(FaultMechanicsTest, CoulombSlipWeakeningLinearRamp) {
    CoulombFriction model;
    model.params.static_friction = 0.6;
    model.params.dynamic_friction = 0.4;
    model.params.slip_weakening_Dc = 0.4;
    model.params.cohesion = 0.0;

    EXPECT_NEAR(model.getFriction(0.0, 0.0, 50e6), 0.6, 1e-12);

    const double slip_mid = 0.2;
    const double mu_mid = 0.6 - (0.6 - 0.4) * (slip_mid / 0.4);
    EXPECT_NEAR(model.getFriction(0.0, slip_mid, 50e6), mu_mid, 1e-12);

    EXPECT_NEAR(model.getFriction(0.0, 0.4, 50e6), 0.4, 1e-12);
    EXPECT_NEAR(model.getFriction(0.0, 1.0, 50e6), 0.4, 1e-12);
}

// Steady-state rate-state: mu_ss = f0 + (a - b) * ln(V / V0)
TEST_F(FaultMechanicsTest, RateStateSteadyStateFriction) {
    RateStateFriction rs(RateStateFriction::EvolutionLaw::AGING);
    rs.params.f0 = 0.6;
    rs.params.a = 0.01;
    rs.params.b = 0.015;
    rs.params.V0 = 1e-6;
    rs.params.Dc = 0.4;

    const double V = 1e-3;
    const double mu_ss = rs.getSteadyStateFriction(V);
    const double expected = rs.params.f0 +
                            (rs.params.a - rs.params.b) * std::log(V / rs.params.V0);
    EXPECT_NEAR(mu_ss, expected, 1e-10);

    // At steady state theta = Dc/V, aging-law friction should match mu_ss
    const double theta_ss = rs.params.Dc / V;
    const double mu_from_state =
        rs.getFriction(V, theta_ss, 50e6);
    EXPECT_NEAR(mu_from_state, mu_ss, 1e-10);
}

// CFF = tau - tau_strength; positive => shear exceeds strength (failure / slipping threshold)
TEST_F(FaultMechanicsTest, CoulombFailureFunctionSignConvention) {
    SeismicFaultModel fault;
    FaultGeometry geom;
    geom.strike = 0.0;
    geom.dip = M_PI / 2.0;   // vertical plane, normal mostly horizontal
    geom.length = 1000.0;
    geom.width = 500.0;
    fault.setGeometry(geom);

    auto coulomb = std::make_unique<CoulombFriction>();
    coulomb->params.static_friction = 0.6;
    coulomb->params.dynamic_friction = 0.4;
    coulomb->params.cohesion = 0.0;
    coulomb->params.slip_weakening_Dc = 1.0;
    fault.setFrictionModel(std::move(coulomb));

    // Convention in FaultGeometry::resolveStress: compression positive for sigma_n
    const double sxx = 100e6, syy = 100e6, szz = 120e6;
    const double sxy = 0.0, sxz = 40e6, syz = 0.0;
    const double pp = 0.0;

    FaultStressState stable = fault.computeStressState(sxx, syy, szz, sxy, sxz, syz, pp);
    EXPECT_LT(stable.CFF, 0.0)
        << "Shear stress should sit below Coulomb strength (stable, CFF < 0)";

    const double sxz_fail = 75e6;
    FaultStressState failing =
        fault.computeStressState(sxx, syy, szz, sxy, sxz_fail, syz, pp);
    EXPECT_GE(failing.CFF, 0.0)
        << "Increased shear should meet or exceed strength (CFF >= 0 promotes failure)";
}

TEST_F(FaultMechanicsTest, FaultStressStateStabilityHelpers) {
    FaultStressState s;
    s.CFF = -1e3;
    s.is_slipping = false;
    EXPECT_TRUE(s.isStable());
    EXPECT_FALSE(s.isCritical());

    s.CFF = 0.0;
    s.is_slipping = false;
    EXPECT_FALSE(s.isStable());
    EXPECT_TRUE(s.isCritical());
}
