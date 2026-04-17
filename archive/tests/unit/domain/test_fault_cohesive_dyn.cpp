/**
 * @file test_fault_cohesive_dyn.cpp
 * @brief Unit tests for PyLith-style dynamic cohesive fault friction (PyLithFault.hpp)
 */

#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>
#include <petscsys.h>

#include "domain/geomechanics/PyLithFault.hpp"

using namespace FSRM;

class FaultCohesiveDynTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    int rank = 0;
};

TEST_F(FaultCohesiveDynTest, SlipWeakeningLinearDecreaseWithSlipState) {
    SlipWeakeningFriction law;
    law.setStaticCoefficient(0.6);
    law.setDynamicCoefficient(0.4);
    law.setCriticalSlipDistance(0.4);

    DynamicFaultState st;
    st.slip_rate = 0.0;
    st.effective_normal = 50e6;
    st.shear_stress = 0.0;

    st.state_variable = 0.0;
    EXPECT_NEAR(law.computeFriction(st), 0.6, 1e-12);

    st.state_variable = 0.2;
    const double mu_mid = 0.6 - (0.6 - 0.4) * (0.2 / 0.4);
    EXPECT_NEAR(law.computeFriction(st), mu_mid, 1e-12);

    st.state_variable = 0.4;
    EXPECT_NEAR(law.computeFriction(st), 0.4, 1e-12);
}

TEST_F(FaultCohesiveDynTest, RateStateAgingStateEvolutionRate) {
    RateStateFrictionAging law;
    law.setCriticalSlipDistance(0.4);
    law.setReferenceSlipRate(1e-6);

    DynamicFaultState st;
    st.effective_normal = 40e6;
    st.shear_stress = 20e6;
    st.state_variable = 0.5;
    st.slip_rate = 1e-3;

    const double V = std::max(std::abs(st.slip_rate), 1e-12);
    const double theta = std::max(st.state_variable, 1e-20);
    const double expected = 1.0 - V * theta / 0.4;
    EXPECT_NEAR(law.computeStateRate(st), expected, 1e-12);
}

TEST_F(FaultCohesiveDynTest, RateStateAgingSteadyStateThetaMatchesFriction) {
    RateStateFrictionAging law;
    law.setDirectEffect(0.01);
    law.setEvolutionEffect(0.015);
    law.setCriticalSlipDistance(0.4);
    law.setReferenceFriction(0.6);
    law.setReferenceSlipRate(1e-6);

    const double V = 1e-4;
    DynamicFaultState st;
    st.slip_rate = V;
    st.state_variable = 0.4 / V;  // theta_ss = Dc / V => dtheta/dt = 0
    st.effective_normal = 30e6;
    st.shear_stress = 15e6;

    EXPECT_NEAR(law.computeStateRate(st), 0.0, 1e-9);

    const double mu_ss = law.computeFriction(st);
    const double mu_rs =
        0.6 + (0.01 - 0.015) * std::log(V / 1e-6);
    EXPECT_NEAR(mu_ss, mu_rs, 1e-9);
}

TEST_F(FaultCohesiveDynTest, DynamicFaultStateLockedVersusSlippingFromResidual) {
    FaultCohesiveDyn fault;
    FaultVertex v;
    v.vertex_id = 0;
    v.vertex_negative = 0;
    v.vertex_positive = 1;
    v.vertex_lagrange = -1;
    v.coords = {0.0, 0.0, 0.0};
    v.normal = {0.0, 0.0, 1.0};
    v.along_strike = {1.0, 0.0, 0.0};
    v.up_dip = {0.0, 1.0, 0.0};
    v.area = 1.0;

    fault.setFaultVertices({v});
    // Static friction gives a fixed strength mu*sigma'_n; slip-weakening uses state_variable
    // as cumulative slip and initialize() seeds a large placeholder theta — use Coulomb-static
    // for a clean stick/slip threshold test.
    auto stat = std::make_unique<StaticFriction>();
    stat->setFrictionCoefficient(0.6);
    fault.setFrictionModel(std::move(stat));
    // Compression-positive effective normal: traction_normal = -50 MPa -> sigma'_n = 50 MPa
    fault.setUniformInitialTraction(20e6, 0.0, -50e6);
    fault.initialize();

    std::vector<double> sol(6, 0.0);  // two vertices, 3 DOFs each

    fault.computeResidual(sol.data(), nullptr);
    EXPECT_TRUE(fault.getState(0).is_locked)
        << "tau=20 MPa below 0.6*50 MPa strength should remain locked";

    fault.setUniformInitialTraction(35e6, 0.0, -50e6);
    fault.initialize();
    fault.computeResidual(sol.data(), nullptr);
    EXPECT_FALSE(fault.getState(0).is_locked)
        << "tau=35 MPa above 0.6*50 MPa should unlock (slipping branch)";
}
