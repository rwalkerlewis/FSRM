/**
 * @file test_viscoelastic_relaxation.cpp
 * @brief Physics validation: viscoelastic stress relaxation
 *
 * Tests the generalized Maxwell body by verifying that the mechanism weight
 * computation produces coefficients that approximate the target Q factor.
 *
 * For a single Maxwell element with relaxation time tau, the Q at the
 * relaxation frequency omega_0 = 1/tau should be approximately 1/(2*Y_0)
 * where Y_0 is the anelastic coefficient.
 *
 * For N mechanisms, the constant-Q approximation should hold over the
 * bandwidth [f_min, f_max]. The test verifies:
 *   1. computeMechanismWeights produces non-negative weights
 *   2. The resulting Q approximation matches the target within tolerance
 *   3. Relaxation times are log-spaced in the correct bandwidth
 */

#include <gtest/gtest.h>
#include <cmath>
#include "numerics/PetscFEViscoelastic.hpp"

using namespace FSRM;

class ViscoelasticRelaxationTest : public ::testing::Test {};

TEST_F(ViscoelasticRelaxationTest, MechanismWeights)
{
    // Test with 3 mechanisms, Q=100, bandwidth [0.1, 10] Hz
    const int N = 3;
    const double f_min = 0.1;
    const double f_max = 10.0;
    const double Q_target = 100.0;

    double tau[5] = {0};
    double delta_mu[5] = {0};

    PetscFEViscoelastic::computeMechanismWeights(N, f_min, f_max, Q_target, tau, delta_mu);

    // Verify relaxation times are positive and log-spaced
    for (int i = 0; i < N; ++i) {
        EXPECT_GT(tau[i], 0.0) << "Relaxation time " << i << " must be positive";
        EXPECT_TRUE(std::isfinite(tau[i])) << "Relaxation time " << i << " must be finite";
    }
    // Log spacing: tau should decrease (higher frequency = shorter time)
    for (int i = 1; i < N; ++i) {
        EXPECT_LT(tau[i], tau[i-1]) << "Relaxation times must decrease (log-spaced frequencies increase)";
    }

    // Verify weights are non-negative
    for (int i = 0; i < N; ++i) {
        EXPECT_GE(delta_mu[i], 0.0) << "Weight " << i << " must be non-negative";
    }

    // Verify the Q approximation: evaluate 1/Q at several frequencies
    // 1/Q(omega) = sum_i Y_i * omega*tau_i / (1 + (omega*tau_i)^2)
    const int n_test = 20;
    double max_Q_error = 0.0;
    for (int j = 0; j < n_test; ++j) {
        double t_frac = static_cast<double>(j) / (n_test - 1);
        double freq = std::exp(std::log(f_min) + t_frac * (std::log(f_max) - std::log(f_min)));
        double omega = 2.0 * M_PI * freq;

        double inv_Q = 0.0;
        for (int i = 0; i < N; ++i) {
            double wt = omega * tau[i];
            inv_Q += delta_mu[i] * wt / (1.0 + wt * wt);
        }

        double Q_approx = (inv_Q > 1e-15) ? 1.0 / inv_Q : 1e15;
        double rel_error = std::abs(Q_approx - Q_target) / Q_target;
        max_Q_error = std::max(max_Q_error, rel_error);
    }

    // Q approximation should be within 30% of target across the bandwidth
    // (3 mechanisms gives a reasonable but not perfect constant-Q approximation)
    EXPECT_LT(max_Q_error, 0.30)
        << "Q approximation error should be less than 30% across bandwidth. "
        << "Max relative error: " << max_Q_error;
}

TEST_F(ViscoelasticRelaxationTest, SingleMechanism)
{
    // For a single mechanism, Q at the relaxation frequency is approximately 1/(2*Y)
    const int N = 1;
    const double f_min = 0.5;
    const double f_max = 2.0;
    const double Q_target = 50.0;

    double tau[5] = {0};
    double delta_mu[5] = {0};

    PetscFEViscoelastic::computeMechanismWeights(N, f_min, f_max, Q_target, tau, delta_mu);

    ASSERT_GT(tau[0], 0.0);
    ASSERT_GT(delta_mu[0], 0.0);

    // At the relaxation frequency omega_0 = 1/tau:
    double omega_0 = 1.0 / tau[0];
    double wt = omega_0 * tau[0];  // = 1.0
    double inv_Q_at_omega0 = delta_mu[0] * wt / (1.0 + wt * wt);
    // = delta_mu[0] * 1.0 / 2.0 = delta_mu[0] / 2
    double Q_at_omega0 = 1.0 / inv_Q_at_omega0;

    // For a single mechanism, Q at the relaxation frequency should be close to target
    // Tolerance is wider because single mechanism cannot achieve constant Q
    EXPECT_NEAR(Q_at_omega0, Q_target, Q_target * 0.5)
        << "Single mechanism Q at relaxation frequency should approximate target";
}

TEST_F(ViscoelasticRelaxationTest, FiveMechanisms)
{
    // Five mechanisms should give a much better constant-Q approximation
    const int N = 5;
    const double f_min = 0.01;
    const double f_max = 100.0;
    const double Q_target = 200.0;

    double tau[5] = {0};
    double delta_mu[5] = {0};

    PetscFEViscoelastic::computeMechanismWeights(N, f_min, f_max, Q_target, tau, delta_mu);

    // Verify all weights are non-negative
    for (int i = 0; i < N; ++i) {
        EXPECT_GE(delta_mu[i], 0.0);
    }

    // Verify Q approximation accuracy (should be within 20% with 5 mechanisms)
    double max_Q_error = 0.0;
    const int n_test = 50;
    for (int j = 0; j < n_test; ++j) {
        double t_frac = static_cast<double>(j) / (n_test - 1);
        double freq = std::exp(std::log(f_min) + t_frac * (std::log(f_max) - std::log(f_min)));
        double omega = 2.0 * M_PI * freq;

        double inv_Q = 0.0;
        for (int i = 0; i < N; ++i) {
            double wt = omega * tau[i];
            inv_Q += delta_mu[i] * wt / (1.0 + wt * wt);
        }

        double Q_approx = (inv_Q > 1e-15) ? 1.0 / inv_Q : 1e15;
        double rel_error = std::abs(Q_approx - Q_target) / Q_target;
        max_Q_error = std::max(max_Q_error, rel_error);
    }

    EXPECT_LT(max_Q_error, 0.20)
        << "5-mechanism Q approximation error should be less than 20% across 4-decade bandwidth";
}
