/**
 * @file test_coulomb_stress_transfer.cpp
 * @brief Unit tests for CoulombStressTransfer (Hooke, fault projection, delta_CFS)
 */

#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>
#include <petscsys.h>

#include "domain/geomechanics/CoulombStressTransfer.hpp"
#include "domain/geomechanics/PyLithFault.hpp"

using namespace FSRM;

class CoulombStressTransferTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }

    // Wrappers for private members (friend access doesn't extend to TEST_F derived classes)
    static void callComputeStressFromStrain(const CoulombStressTransfer& cst,
                                             const double eps[6], double P, double sigma[6]) {
        cst.computeStressFromStrain(eps, P, sigma);
    }

    static std::vector<CoulombStressTransfer::VertexStress>& getCurrent(CoulombStressTransfer& cst) {
        return cst.current_;
    }

    int rank = 0;
};

TEST_F(CoulombStressTransferTest, ComputeStressFromStrainHookeLaw) {
    FaultCohesiveDyn fault;
    FaultVertex v;
    v.normal = {0.0, 0.0, 1.0};
    v.along_strike = {1.0, 0.0, 0.0};
    v.up_dip = {0.0, 1.0, 0.0};
    fault.setFaultVertices({v});

    CoulombStressTransfer cst(PETSC_COMM_WORLD);
    const double lambda = 40e9;
    const double mu = 30e9;
    const double biot = 0.8;
    PetscErrorCode ierr = cst.initialize(nullptr, &fault, lambda, mu, biot);
    ASSERT_EQ(ierr, 0);

    const double exx = 1e-4, eyy = -2e-5, ezz = 3e-5;
    const double exy = 2e-5, exz = -1e-5, eyz = 4e-5;
    const double eps[6] = {exx, eyy, ezz, exy, exz, eyz};
    const double P = 10e6;

    double sigma[6];
    callComputeStressFromStrain(cst, eps, P, sigma);

    const double eps_kk = exx + eyy + ezz;
    EXPECT_NEAR(sigma[0], lambda * eps_kk + 2.0 * mu * exx - biot * P, 1.0);
    EXPECT_NEAR(sigma[1], lambda * eps_kk + 2.0 * mu * eyy - biot * P, 1.0);
    EXPECT_NEAR(sigma[2], lambda * eps_kk + 2.0 * mu * ezz - biot * P, 1.0);
    EXPECT_NEAR(sigma[3], 2.0 * mu * exy, 1.0);
    EXPECT_NEAR(sigma[4], 2.0 * mu * exz, 1.0);
    EXPECT_NEAR(sigma[5], 2.0 * mu * eyz, 1.0);
}

TEST_F(CoulombStressTransferTest, ResolveStressOnFaultNormalAndShear) {
    FaultCohesiveDyn fault;
    FaultVertex v;
    v.normal = {0.0, 0.0, 1.0};
    v.along_strike = {1.0, 0.0, 0.0};
    v.up_dip = {0.0, 1.0, 0.0};
    fault.setFaultVertices({v});

    CoulombStressTransfer cst(PETSC_COMM_WORLD);
    ASSERT_EQ(cst.initialize(nullptr, &fault, 40e9, 30e9, 1.0), 0);

    auto& vs = getCurrent(cst)[0];
    // Pure shear sigma_xz = sigma_zx = 30 MPa gives traction on z-normal: t = (30e6,0,0)
    vs.sigma = {0.0, 0.0, 0.0, 0.0, 30e6, 0.0};
    vs.pressure = 0.0;

    ASSERT_EQ(cst.resolveStressOnFault(), 0);
    EXPECT_NEAR(vs.sigma_n_eff, 0.0, 1.0);
    EXPECT_NEAR(vs.tau, 30e6, 1.0);

    // Uniaxial compression along z on horizontal fault: sigma_n_eff = sigma_zz - alpha*P
    vs.sigma = {0.0, 0.0, -80e6, 0.0, 0.0, 0.0};
    vs.pressure = 10e6;
    ASSERT_EQ(cst.resolveStressOnFault(), 0);
    EXPECT_NEAR(vs.sigma_n_eff, -80e6 - 1.0 * 10e6, 1.0);
    EXPECT_NEAR(vs.tau, 0.0, 1.0);
}

TEST_F(CoulombStressTransferTest, DeltaCFSMatchesTauAndEffectiveNormalChange) {
    FaultCohesiveDyn fault;
    FaultVertex v;
    v.normal = {1.0, 0.0, 0.0};
    v.along_strike = {0.0, 1.0, 0.0};
    v.up_dip = {0.0, 0.0, 1.0};
    fault.setFaultVertices({v});

    CoulombStressTransfer cst(PETSC_COMM_WORLD);
    ASSERT_EQ(cst.initialize(nullptr, &fault, 40e9, 30e9, 1.0), 0);

    getCurrent(cst)[0].tau = 10e6;
    getCurrent(cst)[0].sigma_n_eff = 50e6;
    cst.storeInitialStress();

    getCurrent(cst)[0].tau = 15e6;
    getCurrent(cst)[0].sigma_n_eff = 40e6;

    const double mu_s = 0.6;
    ASSERT_EQ(cst.computeDeltaCFS(mu_s), 0);
    const double delta_tau = 5e6;
    const double delta_sig_eff = -10e6;
    const double expected = delta_tau - mu_s * delta_sig_eff;
    EXPECT_NEAR(getCurrent(cst)[0].delta_cfs, expected, 1.0);
    EXPECT_GT(getCurrent(cst)[0].delta_cfs, 0.0)
        << "Lower effective normal (injection-style) should raise delta_CFS when shear increases";
}

TEST_F(CoulombStressTransferTest, HydrostaticTotalStressZeroShearOnPlane) {
    FaultCohesiveDyn fault;
    FaultVertex v;
    v.normal = {1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)};
    v.along_strike = {0.0, 1.0, 0.0};
    v.up_dip = {-1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)};
    fault.setFaultVertices({v});

    CoulombStressTransfer cst(PETSC_COMM_WORLD);
    ASSERT_EQ(cst.initialize(nullptr, &fault, 40e9, 30e9, 1.0), 0);

    const double p = -60e6;
    auto& vs = getCurrent(cst)[0];
    vs.sigma = {p, p, p, 0.0, 0.0, 0.0};
    vs.pressure = 0.0;

    ASSERT_EQ(cst.resolveStressOnFault(), 0);
    EXPECT_NEAR(vs.tau, 0.0, 1.0);
}
