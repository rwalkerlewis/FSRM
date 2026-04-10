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

/**
 * @brief Verify FEM gradient recovery on a small 3D tet mesh with a known
 *        linear displacement field u_x = a*x, u_y = 0, u_z = 0.
 *
 * For a linear displacement field, the strain is:
 *   eps_xx = a, all others = 0.
 * The stress from Hooke's law with P=0:
 *   sigma_xx = (lambda + 2*mu)*a
 *   sigma_yy = sigma_zz = lambda*a
 *   sigma_xy = sigma_xz = sigma_yz = 0
 *
 * This test creates a DMPlex box, sets the displacement FE, creates a
 * FaultCohesiveDyn with a single vertex whose negative-side vertex is a
 * mesh vertex, projects the linear field, and calls sampleStressAtFaults.
 * It then checks the recovered stress matches the analytical answer.
 */
TEST_F(CoulombStressTransferTest, GradientRecoveryLinearField) {
    // Create a small 3D simplex mesh (2x2x2 box from [0,1]^3)
    DM dm;
    PetscErrorCode ierr;
    const PetscInt dim = 3;
    const PetscInt faces[3] = {2, 2, 2};
    const PetscReal lower[3] = {0.0, 0.0, 0.0};
    const PetscReal upper[3] = {1.0, 1.0, 1.0};

    ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, dim, PETSC_TRUE,
                                faces, lower, upper, nullptr,
                                PETSC_TRUE, 1, PETSC_TRUE, &dm);
    ASSERT_EQ(ierr, 0);

    // Add a vector displacement FE (P1 Lagrange, dim components)
    PetscFE fe;
    ierr = PetscFECreateLagrange(PETSC_COMM_SELF, dim, dim, PETSC_TRUE,
                                  1, PETSC_DETERMINE, &fe);
    ASSERT_EQ(ierr, 0);
    ierr = DMSetField(dm, 0, nullptr, (PetscObject)fe);
    ASSERT_EQ(ierr, 0);
    ierr = DMCreateDS(dm);
    ASSERT_EQ(ierr, 0);
    ierr = PetscFEDestroy(&fe);
    ASSERT_EQ(ierr, 0);

    // Create a global vector and project the linear displacement field
    // u_x = a*x, u_y = 0, u_z = 0
    Vec sol;
    ierr = DMCreateGlobalVector(dm, &sol);
    ASSERT_EQ(ierr, 0);

    const double a = 1.0e-4;  // strain magnitude

    // Project function: u(x) = (a*x, 0, 0)
    auto linearDisp = [](PetscInt dim_in, PetscReal t, const PetscReal x[],
                          PetscInt Nc, PetscScalar *u, void *ctx) -> PetscErrorCode {
        (void)dim_in; (void)t; (void)Nc;
        double *aa = (double *)ctx;
        u[0] = (*aa) * x[0];
        u[1] = 0.0;
        u[2] = 0.0;
        return 0;
    };

    // PetscDS project function
    double a_ctx = a;
    PetscErrorCode (*funcs[1])(PetscInt, PetscReal, const PetscReal[], PetscInt, PetscScalar *, void *) = {linearDisp};
    void *ctxs[1] = {&a_ctx};

    ierr = DMProjectFunction(dm, 0.0, funcs, ctxs, INSERT_ALL_VALUES, sol);
    ASSERT_EQ(ierr, 0);

    // Find a vertex in the interior of the mesh to use as a fault vertex
    PetscInt vStart, vEnd;
    ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);
    ASSERT_EQ(ierr, 0);

    // Pick the first vertex
    PetscInt testVert = vStart;

    // Create a FaultCohesiveDyn with a single vertex
    FaultCohesiveDyn fault;
    FaultVertex fv;
    fv.vertex_negative = testVert;
    fv.vertex_positive = -1;  // not used for gradient recovery
    fv.normal = {0.0, 0.0, 1.0};
    fv.along_strike = {1.0, 0.0, 0.0};
    fv.up_dip = {0.0, 1.0, 0.0};
    fault.setFaultVertices({fv});

    // Initialize CoulombStressTransfer
    const double lambda = 40.0e9;
    const double mu = 30.0e9;
    const double biot = 0.0;  // no pore pressure coupling for this test

    CoulombStressTransfer cst(PETSC_COMM_WORLD);
    ierr = cst.initialize(dm, &fault, lambda, mu, biot);
    ASSERT_EQ(ierr, 0);

    // Sample stress
    ierr = cst.sampleStressAtFaults(sol);
    ASSERT_EQ(ierr, 0);

    // Check recovered stress
    const auto& vs = cst.getCurrentStress()[0];

    // Expected stress for eps_xx = a, all other eps = 0:
    double expected_sxx = (lambda + 2.0 * mu) * a;
    double expected_syy = lambda * a;
    double expected_szz = lambda * a;
    double expected_sxy = 0.0;
    double expected_sxz = 0.0;
    double expected_syz = 0.0;

    // Relative error tolerance: 1e-10
    double scale = std::abs(expected_sxx);
    EXPECT_NEAR(vs.sigma[0], expected_sxx, scale * 1e-10)
        << "sigma_xx should be (lambda + 2*mu)*a";
    EXPECT_NEAR(vs.sigma[1], expected_syy, scale * 1e-10)
        << "sigma_yy should be lambda*a";
    EXPECT_NEAR(vs.sigma[2], expected_szz, scale * 1e-10)
        << "sigma_zz should be lambda*a";
    EXPECT_NEAR(vs.sigma[3], expected_sxy, scale * 1e-10)
        << "sigma_xy should be 0";
    EXPECT_NEAR(vs.sigma[4], expected_sxz, scale * 1e-10)
        << "sigma_xz should be 0";
    EXPECT_NEAR(vs.sigma[5], expected_syz, scale * 1e-10)
        << "sigma_yz should be 0";

    // Clean up
    ierr = VecDestroy(&sol);
    ASSERT_EQ(ierr, 0);
    ierr = DMDestroy(&dm);
    ASSERT_EQ(ierr, 0);
}
