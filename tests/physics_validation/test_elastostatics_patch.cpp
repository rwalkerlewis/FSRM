/**
 * @file test_elastostatics_patch.cpp
 * @brief Patch test: uniform strain should produce exact stress
 *
 * Setup: Unit cube, Q1 hex elements
 * Prescribed displacement: u_x = eps_0 * x, u_y = u_z = 0
 * Check: sigma_xx = (lambda + 2*mu)*eps_0, sigma_yy = sigma_zz = lambda*eps_0
 * This validates the stiffness tensor g3 is correct
 */

#include "physics/PhysicsKernel.hpp"
#include "numerics/PetscFEElasticity.hpp"
#include "core/FSRM.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace FSRM;

class ElastostaticsPatchTest : public ::testing::Test {
protected:
    void SetUp() override {
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    }
    int rank = 0;
};

TEST_F(ElastostaticsPatchTest, F1ProducesCorrectStress) {
    // Test that f1_elastostatics computes stress correctly for uniform strain
    const double E = 10e9;              // 10 GPa Young's modulus
    const double nu = 0.25;              // Poisson's ratio
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));
    const double eps_0 = 0.001;          // 0.1% strain

    // Setup constants array with Simulator's layout
    PetscScalar constants[25] = {0};
    constants[0] = lambda;
    constants[1] = mu;
    constants[2] = 2700;                 // rho_solid

    // Displacement gradient: uniform uniaxial strain eps_xx = eps_0, eps_yy = eps_zz = 0
    // PETSc layout: u_x[uOff_x[f] + c*dim + d] = d(u_c)/dx_d
    PetscScalar u_x[9] = {0};
    u_x[0] = eps_0;                      // du_x/dx = eps_0
    u_x[1] = 0.0;                        // du_x/dy = 0
    u_x[2] = 0.0;                        // du_x/dz = 0
    u_x[3] = 0.0;                        // du_y/dx = 0
    u_x[4] = 0.0;                        // du_y/dy = 0
    u_x[5] = 0.0;                        // du_y/dz = 0
    u_x[6] = 0.0;                        // du_z/dx = 0
    u_x[7] = 0.0;                        // du_z/dy = 0
    u_x[8] = 0.0;                        // du_z/dz = 0

    PetscInt uOff[1] = {0};
    PetscInt uOff_x[1] = {0};
    PetscScalar f1[9] = {0};
    PetscReal x[3] = {0.5, 0.5, 0.5};

    // Call the PetscFE callback directly
    PetscFEElasticity::f1_elastostatics(
        3,          // dim
        1,          // Nf (number of fields)
        0,          // NfAux
        uOff, uOff_x,
        nullptr,    // u (not needed for f1_elastostatics)
        nullptr,    // u_t
        u_x,
        nullptr,    // aOff
        nullptr,    // aOff_x
        nullptr,    // a
        nullptr,    // a_x
        nullptr,    // a_t
        0.0,        // t
        x,
        25,         // numConstants
        constants,
        f1);

    // f1[c*dim+d] = sigma_{cd}
    // Expected: sigma_xx = (lambda+2*mu)*eps_0, sigma_yy = sigma_zz = lambda*eps_0
    const double sigma_xx_expected = (lambda + 2.0*mu) * eps_0;
    const double sigma_yy_expected = lambda * eps_0;
    const double sigma_zz_expected = lambda * eps_0;

    // Check diagonal stress components
    EXPECT_NEAR(PetscRealPart(f1[0]), sigma_xx_expected, 1e-6*E);  // sigma_xx = f1[0*3+0]
    EXPECT_NEAR(PetscRealPart(f1[4]), sigma_yy_expected, 1e-6*E);  // sigma_yy = f1[1*3+1]
    EXPECT_NEAR(PetscRealPart(f1[8]), sigma_zz_expected, 1e-6*E);  // sigma_zz = f1[2*3+2]

    // Check off-diagonal stress components (should be zero for uniform uniaxial strain)
    EXPECT_NEAR(PetscRealPart(f1[1]), 0.0, 1e-12);  // sigma_xy = f1[0*3+1]
    EXPECT_NEAR(PetscRealPart(f1[2]), 0.0, 1e-12);  // sigma_xz = f1[0*3+2]
    EXPECT_NEAR(PetscRealPart(f1[3]), 0.0, 1e-12);  // sigma_yx = f1[1*3+0]
    EXPECT_NEAR(PetscRealPart(f1[5]), 0.0, 1e-12);  // sigma_yz = f1[1*3+2]
    EXPECT_NEAR(PetscRealPart(f1[6]), 0.0, 1e-12);  // sigma_zx = f1[2*3+0]
    EXPECT_NEAR(PetscRealPart(f1[7]), 0.0, 1e-12);  // sigma_zy = f1[2*3+1]

    (void)rank;
}

TEST_F(ElastostaticsPatchTest, G3ProducesCorrectStiffness) {
    // Test that g3_elastostatics produces correct stiffness tensor
    const double E = 10e9;               // 10 GPa Young's modulus
    const double nu = 0.25;               // Poisson's ratio
    const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = E / (2.0 * (1.0 + nu));

    PetscScalar constants[25] = {0};
    constants[0] = lambda;
    constants[1] = mu;
    constants[2] = 2700;                 // rho_solid

    PetscInt uOff[1] = {0};
    PetscInt uOff_x[1] = {0};
    PetscScalar g3[81] = {0};            // 9x9 matrix for 3D elasticity (3 components x 3 directions x 3 components x 3 directions)
    PetscReal x[3] = {0.5, 0.5, 0.5};

    // Call the PetscFE Jacobian callback
    PetscFEElasticity::g3_elastostatics(
        3,          // dim
        1,          // Nf
        0,          // NfAux
        uOff, uOff_x,
        nullptr,    // u
        nullptr,    // u_t
        nullptr,    // u_x
        nullptr,    // aOff
        nullptr,    // aOff_x
        nullptr,    // a
        nullptr,    // a_x
        nullptr,    // a_t
        0.0,        // t
        0.0,        // u_tShift (not used in elastostatics)
        x,
        25,         // numConstants
        constants,
        g3);

    // g3[i*dim*dim*dim + j*dim*dim + d1*dim + d2] = C_{i,d1,j,d2}
    // = d(sigma_{i,d1})/d(du_j/dx_d2)
    //
    // For isotropic elasticity:
    // C_{i,d1,j,d2} = lambda*delta_{i,d1}*delta_{j,d2} + mu*(delta_{i,j}*delta_{d1,d2} + delta_{i,d2}*delta_{d1,j})

    // Check diagonal entries: d(sigma_xx)/d(du_x/dx) = lambda + 2*mu
    // i=0, d1=0, j=0, d2=0: g3[0*27 + 0*9 + 0*3 + 0] = g3[0]
    EXPECT_NEAR(PetscRealPart(g3[0]), lambda + 2*mu, 1e-6*E);

    // Check coupling: d(sigma_yy)/d(du_x/dx) = lambda
    // i=1, d1=1, j=0, d2=0: g3[1*27 + 0*9 + 1*3 + 0] = g3[27 + 3] = g3[30]
    EXPECT_NEAR(PetscRealPart(g3[30]), lambda, 1e-6*E);

    // Check shear: d(sigma_xy)/d(du_x/dy) = mu
    // i=0, d1=1, j=0, d2=1: g3[0*27 + 0*9 + 1*3 + 1] = g3[4]
    EXPECT_NEAR(PetscRealPart(g3[4]), mu, 1e-6*E);

    // Check shear: d(sigma_xy)/d(du_y/dx) = mu
    // i=0, d1=1, j=1, d2=0: g3[0*27 + 1*9 + 1*3 + 0] = g3[12]
    EXPECT_NEAR(PetscRealPart(g3[12]), mu, 1e-6*E);

    (void)rank;
}
