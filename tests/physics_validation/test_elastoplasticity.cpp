/**
 * @file test_elastoplasticity.cpp
 * @brief Physics validation for PetscDS elastoplastic callbacks
 */

#include <gtest/gtest.h>
#include <array>
#include <cmath>

#include "numerics/PetscFEElastoplasticity.hpp"
#include "domain/geomechanics/PlasticityModel.hpp"

using namespace FSRM;

namespace
{

std::array<std::array<double, 6>, 6> makeElasticModulus(double E, double nu)
{
  std::array<std::array<double, 6>, 6> C{};
  const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const double mu = E / (2.0 * (1.0 + nu));

  C[0][0] = lambda + 2.0 * mu;
  C[0][1] = lambda;
  C[0][2] = lambda;
  C[1][0] = lambda;
  C[1][1] = lambda + 2.0 * mu;
  C[1][2] = lambda;
  C[2][0] = lambda;
  C[2][1] = lambda;
  C[2][2] = lambda + 2.0 * mu;

  C[3][3] = mu;
  C[4][4] = mu;
  C[5][5] = mu;

  return C;
}

double yieldFunctionDP(const std::array<double, 6>& stress,
    double cohesion, double friction_angle)
{
  const double I1 = stress[0] + stress[1] + stress[2];
  const double p = I1 / 3.0;
  const double s11 = stress[0] - p;
  const double s22 = stress[1] - p;
  const double s33 = stress[2] - p;
  const double s12 = stress[3];
  const double s13 = stress[4];
  const double s23 = stress[5];
  const double J2 = 0.5 * (s11*s11 + s22*s22 + s33*s33) + s12*s12 + s13*s13 + s23*s23;

  const double sin_phi = std::sin(friction_angle);
  const double cos_phi = std::cos(friction_angle);
  const double alpha = 6.0 * sin_phi / (std::sqrt(3.0) * (3.0 - sin_phi));
  const double kappa = 6.0 * cohesion * cos_phi / (std::sqrt(3.0) * (3.0 - sin_phi));

  return std::sqrt(std::max(0.0, J2)) + alpha * I1 - kappa;
}

} // namespace

TEST(ElastoplasticityTest, ElastoplasticCallbackReturnsNearYieldStress)
{
  const PetscInt dim = 3;
  const PetscInt Nf = 1;
  const PetscInt NfAux = PetscFEElastoplasticity::EP_NUM_AUX_FIELDS;

  PetscInt uOff[1] = {0};
  PetscInt uOff_x[1] = {0};

  PetscInt aOff[PetscFEElastoplasticity::EP_NUM_AUX_FIELDS];
  aOff[PetscFEElastoplasticity::EP_AUX_LAMBDA] = 0;
  aOff[PetscFEElastoplasticity::EP_AUX_MU] = 1;
  aOff[PetscFEElastoplasticity::EP_AUX_RHO] = 2;
  aOff[PetscFEElastoplasticity::EP_AUX_PLASTIC_STRAIN] = 3;
  aOff[PetscFEElastoplasticity::EP_AUX_ACC_PLASTIC_STRAIN] = 9;

  PetscInt aOff_x[PetscFEElastoplasticity::EP_NUM_AUX_FIELDS] = {0, 0, 0, 0, 0};

  // a = [lambda, mu, rho, epsp(6), acc_epsp]
  PetscScalar a[10] = {0.0};
  const double E = 10.0e9;
  const double nu = 0.25;
  const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const double mu = E / (2.0 * (1.0 + nu));
  a[0] = lambda;
  a[1] = mu;
  a[2] = 2650.0;

  PetscScalar u[3] = {0.0, 0.0, 0.0};
  PetscScalar u_t[3] = {0.0, 0.0, 0.0};

  // Gradient with pure shear-type deviatoric loading:
  // eps_xx = 0.01, eps_yy = -0.01, eps_zz = 0 => I1 = 0
  PetscScalar u_x[9] = {0.01, 0.0, 0.0,
                        0.0, -0.01, 0.0,
                        0.0, 0.0, 0.0};

  PetscScalar a_t[10] = {0.0};
  PetscScalar a_x[1] = {0.0};

  const double cohesion = 1.0e6;
  const double friction_angle = 30.0 * M_PI / 180.0;
  PetscScalar constants[PetscFEElastoplasticity::EP_CONST_COUNT] = {0.0};
  constants[PetscFEElastoplasticity::EP_CONST_COHESION] = cohesion;
  constants[PetscFEElastoplasticity::EP_CONST_FRICTION_ANGLE] = friction_angle;
  constants[PetscFEElastoplasticity::EP_CONST_DILATION_ANGLE] = 0.0;
  constants[PetscFEElastoplasticity::EP_CONST_HARDENING_MODULUS] = 0.0;

  PetscReal x[3] = {0.5, 0.5, 0.5};
  PetscScalar f1[9] = {0.0};

  PetscFEElastoplasticity::f1_elastoplastic_aux(dim, Nf, NfAux,
      uOff, uOff_x, u, u_t, u_x,
      aOff, aOff_x, a, a_t, a_x,
      0.0, x, PetscFEElastoplasticity::EP_CONST_COUNT,
      constants, f1);

  std::array<double, 6> stress = {
      PetscRealPart(f1[0]), PetscRealPart(f1[4]), PetscRealPart(f1[8]),
      PetscRealPart(f1[5]), PetscRealPart(f1[2]), PetscRealPart(f1[1])};

  const double F = yieldFunctionDP(stress, cohesion, friction_angle);
  EXPECT_NEAR(F, 0.0, 1.0e-4)
      << "Elastoplastic callback should return stress to DP yield surface";
}

TEST(ElastoplasticityTest, MaterialPointSolveStoresPlasticStateInAuxField)
{
  PlasticityModel model(YieldCriterion::DRUCKER_PRAGER);

  // Use phi=0 so uniaxial compression clearly yields (pressure-insensitive limit)
  DruckerPragerModel::Parameters params;
  params.friction_angle = 0.0;
  params.cohesion = 1.0e6;
  params.dilation_angle = 0.0;
  params.hardening_modulus = 0.0;
  model.setDruckerPragerParameters(params);

  const double E = 10.0e9;
  const double nu = 0.25;
  auto C = makeElasticModulus(E, nu);

  std::array<double, 6> stress = {0, 0, 0, 0, 0, 0};
  std::array<double, 6> strain_inc = {0, 0, -0.01, 0, 0, 0};
  PlasticityState state;

  // Material-point solve (return mapping)
  model.integrateStress(stress, strain_inc, C, state);

  ASSERT_TRUE(state.is_plastic);
  ASSERT_GT(state.accumulated_plastic_strain, 0.0);

  // Emulate storing updated plastic state into auxiliary field
  PetscScalar aux_state[7] = {0.0};  // [epsp(6), acc_eps_p]
  for (int i = 0; i < 6; ++i)
  {
    aux_state[i] = state.plastic_strain[i];
  }
  aux_state[6] = state.accumulated_plastic_strain;

  double epsp_norm = 0.0;
  for (int i = 0; i < 6; ++i)
  {
    epsp_norm += PetscRealPart(aux_state[i]) * PetscRealPart(aux_state[i]);
  }

  EXPECT_GT(epsp_norm, 0.0)
      << "Plastic strain tensor stored in aux field should be nonzero after solve";
  EXPECT_GT(PetscRealPart(aux_state[6]), 0.0)
      << "Accumulated plastic strain in aux field should be nonzero after solve";
}
