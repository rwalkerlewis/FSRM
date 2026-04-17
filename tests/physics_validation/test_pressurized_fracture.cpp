/**
 * @file test_pressurized_fracture.cpp
 * @brief Phase 1 validation for pressurized cohesive fracture callbacks
 */

#include <gtest/gtest.h>
#include <petscsys.h>

#include <cmath>

#include "numerics/PetscFEHydrofrac.hpp"

namespace
{

class PressurizedFractureValidationTest : public ::testing::Test
{
};

TEST_F(PressurizedFractureValidationTest, SneddonApertureIsPositiveAndFinite)
{
  const PetscReal pressure_pa = 8.0e6;
  const PetscReal radius_m = 12.0;
  const PetscReal youngs_modulus_pa = 25.0e9;
  const PetscReal nu = 0.24;

  const PetscReal w_max = FSRM::PetscFEHydrofrac::sneddonPennyCrackMaxAperture(
      pressure_pa, radius_m, youngs_modulus_pa, nu);

  EXPECT_GT(w_max, 0.0);
  EXPECT_LT(w_max, 0.10);
}

TEST_F(PressurizedFractureValidationTest, SneddonProfileDecreasesTowardTip)
{
  const PetscReal pressure_pa = 6.0e6;
  const PetscReal radius_m = 10.0;
  const PetscReal youngs_modulus_pa = 20.0e9;
  const PetscReal nu = 0.25;

  const PetscReal w_max = FSRM::PetscFEHydrofrac::sneddonPennyCrackMaxAperture(
      pressure_pa, radius_m, youngs_modulus_pa, nu);

  auto profile = [w_max, radius_m](PetscReal r) {
    const PetscReal xi = r / radius_m;
    return w_max * std::sqrt(std::max(0.0, 1.0 - xi * xi));
  };

  const PetscReal w_center = profile(0.0);
  const PetscReal w_mid = profile(0.5 * radius_m);
  const PetscReal w_tip = profile(0.95 * radius_m);

  EXPECT_GT(w_center, w_mid);
  EXPECT_GT(w_mid, w_tip);
  EXPECT_GT(w_tip, 0.0);
}

TEST_F(PressurizedFractureValidationTest, PressureBalanceResidualIsZeroAtEquilibrium)
{
  constexpr PetscInt dim = 3;
  constexpr PetscInt Nf = 2;

  // Field layout on cohesive evaluation point: [u] [lambda]
  PetscInt uOff[3] = {0, 2 * dim, 3 * dim};
  PetscScalar u[9] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      -5.0e6, 0.0, 0.0};

  PetscScalar constants[FSRM::PetscFEHydrofrac::HYDROFRAC_CONST_COUNT] = {};
  constants[FSRM::PetscFEHydrofrac::HYDROFRAC_CONST_PRESSURE] = 5.0e6;

  const PetscReal x[3] = {0.0, 0.0, 0.0};
  const PetscReal n[3] = {1.0, 0.0, 0.0};
  PetscScalar f[3] = {1.0, 1.0, 1.0};

  FSRM::PetscFEHydrofrac::f0_lagrange_pressure_balance(
      dim, Nf, 0,
      uOff, nullptr,
      u, nullptr, nullptr,
      nullptr, nullptr, nullptr, nullptr, nullptr,
      0.0, x, n,
      FSRM::PetscFEHydrofrac::HYDROFRAC_CONST_COUNT,
      constants,
      f);

  EXPECT_NEAR(PetscRealPart(f[0]), 0.0, 1e-10);
  EXPECT_NEAR(PetscRealPart(f[1]), 0.0, 1e-10);
  EXPECT_NEAR(PetscRealPart(f[2]), 0.0, 1e-10);
}

TEST_F(PressurizedFractureValidationTest, PressureBalanceResidualSignTracksCompression)
{
  constexpr PetscInt dim = 3;
  PetscInt uOff[3] = {0, 2 * dim, 3 * dim};
  PetscScalar u[9] = {
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      -3.0e6, 0.0, 0.0};

  PetscScalar constants[FSRM::PetscFEHydrofrac::HYDROFRAC_CONST_COUNT] = {};
  constants[FSRM::PetscFEHydrofrac::HYDROFRAC_CONST_PRESSURE] = 5.0e6;

  const PetscReal x[3] = {0.0, 0.0, 0.0};
  const PetscReal n[3] = {1.0, 0.0, 0.0};
  PetscScalar f[3] = {0.0, 0.0, 0.0};

  FSRM::PetscFEHydrofrac::f0_lagrange_pressure_balance(
      dim, 2, 0,
      uOff, nullptr,
      u, nullptr, nullptr,
      nullptr, nullptr, nullptr, nullptr, nullptr,
      0.0, x, n,
      FSRM::PetscFEHydrofrac::HYDROFRAC_CONST_COUNT,
      constants,
      f);

  // Residual > 0 means lambda is less compressive than required by pressure loading.
  EXPECT_GT(PetscRealPart(f[0]), 0.0);
}

} // namespace
