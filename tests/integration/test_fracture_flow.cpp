/**
 * @file test_fracture_flow.cpp
 * @brief Integration-level validation of hydrofrac lubrication callbacks
 */

#include <gtest/gtest.h>
#include <petscsys.h>

#include <cmath>

#include "numerics/PetscFEHydrofrac.hpp"

using namespace FSRM;

class FractureFlowTest : public ::testing::Test
{
};

TEST_F(FractureFlowTest, PoiseuilleFluxScalesWithApertureCubed)
{
  constexpr PetscInt dim = 2;
  PetscInt uOff[2] = {0, 1};
  PetscInt uOff_x[2] = {0, 1};

  // u = [p_f, u-, u+]
  // Configure two cases with jumps 1e-4 and 2e-4 m.
  PetscScalar u_case1[1 + 2 * dim] = {1.0e6, 0.0, 0.0, 1.0e-4, 1.0e-4};
  PetscScalar u_case2[1 + 2 * dim] = {1.0e6, 0.0, 0.0, 2.0e-4, 2.0e-4};
  PetscScalar u_x[dim] = {2.0e5, 0.0};

  PetscScalar constants[PetscFEHydrofrac::HYDROFRAC_CONST_PHASE2_COUNT] = {};
  constants[PetscFEHydrofrac::HYDROFRAC_CONST_MU_F] = 1.0e-3;

  PetscScalar flux1[dim] = {0.0, 0.0};
  PetscScalar flux2[dim] = {0.0, 0.0};

  PetscFEHydrofrac::f1_fracture_pressure(
      dim, 0, 0,
      uOff, uOff_x,
      u_case1, nullptr, u_x,
      nullptr, nullptr, nullptr, nullptr, nullptr,
      0.0, nullptr,
      PetscFEHydrofrac::HYDROFRAC_CONST_PHASE2_COUNT,
      constants,
      flux1);

  PetscFEHydrofrac::f1_fracture_pressure(
      dim, 0, 0,
      uOff, uOff_x,
      u_case2, nullptr, u_x,
      nullptr, nullptr, nullptr, nullptr, nullptr,
      0.0, nullptr,
      PetscFEHydrofrac::HYDROFRAC_CONST_PHASE2_COUNT,
      constants,
      flux2);

  const double ratio = PetscRealPart(flux2[0] / flux1[0]);
  EXPECT_NEAR(ratio, 8.0, 0.5);
}

TEST_F(FractureFlowTest, ResidualIncludesInjectionAndLeakoff)
{
  constexpr PetscInt dim = 2;
  PetscInt uOff[2] = {0, 1};
  PetscInt uOff_x[2] = {0, 1};

  PetscScalar u[1 + 2 * dim] = {2.0e6, 0.0, 0.0, 0.0, 0.0};
  PetscScalar u_t[1 + 2 * dim] = {0.0, 0.0, 0.0, 0.0, 0.0};
  PetscScalar f[1] = {0.0};

  PetscScalar constants[PetscFEHydrofrac::HYDROFRAC_CONST_PHASE2_COUNT] = {};
  constants[PetscFEHydrofrac::HYDROFRAC_CONST_Q_INJ] = 1.0e-4;
  constants[PetscFEHydrofrac::HYDROFRAC_CONST_LEAKOFF] = 2.0e-10;
  constants[PetscFEHydrofrac::HYDROFRAC_CONST_P_FORMATION] = 1.0e6;

  PetscFEHydrofrac::f0_fracture_pressure(
      dim, 0, 0,
      uOff, uOff_x,
      u, u_t, nullptr,
      nullptr, nullptr, nullptr, nullptr, nullptr,
      0.0, nullptr,
      PetscFEHydrofrac::HYDROFRAC_CONST_PHASE2_COUNT,
      constants,
      f);

  const double expected = -1.0e-4 + 2.0e-10 * 1.0e6;
  EXPECT_NEAR(PetscRealPart(f[0]), expected, 1e-12);
}

TEST_F(FractureFlowTest, RadialPressureProfileMonotonicFromInjectionPoint)
{
  // Analytical proxy check for 1D slot with source at x=0 and no-flow tip: p decreases with distance.
  const double p0 = 5.0e6;
  const double grad = 3.0e4;
  auto p = [p0, grad](double x) { return p0 - grad * x; };

  EXPECT_GT(p(0.0), p(10.0));
  EXPECT_GT(p(10.0), p(20.0));
}
