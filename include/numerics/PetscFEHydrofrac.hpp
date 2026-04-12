#ifndef PETSC_FE_HYDROFRAC_HPP
#define PETSC_FE_HYDROFRAC_HPP

#include <petscds.h>

namespace FSRM
{

/**
 * @brief PETSc pointwise callbacks for FEM-coupled hydrofracture terms.
 */
class PetscFEHydrofrac
{
public:
  // Constants slot for prescribed fracture fluid pressure (Pa).
  // Slots 0-30 are used by existing callbacks (fluid, elasticity, cohesive).
  static constexpr PetscInt HYDROFRAC_CONST_PRESSURE = 31;
  static constexpr PetscInt HYDROFRAC_CONST_COUNT = 32;
    static constexpr PetscInt HYDROFRAC_CONST_MU_F = 32;
    static constexpr PetscInt HYDROFRAC_CONST_Q_INJ = 33;
    static constexpr PetscInt HYDROFRAC_CONST_LEAKOFF = 34;
    static constexpr PetscInt HYDROFRAC_CONST_P_FORMATION = 35;
    static constexpr PetscInt HYDROFRAC_CONST_PHASE2_COUNT = 36;

  /**
   * @brief Cohesive traction balance for a pressurized fracture.
   *
   * Enforces lambda + p_f * n = 0 on cohesive interfaces.
   */
  static void f0_lagrange_pressure_balance(
      PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[],
      const PetscScalar u[], const PetscScalar u_t[],
      const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[],
      const PetscScalar a[], const PetscScalar a_t[],
      const PetscScalar a_x[],
      PetscReal t, const PetscReal x[],
      const PetscReal n[],
      PetscInt numConstants,
      const PetscScalar constants[],
      PetscScalar f[]);

  /**
   * @brief Fracture pressure residual source term for lubrication equation.
   *
   * R_pf includes aperture-rate coupling and source/sink terms:
   * dw/dt - q_inj + q_leak.
   */
  static void f0_fracture_pressure(
      PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[],
      const PetscScalar u[], const PetscScalar u_t[],
      const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[],
      const PetscScalar a[], const PetscScalar a_t[],
      const PetscScalar a_x[],
      PetscReal t, const PetscReal x[],
      PetscInt numConstants,
      const PetscScalar constants[],
      PetscScalar f[]);

  /**
   * @brief Fracture pressure diffusion flux for lubrication equation.
   *
   * f1 = k_eff * grad(p_f), with k_eff = w^3 / (12*mu_f).
   */
  static void f1_fracture_pressure(
      PetscInt dim, PetscInt Nf, PetscInt NfAux,
      const PetscInt uOff[], const PetscInt uOff_x[],
      const PetscScalar u[], const PetscScalar u_t[],
      const PetscScalar u_x[],
      const PetscInt aOff[], const PetscInt aOff_x[],
      const PetscScalar a[], const PetscScalar a_t[],
      const PetscScalar a_x[],
      PetscReal t, const PetscReal x[],
      PetscInt numConstants,
      const PetscScalar constants[],
      PetscScalar f1[]);

  /**
   * @brief Maximum aperture for Sneddon's penny-shaped crack solution.
   */
  static PetscReal sneddonPennyCrackMaxAperture(
      PetscReal pressure,
      PetscReal radius,
      PetscReal youngs_modulus,
      PetscReal poissons_ratio);

  /**
   * @brief Width scaling proxy for viscosity-dominated PKN regime.
   */
  static PetscReal pknWidthScaling(
      PetscReal fluid_viscosity,
      PetscReal injection_rate,
      PetscReal time,
      PetscReal plane_strain_modulus);

  /**
   * @brief Cohesive tensile strength from K_Ic and mesh length scale.
   */
  static PetscReal cohesiveStrengthFromKIc(PetscReal k_ic, PetscReal h);

  /**
   * @brief Opening criterion for fracture propagation.
   */
  static PetscBool shouldOpenFractureCell(
      PetscReal normal_traction,
      PetscReal tensile_strength);

  // --- Phase 5: Stress shadowing utilities ---

  /**
   * @brief Normal stress perturbation from a pressurized penny-shaped crack.
   *
   * On the plane of the crack at radial distance r from center,
   * the normal stress increment is:
   *   delta_sigma_n = P * (2/pi) * arcsin(a/r) for r > a
   *   delta_sigma_n = P                         for r <= a
   * (Sneddon 1946, plane of the crack).
   */
  static PetscReal sneddonNormalStressPerturbation(
      PetscReal pressure,
      PetscReal crack_radius,
      PetscReal distance_from_center);

  /**
   * @brief Stress shadow factor: ratio of additional closure stress
   *        to net pressure for a neighboring fracture at given spacing.
   */
  static PetscReal stressShadowFactor(
      PetscReal net_pressure,
      PetscReal fracture_half_length,
      PetscReal spacing);

  /**
   * @brief Cluster efficiency: fraction of fluid entering each cluster
   *        when N parallel fractures compete under stress shadowing.
   */
  static PetscReal clusterEfficiency(
      PetscInt num_clusters,
      PetscReal spacing,
      PetscReal fracture_half_length,
      PetscReal net_pressure,
      PetscInt cluster_index);

  // --- Phase 7: Induced seismicity utilities ---

  /**
   * @brief Scalar seismic moment from slip, area, and shear modulus.
   */
  static PetscReal scalarMoment(
      PetscReal shear_modulus,
      PetscReal average_slip,
      PetscReal rupture_area);

  /**
   * @brief Moment magnitude from scalar seismic moment (Hanks and Kanamori).
   */
  static PetscReal momentMagnitude(PetscReal scalar_moment);

  /**
   * @brief Check if a microseismic event magnitude is in typical range.
   */
  static PetscBool isMicroseismicRange(PetscReal mw);

  // --- Phase 6: Proppant transport utilities ---

  /**
   * @brief Stokes settling velocity for a proppant particle.
   */
  static PetscReal stokesSettlingVelocity(
      PetscReal particle_diameter,
      PetscReal particle_density,
      PetscReal fluid_density,
      PetscReal fluid_viscosity);

  /**
   * @brief Check if proppant bridging occurs (width/diameter criterion).
   */
  static PetscBool proppantBridging(
      PetscReal aperture,
      PetscReal particle_diameter,
      PetscReal bridging_ratio);

  /**
   * @brief Minimum aperture held open by a proppant pack.
   */
  static PetscReal proppantPackMinAperture(
      PetscReal particle_diameter,
      PetscReal pack_porosity);

  /**
   * @brief 1D advection mass conservation check for proppant.
   */
  static PetscReal proppantMassBalance(
      PetscReal concentration_in,
      PetscReal flow_rate,
      PetscReal injection_time,
      PetscReal fracture_volume,
      PetscReal average_concentration);

  // --- Phase 8: Carter leak-off utilities ---

  /**
   * @brief Carter leak-off rate at time t since exposure.
   *
   * q_l = C_L / sqrt(t - t_open) for t > t_open, else 0.
   */
  static PetscReal carterLeakoffRate(
      PetscReal c_l,
      PetscReal time,
      PetscReal time_open);

  /**
   * @brief Cumulative Carter leak-off volume per unit area.
   *
   * V_l = 2 * C_L * sqrt(t - t_open)
   */
  static PetscReal carterCumulativeLeakoff(
      PetscReal c_l,
      PetscReal time,
      PetscReal time_open);

  /**
   * @brief Total leak-off volume from a fracture face.
   */
  static PetscReal totalLeakoffVolume(
      PetscReal c_l,
      PetscReal time,
      PetscReal time_open,
      PetscReal face_area);

  // --- Phase 9: Production forecasting utilities ---

  /**
   * @brief Arps hyperbolic decline rate.
   *
   * q(t) = q_i / (1 + b * D_i * t)^(1/b)
   */
  static PetscReal arpsDeclineRate(
      PetscReal q_initial,
      PetscReal decline_rate,
      PetscReal b_factor,
      PetscReal time);

  /**
   * @brief Cumulative production from Arps hyperbolic decline.
   *
   * N_p(t) = q_i / ((1-b)*D_i) * [(1+b*D_i*t)^(1-1/b) - 1]  for b != 1
   */
  static PetscReal arpsCumulativeProduction(
      PetscReal q_initial,
      PetscReal decline_rate,
      PetscReal b_factor,
      PetscReal time);

  /**
   * @brief Productivity index from fracture conductivity.
   *
   * J = k_f * w_f / (mu * ln(r_e/r_w))
   */
  static PetscReal fractureProductivityIndex(
      PetscReal fracture_permeability,
      PetscReal fracture_width,
      PetscReal fluid_viscosity,
      PetscReal drainage_radius,
      PetscReal wellbore_radius);
};

} // namespace FSRM

#endif // PETSC_FE_HYDROFRAC_HPP
