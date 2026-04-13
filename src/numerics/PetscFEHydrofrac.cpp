#include "numerics/PetscFEHydrofrac.hpp"

#include <petscmath.h>
#include <vector>

namespace FSRM
{

void PetscFEHydrofrac::f0_lagrange_pressure_balance(
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
    PetscScalar f[])
{
  (void)Nf;
  (void)NfAux;
  (void)uOff_x;
  (void)u_t;
  (void)u_x;
  (void)aOff;
  (void)aOff_x;
  (void)a;
  (void)a_t;
  (void)a_x;
  (void)t;
  (void)x;

  const PetscScalar pressure =
      (numConstants > HYDROFRAC_CONST_PRESSURE) ? constants[HYDROFRAC_CONST_PRESSURE] : 0.0;

  // lambda + p_f * n = 0
  for (PetscInt d = 0; d < dim; ++d)
  {
    f[d] = u[uOff[1] + d] + pressure * n[d];
  }
}

void PetscFEHydrofrac::f0_fracture_pressure(
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
    PetscScalar f[])
{
  (void)Nf;
  (void)NfAux;
  (void)uOff_x;
  (void)u;
  (void)u_x;
  (void)aOff;
  (void)aOff_x;
  (void)a;
  (void)a_t;
  (void)a_x;
  (void)t;
  (void)x;

  // Expected layout for coupled cohesive callbacks:
  // field 0: fracture pressure p_f (scalar)
  // field 1: displacement jump pair [u-, u+] packed by PETSc hybrid integration
  PetscScalar dw_dt = 0.0;
  if (u_t && dim > 0 && uOff[1] >= 0) {
    for (PetscInt d = 0; d < dim; ++d) {
      dw_dt += u_t[uOff[1] + dim + d] - u_t[uOff[1] + d];
    }
    dw_dt /= static_cast<PetscScalar>(dim);
  }

  const PetscScalar q_inj =
      (numConstants > HYDROFRAC_CONST_Q_INJ) ? constants[HYDROFRAC_CONST_Q_INJ] : 0.0;
  const PetscScalar c_l =
      (numConstants > HYDROFRAC_CONST_LEAKOFF) ? constants[HYDROFRAC_CONST_LEAKOFF] : 0.0;
  const PetscScalar p_form =
      (numConstants > HYDROFRAC_CONST_P_FORMATION) ? constants[HYDROFRAC_CONST_P_FORMATION] : 0.0;
  const PetscScalar p_frac = u[uOff[0]];

  const PetscScalar q_leak = c_l * (p_frac - p_form);
  f[0] = dw_dt - q_inj + q_leak;
}

void PetscFEHydrofrac::f1_fracture_pressure(
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
    PetscScalar f1[])
{
  (void)Nf;
  (void)NfAux;
  (void)uOff_x;
  (void)u_t;
  (void)aOff;
  (void)aOff_x;
  (void)a;
  (void)a_t;
  (void)a_x;
  (void)t;
  (void)x;

  const PetscScalar mu_f =
      (numConstants > HYDROFRAC_CONST_MU_F) ? constants[HYDROFRAC_CONST_MU_F] : 1.0e-3;

  PetscScalar w = 1.0e-6;
  if (dim > 0 && uOff[1] >= 0) {
    PetscScalar jump_sum = 0.0;
    for (PetscInt d = 0; d < dim; ++d) {
      jump_sum += u[uOff[1] + dim + d] - u[uOff[1] + d];
    }
    w = PetscMax(1.0e-9, PetscAbsScalar(jump_sum / static_cast<PetscScalar>(dim)));
  }

  const PetscScalar k_eff = (w * w * w) / (12.0 * mu_f);
  for (PetscInt d = 0; d < dim; ++d) {
    f1[d] = k_eff * u_x[uOff_x[0] + d];
  }
}

PetscReal PetscFEHydrofrac::sneddonPennyCrackMaxAperture(
    PetscReal pressure,
    PetscReal radius,
    PetscReal youngs_modulus,
    PetscReal poissons_ratio)
{
  if (youngs_modulus <= 0.0 || radius <= 0.0)
  {
    return 0.0;
  }

  return 8.0 * (1.0 - poissons_ratio * poissons_ratio) * pressure * radius /
         (PETSC_PI * youngs_modulus);
}

PetscReal PetscFEHydrofrac::pknWidthScaling(
    PetscReal fluid_viscosity,
    PetscReal injection_rate,
    PetscReal time,
    PetscReal plane_strain_modulus)
{
  if (fluid_viscosity <= 0.0 || injection_rate <= 0.0 || time <= 0.0 || plane_strain_modulus <= 0.0)
  {
    return 0.0;
  }

  const PetscReal arg = fluid_viscosity * injection_rate * time / plane_strain_modulus;
  return PetscPowReal(arg, 0.25);
}

PetscReal PetscFEHydrofrac::cohesiveStrengthFromKIc(PetscReal k_ic, PetscReal h)
{
  if (k_ic <= 0.0 || h <= 0.0)
  {
    return 0.0;
  }
  return k_ic / PetscSqrtReal(PETSC_PI * h / 2.0);
}

PetscBool PetscFEHydrofrac::shouldOpenFractureCell(
    PetscReal normal_traction,
    PetscReal tensile_strength)
{
  return normal_traction > tensile_strength ? PETSC_TRUE : PETSC_FALSE;
}

// --- Phase 5: Stress shadowing ---

PetscReal PetscFEHydrofrac::sneddonNormalStressPerturbation(
    PetscReal pressure,
    PetscReal crack_radius,
    PetscReal distance_from_center)
{
  if (pressure <= 0.0 || crack_radius <= 0.0 || distance_from_center < 0.0)
  {
    return 0.0;
  }
  if (distance_from_center <= crack_radius)
  {
    return pressure;
  }
  // On the crack plane, r > a: delta_sigma = P * (2/pi) * arcsin(a/r)
  return pressure * (2.0 / PETSC_PI) * PetscAsinReal(crack_radius / distance_from_center);
}

PetscReal PetscFEHydrofrac::stressShadowFactor(
    PetscReal net_pressure,
    PetscReal fracture_half_length,
    PetscReal spacing)
{
  if (net_pressure <= 0.0 || fracture_half_length <= 0.0 || spacing <= 0.0)
  {
    return 0.0;
  }
  // Stress shadow factor is the normal stress perturbation at the neighbor
  // divided by the net pressure
  PetscReal delta_sigma = sneddonNormalStressPerturbation(
      net_pressure, fracture_half_length, spacing);
  return delta_sigma / net_pressure;
}

PetscReal PetscFEHydrofrac::clusterEfficiency(
    PetscInt num_clusters,
    PetscReal spacing,
    PetscReal fracture_half_length,
    PetscReal net_pressure,
    PetscInt cluster_index)
{
  if (num_clusters <= 0 || spacing <= 0.0 || fracture_half_length <= 0.0 ||
      net_pressure <= 0.0 || cluster_index < 0 || cluster_index >= num_clusters)
  {
    return 0.0;
  }
  if (num_clusters == 1)
  {
    return 1.0;
  }

  // Compute cumulative stress shadow on each cluster from all others
  std::vector<PetscReal> shadow(num_clusters, 0.0);
  for (PetscInt i = 0; i < num_clusters; ++i)
  {
    for (PetscInt j = 0; j < num_clusters; ++j)
    {
      if (i == j) continue;
      PetscReal dist = PetscAbsReal(static_cast<PetscReal>(i - j)) * spacing;
      shadow[i] += sneddonNormalStressPerturbation(
          net_pressure, fracture_half_length, dist);
    }
  }

  // Efficiency inversely proportional to shadow: e_i = (1/s_i) / sum(1/s_j)
  // Higher shadow means harder to open, lower efficiency
  PetscReal sum_inv = 0.0;
  for (PetscInt i = 0; i < num_clusters; ++i)
  {
    PetscReal effective = net_pressure + shadow[i];
    if (effective > 0.0)
    {
      sum_inv += 1.0 / effective;
    }
  }

  if (sum_inv <= 0.0)
  {
    return 1.0 / static_cast<PetscReal>(num_clusters);
  }

  PetscReal effective_idx = net_pressure + shadow[cluster_index];
  if (effective_idx <= 0.0)
  {
    return 0.0;
  }
  return (1.0 / effective_idx) / sum_inv;
}

// --- Phase 7: Induced seismicity ---

PetscReal PetscFEHydrofrac::scalarMoment(
    PetscReal shear_modulus,
    PetscReal average_slip,
    PetscReal rupture_area)
{
  if (shear_modulus <= 0.0 || average_slip <= 0.0 || rupture_area <= 0.0)
  {
    return 0.0;
  }
  return shear_modulus * average_slip * rupture_area;
}

PetscReal PetscFEHydrofrac::momentMagnitude(PetscReal scalar_moment)
{
  if (scalar_moment <= 0.0)
  {
    return -99.0;
  }
  // Hanks and Kanamori (1979): Mw = (2/3) * log10(M0) - 6.07
  return (2.0 / 3.0) * PetscLog10Real(scalar_moment) - 6.07;
}

PetscBool PetscFEHydrofrac::isMicroseismicRange(PetscReal mw)
{
  // Typical microseismic: -3 < Mw < 1
  return (mw > -3.0 && mw < 1.0) ? PETSC_TRUE : PETSC_FALSE;
}

// --- Phase 6: Proppant transport ---

PetscReal PetscFEHydrofrac::stokesSettlingVelocity(
    PetscReal particle_diameter,
    PetscReal particle_density,
    PetscReal fluid_density,
    PetscReal fluid_viscosity)
{
  if (particle_diameter <= 0.0 || fluid_viscosity <= 0.0)
  {
    return 0.0;
  }
  const PetscReal g = 9.81;
  const PetscReal delta_rho = particle_density - fluid_density;
  // v_s = d^2 * (rho_p - rho_f) * g / (18 * mu)
  return particle_diameter * particle_diameter * delta_rho * g /
         (18.0 * fluid_viscosity);
}

PetscBool PetscFEHydrofrac::proppantBridging(
    PetscReal aperture,
    PetscReal particle_diameter,
    PetscReal bridging_ratio)
{
  if (particle_diameter <= 0.0 || aperture <= 0.0 || bridging_ratio <= 0.0)
  {
    return PETSC_FALSE;
  }
  // Bridging occurs when aperture/diameter < bridging_ratio
  return (aperture / particle_diameter) < bridging_ratio ? PETSC_TRUE : PETSC_FALSE;
}

PetscReal PetscFEHydrofrac::proppantPackMinAperture(
    PetscReal particle_diameter,
    PetscReal pack_porosity)
{
  if (particle_diameter <= 0.0 || pack_porosity < 0.0 || pack_porosity >= 1.0)
  {
    return 0.0;
  }
  // Minimum aperture is one particle diameter (single layer pack)
  return particle_diameter;
}

PetscReal PetscFEHydrofrac::proppantMassBalance(
    PetscReal concentration_in,
    PetscReal flow_rate,
    PetscReal injection_time,
    PetscReal fracture_volume,
    PetscReal average_concentration)
{
  if (flow_rate <= 0.0 || injection_time <= 0.0 || fracture_volume <= 0.0)
  {
    return -1.0;
  }
  // Mass injected = c_in * Q * t
  // Mass in fracture = c_avg * V_f
  // Balance ratio should be ~1.0 for perfect conservation
  const PetscReal mass_in = concentration_in * flow_rate * injection_time;
  const PetscReal mass_frac = average_concentration * fracture_volume;
  if (mass_in <= 0.0)
  {
    return 0.0;
  }
  return mass_frac / mass_in;
}

// --- Phase 8: Carter leak-off ---

PetscReal PetscFEHydrofrac::carterLeakoffRate(
    PetscReal c_l,
    PetscReal time,
    PetscReal time_open)
{
  if (c_l <= 0.0 || time <= time_open)
  {
    return 0.0;
  }
  const PetscReal dt = time - time_open;
  return c_l / PetscSqrtReal(dt);
}

PetscReal PetscFEHydrofrac::carterCumulativeLeakoff(
    PetscReal c_l,
    PetscReal time,
    PetscReal time_open)
{
  if (c_l <= 0.0 || time <= time_open)
  {
    return 0.0;
  }
  const PetscReal dt = time - time_open;
  return 2.0 * c_l * PetscSqrtReal(dt);
}

PetscReal PetscFEHydrofrac::totalLeakoffVolume(
    PetscReal c_l,
    PetscReal time,
    PetscReal time_open,
    PetscReal face_area)
{
  if (face_area <= 0.0)
  {
    return 0.0;
  }
  return carterCumulativeLeakoff(c_l, time, time_open) * face_area;
}

// --- Phase 9: Production forecasting ---

PetscReal PetscFEHydrofrac::arpsDeclineRate(
    PetscReal q_initial,
    PetscReal decline_rate,
    PetscReal b_factor,
    PetscReal time)
{
  if (q_initial <= 0.0 || decline_rate <= 0.0 || time < 0.0)
  {
    return 0.0;
  }
  if (b_factor <= 0.0)
  {
    // Exponential decline: q = q_i * exp(-D_i * t)
    return q_initial * PetscExpReal(-decline_rate * time);
  }
  // Hyperbolic decline: q = q_i / (1 + b*D_i*t)^(1/b)
  const PetscReal denom = 1.0 + b_factor * decline_rate * time;
  if (denom <= 0.0)
  {
    return 0.0;
  }
  return q_initial / PetscPowReal(denom, 1.0 / b_factor);
}

PetscReal PetscFEHydrofrac::arpsCumulativeProduction(
    PetscReal q_initial,
    PetscReal decline_rate,
    PetscReal b_factor,
    PetscReal time)
{
  if (q_initial <= 0.0 || decline_rate <= 0.0 || time <= 0.0)
  {
    return 0.0;
  }
  if (b_factor <= 0.0)
  {
    // Exponential: Np = (q_i / D_i) * (1 - exp(-D_i*t))
    return (q_initial / decline_rate) * (1.0 - PetscExpReal(-decline_rate * time));
  }
  if (PetscAbsReal(b_factor - 1.0) < 1.0e-12)
  {
    // Harmonic: Np = (q_i / D_i) * ln(1 + D_i*t)
    return (q_initial / decline_rate) * PetscLogReal(1.0 + decline_rate * time);
  }
  // Hyperbolic: Np = q_i / ((1-b)*D_i) * [1 - (1+b*D_i*t)^((b-1)/b)]
  const PetscReal exponent = (b_factor - 1.0) / b_factor;
  const PetscReal base = 1.0 + b_factor * decline_rate * time;
  return q_initial / ((1.0 - b_factor) * decline_rate) *
         (1.0 - PetscPowReal(base, exponent));
}

PetscReal PetscFEHydrofrac::fractureProductivityIndex(
    PetscReal fracture_permeability,
    PetscReal fracture_width,
    PetscReal fluid_viscosity,
    PetscReal drainage_radius,
    PetscReal wellbore_radius)
{
  if (fracture_permeability <= 0.0 || fracture_width <= 0.0 ||
      fluid_viscosity <= 0.0 || drainage_radius <= 0.0 ||
      wellbore_radius <= 0.0 || drainage_radius <= wellbore_radius)
  {
    return 0.0;
  }
  return fracture_permeability * fracture_width /
         (fluid_viscosity * PetscLogReal(drainage_radius / wellbore_radius));
}

void PetscFEHydrofrac::f0_lagrange_regularize(
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
    PetscScalar f[])
{
  (void)Nf; (void)NfAux; (void)uOff_x; (void)u_t; (void)u_x;
  (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
  (void)t; (void)x; (void)numConstants; (void)constants;

  // Identity: F_lambda = lambda (makes system non-singular away from fault)
  // Lagrange field is always the last solution field before thermal
  const PetscInt lagr_off = uOff[Nf - 1];
  for (PetscInt d = 0; d < dim; ++d)
  {
    f[d] = u[lagr_off + d];
  }
}

void PetscFEHydrofrac::g0_lagrange_regularize(
    PetscInt dim, PetscInt Nf, PetscInt NfAux,
    const PetscInt uOff[], const PetscInt uOff_x[],
    const PetscScalar u[], const PetscScalar u_t[],
    const PetscScalar u_x[],
    const PetscInt aOff[], const PetscInt aOff_x[],
    const PetscScalar a[], const PetscScalar a_t[],
    const PetscScalar a_x[],
    PetscReal t, PetscReal u_tShift, const PetscReal x[],
    PetscInt numConstants,
    const PetscScalar constants[],
    PetscScalar g0[])
{
  (void)Nf; (void)NfAux; (void)uOff; (void)uOff_x; (void)u; (void)u_t;
  (void)u_x; (void)aOff; (void)aOff_x; (void)a; (void)a_t; (void)a_x;
  (void)t; (void)u_tShift; (void)x; (void)numConstants; (void)constants;

  // Identity Jacobian for the (Lagrange, Lagrange) block
  for (PetscInt d = 0; d < dim; ++d)
  {
    g0[d * dim + d] = 1.0;
  }
}

} // namespace FSRM
