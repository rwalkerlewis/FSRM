/**
 * @file RadialLagrangianOutput.hpp
 * @brief Pass-6 HDF5 / XDMF spatial-profile writer for the 1D radial
 *        Lagrangian shock solver.
 *
 * The Simulator records snapshots of the radial state from
 * RadialLagrangianSolver::getRadialProfile at the configured cadence
 * during the DYNAMIC_PLASTIC + RADIAL_LAGRANGIAN setup phase. After
 * the recording loop completes, writeProfilesHDF5 emits a single
 * near_field_profile.h5 file with one HDF5 group per snapshot, and
 * writeProfilesXDMF emits a sibling near_field_profile.xdmf wrapper
 * so ParaView can read the file as a time series.
 *
 * Schema (frozen for pass-6; consumers depend on it):
 *
 *   /time                 (n_snap,)  double, simulation time [s]
 *   /num_cells            scalar     int, snapshot N
 *   /profiles/<i>/r        (N+1,)    double, face radii [m]
 *   /profiles/<i>/r_cell   (N,)      double, cell-centred radii [m]
 *   /profiles/<i>/v_r      (N+1,)    double, face velocities [m/s]
 *   /profiles/<i>/rho      (N,)      double, density [kg/m^3]
 *   /profiles/<i>/p        (N,)      double, pressure [Pa]
 *   /profiles/<i>/sigma_rr (N,)      double, total radial stress [Pa]
 *   /profiles/<i>/sigma_tt (N,)      double, total hoop stress [Pa]
 *   /profiles/<i>/eps_p    (N,)      double, equivalent plastic strain
 *   /profiles/<i>/damage   (N,)      double, scalar damage
 *   /profiles/<i>/yield_indicator (N,) double, 1.0 if yielded
 *
 * The XDMF references the cell-centred values and exposes them as
 * cell-data on a polyline mesh of N+1 vertices. ParaView's
 * XYChartView can plot the cell arrays vs r_cell as a time-series
 * line chart; the time-keeper drives an animation through the
 * cavity-formation transient.
 */

#ifndef NEAR_FIELD_RADIAL_LAGRANGIAN_OUTPUT_HPP
#define NEAR_FIELD_RADIAL_LAGRANGIAN_OUTPUT_HPP

#include <string>
#include <vector>

#include "domain/explosion/RadialLagrangian.hpp"

namespace FSRM {

/// Write the accumulated radial-profile snapshots to a single
/// near_field_profile.h5 file. Returns true on success. No-op (returns
/// true) if profiles is empty. Rank-0 only; the caller must guard.
bool writeRadialProfilesHDF5(const std::string& path,
                             const std::vector<RadialLagrangianSolver::RadialProfile>& profiles);

/// Write the matching near_field_profile.xdmf wrapper. References the
/// HDF5 file by relative path so ParaView can find it. Returns true on
/// success. No-op (returns true) if profiles is empty.
bool writeRadialProfilesXDMF(const std::string& xdmf_path,
                             const std::string& h5_basename,
                             const std::vector<RadialLagrangianSolver::RadialProfile>& profiles);

} // namespace FSRM

#endif // NEAR_FIELD_RADIAL_LAGRANGIAN_OUTPUT_HPP
