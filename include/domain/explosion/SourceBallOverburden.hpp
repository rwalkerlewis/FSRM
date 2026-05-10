/**
 * @file SourceBallOverburden.hpp
 * @brief Pass-13b (axis-1b physics) asymmetric overburden initial state
 *        for the Source3DBall solver.
 *
 * Pure function (no PETSc dependency, easy to unit test): given a
 * cell-centroid Cartesian coordinate, the host-rock density, the
 * gravitational acceleration g, the source-emplacement depth (positive
 * downward from the free surface, so a 500 m deep shot has source_depth
 * = 500 m), and the at-rest Earth coefficient K_0, return the 6-Voigt
 * stress (xx, yy, zz, xy, xz, yz) at that cell.
 *
 * Sign convention (matches the rest of FSRM): compressive stress is
 * negative.
 *
 * Coordinate convention: the source-ball mesh is centred at the source
 * point, so a cell at z = +5 m relative to the source is 5 m shallower
 * (smaller depth) than the source. Translating to absolute depth d:
 *   d = source_depth - z_centroid_local
 * with z positive upward in the local source-ball frame.
 *
 * Stress tensor:
 *   sigma_zz(d) = -rho * g * d
 *   sigma_xx(d) = sigma_yy(d) = K_0 * sigma_zz(d)
 *   off-diagonal = 0
 *
 * Reference: Hoek & Brown (1980), "Underground excavations in rock"
 * for the at-rest Earth coefficient K_0; Patton (1991) for the cavity-
 * asymmetry consequence under partial K_0 loading.
 */

#ifndef NEAR_FIELD_SOURCE_BALL_OVERBURDEN_HPP
#define NEAR_FIELD_SOURCE_BALL_OVERBURDEN_HPP

#include <array>

namespace FSRM {

/// Configuration for one asymmetric-overburden IC application.
struct SourceBallOverburdenConfig
{
    /// Absolute depth of the source point [m], positive downward from
    /// the free surface. Cell centroid local z is added to give the
    /// per-cell absolute depth.
    double source_depth_m = 500.0;
    /// Host-rock density [kg/m^3].
    double rho_solid = 2700.0;
    /// Gravitational acceleration [m/s^2]. Use 9.81 unless you really
    /// mean otherwise.
    double g = 9.81;
    /// At-rest Earth coefficient. K_0 = 1 -> isotropic (lithostatic).
    /// K_0 < 1 -> compressive vertical, less compressive horizontal,
    /// the realistic crustal case (Hoek & Brown 1980 default 0.5).
    double K_0 = 0.5;
};

/// Compute the per-cell stress tensor for the asymmetric overburden IC.
/// Inputs:
///   centroid_xyz: local Cartesian centroid in the source-ball frame [m].
///   cfg: overburden configuration.
/// Returns the 6-Voigt stress at the centroid (sign: compression < 0).
inline std::array<double, 6> overburdenStressAtCell(
    const std::array<double, 3>& centroid_xyz,
    const SourceBallOverburdenConfig& cfg)
{
    // Local frame: z positive up. Absolute depth d = source_depth - z.
    const double d = cfg.source_depth_m - centroid_xyz[2];
    const double sigma_zz = -cfg.rho_solid * cfg.g * d;
    const double sigma_xx = cfg.K_0 * sigma_zz;
    const double sigma_yy = cfg.K_0 * sigma_zz;
    return {sigma_xx, sigma_yy, sigma_zz, 0.0, 0.0, 0.0};
}

}  // namespace FSRM

#endif  // NEAR_FIELD_SOURCE_BALL_OVERBURDEN_HPP
