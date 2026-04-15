#ifndef VELOCITY_MODEL_READER_HPP
#define VELOCITY_MODEL_READER_HPP

/**
 * @file VelocityModelReader.hpp
 * @brief Read structured binary velocity model files for per-cell material assignment
 *
 * Binary format (little-endian):
 *   Header (32 bytes):
 *     nx, ny, nz       (3 x int32)
 *     x_min, x_max     (2 x float32)
 *     y_min, y_max     (2 x float32)
 *     z_min, z_max     (2 x float32)
 *     padding           (2 x float32, reserved)
 *   Data (nx * ny * nz * 3 * float32):
 *     For each (ix, iy, iz) in C order (iz fastest):
 *       Vp, Vs, rho    (3 x float32)
 *
 * Use case: seismologists distribute 3D velocity models as structured binary
 * files. FSRM reads these and assigns per-cell material properties to the
 * unstructured FEM mesh via trilinear interpolation.
 */

#include <string>
#include <vector>
#include <cstdint>

namespace FSRM {

struct VelocityModel {
    int32_t nx = 0, ny = 0, nz = 0;
    double x_min = 0, x_max = 0;
    double y_min = 0, y_max = 0;
    double z_min = 0, z_max = 0;
    std::vector<float> vp;    // nx * ny * nz
    std::vector<float> vs;    // nx * ny * nz
    std::vector<float> rho;   // nx * ny * nz

    /**
     * @brief Trilinear interpolation of Vp, Vs, rho at point (x, y, z)
     *
     * Points outside the model domain are clamped to the nearest boundary value.
     *
     * @param[in]  x, y, z     Query coordinates
     * @param[out] vp_out      Interpolated P-wave velocity (m/s)
     * @param[out] vs_out      Interpolated S-wave velocity (m/s)
     * @param[out] rho_out     Interpolated density (kg/m^3)
     */
    void interpolate(double x, double y, double z,
                     double &vp_out, double &vs_out, double &rho_out) const;
};

/**
 * @brief Read a binary velocity model file
 * @param filename  Path to the binary file
 * @param model     Output model structure
 * @return 0 on success, nonzero on error
 */
int readVelocityModel(const std::string &filename, VelocityModel &model);

/**
 * @brief Write a binary velocity model file (for test fixtures)
 * @param filename  Path to the output binary file
 * @param model     Model to write
 * @return 0 on success, nonzero on error
 */
int writeVelocityModel(const std::string &filename, const VelocityModel &model);

} // namespace FSRM

#endif // VELOCITY_MODEL_READER_HPP
