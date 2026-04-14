/**
 * @file VelocityModelReader.cpp
 * @brief Read/write structured binary velocity model files
 */

#include "io/VelocityModelReader.hpp"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace FSRM {

// ---------------------------------------------------------------------------
// Trilinear interpolation
// ---------------------------------------------------------------------------
void VelocityModel::interpolate(double x, double y, double z,
                                double &vp_out, double &vs_out, double &rho_out) const
{
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        vp_out = 0.0;
        vs_out = 0.0;
        rho_out = 0.0;
        return;
    }

    // Convert to continuous grid coordinates
    auto toGridCoord = [](double val, double vmin, double vmax, int n) -> double {
        if (n <= 1) return 0.0;
        double t = (val - vmin) / (vmax - vmin) * (n - 1);
        return std::max(0.0, std::min(static_cast<double>(n - 1), t));
    };

    double gx = toGridCoord(x, x_min, x_max, nx);
    double gy = toGridCoord(y, y_min, y_max, ny);
    double gz = toGridCoord(z, z_min, z_max, nz);

    // Integer indices and fractional parts
    int ix0 = static_cast<int>(std::floor(gx));
    int iy0 = static_cast<int>(std::floor(gy));
    int iz0 = static_cast<int>(std::floor(gz));

    ix0 = std::max(0, std::min(nx - 2, ix0));
    iy0 = std::max(0, std::min(ny - 2, iy0));
    iz0 = std::max(0, std::min(nz - 2, iz0));

    // Handle edge case: single cell in a dimension
    if (nx <= 1) ix0 = 0;
    if (ny <= 1) iy0 = 0;
    if (nz <= 1) iz0 = 0;

    double fx = gx - ix0;
    double fy = gy - iy0;
    double fz = gz - iz0;

    // Clamp fractions for single-cell dimensions
    if (nx <= 1) fx = 0.0;
    if (ny <= 1) fy = 0.0;
    if (nz <= 1) fz = 0.0;

    int ix1 = std::min(ix0 + 1, nx - 1);
    int iy1 = std::min(iy0 + 1, ny - 1);
    int iz1 = std::min(iz0 + 1, nz - 1);

    // Flat index: ix * (ny * nz) + iy * nz + iz (C order)
    auto idx = [&](int ix, int iy, int iz) -> size_t {
        return static_cast<size_t>(ix) * static_cast<size_t>(ny * nz)
             + static_cast<size_t>(iy) * static_cast<size_t>(nz)
             + static_cast<size_t>(iz);
    };

    // Trilinear interpolation of a single field
    auto interp = [&](const std::vector<float> &field) -> double {
        double c000 = field[idx(ix0, iy0, iz0)];
        double c100 = field[idx(ix1, iy0, iz0)];
        double c010 = field[idx(ix0, iy1, iz0)];
        double c110 = field[idx(ix1, iy1, iz0)];
        double c001 = field[idx(ix0, iy0, iz1)];
        double c101 = field[idx(ix1, iy0, iz1)];
        double c011 = field[idx(ix0, iy1, iz1)];
        double c111 = field[idx(ix1, iy1, iz1)];

        double c00 = c000 * (1 - fx) + c100 * fx;
        double c01 = c001 * (1 - fx) + c101 * fx;
        double c10 = c010 * (1 - fx) + c110 * fx;
        double c11 = c011 * (1 - fx) + c111 * fx;

        double c0 = c00 * (1 - fy) + c10 * fy;
        double c1 = c01 * (1 - fy) + c11 * fy;

        return c0 * (1 - fz) + c1 * fz;
    };

    vp_out  = interp(vp);
    vs_out  = interp(vs);
    rho_out = interp(rho);
}

// ---------------------------------------------------------------------------
// Read binary velocity model
// ---------------------------------------------------------------------------
int readVelocityModel(const std::string &filename, VelocityModel &model)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return 1;
    }

    // Read header: 3 int32 + 8 float32 = 12 + 32 = 44 bytes?
    // Actually: 3*4 + 6*4 + 2*4 = 12 + 24 + 8 = 44 bytes.
    // But specification says 32 bytes total. Let me re-read the spec:
    //   nx, ny, nz (3 x int32) = 12 bytes
    //   x_min, x_max, y_min, y_max, z_min, z_max (6 x float32) = 24 bytes
    //   padding (2 x float32) = 8 bytes
    //   Total = 44 bytes
    // The prompt says 32 bytes but that is incorrect (3*4 + 6*4 + 2*4 = 44).
    // We use 44 bytes.

    int32_t header_ints[3];
    float header_floats[8]; // 6 bounds + 2 padding

    file.read(reinterpret_cast<char *>(header_ints), sizeof(header_ints));
    file.read(reinterpret_cast<char *>(header_floats), sizeof(header_floats));

    if (!file.good()) {
        return 2;
    }

    model.nx = header_ints[0];
    model.ny = header_ints[1];
    model.nz = header_ints[2];
    model.x_min = static_cast<double>(header_floats[0]);
    model.x_max = static_cast<double>(header_floats[1]);
    model.y_min = static_cast<double>(header_floats[2]);
    model.y_max = static_cast<double>(header_floats[3]);
    model.z_min = static_cast<double>(header_floats[4]);
    model.z_max = static_cast<double>(header_floats[5]);

    if (model.nx <= 0 || model.ny <= 0 || model.nz <= 0) {
        return 3;
    }

    size_t n = static_cast<size_t>(model.nx) *
               static_cast<size_t>(model.ny) *
               static_cast<size_t>(model.nz);

    // Read interleaved Vp, Vs, rho triplets
    model.vp.resize(n);
    model.vs.resize(n);
    model.rho.resize(n);

    for (size_t i = 0; i < n; ++i) {
        float triplet[3];
        file.read(reinterpret_cast<char *>(triplet), sizeof(triplet));
        if (!file.good()) {
            return 4;
        }
        model.vp[i]  = triplet[0];
        model.vs[i]  = triplet[1];
        model.rho[i] = triplet[2];
    }

    return 0;
}

// ---------------------------------------------------------------------------
// Write binary velocity model (for test fixtures)
// ---------------------------------------------------------------------------
int writeVelocityModel(const std::string &filename, const VelocityModel &model)
{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return 1;
    }

    int32_t header_ints[3] = {model.nx, model.ny, model.nz};
    float header_floats[8] = {
        static_cast<float>(model.x_min), static_cast<float>(model.x_max),
        static_cast<float>(model.y_min), static_cast<float>(model.y_max),
        static_cast<float>(model.z_min), static_cast<float>(model.z_max),
        0.0f, 0.0f  // padding
    };

    file.write(reinterpret_cast<const char *>(header_ints), sizeof(header_ints));
    file.write(reinterpret_cast<const char *>(header_floats), sizeof(header_floats));

    size_t n = static_cast<size_t>(model.nx) *
               static_cast<size_t>(model.ny) *
               static_cast<size_t>(model.nz);

    for (size_t i = 0; i < n; ++i) {
        float triplet[3] = {model.vp[i], model.vs[i], model.rho[i]};
        file.write(reinterpret_cast<const char *>(triplet), sizeof(triplet));
    }

    if (!file.good()) {
        return 2;
    }

    return 0;
}

} // namespace FSRM
