/**
 * @file TabulatedDataReader.cpp
 * @brief HDF5 ingestion + log-space bilinear interpolation for the
 *        pass-9 tabulated EOS / opacity patches.
 *
 * Mirrors the pattern used in src/domain/explosion/RadialLagrangianOutput.cpp:
 * the HDF5 C API is included via __has_include and a no-op fallback
 * is provided when HDF5 is absent at compile time. The reader is then
 * effectively unusable in HDF5-less builds, but the rest of the code
 * compiles. At runtime an attempted load() returns false with a clear
 * error_message.
 */

#include "io/TabulatedData/TabulatedDataReader.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>

#if defined(__has_include)
#  if __has_include(<hdf5.h>)
#    include <hdf5.h>
#    define FSRM_HAVE_HDF5_H 1
#  endif
#endif

namespace FSRM {
namespace io {

#ifdef FSRM_HAVE_HDF5_H

namespace {

bool readScalarString(hid_t loc, const std::string& name, std::string& out)
{
    if (!H5Lexists(loc, name.c_str(), H5P_DEFAULT)) return false;
    hid_t dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
    if (dset < 0) return false;
    hid_t dtype = H5Dget_type(dset);
    H5T_class_t cls = H5Tget_class(dtype);
    if (cls != H5T_STRING) {
        H5Tclose(dtype);
        H5Dclose(dset);
        return false;
    }
    if (H5Tis_variable_str(dtype)) {
        // Variable-length string. h5py writes these with UTF-8 char
        // set; read with the file's native dtype to avoid the "no
        // datatype conversion path" error that the simple
        // H5Tcopy(H5T_C_S1) + H5Tset_size(VARIABLE) recipe triggers
        // when the file is UTF-8 encoded but the memory type is
        // ASCII.
        char* str_buf = nullptr;
        herr_t st = H5Dread(dset, dtype, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, &str_buf);
        if (st >= 0 && str_buf) {
            out = str_buf;
            // Reclaim the variable-length buffer that HDF5 allocated.
            hid_t space = H5Dget_space(dset);
            H5Dvlen_reclaim(dtype, space, H5P_DEFAULT, &str_buf);
            H5Sclose(space);
        }
        H5Tclose(dtype);
        H5Dclose(dset);
        return st >= 0;
    } else {
        size_t sz = H5Tget_size(dtype);
        std::vector<char> buf(sz + 1, 0);
        herr_t st = H5Dread(dset, dtype, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, buf.data());
        if (st >= 0) out.assign(buf.data());
        H5Tclose(dtype);
        H5Dclose(dset);
        return st >= 0;
    }
}

bool readScalarDouble(hid_t loc, const std::string& name, double& out)
{
    if (!H5Lexists(loc, name.c_str(), H5P_DEFAULT)) return false;
    hid_t dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
    if (dset < 0) return false;
    herr_t st = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, &out);
    H5Dclose(dset);
    return st >= 0;
}

bool readDataset1D(hid_t loc, const std::string& name,
                   std::vector<double>& out)
{
    if (!H5Lexists(loc, name.c_str(), H5P_DEFAULT)) return false;
    hid_t dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
    if (dset < 0) return false;
    hid_t space = H5Dget_space(dset);
    int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 1) {
        H5Sclose(space);
        H5Dclose(dset);
        return false;
    }
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(space, dims, nullptr);
    out.resize(dims[0]);
    herr_t st = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, out.data());
    H5Sclose(space);
    H5Dclose(dset);
    return st >= 0;
}

bool readDataset2D(hid_t loc, const std::string& name,
                   std::vector<double>& out, hsize_t& dim0, hsize_t& dim1)
{
    if (!H5Lexists(loc, name.c_str(), H5P_DEFAULT)) return false;
    hid_t dset = H5Dopen2(loc, name.c_str(), H5P_DEFAULT);
    if (dset < 0) return false;
    hid_t space = H5Dget_space(dset);
    int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 2) {
        H5Sclose(space);
        H5Dclose(dset);
        return false;
    }
    hsize_t dims[2] = {0, 0};
    H5Sget_simple_extent_dims(space, dims, nullptr);
    out.resize(dims[0] * dims[1]);
    herr_t st = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, out.data());
    H5Sclose(space);
    H5Dclose(dset);
    dim0 = dims[0];
    dim1 = dims[1];
    return st >= 0;
}

TabulatedQuantity parseQuantity(const std::string& s)
{
    if (s == "EOS_PRESSURE") return TabulatedQuantity::EOS_PRESSURE;
    if (s == "EOS_SOUND_SPEED") return TabulatedQuantity::EOS_SOUND_SPEED;
    if (s == "OPACITY_ROSSELAND") return TabulatedQuantity::OPACITY_ROSSELAND;
    if (s == "OPACITY_PLANCK") return TabulatedQuantity::OPACITY_PLANCK;
    return TabulatedQuantity::UNKNOWN;
}

bool isOpacityQuantity(TabulatedQuantity q)
{
    return q == TabulatedQuantity::OPACITY_ROSSELAND ||
           q == TabulatedQuantity::OPACITY_PLANCK;
}

}  // namespace

bool TabulatedDataReader::load(const std::string& path,
                               std::string& error_message)
{
    loaded_ = false;
    error_message.clear();

    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        error_message = "TabulatedDataReader: cannot open " + path;
        return false;
    }

    hid_t meta = H5Gopen2(file, "metadata", H5P_DEFAULT);
    if (meta < 0) {
        error_message = "TabulatedDataReader: missing /metadata in " + path;
        H5Fclose(file);
        return false;
    }

    if (!readScalarString(meta, "medium", meta_.medium)) {
        error_message = "TabulatedDataReader: missing /metadata/medium";
        H5Gclose(meta);
        H5Fclose(file);
        return false;
    }
    if (!readScalarString(meta, "quantity", meta_.quantity_str)) {
        error_message = "TabulatedDataReader: missing /metadata/quantity";
        H5Gclose(meta);
        H5Fclose(file);
        return false;
    }
    meta_.quantity = parseQuantity(meta_.quantity_str);
    if (meta_.quantity == TabulatedQuantity::UNKNOWN) {
        error_message = "TabulatedDataReader: unknown quantity '" +
                        meta_.quantity_str + "'";
        H5Gclose(meta);
        H5Fclose(file);
        return false;
    }
    if (!readScalarString(meta, "source_citation", meta_.source_citation) ||
        meta_.source_citation.empty()) {
        error_message =
            "TabulatedDataReader: source_citation is REQUIRED for " + path +
            " (every shipped table must cite primary literature)";
        H5Gclose(meta);
        H5Fclose(file);
        return false;
    }
    readScalarString(meta, "generation_date", meta_.generation_date);
    readScalarString(meta, "generation_tool", meta_.generation_tool);

    if (!readDataset1D(meta, "rho_axis_kg_per_m3", rho_axis_)) {
        error_message =
            "TabulatedDataReader: missing /metadata/rho_axis_kg_per_m3";
        H5Gclose(meta);
        H5Fclose(file);
        return false;
    }

    const bool is_opacity = isOpacityQuantity(meta_.quantity);
    const std::string other_name = is_opacity ? "T_axis_K" : "e_axis_J_per_kg";
    if (!readDataset1D(meta, other_name, other_axis_)) {
        error_message =
            "TabulatedDataReader: missing /metadata/" + other_name;
        H5Gclose(meta);
        H5Fclose(file);
        return false;
    }

    readScalarDouble(meta, "rho_min", meta_.rho_min);
    readScalarDouble(meta, "rho_max", meta_.rho_max);
    if (is_opacity) {
        readScalarDouble(meta, "T_min", meta_.other_min);
        readScalarDouble(meta, "T_max", meta_.other_max);
    } else {
        readScalarDouble(meta, "e_min", meta_.other_min);
        readScalarDouble(meta, "e_max", meta_.other_max);
    }
    H5Gclose(meta);

    hsize_t d0 = 0, d1 = 0;
    if (!readDataset2D(file, "data", data_, d0, d1)) {
        error_message = "TabulatedDataReader: missing /data dataset";
        H5Fclose(file);
        return false;
    }
    H5Fclose(file);

    n_rho_ = static_cast<int>(rho_axis_.size());
    n_other_ = static_cast<int>(other_axis_.size());
    if (n_rho_ < 2 || n_other_ < 2) {
        error_message =
            "TabulatedDataReader: axes must have at least 2 points each";
        return false;
    }

    // Detect axis ordering of the 2D dataset. RHO_MAJOR if dims =
    // (n_rho, n_other), RHO_MINOR if dims = (n_other, n_rho).
    if (static_cast<hsize_t>(n_rho_) == d0 &&
        static_cast<hsize_t>(n_other_) == d1) {
        ordering_ = AxisOrdering::RHO_MAJOR;
    } else if (static_cast<hsize_t>(n_rho_) == d1 &&
               static_cast<hsize_t>(n_other_) == d0) {
        ordering_ = AxisOrdering::RHO_MINOR;
    } else {
        std::ostringstream m;
        m << "TabulatedDataReader: /data shape (" << d0 << ", " << d1
          << ") does not match axis sizes rho=" << n_rho_
          << ", other=" << n_other_;
        error_message = m.str();
        return false;
    }

    // Sanity-check ascending axes.
    for (int i = 1; i < n_rho_; ++i) {
        if (rho_axis_[i] <= rho_axis_[i - 1]) {
            error_message =
                "TabulatedDataReader: rho_axis must be strictly ascending";
            return false;
        }
    }
    for (int i = 1; i < n_other_; ++i) {
        if (other_axis_[i] <= other_axis_[i - 1]) {
            error_message =
                "TabulatedDataReader: other axis must be strictly ascending";
            return false;
        }
    }

    log_rho_axis_.resize(n_rho_);
    for (int i = 0; i < n_rho_; ++i) log_rho_axis_[i] = logSafe(rho_axis_[i]);
    log_other_axis_.resize(n_other_);
    for (int i = 0; i < n_other_; ++i)
        log_other_axis_[i] = logSafe(other_axis_[i]);

    if (meta_.rho_min == 0.0) meta_.rho_min = rho_axis_.front();
    if (meta_.rho_max == 0.0) meta_.rho_max = rho_axis_.back();
    if (meta_.other_min == 0.0) meta_.other_min = other_axis_.front();
    if (meta_.other_max == 0.0) meta_.other_max = other_axis_.back();

    loaded_ = true;
    oor_warned_low_rho_ = false;
    oor_warned_high_rho_ = false;
    oor_warned_low_other_ = false;
    oor_warned_high_other_ = false;
    return true;
}

#else  // !FSRM_HAVE_HDF5_H

bool TabulatedDataReader::load(const std::string& path,
                               std::string& error_message)
{
    (void)path;
    error_message =
        "TabulatedDataReader: HDF5 not compiled in; tabulated data path "
        "is unavailable in this build";
    loaded_ = false;
    return false;
}

#endif  // FSRM_HAVE_HDF5_H

int TabulatedDataReader::index2(int i_rho, int j_other) const
{
    if (ordering_ == AxisOrdering::RHO_MAJOR) {
        return i_rho * n_other_ + j_other;
    }
    return j_other * n_rho_ + i_rho;
}

double TabulatedDataReader::dataAt(int i_rho, int j_other) const
{
    if (!loaded_) return TABULATED_DATA_OOR_SENTINEL;
    if (i_rho < 0 || i_rho >= n_rho_) return TABULATED_DATA_OOR_SENTINEL;
    if (j_other < 0 || j_other >= n_other_)
        return TABULATED_DATA_OOR_SENTINEL;
    return data_[index2(i_rho, j_other)];
}

void TabulatedDataReader::warnOnce(const char* axis_label, bool& flag,
                                   double queried, double bound,
                                   int cell_index) const
{
    if (flag) return;
    flag = true;
    std::cerr << "TabulatedDataReader[" << meta_.medium << "/"
              << meta_.quantity_str << "]: " << axis_label
              << " out-of-range (queried " << queried << " vs bound "
              << bound << ", cell_index=" << cell_index
              << "). Falling back to analytic model. This warning is "
                 "logged once per (medium, quantity, axis).\n";
}

double TabulatedDataReader::evaluate(double rho, double other,
                                     int cell_index) const
{
    if (!loaded_) return TABULATED_DATA_OOR_SENTINEL;
    if (rho < rho_axis_.front()) {
        warnOnce("rho_low", oor_warned_low_rho_, rho, rho_axis_.front(),
                 cell_index);
        return TABULATED_DATA_OOR_SENTINEL;
    }
    if (rho > rho_axis_.back()) {
        warnOnce("rho_high", oor_warned_high_rho_, rho, rho_axis_.back(),
                 cell_index);
        return TABULATED_DATA_OOR_SENTINEL;
    }
    if (other < other_axis_.front()) {
        warnOnce("other_low", oor_warned_low_other_, other,
                 other_axis_.front(), cell_index);
        return TABULATED_DATA_OOR_SENTINEL;
    }
    if (other > other_axis_.back()) {
        warnOnce("other_high", oor_warned_high_other_, other,
                 other_axis_.back(), cell_index);
        return TABULATED_DATA_OOR_SENTINEL;
    }

    const double lr = logSafe(rho);
    const double lo = logSafe(other);

    auto it_r = std::upper_bound(log_rho_axis_.begin(),
                                 log_rho_axis_.end(), lr);
    int i1 = static_cast<int>(it_r - log_rho_axis_.begin());
    if (i1 <= 0) i1 = 1;
    if (i1 >= n_rho_) i1 = n_rho_ - 1;
    const int i0 = i1 - 1;

    auto it_o = std::upper_bound(log_other_axis_.begin(),
                                 log_other_axis_.end(), lo);
    int j1 = static_cast<int>(it_o - log_other_axis_.begin());
    if (j1 <= 0) j1 = 1;
    if (j1 >= n_other_) j1 = n_other_ - 1;
    const int j0 = j1 - 1;

    const double t_r = (lr - log_rho_axis_[i0]) /
                       (log_rho_axis_[i1] - log_rho_axis_[i0]);
    const double t_o = (lo - log_other_axis_[j0]) /
                       (log_other_axis_[j1] - log_other_axis_[j0]);

    const double f00 = data_[index2(i0, j0)];
    const double f10 = data_[index2(i1, j0)];
    const double f01 = data_[index2(i0, j1)];
    const double f11 = data_[index2(i1, j1)];

    // Bilinear interpolation in log space on the value as well, but
    // only when all four corners are positive. For pressure tables
    // negative values (rare; sound-speed gradients in mixed regime)
    // make a log-value interpolation invalid; fall back to linear in
    // value space in that case.
    const bool all_positive =
        (f00 > 0.0) && (f10 > 0.0) && (f01 > 0.0) && (f11 > 0.0);
    if (all_positive) {
        const double lf00 = std::log(f00);
        const double lf10 = std::log(f10);
        const double lf01 = std::log(f01);
        const double lf11 = std::log(f11);
        const double lv = (1.0 - t_r) * (1.0 - t_o) * lf00 +
                          t_r * (1.0 - t_o) * lf10 +
                          (1.0 - t_r) * t_o * lf01 +
                          t_r * t_o * lf11;
        return std::exp(lv);
    }
    return (1.0 - t_r) * (1.0 - t_o) * f00 +
           t_r * (1.0 - t_o) * f10 +
           (1.0 - t_r) * t_o * f01 +
           t_r * t_o * f11;
}

bool TabulatedDataReader::evaluateDerivatives(double rho, double other,
                                              double& dval_drho,
                                              double& dval_dother,
                                              int cell_index) const
{
    dval_drho = 0.0;
    dval_dother = 0.0;
    if (!loaded_) return false;
    if (rho < rho_axis_.front() || rho > rho_axis_.back() ||
        other < other_axis_.front() || other > other_axis_.back()) {
        return false;
    }

    // Centred finite differences on the bilinear-interpolated value.
    const double drho = std::max(rho * 1.0e-3, 1.0e-12);
    const double dother = std::max(other * 1.0e-3, 1.0e-12);

    const double rp = std::min(rho + drho, rho_axis_.back());
    const double rm = std::max(rho - drho, rho_axis_.front());
    const double op = std::min(other + dother, other_axis_.back());
    const double om = std::max(other - dother, other_axis_.front());

    const double f_rp = evaluate(rp, other, cell_index);
    const double f_rm = evaluate(rm, other, cell_index);
    const double f_op = evaluate(rho, op, cell_index);
    const double f_om = evaluate(rho, om, cell_index);

    if (std::isnan(f_rp) || std::isnan(f_rm) ||
        std::isnan(f_op) || std::isnan(f_om)) {
        return false;
    }
    dval_drho = (f_rp - f_rm) / (rp - rm);
    dval_dother = (f_op - f_om) / (op - om);
    return true;
}

} // namespace io
} // namespace FSRM
