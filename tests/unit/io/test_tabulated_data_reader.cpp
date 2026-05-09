/**
 * @file test_tabulated_data_reader.cpp
 * @brief Pass-9 unit tests for io::TabulatedDataReader.
 *
 * Five gates:
 *   - RoundTripGranite: write a synthetic table, read it back,
 *     bilinear-interpolate, verify round-trip in log space.
 *   - OutOfRangeWarn: out-of-range query returns NaN sentinel; the
 *     reader logs a one-time warning. We capture the (medium,
 *     quantity) signature in the metadata.
 *   - AxisOrderingDetection: tables written rho-major and rho-minor
 *     both parse and produce identical evaluations.
 *   - CitationMetadata: a table without /metadata/source_citation is
 *     refused at load() time with a clear error.
 *   - LogspaceInterp: bilinear interp matches a known analytic
 *     function f(rho, e) = rho * e in log space within 1%.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "io/TabulatedData/TabulatedDataReader.hpp"

#if defined(__has_include)
#  if __has_include(<hdf5.h>)
#    include <hdf5.h>
#    define FSRM_TEST_HAVE_HDF5 1
#  endif
#endif

namespace fs = std::filesystem;
using FSRM::io::AxisOrdering;
using FSRM::io::TabulatedDataReader;
using FSRM::io::TabulatedQuantity;

namespace
{

#ifdef FSRM_TEST_HAVE_HDF5

void writeStringDataset(hid_t loc, const std::string& name,
                        const std::string& value)
{
    hid_t dt = H5Tcopy(H5T_C_S1);
    H5Tset_size(dt, value.size() + 1);
    H5Tset_strpad(dt, H5T_STR_NULLTERM);
    hid_t scalar = H5Screate(H5S_SCALAR);
    hid_t dset = H5Dcreate2(loc, name.c_str(), dt, scalar,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, dt, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.c_str());
    H5Dclose(dset);
    H5Sclose(scalar);
    H5Tclose(dt);
}

void writeDataset1D(hid_t loc, const std::string& name,
                    const std::vector<double>& v)
{
    hsize_t dims[1] = {v.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dset = H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE,
                            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, v.data());
    H5Dclose(dset);
    H5Sclose(space);
}

void writeScalarDouble(hid_t loc, const std::string& name, double v)
{
    hid_t scalar = H5Screate(H5S_SCALAR);
    hid_t dset = H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE,
                            scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, &v);
    H5Dclose(dset);
    H5Sclose(scalar);
}

void writeDataset2D(hid_t loc, const std::string& name,
                    const std::vector<double>& flat,
                    hsize_t d0, hsize_t d1)
{
    hsize_t dims[2] = {d0, d1};
    hid_t space = H5Screate_simple(2, dims, nullptr);
    hid_t dset = H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE,
                            space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, flat.data());
    H5Dclose(dset);
    H5Sclose(space);
}

struct SyntheticEOSTable
{
    std::vector<double> rho_axis;
    std::vector<double> e_axis;
    std::vector<double> data;  // (n_rho, n_e) row-major.

    static SyntheticEOSTable makeAnalytic(int n_rho = 32, int n_e = 32)
    {
        SyntheticEOSTable t;
        t.rho_axis.resize(n_rho);
        t.e_axis.resize(n_e);
        for (int i = 0; i < n_rho; ++i) {
            t.rho_axis[i] = std::pow(10.0, 2.0 + i * (4.0 - 2.0) / (n_rho - 1));
        }
        for (int j = 0; j < n_e; ++j) {
            t.e_axis[j] = std::pow(10.0, 3.0 + j * (9.0 - 3.0) / (n_e - 1));
        }
        t.data.resize(static_cast<size_t>(n_rho) * n_e);
        // Analytic: p(rho, e) = rho * e (linear in both, so log-space
        // bilinear is exact in the log-rho, log-e plane).
        for (int i = 0; i < n_rho; ++i) {
            for (int j = 0; j < n_e; ++j) {
                t.data[static_cast<size_t>(i) * n_e + j] =
                    t.rho_axis[i] * t.e_axis[j];
            }
        }
        return t;
    }
};

bool writeRhoMajorEOSTable(const std::string& path,
                           const SyntheticEOSTable& t,
                           const std::string& citation,
                           const std::string& medium = "GRANITE",
                           bool include_citation = true)
{
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return false;
    hid_t meta = H5Gcreate2(file, "metadata", H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
    writeStringDataset(meta, "medium", medium);
    writeStringDataset(meta, "quantity", "EOS_PRESSURE");
    if (include_citation) {
        writeStringDataset(meta, "source_citation", citation);
    }
    writeStringDataset(meta, "generation_date",
                       "2026-05-09T00:00:00Z");
    writeStringDataset(meta, "generation_tool", "test fixture");
    writeDataset1D(meta, "rho_axis_kg_per_m3", t.rho_axis);
    writeDataset1D(meta, "e_axis_J_per_kg", t.e_axis);
    writeScalarDouble(meta, "rho_min", t.rho_axis.front());
    writeScalarDouble(meta, "rho_max", t.rho_axis.back());
    writeScalarDouble(meta, "e_min", t.e_axis.front());
    writeScalarDouble(meta, "e_max", t.e_axis.back());
    writeDataset2D(file, "data", t.data,
                   t.rho_axis.size(), t.e_axis.size());
    H5Gclose(meta);
    H5Fclose(file);
    return true;
}

bool writeRhoMinorEOSTable(const std::string& path,
                           const SyntheticEOSTable& t,
                           const std::string& citation)
{
    // Same data, transposed so axis ordering is RHO_MINOR
    // (data[i_e, j_rho]).
    SyntheticEOSTable transposed = t;
    transposed.data.assign(t.data.size(), 0.0);
    const int n_rho = static_cast<int>(t.rho_axis.size());
    const int n_e = static_cast<int>(t.e_axis.size());
    for (int i = 0; i < n_rho; ++i) {
        for (int j = 0; j < n_e; ++j) {
            transposed.data[static_cast<size_t>(j) * n_rho + i] =
                t.data[static_cast<size_t>(i) * n_e + j];
        }
    }
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return false;
    hid_t meta = H5Gcreate2(file, "metadata", H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
    writeStringDataset(meta, "medium", "GRANITE");
    writeStringDataset(meta, "quantity", "EOS_PRESSURE");
    writeStringDataset(meta, "source_citation", citation);
    writeStringDataset(meta, "generation_date", "2026-05-09T00:00:00Z");
    writeStringDataset(meta, "generation_tool", "test fixture");
    writeDataset1D(meta, "rho_axis_kg_per_m3", t.rho_axis);
    writeDataset1D(meta, "e_axis_J_per_kg", t.e_axis);
    writeDataset2D(file, "data", transposed.data,
                   t.e_axis.size(), t.rho_axis.size());
    H5Gclose(meta);
    H5Fclose(file);
    return true;
}

#endif  // FSRM_TEST_HAVE_HDF5

}  // namespace


class TabulatedDataReaderTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
#ifndef FSRM_TEST_HAVE_HDF5
        GTEST_SKIP() << "HDF5 not available in this build; "
                        "TabulatedDataReader tests require it.";
#endif
        tmp_dir_ = fs::temp_directory_path() /
                   ("fsrm_tab_test_" + std::to_string(::getpid()));
        fs::create_directories(tmp_dir_);
    }
    void TearDown() override
    {
        if (fs::exists(tmp_dir_)) {
            std::error_code ec;
            fs::remove_all(tmp_dir_, ec);
        }
    }

    fs::path tmp_dir_;
};

TEST_F(TabulatedDataReaderTest, RoundTripGranite)
{
#ifdef FSRM_TEST_HAVE_HDF5
    const auto t = SyntheticEOSTable::makeAnalytic(32, 32);
    const std::string path = (tmp_dir_ / "round_trip.h5").string();
    ASSERT_TRUE(writeRhoMajorEOSTable(path, t,
        "Test citation: Marsh 1980 LASL Hugoniot."));

    TabulatedDataReader reader;
    std::string err;
    ASSERT_TRUE(reader.load(path, err)) << err;
    EXPECT_TRUE(reader.isLoaded());
    EXPECT_EQ(reader.metadata().medium, "GRANITE");
    EXPECT_EQ(reader.metadata().quantity,
              TabulatedQuantity::EOS_PRESSURE);
    EXPECT_EQ(reader.numRhoPoints(), 32);
    EXPECT_EQ(reader.numOtherPoints(), 32);

    // Round-trip on a grid corner: (rho=100, e=1e3) -> p = 1e5
    const double v_corner = reader.evaluate(100.0, 1.0e3);
    EXPECT_NEAR(v_corner, 100.0 * 1.0e3, 1.0e-3 * 100.0 * 1.0e3);

    // Mid-point: (rho=1e3, e=1e6) -> p = 1e9. Log-space bilinear on
    // f(rho, e) = rho*e is exact (within numerical precision) so the
    // tolerance is tight.
    const double v_mid = reader.evaluate(1.0e3, 1.0e6);
    EXPECT_NEAR(v_mid, 1.0e3 * 1.0e6, 1.0e-9 * 1.0e3 * 1.0e6);
#endif
}

TEST_F(TabulatedDataReaderTest, OutOfRangeWarn)
{
#ifdef FSRM_TEST_HAVE_HDF5
    const auto t = SyntheticEOSTable::makeAnalytic(16, 16);
    const std::string path = (tmp_dir_ / "oor.h5").string();
    ASSERT_TRUE(writeRhoMajorEOSTable(path, t, "test"));

    TabulatedDataReader reader;
    std::string err;
    ASSERT_TRUE(reader.load(path, err));

    const double below_rho = reader.evaluate(1.0, 1.0e6);
    EXPECT_TRUE(std::isnan(below_rho));

    const double above_e = reader.evaluate(1.0e3, 1.0e12);
    EXPECT_TRUE(std::isnan(above_e));
#endif
}

TEST_F(TabulatedDataReaderTest, AxisOrderingDetection)
{
#ifdef FSRM_TEST_HAVE_HDF5
    // Use distinct axis sizes (12 != 20) so the on-disk shape
    // unambiguously identifies the ordering. A square table cannot
    // be disambiguated from shape alone; that is a documented
    // limitation of the reader's auto-detection.
    const auto t = SyntheticEOSTable::makeAnalytic(12, 20);
    const std::string p_major = (tmp_dir_ / "major.h5").string();
    const std::string p_minor = (tmp_dir_ / "minor.h5").string();
    ASSERT_TRUE(writeRhoMajorEOSTable(p_major, t, "test"));
    ASSERT_TRUE(writeRhoMinorEOSTable(p_minor, t, "test"));

    TabulatedDataReader r1, r2;
    std::string err;
    ASSERT_TRUE(r1.load(p_major, err)) << err;
    ASSERT_TRUE(r2.load(p_minor, err)) << err;
    EXPECT_EQ(r1.axisOrdering(), AxisOrdering::RHO_MAJOR);
    EXPECT_EQ(r2.axisOrdering(), AxisOrdering::RHO_MINOR);

    // Both should evaluate to the same value.
    const double v1 = r1.evaluate(1.0e3, 1.0e6);
    const double v2 = r2.evaluate(1.0e3, 1.0e6);
    EXPECT_NEAR(v1, v2, 1.0e-9 * std::max(v1, 1.0));
#endif
}

TEST_F(TabulatedDataReaderTest, CitationMetadata)
{
#ifdef FSRM_TEST_HAVE_HDF5
    const auto t = SyntheticEOSTable::makeAnalytic(8, 8);
    const std::string path = (tmp_dir_ / "no_cite.h5").string();
    ASSERT_TRUE(writeRhoMajorEOSTable(path, t, "(unused)",
                                      "GRANITE",
                                      /*include_citation=*/false));

    TabulatedDataReader reader;
    std::string err;
    EXPECT_FALSE(reader.load(path, err));
    EXPECT_NE(err.find("source_citation"), std::string::npos)
        << "Expected error message to mention source_citation, got: "
        << err;
#endif
}

TEST_F(TabulatedDataReaderTest, LogspaceInterp)
{
#ifdef FSRM_TEST_HAVE_HDF5
    // With f(rho, e) = rho * e, log-space bilinear is exact.
    const auto t = SyntheticEOSTable::makeAnalytic(64, 64);
    const std::string path = (tmp_dir_ / "logspace.h5").string();
    ASSERT_TRUE(writeRhoMajorEOSTable(path, t, "test"));

    TabulatedDataReader reader;
    std::string err;
    ASSERT_TRUE(reader.load(path, err));

    // Sample 10 random-but-deterministic interior points.
    const double rho_lo = t.rho_axis.front();
    const double rho_hi = t.rho_axis.back();
    const double e_lo = t.e_axis.front();
    const double e_hi = t.e_axis.back();
    for (int k = 1; k <= 10; ++k) {
        const double fr = (k * 0.0931) - std::floor(k * 0.0931);
        const double fe = (k * 0.137) - std::floor(k * 0.137);
        const double rho = std::pow(10.0,
            std::log10(rho_lo) + fr * (std::log10(rho_hi) - std::log10(rho_lo)));
        const double e = std::pow(10.0,
            std::log10(e_lo) + fe * (std::log10(e_hi) - std::log10(e_lo)));
        const double v = reader.evaluate(rho, e);
        const double expected = rho * e;
        const double rel_err = std::abs(v - expected) / expected;
        EXPECT_LT(rel_err, 0.01) << "rho=" << rho << " e=" << e;
    }
#endif
}
