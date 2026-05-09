/**
 * @file test_sac_reader.cpp
 * @brief Unit tests for the production SAC reader (include/io/SACReader.hpp).
 *
 * Each test writes a minimal SAC binary to a temporary path, reads it
 * back through the production reader, and asserts the parsed fields
 * match. The test-helper writer is intentionally inline (not re-using
 * the FEM seismometer SAC writer) so the unit test exercises only the
 * reader.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "io/SACReader.hpp"

namespace
{

constexpr float SAC_FLOAT_SENTINEL = -12345.0f;
constexpr int SAC_INT_SENTINEL = -12345;
constexpr size_t SAC_HEADER_BYTES = 632;

void writeMinimalSAC(const std::string& path,
                     double delta, double begin_time, int npts,
                     const std::vector<float>& samples,
                     const std::string& kstnm = "TEST",
                     const std::string& kcmpnm = "BHZ",
                     double evla = 31.0, double evlo = -89.0)
{
    std::vector<char> hdr(SAC_HEADER_BYTES, 0);

    // Initialize all floats to sentinel.
    for (int i = 0; i < 70; ++i) {
        std::memcpy(hdr.data() + i * 4, &SAC_FLOAT_SENTINEL, 4);
    }
    // Initialize all ints to sentinel.
    for (int i = 0; i < 40; ++i) {
        std::memcpy(hdr.data() + 280 + i * 4, &SAC_INT_SENTINEL, 4);
    }
    // String section: blank-padded to length 8 (16 for KEVNM); we
    // write blanks for the sentinel-string convention to be cleared.
    std::memset(hdr.data() + 440, ' ', 632 - 440);

    auto setFloat = [&](int idx, double v) {
        const float fv = static_cast<float>(v);
        std::memcpy(hdr.data() + idx * 4, &fv, 4);
    };
    auto setInt = [&](int idx, int v) {
        std::memcpy(hdr.data() + 280 + idx * 4, &v, 4);
    };
    auto setString = [&](int byte_offset, int length,
                         const std::string& s) {
        for (int i = 0; i < length; ++i) {
            hdr[440 + byte_offset + i] =
                (i < static_cast<int>(s.size())) ? s[i] : ' ';
        }
    };

    // Floats: DELTA = 0, B = 5, EVLA = 31, EVLO = 36.
    setFloat(0, delta);
    setFloat(5, begin_time);
    setFloat(35, evla);
    setFloat(36, evlo);
    // Ints: NVHDR = 6 (idx 6), NPTS (idx 9).
    setInt(6, 6);
    setInt(9, npts);
    // Strings: KSTNM at offset 0, KCMPNM at offset 160.
    setString(0, 8, kstnm);
    setString(160, 8, kcmpnm);

    std::ofstream f(path, std::ios::binary);
    ASSERT_TRUE(f.good());
    f.write(hdr.data(), SAC_HEADER_BYTES);
    if (!samples.empty()) {
        f.write(reinterpret_cast<const char*>(samples.data()),
                samples.size() * sizeof(float));
    }
}

std::string tmpPath(const std::string& tag)
{
    return std::string("/tmp/fsrm_sac_") + tag + "_" +
           std::to_string(::getpid()) + ".sac";
}

}  // namespace

TEST(SACReaderTest, ReadDeltaBeginNpts)
{
    const std::string path = tmpPath("delta_b_npts");
    std::vector<float> samples = {1.0f, 2.0f, 3.0f, -1.0f, 0.5f};
    writeMinimalSAC(path, 0.01, 7.5, static_cast<int>(samples.size()),
                    samples);

    FSRM::io::SACTrace trace;
    EXPECT_TRUE(FSRM::io::readSAC(path, trace));
    EXPECT_TRUE(trace.valid);
    EXPECT_NEAR(trace.delta, 0.01, 1e-7);
    EXPECT_NEAR(trace.begin_time, 7.5, 1e-5);
    EXPECT_EQ(trace.npts, 5);
    ASSERT_EQ(trace.samples.size(), 5u);
    EXPECT_FLOAT_EQ(trace.samples[0], 1.0f);
    EXPECT_FLOAT_EQ(trace.samples[3], -1.0f);
    std::remove(path.c_str());
}

TEST(SACReaderTest, ReadStationAndComponent)
{
    const std::string path = tmpPath("station");
    std::vector<float> samples(8, 0.0f);
    writeMinimalSAC(path, 0.05, 0.0, static_cast<int>(samples.size()),
                    samples, "RGSD", "BHZ", 31.135, -89.578);

    FSRM::io::SACTrace trace;
    EXPECT_TRUE(FSRM::io::readSAC(path, trace));
    EXPECT_EQ(trace.getString("KSTNM"), "RGSD");
    EXPECT_EQ(trace.getString("KCMPNM"), "BHZ");
    EXPECT_NEAR(trace.getFloat("EVLA"), 31.135, 1e-3);
    EXPECT_NEAR(trace.getFloat("EVLO"), -89.578, 1e-3);
    std::remove(path.c_str());
}

TEST(SACReaderTest, MissingFileReturnsInvalid)
{
    FSRM::io::SACTrace trace;
    EXPECT_FALSE(FSRM::io::readSAC("/no/such/path.sac", trace));
    EXPECT_FALSE(trace.valid);
}

TEST(SACReaderTest, ResamplePreservesEndpoints)
{
    const std::string path = tmpPath("resample");
    std::vector<float> samples(11, 0.0f);
    for (int i = 0; i < 11; ++i)
        samples[i] = static_cast<float>(i);
    writeMinimalSAC(path, 0.1, 0.0, 11, samples);

    FSRM::io::SACTrace trace;
    ASSERT_TRUE(FSRM::io::readSAC(path, trace));
    auto resampled = FSRM::io::resampleSAC(trace, 0.05);
    EXPECT_NEAR(resampled.delta, 0.05, 1e-7);
    EXPECT_GE(resampled.npts, 21);
    // Sample at t = 0 should match the input start; sample at t = 1.0
    // (idx 20) should equal the input value at t = 1.0 (idx 10) = 10.
    EXPECT_NEAR(resampled.samples[0], 0.0f, 1e-3);
    EXPECT_NEAR(resampled.samples[20], 10.0f, 1e-3);
    std::remove(path.c_str());
}

TEST(SACReaderTest, WindowExtractsRange)
{
    const std::string path = tmpPath("window");
    std::vector<float> samples(101, 0.0f);
    for (int i = 0; i < 101; ++i)
        samples[i] = static_cast<float>(i);
    writeMinimalSAC(path, 0.01, 5.0, 101, samples);

    FSRM::io::SACTrace trace;
    ASSERT_TRUE(FSRM::io::readSAC(path, trace));
    auto windowed = FSRM::io::windowSAC(trace, 5.5, 5.7);
    EXPECT_GE(windowed.npts, 19);
    EXPECT_LE(windowed.npts, 22);
    EXPECT_NEAR(windowed.begin_time, 5.5, 1e-2);
}

TEST(SACReaderTest, DemeanDetrendRemovesLinearTrend)
{
    const std::string path = tmpPath("detrend");
    std::vector<float> samples(101, 0.0f);
    for (int i = 0; i < 101; ++i)
        samples[i] = static_cast<float>(2.0 + 0.05 * i);
    writeMinimalSAC(path, 0.01, 0.0, 101, samples);

    FSRM::io::SACTrace trace;
    ASSERT_TRUE(FSRM::io::readSAC(path, trace));
    auto detrended = FSRM::io::demeanDetrendSAC(trace);
    double sum = 0.0, max_abs = 0.0;
    for (float v : detrended.samples) {
        sum += v;
        max_abs = std::max(max_abs, static_cast<double>(std::abs(v)));
    }
    EXPECT_NEAR(sum / detrended.samples.size(), 0.0, 1e-3);
    EXPECT_LT(max_abs, 1e-3)
        << "Linear trend should be near-perfectly removed";
}

TEST(SACReaderTest, TaperZerosEndpoints)
{
    const std::string path = tmpPath("taper");
    std::vector<float> samples(101, 1.0f);
    writeMinimalSAC(path, 0.01, 0.0, 101, samples);

    FSRM::io::SACTrace trace;
    ASSERT_TRUE(FSRM::io::readSAC(path, trace));
    auto tapered = FSRM::io::taperSAC(trace, 0.1);
    EXPECT_LT(tapered.samples.front(), 0.1f);
    EXPECT_LT(tapered.samples.back(), 0.1f);
    EXPECT_NEAR(tapered.samples[50], 1.0f, 1e-3);
}
