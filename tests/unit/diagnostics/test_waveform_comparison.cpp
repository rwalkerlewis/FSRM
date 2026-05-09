/**
 * @file test_waveform_comparison.cpp
 * @brief Unit tests for the waveform comparison metrics library
 *        (include/diagnostics/WaveformComparison.hpp).
 *
 * Each test constructs synthetic SAC traces in-memory (no file IO),
 * verifies the metric returns the expected value within a documented
 * tolerance.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "diagnostics/WaveformComparison.hpp"
#include "io/SACReader.hpp"

namespace
{

constexpr double PI = 3.14159265358979323846;

FSRM::io::SACTrace makeSineTrace(double freq_hz, double amplitude,
                                 double dt = 0.01, int npts = 1024)
{
    FSRM::io::SACTrace t;
    t.valid = true;
    t.delta = dt;
    t.begin_time = 0.0;
    t.npts = npts;
    t.samples.assign(npts, 0.0f);
    for (int i = 0; i < npts; ++i) {
        const double tt = i * dt;
        t.samples[i] = static_cast<float>(amplitude *
                                          std::sin(2.0 * PI * freq_hz * tt));
    }
    return t;
}

FSRM::io::SACTrace shiftedTrace(const FSRM::io::SACTrace& src, int shift_samples)
{
    FSRM::io::SACTrace out = src;
    out.samples.assign(src.npts, 0.0f);
    for (int i = 0; i < src.npts; ++i) {
        const int j = i - shift_samples;
        if (j >= 0 && j < src.npts) out.samples[i] = src.samples[j];
    }
    return out;
}

}  // namespace

TEST(WaveformComparisonTest, PeakAmplitudeRatioOnSynthetic)
{
    auto a = makeSineTrace(2.0, 1.0);
    auto b = makeSineTrace(2.0, 3.0);
    const double ratio = FSRM::diagnostics::peakAmplitudeRatio(a, b);
    EXPECT_NEAR(ratio, 3.0, 0.05);
}

TEST(WaveformComparisonTest, SpectralAmplitudeRatioOnSynthetic)
{
    auto a = makeSineTrace(2.0, 1.0);
    auto b = makeSineTrace(2.0, 5.0);
    const double ratio =
        FSRM::diagnostics::spectralAmplitudeRatio(a, b, 0.5, 5.0);
    EXPECT_NEAR(ratio, 5.0, 0.5);
}

TEST(WaveformComparisonTest, DominantFrequencyOnSynthetic)
{
    auto t = makeSineTrace(3.0, 1.0, 0.005, 4096);
    const double f = FSRM::diagnostics::dominantFrequency(t);
    EXPECT_NEAR(f, 3.0, 0.5);
}

TEST(WaveformComparisonTest, EnvelopeMisfitZeroForSelfCompare)
{
    auto t = makeSineTrace(2.0, 1.0);
    const double m = FSRM::diagnostics::envelopeMisfit(t, t);
    EXPECT_LT(m, 1e-3);
}

TEST(WaveformComparisonTest, EnvelopeMisfitLargeForUnrelated)
{
    auto a = makeSineTrace(2.0, 1.0);
    auto b = makeSineTrace(7.0, 1.0);
    const double m = FSRM::diagnostics::envelopeMisfit(a, b);
    EXPECT_GT(m, 0.05);
}

TEST(WaveformComparisonTest, CrossCorrelationZeroLagForSelf)
{
    auto t = makeSineTrace(2.0, 1.0);
    double lag = -1.0;
    const double cc = FSRM::diagnostics::crossCorrelation(t, t, 0.5, lag);
    EXPECT_NEAR(cc, 1.0, 0.05);
    EXPECT_NEAR(lag, 0.0, 0.05);
}

TEST(WaveformComparisonTest, CrossCorrelationDetectsKnownShift)
{
    auto a = makeSineTrace(2.0, 1.0, 0.01, 1024);
    auto b = shiftedTrace(a, 10);  // 10 samples = 0.1 s shift
    double lag = -1.0;
    const double cc = FSRM::diagnostics::crossCorrelation(a, b, 0.5, lag);
    EXPECT_GT(cc, 0.5);
    // CC peak should sit near +0.1 s lag.
    EXPECT_NEAR(lag, 0.1, 0.05);
}

TEST(WaveformComparisonTest, FFTRoundTrip)
{
    std::vector<double> re = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> im(8, 0.0);
    auto re0 = re;
    FSRM::diagnostics::fftForward(re, im);
    FSRM::diagnostics::fftInverse(re, im);
    for (size_t i = 0; i < re.size(); ++i) {
        EXPECT_NEAR(re[i], re0[i], 1e-9);
        EXPECT_NEAR(im[i], 0.0, 1e-9);
    }
}
