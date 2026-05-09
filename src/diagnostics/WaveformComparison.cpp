/**
 * @file WaveformComparison.cpp
 * @brief Implementation of the waveform comparison metrics. See
 *        WaveformComparison.hpp for the design rationale.
 */

#include "diagnostics/WaveformComparison.hpp"

#include <algorithm>
#include <cmath>

namespace FSRM {
namespace diagnostics {

namespace
{

constexpr double PI = 3.14159265358979323846;

int nextPowerOfTwo(int n)
{
    int p = 1;
    while (p < n) p *= 2;
    return p;
}

void fft_impl(std::vector<double>& re, std::vector<double>& im, bool inverse)
{
    const int n = static_cast<int>(re.size());
    if (n < 2) return;
    // Bit-reversal permutation.
    int j = 0;
    for (int i = 0; i < n - 1; ++i) {
        if (i < j) { std::swap(re[i], re[j]); std::swap(im[i], im[j]); }
        int m = n >> 1;
        while (j >= m && m > 0) { j -= m; m >>= 1; }
        j += m;
    }
    // Cooley-Tukey butterflies.
    for (int s = 1; s <= 30; ++s) {
        const int m = 1 << s;
        if (m > n) break;
        const int m2 = m >> 1;
        const double theta = (inverse ? 2.0 : -2.0) * PI / m;
        const double wpr = std::cos(theta);
        const double wpi = std::sin(theta);
        for (int k = 0; k < n; k += m) {
            double wr = 1.0, wi = 0.0;
            for (int p = 0; p < m2; ++p) {
                const int t = k + p + m2;
                const int u = k + p;
                const double tr = wr * re[t] - wi * im[t];
                const double ti = wr * im[t] + wi * re[t];
                re[t] = re[u] - tr;
                im[t] = im[u] - ti;
                re[u] += tr;
                im[u] += ti;
                const double wr_new = wr * wpr - wi * wpi;
                wi = wr * wpi + wi * wpr;
                wr = wr_new;
            }
        }
    }
    if (inverse) {
        const double inv = 1.0 / n;
        for (int i = 0; i < n; ++i) { re[i] *= inv; im[i] *= inv; }
    }
}

/// Pad a real-valued time series to the next power-of-two. Output
/// arrays have re = padded, im = 0.
void padToPowerOfTwo(const std::vector<float>& src,
                     std::vector<double>& re, std::vector<double>& im)
{
    const int n_in = static_cast<int>(src.size());
    const int n_pad = nextPowerOfTwo(std::max(n_in, 4));
    re.assign(n_pad, 0.0);
    im.assign(n_pad, 0.0);
    for (int i = 0; i < n_in; ++i) re[i] = static_cast<double>(src[i]);
}

/// Hilbert-transform envelope of a real time series.
std::vector<double> hilbertEnvelope(const std::vector<float>& samples)
{
    std::vector<double> re, im;
    padToPowerOfTwo(samples, re, im);
    const int n = static_cast<int>(re.size());
    fft_impl(re, im, /*inverse=*/false);
    // Construct analytic signal: zero negative frequencies, double
    // positive.
    for (int i = 1; i < n / 2; ++i) {
        re[i] *= 2.0;
        im[i] *= 2.0;
    }
    for (int i = n / 2 + 1; i < n; ++i) {
        re[i] = 0.0;
        im[i] = 0.0;
    }
    fft_impl(re, im, /*inverse=*/true);
    std::vector<double> env(samples.size(), 0.0);
    for (size_t i = 0; i < samples.size(); ++i) {
        env[i] = std::sqrt(re[i] * re[i] + im[i] * im[i]);
    }
    return env;
}

/// Resample / reconcile two traces to a common dt and length. Used
/// by every spectral metric. Returns the common dt; resamples obs and
/// syn into the supplied output vectors.
double conformTraces(const io::SACTrace& obs, const io::SACTrace& syn,
                     std::vector<float>& obs_out, std::vector<float>& syn_out)
{
    if (!obs.valid || !syn.valid) return 0.0;
    const double dt_common = std::min(obs.delta, syn.delta);
    const auto obs_r =
        (std::abs(obs.delta - dt_common) < 1e-12) ? obs
                                                  : io::resampleSAC(obs, dt_common);
    const auto syn_r =
        (std::abs(syn.delta - dt_common) < 1e-12) ? syn
                                                  : io::resampleSAC(syn, dt_common);
    const int n = std::min(obs_r.npts, syn_r.npts);
    obs_out.assign(obs_r.samples.begin(), obs_r.samples.begin() + n);
    syn_out.assign(syn_r.samples.begin(), syn_r.samples.begin() + n);
    return dt_common;
}

double maxAbs(const std::vector<float>& v)
{
    double m = 0.0;
    for (float x : v)
        m = std::max(m, static_cast<double>(std::abs(x)));
    return m;
}

}  // namespace

void fftForward(std::vector<double>& re, std::vector<double>& im)
{
    fft_impl(re, im, /*inverse=*/false);
}

void fftInverse(std::vector<double>& re, std::vector<double>& im)
{
    fft_impl(re, im, /*inverse=*/true);
}

double peakAmplitudeRatio(const io::SACTrace& obs,
                          const io::SACTrace& syn)
{
    if (!obs.valid || !syn.valid) return 0.0;
    const double obs_p = maxAbs(obs.samples);
    const double syn_p = maxAbs(syn.samples);
    if (obs_p <= 0.0) return 0.0;
    return syn_p / obs_p;
}

double spectralAmplitudeRatio(const io::SACTrace& obs,
                              const io::SACTrace& syn,
                              double fmin_hz, double fmax_hz)
{
    std::vector<float> obs_s, syn_s;
    const double dt = conformTraces(obs, syn, obs_s, syn_s);
    if (dt <= 0.0 || obs_s.empty()) return 0.0;
    std::vector<double> re_o, im_o, re_s, im_s;
    padToPowerOfTwo(obs_s, re_o, im_o);
    padToPowerOfTwo(syn_s, re_s, im_s);
    fft_impl(re_o, im_o, false);
    fft_impl(re_s, im_s, false);
    const int n = static_cast<int>(re_o.size());
    const double df = 1.0 / (n * dt);
    double sum_o = 0.0, sum_s = 0.0;
    for (int i = 0; i < n / 2; ++i) {
        const double f = i * df;
        if (f < fmin_hz || f > fmax_hz) continue;
        sum_o += std::sqrt(re_o[i] * re_o[i] + im_o[i] * im_o[i]);
        sum_s += std::sqrt(re_s[i] * re_s[i] + im_s[i] * im_s[i]);
    }
    if (sum_o <= 0.0) return 0.0;
    return sum_s / sum_o;
}

double dominantFrequency(const io::SACTrace& trace)
{
    if (!trace.valid || trace.npts < 4) return 0.0;
    std::vector<double> re, im;
    padToPowerOfTwo(trace.samples, re, im);
    fft_impl(re, im, false);
    const int n = static_cast<int>(re.size());
    const double df = 1.0 / (n * trace.delta);
    std::vector<double> spec(n / 2, 0.0);
    for (int i = 0; i < n / 2; ++i) {
        spec[i] = std::sqrt(re[i] * re[i] + im[i] * im[i]);
    }
    // 5-point boxcar smoothing.
    std::vector<double> sm = spec;
    const int half = 2;
    for (int i = half; i < static_cast<int>(spec.size()) - half; ++i) {
        double s = 0.0;
        for (int k = -half; k <= half; ++k) s += spec[i + k];
        sm[i] = s / (2 * half + 1);
    }
    int idx = 1;  // skip DC
    double best = sm[1];
    for (int i = 2; i < static_cast<int>(sm.size()); ++i) {
        if (sm[i] > best) { best = sm[i]; idx = i; }
    }
    return idx * df;
}

double envelopeMisfit(const io::SACTrace& obs,
                      const io::SACTrace& syn)
{
    std::vector<float> obs_s, syn_s;
    const double dt = conformTraces(obs, syn, obs_s, syn_s);
    if (dt <= 0.0 || obs_s.empty()) return 0.0;
    auto env_o = hilbertEnvelope(obs_s);
    auto env_s = hilbertEnvelope(syn_s);
    double num = 0.0, den = 0.0;
    const size_t n = std::min(env_o.size(), env_s.size());
    for (size_t i = 0; i < n; ++i) {
        const double d = env_s[i] - env_o[i];
        num += d * d;
        den += env_o[i] * env_o[i];
    }
    if (den <= 0.0) return 0.0;
    return std::sqrt(num / den);
}

double crossCorrelation(const io::SACTrace& obs,
                        const io::SACTrace& syn,
                        double max_lag_s, double& out_lag_s)
{
    out_lag_s = 0.0;
    std::vector<float> obs_s, syn_s;
    const double dt = conformTraces(obs, syn, obs_s, syn_s);
    if (dt <= 0.0 || obs_s.empty()) return 0.0;
    const int n = static_cast<int>(obs_s.size());
    const int max_lag = std::min(n - 1,
        static_cast<int>(std::ceil(max_lag_s / dt)));
    // Demean.
    auto demean = [](std::vector<float>& v) {
        double m = 0.0;
        for (float x : v) m += x;
        m /= v.size();
        for (float& x : v) x -= static_cast<float>(m);
    };
    demean(obs_s);
    demean(syn_s);
    auto norm2 = [](const std::vector<float>& v) {
        double s = 0.0;
        for (float x : v) s += static_cast<double>(x) * x;
        return std::sqrt(s);
    };
    const double n_o = norm2(obs_s);
    const double n_s = norm2(syn_s);
    if (n_o <= 0.0 || n_s <= 0.0) return 0.0;
    double best = -1.0e30;
    int best_lag = 0;
    for (int lag = -max_lag; lag <= max_lag; ++lag) {
        double sum = 0.0;
        const int i_lo = std::max(0, -lag);
        const int i_hi = std::min(n, n - lag);
        for (int i = i_lo; i < i_hi; ++i) {
            sum += static_cast<double>(obs_s[i]) *
                   static_cast<double>(syn_s[i + lag]);
        }
        const double cc = sum / (n_o * n_s);
        if (cc > best) { best = cc; best_lag = lag; }
    }
    out_lag_s = best_lag * dt;
    return best;
}

WaveformComparisonResult compareWaveforms(const io::SACTrace& obs,
                                          const io::SACTrace& syn,
                                          double fmin_hz, double fmax_hz,
                                          double max_lag_s)
{
    WaveformComparisonResult r;
    if (!obs.valid || !syn.valid) return r;
    r.peak_amplitude_ratio = peakAmplitudeRatio(obs, syn);
    r.spectral_amplitude_ratio =
        spectralAmplitudeRatio(obs, syn, fmin_hz, fmax_hz);
    r.dominant_freq_obs_hz = dominantFrequency(obs);
    r.dominant_freq_syn_hz = dominantFrequency(syn);
    r.envelope_misfit_l2 = envelopeMisfit(obs, syn);
    r.cross_correlation = crossCorrelation(obs, syn, max_lag_s,
                                            r.cross_correlation_lag_s);
    r.valid = true;
    return r;
}

}  // namespace diagnostics
}  // namespace FSRM
