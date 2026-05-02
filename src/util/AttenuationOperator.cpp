/**
 * @file AttenuationOperator.cpp
 * @brief Frequency-dependent t*(f) attenuation operator (pass-3).
 *
 * Implements applyFrequencyDependentTStar for the DPRK 2017 synthetic
 * mb pipeline. The pass-2 test used a single scalar t* = 1.0 s, which
 * fidelity-document audit identified as already over-attenuating
 * compared to AK135 (~0.6-0.8 s for ~10-deg regional P). Pass-3
 * replaces the single-scalar with a frequency-dependent t*(f) =
 * t_star_ref * (f / f_ref)^(-alpha): low frequencies pass through
 * with less attenuation while the higher-frequency content -- which
 * dominates the band-passed mb measurement -- is more strongly
 * absorbed, recovering the physically motivated AK135 envelope.
 *
 * No third-party FFT dependency is added. The radix-2 Cooley-Tukey
 * transform below handles trace lengths up to 2^20 samples (~10^6),
 * well above the ~10^4 used by the DPRK test. Lengths that are not a
 * power of two are zero-padded to the next power of two before the
 * forward transform; the inverse transform truncates back to the
 * original length.
 *
 * References:
 *   Der, Z. A. and Lees, A. C. (1985), "Methodologies for estimating
 *     t*(f) from short-period body waves", BSSA 75(6).
 *   Choy, G. L. and Boatwright, J. L. (1995), "Global patterns of
 *     radiated seismic energy and apparent stress", JGR 100(B9).
 *   Kennett, B. L. N., Engdahl, E. R., Buland, R. (1995), "Constraints
 *     on seismic velocities in the Earth from traveltimes", GJI 122
 *     (AK135 reference Earth model).
 */

#include "util/AttenuationOperator.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>

namespace FSRM {
namespace util {

namespace {

using cplx = std::complex<double>;

// Round n up to the next power of two (n > 0).
inline std::size_t nextPow2(std::size_t n)
{
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

// In-place radix-2 Cooley-Tukey FFT. invert==true performs the inverse
// transform (with 1/N normalization). Length must be a power of two.
void radix2FFT(std::vector<cplx>& x, bool invert)
{
    const std::size_t n = x.size();
    if (n <= 1) return;

    // Bit-reversal permutation.
    std::size_t j = 0;
    for (std::size_t i = 1; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(x[i], x[j]);
    }

    // Cooley-Tukey butterflies.
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = (invert ? +2.0 : -2.0) * M_PI / static_cast<double>(len);
        const cplx w_n(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            cplx w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; ++k) {
                const cplx u = x[i + k];
                const cplx v = x[i + k + len / 2] * w;
                x[i + k]             = u + v;
                x[i + k + len / 2]   = u - v;
                w *= w_n;
            }
        }
    }

    if (invert) {
        const double inv_n = 1.0 / static_cast<double>(n);
        for (auto& z : x) z *= inv_n;
    }
}

}  // namespace

void applyFrequencyDependentTStar(std::vector<double>& trace,
                                  double dt,
                                  double t_star_ref,
                                  double f_ref,
                                  double alpha)
{
    if (trace.empty() || dt <= 0.0) return;

    // Zero-pad to the next power of two for radix-2 FFT.
    const std::size_t orig_len = trace.size();
    const std::size_t N = nextPow2(orig_len);

    std::vector<cplx> spectrum(N, cplx(0.0, 0.0));
    for (std::size_t i = 0; i < orig_len; ++i) {
        spectrum[i] = cplx(trace[i], 0.0);
    }

    radix2FFT(spectrum, /*invert=*/false);

    // Frequency vector matches NumPy's rfft layout for a real signal:
    // bin k corresponds to frequency k / (N * dt) for k in [0, N/2],
    // and to (k - N) / (N * dt) for k in (N/2, N-1] (negative
    // mirror). Apply the symmetric attenuation envelope in both
    // halves so the inverse FFT returns a real-valued trace.
    const double fs = 1.0 / dt;
    const double df = fs / static_cast<double>(N);

    auto attenuationFactor = [&](double f_hz) -> double {
        const double f_eff = std::max(std::abs(f_hz), df);  // avoid f=0 div-by-zero
        const double t_star_f = t_star_ref *
            std::pow(f_eff / f_ref, -alpha);
        const double atten = std::exp(-M_PI * f_eff * t_star_f);
        return atten;
    };

    for (std::size_t k = 0; k < N; ++k) {
        const double f_hz = (k <= N / 2)
            ? static_cast<double>(k) * df
            : (static_cast<double>(k) - static_cast<double>(N)) * df;
        const double atten = attenuationFactor(f_hz);
        spectrum[k] *= atten;
    }

    if (std::getenv("FSRM_ATTENUATION_DEBUG")) {
        const double f_band = std::sqrt(0.5 * 5.0);
        std::fprintf(stderr,
            "[AttenuationOperator] N=%zu, dt=%.4e, t_star_ref=%.3f, "
            "f_ref=%.3f, alpha=%.3f, atten(%.3f Hz)=%.4e\n",
            N, dt, t_star_ref, f_ref, alpha, f_band,
            attenuationFactor(f_band));
    }

    radix2FFT(spectrum, /*invert=*/true);

    // Truncate back to the original length, taking the real part
    // (imaginary part is round-off only because both spectrum and
    // attenuation are conjugate-symmetric in N).
    for (std::size_t i = 0; i < orig_len; ++i) {
        trace[i] = spectrum[i].real();
    }
}

}  // namespace util
}  // namespace FSRM
