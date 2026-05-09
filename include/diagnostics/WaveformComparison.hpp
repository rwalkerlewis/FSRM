/**
 * @file WaveformComparison.hpp
 * @brief Comparison metrics for the IRIS waveform V&V infrastructure.
 *
 * Pass-8 (axis 1, see docs/HISTORIC_NUCLEAR_ROADMAP.md) introduces a
 * minimal C++ comparison library for matching simulated synthetic
 * seismograms against IRIS-fetched observed traces. The five metrics
 * here gate the new iris_validation CTest label tests; each is
 * individually unit-tested.
 *
 *   peakAmplitudeRatio: |peak(syn)| / |peak(obs)| in a chosen window.
 *   spectralAmplitudeRatio: integrated |F(f)| in [fmin,fmax] ratio.
 *   dominantFrequency: argmax of the smoothed spectrum.
 *   envelopeMisfit: L2 of (env(syn) - env(obs)) / ||env(obs)||.
 *   crossCorrelation: max normalized CC over a lag window.
 *
 * The implementation uses a hand-rolled radix-2 Cooley-Tukey FFT
 * (no external FFT dep). Traces are resampled onto a common dt and
 * zero-padded to the next power-of-two before any spectral metric.
 *
 * State convention. SI throughout. Frequencies in Hz, lags in
 * seconds.
 *
 * References.
 *  - Goldstein, P. and Snoke, A. (2005), "SAC Availability for the
 *    IRIS Community", DMS Electronic Newsletter VII(1).
 *  - Cooley, J. W. and Tukey, J. W. (1965), "An algorithm for the
 *    machine calculation of complex Fourier series", Math. Comp 19,
 *    pp. 297-301.
 *  - Hilbert transform for the analytic signal envelope: Bracewell
 *    (1986), "The Fourier Transform and Its Applications", McGraw-
 *    Hill, ch 11.
 */

#ifndef FSRM_DIAGNOSTICS_WAVEFORM_COMPARISON_HPP
#define FSRM_DIAGNOSTICS_WAVEFORM_COMPARISON_HPP

#include "io/SACReader.hpp"

#include <vector>

namespace FSRM {
namespace diagnostics {

/// Aggregate result returned by compareWaveforms.
struct WaveformComparisonResult
{
    bool valid = false;
    double peak_amplitude_ratio = 0.0;
    double spectral_amplitude_ratio = 0.0;
    double dominant_freq_obs_hz = 0.0;
    double dominant_freq_syn_hz = 0.0;
    double envelope_misfit_l2 = 0.0;
    double cross_correlation = 0.0;
    double cross_correlation_lag_s = 0.0;
};

/// |peak(syn)| / |peak(obs)|.
double peakAmplitudeRatio(const io::SACTrace& obs,
                          const io::SACTrace& syn);

/// Mean integrated |F(f)| over [fmin, fmax]: ratio syn/obs.
double spectralAmplitudeRatio(const io::SACTrace& obs,
                              const io::SACTrace& syn,
                              double fmin_hz, double fmax_hz);

/// Frequency at the maximum of the smoothed amplitude spectrum.
/// Spectrum is smoothed by a 5-point boxcar to suppress numerical
/// noise.
double dominantFrequency(const io::SACTrace& trace);

/// L2 norm of the difference between Hilbert-transform envelopes,
/// normalized by ||env(obs)||.
double envelopeMisfit(const io::SACTrace& obs,
                      const io::SACTrace& syn);

/// Normalized cross-correlation, maximized over lags in
/// [-max_lag_s, +max_lag_s]. out_lag_s receives the lag of the
/// maximum.
double crossCorrelation(const io::SACTrace& obs,
                        const io::SACTrace& syn,
                        double max_lag_s, double& out_lag_s);

/// Aggregate convenience: runs all five metrics, fills the struct.
WaveformComparisonResult compareWaveforms(const io::SACTrace& obs,
                                          const io::SACTrace& syn,
                                          double fmin_hz, double fmax_hz,
                                          double max_lag_s);

/// Hand-rolled radix-2 Cooley-Tukey FFT. Padded to the next
/// power-of-two; the imaginary part of the output is zero on input
/// for real-valued traces. Exposed for tests.
void fftForward(std::vector<double>& re, std::vector<double>& im);
void fftInverse(std::vector<double>& re, std::vector<double>& im);

}  // namespace diagnostics
}  // namespace FSRM

#endif  // FSRM_DIAGNOSTICS_WAVEFORM_COMPARISON_HPP
