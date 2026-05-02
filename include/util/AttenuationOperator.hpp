#ifndef FSRM_UTIL_ATTENUATION_OPERATOR_HPP
#define FSRM_UTIL_ATTENUATION_OPERATOR_HPP

#include <vector>

namespace FSRM {
namespace util {

/**
 * @brief Apply a frequency-dependent t*(f) attenuation in place.
 *
 * Convolves the time-domain trace with the kernel
 *
 *   exp(-pi * f * t*(f))
 *
 * where
 *
 *   t*(f) = t_star_ref * (f / f_ref)^(-alpha)
 *
 * The trace is FFT'd, multiplied by the attenuation envelope, and
 * inverse-FFT'd back to the time domain. The hand-rolled radix-2
 * Cooley-Tukey transform (in src/util/AttenuationOperator.cpp) keeps
 * FSRM dependency-free; an FSRM_ATTENUATION_DEBUG environment
 * variable prints the FFT length and band-centre attenuation factor
 * to stderr for verification runs.
 *
 * Defaults are AK135-consistent for short-period regional P:
 *   t_star_ref = 0.7 s, f_ref = 1.0 Hz, alpha = 0.4
 * (Der and Lees, 1985, BSSA 75, "Methodologies for estimating t*(f)
 * from short-period body waves").
 *
 * @param trace        In/out time-domain samples (uniform spacing).
 * @param dt           Sample interval (s).
 * @param t_star_ref   t* at the reference frequency f_ref.
 * @param f_ref        Reference frequency (Hz). Default 1.0 Hz.
 * @param alpha        Frequency exponent. Default 0.4 from Der & Lees.
 */
void applyFrequencyDependentTStar(std::vector<double>& trace,
                                  double dt,
                                  double t_star_ref,
                                  double f_ref = 1.0,
                                  double alpha = 0.4);

} // namespace util
} // namespace FSRM

#endif // FSRM_UTIL_ATTENUATION_OPERATOR_HPP
