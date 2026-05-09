/**
 * @file SACReader.hpp
 * @brief Production SAC binary reader/writer for the waveform V&V
 *        infrastructure. Supports header keyword extraction beyond
 *        the test-helper minimum (delta / b / npts).
 *
 * SAC binary file format reference: Goldstein & Snoke (2005),
 * "SAC Availability for the IRIS Community", DMS Electronic
 * Newsletter VII(1).
 *
 * Header layout (632 bytes total):
 *   bytes [0,   280):  70 floats   (DELTA, B, E, EVLA, EVLO, ...)
 *   bytes [280, 420):  35 ints     (NVHDR, NPTS, NZYEAR, ...)
 *   bytes [420, 440):   5 logicals
 *   bytes [440, 632):  24 strings  (KSTNM, KCMPNM, KEVNM, KNETWK, ...)
 *
 * The reader produces a SACTrace whose float_headers, int_headers,
 * and string_headers maps are keyed by the canonical SAC keyword
 * names. The pass-8 V&V code consults DELTA, B, NPTS, KSTNM, KCMPNM,
 * KNETWK, KEVNM, EVLA, EVLO, EVDP, MAG, STLA, STLO when present.
 *
 * State convention. SI / SAC convention; sample data is float[npts].
 * The reader does not deconvolve instrument response. The Python
 * refresh tool (tools/waveform_vv/refresh.py) is responsible for
 * removing the instrument response when fetching from IRIS so the
 * cached SAC files are already in displacement (m) at every sample.
 */

#ifndef FSRM_IO_SAC_READER_HPP
#define FSRM_IO_SAC_READER_HPP

#include <map>
#include <string>
#include <vector>

namespace FSRM {
namespace io {

/// Representation of one SAC trace in memory.
struct SACTrace
{
    bool valid = false;          ///< False on any IO or format error.
    bool big_endian = false;     ///< True if file was big-endian.
    double delta = 0.0;          ///< Sample interval [s] (DELTA).
    double begin_time = 0.0;     ///< First-sample time [s] (B).
    int npts = 0;                ///< Number of samples (NPTS).
    std::vector<float> samples;  ///< Sample values (typically displacement).

    /// Header keyword maps. Canonical SAC names (uppercase, no
    /// trailing whitespace). String headers are trimmed of trailing
    /// spaces / nulls. Sentinel values (-12345.0, -12345, "-12345  ")
    /// are removed before insertion so callers can iterate the maps
    /// without sentinel filtering.
    std::map<std::string, double> float_headers;
    std::map<std::string, int> int_headers;
    std::map<std::string, std::string> string_headers;

    /// Convenience accessors with default if missing.
    double getFloat(const std::string& key, double def = 0.0) const;
    int getInt(const std::string& key, int def = 0) const;
    std::string getString(const std::string& key,
                          const std::string& def = "") const;

    /// Time of sample i in absolute seconds since the trace begin.
    double sampleTime(int i) const
    {
        return begin_time + static_cast<double>(i) * delta;
    }
};

/// Read a SAC binary file into an in-memory SACTrace. Returns a
/// trace whose `valid` flag is false on any IO or format error;
/// callers branch on the flag.
bool readSAC(const std::string& path, SACTrace& trace);

/// Resample a trace to a target sample interval by linear
/// interpolation. Used to align observed and synthetic before metric
/// computation. The output trace inherits all headers but updates
/// DELTA and NPTS to match the new sampling.
SACTrace resampleSAC(const SACTrace& trace, double target_delta);

/// Window a trace to [t_begin, t_end] in absolute (B-relative)
/// seconds. The output's begin_time is shifted to t_begin and NPTS
/// is updated.
SACTrace windowSAC(const SACTrace& trace, double t_begin, double t_end);

/// Apply a real-cosine taper to the leading and trailing
/// taper_fraction of the trace samples. Value 0.05 means each end
/// gets a 5 percent cosine taper. No-op for taper_fraction <= 0.
SACTrace taperSAC(const SACTrace& trace, double taper_fraction);

/// Demean and detrend (linear least-squares) the trace samples.
SACTrace demeanDetrendSAC(const SACTrace& trace);

}  // namespace io
}  // namespace FSRM

#endif  // FSRM_IO_SAC_READER_HPP
