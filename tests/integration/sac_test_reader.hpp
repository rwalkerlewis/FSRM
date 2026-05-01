/**
 * @file sac_test_reader.hpp
 * @brief Minimal SAC binary reader for integration tests.
 *
 * Reads enough of the SAC (Seismic Analysis Code) binary file format to
 * extract the per-sample displacement trace, the start time (B), and the
 * sample interval (DELTA). Supports both little-endian and big-endian
 * SAC files via the standard NVHDR == 6 sentinel.
 *
 * The full SAC header is 632 bytes:
 *   70 floats (280 bytes) -> bytes [0, 280)
 *   35 ints   (140 bytes) -> bytes [280, 420)
 *    5 logicals(20 bytes) -> bytes [420, 440)
 *   24 strings(192 bytes) -> bytes [440, 632)
 *
 * Field positions (Goldstein, P. & Snoke, A.,
 * "SAC Availability for the IRIS Community", DMS Electronic Newsletter
 * Volume VII Number 1, March 2005):
 *   DELTA = floats[0]      (sample interval, seconds)
 *   B     = floats[5]      (begin time, seconds)
 *   NVHDR = ints[6]        (always 6)
 *   NPTS  = ints[9]        (number of points)
 *
 * This reader is intentionally header-light -- it does not interpret
 * station coordinates or pick times. Tests that need those should use
 * the simulation config they wrote, since the SAC writer copies them
 * verbatim from the seismometer record.
 */

#ifndef FSRM_SAC_TEST_READER_HPP
#define FSRM_SAC_TEST_READER_HPP

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace FSRM {
namespace test_helpers {

struct SACTrace
{
  bool valid = false;          // True when the trace was loaded successfully
  bool big_endian = false;     // True if file was big-endian
  double delta = 0.0;          // Sample interval (s)
  double begin_time = 0.0;     // First sample time relative to event (s)
  int npts = 0;                // Number of samples
  std::vector<float> samples;  // Sample values (typically displacement)
};

namespace detail {

inline uint32_t swapBytes(uint32_t v)
{
  return ((v & 0x000000FFu) << 24) | ((v & 0x0000FF00u) << 8)
       | ((v & 0x00FF0000u) >> 8)  | ((v & 0xFF000000u) >> 24);
}

inline float readFloat(const char* buf, bool swap)
{
  uint32_t raw;
  std::memcpy(&raw, buf, 4);
  if (swap) raw = swapBytes(raw);
  float f;
  std::memcpy(&f, &raw, 4);
  return f;
}

inline int32_t readInt(const char* buf, bool swap)
{
  uint32_t raw;
  std::memcpy(&raw, buf, 4);
  if (swap) raw = swapBytes(raw);
  int32_t v;
  std::memcpy(&v, &raw, 4);
  return v;
}

}  // namespace detail

/**
 * @brief Read a SAC binary file into an in-memory SACTrace.
 *
 * Returns a SACTrace whose `valid` flag is false on any IO or format
 * error. The caller can branch on `trace.valid` and use the data fields
 * only when valid.
 */
inline SACTrace readSAC(const std::string& path)
{
  SACTrace trace;
  std::ifstream f(path, std::ios::binary);
  if (!f.good()) return trace;

  // Read the whole 632-byte header into a buffer.
  constexpr size_t kHeaderBytes = 632;
  std::vector<char> header(kHeaderBytes);
  f.read(header.data(), kHeaderBytes);
  if (f.gcount() != static_cast<std::streamsize>(kHeaderBytes)) {
    return trace;
  }

  // NVHDR is the 7th int in the int section. Int section starts at byte
  // 280; NVHDR is at offset 280 + 6*4 = 304. Native value should be 6;
  // a swapped value would also be 6 only if endianness matches, so we
  // compare against 6 in both interpretations.
  const char* nvhdr_pos = header.data() + 304;
  int32_t nvhdr_native = detail::readInt(nvhdr_pos, false);
  int32_t nvhdr_swap   = detail::readInt(nvhdr_pos, true);

  bool swap = false;
  if (nvhdr_native == 6) {
    swap = false;
  } else if (nvhdr_swap == 6) {
    swap = true;
  } else {
    return trace;  // Not a recognisable SAC file.
  }

  // Parse the fields we care about.
  trace.delta      = detail::readFloat(header.data() + 0  * 4, swap);
  trace.begin_time = detail::readFloat(header.data() + 5  * 4, swap);
  trace.npts       = detail::readInt(  header.data() + 280 + 9 * 4, swap);

  if (trace.npts <= 0 || trace.delta <= 0.0) {
    return trace;
  }

  trace.samples.resize(static_cast<size_t>(trace.npts));
  std::vector<char> data_buf(static_cast<size_t>(trace.npts) * 4);
  f.read(data_buf.data(), data_buf.size());
  if (f.gcount() != static_cast<std::streamsize>(data_buf.size())) {
    trace.samples.clear();
    return trace;
  }
  for (int i = 0; i < trace.npts; ++i) {
    trace.samples[i] =
        detail::readFloat(data_buf.data() + i * 4, swap);
  }

  trace.big_endian = swap;
  trace.valid = true;
  return trace;
}

/**
 * @brief Peak absolute amplitude of the trace and the time at which it
 * occurs (in seconds since the event origin).
 */
struct PeakResult
{
  double peak_abs = 0.0;
  double peak_time = 0.0;
  double peak_signed = 0.0;
};

inline PeakResult tracePeak(const SACTrace& trace)
{
  PeakResult r;
  if (!trace.valid || trace.samples.empty()) return r;
  for (int i = 0; i < trace.npts; ++i) {
    double v = static_cast<double>(trace.samples[i]);
    if (std::abs(v) > r.peak_abs) {
      r.peak_abs = std::abs(v);
      r.peak_signed = v;
      r.peak_time = trace.begin_time +
                    static_cast<double>(i) * trace.delta;
    }
  }
  return r;
}

/**
 * @brief Sign of the first non-trivial motion in the trace.
 *
 * Returns +1 / -1 / 0. A sample is considered "non-trivial" when its
 * absolute value exceeds `threshold * peak_abs`. Default threshold is
 * 5% of the peak so we ignore numerical noise in the leading window.
 */
inline int firstMotionSign(const SACTrace& trace, double threshold = 0.05)
{
  if (!trace.valid || trace.samples.empty()) return 0;
  PeakResult peak = tracePeak(trace);
  const double cutoff = peak.peak_abs * threshold;
  for (int i = 0; i < trace.npts; ++i) {
    double v = static_cast<double>(trace.samples[i]);
    if (std::abs(v) > cutoff) {
      return (v > 0.0) ? +1 : -1;
    }
  }
  return 0;
}

/**
 * @brief L2 (Euclidean) norm of the trace samples.
 */
inline double traceL2Norm(const SACTrace& trace)
{
  if (!trace.valid) return 0.0;
  double s2 = 0.0;
  for (float v : trace.samples) {
    s2 += static_cast<double>(v) * static_cast<double>(v);
  }
  return std::sqrt(s2);
}

}  // namespace test_helpers
}  // namespace FSRM

#endif  // FSRM_SAC_TEST_READER_HPP
