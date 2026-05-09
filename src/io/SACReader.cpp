/**
 * @file SACReader.cpp
 * @brief Implementation of the production SAC binary reader.
 */

#include "io/SACReader.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <vector>

namespace FSRM {
namespace io {

namespace
{

constexpr size_t kHeaderBytes = 632;
constexpr int kNumFloats = 70;
constexpr int kNumInts = 40;        // 35 ints + 5 logicals (treated as ints)
constexpr int kStringSectionStart = 440;
constexpr float kFloatSentinel = -12345.0f;
constexpr int kIntSentinel = -12345;

uint32_t swapBytes(uint32_t v)
{
    return ((v & 0x000000FFu) << 24) | ((v & 0x0000FF00u) << 8)
         | ((v & 0x00FF0000u) >> 8)  | ((v & 0xFF000000u) >> 24);
}

float readFloat(const char* buf, bool swap)
{
    uint32_t raw;
    std::memcpy(&raw, buf, 4);
    if (swap) raw = swapBytes(raw);
    float f;
    std::memcpy(&f, &raw, 4);
    return f;
}

int32_t readInt(const char* buf, bool swap)
{
    uint32_t raw;
    std::memcpy(&raw, buf, 4);
    if (swap) raw = swapBytes(raw);
    int32_t v;
    std::memcpy(&v, &raw, 4);
    return v;
}

std::string trimTrailing(const std::string& s)
{
    size_t end = s.size();
    while (end > 0 && (s[end - 1] == ' ' || s[end - 1] == '\0'))
        --end;
    return s.substr(0, end);
}

bool isSentinelString(const std::string& s)
{
    // Trimmed sentinel string is "-12345" (the SAC convention).
    const std::string t = trimTrailing(s);
    return t == "-12345" || t.empty();
}

// Canonical SAC header keyword name lists, indexed by their position
// in the float / int / string header sections. From Goldstein & Snoke
// (2005) Table 1.
const char* kFloatNames[kNumFloats] = {
    "DELTA",   "DEPMIN",  "DEPMAX",  "SCALE",   "ODELTA",
    "B",       "E",       "O",       "A",       "INTERNAL1",
    "T0",      "T1",      "T2",      "T3",      "T4",
    "T5",      "T6",      "T7",      "T8",      "T9",
    "F",       "RESP0",   "RESP1",   "RESP2",   "RESP3",
    "RESP4",   "RESP5",   "RESP6",   "RESP7",   "RESP8",
    "RESP9",   "STLA",    "STLO",    "STEL",    "STDP",
    "EVLA",    "EVLO",    "EVEL",    "EVDP",    "MAG",
    "USER0",   "USER1",   "USER2",   "USER3",   "USER4",
    "USER5",   "USER6",   "USER7",   "USER8",   "USER9",
    "DIST",    "AZ",      "BAZ",     "GCARC",   "INTERNAL2",
    "INTERNAL3","DEPMEN", "CMPAZ",   "CMPINC",  "XMINIMUM",
    "XMAXIMUM","YMINIMUM","YMAXIMUM","UNUSED1", "UNUSED2",
    "UNUSED3", "UNUSED4", "UNUSED5", "UNUSED6", "UNUSED7"};

const char* kIntNames[kNumInts] = {
    "NZYEAR", "NZJDAY",  "NZHOUR",  "NZMIN",  "NZSEC",
    "NZMSEC", "NVHDR",   "NORID",   "NEVID",  "NPTS",
    "INTERNAL4","NWFID", "NXSIZE",  "NYSIZE", "UNUSED8",
    "IFTYPE", "IDEP",    "IZTYPE",  "UNUSED9","IINST",
    "ISTREG", "IEVREG",  "IEVTYP",  "IQUAL",  "ISYNTH",
    "IMAGTYP","IMAGSRC", "UNUSED10","UNUSED11","UNUSED12",
    "UNUSED13","UNUSED14","UNUSED15","UNUSED16","UNUSED17",
    "LEVEN",  "LPSPOL",  "LOVROK",  "LCALDA", "UNUSED18"};

// String headers: 24 entries; KEVNM is double-length (16 chars).
struct StringSpec { const char* name; int byte_offset; int length; };
const StringSpec kStringSpecs[] = {
    {"KSTNM",  0,    8}, {"KEVNM",  8,   16}, {"KHOLE",  24,   8},
    {"KO",     32,   8}, {"KA",     40,   8}, {"KT0",    48,   8},
    {"KT1",    56,   8}, {"KT2",    64,   8}, {"KT3",    72,   8},
    {"KT4",    80,   8}, {"KT5",    88,   8}, {"KT6",    96,   8},
    {"KT7",   104,   8}, {"KT8",   112,   8}, {"KT9",   120,   8},
    {"KF",    128,   8}, {"KUSER0",136,   8}, {"KUSER1",144,   8},
    {"KUSER2",152,   8}, {"KCMPNM",160,   8}, {"KNETWK",168,   8},
    {"KDATRD",176,   8}, {"KINST", 184,   8}};
constexpr int kNumStrings = sizeof(kStringSpecs) / sizeof(kStringSpecs[0]);

}  // namespace

double SACTrace::getFloat(const std::string& key, double def) const
{
    auto it = float_headers.find(key);
    return it == float_headers.end() ? def : it->second;
}

int SACTrace::getInt(const std::string& key, int def) const
{
    auto it = int_headers.find(key);
    return it == int_headers.end() ? def : it->second;
}

std::string SACTrace::getString(const std::string& key,
                                const std::string& def) const
{
    auto it = string_headers.find(key);
    return it == string_headers.end() ? def : it->second;
}

bool readSAC(const std::string& path, SACTrace& trace)
{
    trace = SACTrace{};
    std::ifstream f(path, std::ios::binary);
    if (!f.good()) return false;

    std::vector<char> hdr(kHeaderBytes);
    f.read(hdr.data(), kHeaderBytes);
    if (f.gcount() != static_cast<std::streamsize>(kHeaderBytes))
        return false;

    // NVHDR sentinel test for endianness.
    const char* nvhdr_pos = hdr.data() + 280 + 6 * 4;
    int32_t nvhdr_native = readInt(nvhdr_pos, false);
    int32_t nvhdr_swap = readInt(nvhdr_pos, true);

    bool swap = false;
    if (nvhdr_native == 6) {
        swap = false;
    } else if (nvhdr_swap == 6) {
        swap = true;
    } else {
        return false;
    }

    // Read all 70 floats.
    for (int i = 0; i < kNumFloats; ++i) {
        const float v = readFloat(hdr.data() + i * 4, swap);
        if (v == kFloatSentinel) continue;
        trace.float_headers[kFloatNames[i]] = static_cast<double>(v);
    }
    // Read all 40 ints (35 + 5 logicals).
    for (int i = 0; i < kNumInts; ++i) {
        const int32_t v = readInt(hdr.data() + 280 + i * 4, swap);
        if (v == kIntSentinel) continue;
        trace.int_headers[kIntNames[i]] = static_cast<int>(v);
    }
    // Read string section.
    const char* str_section = hdr.data() + kStringSectionStart;
    for (int i = 0; i < kNumStrings; ++i) {
        const StringSpec& spec = kStringSpecs[i];
        std::string raw(str_section + spec.byte_offset, spec.length);
        if (isSentinelString(raw)) continue;
        trace.string_headers[spec.name] = trimTrailing(raw);
    }

    // Pull canonical fields.
    auto df = trace.float_headers.find("DELTA");
    if (df == trace.float_headers.end() || df->second <= 0.0) return false;
    trace.delta = df->second;
    auto bf = trace.float_headers.find("B");
    trace.begin_time = (bf == trace.float_headers.end()) ? 0.0 : bf->second;
    auto npi = trace.int_headers.find("NPTS");
    if (npi == trace.int_headers.end() || npi->second <= 0) return false;
    trace.npts = npi->second;

    // Read sample data.
    trace.samples.resize(static_cast<size_t>(trace.npts));
    std::vector<char> data(static_cast<size_t>(trace.npts) * 4);
    f.read(data.data(), static_cast<std::streamsize>(data.size()));
    if (f.gcount() != static_cast<std::streamsize>(data.size())) {
        trace.samples.clear();
        return false;
    }
    for (int i = 0; i < trace.npts; ++i) {
        trace.samples[i] = readFloat(data.data() + i * 4, swap);
    }
    trace.big_endian = swap;
    trace.valid = true;
    return true;
}

SACTrace resampleSAC(const SACTrace& trace, double target_delta)
{
    SACTrace out = trace;
    if (!trace.valid || target_delta <= 0.0 || trace.npts < 2)
        return out;

    const double t0 = trace.begin_time;
    const double t_end = t0 + (trace.npts - 1) * trace.delta;
    const int n_new = static_cast<int>(
        std::floor((t_end - t0) / target_delta)) + 1;
    if (n_new < 2) return out;

    out.delta = target_delta;
    out.npts = n_new;
    out.samples.assign(n_new, 0.0f);
    out.float_headers["DELTA"] = target_delta;
    out.int_headers["NPTS"] = n_new;
    out.float_headers["E"] = t0 + (n_new - 1) * target_delta;
    for (int i = 0; i < n_new; ++i) {
        const double t = t0 + i * target_delta;
        const double idx = (t - t0) / trace.delta;
        const int idx_lo = static_cast<int>(std::floor(idx));
        const int idx_hi = idx_lo + 1;
        if (idx_lo < 0) {
            out.samples[i] = trace.samples[0];
            continue;
        }
        if (idx_hi >= trace.npts) {
            out.samples[i] = trace.samples[trace.npts - 1];
            continue;
        }
        const double w = idx - idx_lo;
        out.samples[i] = static_cast<float>(
            (1.0 - w) * trace.samples[idx_lo] + w * trace.samples[idx_hi]);
    }
    return out;
}

SACTrace windowSAC(const SACTrace& trace, double t_begin, double t_end)
{
    SACTrace out = trace;
    if (!trace.valid || t_end <= t_begin || trace.npts < 1) {
        out.valid = false;
        return out;
    }
    const double t0 = trace.begin_time;
    const int i_lo = std::max(0,
        static_cast<int>(std::ceil((t_begin - t0) / trace.delta)));
    const int i_hi = std::min(trace.npts - 1,
        static_cast<int>(std::floor((t_end - t0) / trace.delta)));
    if (i_hi < i_lo) {
        out.valid = false;
        return out;
    }
    const int n = i_hi - i_lo + 1;
    out.samples.assign(trace.samples.begin() + i_lo,
                       trace.samples.begin() + i_lo + n);
    out.npts = n;
    out.begin_time = t0 + i_lo * trace.delta;
    out.int_headers["NPTS"] = n;
    out.float_headers["B"] = out.begin_time;
    out.float_headers["E"] = out.begin_time + (n - 1) * trace.delta;
    return out;
}

SACTrace taperSAC(const SACTrace& trace, double taper_fraction)
{
    SACTrace out = trace;
    if (!trace.valid || taper_fraction <= 0.0 || trace.npts < 2) return out;
    const int n = trace.npts;
    const int n_taper = static_cast<int>(std::floor(taper_fraction * n));
    if (n_taper < 1) return out;
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < n_taper; ++i) {
        const double w = 0.5 * (1.0 - std::cos(pi * i / n_taper));
        out.samples[i] = static_cast<float>(out.samples[i] * w);
        out.samples[n - 1 - i] = static_cast<float>(out.samples[n - 1 - i] * w);
    }
    return out;
}

SACTrace demeanDetrendSAC(const SACTrace& trace)
{
    SACTrace out = trace;
    if (!trace.valid || trace.npts < 2) return out;
    const int n = trace.npts;
    double mean_y = 0.0;
    for (int i = 0; i < n; ++i) mean_y += out.samples[i];
    mean_y /= n;
    // Linear LS on y vs i.
    double sum_xy = 0.0, sum_x2 = 0.0;
    const double mean_x = (n - 1) * 0.5;
    for (int i = 0; i < n; ++i) {
        const double dx = i - mean_x;
        sum_xy += dx * (out.samples[i] - mean_y);
        sum_x2 += dx * dx;
    }
    const double slope = (sum_x2 > 0.0) ? sum_xy / sum_x2 : 0.0;
    const double intercept = mean_y - slope * mean_x;
    for (int i = 0; i < n; ++i) {
        out.samples[i] = static_cast<float>(
            out.samples[i] - (slope * i + intercept));
    }
    return out;
}

}  // namespace io
}  // namespace FSRM
