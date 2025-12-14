#include "SeismometerNetwork.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace FSRM {

static constexpr float SAC_UNDEF = -12345.0f;
static constexpr int SAC_UNDEFI = -12345;

SeismometerNetwork::SeismometerNetwork(MPI_Comm comm) : comm_(comm) {
    MPI_Comm_rank(comm_, &rank_);
}

SeismometerNetwork::~SeismometerNetwork() {
    if (interp_) {
        DMInterpolationDestroy(&interp_);
    }
    if (interp_result_) {
        VecDestroy(&interp_result_);
    }
    if (disp_is_) {
        ISDestroy(&disp_is_);
    }
    if (disp_dm_) {
        DMDestroy(&disp_dm_);
    }
}

void SeismometerNetwork::setStations(std::vector<SeismometerSpec> stations) {
    stations_.clear();
    stations_.reserve(stations.size());
    for (auto& s : stations) {
        StationRuntime rt;
        rt.spec = std::move(s);
        stations_.push_back(std::move(rt));
    }
}

bool SeismometerNetwork::startsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::string SeismometerNetwork::trim(const std::string& s) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    size_t b = 0;
    while (b < s.size() && is_space(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && is_space(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::string SeismometerNetwork::upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return s;
}

std::vector<std::string> SeismometerNetwork::split(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        item = trim(item);
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

bool SeismometerNetwork::parseIsoUtc(const std::string& iso, int& year, int& month, int& day,
                                     int& hour, int& min, int& sec) {
    // Minimal ISO8601 parser: YYYY-MM-DDTHH:MM:SSZ
    if (iso.size() < 20) return false;
    char z = iso.back();
    if (z != 'Z' && z != 'z') return false;
    char T = iso[10];
    if (T != 'T' && T != 't') return false;
    try {
        year = std::stoi(iso.substr(0, 4));
        month = std::stoi(iso.substr(5, 2));
        day = std::stoi(iso.substr(8, 2));
        hour = std::stoi(iso.substr(11, 2));
        min = std::stoi(iso.substr(14, 2));
        sec = std::stoi(iso.substr(17, 2));
        return true;
    } catch (...) {
        return false;
    }
}

bool SeismometerNetwork::ymdToYday(int year, int month, int day, int& yday) {
    static const int mdays_norm[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    auto is_leap = [&](int y) {
        if (y % 400 == 0) return true;
        if (y % 100 == 0) return false;
        return (y % 4 == 0);
    };
    if (month < 1 || month > 12) return false;
    int md = mdays_norm[month - 1];
    if (month == 2 && is_leap(year)) md = 29;
    if (day < 1 || day > md) return false;
    int yd = 0;
    for (int m = 1; m < month; ++m) {
        int mlen = mdays_norm[m - 1];
        if (m == 2 && is_leap(year)) mlen = 29;
        yd += mlen;
    }
    yd += day;
    yday = yd;
    return true;
}

void SeismometerNetwork::addSecondsToYdayTime(int& year, int& yday, int& hour, int& min, int& sec,
                                              int& frac10k, double add_seconds) {
    // frac10k = 1e-4 sec units for MiniSEED BTime.
    long long total_10k = static_cast<long long>(std::llround(add_seconds * 10000.0));
    long long base_10k = ((static_cast<long long>(hour) * 3600LL +
                           static_cast<long long>(min) * 60LL +
                           static_cast<long long>(sec)) * 10000LL) +
                         static_cast<long long>(frac10k);
    long long t10k = base_10k + total_10k;
    if (t10k < 0) t10k = 0;

    long long day_10k = 86400LL * 10000LL;
    long long day_offset = t10k / day_10k;
    long long rem = t10k % day_10k;

    // Update time of day
    long long seconds = rem / 10000LL;
    frac10k = static_cast<int>(rem % 10000LL);
    hour = static_cast<int>(seconds / 3600LL);
    seconds %= 3600LL;
    min = static_cast<int>(seconds / 60LL);
    sec = static_cast<int>(seconds % 60LL);

    // Update yday; keep it simple (SAC/MSEED consumers mostly care about continuity).
    yday += static_cast<int>(day_offset);
    while (yday > 366) { yday -= 365; year += 1; }
}

bool SeismometerNetwork::geoToLocalMetersFallback(const GridConfig& grid_cfg,
                                                  const SeismometerSpec& spec,
                                                  double& x_m, double& y_m, double& z_m) {
    // Fallback assumes EPSG:4326 input, local origin given as lon/lat in GRID.local_origin_{x,y}.
    // Produces local tangent-plane meters.
    const double lon0 = grid_cfg.local_origin_x;
    const double lat0 = grid_cfg.local_origin_y;
    const double z0 = grid_cfg.local_origin_z;
    const double R = 6371000.0;
    const double deg = M_PI / 180.0;
    const double dlon = (spec.lon - lon0) * deg;
    const double dlat = (spec.lat - lat0) * deg;
    const double x = dlon * R * std::cos(lat0 * deg);
    const double y = dlat * R;
    const double z = spec.elev - z0;

    x_m = x;
    y_m = y;
    z_m = z;
    return true;
}

PetscErrorCode SeismometerNetwork::initialize(DM dm,
                                              const GridConfig& grid_cfg,
                                              const CoordinateSystemManager* coord_mgr) {
    PetscFunctionBeginUser;
    if (!out_cfg_.enabled || stations_.empty()) PetscFunctionReturn(0);

    dm_ = dm;

    // Determine minimum sample dt for downsampling in MonitorFunction
    min_sample_dt_ = std::numeric_limits<double>::infinity();
    for (const auto& st : stations_) {
        if (st.spec.sample_rate_hz > 0) {
            min_sample_dt_ = std::min(min_sample_dt_, 1.0 / st.spec.sample_rate_hz);
        }
    }
    if (!std::isfinite(min_sample_dt_)) min_sample_dt_ = 0.0;

    // Resolve station model coordinates
    for (auto& st : stations_) {
        const auto& s = st.spec;
        if (s.coord_type == SeismoCoordinateType::MODEL_XYZ) {
            st.xm = s.x; st.ym = s.y; st.zm = s.z;
        } else if (s.coord_type == SeismoCoordinateType::GRID_INDEX) {
            if (grid_cfg.mesh_type != MeshType::CARTESIAN) {
                if (rank_ == 0) {
                    PetscPrintf(comm_, "Warning: station %s uses GRID_INDEX but mesh_type is not CARTESIAN; ignoring.\n",
                                s.station.c_str());
                }
                st.xm = s.x; st.ym = s.y; st.zm = s.z;
            } else {
                const double dx = grid_cfg.Lx / std::max(1, grid_cfg.nx);
                const double dy = grid_cfg.Ly / std::max(1, grid_cfg.ny);
                const double dz = grid_cfg.Lz / std::max(1, grid_cfg.nz);
                if (s.grid_mode == GridIndexMode::NODE) {
                    const double ndx = grid_cfg.Lx / std::max(1, grid_cfg.nx - 1);
                    const double ndy = grid_cfg.Ly / std::max(1, grid_cfg.ny - 1);
                    const double ndz = grid_cfg.Lz / std::max(1, grid_cfg.nz - 1);
                    st.xm = grid_cfg.origin_x + s.i * ndx;
                    st.ym = grid_cfg.origin_y + s.j * ndy;
                    st.zm = grid_cfg.origin_z + s.k * ndz;
                } else {
                    st.xm = grid_cfg.origin_x + (s.i + 0.5) * dx;
                    st.ym = grid_cfg.origin_y + (s.j + 0.5) * dy;
                    st.zm = grid_cfg.origin_z + (s.k + 0.5) * dz;
                }
            }
        } else if (s.coord_type == SeismoCoordinateType::GEOGRAPHIC) {
            st.have_geo = true;
            st.stlo = s.lon;
            st.stla = s.lat;
            st.stel = s.elev;

            bool used_mgr = false;
            if (coord_mgr && coord_mgr->isConfigured()) {
                // If PROJ is not enabled, the internal transformer won't convert degrees->meters,
                // so only use coord_mgr when it is truly projecting (best-effort).
                GeoPoint input(s.lon, s.lat, s.elev);
                GeoPoint model = coord_mgr->toModelCoords(input);
                // Heuristic: if results are still in degree-like range, prefer fallback.
                if (std::abs(model.x) > 1000.0 || std::abs(model.y) > 1000.0) {
                    st.xm = model.x; st.ym = model.y; st.zm = model.z;
                    used_mgr = true;
                }
            }
            if (!used_mgr) {
                geoToLocalMetersFallback(grid_cfg, s, st.xm, st.ym, st.zm);
            }
        }
    }

    // Find displacement field index by name
    PetscInt nfields = 0;
    PetscCall(DMGetNumFields(dm_, &nfields));
    PetscInt disp_field = -1;
    for (PetscInt f = 0; f < nfields; ++f) {
        PetscObject obj = nullptr;
        PetscCall(DMGetField(dm_, f, nullptr, &obj));
        const char* nm = nullptr;
        if (obj) PetscCall(PetscObjectGetName(obj, &nm));
        if (nm && std::string(nm).find("displacement_") == 0) {
            disp_field = f;
            break;
        }
    }
    if (disp_field < 0) {
        if (rank_ == 0) {
            PetscPrintf(comm_, "Warning: seismometers enabled but no displacement_ field exists; disabling.\n");
        }
        out_cfg_.enabled = false;
        PetscFunctionReturn(0);
    }

    // Compute interpolation layout offsets (interp_result_ packs all field components)
    disp_offset_ = 0;
    dof_per_point_ = 0;
    disp_components_ = 3;
    PetscInt running = 0;
    for (PetscInt f = 0; f < nfields; ++f) {
        PetscObject obj = nullptr;
        PetscInt nc = 1;
        PetscCall(DMGetField(dm_, f, nullptr, &obj));
        if (obj) {
            PetscClassId cid = 0;
            PetscCall(PetscObjectGetClassId(obj, &cid));
            if (cid == PETSCFE_CLASSID) {
                PetscCall(PetscFEGetNumComponents((PetscFE)obj, &nc));
            }
        }
        if (f == disp_field) {
            disp_offset_ = running;
            disp_components_ = nc;
        }
        running += nc;
    }
    dof_per_point_ = running;

    // Setup interpolation
    PetscCall(DMInterpolationCreate(comm_, &interp_));
    PetscCall(DMInterpolationSetDim(interp_, 3));

    std::vector<PetscReal> pts;
    pts.reserve(stations_.size() * 3);
    for (const auto& st : stations_) {
        pts.push_back(static_cast<PetscReal>(st.xm));
        pts.push_back(static_cast<PetscReal>(st.ym));
        pts.push_back(static_cast<PetscReal>(st.zm));
    }

    PetscCall(DMInterpolationAddPoints(interp_, static_cast<PetscInt>(stations_.size()), pts.data()));
    PetscCall(DMInterpolationSetUp(interp_, dm_, PETSC_FALSE, PETSC_FALSE));
    PetscCall(DMInterpolationGetVector(interp_, &interp_result_));

    initialized_ = true;
    PetscFunctionReturn(0);
}

PetscErrorCode SeismometerNetwork::sample(double t, Vec U) {
    PetscFunctionBeginUser;
    if (!out_cfg_.enabled || !initialized_ || stations_.empty()) PetscFunctionReturn(0);

    // Lightweight downsampling: record only if enough time passed for the fastest instrument.
    if (min_sample_dt_ > 0.0 && (t - last_sample_time_) < (min_sample_dt_ - 1e-12)) {
        PetscFunctionReturn(0);
    }
    last_sample_time_ = t;

    // Interpolate full solution at station points (then pick displacement components by offset).
    Vec locU = nullptr;
    PetscCall(DMGetLocalVector(dm_, &locU));
    PetscCall(DMGlobalToLocalBegin(dm_, U, INSERT_VALUES, locU));
    PetscCall(DMGlobalToLocalEnd(dm_, U, INSERT_VALUES, locU));
    PetscCall(DMInterpolationEvaluate(interp_, dm_, locU, interp_result_));
    PetscCall(DMRestoreLocalVector(dm_, &locU));

    // Only rank 0 stores and writes traces (avoids duplicates).
    if (rank_ != 0) PetscFunctionReturn(0);

    const PetscScalar* a = nullptr;
    PetscCall(VecGetArrayRead(interp_result_, &a));
    PetscInt nloc = 0;
    PetscCall(VecGetLocalSize(interp_result_, &nloc));
    PetscInt stride = (stations_.empty()) ? 0 : (nloc / static_cast<PetscInt>(stations_.size()));
    if (stride <= 0) stride = dof_per_point_;
    // interp_result_ layout: point-major blocks, with displacement starting at disp_offset_
    for (size_t p = 0; p < stations_.size(); ++p) {
        auto& st = stations_[p];
        st.t.push_back(t);
        const PetscInt base = static_cast<PetscInt>(p) * stride + disp_offset_;
        const double ux = (base + 0 < nloc) ? static_cast<double>(a[base + 0]) : 0.0;
        const double uy = (base + 1 < nloc) ? static_cast<double>(a[base + 1]) : 0.0;
        const double uz = (base + 2 < nloc) ? static_cast<double>(a[base + 2]) : 0.0;
        st.disp.push_back({ux, uy, uz});
    }
    PetscCall(VecRestoreArrayRead(interp_result_, &a));

    PetscFunctionReturn(0);
}

static double medianDt(const std::vector<double>& t) {
    if (t.size() < 2) return 0.0;
    std::vector<double> d;
    d.reserve(t.size() - 1);
    for (size_t i = 1; i < t.size(); ++i) {
        double dt = t[i] - t[i - 1];
        if (dt > 0) d.push_back(dt);
    }
    if (d.empty()) return 0.0;
    std::nth_element(d.begin(), d.begin() + d.size()/2, d.end());
    return d[d.size()/2];
}

std::vector<float> SeismometerNetwork::deriveVelocity(const std::vector<double>& t,
                                                      const std::vector<std::array<double, 3>>& disp,
                                                      int comp) {
    std::vector<float> v;
    v.resize(disp.size(), 0.0f);
    if (disp.size() < 2) return v;
    for (size_t i = 1; i < disp.size(); ++i) {
        double dt = t[i] - t[i - 1];
        if (dt <= 0) continue;
        v[i] = static_cast<float>((disp[i][comp] - disp[i - 1][comp]) / dt);
    }
    // Backfill first
    v[0] = v[1];
    return v;
}

std::vector<float> SeismometerNetwork::deriveAcceleration(const std::vector<double>& t,
                                                          const std::vector<std::array<double, 3>>& disp,
                                                          int comp) {
    auto v = deriveVelocity(t, disp, comp);
    std::vector<float> a(v.size(), 0.0f);
    if (v.size() < 2) return a;
    for (size_t i = 1; i < v.size(); ++i) {
        double dt = t[i] - t[i - 1];
        if (dt <= 0) continue;
        a[i] = static_cast<float>((v[i] - v[i - 1]) / dt);
    }
    a[0] = a[1];
    return a;
}

void SeismometerNetwork::applySimpleHP(const std::vector<float>& x, double fs, double fc, std::vector<float>& y) {
    // First-order high-pass (discrete RC)
    y.assign(x.size(), 0.0f);
    if (x.empty() || fs <= 0 || fc <= 0) { y = x; return; }
    double dt = 1.0 / fs;
    double RC = 1.0 / (2.0 * M_PI * fc);
    double alpha = RC / (RC + dt);
    float yprev = 0.0f;
    float xprev = x[0];
    for (size_t i = 0; i < x.size(); ++i) {
        float xi = x[i];
        float yi = static_cast<float>(alpha * (yprev + xi - xprev));
        y[i] = yi;
        yprev = yi;
        xprev = xi;
    }
}

void SeismometerNetwork::applySimpleLP(const std::vector<float>& x, double fs, double fc, std::vector<float>& y) {
    // First-order low-pass (discrete RC)
    y.assign(x.size(), 0.0f);
    if (x.empty() || fs <= 0 || fc <= 0) { y = x; return; }
    double dt = 1.0 / fs;
    double RC = 1.0 / (2.0 * M_PI * fc);
    double alpha = dt / (RC + dt);
    float yprev = x[0];
    for (size_t i = 0; i < x.size(); ++i) {
        float xi = x[i];
        float yi = static_cast<float>(yprev + alpha * (xi - yprev));
        y[i] = yi;
        yprev = yi;
    }
}

void SeismometerNetwork::applyInstrument(const InstrumentPerformance& inst,
                                         double sample_rate_hz,
                                         std::vector<float>& data) {
    if (data.empty()) return;

    // Noise (in recorded units)
    if (inst.noise_std > 0.0) {
        std::mt19937_64 rng(1234567);
        std::normal_distribution<double> n(0.0, inst.noise_std);
        for (auto& v : data) v = static_cast<float>(v + n(rng));
    }

    // Simple band-limiting
    if (inst.highpass_corner_hz > 0.0) {
        std::vector<float> y;
        applySimpleHP(data, sample_rate_hz, inst.highpass_corner_hz, y);
        data.swap(y);
    }
    if (inst.lowpass_corner_hz > 0.0) {
        std::vector<float> y;
        applySimpleLP(data, sample_rate_hz, inst.lowpass_corner_hz, y);
        data.swap(y);
    }

    // Gain (counts-per-unit or pre-ADC scaling)
    if (inst.gain != 1.0) {
        for (auto& v : data) v = static_cast<float>(v * inst.gain);
    }

    // Clip in recorded units (after gain, before ADC)
    if (inst.clip > 0.0) {
        const float c = static_cast<float>(inst.clip);
        for (auto& v : data) {
            if (v > c) v = c;
            if (v < -c) v = -c;
        }
    }
}

void SeismometerNetwork::quantizeToCounts(const InstrumentPerformance& inst,
                                          std::vector<float>& data) {
    if (inst.adc_bits <= 0) return;
    if (inst.full_scale <= 0.0) return;
    const double max_count = static_cast<double>((1LL << (inst.adc_bits - 1)) - 1LL);
    const double fs = inst.full_scale;
    for (auto& v : data) {
        double c = (static_cast<double>(v) / fs) * max_count;
        if (c > max_count) c = max_count;
        if (c < -max_count) c = -max_count;
        v = static_cast<float>(std::llround(c));
    }
}

// SAC header layout (binary) – enough to be readable by ObsPy/SAC.
struct SacHeader {
    float f[70];
    int32_t i[40];
    char c[192];
};

static void sacSetString(char* dest, size_t n, const std::string& s) {
    std::string t = s;
    if (t.size() > n) t = t.substr(0, n);
    std::memset(dest, ' ', n);
    std::memcpy(dest, t.data(), t.size());
}

void SeismometerNetwork::writeSAC(const std::string& filename,
                                  const std::string& net, const std::string& sta,
                                  const std::string& loc, const std::string& chan,
                                  const StationRuntime& st,
                                  const std::vector<float>& data,
                                  double delta,
                                  const SeismometerOutputConfig& cfg) {
    (void)cfg;
    SacHeader h{};
    for (float& v : h.f) v = SAC_UNDEF;
    for (int32_t& v : h.i) v = SAC_UNDEFI;
    std::memset(h.c, ' ', sizeof(h.c));

    // Common SAC fields
    h.f[0] = static_cast<float>(delta);          // DELTA
    h.f[5] = 0.0f;                               // B
    h.f[6] = static_cast<float>((data.size() > 0) ? (delta * (data.size() - 1)) : 0.0); // E
    h.f[31] = st.have_geo ? static_cast<float>(st.stla) : SAC_UNDEF; // STLA
    h.f[32] = st.have_geo ? static_cast<float>(st.stlo) : SAC_UNDEF; // STLO
    h.f[33] = st.have_geo ? static_cast<float>(st.stel) : SAC_UNDEF; // STEL

    // Stats
    float depmin = std::numeric_limits<float>::infinity();
    float depmax = -std::numeric_limits<float>::infinity();
    double sum = 0.0;
    for (float v : data) {
        depmin = std::min(depmin, v);
        depmax = std::max(depmax, v);
        sum += v;
    }
    h.f[1] = (data.empty() ? SAC_UNDEF : depmin); // DEPMIN
    h.f[2] = (data.empty() ? SAC_UNDEF : depmax); // DEPMAX
    h.f[56] = (data.empty() ? SAC_UNDEF : static_cast<float>(sum / data.size())); // DEPMEN

    // Ints
    h.i[9]  = static_cast<int32_t>(data.size());  // NPTS
    h.i[6]  = 6;                                  // NVHDR
    h.i[15] = 1;                                  // IFTYPE = ITIME
    h.i[35] = 1;                                  // LEVEN

    // Reference time: leave unset (many tools accept), but set NZ* when possible.
    int y, mo, d, hh, mm, ss;
    if (parseIsoUtc(cfg.start_time_utc, y, mo, d, hh, mm, ss)) {
        int yday = 1;
        ymdToYday(y, mo, d, yday);
        h.i[0] = y;                // NZYEAR
        h.i[1] = yday;             // NZJDAY
        h.i[2] = hh;               // NZHOUR
        h.i[3] = mm;               // NZMIN
        h.i[4] = ss;               // NZSEC
        h.i[5] = 0;                // NZMSEC
    }

    // Strings: offsets per SAC standard
    // c[0..7]   KSTNM (8)
    // c[160..167] KCMPNM (8)
    // c[168..175] KNETWK (8)
    // c[176..183] KHOLE (8)
    sacSetString(&h.c[0], 8, sta);
    sacSetString(&h.c[160], 8, chan);
    sacSetString(&h.c[168], 8, net);
    sacSetString(&h.c[176], 8, loc);

    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) return;
    out.write(reinterpret_cast<const char*>(&h), sizeof(h));
    out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(float)));
}

// MiniSEED v2: write 512-byte records with encoding 4 (float32), blockette 1000.
static void write_be16(std::ostream& os, int16_t v) {
    uint16_t u = static_cast<uint16_t>(v);
    char b[2] = {static_cast<char>((u >> 8) & 0xFF), static_cast<char>(u & 0xFF)};
    os.write(b, 2);
}
static void write_be32(std::ostream& os, int32_t v) {
    uint32_t u = static_cast<uint32_t>(v);
    char b[4] = {
        static_cast<char>((u >> 24) & 0xFF),
        static_cast<char>((u >> 16) & 0xFF),
        static_cast<char>((u >> 8) & 0xFF),
        static_cast<char>(u & 0xFF)
    };
    os.write(b, 4);
}
static void write_befloat(std::ostream& os, float vf) {
    static_assert(sizeof(float) == 4, "float must be 4 bytes");
    uint32_t u;
    std::memcpy(&u, &vf, sizeof(u));
    write_be32(os, static_cast<int32_t>(u));
}

static void seedFixedString(char* dest, size_t n, const std::string& s) {
    std::memset(dest, ' ', n);
    std::string t = s;
    if (t.size() > n) t = t.substr(0, n);
    std::memcpy(dest, t.data(), t.size());
}

static void computeSeedSampRate(double sr, int16_t& factor, int16_t& mult) {
    // Simple representation: prefer integer factor.
    if (sr <= 0) { factor = 0; mult = 0; return; }
    double r = sr;
    int ir = static_cast<int>(std::llround(r));
    if (std::abs(r - ir) < 1e-6 && ir >= 1 && ir <= 32767) {
        factor = static_cast<int16_t>(ir);
        mult = 1;
        return;
    }
    // Fallback: represent dt as reciprocal factor (negative factor convention).
    double dt = 1.0 / sr;
    int idt = static_cast<int>(std::llround(dt));
    if (std::abs(dt - idt) < 1e-6 && idt >= 1 && idt <= 32767) {
        factor = -static_cast<int16_t>(idt);
        mult = 1;
        return;
    }
    factor = static_cast<int16_t>(std::min(32767, std::max(1, ir)));
    mult = 1;
}

void SeismometerNetwork::writeMiniSEED(const std::string& filename,
                                       const std::string& net, const std::string& sta,
                                       const std::string& loc, const std::string& chan,
                                       const StationRuntime& st,
                                       const std::vector<float>& data,
                                       double sample_rate_hz,
                                       const SeismometerOutputConfig& cfg) {
    // Base time
    int by, bmo, bd, bh, bmin, bs;
    if (!parseIsoUtc(cfg.start_time_utc, by, bmo, bd, bh, bmin, bs)) {
        by = 1970; bmo = 1; bd = 1; bh = 0; bmin = 0; bs = 0;
    }
    int yday = 1;
    ymdToYday(by, bmo, bd, yday);

    // Determine dt from sample rate (seconds)
    const double dt = (sample_rate_hz > 0.0) ? (1.0 / sample_rate_hz) : 0.0;

    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) return;

    const int record_len = 512;
    const int header_len = 56; // 48 FSDH + 8 blockette 1000
    const int bytes_per_sample = 4;
    const int max_samps = (record_len - header_len) / bytes_per_sample;

    int16_t sfactor = 0, smult = 0;
    computeSeedSampRate(sample_rate_hz, sfactor, smult);

    int seq = 1;
    size_t idx = 0;
    while (idx < data.size()) {
        size_t n = std::min(static_cast<size_t>(max_samps), data.size() - idx);

        // Compute record start time = base + t0 + idx*dt (approx)
        int ry = by, rj = yday, rh = bh, rmin = bmin, rs = bs, frac10k = 0;
        double add = (st.t.empty() ? 0.0 : st.t.front()) + (idx * dt);
        addSecondsToYdayTime(ry, rj, rh, rmin, rs, frac10k, add);

        // ---- Fixed Section of Data Header (48 bytes) ----
        // SEED sequence is 6 ASCII digits (no NUL). Wrap at 1,000,000.
        int seq6 = seq % 1000000;
        char seqno[6];
        for (int p = 5; p >= 0; --p) {
            seqno[p] = static_cast<char>('0' + (seq6 % 10));
            seq6 /= 10;
        }
        out.write(seqno, 6);
        out.put('D'); // data quality
        out.put(' '); // reserved

        char sta5[5], loc2[2], chan3[3], net2[2];
        seedFixedString(sta5, 5, sta);
        seedFixedString(loc2, 2, loc);
        seedFixedString(chan3, 3, chan);
        seedFixedString(net2, 2, net);
        out.write(sta5, 5);
        out.write(loc2, 2);
        out.write(chan3, 3);
        out.write(net2, 2);

        // BTime (10 bytes), big-endian fields
        write_be16(out, static_cast<int16_t>(ry));
        write_be16(out, static_cast<int16_t>(rj));
        out.put(static_cast<char>(rh));
        out.put(static_cast<char>(rmin));
        out.put(static_cast<char>(rs));
        out.put(static_cast<char>(0)); // unused
        write_be16(out, static_cast<int16_t>(frac10k)); // 1e-4 sec

        write_be16(out, static_cast<int16_t>(n));       // number of samples
        write_be16(out, sfactor);                        // sample rate factor
        write_be16(out, smult);                          // sample rate multiplier
        out.put(static_cast<char>(0));                   // activity flags
        out.put(static_cast<char>(0));                   // IO flags
        out.put(static_cast<char>(0));                   // data quality flags
        out.put(static_cast<char>(1));                   // num blockettes to follow
        write_be32(out, 0);                              // time correction
        write_be16(out, static_cast<int16_t>(header_len)); // begin data offset
        write_be16(out, static_cast<int16_t>(48));       // first blockette offset

        // ---- Blockette 1000 (8 bytes) ----
        write_be16(out, 1000);                           // blockette type
        write_be16(out, 0);                              // next blockette offset
        out.put(static_cast<char>(4));                   // encoding: 4 = IEEE float32
        out.put(static_cast<char>(1));                   // word order: 1 = big endian
        out.put(static_cast<char>(9));                   // record length: 2^9 = 512
        out.put(static_cast<char>(0));                   // reserved

        // ---- Data ----
        for (size_t k = 0; k < n; ++k) {
            write_befloat(out, data[idx + k]);
        }

        // Pad record
        const int written = header_len + static_cast<int>(n) * bytes_per_sample;
        const int pad = record_len - written;
        for (int p = 0; p < pad; ++p) out.put('\0');

        idx += n;
        seq += 1;
    }
}

PetscErrorCode SeismometerNetwork::finalizeAndWrite() const {
    PetscFunctionBeginUser;
    if (!out_cfg_.enabled || rank_ != 0 || stations_.empty()) PetscFunctionReturn(0);

    // Ensure output directory exists
    try {
        std::filesystem::create_directories(out_cfg_.output_dir);
    } catch (...) {
        // Ignore; file opens will fail if directory is missing/unwritable.
    }

    for (const auto& st : stations_) {
        if (st.t.size() < 2) continue;

        const double dt_med = medianDt(st.t);
        const double sr = (dt_med > 0) ? (1.0 / dt_med) : st.spec.sample_rate_hz;
        const double delta = (sr > 0) ? (1.0 / sr) : 0.0;

        for (int comp = 0; comp < 3; ++comp) {
            std::vector<float> series;
            if (st.spec.quantity == SeismoQuantity::DISPLACEMENT) {
                series.resize(st.disp.size());
                for (size_t i = 0; i < st.disp.size(); ++i) series[i] = static_cast<float>(st.disp[i][comp]);
            } else if (st.spec.quantity == SeismoQuantity::ACCELERATION) {
                series = deriveAcceleration(st.t, st.disp, comp);
            } else {
                series = deriveVelocity(st.t, st.disp, comp);
            }

            // Apply instrument performance and digitizer
            applyInstrument(st.spec.instrument, sr, series);
            quantizeToCounts(st.spec.instrument, series);

            const std::string chan = st.spec.channels[comp];
            const std::string base =
                out_cfg_.output_dir + "/" +
                st.spec.network + "." + st.spec.station + "." + st.spec.location + "." + chan;

            if (out_cfg_.write_sac) {
                writeSAC(base + ".sac", st.spec.network, st.spec.station, st.spec.location, chan, st, series, delta, out_cfg_);
            }
            if (out_cfg_.write_mseed) {
                writeMiniSEED(base + ".mseed", st.spec.network, st.spec.station, st.spec.location, chan, st, series, sr, out_cfg_);
            }
        }
    }

    PetscFunctionReturn(0);
}

} // namespace FSRM

