/**
 * @file SyntheticSeismogram.cpp
 * @brief Implementation of regional synthetic seismogram generation.
 *
 * Translates the Python propagation module into C++ for high-throughput
 * waveform synthesis (e.g. inversion grids, Monte-Carlo uncertainty).
 */

#include "domain/seismic/SyntheticSeismogram.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

// Minimal in-place real-FFT: Cooley-Tukey radix-2 DIT.
// For production builds this should be replaced by FFTW / MKL; this keeps the
// translation self-contained and dependency-free.
namespace
{

// ─────────────────────────────────────────────────────────────────────────────
// Bit-reverse permutation
// ─────────────────────────────────────────────────────────────────────────────
static void bitReverse(std::vector<std::complex<double>>& a)
{
  int n = static_cast<int>(a.size());
  for (int i = 1, j = 0; i < n; ++i)
  {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j)
      std::swap(a[i], a[j]);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// In-place complex FFT (forward = true) / IFFT (forward = false)
// ─────────────────────────────────────────────────────────────────────────────
static void fft(std::vector<std::complex<double>>& a, bool forward)
{
  int n = static_cast<int>(a.size());
  if (n == 0)
    return;
  bitReverse(a);

  for (int len = 2; len <= n; len <<= 1)
  {
    double ang = 2.0 * M_PI / len * (forward ? -1.0 : 1.0);
    std::complex<double> wlen(std::cos(ang), std::sin(ang));
    for (int i = 0; i < n; i += len)
    {
      std::complex<double> w(1.0, 0.0);
      for (int j = 0; j < len / 2; ++j)
      {
        std::complex<double> u = a[i + j];
        std::complex<double> v = a[i + j + len / 2] * w;
        a[i + j]             = u + v;
        a[i + j + len / 2]   = u - v;
        w *= wlen;
      }
    }
  }

  if (!forward)
  {
    for (auto& x : a)
      x /= static_cast<double>(n);
  }
}

// Next power of two >= n.
static int nextPow2(int n)
{
  int p = 1;
  while (p < n)
    p <<= 1;
  return p;
}

// ─────────────────────────────────────────────────────────────────────────────
// Real-to-complex forward FFT.  Returns (N/2+1) complex bins.
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<std::complex<double>> rfft(const std::vector<double>& x, int N)
{
  std::vector<std::complex<double>> c(N);
  int nx = static_cast<int>(x.size());
  for (int i = 0; i < N; ++i)
    c[i] = (i < nx) ? std::complex<double>(x[i], 0.0)
                     : std::complex<double>(0.0, 0.0);
  fft(c, true);

  int M = N / 2 + 1;
  std::vector<std::complex<double>> out(M);
  for (int i = 0; i < M; ++i)
    out[i] = c[i];
  return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Complex-to-real inverse FFT.  Input has (N/2+1) bins.
// ─────────────────────────────────────────────────────────────────────────────
static std::vector<double> irfft(const std::vector<std::complex<double>>& X, int N)
{
  std::vector<std::complex<double>> c(N);
  int M = static_cast<int>(X.size());
  for (int i = 0; i < M && i < N; ++i)
    c[i] = X[i];
  // Hermitian symmetry
  for (int i = M; i < N; ++i)
    c[i] = std::conj(c[N - i]);
  fft(c, false);

  std::vector<double> out(N);
  for (int i = 0; i < N; ++i)
    out[i] = c[i].real();
  return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Mueller-Murphy source spectrum (displacement).
// ─────────────────────────────────────────────────────────────────────────────
static double muellerMurphySpectrum(double f, double M0, double fc,
                                    double overshoot = 1.1)
{
  return overshoot * M0 / (1.0 + (f / fc) * (f / fc));
}

// ─────────────────────────────────────────────────────────────────────────────
// Patton corner frequency.
// ─────────────────────────────────────────────────────────────────────────────
static double pattonCornerFreq(double yield_kt)
{
  return 2.5 * std::pow(yield_kt, -1.0 / 3.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Simple xorshift32 PRNG (deterministic, fast).
// ─────────────────────────────────────────────────────────────────────────────
struct Xorshift32
{
  uint32_t state;
  explicit Xorshift32(uint32_t seed) : state(seed ? seed : 1u) {}
  uint32_t next()
  {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
  }
  // Uniform in [-1, 1)
  double uniform()
  {
    return static_cast<double>(next()) / 2147483648.0 - 1.0;
  }
  // Approximate normal (Box-Muller).
  double randn()
  {
    double u1, u2;
    do { u1 = (next() + 0.5) / 4294967296.0; } while (u1 <= 0.0);
    u2 = (next() + 0.5) / 4294967296.0;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
  }
};

} // anonymous namespace

// ============================================================================
// FSRM namespace implementation
// ============================================================================
namespace FSRM
{

// --------------------------------------------------------------------------
// VelocityModel1D
// --------------------------------------------------------------------------

VelocityModel1D::VelocityModel1D(const std::vector<VelocityLayer>& layers,
                                 double moho_depth_km,
                                 double pn_velocity,
                                 const std::string& model_name)
  : layers_(layers),
    moho_depth_(moho_depth_km),
    pn_velocity_(pn_velocity),
    name_(model_name)
{
  computeAverages();
}

void VelocityModel1D::computeAverages()
{
  if (layers_.empty())
    return;

  // Average all crustal layers (exclude deepest = mantle).
  std::size_t n_crust = layers_.size() > 1 ? layers_.size() - 1 : 1;
  double sum_vp = 0.0, sum_vs = 0.0;
  for (std::size_t i = 0; i < n_crust; ++i)
  {
    sum_vp += layers_[i].vp_km_s;
    sum_vs += layers_[i].vs_km_s;
  }
  avg_crustal_vp_ = sum_vp / static_cast<double>(n_crust);
  avg_crustal_vs_ = sum_vs / static_cast<double>(n_crust);
}

std::size_t VelocityModel1D::layerIndex(double z_km) const
{
  for (std::size_t i = layers_.size(); i-- > 0; )
  {
    if (z_km >= layers_[i].depth_top_km)
      return i;
  }
  return 0;
}

void VelocityModel1D::get(double z_km,
                           double& vp, double& vs, double& rho,
                           double& qp, double& qs) const
{
  if (layers_.empty())
  {
    vp = 6.0; vs = 3.5; rho = 2700; qp = 500; qs = 250;
    return;
  }
  auto idx = layerIndex(z_km);
  const auto& L = layers_[idx];
  vp  = L.vp_km_s;
  vs  = L.vs_km_s;
  rho = L.rho_kg_m3;
  qp  = L.Qp;
  qs  = L.Qs;
}

double VelocityModel1D::snVelocity() const
{
  if (layers_.empty())
    return 4.5;
  return layers_.back().vs_km_s;
}

// ---- Built-in models -----------------------------------------------------

VelocityModel1D VelocityModel1D::punggyeRi()
{
  return VelocityModel1D(
    {
      { 0.00, 3.50, 2.00, 2100,   80,   40},
      { 0.10, 5.60, 3.25, 2700,  350,  175},
      { 5.00, 6.10, 3.55, 2800,  500,  250},
      {18.00, 6.50, 3.75, 2900,  600,  300},
      {28.00, 7.00, 4.05, 3100,  700,  350},
      {33.00, 8.05, 4.55, 3350, 1000,  500},
    },
    33.0, 8.05, "punggye_ri");
}

VelocityModel1D VelocityModel1D::lopNor()
{
  return VelocityModel1D(
    {
      { 0.00, 4.00, 2.20, 2300,  100,   50},
      { 0.50, 5.80, 3.35, 2650,  300,  150},
      { 8.00, 6.20, 3.58, 2800,  500,  250},
      {20.00, 6.60, 3.81, 2950,  600,  300},
      {35.00, 7.00, 4.05, 3100,  700,  350},
      {48.00, 8.10, 4.60, 3400, 1000,  500},
    },
    48.0, 8.10, "lop_nor");
}

VelocityModel1D VelocityModel1D::nts()
{
  return VelocityModel1D(
    {
      { 0.00, 3.50, 2.00, 2200,  100,   50},
      { 0.50, 5.00, 2.90, 2500,  200,  100},
      { 3.00, 6.00, 3.50, 2700,  400,  200},
      {15.00, 6.30, 3.65, 2850,  500,  250},
      {25.00, 6.80, 3.95, 3000,  650,  325},
      {30.00, 7.80, 4.45, 3300,  900,  450},
    },
    30.0, 7.80, "nts");
}

VelocityModel1D VelocityModel1D::genericGranite()
{
  return VelocityModel1D(
    {
      { 0.00, 5.80, 3.36, 2700,  300,  150},
      { 5.00, 6.10, 3.55, 2800,  500,  250},
      {15.00, 6.50, 3.75, 2900,  600,  300},
      {25.00, 7.00, 4.05, 3100,  700,  350},
      {35.00, 8.10, 4.55, 3350, 1000,  500},
    },
    35.0, 8.10, "generic");
}

VelocityModel1D VelocityModel1D::byName(const std::string& model_name)
{
  if (model_name == "punggye_ri")
    return punggyeRi();
  if (model_name == "lop_nor")
    return lopNor();
  if (model_name == "nts")
    return nts();
  if (model_name == "generic")
    return genericGranite();
  throw std::invalid_argument(
    "Unknown velocity model: " + model_name +
    ". Available: punggye_ri, lop_nor, nts, generic");
}

// --------------------------------------------------------------------------
// Free-function propagation helpers
// --------------------------------------------------------------------------

std::map<std::string, PhaseInfo>
regionalTravelTimes(double dist_km, const VelocityModel1D& model)
{
  double moho     = model.mohoDepth();
  double pn       = model.pnVelocity();
  double vp_crust = model.avgCrustalVp();
  double vs_crust = model.avgCrustalVs();

  // Pn: head wave refracted along Moho
  double cos_ic_p = vp_crust / pn;
  double t_pn = 2.0 * moho * cos_ic_p / vp_crust + dist_km / pn;

  // Pg: direct crustal P
  double t_pg = dist_km / vp_crust;

  // Sn: S refraction along Moho
  double sn_vel   = model.snVelocity();
  double cos_ic_s = std::min(vs_crust / sn_vel, 0.999);
  double t_sn = 2.0 * moho * cos_ic_s / vs_crust + dist_km / sn_vel;

  // Lg: crustal guided S / surface wave
  double lg_vel = 3.5;
  double t_lg  = dist_km / lg_vel;

  return {
    {"Pn", {t_pn, pn,       "P"}},
    {"Pg", {t_pg, vp_crust, "P"}},
    {"Sn", {t_sn, sn_vel,   "S"}},
    {"Lg", {t_lg, lg_vel,   "S"}},
  };
}

double geometricSpreading(double dist_km, const std::string& wave_type)
{
  double r = std::max(dist_km, 0.1);
  if (wave_type == "surface")
    return 1.0 / std::sqrt(r);
  return 1.0 / r;
}

double anelasticAttenuation(double freq, double travel_time_s, double Q)
{
  return std::exp(-M_PI * freq * travel_time_s / Q);
}

void depthPhaseDelays(double depth_km, const VelocityModel1D& model,
                      double& pP_delay, double& sP_delay)
{
  double vp, vs, rho, qp, qs;
  model.get(depth_km, vp, vs, rho, qp, qs);
  pP_delay = 2.0 * depth_km / vp;
  sP_delay = depth_km / vs + depth_km / vp;
}

// --------------------------------------------------------------------------
// SyntheticSeismogramGenerator
// --------------------------------------------------------------------------

SyntheticSeismogramGenerator::SyntheticSeismogramGenerator(
  const VelocityModel1D& model)
  : model_(model)
{
}

void SyntheticSeismogramGenerator::applyPhase(
  const std::string& phase_name,
  const PhaseInfo& info,
  const SyntheticConfig& cfg,
  int npts,
  const std::vector<double>& freq,
  const std::vector<double>& bandpass,
  double M0,
  double fc,
  double pP_delay,
  double sP_delay,
  uint32_t seed,
  std::vector<double>& vel) const
{
  int M = static_cast<int>(freq.size());  // N/2+1
  int N = nextPow2(npts);
  double dt = cfg.dt;

  // Phase-specific parameters
  static const std::map<std::string, double> phase_weights =
    {{"Pn", 1.0}, {"Pg", 0.6}, {"Sn", 0.35}, {"Lg", 2.5}};
  static const std::map<std::string, double> q_base =
    {{"Pn", 800}, {"Pg", 400}, {"Sn", 300}, {"Lg", 250}};
  static const std::map<std::string, std::string> spreading_type =
    {{"Pn", "body"}, {"Pg", "body"}, {"Sn", "body"}, {"Lg", "surface"}};
  static const std::map<std::string, double> env_duration =
    {{"Pn", 12.0}, {"Pg", 18.0}, {"Sn", 22.0}, {"Lg", 80.0}};
  static const std::map<std::string, double> rise_time =
    {{"Pn", 0.8}, {"Pg", 1.2}, {"Sn", 2.5}, {"Lg", 4.0}};
  static const std::map<std::string, double> coda_Q =
    {{"Pn", 600}, {"Pg", 400}, {"Sn", 300}, {"Lg", 250}};

  double t_arr = info.travel_time_s;
  if (t_arr < 0 || t_arr >= cfg.duration)
    return;

  auto lookup = [](const std::map<std::string, double>& m,
                   const std::string& k, double def) -> double {
    auto it = m.find(k);
    return it != m.end() ? it->second : def;
  };
  auto lookup_s = [](const std::map<std::string, std::string>& m,
                     const std::string& k,
                     const std::string& def) -> std::string {
    auto it = m.find(k);
    return it != m.end() ? it->second : def;
  };

  double weight = lookup(phase_weights, phase_name, 1.0);
  double Q = lookup(q_base, phase_name, 400.0) *
             (1.0 + 0.15 * std::log10(std::max(cfg.dist_km, 100.0) / 100.0));
  double spr = geometricSpreading(
    cfg.dist_km, lookup_s(spreading_type, phase_name, "body"));

  // Build spectral representation of this phase
  std::vector<std::complex<double>> spec(M);
  for (int i = 0; i < M; ++i)
  {
    double f = freq[i];
    double src = muellerMurphySpectrum(f, M0, fc);
    double att = anelasticAttenuation(f, t_arr, Q);
    double amp = src * att * spr * weight;
    double omega = 2.0 * M_PI * f;

    // Displacement → velocity (* omega), apply bandpass
    double spec_vel = amp * omega * bandpass[i];

    // Phase shift for arrival time
    double angle = -omega * t_arr;
    spec[i] = std::complex<double>(spec_vel * std::cos(angle),
                                   spec_vel * std::sin(angle));
  }

  // For P phases: add depth phases (pP, sP)
  std::vector<double> wavelet;
  if (info.wave_type == "P")
  {
    // pP reflection coefficient ~ -0.85, sP ~ +0.45
    std::vector<std::complex<double>> pP_spec(M), sP_spec(M);
    for (int i = 0; i < M; ++i)
    {
      double omega = 2.0 * M_PI * freq[i];
      double angle_pP = -omega * (t_arr + pP_delay);
      double angle_sP = -omega * (t_arr + sP_delay);
      double f = freq[i];
      double src = muellerMurphySpectrum(f, M0, fc);
      double att = anelasticAttenuation(f, t_arr, Q);
      double s_vel = src * att * spr * weight * omega * bandpass[i];

      pP_spec[i] = std::complex<double>(
        s_vel * (-0.85) * std::cos(angle_pP),
        s_vel * (-0.85) * std::sin(angle_pP));
      sP_spec[i] = std::complex<double>(
        s_vel * 0.45 * std::cos(angle_sP),
        s_vel * 0.45 * std::sin(angle_sP));
    }

    // Sum direct + pP + sP
    std::vector<std::complex<double>> total(M);
    for (int i = 0; i < M; ++i)
      total[i] = spec[i] + pP_spec[i] + sP_spec[i];
    wavelet = irfft(total, N);
  }
  else
  {
    // S phases: random phase perturbation
    Xorshift32 rng(seed ^ static_cast<uint32_t>(
      std::hash<std::string>{}(phase_name)));
    for (int i = 0; i < M; ++i)
    {
      double phi = rng.uniform() * M_PI / 4.0;
      std::complex<double> rp(std::cos(phi), std::sin(phi));
      spec[i] *= rp;
    }
    wavelet = irfft(spec, N);
  }

  // Truncate to npts
  wavelet.resize(npts, 0.0);

  // Envelope shaping
  double ed = lookup(env_duration, phase_name, 15.0);
  double rt = lookup(rise_time, phase_name, 2.0);
  double cQ = lookup(coda_Q, phase_name, 300.0);
  double f_mean = 0.5 * (cfg.fmin + cfg.fmax);
  int idx_arr = static_cast<int>(t_arr / dt);

  for (int j = 0; j < npts; ++j)
  {
    double t_rel = (j - idx_arr) * dt;
    if (t_rel < -rt)
    {
      wavelet[j] = 0.0;
    }
    else if (t_rel < 0.0)
    {
      wavelet[j] *= std::exp(-0.5 * (t_rel / (rt * 0.3)) * (t_rel / (rt * 0.3)));
    }
    else if (t_rel < rt)
    {
      double ramp = std::min(1.0, t_rel / (rt * 0.6));
      wavelet[j] *= ramp;
    }
    else if (t_rel >= ed * 0.3)
    {
      wavelet[j] *= std::exp(-M_PI * f_mean * (t_rel - ed * 0.3) / cQ);
    }
  }

  // Coda scattering noise
  Xorshift32 coda_rng(seed ^ static_cast<uint32_t>(
    std::hash<std::string>{}(phase_name + "coda")));
  std::vector<double> coda(npts);
  for (int j = 0; j < npts; ++j)
    coda[j] = coda_rng.randn();

  // Simple bandpass on coda (2nd-order Butterworth approximation via
  // frequency-domain windowing avoids pulling in a filter library).
  {
    int Nc = nextPow2(npts);
    auto C = rfft(coda, Nc);
    double df = 1.0 / (Nc * dt);
    int Mc = Nc / 2 + 1;
    for (int i = 0; i < Mc; ++i)
    {
      double f = i * df;
      double lo = f / cfg.fmin;
      double hi = f / cfg.fmax;
      // 2nd-order Butterworth magnitude
      double H = 1.0 / std::sqrt(1.0 + 1.0 / (lo * lo * lo * lo)) *
                 1.0 / std::sqrt(1.0 + hi * hi * hi * hi);
      C[i] *= H;
    }
    coda = irfft(C, Nc);
    coda.resize(npts, 0.0);
  }

  // Envelope coda
  for (int j = 0; j < npts; ++j)
  {
    double t_rel = (j - idx_arr) * dt;
    if (t_rel < rt * 1.5)
    {
      coda[j] = 0.0;
    }
    else if (t_rel < ed * 0.5)
    {
      coda[j] *= (1.0 - std::exp(-(t_rel - rt * 1.5) / (ed * 0.15))) * 0.35;
    }
    else if (t_rel < ed)
    {
      coda[j] *= std::exp(-(t_rel - ed * 0.5) / (ed * 0.4)) * 0.35;
    }
    else
    {
      coda[j] *= std::exp(-(t_rel - ed) / (ed * 0.25)) * 0.1;
    }
  }

  // Scale coda relative to wavelet peak
  double peak_w = 0.0;
  for (int j = 0; j < npts; ++j)
    peak_w = std::max(peak_w, std::abs(wavelet[j]));
  if (peak_w > 0.0)
  {
    for (int j = 0; j < npts; ++j)
      wavelet[j] += coda[j] * peak_w * 0.25;
  }

  // Accumulate into output
  for (int j = 0; j < npts; ++j)
    vel[j] += wavelet[j];
}

SyntheticResult
SyntheticSeismogramGenerator::generate(const SyntheticConfig& cfg) const
{
  int npts = static_cast<int>(cfg.duration / cfg.dt);
  int N = nextPow2(npts);  // FFT length
  int M = N / 2 + 1;       // Number of frequency bins

  // Frequency axis
  std::vector<double> freq(M);
  double df = 1.0 / (N * cfg.dt);
  for (int i = 0; i < M; ++i)
    freq[i] = std::max(i * df, 1e-10);

  // Cosine-taper bandpass
  std::vector<double> bp(M, 0.0);
  for (int i = 0; i < M; ++i)
  {
    double f = freq[i];
    if (f >= cfg.fmin && f <= cfg.fmax)
    {
      bp[i] = 1.0;
    }
    if (f >= cfg.fmin * 0.5 && f < cfg.fmin)
    {
      bp[i] = 0.5 * (1.0 - std::cos(
        M_PI * (f - cfg.fmin * 0.5) / (cfg.fmin * 0.5)));
    }
    if (f > cfg.fmax && f <= cfg.fmax * 1.25)
    {
      bp[i] = 0.5 * (1.0 + std::cos(
        M_PI * (f - cfg.fmax) / (cfg.fmax * 0.25)));
    }
  }

  // Source parameters
  double fc = pattonCornerFreq(cfg.yield_kt);
  double M0_coupled = std::pow(10.0, 17.0 + std::log10(
    std::max(cfg.yield_kt, 1e-6)));
  double M0 = M0_coupled / cfg.decoupling_factor;

  // Depth-phase delays
  double pP_delay = 0.0, sP_delay = 0.0;
  depthPhaseDelays(cfg.depth_km, model_, pP_delay, sP_delay);

  // Regional phases
  auto phases = regionalTravelTimes(cfg.dist_km, model_);

  // Seed
  uint32_t seed = cfg.random_seed;
  if (seed == 0)
    seed = static_cast<uint32_t>(cfg.dist_km * 1000.0) ^ 0xDEADBEEFu;

  // Accumulate velocity waveform
  std::vector<double> vel(npts, 0.0);

  for (auto& [ph_name, ph_info] : phases)
  {
    applyPhase(ph_name, ph_info, cfg, npts, freq, bp,
               M0, fc, pP_delay, sP_delay, seed, vel);
  }

  // Build result
  SyntheticResult result;
  result.npts = npts;
  result.dt = cfg.dt;
  result.time.resize(npts);
  for (int i = 0; i < npts; ++i)
    result.time[i] = i * cfg.dt;
  result.velocity = std::move(vel);

  return result;
}

std::vector<SyntheticResult>
SyntheticSeismogramGenerator::generateBatch(
  const std::vector<double>& distances_km,
  const SyntheticConfig& base_cfg) const
{
  std::vector<SyntheticResult> results;
  results.reserve(distances_km.size());

  SyntheticConfig cfg = base_cfg;
  for (double d : distances_km)
  {
    cfg.dist_km = d;
    results.push_back(generate(cfg));
  }

  return results;
}

} // namespace FSRM
