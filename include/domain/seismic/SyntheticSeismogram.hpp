/**
 * @file SyntheticSeismogram.hpp
 * @brief Regional synthetic seismogram generation for underground explosions
 *
 * Provides 1-D layered velocity models and a frequency-domain phase-summation
 * engine that generates velocity seismograms at regional distances
 * (100–2000 km).  The source is parameterised through the Mueller-Murphy
 * Reduced Displacement Potential model.
 *
 * Physics:
 *   - Regional phases: Pn, Pg, Sn, Lg (with pP / sP depth phases on P)
 *   - Geometric spreading: 1/r (body) and 1/sqrt(r) (surface)
 *   - Anelastic attenuation: exp(-pi f t / Q)
 *   - Mueller-Murphy omega-squared source spectrum
 *   - Frequency-domain Butterworth bandpass
 *   - Envelope shaping and frequency-dependent coda scattering
 *
 * References:
 *   Mueller & Murphy (1971) — Seismic characteristics of underground
 *       nuclear detonations
 *   Patton (1988) — Corner frequency scaling
 *
 * @author FSRM Development Team
 */

#ifndef SYNTHETIC_SEISMOGRAM_HPP
#define SYNTHETIC_SEISMOGRAM_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace FSRM
{

// ============================================================================
// 1-D Layered Velocity Model
// ============================================================================

/**
 * @brief A single layer in a 1-D velocity model.
 */
struct VelocityLayer
{
  double depth_top_km;   ///< Depth to top of layer (km)
  double vp_km_s;        ///< P-wave velocity (km/s)
  double vs_km_s;        ///< S-wave velocity (km/s)
  double rho_kg_m3;      ///< Density (kg/m^3)
  double Qp;             ///< P-wave quality factor
  double Qs;             ///< S-wave quality factor
};

/**
 * @brief 1-D layered velocity model for regional propagation.
 *
 * Stores a stack of constant-velocity layers and derived crustal averages.
 */
class VelocityModel1D
{
public:
  VelocityModel1D() = default;

  /**
   * @brief Construct from a layer table.
   *
   * @param layers         Ordered layer list (shallow → deep)
   * @param moho_depth_km  Moho depth (km)
   * @param pn_velocity    Pn phase velocity (km/s)
   * @param model_name     Human-readable identifier
   */
  VelocityModel1D(const std::vector<VelocityLayer>& layers,
                  double moho_depth_km,
                  double pn_velocity,
                  const std::string& model_name = "custom");

  /// Name of the model.
  const std::string& name() const { return name_; }

  /// Number of layers.
  std::size_t size() const { return layers_.size(); }

  /// Moho depth (km).
  double mohoDepth() const { return moho_depth_; }

  /// Pn velocity (km/s).
  double pnVelocity() const { return pn_velocity_; }

  /// Average crustal P-wave velocity (km/s).
  double avgCrustalVp() const { return avg_crustal_vp_; }

  /// Average crustal S-wave velocity (km/s).
  double avgCrustalVs() const { return avg_crustal_vs_; }

  /**
   * @brief Get elastic properties at depth.
   *
   * @param z_km  Depth (km)
   * @param vp    P-velocity (km/s)
   * @param vs    S-velocity (km/s)
   * @param rho   Density (kg/m^3)
   * @param qp    P quality factor
   * @param qs    S quality factor
   */
  void get(double z_km,
           double& vp, double& vs, double& rho,
           double& qp, double& qs) const;

  /// Convenience: Sn velocity (km/s) — deepest layer S velocity.
  double snVelocity() const;

  // ---- Built-in models ---------------------------------------------------
  static VelocityModel1D punggyeRi();
  static VelocityModel1D lopNor();
  static VelocityModel1D nts();
  static VelocityModel1D genericGranite();

  /// Retrieve a built-in model by name (punggye_ri, lop_nor, nts, generic).
  static VelocityModel1D byName(const std::string& model_name);

private:
  std::vector<VelocityLayer> layers_;
  double moho_depth_  = 35.0;
  double pn_velocity_ = 8.0;
  std::string name_   = "base";

  double avg_crustal_vp_ = 6.0;
  double avg_crustal_vs_ = 3.5;

  void computeAverages();
  std::size_t layerIndex(double z_km) const;
};

// ============================================================================
// Travel-time and propagation helpers
// ============================================================================

/**
 * @brief Descriptor for one seismic phase.
 */
struct PhaseInfo
{
  double travel_time_s;   ///< First-arrival time (s)
  double velocity_km_s;   ///< Group / apparent velocity (km/s)
  std::string wave_type;  ///< "P" or "S"
};

/**
 * @brief Compute regional travel times for Pn, Pg, Sn, Lg.
 */
std::map<std::string, PhaseInfo>
regionalTravelTimes(double dist_km, const VelocityModel1D& model);

/// Geometric spreading factor (body → 1/r, surface → 1/sqrt(r)).
double geometricSpreading(double dist_km, const std::string& wave_type);

/// Anelastic attenuation: exp(-pi * f * t / Q).
double anelasticAttenuation(double freq, double travel_time_s, double Q);

/// Depth-phase delays (pP, sP) relative to direct P.
void depthPhaseDelays(double depth_km, const VelocityModel1D& model,
                      double& pP_delay, double& sP_delay);

// ============================================================================
// Synthetic-seismogram generator
// ============================================================================

/**
 * @brief Configuration for generating a single synthetic seismogram.
 */
struct SyntheticConfig
{
  double dist_km            = 500.0;   ///< Epicentral distance (km)
  double yield_kt           = 100.0;   ///< Yield (kilotons)
  double decoupling_factor  = 1.0;     ///< 1 = fully coupled
  double depth_km           = 0.76;    ///< Source depth (km)
  double dt                 = 0.05;    ///< Sample interval (s)
  double duration           = 500.0;   ///< Record length (s)
  double fmin               = 0.5;     ///< Low bandpass corner (Hz)
  double fmax               = 8.0;     ///< High bandpass corner (Hz)
  uint32_t random_seed      = 0;       ///< 0 = derive from distance
};

/**
 * @brief Result of a synthetic seismogram computation.
 */
struct SyntheticResult
{
  std::vector<double> time;      ///< Time axis (s)
  std::vector<double> velocity;  ///< Velocity waveform
  int npts = 0;                  ///< Number of samples
  double dt = 0.0;               ///< Sample interval (s)
};

/**
 * @brief Generate regional synthetic velocity seismograms.
 *
 * Orchestrates Mueller-Murphy source, regional phase travel times,
 * geometric spreading, Q attenuation, depth phases, and coda scattering
 * to produce a complete far-field velocity seismogram.
 *
 * Typical usage:
 * @code
 *   auto model = VelocityModel1D::punggyeRi();
 *   SyntheticSeismogramGenerator gen(model);
 *
 *   SyntheticConfig cfg;
 *   cfg.dist_km  = 370.0;
 *   cfg.yield_kt = 250.0;
 *
 *   auto result = gen.generate(cfg);
 * @endcode
 */
class SyntheticSeismogramGenerator
{
public:
  explicit SyntheticSeismogramGenerator(const VelocityModel1D& model);

  /**
   * @brief Generate one synthetic seismogram.
   *
   * @param cfg  Configuration (distance, yield, filters, …)
   * @return     Time series and metadata.
   */
  SyntheticResult generate(const SyntheticConfig& cfg) const;

  /**
   * @brief Generate seismograms at multiple distances.
   *
   * @param distances_km  Vector of epicentral distances
   * @param cfg           Base configuration (dist_km is overridden)
   * @return              One result per distance.
   */
  std::vector<SyntheticResult>
  generateBatch(const std::vector<double>& distances_km,
                const SyntheticConfig& cfg) const;

  /// Replace the velocity model.
  void setModel(const VelocityModel1D& model) { model_ = model; }

  /// Read-only access to current model.
  const VelocityModel1D& model() const { return model_; }

private:
  VelocityModel1D model_;

  // Internal helpers
  void applyPhase(const std::string& phase_name,
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
                  std::vector<double>& vel) const;
};

} // namespace FSRM

#endif // SYNTHETIC_SEISMOGRAM_HPP
