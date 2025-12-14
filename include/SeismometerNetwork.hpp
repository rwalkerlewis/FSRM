#ifndef SEISMOMETER_NETWORK_HPP
#define SEISMOMETER_NETWORK_HPP

#include "FSRM.hpp"
#include "CoordinateSystem.hpp"

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace FSRM {

enum class SeismoQuantity {
    DISPLACEMENT,
    VELOCITY,
    ACCELERATION
};

enum class SeismoCoordinateType {
    MODEL_XYZ,     // x,y,z in model coordinates (meters)
    GRID_INDEX,    // i,j,k in structured grid indices
    GEOGRAPHIC     // lon,lat,elev in degrees/meters (input CRS, typically EPSG:4326)
};

enum class GridIndexMode {
    CELL_CENTER,
    NODE
};

struct InstrumentPerformance {
    // Simple band-limited instrument model applied in the recorded quantity's units
    double highpass_corner_hz = 0.0;   // 0 => off
    double lowpass_corner_hz = 0.0;    // 0 => off
    double gain = 1.0;                // counts per unit (if adc_bits==0), else applied before ADC
    double noise_std = 0.0;           // Gaussian white noise std (in recorded units)

    // Digitizer/ADC (optional)
    int adc_bits = 0;                 // 0 => disabled
    double full_scale = 0.0;          // peak full-scale amplitude in recorded units (required if adc_bits>0)
    double clip = 0.0;                // optional hard clip in recorded units (0 => no extra clip)
};

struct SeismometerSpec {
    // SEED metadata
    std::string network = "XX";
    std::string station = "STAT";
    std::string location = "00";
    std::array<std::string, 3> channels = {"BHN", "BHE", "BHZ"};

    SeismoQuantity quantity = SeismoQuantity::VELOCITY;
    double sample_rate_hz = 100.0;

    // Placement
    SeismoCoordinateType coord_type = SeismoCoordinateType::MODEL_XYZ;
    GridIndexMode grid_mode = GridIndexMode::CELL_CENTER;

    // MODEL_XYZ
    double x = 0.0, y = 0.0, z = 0.0;

    // GRID_INDEX
    int i = 0, j = 0, k = 0;

    // GEOGRAPHIC
    double lon = 0.0, lat = 0.0, elev = 0.0;

    // Instrument performance
    InstrumentPerformance instrument;
};

struct SeismometerOutputConfig {
    bool enabled = false;
    std::string output_dir = "output/seismometers";
    bool write_sac = true;
    bool write_mseed = true;

    // Absolute reference time for data formats (ISO-8601 UTC, e.g. 1992-09-23T15:00:00Z)
    // Simulation time t (seconds) is added to this base.
    std::string start_time_utc = "1970-01-01T00:00:00Z";
};

class SeismometerNetwork {
public:
    explicit SeismometerNetwork(MPI_Comm comm);
    ~SeismometerNetwork();

    SeismometerNetwork(const SeismometerNetwork&) = delete;
    SeismometerNetwork& operator=(const SeismometerNetwork&) = delete;

    void setOutputConfig(const SeismometerOutputConfig& cfg) { out_cfg_ = cfg; }
    void setStations(std::vector<SeismometerSpec> stations);

    // Must be called after DM + fields exist (setupFields), before sampling.
    PetscErrorCode initialize(DM dm,
                             const GridConfig& grid_cfg,
                             const CoordinateSystemManager* coord_mgr);

    // Called from TS monitor at each step; records if sampling cadence permits.
    PetscErrorCode sample(double t, Vec U);

    // Call once after simulation to write files.
    PetscErrorCode finalizeAndWrite() const;

private:
    struct StationRuntime {
        SeismometerSpec spec;

        // Final model coordinates (meters)
        double xm = 0.0, ym = 0.0, zm = 0.0;

        // Optional geographic (degrees, meters)
        bool have_geo = false;
        double stla = 0.0, stlo = 0.0, stel = 0.0;

        // Raw sampled time series (simulation times, seconds)
        std::vector<double> t;

        // Raw sampled displacement (meters)
        std::vector<std::array<double, 3>> disp;
    };

    // Helpers
    static bool startsWith(const std::string& s, const std::string& prefix);
    static std::string trim(const std::string& s);
    static std::string upper(std::string s);
    static std::vector<std::string> split(const std::string& s, char delim);

    static bool parseIsoUtc(const std::string& iso, int& year, int& month, int& day,
                            int& hour, int& min, int& sec);
    static bool ymdToYday(int year, int month, int day, int& yday);
    static void addSecondsToYdayTime(int& year, int& yday, int& hour, int& min, int& sec,
                                     int& frac10k, double add_seconds);

    static void applySimpleHP(const std::vector<float>& x, double fs, double fc, std::vector<float>& y);
    static void applySimpleLP(const std::vector<float>& x, double fs, double fc, std::vector<float>& y);

    static void writeSAC(const std::string& filename,
                         const std::string& net, const std::string& sta,
                         const std::string& loc, const std::string& chan,
                         const StationRuntime& st,
                         const std::vector<float>& data,
                         double delta,
                         const SeismometerOutputConfig& cfg);

    static void writeMiniSEED(const std::string& filename,
                              const std::string& net, const std::string& sta,
                              const std::string& loc, const std::string& chan,
                              const StationRuntime& st,
                              const std::vector<float>& data,
                              double sample_rate_hz,
                              const SeismometerOutputConfig& cfg);

    static std::vector<float> deriveVelocity(const std::vector<double>& t,
                                             const std::vector<std::array<double, 3>>& disp,
                                             int comp);
    static std::vector<float> deriveAcceleration(const std::vector<double>& t,
                                                 const std::vector<std::array<double, 3>>& disp,
                                                 int comp);

    static void applyInstrument(const InstrumentPerformance& inst,
                                double sample_rate_hz,
                                std::vector<float>& data);

    static void quantizeToCounts(const InstrumentPerformance& inst,
                                 std::vector<float>& data);

    // Coordinate conversion for lon/lat when PROJ isn't available.
    static bool geoToLocalMetersFallback(const GridConfig& grid_cfg,
                                         const SeismometerSpec& spec,
                                         double& x_m, double& y_m, double& z_m);

private:
    MPI_Comm comm_;
    int rank_ = 0;

    SeismometerOutputConfig out_cfg_;
    std::vector<StationRuntime> stations_;

    // PETSc interpolation setup (displacement field subDM)
    DM dm_ = nullptr;
    DM disp_dm_ = nullptr;
    IS disp_is_ = nullptr;
    DMInterpolationInfo interp_ = nullptr;
    Vec interp_result_ = nullptr;

    // Sampling control
    double last_sample_time_ = -1.0e300;
    double min_sample_dt_ = 0.0;
    bool initialized_ = false;
};

} // namespace FSRM

#endif // SEISMOMETER_NETWORK_HPP

