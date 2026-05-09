/**
 * @file TabulatedDataReader.hpp
 * @brief Pass-9 (axis 1) tabulated data reader for the EOS and opacity
 *        patches. Loads the binary HDF5 layout described in
 *        docs/EXPLOSION_IMPACT_PHYSICS.md "Pass-9 tabulated data patches"
 *        and tools/tabulated_data/README.md and exposes log-space
 *        bilinear interpolation with explicit out-of-range handling.
 *
 * The on-disk layout is shared between EOS and opacity tables:
 *   /metadata/medium             string  ("GRANITE"|"SALT"|"TUFF"|"ALLUVIUM")
 *   /metadata/quantity           string  ("EOS_PRESSURE"|"EOS_SOUND_SPEED"|
 *                                          "OPACITY_ROSSELAND"|"OPACITY_PLANCK")
 *   /metadata/rho_axis_kg_per_m3 1D dataset (log-spaced, ascending)
 *   /metadata/e_axis_J_per_kg    1D dataset (EOS only)
 *   /metadata/T_axis_K           1D dataset (opacity only)
 *   /metadata/rho_min, rho_max   scalars
 *   /metadata/e_min, e_max       scalars (or T_min, T_max)
 *   /metadata/source_citation    string (one-line literature reference,
 *                                        REQUIRED; readers refuse a table
 *                                        without it).
 *   /metadata/generation_date    ISO 8601 string
 *   /metadata/generation_tool    string ("scripts/generate_aneos_table.py"
 *                                        with version, etc.)
 *   /data                        2D dataset (n_rho x n_other)
 *
 * Out-of-range queries do not silently extrapolate. The reader returns
 * a sentinel (NaN) and the caller (the EOS / opacity dispatch) is
 * expected to fall back to the underlying analytic model. A one-time
 * warning per (medium, quantity, axis) pair is logged to stderr.
 *
 * State convention: SI units throughout. Density rho [kg/m^3], specific
 * internal energy e [J/kg], temperature T [K], pressure p [Pa], opacity
 * kappa [m^2/kg].
 *
 * References.
 *  - Thompson, S. L. and Lauson, H. S. (1972), "Improvements in the
 *    Chart D radiation-hydrodynamic CDC 6600 computer code describing
 *    a thermonuclear weapon", Sandia SC-RR-71-0714 (the canonical
 *    description of the ANEOS-style EOS that the table generator
 *    re-implements).
 *  - Melosh, H. J. (1989), "Impact Cratering: A Geologic Process",
 *    Oxford University Press, sec A2.2 (Tillotson-form EOS, parameter
 *    sets for granite, basalt, NaCl).
 *  - Zel'dovich, Y. B. and Raizer, Y. P. (1967), "Physics of Shock
 *    Waves and High-Temperature Hydrodynamic Phenomena", vol I, ch X
 *    sec 7-10 (Kramers' free-free opacity, plasma-regime corrections;
 *    target literature for the opacity table generator).
 *  - Mihalas, D. and Mihalas, B. W. (1984), "Foundations of Radiation
 *    Hydrodynamics", Oxford University Press, sec 82-83 (Rosseland
 *    and Planck mean averaging).
 */

#ifndef FSRM_IO_TABULATED_DATA_READER_HPP
#define FSRM_IO_TABULATED_DATA_READER_HPP

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace FSRM {
namespace io {

/// Sentinel returned by the reader on out-of-range queries. The caller
/// must check std::isnan(result) and fall back to the analytic model.
constexpr double TABULATED_DATA_OOR_SENTINEL =
    std::numeric_limits<double>::quiet_NaN();

/// Quantity stored in the table. The reader does not enforce units on
/// the stored values; the consumer (EOS / opacity dispatch) is
/// responsible for reading the metadata and confirming the quantity it
/// expects.
enum class TabulatedQuantity
{
    EOS_PRESSURE,        ///< Pressure p(rho, e) in Pa.
    EOS_SOUND_SPEED,     ///< Sound speed c(rho, e) in m/s.
    OPACITY_ROSSELAND,   ///< Rosseland mean opacity kappa_R(rho, T) in m^2/kg.
    OPACITY_PLANCK,      ///< Planck mean opacity kappa_P(rho, T) in m^2/kg.
    UNKNOWN
};

/// Axis ordering of the 2D data dataset. Different generators may
/// write rho-major or rho-minor; the reader detects from the dataset
/// shape and the axis lengths in /metadata.
enum class AxisOrdering
{
    RHO_MAJOR,   ///< data[i_rho, j_other], n_rho rows by n_other cols.
    RHO_MINOR,   ///< data[i_other, j_rho], n_other rows by n_rho cols.
};

struct TabulatedTableMetadata
{
    std::string medium;            ///< "GRANITE" | "SALT" | "TUFF" | "ALLUVIUM"
    std::string quantity_str;      ///< Verbatim quantity string from /metadata.
    TabulatedQuantity quantity = TabulatedQuantity::UNKNOWN;
    std::string source_citation;   ///< One-line literature reference. Required.
    std::string generation_date;   ///< ISO 8601.
    std::string generation_tool;   ///< Tool name + version.
    double rho_min = 0.0;
    double rho_max = 0.0;
    double other_min = 0.0;        ///< e_min for EOS, T_min for opacity.
    double other_max = 0.0;
};

/**
 * @brief HDF5 tabulated-data reader. Bilinear interpolation in
 *        log(rho), log(other) space. Constructed against a file path,
 *        loads metadata + axes + the 2D data dataset, then serves
 *        evaluate(rho, other) calls in O(log n) per axis (binary search)
 *        plus O(1) interpolation.
 *
 * The class is intentionally simple and stateful only after a
 * successful load(). Failed loads leave the reader in the
 * "not-loaded" state; isLoaded() returns false; evaluate() returns
 * the OOR sentinel.
 */
class TabulatedDataReader
{
public:
    TabulatedDataReader() = default;

    /// Load a table from an HDF5 file. Returns true on success.
    /// On failure (file missing, schema mismatch, missing
    /// source_citation), populates error_message and returns false.
    /// Idempotent: a second successful load replaces the first.
    bool load(const std::string& path, std::string& error_message);

    /// Whether a successful load has been performed.
    bool isLoaded() const { return loaded_; }

    /// Bilinear interpolation in log space. For an EOS table the
    /// "other" argument is e [J/kg]; for an opacity table it is
    /// T [K]. Out-of-range returns TABULATED_DATA_OOR_SENTINEL and
    /// logs a one-time warning identifying the (medium, quantity,
    /// axis) and the excursion. Below floor, the reader does not
    /// extrapolate either; the caller falls back.
    ///
    /// The optional cell_index parameter is forwarded into the warning
    /// message to help users localise the diagnostic; -1 means
    /// unspecified.
    double evaluate(double rho, double other, int cell_index = -1) const;

    /// Bilinear interpolation derivatives in log space. Returns
    /// dp_drho and dp_dother (or dkappa_dT etc.) at (rho, other).
    /// On out-of-range both derivatives are set to 0.0 and the
    /// function returns false.
    bool evaluateDerivatives(double rho, double other,
                             double& dval_drho, double& dval_dother,
                             int cell_index = -1) const;

    const TabulatedTableMetadata& metadata() const { return meta_; }

    /// Diagnostic accessors for the unit / round-trip tests.
    int numRhoPoints() const { return static_cast<int>(rho_axis_.size()); }
    int numOtherPoints() const { return static_cast<int>(other_axis_.size()); }
    AxisOrdering axisOrdering() const { return ordering_; }

    /// Lookup helpers for tests and diagnostics. Indexing is
    /// (i_rho, j_other) regardless of the on-disk axis ordering.
    double dataAt(int i_rho, int j_other) const;

    /// Internal log-space helper exposed for tests.
    static double logSafe(double x)
    {
        return std::log(x > 1.0e-300 ? x : 1.0e-300);
    }

private:
    bool loaded_ = false;
    TabulatedTableMetadata meta_;
    AxisOrdering ordering_ = AxisOrdering::RHO_MAJOR;
    std::vector<double> rho_axis_;
    std::vector<double> log_rho_axis_;
    std::vector<double> other_axis_;
    std::vector<double> log_other_axis_;
    std::vector<double> data_;       ///< Row-major in (rho, other) ordering.
    int n_rho_ = 0;
    int n_other_ = 0;

    mutable bool oor_warned_low_rho_ = false;
    mutable bool oor_warned_high_rho_ = false;
    mutable bool oor_warned_low_other_ = false;
    mutable bool oor_warned_high_other_ = false;

    /// Log a one-time stderr warning. Mutates the matching
    /// oor_warned_* member.
    void warnOnce(const char* axis_label, bool& flag,
                  double queried, double bound,
                  int cell_index) const;

    /// Map (i_rho, j_other) to the 1D index in data_, given the
    /// stored axis ordering.
    int index2(int i_rho, int j_other) const;
};

} // namespace io
} // namespace FSRM

#endif // FSRM_IO_TABULATED_DATA_READER_HPP
