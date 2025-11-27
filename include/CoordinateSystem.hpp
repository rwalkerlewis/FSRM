#ifndef COORDINATE_SYSTEM_HPP
#define COORDINATE_SYSTEM_HPP

#include <string>
#include <vector>
#include <memory>
#include <array>
#include <cmath>

// Forward declaration for PROJ types (avoid including proj.h in header)
struct PJ;
struct PJ_CONTEXT;

namespace FSRM {

/**
 * @brief 3D point with optional elevation
 */
struct GeoPoint {
    double x;       ///< X coordinate (easting/longitude)
    double y;       ///< Y coordinate (northing/latitude)
    double z;       ///< Z coordinate (elevation/depth)
    
    GeoPoint() : x(0), y(0), z(0) {}
    GeoPoint(double xx, double yy, double zz = 0) : x(xx), y(yy), z(zz) {}
};

/**
 * @brief Bounding box in geographic or projected coordinates
 */
struct GeoBounds {
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
    
    GeoBounds() : x_min(0), x_max(0), y_min(0), y_max(0), z_min(0), z_max(0) {}
    
    double getLx() const { return x_max - x_min; }
    double getLy() const { return y_max - y_min; }
    double getLz() const { return z_max - z_min; }
    
    GeoPoint getCenter() const {
        return GeoPoint((x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2);
    }
};

/**
 * @brief Coordinate Reference System definition
 */
struct CRSDefinition {
    std::string epsg_code;      ///< EPSG code (e.g., "EPSG:4326")
    std::string proj_string;    ///< PROJ string (alternative to EPSG)
    std::string wkt;            ///< WKT definition (alternative)
    std::string name;           ///< Human-readable name
    
    bool isGeographic() const;  ///< True if geographic (lat/lon)
    bool isProjected() const;   ///< True if projected (meters/feet)
    
    CRSDefinition() = default;
    CRSDefinition(const std::string& epsg) : epsg_code(epsg) {}
};

/**
 * @brief Coordinate transformation between CRS using PROJ library
 * 
 * This class provides coordinate transformations between different coordinate
 * reference systems using the PROJ library. It supports EPSG codes, PROJ strings,
 * and WKT definitions.
 * 
 * Usage:
 * @code
 * // Transform from WGS84 to UTM Zone 10N
 * CoordinateTransformer transformer;
 * transformer.setSourceCRS("EPSG:4326");
 * transformer.setTargetCRS("EPSG:32610");
 * 
 * GeoPoint pt_wgs84(-122.4194, 37.7749);  // San Francisco
 * GeoPoint pt_utm = transformer.transform(pt_wgs84);
 * @endcode
 */
class CoordinateTransformer {
public:
    CoordinateTransformer();
    ~CoordinateTransformer();
    
    // Disable copy (PROJ handles are not copyable)
    CoordinateTransformer(const CoordinateTransformer&) = delete;
    CoordinateTransformer& operator=(const CoordinateTransformer&) = delete;
    
    // Move semantics
    CoordinateTransformer(CoordinateTransformer&& other) noexcept;
    CoordinateTransformer& operator=(CoordinateTransformer&& other) noexcept;
    
    // =========================================================================
    // CRS Configuration
    // =========================================================================
    
    /**
     * @brief Set source CRS by EPSG code
     * @param epsg EPSG code (e.g., "EPSG:4326" or "4326")
     * @return true if CRS was set successfully
     */
    bool setSourceCRS(const std::string& epsg);
    
    /**
     * @brief Set target CRS by EPSG code
     * @param epsg EPSG code (e.g., "EPSG:32610" or "32610")
     * @return true if CRS was set successfully
     */
    bool setTargetCRS(const std::string& epsg);
    
    /**
     * @brief Set source CRS by PROJ string
     * @param proj_string PROJ definition string
     * @return true if CRS was set successfully
     */
    bool setSourceCRSFromProj(const std::string& proj_string);
    
    /**
     * @brief Set target CRS by PROJ string
     * @param proj_string PROJ definition string
     * @return true if CRS was set successfully
     */
    bool setTargetCRSFromProj(const std::string& proj_string);
    
    /**
     * @brief Set source CRS by WKT string
     * @param wkt WKT definition string
     * @return true if CRS was set successfully
     */
    bool setSourceCRSFromWKT(const std::string& wkt);
    
    /**
     * @brief Set target CRS by WKT string
     * @param wkt WKT definition string
     * @return true if CRS was set successfully
     */
    bool setTargetCRSFromWKT(const std::string& wkt);
    
    // =========================================================================
    // Transformation Methods
    // =========================================================================
    
    /**
     * @brief Initialize transformation pipeline
     * @return true if transformation is ready
     */
    bool initialize();
    
    /**
     * @brief Transform a single point
     * @param point Input point in source CRS
     * @return Transformed point in target CRS
     */
    GeoPoint transform(const GeoPoint& point) const;
    
    /**
     * @brief Transform multiple points (in-place)
     * @param points Vector of points to transform
     * @return true if all points transformed successfully
     */
    bool transform(std::vector<GeoPoint>& points) const;
    
    /**
     * @brief Transform coordinate arrays
     * @param x X coordinates (modified in place)
     * @param y Y coordinates (modified in place)
     * @param z Z coordinates (modified in place, can be nullptr)
     * @param n Number of points
     * @return true if transformation succeeded
     */
    bool transform(double* x, double* y, double* z, size_t n) const;
    
    /**
     * @brief Inverse transform (target to source)
     * @param point Input point in target CRS
     * @return Transformed point in source CRS
     */
    GeoPoint inverseTransform(const GeoPoint& point) const;
    
    // =========================================================================
    // Query Methods
    // =========================================================================
    
    /**
     * @brief Check if transformation is valid
     */
    bool isValid() const { return is_valid_; }
    
    /**
     * @brief Get source CRS definition
     */
    const CRSDefinition& getSourceCRS() const { return source_crs_; }
    
    /**
     * @brief Get target CRS definition
     */
    const CRSDefinition& getTargetCRS() const { return target_crs_; }
    
    /**
     * @brief Get last error message
     */
    const std::string& getLastError() const { return last_error_; }
    
    /**
     * @brief Check if PROJ library is available
     * @return true if PROJ is available
     */
    static bool isProjAvailable();
    
    /**
     * @brief Get PROJ library version
     */
    static std::string getProjVersion();
    
private:
    CRSDefinition source_crs_;
    CRSDefinition target_crs_;
    
    // PROJ handles (opaque pointers)
    PJ_CONTEXT* ctx_;
    PJ* transform_;
    
    bool is_valid_;
    mutable std::string last_error_;
    
    void cleanup();
    std::string normalizeEPSG(const std::string& epsg) const;
};

/**
 * @brief Common CRS definitions
 */
namespace CRS {
    // Geographic CRS
    const std::string WGS84 = "EPSG:4326";           ///< WGS 84 (GPS standard)
    const std::string NAD83 = "EPSG:4269";           ///< NAD83 (North America)
    const std::string NAD27 = "EPSG:4267";           ///< NAD27 (legacy North America)
    const std::string ETRS89 = "EPSG:4258";          ///< ETRS89 (Europe)
    
    // UTM Zones (example - WGS84 based)
    const std::string UTM10N = "EPSG:32610";         ///< UTM Zone 10N (California)
    const std::string UTM11N = "EPSG:32611";         ///< UTM Zone 11N (Nevada)
    const std::string UTM12N = "EPSG:32612";         ///< UTM Zone 12N (Arizona)
    const std::string UTM13N = "EPSG:32613";         ///< UTM Zone 13N (Colorado)
    const std::string UTM14N = "EPSG:32614";         ///< UTM Zone 14N (Texas)
    const std::string UTM15N = "EPSG:32615";         ///< UTM Zone 15N (Louisiana)
    
    // State Plane (NAD83)
    const std::string SPCS_TX_C = "EPSG:2277";       ///< Texas Central
    const std::string SPCS_CA_3 = "EPSG:2227";       ///< California Zone 3
    
    // Web Mercator (common for web maps)
    const std::string WEB_MERCATOR = "EPSG:3857";    ///< Web Mercator
    
    /**
     * @brief Get UTM zone EPSG code for WGS84
     * @param zone Zone number (1-60)
     * @param north True for northern hemisphere
     * @return EPSG code string
     */
    inline std::string getUTMZone(int zone, bool north = true) {
        int base = north ? 32600 : 32700;
        return "EPSG:" + std::to_string(base + zone);
    }
    
    /**
     * @brief Calculate UTM zone from longitude
     * @param longitude Longitude in degrees (-180 to 180)
     * @return UTM zone number (1-60)
     */
    inline int calculateUTMZone(double longitude) {
        return static_cast<int>((longitude + 180) / 6) + 1;
    }
}

/**
 * @brief Coordinate system manager for a simulation
 * 
 * Manages coordinate transformations for a simulation, converting between
 * the user's input CRS and the internal model CRS (typically local meters).
 */
class CoordinateSystemManager {
public:
    CoordinateSystemManager();
    ~CoordinateSystemManager() = default;
    
    // =========================================================================
    // Configuration
    // =========================================================================
    
    /**
     * @brief Set the input coordinate system
     * @param epsg EPSG code for input coordinates
     */
    bool setInputCRS(const std::string& epsg);
    
    /**
     * @brief Set the model coordinate system (internal)
     * @param epsg EPSG code for model coordinates
     */
    bool setModelCRS(const std::string& epsg);
    
    /**
     * @brief Set local origin for model coordinates
     * @param origin Origin point in input CRS
     */
    void setLocalOrigin(const GeoPoint& origin);
    
    /**
     * @brief Use auto-detected local CRS (UTM based on origin)
     */
    bool useAutoLocalCRS();
    
    // =========================================================================
    // Coordinate Transformation
    // =========================================================================
    
    /**
     * @brief Transform point from input CRS to model CRS
     */
    GeoPoint toModelCoords(const GeoPoint& input) const;
    
    /**
     * @brief Transform point from model CRS to input CRS
     */
    GeoPoint toInputCoords(const GeoPoint& model) const;
    
    /**
     * @brief Transform points array from input to model CRS
     */
    void toModelCoords(std::vector<GeoPoint>& points) const;
    
    /**
     * @brief Transform points array from model to input CRS
     */
    void toInputCoords(std::vector<GeoPoint>& points) const;
    
    /**
     * @brief Apply local origin offset (after CRS transformation)
     */
    GeoPoint applyLocalOrigin(const GeoPoint& pt) const;
    
    /**
     * @brief Remove local origin offset (before inverse transformation)
     */
    GeoPoint removeLocalOrigin(const GeoPoint& pt) const;
    
    // =========================================================================
    // Query Methods
    // =========================================================================
    
    /**
     * @brief Check if coordinate transformation is configured
     */
    bool isConfigured() const;
    
    /**
     * @brief Get the input CRS
     */
    const CRSDefinition& getInputCRS() const { return input_crs_; }
    
    /**
     * @brief Get the model CRS
     */
    const CRSDefinition& getModelCRS() const { return model_crs_; }
    
    /**
     * @brief Get the local origin
     */
    const GeoPoint& getLocalOrigin() const { return local_origin_; }
    
    /**
     * @brief Check if using local origin
     */
    bool hasLocalOrigin() const { return has_local_origin_; }
    
private:
    CRSDefinition input_crs_;
    CRSDefinition model_crs_;
    GeoPoint local_origin_;
    bool has_local_origin_;
    
    std::unique_ptr<CoordinateTransformer> forward_transform_;
    std::unique_ptr<CoordinateTransformer> inverse_transform_;
};

/**
 * @brief Geodetic calculation utilities
 */
namespace Geodetic {
    /**
     * @brief Calculate distance between two points on WGS84 ellipsoid (Vincenty's formula)
     * @param lat1, lon1 First point in degrees
     * @param lat2, lon2 Second point in degrees
     * @return Distance in meters
     */
    double vincentyDistance(double lat1, double lon1, double lat2, double lon2);
    
    /**
     * @brief Calculate distance using Haversine formula (spherical approximation)
     * @param lat1, lon1 First point in degrees
     * @param lat2, lon2 Second point in degrees
     * @return Distance in meters
     */
    double haversineDistance(double lat1, double lon1, double lat2, double lon2);
    
    /**
     * @brief Convert degrees to radians
     */
    inline double deg2rad(double degrees) {
        return degrees * M_PI / 180.0;
    }
    
    /**
     * @brief Convert radians to degrees
     */
    inline double rad2deg(double radians) {
        return radians * 180.0 / M_PI;
    }
    
    /**
     * @brief Calculate destination point given start, bearing, and distance
     * @param lat, lon Starting point in degrees
     * @param bearing Bearing in degrees (0 = north, 90 = east)
     * @param distance Distance in meters
     * @return Destination point
     */
    GeoPoint destinationPoint(double lat, double lon, double bearing, double distance);
    
    /**
     * @brief Calculate bearing from point 1 to point 2
     * @return Bearing in degrees
     */
    double calculateBearing(double lat1, double lon1, double lat2, double lon2);
}

} // namespace FSRM

#endif // COORDINATE_SYSTEM_HPP
