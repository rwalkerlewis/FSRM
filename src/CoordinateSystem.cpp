#include "CoordinateSystem.hpp"
#include <cmath>
#include <sstream>
#include <iostream>

// Conditionally include PROJ if available
#ifdef HAVE_PROJ
#include <proj.h>
#define PROJ_AVAILABLE 1
#else
#define PROJ_AVAILABLE 0
#endif

namespace FSRM {

// ============================================================================
// CRSDefinition Implementation
// ============================================================================

bool CRSDefinition::isGeographic() const {
    // Common geographic CRS codes
    if (epsg_code == "EPSG:4326" || epsg_code == "4326") return true;
    if (epsg_code == "EPSG:4269" || epsg_code == "4269") return true;
    if (epsg_code == "EPSG:4267" || epsg_code == "4267") return true;
    if (epsg_code == "EPSG:4258" || epsg_code == "4258") return true;
    
    // Check PROJ string for geographic keywords
    if (proj_string.find("+proj=longlat") != std::string::npos) return true;
    if (proj_string.find("+proj=latlong") != std::string::npos) return true;
    
    return false;
}

bool CRSDefinition::isProjected() const {
    return !isGeographic() && !epsg_code.empty();
}

// ============================================================================
// CoordinateTransformer Implementation
// ============================================================================

CoordinateTransformer::CoordinateTransformer() 
    : ctx_(nullptr), transform_(nullptr), is_valid_(false) {
#if PROJ_AVAILABLE
    ctx_ = proj_context_create();
#endif
}

CoordinateTransformer::~CoordinateTransformer() {
    cleanup();
}

CoordinateTransformer::CoordinateTransformer(CoordinateTransformer&& other) noexcept
    : source_crs_(std::move(other.source_crs_)),
      target_crs_(std::move(other.target_crs_)),
      ctx_(other.ctx_),
      transform_(other.transform_),
      is_valid_(other.is_valid_),
      last_error_(std::move(other.last_error_)) {
    other.ctx_ = nullptr;
    other.transform_ = nullptr;
    other.is_valid_ = false;
}

CoordinateTransformer& CoordinateTransformer::operator=(CoordinateTransformer&& other) noexcept {
    if (this != &other) {
        cleanup();
        source_crs_ = std::move(other.source_crs_);
        target_crs_ = std::move(other.target_crs_);
        ctx_ = other.ctx_;
        transform_ = other.transform_;
        is_valid_ = other.is_valid_;
        last_error_ = std::move(other.last_error_);
        other.ctx_ = nullptr;
        other.transform_ = nullptr;
        other.is_valid_ = false;
    }
    return *this;
}

void CoordinateTransformer::cleanup() {
#if PROJ_AVAILABLE
    if (transform_) {
        proj_destroy(transform_);
        transform_ = nullptr;
    }
    if (ctx_) {
        proj_context_destroy(ctx_);
        ctx_ = nullptr;
    }
#endif
    is_valid_ = false;
}

std::string CoordinateTransformer::normalizeEPSG(const std::string& epsg) const {
    // If already has EPSG: prefix, return as is
    if (epsg.substr(0, 5) == "EPSG:") {
        return epsg;
    }
    // If it's just a number, add prefix
    try {
        std::stoi(epsg);
        return "EPSG:" + epsg;
    } catch (...) {
        // Not a number, return as is (might be a PROJ string)
        return epsg;
    }
}

bool CoordinateTransformer::setSourceCRS(const std::string& epsg) {
    source_crs_.epsg_code = normalizeEPSG(epsg);
    is_valid_ = false;
    return true;
}

bool CoordinateTransformer::setTargetCRS(const std::string& epsg) {
    target_crs_.epsg_code = normalizeEPSG(epsg);
    is_valid_ = false;
    return true;
}

bool CoordinateTransformer::setSourceCRSFromProj(const std::string& proj_string) {
    source_crs_.proj_string = proj_string;
    source_crs_.epsg_code.clear();
    is_valid_ = false;
    return true;
}

bool CoordinateTransformer::setTargetCRSFromProj(const std::string& proj_string) {
    target_crs_.proj_string = proj_string;
    target_crs_.epsg_code.clear();
    is_valid_ = false;
    return true;
}

bool CoordinateTransformer::setSourceCRSFromWKT(const std::string& wkt) {
    source_crs_.wkt = wkt;
    source_crs_.epsg_code.clear();
    source_crs_.proj_string.clear();
    is_valid_ = false;
    return true;
}

bool CoordinateTransformer::setTargetCRSFromWKT(const std::string& wkt) {
    target_crs_.wkt = wkt;
    target_crs_.epsg_code.clear();
    target_crs_.proj_string.clear();
    is_valid_ = false;
    return true;
}

bool CoordinateTransformer::initialize() {
#if PROJ_AVAILABLE
    if (transform_) {
        proj_destroy(transform_);
        transform_ = nullptr;
    }
    
    // Determine source CRS string
    std::string src_str;
    if (!source_crs_.epsg_code.empty()) {
        src_str = source_crs_.epsg_code;
    } else if (!source_crs_.proj_string.empty()) {
        src_str = source_crs_.proj_string;
    } else if (!source_crs_.wkt.empty()) {
        src_str = source_crs_.wkt;
    } else {
        last_error_ = "Source CRS not specified";
        return false;
    }
    
    // Determine target CRS string
    std::string tgt_str;
    if (!target_crs_.epsg_code.empty()) {
        tgt_str = target_crs_.epsg_code;
    } else if (!target_crs_.proj_string.empty()) {
        tgt_str = target_crs_.proj_string;
    } else if (!target_crs_.wkt.empty()) {
        tgt_str = target_crs_.wkt;
    } else {
        last_error_ = "Target CRS not specified";
        return false;
    }
    
    // Create transformation
    transform_ = proj_create_crs_to_crs(ctx_, src_str.c_str(), tgt_str.c_str(), nullptr);
    
    if (!transform_) {
        int err = proj_context_errno(ctx_);
        last_error_ = std::string("Failed to create transformation: ") + 
                     proj_errno_string(err);
        return false;
    }
    
    // Normalize for longitude/latitude ordering
    PJ* norm = proj_normalize_for_visualization(ctx_, transform_);
    if (norm) {
        proj_destroy(transform_);
        transform_ = norm;
    }
    
    is_valid_ = true;
    return true;
#else
    last_error_ = "PROJ library not available. Coordinate transformations disabled.";
    
    // Check if source and target are the same (identity transform)
    if (source_crs_.epsg_code == target_crs_.epsg_code &&
        source_crs_.proj_string == target_crs_.proj_string) {
        is_valid_ = true;
        return true;
    }
    
    return false;
#endif
}

GeoPoint CoordinateTransformer::transform(const GeoPoint& point) const {
#if PROJ_AVAILABLE
    if (!is_valid_ || !transform_) {
        last_error_ = "Transformation not initialized";
        return point;
    }
    
    PJ_COORD in, out;
    in.xyz.x = point.x;
    in.xyz.y = point.y;
    in.xyz.z = point.z;
    
    out = proj_trans(transform_, PJ_FWD, in);
    
    if (proj_errno(transform_)) {
        last_error_ = proj_errno_string(proj_errno(transform_));
        return point;
    }
    
    return GeoPoint(out.xyz.x, out.xyz.y, out.xyz.z);
#else
    // Identity transform if PROJ not available
    return point;
#endif
}

bool CoordinateTransformer::transform(std::vector<GeoPoint>& points) const {
    bool success = true;
    for (auto& pt : points) {
        GeoPoint result = transform(pt);
        if (!last_error_.empty()) {
            success = false;
        }
        pt = result;
    }
    return success;
}

bool CoordinateTransformer::transform(double* x, double* y, double* z, size_t n) const {
#if PROJ_AVAILABLE
    if (!is_valid_ || !transform_) {
        last_error_ = "Transformation not initialized";
        return false;
    }
    
    // Transform in batches for efficiency
    int result = proj_trans_generic(
        transform_, PJ_FWD,
        x, sizeof(double), n,
        y, sizeof(double), n,
        z, z ? sizeof(double) : 0, z ? n : 0,
        nullptr, 0, 0
    );
    
    if (result != static_cast<int>(n)) {
        last_error_ = "Not all points transformed successfully";
        return false;
    }
    
    return true;
#else
    // Identity transform if PROJ not available
    (void)x; (void)y; (void)z; (void)n;
    return true;
#endif
}

GeoPoint CoordinateTransformer::inverseTransform(const GeoPoint& point) const {
#if PROJ_AVAILABLE
    if (!is_valid_ || !transform_) {
        last_error_ = "Transformation not initialized";
        return point;
    }
    
    PJ_COORD in, out;
    in.xyz.x = point.x;
    in.xyz.y = point.y;
    in.xyz.z = point.z;
    
    out = proj_trans(transform_, PJ_INV, in);
    
    if (proj_errno(transform_)) {
        last_error_ = proj_errno_string(proj_errno(transform_));
        return point;
    }
    
    return GeoPoint(out.xyz.x, out.xyz.y, out.xyz.z);
#else
    return point;
#endif
}

bool CoordinateTransformer::isProjAvailable() {
#if PROJ_AVAILABLE
    return true;
#else
    return false;
#endif
}

std::string CoordinateTransformer::getProjVersion() {
#if PROJ_AVAILABLE
    return std::string(proj_info().version);
#else
    return "Not available";
#endif
}

// ============================================================================
// CoordinateSystemManager Implementation
// ============================================================================

CoordinateSystemManager::CoordinateSystemManager()
    : has_local_origin_(false) {
}

bool CoordinateSystemManager::setInputCRS(const std::string& epsg) {
    input_crs_.epsg_code = epsg;
    
    // Recreate transformers
    forward_transform_ = std::make_unique<CoordinateTransformer>();
    inverse_transform_ = std::make_unique<CoordinateTransformer>();
    
    forward_transform_->setSourceCRS(epsg);
    inverse_transform_->setTargetCRS(epsg);
    
    return true;
}

bool CoordinateSystemManager::setModelCRS(const std::string& epsg) {
    model_crs_.epsg_code = epsg;
    
    if (forward_transform_) {
        forward_transform_->setTargetCRS(epsg);
    }
    if (inverse_transform_) {
        inverse_transform_->setSourceCRS(epsg);
    }
    
    return true;
}

void CoordinateSystemManager::setLocalOrigin(const GeoPoint& origin) {
    local_origin_ = origin;
    has_local_origin_ = true;
}

bool CoordinateSystemManager::useAutoLocalCRS() {
    if (input_crs_.epsg_code.empty()) {
        return false;
    }
    
    // If input is geographic, auto-select UTM zone based on origin
    if (input_crs_.isGeographic() && has_local_origin_) {
        int zone = CRS::calculateUTMZone(local_origin_.x);
        bool north = (local_origin_.y >= 0);
        
        model_crs_.epsg_code = CRS::getUTMZone(zone, north);
        
        if (forward_transform_) {
            forward_transform_->setTargetCRS(model_crs_.epsg_code);
        }
        if (inverse_transform_) {
            inverse_transform_->setSourceCRS(model_crs_.epsg_code);
        }
        
        return true;
    }
    
    return false;
}

GeoPoint CoordinateSystemManager::toModelCoords(const GeoPoint& input) const {
    GeoPoint result = input;
    
    // Apply CRS transformation
    if (forward_transform_ && forward_transform_->isValid()) {
        result = forward_transform_->transform(result);
    } else if (forward_transform_ && !input_crs_.epsg_code.empty() && !model_crs_.epsg_code.empty()) {
        // Try to initialize
        const_cast<CoordinateTransformer*>(forward_transform_.get())->initialize();
        if (forward_transform_->isValid()) {
            result = forward_transform_->transform(result);
        }
    }
    
    // Apply local origin
    if (has_local_origin_) {
        result = applyLocalOrigin(result);
    }
    
    return result;
}

GeoPoint CoordinateSystemManager::toInputCoords(const GeoPoint& model) const {
    GeoPoint result = model;
    
    // Remove local origin
    if (has_local_origin_) {
        result = removeLocalOrigin(result);
    }
    
    // Apply inverse CRS transformation
    if (inverse_transform_ && inverse_transform_->isValid()) {
        result = inverse_transform_->transform(result);
    } else if (inverse_transform_ && !input_crs_.epsg_code.empty() && !model_crs_.epsg_code.empty()) {
        // Try to initialize
        const_cast<CoordinateTransformer*>(inverse_transform_.get())->initialize();
        if (inverse_transform_->isValid()) {
            result = inverse_transform_->transform(result);
        }
    }
    
    return result;
}

void CoordinateSystemManager::toModelCoords(std::vector<GeoPoint>& points) const {
    for (auto& pt : points) {
        pt = toModelCoords(pt);
    }
}

void CoordinateSystemManager::toInputCoords(std::vector<GeoPoint>& points) const {
    for (auto& pt : points) {
        pt = toInputCoords(pt);
    }
}

GeoPoint CoordinateSystemManager::applyLocalOrigin(const GeoPoint& pt) const {
    // Transform origin to model CRS first
    GeoPoint origin_model = local_origin_;
    if (forward_transform_ && forward_transform_->isValid()) {
        origin_model = forward_transform_->transform(local_origin_);
    }
    
    return GeoPoint(
        pt.x - origin_model.x,
        pt.y - origin_model.y,
        pt.z - origin_model.z
    );
}

GeoPoint CoordinateSystemManager::removeLocalOrigin(const GeoPoint& pt) const {
    // Transform origin to model CRS first
    GeoPoint origin_model = local_origin_;
    if (forward_transform_ && forward_transform_->isValid()) {
        origin_model = forward_transform_->transform(local_origin_);
    }
    
    return GeoPoint(
        pt.x + origin_model.x,
        pt.y + origin_model.y,
        pt.z + origin_model.z
    );
}

bool CoordinateSystemManager::isConfigured() const {
    return !input_crs_.epsg_code.empty() && !model_crs_.epsg_code.empty();
}

// ============================================================================
// Geodetic Utilities Implementation
// ============================================================================

namespace Geodetic {

// WGS84 ellipsoid constants
constexpr double WGS84_A = 6378137.0;          // Semi-major axis
constexpr double WGS84_F = 1.0 / 298.257223563; // Flattening
constexpr double WGS84_B = WGS84_A * (1 - WGS84_F); // Semi-minor axis
constexpr double EARTH_RADIUS = 6371000.0;      // Mean radius for spherical calculations

double vincentyDistance(double lat1, double lon1, double lat2, double lon2) {
    // Convert to radians
    double phi1 = deg2rad(lat1);
    double phi2 = deg2rad(lat2);
    double L = deg2rad(lon2 - lon1);
    
    double U1 = std::atan((1 - WGS84_F) * std::tan(phi1));
    double U2 = std::atan((1 - WGS84_F) * std::tan(phi2));
    
    double sinU1 = std::sin(U1), cosU1 = std::cos(U1);
    double sinU2 = std::sin(U2), cosU2 = std::cos(U2);
    
    double lambda = L;
    double lambda_prev;
    int iterations = 0;
    const int max_iterations = 100;
    const double tolerance = 1e-12;
    
    double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, cos2SigmaM;
    
    do {
        double sinLambda = std::sin(lambda);
        double cosLambda = std::cos(lambda);
        
        sinSigma = std::sqrt(
            std::pow(cosU2 * sinLambda, 2) +
            std::pow(cosU1 * sinU2 - sinU1 * cosU2 * cosLambda, 2)
        );
        
        if (sinSigma == 0) return 0;  // Coincident points
        
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma = std::atan2(sinSigma, cosSigma);
        
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1 - sinAlpha * sinAlpha;
        
        cos2SigmaM = (cosSqAlpha == 0) ? 0 : cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
        
        double C = WGS84_F / 16 * cosSqAlpha * (4 + WGS84_F * (4 - 3 * cosSqAlpha));
        
        lambda_prev = lambda;
        lambda = L + (1 - C) * WGS84_F * sinAlpha * (
            sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM))
        );
        
    } while (std::abs(lambda - lambda_prev) > tolerance && ++iterations < max_iterations);
    
    if (iterations >= max_iterations) {
        // Fall back to haversine if Vincenty doesn't converge
        return haversineDistance(lat1, lon1, lat2, lon2);
    }
    
    double uSq = cosSqAlpha * (WGS84_A * WGS84_A - WGS84_B * WGS84_B) / (WGS84_B * WGS84_B);
    double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    
    double deltaSigma = B * sinSigma * (
        cos2SigmaM + B / 4 * (
            cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
            B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)
        )
    );
    
    return WGS84_B * A * (sigma - deltaSigma);
}

double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
    double phi1 = deg2rad(lat1);
    double phi2 = deg2rad(lat2);
    double dPhi = deg2rad(lat2 - lat1);
    double dLambda = deg2rad(lon2 - lon1);
    
    double a = std::sin(dPhi / 2) * std::sin(dPhi / 2) +
               std::cos(phi1) * std::cos(phi2) *
               std::sin(dLambda / 2) * std::sin(dLambda / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    
    return EARTH_RADIUS * c;
}

GeoPoint destinationPoint(double lat, double lon, double bearing, double distance) {
    double phi1 = deg2rad(lat);
    double lambda1 = deg2rad(lon);
    double theta = deg2rad(bearing);
    double delta = distance / EARTH_RADIUS;
    
    double phi2 = std::asin(
        std::sin(phi1) * std::cos(delta) +
        std::cos(phi1) * std::sin(delta) * std::cos(theta)
    );
    
    double lambda2 = lambda1 + std::atan2(
        std::sin(theta) * std::sin(delta) * std::cos(phi1),
        std::cos(delta) - std::sin(phi1) * std::sin(phi2)
    );
    
    // Normalize longitude to -180 to 180
    double lon2 = rad2deg(lambda2);
    while (lon2 > 180) lon2 -= 360;
    while (lon2 < -180) lon2 += 360;
    
    return GeoPoint(lon2, rad2deg(phi2));
}

double calculateBearing(double lat1, double lon1, double lat2, double lon2) {
    double phi1 = deg2rad(lat1);
    double phi2 = deg2rad(lat2);
    double dLambda = deg2rad(lon2 - lon1);
    
    double x = std::cos(phi2) * std::sin(dLambda);
    double y = std::cos(phi1) * std::sin(phi2) -
               std::sin(phi1) * std::cos(phi2) * std::cos(dLambda);
    
    double theta = std::atan2(x, y);
    double bearing = rad2deg(theta);
    
    // Normalize to 0-360
    while (bearing < 0) bearing += 360;
    while (bearing >= 360) bearing -= 360;
    
    return bearing;
}

} // namespace Geodetic

} // namespace FSRM
