#ifndef STRESS_SHADOWING_HPP
#define STRESS_SHADOWING_HPP

/**
 * @file StressShadowing.hpp
 * @brief Stress shadowing and inter-fracture mechanical interactions
 * 
 * ResFrac-equivalent capabilities:
 * - Stress perturbation from hydraulic fractures
 * - Inter-cluster stress interference
 * - Optimal cluster/stage spacing
 * - Fracture reorientation near wellbore
 * - Parent-child well stress interactions
 * - 3D stress field computation
 */

#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <cmath>

namespace FSRM {

/**
 * @brief 3D stress tensor
 */
struct Stress3D {
    double sxx, syy, szz;           ///< Normal stresses
    double sxy, sxz, syz;           ///< Shear stresses
    
    Stress3D() : sxx(0), syy(0), szz(0), sxy(0), sxz(0), syz(0) {}
    
    Stress3D(double xx, double yy, double zz, 
             double xy, double xz, double yz)
        : sxx(xx), syy(yy), szz(zz), sxy(xy), sxz(xz), syz(yz) {}
    
    // Add two stress tensors
    Stress3D operator+(const Stress3D& other) const;
    Stress3D& operator+=(const Stress3D& other);
    
    // Multiply by scalar
    Stress3D operator*(double s) const;
    
    // Principal stresses
    std::array<double, 3> principalStresses() const;
    std::array<std::array<double, 3>, 3> principalDirections() const;
    
    // Invariants
    double I1() const { return sxx + syy + szz; }
    double I2() const;
    double I3() const;
    
    // Mean and deviatoric
    double meanStress() const { return I1() / 3.0; }
    Stress3D deviatoricStress() const;
    double vonMises() const;
};

/**
 * @brief Point in 3D space
 */
struct Point3D {
    double x, y, z;
    
    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    double distanceTo(const Point3D& other) const;
    Point3D operator+(const Point3D& other) const;
    Point3D operator-(const Point3D& other) const;
    Point3D operator*(double s) const;
    double dot(const Point3D& other) const;
    Point3D cross(const Point3D& other) const;
    double norm() const;
    Point3D normalized() const;
};

/**
 * @brief Fracture geometry for stress calculation
 */
struct FractureGeometry {
    Point3D center;                 ///< Fracture center
    Point3D normal;                 ///< Fracture normal (opening direction)
    Point3D strike;                 ///< Strike direction
    Point3D dip;                    ///< Dip direction
    
    double half_length;             ///< Half-length in strike direction (m)
    double half_height;             ///< Half-height in dip direction (m)
    double average_width;           ///< Average aperture (m)
    
    // Net pressure (fluid pressure - closure stress)
    double net_pressure;            ///< Net pressure (Pa)
    
    // Discretization for numerical integration
    int num_elements_length;        ///< Elements along length
    int num_elements_height;        ///< Elements along height
    
    FractureGeometry();
    
    // Check if point is inside fracture
    bool containsPoint(const Point3D& p) const;
    
    // Get fracture area
    double getArea() const { return 4.0 * half_length * half_height; }
    
    // Get fracture volume
    double getVolume() const { return getArea() * average_width; }
};

/**
 * @brief Sneddon solution for penny-shaped crack
 */
class SneddonSolution {
public:
    SneddonSolution();
    
    /**
     * @brief Calculate stress perturbation from penny-shaped fracture
     * 
     * @param obs_point Observation point
     * @param frac_center Fracture center
     * @param frac_radius Fracture radius
     * @param frac_normal Fracture normal direction
     * @param net_pressure Net pressure in fracture
     * @param youngs_modulus Rock Young's modulus (Pa)
     * @param poisson_ratio Rock Poisson's ratio
     * @return Stress tensor perturbation
     */
    Stress3D calculateStress(const Point3D& obs_point,
                             const Point3D& frac_center,
                             double frac_radius,
                             const Point3D& frac_normal,
                             double net_pressure,
                             double youngs_modulus,
                             double poisson_ratio) const;
    
    /**
     * @brief Calculate displacement from penny-shaped fracture
     * 
     * @return Displacement vector
     */
    Point3D calculateDisplacement(const Point3D& obs_point,
                                   const Point3D& frac_center,
                                   double frac_radius,
                                   const Point3D& frac_normal,
                                   double net_pressure,
                                   double youngs_modulus,
                                   double poisson_ratio) const;
                                   
private:
    // Elliptic integrals
    double completeEllipticK(double m) const;
    double completeEllipticE(double m) const;
};

/**
 * @brief Rectangular fracture stress solution (DDM-based)
 */
class RectangularFractureSolution {
public:
    RectangularFractureSolution();
    
    /**
     * @brief Calculate stress from rectangular fracture
     * 
     * Uses displacement discontinuity method (DDM).
     * 
     * @param obs_point Observation point
     * @param fracture Fracture geometry
     * @param youngs_modulus Rock Young's modulus (Pa)
     * @param poisson_ratio Rock Poisson's ratio
     * @return Stress tensor perturbation
     */
    Stress3D calculateStress(const Point3D& obs_point,
                             const FractureGeometry& fracture,
                             double youngs_modulus,
                             double poisson_ratio) const;
    
    /**
     * @brief Calculate displacement from rectangular fracture
     */
    Point3D calculateDisplacement(const Point3D& obs_point,
                                   const FractureGeometry& fracture,
                                   double youngs_modulus,
                                   double poisson_ratio) const;
    
    /**
     * @brief Set numerical integration order
     */
    void setIntegrationOrder(int n) { integration_order_ = n; }
    
private:
    int integration_order_;
    
    // Okada (1985) dislocation solutions
    Stress3D okadaSolution(const Point3D& p, double length, double width,
                           double dip, double displacement,
                           double youngs_modulus, double poisson_ratio) const;
};

/**
 * @brief Stress shadow calculator for multiple fractures
 */
class StressShadowCalculator {
public:
    StressShadowCalculator();
    
    // Set rock properties
    void setRockProperties(double youngs_modulus, double poisson_ratio);
    
    // Set in-situ stress
    void setInSituStress(const Stress3D& stress);
    void setInSituStressGradients(double shmin_grad, double shmax_grad, 
                                   double sv_grad, double pp_grad);
    void setReferenceDepth(double depth) { ref_depth_ = depth; }
    
    // Add fractures
    void addFracture(const FractureGeometry& fracture);
    void clearFractures();
    
    /**
     * @brief Calculate total stress at point
     * 
     * Includes in-situ stress + perturbations from all fractures.
     * 
     * @param point Observation point
     * @return Total stress tensor
     */
    Stress3D calculateTotalStress(const Point3D& point) const;
    
    /**
     * @brief Calculate stress perturbation only
     * 
     * @param point Observation point
     * @return Stress perturbation from fractures
     */
    Stress3D calculatePerturbation(const Point3D& point) const;
    
    /**
     * @brief Calculate minimum horizontal stress change
     * 
     * @param point Observation point
     * @return Delta Shmin (Pa)
     */
    double calculateShminChange(const Point3D& point) const;
    
    /**
     * @brief Calculate stress shadow profile along line
     * 
     * @param start Start point
     * @param end End point
     * @param num_points Number of sample points
     * @return Stress at each point
     */
    std::vector<Stress3D> calculateProfileAlongLine(
        const Point3D& start,
        const Point3D& end,
        int num_points) const;
    
    /**
     * @brief Calculate stress shadow on 2D grid
     * 
     * @param corner1 First corner of grid
     * @param corner2 Opposite corner
     * @param nx Grid points in x
     * @param ny Grid points in y
     * @return 2D array of stress values
     */
    std::vector<std::vector<Stress3D>> calculateStressField(
        const Point3D& corner1,
        const Point3D& corner2,
        int nx, int ny) const;
    
    /**
     * @brief Calculate optimal fracture spacing
     * 
     * Find spacing where stress shadow doesn't significantly
     * affect fracture initiation.
     * 
     * @param fracture_template Template fracture geometry
     * @param max_stress_increase Maximum acceptable Shmin increase (Pa)
     * @return Optimal spacing (m)
     */
    double calculateOptimalSpacing(const FractureGeometry& fracture_template,
                                    double max_stress_increase) const;
    
    /**
     * @brief Calculate cluster efficiency due to stress shadow
     * 
     * @param cluster_positions Positions of clusters
     * @param pumping_pressure Applied pressure
     * @return Expected efficiency for each cluster
     */
    std::vector<double> calculateClusterEfficiency(
        const std::vector<Point3D>& cluster_positions,
        double pumping_pressure) const;
    
    // Output
    void writeVTKStressField(const std::string& filename,
                              const Point3D& corner1,
                              const Point3D& corner2,
                              int nx, int ny, int nz) const;
    
private:
    double youngs_modulus_;
    double poisson_ratio_;
    
    Stress3D in_situ_stress_;
    double shmin_grad_, shmax_grad_, sv_grad_, pp_grad_;
    double ref_depth_;
    
    std::vector<FractureGeometry> fractures_;
    
    // Solution methods
    RectangularFractureSolution rect_solution_;
    SneddonSolution sneddon_solution_;
    
    // Use rectangular or circular solution based on aspect ratio
    Stress3D calculateSingleFracturePerturbation(
        const Point3D& point,
        const FractureGeometry& fracture) const;
    
    // Get in-situ stress at depth
    Stress3D getInSituStressAtDepth(double depth) const;
};

/**
 * @brief Fracture reorientation predictor
 */
class FractureReorientation {
public:
    FractureReorientation();
    
    /**
     * @brief Calculate expected fracture orientation
     * 
     * Fractures may reorient near wellbore due to perforation
     * direction, near-wellbore stress, and existing fractures.
     * 
     * @param position Position in reservoir
     * @param wellbore_axis Wellbore axis direction
     * @param perforation_direction Perforation orientation
     * @param stress_calculator Stress shadow calculator
     * @return Expected fracture normal direction
     */
    Point3D predictFractureOrientation(
        const Point3D& position,
        const Point3D& wellbore_axis,
        const Point3D& perforation_direction,
        const StressShadowCalculator& stress_calculator) const;
    
    /**
     * @brief Calculate tortuosity factor
     * 
     * Higher tortuosity reduces effective conductivity.
     * 
     * @param initial_direction Initial fracture direction
     * @param final_direction Far-field direction
     * @param distance Distance over which reorientation occurs
     * @return Tortuosity factor (1.0 = no tortuosity)
     */
    double calculateTortuosity(const Point3D& initial_direction,
                                const Point3D& final_direction,
                                double distance) const;
    
    /**
     * @brief Check if fracture will turn towards or away from existing
     * 
     * @param new_frac_position New fracture position
     * @param existing_fracs Existing fracture geometries
     * @return Turn direction (positive = towards, negative = away)
     */
    double checkFractureTurning(
        const Point3D& new_frac_position,
        const std::vector<FractureGeometry>& existing_fracs) const;
    
    // Settings
    void setNearWellboreRadius(double r) { near_wellbore_radius_ = r; }
    
private:
    double near_wellbore_radius_;
};

/**
 * @brief Parent-child well stress interaction calculator
 */
class ParentChildInteraction {
public:
    ParentChildInteraction();
    
    // Set parent well fractures
    void setParentFractures(const std::vector<FractureGeometry>& fractures);
    
    // Set depletion (pressure drop from parent production)
    void setDepletion(double pressure_drop, double time_since_parent);
    
    // Set rock properties
    void setRockProperties(double E, double nu, double biot_coeff);
    
    /**
     * @brief Calculate stress change at child well due to parent
     * 
     * Includes both stress shadow and poroelastic effects from depletion.
     * 
     * @param child_position Position along child wellbore
     * @return Total stress change tensor
     */
    Stress3D calculateStressAtChild(const Point3D& child_position) const;
    
    /**
     * @brief Calculate frac hit risk
     * 
     * Determines probability that child fracture will connect to parent.
     * 
     * @param child_cluster_position Child cluster position
     * @param child_frac_azimuth Expected child fracture direction
     * @param pumping_pressure Child treatment pressure
     * @return Risk score (0-1) and distance to nearest parent frac
     */
    std::pair<double, double> calculateFracHitRisk(
        const Point3D& child_cluster_position,
        double child_frac_azimuth,
        double pumping_pressure) const;
    
    /**
     * @brief Calculate optimal child well placement
     * 
     * Find child well trajectory that minimizes frac hit risk
     * while maximizing drainage area.
     * 
     * @param parent_heel Parent well heel position
     * @param parent_toe Parent well toe position
     * @param target_spacing Target well spacing
     * @return Recommended child trajectory points
     */
    std::vector<Point3D> optimizeChildWellPlacement(
        const Point3D& parent_heel,
        const Point3D& parent_toe,
        double target_spacing) const;
    
    /**
     * @brief Calculate depleted stress zone
     * 
     * Returns the extent of stress alteration around parent fractures.
     * 
     * @return Approximate radius of stress alteration (m)
     */
    double calculateDepletedZoneRadius() const;
    
private:
    std::vector<FractureGeometry> parent_fractures_;
    double depletion_pressure_;
    double depletion_time_;
    
    double youngs_modulus_;
    double poisson_ratio_;
    double biot_coefficient_;
    
    StressShadowCalculator stress_calc_;
    
    // Poroelastic stress change from depletion
    Stress3D calculatePoroelasticStress(const Point3D& position) const;
};

/**
 * @brief Helper function to convert between coordinate systems
 */
Point3D localToGlobal(const Point3D& local, 
                       const Point3D& origin,
                       const Point3D& x_axis,
                       const Point3D& y_axis,
                       const Point3D& z_axis);

Point3D globalToLocal(const Point3D& global,
                       const Point3D& origin,
                       const Point3D& x_axis,
                       const Point3D& y_axis,
                       const Point3D& z_axis);

} // namespace FSRM

#endif // STRESS_SHADOWING_HPP
