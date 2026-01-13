#ifndef FRACTURE_MODEL_HPP
#define FRACTURE_MODEL_HPP

/**
 * @file FractureModel.hpp
 * @brief Core fracture model classes for reservoir simulation
 * 
 * This file provides the base FractureModel class and derived classes for:
 * - NaturalFractureNetwork: Natural/existing fracture systems
 * - HydraulicFractureModel: Induced hydraulic fractures (PKN, KGD, P3D)
 * - FaultModel: Fault mechanics with rate-state friction
 * 
 * For advanced discrete fracture network (DFN) algorithms including:
 * - Stochastic DFN generation with Fisher/Bingham/power-law distributions
 * - Graph-based connectivity analysis and percolation
 * - MINC and EDFM upscaling methods
 * - Oda's crack tensor and Snow's permeability models
 * - LEFM-based fracture propagation
 * 
 * See FractureNetwork.hpp for the comprehensive FractureNetwork class
 * and related algorithms.
 */

#include "FSRM.hpp"
#include <vector>
#include <memory>

namespace FSRM {

// Forward declarations for advanced fracture network algorithms
// See FractureNetwork.hpp for full definitions
class FractureNetwork;
class FractureSet;
class FractureConnectivityGraph;
class FracturePropagationModel;
class FractureFlowSimulator;
struct DiscreteFracture;
struct FractureIntensity;

// Base fracture model
class FractureModel {
public:
    FractureModel(FractureType type);
    virtual ~FractureModel() = default;
    
    // Fracture geometry
    virtual void setGeometry(const std::vector<double>& coords) = 0;
    virtual void updateGeometry(double dt) = 0;
    
    // Fracture properties
    void setAperture(double aperture);
    void setPermeability(double perm);
    void setCompressibility(double comp);
    
    // Coupling with reservoir
    virtual void computeFluidExchange(Vec U, DM dm, 
                                     std::vector<double>& exchange_rates) const = 0;
    
    // Contribution to governing equations
    virtual void contributeToResidual(Vec F, Vec U, DM dm) const = 0;
    virtual void contributeToJacobian(Mat J, Vec U, DM dm) const = 0;
    
    FractureType getType() const { return frac_type; }
    
protected:
    FractureType frac_type;
    double aperture;
    double permeability;
    double compressibility;
    std::vector<double> coordinates;
};

// Natural fracture network
class NaturalFractureNetwork : public FractureModel {
public:
    NaturalFractureNetwork();
    
    void setGeometry(const std::vector<double>& coords) override;
    void updateGeometry(double dt) override;
    
    // Discrete Fracture Network (DFN) approach
    void addFracture(const std::vector<double>& points, double aperture);
    void generateStochasticNetwork(int num_fractures, double mean_length,
                                  double std_length, int seed);
    
    /**
     * @brief Import fractures from advanced FractureNetwork object
     * 
     * Converts DiscreteFracture objects from the FractureNetwork class
     * into the internal Fracture representation used by this model.
     * This allows using advanced DFN generation algorithms from
     * FractureNetwork.hpp with the reservoir simulation coupling.
     * 
     * @param network Pointer to FractureNetwork with generated fractures
     * @see FractureNetwork::generate() for DFN generation
     */
    void importFromFractureNetwork(const FractureNetwork* network);
    
    // Dual porosity/permeability model
    void enableDualPorosity(bool enable);
    void setShapeFactorModel(const std::string& model); // Warren-Root, Kazemi, etc.
    void setMatrixBlockSize(double size);
    
    // Fracture flow
    void computeFluidExchange(Vec U, DM dm,
                             std::vector<double>& exchange_rates) const override;
    
    void contributeToResidual(Vec F, Vec U, DM dm) const override;
    void contributeToJacobian(Mat J, Vec U, DM dm) const override;
    
private:
    struct Fracture {
        std::vector<double> points;
        double aperture;
        double permeability;
        int orientation; // 0=x, 1=y, 2=z dominant
    };
    
    std::vector<Fracture> fractures;
    bool use_dual_porosity;
    double shape_factor;
    double matrix_block_size;
    
    double computeShapeFactor(const std::string& model) const;
};

// Hydraulic fracture model
class HydraulicFractureModel : public FractureModel {
public:
    HydraulicFractureModel();
    
    void setGeometry(const std::vector<double>& coords) override;
    void updateGeometry(double dt) override;
    
    // PKN, KGD, or P3D model
    void setFractureModel(const std::string& model);
    
    // Fracture propagation
    void setPropagationCriteria(double Kc, double sigma_min);
    void computePropagationDirection(const Vec stress_field);
    void updateFractureGeometry(double pressure, double dt);
    
    // Proppant placement
    void enableProppantTransport(bool enable);
    void setProppantProperties(double diameter, double density, double concentration);
    void computeProppantDistribution(const Vec velocity_field);
    
    // Leak-off
    void enableLeakoff(bool enable);
    void setLeakoffCoefficient(double C_L);
    void computeLeakoff(Vec U, DM dm, std::vector<double>& leak_rates) const;
    
    void computeFluidExchange(Vec U, DM dm,
                             std::vector<double>& exchange_rates) const override;
    
    void contributeToResidual(Vec F, Vec U, DM dm) const override;
    void contributeToJacobian(Mat J, Vec U, DM dm) const override;
    
    // Fracture closure
    void checkClosure();
    void computeEffectivePermeability(const std::vector<double>& proppant_conc,
                                     double& k_eff) const;
    
private:
    std::string fracture_model; // PKN, KGD, P3D, planar3D
    
    // Geometry
    double length;
    double height;
    double width;
    std::vector<double> fracture_center;
    std::vector<double> fracture_normal;
    
    // Propagation
    double fracture_toughness;
    double min_horizontal_stress;
    double propagation_velocity;
    
    // Proppant
    bool transport_proppant;
    double proppant_diameter;
    double proppant_density;
    double proppant_concentration;
    std::vector<double> proppant_distribution;
    
    // Leak-off
    bool enable_leakoff_model;
    double leakoff_coefficient;
    
    // Formation/rock properties for fracture mechanics
    double youngs_modulus;        // Formation Young's modulus [Pa]
    double poisson_ratio;         // Formation Poisson's ratio
    double fluid_density;         // Fracturing fluid density [kg/m³]
    double fluid_viscosity;       // Fracturing fluid viscosity [Pa·s]
    double proppant_pack_porosity; // Proppant pack porosity
    
    // Propagation models
    void propagatePKN(double pressure, double dt);
    void propagateKGD(double pressure, double dt);
    void propagateP3D(double pressure, double dt);
    
public:
    // Setters for formation properties
    void setFormationProperties(double E, double nu);
    void setFluidProperties(double rho, double mu);
    void setProppantPackPorosity(double phi);
};

// Fault model
class FaultModel {
public:
    FaultModel();
    
    // Geometry
    void setFaultPlane(const std::vector<double>& strike,
                      const std::vector<double>& dip,
                      double length, double width);
    
    // Mechanical properties
    void setFrictionCoefficient(double static_mu, double dynamic_mu);
    void setCohesion(double cohesion);
    void setDilationAngle(double angle);
    
    // Fault slip
    enum class SlipMode {
        NO_SLIP,
        STICK_SLIP,
        STABLE_SLIP,
        SEISMIC
    };
    
    void computeSlip(const Vec stress_field, Vec displacement_field);
    void checkSlipCriteria(double shear_stress, double normal_stress,
                          SlipMode& mode) const;
    
    // Rate-and-state friction
    void enableRateStateFriction(bool enable);
    void setRateStateParameters(double a, double b, double Dc);
    void updateStateVariable(double slip_rate, double dt);
    
    // Fluid pressure effects
    void computeEffectiveStress(double pore_pressure, double& eff_normal_stress);
    
    // Seismicity
    void computeMomentMagnitude(double slip_area, double slip_amount, 
                               double& magnitude) const;
    
    // Permeability evolution
    void updatePermeability(double slip_amount);
    double getPermeability() const { return permeability; }
    
private:
    // Geometry
    std::vector<double> fault_strike;
    std::vector<double> fault_dip;
    double fault_length;
    double fault_width;
    
    // Mechanical properties
    double static_friction;
    double dynamic_friction;
    double cohesion;
    double dilation_angle;
    
    // Rate-state friction
    bool use_rate_state;
    double a_parameter;
    double b_parameter;
    double Dc_parameter;
    double state_variable;
    
    // Current state
    double permeability;
    double cumulative_slip;
    mutable SlipMode current_slip_mode;
};

} // namespace FSRM

#endif // FRACTURE_MODEL_HPP
