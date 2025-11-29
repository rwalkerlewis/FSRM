#ifndef PROPPANT_TRANSPORT_HPP
#define PROPPANT_TRANSPORT_HPP

/**
 * @file ProppantTransport.hpp
 * @brief Advanced 2D proppant transport in hydraulic fractures
 * 
 * ResFrac-equivalent capabilities:
 * - 2D proppant concentration field in fracture
 * - Settling and convection
 * - Screenout prediction
 * - Bridging mechanics
 * - Proppant crushing and embedment
 * - Multiple proppant types/sizes
 * - Bank formation and dune transport
 * - Flowback prediction
 */

#include <vector>
#include <memory>
#include <array>
#include <string>
#include <map>
#include <functional>
#include <cmath>

namespace FSRM {

/**
 * @brief Proppant particle properties
 */
struct ProppantParticle {
    std::string name;                   ///< Particle type name
    double diameter;                    ///< Mean diameter (m)
    double density;                     ///< Particle density (kg/m³)
    double sphericity;                  ///< Sphericity factor (0-1)
    double roundness;                   ///< Roundness factor (0-1)
    
    // Size distribution (optional)
    std::vector<double> size_bins;      ///< Size bin edges (m)
    std::vector<double> mass_fractions; ///< Mass fraction in each bin
    
    // Strength properties
    double crush_strength;              ///< K-value crush strength (Pa)
    double embedment_strength;          ///< Formation hardness for embedment
    
    // Permeability correlation
    // k = k_ref * (σ/σ_ref)^(-α) where σ is closure stress
    double k_ref;                       ///< Reference permeability (m²)
    double stress_ref;                  ///< Reference stress (Pa)
    double stress_exponent;             ///< Stress sensitivity exponent
    
    // Non-Darcy coefficient
    double beta_factor;                 ///< Forchheimer beta (1/m)
    
    ProppantParticle();
    
    // Get settling velocity in given fluid
    double getSettlingVelocity(double fluid_viscosity, 
                               double fluid_density) const;
    
    // Get permeability at given stress
    double getPermeability(double closure_stress) const;
    
    // Get pack porosity
    double getPackPorosity() const;
};

/**
 * @brief Slurry (fluid + proppant) properties
 */
class SlurryModel {
public:
    SlurryModel();
    
    // Set base fluid properties
    void setFluidProperties(double density, double viscosity, 
                            double n_prime, double k_prime);
    
    // Add proppant
    void setProppant(const ProppantParticle& proppant);
    
    /**
     * @brief Calculate slurry density
     * 
     * @param proppant_concentration Proppant mass concentration (kg/m³)
     * @return Slurry density (kg/m³)
     */
    double getSlurryDensity(double proppant_concentration) const;
    
    /**
     * @brief Calculate slurry viscosity
     * 
     * Uses various correlations for concentrated suspensions.
     * 
     * @param proppant_concentration Proppant concentration (kg/m³)
     * @param shear_rate Shear rate (1/s)
     * @return Effective viscosity (Pa·s)
     */
    double getSlurryViscosity(double proppant_concentration,
                              double shear_rate) const;
    
    /**
     * @brief Calculate volume fraction from mass concentration
     * 
     * @param mass_concentration Mass concentration (kg/m³ slurry)
     * @return Volume fraction (dimensionless)
     */
    double getVolumeFraction(double mass_concentration) const;
    
    /**
     * @brief Get maximum packing fraction
     */
    double getMaxPackingFraction() const { return max_packing_; }
    
    // Viscosity model selection
    enum class ViscosityModel {
        EINSTEIN,                       ///< Einstein dilute suspension
        KRIEGER_DOUGHERTY,             ///< Krieger-Dougherty
        QUEMADA,                        ///< Quemada model
        THOMAS                          ///< Thomas correlation
    };
    
    void setViscosityModel(ViscosityModel model) { visc_model_ = model; }
    
private:
    double fluid_density_;
    double fluid_viscosity_;
    double n_prime_, k_prime_;
    
    ProppantParticle proppant_;
    ViscosityModel visc_model_;
    double max_packing_;                ///< Maximum packing fraction (~0.65)
    double intrinsic_viscosity_;        ///< [η] for Einstein model (~2.5)
};

/**
 * @brief 2D proppant transport solver in fracture
 */
class Proppant2DTransport {
public:
    Proppant2DTransport();
    
    /**
     * @brief Initialize fracture grid
     * 
     * @param length Fracture half-length (m)
     * @param height Fracture height (m)
     * @param nx Grid points along length
     * @param nz Grid points along height
     */
    void initializeGrid(double length, double height, int nx, int nz);
    
    /**
     * @brief Set fracture width field
     * 
     * @param width 2D array of width values (m)
     */
    void setWidthField(const std::vector<std::vector<double>>& width);
    
    /**
     * @brief Set velocity field in fracture
     * 
     * @param vx Velocity in length direction (m/s)
     * @param vz Velocity in height direction (m/s)
     */
    void setVelocityField(const std::vector<std::vector<double>>& vx,
                          const std::vector<std::vector<double>>& vz);
    
    // Set proppant and fluid
    void setProppant(const ProppantParticle& proppant);
    void setSlurryModel(const SlurryModel& slurry);
    
    /**
     * @brief Set inlet boundary condition
     * 
     * @param inlet_concentration Inlet proppant concentration (kg/m³)
     * @param inlet_rate Inlet slurry rate per unit height (m²/s)
     */
    void setInletCondition(double inlet_concentration, double inlet_rate);
    
    /**
     * @brief Advance solution one time step
     * 
     * @param dt Time step (s)
     * @return True if converged, false if screenout occurred
     */
    bool step(double dt);
    
    /**
     * @brief Get proppant concentration field
     * 
     * @return 2D concentration field (kg/m³)
     */
    const std::vector<std::vector<double>>& getConcentration() const {
        return concentration_;
    }
    
    /**
     * @brief Get proppant bank height
     * 
     * @return Height of proppant bank at each x position (m)
     */
    std::vector<double> getBankHeight() const;
    
    /**
     * @brief Get proppant coverage (effective propped area)
     * 
     * @return Fraction of fracture area with significant proppant
     */
    double getProppantCoverage() const;
    
    /**
     * @brief Check for screenout condition
     * 
     * @return True if screenout detected
     */
    bool isScreenout() const;
    
    /**
     * @brief Get screenout location
     * 
     * @return (x, z) position of screenout
     */
    std::pair<double, double> getScreenoutLocation() const;
    
    /**
     * @brief Calculate average conductivity
     * 
     * @param closure_stress Closure stress (Pa)
     * @return Area-weighted average conductivity (mD·ft)
     */
    double getAverageConductivity(double closure_stress) const;
    
    /**
     * @brief Get conductivity field
     * 
     * @param closure_stress Closure stress (Pa)
     * @return 2D conductivity field (m²·m)
     */
    std::vector<std::vector<double>> getConductivityField(
        double closure_stress) const;
    
    // Numerical settings
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    void setMinWidth(double w) { min_width_ = w; }
    void setScreenoutThreshold(double threshold) { screenout_threshold_ = threshold; }
    
    // Output
    void writeVTK(const std::string& filename) const;
    void writeConcentrationProfile(std::ostream& os) const;
    
private:
    // Grid
    double length_, height_;
    int nx_, nz_;
    double dx_, dz_;
    
    // Fields
    std::vector<std::vector<double>> width_;
    std::vector<std::vector<double>> vx_, vz_;
    std::vector<std::vector<double>> concentration_;
    std::vector<std::vector<double>> bank_indicator_;  ///< 1.0 if banked, 0.0 otherwise
    
    // Models
    ProppantParticle proppant_;
    SlurryModel slurry_;
    
    // Boundary conditions
    double inlet_concentration_;
    double inlet_rate_;
    
    // Numerical parameters
    int max_iterations_;
    double tolerance_;
    double min_width_;
    double screenout_threshold_;
    
    // State
    bool screenout_occurred_;
    int screenout_i_, screenout_j_;
    
    // Solver components
    void computeSettlingVelocity(std::vector<std::vector<double>>& v_settle);
    void computeHinderedSettling(double local_conc, double& settling_factor);
    void advectionStep(double dt);
    void settlingStep(double dt);
    void bankFormation(double dt);
    void checkScreenout();
    
    // Flux limiters
    double minmod(double a, double b) const;
    double superbee(double a, double b) const;
};

/**
 * @brief Proppant bridging model
 */
class ProppantBridging {
public:
    ProppantBridging();
    
    /**
     * @brief Check if bridging occurs
     * 
     * Bridging typically occurs when width < N * diameter
     * where N depends on proppant concentration and sphericity.
     * 
     * @param width Fracture width (m)
     * @param proppant_diameter Proppant diameter (m)
     * @param concentration Proppant concentration (kg/m³)
     * @param sphericity Particle sphericity
     * @return True if bridging occurs
     */
    bool checkBridging(double width, double proppant_diameter,
                       double concentration, double sphericity) const;
    
    /**
     * @brief Calculate bridging width ratio
     * 
     * @param concentration Volume concentration
     * @param sphericity Particle sphericity
     * @return Critical width/diameter ratio
     */
    double getBridgingRatio(double concentration, double sphericity) const;
    
    /**
     * @brief Calculate pressure drop due to bridging
     * 
     * @param width Fracture width (m)
     * @param proppant_diameter Particle diameter (m)
     * @param concentration Concentration (kg/m³)
     * @param flow_rate Flow rate (m³/s)
     * @param fluid_viscosity Fluid viscosity (Pa·s)
     * @return Additional pressure drop (Pa)
     */
    double calculateBridgingPressureDrop(double width, 
                                          double proppant_diameter,
                                          double concentration,
                                          double flow_rate,
                                          double fluid_viscosity) const;
    
    // Model parameters
    void setBridgingMultiplier(double mult) { bridging_multiplier_ = mult; }
    void setMinBridgingConcentration(double c) { min_bridging_conc_ = c; }
    
private:
    double bridging_multiplier_;        ///< Adjusts bridging criterion
    double min_bridging_conc_;          ///< Minimum concentration for bridging
};

/**
 * @brief Proppant crushing and embedment model
 */
class ProppantDamage {
public:
    ProppantDamage();
    
    /**
     * @brief Calculate proppant crushing
     * 
     * @param closure_stress Effective closure stress (Pa)
     * @param proppant Proppant properties
     * @return Fraction of proppant crushed (0-1)
     */
    double calculateCrushing(double closure_stress,
                              const ProppantParticle& proppant) const;
    
    /**
     * @brief Calculate embedment depth
     * 
     * @param closure_stress Effective closure stress (Pa)
     * @param proppant_diameter Proppant diameter (m)
     * @param formation_hardness Formation Brinell hardness (Pa)
     * @return Embedment depth per side (m)
     */
    double calculateEmbedment(double closure_stress,
                               double proppant_diameter,
                               double formation_hardness) const;
    
    /**
     * @brief Calculate effective width after embedment
     * 
     * @param initial_width Initial propped width (m)
     * @param closure_stress Closure stress (Pa)
     * @param proppant_diameter Proppant diameter (m)
     * @param formation_hardness Formation hardness (Pa)
     * @return Effective width (m)
     */
    double getEffectiveWidth(double initial_width,
                              double closure_stress,
                              double proppant_diameter,
                              double formation_hardness) const;
    
    /**
     * @brief Calculate conductivity reduction factor
     * 
     * Combines crushing, embedment, and fines migration effects.
     * 
     * @param closure_stress Closure stress (Pa)
     * @param proppant Proppant properties
     * @param formation_hardness Formation hardness (Pa)
     * @return Conductivity multiplier (0-1)
     */
    double getConductivityMultiplier(double closure_stress,
                                      const ProppantParticle& proppant,
                                      double formation_hardness) const;
    
    // Model calibration
    void setCrushingExponent(double n) { crushing_exponent_ = n; }
    void setEmbedmentCoefficient(double c) { embedment_coeff_ = c; }
    void setFinesMultiplier(double m) { fines_multiplier_ = m; }
    
private:
    double crushing_exponent_;
    double embedment_coeff_;
    double fines_multiplier_;
};

/**
 * @brief Multiple proppant types transport
 */
class MultiProppantTransport {
public:
    MultiProppantTransport();
    
    /**
     * @brief Initialize with fracture geometry
     */
    void initializeGrid(double length, double height, int nx, int nz);
    
    /**
     * @brief Add proppant type
     * 
     * @param name Proppant name
     * @param proppant Proppant properties
     */
    void addProppantType(const std::string& name, 
                          const ProppantParticle& proppant);
    
    /**
     * @brief Set pumping schedule for proppants
     * 
     * @param schedule Vector of (time, proppant_name, concentration)
     */
    void setPumpingSchedule(
        const std::vector<std::tuple<double, std::string, double>>& schedule);
    
    /**
     * @brief Advance simulation
     * 
     * @param dt Time step (s)
     * @param current_time Current simulation time (s)
     * @return Success flag
     */
    bool step(double dt, double current_time);
    
    /**
     * @brief Get total concentration field
     */
    std::vector<std::vector<double>> getTotalConcentration() const;
    
    /**
     * @brief Get concentration field for specific proppant
     */
    const std::vector<std::vector<double>>& getConcentration(
        const std::string& name) const;
    
    /**
     * @brief Calculate effective conductivity with all proppants
     */
    std::vector<std::vector<double>> getEffectiveConductivity(
        double closure_stress) const;
    
private:
    std::map<std::string, ProppantParticle> proppants_;
    std::map<std::string, std::vector<std::vector<double>>> concentrations_;
    std::map<std::string, Proppant2DTransport> transporters_;
    
    std::vector<std::tuple<double, std::string, double>> pump_schedule_;
    
    double length_, height_;
    int nx_, nz_;
};

/**
 * @brief Proppant flowback predictor
 */
class ProppantFlowback {
public:
    ProppantFlowback();
    
    /**
     * @brief Calculate critical flowback rate
     * 
     * Rate above which proppant will flow back with produced fluids.
     * 
     * @param proppant Proppant properties
     * @param fracture_width Propped width (m)
     * @param closure_stress Closure stress (Pa)
     * @param fluid_density Produced fluid density (kg/m³)
     * @param fluid_viscosity Produced fluid viscosity (Pa·s)
     * @return Critical production rate (m³/s)
     */
    double getCriticalFlowbackRate(const ProppantParticle& proppant,
                                    double fracture_width,
                                    double closure_stress,
                                    double fluid_density,
                                    double fluid_viscosity) const;
    
    /**
     * @brief Estimate proppant production
     * 
     * @param production_rate Current production rate (m³/s)
     * @param critical_rate Critical flowback rate (m³/s)
     * @param fracture_area Fracture area (m²)
     * @param proppant_volume Proppant volume in fracture (m³)
     * @param dt Time step (s)
     * @return Mass of proppant produced (kg)
     */
    double estimateProppantProduction(double production_rate,
                                       double critical_rate,
                                       double fracture_area,
                                       double proppant_volume,
                                       double dt) const;
    
    /**
     * @brief Calculate resin-coated proppant stability
     * 
     * @param proppant Proppant (must be resin-coated)
     * @param temperature Downhole temperature (K)
     * @param time Exposure time (s)
     * @return Bond strength development (0-1)
     */
    double getResinCureLevel(const ProppantParticle& proppant,
                              double temperature,
                              double time) const;
    
    // Settings
    void setFlowbackMultiplier(double m) { flowback_mult_ = m; }
    
private:
    double flowback_mult_;
};

/**
 * @brief Parse proppant database from file
 */
std::map<std::string, ProppantParticle> loadProppantDatabase(
    const std::string& filename);

} // namespace FSRM

#endif // PROPPANT_TRANSPORT_HPP
