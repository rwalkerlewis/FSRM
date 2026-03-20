#ifndef FLUID_MODEL_HPP
#define FLUID_MODEL_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cmath>

namespace FSRM {

/**
 * @brief Enumeration of supported fluid model types
 */
enum class FluidType {
    SINGLE_PHASE,       // Single-phase incompressible/compressible
    BLACK_OIL,          // Three-phase black oil (oil, water, gas)
    COMPOSITIONAL,      // Multi-component with flash calculations
    DEAD_OIL,           // Oil without dissolved gas
    DRY_GAS,            // Gas without condensate
    WET_GAS,            // Gas with condensate
    BRINE,              // Saline water with salinity effects
    CO2,                // CO2 with supercritical behavior
    THERMAL             // Temperature-dependent properties
};

/**
 * @brief PVT correlation types for black oil models
 */
enum class PVTCorrelation {
    STANDING,           // Standing correlations
    VASQUEZ_BEGGS,      // Vasquez-Beggs correlations
    GLASO,              // Glaso correlations
    PETROSKY_FARSHAD,   // Petrosky-Farshad correlations
    AL_MARHOUN,         // Al-Marhoun correlations
    TABLE              // User-provided tables
};

/**
 * @brief Equation of State types for compositional models
 */
enum class EOSType {
    PENG_ROBINSON,      // Peng-Robinson EOS
    SRK,                // Soave-Redlich-Kwong EOS
    PR78,               // Peng-Robinson (1978) modification
    ZUDKEVITCH_JOFFE,   // Z-J modification
    IDEAL_GAS           // Ideal gas (for testing)
};

/**
 * @brief Base class for all fluid models
 * 
 * This provides a common interface for fluid property calculations
 * that can be configured entirely from text configuration files.
 */
class FluidModelBase {
public:
    FluidModelBase(FluidType type) : fluid_type(type) {}
    virtual ~FluidModelBase() = default;
    
    // Pure virtual methods for property calculations
    virtual double getDensity(double P, double T) const = 0;
    virtual double getViscosity(double P, double T) const = 0;
    virtual double getCompressibility(double P, double T) const = 0;
    
    // Configuration from key-value pairs
    virtual void configure(const std::map<std::string, std::string>& config) = 0;
    
    FluidType getType() const { return fluid_type; }
    
protected:
    FluidType fluid_type;
};

/**
 * @brief Single-phase fluid model
 * 
 * Supports incompressible and compressible single-phase fluids.
 * All properties configurable from file.
 */
class SinglePhaseFluid : public FluidModelBase {
public:
    SinglePhaseFluid();
    
    double getDensity(double P, double T) const override;
    double getViscosity(double P, double T) const override;
    double getCompressibility(double P, double T) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Setters for programmatic configuration
    void setDensity(double rho_ref, double P_ref = 1e5, double T_ref = 293.15);
    void setViscosity(double mu_ref, double P_ref = 1e5, double T_ref = 293.15);
    void setCompressibility(double c_total);
    void setThermalExpansion(double beta);
    void setViscosityModel(const std::string& model);
    
private:
    double density_ref;         // Reference density (kg/m³)
    double P_reference;         // Reference pressure (Pa)
    double T_reference;         // Reference temperature (K)
    double viscosity_ref;       // Reference viscosity (Pa·s)
    double compressibility;     // Total compressibility (1/Pa)
    double thermal_expansion;   // Thermal expansion coefficient (1/K)
    double viscosity_P_coeff;   // Viscosity pressure coefficient
    double viscosity_T_coeff;   // Viscosity temperature coefficient
    std::string viscosity_model; // "constant", "exponential", "arrhenius"
};

/**
 * @brief Black oil fluid model (three-phase: oil, water, gas)
 * 
 * Complete implementation with multiple PVT correlations.
 * Supports solution gas-oil ratio, formation volume factors,
 * relative permeabilities, and capillary pressures.
 */
class BlackOilFluid : public FluidModelBase {
public:
    BlackOilFluid();
    
    // Basic interface (averages over phases)
    double getDensity(double P, double T) const override;
    double getViscosity(double P, double T) const override;
    double getCompressibility(double P, double T) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Phase-specific properties
    double getOilDensity(double P, double Rs) const;
    double getGasDensity(double P) const;
    double getWaterDensity(double P) const;
    
    double getOilViscosity(double P, double Rs) const;
    double getGasViscosity(double P) const;
    double getWaterViscosity(double P) const;
    
    double getOilFVF(double P, double Rs) const;  // Bo
    double getGasFVF(double P) const;             // Bg
    double getWaterFVF(double P) const;           // Bw
    
    double getSolutionGOR(double P) const;        // Rs
    double getBubblePointPressure() const;        // Pb
    
    // Relative permeabilities (Corey model)
    double getKrw(double Sw) const;
    double getKro(double So, double Sg) const;
    double getKrg(double Sg) const;
    
    // Capillary pressures
    double getPcow(double Sw) const;
    double getPcog(double Sg) const;
    
    // Setters
    void setPVTCorrelation(PVTCorrelation corr);
    void setOilProperties(double rho_std, double mu_dead, double API);
    void setGasProperties(double rho_std, double gamma_g);
    void setWaterProperties(double rho_std, double mu, double salinity);
    void setSolutionGOR(double Rs_max, double Pb);
    void setCoreyParameters(double Swc, double Sor, double Sgc, 
                           double nw, double no, double ng,
                           double krw_max, double kro_max, double krg_max);
    
private:
    // PVT correlation type
    PVTCorrelation pvt_correlation;
    
    // Reference conditions (standard conditions)
    double P_std;               // Standard pressure (Pa)
    double T_std;               // Standard temperature (K)
    
    // Oil properties
    double oil_density_std;     // Oil density at SC (kg/m³)
    double oil_viscosity_dead;  // Dead oil viscosity (Pa·s)
    double oil_API;             // API gravity
    double oil_compressibility; // Oil compressibility (1/Pa)
    
    // Gas properties  
    double gas_density_std;     // Gas density at SC (kg/m³)
    double gas_gravity;         // Gas specific gravity (air=1)
    double gas_compressibility; // Gas compressibility (1/Pa)
    
    // Water properties
    double water_density_std;   // Water density at SC (kg/m³)
    double water_viscosity;     // Water viscosity (Pa·s)
    double water_salinity;      // Salinity (ppm or mass fraction)
    double water_compressibility;
    
    // Solution gas properties
    double Rs_max;              // Maximum solution GOR (scf/stb or sm³/sm³)
    double Pb;                  // Bubble point pressure (Pa)
    
    // Corey relative permeability parameters
    double Swc;                 // Connate water saturation
    double Sor;                 // Residual oil saturation
    double Sgc;                 // Critical gas saturation
    double nw, no, ng;          // Corey exponents
    double krw_max, kro_max, krg_max;  // Endpoint relative perms
    
    // Capillary pressure parameters
    double Pc_entry;            // Entry pressure (Pa)
    double lambda_pc;           // Pore size distribution index
    
    // Helper functions for correlations
    double standingRs(double P) const;
    double standingBo(double P, double Rs) const;
    double standingMuoLive(double P, double Rs) const;
    
    double vazquezBeggsRs(double P) const;
    double vazquezBeggsBo(double P, double Rs) const;
    double vazquesBeggsViscosity(double P, double Rs) const;
};

/**
 * @brief Compositional fluid model with equation of state
 * 
 * Multi-component fluid with flash calculations using
 * Peng-Robinson or SRK equations of state.
 */
class CompositionalFluid : public FluidModelBase {
public:
    CompositionalFluid(int num_components = 3);
    
    double getDensity(double P, double T) const override;
    double getViscosity(double P, double T) const override;
    double getCompressibility(double P, double T) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Component management
    void addComponent(const std::string& name, double Mw, double Tc, 
                     double Pc, double omega, double Vc = 0.0);
    void setComposition(const std::vector<double>& z);
    
    // Flash calculation
    struct FlashResult {
        std::vector<double> x;      // Liquid composition
        std::vector<double> y;      // Vapor composition
        double L;                   // Liquid mole fraction
        double V;                   // Vapor mole fraction
        double rho_L;               // Liquid density
        double rho_V;               // Vapor density
        bool converged;
    };
    
    FlashResult flash(double P, double T) const;
    FlashResult flash(double P, double T, const std::vector<double>& z) const;
    
    // Phase properties
    double getLiquidDensity(double P, double T, const std::vector<double>& x) const;
    double getVaporDensity(double P, double T, const std::vector<double>& y) const;
    double getLiquidViscosity(double P, double T, const std::vector<double>& x) const;
    double getVaporViscosity(double P, double T, const std::vector<double>& y) const;
    
    // EOS calculations
    void setEOS(EOSType type);
    double getCompressibilityFactor(double P, double T, const std::vector<double>& comp, 
                                    bool vapor) const;
    std::vector<double> getFugacityCoefficients(double P, double T, 
                                                const std::vector<double>& comp,
                                                bool vapor) const;
    
private:
    int nc;  // Number of components
    EOSType eos_type;
    
    // Component properties
    std::vector<std::string> component_names;
    std::vector<double> Mw;     // Molecular weights (g/mol)
    std::vector<double> Tc;     // Critical temperatures (K)
    std::vector<double> Pc;     // Critical pressures (Pa)
    std::vector<double> omega;  // Acentric factors
    std::vector<double> Vc;     // Critical volumes (m³/mol)
    
    // Binary interaction parameters (nc x nc matrix)
    std::vector<std::vector<double>> kij;
    
    // Current overall composition
    std::vector<double> z_global;
    
    // Cached flash results
    mutable FlashResult last_flash;
    mutable double last_P, last_T;
    
    // EOS helper functions
    double getPRa(int i, double T) const;
    double getPRb(int i) const;
    double getMixturea(const std::vector<double>& comp, double T) const;
    double getMixtureb(const std::vector<double>& comp) const;
    
    // Cubic equation solver
    std::vector<double> solveCubic(double A, double B) const;
};

/**
 * @brief Brine fluid model with salinity effects
 */
class BrineFluid : public FluidModelBase {
public:
    BrineFluid();
    
    double getDensity(double P, double T) const override;
    double getViscosity(double P, double T) const override;
    double getCompressibility(double P, double T) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    void setSalinity(double S);  // kg/kg or ppm
    void setIonicStrength(double I);
    
    // Specific brine correlations
    double getActivityCoefficient() const;
    double getOsmoticPressure(double T) const;
    
private:
    double salinity;            // Mass fraction
    double ionic_strength;
    double pure_water_density;
    double pure_water_viscosity;
};

/**
 * @brief CO2 fluid model with supercritical behavior
 */
class CO2Fluid : public FluidModelBase {
public:
    CO2Fluid();
    
    double getDensity(double P, double T) const override;
    double getViscosity(double P, double T) const override;
    double getCompressibility(double P, double T) const override;
    
    void configure(const std::map<std::string, std::string>& config) override;
    
    // Phase identification
    bool isSupercritical(double P, double T) const;
    bool isLiquid(double P, double T) const;
    bool isGas(double P, double T) const;
    
    // Solubility
    double getSolubilityInWater(double P, double T) const;
    double getSolubilityInOil(double P, double T) const;
    
private:
    static constexpr double Tc_CO2 = 304.13;  // K
    static constexpr double Pc_CO2 = 7.38e6;  // Pa
    static constexpr double rhoc_CO2 = 467.6; // kg/m³
    
    // Span-Wagner EOS coefficients (simplified)
    double getReducedDensity(double P, double T) const;
};

/**
 * @brief Factory function to create fluid models from configuration
 */
std::unique_ptr<FluidModelBase> createFluidModel(
    const std::string& type_str,
    const std::map<std::string, std::string>& config);

/**
 * @brief Parse fluid type from string
 */
FluidType parseFluidType(const std::string& type_str);

} // namespace FSRM

#endif // FLUID_MODEL_HPP
