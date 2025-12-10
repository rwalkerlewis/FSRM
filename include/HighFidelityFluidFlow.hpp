#ifndef HIGH_FIDELITY_FLUID_FLOW_HPP
#define HIGH_FIDELITY_FLUID_FLOW_HPP

/**
 * @file HighFidelityFluidFlow.hpp
 * @brief High-fidelity fluid flow models for complex reservoir simulation
 * 
 * This file provides advanced fluid flow formulations that extend the basic
 * Darcy flow capabilities in FluidModel.hpp and TwoPhaseFlow.hpp. These are
 * optional high-fidelity extensions that can be enabled when greater physical
 * accuracy is required.
 * 
 * Features:
 * - Non-Darcy flow (Forchheimer, turbulent, inertial corrections)
 * - Dual-porosity/dual-permeability for fractured reservoirs
 * - Multi-continuum models (matrix-fracture-vugs)
 * - Non-isothermal multiphase flow with phase change
 * - Miscible/near-miscible displacement
 * - Dynamic relative permeability with rate effects
 * - Compositional effects with EOS coupling
 * 
 * All models integrate with existing FluidModel, RelativePermeability, and
 * TwoPhaseFlow classes without replacing them.
 * 
 * @note These are high-fidelity options - basic simulations should use
 *       the simpler models in FluidModel.hpp
 */

#include "FluidModel.hpp"
#include "RelativePermeability.hpp"
#include "TwoPhaseFlow.hpp"
#include "MaterialModel.hpp"
#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <map>
#include <cmath>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Non-Darcy Flow Models
// =============================================================================

/**
 * @brief Non-Darcy flow correction types
 * 
 * At high velocities, deviations from Darcy's law become significant.
 * These models capture inertial and turbulent effects.
 */
enum class NonDarcyModel {
    DARCY_ONLY,          ///< Standard Darcy flow (no correction)
    FORCHHEIMER,         ///< Forchheimer inertial correction: -∇P = μ/k·v + β·ρ·v²
    BARREE_CONWAY,       ///< Barree-Conway model for proppant packs
    ERGUN,               ///< Ergun equation for packed beds
    TURBULENT_DARCY,     ///< Turbulent flow with effective permeability
    KLINKENBERG,         ///< Gas slippage effect at low pressures
    KNUDSEN              ///< Knudsen diffusion for tight formations
};

/**
 * @brief Non-Darcy flow coefficient models
 */
enum class BetaCoefficientModel {
    CONSTANT,            ///< User-specified constant β
    COOKE,               ///< Cooke (1973) correlation
    GEERTSMA,            ///< Geertsma (1974) correlation
    FREDERICK_GRAVES,    ///< Frederick & Graves correlation
    EVANS_CIVAN,         ///< Evans & Civan correlation
    THAUVIN_MOHANTY      ///< Thauvin & Mohanty (1998)
};

/**
 * @brief Forchheimer non-Darcy flow model
 * 
 * Extends Darcy's law with an inertial term:
 * -∇P = (μ/k)·v + β·ρ·|v|·v
 * 
 * where β is the non-Darcy flow coefficient (1/m).
 * 
 * This is critical for:
 * - High-rate wells
 * - Near-wellbore flow
 * - Flow through proppant packs
 * - Fracture flow at high rates
 */
class ForchheimerFlow {
public:
    struct Parameters {
        double beta = 1e5;              ///< Non-Darcy coefficient [1/m]
        double min_permeability = 1e-20; ///< Minimum permeability [m²]
        double critical_Re = 10.0;       ///< Critical Reynolds number
        bool use_multiphase = true;      ///< Apply to multiphase flow
        BetaCoefficientModel beta_model = BetaCoefficientModel::GEERTSMA;
        
        Parameters() = default;
    };
    
    ForchheimerFlow();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate non-Darcy corrected velocity
     * 
     * Solves the Forchheimer equation implicitly for velocity given
     * pressure gradient.
     * 
     * @param grad_P Pressure gradient magnitude [Pa/m]
     * @param k Permeability [m²]
     * @param mu Viscosity [Pa·s]
     * @param rho Density [kg/m³]
     * @return Darcy velocity [m/s]
     */
    double calculateVelocity(double grad_P, double k, double mu, double rho) const;
    
    /**
     * @brief Calculate effective permeability reduction
     * 
     * Returns an apparent permeability that gives the same flow rate
     * when used in Darcy's law: k_app = k / (1 + Fo)
     * 
     * @param k Intrinsic permeability [m²]
     * @param mu Viscosity [Pa·s]
     * @param rho Density [kg/m³]
     * @param velocity Darcy velocity [m/s]
     * @return Apparent permeability [m²]
     */
    double getApparentPermeability(double k, double mu, double rho, double velocity) const;
    
    /**
     * @brief Calculate Forchheimer number
     * 
     * Fo = β·ρ·k·|v|/μ
     * 
     * Fo >> 1 indicates strong non-Darcy effects.
     */
    double getForchheimerNumber(double k, double mu, double rho, double velocity) const;
    
    /**
     * @brief Calculate beta coefficient from correlations
     * 
     * @param k Permeability [m²]
     * @param phi Porosity [-]
     * @return Beta coefficient [1/m]
     */
    double calculateBeta(double k, double phi) const;
    
    /**
     * @brief Multiphase extension
     * 
     * For multiphase flow, applies correction to total mobility:
     * λ_t = Σ(kr_i/μ_i) → λ_t_eff accounting for phase velocities
     */
    double getEffectiveMobility(double k, 
                                const std::vector<double>& kr,
                                const std::vector<double>& mu,
                                const std::vector<double>& rho,
                                const std::vector<double>& S,
                                double grad_P) const;
    
    /**
     * @brief Calculate Reynolds number for flow regime
     */
    double getReynoldsNumber(double k, double phi, double rho, double mu, double velocity) const;
    
    /**
     * @brief Check if non-Darcy effects are significant
     */
    bool isNonDarcySignificant(double k, double mu, double rho, double velocity,
                               double threshold = 0.1) const;
    
private:
    Parameters params_;
    
    // Newton iteration for velocity from Forchheimer equation
    double solveForchheimer(double grad_P, double k, double mu, double rho) const;
    
    // Correlation functions for beta
    double betaCooke(double k) const;
    double betaGeertsma(double k, double phi) const;
    double betaFrederickGraves(double k) const;
    double betaEvansCivan(double k, double phi) const;
    double betaThauvinMohanty(double k, double phi) const;
};

/**
 * @brief Klinkenberg gas slippage model
 * 
 * For gas flow in tight formations, mean free path effects cause
 * apparent permeability to exceed intrinsic permeability:
 * k_app = k_∞ * (1 + b/P)
 * 
 * where b is the Klinkenberg slip factor.
 */
class KlinkenbergCorrection {
public:
    struct Parameters {
        double slip_factor = 1e5;        ///< Klinkenberg slip factor b [Pa]
        double reference_pressure = 1e5; ///< Reference pressure [Pa]
        bool use_florence_model = false; ///< Use Florence (2007) model
        double pore_diameter = 1e-8;     ///< Mean pore diameter [m] (for Florence)
    };
    
    KlinkenbergCorrection();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate apparent permeability with slip correction
     * 
     * @param k_infinity Intrinsic permeability [m²]
     * @param P Pressure [Pa]
     * @param T Temperature [K] (optional for Florence model)
     * @param Mw Molecular weight [kg/mol] (optional)
     * @return Apparent permeability [m²]
     */
    double getApparentPermeability(double k_infinity, double P, 
                                   double T = 293.15, double Mw = 0.028) const;
    
    /**
     * @brief Calculate slip factor from pore geometry
     * 
     * Uses kinetic theory of gases:
     * b = 4*c*λ*P / r
     * 
     * where λ is mean free path and c is a constant.
     */
    double calculateSlipFactor(double pore_radius, double P, double T, double Mw) const;
    
    /**
     * @brief Calculate Knudsen number
     * 
     * Kn = λ/d where λ is mean free path and d is pore diameter.
     * Kn > 0.01 indicates slip flow regime.
     */
    double getKnudsenNumber(double P, double T, double Mw, double pore_diameter) const;
    
private:
    Parameters params_;
    
    double meanFreePath(double P, double T, double Mw) const;
};


// =============================================================================
// Dual-Porosity / Dual-Permeability Models
// =============================================================================

/**
 * @brief Transfer function type for matrix-fracture flow
 */
enum class TransferFunctionType {
    KAZEMI,              ///< Kazemi et al. (1976) shape factor
    WARREN_ROOT,         ///< Warren & Root (1963) classic
    COATS,               ///< Coats (1989) for multiphase
    GILMAN_KAZEMI,       ///< Gilman & Kazemi (1983)
    LEMONNIER_BOURBIAUX, ///< Lemonnier & Bourbiaux gravity term
    QUANDALLE_SABATHIER, ///< Quandalle & Sabathier (1989)
    PRUESS_NARASIMHAN,   ///< MINC-style nested continua
    SUBDOMAIN            ///< Subdomain/MINC approach
};

/**
 * @brief Shape factor calculation method
 */
enum class ShapeFactorModel {
    KAZEMI,              ///< σ = 4(1/Lx² + 1/Ly² + 1/Lz²)
    COATS,               ///< σ = 8/L² (single block)
    GILMAN,              ///< σ = 4(n² + 1)/L²
    LEMONNIER,           ///< Gravity-corrected
    UEDA,                ///< Ueda et al. shape factor
    CUSTOMIZABLE         ///< User-defined
};

/**
 * @brief Dual-porosity/dual-permeability model
 * 
 * Models flow in fractured reservoirs with:
 * - Matrix continuum (storage, low permeability)
 * - Fracture continuum (flow paths, high permeability)
 * 
 * Can operate in:
 * - Dual-porosity (matrix-fracture transfer, no matrix-matrix flow)
 * - Dual-permeability (both continua flow, with exchange)
 * 
 * Governing equations:
 * Matrix: ∂(φ_m·ρ)/∂t - ∇·(ρ·k_m/μ·∇P_m) + T_mf = 0
 * Fracture: ∂(φ_f·ρ)/∂t - ∇·(ρ·k_f/μ·∇P_f) - T_mf + Q = 0
 * 
 * Transfer: T_mf = σ·k_m/μ·(P_m - P_f)
 */
class DualPorosityModel {
public:
    struct MatrixProperties {
        double porosity = 0.2;           ///< Matrix porosity [-]
        double permeability = 1e-18;     ///< Matrix permeability [m²]
        double compressibility = 1e-9;   ///< Matrix compressibility [1/Pa]
        double block_size_x = 1.0;       ///< Matrix block dimension X [m]
        double block_size_y = 1.0;       ///< Matrix block dimension Y [m]
        double block_size_z = 1.0;       ///< Matrix block dimension Z [m]
    };
    
    struct FractureProperties {
        double porosity = 0.01;          ///< Fracture porosity [-]
        double permeability = 1e-12;     ///< Fracture permeability [m²]
        double compressibility = 1e-8;   ///< Fracture compressibility [1/Pa]
        double spacing = 1.0;            ///< Fracture spacing [m]
        double aperture = 1e-4;          ///< Fracture aperture [m]
    };
    
    struct Parameters {
        bool dual_permeability = false;  ///< Allow matrix-matrix flow
        TransferFunctionType transfer_type = TransferFunctionType::KAZEMI;
        ShapeFactorModel shape_model = ShapeFactorModel::KAZEMI;
        double shape_factor_multiplier = 1.0;
        bool include_gravity = true;     ///< Include gravity in transfer
        bool include_capillary = true;   ///< Include capillary in transfer
        int num_minc_shells = 5;         ///< MINC shells (if using MINC)
    };
    
    DualPorosityModel();
    
    void setMatrixProperties(const MatrixProperties& props);
    void setFractureProperties(const FractureProperties& props);
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate shape factor for matrix-fracture transfer
     * 
     * @return Shape factor σ [1/m²]
     */
    double calculateShapeFactor() const;
    
    /**
     * @brief Calculate single-phase transfer rate
     * 
     * T_mf = σ · (k_m/μ) · (P_m - P_f)
     * 
     * @param P_matrix Matrix pressure [Pa]
     * @param P_fracture Fracture pressure [Pa]
     * @param mu Viscosity [Pa·s]
     * @return Transfer rate [m³/s/m³ bulk]
     */
    double calculateTransferRate(double P_matrix, double P_fracture, double mu) const;
    
    /**
     * @brief Calculate multiphase transfer rate for a phase
     * 
     * Uses upstream weighting for relative permeability.
     * 
     * @param P_m Matrix phase pressure [Pa]
     * @param P_f Fracture phase pressure [Pa]
     * @param kr_m Matrix relative permeability [-]
     * @param kr_f Fracture relative permeability [-]
     * @param mu Phase viscosity [Pa·s]
     * @param rho Phase density [kg/m³]
     * @param depth_diff Height difference [m] (for gravity)
     * @return Phase transfer rate [m³/s/m³ bulk]
     */
    double calculatePhaseTransfer(double P_m, double P_f,
                                  double kr_m, double kr_f,
                                  double mu, double rho, double depth_diff = 0.0) const;
    
    /**
     * @brief Calculate capillary-driven imbibition transfer
     * 
     * Counter-current flow driven by capillary pressure difference.
     */
    double calculateCapillaryTransfer(double Sw_m, double Sw_f,
                                      double Pc_m, double Pc_f,
                                      double mu_w, double mu_o) const;
    
    /**
     * @brief MINC (Multiple INteracting Continua) transfer
     * 
     * Discretizes matrix into nested shells for transient drainage.
     * More accurate than single-transfer-function for early time.
     */
    struct MINCState {
        std::vector<double> shell_pressures;
        std::vector<double> shell_saturations;
        std::vector<double> shell_volumes;
    };
    
    void initializeMINC(MINCState& state, double P_init, double S_init) const;
    void updateMINC(MINCState& state, double P_fracture, double S_fracture,
                    double mu, double dt) const;
    double getMINCTransferRate(const MINCState& state, double P_fracture, double mu) const;
    
    // Accessors
    const MatrixProperties& getMatrixProperties() const { return matrix_; }
    const FractureProperties& getFractureProperties() const { return fracture_; }
    const Parameters& getParameters() const { return params_; }
    
private:
    MatrixProperties matrix_;
    FractureProperties fracture_;
    Parameters params_;
    double shape_factor_;
    
    void updateShapeFactor();
    
    // Shape factor calculations
    double shapeFactor_Kazemi() const;
    double shapeFactor_Coats() const;
    double shapeFactor_Gilman() const;
    double shapeFactor_Lemonnier() const;
};

/**
 * @brief Triple-porosity / Triple-permeability model
 * 
 * Extends dual-porosity for reservoirs with three distinct continua:
 * - Matrix (primary storage)
 * - Natural fractures (secondary flow)
 * - Vugs/Karst (tertiary features)
 * 
 * Common in carbonate reservoirs with dissolution features.
 */
class TriplePorosityModel {
public:
    struct VugProperties {
        double porosity = 0.05;          ///< Vug porosity [-]
        double permeability = 1e-10;     ///< Vug permeability [m²]
        double compressibility = 1e-7;   ///< Vug compressibility [1/Pa]
        double connectivity = 0.5;       ///< Vug connectivity factor [-]
    };
    
    TriplePorosityModel();
    
    void setVugProperties(const VugProperties& props);
    
    // Inherits matrix/fracture from DualPorosityModel
    void setDualPorosityBase(const DualPorosityModel& dp);
    
    /**
     * @brief Calculate all inter-continuum transfer rates
     * 
     * Returns rates for: matrix-fracture, fracture-vug, matrix-vug
     */
    std::array<double, 3> calculateAllTransfers(double P_m, double P_f, double P_v,
                                                 double mu) const;
    
private:
    DualPorosityModel dual_porosity_;
    VugProperties vugs_;
    double shape_factor_mv_;  // Matrix-vug
    double shape_factor_fv_;  // Fracture-vug
};


// =============================================================================
// Non-Isothermal Multiphase Flow
// =============================================================================

/**
 * @brief Non-isothermal multiphase flow model for reservoir simulations
 * 
 * Extends standard multiphase flow with thermal effects:
 * - Temperature-dependent fluid properties
 * - Heat transport (conduction + convection)
 * - Phase change (vaporization/condensation)
 * - Thermal expansion effects
 * 
 * Couples with TwoPhaseFlow and FluidModel for property calculations.
 * 
 * @note Relationship to other thermal modules:
 * 
 * - **ThermalPressurization.hpp**: Fault-focused thermal physics for earthquake
 *   mechanics. Models frictional heating and 1D diffusion across narrow fault
 *   shear zones. Use for dynamic rupture simulations.
 * 
 * - **PhysicsKernel.hpp (ThermalKernel)**: Basic FE thermal kernel for heat
 *   conduction/convection. This class extends those capabilities with
 *   multiphase-specific effects.
 * 
 * - **AdvancedCoupledPhysics.hpp (THMCoupling)**: Coupling framework that can
 *   use this class's thermal calculations for reservoir THM simulations.
 * 
 * This class is specifically designed for reservoir/subsurface flow with
 * multiple fluid phases. For geothermal or EOR applications.
 */
class NonIsothermalMultiphaseFlow {
public:
    struct ThermalProperties {
        // Rock thermal properties
        double rock_thermal_conductivity = 2.5;   ///< k_rock [W/(m·K)]
        double rock_specific_heat = 900.0;        ///< c_rock [J/(kg·K)]
        double rock_density = 2500.0;             ///< ρ_rock [kg/m³]
        
        // Fluid thermal properties (reference values)
        double water_thermal_conductivity = 0.6;  ///< k_w [W/(m·K)]
        double water_specific_heat = 4186.0;      ///< c_w [J/(kg·K)]
        double oil_thermal_conductivity = 0.15;   ///< k_o [W/(m·K)]
        double oil_specific_heat = 2000.0;        ///< c_o [J/(kg·K)]
        double gas_thermal_conductivity = 0.03;   ///< k_g [W/(m·K)]
        double gas_specific_heat = 2200.0;        ///< c_g [J/(kg·K)]
        
        // Thermal expansion
        double water_thermal_expansion = 2.1e-4;  ///< β_w [1/K]
        double oil_thermal_expansion = 7e-4;      ///< β_o [1/K]
        
        // Phase change
        double latent_heat_vaporization = 2.26e6; ///< L_v [J/kg]
        double boiling_point = 373.15;            ///< T_boil at 1 atm [K]
    };
    
    struct Parameters {
        bool include_conduction = true;
        bool include_convection = true;
        bool include_phase_change = false;
        bool include_joule_thomson = false;      ///< J-T cooling effect
        bool include_viscous_dissipation = false; ///< Frictional heating
        double reference_temperature = 293.15;    ///< T_ref [K]
    };
    
    NonIsothermalMultiphaseFlow();
    
    void setThermalProperties(const ThermalProperties& props);
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate effective thermal conductivity
     * 
     * Mixture model considering porosity and saturation:
     * k_eff = (1-φ)·k_rock + φ·(Sw·k_w + So·k_o + Sg·k_g)
     */
    double getEffectiveConductivity(double phi, double Sw, double So, double Sg) const;
    
    /**
     * @brief Calculate effective heat capacity
     * 
     * (ρc)_eff = (1-φ)·ρ_rock·c_rock + φ·Σ(S_i·ρ_i·c_i)
     */
    double getEffectiveHeatCapacity(double phi, double Sw, double So, double Sg,
                                    double rho_w, double rho_o, double rho_g) const;
    
    /**
     * @brief Calculate heat flux (conductive + convective)
     * 
     * q_heat = -k_eff·∇T + Σ(ρ_i·h_i·v_i)
     * 
     * @param grad_T Temperature gradient [K/m]
     * @param velocities Phase Darcy velocities [m/s]
     * @param rho Phase densities [kg/m³]
     * @param T Temperature [K]
     * @return Heat flux [W/m²]
     */
    std::array<double, 3> calculateHeatFlux(double phi, double Sw, double So, double Sg,
                                            const std::array<double, 3>& grad_T,
                                            const std::array<double, 3>& vel_w,
                                            const std::array<double, 3>& vel_o,
                                            const std::array<double, 3>& vel_g,
                                            double rho_w, double rho_o, double rho_g,
                                            double T) const;
    
    /**
     * @brief Calculate temperature-dependent viscosity
     * 
     * Uses Andrade equation: μ = A·exp(B/T)
     */
    double getViscosity(double mu_ref, double T_ref, double T, double activation_energy) const;
    
    /**
     * @brief Calculate temperature-dependent density
     * 
     * ρ(T) = ρ_ref / (1 + β·(T - T_ref))
     */
    double getDensity(double rho_ref, double T_ref, double T, double beta) const;
    
    /**
     * @brief Calculate phase change rate (evaporation/condensation)
     * 
     * Returns mass transfer rate between liquid and vapor phases.
     * 
     * @param P Pressure [Pa]
     * @param T Temperature [K]
     * @param Sw Water saturation [-]
     * @param Sg Gas saturation [-]
     * @return Mass transfer rate [kg/(m³·s)] (positive = evaporation)
     */
    double calculatePhaseChangeRate(double P, double T, double Sw, double Sg) const;
    
    /**
     * @brief Calculate Joule-Thomson coefficient
     * 
     * μ_JT = -(1/ρ·c_p)·[T·(∂V/∂T)_P - V]
     */
    double getJouleThomsonCoefficient(double P, double T, double rho, double cp,
                                      double thermal_expansion) const;
    
private:
    ThermalProperties thermal_;
    Parameters params_;
    
    double calculateEnthalpy(double T, double phase_cp, double latent = 0.0) const;
    double saturationPressure(double T) const;  // Clausius-Clapeyron
};


// =============================================================================
// Dynamic/Rate-Dependent Relative Permeability
// =============================================================================

/**
 * @brief Dynamic relative permeability model
 * 
 * Standard relative permeability assumes equilibrium saturation distribution.
 * At high flow rates or with viscous fingering, this assumption breaks down.
 * 
 * This model accounts for:
 * - Capillary number effects (low Ca = capillary dominated)
 * - Velocity-dependent endpoint shifts
 * - Transient drainage/imbibition
 * - Viscous fingering corrections
 */
class DynamicRelativePermeability {
public:
    struct Parameters {
        // Capillary number effects
        bool include_capillary_number = true;
        double critical_capillary_number = 1e-5;  ///< Ca_c for desaturation
        double Sor_at_high_Ca = 0.05;             ///< Residual oil at high Ca
        double Swc_at_high_Ca = 0.05;             ///< Connate water at high Ca
        
        // CDC (Capillary Desaturation Curve) parameters
        double cdc_exponent = 2.0;                ///< Exponent in CDC model
        double cdc_transition_width = 1.0;        ///< Log-decades for transition
        
        // Viscous fingering
        bool include_fingering = false;
        double fingering_exponent = 0.5;          ///< Todd-Longstaff mixing parameter
        double mobility_ratio_threshold = 10.0;   ///< Above this, apply correction
        
        // Rate effects
        bool include_rate_effects = false;
        double rate_exponent = 0.2;               ///< Sensitivity to velocity
        double reference_velocity = 1e-5;         ///< Reference Darcy velocity [m/s]
    };
    
    DynamicRelativePermeability();
    
    void setParameters(const Parameters& params);
    void setBaseRelPerm(std::shared_ptr<ThreePhaseRelPerm> base_relperm);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate capillary number
     * 
     * Ca = μ·v / σ
     * 
     * @param mu Viscosity [Pa·s]
     * @param v Darcy velocity [m/s]
     * @param sigma Interfacial tension [N/m]
     * @return Capillary number [-]
     */
    double calculateCapillaryNumber(double mu, double v, double sigma) const;
    
    /**
     * @brief Get modified endpoints from CDC
     * 
     * Residual saturations decrease at high Ca:
     * Sor(Ca) = Sor_0 · [1 - exp(-(Ca/Ca_c)^n)]
     */
    void getModifiedEndpoints(double Ca, double& Sor, double& Swc) const;
    
    /**
     * @brief Calculate rate-adjusted relative permeability
     * 
     * @param Sw Water saturation
     * @param So Oil saturation  
     * @param Sg Gas saturation
     * @param velocity Darcy velocity [m/s]
     * @param mu_w Water viscosity [Pa·s]
     * @param mu_o Oil viscosity [Pa·s]
     * @param sigma Interfacial tension [N/m]
     * @return {kr_w, kr_o, kr_g}
     */
    std::array<double, 3> getRelPerm(double Sw, double So, double Sg,
                                      double velocity, double mu_w, double mu_o,
                                      double sigma) const;
    
    /**
     * @brief Todd-Longstaff viscous fingering correction
     * 
     * Effective viscosity in fingered zone:
     * μ_eff = μ_o^(1-ω) · μ_mix^ω
     * 
     * where ω is the mixing parameter (0 = no mixing, 1 = full mixing).
     */
    double getEffectiveViscosity(double mu_o, double mu_w, double fw, 
                                  double mixing_param) const;
    
    /**
     * @brief Calculate mobility ratio
     * 
     * M = (kr_w/μ_w) / (kr_o/μ_o)
     */
    double getMobilityRatio(double kr_w, double kr_o, double mu_w, double mu_o) const;
    
private:
    Parameters params_;
    std::shared_ptr<ThreePhaseRelPerm> base_relperm_;
    
    void applyFingeringCorrection(double M, std::array<double, 3>& kr) const;
};


// =============================================================================
// Miscible/Near-Miscible Flow
// =============================================================================

/**
 * @brief Miscibility state enumeration
 */
enum class MiscibilityState {
    IMMISCIBLE,          ///< Complete phase separation
    FIRST_CONTACT,       ///< FCM - immediate miscibility
    MULTI_CONTACT,       ///< MCM - develops miscibility
    NEAR_MISCIBLE        ///< Close to MMP but not miscible
};

/**
 * @brief Miscible displacement model
 * 
 * Models EOR processes where injected fluid mixes with reservoir oil:
 * - CO2 flooding
 * - Hydrocarbon gas injection
 * - Solvent floods
 * 
 * Key physics:
 * - Mixing zone development
 * - Effective properties in transition
 * - Minimum miscibility pressure (MMP)
 * - Dispersion effects
 */
class MiscibleFlowModel {
public:
    struct Parameters {
        // Miscibility parameters
        double MMP = 15e6;                        ///< Minimum miscibility pressure [Pa]
        double miscibility_transition_width = 2e6; ///< Width of transition region [Pa]
        MiscibilityState state = MiscibilityState::IMMISCIBLE;
        
        // Mixing parameters
        double longitudinal_dispersivity = 0.1;   ///< α_L [m]
        double transverse_dispersivity = 0.01;    ///< α_T [m]
        double molecular_diffusion = 1e-9;        ///< D_m [m²/s]
        
        // Todd-Longstaff
        double mixing_parameter = 0.67;           ///< ω (omega)
        
        // Solubility
        double solvent_solubility = 0.3;          ///< Max solvent in oil [-]
    };
    
    MiscibleFlowModel();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Determine miscibility state from pressure
     */
    MiscibilityState getMiscibilityState(double P) const;
    
    /**
     * @brief Calculate miscibility factor (0 = immiscible, 1 = miscible)
     * 
     * F_misc = (P - MMP_lower) / (MMP - MMP_lower)
     */
    double getMiscibilityFactor(double P) const;
    
    /**
     * @brief Calculate effective properties in mixing zone
     * 
     * Using Todd-Longstaff model:
     * ρ_eff = (ρ_s·S_s + ρ_o·S_o) / (S_s + S_o)  for fully mixed
     */
    void getEffectiveProperties(double S_solvent, double S_oil, double P, double T,
                                double rho_s, double rho_o, double mu_s, double mu_o,
                                double& rho_eff, double& mu_eff) const;
    
    /**
     * @brief Calculate dispersion coefficient
     * 
     * D_eff = D_m + α_L·|v| + α_T·|v_perp|
     */
    double getDispersionCoefficient(double velocity, bool longitudinal = true) const;
    
    /**
     * @brief Calculate relative permeability with miscibility
     * 
     * Interpolates between immiscible and miscible curves:
     * kr = F_misc·kr_misc + (1-F_misc)·kr_immisc
     */
    std::array<double, 3> getRelPermWithMiscibility(double Sw, double So, double Sg,
                                                     double P,
                                                     const std::array<double, 3>& kr_immisc) const;
    
private:
    Parameters params_;
};


// =============================================================================
// High-Fidelity Flow Solver Configuration
// =============================================================================

/**
 * @brief Configuration for high-fidelity flow solver
 * 
 * Aggregates all high-fidelity options into a single configuration
 * that can be parsed from config files.
 */
struct HighFidelityFlowConfig {
    // Master enable flags
    bool enable_high_fidelity = false;           ///< Master switch
    
    // Non-Darcy
    bool enable_non_darcy = false;
    NonDarcyModel non_darcy_model = NonDarcyModel::DARCY_ONLY;
    ForchheimerFlow::Parameters forchheimer_params;
    
    // Dual-porosity
    bool enable_dual_porosity = false;
    bool enable_dual_permeability = false;
    DualPorosityModel::Parameters dual_porosity_params;
    DualPorosityModel::MatrixProperties matrix_props;
    DualPorosityModel::FractureProperties fracture_props;
    
    // Triple-porosity (vugs)
    bool enable_triple_porosity = false;
    TriplePorosityModel::VugProperties vug_props;
    
    // Thermal
    bool enable_non_isothermal = false;
    NonIsothermalMultiphaseFlow::ThermalProperties thermal_props;
    NonIsothermalMultiphaseFlow::Parameters thermal_params;
    
    // Dynamic rel perm
    bool enable_dynamic_relperm = false;
    DynamicRelativePermeability::Parameters dynamic_relperm_params;
    
    // Miscible flow
    bool enable_miscible = false;
    MiscibleFlowModel::Parameters miscible_params;
    
    // Klinkenberg
    bool enable_klinkenberg = false;
    KlinkenbergCorrection::Parameters klinkenberg_params;
    
    // Parse from config file
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Validate configuration
    bool validate(std::string& error_msg) const;
};

/**
 * @brief Factory function to create high-fidelity flow model from config
 */
std::unique_ptr<ForchheimerFlow> createNonDarcyModel(
    const std::map<std::string, std::string>& config);

std::unique_ptr<DualPorosityModel> createDualPorosityModel(
    const std::map<std::string, std::string>& config);

std::unique_ptr<NonIsothermalMultiphaseFlow> createThermalFlowModel(
    const std::map<std::string, std::string>& config);

} // namespace HighFidelity
} // namespace FSRM

#endif // HIGH_FIDELITY_FLUID_FLOW_HPP
