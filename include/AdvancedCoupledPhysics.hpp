#ifndef ADVANCED_COUPLED_PHYSICS_HPP
#define ADVANCED_COUPLED_PHYSICS_HPP

/**
 * @file AdvancedCoupledPhysics.hpp
 * @brief Advanced multi-physics coupling for high-fidelity simulations
 * 
 * This file provides sophisticated coupling mechanisms that extend the basic
 * poroelastic coupling in PoroelasticSolver.hpp. These are optional high-fidelity
 * extensions for complex coupled simulations.
 * 
 * Features:
 * - Full Biot dynamic equations with all three wave types
 * - Thermo-Hydro-Mechanical (THM) coupling
 * - Thermo-Hydro-Mechanical-Chemical (THMC) coupling
 * - Multi-scale homogenization approaches
 * - Non-linear Biot parameters (pressure/stress-dependent)
 * - Unsaturated soil mechanics (Richards + Biot)
 * - Reactive transport coupling
 * - Fluid-Structure Interaction (FSI) for fractures
 * 
 * Integrates with existing:
 * - PoroelasticSolver.hpp for basic poroelastic coupling
 * - FluidModel.hpp for fluid properties
 * - MaterialModel.hpp for solid properties
 * - PhysicsKernel.hpp for kernel implementations
 * 
 * @note These are advanced coupling options - basic simulations should use
 *       the simpler formulations in PoroelasticSolver.hpp
 */

#include "PoroelasticSolver.hpp"
#include "FluidModel.hpp"
#include "MaterialModel.hpp"
#include "PhysicsKernel.hpp"
#include "ThermalPressurization.hpp"  // For fault-focused thermal pressurization
#include "HighFidelityFluidFlow.hpp"
#include "HighFidelityGeomechanics.hpp"
#include <vector>
#include <memory>
#include <array>
#include <functional>
#include <map>

namespace FSRM {
namespace HighFidelity {

// =============================================================================
// Full Biot Dynamic Poroelasticity
// =============================================================================

/**
 * @brief Biot wave types
 */
enum class BiotWaveType {
    FAST_P,              ///< Fast compressional (frame dominated)
    SLOW_P,              ///< Slow compressional (Biot type II, fluid dominated)
    SHEAR                ///< Shear wave
};

/**
 * @brief Biot formulation type
 */
enum class BiotFormulation {
    LOW_FREQUENCY,       ///< Quasi-static fluid (standard Biot)
    HIGH_FREQUENCY,      ///< Full dynamic including relative acceleration
    U_P,                 ///< u-p formulation (displacement-pressure)
    U_W,                 ///< u-w formulation (solid-fluid displacements)
    U_P_W                ///< u-p-w mixed formulation
};

/**
 * @brief Full Biot dynamic poroelasticity model
 * 
 * Complete formulation of Biot's equations for wave propagation in
 * fluid-saturated porous media:
 * 
 * Solid momentum: ρ·ü - ∇·σ' + α·∇p + ρ_f·ẅ = f_s
 * Fluid momentum: ρ_f·ü + ρ_a·ẅ + (μ/k)·ẇ + ∇p = f_f
 * Mass balance: ∂(α·∇·u + p/M)/∂t + ∇·ẇ = q
 * 
 * where:
 * - u = solid displacement
 * - w = relative fluid displacement (= φ·(U_f - u))
 * - p = pore pressure
 * - σ' = effective stress
 * - α = Biot coefficient
 * - M = Biot modulus
 * - ρ_a = added mass (tortuosity effect)
 * 
 * Supports three wave types:
 * - Fast P-wave (skeleton dominated, ~V_p1)
 * - Slow P-wave (Biot type II, ~V_p2, highly attenuated)
 * - Shear wave (~V_s)
 */
class FullBiotDynamics {
public:
    struct FluidProperties {
        double density = 1000.0;             ///< ρ_f [kg/m³]
        double viscosity = 0.001;            ///< μ [Pa·s]
        double bulk_modulus = 2.2e9;         ///< K_f [Pa]
    };
    
    struct SolidFrameProperties {
        double density = 2650.0;             ///< ρ_s (grain density) [kg/m³]
        double bulk_modulus = 36e9;          ///< K_s (grain bulk modulus) [Pa]
        double frame_bulk_modulus = 10e9;    ///< K_fr (frame/drained) [Pa]
        double shear_modulus = 8e9;          ///< G [Pa]
    };
    
    struct PorousMediaProperties {
        double porosity = 0.2;               ///< φ [-]
        double permeability = 1e-13;         ///< k [m²]
        double tortuosity = 2.0;             ///< τ [-]
        double biot_coefficient = 0.0;       ///< α (computed if 0)
        double biot_modulus = 0.0;           ///< M (computed if 0)
    };
    
    struct Parameters {
        BiotFormulation formulation = BiotFormulation::U_P;
        
        // Frequency regime
        double characteristic_frequency = 1.0; ///< f_c [Hz] transition frequency
        bool use_high_frequency_correction = true;
        
        // Damping
        double viscous_damping = 1.0;          ///< Viscous coupling multiplier
        bool include_inertial_coupling = true;  ///< Include added mass
        
        // Non-linear effects
        bool nonlinear_biot_coefficient = false;
        bool pressure_dependent_permeability = false;
        
        // Numerical
        double mass_matrix_lumping = 0.0;      ///< 0 = consistent, 1 = lumped
    };
    
    FullBiotDynamics();
    
    void setFluidProperties(const FluidProperties& props);
    void setSolidProperties(const SolidFrameProperties& props);
    void setPorousMediaProperties(const PorousMediaProperties& props);
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate Biot parameters from constituents
     * 
     * α = 1 - K_fr/K_s
     * M = [(α - φ)/K_s + φ/K_f]^(-1)
     */
    void computeBiotParameters();
    
    /**
     * @brief Get computed Biot coefficient
     */
    double getBiotCoefficient() const { return biot_alpha_; }
    
    /**
     * @brief Get computed Biot modulus
     */
    double getBiotModulus() const { return biot_M_; }
    
    /**
     * @brief Calculate wave velocities
     * 
     * Solves dispersion relation for three body waves.
     * 
     * @param frequency Frequency [Hz]
     * @return {V_p_fast, V_p_slow, V_s}
     */
    std::array<double, 3> waveVelocities(double frequency = 0.0) const;
    
    /**
     * @brief Calculate wave attenuations (Q^-1)
     * 
     * @param frequency Frequency [Hz]
     * @return {Q_fast^-1, Q_slow^-1, Q_s^-1}
     */
    std::array<double, 3> waveAttenuations(double frequency) const;
    
    /**
     * @brief Calculate characteristic frequency
     * 
     * f_c = φ·μ / (2π·ρ_f·k·τ)
     * 
     * Below f_c: viscosity-dominated (low-frequency Biot)
     * Above f_c: inertia-dominated (high-frequency Biot)
     */
    double characteristicFrequency() const;
    
    /**
     * @brief Get mass matrices for FE assembly
     * 
     * For u-p formulation:
     * [M_uu  M_up] [ü]   [K_uu  K_up] [u]   [f_u]
     * [M_pu  M_pp] [p̈] + [K_pu  K_pp] [p] = [f_p]
     * 
     * @param[out] M_uu Solid-solid mass
     * @param[out] M_up Solid-fluid coupling mass
     * @param[out] M_pp Fluid-fluid mass
     */
    void getMassMatrices(double& M_uu, double& M_up, double& M_pp) const;
    
    /**
     * @brief Get stiffness/coupling matrices
     * 
     * @param[out] K_uu Solid stiffness
     * @param[out] K_up Coupling
     * @param[out] K_pp Fluid compressibility
     */
    void getStiffnessMatrices(std::array<std::array<double, 6>, 6>& K_uu,
                              std::array<double, 6>& K_up,
                              double& K_pp) const;
    
    /**
     * @brief Calculate residual for u-p formulation
     * 
     * For implicit time integration.
     */
    void calculateResidual(const std::array<double, 3>& u,
                          const std::array<double, 3>& u_dot,
                          const std::array<double, 3>& u_ddot,
                          double p, double p_dot,
                          const std::array<double, 6>& grad_u,
                          const std::array<double, 3>& grad_p,
                          std::array<double, 3>& R_u,
                          double& R_p) const;
    
    /**
     * @brief Calculate Jacobian for Newton iteration
     */
    void calculateJacobian(const std::array<double, 3>& u,
                          double p,
                          const std::array<double, 6>& grad_u,
                          const std::array<double, 3>& grad_p,
                          double dt, double theta,
                          std::array<std::array<double, 4>, 4>& J) const;
    
    /**
     * @brief Non-linear Biot coefficient
     * 
     * α(σ) = α_0 + α_1·exp(-σ_eff/σ_ref)
     */
    double nonlinearBiotCoefficient(double effective_stress) const;
    
    /**
     * @brief Pressure-dependent permeability
     * 
     * k(p) = k_0·exp(α_k·(p - p_ref))
     */
    double pressureDependentPermeability(double pressure) const;
    
private:
    FluidProperties fluid_;
    SolidFrameProperties solid_;
    PorousMediaProperties porous_;
    Parameters params_;
    
    // Computed Biot parameters
    double biot_alpha_;
    double biot_M_;
    double added_mass_;  // ρ_a = (τ - 1)·ρ_f / φ
    
    // Effective densities
    double rho_bulk_;    // (1-φ)·ρ_s + φ·ρ_f
    double rho_11_;      // Solid inertia coefficient
    double rho_12_;      // Coupling coefficient
    double rho_22_;      // Fluid inertia coefficient
    
    void computeEffectiveDensities();
};


// =============================================================================
// Thermo-Hydro-Mechanical (THM) Coupling
// =============================================================================

/**
 * @brief THM coupling type
 */
enum class THMCouplingType {
    SEQUENTIAL,          ///< Staggered solve (T → H → M)
    ITERATIVE,           ///< Iterative coupling with convergence check
    MONOLITHIC,          ///< Fully coupled Newton
    OPERATOR_SPLIT       ///< Operator splitting with sub-cycling
};

/**
 * @brief Thermal coupling effects for THM simulations
 * 
 * @note For fault-specific thermal pressurization (frictional heating during
 *       dynamic rupture), use ThermalPressurization class from 
 *       ThermalPressurization.hpp instead. This struct is for general THM
 *       coupling in reservoir/geomechanics simulations.
 * 
 * Relationship to existing thermal modules:
 * - ThermalPressurization.hpp: Fault-focused, models 1D diffusion across shear
 *   zone with rate-state friction coupling. Use for earthquake simulations.
 * - PhysicsKernel.hpp (ThermalKernel): Basic heat conduction/convection FE kernel.
 *   THMCoupling uses this as its thermal solver.
 * - HighFidelityFluidFlow.hpp (NonIsothermalMultiphaseFlow): Multiphase reservoir
 *   thermal effects. Can be used with THMCoupling for reservoir simulations.
 */
struct ThermalCouplingEffects {
    // T → M (Thermal to Mechanical)
    double thermal_expansion_solid = 1e-5;    ///< α_T [1/K]
    double thermal_softening_modulus = 0.0;   ///< dE/dT [Pa/K]
    
    // T → H (Thermal to Hydraulic)
    double thermal_expansion_fluid = 2.1e-4;  ///< β_T [1/K]
    double viscosity_temperature_coeff = 0.02; ///< dln(μ)/dT
    
    // T → (H,M) (Thermal pressurization coefficient)
    // Note: For fault-zone thermal pressurization, use ThermalPressurization class
    // which provides full diffusion equations. This Lambda is for bulk THM coupling.
    double Lambda = 0.1e6;                    ///< ∂p/∂T at const vol [Pa/K]
    
    // Frictional heating (for bulk plastic dissipation)
    // Note: For fault frictional heating, use ThermalPressurization class
    double heat_partition = 0.5;              ///< Fraction to solid
    
    // Optional: Use ThermalPressurization for detailed fault physics
    bool use_fault_thermal_pressurization = false;
    std::shared_ptr<ThermalPressurization> fault_tp_model = nullptr;
};

/**
 * @brief Hydraulic coupling effects
 */
struct HydraulicCouplingEffects {
    // H → M (Pore pressure effects)
    double biot_coefficient = 0.7;            ///< α
    double biot_modulus = 10e9;               ///< M [Pa]
    
    // H → T (Advective heat transport)
    bool include_advection = true;
    
    // H → M (Effective stress)
    bool use_effective_stress = true;
    double pore_pressure_coefficient = 1.0;   ///< Terzaghi = 1, Biot = α
};

/**
 * @brief Mechanical coupling effects
 */
struct MechanicalCouplingEffects {
    // M → H (Volumetric strain effect on storage)
    double storage_compression_coefficient = 1e-9; ///< [1/Pa]
    
    // M → H (Permeability evolution)
    bool stress_dependent_permeability = false;
    double permeability_stress_coefficient = 1e-8; ///< dk/dσ
    double permeability_strain_coefficient = 1.0;  ///< Kozeny-Carman exponent
    
    // M → T (Mechanical dissipation)
    bool include_plastic_heating = false;
    double plastic_work_fraction = 0.9;       ///< Taylor-Quinney coefficient
};

/**
 * @brief Thermo-Hydro-Mechanical coupled model
 * 
 * Full THM coupling with:
 * 
 * Energy equation (T):
 * (ρc)_eff·∂T/∂t - ∇·(k_eff·∇T) + ρ_f·c_f·v·∇T = Q_h
 * 
 * Mass balance (H):
 * S·∂p/∂t + α·∂(∇·u)/∂t - ∇·(k/μ·∇p) = Q_f + Λ·∂T/∂t
 * 
 * Momentum balance (M):
 * ∇·σ' - α·∇p - 3K·α_T·∇T + ρ·g = 0
 * 
 * where:
 * - T = temperature
 * - p = pore pressure  
 * - u = displacement
 * - σ' = effective stress
 * - Λ = thermal pressurization coefficient
 */
class THMCoupling {
public:
    struct Parameters {
        THMCouplingType coupling_type = THMCouplingType::ITERATIVE;
        
        // Coupling effects
        ThermalCouplingEffects thermal_effects;
        HydraulicCouplingEffects hydraulic_effects;
        MechanicalCouplingEffects mechanical_effects;
        
        // Iterative coupling parameters
        int max_coupling_iterations = 20;
        double coupling_tolerance = 1e-6;
        double relaxation_factor = 1.0;       ///< Under-relaxation (< 1)
        
        // Time stepping
        bool adaptive_subcycling = true;
        double thermal_subcycle_factor = 10.0;
        double hydraulic_subcycle_factor = 1.0;
    };
    
    THMCoupling();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Set base models for each physics
     * 
     * THM coupling uses existing model implementations.
     */
    void setThermalModel(std::shared_ptr<FluidModelBase> fluid_thermal);
    void setHydraulicModel(std::shared_ptr<FluidModelBase> fluid);
    void setMechanicalModel(std::shared_ptr<MaterialModelBase> solid);
    
    /**
     * @brief Calculate effective thermal conductivity
     * 
     * k_eff = (1-φ)·k_s + φ·S·k_w + φ·(1-S)·k_g
     */
    double effectiveThermalConductivity(double phi, double Sw,
                                        double k_solid, double k_water,
                                        double k_gas) const;
    
    /**
     * @brief Calculate effective heat capacity
     * 
     * (ρc)_eff = (1-φ)·ρ_s·c_s + φ·S·ρ_w·c_w + φ·(1-S)·ρ_g·c_g
     */
    double effectiveHeatCapacity(double phi, double Sw,
                                 double rho_s, double c_s,
                                 double rho_w, double c_w,
                                 double rho_g, double c_g) const;
    
    /**
     * @brief Calculate thermal stress
     * 
     * σ_T = -3·K·α_T·(T - T_ref)·I
     */
    std::array<double, 6> thermalStress(double K, double alpha_T,
                                        double T, double T_ref) const;
    
    /**
     * @brief Calculate thermal pressurization rate
     * 
     * For undrained conditions:
     * ∂p/∂t = Λ·∂T/∂t
     */
    double thermalPressurizationRate(double dT_dt) const;
    
    /**
     * @brief Calculate viscosity at temperature
     * 
     * μ(T) = μ_ref·exp(-E_a/(R·T))
     */
    double viscosityAtTemperature(double mu_ref, double T_ref,
                                  double T, double E_activation) const;
    
    /**
     * @brief Perform sequential THM solve step
     * 
     * 1. Solve thermal (fixed H, M)
     * 2. Solve hydraulic (fixed T, M)
     * 3. Solve mechanical (fixed T, H)
     */
    bool sequentialStep(double dt,
                       std::vector<double>& T,
                       std::vector<double>& P,
                       std::vector<std::array<double, 3>>& U);
    
    /**
     * @brief Perform iterative THM solve
     * 
     * Iterates sequential solve until convergence.
     */
    bool iterativeStep(double dt,
                      std::vector<double>& T,
                      std::vector<double>& P,
                      std::vector<std::array<double, 3>>& U);
    
    /**
     * @brief Check coupling convergence
     */
    bool checkConvergence(const std::vector<double>& T_old,
                         const std::vector<double>& T_new,
                         const std::vector<double>& P_old,
                         const std::vector<double>& P_new) const;
    
    /**
     * @brief Calculate coupling matrix for monolithic solve
     * 
     * Returns blocks of the full THM Jacobian matrix.
     */
    void couplingMatrix(double phi, double k, double K, double G,
                       double alpha, double M, double alpha_T,
                       std::array<std::array<double, 3>, 3>& J_TT,
                       std::array<std::array<double, 3>, 3>& J_TH,
                       std::array<std::array<double, 3>, 3>& J_TM,
                       std::array<std::array<double, 3>, 3>& J_HT,
                       std::array<std::array<double, 3>, 3>& J_HH,
                       std::array<std::array<double, 3>, 3>& J_HM,
                       std::array<std::array<double, 3>, 3>& J_MT,
                       std::array<std::array<double, 3>, 3>& J_MH,
                       std::array<std::array<double, 6>, 6>& J_MM) const;
    
    /**
     * @brief Enable fault thermal pressurization
     * 
     * For fault zones, use the ThermalPressurization model from
     * ThermalPressurization.hpp which solves 1D diffusion equations
     * across the shear zone. This provides more accurate physics than
     * the bulk Lambda coefficient for localized shearing.
     * 
     * @param tp ThermalPressurization instance (will be stored in params)
     */
    void enableFaultThermalPressurization(std::shared_ptr<ThermalPressurization> tp);
    
    /**
     * @brief Get fault thermal pressurization contribution
     * 
     * If fault TP is enabled, returns additional pore pressure change
     * from the detailed 1D diffusion model.
     * 
     * @param state Current TP state at fault
     * @param dt Time step
     * @return Additional delta_p from fault TP
     */
    double getFaultTPContribution(ThermalPressurization::State& state, double dt) const;
    
private:
    Parameters params_;
    
    std::shared_ptr<FluidModelBase> thermal_model_;
    std::shared_ptr<FluidModelBase> hydraulic_model_;
    std::shared_ptr<MaterialModelBase> mechanical_model_;
    
    // Coupling flags
    bool T_to_H_coupled_ = true;
    bool T_to_M_coupled_ = true;
    bool H_to_T_coupled_ = true;
    bool H_to_M_coupled_ = true;
    bool M_to_T_coupled_ = false;
    bool M_to_H_coupled_ = true;
    
    // Fault thermal pressurization integration
    // Uses ThermalPressurization from ThermalPressurization.hpp for detailed
    // fault-zone physics (1D diffusion equations across shear zone)
    bool use_fault_tp_ = false;
    std::shared_ptr<ThermalPressurization> fault_tp_model_ = nullptr;
};


// =============================================================================
// THMC (Thermo-Hydro-Mechanical-Chemical) Coupling
// =============================================================================

/**
 * @brief Chemical reaction type
 */
enum class ReactionType {
    EQUILIBRIUM,         ///< Instantaneous equilibrium
    KINETIC,             ///< Rate-limited kinetics
    SURFACE,             ///< Surface reactions
    PRECIPITATION,       ///< Mineral precipitation
    DISSOLUTION,         ///< Mineral dissolution
    SORPTION             ///< Adsorption/desorption
};

/**
 * @brief Chemical species
 */
struct ChemicalSpecies {
    std::string name;
    double molar_mass;           ///< [kg/mol]
    double diffusion_coefficient; ///< [m²/s]
    double initial_concentration; ///< [mol/m³]
    bool is_primary = true;       ///< Primary or secondary species
};

/**
 * @brief Chemical reaction
 */
struct ChemicalReaction {
    std::string name;
    ReactionType type;
    
    // Stoichiometry (species index → coefficient)
    std::map<int, double> reactants;
    std::map<int, double> products;
    
    // Kinetic parameters
    double rate_constant = 1e-10;     ///< k [mol/(m²·s)]
    double activation_energy = 60e3;  ///< E_a [J/mol]
    double equilibrium_constant = 1.0; ///< K_eq
    
    // Surface area (for mineral reactions)
    double specific_surface_area = 100.0; ///< [m²/m³]
};

/**
 * @brief Chemical effects on other physics
 */
struct ChemicalEffects {
    // C → H (Chemistry to Hydraulic)
    bool porosity_change = true;
    bool permeability_change = true;
    double porosity_reaction_coefficient = 0.1;  // Δφ per mole reacted
    
    // C → M (Chemistry to Mechanical)
    bool weakening = false;
    double strength_reduction_factor = 0.01;     // Per unit concentration
    
    // C → T (Chemistry to Thermal)
    bool reaction_enthalpy = true;
    double heat_of_reaction = 0.0;               // [J/mol]
};

/**
 * @brief THMC fully coupled model
 * 
 * Extends THM with chemical transport and reactions:
 * 
 * Transport equation for species i:
 * ∂(φ·C_i)/∂t + ∇·(v·C_i) - ∇·(D_eff·∇C_i) = R_i
 * 
 * where:
 * - C_i = concentration of species i
 * - v = Darcy velocity
 * - D_eff = effective diffusion/dispersion
 * - R_i = reaction source/sink term
 * 
 * Chemical effects feedback to THM through:
 * - Porosity changes (dissolution/precipitation)
 * - Permeability changes (clogging/enhancement)
 * - Mechanical weakening (chemical softening)
 */
class THMCCoupling : public THMCoupling {
public:
    struct Parameters {
        // Inherit THM parameters
        THMCoupling::Parameters thm_params;
        
        // Chemical parameters
        std::vector<ChemicalSpecies> species;
        std::vector<ChemicalReaction> reactions;
        ChemicalEffects chemical_effects;
        
        // Transport parameters
        double longitudinal_dispersivity = 0.1;   ///< α_L [m]
        double transverse_dispersivity = 0.01;    ///< α_T [m]
        double tortuosity = 2.0;
        
        // Numerical
        bool operator_split_chemistry = true;
        int chemistry_substeps = 10;
    };
    
    THMCCoupling();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Add chemical species
     */
    void addSpecies(const ChemicalSpecies& species);
    
    /**
     * @brief Add chemical reaction
     */
    void addReaction(const ChemicalReaction& reaction);
    
    /**
     * @brief Calculate effective dispersion tensor
     * 
     * D_ij = α_T·|v|·δ_ij + (α_L - α_T)·v_i·v_j/|v| + D_m/τ
     */
    std::array<std::array<double, 3>, 3> dispersionTensor(
        const std::array<double, 3>& velocity) const;
    
    /**
     * @brief Calculate reaction rates for all species
     * 
     * Uses kinetic rate laws or equilibrium constraints.
     */
    std::vector<double> calculateReactionRates(
        const std::vector<double>& concentrations,
        double T, double pH) const;
    
    /**
     * @brief Calculate porosity change from reactions
     * 
     * Δφ = -Σ(V_m · r_m)
     * 
     * where V_m is molar volume and r_m is reaction rate.
     */
    double porosityChange(const std::vector<double>& reaction_rates,
                         double dt) const;
    
    /**
     * @brief Calculate permeability change
     * 
     * Kozeny-Carman: k/k_0 = (φ/φ_0)^3 · ((1-φ_0)/(1-φ))^2
     */
    double permeabilityChange(double phi, double phi_0, double k_0) const;
    
    /**
     * @brief Solve reactive transport step
     * 
     * 1. Transport (advection-dispersion)
     * 2. Chemistry (speciation + kinetics)
     */
    void reactiveTransportStep(double dt,
                              const std::array<double, 3>& velocity,
                              std::vector<std::vector<double>>& concentrations,
                              std::vector<double>& porosity);
    
    /**
     * @brief Solve speciation (equilibrium chemistry)
     * 
     * Newton iteration for equilibrium concentrations.
     */
    void solveSpeciation(std::vector<double>& concentrations,
                        double T, double ionic_strength) const;
    
private:
    Parameters params_thmc_;
    
    std::vector<ChemicalSpecies> species_;
    std::vector<ChemicalReaction> reactions_;
    
    // Equilibrium constants (temperature-dependent)
    double equilibriumConstant(const ChemicalReaction& rxn, double T) const;
    
    // Activity coefficients
    double activityCoefficient(int species_idx, double ionic_strength) const;
};


// =============================================================================
// Unsaturated Flow Coupling (Richards + Biot)
// =============================================================================

/**
 * @brief Unsaturated soil constitutive model
 */
enum class UnsaturatedModel {
    VAN_GENUCHTEN,       ///< Van Genuchten SWRC
    BROOKS_COREY,        ///< Brooks-Corey SWRC
    FREDLUND_XING,       ///< Fredlund-Xing model
    MODIFIED_CAMCLAY     ///< Barcelona Basic Model (BBM)
};

/**
 * @brief Unsaturated poromechanics model
 * 
 * Extends Biot's theory to partially saturated conditions:
 * 
 * Bishop's effective stress:
 * σ' = σ - [S·p_w + (1-S)·p_g]·I = σ - p̄·I
 * 
 * Or with χ parameter:
 * σ' = σ - χ·p_w - (1-χ)·p_g
 * 
 * where χ is often taken as S (saturation).
 * 
 * Mass balance includes:
 * - Water: ∂(φ·S·ρ_w)/∂t + ∇·(ρ_w·q_w) = Q_w
 * - Air: ∂(φ·(1-S)·ρ_a)/∂t + ∇·(ρ_a·q_a) = Q_a
 * 
 * with capillary pressure: p_c = p_g - p_w = f(S)
 */
class UnsaturatedCoupling {
public:
    struct SoilWaterCharacteristics {
        UnsaturatedModel model = UnsaturatedModel::VAN_GENUCHTEN;
        
        // Van Genuchten parameters
        double alpha = 0.01;              ///< α [1/Pa]
        double n = 1.5;                   ///< n [-]
        double m = 0.333;                 ///< m = 1 - 1/n
        double S_residual = 0.1;          ///< S_r [-]
        double S_saturated = 1.0;         ///< S_s [-]
        
        // Brooks-Corey parameters
        double air_entry_pressure = 1e4;  ///< P_b [Pa]
        double lambda = 0.5;              ///< Pore size distribution
        
        // Relative permeability
        double kr_model_exponent = 0.5;   ///< L in Mualem model
    };
    
    struct Parameters {
        SoilWaterCharacteristics swrc;
        
        // Coupling type
        bool use_bishop_stress = true;
        bool include_suction_hardening = false;
        
        // Two-phase flow
        bool include_air_flow = false;    ///< Include air phase explicitly
        double air_compressibility = 1e-5; ///< [1/Pa]
        
        // Hysteresis
        bool include_hysteresis = false;
    };
    
    UnsaturatedCoupling();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate saturation from capillary pressure (SWRC)
     * 
     * Van Genuchten: S = S_r + (S_s - S_r) / [1 + (α·p_c)^n]^m
     */
    double saturation(double capillary_pressure) const;
    
    /**
     * @brief Calculate capillary pressure from saturation
     * 
     * Inverse of SWRC.
     */
    double capillaryPressure(double saturation) const;
    
    /**
     * @brief Calculate specific moisture capacity
     * 
     * C = ∂S/∂p_c
     */
    double moistureCapacity(double capillary_pressure) const;
    
    /**
     * @brief Calculate relative permeability (Mualem model)
     * 
     * k_r = S_e^L · [1 - (1 - S_e^(1/m))^m]^2
     */
    double relativePermeability(double saturation) const;
    
    /**
     * @brief Calculate Bishop's χ parameter
     * 
     * Various formulations:
     * - χ = S (degree of saturation)
     * - χ = S_e (effective saturation)
     * - χ = f(S, p_c) (more general)
     */
    double bishopChi(double saturation, double capillary_pressure) const;
    
    /**
     * @brief Calculate equivalent pore pressure
     * 
     * p̄ = χ·p_w + (1-χ)·p_g
     */
    double equivalentPorePressure(double p_water, double p_gas,
                                  double saturation) const;
    
    /**
     * @brief Calculate effective stress
     */
    std::array<double, 6> effectiveStress(const std::array<double, 6>& total_stress,
                                          double p_water, double p_gas,
                                          double saturation) const;
    
    /**
     * @brief Richards equation residual
     * 
     * ∂θ/∂t - ∇·(K·∇(ψ + z)) = 0
     * 
     * where θ = φ·S, ψ = p_w/ρ_w/g (pressure head)
     */
    double richardsResidual(double theta, double theta_dot,
                           double K, const std::array<double, 3>& grad_psi,
                           const std::array<double, 3>& grad_z) const;
    
    /**
     * @brief BBM (Barcelona Basic Model) yield surface
     * 
     * Extends Cam-Clay for unsaturated soils:
     * f = q² - M²·(p + p_s)·(p_0(s) - p) = 0
     * 
     * where s = suction, p_s = cohesion from suction.
     */
    double bbmYieldSurface(double p, double q, double suction,
                          double p0_sat, double M) const;
    
    /**
     * @brief Loading-collapse (LC) curve
     * 
     * p_0(s) = p_c * [(p_0_star/p_c)^((lambda_0 - kappa)/(lambda_s - kappa))]
     */
    double loadingCollapseCurve(double suction, double p0_star,
                               double p_c, double lambda_0,
                               double lambda_s, double kappa) const;
    
private:
    Parameters params_;
    
    double effectiveSaturation(double S) const;
};


// =============================================================================
// Multi-Scale Coupling
// =============================================================================

/**
 * @brief Multi-scale method type
 */
enum class MultiScaleMethod {
    FE2,                 ///< FE² computational homogenization
    HMM,                 ///< Heterogeneous Multiscale Method
    VARIATIONAL,         ///< Variational multiscale
    ASYMPTOTIC           ///< Asymptotic homogenization
};

/**
 * @brief Multi-scale coupling for heterogeneous media
 * 
 * Bridges microscale (pore/grain scale) to macroscale (continuum):
 * 
 * FE² approach:
 * - Macro: Standard FE with effective properties
 * - Micro: RVE solution at each integration point
 * - Coupling: Homogenization of micro response
 * 
 * Asymptotic homogenization:
 * - Solves cell problems on periodic RVE
 * - Provides effective tensors (stiffness, permeability)
 */
class MultiScaleCoupling {
public:
    struct MicroStructure {
        // RVE dimensions
        double size_x = 0.001;           ///< RVE size X [m]
        double size_y = 0.001;           ///< RVE size Y [m]
        double size_z = 0.001;           ///< RVE size Z [m]
        
        // Micro mesh (or analytical description)
        bool use_analytical = true;
        std::string micro_mesh_file;
        
        // Analytical micro-structure
        double inclusion_fraction = 0.3;  ///< Volume fraction
        double inclusion_radius = 0.0001; ///< For spherical inclusions [m]
        
        // Constituent properties
        double matrix_modulus = 10e9;
        double inclusion_modulus = 50e9;
        double matrix_permeability = 1e-15;
        double inclusion_permeability = 1e-18;
    };
    
    struct Parameters {
        MultiScaleMethod method = MultiScaleMethod::ASYMPTOTIC;
        
        MicroStructure micro;
        
        // FE² parameters
        int micro_elements_per_direction = 10;
        bool use_periodic_bc = true;
        bool use_reduced_integration = false;
        
        // Caching
        bool cache_micro_solutions = true;
        double cache_tolerance = 1e-4;
    };
    
    MultiScaleCoupling();
    
    void setParameters(const Parameters& params);
    void configure(const std::map<std::string, std::string>& config);
    
    /**
     * @brief Calculate effective stiffness tensor
     * 
     * Solves 6 cell problems (3D) or 3 cell problems (2D)
     * for effective elastic moduli.
     */
    std::array<std::array<double, 6>, 6> effectiveStiffness() const;
    
    /**
     * @brief Calculate effective permeability tensor
     * 
     * Solves 3 cell problems for effective permeability.
     */
    std::array<std::array<double, 3>, 3> effectivePermeability() const;
    
    /**
     * @brief Calculate effective thermal conductivity
     */
    std::array<std::array<double, 3>, 3> effectiveThermalConductivity() const;
    
    /**
     * @brief Calculate effective Biot parameters
     * 
     * Homogenization of poroelastic coupling.
     */
    void effectiveBiotParameters(double& alpha_eff, double& M_eff) const;
    
    /**
     * @brief Solve micro problem for given macro strain
     * 
     * FE² approach: Returns homogenized stress from RVE.
     */
    std::array<double, 6> solveMicro(const std::array<double, 6>& macro_strain,
                                     double macro_pressure) const;
    
    /**
     * @brief Get consistent tangent from micro solution
     * 
     * Numerical differentiation or analytical for linear micro.
     */
    std::array<std::array<double, 6>, 6> microTangent(
        const std::array<double, 6>& macro_strain) const;
    
    /**
     * @brief Hashin-Shtrikman bounds for composite
     * 
     * Provides bounds on effective moduli.
     */
    void hashinShtrikmanBounds(double f, double K1, double G1,
                               double K2, double G2,
                               double& K_lower, double& K_upper,
                               double& G_lower, double& G_upper) const;
    
    /**
     * @brief Self-consistent estimate
     * 
     * Eshelby-based effective medium approximation.
     */
    void selfConsistentEstimate(double f, double K1, double G1,
                                double K2, double G2,
                                double& K_eff, double& G_eff) const;
    
private:
    Parameters params_;
    
    // Cached effective properties
    mutable bool properties_computed_ = false;
    mutable std::array<std::array<double, 6>, 6> C_eff_;
    mutable std::array<std::array<double, 3>, 3> k_eff_;
    
    void computeEffectiveProperties() const;
    
    // Cell problem solvers
    std::array<double, 6> solveCellProblem(int load_case) const;
};


// =============================================================================
// Advanced Coupled Physics Configuration
// =============================================================================

/**
 * @brief Configuration for advanced coupled physics
 */
struct AdvancedCouplingConfig {
    // Master enable
    bool enable_advanced_coupling = false;
    
    // Full Biot dynamics
    bool enable_full_biot = false;
    FullBiotDynamics::Parameters biot_params;
    FullBiotDynamics::FluidProperties biot_fluid;
    FullBiotDynamics::SolidFrameProperties biot_solid;
    FullBiotDynamics::PorousMediaProperties biot_porous;
    
    // THM coupling
    bool enable_thm = false;
    THMCoupling::Parameters thm_params;
    
    // THMC coupling
    bool enable_thmc = false;
    THMCCoupling::Parameters thmc_params;
    
    // Unsaturated coupling
    bool enable_unsaturated = false;
    UnsaturatedCoupling::Parameters unsaturated_params;
    
    // Multi-scale
    bool enable_multiscale = false;
    MultiScaleCoupling::Parameters multiscale_params;
    
    // Parse from config
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Validate
    bool validate(std::string& error_msg) const;
};

/**
 * @brief Factory functions for advanced coupling
 */
std::unique_ptr<FullBiotDynamics> createFullBiotModel(
    const std::map<std::string, std::string>& config);

std::unique_ptr<THMCoupling> createTHMCoupling(
    const std::map<std::string, std::string>& config);

std::unique_ptr<THMCCoupling> createTHMCCoupling(
    const std::map<std::string, std::string>& config);

std::unique_ptr<UnsaturatedCoupling> createUnsaturatedCoupling(
    const std::map<std::string, std::string>& config);

std::unique_ptr<MultiScaleCoupling> createMultiScaleCoupling(
    const std::map<std::string, std::string>& config);

} // namespace HighFidelity
} // namespace FSRM

#endif // ADVANCED_COUPLED_PHYSICS_HPP
