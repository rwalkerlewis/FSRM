/**
 * @file ExplosionDamageKernels.hpp
 * @brief Physics kernels for explosion sources, damage evolution, and crater formation
 * 
 * Implements full physics implementations for:
 * - Explosion source models (nuclear, chemical, volcanic)
 * - Near-field damage evolution (cavity, crushed zone, fractured zone)
 * - Hydrodynamic flow (Euler equations with EOS)
 * - Crater formation (excavation flow, ejecta)
 * 
 * All calculations use SI base units (m, kg, s).
 */

#ifndef EXPLOSION_DAMAGE_KERNELS_HPP
#define EXPLOSION_DAMAGE_KERNELS_HPP

#include "PhysicsKernel.hpp"
#include "ExplosionImpactPhysics.hpp"
#include "UnitSystem.hpp"
#include <memory>
#include <functional>
#include <cmath>

namespace FSRM {

// =============================================================================
// Explosion Source Kernel
// =============================================================================

/**
 * @brief Physics kernel for explosion source models
 * 
 * Provides seismic source terms for underground explosions following
 * the Mueller-Murphy reduced displacement potential model.
 * 
 * Governing Equation (Moment Tensor Source):
 *   Mij(t) = M₀ · S(t) · δij
 * 
 * where M₀ is the scalar moment and S(t) is the source time function.
 * 
 * For underground explosions, the source is predominantly isotropic
 * with possible CLVD (compensated linear vector dipole) component
 * from tectonic release.
 */
class ExplosionSourceKernel : public PhysicsKernel {
public:
    ExplosionSourceKernel();
    ~ExplosionSourceKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 6; }  // Moment tensor
    
    // Configuration
    void setExplosionParameters(double yield_kt, double depth_m);
    void setMediumProperties(double rho, double vp, double vs);
    void setSourceLocation(double x, double y, double z);
    void setTectonicRelease(bool enable, double orientation, double fraction);
    
    // Query
    double getScalarMoment() const;
    double getCornerFrequency() const;
    double getMagnitude() const;  // mb or Ms
    
    // Source time function
    double sourceTimeFunction(double t) const;
    void getMomentTensor(double t, double Mij[6]) const;  // Voigt notation
    
    // Generate equivalent point source force
    void getPointForce(const PetscReal x[], double t, PetscScalar f[3]) const;
    
private:
    // Source parameters
    double yield_kt_;
    double depth_;
    double source_x_, source_y_, source_z_;
    
    // Medium properties
    double density_;
    double p_velocity_;
    double s_velocity_;
    
    // Derived quantities
    double scalar_moment_;
    double corner_frequency_;
    double rise_time_;
    double overshoot_;
    
    // Tectonic release
    bool tectonic_release_;
    double tectonic_orientation_;
    double tectonic_fraction_;
    
    // Internal source model
    std::unique_ptr<MuellerMurphySource> source_model_;
    
    void computeDerivedQuantities();
};

// =============================================================================
// Near-Field Damage Kernel
// =============================================================================

/**
 * @brief Physics kernel for near-field damage zone evolution
 * 
 * Models the evolution of damage zones around an explosion cavity:
 * - Cavity zone (D = 1): Complete destruction
 * - Crushed zone: Pervasive fracturing, pore collapse
 * - Fractured zone: Discrete fractures, enhanced permeability
 * 
 * Governing Equations:
 * 
 * Damage evolution:
 *   ∂D/∂t = (1-D) · f(ε, σ, ε̇)   for D < 1
 * 
 * Permeability enhancement:
 *   k = k₀ · (1 + α·D)^n
 * 
 * Porosity change:
 *   φ = φ₀ · (1 + β·D)
 */
class NearFieldDamageKernel : public PhysicsKernel {
public:
    NearFieldDamageKernel();
    ~NearFieldDamageKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 1; }  // Damage D
    
    // Configuration
    void setExplosionSource(double yield_kt, double x, double y, double z);
    void setMaterialProperties(double tensile_strength, double compressive_strength);
    void setDamageParameters(double rate_coeff, double strain_threshold);
    void setPermeabilityModel(double alpha, double n);
    
    // Query damage zones
    double getCavityRadius() const;
    double getCrushedZoneRadius() const;
    double getFracturedZoneRadius() const;
    
    // Permeability from damage
    double getEnhancedPermeability(double D, double k0) const;
    
    // Zone identification
    enum class DamageZone { INTACT, FRACTURED, CRUSHED, CAVITY };
    DamageZone classifyZone(const PetscReal x[]) const;
    DamageZone classifyZoneFromDamage(double D) const;
    
private:
    // Source
    double yield_kt_;
    double source_x_, source_y_, source_z_;
    
    // Material
    double tensile_strength_;
    double compressive_strength_;
    
    // Damage model
    double damage_rate_coeff_;
    double strain_threshold_;
    
    // Permeability model
    double perm_alpha_;
    double perm_exponent_;
    
    // Zone radii (computed)
    double cavity_radius_;
    double crushed_radius_;
    double fractured_radius_;
    
    void computeZoneRadii();
    double damageRate(double D, double strain, double strain_rate) const;
};

// =============================================================================
// Hydrodynamic Kernel (Full Implementation)
// =============================================================================

/**
 * @brief Full implementation of hydrodynamic (Euler) equations
 * 
 * Implements compressible Euler equations for shock propagation:
 * 
 * Conservation form:
 *   ∂U/∂t + ∇·F(U) = S
 * 
 * where U = (ρ, ρv, E) is the conservative state vector.
 * 
 * Supports multiple equations of state:
 * - Ideal gas: p = (γ-1)ρe
 * - Mie-Grüneisen: p = p_H + Γρ(e - e_H)
 * - Tillotson: For impact/vaporization
 * - JWL: For detonation products
 */
class HydrodynamicKernel : public PhysicsKernel {
public:
    /// Equation of state types
    enum class EOSType {
        IDEAL_GAS,
        STIFFENED_GAS,
        MIE_GRUNEISEN,
        TILLOTSON,
        JWL,
        TABULATED
    };
    
    HydrodynamicKernel();
    ~HydrodynamicKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 3; }  // ρ, ρv, E
    int getNumComponents(int field) const override {
        if (field == 0) return 1;       // density
        if (field == 1) return 3;       // momentum
        return 1;                        // energy
    }
    
    // Configuration
    void setEOS(EOSType type);
    void setIdealGasParameters(double gamma);
    void setMieGruneisenParameters(double rho0, double c0, double s, double Gamma0);
    void setTillotsonParameters(double rho0, double A, double B, double E0,
                                double a, double b, double alpha, double beta);
    void setJWLParameters(double A, double B, double R1, double R2, double omega);
    
    // Artificial viscosity
    void setArtificialViscosity(double C_q, double C_l);
    
    // Source term (for explosions)
    void setExplosionSource(double yield_kt, double x, double y, double z, double t0);
    
    // Gravity
    void setGravity(double gx, double gy, double gz);
    
    // EOS evaluation
    double pressure(double rho, double e) const;
    double soundSpeed(double rho, double p) const;
    double temperature(double rho, double e) const;
    
    // Riemann solver flux
    void riemannFlux(const double UL[5], const double UR[5], 
                    const double n[3], double flux[5]) const;
    
private:
    EOSType eos_type_;
    
    // Ideal gas
    double gamma_;
    
    // Mie-Grüneisen
    double rho0_mg_;
    double c0_mg_;
    double s_mg_;
    double Gamma0_mg_;
    
    // Tillotson
    double rho0_t_;
    double A_t_, B_t_, E0_t_;
    double a_t_, b_t_, alpha_t_, beta_t_;
    
    // JWL
    double A_jwl_, B_jwl_, R1_jwl_, R2_jwl_, omega_jwl_;
    
    // Artificial viscosity
    double C_q_;  // Quadratic coefficient
    double C_l_;  // Linear coefficient
    
    // Gravity
    double gravity_[3];
    
    // Explosion source
    bool has_source_;
    double source_yield_kt_;
    double source_x_, source_y_, source_z_;
    double source_t0_;
    
    // Internal methods
    double pressureIdealGas(double rho, double e) const;
    double pressureMieGruneisen(double rho, double e) const;
    double pressureTillotson(double rho, double e) const;
    double pressureJWL(double rho, double e) const;
    
    void HLLCFlux(const double UL[5], const double UR[5], 
                 const double n[3], double flux[5]) const;
};

// =============================================================================
// Crater Formation Kernel
// =============================================================================

/**
 * @brief Physics kernel for impact crater formation
 * 
 * Models the excavation flow and crater formation following impact:
 * - Transient crater formation (excavation)
 * - Crater modification (collapse for large craters)
 * - Ejecta dynamics and deposition
 * 
 * Uses Z-model for excavation flow field:
 *   v_r = A/r² · (1 - z/z_c)^Z
 * 
 * Scaling laws from Holsapple-Schmidt pi-group analysis.
 */
class CraterFormationKernel : public PhysicsKernel {
public:
    CraterFormationKernel();
    ~CraterFormationKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 2; }  // Displacement, Velocity
    int getNumComponents(int field) const override { return 3; }
    
    // Configuration
    void setImpactor(double diameter, double velocity, double density, double angle);
    void setTarget(double density, double strength, double gravity);
    void setZModelParameters(double Z, double flow_depth);
    
    // Query results
    double getTransientCraterDiameter() const;
    double getTransientCraterDepth() const;
    double getFinalCraterDiameter() const;
    double getFinalCraterDepth() const;
    double getCraterFormationTime() const;
    bool isComplexCrater() const;
    
    // Excavation flow field
    void getExcavationVelocity(const PetscReal x[], double t, 
                               PetscScalar v[3]) const;
    
    // Ejecta
    double getEjectaThickness(double r) const;
    void getEjectaVelocity(double r, double& v, double& angle) const;
    
    // Crater profile at time t
    double getCraterDepth(double r, double t) const;
    
private:
    // Impactor
    double impactor_diameter_;
    double impactor_velocity_;
    double impactor_density_;
    double impact_angle_;
    
    // Target
    double target_density_;
    double target_strength_;
    double gravity_;
    
    // Z-model
    double Z_exponent_;
    double flow_center_depth_;
    
    // Computed crater properties
    double transient_diameter_;
    double transient_depth_;
    double final_diameter_;
    double final_depth_;
    double formation_time_;
    
    // Scaling model
    std::unique_ptr<CraterScalingModel> scaling_model_;
    std::unique_ptr<ZModelExcavation> z_model_;
    
    void computeCraterDimensions();
};

// =============================================================================
// Chemical Reaction Kernel
// =============================================================================

/**
 * @brief Physics kernel for reactive transport
 * 
 * Implements advection-diffusion-reaction equations:
 *   φ ∂Cᵢ/∂t + ∇·(vCᵢ) - ∇·(D∇Cᵢ) = Rᵢ(C₁,...,Cₙ)
 * 
 * Supports:
 * - Homogeneous reactions (aqueous phase)
 * - Heterogeneous reactions (mineral dissolution/precipitation)
 * - Equilibrium and kinetic reactions
 * - Temperature-dependent rate laws
 */
class ChemicalReactionKernel : public PhysicsKernel {
public:
    /// Reaction types
    enum class ReactionType {
        EQUILIBRIUM,       // Fast equilibrium
        KINETIC,           // Rate-limited kinetic
        TST,               // Transition state theory
        ARRHENIUS          // Temperature-dependent Arrhenius
    };
    
    ChemicalReactionKernel(int num_species);
    ~ChemicalReactionKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return n_species_; }
    
    // Species configuration
    void addSpecies(const std::string& name, double diffusivity, double charge);
    
    // Reaction configuration
    void addReaction(const std::string& name, ReactionType type,
                    const std::vector<int>& reactants,
                    const std::vector<int>& products,
                    const std::vector<double>& stoichiometry);
    void setReactionRate(int reaction_id, double k_forward, double k_backward);
    void setActivationEnergy(int reaction_id, double Ea);
    
    // Mineral reactions
    void addMineralReaction(const std::string& mineral,
                           const std::vector<int>& species,
                           const std::vector<double>& stoichiometry,
                           double k_rate, double surface_area);
    
    // Porosity-permeability coupling
    void enablePorosityFeedback(bool enable);
    double computePorosityChange(const PetscScalar C[]) const;
    double computePermeabilityChange(double dphi) const;
    
    // Query
    double getReactionRate(int reaction_id, const PetscScalar C[], 
                          double T = 298.15) const;
    
private:
    int n_species_;
    int n_reactions_;
    
    // Species properties
    std::vector<std::string> species_names_;
    std::vector<double> diffusivities_;
    std::vector<double> charges_;
    
    // Reaction definitions
    struct Reaction {
        std::string name;
        ReactionType type;
        std::vector<int> reactant_ids;
        std::vector<int> product_ids;
        std::vector<double> stoichiometry;
        double k_forward;
        double k_backward;
        double activation_energy;
        bool is_mineral;
        double mineral_surface_area;
    };
    std::vector<Reaction> reactions_;
    
    // Porosity feedback
    bool porosity_feedback_;
    double initial_porosity_;
    
    // Internal methods
    double equilibriumReactionRate(const Reaction& rxn, const PetscScalar C[]) const;
    double kineticReactionRate(const Reaction& rxn, const PetscScalar C[], double T) const;
    double mineralReactionRate(const Reaction& rxn, const PetscScalar C[], double T) const;
};

// =============================================================================
// Thermal Radiation Kernel
// =============================================================================

/**
 * @brief Physics kernel for thermal radiation transport
 * 
 * Implements thermal radiation from explosions/fireballs:
 *   (1/c) ∂I/∂t + n̂·∇I + κI = j
 * 
 * or simplified thermal pulse model:
 *   Q(r, t) = Q₀ · τ(r) · f(t)
 * 
 * where τ is atmospheric transmittance and f(t) is time profile.
 */
class ThermalRadiationKernel : public PhysicsKernel {
public:
    ThermalRadiationKernel();
    ~ThermalRadiationKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 1; }  // Thermal fluence
    
    // Configuration
    void setFireball(double yield_kt, double burst_height);
    void setAtmosphericTransmittance(std::function<double(double)> tau);
    
    // Query
    double getThermalFluence(double r, double t) const;
    double getPeakFluence(double r) const;
    double getThermalPulseWidth() const;
    
private:
    double yield_kt_;
    double burst_height_;
    std::function<double(double)> transmittance_;
    
    std::unique_ptr<FireballModel> fireball_;
};

// =============================================================================
// EMP Kernel
// =============================================================================

/**
 * @brief Physics kernel for electromagnetic pulse
 * 
 * Implements E1, E2, E3 EMP components from high-altitude bursts:
 * 
 * Maxwell's equations with source current J_c (Compton):
 *   ∂E/∂t = c²∇×B - J_c/ε₀
 *   ∂B/∂t = -∇×E
 */
class EMPKernel : public PhysicsKernel {
public:
    EMPKernel();
    ~EMPKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 2; }  // E, B
    int getNumComponents(int field) const override { return 3; }
    
    // Configuration
    void setBurst(double yield_kt, double altitude, double lat, double lon);
    void setGeomagneticField(double B0, double dip_angle);
    
    // Query
    void getElectricField(const PetscReal x[], double t, PetscScalar E[3]) const;
    void getMagneticField(const PetscReal x[], double t, PetscScalar B[3]) const;
    double getE1Peak() const;
    double getE2Peak() const;
    double getE3Peak() const;
    
private:
    double yield_kt_;
    double altitude_;
    double latitude_, longitude_;
    double geomagnetic_field_;
    double geomagnetic_dip_;
    
    std::unique_ptr<EMPModel> emp_model_;
};

// =============================================================================
// Fallout Kernel
// =============================================================================

/**
 * @brief Physics kernel for radioactive fallout transport
 * 
 * Advection-diffusion with radioactive decay:
 *   ∂A/∂t + ∇·(vA) - ∇·(K∇A) + v_s ∂A/∂z = -λA + S
 * 
 * where A is activity, v is wind, v_s is settling velocity,
 * λ is decay constant, S is source.
 */
class FalloutKernel : public PhysicsKernel {
public:
    FalloutKernel();
    ~FalloutKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 1; }  // Activity
    
    // Configuration
    void setBurst(double yield_kt, double fission_yield_kt, 
                 double x, double y, double height);
    void setWind(std::function<void(double z, double& vx, double& vy)> wind);
    void setDiffusivity(double Kh, double Kv);
    
    // Query
    double getGroundDeposition(double x, double y, double t) const;
    double getDoseRate(double x, double y, double t) const;
    
private:
    double yield_kt_;
    double fission_yield_kt_;
    double source_x_, source_y_, source_height_;
    
    std::function<void(double, double&, double&)> wind_profile_;
    double K_horizontal_;
    double K_vertical_;
    
    std::unique_ptr<FalloutModel> fallout_model_;
};

// =============================================================================
// Tsunami Kernel (Full Implementation)
// =============================================================================

/**
 * @brief Full implementation of tsunami (shallow water) physics
 * 
 * Nonlinear shallow water equations with:
 * - Variable bathymetry
 * - Wetting/drying (inundation)
 * - Bottom friction (Manning/Chezy)
 * - Coriolis force
 * - Seafloor deformation source
 */
class TsunamiKernel : public PhysicsKernel {
public:
    TsunamiKernel();
    ~TsunamiKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 2; }  // η, (hu, hv)
    int getNumComponents(int field) const override { 
        return field == 0 ? 1 : 2; 
    }
    
    // Configuration
    void setBathymetry(std::function<double(double, double)> depth);
    void setManning(double n);
    void setCoriolis(double latitude);
    void setMinDepth(double h_min);
    
    // Seafloor deformation source
    void setSeafloorDeformation(std::function<double(double, double, double)> dz);
    void setOkadaSource(double lat, double lon, double depth, 
                       double strike, double dip, double rake,
                       double length, double width, double slip);
    
    // Query
    double getWaterElevation(double x, double y) const;
    double getWaveHeight(double x, double y) const;
    double getFlowSpeed(double x, double y) const;
    bool isInundated(double x, double y) const;
    
private:
    std::function<double(double, double)> bathymetry_;
    std::function<double(double, double, double)> seafloor_deformation_;
    double manning_n_;
    double gravity_;
    double coriolis_f_;
    double min_depth_;
    
    // Okada source parameters
    bool use_okada_;
    double okada_params_[10];
    
    // Helper methods
    double computeFlux(double hL, double hR, double uL, double uR, double g) const;
    void wettingDrying(double& h, double& hu, double& hv) const;
};

// =============================================================================
// Surface Deformation Kernel
// =============================================================================

/**
 * @brief Physics kernel for ground surface deformation
 * 
 * Tracks vertical surface displacement from subsurface processes:
 * - Pore pressure changes (poroelastic)
 * - Fault slip
 * - Magma intrusion
 * - Fluid extraction/injection
 * 
 * Provides:
 * - LOS (line-of-sight) displacement for InSAR
 * - Tilt and strain at surface
 */
class SurfaceDeformationKernel : public PhysicsKernel {
public:
    SurfaceDeformationKernel();
    ~SurfaceDeformationKernel() override = default;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 3; }  // uz, tilt_x, tilt_y
    
    // Configuration
    void setSubsurfaceSource(std::function<void(double, double, double, 
                                                double&, double&, double&)> src);
    void setInSARParameters(double incidence, double azimuth);
    
    // Query
    double getVerticalDisplacement(double x, double y) const;
    double getLOSDisplacement(double x, double y) const;
    void getTilt(double x, double y, double& tilt_x, double& tilt_y) const;
    
    // Okada (Mogi, etc.) analytical sources
    void addOkadaFault(double lat, double lon, double depth,
                      double strike, double dip, double rake,
                      double length, double width, double slip);
    void addMogiSource(double x, double y, double depth, double volume_change);
    
private:
    std::function<void(double, double, double, double&, double&, double&)> 
        subsurface_source_;
    
    double insar_incidence_;
    double insar_azimuth_;
    
    // Analytical sources
    struct MogiSource { double x, y, depth, dV; };
    struct OkadaFault { double params[10]; };
    std::vector<MogiSource> mogi_sources_;
    std::vector<OkadaFault> okada_faults_;
    
    void computeOkadaDisplacement(const OkadaFault& fault, 
                                  double x, double y,
                                  double& ux, double& uy, double& uz) const;
    void computeMogiDisplacement(const MogiSource& src,
                                double x, double y,
                                double& ux, double& uy, double& uz) const;
};

// =============================================================================
// Factory Functions
// =============================================================================

/**
 * @brief Create physics kernel by type
 */
std::shared_ptr<PhysicsKernel> createPhysicsKernel(PhysicsType type);

/**
 * @brief Check if kernel type is implemented
 */
bool isKernelImplemented(PhysicsType type);

/**
 * @brief Get all implemented kernel types
 */
std::vector<PhysicsType> getImplementedKernels();

} // namespace FSRM

#endif // EXPLOSION_DAMAGE_KERNELS_HPP
