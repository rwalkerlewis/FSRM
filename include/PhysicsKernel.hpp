#ifndef PHYSICS_KERNEL_HPP
#define PHYSICS_KERNEL_HPP

#include "FSRM.hpp"
#include <petscfe.h>
#include <functional>

namespace FSRM {

// =============================================================================
// Execution Policy System
// =============================================================================
// 
// The execution policy system allows physics kernels to run on different 
// hardware (CPU/GPU) without code changes. Each kernel can be tagged with
// execution traits that determine how it should be dispatched.

/**
 * @brief Execution backend enumeration
 * 
 * Specifies where kernel computation should occur.
 */
enum class ExecutionBackend {
    CPU,           ///< Standard CPU execution
    CUDA,          ///< NVIDIA CUDA GPU execution
    HIP,           ///< AMD ROCm/HIP GPU execution
    SYCL,          ///< SYCL (Intel/generic GPU) execution
    AUTOMATIC      ///< Auto-select based on availability
};

/**
 * @brief Execution traits for physics kernels
 * 
 * Defines capabilities and requirements for kernel execution.
 */
struct ExecutionTraits {
    bool supports_gpu = false;            ///< Kernel has GPU implementation
    bool requires_gpu = false;            ///< Kernel must run on GPU
    bool supports_batched = true;         ///< Can process multiple elements at once
    bool thread_safe = true;              ///< Safe for OpenMP parallelism
    size_t min_elements_for_gpu = 1000;   ///< Minimum elements to benefit from GPU
    size_t preferred_block_size = 256;    ///< Preferred GPU block size
    
    /// Default CPU-only traits
    static ExecutionTraits CPUOnly() {
        return ExecutionTraits{false, false, true, true, 0, 256};
    }
    
    /// GPU-accelerated traits
    static ExecutionTraits GPUAccelerated(size_t min_elements = 1000) {
        return ExecutionTraits{true, false, true, true, min_elements, 256};
    }
};

/**
 * @brief Kernel capability flags
 * 
 * Bit flags indicating what a kernel can compute.
 */
enum class KernelCapability : uint32_t {
    NONE              = 0,
    RESIDUAL          = 1 << 0,    ///< Can compute residual F(u)
    JACOBIAN          = 1 << 1,    ///< Can compute Jacobian dF/du
    BOUNDARY_INTEGRAL = 1 << 2,    ///< Can compute boundary terms
    POINT_SOURCE      = 1 << 3,    ///< Can add point source terms
    GPU_RESIDUAL      = 1 << 4,    ///< GPU-accelerated residual
    GPU_JACOBIAN      = 1 << 5,    ///< GPU-accelerated Jacobian
    GPU_MATRIX_FREE   = 1 << 6,    ///< GPU matrix-free Jacobian-vector product
    ALL_CPU           = RESIDUAL | JACOBIAN | BOUNDARY_INTEGRAL | POINT_SOURCE,
    ALL_GPU           = GPU_RESIDUAL | GPU_JACOBIAN | GPU_MATRIX_FREE
};

inline KernelCapability operator|(KernelCapability a, KernelCapability b) {
    return static_cast<KernelCapability>(static_cast<uint32_t>(a) | static_cast<uint32_t>(b));
}

inline bool hasCapability(KernelCapability caps, KernelCapability test) {
    return (static_cast<uint32_t>(caps) & static_cast<uint32_t>(test)) != 0;
}

// =============================================================================
// Physics Kernel Base Class
// =============================================================================

/**
 * @brief Base class for all physics kernels
 * 
 * PhysicsKernel provides the abstract interface for finite element kernel 
 * evaluation. All physics implementations (flow, mechanics, thermal, etc.)
 * derive from this class and implement the residual and Jacobian methods.
 * 
 * ## Execution Model
 * 
 * Kernels are evaluated at quadrature points during finite element assembly.
 * The PETSc framework calls residual() and jacobian() for each quadrature 
 * point, providing local solution values and gradients.
 * 
 * ## GPU Support
 * 
 * GPU-accelerated kernels inherit from both PhysicsKernel and provide
 * GPU-specific implementations. The execution policy system handles
 * dispatch between CPU and GPU backends.
 * 
 * ## Thread Safety
 * 
 * Kernel methods must be thread-safe for OpenMP parallelization.
 * Instance state should be read-only during residual/jacobian evaluation.
 */
class PhysicsKernel {
public:
    /**
     * @brief Construct kernel with specified physics type
     * @param type The physics type (flow, mechanics, etc.)
     */
    explicit PhysicsKernel(PhysicsType type) 
        : physics_type(type), 
          execution_traits(ExecutionTraits::CPUOnly()),
          capabilities(KernelCapability::ALL_CPU) {}
    
    virtual ~PhysicsKernel() = default;
    
    // =========================================================================
    // Core Interface (must be implemented by derived classes)
    // =========================================================================
    
    /**
     * @brief Setup kernel with PETSc finite element infrastructure
     * @param dm PETSc distributed mesh
     * @param fe PETSc finite element space
     * @return PETSc error code
     */
    virtual PetscErrorCode setup(DM dm, PetscFE fe) = 0;
    
    /**
     * @brief Compute residual contribution at a quadrature point
     * 
     * Evaluates F(u) for the nonlinear system F(u) = 0.
     * 
     * @param u Solution values at quadrature point
     * @param u_t Time derivatives of solution (for transient problems)
     * @param u_x Spatial gradients of solution (∇u)
     * @param a Auxiliary data (shift parameters from time integrator)
     * @param x Physical coordinates of quadrature point
     * @param[out] f Residual contribution (output)
     * 
     * @note For transient problems, the shift parameter a[0] relates the
     *       time derivative to the solution: du/dt ≈ a[0] * (u - u_old)
     */
    virtual void residual(const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscScalar a[],
                         const PetscReal x[], PetscScalar f[]) = 0;
    
    /**
     * @brief Compute Jacobian contribution at a quadrature point
     * 
     * Evaluates dF/du for Newton iteration.
     * 
     * @param u Solution values at quadrature point
     * @param u_t Time derivatives of solution
     * @param u_x Spatial gradients of solution
     * @param a Auxiliary data (shift = a[0])
     * @param x Physical coordinates
     * @param[out] J Jacobian contribution (output)
     * 
     * @note The Jacobian combines:
     *       - g0: ∂f/∂u (point terms)
     *       - g1: ∂f/∂(∇u) (flux gradient)
     *       - g2: ∂(flux)/∂u (flux solution)
     *       - g3: ∂(flux)/∂(∇u) (flux gradient-gradient)
     * 
     * @note For transient: J += a[0] * ∂f/∂(∂u/∂t) (mass matrix contribution)
     */
    virtual void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscScalar a[],
                         const PetscReal x[], PetscScalar J[]) = 0;
    
    /**
     * @brief Get number of solution fields
     * @return Number of fields (e.g., 1 for pressure, 2 for u+p, etc.)
     */
    virtual int getNumFields() const = 0;
    
    /**
     * @brief Get number of components for a field
     * @param field Field index (0-based)
     * @return Number of components (1 for scalar, 3 for vector in 3D)
     */
    virtual int getNumComponents(int field) const = 0;
    
    // =========================================================================
    // Query Methods
    // =========================================================================
    
    /// Get physics type
    PhysicsType getType() const { return physics_type; }
    
    /// Get physics type name as string
    virtual const char* getTypeName() const {
        switch (physics_type) {
            // Core reservoir
            case PhysicsType::FLUID_FLOW: return "FluidFlow";
            case PhysicsType::GEOMECHANICS: return "Geomechanics";
            case PhysicsType::THERMAL: return "Thermal";
            case PhysicsType::PARTICLE_TRANSPORT: return "ParticleTransport";
            case PhysicsType::FRACTURE_PROPAGATION: return "FracturePropagation";
            case PhysicsType::TIDAL_FORCES: return "TidalForces";
            case PhysicsType::CHEMICAL_REACTION: return "ChemicalReaction";
            // Dynamic waves
            case PhysicsType::ELASTODYNAMICS: return "Elastodynamics";
            case PhysicsType::POROELASTODYNAMICS: return "Poroelastodynamics";
            // Explosion/impact
            case PhysicsType::EXPLOSION_SOURCE: return "ExplosionSource";
            case PhysicsType::NEAR_FIELD_DAMAGE: return "NearFieldDamage";
            case PhysicsType::HYDRODYNAMIC: return "Hydrodynamic";
            case PhysicsType::CRATER_FORMATION: return "CraterFormation";
            // Atmospheric
            case PhysicsType::ATMOSPHERIC_BLAST: return "AtmosphericBlast";
            case PhysicsType::ATMOSPHERIC_ACOUSTIC: return "AtmosphericAcoustic";
            case PhysicsType::INFRASOUND: return "Infrasound";
            case PhysicsType::THERMAL_RADIATION: return "ThermalRadiation";
            case PhysicsType::EMP: return "ElectromagneticPulse";
            case PhysicsType::FALLOUT: return "Fallout";
            // Surface/water
            case PhysicsType::TSUNAMI: return "Tsunami";
            case PhysicsType::SURFACE_DEFORMATION: return "SurfaceDeformation";
            // Hydrodynamics and Ocean Physics
            case PhysicsType::OCEAN_CIRCULATION: return "OceanCirculation";
            case PhysicsType::COASTAL_HYDRODYNAMICS: return "CoastalHydrodynamics";
            case PhysicsType::ESTUARINE_DYNAMICS: return "EstuarineDynamics";
            case PhysicsType::STORM_SURGE: return "StormSurge";
            case PhysicsType::OCEAN_ACOUSTICS: return "OceanAcoustics";
            case PhysicsType::SEDIMENT_TRANSPORT: return "SedimentTransport";
            case PhysicsType::WAVE_DYNAMICS: return "WaveDynamics";
            case PhysicsType::INTERNAL_WAVES: return "InternalWaves";
            case PhysicsType::OCEAN_MIXING: return "OceanMixing";
            case PhysicsType::THERMOHALINE: return "Thermohaline";
            default: return "Unknown";
        }
    }
    
    /// Get total degrees of freedom per point
    virtual int getTotalDOF() const {
        int total = 0;
        for (int f = 0; f < getNumFields(); ++f) {
            total += getNumComponents(f);
        }
        return total;
    }
    
    // =========================================================================
    // Execution Policy Methods
    // =========================================================================
    
    /// Get execution traits
    const ExecutionTraits& getExecutionTraits() const { return execution_traits; }
    
    /// Check if kernel supports GPU execution
    bool supportsGPU() const { return execution_traits.supports_gpu; }
    
    /// Get kernel capabilities
    KernelCapability getCapabilities() const { return capabilities; }
    
    /// Check if kernel has a specific capability
    bool hasCapability(KernelCapability cap) const {
        return FSRM::hasCapability(capabilities, cap);
    }
    
    /// Get recommended minimum problem size for GPU benefit
    size_t getMinGPUElements() const { return execution_traits.min_elements_for_gpu; }
    
    /// Determine best execution backend for given problem size
    virtual ExecutionBackend selectBackend(size_t n_elements) const {
        if (!execution_traits.supports_gpu) {
            return ExecutionBackend::CPU;
        }
        if (execution_traits.requires_gpu) {
            return ExecutionBackend::CUDA;  // Or detect available
        }
        if (n_elements >= execution_traits.min_elements_for_gpu) {
            return ExecutionBackend::AUTOMATIC;  // Let runtime decide
        }
        return ExecutionBackend::CPU;
    }
    
protected:
    PhysicsType physics_type;
    ExecutionTraits execution_traits;
    KernelCapability capabilities;
    
    /// Set execution traits (for derived classes)
    void setExecutionTraits(const ExecutionTraits& traits) {
        execution_traits = traits;
    }
    
    /// Set capabilities (for derived classes)
    void setCapabilities(KernelCapability caps) {
        capabilities = caps;
    }
};

// Single phase flow kernel (Darcy flow)
class SinglePhaseFlowKernel : public PhysicsKernel {
public:
    SinglePhaseFlowKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 1; }
    
    void setProperties(double phi, double k, double ct, double mu, double rho);
    
private:
    double porosity;
    double permeability;
    double compressibility;
    double viscosity;
    double density;
};

// Black oil model kernel (oil, gas, water with dissolution)
class BlackOilKernel : public PhysicsKernel {
public:
    BlackOilKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 3; } // Po, Sw, Sg
    int getNumComponents(int field) const override { return 1; }
    
    void setFluidProperties(const FluidProperties& props);
    void setRockProperties(double phi, double kx, double ky, double kz);
    
    // PVT functions
    double oilDensity(double P, double Rs) const;
    double gasDensity(double P) const;
    double waterDensity(double P) const;
    double oilViscosity(double P, double Rs) const;
    double gasViscosity(double P) const;
    double waterViscosity(double P) const;
    double solutionGOR(double P) const;
    
    // Relative permeability
    double krw(double Sw) const;
    double kro(double So, double Sg) const;
    double krg(double Sg) const;
    
    // Capillary pressure
    double Pcow(double Sw) const;
    double Pcog(double Sg) const;
    
private:
    FluidProperties fluid_props;
    double porosity;
    double perm_x, perm_y, perm_z;
};

// Compositional flow kernel (multi-component with EOS)
class CompositionalKernel : public PhysicsKernel {
public:
    CompositionalKernel(int num_components);
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return nc + 1; } // P + nc-1 compositions
    int getNumComponents(int field) const override { return 1; }
    
    void setComponentProperties(const std::vector<double>& mw,
                               const std::vector<double>& Tc,
                               const std::vector<double>& Pc,
                               const std::vector<double>& omega);
    
    // Equation of state (Peng-Robinson)
    void flashCalculation(double P, double T, const std::vector<double>& z,
                         std::vector<double>& x, std::vector<double>& y,
                         double& S_L, double& S_V) const;
    
    double fugacityCoefficient(double P, double T, double z_factor,
                              const std::vector<double>& comp, int i) const;
    
private:
    int nc; // number of components
    std::vector<double> component_mw;
    std::vector<double> component_Tc;
    std::vector<double> component_Pc;
    std::vector<double> component_omega;
    double porosity;
    double perm_x, perm_y, perm_z;
};

// Geomechanics kernel (linear elasticity or viscoelasticity)
class GeomechanicsKernel : public PhysicsKernel {
public:
    GeomechanicsKernel(SolidModelType model_type);
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 3; } // ux, uy, uz
    
    void setMaterialProperties(double E, double nu, double rho);
    void setViscoelasticProperties(double tau, double eta);
    void setPoroelasticCoupling(double alpha, double M);
    
private:
    SolidModelType model_type;
    double youngs_modulus;
    double poisson_ratio;
    double density;
    double relaxation_time;  // for viscoelasticity
    double viscosity;         // for viscoelasticity
    double biot_coefficient;
    double biot_modulus;
};

// Thermal kernel (heat conduction and convection)
class ThermalKernel : public PhysicsKernel {
public:
    ThermalKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; } // Temperature
    int getNumComponents(int field) const override { return 1; }
    
    void setThermalProperties(double k, double rho, double cp);
    void setFluidThermalProperties(double k_f, double rho_f, double cp_f);
    
private:
    double thermal_conductivity_solid;
    double density_solid;
    double heat_capacity_solid;
    double thermal_conductivity_fluid;
    double density_fluid;
    double heat_capacity_fluid;
    double porosity;
};

// Particle transport kernel (proppant, tracers, etc.)
class ParticleTransportKernel : public PhysicsKernel {
public:
    ParticleTransportKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; } // Concentration
    int getNumComponents(int field) const override { return 1; }
    
    void setParticleProperties(double diameter, double density, double diffusivity);
    void enableGravitationalSettling(bool enable);
    void enableBridging(bool enable);
    void setDispersivity(double alpha_L, double alpha_T);
    void setMediumProperties(double phi, double k);
    
private:
    double particle_diameter;
    double particle_density;
    double diffusivity;
    bool gravity_settling;
    bool enable_bridging;
    double porosity;
    double permeability;
    double longitudinal_dispersivity;   // Longitudinal dispersivity [m]
    double transverse_dispersivity;     // Transverse dispersivity [m]
};

// Fracture propagation kernel (cohesive zone model)
class FracturePropagationKernel : public PhysicsKernel {
public:
    FracturePropagationKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; } // Fracture width
    int getNumComponents(int field) const override { return 1; }
    
    void setFractureProperties(double Kc, double Gc, double sigma_c);
    
private:
    double fracture_toughness;
    double fracture_energy;
    double critical_stress;
};

// Tidal forces kernel
class TidalForcesKernel : public PhysicsKernel {
public:
    TidalForcesKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 0; } // Body force addition
    int getNumComponents(int field) const override { return 0; }
    
    void setLocationAndTime(double lat, double lon, double time);
    void computeTidalStress(const PetscReal x[], PetscScalar stress[]);
    
private:
    double latitude;
    double longitude;
    double current_time;
};

// Elastodynamics kernel (wave propagation with inertia)
class ElastodynamicsKernel : public PhysicsKernel {
public:
    ElastodynamicsKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 3; } // ux, uy, uz
    
    void setMaterialProperties(double E, double nu, double rho);
    void setWaveProperties(double vp, double vs, double Q);
    void setDamping(double alpha, double beta);
    void setStaticTriggeringMode(bool enable, double threshold, double duration);
    
    // Check if static stress exceeds threshold
    bool checkStaticTrigger(const PetscScalar stress[], double current_time);
    bool isInDynamicEvent(double current_time) const;
    
private:
    double youngs_modulus;
    double poisson_ratio;
    double density;
    double p_wave_velocity;
    double s_wave_velocity;
    double quality_factor;
    double damping_alpha;  // Rayleigh mass damping
    double damping_beta;   // Rayleigh stiffness damping
    
    // Static triggering parameters
    bool enable_static_triggering;
    double trigger_threshold;
    double event_duration;
    double trigger_time;
    bool event_active;
};

// Poroelastodynamics kernel (Biot's equations with dynamics)
class PoroelastodynamicsKernel : public PhysicsKernel {
public:
    PoroelastodynamicsKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 2; } // Displacement + Pressure
    int getNumComponents(int field) const override { return (field == 0 ? 3 : 1); }
    
    void setMaterialProperties(double E, double nu, double rho, double phi);
    void setFluidProperties(double rho_f, double mu, double K_f);
    void setBiotParameters(double alpha, double M);
    void setWaveProperties(double vp_fast, double vs, double vp_slow);
    void setDamping(double alpha, double beta);
    void setStaticTriggeringMode(bool enable, double threshold, double duration);
    
    // Permeability dynamics
    void enableDynamicPermeabilityChange(bool enable, double strain_coeff, 
                                         double stress_coeff, double recovery_time);
    double computePermeabilityChange(const PetscScalar u_x[], const PetscScalar stress[], 
                                     double k_initial, double dt);
    
    bool checkStaticTrigger(const PetscScalar stress[], double current_time);
    bool isInDynamicEvent(double current_time) const;
    
private:
    // Solid properties
    double youngs_modulus;
    double poisson_ratio;
    double density_solid;
    double porosity;
    
    // Fluid properties
    double density_fluid;
    double viscosity;
    double bulk_modulus_fluid;
    
    // Biot parameters
    double biot_coefficient;
    double biot_modulus;
    
    // Wave properties
    double p_wave_fast;   // Fast P-wave (Biot wave I)
    double s_wave;        // Shear wave
    double p_wave_slow;   // Slow P-wave (Biot wave II)
    double damping_alpha;
    double damping_beta;
    
    // Permeability
    double permeability;
    double permeability_initial;
    
    // Permeability dynamics
    bool enable_dynamic_k;
    double k_strain_coeff;
    double k_stress_coeff;
    double k_recovery_time;
    
    // Static triggering
    bool enable_static_triggering;
    double trigger_threshold;
    double event_duration;
    double trigger_time;
    bool event_active;
};

// Dynamic permeability change model (can be used by other kernels)
class DynamicPermeabilityModel {
public:
    DynamicPermeabilityModel();
    
    void setParameters(double k0, double strain_coeff, double stress_coeff,
                      double recovery_time, double k_min, double k_max);
    
    // Compute instantaneous permeability change from wave passage
    double computeInstantaneousChange(double strain_amplitude, double stress_amplitude);
    
    // Compute time-dependent permeability evolution
    double computeTimeEvolution(double k_current, double k0, double dt);
    
    // Mechanistic models for permeability change
    double strainModel(double volumetric_strain);
    double stressModel(double effective_stress);
    double shearDilationModel(double shear_strain);
    double crackOpeningModel(double normal_stress);
    
private:
    double k_initial;
    double strain_coefficient;
    double stress_coefficient;
    double recovery_tau;
    double k_minimum;
    double k_maximum;
};

} // namespace FSRM

#endif // PHYSICS_KERNEL_HPP
