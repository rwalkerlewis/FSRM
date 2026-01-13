/**
 * @file PhysicsKernel_GPU_Extended.hpp
 * @brief Extended GPU implementations for physics kernels
 * 
 * This file provides GPU-accelerated implementations for:
 * - Black Oil flow
 * - Thermal transport
 * - Geomechanics (static and dynamic)
 * - Hydrodynamic (Euler equations)
 * - Atmospheric blast
 * - Infrasound propagation
 * - Tsunami (shallow water)
 * 
 * All kernels inherit from the CPU base class and use CRTP pattern
 * for efficient GPU dispatch.
 */

#ifndef PHYSICS_KERNEL_GPU_EXTENDED_HPP
#define PHYSICS_KERNEL_GPU_EXTENDED_HPP

#include "PhysicsKernel.hpp"
#include "PhysicsKernel_GPU.hpp"
#include "AtmosphericInfrasound.hpp"
#include "ExplosionImpactPhysics.hpp"

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#ifdef USE_HIP
#include <hip/hip_runtime.h>
#endif

namespace FSRM {

// =============================================================================
// Black Oil GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated Black Oil kernel
 * 
 * Implements three-phase flow (oil, water, gas) with:
 * - PVT property calculations on GPU
 * - Relative permeability evaluation
 * - Capillary pressure
 * - Solution gas/oil ratio
 * 
 * ## GPU Operations:
 * - Phase property evaluation (Bo, Bg, Bw, Rs, viscosities)
 * - Relative permeability (Corey/tabulated)
 * - Phase mobility computation
 * - Accumulation and flux terms
 * 
 * ## Performance:
 * - 3-8× speedup over CPU for >50K cells
 * - Benefits from coalesced memory access in PVT lookups
 */
class BlackOilKernelGPU : public GPUKernelBase<BlackOilKernelGPU, BlackOilKernel> {
public:
    BlackOilKernelGPU();
    ~BlackOilKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    // GPU memory management
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computePVTPropertiesGPU();
    void computeRelPermGPU();
    void computeMobilitiesGPU();
    void computeAccumulationGPU(double dt);
    void computeFluxGPU();
    
private:
#ifdef USE_CUDA
    // Solution fields
    double* d_pressure_oil;
    double* d_saturation_water;
    double* d_saturation_gas;
    
    // PVT properties
    double* d_Bo;              // Oil formation volume factor
    double* d_Bg;              // Gas formation volume factor
    double* d_Bw;              // Water formation volume factor
    double* d_Rs;              // Solution gas-oil ratio
    double* d_viscosity_oil;
    double* d_viscosity_gas;
    double* d_viscosity_water;
    double* d_density_oil;
    double* d_density_gas;
    double* d_density_water;
    
    // Relative permeability
    double* d_krw;
    double* d_kro;
    double* d_krg;
    
    // Mobilities
    double* d_lambda_oil;
    double* d_lambda_gas;
    double* d_lambda_water;
    double* d_lambda_total;
    
    // Rock properties
    double* d_porosity;
    double* d_permeability;
    
    // PVT table data (device copies)
    double* d_pvt_pressure;
    double* d_pvt_Bo;
    double* d_pvt_Bg;
    double* d_pvt_Rs;
    int pvt_table_size;
#endif
};

// =============================================================================
// Thermal GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated thermal kernel
 * 
 * Implements heat transport with conduction and convection:
 *   (ρcp)_eff ∂T/∂t = ∇·(k_eff ∇T) + ρf cp_f v·∇T + Q
 * 
 * ## GPU Operations:
 * - Effective property computation
 * - Conduction flux (parallel over faces)
 * - Convection flux (parallel over faces)
 * - Source term integration
 * 
 * ## Performance:
 * - 2-5× speedup for conduction-dominated
 * - 4-10× speedup for convection-dominated
 */
class ThermalKernelGPU : public GPUKernelBase<ThermalKernelGPU, ThermalKernel> {
public:
    ThermalKernelGPU();
    ~ThermalKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeEffectivePropertiesGPU();
    void computeConductionFluxGPU();
    void computeConvectionFluxGPU();
    void applySourceTermsGPU();
    
private:
#ifdef USE_CUDA
    double* d_temperature;
    double* d_temperature_old;
    double* d_heat_flux;
    
    // Effective properties
    double* d_rho_cp_eff;
    double* d_k_eff;
    
    // Material properties
    double* d_k_solid;
    double* d_k_fluid;
    double* d_rho_solid;
    double* d_rho_fluid;
    double* d_cp_solid;
    double* d_cp_fluid;
    double* d_porosity;
    
    // Velocity field (for convection)
    double* d_velocity;
    
    // Heat sources
    double* d_source;
#endif
};

// =============================================================================
// Geomechanics GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated geomechanics kernel
 * 
 * Implements static/quasi-static elasticity with optional plasticity:
 *   ∇·σ + ρg = 0  (equilibrium)
 *   σ = C : ε - α p I  (constitutive, with pore pressure)
 * 
 * ## GPU Operations:
 * - Strain computation from displacement gradient
 * - Stress computation (elastic or elastoplastic)
 * - Pore pressure coupling (Biot)
 * - Residual force computation
 * 
 * ## Performance:
 * - 3-10× speedup for >100K cells
 * - Excellent GPU utilization due to local computations
 */
class GeomechanicsKernelGPU : public GPUKernelBase<GeomechanicsKernelGPU, GeomechanicsKernel> {
public:
    GeomechanicsKernelGPU();
    ~GeomechanicsKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeStrainGPU();
    void computeStressGPU();
    void computeForceGPU();
    void applyBiotCouplingGPU(const double* pressure);
    
    // Plasticity (Drucker-Prager or Mohr-Coulomb)
    void computePlasticCorrectionGPU();
    bool checkYieldGPU();
    
private:
#ifdef USE_CUDA
    double* d_displacement;
    double* d_strain;           // Voigt notation: εxx, εyy, εzz, γxy, γyz, γxz
    double* d_stress;           // Voigt notation: σxx, σyy, σzz, τxy, τyz, τxz
    double* d_force;
    
    // Material properties
    double* d_youngs_modulus;
    double* d_poisson_ratio;
    double* d_density;
    double* d_biot_coefficient;
    
    // Plasticity state
    double* d_plastic_strain;
    double* d_yield_state;       // 0 = elastic, 1 = yielded
    double* d_cohesion;
    double* d_friction_angle;
    
    // Pore pressure (coupling)
    double* d_pressure;
    
    // Mesh
    int* d_connectivity;
#endif
};

// =============================================================================
// Hydrodynamic GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated hydrodynamic kernel (Euler equations)
 * 
 * Implements compressible Euler equations for shock propagation:
 *   ∂ρ/∂t + ∇·(ρv) = 0
 *   ∂(ρv)/∂t + ∇·(ρv⊗v) + ∇p = ρg
 *   ∂E/∂t + ∇·((E+p)v) = ρv·g + Q
 * 
 * Supports equations of state:
 * - Ideal gas
 * - Mie-Grüneisen (condensed matter)
 * - Tillotson (impact/explosion)
 * - JWL (detonation products)
 * 
 * ## GPU Operations:
 * - Flux computation (HLLC or Roe)
 * - EOS evaluation
 * - Artificial viscosity
 * - Time integration (RK or MUSCL)
 * 
 * ## Performance:
 * - 10-50× speedup for shock problems
 * - Ideal for GPU due to local, stencil-based operations
 */
class HydrodynamicKernelGPU : public GPUKernelBase<HydrodynamicKernelGPU, HydrodynamicKernel> {
public:
    HydrodynamicKernelGPU();
    ~HydrodynamicKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeFluxGPU();              // HLLC Riemann solver
    void computeEOSGPU();               // Equation of state
    void applyArtificialViscosityGPU(); // Shock capturing
    void updateConservativeGPU(double dt);
    
    // Equation of state selection
    enum class EOSType { IDEAL_GAS, MIE_GRUNEISEN, TILLOTSON, JWL };
    void setEOS(EOSType type);
    
private:
#ifdef USE_CUDA
    // Conservative variables
    double* d_density;
    double* d_momentum_x;
    double* d_momentum_y;
    double* d_momentum_z;
    double* d_energy;
    
    // Primitive variables
    double* d_velocity_x;
    double* d_velocity_y;
    double* d_velocity_z;
    double* d_pressure;
    double* d_temperature;
    double* d_sound_speed;
    
    // Fluxes
    double* d_flux_density;
    double* d_flux_momentum_x;
    double* d_flux_momentum_y;
    double* d_flux_momentum_z;
    double* d_flux_energy;
    
    // Artificial viscosity
    double* d_q_viscosity;
    
    // EOS parameters
    EOSType eos_type;
    double* d_eos_params;
    
    // Mesh
    double* d_cell_volume;
    double* d_face_area;
    int* d_face_cells;
#endif
};

// =============================================================================
// Atmospheric Blast GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated atmospheric blast kernel
 * 
 * Implements atmospheric blast wave propagation with:
 * - Stratified atmosphere (density, temperature profiles)
 * - Wind effects
 * - Ground reflection and Mach stem
 * - Spherical/cylindrical divergence
 * 
 * ## GPU Operations:
 * - Compressible Euler with source terms
 * - Atmospheric profile interpolation
 * - Ground boundary condition
 * - Thermal precursor (optional)
 * 
 * ## Performance:
 * - 5-20× speedup for 2D/3D blast calculations
 */
class AtmosphericBlastKernelGPU : public GPUKernelBase<AtmosphericBlastKernelGPU, AtmosphericBlastKernel> {
public:
    AtmosphericBlastKernelGPU();
    ~AtmosphericBlastKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeAtmosphericStateGPU();
    void computeBlastFluxGPU();
    void applyGroundReflectionGPU();
    void computeSourceTermsGPU(double t);
    
private:
#ifdef USE_CUDA
    // State variables
    double* d_density;
    double* d_momentum;        // 3 components
    double* d_energy;
    
    // Atmospheric profile (device copy)
    double* d_atm_altitude;
    double* d_atm_density;
    double* d_atm_temperature;
    double* d_atm_pressure;
    double* d_wind_east;
    double* d_wind_north;
    int atm_profile_size;
    
    // Ground reflection
    double* d_ground_elevation;
    
    // Source
    double source_x, source_y, source_z;
    double yield_kt;
    double burst_time;
#endif
};

// =============================================================================
// Infrasound GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated infrasound kernel
 * 
 * Implements linearized acoustic wave equation for infrasound:
 *   ∂p'/∂t + ρc² ∇·v' + v₀·∇p' = S
 *   ρ ∂v'/∂t + ∇p' = 0
 * 
 * with atmospheric refraction and attenuation.
 * 
 * ## GPU Operations:
 * - Acoustic propagation (split-step or FD)
 * - Atmospheric profile interpolation
 * - Absorption coefficient computation
 * - Source injection
 * 
 * ## Performance:
 * - 5-30× speedup for PE method
 * - Ray tracing less amenable to GPU (adaptive)
 */
class InfrasoundKernelGPU : public GPUKernelBase<InfrasoundKernelGPU, InfrasoundKernel> {
public:
    InfrasoundKernelGPU();
    ~InfrasoundKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeAcousticFluxGPU();
    void applyAttenuationGPU(double frequency);
    void injectSourceGPU(double t);
    void computeEffectiveSoundSpeedGPU();
    
    // Parabolic equation on GPU
    void propagatePE_GPU(double frequency, double dr);
    
private:
#ifdef USE_CUDA
    // Acoustic field
    double* d_pressure;
    double* d_velocity_x;
    double* d_velocity_y;
    double* d_velocity_z;
    
    // Complex field for PE (real and imaginary)
    double* d_field_real;
    double* d_field_imag;
    
    // Atmospheric properties
    double* d_sound_speed;
    double* d_density;
    double* d_wind_u;
    double* d_wind_v;
    double* d_effective_c;    // c_eff = c + v·n
    
    // Absorption
    double* d_absorption_coeff;
    
    // Topography
    double* d_terrain_elevation;
    
    // Source
    double source_x, source_y, source_z;
    double source_strength;
    double source_frequency;
#endif
};

// =============================================================================
// Tsunami GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated tsunami (shallow water) kernel
 * 
 * Implements nonlinear shallow water equations:
 *   ∂η/∂t + ∇·((h+η)v) = 0
 *   ∂v/∂t + v·∇v + g∇η = -Cᵢ|v|v/(h+η) + f×v
 * 
 * where η is surface elevation, h is water depth, v is velocity.
 * 
 * ## GPU Operations:
 * - Flux computation (HLL or central)
 * - Wetting/drying treatment
 * - Coriolis force
 * - Bottom friction
 * 
 * ## Performance:
 * - 10-50× speedup for large ocean domains
 * - Excellent scaling due to local stencil operations
 */
class TsunamiKernelGPU : public GPUKernelBase<TsunamiKernelGPU, TsunamiKernel> {
public:
    TsunamiKernelGPU();
    ~TsunamiKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeFluxGPU();
    void applyBottomFrictionGPU();
    void applyCoriolisGPU();
    void handleWettingDryingGPU();
    
    // Seafloor deformation source
    void applySeafloorDeformationGPU(const double* dz, double t);
    
private:
#ifdef USE_CUDA
    // State variables
    double* d_eta;             // Surface elevation
    double* d_hu;              // x-momentum (depth-averaged)
    double* d_hv;              // y-momentum (depth-averaged)
    
    // Derived
    double* d_H;               // Total water depth (h + η)
    double* d_u;               // x-velocity
    double* d_v;               // y-velocity
    double* d_speed;           // |v|
    
    // Bathymetry
    double* d_depth;           // Still water depth h
    double* d_manning;         // Manning roughness coefficient
    
    // Coriolis
    double* d_coriolis;        // f = 2Ω sin(φ)
    
    // Fluxes
    double* d_flux_eta;
    double* d_flux_hu;
    double* d_flux_hv;
    
    // Wetting/drying
    double* d_wet_mask;        // 1 = wet, 0 = dry
    double min_depth;          // Minimum depth for wet cell
#endif
};

// =============================================================================
// Particle Transport GPU Kernel
// =============================================================================

/**
 * @brief GPU-accelerated particle transport kernel
 * 
 * Implements advection-diffusion with settling:
 *   φ ∂C/∂t + ∇·(vC) - ∇·(D∇C) + vs ∂C/∂z = 0
 * 
 * ## GPU Operations:
 * - Advection (upwind or TVD)
 * - Diffusion/dispersion
 * - Gravitational settling
 * - Deposition/resuspension
 */
class ParticleTransportKernelGPU : public GPUKernelBase<ParticleTransportKernelGPU, ParticleTransportKernel> {
public:
    ParticleTransportKernelGPU();
    ~ParticleTransportKernelGPU() override;
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    void allocateGPUMemory(int n_cells);
    void freeGPUMemory();
    
    // GPU batch operations
    void computeAdvectionFluxGPU();
    void computeDiffusionFluxGPU();
    void computeSettlingGPU();
    void applyDepositionGPU();
    
private:
#ifdef USE_CUDA
    double* d_concentration;
    double* d_concentration_old;
    double* d_velocity;        // 3 components
    double* d_diffusivity;
    double* d_settling_velocity;
    double* d_porosity;
    double* d_deposition_rate;
    
    // Dispersion tensor
    double* d_dispersivity_L;
    double* d_dispersivity_T;
#endif
};

// =============================================================================
// GPU Kernel Factory (Extended)
// =============================================================================

/**
 * @brief Extended factory for creating GPU kernels
 * 
 * Supports all implemented GPU kernels including extended physics.
 */
std::shared_ptr<PhysicsKernel> createGPUKernelExtended(PhysicsType type);

/**
 * @brief Check if extended GPU kernel is available
 */
bool hasExtendedGPUKernel(PhysicsType type);

/**
 * @brief Get list of all available GPU kernel types
 */
std::vector<PhysicsType> getAvailableGPUKernels();

/**
 * @brief Benchmark all available GPU kernels
 */
void benchmarkAllGPUKernels(int n_cells, int n_iterations);

} // namespace FSRM

#endif // PHYSICS_KERNEL_GPU_EXTENDED_HPP
