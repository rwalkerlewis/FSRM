/**
 * @file GPUPhysicsKernels.cuh
 * @brief Extended CUDA kernels for radiation, atmospheric, and explosion physics
 * 
 * Provides GPU-accelerated implementations for:
 * - Radiation transport and fallout
 * - Atmospheric dynamics
 * - Explosion and blast wave physics
 * - Fluid instabilities (RT, KH, RM)
 */

#ifndef GPU_PHYSICS_KERNELS_CUH
#define GPU_PHYSICS_KERNELS_CUH

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cuda.h>
#include <cufft.h>

namespace FSRM {
namespace GPU {
namespace Kernels {

// Block sizes
constexpr int BLOCK_1D = 256;
constexpr int BLOCK_2D = 16;
constexpr int BLOCK_3D = 8;

// =============================================================================
// Atmospheric Model Kernels
// =============================================================================

/**
 * @brief Compute US Standard Atmosphere 1976 profiles
 */
__global__ void computeUSStandard1976(
    const double* altitude,       // Input: geometric altitude (m)
    int n,
    double* temperature,          // Output: K
    double* pressure,             // Output: Pa
    double* density               // Output: kg/m³
);

/**
 * @brief Compute NRLMSISE-00 model (simplified)
 */
__global__ void computeNRLMSISE(
    const double* altitude,       // Input: altitude (m)
    int n,
    double f107,                  // Solar flux
    double ap,                    // Geomagnetic index
    int day_of_year,
    double latitude,
    double* temperature,
    double* density,
    double* n2_density,           // Molecular nitrogen
    double* o2_density,           // Molecular oxygen
    double* o_density             // Atomic oxygen
);

/**
 * @brief Interpolate wind profile
 */
__global__ void interpolateWindProfile(
    const double* query_alt,      // Query altitudes
    int n_query,
    const double* profile_alt,    // Profile altitudes
    const double* wind_u,         // East wind component
    const double* wind_v,         // North wind component
    int n_profile,
    double* wind_u_out,
    double* wind_v_out
);

/**
 * @brief Compute sound speed with humidity correction
 */
__global__ void computeSoundSpeedHumid(
    const double* temperature,
    const double* relative_humidity,
    const double* pressure,
    int n,
    double* sound_speed
);

/**
 * @brief Compute atmospheric refraction for acoustic/infrasound
 */
__global__ void computeAcousticRefraction(
    const double* sound_speed,    // Sound speed profile
    const double* wind_u,         // Wind components
    const double* wind_v,
    const double* ray_dx,         // Ray direction
    const double* ray_dy,
    const double* ray_dz,
    int n,
    double* effective_c,          // Effective sound speed
    double* refraction_angle      // Ray bending
);

/**
 * @brief Apply atmospheric absorption to acoustic signal
 */
__global__ void applyAtmosphericAbsorption(
    double* amplitude,            // Signal amplitude
    const double* distance,       // Propagation distance
    const double* frequency,      // Frequency array
    const double* temperature,
    const double* pressure,
    const double* humidity,
    int n_points,
    int n_freq
);

// =============================================================================
// Radiation Physics Kernels
// =============================================================================

/**
 * @brief Compute Planck blackbody spectral radiance
 */
__global__ void planckRadiance(
    const double* temperature,    // K
    const double* wavelength,     // m
    int n_temp,
    int n_wave,
    double* radiance              // W/(m²·sr·m)
);

/**
 * @brief Stefan-Boltzmann total radiance
 */
__global__ void stefanBoltzmannRadiance(
    const double* temperature,
    int n,
    double emissivity,
    double* radiance              // W/m²
);

/**
 * @brief Compute absorption in participating medium
 */
__global__ void computeAbsorption(
    const double* incident_intensity,
    const double* absorption_coeff,
    const double* path_length,
    int n,
    double* transmitted_intensity
);

/**
 * @brief Radiative heat transfer view factor calculation
 */
__global__ void computeViewFactors(
    const double* positions,      // Surface center positions [n, 3]
    const double* normals,        // Surface normals [n, 3]
    const double* areas,          // Surface areas [n]
    int n_surfaces,
    double* view_factors          // Output [n, n]
);

/**
 * @brief Discrete ordinates radiative transfer step
 */
__global__ void discreteOrdinatesStep(
    double* intensity,            // Radiation intensity [nx, ny, n_angles]
    const double* source,         // Source term
    const double* sigma_a,        // Absorption coefficient
    const double* sigma_s,        // Scattering coefficient
    const double* weights,        // Quadrature weights [n_angles]
    const double* mu,             // Direction cosines [n_angles, 3]
    double dx, double dy,
    int nx, int ny, int n_angles
);

// =============================================================================
// Nuclear/Radiation Fallout Kernels
// =============================================================================

/**
 * @brief Compute radioactive decay for multiple species
 */
__global__ void radioactiveDecay(
    double* activity,             // Current activity [n_species]
    const double* half_life,      // Half-life in seconds
    double dt,
    int n_species
);

/**
 * @brief Compute decay chain (parent -> daughter)
 */
__global__ void decayChain(
    double* parent_activity,
    double* daughter_activity,
    const double* parent_half_life,
    const double* daughter_half_life,
    const double* branching_ratio,
    double dt,
    int n_chains
);

/**
 * @brief Compute gamma dose rate from point sources
 */
__global__ void gammaDoseRate(
    const double* source_positions,   // [n_sources, 3]
    const double* source_activity,    // Bq
    const double* dose_conversion,    // Sv·m²/(Bq·s)
    const double* observer_positions, // [n_observers, 3]
    int n_sources,
    int n_observers,
    double* dose_rate                 // Sv/s
);

/**
 * @brief Compute fallout particle settling velocity
 */
__global__ void settlingVelocity(
    const double* particle_diameter,  // m
    const double* particle_density,   // kg/m³
    const double* air_density,        // kg/m³
    const double* air_viscosity,      // Pa·s
    int n_particles,
    double* v_settle                  // m/s
);

/**
 * @brief Advect fallout particles with wind
 */
__global__ void advectFalloutParticles(
    double* pos_x, double* pos_y, double* pos_z,
    const double* v_settle,
    const double* wind_u,             // Wind interpolated at particle positions
    const double* wind_v,
    const double* wind_w,
    const double* diffusivity,
    int n_particles,
    double dt,
    unsigned long seed                // For random diffusion
);

/**
 * @brief Deposit fallout on ground
 */
__global__ void depositFallout(
    const double* pos_x, const double* pos_y, const double* pos_z,
    const double* activity,
    const int* deposited_flag,        // 1 if particle reached ground
    double* ground_deposition,        // [nx, ny] ground activity
    double dx, double dy,
    double x_origin, double y_origin,
    int nx, int ny,
    int n_particles
);

/**
 * @brief Compute total dose from time-integrated exposure
 */
__global__ void integratedDose(
    const double* dose_rate,          // Time series [n_times]
    const double* times,
    int n_times,
    double shielding_factor,          // 0-1
    double* total_dose                // Output
);

// =============================================================================
// Explosion/Blast Wave Kernels
// =============================================================================

/**
 * @brief Sedov-Taylor blast wave solution
 */
__global__ void sedovTaylorBlast(
    double energy,                // Total energy (J)
    double rho0,                  // Ambient density (kg/m³)
    double gamma,                 // Ratio of specific heats
    double time,                  // Time since detonation
    const double* r,              // Radial distances
    int n,
    double* pressure,
    double* density,
    double* velocity
);

/**
 * @brief Nuclear fireball evolution
 */
__global__ void nuclearFireball(
    double yield_kt,              // Yield in kilotons
    double altitude,              // Burst altitude (m)
    double time,                  // Time since detonation
    double* radius,               // Fireball radius
    double* temperature,          // Surface temperature
    double* rise_velocity,        // Buoyant rise velocity
    double* luminosity            // Radiated power
);

/**
 * @brief Thermal radiation from fireball
 */
__global__ void fireballThermalFlux(
    double fireball_radius,
    double fireball_temperature,
    const double* observer_positions,  // [n, 3]
    const double* fireball_position,   // [3]
    int n_observers,
    double* thermal_flux               // W/m²
);

/**
 * @brief Blast overpressure at distance
 */
__global__ void blastOverpressure(
    double yield_kt,
    double altitude,
    const double* observer_positions,  // [n, 3]
    const double* burst_position,      // [3]
    int n_observers,
    double* peak_overpressure,         // Pa
    double* arrival_time,              // s
    double* positive_duration          // s
);

/**
 * @brief Ground reflection and Mach stem
 */
__global__ void machStem(
    double incident_overpressure,
    double incident_angle,        // From vertical
    double ground_reflection_coeff,
    double* reflected_overpressure,
    double* mach_stem_height,
    double* mach_overpressure
);

/**
 * @brief Mushroom cloud dynamics step
 */
__global__ void mushroomCloudStep(
    double* cloud_density,        // [nx, ny, nz]
    double* cloud_temperature,
    double* velocity_u,
    double* velocity_v,
    double* velocity_w,
    const double* ambient_density,
    const double* ambient_temperature,
    double g,
    double entrainment_coeff,
    double dx, double dy, double dz,
    double dt,
    int nx, int ny, int nz
);

// =============================================================================
// Fluid Instability Kernels
// =============================================================================

/**
 * @brief Rayleigh-Taylor instability growth
 */
__global__ void rayleighTaylorGrowth(
    double* interface_height,     // z position of interface [nx, ny]
    double atwood_number,         // (ρ2 - ρ1)/(ρ2 + ρ1)
    double g,                     // Acceleration
    double dt,
    int nx, int ny
);

/**
 * @brief Rayleigh-Taylor mixing zone width
 */
__global__ void rtMixingWidth(
    const double* density,        // [nx, ny, nz]
    double rho_light,
    double rho_heavy,
    int nx, int ny, int nz,
    double* mixing_width          // Output
);

/**
 * @brief Kelvin-Helmholtz instability
 */
__global__ void kelvinHelmholtzStep(
    double* vorticity,            // [nx, ny]
    const double* velocity_u,
    const double* velocity_v,
    double density_ratio,
    double dt,
    double dx, double dy,
    int nx, int ny
);

/**
 * @brief Richtmyer-Meshkov instability (shock-interface)
 */
__global__ void richtmyerMeshkovGrowth(
    double* interface_amplitude,  // Mode amplitudes
    double atwood_number,
    double shock_mach,
    double delta_v,               // Velocity jump from shock
    double time,
    int n_modes
);

/**
 * @brief Buoyancy force computation
 */
__global__ void computeBuoyancy(
    const double* density,
    double reference_density,
    double g,
    int nx, int ny, int nz,
    double* buoyancy_z            // Vertical buoyancy force
);

// =============================================================================
// Infrasound Propagation Kernels
// =============================================================================

/**
 * @brief Parabolic equation (PE) propagation step
 */
__global__ void parabolicEquationStep(
    double* field_real,           // Real part of pressure field
    double* field_imag,           // Imaginary part
    const double* n_squared,      // Refractive index squared
    double k0,                    // Reference wavenumber
    double dr,                    // Range step
    double dz,                    // Vertical step
    int nr, int nz
);

/**
 * @brief Wide-angle PE propagation
 */
__global__ void wideAnglePE(
    double* field_real,
    double* field_imag,
    const double* n_squared,
    const double* absorption,
    double k0,
    double dr, double dz,
    int nr, int nz
);

/**
 * @brief Apply ground boundary condition
 */
__global__ void groundBoundaryPE(
    double* field_real,
    double* field_imag,
    double ground_impedance_real,
    double ground_impedance_imag,
    int nz
);

/**
 * @brief Compute transmission loss
 */
__global__ void transmissionLoss(
    const double* field_real,
    const double* field_imag,
    double source_level,
    int n,
    double* TL                    // dB re 1 m
);

// =============================================================================
// Helper Device Functions
// =============================================================================

/**
 * @brief Get global thread ID for 1D kernel
 */
__device__ __forceinline__ int getGlobalId1D() {
    return blockIdx.x * blockDim.x + threadIdx.x;
}

/**
 * @brief Get global thread ID for 2D kernel
 */
__device__ __forceinline__ void getGlobalId2D(int& i, int& j) {
    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
}

/**
 * @brief Get global thread ID for 3D kernel
 */
__device__ __forceinline__ void getGlobalId3D(int& i, int& j, int& k) {
    i = blockIdx.x * blockDim.x + threadIdx.x;
    j = blockIdx.y * blockDim.y + threadIdx.y;
    k = blockIdx.z * blockDim.z + threadIdx.z;
}

/**
 * @brief Linear interpolation
 */
__device__ __forceinline__ double lerp(double a, double b, double t) {
    return a + t * (b - a);
}

/**
 * @brief Bilinear interpolation
 */
__device__ double bilerp(
    double v00, double v01, double v10, double v11,
    double tx, double ty
);

/**
 * @brief Trilinear interpolation
 */
__device__ double trilerp(
    double v000, double v001, double v010, double v011,
    double v100, double v101, double v110, double v111,
    double tx, double ty, double tz
);

/**
 * @brief Complex multiplication
 */
__device__ __forceinline__ void cmul(
    double ar, double ai, double br, double bi,
    double& cr, double& ci
) {
    cr = ar * br - ai * bi;
    ci = ar * bi + ai * br;
}

/**
 * @brief Atomic add for double precision
 */
__device__ double atomicAddDouble(double* address, double val);

/**
 * @brief Compute Euclidean distance
 */
__device__ __forceinline__ double distance3D(
    double x1, double y1, double z1,
    double x2, double y2, double z2
) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

// =============================================================================
// Kernel Launch Helpers
// =============================================================================

/**
 * @brief Get optimal block size for 1D kernel
 */
inline dim3 getBlocks1D(int n, int block_size = BLOCK_1D) {
    return dim3((n + block_size - 1) / block_size);
}

/**
 * @brief Get optimal grid for 2D kernel
 */
inline dim3 getBlocks2D(int nx, int ny, int block_size = BLOCK_2D) {
    return dim3((nx + block_size - 1) / block_size,
                (ny + block_size - 1) / block_size);
}

/**
 * @brief Get optimal grid for 3D kernel
 */
inline dim3 getBlocks3D(int nx, int ny, int nz, int block_size = BLOCK_3D) {
    return dim3((nx + block_size - 1) / block_size,
                (ny + block_size - 1) / block_size,
                (nz + block_size - 1) / block_size);
}

} // namespace Kernels
} // namespace GPU
} // namespace FSRM

#endif // USE_CUDA

#endif // GPU_PHYSICS_KERNELS_CUH
